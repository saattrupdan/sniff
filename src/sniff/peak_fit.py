"""Measured-shape peak fitting for overlapping PTR-MS channels.

The functions in this module are deliberately independent of HDF5. They learn a
normalised profile from an average spectrum, fit a small overlapping group, and return
a projection that can be applied efficiently to every cycle by :mod:`sniff.ptrms`.
"""

from __future__ import annotations

import numpy as np

_PROFILE_X = np.linspace(-4.5, 4.5, 181)


def estimate_empirical_profile(spectrum, centres, sigmas):
    """Estimate a robust normalised line shape from isolated target peaks.

    Args:
        spectrum:
            Average spectrum indexed by timebin.
        centres:
            Timebin centres for isolated candidate peaks.
        sigmas:
            Corresponding physical Gaussian sigma estimates in timebins.

    Returns:
        A dictionary containing the profile and diagnostics. ``usable`` is false when
        fewer than three clean peaks support the estimate.
    """
    spectrum = np.asarray(spectrum, dtype=np.float64)
    profiles = []
    for centre, sigma in zip(centres, sigmas):
        if not np.isfinite(centre) or not np.isfinite(sigma) or sigma <= 0:
            continue
        x = centre + _PROFILE_X * sigma
        lo = max(0, int(np.floor(x[0])))
        hi = min(len(spectrum), int(np.ceil(x[-1])) + 1)
        if hi - lo < 7:
            continue
        bins = np.arange(lo, hi, dtype=np.float64)
        values = np.interp(x, bins, spectrum[lo:hi], left=np.nan, right=np.nan)
        if not np.all(np.isfinite(values)):
            continue
        edge = np.concatenate((values[:20], values[-20:]))
        baseline = float(np.median(edge))
        signal = values - baseline
        height = float(np.max(signal))
        noise = float(np.median(np.abs(edge - baseline))) * 1.4826
        if height <= 0 or height < max(8.0 * noise, 1.0):
            continue
        signal = np.clip(signal / height, 0.0, None)
        if signal[len(signal) // 2] < 0.35:
            continue
        profiles.append(signal)

    if len(profiles) < 3:
        return {
            "usable": False,
            "reason": "fewer than three clean isolated peaks",
            "n_reference_peaks": len(profiles),
            "x": _PROFILE_X,
            "y": np.exp(-0.5 * _PROFILE_X**2),
        }

    profile = np.median(np.asarray(profiles), axis=0)
    maximum = float(np.max(profile))
    if maximum <= 0:
        return {
            "usable": False,
            "reason": "the empirical profile has no positive signal",
            "n_reference_peaks": len(profiles),
            "x": _PROFILE_X,
            "y": np.exp(-0.5 * _PROFILE_X**2),
        }
    profile = np.clip(profile / maximum, 0.0, None)
    return {
        "usable": True,
        "reason": None,
        "n_reference_peaks": len(profiles),
        "x": _PROFILE_X,
        "y": profile,
    }


def solve_nonnegative(design, values):
    """Solve a small non-negative least-squares problem without SciPy.

    Negative coefficients are removed from the active set and the remaining columns
    are refitted. Unlike clipping an unconstrained solution, this recomputes the other
    amplitudes for the small, mildly correlated peak groups used here.
    """
    design = np.asarray(design, dtype=np.float64)
    values = np.asarray(values, dtype=np.float64)
    active = list(range(design.shape[1]))
    result = np.zeros(design.shape[1], dtype=np.float64)
    while active:
        fitted, _, _, _ = np.linalg.lstsq(design[:, active], values, rcond=None)
        negative = [active[i] for i, value in enumerate(fitted) if value < 0]
        if not negative:
            result[active] = fitted
            break
        worst = min(negative, key=lambda index: fitted[active.index(index)])
        active.remove(worst)
    return result


def fit_group_design(
    spectrum,
    centres,
    sigmas,
    profile,
    *,
    max_shift_sigma=0.6,
    width_limits=(0.75, 1.35),
):
    """Fit a shared shift and width for one overlapping peak group.

    The expensive shape search uses the run or interval average spectrum. The returned
    projection then solves amplitudes for all cycles with one matrix multiplication.
    """
    spectrum = np.asarray(spectrum, dtype=np.float64)
    centres = np.asarray(centres, dtype=np.float64)
    sigmas = np.asarray(sigmas, dtype=np.float64)
    margin = 6.0 * float(np.max(sigmas))
    lo = max(0, int(np.floor(np.min(centres) - margin)))
    hi = min(len(spectrum), int(np.ceil(np.max(centres) + margin)) + 1)
    x = np.arange(lo, hi, dtype=np.float64)
    values = spectrum[lo:hi]
    if len(x) <= len(centres) + 1:
        return _failed_fit(lo, hi, len(centres), "too few bins for the fitted group")

    profile_x = np.asarray(profile["x"], dtype=np.float64)
    profile_y = np.asarray(profile["y"], dtype=np.float64)
    scale_sigma = float(np.median(sigmas))
    shifts = np.linspace(
        -max_shift_sigma * scale_sigma,
        max_shift_sigma * scale_sigma,
        9,
    )
    widths = np.linspace(float(width_limits[0]), float(width_limits[1]), 9)
    best = None
    for shift in shifts:
        for width in widths:
            components = _components(
                x, centres + shift, sigmas * width, profile_x, profile_y
            )
            design = np.column_stack((components, np.ones(len(x))))
            if (
                not np.isfinite(design).all()
                or float(np.max(np.abs(design))) > 1e6
            ):
                continue
            singular = np.linalg.svd(design, compute_uv=False)
            if (
                np.count_nonzero(singular > singular[0] * 1e-10)
                != design.shape[1]
                or singular[-1] <= 0
                or singular[0] / singular[-1] > 1e8
            ):
                continue
            value_scale = max(float(np.max(np.abs(values))), 1.0)
            scaled_values = values / value_scale
            with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
                coeff = solve_nonnegative(design, scaled_values)
            if (
                not np.isfinite(coeff).all()
                or float(np.max(np.abs(coeff))) > 1e12
            ):
                continue
            with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
                residual = scaled_values - design @ coeff
                score = float(np.dot(residual, residual))
            if np.isfinite(score) and (best is None or score < best[0]):
                best = (score, shift, width, components, design, coeff, residual)

    if best is None:
        return _failed_fit(lo, hi, len(centres), "group design is rank-deficient")
    _, shift, width, components, design, coeff, residual = best
    component_design = design[:, :-1]
    singular = np.linalg.svd(component_design, compute_uv=False)
    rank = int(np.linalg.matrix_rank(component_design))
    condition = (
        float(singular[0] / singular[-1])
        if singular.size and singular[-1] > 0
        else float("inf")
    )
    scaled_values = values / max(float(np.max(np.abs(values))), 1.0)
    scale = float(np.linalg.norm(scaled_values - coeff[-1]))
    relative_residual = float(np.linalg.norm(residual) / scale) if scale > 0 else 0.0
    correlation = _max_correlation(component_design)
    reliable = (
        rank == len(centres)
        and condition <= 100.0
        and correlation <= 0.995
        and relative_residual <= 0.35
    )
    if reliable:
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            projection = np.linalg.pinv(component_design)
        if not np.isfinite(projection).all():
            reliable = False
    if not reliable:
        projection = np.full(
            (len(centres), component_design.shape[0]),
            np.nan,
            dtype=np.float64,
        )
    status = "reliable" if reliable else "unresolved"

    # Baseline is estimated per cycle from edge bins before this projection. Unreliable
    # designs never receive a usable pseudoinverse.
    return {
        "usable": True,
        "lo": lo,
        "hi": hi,
        "components": components,
        "projection": projection,
        "shift_tb": float(shift),
        "width_scale": float(width),
        "rank": rank,
        "condition": condition,
        "correlation": correlation,
        "relative_residual": relative_residual,
        "baseline": float(coeff[-1]),
        "status": status,
        "reason": None if reliable else "components are not independently identifiable",
    }


def apply_group_design(chunk, fit):
    """Apply a fitted group design to a block of spectra."""
    values = np.asarray(chunk[:, fit["lo"] : fit["hi"]], dtype=np.float64)
    edge_count = max(1, min(8, values.shape[1] // 4))
    edges = np.concatenate((values[:, :edge_count], values[:, -edge_count:]), axis=1)
    baseline = np.median(edges, axis=1)
    centred = np.clip(values - baseline[:, None], 0.0, None)
    scale = np.maximum(np.max(centred, axis=1), 1.0)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        unconstrained = (centred / scale[:, None]) @ fit["projection"].T
        unconstrained *= scale[:, None]
    unconstrained[~np.isfinite(unconstrained)] = np.nan
    if np.all(unconstrained >= 0):
        return unconstrained
    result = np.empty_like(unconstrained)
    for row, values_row in enumerate(centred):
        row_scale = max(float(np.max(values_row)), 1.0)
        solved = solve_nonnegative(fit["components"], values_row / row_scale)
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            result[row] = solved * row_scale
    result[~np.isfinite(result)] = np.nan
    return result


def _components(x, centres, sigmas, profile_x, profile_y):
    columns = []
    for centre, sigma in zip(centres, sigmas):
        normalised = (x - centre) / sigma
        columns.append(np.interp(normalised, profile_x, profile_y, left=0.0, right=0.0))
    return np.column_stack(columns)


def _max_correlation(design):
    if design.shape[1] < 2:
        return 0.0
    normalised = design / np.maximum(np.linalg.norm(design, axis=0), 1e-15)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        correlation = normalised.T @ normalised
    if not np.isfinite(correlation).all():
        return float("inf")
    correlation[np.diag_indices_from(correlation)] = 0.0
    return float(np.max(np.abs(correlation)))


def _failed_fit(lo, hi, count, reason):
    return {
        "usable": False,
        "lo": lo,
        "hi": hi,
        "components": np.empty((max(0, hi - lo), count)),
        "projection": None,
        "shift_tb": 0.0,
        "width_scale": 1.0,
        "rank": 0,
        "condition": float("inf"),
        "correlation": 1.0,
        "relative_residual": float("inf"),
        "baseline": None,
        "status": "fallback",
        "reason": reason,
    }
