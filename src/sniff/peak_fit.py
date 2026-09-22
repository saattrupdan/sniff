"""Measured-shape peak fitting for overlapping PTR-MS channels.

The functions in this module are deliberately independent of HDF5. They learn a
normalised profile from an average spectrum, fit a small overlapping group, and return
a projection that can be applied efficiently to every cycle by :mod:`sniff.ptrms`.
"""

from __future__ import annotations

import numpy as np

_PROFILE_X = np.linspace(-4.5, 4.5, 181)


class PeakFitCancelled(RuntimeError):
    """Raised when a caller cancels a whole-run overlap fit."""


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
            if not np.isfinite(design).all() or float(np.max(np.abs(design))) > 1e6:
                continue
            singular = np.linalg.svd(design, compute_uv=False)
            if (
                np.count_nonzero(singular > singular[0] * 1e-10) != design.shape[1]
                or singular[-1] <= 0
                or singular[0] / singular[-1] > 1e8
            ):
                continue
            value_scale = max(float(np.max(np.abs(values))), 1.0)
            scaled_values = values / value_scale
            with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
                coeff = solve_nonnegative(design, scaled_values)
            if not np.isfinite(coeff).all() or float(np.max(np.abs(coeff))) > 1e12:
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


def fit_global_group(
    spectra,
    fit,
    *,
    breaks=(),
    lambda_grid=(0.001, 0.01, 0.1, 1.0),
    shape=None,
    should_stop=None,
    progress=None,
):
    """Fit non-negative, temporally regularised traces for one whole-run group.

    Shared centres and widths come from :func:`fit_group_design`. Every validation fold
    refits those shape parameters without its cycle block or hidden spectral bins. A
    complete decomposition must also outperform every independently refitted model with
    one component removed.
    """
    _check_stop(should_stop)
    spectra = np.asarray(spectra)
    components = np.asarray(fit["components"], dtype=np.float64)
    if spectra.ndim != 2 or spectra.shape[1] != components.shape[0]:
        return _failed_global(fit, "whole-run spectral window has incompatible shape")
    if fit.get("status") != "reliable" or not np.isfinite(spectra).all():
        return _failed_global(
            fit,
            fit.get("reason") or "shared peak shape is not independently identifiable",
        )
    centred, noise = _baseline_correct(spectra)
    design_scale = max(float(np.linalg.norm(components, ord=2) ** 2), 1e-12)
    lambda_multipliers = [float(value) for value in lambda_grid if value > 0]
    if not lambda_multipliers:
        return _failed_global(fit, "temporal regularisation grid is empty")

    folds = _validation_folds(spectra, components, breaks, shape, should_stop)
    _report_progress(progress, 0.2)
    cv = _cross_validated_errors(folds, lambda_multipliers, should_stop)
    if not cv:
        return _failed_global(fit, "held-out spectral prediction is unavailable")
    minimum = min(item["mean_error"] for item in cv)
    best = min(cv, key=lambda item: item["mean_error"])
    eligible = [
        item for item in cv if item["mean_error"] <= minimum + best["standard_error"]
    ]
    selected = max(eligible, key=lambda item: item["lambda_multiplier"])
    _report_progress(progress, 0.4)
    selected_lambda = selected["lambda_multiplier"] * design_scale
    amplitudes = _solve_temporal_nonnegative(
        centred, components, selected_lambda, breaks, should_stop
    )

    reduced = []
    if components.shape[1] > 1:
        for removed in range(components.shape[1]):
            _check_stop(should_stop)
            indexes = [
                index for index in range(components.shape[1]) if index != removed
            ]
            reduced_folds = _validation_folds(
                spectra,
                components[:, indexes],
                breaks,
                shape,
                should_stop,
                component_indexes=indexes,
            )
            reduced_cv = _cross_validated_errors(
                reduced_folds, lambda_multipliers, should_stop
            )
            if reduced_cv:
                reduced.append(min(item["mean_error"] for item in reduced_cv))
            _report_progress(progress, 0.4 + 0.4 * (removed + 1) / components.shape[1])
    best_reduced = min(reduced) if reduced else float("inf")
    improvement = (
        (best_reduced - selected["mean_error"]) / best_reduced
        if np.isfinite(best_reduced) and best_reduced > 0
        else 0.0
    )
    threshold = max(
        3.0 * max(float(np.median(noise)), 1e-12),
        0.01 * max(float(np.nanmax(amplitudes)), 1e-12),
    )
    active_counts = np.count_nonzero(amplitudes > threshold, axis=0)
    minimum_active = max(10, int(np.ceil(0.01 * spectra.shape[0])))
    shifts = np.asarray([fold["shift"] for fold in folds], dtype=np.float64)
    widths = np.asarray([fold["width"] for fold in folds], dtype=np.float64)
    shift_span = float(np.ptp(shifts)) if shifts.size else float("inf")
    width_span = float(np.ptp(widths)) if widths.size else float("inf")
    shift_limit = 0.2 * float(np.median(shape["sigmas"])) if shape else float("inf")
    failed_gates = []
    if selected["mean_error"] > 0.35:
        failed_gates.append("held-out relative RMSE exceeds 0.35")
    if improvement < 0.05:
        failed_gates.append("full model does not improve held-out error by 5%")
    if shape and shift_span > shift_limit:
        failed_gates.append("held-out centre shift varies by more than 0.2 sigma")
    if shape and width_span > 0.15:
        failed_gates.append("held-out width scale varies by more than 0.15")
    if np.any(active_counts < minimum_active):
        failed_gates.append("one or more components lack sufficient active cycles")
    if not np.isfinite(amplitudes).all():
        failed_gates.append("global amplitudes contain non-finite values")
    reliable = not failed_gates
    result = dict(fit)
    result.update(
        {
            "model": "joint-temporal-v2",
            "status": "reliable" if reliable else "unresolved",
            "reason": None if reliable else "; ".join(failed_gates),
            "amplitudes": (
                amplitudes
                if reliable
                else np.full_like(amplitudes, np.nan, dtype=np.float64)
            ),
            "selected_lambda": selected_lambda,
            "selected_lambda_multiplier": selected["lambda_multiplier"],
            "lambda_grid": [value * design_scale for value in lambda_multipliers],
            "held_out_relative_rmse": selected["mean_error"],
            "held_out_standard_error": selected["standard_error"],
            "best_reduced_relative_rmse": (
                best_reduced if np.isfinite(best_reduced) else None
            ),
            "held_out_improvement": improvement,
            "held_out_shift_span": shift_span if shape else None,
            "held_out_width_span": width_span if shape else None,
            "active_cycle_counts": active_counts.tolist(),
            "minimum_active_cycles": minimum_active,
            "failed_gates": failed_gates,
        }
    )
    _report_progress(progress, 1.0)
    return result


def _report_progress(progress, value):
    if progress is not None:
        progress(max(0.0, min(1.0, float(value))))


def _baseline_correct(spectra):
    edge_count = max(1, min(8, spectra.shape[1] // 4))
    edges = np.concatenate((spectra[:, :edge_count], spectra[:, -edge_count:]), axis=1)
    baseline = np.median(edges, axis=1)
    centred = np.clip(spectra - baseline[:, None], 0.0, None)
    deviation = np.median(np.abs(edges - baseline[:, None]), axis=1)
    return centred, 1.4826 * deviation


def _validation_folds(
    values,
    components,
    breaks,
    shape,
    should_stop,
    *,
    component_indexes=None,
):
    n_cycles, n_bins = values.shape
    if n_cycles < 20 or n_bins < 4:
        return []
    blocks = [block for block in np.array_split(np.arange(n_cycles), 8) if len(block)]
    slices = [slice(int(block[0]), int(block[-1]) + 1) for block in blocks]
    bins = np.arange(n_bins)
    total = np.sum(values, axis=0, dtype=np.float64)
    edge_count = max(1, min(8, n_bins // 4))
    edge_bins = np.zeros(n_bins, dtype=bool)
    edge_bins[:edge_count] = True
    edge_bins[-edge_count:] = True
    folds = []
    for block in slices:
        block_count = block.stop - block.start
        raw_training_mean = (
            total - np.sum(values[block], axis=0, dtype=np.float64)
        ) / (n_cycles - block_count)
        for parity in (0, 1):
            _check_stop(should_stop)
            observed = bins % 2 == parity
            hidden = ~observed
            if np.count_nonzero(observed) <= components.shape[1] or not np.any(hidden):
                continue
            baseline_bins = observed & edge_bins
            if not np.any(baseline_bins):
                continue
            baselines = np.median(values[:, baseline_bins], axis=1)
            training_baseline = (
                float(np.sum(baselines)) - float(np.sum(baselines[block]))
            ) / (n_cycles - block_count)
            training_mean = np.clip(raw_training_mean - training_baseline, 0.0, None)
            fold_values = np.clip(values[block] - baselines[block, None], 0.0, None)
            shift = 0.0
            width = 1.0
            fold_components = components
            if shape is not None:
                fitted = _fit_masked_shape(
                    training_mean, observed, shape, component_indexes
                )
                if fitted is None:
                    continue
                fold_components, shift, width = fitted
            local_breaks = [
                int(boundary - block.start)
                for boundary in breaks
                if block.start < boundary < block.stop
            ]
            folds.append(
                {
                    "values": fold_values,
                    "observed": observed,
                    "hidden": hidden,
                    "components": fold_components,
                    "breaks": local_breaks,
                    "shift": shift,
                    "width": width,
                }
            )
    return folds


def _fit_masked_shape(values, observed, shape, component_indexes):
    centres = np.asarray(shape["centres"], dtype=np.float64)
    sigmas = np.asarray(shape["sigmas"], dtype=np.float64)
    if component_indexes is not None:
        centres = centres[component_indexes]
        sigmas = sigmas[component_indexes]
    x = np.asarray(shape["x"], dtype=np.float64)
    profile_x = np.asarray(shape["profile_x"], dtype=np.float64)
    profile_y = np.asarray(shape["profile_y"], dtype=np.float64)
    scale_sigma = float(np.median(sigmas))
    shifts = np.linspace(-0.6 * scale_sigma, 0.6 * scale_sigma, 9)
    widths = np.linspace(0.75, 1.35, 9)
    value_scale = max(float(np.max(np.abs(values[observed]))), 1.0)
    target = values[observed] / value_scale
    best = None
    for shift in shifts:
        for width in widths:
            components = _components(
                x, centres + shift, sigmas * width, profile_x, profile_y
            )
            design = np.column_stack(
                (components[observed], np.ones(np.count_nonzero(observed)))
            )
            if np.linalg.matrix_rank(design) != design.shape[1]:
                continue
            with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
                coefficients = solve_nonnegative(design, target)
                residual = target - design @ coefficients
                score = float(np.dot(residual, residual))
            if np.isfinite(score) and (best is None or score < best[0]):
                best = (score, components, float(shift), float(width))
    return None if best is None else best[1:]


def _cross_validated_errors(folds, lambda_multipliers, should_stop):
    results = []
    for multiplier in lambda_multipliers:
        fold_errors = []
        for fold in folds:
            _check_stop(should_stop)
            observed = fold["observed"]
            hidden = fold["hidden"]
            components = fold["components"]
            fold_scale = max(
                float(np.linalg.norm(components[observed], ord=2) ** 2), 1e-12
            )
            amplitudes = _solve_temporal_nonnegative(
                fold["values"][:, observed],
                components[observed],
                multiplier * fold_scale,
                fold["breaks"],
                should_stop,
            )
            with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
                predicted = amplitudes @ components[hidden].T
            actual = fold["values"][:, hidden]
            scale = max(float(np.linalg.norm(actual)), 1e-12)
            error = float(np.linalg.norm(actual - predicted) / scale)
            if np.isfinite(error):
                fold_errors.append(error)
        if fold_errors:
            spread = float(np.std(fold_errors, ddof=1)) if len(fold_errors) > 1 else 0.0
            results.append(
                {
                    "lambda_multiplier": multiplier,
                    "mean_error": float(np.mean(fold_errors)),
                    "standard_error": spread / np.sqrt(len(fold_errors)),
                    "fold_count": len(fold_errors),
                }
            )
    return results


def _solve_temporal_nonnegative(
    values, components, smoothing, breaks, should_stop=None
):
    values = np.asarray(values, dtype=np.float64)
    components = np.asarray(components, dtype=np.float64)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        gram = components.T @ components
        cross = values @ components
        initial = values @ np.linalg.pinv(components).T
    amplitudes = np.clip(initial, 0.0, None)
    edge_weights = np.ones(max(0, values.shape[0] - 1), dtype=np.float64)
    for boundary in breaks:
        boundary = int(boundary)
        if 0 < boundary < values.shape[0]:
            edge_weights[boundary - 1] = 0.0
    lipschitz = max(float(np.linalg.norm(gram, ord=2)) + 4.0 * smoothing, 1e-12)
    for iteration in range(200):
        if iteration % 20 == 0:
            _check_stop(should_stop)
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            gradient = amplitudes @ gram - cross
        if not np.isfinite(gradient).all():
            return np.full_like(amplitudes, np.nan)
        if len(edge_weights):
            delta = (amplitudes[1:] - amplitudes[:-1]) * edge_weights[:, None]
            gradient[:-1] -= smoothing * delta
            gradient[1:] += smoothing * delta
        updated = np.clip(amplitudes - gradient / lipschitz, 0.0, None)
        change = float(np.linalg.norm(updated - amplitudes))
        scale = max(float(np.linalg.norm(amplitudes)), 1.0)
        amplitudes = updated
        if change / scale <= 1e-7:
            break
    amplitudes[~np.isfinite(amplitudes)] = np.nan
    return amplitudes


def _check_stop(should_stop):
    if should_stop is not None and should_stop():
        raise PeakFitCancelled("the overlap fit was cancelled")


def _failed_global(fit, reason):
    result = dict(fit)
    n_components = int(np.asarray(fit["components"]).shape[1])
    result.update(
        {
            "model": "joint-temporal-v2",
            "status": "unresolved",
            "reason": reason,
            "amplitudes": np.full((0, n_components), np.nan),
            "selected_lambda": None,
            "lambda_grid": [],
            "held_out_relative_rmse": None,
            "held_out_standard_error": None,
            "best_reduced_relative_rmse": None,
            "held_out_improvement": None,
            "active_cycle_counts": [],
            "minimum_active_cycles": None,
            "failed_gates": [reason],
        }
    )
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
