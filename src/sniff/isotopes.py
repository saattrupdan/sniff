"""Natural-isotope envelopes and guarded PTR-MS spillover correction."""

from __future__ import annotations

import re

import numpy as np

PROTON = 1.007276
MODEL_VERSION = "natural-abundance-v1"
ENVELOPE_MODEL_VERSION = "formula-envelope-v2"
_FORMULA_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")

# Monoisotopic masses and representative terrestrial abundances. Values are kept in
# source so the model is reviewable and frozen with every Sniff release. The model is
# based on IUPAC conventional isotope abundances and exact isotope masses (2016 table).
ISOTOPES = {
    "H": ((1.00782503223, 0.999885), (2.01410177812, 0.000115)),
    "C": ((12.0, 0.9893), (13.00335483507, 0.0107)),
    "N": ((14.00307400443, 0.99636), (15.00010889888, 0.00364)),
    "O": (
        (15.99491461957, 0.99757),
        (16.99913175650, 0.00038),
        (17.99915961286, 0.00205),
    ),
    "F": ((18.99840316273, 1.0),),
    "P": ((30.97376199842, 1.0),),
    "S": (
        (31.97207117440, 0.9499),
        (32.97145890980, 0.0075),
        (33.96786700400, 0.0425),
        (35.96708071000, 0.0001),
    ),
    "Cl": ((34.968852682, 0.7576), (36.965902602, 0.2424)),
    "Br": ((78.9183376, 0.5069), (80.9162897, 0.4931)),
    "Si": (
        (27.97692653465, 0.92223),
        (28.97649466490, 0.04685),
        (29.97377013600, 0.03092),
    ),
    "I": ((126.9044719, 1.0),),
}


def parse_formula(formula):
    """Parse a neutral molecular formula into element counts.

    Raises:
        ValueError:
            If the formula is malformed or contains an unsupported element.
    """
    if not isinstance(formula, str) or not formula:
        raise ValueError("formula must be a non-empty string")
    counts = {}
    position = 0
    for match in _FORMULA_TOKEN.finditer(formula):
        if match.start() != position:
            raise ValueError(f"malformed formula: {formula}")
        element, raw_count = match.groups()
        if element not in ISOTOPES:
            raise ValueError(f"unsupported isotope element: {element}")
        count = int(raw_count or "1")
        if count <= 0:
            raise ValueError(f"invalid atom count for {element}")
        counts[element] = counts.get(element, 0) + count
        position = match.end()
    if position != len(formula) or not counts:
        raise ValueError(f"malformed formula: {formula}")
    return counts


def formula_isotope_model(formula, protonated=True):
    """Return the natural M, M+1 and M+2 model for an assigned formula."""
    counts = parse_formula(formula)
    ion_counts = dict(counts)
    if protonated:
        # Preserve Sniff's established [M+H]+ isotope convention. The contribution
        # from deuterium on this additional hydrogen is tiny but versioned here.
        ion_counts["H"] = ion_counts.get("H", 0) + 1

    probability = np.zeros(3, dtype=np.float64)
    shift_sum = np.zeros(3, dtype=np.float64)
    probability[0] = 1.0
    for element, count in ion_counts.items():
        isotopes = ISOTOPES[element]
        mono_mass = isotopes[0][0]
        for _ in range(count):
            next_probability = np.zeros(3, dtype=np.float64)
            next_shift_sum = np.zeros(3, dtype=np.float64)
            for current_order in range(3):
                if probability[current_order] <= 0:
                    continue
                for mass, abundance in isotopes:
                    delta = mass - mono_mass
                    isotope_order = int(round(delta))
                    order = current_order + isotope_order
                    if order > 2:
                        continue
                    next_probability[order] += probability[current_order] * abundance
                    next_shift_sum[order] += abundance * (
                        shift_sum[current_order] + probability[current_order] * delta
                    )
            probability, shift_sum = next_probability, next_shift_sum

    neutral_mass = sum(
        ISOTOPES[element][0][0] * count for element, count in counts.items()
    )
    parent_mz = neutral_mass + PROTON
    f0 = float(probability[0])
    channels = []
    for order in (1, 2):
        fraction = float(probability[order])
        shift = float(shift_sum[order] / fraction) if fraction > 0 else float(order)
        channels.append(
            {
                "order": order,
                "mz": parent_mz + shift,
                "shift": shift,
                "fraction": fraction,
                "ratio": fraction / f0 if f0 > 0 else 0.0,
            }
        )
    return {
        "version": MODEL_VERSION,
        "formula": formula,
        "parent_mz": parent_mz,
        "monoisotopic_fraction": f0,
        "channels": channels,
    }


def isotope_ratios(counts, protonated=True):
    """Return M+1/M and M+2/M for element counts."""
    formula = "".join(
        element + (str(count) if count != 1 else "")
        for element, count in counts.items()
        if count
    )
    model = formula_isotope_model(formula, protonated=protonated)
    return tuple(channel["ratio"] for channel in model["channels"])


def build_isotope_plan(peaks, R_phys=2400.0, model="formula-v1"):
    """Derive auxiliary channels and shared parent/isotope observations.

    Only explicitly assigned formulas participate. Derived channels remain runtime
    products and therefore cannot become accidental analyte rows or rate-constant keys.
    """
    parent_masses = [float(peak["mz"]) for peak in peaks]
    parent_specs = []
    channel_specs = []
    warnings = []
    for parent_index, peak in enumerate(peaks):
        formula = peak.get("formula")
        if not formula:
            continue
        parent_mz = float(peak["mz"])
        try:
            formula_model = formula_isotope_model(formula)
        except ValueError as exc:
            warnings.append({"mz": parent_mz, "formula": formula, "reason": str(exc)})
            continue
        parent_spec = {
            "parent_index": parent_index,
            "mz": parent_mz,
            "formula": formula,
            "version": formula_model["version"],
            "monoisotopic_fraction": formula_model["monoisotopic_fraction"],
            "channels": [],
        }
        parent_specs.append(parent_spec)
        for channel in formula_model["channels"]:
            expected = parent_mz + channel["shift"]
            tolerance = max(expected / float(R_phys), 0.002)
            candidates = [
                (abs(mass - expected), index, mass)
                for index, mass in enumerate(parent_masses)
                if index != parent_index and abs(mass - expected) <= tolerance
            ]
            shared_parent = min(candidates)[2] if candidates else None
            channel_spec = {
                **channel,
                "mz": expected,
                "tolerance": tolerance,
                "overlaps_parent_mz": shared_parent,
                "parent": parent_spec,
            }
            channel_specs.append(channel_spec)
            parent_spec["channels"].append(channel_spec)

    auxiliary = [item for item in channel_specs if item["overlaps_parent_mz"] is None]
    _assign_channel_observations(auxiliary)
    extraction_masses = list(parent_masses)
    parents = []
    for parent_spec in parent_specs:
        channels = []
        for channel in parent_spec["channels"]:
            observed_mz = channel.get("observation_mz")
            if channel["overlaps_parent_mz"] is not None:
                observed_mz = channel["overlaps_parent_mz"]
            extraction_masses.append(observed_mz)
            channels.append(
                {
                    key: value
                    for key, value in channel.items()
                    if key not in {"parent", "tolerance"}
                }
                | {"observation_mz": observed_mz}
            )
        parents.append(
            {
                key: value
                for key, value in parent_spec.items()
                if key not in {"parent_index", "channels"}
            }
            | {"channels": channels}
        )
    return {
        "version": model,
        "isotope_model": MODEL_VERSION,
        "parents": parents,
        "analyte_masses": parent_masses,
        "extraction_masses": sorted(set(extraction_masses)),
        "warnings": warnings,
    }


def _assign_channel_observations(channels):
    remaining = sorted(channels, key=lambda item: (item["mz"], item["order"]))
    groups = []
    for channel in remaining:
        if not groups:
            groups.append([channel])
            continue
        previous = groups[-1][-1]
        if channel["mz"] - previous["mz"] <= max(
            channel["tolerance"], previous["tolerance"]
        ):
            groups[-1].append(channel)
        else:
            groups.append([channel])
    for group in groups:
        observation = float(np.mean([item["mz"] for item in group]))
        for channel in group:
            channel["observation_mz"] = observation


def correct_parent_signals(corrected, plan, *, abundance_basis="unknown"):
    """Fit or subtract formula isotope envelopes in corrected signal space."""
    if plan.get("version") == ENVELOPE_MODEL_VERSION:
        return _fit_parent_envelopes(
            corrected,
            plan,
            abundance_basis=abundance_basis,
        )
    return _correct_parent_signals_legacy(
        corrected,
        plan,
        abundance_basis=abundance_basis,
    )


def _fit_parent_envelopes(corrected, plan, *, abundance_basis):
    parent_models = list(plan.get("parents", []))
    if not parent_models:
        return {}, []
    components = _envelope_components(parent_models)
    output = {}
    diagnostics = []
    for component_index, parent_indices in enumerate(components, start=1):
        parents = [parent_models[index] for index in parent_indices]
        observations = sorted(
            {
                float(parent["mz"])
                for parent in parents
            }
            | {
                float(channel["observation_mz"])
                for parent in parents
                for channel in parent["channels"]
            }
        )
        design = np.zeros((len(observations), len(parents)), dtype=np.float64)
        observation_index = {mass: index for index, mass in enumerate(observations)}
        for column, parent in enumerate(parents):
            design[observation_index[float(parent["mz"])], column] += 1.0
            for channel in parent["channels"]:
                row = observation_index[float(channel["observation_mz"])]
                design[row, column] += float(channel["ratio"])
        available = all(mass in corrected for mass in observations)
        ncycles = (
            len(np.asarray(corrected[observations[0]]))
            if available and observations
            else 0
        )
        amplitudes = np.full((len(parents), ncycles), np.nan, dtype=np.float64)
        standard_errors = np.full_like(amplitudes, np.nan)
        residuals = np.full(ncycles, np.nan, dtype=np.float64)
        conditions = []
        if available and np.linalg.matrix_rank(design) == len(parents):
            values = np.vstack(
                [np.asarray(corrected[mass], dtype=np.float64) for mass in observations]
            )
            for cycle in range(ncycles):
                observed = values[:, cycle]
                if not np.isfinite(observed).all():
                    continue
                variance = np.maximum(np.abs(observed), 1.0)
                whiten = np.sqrt(variance)
                weighted_design = design / whiten[:, None]
                weighted_observed = observed / whiten
                with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
                    condition = float(np.linalg.cond(weighted_design))
                conditions.append(condition)
                if not np.isfinite(condition) or condition > 1e6:
                    continue
                estimate, active = _solve_weighted_nonnegative(
                    weighted_design,
                    weighted_observed,
                )
                if not active or not np.isfinite(estimate).all():
                    continue
                with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
                    fitted = design @ estimate
                    scaled_residual = (observed - fitted) / whiten
                if not np.isfinite(scaled_residual).all():
                    continue
                degrees = max(len(observations) - len(active), 1)
                reduced_error = max(
                    float(np.dot(scaled_residual, scaled_residual) / degrees),
                    1.0,
                )
                with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
                    pseudo_inverse = np.linalg.pinv(weighted_design[:, active])
                if not np.isfinite(pseudo_inverse).all():
                    continue
                active_error = np.hypot.reduce(pseudo_inverse, axis=1)
                active_error *= np.sqrt(reduced_error)
                error = np.full(len(parents), np.nan, dtype=np.float64)
                error[active] = active_error
                amplitudes[:, cycle] = estimate
                standard_errors[:, cycle] = error
                residuals[cycle] = float(np.sqrt(np.mean(scaled_residual**2)))
        finite_residuals = residuals[np.isfinite(residuals)]
        component_condition = max(conditions) if conditions else None
        component_residual = (
            float(np.median(finite_residuals)) if finite_residuals.size else None
        )
        for row, parent in enumerate(parents):
            fitted = amplitudes[row]
            error = standard_errors[row]
            valid = np.isfinite(fitted) & np.isfinite(error) & (fitted > 0)
            relative = error[valid] / fitted[valid]
            relative_uncertainty = (
                float(np.median(relative)) if relative.size else None
            )
            enough_cycles = int(np.count_nonzero(valid)) >= max(1, ncycles // 2)
            reliable = (
                available
                and enough_cycles
                and component_condition is not None
                and component_condition <= 1e6
                and component_residual is not None
                and component_residual <= 3.0
                and relative_uncertainty is not None
                and relative_uncertainty <= 0.5
            )
            if reliable:
                trace = fitted.copy()
                if abundance_basis == "total":
                    fraction = float(parent["monoisotopic_fraction"])
                    if fraction > 0:
                        trace[np.isfinite(trace)] /= fraction
                status = "envelope-fitted"
            else:
                trace = np.full(ncycles, np.nan, dtype=np.float64)
                status = _envelope_withheld_status(
                    available=available,
                    enough_cycles=enough_cycles,
                    condition=component_condition,
                    residual=component_residual,
                    relative_uncertainty=relative_uncertainty,
                )
            mass = float(parent["mz"])
            output[mass] = trace
            diagnostics.append(
                {
                    "mz": mass,
                    "formula": parent["formula"],
                    "status": status,
                    "model": ENVELOPE_MODEL_VERSION,
                    "component": component_index,
                    "component_parents": [float(item["mz"]) for item in parents],
                    "observations": observations,
                    "rank": int(np.linalg.matrix_rank(design)),
                    "condition": component_condition,
                    "median_weighted_residual": component_residual,
                    "fitted_cycles": int(np.count_nonzero(valid)),
                    "total_cycles": ncycles,
                    "median_relative_uncertainty": relative_uncertainty,
                    "uncertainty_method": "weighted-design-pseudoinverse-v1",
                    "abundance_basis": abundance_basis,
                    "abundance_applied": reliable and abundance_basis == "total",
                }
            )
    return output, diagnostics


def _solve_weighted_nonnegative(design, observed):
    active = list(range(design.shape[1]))
    estimate = np.zeros(design.shape[1], dtype=np.float64)
    while active:
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            fitted, _, rank, _ = np.linalg.lstsq(
                design[:, active],
                observed,
                rcond=None,
            )
        if rank != len(active) or not np.isfinite(fitted).all():
            return estimate, []
        negative = [index for index, value in enumerate(fitted) if value < 0]
        if not negative:
            estimate[active] = fitted
            return estimate, active
        worst = min(negative, key=lambda index: fitted[index])
        del active[worst]
    return estimate, []


def _envelope_components(parents):
    observations = []
    for parent in parents:
        observations.append(
            {float(parent["mz"])}
            | {
                float(channel["observation_mz"])
                for channel in parent["channels"]
            }
        )
    remaining = set(range(len(parents)))
    components = []
    while remaining:
        component = {remaining.pop()}
        changed = True
        while changed:
            changed = False
            occupied = set().union(*(observations[index] for index in component))
            linked = {
                index
                for index in remaining
                if observations[index] & occupied
            }
            if linked:
                component.update(linked)
                remaining.difference_update(linked)
                changed = True
        components.append(sorted(component))
    return components


def _envelope_withheld_status(
    *, available, enough_cycles, condition, residual, relative_uncertainty
):
    if not available:
        return "withheld-missing-envelope-channel"
    if condition is None or condition > 1e6:
        return "withheld-ill-conditioned-envelope"
    if not enough_cycles:
        return "withheld-insufficient-envelope-cycles"
    if residual is None or residual > 3.0:
        return "withheld-envelope-residual"
    if relative_uncertainty is None or relative_uncertainty > 0.5:
        return "withheld-envelope-uncertainty"
    return "withheld-envelope-fit"


def _correct_parent_signals_legacy(corrected, plan, *, abundance_basis="unknown"):
    """Remove identifiable isotope spillover from parent corrected traces.

    Args:
        corrected:
            Mapping from extracted m/z to transmission-corrected per-cycle traces.
        plan:
            Output of :func:`build_isotope_plan`.
        abundance_basis:
            ``total`` applies monoisotopic abundance scaling. ``calibrated`` and
            ``unknown`` retain the calibrated parent-signal basis.

    Returns:
        Tuple ``(parent traces, diagnostics)``.
    """
    net = {
        mass: np.asarray(values, dtype=np.float64).copy()
        for mass, values in corrected.items()
    }
    diagnostics = []
    parent_models = sorted(plan.get("parents", []), key=lambda item: item["mz"])
    for parent in parent_models:
        mass = parent["mz"]
        if mass not in net:
            continue
        original = net[mass].copy()
        contributions = []
        status = "not-needed"
        for source in parent_models:
            if source["mz"] >= mass or source["mz"] not in net:
                continue
            for channel in source["channels"]:
                if channel.get("overlaps_parent_mz") != mass:
                    continue
                source_signal = net[source["mz"]]
                finite_inputs = np.isfinite(net[mass]) & np.isfinite(source_signal)
                if not finite_inputs.any():
                    net[mass] = np.full_like(original, np.nan)
                    status = "withheld-unavailable-source-or-target"
                    contributions.append(
                        {
                            "source_mz": source["mz"],
                            "order": channel["order"],
                            "ratio": channel["ratio"],
                            "applied": False,
                        }
                    )
                    break
                expected_ratio = float(channel["ratio"])
                observed_ratio = float(
                    np.median(net[mass][finite_inputs] / source_signal[finite_inputs])
                )
                tolerance = max(
                    0.015 if channel["order"] == 1 else 0.008,
                    (0.50 if channel["order"] == 1 else 0.60) * expected_ratio,
                )
                if observed_ratio < max(0.0, expected_ratio - tolerance):
                    net[mass] = np.full_like(original, np.nan)
                    status = "withheld-isotope-ratio-deficit"
                    contributions.append(
                        {
                            "source_mz": source["mz"],
                            "order": channel["order"],
                            "ratio": expected_ratio,
                            "observed_ratio": observed_ratio,
                            "applied": False,
                        }
                    )
                    break
                contribution = source_signal * expected_ratio
                candidate = net[mass] - contribution
                candidate[~finite_inputs] = np.nan
                materially_negative = finite_inputs & (
                    candidate < -0.02 * np.maximum(np.abs(original), 1.0)
                )
                if materially_negative.any():
                    net[mass] = np.full_like(original, np.nan)
                    status = "withheld-negative-after-subtraction"
                    contributions.append(
                        {
                            "source_mz": source["mz"],
                            "order": channel["order"],
                            "ratio": channel["ratio"],
                            "applied": False,
                        }
                    )
                    break
                net[mass] = np.maximum(candidate, 0.0)
                status = "spillover-corrected"
                contributions.append(
                    {
                        "source_mz": source["mz"],
                        "order": channel["order"],
                        "ratio": channel["ratio"],
                        "applied": True,
                    }
                )
            if status.startswith("withheld-"):
                break
        channel_reports = []
        for channel in parent["channels"]:
            observation = channel["observation_mz"]
            observed_ratio = None
            ratio_status = "missing"
            if observation in corrected:
                denominator = corrected[mass]
                numerator = corrected[observation]
                valid = (
                    np.isfinite(denominator)
                    & np.isfinite(numerator)
                    & (denominator > 0)
                )
                if valid.any():
                    observed_ratio = float(
                        np.median(numerator[valid] / denominator[valid])
                    )
                    expected = float(channel["ratio"])
                    tolerance = max(
                        0.015 if channel["order"] == 1 else 0.008,
                        (0.50 if channel["order"] == 1 else 0.60) * expected,
                    )
                    if observed_ratio < max(0.0, expected - tolerance):
                        ratio_status = "deficient"
                    elif observed_ratio > expected + tolerance:
                        ratio_status = "excess-or-overlap"
                    else:
                        ratio_status = "consistent"
            channel_reports.append(
                {
                    "order": channel["order"],
                    "mz": channel["mz"],
                    "expected_ratio": channel["ratio"],
                    "observed_ratio": observed_ratio,
                    "status": ratio_status,
                    "overlaps_parent_mz": channel.get("overlaps_parent_mz"),
                }
            )
        diagnostics.append(
            {
                "mz": mass,
                "formula": parent["formula"],
                "status": status,
                "contributions": contributions,
                "channels": channel_reports,
                "monoisotopic_fraction": parent["monoisotopic_fraction"],
                "abundance_basis": abundance_basis,
                "abundance_applied": False,
            }
        )

    # Spillover ratios are relative to monoisotopic parent signal. Keep every source in
    # that basis until the complete low-to-high correction graph has been solved, then
    # scale final parent signals. Scaling inside the loop inflates later subtraction.
    if abundance_basis == "total":
        diagnostic_by_mass = {item["mz"]: item for item in diagnostics}
        for parent in parent_models:
            mass = parent["mz"]
            if mass not in net:
                continue
            finite = np.isfinite(net[mass])
            f0 = float(parent["monoisotopic_fraction"])
            if f0 > 0 and np.any(finite):
                net[mass][finite] /= f0
                diagnostic_by_mass[mass]["abundance_applied"] = True

    parent_masses = {parent["mz"] for parent in parent_models}
    return {mass: net[mass] for mass in parent_masses if mass in net}, diagnostics


def _nearest(values, target):
    if not values:
        return None
    return min(values, key=lambda value: abs(value - target))
