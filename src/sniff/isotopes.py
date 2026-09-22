"""Natural-isotope envelopes and guarded PTR-MS spillover correction."""

from __future__ import annotations

import re

import numpy as np

PROTON = 1.007276
MODEL_VERSION = "natural-abundance-v1"
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


def build_isotope_plan(peaks, R_phys=2400.0):
    """Derive auxiliary channels and shared parent/isotope observations.

    Only explicitly assigned formulas participate. Derived channels remain runtime
    products and therefore cannot become accidental analyte rows or rate-constant keys.
    """
    parent_masses = [float(peak["mz"]) for peak in peaks]
    extraction_masses = list(parent_masses)
    parents = []
    warnings = []
    for peak in peaks:
        formula = peak.get("formula")
        if not formula:
            continue
        parent_mz = float(peak["mz"])
        try:
            model = formula_isotope_model(formula)
        except ValueError as exc:
            warnings.append({"mz": parent_mz, "formula": formula, "reason": str(exc)})
            continue
        channels = []
        for channel in model["channels"]:
            expected = parent_mz + channel["shift"]
            tolerance = max(expected / float(R_phys), 0.002)
            existing = _nearest(extraction_masses, expected)
            shared_parent = _nearest(parent_masses, expected)
            if existing is not None and abs(existing - expected) <= tolerance:
                observed_mz = existing
            else:
                observed_mz = expected
                extraction_masses.append(observed_mz)
            overlaps_parent = (
                shared_parent
                if shared_parent is not None
                and shared_parent != parent_mz
                and abs(shared_parent - expected) <= tolerance
                else None
            )
            channels.append(
                {
                    **channel,
                    "mz": expected,
                    "observation_mz": observed_mz,
                    "overlaps_parent_mz": overlaps_parent,
                }
            )
        parents.append(
            {
                "mz": parent_mz,
                "formula": formula,
                "version": model["version"],
                "monoisotopic_fraction": model["monoisotopic_fraction"],
                "channels": channels,
            }
        )
    return {
        "version": MODEL_VERSION,
        "parents": parents,
        "analyte_masses": parent_masses,
        "extraction_masses": sorted(set(extraction_masses)),
        "warnings": warnings,
    }


def correct_parent_signals(corrected, plan, *, abundance_basis="unknown"):
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
