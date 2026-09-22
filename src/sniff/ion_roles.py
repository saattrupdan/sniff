"""Conservative alternative-ion and charge-state candidate hypotheses."""

from __future__ import annotations

import math

import numpy as np

from . import formula_id, isotopes

MODEL = "ptr-ion-hypotheses-v1"
WATER_MASS = formula_id.formula_mass({"H": 2, "O": 1})
HYDROGEN_MASS = formula_id.MONO["H"]
ELECTRON_MASS = 0.000548579909
MIN_LEVEL_CORRELATION = 0.80
MIN_CHANGE_CORRELATION = 0.50


def annotate_ion_candidates(
    peaks,
    compound_catalogue,
    context,
    *,
    traces=None,
    drift=1.0,
    max_candidates=5,
):
    """Attach formula/compound proposals under explicit non-default ion chemistry.

    Direct ``[M+H]+`` candidates remain in ``candidates``. Alternative hypotheses are
    deliberately separate and never assignment-eligible: exact mass can support their
    formula arithmetic, but it cannot establish that the proposed ion pathway occurred.
    """
    traces = traces or {}
    reagent = str(context.get("reagent") or "").upper()
    context_available = context.get("status") == "available"
    h3o_mode = context_available and ("H3O" in reagent or "H₃O" in reagent)
    no_mode = context_available and "NO+" in reagent
    o2_mode = context_available and ("O2+" in reagent or "O₂+" in reagent)
    charge_evidence = _charge_state_evidence(peaks, traces)

    for index, peak in enumerate(peaks):
        peak.pop("ion_candidates", None)
        if context_available and not h3o_mode and (no_mode or o2_mode):
            for candidate in peak.get("candidates") or []:
                candidate["assignment_eligible"] = False
                candidate["reagent_compatible"] = False
                candidate["assignment_limitation"] = (
                    "a direct [M+H]+ assignment is incompatible with the validated "
                    f"active reagent {reagent}"
                )
        hypotheses = []
        if h3o_mode:
            hypotheses.extend(
                [
                    {
                        "kind": "dehydrated-protonated",
                        "notation": "[M+H-H2O]+",
                        "offset": formula_id.PROTON - WATER_MASS,
                        "charge": 1,
                        "reason": (
                            "water-loss product-ion proposal in the measured H3O+ "
                            "reaction context"
                        ),
                    },
                    {
                        "kind": "hydrated-protonated",
                        "notation": "[M+H+H2O]+",
                        "offset": formula_id.PROTON + WATER_MASS,
                        "charge": 1,
                        "reason": (
                            "hydrated product-ion proposal in the measured H3O+ "
                            "reaction context"
                        ),
                    },
                ]
            )
        if o2_mode:
            hypotheses.append(
                {
                    "kind": "charge-transfer",
                    "notation": "M+",
                    "offset": -ELECTRON_MASS,
                    "charge": 1,
                    "reason": "the validated active reagent is O2+ in this run",
                }
            )
        if no_mode:
            hypotheses.append(
                {
                    "kind": "hydride-abstraction",
                    "notation": "[M-H]+",
                    "offset": -HYDROGEN_MASS - ELECTRON_MASS,
                    "charge": 1,
                    "reason": "the validated active reagent is NO+ in this run",
                }
            )
        if index in charge_evidence:
            evidence = charge_evidence[index]
            charge = int(evidence["charge"])
            hypotheses.append(
                {
                    "kind": "multiply-charged-protonated",
                    "notation": f"[M+{charge}H]{charge}+",
                    "offset": charge * formula_id.PROTON,
                    "charge": charge,
                    "reason": (
                        f"two resolved isotope satellites support charge {charge}"
                    ),
                    "charge_evidence": evidence,
                }
            )

        candidates = []
        for hypothesis in hypotheses:
            candidates.extend(
                _candidates_for_hypothesis(
                    peak,
                    hypothesis,
                    compound_catalogue,
                    drift=drift,
                    max_candidates=max_candidates,
                )
            )
        unique = {}
        for candidate in candidates:
            key = (candidate["ion_kind"], candidate["formula"])
            current = unique.get(key)
            if current is None or _candidate_order(candidate) < _candidate_order(current):
                unique[key] = candidate
        retained = sorted(unique.values(), key=_candidate_order)[:max_candidates]
        if retained:
            peak["ion_candidates"] = retained
    return peaks


def candidate_summaries(peak):
    """Return deduplicated formula/name summaries from every candidate route."""
    summaries = []
    seen = set()
    for candidate in [
        *(peak.get("candidates") or []),
        *(peak.get("ion_candidates") or []),
    ]:
        formula = candidate.get("formula")
        names = _candidate_names(candidate)
        key = (formula, tuple(names))
        if formula and key not in seen:
            summaries.append(
                {
                    "formula": formula,
                    "names": names,
                    **(
                        {"ion_notation": candidate["ion_notation"]}
                        if candidate.get("ion_notation")
                        else {}
                    ),
                    **(
                        {"iso_pred": candidate["iso_pred"]}
                        if candidate.get("iso_pred")
                        else {}
                    ),
                }
            )
            seen.add(key)
    return summaries


def coverage_summary(peaks):
    """Return candidate coverage over canonical unresolved-peak components."""
    category_order = (
        "direct_protonated_formula",
        "alternative_ion_formula",
        "inherited_parent_candidate",
        "authored_assignment",
        "interpreted_noncompound",
        "unknown",
    )
    categories = {category: 0 for category in category_order}
    named_compound = 0
    candidate_weight = 0.0
    total_weight = 0.0
    measured_weighting = bool(peaks) and all(
        "role_trace_status" in peak for peak in peaks
    )
    weighted_available = measured_weighting
    components = _canonical_components(peaks)
    for indexes in components:
        component_peaks = [peaks[index] for index in indexes]
        member_categories = [_peak_category(peak) for peak in component_peaks]
        category = min(member_categories, key=category_order.index)
        categories[category] += 1
        named = any(_peak_has_name(peak) for peak in component_peaks)
        named_compound += int(named)
        valid_signals = []
        component_signal_invalid = False
        for peak in component_peaks:
            if (
                peak.get("role_trace_status") != "available"
                or peak.get("role_signal") is None
            ):
                continue
            signal = float(peak["role_signal"])
            if math.isfinite(signal):
                valid_signals.append(max(0.0, signal))
            else:
                component_signal_invalid = True
        if measured_weighting:
            if valid_signals and not component_signal_invalid:
                weight = max(valid_signals)
            else:
                weighted_available = False
                weight = 0.0
        else:
            weight = 0.0
        total_weight += weight
        if category in category_order[:3]:
            candidate_weight += weight
    total_components = len(components)
    with_candidate = sum(categories[category] for category in category_order[:3])
    raw_with_candidate = sum(
        _peak_category(peak) in category_order[:3] for peak in peaks
    )
    signal_source = (
        "mean extracted Raw signal"
        if measured_weighting and weighted_available
        else "unavailable because at least one canonical component has no separable Raw trace"
        if measured_weighting
        else "unavailable because extracted Raw trace status is missing"
    )
    return {
        "model": MODEL,
        "denominator": "canonical unresolved-peak components",
        "detected_peaks": len(peaks),
        "total_components": total_components,
        "with_formula_or_linked_candidate": with_candidate,
        "candidate_coverage_percent": round(
            100.0 * with_candidate / total_components, 1
        )
        if total_components
        else 0.0,
        "with_named_compound_candidate": named_compound,
        "named_compound_coverage_percent": round(
            100.0 * named_compound / total_components, 1
        )
        if total_components
        else 0.0,
        "signal_weighted_candidate_coverage_percent": (
            round(100.0 * candidate_weight / total_weight, 1)
            if weighted_available and total_weight
            else None
        ),
        "signal_weight_source": signal_source,
        "peak_level_candidate_count": raw_with_candidate,
        "categories": categories,
        "limitation": (
            "candidate coverage counts proposals, not unique identities or calibrated "
            "identification confidence"
        ),
    }


def _peak_category(peak):
    interpretations = peak.get("interpretation_candidates") or []
    if any(
        item.get("kind") == "isotope" and item.get("candidate_formulas")
        for item in interpretations
    ):
        return "inherited_parent_candidate"
    if peak.get("candidates"):
        return "direct_protonated_formula"
    if peak.get("ion_candidates"):
        return "alternative_ion_formula"
    if any(item.get("candidate_formulas") for item in interpretations):
        return "inherited_parent_candidate"
    if any(item.get("kind") == "authored" for item in interpretations):
        return "authored_assignment"
    if any(
        not str(item.get("kind", "")).startswith("unresolved")
        for item in interpretations
    ):
        return "interpreted_noncompound"
    return "unknown"


def _peak_has_name(peak):
    candidates = [
        *(peak.get("candidates") or []),
        *(peak.get("ion_candidates") or []),
    ]
    return any(_candidate_names(candidate) for candidate in candidates) or any(
        item.get("compound_candidates")
        for item in peak.get("interpretation_candidates") or []
    )


def _canonical_components(peaks):
    parents = list(range(len(peaks)))

    def root(index):
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    def join(left, right):
        left, right = root(left), root(right)
        if left != right:
            parents[right] = left

    masses = [float(peak.get("apex", peak.get("mz", 0.0))) for peak in peaks]
    for index, peak in enumerate(peaks):
        overlap = peak.get("overlap") or {}
        if overlap.get("level") != "unresolved" or overlap.get("neighbor") is None:
            continue
        neighbour = float(overlap["neighbor"])
        other = min(range(len(peaks)), key=lambda item: abs(masses[item] - neighbour))
        if other != index and abs(masses[other] - neighbour) <= 0.001:
            join(index, other)
    grouped = {}
    for index in range(len(peaks)):
        grouped.setdefault(root(index), []).append(index)
    return list(grouped.values())


def _candidates_for_hypothesis(
    peak,
    hypothesis,
    compound_catalogue,
    *,
    drift,
    max_candidates,
):
    observed_mz = float(peak["mz"]) / float(drift)
    charge = int(hypothesis["charge"])
    offset = float(hypothesis["offset"])
    neutral_mass = charge * observed_mz - offset
    if neutral_mass <= 0:
        return []
    tolerance = peak.get("formula_tolerance") or {}
    strict_ppm = float(tolerance.get("ppm", 10.0))
    proposal_ppm = float(tolerance.get("proposal_ppm", 200.0))
    synthetic_mz = neutral_mass + formula_id.PROTON
    ppm_scale = observed_mz * charge / synthetic_mz
    extra_formulas = []
    if (
        neutral_mass > formula_id.MAX_ENUMERATED_NEUTRAL_MASS
        and hasattr(compound_catalogue, "formulas_in_mass_range")
    ):
        neutral_tolerance = observed_mz * charge * proposal_ppm / 1e6
        extra_formulas = [
            item["formula"]
            for item in compound_catalogue.formulas_in_mass_range(
                neutral_mass,
                neutral_tolerance,
            )
        ]
    scored = formula_id.score_peak(
        synthetic_mz,
        1.0,
        obs_ratios=None,
        tolerance_ppm=strict_ppm * ppm_scale,
        mass_sigma_ppm=(
            float(tolerance.get("score_sigma_ppm", 4.0)) * ppm_scale
        ),
        proposal_tolerance_ppm=proposal_ppm * ppm_scale,
        extra_formulas=extra_formulas,
        max_candidates=max(250, max_candidates),
    )
    output = []
    for candidate in scored:
        counts = isotopes.parse_formula(candidate["formula"])
        if hypothesis["kind"] == "dehydrated-protonated" and (
            counts.get("O", 0) < 1 or counts.get("H", 0) < 2
        ):
            continue
        theoretical_mz = (
            formula_id.formula_mass(counts) + offset
        ) / charge
        delta_ppm = (observed_mz - theoretical_mz) / theoretical_mz * 1e6
        if abs(delta_ppm) > proposal_ppm:
            continue
        mass_eligible = abs(delta_ppm) <= strict_ppm
        enriched = dict(candidate)
        enriched.update(
            {
                "model": MODEL,
                "ion_kind": hypothesis["kind"],
                "ion_notation": hypothesis["notation"],
                "charge": charge,
                "neutral_mass": round(formula_id.formula_mass(counts), 6),
                "ion_mz": round(theoretical_mz, 4),
                "delta_mDa": round((observed_mz - theoretical_mz) * 1000.0, 1),
                "delta_ppm": round(delta_ppm, 2),
                "mass_match": (
                    "alternative-ion-within-run-tolerance"
                    if mass_eligible
                    else "alternative-ion-broad-proposal"
                ),
                "mass_eligible": mass_eligible,
                "assignment_eligible": False,
                "pathway_reason": hypothesis["reason"],
                "limitation": (
                    "exact mass supports this formula/ion arithmetic but does not "
                    "establish that the ion pathway occurred"
                ),
            }
        )
        if hypothesis.get("charge_evidence"):
            enriched["charge_evidence"] = hypothesis["charge_evidence"]
        output.append(enriched)
    retained = sorted(output, key=_candidate_order)[:max_candidates]
    return compound_catalogue.enrich_candidates(retained)


def _candidate_order(candidate):
    return (
        not bool(candidate.get("mass_eligible")),
        abs(float(candidate.get("delta_ppm", math.inf))),
        candidate.get("ion_kind", ""),
        candidate.get("formula", ""),
    )


def _candidate_names(candidate):
    names = []
    seen = set()
    for value in [
        candidate.get("preferred_name"),
        candidate.get("name"),
        *(candidate.get("names") or []),
        *(entry.get("name") for entry in candidate.get("catalogue") or []),
    ]:
        value = " ".join(str(value or "").split())
        key = value.casefold()
        if value and key not in seen:
            names.append(value)
            seen.add(key)
    return names


def _charge_state_evidence(peaks, traces):
    evidence = {}
    masses = [float(peak["mz"]) for peak in peaks]
    for parent_index, peak in enumerate(peaks):
        parent_mz = masses[parent_index]
        parent_height = max(0.0, float(peak.get("height", 0.0)))
        tolerance = max(
            0.002,
            2.0
            * float((peak.get("formula_tolerance") or {}).get("mDa", 1.0))
            / 1000.0,
        )
        for charge in (2, 3):
            satellites = []
            valid = True
            for order in (1, 2):
                expected = parent_mz + order * formula_id.DM1 / charge
                index = min(
                    range(len(peaks)), key=lambda item: abs(masses[item] - expected)
                )
                residual = masses[index] - expected
                height = max(0.0, float(peaks[index].get("height", 0.0)))
                if (
                    index == parent_index
                    or abs(residual) > tolerance
                    or height <= 0
                    or height >= parent_height
                ):
                    valid = False
                    break
                correlation = _trace_correlation(
                    traces.get(parent_mz), traces.get(masses[index])
                )
                if correlation is None or not correlation["supported"]:
                    valid = False
                    break
                satellites.append(
                    {
                        "order": order,
                        "mz": round(masses[index], 4),
                        "spacing_residual_mDa": round(residual * 1000.0, 2),
                        **correlation,
                    }
                )
            if valid:
                evidence[parent_index] = {
                    "model": MODEL,
                    "charge": charge,
                    "spacing_da": round(formula_id.DM1 / charge, 6),
                    "satellites": satellites,
                    "limitation": (
                        "isotope spacing supports charge state but does not establish "
                        "a neutral structure"
                    ),
                }
                break
    return evidence


def _trace_correlation(left, right):
    if left is None or right is None:
        return None
    left = np.asarray(left, dtype=np.float64)
    right = np.asarray(right, dtype=np.float64)
    keep = np.isfinite(left) & np.isfinite(right)
    left, right = left[keep], right[keep]
    if left.size < 8 or np.ptp(left) <= 0 or np.ptp(right) <= 0:
        return None
    level = float(np.corrcoef(left, right)[0, 1])
    change = float(np.corrcoef(np.diff(left), np.diff(right))[0, 1])
    if not np.isfinite([level, change]).all():
        return None
    return {
        "level_correlation": round(level, 3),
        "change_correlation": round(change, 3),
        "supported": (
            level >= MIN_LEVEL_CORRELATION and change >= MIN_CHANGE_CORRELATION
        ),
    }
