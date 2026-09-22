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
        if peak.get("candidates"):
            continue
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
    """Return non-overlapping candidate/interpretation coverage counts."""
    direct_formula = 0
    alternative_formula = 0
    inherited_formula = 0
    named_compound = 0
    interpreted_noncompound = 0
    unknown = 0
    for peak in peaks:
        direct = peak.get("candidates") or []
        alternative = peak.get("ion_candidates") or []
        interpretations = peak.get("interpretation_candidates") or []
        inherited = any(item.get("candidate_formulas") for item in interpretations)
        named = any(_candidate_names(candidate) for candidate in [*direct, *alternative])
        named = named or any(
            item.get("compound_candidates") for item in interpretations
        )
        if direct:
            direct_formula += 1
        elif alternative:
            alternative_formula += 1
        elif inherited:
            inherited_formula += 1
        elif any(item.get("kind") != "unresolved" for item in interpretations):
            interpreted_noncompound += 1
        else:
            unknown += 1
        if named:
            named_compound += 1
    total = len(peaks)
    with_candidate = direct_formula + alternative_formula + inherited_formula
    return {
        "model": MODEL,
        "total_peaks": total,
        "with_formula_or_linked_candidate": with_candidate,
        "candidate_coverage_percent": round(100.0 * with_candidate / total, 1)
        if total
        else 0.0,
        "with_named_compound_candidate": named_compound,
        "named_compound_coverage_percent": round(100.0 * named_compound / total, 1)
        if total
        else 0.0,
        "categories": {
            "direct_protonated_formula": direct_formula,
            "alternative_ion_formula": alternative_formula,
            "inherited_parent_candidate": inherited_formula,
            "interpreted_noncompound": interpreted_noncompound,
            "unknown": unknown,
        },
        "limitation": (
            "candidate coverage counts proposals, not unique identities or calibrated "
            "identification confidence"
        ),
    }


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
    scored = formula_id.score_peak(
        synthetic_mz,
        1.0,
        obs_ratios=None,
        tolerance_ppm=strict_ppm * ppm_scale,
        mass_sigma_ppm=(
            float(tolerance.get("score_sigma_ppm", 4.0)) * ppm_scale
        ),
        proposal_tolerance_ppm=proposal_ppm * ppm_scale,
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
