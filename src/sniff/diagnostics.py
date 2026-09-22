"""Read-only component-level diagnostics for curated PTR-MS reviews."""

from collections import Counter

from . import catalogue, ion_roles

_REASON_MODEL = "unresolved-component-reasons-v1"
_CATEGORY_ORDER = (
    "direct_protonated_formula",
    "alternative_ion_formula",
    "inherited_parent_candidate",
    "authored_assignment",
    "interpreted_noncompound",
    "unknown",
)


def component_report(peaks, *, include_resolved=False):
    """Explain why each canonical component is or is not chemically supported."""
    compounds = catalogue.CompoundCatalogue()
    records = []
    for component_id, indexes in enumerate(
        ion_roles._canonical_components(peaks), start=1
    ):
        members = [peaks[index] for index in indexes]
        categories = [ion_roles._peak_category(peak) for peak in members]
        category = min(categories, key=_CATEGORY_ORDER.index)
        if category != "unknown" and not include_resolved:
            continue
        representative = max(
            members,
            key=lambda peak: float(
                peak.get("role_signal")
                or peak.get("abundance")
                or peak.get("height")
                or 0.0
            ),
        )
        records.append(
            _component_record(
                component_id,
                indexes,
                members,
                representative,
                category,
                compounds,
            )
        )
    reasons = Counter(record["primary_reason"] for record in records)
    constraints = Counter(
        constraint
        for record in records
        for constraint in record["additional_constraints"]
    )
    return {
        "model": _REASON_MODEL,
        "denominator": "canonical unresolved-peak components",
        "component_count": len(records),
        "reason_counts": dict(sorted(reasons.items())),
        "constraint_counts": dict(sorted(constraints.items())),
        "components": records,
        "limitations": [
            "A diagnostic nearest formula outside the 200 ppm proposal gate is not a candidate.",
            "Formula and database matches cannot distinguish structural isomers.",
            "Unsupported ion pathways remain unknown rather than being inferred from nominal mass.",
        ],
    }


def compare_baseline_unknowns(baseline_peaks, current_peaks):
    """Explain outcomes for the unresolved cohort in a prior review payload."""
    baseline = component_report(baseline_peaks, include_resolved=True)
    current = component_report(current_peaks, include_resolved=True)
    current_by_indexes = {
        tuple(record["peak_indexes"]): record for record in current["components"]
    }
    outcomes = []
    newly_unknown = []
    transition_counts = Counter()
    for before in baseline["components"]:
        indexes = tuple(before["peak_indexes"])
        after = current_by_indexes.get(indexes)
        if after is None:
            raise ValueError(
                "baseline and current reviews do not have matching canonical components"
            )
        if before["category"] == "unknown":
            outcomes.append({"baseline": before, "current": after})
            transition_counts[after["category"]] += 1
        elif after["category"] == "unknown":
            newly_unknown.append(after["component_id"])
    return {
        "model": "baseline-unresolved-outcomes-v1",
        "baseline_component_count": len(outcomes),
        "transition_counts": dict(sorted(transition_counts.items())),
        "gained_support": sum(
            count
            for category, count in transition_counts.items()
            if category != "unknown"
        ),
        "remaining_unknown": transition_counts["unknown"],
        "newly_unknown_component_ids": newly_unknown,
        "components": outcomes,
        "limitation": (
            "The baseline cohort preserves the prior review's interpretations for "
            "comparison; those historical roles are not re-endorsed by current gates."
        ),
    }


def _component_record(
    component_id,
    indexes,
    members,
    representative,
    category,
    compounds,
):
    tolerance = representative.get("formula_tolerance") or {}
    measured_mz = float(representative.get("apex", representative["mz"]))
    nearest = (
        _nearest_formula(compounds, measured_mz, tolerance)
        if category == "unknown"
        else None
    )
    supporting_evidence = _supporting_evidence(members, category)
    roles = [
        (member.get("ion_role") or {}).get("kind")
        for member in members
        if member.get("ion_role")
    ]
    trace_statuses = sorted(
        {str(member.get("role_trace_status") or "not-reported") for member in members}
    )
    overlap_levels = sorted(
        {
            str((member.get("overlap") or {}).get("level"))
            for member in members
            if (member.get("overlap") or {}).get("level")
        }
    )
    additional = []
    if "unresolved" in trace_statuses or "unresolved-overlap" in roles:
        additional.append("inseparable_overlap")
    elif overlap_levels:
        additional.append("overlap_deconvolved_or_nearby")
    range_failures = _range_fit_failures(members)
    if range_failures:
        additional.append("one_or_more_interval_fits_withheld")
    background = representative.get("background_evidence")
    if background:
        additional.append("background_behaviour_" + background["status"])
    artefact_evidence = [
        member["artifact_evidence"]
        for member in members
        if member.get("artifact_evidence")
    ]
    if any(item.get("status") == "supporting" for item in artefact_evidence):
        additional.append("detector_echo_pattern_without_full_response_support")
    if measured_mz < 40.0:
        additional.append("low_mass_operational_or_fragment_region")

    resolved_reasons = {
        "direct_protonated_formula": "direct_formula_proposal_within_200_ppm",
        "alternative_ion_formula": "alternative_ion_hypothesis_supported",
        "inherited_parent_candidate": "parent_formula_inherited",
        "authored_assignment": "authored_assignment_preserved",
        "interpreted_noncompound": "noncompound_role_supported",
    }
    if category != "unknown":
        primary_reason = resolved_reasons[category]
        next_evidence = (
            "This support resolves the diagnostic category, not chemical identity. "
            "Formulae, ion pathways and database names remain reviewer evidence."
        )
    elif nearest is None:
        primary_reason = "no_plausible_formula_within_2000_ppm"
        next_evidence = (
            "A candidate needs a supported non-[M+H]+ pathway, a linked parent/fragment "
            "relationship, or an external standard; widening the mass gate is unsafe."
        )
    else:
        primary_reason = "nearest_formula_outside_200_ppm_proposal_gate"
        next_evidence = (
            "Independent recalibration or a standards run would be needed to move this "
            "formula inside the existing gate; the diagnostic match cannot be promoted."
        )
    if "inseparable_overlap" in additional:
        next_evidence += (
            " A higher-resolution measurement or independently constrained peak shape "
            "is also required for a separable trace."
        )

    return {
        "component_id": component_id,
        "peak_indexes": indexes,
        "category": category,
        "representative_mz": round(measured_mz, 6),
        "member_mz": [
            round(float(member.get("apex", member["mz"])), 6) for member in members
        ],
        "labels": [str(member.get("label") or "") for member in members],
        "primary_reason": primary_reason,
        "additional_constraints": additional,
        "nearest_formula_diagnostic": nearest,
        "supporting_evidence": supporting_evidence,
        "route_checks": {
            "direct_protonated_formula": (
                "available"
                if any(member.get("candidates") for member in members)
                else "none within proposal gate"
            ),
            "alternative_ion": (
                "none supported"
                if not any(member.get("ion_candidates") for member in members)
                else "available"
            ),
            "fragment_link": (
                "none supported"
                if not any(member.get("fragmentation_links") for member in members)
                else "available"
            ),
            "isotope_parent": (
                "none supported" if "isotope" not in roles else "available"
            ),
            "artefact_or_reagent": (
                "not supported" if not roles else sorted(set(roles))
            ),
        },
        "trace_statuses": trace_statuses,
        "overlap_levels": overlap_levels,
        "interval_fit_failures": range_failures,
        "background_evidence": background,
        "artefact_evidence": artefact_evidence,
        "next_evidence_needed": next_evidence,
    }


def _supporting_evidence(members, category):
    """Return compact evidence for the component's winning category."""
    if category == "direct_protonated_formula":
        return [
            {
                key: candidate[key]
                for key in (
                    "formula",
                    "delta_ppm",
                    "mass_match",
                    "assignment_eligible",
                )
                if key in candidate
            }
            for member in members
            for candidate in (member.get("candidates") or [])[:3]
        ]
    if category == "alternative_ion_formula":
        return [
            {
                key: candidate[key]
                for key in (
                    "formula",
                    "ion_notation",
                    "pathway",
                    "delta_ppm",
                    "assignment_eligible",
                )
                if key in candidate
            }
            for member in members
            for candidate in (member.get("ion_candidates") or [])[:3]
        ]
    if category == "inherited_parent_candidate":
        return [
            {
                "kind": item.get("kind"),
                "parent_mz": item.get("parent_mz"),
                "candidate_formulas": item.get("candidate_formulas"),
            }
            for member in members
            for item in member.get("interpretation_candidates") or []
            if item.get("candidate_formulas")
        ]
    return [
        {
            key: item[key]
            for key in ("kind", "label", "reason", "parent_mz")
            if key in item
        }
        for member in members
        for item in member.get("interpretation_candidates") or []
        if item.get("kind") and not str(item["kind"]).startswith("unresolved")
    ]


def _nearest_formula(compounds, mz, tolerance):
    candidates = compounds.score_peak(
        mz,
        tolerance_ppm=float(tolerance.get("ppm") or 10.0),
        proposal_tolerance_ppm=2000.0,
        max_candidates=1,
    )
    if not candidates:
        return None
    candidate = candidates[0]
    if abs(float(candidate["delta_ppm"])) <= float(
        tolerance.get("proposal_ppm") or 200.0
    ):
        return None
    return {
        "formula": candidate["formula"],
        "ion_mz": candidate["ion_mz"],
        "delta_ppm": candidate["delta_ppm"],
        "delta_mDa": candidate["delta_mDa"],
        "status": "diagnostic-only-outside-proposal-gate",
    }


def _range_fit_failures(members):
    failures = set()
    for member in members:
        for label, fit in ((member.get("fit") or {}).get("ranges") or {}).items():
            if fit.get("status") != "reliable":
                failures.add(f"{label}: {fit.get('reason') or fit.get('status')}")
    return sorted(failures)
