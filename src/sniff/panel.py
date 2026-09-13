"""Conservative adaptation of a curated peak table to a new PTR-MS file."""

from __future__ import annotations

import copy

_IDENTITY_FIELDS = (
    "label",
    "formula",
    "k",
    "k_estimated",
    "flags",
    "isotope",
)


def adapt_peak_table(template_peaks, detected_peaks, R_phys=2400.0):
    """Match curated identities to a comprehensive target-file detection panel.

    The target detections own every file-derived mass. A template contributes identity
    metadata only when exactly one nearest unused detection is within a conservative
    resolution-aware tolerance.
    """
    detected = [copy.deepcopy(peak) for peak in detected_peaks]
    unused = set(range(len(detected)))
    adapted = []
    matches = []
    missing = []

    for source in sorted(template_peaks or [], key=lambda peak: float(peak["mz"])):
        source_mass = float(source["mz"])
        # The corrected mass axis is calibrated by two internal references. Identity
        # transfer therefore uses mass accuracy, not peak width: a resolution-sized
        # window can contain a different compound when the source target is absent.
        tolerance = 0.015
        candidates = [
            index
            for index in unused
            if abs(float(detected[index]["mz"]) - source_mass) <= tolerance
        ]
        if not candidates:
            missing.append(
                {
                    "mz": source_mass,
                    "label": source.get("label") or source.get("formula") or "",
                    "reason": "no credible target-file peak within tolerance",
                }
            )
            continue
        candidates.sort(
            key=lambda index: abs(float(detected[index]["mz"]) - source_mass)
        )
        best = candidates[0]
        if len(candidates) > 1:
            first_error = abs(float(detected[candidates[0]]["mz"]) - source_mass)
            second_error = abs(float(detected[candidates[1]]["mz"]) - source_mass)
            if second_error - first_error < max(0.002, 0.2 * tolerance):
                missing.append(
                    {
                        "mz": source_mass,
                        "label": source.get("label") or source.get("formula") or "",
                        "reason": "ambiguous target-file match",
                    }
                )
                continue
        target = detected[best]
        # A transferred identity is authoritative as a unit. Do not combine a curated
        # label-only source with an unrelated automatic formula from the target run.
        for field in _IDENTITY_FIELDS:
            target.pop(field, None)
        for field in _IDENTITY_FIELDS:
            if field in source:
                target[field] = copy.deepcopy(source[field])
        adapted.append(target)
        unused.remove(best)
        matches.append(
            {
                "source_mz": source_mass,
                "target_mz": float(target["mz"]),
                "shift_mDa": round((float(target["mz"]) - source_mass) * 1000.0, 2),
            }
        )

    new_peaks = [detected[index] for index in sorted(unused)]
    used_formulas = {
        str(peak.get("formula") or "").strip().upper()
        for peak in adapted
        if peak.get("formula")
    }
    used_labels = {
        str(peak.get("label") or "").strip().casefold()
        for peak in adapted
        if peak.get("label")
    }
    duplicates_suppressed = []
    for peak in new_peaks:
        formula = str(peak.get("formula") or "").strip().upper()
        label = str(peak.get("label") or "").strip().casefold()
        if (formula and formula in used_formulas) or (label and label in used_labels):
            duplicates_suppressed.append(float(peak["mz"]))
            peak.pop("formula", None)
            peak["label"] = ""
            continue
        if formula:
            used_formulas.add(formula)
        if label:
            used_labels.add(label)
    adapted.extend(new_peaks)
    adapted.sort(key=lambda peak: float(peak["mz"]))
    return adapted, {
        "matched": matches,
        "missing": missing,
        "new_target_peaks": [float(peak["mz"]) for peak in new_peaks],
        "automatic_duplicates_suppressed": duplicates_suppressed,
        "n_matched": len(matches),
        "n_missing": len(missing),
        "n_new": len(new_peaks),
        "n_automatic_duplicates_suppressed": len(duplicates_suppressed),
    }
