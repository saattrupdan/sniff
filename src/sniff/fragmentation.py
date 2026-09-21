"""Conservative PTR fragmentation evidence from library pathways and run traces.

Fragmentation evidence is secondary to exact mass.  It can reorder formula proposals
that already exist, but it never creates a formula candidate or changes mass-based
assignment eligibility.
"""

import math

import numpy as np

MODEL = "ptr-library-covariation-v1"
MAX_EN_DIFFERENCE_TD = 20.0
SUPPORT_FACTOR = 1.35
MAX_VALIDATED_RESULTS = 5
MAX_BROAD_RESULTS = 2


def _text(value):
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace").strip()
    return str(value).strip()


def reaction_context(h5):
    """Return measured reagent and E/N context without inventing missing values."""
    context = {
        "model": MODEL,
        "source": "AddTraces/PTR-Reaction",
        "reagent": None,
        "primary_ion_index": None,
        "e_n_td": None,
        "e_n_range_td": None,
        "status": "unavailable",
        "reason": "PTR reaction metadata are unavailable",
    }
    info_path = "AddTraces/PTR-Reaction/Info"
    data_path = "AddTraces/PTR-Reaction/Data"
    if info_path not in h5 or data_path not in h5:
        return context
    info = np.asarray(h5[info_path][...])
    data = np.asarray(h5[data_path][...], dtype=np.float64)
    if info.ndim < 2 or not info.shape[0] or data.ndim != 2:
        return context
    names = [_text(value) for value in info[0]]

    def column(name):
        if name not in names:
            return None
        values = data[:, names.index(name)]
        values = values[np.isfinite(values)]
        return values if values.size else None

    e_n = column("E/N_Act")
    primary_index = column("PrimionIdx")
    e_n_stable = False
    if e_n is not None:
        context["e_n_td"] = round(float(np.median(e_n)), 3)
        context["e_n_range_td"] = [
            round(float(np.min(e_n)), 3),
            round(float(np.max(e_n)), 3),
        ]
        e_n_stable = float(np.max(e_n) - np.min(e_n)) <= 10.0
    primary_stable = False
    if primary_index is not None:
        rounded = np.rint(primary_index).astype(np.int64)
        values = np.unique(rounded)
        primary_stable = values.size == 1
        index = int(values[0]) if primary_stable else None
        context["primary_ion_index"] = index
        descriptions_path = "PTR-PrimaryIons/Descriptions"
        if descriptions_path in h5:
            descriptions = np.asarray(h5[descriptions_path][...]).reshape(-1)
            if index is not None and 0 <= index < descriptions.size:
                context["reagent"] = _text(descriptions[index]) or None
    if not primary_stable:
        context["reason"] = "the active primary ion changes during the run"
    elif context["reagent"] is None:
        context["reason"] = "the active primary ion could not be identified"
    elif context["e_n_td"] is None:
        context["reason"] = "the run E/N was not recorded"
    elif not e_n_stable:
        context["reason"] = "the run E/N varies by more than 10 Td"
    else:
        context["status"] = "available"
        context["reason"] = None
    return context


def _block_means(values, maximum=256):
    values = np.asarray(values, dtype=np.float64)
    if values.size <= maximum:
        return values
    block = int(math.ceil(values.size / float(maximum)))
    usable = values.size // block * block
    if usable < block:
        return values
    means = np.nanmean(values[:usable].reshape(-1, block), axis=1)
    if usable < values.size:
        means = np.concatenate((means, [np.nanmean(values[usable:])]))
    return means


def _correlation(left, right):
    left = np.asarray(left, dtype=np.float64)
    right = np.asarray(right, dtype=np.float64)
    if left.shape != right.shape:
        return None
    finite = np.isfinite(left) & np.isfinite(right) & (left >= 0) & (right >= 0)
    if int(np.count_nonzero(finite)) < 20:
        return None
    left = _block_means(left[finite])
    right = _block_means(right[finite])
    finite = np.isfinite(left) & np.isfinite(right)
    left = np.log1p(left[finite])
    right = np.log1p(right[finite])
    if left.size < 12:
        return None
    left_dynamic = float(np.quantile(left, 0.9) - np.quantile(left, 0.1))
    right_dynamic = float(np.quantile(right, 0.9) - np.quantile(right, 0.1))
    if left_dynamic < 0.15 or right_dynamic < 0.15:
        return None

    def coefficient(x, y):
        if x.size < 3 or float(np.std(x)) <= 1e-12 or float(np.std(y)) <= 1e-12:
            return None
        value = float(np.corrcoef(x, y)[0, 1])
        return value if math.isfinite(value) else None

    raw = coefficient(left, right)
    changes = coefficient(np.diff(left), np.diff(right))
    if raw is None or changes is None:
        return None
    supported = raw >= 0.80 and changes >= 0.35
    return {
        "n_blocks": int(left.size),
        "log_dynamic_range": [round(left_dynamic, 3), round(right_dynamic, 3)],
        "level_correlation": round(raw, 3),
        "change_correlation": round(changes, 3),
        "supported": supported,
    }


def _artifact_peak(peak):
    flags = " ".join(str(value).lower() for value in peak.get("likely_artifact", []))
    return any(word in flags for word in ("noise", "ringing", "saturation"))


def _nearest_peak(peaks, expected_mz, parent_index, r_phys):
    tolerance = max(0.015, min(0.05, float(expected_mz) / float(r_phys) * 0.75))
    options = []
    for index, peak in enumerate(peaks):
        if index == parent_index or _artifact_peak(peak):
            continue
        observed = float(peak.get("apex", peak["mz"]))
        error = abs(observed - float(expected_mz))
        if error <= tolerance:
            options.append((error, index, observed))
    return min(options) if options else None


def _profiles_by_formula(rate_table):
    return {
        item.get("formula"): list(item.get("fragmentation_profiles") or [])
        for item in (rate_table or {}).get("compounds", [])
        if item.get("formula")
    }


def _applicable_profiles(profiles, context):
    if context.get("status") != "available" or context.get("reagent") != "H3O+":
        return []
    run_e_n = float(context["e_n_td"])
    return [
        profile
        for profile in profiles
        if profile.get("reagent") == "H3O+"
        and profile.get("e_n_td") is not None
        and abs(float(profile["e_n_td"]) - run_e_n) <= MAX_EN_DIFFERENCE_TD
    ]


def apply_fragmentation_evidence(
    peaks,
    traces,
    rate_table,
    context,
    *,
    r_phys=2400.0,
    max_validated=MAX_VALIDATED_RESULTS,
    max_broad=MAX_BROAD_RESULTS,
):
    """Attach evidence, rerank existing candidates, and link supported fragments."""
    profiles_by_formula = _profiles_by_formula(rate_table)
    links = {}
    for parent_index, peak in enumerate(peaks):
        parent_mz = float(peak["mz"])
        parent_trace = traces.get(parent_mz)
        candidates = list(peak.get("candidates") or [])
        for candidate in candidates:
            profiles = profiles_by_formula.get(candidate.get("formula"), [])
            applicable = _applicable_profiles(profiles, context)
            matches = []
            for profile in applicable:
                for product in profile.get("products") or []:
                    if product.get("pathway_type", "fragment") != "fragment":
                        continue
                    expected_mz = product.get("mz")
                    if expected_mz is None:
                        continue
                    nearest = _nearest_peak(peaks, expected_mz, parent_index, r_phys)
                    if nearest is None or parent_trace is None:
                        continue
                    error, fragment_index, observed_mz = nearest
                    fragment_mz = float(peaks[fragment_index]["mz"])
                    fragment_trace = traces.get(fragment_mz)
                    if fragment_trace is None:
                        continue
                    correlation = _correlation(parent_trace, fragment_trace)
                    if correlation is None:
                        continue
                    match = {
                        "profile_id": profile.get("profile_id"),
                        "compound": profile.get("name"),
                        "expected_mz": round(float(expected_mz), 4),
                        "observed_mz": round(observed_mz, 4),
                        "error_mDa": round(error * 1000.0, 1),
                        "yield_percent": product.get("yield_percent"),
                        "library_e_n_td": profile.get("e_n_td"),
                        "reference": profile.get("reference"),
                        "doi": profile.get("doi"),
                        **correlation,
                    }
                    matches.append(match)
                    if correlation["supported"]:
                        links.setdefault(fragment_index, []).append(
                            {
                                "model": MODEL,
                                "parent_mz": round(
                                    float(peak.get("apex", peak["mz"])), 4
                                ),
                                "parent_formula": candidate.get("formula"),
                                "candidate_name": profile.get("name"),
                                "expected_fragment_mz": round(float(expected_mz), 4),
                                "observed_mz": round(observed_mz, 4),
                                "yield_percent": product.get("yield_percent"),
                                "level_correlation": correlation["level_correlation"],
                                "change_correlation": correlation["change_correlation"],
                                "reference": profile.get("reference"),
                                "doi": profile.get("doi"),
                            }
                        )
            supported = [match for match in matches if match["supported"]]
            if not profiles:
                status = "unavailable"
                reason = "no PTR Library fragmentation profile for this formula"
            elif not applicable:
                status = "unavailable"
                reason = "no profile matches the measured reagent and E/N context"
            elif supported:
                status = "support"
                reason = (
                    f"{len(supported)} library fragment channel(s) co-vary with the "
                    "candidate parent"
                )
            elif matches:
                status = "inconclusive"
                reason = "matched fragment channels do not show specific co-variation"
            else:
                status = "inconclusive"
                reason = "no suitable detected fragment trace could test this profile"
            candidate["fragmentation_evidence"] = {
                "model": MODEL,
                "status": status,
                "reason": reason,
                "run_reagent": context.get("reagent"),
                "run_e_n_td": context.get("e_n_td"),
                "applicable_profiles": len(applicable),
                "matches": sorted(
                    matches,
                    key=lambda item: (
                        not item["supported"],
                        -item["level_correlation"],
                        item["error_mDa"],
                    ),
                )[:5],
                "limitation": (
                    "co-variation supports a known pathway but is not MS/MS evidence "
                    "and may reflect shared sample timing"
                ),
            }
            candidate["fragmentation_factor"] = SUPPORT_FACTOR if supported else 1.0
            candidate["combined_score"] = float(
                candidate.get("probability", candidate.get("score", 0.0))
            ) * candidate["fragmentation_factor"]

        candidates.sort(
            key=lambda item: (
                bool(item.get("assignment_eligible", True)),
                float(item.get("combined_score", 0.0)),
                float(item.get("score", 0.0)),
            ),
            reverse=True,
        )
        eligible = [item for item in candidates if item.get("assignment_eligible", True)]
        broad = [item for item in candidates if not item.get("assignment_eligible", True)]
        retained = eligible[:max_validated]
        retained.extend(broad[:max_broad] if eligible else broad[:max_validated])
        total = sum(float(item.get("combined_score", 0.0)) for item in retained) or 1.0
        for item in retained:
            item["probability"] = round(float(item["combined_score"]) / total, 3)
            item["combined_score"] = round(float(item["combined_score"]), 6)
        peak["candidates"] = retained
        if retained:
            peak["id_confidence"] = retained[0]["probability"]
            peak["id_ambiguous"] = bool(
                retained[0]["probability"] < 0.6
                or (
                    len(retained) > 1
                    and retained[0]["probability"] - retained[1]["probability"] < 0.2
                )
            )
        else:
            peak.pop("id_confidence", None)
            peak.pop("id_ambiguous", None)

    for fragment_index, fragment_links in links.items():
        unique = {}
        for link in fragment_links:
            key = (link["parent_mz"], link["parent_formula"], link["expected_fragment_mz"])
            unique[key] = link
        peaks[fragment_index]["fragmentation_links"] = sorted(
            unique.values(),
            key=lambda item: (
                -item["change_correlation"],
                -item["level_correlation"],
                item["parent_mz"],
            ),
        )[:5]
    return peaks
