"""Run-wide detector-response evidence for PTR-MS ringing echoes."""

from __future__ import annotations

import numpy as np


def response_metrics(parent, child, *, n_folds=8):
    """Return blocked held-out proportional-response diagnostics for two traces.

    An interval-specific intercept absorbs unrelated constant background, while the
    non-negative gain must predict the echo's variation. Contiguous folds avoid
    treating adjacent cycles as independent evidence.
    """
    parent = np.asarray(parent, dtype=np.float64)
    child = np.asarray(child, dtype=np.float64)
    if parent.shape != child.shape or parent.ndim != 1:
        return _failed_metrics("trace shapes are incompatible")
    finite = np.isfinite(parent) & np.isfinite(child) & (parent >= 0) & (child >= 0)
    indexes = np.flatnonzero(finite)
    if indexes.size < 40:
        return _failed_metrics("fewer than 40 finite non-negative cycles")

    folds = [
        fold
        for fold in np.array_split(indexes, min(n_folds, indexes.size // 5))
        if len(fold)
    ]
    gains = []
    errors = []
    for held_out in folds:
        train = np.setdiff1d(indexes, held_out, assume_unique=True)
        parent_mean = float(np.mean(parent[train]))
        child_mean = float(np.mean(child[train]))
        parent_centred = parent[train] - parent_mean
        denominator = float(np.dot(parent_centred, parent_centred))
        if denominator <= 0 or not np.isfinite(denominator):
            continue
        gain = max(
            0.0,
            float(np.dot(parent_centred, child[train] - child_mean) / denominator),
        )
        intercept = child_mean - gain * parent_mean
        predicted = intercept + gain * parent[held_out]
        held_out_centre = float(np.median(child[held_out]))
        scale = max(
            float(np.linalg.norm(child[held_out] - held_out_centre)),
            0.05 * float(np.linalg.norm(child[held_out])),
            1e-12,
        )
        error = float(np.linalg.norm(child[held_out] - predicted) / scale)
        if np.isfinite(gain) and np.isfinite(error):
            gains.append(gain)
            errors.append(error)
    if len(gains) < 4:
        return _failed_metrics("fewer than four usable held-out folds")

    gains_array = np.asarray(gains)
    gain = float(np.median(gains_array))
    gain_mad = float(np.median(np.abs(gains_array - gain)))
    gain_relative_mad = gain_mad / max(gain, 1e-12)
    return {
        "usable": bool(gain > 0),
        "reason": None if gain > 0 else "non-positive detector response gain",
        "gain": gain,
        "gain_relative_mad": gain_relative_mad,
        "held_out_nrmse": float(np.median(errors)),
        "fold_count": len(gains),
    }


def response_families(pairs):
    """Return connected detector-family IDs for parent/child response pairs."""
    adjacency = {}
    for pair in pairs:
        parent = pair["parent_id"]
        child = pair["child_id"]
        adjacency.setdefault(parent, set()).add(child)
        adjacency.setdefault(child, set()).add(parent)
    families = {}
    for node in sorted(adjacency):
        if node in families:
            continue
        pending = [node]
        members = set()
        while pending:
            current = pending.pop()
            if current in members:
                continue
            members.add(current)
            pending.extend(adjacency.get(current, ()))
        family_id = min(members)
        for member in members:
            families[member] = family_id
    return families


def fit_response_model(reference_pairs):
    """Fit robust delay and response-gain distributions across independent families."""
    by_family = {}
    for pair in reference_pairs:
        metrics = pair["metrics"]
        if (
            not metrics.get("usable")
            or metrics["held_out_nrmse"] > 0.35
            or metrics["gain_relative_mad"] > 0.30
        ):
            continue
        family_id = pair.get("family_id", pair["parent_id"])
        previous = by_family.get(family_id)
        if (
            previous is None
            or pair["metrics"]["held_out_nrmse"] < previous["metrics"]["held_out_nrmse"]
        ):
            by_family[family_id] = pair
    usable = list(by_family.values())
    family_ids = set(by_family)
    if len(family_ids) < 3:
        return _failed_model("fewer than three independent response parents")

    delays = np.asarray([pair["delay"] for pair in usable], dtype=np.float64)
    gains = np.asarray([pair["metrics"]["gain"] for pair in usable], dtype=np.float64)
    errors = np.asarray(
        [pair["metrics"]["held_out_nrmse"] for pair in usable], dtype=np.float64
    )
    delay = float(np.median(delays))
    delay_mad = float(np.median(np.abs(delays - delay)))
    gain = float(np.median(gains))
    gain_mad = float(np.median(np.abs(gains - gain)))
    gain_relative_mad = gain_mad / max(gain, 1e-12)
    # A finite floor avoids treating a tiny synthetic MAD as exact detector physics.
    half_width = max(3.0 * gain_mad, 0.25 * gain)
    gain_interval = [max(0.0, gain - half_width), gain + half_width]
    held_out_nrmse = float(np.median(errors))
    reliable = bool(
        delay_mad <= 6.0
        and gain > 0
        and gain_relative_mad <= 0.30
        and held_out_nrmse <= 0.35
    )
    failed = []
    if delay_mad > 6.0:
        failed.append("delay dispersion exceeds 6 timebins")
    if gain <= 0:
        failed.append("response gain is non-positive")
    if gain_relative_mad > 0.30:
        failed.append("response gain varies by more than 30%")
    if held_out_nrmse > 0.35:
        failed.append("reference held-out error exceeds 0.35")
    return {
        "usable": True,
        "reliable": reliable,
        "reason": None if reliable else "; ".join(failed),
        "parent_count": len(family_ids),
        "pair_count": len(usable),
        "delay_timebins": delay,
        "delay_mad_timebins": delay_mad,
        "gain": gain,
        "gain_interval": gain_interval,
        "gain_relative_mad": gain_relative_mad,
        "held_out_nrmse": held_out_nrmse,
    }


def score_response(metrics, model):
    """Apply fixed conservative gates to one candidate response."""
    failed = []
    if not model.get("reliable"):
        failed.append("run-wide response model is not reliable")
    if not metrics.get("usable"):
        failed.append(metrics.get("reason") or "candidate response is unavailable")
    else:
        if metrics["held_out_nrmse"] > 0.35:
            failed.append("candidate held-out error exceeds 0.35")
        if metrics["gain_relative_mad"] > 0.30:
            failed.append("candidate gain varies by more than 30% across folds")
        interval = model.get("gain_interval")
        low, high = interval or [float("inf"), -float("inf")]
        if not low <= metrics["gain"] <= high:
            failed.append("candidate gain is outside the run-wide prediction interval")
    return {"supported": not failed, "failed_gates": failed}


def _failed_metrics(reason):
    return {
        "usable": False,
        "reason": reason,
        "gain": None,
        "gain_relative_mad": None,
        "held_out_nrmse": None,
        "fold_count": 0,
    }


def _failed_model(reason):
    return {
        "usable": False,
        "reliable": False,
        "reason": reason,
        "parent_count": 0,
        "pair_count": 0,
        "delay_timebins": None,
        "delay_mad_timebins": None,
        "gain": None,
        "gain_interval": None,
        "gain_relative_mad": None,
        "held_out_nrmse": None,
    }
