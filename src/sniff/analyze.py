#!/usr/bin/env python3
"""PTR-MS analysis CLI — open-source reprocessor for IONICON IoniTOF HDF5 files.

Designed to be driven by an agent, not a human. From a checkout, invoke it with
`uv run sniff <subcommand>` — do not
read this source; every operation is a subcommand and every value is in its JSON
output.

Commands (all discovery output is JSON on stdout):
  inspect    File metadata, calibration, transmission, concentration-K, molar volume.
  peaks      Peak detection -> {mz, height, neutral_mass, suggested_label,
             top_candidate, [likely_artifact]} (compact; --full for all candidates).
  segments   Time-segment detection -> stable plateaus (samples vs background).
  analyze    Full pipeline (from a config) -> PTR-MS-Viewer-style results CSV.
  viz        Review app for an existing peak list + ranges (live-save to a config
             with --serve, or a standalone HTML with --html). Does NOT detect.
  calibrate  Fit the concentration constant K to a reference Viewer CSV.
  compare    Error stats of a results CSV vs a reference Viewer CSV.
  rates      Browse the bundled proton-transfer rate-constant table.

There is deliberately NO one-shot command. The agent detects with `peaks`/`segments`,
applies its own chemistry + curation judgment, writes a config, and only then reviews or
analyses it — so `viz`/`analyze` always operate on the best solution, not a mechanical guess.

Flows:
  primary — agent curates, then an expert confirms in the browser:
    1. sniff peaks FILE        -> pick assignments from each peak's `candidates`
    2. sniff segments FILE     -> sample_01/background_01 ranges; never ask names
    3. write analysis-config.json (the curated peaks + ranges)
    4. sniff viz FILE --config analysis-config.json --out results.csv
       (serves the browser app; clicking 'Done' writes results.csv itself)
  no review (the same curated config, straight to CSV):
       sniff analyze FILE --config analysis-config.json --include-cycle-rows --out results.csv
  quick deterministic fallback (no agent judgment, no browser — detect + quantify only):
       sniff analyze FILE --auto-peaks --auto-segments --include-cycle-rows --out results.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import secrets
import sys
from contextlib import ExitStack

import h5py
import numpy as np

from . import catalogue, formula_id, fragmentation, ptrms

_ANALYSIS_DEFAULTS = {
    "R": 1200.0,
    "R_phys": 2400.0,
    "K": None,
    "molar_volume": None,
    "primary_mz": 21.022,
    "kinetic": False,
    "k_anchor": ptrms.K_ANCHOR_DEFAULT,
    "humidity_correct": False,
    "humidity_p": 1.0,
    "humidity_ref": None,
    "whole_run_windows": False,
    "peak_fit": "gaussian-v1",
    "isotope_mode": "off",
    "isotope_abundance_basis": "unknown",
}

_X_AXIS_UNITS = ("cycle", "relative", "absolute")


def resolve_x_axis_unit(config=None, args=None):
    """Resolve the viz x-axis unit with CLI > config > cycle precedence."""
    config = config or {}
    curated = config.get("viz")
    if curated is not None and not isinstance(curated, dict):
        raise ValueError("viz must be a JSON object")
    curated = curated or {}
    cli_value = getattr(args, "x_axis_unit", None) if args is not None else None
    value = cli_value if cli_value is not None else curated.get("x_axis_unit", "cycle")
    if value not in _X_AXIS_UNITS:
        raise ValueError("viz.x_axis_unit must be one of: " + ", ".join(_X_AXIS_UNITS))
    return value


def _load_config(args):
    """Load the complete config once, retaining fields unknown to this CLI."""
    if not getattr(args, "config", None):
        return {}
    with open(args.config, encoding="utf-8") as fh:
        return json.load(fh)


def _migrate_loaded_config(config, args, mass_axis):
    """Migrate and scrub a file-backed config after successful calibration."""
    if not getattr(args, "config", None) or not config:
        return config
    migrated, changed = ptrms.migrate_config_mass_axis(config, mass_axis)
    migrated, retired_fields_removed = _scrub_retired_review_fields(migrated)
    if changed or retired_fields_removed:
        with open(args.config, "w", encoding="utf-8") as fh:
            json.dump(migrated, fh, indent=2)
    return migrated


def _scrub_retired_review_fields(config):
    """Remove retired checklist fields without discarding unknown config data."""
    cleaned = dict(config)
    changed = "checklist" in cleaned
    cleaned.pop("checklist", None)
    review = cleaned.get("review")
    if isinstance(review, dict) and "checklist" in review:
        review = dict(review)
        review.pop("checklist", None)
        changed = True
        if review:
            cleaned["review"] = review
        else:
            cleaned.pop("review", None)
    return cleaned, changed


def resolve_analysis_settings(config=None, args=None):
    """Resolve analysis settings with CLI > curated config > legacy defaults.

    ``None`` is the argparse sentinel for an omitted option.  Keeping this
    distinction here prevents a parser default from silently replacing a
    curated value when a config is used by either command.
    """
    config = config or {}
    curated = config.get("analyze") or {}
    settings = {}
    sources = {}
    for key, default in _ANALYSIS_DEFAULTS.items():
        cli_value = None
        if args is not None:
            cli_value = getattr(args, key, None)
            if key == "whole_run_windows":
                no_per = getattr(args, "no_per_interval", None)
                if no_per is not None:
                    cli_value = bool(no_per)
        if cli_value is not None:
            settings[key] = cli_value
            sources[key] = "cli"
        elif key in curated:
            settings[key] = curated[key]
            sources[key] = "config.analyze"
        else:
            settings[key] = default
            sources[key] = "legacy default"
    if settings["peak_fit"] not in ("gaussian-v1", "empirical-v1"):
        raise ValueError("peak_fit must be gaussian-v1 or empirical-v1")
    if settings["isotope_mode"] not in ("off", "formula-v1"):
        raise ValueError("isotope_mode must be off or formula-v1")
    if settings["isotope_abundance_basis"] not in (
        "unknown",
        "calibrated",
        "total",
    ):
        raise ValueError(
            "isotope_abundance_basis must be unknown, calibrated, or total"
        )
    settings["sources"] = sources
    settings["per_interval_windows"] = not bool(settings["whole_run_windows"])
    return settings


def _effective_sources(settings, humidity_ref, molar_volume_source=None):
    """Describe the runtime provenance of derived and configured values."""
    sources = dict(settings["sources"])
    if settings["K"] is None:
        sources["K"] = "file acquisition calibration"
    if settings["molar_volume"] is None:
        sources["molar_volume"] = molar_volume_source or "file drift temperature"
    if settings["humidity_ref"] is None:
        sources["humidity_ref"] = (
            "run median" if humidity_ref is not None else "unavailable"
        )
    return sources


def _humidity_proxy_label(primary_mz):
    """Describe the humidity ratio with the effective denominator."""
    return f"m/z 37 / m/z {primary_mz:g} (water-cluster ratio; m/z 19 saturates)"


def _peak_windows(peaks):
    """Return curated asymmetric half-widths for peaks that specify windows."""

    def _winlr(peak):
        window = peak["window"]
        if isinstance(window, dict):
            return float(window["left"]), float(window["right"])
        width = float(window) / 2.0
        return width, width

    return {float(p["mz"]): _winlr(p) for p in peaks if p.get("window")}


def _emit(obj, raw=True):
    json.dump(obj, sys.stdout, indent=None if raw else 2)
    sys.stdout.write("\n")


def detect_peaks(
    f,
    min_rel_height=1e-3,
    max_peaks=300,
    mz_min=15.0,
    mz_max=None,
    R_phys=2400.0,
    noise_sigma=6.0,
    mass_axis=None,
):
    """Untargeted peak detection on the average spectrum.

    Local maxima above a height threshold, then merged if closer than one
    instrument linewidth (FWHM = m/R_phys), keeping the tallest — two maxima closer
    than the resolution are one peak (ripples on a shared apex), never two. Each kept
    peak also gets a `prominence`: how far its apex rises above the local baseline
    within ~1 linewidth each side. A real peak rises from near-zero baseline, so its
    prominence ≈ its height; a ripple or shoulder riding a taller peak's flank barely
    dips, so its prominence is tiny — this lets `annotate_peaks` flag noise combs
    (dozens of ~20 cps satellites around an intense ion) without deleting anything.

    The threshold is the LARGER of a relative floor (amax*min_rel_height) and an
    absolute noise floor (median + noise_sigma * robust_sigma). The relative floor
    dominates on real spectra; the noise floor stops a low-count/blank file — where
    amax itself is noise — from returning hundreds of spurious maxima. Non-finite
    bins (rare file corruption) are treated as zero rather than poisoning amax."""
    if mass_axis is None:
        mass_axis = ptrms.load_mass_axis(f)
    else:
        ptrms.validate_mass_axis(mass_axis)
    a, b = mass_axis.a, mass_axis.b
    avg = np.asarray(f["SPECdata/AverageSpec"][:], dtype=np.float64)
    avg = np.where(np.isfinite(avg), avg, 0.0)
    amax = float(avg.max())
    med = float(np.median(avg))
    mad = float(np.median(np.abs(avg - med)))
    sigma = 1.4826 * mad if mad > 0 else (float(avg.std()) or 1e-9)
    thr = max(amax * min_rel_height, med + noise_sigma * sigma)
    hi = (avg[1:-1] > avg[:-2]) & (avg[1:-1] >= avg[2:]) & (avg[1:-1] > thr)
    idx = np.where(hi)[0] + 1
    mz = ptrms.tb_to_m(idx, a, b, mass_axis)
    keep = mz >= mz_min
    if mz_max:
        keep &= mz <= mz_max
    idx, mz = idx[keep], mz[keep]
    # merge maxima within one FWHM (m/R_phys), keeping the tallest
    order_m = np.argsort(mz)
    idx, mz = idx[order_m], mz[order_m]
    merged = []  # (mz, height, bin)
    for k in range(len(mz)):
        h = float(avg[idx[k]])
        if merged and mz[k] - merged[-1][0] < mz[k] / R_phys:
            if h > merged[-1][1]:
                merged[-1] = (float(mz[k]), h, int(idx[k]))
        else:
            merged.append((float(mz[k]), h, int(idx[k])))
    merged.sort(key=lambda t: t[1], reverse=True)
    merged = merged[:max_peaks]

    def prominence(i):
        # half-window ≈ 1 FWHM in timebins: d(bin)/d(mz)=a/(2√mz), FWHM_mz=mz/R_phys
        m = ptrms.tb_to_m(i, a, b, mass_axis)
        w = max(2, round(a * np.sqrt(m) / (2 * R_phys)))
        lo, hiw = max(0, i - w), min(len(avg), i + w + 1)
        base = max(float(avg[lo : i + 1].min()), float(avg[i:hiw].min()))
        return float(avg[i]) - base

    peaks = [
        {
            "mz": round(m, 4),
            "height": round(h, 1),
            "rel_height": round(h / amax, 5) if amax > 0 else 0.0,
            "prominence": round(prominence(ib), 1),
        }
        for m, h, ib in merged
    ]
    peaks.sort(key=lambda p: p["mz"])
    return peaks


# reagent (primary) ions: at least one dominates every real PTR spectrum —
# H3O+ isotope at m/z 21 (H3O+ mode), or NO+/O2+ in switched-reagent modes.
_PRIMARY_MZ = (21.022, 30.994, 31.989, 33.994, 37.033)


def assess_signal(f, avg=None, a=None, b=None, mass_axis=None):
    """Judge whether a file holds real measurement signal or is a blank capture.

    Every real PTR run is dominated by its reagent (primary) ion, which towers over
    the spectral noise. A blank / no-beam / aborted acquisition has no ion beam, so
    even the primary ion sits at the noise floor and every apparent 'peak' is just
    Poisson noise — detecting analytes from it is meaningless. Robustly compares the
    strongest primary-ion region against the spectrum's median + robust sigma.
    Returns dict(signal_present, primary_snr, reason)."""
    if avg is None:
        avg = f["SPECdata/AverageSpec"][:]
    avg = np.asarray(avg, dtype=np.float64)
    finite = avg[np.isfinite(avg)]
    if finite.size == 0:
        return {
            "signal_present": False,
            "primary_snr": 0.0,
            "reason": "average spectrum is entirely non-finite (corrupt file)",
        }
    med = float(np.median(finite))
    mad = float(np.median(np.abs(finite - med)))
    sigma = 1.4826 * mad if mad > 0 else (float(finite.std()) or 1e-9)
    if mass_axis is None:
        mass_axis = ptrms.load_mass_axis(f)
    if mass_axis is not None:
        ptrms.validate_mass_axis(mass_axis)
    if a is None:
        a, b = mass_axis.a, mass_axis.b

    def local_max(mz):
        tb = int(ptrms.m_to_tb(mz, a, b, mass_axis))
        w = avg[max(0, tb - 30) : tb + 30]
        w = w[np.isfinite(w)]
        return float(w.max()) if w.size else med

    primary_snr = max((local_max(mz) - med) / sigma for mz in _PRIMARY_MZ)
    present = primary_snr >= 20.0
    reason = (
        ""
        if present
        else (
            "no reagent (primary) ion detectable above the spectral noise "
            f"(primary-ion S/N {primary_snr:.1f} < 20) — this file appears to be a blank / no-beam / "
            "aborted acquisition, not a measurement, so no analyte peaks can be "
            "extracted from it."
        )
    )
    return {
        "signal_present": present,
        "primary_snr": round(primary_snr, 1),
        "reason": reason,
    }


# ----------------------------- commands -----------------------------
def cmd_inspect(args):
    with h5py.File(args.h5, "r") as f:
        mass_axis = ptrms.load_mass_axis(f)
        a, b = mass_axis.a, mass_axis.b
        tm, tf = ptrms.load_transmission(f)
        ncyc = int(f["SPECdata/Intensities"].shape[0])
        dur = ptrms.spec_duration_s(f)
        created = ""
        try:
            created = f.attrs["FileCreatedTimeSTR_LOCAL"][0].decode("latin-1")
        except (AttributeError, IndexError, KeyError, OSError, TypeError):
            created = ""
        _emit(
            {
                "file": args.h5,
                "instrument": _attr(f, "InstrumentType"),
                "created_local": created,
                "n_cycles": ncyc,
                "cycle_duration_s": dur,
                "duration_min": round(ncyc * dur / 60, 1),
                "n_spectrum_bins": int(f["SPECdata/Intensities"].shape[1]),
                "mass_cal": {
                    "model": "timebin = a*sqrt(m_file) + b",
                    "a": a,
                    "b": b,
                },
                "mass_axis_calibration": mass_axis.to_dict(),
                "fragmentation_context": fragmentation.reaction_context(f),
                "transmission_available": ptrms.has_transmission(f),
                "transmission_masses": [round(x, 3) for x in tm.tolist()],
                "transmission_factors": [round(x, 4) for x in tf.tolist()],
                "concentration_K_from_file": ptrms.derive_K(
                    f, ptrms.extract_primary(f, mass_axis=mass_axis)
                ),
                "molar_volume_L_per_mol": round(
                    ptrms.derive_molar_volume_info(f)[0], 3
                ),
                "molar_volume_source": ptrms.derive_molar_volume_info(f)[1],
                "has_precomputed_traces": "TRACEdata/TraceConcentration" in f,
            },
            args.raw,
        )


def _attr(f, key):
    try:
        v = f.attrs[key]
        v = v[0] if hasattr(v, "__len__") and not isinstance(v, (bytes, str)) else v
        return v.decode("latin-1") if isinstance(v, bytes) else v
    except (AttributeError, IndexError, KeyError, OSError, TypeError, ValueError):
        return None


_REAGENT_MZ = {
    19.018: "H3O+ primary",
    21.022: "H3O+ (18O) isotope",
    37.033: "H3O+·H2O cluster (operational calibration water)",
    55.039: "H3O+·(H2O)2 cluster",
    73.049: "H3O+·(H2O)3 cluster",
    31.989: "O2+",
    32.997: "O2+ (17O)",
    33.994: "O2+ (18O)",
    29.997: "NO+",
    30.994: "NO+ (15N) isotope",
}
_REAGENT_TOL_DA = 0.012


def annotate_peaks(
    peaks,
    avgspec=None,
    a=None,
    b=None,
    R=1200.0,
    R_phys=2400.0,
    elements=None,
    mass_axis=None,
    compounds_of_interest=None,
    assign_all_library=False,
    candidate_pool_size=5,
):
    """Enrich detected peaks with candidate FORMULA assignments (scored by mass +
    isotope pattern + plausibility) and artifact flags, so the agent/expert can
    pick assignments without scripting mass-matching.

    Candidates come from offline formula enumeration (`formula_id`), not a fixed
    short list, so near-isobars are disambiguated by their measured 13C(M+1) and
    heteroatom(M+2, e.g. S/Cl) isotope ratios rather than by "nearest mass". When
    `avgspec`+`a`+`b` are supplied the isotope ratios are measured from the average
    spectrum; without them the ranking falls back to mass + plausibility only.

    Returns (drift, annotated_peaks). On an accepted affine mass axis, `drift` is
    exactly 1 because the correction has already moved every peak. Each candidate's
    `delta_mDa` is the exact-mass residual after that handling, plus
    predicted/observed isotope ratios and a normalised candidate score/share
    (`probability`). It is not a
    calibrated identification probability. Automatic names and formulas are editable
    best guesses when a review workflow requests maximum-coverage defaults."""
    if mass_axis is not None:
        ptrms.validate_mass_axis(mass_axis)
    tbl = ptrms.load_rate_constants()
    comps = tbl["compounds"] if tbl else []
    ratios = []
    for p in peaks:
        near = [c for c in comps if abs(c["mz"] - p["mz"]) < 0.08]
        if len(near) == 1:
            ratios.append(p["mz"] / near[0]["mz"])
    drift = (
        1.0
        if mass_axis is not None and mass_axis.applied
        else (float(np.median(ratios)) if ratios else 1.0)
    )

    have_spec = avgspec is not None and a is not None and b is not None
    compound_catalogue = catalogue.CompoundCatalogue()

    def obs_ratios(mz):
        if not have_spec:
            return None

        def wsum(center):
            wl, wr = ptrms.peak_window(center, a, b, R, mass_axis)
            lo, hi = max(0, wl), min(len(avgspec), wr)
            return float(avgspec[lo:hi].sum()) if hi > lo else 0.0

        i0 = wsum(mz)
        if i0 <= 0:
            return None
        return (wsum(mz + formula_id.DM1) / i0, wsum(mz + formula_id.DM2) / i0)

    all_mz = sorted(q["mz"] for q in peaks)
    reagent_labels = {}
    for reagent_mz, reagent_name in _REAGENT_MZ.items():
        if not peaks:
            break
        index, nearest = min(
            enumerate(peaks),
            key=lambda item: abs(item[1]["mz"] - reagent_mz * drift),
        )
        difference = abs(nearest["mz"] - reagent_mz * drift)
        if difference <= _REAGENT_TOL_DA:
            reagent_labels[index] = reagent_name

    def nearest_other(mz):
        best = None
        for x in all_mz:
            if x == mz:
                continue
            if best is None or abs(x - mz) < abs(best - mz):
                best = x
        return best

    out = []
    for peak_index, p in enumerate(peaks):
        mz, h = p["mz"], p.get("height", 0.0)
        e = dict(p)
        e["neutral_mass"] = round(mz - ptrms.PROTON, 4)
        observed_isotopes = obs_ratios(mz)
        tolerance = (
            ptrms.formula_assignment_tolerance(mass_axis, mz)
            if mass_axis is not None
            else {
                "model": "no-run-calibration-fallback",
                "source": "no run calibration was supplied",
                "status": "fallback",
                "reason": "independent calibration residuals are unavailable",
                "ppm": 10.0,
                "mDa": 10.0 * mz / 1000.0,
                "proposal_ppm": ptrms.FORMULA_PROPOSAL_TOLERANCE_PPM,
                "proposal_mDa": (
                    ptrms.FORMULA_PROPOSAL_TOLERANCE_PPM * mz / 1000.0
                ),
                "score_sigma_ppm": 4.0,
                "candidate_generation_allowed": True,
                "automatic_assignment_allowed": False,
            }
        )
        e["formula_tolerance"] = {
            **tolerance,
            "ppm": round(tolerance["ppm"], 3),
            "mDa": round(tolerance["mDa"], 3),
            "proposal_ppm": round(tolerance["proposal_ppm"], 3),
            "proposal_mDa": round(tolerance["proposal_mDa"], 3),
            "score_sigma_ppm": round(tolerance["score_sigma_ppm"], 3),
        }
        local_candidates = (
            formula_id.score_peak(
                mz,
                drift,
                obs_ratios=observed_isotopes,
                elements=elements,
                compounds_of_interest=compounds_of_interest,
                tolerance_ppm=tolerance["ppm"],
                mass_sigma_ppm=tolerance["score_sigma_ppm"],
                proposal_tolerance_ppm=tolerance["proposal_ppm"],
                max_candidates=candidate_pool_size,
            )
            if tolerance["candidate_generation_allowed"]
            else []
        )
        cands = (
            compound_catalogue.score_peak(
                mz,
                drift=drift,
                candidates=local_candidates,
                obs_ratios=observed_isotopes,
                compounds_of_interest=compounds_of_interest,
                elements=elements,
                tolerance_ppm=tolerance["ppm"],
                mass_sigma_ppm=tolerance["score_sigma_ppm"],
                proposal_tolerance_ppm=tolerance["proposal_ppm"],
                max_candidates=candidate_pool_size,
            )
            if tolerance["candidate_generation_allowed"]
            else []
        )
        if not tolerance["automatic_assignment_allowed"]:
            for candidate in cands:
                candidate["assignment_eligible"] = False
                candidate["mass_match"] = "calibration-unvalidated-proposal"
        e["candidates"] = cands
        # normalized top-candidate score / near-isobar ambiguity, surfaced explicitly;
        # conservative assignment gates below require multiple candidates
        if cands:
            e["id_confidence"] = cands[0]["probability"]
            top2 = (
                len(cands) > 1
                and cands[0]["probability"] - cands[1]["probability"] < 0.2
            )
            if cands[0]["probability"] < 0.6 or top2:
                e["id_ambiguous"] = [
                    {
                        "formula": c["formula"],
                        "name": c["name"],
                        "probability": c["probability"],
                    }
                    for c in cands[:3]
                    if c["probability"] >= 0.05
                ]
        # spectral overlap with a neighbouring peak (affects quantification)
        nb = nearest_other(mz)
        if nb is not None:
            sep = abs(nb - mz)
            if sep < mz / R_phys * 1.5:  # within ~1.5 physical FWHM
                e["overlap"] = {
                    "neighbor": round(nb, 4),
                    "sep_mDa": round(sep * 1000, 1),
                    "level": "unresolved",
                    "note": "closer than the instrument resolution — "
                    "Raw is unreliable even after deconvolution",
                }
            elif sep < 0.20:
                e["overlap"] = {
                    "neighbor": round(nb, 4),
                    "sep_mDa": round(sep * 1000, 1),
                    "level": "deconvolved",
                    "note": "overlaps a neighbour; the configured analysis model "
                    "must establish whether independent deconvolution is reliable",
                }
        flags = []
        if peak_index in reagent_labels:
            reagent_name = reagent_labels[peak_index]
            isobaric_water_cluster = reagent_name.startswith("H3O+·(")
            if not (isobaric_water_cluster and cands):
                flags.append("reagent/cluster: " + reagent_name)
        for q in peaks:
            if 0.008 < mz - q["mz"] < 0.4 and q.get("height", 0) > 20 * max(h, 1):
                flags.append(f"possible tail/ringing of taller m/z {q['mz']:.3f}")
                break
        # low-prominence noise: an apex that barely rises above its local baseline
        # is a ripple/shoulder on a taller peak's flank, not a resolved peak. Catches
        # the dense combs of ~20 cps satellites around intense ions (e.g. TOF
        # ringing) that the tail/ringing test misses when the tall parent is >0.4 Da
        # away. Requires BOTH a small absolute rise and a small fraction of its own
        # height, so a genuinely isolated small peak (which rises from ~0) is kept.
        prom = p.get("prominence")
        if prom is not None and prom < 10.0 and prom < 0.2 * max(h, 1):
            flags.append(
                f"low prominence ({prom:.1f} cps above local baseline) — likely a "
                "noise ripple / shoulder of a nearby taller peak"
            )
        # H3O+ reagent saturation skirt: the primary ion at m/z ~19 saturates the
        # detector (the run normally normalises on its configured primary
        # isotope), and its
        # ringing throws a skirt of strong, sharp satellites between the primary and
        # that isotope. No real H3O+-chemistry analyte lives at m/z ~19.05–20.95
        # (ammonia at 18.03 sits BELOW the primary and is untouched), so peaks there
        # — which are high-prominence and thus escape the ripple test — are flagged
        # as reagent-region artifacts, not analytes.
        if 19.05 < mz < 20.95 and not any(
            fl.startswith("reagent/cluster") for fl in flags
        ):
            flags.append(
                "H3O+ primary saturation region (m/z 19–21) — detector "
                "ringing of the saturated reagent ion, not an analyte"
            )
        if flags:
            e["likely_artifact"] = flags
        out.append(e)
    interpret_peak_roles(out, drift=drift, R_phys=R_phys)
    _assign_suggested_identities(out, assign_all_library=assign_all_library)
    return drift, out


def apply_run_fragmentation_evidence(
    f,
    peaks,
    *,
    mass_axis=None,
    R=1200.0,
    R_phys=2400.0,
    drift=1.0,
    progress=None,
    should_stop=None,
):
    """Rerank existing formula candidates from measured PTR fragment co-variation."""
    context = fragmentation.reaction_context(f)
    masses = [float(peak["mz"]) for peak in peaks]
    if not masses or "SPECdata/Intensities" not in f:
        return context
    traces, _ = ptrms.extract_traces(
        f,
        masses,
        R=R,
        R_phys=R_phys,
        mass_axis=mass_axis,
        progress=progress,
        should_stop=should_stop,
    )
    fragmentation.apply_fragmentation_evidence(
        peaks,
        {mass: traces[mass][0] for mass in masses},
        ptrms.load_rate_constants(),
        context,
        r_phys=R_phys,
    )
    interpret_peak_roles(peaks, drift=drift, R_phys=R_phys)
    return context


def interpret_peak_roles(peaks, *, drift=1.0, R_phys=2400.0):
    """Attach evidence-backed non-analyte or unresolved interpretations.

    Formula candidates describe protonated neutral analytes. This separate list covers
    peaks for which assigning a neutral compound would be misleading: known reagent
    ions, detector artefacts, likely isotope channels, authored identities outside the
    current candidate window, and unresolved real ions. A possible isotope requires a
    stronger measured parent at the expected spacing; charge states and fragments are
    intentionally not inferred from mass alone.
    """
    if not peaks:
        return peaks

    def measured_mz(peak):
        return float(peak.get("apex", peak["mz"])) / float(drift)

    def abundance(peak):
        value = peak.get("abundance", peak.get("height", 0.0))
        return 0.0 if value is None else max(0.0, float(value))

    masses = [measured_mz(peak) for peak in peaks]
    heights = [abundance(peak) for peak in peaks]
    reagent_by_index = {}
    for marker_mz, marker_name in _REAGENT_MZ.items():
        nearest_index = min(
            range(len(peaks)), key=lambda candidate: abs(masses[candidate] - marker_mz)
        )
        if abs(masses[nearest_index] - marker_mz) <= _REAGENT_TOL_DA:
            reagent_by_index[nearest_index] = marker_name
    for index, peak in enumerate(peaks):
        interpretations = []
        validated_candidates = [
            candidate
            for candidate in peak.get("candidates") or []
            if candidate.get("assignment_eligible", True)
        ]
        flags = list(peak.get("likely_artifact") or [])
        reagent = next(
            (
                flag.split(": ", 1)[1]
                for flag in flags
                if flag.startswith("reagent/cluster: ")
            ),
            None,
        )
        if reagent is None:
            reagent = reagent_by_index.get(index)
        isobaric_water_cluster = bool(reagent and reagent.startswith("H3O+·("))
        if isobaric_water_cluster and validated_candidates:
            reagent = None
        if reagent:
            interpretations.append(
                {
                    "kind": "reagent",
                    "label": reagent,
                    "source": "known reagent-ion exact mass",
                    "evidence": [
                        f"nearest reagent marker within {_REAGENT_TOL_DA * 1000:.0f} mDa"
                    ],
                    "exclude_from_analyte_assignment": True,
                }
            )
        noise_flags = [flag for flag in flags if _is_noise_artifact([flag])]
        if noise_flags:
            interpretations.append(
                {
                    "kind": "artifact",
                    "label": "detector artefact / noise",
                    "source": "peak-shape and reagent-region diagnostics",
                    "evidence": noise_flags,
                    "exclude_from_analyte_assignment": True,
                }
            )

        if (
            not interpretations
            and not validated_candidates
            and peak.get("fragmentation_links")
        ):
            link = peak["fragmentation_links"][0]
            identity = link.get("candidate_name") or link.get("parent_formula")
            interpretations.append(
                {
                    "kind": "fragment",
                    "label": (
                        f"possible fragment of {identity} at parent m/z "
                        f"{float(link['parent_mz']):.4f}"
                    ),
                    "source": "PTR Library pathway and measured temporal co-variation",
                    "related_mz": link["parent_mz"],
                    "evidence": [
                        f"library fragment m/z {float(link['expected_fragment_mz']):.4f}",
                        (
                            "level/change correlations "
                            f"{float(link['level_correlation']):.2f}/"
                            f"{float(link['change_correlation']):.2f}"
                        ),
                        "supporting evidence only; not MS/MS proof",
                    ],
                    "exclude_from_analyte_assignment": True,
                }
            )

        if not interpretations and not validated_candidates:
            isotope_options = []
            for order, spacing in ((1, formula_id.DM1), (2, formula_id.DM2)):
                for parent_index, parent_mass in enumerate(masses):
                    parent_candidates = [
                        candidate
                        for candidate in peaks[parent_index].get("candidates") or []
                        if candidate.get("assignment_eligible", True)
                    ]
                    if (
                        parent_index == index
                        or heights[parent_index] <= heights[index]
                        or not parent_candidates
                    ):
                        continue
                    residual = masses[index] - parent_mass - spacing
                    predicted = float(
                        (parent_candidates[0].get("iso_pred") or [0.0, 0.0])[order - 1]
                    )
                    ratio = heights[index] / heights[parent_index]
                    compatible = predicted > 0 and ratio <= max(
                        predicted * 3.0, predicted + 0.02
                    )
                    if abs(residual) <= _REAGENT_TOL_DA and compatible:
                        isotope_options.append(
                            (
                                abs(residual),
                                -heights[parent_index],
                                parent_mass,
                                parent_index,
                                order,
                                residual,
                                ratio,
                                predicted,
                            )
                        )
            if isotope_options:
                (
                    _,
                    _,
                    parent_mass,
                    _parent_index,
                    order,
                    residual,
                    ratio,
                    predicted,
                ) = min(isotope_options)
                interpretations.append(
                    {
                        "kind": "isotope",
                        "label": f"possible M+{order} isotope of m/z {parent_mass:.4f}",
                        "source": "measured spacing and predicted isotope envelope",
                        "related_mz": round(parent_mass, 4),
                        "isotope_order": order,
                        "evidence": [
                            f"spacing residual {residual * 1000:+.1f} mDa",
                            f"observed ratio {ratio:.3g}; predicted {predicted:.3g}",
                        ],
                        "exclude_from_analyte_assignment": True,
                    }
                )

        if not interpretations and not validated_candidates:
            authored_formula = peak.get("formula") or peak.get("suggested_formula")
            if authored_formula:
                interpretations.append(
                    {
                        "kind": "authored",
                        "label": str(peak.get("label") or authored_formula),
                        "formula": str(authored_formula),
                        "source": "saved review assignment",
                        "evidence": [
                            "preserved even though it is outside the current generated "
                            "candidate set"
                        ],
                        "exclude_from_analyte_assignment": False,
                    }
                )
            else:
                tolerance = peak.get("formula_tolerance") or {}
                if tolerance.get("candidate_generation_allowed") is False:
                    evidence = [
                        "formula generation withheld because run calibration does not "
                        "support a 10 ppm assignment radius",
                        str(tolerance.get("reason") or "calibration residuals degraded"),
                    ]
                elif peak.get("candidates"):
                    evidence = [
                        f"{len(peak['candidates'])} broad formula proposal(s) are "
                        f"available within {float(tolerance['proposal_ppm']):.0f} ppm, "
                        f"but none fit the run-validated {float(tolerance['ppm']):.1f} "
                        "ppm assignment radius"
                    ]
                elif tolerance.get("ppm") is not None:
                    evidence = [
                        "no plausible protonated-neutral formula proposal fits within "
                        f"{float(tolerance['proposal_ppm']):.0f} ppm "
                        f"({float(tolerance['proposal_mDa']):.2f} mDa here)"
                    ]
                else:
                    evidence = [
                        "no plausible protonated-neutral formula fits the current "
                        "exact-mass tolerance"
                    ]
                if peak.get("overlap"):
                    evidence.append("the peak overlaps a neighbouring channel")
                interpretations.append(
                    {
                        "kind": "unresolved",
                        "label": f"unresolved ion at m/z {masses[index]:.4f}",
                        "source": "strict local formula search",
                        "evidence": evidence,
                        "exclude_from_analyte_assignment": True,
                    }
                )
        if interpretations:
            peak["interpretation_candidates"] = interpretations
        else:
            peak.pop("interpretation_candidates", None)
    return peaks


def _assign_suggested_identities(peaks, *, assign_all_library=False):
    """Choose editable identity defaults without assigning one formula family twice.

    The matching favours each peak's ranked candidate order while maximising the number
    of distinct assignments. Compounds of interest have already received their score
    prior and provide the preferred isomer label. In maximum-coverage review mode, a
    candidate without a library name uses its formula as the best available label.
    Reagent ions remain outside the analyte matching.
    """
    options = {}
    for index, peak in enumerate(peaks):
        existing_label = peak.get("suggested_label")
        existing_formula = peak.pop("suggested_formula", None)
        existing_rank = peak.pop("suggested_candidate_rank", None)
        reagent = next(
            (
                flag.split(": ", 1)[1]
                for flag in peak.get("likely_artifact", [])
                if flag.startswith("reagent/cluster: ")
            ),
            None,
        )
        if reagent:
            peak["suggested_label"] = reagent
            continue
        interpretation = next(
            (
                item
                for item in peak.get("interpretation_candidates", [])
                if item.get("exclude_from_analyte_assignment")
            ),
            None,
        )
        if interpretation:
            peak["suggested_label"] = interpretation["label"]
            continue
        tolerance = peak.get("formula_tolerance") or {}
        if tolerance.get("automatic_assignment_allowed") is False:
            peak["suggested_label"] = existing_label or f"unknown m/z {peak['mz']:.3f}"
            if existing_formula:
                peak["suggested_formula"] = existing_formula
            if existing_rank:
                peak["suggested_candidate_rank"] = existing_rank
            continue
        candidates = peak.get("candidates") or []
        if not candidates and existing_label:
            peak["suggested_label"] = existing_label
            if existing_formula:
                peak["suggested_formula"] = existing_formula
            if existing_rank:
                peak["suggested_candidate_rank"] = existing_rank
            continue
        peak["suggested_label"] = f"unknown m/z {peak['mz']:.3f}"
        identity_options = []
        for rank, candidate in enumerate(candidates, start=1):
            if not candidate.get("assignment_eligible", True):
                continue
            formula = candidate.get("formula")
            compound_name = candidate.get("preferred_name") or candidate.get("name")
            label = compound_name or (formula if assign_all_library else None)
            eligible = assign_all_library or (
                rank == 1
                and compound_name
                and peak.get("id_confidence", 0) >= 0.6
            )
            if label and formula and eligible:
                identity_options.append(
                    {
                        "formula": formula,
                        "label": label,
                        "rank": rank,
                        "share": candidate.get("probability", 0.0),
                        "interest": bool(candidate.get("interest_matches")),
                    }
                )
        if identity_options:
            options[index] = identity_options

    assignments = optimal_identity_assignments(peaks, options)
    owners = {
        option["formula"].upper(): index for index, option in assignments.items()
    }

    for index, option in assignments.items():
        peaks[index]["suggested_label"] = option["label"]
        peaks[index]["suggested_formula"] = option["formula"]
        peaks[index]["suggested_candidate_rank"] = option["rank"]

    fallback_order = sorted(
        range(len(peaks)),
        key=lambda index: (
            -float(peaks[index].get("id_confidence", 0.0)),
            -float(peaks[index].get("height", 0.0)),
            float(peaks[index]["mz"]),
        ),
    )
    for index in fallback_order:
        peak = peaks[index]
        if index in assignments or not peak["suggested_label"].startswith("unknown m/z"):
            continue
        candidates = peak.get("candidates") or []
        top = candidates[0] if candidates else None
        if (
            top
            and top.get("assignment_eligible", True)
            and top.get("formula")
            and not top.get("name")
            and top["formula"].upper() not in owners
            and all(ch in "CHNO0123456789" for ch in top["formula"])
            and top["formula"].find("C") == 0
            and not peak.get("id_ambiguous")
            and len(candidates) >= 2
            and peak.get("id_confidence", 0) >= 0.9
        ):
            peak["suggested_label"] = top["formula"]
            peak["suggested_formula"] = top["formula"]
            peak["suggested_candidate_rank"] = 1
            owners[top["formula"].upper()] = index


def optimal_identity_assignments(peaks, options):
    """Maximise assignment count, then total candidate quality, deterministically."""
    if not options:
        return {}

    peak_indices = sorted(options)
    formulas = sorted(
        {
            option["formula"].upper()
            for peak_options in options.values()
            for option in peak_options
        }
    )
    n_peaks = len(peak_indices)
    n_formulas = len(formulas)
    strongest = sorted(
        peak_indices,
        key=lambda index: (
            -int(options[index][0]["interest"]),
            -float(options[index][0]["share"]),
            -float(
                peaks[index].get("height", peaks[index].get("abundance", 0.0))
            ),
            float(peaks[index]["mz"]),
        ),
    )
    peak_priority = {
        index: n_peaks - rank for rank, index in enumerate(strongest)
    }
    option_maps = {
        index: {option["formula"].upper(): option for option in options[index]}
        for index in peak_indices
    }
    max_edge_quality = 1_000_000_401 + 6 * n_peaks
    cardinality_bonus = (n_peaks + 1) * max_edge_quality
    invalid_cost = cardinality_bonus
    costs = []
    for index in peak_indices:
        row = []
        for formula in formulas:
            option = option_maps[index].get(formula)
            if option is None:
                row.append(invalid_cost)
                continue
            quality = round(float(option["share"]) * 1_000_000_000)
            rank_preference = 6 - min(5, int(option["rank"]))
            quality += (rank_preference - 1) * 100
            quality += int(option["interest"]) * (n_peaks + 1)
            quality += peak_priority[index] * rank_preference
            row.append(-cardinality_bonus - quality)
        row.extend([0] * n_peaks)
        costs.append(row)

    # One private dummy column per peak permits an unmatched result. The cardinality
    # bonus makes every valid edge preferable and dominates all quality differences,
    # so the rectangular Hungarian algorithm optimises quality only after coverage.
    n_columns = n_formulas + n_peaks
    row_potential = [0] * (n_peaks + 1)
    column_potential = [0] * (n_columns + 1)
    column_match = [0] * (n_columns + 1)
    predecessor = [0] * (n_columns + 1)
    infinity = 10**30
    for row_index in range(1, n_peaks + 1):
        column_match[0] = row_index
        column = 0
        min_cost = [infinity] * (n_columns + 1)
        used = [False] * (n_columns + 1)
        while True:
            used[column] = True
            active_row = column_match[column]
            delta = infinity
            next_column = 0
            for candidate_column in range(1, n_columns + 1):
                if used[candidate_column]:
                    continue
                reduced = (
                    costs[active_row - 1][candidate_column - 1]
                    - row_potential[active_row]
                    - column_potential[candidate_column]
                )
                if reduced < min_cost[candidate_column]:
                    min_cost[candidate_column] = reduced
                    predecessor[candidate_column] = column
                if min_cost[candidate_column] < delta:
                    delta = min_cost[candidate_column]
                    next_column = candidate_column
            for candidate_column in range(n_columns + 1):
                if used[candidate_column]:
                    row_potential[column_match[candidate_column]] += delta
                    column_potential[candidate_column] -= delta
                else:
                    min_cost[candidate_column] -= delta
            column = next_column
            if column_match[column] == 0:
                break
        while True:
            next_column = predecessor[column]
            column_match[column] = column_match[next_column]
            column = next_column
            if column == 0:
                break

    assignments = {}
    for column in range(1, n_formulas + 1):
        matched_row = column_match[column]
        if matched_row == 0:
            continue
        peak_index = peak_indices[matched_row - 1]
        option = option_maps[peak_index].get(formulas[column - 1])
        if option is not None:
            assignments[peak_index] = option
    return assignments


def _is_noise_artifact(flags):
    """True if a peak's likely_artifact flags mark it as instrument NOISE — a ringing
    comb, a low-prominence ripple/shoulder, a tail, or the reagent saturation-region
    skirt. These are never analytes and are dropped from the default `peaks` menu and
    the --auto-peaks panel. A plain reagent/cluster diagnostic ion is NOT noise: it is
    a real ion, kept but labelled, so it stays visible."""
    return any(
        ("ringing" in x)
        or ("low prominence" in x)
        or ("tail" in x)
        or ("saturation region" in x)
        for x in (flags or [])
    )


def _compact_peak(e):
    """Trim a fully-annotated peak to the fields needed for curation: the top
    candidate summary + flags, dropping the per-candidate isotope arrays and the
    long tail of low-probability formulas. `sniff peaks --full` keeps everything."""
    cands = e.get("candidates") or []
    top = cands[0] if cands else None
    out = {
        "mz": e["mz"],
        "height": e.get("height"),
        "rel_height": e.get("rel_height"),
        "prominence": e.get("prominence"),
        "neutral_mass": e.get("neutral_mass"),
        "formula_tolerance": e.get("formula_tolerance"),
        "suggested_label": e.get("suggested_label"),
    }
    if e.get("suggested_formula"):
        out["suggested_formula"] = e["suggested_formula"]
    if e.get("suggested_candidate_rank"):
        out["suggested_candidate_rank"] = e["suggested_candidate_rank"]
    if top:
        out["top_candidate"] = {
            "formula": top.get("formula"),
            "name": top.get("name"),
            "delta_mDa": top.get("delta_mDa"),
            "delta_ppm": top.get("delta_ppm"),
            "mass_match": top.get("mass_match"),
            "assignment_eligible": top.get("assignment_eligible"),
            "k": top.get("k"),
            "k_estimated": top.get("k_estimated"),
            "fragmentation_evidence": top.get("fragmentation_evidence"),
        }
    if "id_confidence" in e:
        out["id_confidence"] = e["id_confidence"]
    if "id_ambiguous" in e:
        out["id_ambiguous"] = e["id_ambiguous"]
    if "overlap" in e:  # keep the facts, drop the prose
        o = e["overlap"]  # (explained once in the header note)
        out["overlap"] = {
            "neighbor": o.get("neighbor"),
            "sep_mDa": o.get("sep_mDa"),
            "level": o.get("level"),
        }
    if "likely_artifact" in e:
        out["likely_artifact"] = e["likely_artifact"]
    if "interpretation_candidates" in e:
        out["interpretation_candidates"] = e["interpretation_candidates"]
    if "fragmentation_links" in e:
        out["fragmentation_links"] = e["fragmentation_links"]
    return out


def cmd_peaks(args):
    with h5py.File(args.h5, "r") as f:
        mass_axis = ptrms.load_mass_axis(f)
        a, b = mass_axis.a, mass_axis.b
        avg = np.where(
            np.isfinite(f["SPECdata/AverageSpec"][:]), f["SPECdata/AverageSpec"][:], 0.0
        )
        sig = assess_signal(f, avg=avg, a=a, b=b, mass_axis=mass_axis)
        if not sig["signal_present"]:
            _emit(
                {
                    "n_peaks": 0,
                    "signal_present": False,
                    "primary_ion_snr": sig["primary_snr"],
                    "mass_axis_calibration": mass_axis.to_dict(),
                    "note": "No significant signal. "
                    + sig["reason"]
                    + " Report this file as a blank/no-beam capture — do not "
                    "fabricate an analyte list from the noise.",
                    "peaks": [],
                },
                args.raw,
            )
            return
        R_phys = getattr(args, "R_phys", None) or 2400.0
        peaks = detect_peaks(
            f,
            args.min_height,
            args.max_peaks,
            args.mz_min,
            args.mz_max,
            R_phys=R_phys,
            mass_axis=mass_axis,
        )
        drift, peaks = annotate_peaks(
            peaks,
            avgspec=avg,
            a=a,
            b=b,
            R_phys=R_phys,
            mass_axis=mass_axis,
            candidate_pool_size=20,
        )
        fragmentation_context = apply_run_fragmentation_evidence(
            f,
            peaks,
            mass_axis=mass_axis,
            R_phys=R_phys,
            drift=drift,
        )
    # By default the menu excludes instrument-noise artifacts (ringing combs,
    # low-prominence ripples, reagent saturation-region skirt) so that copying the
    # list straight into a config can't ship a noise comb; reagent/cluster diagnostic
    # ions stay (labelled). --include-artifacts shows the raw list with every flag.
    include_art = getattr(args, "include_artifacts", False)
    n_noise = sum(1 for p in peaks if _is_noise_artifact(p.get("likely_artifact")))
    if not include_art:
        peaks = [p for p in peaks if not _is_noise_artifact(p.get("likely_artifact"))]
    _assign_suggested_identities(peaks)
    n_amb = sum(1 for p in peaks if p.get("id_ambiguous"))
    n_ovl = sum(1 for p in peaks if p.get("overlap"))
    n_with_formula_proposals = sum(bool(p.get("candidates")) for p in peaks)
    n_mass_validated = sum(
        any(candidate.get("assignment_eligible", True) for candidate in p.get("candidates", []))
        for p in peaks
    )
    # near-duplicate 'peak intervals': windows almost on top of each other
    dup_pairs = _window_overlap_pairs(peaks)
    full = getattr(args, "full", False)
    if full:
        note = (
            "Each peak lists broad candidate FORMULA proposals (up to 200 ppm) "
            "ranked by `probability` (combining exact-mass error, the measured vs "
            "predicted 13C(M+1)/heteroatom(M+2) isotope ratios, plausibility, and "
            "secondary known-pathway fragment co-variation where available). "
            "`mass_match` and `assignment_eligible` distinguish proposals inside the "
            "run-validated 5–10 ppm radius; only those may be assigned automatically. "
            "Use the complete evidence, not nearest-mass, to resolve isobars. "
            "`id_confidence` "
            "is the top candidate's normalized score/share; a sole candidate is "
            "not a 100% confidence estimate. `id_ambiguous` lists the "
            "close rivals when the call is not clear-cut; `overlap` flags a "
            "neighbouring peak whose spectral overlap adds quantification "
            "uncertainty (unresolved = worse than deconvolved). `name`/`k` are "
            "filled when the formula is in the rate table (else k_estimated). "
            "`iso_pred` vs `iso_obs` = predicted vs observed (M+1,M+2)/M. "
            "`fragmentation_evidence` reports condition-matched PTR Library pathways "
            "and temporal co-variation; it can reorder existing candidates but never "
            "create one or override `assignment_eligible`. `suggested_label` is an "
            "editable, globally unique high-ranked default. "
            "`interpretation_candidates` explains reagent, isotope, artefact, authored "
            "or unresolved channels when a neutral compound assignment would mislead. "
            "`prominence` is the apex's rise above local baseline in cps (real "
            "peak ≈ height; noise "
            "ripple ≈ 0). neutral_mass = mz − proton."
        )
        out_peaks = peaks
    else:
        note = (
            "Compact view (default). Each peak: `suggested_label` (a ready-to-use "
            "editable label; automatic library compounds are globally unique), "
            "`top_candidate` (best broad formula/name/mass-error proposal, chosen by "
            "isotope pattern, plausibility, and applicable fragmentation evidence, NOT "
            "nearest-mass; its `mass_match` shows "
            "whether it is inside the run-validated radius), `id_confidence` (a normalized "
            "top-candidate score used by conservative gates, not a calibrated "
            "probability), and, when "
            "relevant, `id_ambiguous` (close rivals), `overlap` (quantification "
            "uncertainty), and `interpretation_candidates` (evidence-backed reagent, "
            "isotope, artefact or unresolved-ion explanations). "
            "`prominence` is the apex's rise above its local baseline in cps: a "
            "real peak's ≈ its height, a noise ripple/shoulder's is near 0. "
            "neutral_mass = mz − proton. Pass `--full` for every candidate + the "
            "isotope arrays."
        )
        out_peaks = [_compact_peak(p) for p in peaks]
    note += (
        (
            f" This list is ALREADY cleaned: {n_noise} instrument-noise peaks (ringing "
            "combs, low-prominence ripples, reagent saturation-region skirt) were "
            "dropped — pass --include-artifacts to see them. Any peak here is safe "
            "to quantify."
        )
        if (n_noise and not include_art)
        else ""
    )
    if dup_pairs:
        examples = ", ".join(
            f"m/z {x['mz']:.4f}≈{y['mz']:.4f}" for x, y, _ in dup_pairs[:4]
        )
        note += (
            f" WARNING: {len(dup_pairs)} pair(s) of peaks have integration windows "
            "that almost coincide (>60% overlap) — e.g. "
            f"{examples}. These double-count the same signal; keep only one m/z "
            "from each pair in your config."
        )
    _emit(
        {
            "n_peaks": len(peaks),
            "n_noise_dropped": (0 if include_art else n_noise),
            "mass_drift": round(drift, 6),
            "mass_axis_calibration": mass_axis.to_dict(),
            "fragmentation_context": fragmentation_context,
            "n_with_formula_proposals": n_with_formula_proposals,
            "n_mass_validated": n_mass_validated,
            "n_provisional_only": n_with_formula_proposals - n_mass_validated,
            "n_without_formula_proposals": len(peaks) - n_with_formula_proposals,
            "n_ambiguous": n_amb,
            "n_overlapping": n_ovl,
            "n_window_overlap_pairs": len(dup_pairs),
            "note": note,
            "peaks": out_peaks,
        },
        args.raw,
    )


def cmd_segments(args):
    with h5py.File(args.h5, "r") as f:
        ncyc = int(f["SPECdata/Intensities"].shape[0])
        mass_axis = ptrms.load_mass_axis(f)
        sig = assess_signal(f, mass_axis=mass_axis)
        D = ptrms.build_discriminator(f, mass_axis=mass_axis)
        segs = ptrms.detect_segments(
            f,
            discriminator=D,
            min_duration=args.min_duration,
            grad_thr=args.grad_thr,
            high_ratio=args.high_ratio,
        )
        # always consolidate fragmented backgrounds; merge samples only on evidence
        segs = ptrms.merge_adjacent_segments(
            segs,
            high_gap=args.merge_high_gap,
            low_gap=200,
            discriminator=D,
            baseline=ptrms.discriminator_baseline(D),
            cap=ptrms.merge_gap_cap(f),
        )
    note = (
        "class 'high' = elevated signal (likely a sample); 'low' = "
        "background or pre-run setup. Final outputs use chronological "
        "sample_01/background_01 labels; do not ask for sample names. "
        "merged_segments > 1 marks plateaus joined across a short unclassified "
        "gap that never left its neighbours' level; merged_gaps gives each "
        "gap's length, its level range and why it merged."
    )
    gap_note = ptrms.merge_gaps_note(
        gap for s in segs for gap in s.get("merged_gaps", [])
    )
    if gap_note:
        note += " " + gap_note[0].upper() + gap_note[1:] + "."
    out = {"n_segments": len(segs), "note": note, "segments": segs}
    if not sig["signal_present"]:
        out["signal_present"] = False
        out["warning"] = (
            "No significant signal — "
            + sig["reason"]
            + " Segmentation is not meaningful for a blank file."
        )
    elif not segs:
        out["warning"] = (
            f"No stable plateaus found (file has {ncyc} cycles). If the run is very "
            "short, analyse the whole file as one interval (omit ranges); otherwise "
            "loosen --min-duration / --grad-thr."
        )
    _emit(out, args.raw)


def _window_overlap_pairs(peaks, R=1200.0, thresh=0.6):
    """Pairs of peaks whose default integration windows (±mz/2R) overlap by more
    than `thresh` of the narrower window — i.e. two 'peak intervals' that sit almost
    on top of each other. Such near-duplicates double-count the same signal, so the
    caller either merges them (auto path) or warns (curation)."""
    ps = sorted(peaks, key=lambda p: p["mz"])
    pairs = []
    for i in range(len(ps) - 1):
        m1, m2 = ps[i]["mz"], ps[i + 1]["mz"]
        hw1, hw2 = m1 / (2 * R), m2 / (2 * R)
        ov = (m1 + hw1) - (m2 - hw2)
        if ov > 0 and ov / min(2 * hw1, 2 * hw2) > thresh:
            pairs.append((ps[i], ps[i + 1], ov / min(2 * hw1, 2 * hw2)))
    return pairs


def _merge_overlapping_windows(peaks, R=1200.0, thresh=0.6):
    """Collapse peaks whose integration windows overlap by > `thresh`, keeping the
    taller of each pair — so an untargeted export never ships two nearly-coincident
    'peak intervals' for what the instrument cannot resolve as two ions."""
    ps = sorted(peaks, key=lambda p: p["mz"])
    out = []
    for p in ps:
        if out:
            q = out[-1]
            hwq, hwp = q["mz"] / (2 * R), p["mz"] / (2 * R)
            ov = (q["mz"] + hwq) - (p["mz"] - hwp)
            if ov > 0 and ov / min(2 * hwq, 2 * hwp) > thresh:
                if p.get("height", 0) > q.get("height", 0):
                    out[-1] = p  # keep the taller apex
                continue
        out.append(p)
    return out


def auto_peaks(
    f,
    *,
    min_height=1e-3,
    max_peaks=300,
    mz_min=15.0,
    mz_max=None,
    R=None,
    R_phys=None,
    mass_axis=None,
    compounds_of_interest=None,
    assign_all_library=False,
):
    """Deterministic peak panel for a file: detect, annotate, drop instrument-noise
    artifacts, collapse windows that coincide, and carry the suggested name.

    Reagent/cluster diagnostic ions are kept (labelled) — they are real ions and an
    untargeted panel usually wants them. Returns [] for a blank/no-beam file, which
    must never be turned into a fabricated analyte list."""
    R = 1200.0 if R is None else R
    R_phys = 2400.0 if R_phys is None else R_phys
    if mass_axis is None:
        mass_axis = ptrms.load_mass_axis(f)
    a, b = mass_axis.a, mass_axis.b
    avg = np.where(
        np.isfinite(f["SPECdata/AverageSpec"][:]), f["SPECdata/AverageSpec"][:], 0.0
    )
    if not assess_signal(
        f, avg=avg, a=a, b=b, mass_axis=mass_axis
    )["signal_present"]:
        return []
    peaks = detect_peaks(
        f,
        min_height,
        max_peaks,
        mz_min,
        mz_max,
        R_phys=R_phys,
        mass_axis=mass_axis,
    )
    drift, peaks = annotate_peaks(
        peaks,
        avgspec=avg,
        a=a,
        b=b,
        R=R,
        R_phys=R_phys,
        mass_axis=mass_axis,
        compounds_of_interest=compounds_of_interest,
        assign_all_library=assign_all_library,
    )
    peaks = [p for p in peaks if not _is_noise_artifact(p.get("likely_artifact"))]
    peaks = _merge_overlapping_windows(peaks, R=R)
    _assign_suggested_identities(
        peaks, assign_all_library=assign_all_library
    )
    out = []
    for p in peaks:
        sl = p.get("suggested_label", "") or ""
        # a real name -> label the channel; an "unknown m/z X" -> leave blank so the
        # CSV shows a clean `m<mz>` (the mass is already the variable name)
        o = {"mz": p["mz"], "label": "" if sl.startswith("unknown m/z") else sl}
        if p.get("suggested_formula"):
            o["formula"] = p["suggested_formula"]
        out.append(o)
    return out


def auto_ranges(
    f,
    *,
    min_duration=30,
    grad_thr=0.02,
    high_gap=None,
    low_gap=200,
    mass_axis=None,
):
    """Deterministic interval list: stable plateaus, consolidated, and named
    `sample_NN` / `background_NN` in chronological order. The name carries the class,
    which is what the analysis blanks against, so these labels are load-bearing.

    `high_gap=None` (the default) lets the merge judge each gap from the signal
    itself across ~60 s of acquisition; a number forces that cycle cap instead and 0
    never merges high plateaus. Merged gaps carry their provenance on the range as
    `merged_gaps` for JSON and config consumers; the review UI does not show a separate
    merge summary."""
    D = ptrms.build_discriminator(f, mass_axis=mass_axis)
    segs = ptrms.detect_segments(
        f, discriminator=D, min_duration=min_duration, grad_thr=grad_thr
    )
    segs = ptrms.merge_adjacent_segments(
        segs,
        high_gap=high_gap,
        low_gap=low_gap,
        discriminator=D,
        baseline=ptrms.discriminator_baseline(D),
        cap=ptrms.merge_gap_cap(f),
    )
    out = []
    counts = {"high": 0, "low": 0}
    for s in segs:
        kind = s["class"]
        counts[kind] += 1
        prefix = "sample" if kind == "high" else "background"
        entry = {
            "label": f"{prefix}_{counts[kind]:02d}",
            "start": s["start_cycle"],
            "end": s["end_cycle"],
            "unit": "cycle",
        }
        if s.get("merged_gaps"):
            entry["merged_gaps"] = s["merged_gaps"]
        out.append(entry)
    return out


def auto_ranges_note(ranges):
    """One line on what `auto_ranges` joined up, or "" when it joined nothing.

    The provenance rides on the ranges, so this is how the app mode — where nobody
    ever sees a command line — explains a merge instead of silently making one."""
    return ptrms.merge_gaps_note(
        gap for r in ranges or [] for gap in (r.get("merged_gaps") or [])
    )


def _auto_peaks(
    f,
    args,
    R=None,
    R_phys=None,
    mass_axis=None,
    compounds_of_interest=None,
):
    """Peaks for the --auto-peaks fallback: detect, annotate, DROP instrument-noise
    artifacts (ringing combs / low-prominence ripples), and carry each peak's
    `suggested_label` so a headless one-shot `analyze` yields a clean, labelled panel
    instead of the raw local-maxima list. Reagent/cluster ions are kept (labelled) —
    they are real ions, and an untargeted export usually wants them. Hand-curating a
    config still gives finer chemistry and segment judgment; this is a safe default,
    not a substitute for it."""
    return auto_peaks(
        f,
        min_height=args.min_height,
        max_peaks=args.max_peaks,
        mz_min=args.mz_min,
        mz_max=args.mz_max,
        R=R,
        R_phys=R_phys,
        mass_axis=mass_axis,
        compounds_of_interest=compounds_of_interest,
    )


def _load_peaks(args, f, settings=None, config=None, mass_axis=None):
    if args.peaks_json:
        return json.loads(args.peaks_json)
    cfg = config
    if args.config:
        if cfg is None:
            with open(args.config, encoding="utf-8") as fh:
                cfg = json.load(fh)
        if cfg.get("peaks"):
            return cfg["peaks"]
    if getattr(args, "auto_peaks", False):
        if settings is None:
            settings = resolve_analysis_settings(_load_config(args), args)
        return _auto_peaks(
            f,
            args,
            R=settings["R"],
            R_phys=settings["R_phys"],
            mass_axis=mass_axis,
            compounds_of_interest=(cfg or {}).get("compounds_of_interest"),
        )
    return None


def _load_ranges(args, f, mass_axis=None):
    if args.ranges_json:
        return json.loads(args.ranges_json)
    if args.config:
        with open(args.config, encoding="utf-8") as fh:
            cfg = json.load(fh)
        if cfg.get("ranges"):
            return cfg["ranges"]
    if getattr(args, "auto_segments", False):
        return auto_ranges(
            f,
            high_gap=getattr(args, "merge_high_gap", None),
            mass_axis=mass_axis,
        )
    return None


def _resolve_ranges(f, ranges_cfg):
    ncyc = int(f["SPECdata/Intensities"].shape[0])
    if not ranges_cfg:
        return {"All": (1, ncyc)}
    dur = ptrms.spec_duration_s(f)
    out = {}
    for r in ranges_cfg:
        if r.get("unit", "cycle") == "second":
            lo = max(1, round(r["start"] / dur) + 1)
            hi = min(ncyc, round(r["end"] / dur) + 1)
        else:
            lo, hi = max(1, int(r["start"])), min(ncyc, int(r["end"]))
        out[r["label"]] = (lo, hi)
    return out


def cmd_analyze(args):
    config = _load_config(args)
    with h5py.File(args.h5, "r") as f:
        mass_axis = ptrms.load_mass_axis(f)
        config = _migrate_loaded_config(config, args, mass_axis)
        settings = resolve_analysis_settings(config, args)
        peaks = _load_peaks(
            args, f, settings=settings, config=config, mass_axis=mass_axis
        )
        if not peaks:
            # distinguish a genuinely blank file from a missing peak list
            if getattr(args, "auto_peaks", False):
                sig = assess_signal(f, mass_axis=mass_axis)
                if not sig["signal_present"]:
                    _emit(
                        {
                            "out": None,
                            "n_rows": 0,
                            "n_peaks": 0,
                            "signal_present": False,
                            "primary_ion_snr": sig["primary_snr"],
                            "mass_axis_calibration": mass_axis.to_dict(),
                            "note": "No output written. "
                            + sig["reason"]
                            + " Report this file as a blank/no-beam capture.",
                        },
                        args.raw,
                    )
                    return
                sys.exit(
                    "No analyte peaks cleared the noise threshold (file has "
                    "signal but no resolvable peaks). Inspect with `sniff peaks`."
                )
            sys.exit(
                "No peaks. Pass --peaks-json '[{\"mz\":..}]', --config, or --auto-peaks."
            )
        masses = [float(p["mz"]) for p in peaks]
        labels = {
            float(p["mz"]): formula_id.identity_label(p.get("label"), p.get("formula"))
            for p in peaks
        }
        ranges = _resolve_ranges(f, _load_ranges(args, f, mass_axis=mass_axis))

        R = settings["R"]
        R_phys = settings["R_phys"]
        primary_mz = settings["primary_mz"]

        # resolve rate constants once (for kinetic correction and/or humid flags)
        resolved = ptrms.resolve_k(peaks, ptrms.load_rate_constants())
        k_map = resolved if settings["kinetic"] else None
        humid_masses = {
            m for m, info in resolved.items() if "humid" in info.get("flags", [])
        }

        # humidity proxy (per-cycle water-cluster ratio) — always computed if any
        # humid compound is present, so it can be reported as a diagnostic
        hum_ratio = (
            ptrms.water_cluster_ratio(
                f, primary_mz=primary_mz, R=R, mass_axis=mass_axis
            )
            if humid_masses
            else None
        )

        # per-interval windows (default on): re-centre each isolated peak's window
        # on every interval's own spectrum. Disable with --no-per-interval to get
        # one whole-run window per compound (the pre-2026-08 behaviour).
        real_ranges = not (len(ranges) == 1 and "All" in ranges)
        per_range = (
            ranges if (real_ranges and settings["per_interval_windows"]) else None
        )
        isotope_plan = (
            ptrms.isotopes.build_isotope_plan(peaks, R_phys=R_phys)
            if settings["isotope_mode"] == "formula-v1"
            else None
        )
        extraction_masses = (
            isotope_plan["extraction_masses"] if isotope_plan is not None else masses
        )
        fit_diagnostics = {}
        traces, _ = ptrms.extract_traces(
            f,
            extraction_masses,
            R=R,
            R_phys=R_phys,
            windows=_peak_windows(peaks) or None,
            per_range=per_range,
            mass_axis=mass_axis,
            peak_fit_model=settings["peak_fit"],
            fit_diagnostics=fit_diagnostics,
        )
        rows, params = ptrms.quantify(
            traces,
            f,
            ranges,
            K=settings["K"],
            primary_mz=primary_mz,
            molar_volume=settings["molar_volume"],
            R_used=R,
            k_map=k_map,
            k_anchor=settings["k_anchor"],
            humid_masses=(humid_masses if settings["humidity_correct"] else None),
            humidity_ratio=hum_ratio,
            humidity_ref=settings["humidity_ref"],
            humidity_p=settings["humidity_p"],
            mass_axis=mass_axis,
            isotope_plan=isotope_plan,
            isotope_abundance_basis=settings["isotope_abundance_basis"],
        )
        params["peak_fit"] = fit_diagnostics
        apexes = {m: traces[m][1] for m in masses}
        humidity_ref = settings["humidity_ref"]
        if humidity_ref is None and hum_ratio is not None:
            good = np.isfinite(hum_ratio) & (hum_ratio > 0)
            humidity_ref = float(np.median(hum_ratio[good])) if good.any() else None
        sources = _effective_sources(
            settings, humidity_ref, params.get("molar_volume_source")
        )
        params.update(
            {
                "R_phys": R_phys,
                "whole_run_windows": settings["whole_run_windows"],
                "per_interval_windows": settings["per_interval_windows"],
                "humidity_correct": settings["humidity_correct"],
                "humidity_p": settings["humidity_p"],
                "humidity_ref": humidity_ref,
                "humidity_ref_source": sources["humidity_ref"],
                "molar_volume_source": sources["molar_volume"],
                "sources": sources,
                "mass_axis_calibration": mass_axis.to_dict(),
            }
        )

        # per-range humidity proxy + cross-range spread diagnostic
        humidity_report = None
        if humid_masses and hum_ratio is not None:
            per_range = {}
            for label, (lo, hi) in ranges.items():
                seg = hum_ratio[lo - 1 : hi]
                seg = seg[np.isfinite(seg)]
                per_range[label] = round(float(seg.mean()), 5) if seg.size else None
            vals = [v for v in per_range.values() if v]
            spread = (max(vals) - min(vals)) / (sum(vals) / len(vals)) if vals else 0.0
            humidity_report = {
                "humid_compounds": sorted(f"{m:.3f}" for m in humid_masses),
                "proxy": _humidity_proxy_label(primary_mz),
                "per_range": per_range,
                "cross_range_spread_pct": round(100 * spread, 1),
                "corrected": params.get("humidity_corrected", False),
                "p": params.get("humidity_p"),
                "reference_ratio": params.get("humidity_ref"),
            }
            if spread > 0.1 and not params.get("humidity_corrected"):
                humidity_report["warning"] = (
                    f"Humidity varies {100 * spread:.0f}% across ranges — relative "
                    "concentrations of the humid compounds are confounded. Add "
                    "--humidity-correct (needs a calibrated --humidity-p for accuracy)."
                )

    _write_csv(
        args.out,
        args.h5,
        rows,
        labels,
        args.sep,
        ranges=ranges,
        include_cycle_rows=args.include_cycle_rows,
    )

    # Peaks beyond the shared correction may have snapped to a neighbour, be absent,
    # or be mis-assigned. On fallback, retain the legacy median-relative diagnostic.
    warn = []
    rel = np.array([apexes[m] / m for m in masses])
    drift = 1.0 if mass_axis.applied else float(np.median(rel))
    for m in masses:
        resid_da = (
            apexes[m] - m
            if mass_axis.applied
            else (apexes[m] / m - drift) * m
        )
        if abs(resid_da) > 0.03:
            warn.append(
                f"m{m:.3f}: apex {apexes[m]:.4f} deviates "
                f"{resid_da:+.3f} Da beyond the run's mass-axis correction "
                f"(check assignment / possible peak overlap)"
            )
    note = None
    if settings["K"] is None and params.get("concentration_available"):
        note = (
            "Concentration uses K derived from the file's own calibration; "
            "absolute scale may differ from a specific PTR-MS Viewer project. "
            "Run `calibrate` against a reference CSV, or pass --K, to match exactly."
        )
    if not params.get("concentration_available"):
        note = "No primary-ion/pre-computed data: Conc columns are NaN. Pass --K and ensure a primary-ion peak exists."
    trans_note = None
    if not params.get("transmission_available", True):
        trans_note = (
            "This file carries no transmission curve, so unit transmission "
            "was assumed: Corrected == Raw. Absolute Corrected/Conc values "
            "are uncalibrated for mass-dependent transmission."
        )
        note = (note + " " + trans_note) if note else trans_note

    # per-compound kinetic reporting + humidity flags
    kinetic_info = None
    if k_map is not None:
        used, estimated, missing, humid = {}, [], [], []
        for m in masses:
            info = k_map.get(m, {})
            if info.get("k") and not info.get("k_estimated"):
                used[f"{m:.3f}"] = {"k": info["k"], "source": info["source"]}
                if "humid" in info.get("flags", []):
                    humid.append(f"{m:.3f}")
            elif info.get("k"):  # k exists but is estimated -> kept on shared K
                estimated.append(
                    f"{m:.3f}" + (f" ({info['source']})" if info.get("source") else "")
                )
            else:
                missing.append(
                    f"{m:.3f}" + (f" ({info['source']})" if info.get("source") else "")
                )
        kinetic_info = {
            "k_anchor": settings["k_anchor"],
            "resolved": used,
            "estimated_shared_K": estimated,
            "no_k": missing,
        }
        if humid:
            kinetic_info["humidity_warning"] = (
                "These masses have proton affinity near water (HCN/formaldehyde/"
                "H2S/acids/ammonia): a fixed k is unreliable — sensitivity is "
                "humidity/temperature dependent. Use a dedicated standard/humidity "
                f"model for: {humid}"
            )

    # sample-vs-background diagnostic: flag channels that behave like instrument
    # background / contamination rather than analytes — higher in backgrounds than
    # samples (S/B < 1) and/or drifting monotonically across the run. Lets the
    # agent relabel/drop them (e.g. an `unknown m/z 331` that is really background)
    # instead of shipping a bare "unknown" channel. Needs sample_/background_ labels.
    background_report = None
    samp_labels = [l for l in ranges if l.startswith("sample")]
    bg_labels = [l for l in ranges if l.startswith("background")]
    if samp_labels and bg_labels:
        by_mass = {}
        for r in rows:
            by_mass.setdefault(r["mass"], {})[r["range"]] = r["raw"]["Average"]
        flagged = {}
        for m in masses:
            per = by_mass.get(m, {})
            s = [per[l] for l in samp_labels if l in per]
            bg = [per[l] for l in bg_labels if l in per]
            if not s or not bg:
                continue
            smean, bmean = sum(s) / len(s), sum(bg) / len(bg)
            sb = (smean / bmean) if bmean else float("inf")
            bg_series = [per[l] for l in sorted(bg_labels) if l in per]
            trend = (
                (bg_series[-1] / bg_series[0])
                if len(bg_series) >= 2 and bg_series[0]
                else None
            )
            if sb < 0.9:  # not elevated in samples -> background-like
                flagged[f"{m:.3f}"] = {
                    "label": labels.get(m, ""),
                    "S_over_B": round(sb, 2),
                    "bg_trend_last_over_first": round(trend, 2) if trend else None,
                }
        background_report = {
            "metric": "mean Raw over sample_* vs background_* ranges (S/B); "
            "channels with S/B < 0.9 are flagged below; "
            "bg_trend = last/first background range (>1 = rising across run)",
            "n_samples": len(samp_labels),
            "n_backgrounds": len(bg_labels),
            "background_like": flagged,
        }
        if flagged:
            background_report["warning"] = (
                f"{len(flagged)} channel(s) are higher in backgrounds than samples "
                "(S/B < 0.9) — likely instrument background/contamination, not breath "
                "analytes (real VOCs have S/B >> 1). Relabel these as 'background m/z ...' "
                "or drop them from an analyte panel. Scrutinise unidentified/high-m/z "
                "peaks first; reagent/cluster diagnostic ions flagging here is expected."
            )

    n_cycle_rows = len(ranges) if args.include_cycle_rows else 0
    _emit(
        {
            "out": args.out,
            "n_rows": len(rows) + n_cycle_rows,
            "n_quant_rows": len(rows),
            "n_cycle_rows": n_cycle_rows,
            "n_peaks": len(masses),
            "n_ranges": len(ranges),
            "measured_apexes": {f"{m:.3f}": round(apexes[m], 4) for m in masses},
            "params": params,
            "kinetic": kinetic_info,
            "humidity": humidity_report,
            "background": background_report,
            "apex_warnings": warn,
            "note": note,
        },
        args.raw,
    )


def cmd_viz(args):
    """Review app for an EXISTING peak list + time ranges. `viz` does NOT detect
    peaks or segments — build those with `peaks`/`segments`, curate them into a
    config, and pass it via --config (or --peaks-json/--ranges-json).

    Two modes:
      * serve (default): run a localhost server, open the browser, and LIVE-SAVE
        every edit to the --config file. When the expert clicks 'Done' it runs the
        full-precision analysis and writes the results CSV (--out). Blocks until
        Done or --timeout. Because it blocks on the browser, run it backgrounded.
      * --html review.html: write a standalone, portable HTML file instead (no
        server, no CSV; edits exported via the page's Download button)."""
    from . import viz

    config = _load_config(args)
    x_axis_unit = resolve_x_axis_unit(config, args)
    with h5py.File(args.h5, "r") as f:
        mass_axis = ptrms.load_mass_axis(f)
        config = _migrate_loaded_config(config, args, mass_axis)
        settings = resolve_analysis_settings(config, args)
        peaks = _load_peaks(
            args, f, settings=settings, config=config, mass_axis=mass_axis
        )
        ranges_cfg = _load_ranges(args, f, mass_axis=mass_axis)
        if not peaks or not ranges_cfg:
            sys.exit(
                "viz needs an explicit peak list AND time ranges — it does not "
                "detect them. Pass --config with 'peaks' and 'ranges' (or "
                "--peaks-json/--ranges-json). Build them with `sniff peaks` and "
                "`sniff segments`, then curate into the config."
            )
        R = settings["R"]
        R_phys = settings["R_phys"]
        # Large files take ~30-90 s to load and pre-compute traces BEFORE the server
        # starts. Announce it so a watching agent waits for "review app running at …"
        # (below) rather than polling the port — which refuses until this finishes.
        print(
            "sniff: preparing the review (loading the file + computing traces; large "
            "files take ~30-90 s) — the URL is printed when it's ready…",
            file=sys.stderr,
            flush=True,
        )
        data = viz.build_viz_data(
            f,
            peaks,
            ranges_cfg,
            R=R,
            R_phys=R_phys,
            primary_mz=settings["primary_mz"],
            K=settings["K"],
            molar_volume=settings["molar_volume"],
            analysis_settings=settings,
            config_base=config,
            x_axis_unit=x_axis_unit,
            merge_note=config.get("merge_note") or "",
            mass_axis=mass_axis,
        )

    serve_mode = args.serve if args.serve is not None else (not args.html)
    if serve_mode:
        cfg_path = args.config or args.save_config
        if not cfg_path:
            sys.exit(
                "serve mode needs a config path to save to: pass --config PATH "
                "(the source) or --save-config PATH."
            )
        if not os.path.exists(cfg_path):
            initial = dict(config)
            initial.update({"peaks": peaks, "ranges": ranges_cfg})
            initial["mass_axis_domain"] = ptrms.MASS_AXIS_CONFIG_DOMAIN
            initial["mass_axis_version"] = ptrms.MASS_AXIS_CONFIG_VERSION
            initial["analyze"] = {
                **(config.get("analyze") or {}),
                **{key: settings[key] for key in _ANALYSIS_DEFAULTS},
            }
            initial["viz"] = {**(config.get("viz") or {}), "x_axis_unit": x_axis_unit}
            with open(cfg_path, "w", encoding="utf-8") as fh:
                json.dump(initial, fh, indent=2)
        review_token = secrets.token_urlsafe(24)
        html = viz.render_html(
            data,
            config_path=cfg_path,
            page_token=review_token,
        )
        run = lambda cfg: analyze_config_to_csv(
            args.h5, cfg, args.out, args.sep, args.include_cycle_rows
        )
        spec_fn = lambda lo, hi: interval_spectrum(args.h5, lo, hi)
        def peak_preview_fn(lo, hi):
            with h5py.File(args.h5, "r") as source:
                axis = ptrms.load_mass_axis(source)
                return viz.preview_peak(
                    source,
                    lo,
                    hi,
                    R=settings["R"],
                    mass_axis=axis,
                    compounds_of_interest=config.get("compounds_of_interest"),
                )

        final, finished, summary = viz.serve(
            html,
            cfg_path,
            port=args.port,
            timeout=args.timeout,
            open_browser=not args.no_open,
            run_analysis=run,
            spectrum_fn=spec_fn,
            peak_preview_fn=peak_preview_fn,
        )
        if final is not None:
            cfg = final
        else:
            with open(cfg_path, encoding="utf-8") as fh:
                cfg = json.load(fh)
        if summary is None:
            summary = run(cfg)
        summary.update(
            {
                "mode": "served",
                "config": cfg_path,
                "review_finished": finished,
                "note": (
                    "Done: the expert's review was saved to the config and the "
                    "full-precision analysis was written to the CSV."
                    if finished
                    else "Review timed out; the analysis ran on the last auto-saved "
                    "config. Re-open with the same command to continue editing."
                ),
            }
        )
        _emit(summary, args.raw)
    else:
        html = viz.render_html(data)
        with open(args.html, "w", encoding="utf-8") as fh:
            fh.write(html)
        _emit(
            {
                "out": args.html,
                "served": False,
                "n_peaks": len(data["peaks"]),
                "n_ranges": len(data["ranges"]),
                "n_cycles": data["meta"]["ncyc"],
                "concentration_available": data["meta"]["concentration_available"],
                "note": "Standalone portable review app written. Open in a browser to "
                "sanity-check/tweak and Download config.json to hand back for "
                "`sniff analyze`. For a live-saving session that also writes the "
                "CSV on Done, drop --html and pass --config cfg.json.",
            },
            args.raw,
        )


def cmd_app(args):
    """Run the persistent local review app, blocking until interrupted.

    This is the mode for a user at a desk rather than an agent: the server stays up
    between files, opening a file loads the config saved beside it or runs the
    deterministic pipeline to make one, and exporting writes the CSV without shutting
    anything down.

    The app normally opens in a browser tab. Frozen macOS and Windows bundles open in a
    desktop window instead, because a double-clicked app has no terminal to read a URL
    out of. Frozen Linux bundles deliberately retain the browser for portability;
    `--window` asks for a desktop window explicitly, and `--no-browser` declines both.
    """
    from . import app as app_mode

    if args.h5 and not os.path.isfile(args.h5):
        raise SystemExit(f"sniff: file not found: {args.h5}")
    agent = args.agent or os.environ.get("SNIFF_AGENT_URL")
    app_mode.serve_app(
        port=args.port,
        open_browser=not args.no_browser,
        agent_url=agent,
        agent_timeout=args.agent_timeout,
        initial=args.h5,
        window=_wants_window(args),
    )
    return 0


def _wants_window(args) -> bool:
    """Window, browser tab, or neither: the precedence, in one place.

    An explicit `--window` wins. Otherwise packaged macOS and Windows bundles default
    to the window unless `--no-browser` says nothing is to be opened. Linux deliberately
    keeps the browser-based surface so its package does not depend on a particular GUI
    toolkit or web renderer. Anywhere else the browser is left exactly as it was.
    """
    if getattr(args, "window", False):
        return True
    desktop_bundle = bool(getattr(sys, "frozen", False)) and sys.platform in (
        "darwin",
        "win32",
    )
    return desktop_bundle and not args.no_browser


def cmd_rates(args):
    """Look up / list proton-transfer rate constants."""
    tbl = ptrms.load_rate_constants()
    if not tbl:
        sys.exit("bundled rate_constants.json could not be loaded.")
    comps = tbl["compounds"]
    q = args.query
    if q:
        try:
            target = float(q)
            comps = [c for c in comps if abs(c["mz"] - target) < 0.3]
        except ValueError:
            ql = q.lower()
            comps = [
                c
                for c in comps
                if ql in c["name"].lower()
                or ql in c["formula"].lower()
                or any(ql in n.lower() for n in c.get("isomers", []))
            ]
    _emit(
        {
            "units": "1e-9 cm3/s",
            "n": len(comps),
            "source": tbl.get("_source", ""),
            "compounds": comps,
        },
        args.raw,
    )


def cmd_calibrate(args):
    """Fit the concentration constant K against a reference Viewer CSV."""
    config = _load_config(args)
    ref = _parse_viewer_csv(args.reference)
    ref_conc = {k: v["con"] for k, v in ref.items()}
    with h5py.File(args.h5, "r") as f:
        mass_axis = ptrms.load_mass_axis(f)
        config = _migrate_loaded_config(config, args, mass_axis)
        settings = resolve_analysis_settings(config, args)
        peaks = _load_peaks(
            args, f, settings=settings, config=config, mass_axis=mass_axis
        )
        # default: calibrate on whatever masses appear in the reference
        if not peaks:
            masses = sorted({mz for (mz, _) in ref_conc})
        else:
            masses = [float(p["mz"]) for p in peaks]
        ranges = _resolve_ranges(f, _load_ranges(args, f, mass_axis=mass_axis))
        if len(ranges) == 1 and "All" in ranges:
            # derive ranges from the reference's own labels via its Cycle rows
            ranges = _ranges_from_reference(args.reference)
        R = settings["R"]
        R_phys = settings["R_phys"]
        primary_mz = settings["primary_mz"]
        traces, _ = ptrms.extract_traces(
            f,
            masses,
            R=R,
            R_phys=R_phys,
            mass_axis=mass_axis,
            peak_fit_model=settings["peak_fit"],
        )
        K, resid, n = ptrms.calibrate_K(
            f,
            traces,
            ref_conc,
            ranges,
            primary_mz=primary_mz,
            R_used=R,
            mass_axis=mass_axis,
        )
        K_file = ptrms.derive_K(
            f,
            ptrms.extract_primary(f, primary_mz, R, mass_axis=mass_axis),
        )
    _emit(
        {
            "K_calibrated": K,
            "K_from_file": K_file,
            "calibration_points": n,
            "residual_median_pct": resid,
            "usage": f"pass --K {K} to analyze to match this reference",
        },
        args.raw,
    )


def _ranges_from_reference(path):
    """Recover {label:(lo,hi)} from a reference CSV's own Cycle variable rows."""
    out = {}
    with open(path, encoding="utf-8-sig") as fh:
        for r in csv.reader(fh, delimiter=";"):
            if len(r) < 6 or r[1].strip() != "Cycle":
                continue

            def n(x):
                return int(float(x.replace(",", ".")))

            out[r[2].strip()] = (n(r[4]), n(r[3]))  # (min, max)
    return out


def _write_csv(path, src, rows, labels, sep, ranges=None, include_cycle_rows=False):
    header = [
        "File",
        "Variable",
        "Range",
        "Max(Raw)",
        "Min(Raw)",
        "Average(Raw)",
        "Deviation(Raw)",
        "Max(Corrected)",
        "Min(Corrected)",
        "Average(Corrected)",
        "Deviation(Corrected)",
        "Max(Conc)",
        "Min(Conc)",
        "Average(Conc)",
        "Deviation(Conc)",
        "Max(Conc [ug])",
        "Min(Conc [ug])",
        "Average(Conc [ug])",
        "Deviation(Conc [ug])",
    ]

    def fmt(v):
        s = f"{v:.6f}"
        return s.replace(".", ",") if sep == ";" else s

    with ExitStack() as stack:
        fh = (
            sys.stdout
            if path == "-"
            else stack.enter_context(open(path, "w", newline="", encoding="utf-8-sig"))
        )
        w = csv.writer(fh, delimiter=sep)
        w.writerow(header)
        for r in rows:
            m = r["mass"]
            lbl = labels.get(m, "")
            var = f"m{m:.3f}".replace(".", ",") + (f" ({lbl})" if lbl else "")
            row = [src, var, r["range"]]
            for q in ("raw", "cor", "con", "ug"):
                s = r[q]
                row += [
                    fmt(s["Max"]),
                    fmt(s["Min"]),
                    fmt(s["Average"]),
                    fmt(s["Deviation"]),
                ]
            w.writerow(row)
        if include_cycle_rows:
            for label, (lo, hi) in (ranges or {}).items():
                cycles = np.arange(lo, hi + 1, dtype=float)
                deviation = float(np.std(cycles, ddof=1)) if cycles.size > 1 else 0.0
                w.writerow(
                    [
                        src,
                        "Cycle",
                        label,
                        str(hi),
                        str(lo),
                        fmt(float(cycles.mean())),
                        fmt(deviation),
                        *("" for _ in range(12)),
                    ]
                )


def interval_spectrum(h5_path, lo, hi, block=512):
    """Average mass spectrum over cycles [lo, hi] (1-based inclusive) as a list of
    ints (index = timebin), for the viz app's per-interval spectrum view. Streams
    in cycle blocks so memory stays bounded regardless of interval length."""
    with h5py.File(h5_path, "r") as f:
        inten = f["SPECdata/Intensities"]
        ncyc, nbin = inten.shape
        lo = max(1, int(lo))
        hi = min(int(ncyc), int(hi))
        if hi < lo:
            lo, hi = hi, lo
        acc = np.zeros(nbin, dtype=np.float64)
        n = 0
        for i in range(lo - 1, hi, block):  # 0-based half-open
            j = min(i + block, hi)
            acc += np.asarray(inten[i:j, :], dtype=np.float64).sum(axis=0)
            n += j - i
        avg = acc / max(1, n)
    return [round(x) for x in avg]


def analyze_config_to_csv(h5_path, config, out, sep=";", include_cycle_rows=True):
    """Run the full analyze pipeline from a viz/exported config dict -> results CSV.

    Honours `config['analyze']` settings (R, K, molar_volume, kinetic, k_anchor,
    humidity_*) when present. Used by `viz` after an interactive review and as the
    shared quantify path. Returns a small JSON summary."""
    with h5py.File(h5_path, "r") as f:
        mass_axis = ptrms.load_mass_axis(f)
        config, _ = ptrms.migrate_config_mass_axis(config, mass_axis)
        settings = resolve_analysis_settings(config)
        peaks = config["peaks"]
        ranges_cfg = config.get("ranges") or []
        masses = [float(p["mz"]) for p in peaks]
        labels = {
            float(p["mz"]): formula_id.identity_label(p.get("label"), p.get("formula"))
            for p in peaks
        }
        ranges = _resolve_ranges(f, ranges_cfg)
        R = settings["R"]
        R_phys = settings["R_phys"]
        primary_mz = settings["primary_mz"]
        resolved = ptrms.resolve_k(peaks, ptrms.load_rate_constants())
        k_map = resolved if settings["kinetic"] else None
        humid_masses = {
            m for m, info in resolved.items() if "humid" in info.get("flags", [])
        }
        hum_ratio = (
            ptrms.water_cluster_ratio(
                f, primary_mz=primary_mz, R=R, mass_axis=mass_axis
            )
            if humid_masses
            else None
        )

        # per-peak integration-window overrides: `window` is either a full-width
        # number (symmetric) or {"left":hwL,"right":hwR} half-widths (asymmetric)
        def _winlr(p):
            w = p["window"]
            if isinstance(w, dict):
                return (float(w["left"]), float(w["right"]))
            return (float(w) / 2.0, float(w) / 2.0)

        windows = {float(p["mz"]): _winlr(p) for p in peaks if p.get("window")}
        # per-interval windows: each interval integrates each isolated peak with an
        # apex/window re-centred on that interval's own spectrum (peaks drift). On
        # by default when real intervals exist; matches the viz per-interval review.
        per_range = (
            ranges if (ranges_cfg and settings["per_interval_windows"]) else None
        )
        isotope_plan = (
            ptrms.isotopes.build_isotope_plan(peaks, R_phys=R_phys)
            if settings["isotope_mode"] == "formula-v1"
            else None
        )
        extraction_masses = (
            isotope_plan["extraction_masses"] if isotope_plan is not None else masses
        )
        fit_diagnostics = {}
        traces, _ = ptrms.extract_traces(
            f,
            extraction_masses,
            R=R,
            R_phys=R_phys,
            windows=windows or None,
            per_range=per_range,
            mass_axis=mass_axis,
            peak_fit_model=settings["peak_fit"],
            fit_diagnostics=fit_diagnostics,
        )
        rows, params = ptrms.quantify(
            traces,
            f,
            ranges,
            K=settings["K"],
            primary_mz=primary_mz,
            molar_volume=settings["molar_volume"],
            R_used=R,
            k_map=k_map,
            k_anchor=settings["k_anchor"],
            humid_masses=(humid_masses if settings["humidity_correct"] else None),
            humidity_ratio=hum_ratio,
            humidity_ref=settings["humidity_ref"],
            humidity_p=settings["humidity_p"],
            mass_axis=mass_axis,
            isotope_plan=isotope_plan,
            isotope_abundance_basis=settings["isotope_abundance_basis"],
        )
        params["peak_fit"] = fit_diagnostics
        humidity_ref = settings["humidity_ref"]
        if humidity_ref is None and hum_ratio is not None:
            good = np.isfinite(hum_ratio) & (hum_ratio > 0)
            humidity_ref = float(np.median(hum_ratio[good])) if good.any() else None
        sources = _effective_sources(
            settings, humidity_ref, params.get("molar_volume_source")
        )
        params.update(
            {
                "R_phys": R_phys,
                "whole_run_windows": settings["whole_run_windows"],
                "per_interval_windows": settings["per_interval_windows"],
                "humidity_correct": settings["humidity_correct"],
                "humidity_p": settings["humidity_p"],
                "humidity_ref": humidity_ref,
                "humidity_ref_source": sources["humidity_ref"],
                "molar_volume_source": sources["molar_volume"],
                "sources": sources,
                "mass_axis_calibration": mass_axis.to_dict(),
            }
        )
    _write_csv(
        out,
        h5_path,
        rows,
        labels,
        sep,
        ranges=ranges,
        include_cycle_rows=include_cycle_rows,
    )
    n_cycle = len(ranges) if include_cycle_rows else 0
    return {
        "out": out,
        "n_rows": len(rows) + n_cycle,
        "n_peaks": len(masses),
        "n_ranges": len(ranges),
        "params": params,
        "K": params.get("K"),
        "molar_volume": params.get("molar_volume"),
        "primary_mz": params.get("primary_mz"),
        "kinetic": params.get("kinetic"),
        "concentration_available": params.get("concentration_available"),
    }


def _parse_viewer_csv(path):
    def num(s):
        return float(s.replace(",", ".")) if s.strip() else float("nan")

    out = {}
    with open(path, encoding="utf-8-sig") as fh:
        rd = csv.reader(fh, delimiter=";")
        next(rd, None)
        for row in rd:
            if len(row) < 18 or row[1].strip() == "Cycle":
                continue
            mstr = row[1].split()[0].lstrip("m").replace(",", ".")
            try:
                mz = round(float(mstr), 3)
            except ValueError:
                continue
            out[(mz, row[2].strip())] = {
                "raw": num(row[5]),
                "cor": num(row[9]),
                "con": num(row[13]),
                "ug": num(row[17]),
            }
    return out


def cmd_compare(args):
    mine, ref = _parse_viewer_csv(args.mine), _parse_viewer_csv(args.reference)
    errs = {"raw": [], "cor": [], "con": [], "ug": []}
    per_mass, n = {}, 0
    for key, rr in ref.items():
        if key not in mine:
            continue
        n += 1
        mm = mine[key]
        for q, values in errs.items():
            if rr[q] and np.isfinite(rr[q]) and np.isfinite(mm[q]):
                e = abs(100 * (mm[q] - rr[q]) / rr[q])
                values.append(e)
                if q == "raw":
                    per_mass.setdefault(key[0], []).append(e)
    if n == 0:
        sys.exit("No overlapping (mass, range) rows between the two files.")
    summary = {"matched_rows": n}
    for q in ("raw", "cor", "con", "ug"):
        a = np.array(errs[q]) if errs[q] else np.array([np.nan])
        summary[q] = {
            "median_pct": round(float(np.nanmedian(a)), 2),
            "mean_pct": round(float(np.nanmean(a)), 2),
            "p90_pct": round(float(np.nanpercentile(a, 90)), 2),
            "max_pct": round(float(np.nanmax(a)), 2),
        }
    if args.per_mass:
        summary["per_mass_raw_median_pct"] = {
            f"{m:.3f}": round(float(np.median(per_mass[m])), 2)
            for m in sorted(per_mass)
        }
    _emit(summary, args.raw)


def main():
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument(
        "--pretty",
        dest="raw",
        action="store_false",
        help="Pretty-print JSON (default compact)",
    )

    p = argparse.ArgumentParser(
        prog="sniff",
        description=__doc__,
        parents=[common],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    pi = sub.add_parser(
        "inspect", parents=[common], help="File metadata & calibration (JSON)"
    )
    pi.add_argument("h5")
    pi.set_defaults(func=cmd_inspect)

    pp = sub.add_parser(
        "peaks", parents=[common], help="Detect peaks (JSON) for chemistry assignment"
    )
    pp.add_argument("h5")
    pp.add_argument(
        "--full",
        action="store_true",
        help="Emit every candidate formula + isotope arrays per peak "
        "(default is a compact top-candidate view)",
    )
    pp.add_argument(
        "--include-artifacts",
        action="store_true",
        help="Include instrument-noise peaks (ringing combs, "
        "low-prominence ripples, reagent saturation skirt) that the "
        "default view drops",
    )
    pp.add_argument(
        "--min-height",
        type=float,
        default=1e-3,
        help="Threshold as fraction of tallest peak (default 1e-3)",
    )
    pp.add_argument("--max-peaks", type=int, default=300)
    pp.add_argument("--mz-min", type=float, default=15.0)
    pp.add_argument("--mz-max", type=float, default=None)
    pp.add_argument(
        "--R-phys",
        dest="R_phys",
        type=float,
        default=None,
        help="Physical peak resolution for detection and merging (default 2400)",
    )
    pp.set_defaults(func=cmd_peaks)

    ps = sub.add_parser(
        "segments", parents=[common], help="Detect time segments (JSON) for labelling"
    )
    ps.add_argument("h5")
    ps.add_argument(
        "--min-duration", type=int, default=30, help="Min cycles per segment"
    )
    ps.add_argument(
        "--grad-thr",
        type=float,
        default=0.02,
        help="Log-signal gradient threshold for stability",
    )
    ps.add_argument(
        "--high-ratio",
        type=float,
        default=3.0,
        help="x-baseline above which a segment is 'high' (sample)",
    )
    ps.add_argument(
        "--merge-high-gap",
        type=int,
        default=None,
        metavar="CYCLES",
        help="Force the gap cap when joining consecutive high plateaus. By default "
        "every gap is judged on evidence — it merges only if it stayed within its "
        "neighbours' level and covered under ~60 s of acquisition. 0 never merges "
        "high plateaus",
    )
    ps.set_defaults(func=cmd_segments)

    pa = sub.add_parser("analyze", parents=[common], help="Run pipeline -> results CSV")
    pa.add_argument("h5")
    pa.add_argument(
        "--peaks-json", help="Inline JSON: [{'mz':.., 'label':.., 'formula':..}]"
    )
    pa.add_argument(
        "--ranges-json",
        help="Inline JSON: [{'label':.., 'start':.., 'end':.., 'unit':'cycle|second'}]",
    )
    pa.add_argument(
        "--config", help="JSON file with 'peaks'/'ranges' (alternative to inline)"
    )
    pa.add_argument(
        "--auto-peaks", action="store_true", help="Auto-detect peaks if none given"
    )
    pa.add_argument(
        "--auto-segments",
        action="store_true",
        help="Auto-detect segments if no ranges given (generic labels)",
    )
    pa.add_argument(
        "--merge-high-gap",
        type=int,
        default=None,
        metavar="CYCLES",
        help="With --auto-segments, force the gap cap when joining consecutive "
        "high plateaus (default: judge each gap on the signal itself, across "
        "~60 s of acquisition; 0 never merges high plateaus)",
    )
    pa.add_argument(
        "--include-cycle-rows",
        action="store_true",
        help="Append Viewer-style Cycle rows with each range's boundaries",
    )
    pa.add_argument("--out", default="-", help="Output CSV path (default stdout)")
    pa.add_argument(
        "--sep",
        default=";",
        help="Delimiter (default ';' + comma decimals = "
        "PTR-MS Viewer format; pass ',' for a standard ','-delimited, "
        "dot-decimal CSV for other tools)",
    )
    pa.add_argument("--min-height", type=float, default=1e-3)
    pa.add_argument("--max-peaks", type=int, default=300)
    pa.add_argument("--mz-min", type=float, default=15.0)
    pa.add_argument("--mz-max", type=float, default=None)
    pa.add_argument(
        "--no-per-interval",
        dest="no_per_interval",
        action="store_true",
        default=None,
        help="Use one whole-run integration window per compound instead of "
        "re-centring each peak's window on every interval's own spectrum "
        "(per-interval is the default; peaks drift between intervals)",
    )
    pa.add_argument(
        "--per-interval",
        dest="no_per_interval",
        action="store_false",
        help="Explicitly use isolated per-interval windows",
    )
    pa.add_argument(
        "--R",
        type=float,
        default=None,
        help="Integration-window resolution (default 1200)",
    )
    pa.add_argument(
        "--R-phys",
        dest="R_phys",
        type=float,
        default=None,
        help="Physical peak resolution for deconvolution (default 2400)",
    )
    pa.add_argument(
        "--K",
        type=float,
        help="Concentration constant (Conc=Corrected*K/primary). "
        "Default: derived from file; use `calibrate` to match a Viewer project.",
    )
    pa.add_argument(
        "--primary-mz",
        type=float,
        default=None,
        help="Primary-ion m/z for normalisation (default 21.022, H3(18O)+)",
    )
    pa.add_argument(
        "--kinetic",
        dest="kinetic",
        action="store_true",
        default=None,
        help="Apply per-compound rate-constant (k) correction for physically "
        "resolved sensitivities (looks up k by peak 'k'/'formula'/m/z). "
        "Diverges from a single-k reference but is more accurate.",
    )
    pa.add_argument(
        "--no-kinetic",
        dest="kinetic",
        action="store_false",
        help="Disable config kinetic correction",
    )
    pa.add_argument(
        "--k-anchor",
        type=float,
        default=None,
        help="Rate constant (1e-9 cm3/s) the baseline K assumes (default 2.0)",
    )
    pa.add_argument(
        "--humidity-correct",
        dest="humidity_correct",
        action="store_true",
        default=None,
        help="Humidity-correct near-thermoneutral compounds (HCN etc.) using "
        "the per-cycle water-cluster ratio. Needs a calibrated --humidity-p.",
    )
    pa.add_argument(
        "--no-humidity-correct",
        dest="humidity_correct",
        action="store_false",
        help="Disable config humidity correction",
    )
    pa.add_argument(
        "--humidity-p",
        type=float,
        default=None,
        help="Humidity exponent in [0,1]: 0=off, 1=equilibrium upper bound "
        "(default 1.0). Calibrate from a standard at >=2 humidities.",
    )
    pa.add_argument(
        "--humidity-ref",
        type=float,
        default=None,
        help="Reference water-cluster ratio to normalise to (default: run median)",
    )
    pa.add_argument(
        "--molar-volume",
        type=float,
        default=None,
        help="Molar volume L/mol (else from drift temperature)",
    )
    pa.set_defaults(func=cmd_analyze)

    pv = sub.add_parser(
        "viz",
        parents=[common],
        help="Review app for an existing peak list + ranges (live-save "
        "to a config with --serve, or a standalone HTML with --html)",
    )
    pv.add_argument("h5")
    pv.add_argument(
        "--config",
        help="JSON file with 'peaks'/'ranges' (source; live-save target when serving)",
    )
    pv.add_argument("--peaks-json", help="Inline peaks JSON (alternative to --config)")
    pv.add_argument(
        "--ranges-json", help="Inline ranges JSON (alternative to --config)"
    )
    pv.add_argument(
        "--x-axis-unit",
        choices=_X_AXIS_UNITS,
        default=None,
        help="Time-trace x-axis: cycle, relative, or absolute (default cycle)",
    )
    pv.add_argument(
        "--save-config",
        help="Config path to live-save to when serving without --config",
    )
    pv.add_argument(
        "--serve",
        dest="serve",
        action="store_true",
        default=None,
        help="Serve on localhost, live-save edits, and run analysis on 'Done' (the default)",
    )
    pv.add_argument(
        "--html", help="Instead of serving, write a standalone portable HTML file here"
    )
    pv.add_argument(
        "--out",
        default="results.csv",
        help="Results CSV written when the expert clicks 'Done' (default results.csv)",
    )
    pv.add_argument(
        "--sep", default=";", help="CSV delimiter (default ';' with comma decimals)"
    )
    pv.add_argument(
        "--include-cycle-rows",
        dest="include_cycle_rows",
        action="store_true",
        default=True,
        help="Append Viewer-style Cycle rows (default on)",
    )
    pv.add_argument("--no-cycle-rows", dest="include_cycle_rows", action="store_false")
    pv.add_argument(
        "--port",
        type=int,
        default=8765,
        help="Server port (default 8765; scans upward if busy)",
    )
    pv.add_argument(
        "--timeout",
        type=int,
        default=None,
        help="Seconds to wait for 'Done' (default: indefinitely)",
    )
    pv.add_argument(
        "--no-open", action="store_true", help="Do not auto-open the browser"
    )
    pv.add_argument("--R", type=float, default=None)
    pv.add_argument("--R-phys", dest="R_phys", type=float, default=None)
    pv.add_argument(
        "--K",
        type=float,
        default=None,
        help="Initial concentration constant (default from file)",
    )
    pv.add_argument("--primary-mz", type=float, default=None)
    pv.add_argument("--kinetic", dest="kinetic", action="store_true", default=None)
    pv.add_argument("--no-kinetic", dest="kinetic", action="store_false")
    pv.add_argument("--k-anchor", type=float, default=None)
    pv.add_argument(
        "--humidity-correct", dest="humidity_correct", action="store_true", default=None
    )
    pv.add_argument(
        "--no-humidity-correct", dest="humidity_correct", action="store_false"
    )
    pv.add_argument("--humidity-p", type=float, default=None)
    pv.add_argument("--humidity-ref", type=float, default=None)
    pv.add_argument(
        "--no-per-interval", dest="no_per_interval", action="store_true", default=None
    )
    pv.add_argument("--per-interval", dest="no_per_interval", action="store_false")
    pv.add_argument("--molar-volume", type=float, default=None)
    pv.set_defaults(func=cmd_viz)

    pap = sub.add_parser(
        "app",
        help="Run the persistent review app; open files from its own start screen",
    )
    pap.add_argument("h5", nargs="?", help="Open this file immediately (optional)")
    pap.add_argument(
        "--port",
        type=int,
        default=8765,
        help="Port for the localhost server (default 8765)",
    )
    pap.add_argument(
        "--agent",
        default=None,
        help="URL to post a freshly detected config to for curation; the "
        "deterministic config is kept if the endpoint fails (or set SNIFF_AGENT_URL)",
    )
    pap.add_argument(
        "--agent-timeout",
        dest="agent_timeout",
        type=float,
        default=300.0,
        help="Seconds to wait for the agent endpoint (default 300)",
    )
    pap.add_argument(
        "--no-browser", action="store_true", help="Do not open a browser window"
    )
    pap.add_argument(
        "--window",
        action="store_true",
        help="Open the app in its own desktop window instead of a browser tab "
        "(needs the 'desktop' extra; a packaged bundle uses it by default)",
    )
    pap.set_defaults(func=cmd_app)

    pr = sub.add_parser(
        "rates",
        parents=[common],
        help="Look up proton-transfer rate constants (k) by name/formula/mz",
    )
    pr.add_argument(
        "query", nargs="?", help="Substring of name/formula, or an m/z number"
    )
    pr.set_defaults(func=cmd_rates)

    pk = sub.add_parser(
        "calibrate",
        parents=[common],
        help="Fit concentration constant K to a reference Viewer CSV",
    )
    pk.add_argument("h5")
    pk.add_argument(
        "reference", help="Reference PTR-MS Viewer CSV with known concentrations"
    )
    pk.add_argument(
        "--peaks-json",
        help="Restrict calibration to these peaks (else use reference's)",
    )
    pk.add_argument(
        "--ranges-json", help="Ranges (else recovered from the reference's Cycle rows)"
    )
    pk.add_argument("--config")
    pk.add_argument("--auto-peaks", action="store_true")
    pk.add_argument("--auto-segments", action="store_true")
    pk.add_argument("--min-height", type=float, default=1e-3)
    pk.add_argument("--max-peaks", type=int, default=300)
    pk.add_argument("--mz-min", type=float, default=15.0)
    pk.add_argument("--mz-max", type=float, default=None)
    pk.add_argument(
        "--primary-mz",
        type=float,
        default=None,
        help="Primary-ion m/z for normalisation (default 21.022, H3(18O)+)",
    )
    pk.add_argument(
        "--R",
        type=float,
        default=None,
        help="Integration-window resolution (default 1200)",
    )
    pk.add_argument(
        "--R-phys",
        dest="R_phys",
        type=float,
        default=None,
        help="Physical peak resolution for deconvolution (default 2400)",
    )
    pk.set_defaults(func=cmd_calibrate)

    pc = sub.add_parser(
        "compare", parents=[common], help="Compare results CSV vs reference Viewer CSV"
    )
    pc.add_argument("mine")
    pc.add_argument("reference")
    pc.add_argument("--per-mass", action="store_true")
    pc.set_defaults(func=cmd_compare)

    args = p.parse_args()
    try:
        args.func(args)
    except ptrms.MassCalibrationError as exc:
        _emit(
            {
                "error": "mass_calibration_failed",
                "message": str(exc),
                "mass_axis_calibration": exc.diagnostics,
            },
            getattr(args, "raw", True),
        )
        raise SystemExit(2)


if __name__ == "__main__":
    main()
