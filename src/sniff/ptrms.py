"""Generalizable PTR-MS reprocessing pipeline (open-source replacement for PTR-MS Viewer).

Everything instrument-specific is read from the HDF5 file itself:
  - mass calibration      <- CALdata/Mapping  (timebin = a*sqrt(m) + b)
  - transmission curve    <- PTR-Transmission
  - concentration factor  <- derived from pre-computed TRACEdata (Conc/Corrected)
  - molar volume Vm       <- drift temperature (AddTraces/PTR-Reaction)

User inputs (experiment-specific, not in the raw file):
  - target peak list (m/z of product ions to quantify)
  - time ranges (labelled cycle windows)
"""

import copy
import json
from importlib import resources

import numpy as np

from . import isotopes, peak_fit

PROTON = 1.007276

# Date.toISOString() switches to expanded years outside this UTC interval.
_JS_NORMAL_YEAR_0000_START_S = -62167219200.0
_JS_NORMAL_YEAR_10000_START_S = 253402300800.0

K_ANCHOR_DEFAULT = 2.0  # 1e-9 cm3/s: the single k a non-kinetic calibration assumes

# Run-internal reference ions. For authoritative multi-point Mapping files their
# centres are stability diagnostics only. Older two-point/Spectrum calibrations retain
# the separate affine correction in the mass domain.
INTERNAL_MASS_ANCHORS = (
    ("water_cluster", 37.033),
    ("iodobenzene", 204.951),
)
INTERNAL_MASS_CORRECTION_MODEL = "m_corrected = scale*m_file + offset"
MAPPING_AUTHORITY_MODEL = "authoritative-multi-point-mapping-v1"
MASS_AXIS_ALGORITHM_VERSION = 2
FILE_MASS_CALIBRATION_MODEL = "timebin = a*sqrt(m_file) + b"
INTERNAL_ANCHOR_SEARCH_DA = 0.20
INTERNAL_ANCHOR_MIN_PROMINENCE = 8.0
INTERNAL_ANCHOR_MIN_SNR = 8.0
INTERNAL_ANCHOR_AMBIGUITY_RATIO = 0.50
INTERNAL_ANCHOR_BLOCKS = 8
INTERNAL_ANCHOR_MIN_PERSISTENCE = 0.625
INTERNAL_ANCHOR_PROXIMITY_DA = 0.025
INTERNAL_MASS_SCALE_LIMIT = 0.005
INTERNAL_MASS_OFFSET_LIMIT_DA = 0.25

FORMULA_TOLERANCE_MODEL = "run-calibration-residuals-v1"
FORMULA_TOLERANCE_FLOOR_PPM = 5.0
FORMULA_TOLERANCE_MAX_PPM = 10.0
FORMULA_TOLERANCE_MARGIN_PPM = 2.0
FORMULA_PROPOSAL_TOLERANCE_PPM = 200.0
FORMULA_SCORE_SIGMA_FLOOR_PPM = 2.0


MASS_AXIS_CONFIG_DOMAIN = "corrected"
MASS_AXIS_CONFIG_VERSION = 2


class MassCalibrationError(ValueError):
    """Raised when a run cannot establish a trustworthy mass axis.

    Attributes:
        diagnostics:
            JSON-serialisable calibration evidence, including failed anchors.
    """

    def __init__(self, message, diagnostics):
        self.diagnostics = diagnostics
        super().__init__(f"mass calibration failed: {message}")


class MassAxisCalibration:
    """File timebin calibration plus an affine mass-domain correction."""

    def __init__(self, a, b, scale=1.0, offset=0.0, diagnostics=None):
        self.a = float(a)
        self.b = float(b)
        self.scale = float(scale)
        self.offset = float(offset)
        self.diagnostics = diagnostics or {}

    @property
    def applied(self):
        """Whether the Mapping or fallback internal calibration was accepted."""
        return bool(self.diagnostics.get("applied", False))

    def corrected_to_file(self, mass):
        """Map a corrected m/z value back onto the file's baseline mass domain."""
        return (np.asarray(mass) - self.offset) / self.scale

    def file_to_corrected(self, mass):
        """Map a baseline file m/z value onto the corrected mass domain."""
        return self.scale * np.asarray(mass) + self.offset

    def m_to_tb(self, mass):
        """Convert corrected m/z to a file timebin."""
        return self.a * np.sqrt(self.corrected_to_file(mass)) + self.b

    def tb_to_m(self, timebin):
        """Convert a file timebin to corrected m/z."""
        file_mass = ((np.asarray(timebin) - self.b) / self.a) ** 2
        return self.file_to_corrected(file_mass)

    def to_dict(self):
        """Return deterministic JSON-safe calibration provenance."""
        return dict(self.diagnostics)


def formula_assignment_tolerance(mass_axis, mz):
    """Return the run-validated ppm exact-mass tolerance at one ion m/z.

    A ppm radius naturally expands in mDa with m/z. Three or more independent file
    calibration references validate a 5-10 ppm run-specific radius. Files without
    independent residuals retain conservative review candidates at 10 ppm. When a
    run exceeds 10 ppm, all broad candidates remain reviewer proposals only.
    """
    validate_mass_axis(mass_axis)
    ion_mz = float(mz)
    if not np.isfinite(ion_mz) or ion_mz <= 0:
        raise ValueError("formula-assignment m/z must be finite and positive")
    model = mass_axis.diagnostics.get("formula_assignment_tolerance")
    if not isinstance(model, dict) or model.get("model") != FORMULA_TOLERANCE_MODEL:
        model = {
            "model": "legacy-conservative-fallback",
            "source": "saved calibration lacks independent residual evidence",
            "status": "fallback",
            "tolerance_ppm": FORMULA_TOLERANCE_MAX_PPM,
            "score_sigma_ppm": FORMULA_TOLERANCE_MAX_PPM / 2.5,
            "candidate_generation_allowed": True,
            "automatic_assignment_allowed": False,
        }
    tolerance_ppm = float(model["tolerance_ppm"])
    return {
        "model": model["model"],
        "source": model["source"],
        "status": model["status"],
        "reason": model.get("reason"),
        "ppm": tolerance_ppm,
        "mDa": tolerance_ppm * ion_mz / 1000.0,
        "proposal_ppm": FORMULA_PROPOSAL_TOLERANCE_PPM,
        "proposal_mDa": FORMULA_PROPOSAL_TOLERANCE_PPM * ion_mz / 1000.0,
        "score_sigma_ppm": float(model["score_sigma_ppm"]),
        "candidate_generation_allowed": bool(
            model["candidate_generation_allowed"]
        ),
        "automatic_assignment_allowed": bool(
            model["automatic_assignment_allowed"]
        ),
    }


def _derive_reference_stability(anchors, scale, offset):
    """Summarise temporal reference movement without treating it as mass accuracy."""
    output = []
    for anchor in anchors:
        target = float(anchor["target_mz"])
        persistence = anchor.get("persistence") or {}
        centres = [
            scale * float(centre) + offset
            for status, centre in zip(
                persistence.get("statuses", []),
                persistence.get("block_centres_file_mz", []),
            )
            if status == "accepted" and centre is not None
        ]
        if not centres:
            continue
        reference = float(np.median(centres))
        values = (np.asarray(centres, dtype=np.float64) - reference) / reference * 1e6
        median = float(np.median(values))
        output.append(
            {
                "name": anchor["name"],
                "mz": target,
                "n_blocks": len(centres),
                "median_ppm": median,
                "robust_sigma_ppm": float(
                    1.4826 * np.median(np.abs(values - median))
                ),
                "max_abs_ppm": float(np.max(np.abs(values))),
            }
        )
    return {
        "source": "accepted internal-reference block centres",
        "interpretation": "temporal stability diagnostic, not mass-accuracy evidence",
        "references": output,
    }


def _derive_formula_tolerance_model(f, a, b):
    """Validate a 5-10 ppm radius against independent Mapping residuals."""
    try:
        mapping = np.asarray(f["CALdata/Mapping"][:], dtype=np.float64)
    except (KeyError, OSError, TypeError, ValueError):
        mapping = None
    points = []
    if (
        mapping is not None
        and mapping.ndim == 2
        and mapping.shape[0] >= 3
        and mapping.shape[1] == 2
        and np.isfinite(mapping).all()
        and (mapping > 0).all()
    ):
        mapping = mapping[np.argsort(mapping[:, 0])]
        masses = mapping[:, 0]
        timebins = mapping[:, 1]
        if np.all(np.diff(masses) > 0) and np.all(np.diff(timebins) > 0):
            fitted_masses = ((timebins - b) / a) ** 2
            residuals = (fitted_masses - masses) / masses * 1e6
            if np.isfinite(residuals).all():
                points = [
                    {"mz": float(mass), "residual_ppm": float(residual)}
                    for mass, residual in zip(masses, residuals)
                ]
    if not points:
        return {
            "model": FORMULA_TOLERANCE_MODEL,
            "source": "no independent multi-point Mapping residuals",
            "status": "fallback",
            "reason": "fewer than three usable independent calibration references",
            "mass_error_convention": "1e6 * (observed - theoretical) / theoretical",
            "minimum_ppm": FORMULA_TOLERANCE_FLOOR_PPM,
            "maximum_ppm": FORMULA_TOLERANCE_MAX_PPM,
            "safety_margin_ppm": FORMULA_TOLERANCE_MARGIN_PPM,
            "tolerance_ppm": FORMULA_TOLERANCE_MAX_PPM,
            "score_sigma_ppm": FORMULA_TOLERANCE_MAX_PPM / 2.5,
            "candidate_generation_allowed": True,
            "automatic_assignment_allowed": False,
            "calibration_points": [],
        }
    absolute = np.abs([point["residual_ppm"] for point in points])
    q95_abs = float(np.percentile(absolute, 95.0))
    accepted = q95_abs <= FORMULA_TOLERANCE_MAX_PPM
    tolerance_ppm = (
        max(
            FORMULA_TOLERANCE_FLOOR_PPM,
            min(FORMULA_TOLERANCE_MAX_PPM, q95_abs + FORMULA_TOLERANCE_MARGIN_PPM),
        )
        if accepted
        else FORMULA_TOLERANCE_MAX_PPM
    )
    return {
        "model": FORMULA_TOLERANCE_MODEL,
        "source": "independent CALdata/Mapping fit residuals",
        "status": "accepted" if accepted else "degraded",
        "reason": (
            None
            if accepted
            else f"95th-percentile calibration residual {q95_abs:.2f} ppm exceeds "
            f"the {FORMULA_TOLERANCE_MAX_PPM:.0f} ppm assignment limit"
        ),
        "mass_error_convention": "1e6 * (observed - theoretical) / theoretical",
        "minimum_ppm": FORMULA_TOLERANCE_FLOOR_PPM,
        "maximum_ppm": FORMULA_TOLERANCE_MAX_PPM,
        "safety_margin_ppm": FORMULA_TOLERANCE_MARGIN_PPM,
        "q95_abs_ppm": q95_abs,
        "tolerance_ppm": tolerance_ppm,
        "score_sigma_ppm": max(
            FORMULA_SCORE_SIGMA_FLOOR_PPM, tolerance_ppm / 2.5
        ),
        "candidate_generation_allowed": True,
        "automatic_assignment_allowed": accepted,
        "calibration_points": points,
    }


def validate_mass_axis(mass_axis):
    """Validate the complete evidence for an accepted mass calibration.

    This is the trust boundary for axes passed between high-level operations. The
    ``applied`` flag is not evidence by itself: either authoritative multi-point Mapping
    evidence or the legacy two-reference affine evidence must agree with the numerical
    coefficients.

    Raises:
        MassCalibrationError:
            If the object or any part of its diagnostic evidence is malformed or
            contradictory.
    """
    diagnostics = getattr(mass_axis, "diagnostics", None)
    if not isinstance(mass_axis, MassAxisCalibration) or not isinstance(
        diagnostics, dict
    ):
        raise MassCalibrationError(
            "caller-supplied mass axis is not an internal calibration",
            {"applied": False, "fallback_reason": "invalid mass-axis object"},
        )

    def fail(reason):
        raise MassCalibrationError(reason, dict(diagnostics))

    def numeric(value):
        return isinstance(
            value, (int, float, np.integer, np.floating)
        ) and not isinstance(value, bool)

    try:
        values = np.asarray(
            [mass_axis.a, mass_axis.b, mass_axis.scale, mass_axis.offset],
            dtype=np.float64,
        )
    except (OverflowError, TypeError, ValueError):
        fail("caller-supplied mass axis has non-numeric coefficients")
    if not np.isfinite(values).all() or mass_axis.a <= 0 or mass_axis.scale <= 0:
        fail("caller-supplied mass axis has non-finite or non-physical coefficients")
    if diagnostics.get("applied") is not True:
        fail("caller-supplied mass axis is not an applied internal calibration")
    if diagnostics.get("fallback_reason") is not None:
        fail("caller-supplied mass axis retains a fallback reason")
    if diagnostics.get("model") not in {
        INTERNAL_MASS_CORRECTION_MODEL,
        MAPPING_AUTHORITY_MODEL,
    }:
        fail("caller-supplied mass axis has an unsupported correction model")

    try:
        diag_scale_value = diagnostics["scale"]
        diag_offset_value = diagnostics["offset_da"]
        file_cal = diagnostics["file_calibration"]
        if (
            not isinstance(file_cal, dict)
            or file_cal.get("model") != FILE_MASS_CALIBRATION_MODEL
        ):
            raise TypeError
        file_a_value = file_cal["a"]
        file_b_value = file_cal["b"]
        if not all(
            numeric(value)
            for value in (
                diag_scale_value,
                diag_offset_value,
                file_a_value,
                file_b_value,
            )
        ):
            raise TypeError
        diag_scale = float(diag_scale_value)
        diag_offset = float(diag_offset_value)
        file_a = float(file_a_value)
        file_b = float(file_b_value)
    except (KeyError, OverflowError, TypeError, ValueError):
        fail("caller-supplied mass axis has incomplete calibration diagnostics")
    if not np.isfinite([diag_scale, diag_offset, file_a, file_b]).all():
        fail("caller-supplied mass axis has non-finite calibration diagnostics")
    if file_a <= 0:
        fail("caller-supplied mass axis has non-physical file calibration")
    if (
        abs(mass_axis.scale - 1.0) > INTERNAL_MASS_SCALE_LIMIT
        or abs(mass_axis.offset) > INTERNAL_MASS_OFFSET_LIMIT_DA
    ):
        fail("caller-supplied mass axis has an implausible correction")
    if not np.allclose(
        [mass_axis.scale, mass_axis.offset, mass_axis.a, mass_axis.b],
        [diag_scale, diag_offset, file_a, file_b],
        rtol=0,
        atol=1e-12,
    ):
        fail("caller-supplied mass axis has contradictory calibration diagnostics")

    if diagnostics.get("model") == MAPPING_AUTHORITY_MODEL:
        try:
            mapping = diagnostics["mapping_calibration"]
            points = mapping["points"]
            if (
                diagnostics.get("algorithm_version") != MASS_AXIS_ALGORITHM_VERSION
                or diagnostics.get("authority") != "CALdata/Mapping"
                or diagnostics.get("mass_domain_correction_applied") is not False
                or mass_axis.scale != 1.0
                or mass_axis.offset != 0.0
                or mapping.get("source") != "CALdata/Mapping"
                or mapping.get("n_points") != len(points)
                or len(points) < 3
            ):
                raise TypeError
            masses = np.asarray([float(point["mz"]) for point in points])
            timebins = np.asarray([float(point["timebin"]) for point in points])
            residuals = np.asarray(
                [float(point["residual_ppm"]) for point in points]
            )
            expected = (((timebins - file_b) / file_a) ** 2 - masses) / masses * 1e6
            tolerance = diagnostics["formula_assignment_tolerance"]
            tolerance_points = tolerance["calibration_points"]
            tolerance_masses = np.asarray(
                [float(point["mz"]) for point in tolerance_points]
            )
            tolerance_residuals = np.asarray(
                [float(point["residual_ppm"]) for point in tolerance_points]
            )
            q95_abs = float(np.percentile(np.abs(residuals), 95.0))
            accepted = q95_abs <= FORMULA_TOLERANCE_MAX_PPM
            expected_tolerance = (
                max(
                    FORMULA_TOLERANCE_FLOOR_PPM,
                    min(
                        FORMULA_TOLERANCE_MAX_PPM,
                        q95_abs + FORMULA_TOLERANCE_MARGIN_PPM,
                    ),
                )
                if accepted
                else FORMULA_TOLERANCE_MAX_PPM
            )
            if (
                not np.isfinite([*masses, *timebins, *residuals]).all()
                or not (masses > 0).all()
                or not (timebins > 0).all()
                or not np.all(np.diff(masses) > 0)
                or not np.all(np.diff(timebins) > 0)
                or not np.allclose(residuals, expected, rtol=0, atol=1e-9)
                or np.any(
                    np.abs(residuals) > MAPPING_MAX_RELATIVE_MASS_ERROR * 1e6
                )
                or tolerance.get("model") != FORMULA_TOLERANCE_MODEL
                or tolerance.get("status")
                != ("accepted" if accepted else "degraded")
                or not np.isclose(
                    float(tolerance["q95_abs_ppm"]), q95_abs, rtol=0, atol=1e-12
                )
                or not np.isclose(
                    float(tolerance["tolerance_ppm"]),
                    expected_tolerance,
                    rtol=0,
                    atol=1e-12,
                )
                or not np.isclose(
                    float(tolerance["score_sigma_ppm"]),
                    max(
                        FORMULA_SCORE_SIGMA_FLOOR_PPM,
                        expected_tolerance / 2.5,
                    ),
                    rtol=0,
                    atol=1e-12,
                )
                or tolerance.get("candidate_generation_allowed") is not True
                or tolerance.get("automatic_assignment_allowed") is not accepted
                or not np.allclose(tolerance_masses, masses, rtol=0, atol=1e-9)
                or not np.allclose(
                    tolerance_residuals, residuals, rtol=0, atol=1e-9
                )
            ):
                raise TypeError
        except (KeyError, OverflowError, TypeError, ValueError):
            fail("authoritative Mapping calibration evidence is contradictory")
        return mass_axis

    required = {name: float(target) for name, target in INTERNAL_MASS_ANCHORS}
    anchors = diagnostics.get("anchors")
    if not isinstance(anchors, list) or len(anchors) != len(required):
        fail("caller-supplied mass axis must contain exactly the required anchors")
    if any(not isinstance(anchor, dict) for anchor in anchors):
        fail("caller-supplied mass axis contains malformed anchor evidence")
    names = [anchor.get("name") for anchor in anchors]
    if (
        any(not isinstance(name, str) for name in names)
        or len(set(names)) != len(names)
        or set(names) != set(required)
    ):
        fail("caller-supplied mass axis contains contradictory anchor names")

    min_blocks = 2
    for anchor in anchors:
        name = anchor["name"]
        target = required[name]
        try:
            target_raw = anchor["target_mz"]
            observed_raw = anchor["observed_file_mz"]
            corrected_raw = anchor["corrected_mz"]
            timebin_raw = anchor["timebin"]
            prominence_raw = anchor["prominence"]
            snr_raw = anchor["snr"]
            persistence = anchor["persistence"]
            if not all(
                numeric(value)
                for value in (
                    target_raw,
                    observed_raw,
                    corrected_raw,
                    timebin_raw,
                    prominence_raw,
                    snr_raw,
                )
            ):
                raise TypeError
            target_value = float(target_raw)
            observed = float(observed_raw)
            corrected = float(corrected_raw)
            timebin = float(timebin_raw)
            prominence = float(prominence_raw)
            snr = float(snr_raw)
        except (KeyError, OverflowError, TypeError, ValueError):
            fail(f"{name} anchor evidence is incomplete")
        if (
            anchor.get("status") != "accepted"
            or anchor.get("reason") != ""
            or not np.isfinite(
                [target_value, observed, corrected, timebin, prominence, snr]
            ).all()
            or not np.isclose(target_value, target, rtol=0, atol=1e-12)
            or observed <= 0
            or (
                INTERNAL_ANCHOR_SEARCH_DA - abs(observed - target)
                < INTERNAL_ANCHOR_PROXIMITY_DA
            )
            or not np.isclose(
                timebin,
                mass_axis.a * np.sqrt(observed) + mass_axis.b,
                rtol=0,
                atol=1e-6,
            )
            or prominence < INTERNAL_ANCHOR_MIN_PROMINENCE
            or snr < INTERNAL_ANCHOR_MIN_SNR
            or not np.isclose(
                corrected,
                mass_axis.scale * observed + mass_axis.offset,
                rtol=0,
                atol=1e-9,
            )
            or not np.isclose(corrected, target, rtol=0, atol=1e-9)
        ):
            fail(f"{name} anchor is not a valid corrected endpoint")
        if (
            not isinstance(persistence, dict)
            or persistence.get("available") is not True
        ):
            fail(f"{name} anchor persistence is unavailable")
        if persistence.get("checked", True) is not True:
            fail(f"{name} anchor persistence was not checked")
        try:
            blocks_value = persistence["blocks"]
            accepted_value = persistence["accepted_blocks"]
            fraction_value = persistence["fraction"]
            statuses = persistence["statuses"]
            if (
                isinstance(blocks_value, bool)
                or not isinstance(blocks_value, (int, np.integer))
                or isinstance(accepted_value, bool)
                or not isinstance(accepted_value, (int, np.integer))
                or isinstance(fraction_value, bool)
                or not isinstance(fraction_value, (float, int, np.floating, np.integer))
            ):
                raise TypeError
            blocks = int(blocks_value)
            accepted_blocks = int(accepted_value)
            fraction = float(fraction_value)
        except (KeyError, OverflowError, TypeError, ValueError):
            fail(f"{name} anchor persistence evidence is incomplete")
        if (
            blocks < min_blocks
            or blocks > INTERNAL_ANCHOR_BLOCKS
            or accepted_blocks < 0
            or accepted_blocks > blocks
            or not isinstance(statuses, list)
            or len(statuses) != blocks
            or any(
                not isinstance(status, str)
                or status
                not in {"accepted", "missing", "weak", "ambiguous", "implausible"}
                for status in statuses
            )
            or statuses.count("accepted") != accepted_blocks
            or not np.isfinite(fraction)
            or not np.isclose(fraction, accepted_blocks / blocks, rtol=0, atol=1e-12)
            or fraction < INTERNAL_ANCHOR_MIN_PERSISTENCE
        ):
            fail(f"{name} anchor persistence evidence is contradictory")

    tolerance_model = diagnostics.get("formula_assignment_tolerance")
    stability = diagnostics.get("reference_stability")
    if tolerance_model is not None or stability is not None:
        try:
            for anchor in anchors:
                persistence = anchor["persistence"]
                statuses = persistence["statuses"]
                centres = persistence["block_centres_file_mz"]
                if not isinstance(centres, list) or len(centres) != len(statuses):
                    raise TypeError
                for status, centre in zip(statuses, centres):
                    if status == "accepted":
                        if not numeric(centre) or not np.isfinite(centre) or centre <= 0:
                            raise TypeError
                    elif centre is not None:
                        raise TypeError
            expected_stability = _derive_reference_stability(
                anchors, mass_axis.scale, mass_axis.offset
            )
        except (KeyError, OverflowError, TypeError, ValueError):
            fail("formula-assignment calibration evidence is malformed")
        if stability is not None and stability != expected_stability:
            fail("reference-stability evidence is contradictory")
    if tolerance_model is not None:
        try:
            if tolerance_model["model"] != FORMULA_TOLERANCE_MODEL:
                raise TypeError
            status = tolerance_model["status"]
            tolerance_ppm = float(tolerance_model["tolerance_ppm"])
            sigma_ppm = float(tolerance_model["score_sigma_ppm"])
            points = tolerance_model["calibration_points"]
            if (
                status not in {"accepted", "degraded", "fallback"}
                or not isinstance(points, list)
                or not np.isfinite([tolerance_ppm, sigma_ppm]).all()
                or not FORMULA_TOLERANCE_FLOOR_PPM
                <= tolerance_ppm
                <= FORMULA_TOLERANCE_MAX_PPM
                or sigma_ppm <= 0
            ):
                raise TypeError
            if status == "fallback":
                valid = (
                    not points
                    and tolerance_ppm == FORMULA_TOLERANCE_MAX_PPM
                    and sigma_ppm
                    == FORMULA_TOLERANCE_MAX_PPM / 2.5
                    and tolerance_model["candidate_generation_allowed"] is True
                    and tolerance_model["automatic_assignment_allowed"] is False
                )
            else:
                residuals = np.asarray(
                    [float(point["residual_ppm"]) for point in points],
                    dtype=np.float64,
                )
                masses = np.asarray(
                    [float(point["mz"]) for point in points], dtype=np.float64
                )
                q95_abs = float(np.percentile(np.abs(residuals), 95.0))
                accepted = q95_abs <= FORMULA_TOLERANCE_MAX_PPM
                expected_tolerance = (
                    max(
                        FORMULA_TOLERANCE_FLOOR_PPM,
                        min(
                            FORMULA_TOLERANCE_MAX_PPM,
                            q95_abs + FORMULA_TOLERANCE_MARGIN_PPM,
                        ),
                    )
                    if accepted
                    else FORMULA_TOLERANCE_MAX_PPM
                )
                valid = (
                    len(points) >= 3
                    and np.isfinite(residuals).all()
                    and np.isfinite(masses).all()
                    and (masses > 0).all()
                    and np.all(np.diff(masses) > 0)
                    and np.isclose(
                        float(tolerance_model["q95_abs_ppm"]),
                        q95_abs,
                        rtol=0,
                        atol=1e-12,
                    )
                    and np.isclose(tolerance_ppm, expected_tolerance, rtol=0, atol=1e-12)
                    and np.isclose(
                        sigma_ppm,
                        max(FORMULA_SCORE_SIGMA_FLOOR_PPM, tolerance_ppm / 2.5),
                        rtol=0,
                        atol=1e-12,
                    )
                    and (status == "accepted") == accepted
                    and tolerance_model["candidate_generation_allowed"] is True
                    and tolerance_model["automatic_assignment_allowed"] is accepted
                )
            if not valid:
                raise TypeError
        except (KeyError, OverflowError, TypeError, ValueError):
            fail("formula-assignment tolerance evidence is contradictory")
    return mass_axis


def mass_axis_from_dict(diagnostics):
    """Rebuild and validate an internal calibration from saved diagnostics."""
    if not isinstance(diagnostics, dict):
        raise MassCalibrationError(
            "saved mass-axis calibration is not a mapping",
            {"applied": False, "fallback_reason": "invalid saved calibration"},
        )
    try:
        file_calibration = diagnostics["file_calibration"]
        calibration = MassAxisCalibration(
            file_calibration["a"],
            file_calibration["b"],
            scale=diagnostics["scale"],
            offset=diagnostics["offset_da"],
            diagnostics=copy.deepcopy(diagnostics),
        )
    except (KeyError, OverflowError, TypeError, ValueError) as exc:
        raise MassCalibrationError(
            "saved mass-axis calibration has incomplete coefficients",
            dict(diagnostics),
        ) from exc
    return validate_mass_axis(calibration)


def _historical_axis_for_migration(config, mass_axis):
    """Reconstruct the version-1 affine axis without applying version-2 rules."""
    diagnostics = config.get("mass_axis_calibration")
    if isinstance(diagnostics, dict):
        try:
            file_calibration = diagnostics["file_calibration"]
            scale = float(diagnostics["scale"])
            offset = float(diagnostics["offset_da"])
            a = float(file_calibration["a"])
            b = float(file_calibration["b"])
            anchors = diagnostics["anchors"]
            required = {name: float(target) for name, target in INTERNAL_MASS_ANCHORS}
            if (
                diagnostics.get("model") != INTERNAL_MASS_CORRECTION_MODEL
                or diagnostics.get("applied") is not True
                or file_calibration.get("model") != FILE_MASS_CALIBRATION_MODEL
                or not np.isfinite([a, b, scale, offset]).all()
                or a <= 0
                or scale <= 0
                or abs(scale - 1.0) > INTERNAL_MASS_SCALE_LIMIT
                or abs(offset) > INTERNAL_MASS_OFFSET_LIMIT_DA
                or not np.allclose(
                    [a, b], [mass_axis.a, mass_axis.b], rtol=0, atol=1e-9
                )
                or not isinstance(anchors, list)
                or {anchor.get("name") for anchor in anchors} != set(required)
            ):
                raise TypeError
            for anchor in anchors:
                observed = float(anchor["observed_file_mz"])
                target = required[anchor["name"]]
                if (
                    anchor.get("status") != "accepted"
                    or not np.isclose(
                        scale * observed + offset, target, rtol=0, atol=1e-8
                    )
                ):
                    raise TypeError
            return MassAxisCalibration(a, b, scale=scale, offset=offset)
        except (KeyError, OverflowError, TypeError, ValueError):
            pass

    anchors = mass_axis.diagnostics.get("anchors")
    required = {name: float(target) for name, target in INTERNAL_MASS_ANCHORS}
    try:
        if (
            not isinstance(anchors, list)
            or {anchor.get("name") for anchor in anchors} != set(required)
        ):
            raise TypeError
        observed = {}
        for anchor in anchors:
            persistence = anchor.get("persistence") or {}
            if (
                anchor.get("status") != "accepted"
                or persistence.get("available") is not True
                or float(persistence.get("fraction", 0.0))
                < INTERNAL_ANCHOR_MIN_PERSISTENCE
            ):
                raise TypeError
            observed[anchor["name"]] = float(anchor["observed_file_mz"])
        observed_lo = observed[INTERNAL_MASS_ANCHORS[0][0]]
        observed_hi = observed[INTERNAL_MASS_ANCHORS[1][0]]
        target_lo = INTERNAL_MASS_ANCHORS[0][1]
        target_hi = INTERNAL_MASS_ANCHORS[1][1]
        scale = (target_hi - target_lo) / (observed_hi - observed_lo)
        offset = target_lo - scale * observed_lo
        if (
            not np.isfinite([scale, offset]).all()
            or scale <= 0
            or abs(scale - 1.0) > INTERNAL_MASS_SCALE_LIMIT
            or abs(offset) > INTERNAL_MASS_OFFSET_LIMIT_DA
        ):
            raise TypeError
    except (KeyError, OverflowError, TypeError, ValueError, ZeroDivisionError) as exc:
        raise ValueError(
            "version-1 corrected config lacks reconstructable historical mass-axis "
            "evidence"
        ) from exc
    return MassAxisCalibration(
        mass_axis.a,
        mass_axis.b,
        scale=scale,
        offset=offset,
    )


def migrate_config_mass_axis(config, mass_axis):
    """Migrate an unmarked or version-1 config onto the current mass axis.

    The migration is deliberately limited to the documented config schema. Unknown
    fields are copied unchanged so an agent's provenance and review notes survive.
    A marked config is returned untouched, making repeated opens a no-op.

    Returns:
        A ``(config, migrated)`` pair. ``migrated`` is true when the marker was added.

    Raises:
        ValueError:
            If a config carries an unsupported mass-axis marker.
    """
    validate_mass_axis(mass_axis)
    result = copy.deepcopy(config)
    domain = result.get("mass_axis_domain")
    version = result.get("mass_axis_version")
    old_axis = None
    if domain is not None or version is not None:
        if domain != MASS_AXIS_CONFIG_DOMAIN or version not in {
            1,
            MASS_AXIS_CONFIG_VERSION,
        }:
            raise ValueError(
                "unsupported config mass axis marker: "
                f"domain={domain!r}, version={version!r}"
            )
        if version == MASS_AXIS_CONFIG_VERSION:
            return result, False
        old_axis = _historical_axis_for_migration(result, mass_axis)

    def corrected(value):
        file_mass = (
            old_axis.corrected_to_file(float(value))
            if old_axis is not None
            else float(value)
        )
        return float(mass_axis.file_to_corrected(file_mass))

    def width(value):
        old_scale = old_axis.scale if old_axis is not None else 1.0
        return float(value) * mass_axis.scale / old_scale

    for peak in result.get("peaks", []):
        if not isinstance(peak, dict):
            continue
        if "mz" in peak:
            peak["mz"] = corrected(peak["mz"])
        # These fields are emitted by older review payloads and are absolute m/z
        # coordinates, unlike cycle ranges and labels.
        for key in ("apex", "mass", "center"):
            if key in peak:
                peak[key] = corrected(peak[key])
        if "window" in peak:
            window = peak["window"]
            if isinstance(window, dict):
                for key in ("left", "right"):
                    if key in window:
                        window[key] = width(window[key])
            elif isinstance(window, (int, float)):
                peak["window"] = width(window)
        for key in ("win_l", "win_r"):
            if key in peak:
                peak[key] = width(peak[key])

    analysis = result.get("analyze")
    if isinstance(analysis, dict) and "primary_mz" in analysis:
        analysis["primary_mz"] = corrected(analysis["primary_mz"])
    result["mass_axis_domain"] = MASS_AXIS_CONFIG_DOMAIN
    result["mass_axis_version"] = MASS_AXIS_CONFIG_VERSION
    result["mass_axis_calibration"] = mass_axis.to_dict()
    return result, True


class AnalysisCancelled(Exception):
    """Raised when a caller's ``should_stop`` callback asks for a pass to stop.

    It is a user action, not a failure: an open that raises it leaves nothing
    behind and says nothing about the file."""


# 100 ppm: a deliberately generous corruption/model-consistency ceiling, not an
# accuracy claim; the real Data_10_26_33 fixture is about 8 ppm.
MAPPING_MAX_RELATIVE_MASS_ERROR = 100e-6


# ---------- per-compound rate constants (kinetic sensitivity) ----------
def load_rate_constants(path=None):
    """Load the bundled proton-transfer rate-constant table (or None).

    The table is loaded from package resources when no path is supplied. An
    explicit path remains available for custom or regenerated tables.
    """
    try:
        if path is None:
            resource = (
                resources.files("sniff")
                .joinpath("reference")
                .joinpath("rate_constants.json")
            )
            with resource.open("r", encoding="utf-8") as fh:
                return json.load(fh)
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, TypeError, ValueError):
        return None


def resolve_k(peaks, rate_table, mz_tol=0.03):
    """Determine each peak's rate constant k (in 1e-9 cm3/s units) and its source.

    Rate-constant priority is an explicit `k` on the peak > exact `formula` match >
    a unique m/z match in the table. Chemical flags are resolved independently from
    an exact formula match, otherwise a unique m/z match, without replacing an
    explicit k or its provenance. Returns {mz: {k, source, flags, k_estimated}}.
    `flags` carries 'humid' (needs humidity handling) / 'frag' (fragments) from the
    table.
    """
    by_formula, by_mz = {}, {}
    if rate_table:
        for c in rate_table.get("compounds", []):
            by_formula[c["formula"].upper()] = c
            by_mz.setdefault(round(c["mz"], 1), []).append(c)
    out = {}
    for p in peaks:
        mz = float(p["mz"])
        formula_match = (
            by_formula.get(p["formula"].upper()) if p.get("formula") else None
        )
        mz_matches = [
            c for c in by_mz.get(round(mz, 1), []) if abs(c["mz"] - mz) < mz_tol
        ]
        flag_match = formula_match or (mz_matches[0] if len(mz_matches) == 1 else None)
        flags = list(flag_match.get("flags", [])) if flag_match else []

        k, src, kest = None, None, False
        if p.get("k") is not None:
            k = float(p["k"])
            k = k / 1e-9 if k < 1e-6 else k  # accept SI or 1e-9 units
            src = "explicit"
            kest = bool(p.get("k_estimated", False))
        elif formula_match:
            k, src = formula_match["k"], "formula:" + formula_match["name"]
            kest = bool(formula_match.get("k_estimated", False))
        elif len(mz_matches) == 1:
            k, src = mz_matches[0]["k"], "mz:" + mz_matches[0]["name"]
            kest = bool(mz_matches[0].get("k_estimated", False))
        elif len(mz_matches) > 1:
            src = "ambiguous:" + ",".join(c["name"] for c in mz_matches)
        out[mz] = {"k": k, "source": src, "flags": flags, "k_estimated": kest}
    return out


# ---------- calibration read from file ----------
def load_mass_cal(f):
    """Mass calibration coefficients for ``timebin = a*sqrt(m) + b``.

    ``CALdata/Mapping`` is preferred when it contains a valid set of anchors.
    Exactly two anchors retain the original scalar calculation; additional anchors
    are sorted and fit by least squares after strict physical and numerical
    validation, including a generous reconstructed-mass residual sanity check.
    Raw-acquisition exports omit Mapping but store per-cycle ``[a, b]``
    coefficients directly in ``CALdata/Spectrum``; fall back to their median across
    cycles (robust to drift and zero/blank rows).
    """
    if "CALdata/Mapping" in f:
        try:
            raw_mapping = np.asarray(f["CALdata/Mapping"][:])
        except (TypeError, ValueError, OSError):
            raw_mapping = None

        if raw_mapping is not None and raw_mapping.shape == (2, 2):
            try:
                (m1, tb1), (m2, tb2) = raw_mapping
                if np.isfinite(raw_mapping).all() and m1 > 0 and m2 > 0 and m1 != m2:
                    # This is intentionally the original expression on the raw
                    # scalar dtype.  Some valid files use float32 and changing the
                    # order or precision changes their legacy calibration exactly.
                    a = (tb2 - tb1) / (np.sqrt(m2) - np.sqrt(m1))
                    b = tb1 - a * np.sqrt(m1)
                    if np.isfinite(a) and np.isfinite(b) and a > 0:
                        return float(a), float(b)
            except (TypeError, ValueError, FloatingPointError):
                pass

        if (
            raw_mapping is not None
            and raw_mapping.ndim == 2
            and raw_mapping.shape[1] == 2
            and raw_mapping.shape[0] > 2
        ):
            try:
                # Do not filter individual rows: one bad anchor invalidates the
                # Mapping rather than making an inconsistent fit look plausible.
                if np.isfinite(raw_mapping).all() and (raw_mapping > 0).all():
                    anchors = np.asarray(raw_mapping, dtype=np.float64)
                    anchors = anchors[np.argsort(anchors[:, 0])]
                    masses = anchors[:, 0]
                    timebins = anchors[:, 1]
                    if np.all(np.diff(masses) > 0) and np.all(np.diff(timebins) > 0):
                        sqrt_masses = np.sqrt(masses)
                        design = np.column_stack(
                            (sqrt_masses, np.ones_like(sqrt_masses))
                        )
                        if np.linalg.matrix_rank(design) == 2:
                            condition = np.linalg.cond(design)
                            # Limit round-off amplification to sqrt(epsilon), a
                            # standard useful-digit criterion for a float64 fit.
                            condition_limit = 1.0 / np.sqrt(np.finfo(np.float64).eps)
                            if np.isfinite(condition) and condition <= condition_limit:
                                a, b = np.linalg.lstsq(design, timebins, rcond=None)[0]
                                if np.isfinite(a) and np.isfinite(b) and a > 0:
                                    with np.errstate(over="ignore", invalid="ignore"):
                                        inferred_masses = ((timebins - b) / a) ** 2
                                        relative_mass_errors = np.abs(
                                            (inferred_masses - masses) / masses
                                        )
                                    if np.isfinite(
                                        relative_mass_errors
                                    ).all() and np.all(
                                        relative_mass_errors
                                        <= MAPPING_MAX_RELATIVE_MASS_ERROR
                                    ):
                                        return float(a), float(b)
            except (TypeError, ValueError, FloatingPointError, np.linalg.LinAlgError):
                pass

    if "CALdata/Spectrum" in f:
        try:
            sp = np.asarray(f["CALdata/Spectrum"][:], dtype=np.float64)
        except (TypeError, ValueError, OSError):
            sp = None
        if sp is not None and sp.ndim == 2 and sp.shape[1] == 2 and sp.shape[0] > 0:
            good = np.isfinite(sp).all(axis=1) & (sp[:, 0] > 0)
            if good.any():
                a, b = np.median(sp[good], axis=0)
                if np.isfinite(a) and np.isfinite(b) and a > 0:
                    return float(a), float(b)
    raise ValueError(
        "no mass calibration in file (neither CALdata/Mapping nor a usable "
        "CALdata/Spectrum)"
    )


def _authoritative_mapping_evidence(f, a, b):
    """Return validated multi-point Mapping evidence matching ``a,b``, or None."""
    if "CALdata/Mapping" not in f:
        return None
    try:
        mapping = np.asarray(f["CALdata/Mapping"][:], dtype=np.float64)
    except (OSError, TypeError, ValueError):
        return None
    if (
        mapping.ndim != 2
        or mapping.shape[0] < 3
        or mapping.shape[1] != 2
        or not np.isfinite(mapping).all()
        or not (mapping > 0).all()
    ):
        return None
    mapping = mapping[np.argsort(mapping[:, 0])]
    masses = mapping[:, 0]
    timebins = mapping[:, 1]
    if not np.all(np.diff(masses) > 0) or not np.all(np.diff(timebins) > 0):
        return None
    inferred = ((timebins - float(b)) / float(a)) ** 2
    residuals = (inferred - masses) / masses * 1e6
    if (
        not np.isfinite(residuals).all()
        or np.any(np.abs(residuals) > MAPPING_MAX_RELATIVE_MASS_ERROR * 1e6)
    ):
        return None
    return {
        "source": "CALdata/Mapping",
        "n_points": int(len(masses)),
        "points": [
            {
                "mz": float(mass),
                "timebin": float(timebin),
                "residual_ppm": float(residual),
            }
            for mass, timebin, residual in zip(masses, timebins, residuals)
        ],
    }


def load_mass_axis(f, *, progress=None, should_stop=None):
    """Build the mass axis from authoritative Mapping or fallback references.

    Args:
        progress (optional): Callback receiving a monotonic 0..1 calibration fraction.
        should_stop (optional): Callback polled between raw-cycle blocks.

    A valid Mapping with three or more points is authoritative and receives no second
    mass-domain translation or scaling. Water-cluster and iodobenzene movement remains
    diagnostic on that path. Exactly two Mapping points or Spectrum fallback retain the
    mandatory two-reference ``m_corrected = scale*m_file + offset`` correction.
    """
    if progress is not None:
        progress(0.0)
    if should_stop is not None and should_stop():
        raise AnalysisCancelled("the analysis was cancelled")
    try:
        a, b = load_mass_cal(f)
    except ValueError as exc:
        diagnostics = {
            "model": INTERNAL_MASS_CORRECTION_MODEL,
            "applied": False,
            "scale": 1.0,
            "offset_da": 0.0,
            "fallback_reason": str(exc),
            "anchors": [
                {
                    "name": name,
                    "target_mz": target,
                    "status": "unavailable",
                    "reason": "file mass calibration is unavailable",
                }
                for name, target in INTERNAL_MASS_ANCHORS
            ],
        }
        raise MassCalibrationError(str(exc), diagnostics) from exc
    mapping_evidence = _authoritative_mapping_evidence(f, a, b)
    base = {
        "model": INTERNAL_MASS_CORRECTION_MODEL,
        "algorithm_version": MASS_AXIS_ALGORITHM_VERSION,
        "authority": (
            "CALdata/Mapping" if mapping_evidence is not None else "internal-affine"
        ),
        "applied": False,
        "scale": 1.0,
        "offset_da": 0.0,
        "fallback_reason": None,
        "file_calibration": {
            "model": FILE_MASS_CALIBRATION_MODEL,
            "a": a,
            "b": b,
        },
        "anchors": [
            {
                "name": name,
                "target_mz": float(target_mz),
                "status": "unavailable",
                "reason": "anchor search did not run",
            }
            for name, target_mz in INTERNAL_MASS_ANCHORS
        ],
    }
    try:
        raw = np.asarray(f["SPECdata/AverageSpec"][:], dtype=np.float64)
        if progress is not None:
            progress(0.1)
    except (KeyError, OSError, TypeError, ValueError) as exc:
        _mass_axis_failure(
            a, b, base, f"average spectrum is unavailable or malformed: {exc}"
        )
    if raw.ndim != 1 or raw.size < 3:
        _mass_axis_failure(
            a,
            b,
            base,
            "average spectrum is malformed: expected a one-dimensional array with "
            "at least three bins",
        )
    finite = np.isfinite(raw)
    if not finite.any():
        _mass_axis_failure(
            a, b, base, "average spectrum is malformed: it has no finite bins"
        )
    avg = np.where(finite & (raw > 0), raw, 0.0)

    anchors = [
        _detect_internal_anchor(
            avg=avg,
            a=a,
            b=b,
            name=name,
            target_mz=target_mz,
        )
        for name, target_mz in INTERNAL_MASS_ANCHORS
    ]
    _add_anchor_persistence(
        f,
        avg,
        a,
        b,
        anchors,
        progress=progress,
        should_stop=should_stop,
    )
    base["anchors"] = anchors
    if mapping_evidence is not None:
        base.update(
            {
                "model": MAPPING_AUTHORITY_MODEL,
                "applied": True,
                "scale": 1.0,
                "offset_da": 0.0,
                "fallback_reason": None,
                "mass_domain_correction_applied": False,
                "mapping_calibration": mapping_evidence,
                "reference_stability": _derive_reference_stability(anchors, 1.0, 0.0),
                "formula_assignment_tolerance": _derive_formula_tolerance_model(f, a, b),
            }
        )
        calibration = MassAxisCalibration(a, b, diagnostics=base)
        validate_mass_axis(calibration)
        if progress is not None:
            progress(1.0)
        return calibration

    failures = [
        anchor
        for anchor in anchors
        if anchor["status"] != "accepted"
        or anchor.get("persistence", {}).get("available") is not True
    ]
    if failures:
        reasons = "; ".join(
            f"{anchor['name']} anchor {anchor['status']}: "
            f"{anchor['reason'] if anchor['status'] != 'accepted' else anchor.get('persistence', {}).get('reason', anchor['reason'])}"
            for anchor in failures
        )
        _mass_axis_failure(a, b, base, reasons)

    observed_lo = anchors[0]["observed_file_mz"]
    observed_hi = anchors[1]["observed_file_mz"]
    target_lo = anchors[0]["target_mz"]
    target_hi = anchors[1]["target_mz"]
    if observed_hi <= observed_lo:
        _mass_axis_failure(
            a,
            b,
            base,
            "internal anchors are physically implausible: mass order reversed",
        )
    scale = (target_hi - target_lo) / (observed_hi - observed_lo)
    offset = target_lo - scale * observed_lo
    if not np.isfinite([scale, offset]).all() or scale <= 0:
        _mass_axis_failure(
            a, b, base, "internal affine correction is non-finite or non-monotonic"
        )
    if abs(scale - 1.0) > INTERNAL_MASS_SCALE_LIMIT:
        _mass_axis_failure(
            a,
            b,
            base,
            "internal affine correction is implausible: "
            f"scale {scale:.9f} differs from unity by more than "
            f"{INTERNAL_MASS_SCALE_LIMIT * 100:.2f}%",
        )
    if abs(offset) > INTERNAL_MASS_OFFSET_LIMIT_DA:
        _mass_axis_failure(
            a,
            b,
            base,
            "internal affine correction is implausible: "
            f"offset {offset:+.6f} Da exceeds "
            f"{INTERNAL_MASS_OFFSET_LIMIT_DA:.2f} Da",
        )

    base.update(
        {
            "applied": True,
            "scale": float(scale),
            "offset_da": float(offset),
            "fallback_reason": None,
            "mass_domain_correction_applied": True,
        }
    )
    for anchor in anchors:
        anchor["corrected_mz"] = float(scale * anchor["observed_file_mz"] + offset)
    base["reference_stability"] = _derive_reference_stability(
        anchors, scale, offset
    )
    base["formula_assignment_tolerance"] = _derive_formula_tolerance_model(f, a, b)
    calibration = MassAxisCalibration(
        a, b, scale=scale, offset=offset, diagnostics=base
    )
    validate_mass_axis(calibration)
    if progress is not None:
        progress(1.0)
    return calibration


def _mass_axis_failure(a, b, diagnostics, reason):
    diagnostics.update(
        {
            "applied": False,
            "scale": 1.0,
            "offset_da": 0.0,
            "fallback_reason": reason,
        }
    )
    raise MassCalibrationError(reason, diagnostics)


def _add_anchor_persistence(f, avg, a, b, anchors, *, progress=None, should_stop=None):
    """Check both anchors in one deterministic raw-cycle persistence scan."""
    dataset = f.get("SPECdata/Intensities")
    reason = None
    if dataset is None or len(getattr(dataset, "shape", ())) != 2:
        reason = "raw cycle spectra are unavailable or malformed"
    elif int(dataset.shape[1]) != int(avg.size):
        reason = "raw cycle spectra do not match the average spectrum"
    elif int(dataset.shape[0]) < 2:
        reason = "fewer than two raw cycle spectra are available"
    if reason is not None:
        for anchor in anchors:
            anchor["persistence"] = {"available": False, "reason": reason}
        return

    ncyc = int(dataset.shape[0])
    block_count = min(INTERNAL_ANCHOR_BLOCKS, ncyc)
    edges = np.linspace(0, ncyc, block_count + 1, dtype=int)
    statuses = {anchor["name"]: [] for anchor in anchors}
    block_centres = {anchor["name"]: [] for anchor in anchors}
    windows = {}
    for anchor in anchors:
        target = float(anchor["target_mz"])
        tlo = max(1, int(np.floor(m_to_tb(target - INTERNAL_ANCHOR_SEARCH_DA, a, b))))
        thi = min(
            int(dataset.shape[1]) - 2,
            int(np.ceil(m_to_tb(target + INTERNAL_ANCHOR_SEARCH_DA, a, b))),
        )
        windows[anchor["name"]] = (tlo, thi)
        if thi <= tlo and anchor["status"] == "accepted":
            anchor["status"] = "implausible"
            anchor["reason"] = "raw persistence window falls outside the spectrum"

    usable_totals = {anchor["name"]: 0 for anchor in anchors}
    for block_index, (start, end) in enumerate(zip(edges[:-1], edges[1:]), 1):
        if should_stop is not None and should_stop():
            raise AnalysisCancelled("the analysis was cancelled")
        for anchor in anchors:
            name = anchor["name"]
            if anchor["status"] != "accepted":
                statuses[name].append("missing")
                block_centres[name].append(None)
                continue
            tlo, thi = windows[name]
            if thi <= tlo:
                statuses[name].append("missing")
                block_centres[name].append(None)
                continue
            try:
                window = np.asarray(dataset[start:end, tlo : thi + 1], dtype=np.float64)
            except (OSError, TypeError, ValueError):
                reason = f"{name} anchor persistence window is unavailable or malformed"
                break
            finite = np.isfinite(window)
            usable_rows = (finite & (window > 0)).any(axis=1)
            usable_totals[name] += int(usable_rows.sum())
            if not usable_rows.any():
                statuses[name].append("missing")
                block_centres[name].append(None)
                continue

            usable = window[usable_rows]
            finite_usable = np.isfinite(usable)
            counts = finite_usable.sum(axis=0)
            totals = np.where(finite_usable, usable, 0.0).sum(axis=0)
            block_mean = np.divide(
                totals,
                counts,
                out=np.zeros_like(totals),
                where=counts > 0,
            )
            block = np.zeros(int(dataset.shape[1]), dtype=np.float64)
            block[tlo : thi + 1] = block_mean
            candidate = _detect_internal_anchor(
                avg=block,
                a=a,
                b=b,
                name=name,
                target_mz=float(anchor["target_mz"]),
            )
            statuses[name].append(candidate["status"])
            block_centres[name].append(
                candidate.get("observed_file_mz")
                if candidate["status"] == "accepted"
                else None
            )
        if reason is not None:
            break
        if progress is not None:
            progress(0.1 + 0.9 * block_index / block_count)

    if reason is not None:
        for anchor in anchors:
            anchor["persistence"] = {
                "available": False,
                "reason": reason,
            }
        return
    for anchor in anchors:
        name = anchor["name"]
        if usable_totals[name] < 2:
            anchor["persistence"] = {
                "available": False,
                "reason": (
                    f"{name} anchor persistence window has fewer than two usable "
                    f"raw cycles (found {usable_totals[name]}); requires at least 2"
                ),
            }
            continue
        block_statuses = statuses[name]
        accepted = block_statuses.count("accepted")
        fraction = accepted / block_count
        anchor["persistence"] = {
            "available": True,
            "blocks": block_count,
            "accepted_blocks": accepted,
            "fraction": float(fraction),
            "statuses": block_statuses,
            "block_centres_file_mz": block_centres[name],
        }
        if (
            anchor["status"] == "accepted"
            and fraction < INTERNAL_ANCHOR_MIN_PERSISTENCE
        ):
            anchor["status"] = "ambiguous"
            anchor["reason"] = (
                f"raw-spectrum persistence is {fraction:.0%} "
                f"({accepted}/{block_count} blocks); requires at least "
                f"{INTERNAL_ANCHOR_MIN_PERSISTENCE:.0%}"
            )


def _detect_internal_anchor(avg, a, b, name, target_mz):
    tlo = max(1, int(np.floor(m_to_tb(target_mz - INTERNAL_ANCHOR_SEARCH_DA, a, b))))
    thi = min(
        avg.size - 2,
        int(np.ceil(m_to_tb(target_mz + INTERNAL_ANCHOR_SEARCH_DA, a, b))),
    )
    result = {
        "name": name,
        "target_mz": float(target_mz),
        "status": "missing",
        "reason": "search window falls outside the recorded spectrum",
        "observed_file_mz": None,
        "corrected_mz": None,
        "timebin": None,
        "prominence": None,
        "snr": None,
    }
    if thi <= tlo:
        return result

    section = avg[tlo : thi + 1]
    maxima = (
        np.where((section[1:-1] > section[:-2]) & (section[1:-1] >= section[2:]))[0]
        + tlo
        + 1
    )
    if maxima.size == 0:
        result["reason"] = (
            f"no local maximum within +-{INTERNAL_ANCHOR_SEARCH_DA:.2f} Da"
        )
        return result

    low = section[section <= np.percentile(section, 70.0)]
    baseline = float(np.median(low)) if low.size else 0.0
    mad = float(np.median(np.abs(low - baseline))) if low.size else 0.0
    sigma = max(1.4826 * mad, np.sqrt(max(baseline, 1.0)), 1e-9)
    fwhm_tb = max(2.0, a * np.sqrt(target_mz) / (2.0 * 2400.0))

    # Multiple sampled maxima inside one physical linewidth are one peak, not
    # competing anchors. Keep the tallest representative before ambiguity checks.
    ordered = sorted(maxima.tolist(), key=lambda idx: float(avg[idx]), reverse=True)
    representatives = []
    for idx in ordered:
        if all(abs(idx - kept) >= fwhm_tb for kept in representatives):
            representatives.append(idx)

    candidates = []
    for idx in representatives:
        radius = max(3, int(np.ceil(2.5 * fwhm_tb)))
        left = avg[max(tlo, idx - radius) : idx + 1]
        right = avg[idx : min(thi, idx + radius) + 1]
        local_floor = max(float(left.min()), float(right.min()))
        prominence = max(0.0, float(avg[idx]) - local_floor)
        snr = prominence / sigma
        candidates.append((idx, prominence, snr))
    candidates.sort(
        key=lambda item: (item[1], -abs(item[0] - m_to_tb(target_mz, a, b))),
        reverse=True,
    )
    idx, prominence, snr = candidates[0]
    result.update({"prominence": float(prominence), "snr": float(snr)})
    if prominence < INTERNAL_ANCHOR_MIN_PROMINENCE or snr < INTERNAL_ANCHOR_MIN_SNR:
        result.update(
            {
                "status": "weak",
                "reason": (
                    f"best local maximum has prominence {prominence:.2f} cps and "
                    f"S/N {snr:.2f}; requires at least "
                    f"{INTERNAL_ANCHOR_MIN_PROMINENCE:.1f} cps and "
                    f"S/N {INTERNAL_ANCHOR_MIN_SNR:.1f}"
                ),
            }
        )
        return result
    strong = [
        item
        for item in candidates[1:]
        if item[1] >= INTERNAL_ANCHOR_MIN_PROMINENCE
        and item[2] >= INTERNAL_ANCHOR_MIN_SNR
        and item[1] >= INTERNAL_ANCHOR_AMBIGUITY_RATIO * prominence
    ]
    if strong:
        result.update(
            {
                "status": "ambiguous",
                "reason": (
                    "multiple resolved maxima pass the anchor checks; the second "
                    f"has {strong[0][1] / prominence:.0%} of the leading prominence"
                ),
            }
        )
        return result

    centre_tb = _sub_bin_centre(avg, idx)
    observed = float(tb_to_m(centre_tb, a, b))
    if (
        not np.isfinite(observed)
        or abs(observed - target_mz) > INTERNAL_ANCHOR_SEARCH_DA
    ):
        result.update(
            {
                "status": "implausible",
                "reason": "sub-bin centre lies outside the allowed search window",
            }
        )
        return result
    if (
        min(
            observed - (target_mz - INTERNAL_ANCHOR_SEARCH_DA),
            (target_mz + INTERNAL_ANCHOR_SEARCH_DA) - observed,
        )
        < INTERNAL_ANCHOR_PROXIMITY_DA
    ):
        result.update(
            {
                "status": "ambiguous",
                "reason": "candidate is too close to the search-window boundary",
            }
        )
        return result
    result.update(
        {
            "status": "accepted",
            "reason": "",
            "observed_file_mz": observed,
            "timebin": float(centre_tb),
        }
    )
    return result


def _sub_bin_centre(spectrum, index):
    left, centre, right = (float(x) for x in spectrum[index - 1 : index + 2])
    curvature = left - 2.0 * centre + right
    if curvature >= 0 or not np.isfinite(curvature):
        return float(index)
    delta = 0.5 * (left - right) / curvature
    return float(index + np.clip(delta, -0.5, 0.5))


def m_to_tb(m, a, b, mass_axis=None):
    if mass_axis is not None:
        validate_mass_axis(mass_axis)
        return mass_axis.m_to_tb(m)
    return a * np.sqrt(m) + b


def tb_to_m(tb, a, b, mass_axis=None):
    if mass_axis is not None:
        validate_mass_axis(mass_axis)
        return mass_axis.tb_to_m(tb)
    return ((tb - b) / a) ** 2


def has_transmission(f):
    """True if the file carries a real transmission curve (some raw exports omit it)."""
    if "PTR-Transmission/Masses_Factors" not in f:
        return False
    try:
        mf = f["PTR-Transmission/Masses_Factors"][:]
        return bool((mf[0, 0, :] > 0).any())
    except (IndexError, KeyError, OSError, TypeError, ValueError):
        return False


def load_transmission(f):
    """Transmission curve (sorted m/z -> relative transmission).

    Raw-acquisition exports may omit PTR-Transmission entirely; fall back to flat
    unity transmission so Corrected == Raw (surfaced via has_transmission / the
    ``transmission_available`` flag in analysis output) instead of crashing."""
    if "PTR-Transmission/Masses_Factors" in f:
        mf = f["PTR-Transmission/Masses_Factors"][:]
        tm, tf = mf[0, 0, :], mf[0, 1, :]
        keep = tm > 0
        tm, tf = tm[keep], tf[keep]
        if tm.size:
            o = np.argsort(tm)
            return tm[o], tf[o]
    # no transmission table in file: unit transmission across the full m/z range
    return np.array([1.0, 1000.0]), np.array([1.0, 1.0])


def derive_sensitivity_percycle(f, min_corrected=1000.0):
    """Per-cycle ppb-per-corrected-cps from the file's own pre-computed traces.

    The instrument's concentration model folds primary-ion normalisation into the
    Conc/Corrected ratio, so this ratio is constant across masses within a cycle
    but drifts over the run. Returning it per cycle tracks that drift. Returns a
    1-D array of length n_cycles, or None if the file has no pre-computed traces.
    """
    if "TRACEdata/TraceConcentration" not in f:
        return None
    cor = np.asarray(f["TRACEdata/TraceCorrected"][:], dtype=np.float64)
    con = np.asarray(f["TRACEdata/TraceConcentration"][:], dtype=np.float64)
    ncyc = cor.shape[0]
    s = np.full(ncyc, np.nan)
    for i in range(ncyc):
        c = cor[i]
        k = con[i]
        m = (c > min_corrected) & np.isfinite(c) & np.isfinite(k) & (k > 0)
        if m.any():
            s[i] = np.median(k[m] / c[m])
    # fill any gaps with the global median
    good = np.isfinite(s)
    if not good.any():
        return None
    s[~good] = np.median(s[good])
    return s


def extract_primary(f, primary_mz=21.022, R=1200.0, block=400, mass_axis=None):
    """Per-cycle primary-ion (reagent-ion) signal used to normalise concentration.

    In H3O+ mode the primary ion is monitored via its H3(18O)+ isotope at m/z ~21
    (m/z 19 saturates the detector). Uses the pre-computed TraceRaw column nearest
    primary_mz when available (fast, and what the instrument/Viewer normalise to);
    otherwise integrates the peak from the raw spectra. Returns a length-n_cycles
    array, or None if it cannot be obtained."""
    if mass_axis is not None:
        validate_mass_axis(mass_axis)
    if "TRACEdata/TraceRaw" in f and "TRACEdata/TraceInfo" in f:
        ti = f["TRACEdata/TraceInfo"][:]
        centers = np.array([float(ti[2, c]) for c in range(ti.shape[1])])
        if mass_axis is None:
            mass_axis = load_mass_axis(f)
        centers = mass_axis.file_to_corrected(centers)
        j = int(np.argmin(np.abs(centers - primary_mz)))
        if abs(centers[j] - primary_mz) < 0.1:
            return np.asarray(f["TRACEdata/TraceRaw"][:, j], dtype=np.float64)
    try:
        (traces, _) = extract_traces(
            f, [primary_mz], R=R, block=block, mass_axis=mass_axis
        )
        return traces[primary_mz][0]
    except (IndexError, KeyError, OSError, TypeError, ValueError):
        return None


def water_cluster_ratio(
    f, cluster_mz=37.033, primary_mz=21.022, R=1200.0, mass_axis=None
):
    """Per-cycle humidity proxy X(t) = I(first water cluster) / I(primary isotope).

    The standard PTR-MS humidity measure is I(H3O+.H2O)/I(H3O+) = m/z 37 / m/z 19,
    but m/z 19 is usually saturated/blanked, so this uses the configured primary m/z
    (m/z 21.022 by default) as the denominator instead. A constant isotope factor
    cancels once the ratio is normalised to a reference. Returns a length-n_cycles
    array, or None."""
    if mass_axis is not None:
        validate_mass_axis(mass_axis)

    def get(mz):
        if "TRACEdata/TraceRaw" in f and "TRACEdata/TraceInfo" in f:
            ti = f["TRACEdata/TraceInfo"][:]
            centers = np.array([float(ti[2, c]) for c in range(ti.shape[1])])
            axis = mass_axis or load_mass_axis(f)
            centers = axis.file_to_corrected(centers)
            j = int(np.argmin(np.abs(centers - mz)))
            if abs(centers[j] - mz) < 0.1:
                return np.asarray(f["TRACEdata/TraceRaw"][:, j], dtype=np.float64)
        try:
            return extract_traces(f, [mz], R=R, mass_axis=mass_axis)[0][mz][0]
        except (IndexError, KeyError, OSError, TypeError, ValueError):
            return None

    c = get(cluster_mz)
    p = get(primary_mz)
    if c is None or p is None:
        return None
    x = np.full_like(p, np.nan)
    ok = p > 0
    x[ok] = c[ok] / p[ok]
    return x


def humidity_factor(ratio, ref, p):
    """Per-cycle correction for near-thermoneutral (humid-flagged) compounds.

    In the equilibrium limit the proton transfer is reversible and sensitivity
    scales as 1/[H2O], so signal must be multiplied by (X/X_ref)**p to normalise to
    a reference humidity X_ref. p in [0,1]: p=0 no correction (kinetic limit),
    p=1 full equilibrium (upper bound). Calibrate p from a standard at >=2
    humidities; without that it is approximate and only makes RELATIVE comparisons
    at differing humidity valid, not absolute values."""
    f = np.ones_like(ratio)
    good = np.isfinite(ratio) & (ratio > 0) & (ref > 0)
    f[good] = (ratio[good] / ref) ** p
    return f


def derive_K(f, primary, min_corrected=1000.0):
    """Concentration calibration constant K for Conc = Corrected * K / primary.

    Derived from the file's own pre-computed concentration so the default output
    reproduces the instrument's concentration. K = median over cycles of
    sensitivity(t) * primary(t), where sensitivity(t) = TraceConc/TraceCorrected.
    Returns None if the file has no pre-computed concentration data.

    NOTE: this is the *acquisition* calibration. A specific PTR-MS Viewer project
    may use a different absolute K (its own sensitivity setting); use `calibrate`
    against a reference, or pass K explicitly, to match that exactly."""
    s = derive_sensitivity_percycle(f, min_corrected=min_corrected)
    if s is None or primary is None:
        return None
    prod = s * primary
    good = np.isfinite(prod) & (primary > 0)
    if not good.any():
        return None
    return float(np.median(prod[good]))


def derive_molar_volume_info(f):
    """Return ``(Vm, source)`` for the file's drift-temperature calibration.

    A missing drift-temperature trace is deliberately distinguishable from a real
    25 °C measurement: both values are plausible, but only the former is a fallback.
    """
    try:
        data = f["AddTraces/PTR-Reaction/Data"]
        info = f["AddTraces/PTR-Reaction/Info"][0]
        names = [x.decode("latin-1").strip() for x in info]
        ti = names.index("T-Drift_Act")
        T_C = float(np.nanmean(data[:, ti]))
        if not np.isfinite(T_C):
            raise ValueError("drift temperature is not finite")
        return (22.414 * (T_C + 273.15) / 273.15, "file drift temperature")
    except (IndexError, KeyError, OSError, TypeError, ValueError):
        return (24.465, "25 °C fallback (drift metadata unavailable)")


def derive_molar_volume(f):
    """Vm [L/mol] at drift temperature, or the documented 25 °C fallback."""
    return derive_molar_volume_info(f)[0]


# ---------- peak extraction ----------
def find_apex(avgspec, a, b, target_m, tol=0.15, mass_axis=None):
    tlo = max(0, int(m_to_tb(target_m - tol, a, b, mass_axis)))
    thi = min(len(avgspec), int(m_to_tb(target_m + tol, a, b, mass_axis)))
    if thi <= tlo:
        return target_m, tlo, thi
    apex_tb = tlo + int(np.argmax(avgspec[tlo:thi]))
    return tb_to_m(apex_tb, a, b, mass_axis), tlo, thi


def refine_apex_local(avgspec, a, b, apex0, tol=0.035, mass_axis=None):
    """Re-centre a peak on THIS spectrum's real maximum near a known apex.

    Peak positions drift between time intervals (mass-cal drift; a compound may be
    absent in a background). Given a canonical apex, return the local maximum
    within +-tol only if it is a genuine interior peak that clears the noise —
    otherwise None, meaning keep the canonical position (so a background where the
    compound is absent does NOT chase an unrelated neighbour). Mirrors the viz
    per-interval overlay exactly so the delivered CSV matches what was reviewed."""
    lo = max(0, int(np.floor(m_to_tb(apex0 - tol, a, b, mass_axis))))
    hi = min(
        len(avgspec) - 1,
        int(np.ceil(m_to_tb(apex0 + tol, a, b, mass_axis))),
    )
    if hi - lo < 2:
        return None
    bi = lo + int(np.argmax(avgspec[lo : hi + 1]))
    bv = float(avgspec[bi])
    if bi <= lo or bi >= hi:  # max at an edge -> monotonic climb, no clear peak
        return None
    floor = max(float(avgspec[lo]), float(avgspec[hi]))
    if bv >= 3 and bv >= 1.25 * floor:
        return tb_to_m(bi, a, b, mass_axis)
    return None


def peak_window(apex_m, a, b, R, mass_axis=None):
    hw = apex_m / (2 * R)
    wl = int(np.floor(m_to_tb(apex_m - hw, a, b, mass_axis)))
    wr = int(np.ceil(m_to_tb(apex_m + hw, a, b, mass_axis)))
    return wl, wr


def peak_window_lr(apex_m, a, b, hwL, hwR, mass_axis=None):
    """Integration window from explicit left/right half-widths in m/z (per-peak,
    possibly asymmetric — the window need not be centred on the apex)."""
    wl = int(np.floor(m_to_tb(apex_m - hwL, a, b, mass_axis)))
    wr = int(np.ceil(m_to_tb(apex_m + hwR, a, b, mass_axis)))
    return wl, wr


def _hw_for(m, apex_m, R, windows):
    """(left, right) half-widths in m/z for a peak: an explicit per-peak override
    (a scalar for symmetric, or a (left,right) pair), else R-derived symmetric."""
    if windows and m in windows and windows[m]:
        w = windows[m]
        if isinstance(w, (tuple, list)):
            return float(w[0]), float(w[1])
        return float(w), float(w)
    hw = apex_m / (2 * R)
    return hw, hw


def _cluster(masses, gap=0.20):
    """Group masses whose neighbours are closer than `gap` (needs deconvolution)."""
    ms = sorted(masses)
    groups, cur = [], [ms[0]]
    for m in ms[1:]:
        if m - cur[-1] < gap:
            cur.append(m)
        else:
            groups.append(cur)
            cur = [m]
    groups.append(cur)
    return groups


def _sigma_tb(mu_m, a, R_phys, mass_axis=None):
    """Gaussian sigma in timebins for a peak at mu_m given physical resolution."""
    sigma_m = mu_m / (2.3548 * R_phys)
    if mass_axis is None:
        dtb_dm = a / (2 * np.sqrt(mu_m))
    else:
        file_mass = mass_axis.corrected_to_file(mu_m)
        dtb_dm = a / (2 * mass_axis.scale * np.sqrt(file_mass))
    return sigma_m * dtb_dm


def _cluster_design(
    centers_m,
    a,
    b,
    R=1200.0,
    R_phys=2400.0,
    windows=None,
    nbin=None,
    mass_axis=None,
):
    """Precompute the Gaussian-unmixing design for one cluster of overlapping peaks.

    Returns (tlo, thi, P, norm): the timebin span to read, the projection matrix P
    (cluster amplitudes A = Y[:, tlo:thi] @ P), and per-peak normalisation to the
    window-sum scale so deconvolved peaks share the isolated-peak Raw scale. Split
    out from deconvolve_cluster so many clusters can be applied in a single shared
    streaming pass (a gzip-compressed file decompresses whole rows, so a per-cluster
    pass would re-decompress the entire dataset once per cluster)."""
    centers_tb = np.array([m_to_tb(m, a, b, mass_axis) for m in centers_m])
    sig_tb = np.array([_sigma_tb(m, a, R_phys, mass_axis=mass_axis) for m in centers_m])
    tlo = int(np.floor(centers_tb.min() - 6 * sig_tb.max()))
    thi = int(np.ceil(centers_tb.max() + 6 * sig_tb.max()))
    tlo = max(0, tlo)
    if nbin is not None:
        thi = min(nbin, thi)
    x = np.arange(tlo, thi)
    # design matrix G (n_tb x K), unit-height Gaussians
    G = np.exp(-0.5 * ((x[:, None] - centers_tb[None, :]) / sig_tb[None, :]) ** 2)
    P = G @ np.linalg.inv(G.T @ G)  # n_tb x K : A = Y @ P
    # normalisation: window-sum of each peak's own unit Gaussian (matches isolated)
    norm = np.zeros(len(centers_m))
    for k, m in enumerate(centers_m):
        hwL, hwR = _hw_for(m, m, R, windows)
        wl, wr = peak_window_lr(m, a, b, hwL, hwR, mass_axis)
        xx = np.arange(wl, wr)
        norm[k] = np.exp(-0.5 * ((xx - centers_tb[k]) / sig_tb[k]) ** 2).sum()
    return tlo, thi, P, norm


def deconvolve_cluster(
    f,
    centers_m,
    a,
    b,
    R=1200.0,
    R_phys=2400.0,
    block=400,
    windows=None,
    mass_axis=None,
):
    """Separate overlapping peaks by vectorised linear least-squares Gaussian
    unmixing. Returns dict center_m -> raw_trace (scaled to match the window-sum
    definition so isolated and deconvolved peaks share one Raw scale)."""
    if mass_axis is None:
        mass_axis = load_mass_axis(f)
    else:
        validate_mass_axis(mass_axis)
    a, b = mass_axis.a, mass_axis.b
    inten = f["SPECdata/Intensities"]
    ncyc = inten.shape[0]
    tlo, thi, P, norm = _cluster_design(
        centers_m,
        a,
        b,
        R,
        R_phys,
        windows,
        nbin=inten.shape[1],
        mass_axis=mass_axis,
    )
    traces = np.empty((ncyc, len(centers_m)))
    for i in range(0, ncyc, block):
        j = min(i + block, ncyc)
        Y = np.asarray(inten[i:j, tlo:thi], dtype=np.float64)
        Y[~np.isfinite(Y)] = 0.0
        A = Y @ P  # (j-i) x K amplitudes
        np.clip(A, 0, None, out=A)
        traces[i:j, :] = A * norm[None, :]
    return {m: traces[:, k] for k, m in enumerate(centers_m)}


def _profile_report(profile):
    if profile is None:
        return {"usable": False, "reason": "legacy model selected"}
    return {
        "usable": bool(profile["usable"]),
        "reason": profile["reason"],
        "n_reference_peaks": int(profile["n_reference_peaks"]),
    }


def _fit_report(masses, fitted, n_reference_peaks):
    def _finite(value):
        return round(float(value), 6) if np.isfinite(value) else None

    return {
        "masses": [float(mass) for mass in masses],
        "method": "empirical-v1",
        "status": fitted["status"],
        "reason": fitted["reason"],
        "n_reference_peaks": int(n_reference_peaks),
        "shift_timebins": round(float(fitted["shift_tb"]), 4),
        "width_scale": round(float(fitted["width_scale"]), 4),
        "rank": int(fitted["rank"]),
        "condition": _finite(fitted["condition"]),
        "component_correlation": _finite(fitted["correlation"]),
        "relative_residual": _finite(fitted["relative_residual"]),
    }


def extract_traces(
    f,
    target_masses,
    R=1200.0,
    R_phys=2400.0,
    block=400,
    refine_tol=0.02,
    cluster_gap=0.20,
    windows=None,
    per_range=None,
    range_refine_tol=0.035,
    *,
    mass_axis=None,
    progress=None,
    should_stop=None,
    peak_fit_model="gaussian-v1",
    fit_diagnostics=None,
):
    """Return dict m -> (raw_trace[ncycles], apex_m). One streaming pass.

    Peak centring is two-stage: (1) the run's accepted two-anchor affine mass-axis
    correction aligns every timebin, independently of the selected panel; (2) a tight
    local apex search (+-refine_tol) around the corrected position snaps isolated
    peaks to their exact centres. This lets closely-spaced peaks be resolved without
    applying the former target-panel-derived multiplicative drift a second time.

    windows: optional {target_mass: half_width_m} to override the R-derived
    integration window for specific peaks (from an expert's viz adjustment).

    per_range: optional {label: (lo, hi)} (1-based inclusive cycles). When given,
    each interval's cycles are re-integrated with each isolated peak's apex/window
    RE-CENTRED on that interval's own average spectrum (peaks drift between
    intervals). apex_m in the return stays the whole-run value (transmission moves
    <0.1% over the drift); only the per-cycle window changes. Clustered peaks keep
    their whole-run deconvolved trace.

    progress: optional callback with cycles consumed / cycles this pass reads, after
    each block of either pass. Reads are the honest axis: on the 2 GB / 20,725-cycle
    run the streaming pass takes 14.6 s and the interval re-centring 13.9 s of the
    28.6 s they share, so a bar driven by the first pass alone would sit at 100 %
    through the second. It is 89 % of the ~33 s an open costs in total.
    should_stop: optional callback polled once per block in both passes; when it
    returns true the pass raises AnalysisCancelled. ``peak_fit_model`` may be the
    historical ``gaussian-v1`` or measured-shape ``empirical-v1``. Diagnostics are
    written into the optional mutable ``fit_diagnostics`` mapping. All additions
    default to the historical arithmetic.
    """
    if mass_axis is None:
        mass_axis = load_mass_axis(f)
    else:
        validate_mass_axis(mass_axis)
    a, b = mass_axis.a, mass_axis.b
    inten = f["SPECdata/Intensities"]
    ncyc = inten.shape[0]
    avg = np.asarray(f["SPECdata/AverageSpec"][:], dtype=np.float64)
    avg = np.where(np.isfinite(avg), avg, 0.0)  # tolerate rare corrupt bins
    nbin = avg.shape[0]

    # The run-wide correction is already encoded in ``mass_axis``. Applying the
    # old target-panel-derived multiplicative drift here as well would scale the
    # axis twice and make the result depend on which compounds were selected.

    # isolated peaks -> window-sum; clustered peaks -> Gaussian deconvolution
    groups = _cluster(target_masses, gap=cluster_gap)
    isolated = [g[0] for g in groups if len(g) == 1]
    clusters = [g for g in groups if len(g) > 1]

    # Isolated peaks retain a tight local snap after the shared axis correction.
    # Clustered model centres always stay on the common mass axis.
    apexes = {}
    search_tol = refine_tol
    for m in isolated:
        apex_m, _, _ = find_apex(avg, a, b, m, tol=search_tol, mass_axis=mass_axis)
        apexes[m] = apex_m
    for g in clusters:
        for m in g:
            apexes[m] = m

    # per-range average-spectrum accumulators, filled during the isolated pass so
    # interval re-centring costs no extra read of the whole run
    want_ranges = {
        lbl: (int(lo), int(hi))
        for lbl, (lo, hi) in (per_range or {}).items()
        if int(hi) >= int(lo)
    }
    rsum = {lbl: np.zeros(nbin, dtype=np.float64) for lbl in want_ranges}
    rcnt = {lbl: 0 for lbl in want_ranges}

    # ONE shared streaming pass computes every isolated window-sum, every cluster's
    # deconvolved amplitudes, and the per-range average-spectrum accumulators. A
    # gzip-compressed file decompresses whole rows, so reading the file once here —
    # instead of once per cluster (deconvolve_cluster) — turns an O(n_clusters) set
    # of full-file decompression passes into a single pass (the dominant cost on
    # dense breath spectra: ~15 clusters was ~15x slower).
    win_tb = {
        m: peak_window_lr(
            apexes[m],
            a,
            b,
            *_hw_for(m, apexes[m], R, windows),
            mass_axis=mass_axis,
        )
        for m in isolated
    }
    iso_buf = {m: np.empty(ncyc) for m in isolated}
    empirical_profile = None
    if peak_fit_model == "empirical-v1":
        isolated_tb = [m_to_tb(apexes[m], a, b, mass_axis) for m in isolated]
        isolated_sigma = [
            _sigma_tb(apexes[m], a, R_phys, mass_axis=mass_axis) for m in isolated
        ]
        empirical_profile = peak_fit.estimate_empirical_profile(
            avg, isolated_tb, isolated_sigma
        )
    elif peak_fit_model != "gaussian-v1":
        raise ValueError(f"unknown peak-fit model: {peak_fit_model}")

    cluster_design = []
    cluster_reports = []
    for g in clusters:
        caps = [apexes[m] for m in g]
        design = None
        if empirical_profile is not None:
            centers_tb = np.array([m_to_tb(m, a, b, mass_axis) for m in caps])
            sigmas_tb = np.array(
                [_sigma_tb(m, a, R_phys, mass_axis=mass_axis) for m in caps]
            )
            fitted = peak_fit.fit_group_design(
                avg, centers_tb, sigmas_tb, empirical_profile
            )
            if fitted["usable"]:
                shifted_tb = centers_tb + fitted["shift_tb"]
                shifted_m = [tb_to_m(value, a, b, mass_axis) for value in shifted_tb]
                for m, measured in zip(g, shifted_m):
                    apexes[m] = measured
                norm = np.zeros(len(g))
                for k, (m, measured) in enumerate(zip(g, shifted_m)):
                    hw_l, hw_r = _hw_for(m, measured, R, windows)
                    wl, wr = peak_window_lr(
                        measured, a, b, hw_l, hw_r, mass_axis
                    )
                    left = max(0, wl - fitted["lo"])
                    right = min(fitted["components"].shape[0], wr - fitted["lo"])
                    norm[k] = fitted["components"][left:right, k].sum()
                design = {
                    "method": "empirical-v1",
                    "fit": fitted,
                    "norm": norm,
                }
                report = _fit_report(
                    g, fitted, empirical_profile["n_reference_peaks"]
                )
                if not empirical_profile["usable"]:
                    report["method"] = "gaussian-fallback-v1"
                    report["fallback_reason"] = empirical_profile["reason"]
                cluster_reports.append(report)
            else:
                design = {"method": "unresolved"}
                cluster_reports.append(
                    {
                        "masses": [float(m) for m in g],
                        "method": "gaussian-fallback-v1",
                        "status": "unresolved",
                        "reason": fitted["reason"],
                        "fallback_reason": empirical_profile["reason"],
                    }
                )
        if design is None:
            apex_hw = {ap: _hw_for(m, ap, R, windows) for m, ap in zip(g, caps)}
            gaussian = _cluster_design(
                caps,
                a,
                b,
                R,
                R_phys,
                apex_hw,
                nbin=nbin,
                mass_axis=mass_axis,
            )
            design = {"method": "gaussian-v1", "design": gaussian}
            cluster_reports.append(
                {
                    "masses": [float(m) for m in g],
                    "method": "gaussian-v1",
                    "status": "legacy",
                    "reason": "legacy model selected",
                }
            )
        cluster_design.append(design)
    cluster_buf = [np.empty((ncyc, len(g))) for g in clusters]
    if fit_diagnostics is not None:
        fit_diagnostics.clear()
        fit_diagnostics.update(
            {
                "model": peak_fit_model,
                "profile": _profile_report(empirical_profile),
                "clusters": cluster_reports,
            }
        )

    # Cycles the two passes read between them: every cycle once here, plus every
    # cycle the interval re-centring below reads again. On the 2 GB fixture with its
    # 23 intervals that is 20,725 + 19,523 = 40,248, and the two passes take 14.6 s
    # and 13.9 s — so reporting only this pass would leave the bar sitting at 100 %
    # for the whole of the second one. The pass is one axis, not two. With no
    # intervals to re-centre the span is just ncyc, which is the plain axis.
    span = ncyc + sum(hi - lo + 1 for lo, hi in want_ranges.values())

    for i in range(0, ncyc, block):
        if should_stop is not None and should_stop():
            raise AnalysisCancelled("the analysis was cancelled")
        j = min(i + block, ncyc)
        chunk = np.asarray(inten[i:j, :], dtype=np.float64)
        chunk[~np.isfinite(chunk)] = 0.0
        for m in isolated:
            wl, wr = win_tb[m]
            iso_buf[m][i:j] = chunk[:, wl:wr].sum(axis=1)
        for ci, design in enumerate(cluster_design):
            if design["method"] == "empirical-v1":
                fitted = design["fit"]
                if fitted["status"] == "reliable":
                    amplitudes = peak_fit.apply_group_design(chunk, fitted)
                    cluster_buf[ci][i:j, :] = amplitudes * design["norm"][None, :]
                else:
                    cluster_buf[ci][i:j, :] = np.nan
            elif design["method"] == "unresolved":
                cluster_buf[ci][i:j, :] = np.nan
            else:
                tlo, thi, projection, norm = design["design"]
                amplitudes = chunk[:, tlo:thi] @ projection
                np.clip(amplitudes, 0, None, out=amplitudes)
                cluster_buf[ci][i:j, :] = amplitudes * norm[None, :]
        for lbl, (lo, hi) in want_ranges.items():  # cycles are 1-based inclusive
            c0, c1 = max(i, lo - 1), min(j, hi)
            if c1 > c0:
                rsum[lbl] += chunk[c0 - i : c1 - i, :].sum(axis=0)
                rcnt[lbl] += c1 - c0
        if progress is not None:
            progress(min(j, span) / span)

    traces = {}
    for m in isolated:
        traces[m] = iso_buf[m]
    for g, buf in zip(clusters, cluster_buf):
        for k, m in enumerate(g):
            traces[m] = buf[:, k]

    # second pass over ONLY each interval's cycles: re-centre each isolated peak on
    # that interval's average spectrum and overwrite those cycles (clustered peaks
    # keep their whole-run deconvolved trace)
    read = ncyc  # cycles consumed over both passes: the one progress axis
    for lbl, (lo, hi) in want_ranges.items():
        if not rcnt[lbl]:
            continue
        avg_r = rsum[lbl] / rcnt[lbl]
        rwin = {}
        for m in isolated:
            if windows and m in windows:
                continue  # hand-placed window: keep it everywhere (matches viz winManual)
            ap = refine_apex_local(
                avg_r,
                a,
                b,
                apexes[m],
                tol=range_refine_tol,
                mass_axis=mass_axis,
            )
            if ap is None:
                continue  # no clear interval peak -> keep whole-run window
            rwin[m] = peak_window_lr(
                ap,
                a,
                b,
                *_hw_for(m, ap, R, windows),
                mass_axis=mass_axis,
            )
        range_cluster_design = []
        if empirical_profile is not None:
            for cluster_index, group in enumerate(clusters):
                centres_tb = np.array(
                    [m_to_tb(apexes[m], a, b, mass_axis) for m in group]
                )
                sigmas_tb = np.array(
                    [
                        _sigma_tb(apexes[m], a, R_phys, mass_axis=mass_axis)
                        for m in group
                    ]
                )
                fitted = peak_fit.fit_group_design(
                    avg_r, centres_tb, sigmas_tb, empirical_profile
                )
                if not fitted["usable"]:
                    continue
                norm = np.zeros(len(group))
                shifted_tb = centres_tb + fitted["shift_tb"]
                for k, (m, centre_tb) in enumerate(zip(group, shifted_tb)):
                    measured = tb_to_m(centre_tb, a, b, mass_axis)
                    hw_l, hw_r = _hw_for(m, measured, R, windows)
                    window_lo, window_hi = peak_window_lr(
                        measured, a, b, hw_l, hw_r, mass_axis
                    )
                    left = max(0, window_lo - fitted["lo"])
                    right = min(
                        fitted["components"].shape[0],
                        window_hi - fitted["lo"],
                    )
                    norm[k] = fitted["components"][left:right, k].sum()
                range_cluster_design.append(
                    (cluster_index, {"fit": fitted, "norm": norm})
                )
                range_report = _fit_report(
                    group,
                    fitted,
                    empirical_profile["n_reference_peaks"],
                )
                if not empirical_profile["usable"]:
                    range_report["method"] = "gaussian-fallback-v1"
                    range_report["fallback_reason"] = empirical_profile["reason"]
                cluster_reports[cluster_index].setdefault("ranges", {})[lbl] = (
                    range_report
                )
        if not rwin and not range_cluster_design:
            continue
        for i in range(lo - 1, hi, block):
            if should_stop is not None and should_stop():
                raise AnalysisCancelled("the analysis was cancelled")
            j = min(i + block, hi)
            chunk = np.asarray(inten[i:j, :], dtype=np.float64)
            chunk[~np.isfinite(chunk)] = 0.0
            for m, (wl, wr) in rwin.items():
                traces[m][i:j] = chunk[:, wl:wr].sum(axis=1)
            for cluster_index, design in range_cluster_design:
                fitted = design["fit"]
                group = clusters[cluster_index]
                if fitted["status"] == "reliable":
                    amplitudes = peak_fit.apply_group_design(chunk, fitted)
                    values = amplitudes * design["norm"][None, :]
                else:
                    values = np.full((j - i, len(group)), np.nan)
                for k, mass in enumerate(group):
                    traces[mass][i:j] = values[:, k]
            if progress is not None:
                read += j - i
                progress(min(read, span) / span)

    # An interval with nothing to re-centre was never re-read, so the count can stop
    # short of the span it planned. The pass has still ended here.
    if progress is not None and read < span:
        progress(1.0)

    return {m: (traces[m], apexes[m]) for m in target_masses}, (a, b)


# ---------- automatic segmentation ----------

# Gap merging. A plateau can be split by a wobble in one sample rather than by a
# real change of phase, and a cycle count cannot tell those apart across
# instruments, so the limit is a wall-clock window and the decision is made from
# the discriminator itself (see merge_adjacent_segments).
MERGE_BAND = 2.0  # how far a gap may leave its neighbours' level and still merge
MERGE_GAP_WINDOW_S = 60.0  # merged gaps cover at most this much acquisition time
MERGE_MIN_GAP_CYCLES = 30  # ... but never fewer cycles than this
# the historical length limits, kept for callers with segments but no signal
MERGE_HIGH_GAP_DEFAULT = 60
MERGE_LOW_GAP_DEFAULT = 200
MERGE_REASON_HELD = "level held"
MERGE_REASON_FELL = "fell to baseline"
MERGE_REASON_ADJACENT = "adjacent"
MERGE_REASON_LENGTH = "length only"


def build_discriminator(f, mz_lo=40.0, mz_hi=200.0, block=400, mass_axis=None):
    """Per-cycle composite VOC signal, ~1 at background and high during samples.

    Fast path uses the pre-computed TraceRaw (normalising each strong VOC trace to
    its own baseline so no single peak dominates). Fallback streams the raw spectra
    and sums a VOC m/z band. Returns a 1-D array length n_cycles."""
    if mass_axis is not None:
        validate_mass_axis(mass_axis)
    if "TRACEdata/TraceRaw" in f and "TRACEdata/TraceInfo" in f:
        ti = f["TRACEdata/TraceInfo"][:]
        centers = np.array([float(ti[2, c]) for c in range(ti.shape[1])])
        if mass_axis is None:
            mass_axis = load_mass_axis(f)
        centers = mass_axis.file_to_corrected(centers)
        band = np.where((centers >= mz_lo) & (centers <= mz_hi))[0]
        R = np.asarray(f["TRACEdata/TraceRaw"][:, band], dtype=np.float64)
        med = np.median(R, axis=0)
        pos = med[med > 0]
        if pos.size:
            strong = med > np.percentile(pos, 70)
            if strong.any():
                Rs = R[:, strong] / med[strong]
                return Rs.mean(axis=1)
    # fallback: total ion current in a VOC timebin band, streamed
    if mass_axis is None:
        mass_axis = load_mass_axis(f)
    a, b = mass_axis.a, mass_axis.b
    inten = f["SPECdata/Intensities"]
    ncyc = inten.shape[0]
    tlo = max(0, int(m_to_tb(mz_lo, a, b, mass_axis)))
    thi = min(inten.shape[1], int(m_to_tb(mz_hi, a, b, mass_axis)))
    tic = np.empty(ncyc)
    for i in range(0, ncyc, block):
        j = min(i + block, ncyc)
        chunk = np.asarray(inten[i:j, tlo:thi], dtype=np.float64)
        chunk[~np.isfinite(chunk)] = 0.0
        tic[i:j] = chunk.sum(axis=1)
    base = np.median(tic[tic > 0]) or 1.0
    return tic / base


def detect_segments(
    f,
    discriminator=None,
    min_duration=30,
    trim=8,
    grad_thr=0.02,
    smooth=9,
    high_ratio=3.0,
):
    """Detect stable measurement plateaus (candidate time ranges).

    Works in log space so the large sample/background dynamic range is handled by
    relative changes. A cycle is 'stable' where the smoothed log-signal gradient is
    small; runs of stable cycles longer than min_duration (edge-trimmed) become
    segments. Each is classified 'high' (elevated / sample) or 'low' (background or
    setup) relative to the run baseline. The caller assigns meaningful labels.

    Returns list of dicts: start_cycle/end_cycle (1-based inclusive),
    start_s/end_s, n_cycles, level (x baseline), class."""
    D = build_discriminator(f) if discriminator is None else discriminator
    ncyc = len(D)
    dur = spec_duration_s(f)
    L = np.log10(np.clip(D, 1e-3, None))
    if smooth > 1:
        L = np.convolve(L, np.ones(smooth) / smooth, mode="same")
    g = np.abs(np.gradient(L))
    stable = g < grad_thr
    baseline = discriminator_baseline(D)

    segs = []
    i = 0
    while i < ncyc:
        if not stable[i]:
            i += 1
            continue
        j = i
        while j < ncyc and stable[j]:
            j += 1
        if j - i >= min_duration:
            lo, hi = i + trim, j - trim  # 0-based, trimmed
            if hi - lo >= 15:
                level = float(D[lo:hi].mean() / baseline)
                segs.append(
                    dict(
                        start_cycle=lo + 1,
                        end_cycle=hi,
                        n_cycles=hi - lo,
                        start_s=round(lo * dur, 1),
                        end_s=round((hi - 1) * dur, 1),
                        level=round(level, 2),
                        **{"class": "high" if level >= high_ratio else "low"},
                    )
                )
        i = j
    return segs


def discriminator_baseline(D):
    """The level reference every segment level is measured against.

    The 20th percentile of the discriminator — the run's own background — or 1.0
    when that is degenerate (an empty or wholly non-finite trace). detect_segments
    and merge_adjacent_segments both take their levels from here, which is what makes
    a gap's level directly comparable to a plateau's."""
    values = np.asarray(D, dtype=np.float64)
    if values.size == 0:
        return 1.0
    value = float(np.percentile(values, 20))
    return value if np.isfinite(value) and value > 0 else 1.0


def merge_gap_cap(f, window_s=MERGE_GAP_WINDOW_S):
    """How many cycles a merged gap may cover, from the file's own cycle time.

    ``window_s`` seconds of acquisition, floored at ``MERGE_MIN_GAP_CYCLES`` so a fast
    file is never capped below a few dozen cycles. This is what makes the merge verdict
    independent of the acquisition speed: the same 40-second wobble is 40 cycles at
    1 s/cycle and 8 at 5 s/cycle, and both are one sample either way."""
    return max(MERGE_MIN_GAP_CYCLES, int(round(window_s / spec_duration_s(f))))


def _gap_record(
    previous,
    current,
    gap,
    D,
    baseline,
    band,
    evidence,
    previous_level=None,
    lower_bound=True,
):
    """Provenance for one candidate gap, or None when the gap refuses the merge.

    ``previous_level`` is the level of the plateau the gap actually abuts. A merged
    interval carries a blended level, and testing a decay against that average rather
    than against its last plateau makes the verdict depend on what was merged first.

    The gap's cycles are ``prev.end_cycle + 1 … cur.start_cycle - 1`` (1-based
    inclusive), i.e. ``D[prev.end_cycle : cur.start_cycle - 1]`` in 0-based slice
    form — the same convention detect_segments writes start_cycle/end_cycle with.
    """
    if not evidence:
        return {
            "cycles": gap,
            "min_level": None,
            "max_level": None,
            "reason": MERGE_REASON_LENGTH,
        }
    if gap <= 0:
        return {
            "cycles": gap,
            "min_level": None,
            "max_level": None,
            "reason": MERGE_REASON_ADJACENT,
        }
    lo = max(0, min(int(previous["end_cycle"]), D.size))
    hi = max(lo, min(int(current["start_cycle"]) - 1, D.size))
    window = D[lo:hi]
    if window.size == 0:
        return None  # the gap lies outside the signal: no evidence either way
    low = float(np.min(window)) / baseline
    high = float(np.max(window)) / baseline
    if not (np.isfinite(low) and np.isfinite(high)):
        return None  # cycles we cannot measure are not evidence of a wobble
    levels = [
        float(
            previous_level
            if previous_level is not None
            else previous.get("level") or 0.0
        ),
        float(current.get("level") or 0.0),
    ]
    # the gap must never have left the phase its neighbours are in: no excursion out
    # of the phase, and for a sample no collapse back to the baseline either. The
    # collapse test only means something for a sample: a background cannot fall out of
    # itself, so a dropout toward zero is still background and must not split a blank
    # into two shorter ones, which is the one thing a good blank is needed for.
    if (lower_bound and low < min(levels) / band) or high > max(levels) * band:
        return None
    # the reason says how far the gap came down, measured against the background
    # itself: at least band x baseline it held a level of its own (a wobble), below
    # that it sat at the run's background and only stayed mergeable inside the band
    reason = MERGE_REASON_HELD if low >= band else MERGE_REASON_FELL
    return {
        "cycles": gap,
        "min_level": round(low, 2),
        "max_level": round(high, 2),
        "reason": reason,
    }


def merge_adjacent_segments(
    segments,
    high_gap=MERGE_HIGH_GAP_DEFAULT,
    low_gap=MERGE_LOW_GAP_DEFAULT,
    *,
    discriminator=None,
    baseline=None,
    band=MERGE_BAND,
    cap=None,
):
    """Merge consecutive same-class plateaus separated only by a short transition.

    Plateau detection splits one physical period into several entries whenever the
    signal briefly wobbles — a sample that momentarily dips, or (very commonly) a
    long background/setup phase broken by transients into a run of small pieces plus
    slivers. Merge adjacent entries of the SAME class whose unclassified gap passes
    the test for that class; an opposite-class plateau between them is always a hard
    boundary (samples never merge across a background and vice-versa).

    Two tests are available, and every merged gap says which one ran:

    * **evidence** — given whenever the caller has the file, i.e. every real path.
      With ``gap`` the cycles strictly between the two plateaus, ``level`` the same
      x-baseline level detect_segments reports, and ``baseline`` the run's own
      background, a gap merges only when ALL of

          gap_cycles <= cap
          max(D[gap]) / baseline <= max(prev.level, cur.level) * band
          min(D[gap]) / baseline >= min(prev.level, cur.level) / band  (samples only)

      The upper test belongs to both classes: signal that left the phase is a
      boundary. The lower one is a sample test only, because a background cannot fall
      out of itself — a dropout toward zero is still the same blank, and splitting it
      would cost the longer reference interval a good blank exists to provide. A
      sample that came back down to the baseline, on the other hand, did end.

      A gap is judged against the plateau it abuts: in a chain of merges that plateau,
      not the running average of everything merged so far, or the verdict would depend
      on which plateau happened to come first. ``cap`` is the wall-clock-derived length
      limit (merge_gap_cap); ``high_gap`` / ``low_gap`` override it with a forced cycle
      count for high / low runs, and a class capped at 0 never merges. A gap holding a
      non-finite cycle is refused.
    * **length** — no ``discriminator``: the historical rule, merge when the gap
      ``<= high_gap``/``low_gap``. Kept so callers holding segments but no signal
      still get the behaviour they had.

    A class capped at 0 never merges, whatever the path.

    Returns copies. Each merged gap is appended to ``merged_gaps`` as
    ``{cycles, min_level, max_level, reason}``: ``level held`` means the gap kept a
    level of its own (at least ``band`` x baseline — a wobble in one sample),
    ``fell to baseline`` means it came back down to the run's background yet stayed
    inside the band, ``adjacent`` means nothing separated them, and ``length only``
    is the legacy path, which knows nothing about the level."""
    evidence = discriminator is not None
    if evidence:
        D = np.asarray(discriminator, dtype=np.float64)
        floor = discriminator_baseline(D) if baseline is None else float(baseline)
        if not np.isfinite(floor) or floor <= 0:
            floor = 1.0
        if high_gap is None:
            # 'no cap asked for' means whatever the cycle time implies
            high_limit = MERGE_HIGH_GAP_DEFAULT if cap is None else int(cap)
        else:
            high_limit = int(high_gap)
        limits = {
            "high": max(0, high_limit),
            "low": max(0, int(low_gap if low_gap is not None else 0)),
        }
    else:
        D = np.zeros(0)
        floor = 1.0
        high_limit = MERGE_HIGH_GAP_DEFAULT if high_gap is None else int(high_gap)
        low_limit = MERGE_LOW_GAP_DEFAULT if low_gap is None else int(low_gap)
        limits = {
            "high": max(0, high_limit),
            "low": max(0, low_limit),
        }
    if not any(limits.values()):
        return [dict(segment) for segment in segments]

    merged = []
    edges = []  # the level each merged interval ends with, which the next gap abuts
    stable = []  # plateau cycles only, so a level never counts the gaps it spans
    for segment in segments:
        current = dict(segment)
        current.setdefault("merged_segments", 1)
        current.setdefault("merged_gaps", [])
        if merged:
            previous = merged[-1]
            cls = current.get("class")
            limit = limits.get(cls, 0)
            gap = current["start_cycle"] - previous["end_cycle"] - 1
            record = None
            if previous.get("class") == cls and 0 <= gap <= limit and limit > 0:
                record = _gap_record(
                    previous,
                    current,
                    gap,
                    D,
                    floor,
                    band,
                    evidence,
                    previous_level=edges[-1],
                    lower_bound=cls == "high",
                )
            if record is not None:
                previous_cycles = stable[-1]
                current_cycles = current["n_cycles"]
                stable_cycles = previous_cycles + current_cycles
                previous["end_cycle"] = current["end_cycle"]
                previous["end_s"] = current["end_s"]
                previous["n_cycles"] = (
                    previous["end_cycle"] - previous["start_cycle"] + 1
                )
                previous["level"] = round(
                    (
                        previous["level"] * previous_cycles
                        + current["level"] * current_cycles
                    )
                    / stable_cycles,
                    2,
                )
                previous["merged_segments"] += current["merged_segments"]
                previous["merged_gaps"].append(record)
                previous["merged_gaps"].extend(current["merged_gaps"])
                edges[-1] = float(current.get("level") or 0.0)
                stable[-1] = stable_cycles
                continue
        merged.append(current)
        edges.append(float(current.get("level") or 0.0))
        stable.append(int(current["n_cycles"]))
    return merged


def merge_gaps_note(gaps):
    """Say in one line what a merge did, e.g. `joined 2 wobbles, level held
    (<= 28 cycles)` — the reviewer in the browser has no command line to ask with."""
    records = [gap for gap in gaps if isinstance(gap, dict)]
    if not records:
        return ""
    reasons = [str(gap.get("reason") or "no evidence recorded") for gap in records]
    distinct = list(dict.fromkeys(reasons))
    count = len(records)
    if len(distinct) == 1:
        why = distinct[0]
        noun = "wobble" if why == MERGE_REASON_HELD else "gap"
    else:
        why = ", ".join(f"{reasons.count(r)} {r}" for r in distinct)
        noun = "gap"
    widest = max(int(gap.get("cycles") or 0) for gap in records)
    return (
        f"joined {count} {noun}{'' if count == 1 else 's'}, {why} "
        f"(\u2264 {widest} cycles)"
    )


def merge_adjacent_high_segments(segments, max_gap_cycles=MERGE_HIGH_GAP_DEFAULT):
    """Back-compat: merge only high plateaus (see merge_adjacent_segments)."""
    return merge_adjacent_segments(segments, high_gap=max_gap_cycles, low_gap=0)


def spec_duration_s(f):
    try:
        duration = float(f.attrs["Single Spec Duration (ms)"][0]) / 1000.0
    except (IndexError, KeyError, OSError, TypeError, ValueError, OverflowError):
        return 1.0
    return duration if np.isfinite(duration) and duration > 0 else 1.0


def load_pc_times(f):
    """Return validated per-cycle PC Unix timestamps, or ``None``.

    PCTime is optional in IoniTOF exports.  It is only suitable for an axis when
    it contains exactly one finite, strictly increasing timestamp per spectrum.
    In particular, malformed timestamps must not be replaced with plausible
    absolute dates.
    """
    try:
        ncyc = int(f["SPECdata/Intensities"].shape[0])
        ds = f["SPECdata/PCTime"]
        if ds.ndim == 2 and ds.shape[1] == 1:
            values = np.asarray(ds[:, 0], dtype=np.float64)
        elif ds.ndim == 1:
            values = np.asarray(ds[:], dtype=np.float64)
        else:
            return None
    except (KeyError, OSError, TypeError, ValueError):
        return None
    if values.shape != (ncyc,) or not np.all(np.isfinite(values)):
        return None
    if ncyc > 1 and not np.all(np.diff(values) > 0):
        return None
    return values


def load_pc_timezone_offset(f):
    """Return the lab-PC UTC offset in seconds, or zero when unavailable.

    IONICON stores ``UTC_Offset`` as a root attribute.  It is metadata for
    displaying the PC Unix timestamps in the lab's wall-clock time; it does not
    affect relative elapsed time.
    """
    try:
        values = np.asarray(f.attrs["UTC_Offset"], dtype=np.float64).reshape(-1)
    except (KeyError, TypeError, ValueError):
        return 0.0
    if values.size != 1 or not np.isfinite(values[0]):
        return 0.0
    return float(values[0])


def viz_x_axis_data(f):
    """Return finite, browser-safe cycle, relative, and local-time axis data.

    PCTime is useful for the relative axis even when it is outside JavaScript's
    Date range.  In that case only the absolute lab-local axis is disabled.  A bad
    duration must not leak NaN, infinity, or a non-increasing domain into the
    browser, so the documented one-second fallback is used instead.
    """
    ncyc = int(f["SPECdata/Intensities"].shape[0])
    pctimes = load_pc_times(f)
    cycles = np.arange(ncyc, dtype=np.float64) + 1.0

    duration = spec_duration_s(f)
    if not np.isfinite(duration) or duration <= 0:
        duration = 1.0
    fallback = np.arange(ncyc, dtype=np.float64) * duration
    if not np.all(np.isfinite(fallback)) or (
        ncyc > 1 and not np.all(np.diff(fallback) > 0)
    ):
        fallback = np.arange(ncyc, dtype=np.float64)

    relative = fallback
    absolute = None
    timezone_offset = load_pc_timezone_offset(f)
    if pctimes is not None:
        candidate = pctimes - pctimes[0]
        if np.all(np.isfinite(candidate)) and (
            ncyc <= 1 or np.all(np.diff(candidate) > 0)
        ):
            relative = candidate
        # Keep absolute dates only where Date.toISOString() emits a normal
        # four-digit year.  Expanded years use a leading sign and break the
        # browser's current tick formatting.
        local_times = pctimes + timezone_offset
        if np.all(
            (local_times >= _JS_NORMAL_YEAR_0000_START_S)
            & (local_times < _JS_NORMAL_YEAR_10000_START_S)
        ):
            absolute = local_times

    def serialise_axis(values):
        # Do not round Unix seconds: at ordinary acquisition dates a millisecond
        # is already close to the precision of a JavaScript Number, and rounding
        # here would collapse valid adjacent sub-millisecond spectra.
        serialised = [float(x) for x in values]
        if not np.all(np.isfinite(serialised)) or (
            len(serialised) > 1 and not np.all(np.diff(serialised) > 0)
        ):
            return None
        return serialised

    relative_axis = serialise_axis(relative)
    absolute_axis = None if absolute is None else serialise_axis(absolute)
    return {
        "relative": relative_axis,
        "absolute": absolute_axis,
        "absolute_available": absolute_axis is not None,
        "absolute_offset_s": timezone_offset,
        "cycle": [int(x) for x in cycles],
    }


# ---------- quantification ----------
def stats(x):
    return {
        "Max": float(x.max()),
        "Min": float(x.min()),
        "Average": float(x.mean()),
        "Deviation": float(x.std(ddof=1)),
    }


def quantify(
    traces,
    f,
    ranges,
    K=None,
    primary=None,
    primary_mz=21.022,
    molar_volume=None,
    R_used=1200.0,
    k_map=None,
    k_anchor=K_ANCHOR_DEFAULT,
    humid_masses=None,
    humidity_ratio=None,
    humidity_ref=None,
    humidity_p=1.0,
    mass_axis=None,
    isotope_plan=None,
    isotope_abundance_basis="unknown",
):
    """Turn raw traces into Corrected / Conc / Conc[ug] and per-range statistics.

    Concentration uses the standard primary-ion-normalised model
        Conc[ppb] = Corrected * K / I_primary(t) * (k_anchor / k_compound)
    The last factor is the optional per-compound kinetic correction: without it
    (k_map None) every compound shares one effective rate constant and the output
    reproduces a single-sensitivity reference; with it, each compound is scaled by
    its own proton-transfer rate constant, which is physically more accurate.

    K:       None -> derived from the file's own pre-computed concentration; else
             a float (from `calibrate` / a standard).
    primary: per-cycle primary-ion signal; extracted from the file if None.
    k_map:   {mz: {'k': value_in_1e-9, ...}} from resolve_k(); None disables the
             kinetic correction.
    k_anchor: the single rate constant (1e-9 units) the baseline K assumes.
    """
    if mass_axis is not None:
        validate_mass_axis(mass_axis)
    tm, tf = load_transmission(f)
    if primary is None:
        primary = extract_primary(
            f, primary_mz=primary_mz, R=R_used, mass_axis=mass_axis
        )
    if molar_volume is None:
        molar_volume, molar_volume_source = derive_molar_volume_info(f)
    else:
        molar_volume_source = "configured"
    if K is None:
        K = derive_K(f, primary)

    # guard the divide; where primary is missing/zero, concentration is undefined
    if primary is not None and K is not None:
        pos = primary > 0
        norm = np.zeros_like(primary)
        norm[pos] = K / primary[pos]
    else:
        norm = None

    # humidity correction (per-cycle) for near-thermoneutral compounds
    humid_masses = set(humid_masses or [])
    humid_applied = False
    hfac = None
    if humid_masses and humidity_ratio is not None:
        if humidity_ref is None:
            good = np.isfinite(humidity_ratio) & (humidity_ratio > 0)
            humidity_ref = (
                float(np.median(humidity_ratio[good])) if good.any() else None
            )
        if humidity_ref:
            hfac = humidity_factor(humidity_ratio, humidity_ref, humidity_p)
            humid_applied = True

    corrected = {}
    for mass, (raw_trace, apex_mass) in traces.items():
        transmission = float(np.interp(apex_mass, tm, tf))
        corrected[mass] = raw_trace / transmission
    isotope_diagnostics = []
    if isotope_plan is not None:
        net_corrected, isotope_diagnostics = isotopes.correct_parent_signals(
            corrected,
            isotope_plan,
            abundance_basis=isotope_abundance_basis,
        )
        parent_masses = {
            float(mass) for mass in isotope_plan.get("analyte_masses", [])
        }
    else:
        net_corrected = corrected
        parent_masses = set(traces)

    rows = []
    for m, (raw, apex_m) in traces.items():
        if m not in parent_masses:
            continue
        T = float(np.interp(apex_m, tm, tf))
        cor = corrected[m]
        quantitative_signal = net_corrected.get(m, cor)
        kfac = 1.0
        # hybrid kinetic: only scale by a compound's own k when that k is a
        # measured value; compounds with an estimated k stay on the shared K.
        if k_map and k_map.get(m, {}).get("k") and not k_map[m].get("k_estimated"):
            kfac = k_anchor / float(k_map[m]["k"])
        if norm is not None:
            con = quantitative_signal * norm * kfac
            if hfac is not None and m in humid_masses:
                con = con * hfac
            ug = con * (m - PROTON) / molar_volume
        else:
            con = np.full_like(cor, np.nan)
            ug = con
        for label, (lo, hi) in ranges.items():
            s = slice(lo - 1, hi)  # 1-based inclusive cycle window
            rows.append(
                {
                    "mass": m,
                    "apex": apex_m,
                    "range": label,
                    "transmission": T,
                    "raw": stats(raw[s]),
                    "cor": stats(cor[s]),
                    "con": stats(con[s]),
                    "ug": stats(ug[s]),
                }
            )
    return rows, {
        "K": K,
        "molar_volume": molar_volume,
        "R": R_used,
        "primary_mz": primary_mz,
        "kinetic": k_map is not None,
        "k_anchor": k_anchor,
        "concentration_available": norm is not None,
        "transmission_available": has_transmission(f),
        "humidity_corrected": humid_applied,
        "humidity_ref": humidity_ref,
        "humidity_p": humidity_p if humid_applied else None,
        "molar_volume_source": molar_volume_source,
        "isotopes": {
            "enabled": isotope_plan is not None,
            "model": isotope_plan.get("version") if isotope_plan is not None else None,
            "abundance_basis": isotope_abundance_basis,
            "warnings": isotope_plan.get("warnings", [])
            if isotope_plan is not None
            else [],
            "corrections": isotope_diagnostics,
        },
    }


def calibrate_K(
    f,
    traces,
    ref_rows,
    ranges,
    primary=None,
    primary_mz=21.022,
    R_used=1200.0,
    mass_axis=None,
):
    """Fit the concentration constant K so output matches a reference.

    ref_rows: {(round(mz,3), range_label): reference_conc_ppb}. Returns
    (K, residual_median_pct, n_points). K is the median of
    ref_conc * I_primary / Corrected over all matched reference points."""
    if mass_axis is not None:
        validate_mass_axis(mass_axis)
    tm, tf = load_transmission(f)
    if primary is None:
        primary = extract_primary(
            f, primary_mz=primary_mz, R=R_used, mass_axis=mass_axis
        )
    ks, mine_cor = [], {}
    for m, (raw, apex_m) in traces.items():
        T = float(np.interp(apex_m, tm, tf))
        cor = raw / T
        for label, (lo, hi) in ranges.items():
            mine_cor[(round(m, 3), label)] = (
                cor[lo - 1 : hi].mean(),
                primary[lo - 1 : hi].mean(),
            )
    for key, ref_c in ref_rows.items():
        if key in mine_cor and ref_c:
            c, p = mine_cor[key]
            if c > 0 and p > 0:
                ks.append(ref_c * p / c)
    if not ks:
        return None, None, 0
    K = float(np.median(ks))
    resid = [abs(100 * (K - k) / k) for k in ks]
    return K, float(np.median(resid)), len(ks)
