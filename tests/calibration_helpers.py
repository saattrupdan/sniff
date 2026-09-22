"""Reusable mass-calibration fixtures."""

from sniff import ptrms


def identity_mass_axis(a=10.0, b=1.0):
    """Return a deliberately explicit identity correction for tiny fixtures."""
    anchors = [
        {
            "name": "water_cluster",
            "target_mz": 37.028405,
            "status": "accepted",
            "reason": "",
            "observed_file_mz": 37.028405,
            "corrected_mz": 37.028405,
            "timebin": a * 37.028405**0.5 + b,
            "prominence": 100.0,
            "snr": 100.0,
            "persistence": {
                "available": True,
                "blocks": 8,
                "accepted_blocks": 8,
                "fraction": 1.0,
                "statuses": ["accepted"] * 8,
                "block_centres_file_mz": [37.028405] * 8,
            },
        },
        {
            "name": "iodobenzene",
            "target_mz": 203.942993,
            "status": "accepted",
            "reason": "",
            "observed_file_mz": 203.942993,
            "corrected_mz": 203.942993,
            "timebin": a * 203.942993**0.5 + b,
            "prominence": 100.0,
            "snr": 100.0,
            "persistence": {
                "available": True,
                "blocks": 8,
                "accepted_blocks": 8,
                "fraction": 1.0,
                "statuses": ["accepted"] * 8,
                "block_centres_file_mz": [203.942993] * 8,
            },
        },
    ]
    diagnostics = {
        "model": "m_corrected = scale*m_file + offset",
        "applied": True,
        "scale": 1.0,
        "offset_da": 0.0,
        "fallback_reason": None,
        "file_calibration": {
            "model": "timebin = a*sqrt(m_file) + b",
            "a": float(a),
            "b": float(b),
        },
        "anchors": anchors,
        "formula_assignment_tolerance": {
            "model": ptrms.FORMULA_TOLERANCE_MODEL,
            "source": "synthetic exact calibration references",
            "status": "accepted",
            "reason": None,
            "mass_error_convention": (
                "1e6 * (observed - theoretical) / theoretical"
            ),
            "minimum_ppm": ptrms.FORMULA_TOLERANCE_FLOOR_PPM,
            "maximum_ppm": ptrms.FORMULA_TOLERANCE_MAX_PPM,
            "q95_abs_ppm": 0.0,
            "tolerance_ppm": ptrms.FORMULA_TOLERANCE_FLOOR_PPM,
            "score_sigma_ppm": ptrms.FORMULA_SCORE_SIGMA_FLOOR_PPM,
            "candidate_generation_allowed": True,
            "automatic_assignment_allowed": True,
            "calibration_points": [
                {"mz": 37.028405, "residual_ppm": 0.0},
                {"mz": 100.0, "residual_ppm": 0.0},
                {"mz": 203.942993, "residual_ppm": 0.0},
            ],
        },
    }
    return ptrms.MassAxisCalibration(a, b, diagnostics=diagnostics)
