"""Required internal calibration and legacy config migration tests."""

import json
from argparse import Namespace

import h5py
import numpy as np
import pytest

from sniff import analyze, ptrms

A = 1000.0
B = 0.0
SCALE = 1.0007
OFFSET = -0.020
NBIN = 15000


def _observed(corrected):
    return (corrected - OFFSET) / SCALE


def _spectrum(include_water=True, include_iodobenzene=True):
    spectrum = np.full(NBIN, 4.0, dtype=np.float64)
    for mass, height in (
        (_observed(37.028405), 1200.0),
        (_observed(203.942993), 1000.0),
    ):
        if mass == _observed(37.028405) and not include_water:
            continue
        if mass == _observed(203.942993) and not include_iodobenzene:
            continue
        centre = A * np.sqrt(mass) + B
        bins = np.arange(int(centre) - 5, int(centre) + 7, dtype=np.float64)
        shape = np.maximum(0.0, height * (1.0 - ((bins - centre) / 3.5) ** 2))
        spectrum[int(centre) - 5 : int(centre) + 7] += shape
    return spectrum


def _file(spectrum):
    handle = h5py.File("calibration", "w", driver="core", backing_store=False)
    handle.create_dataset("CALdata/Spectrum", data=np.array([[A, B]]))
    handle.create_dataset("SPECdata/AverageSpec", data=spectrum)
    handle.create_dataset("SPECdata/Intensities", data=np.vstack([spectrum] * 8))
    return handle


def valid_axis(scale=SCALE, offset=OFFSET):
    return ptrms.MassAxisCalibration(
        A,
        B,
        scale=scale,
        offset=offset,
        diagnostics={
            "model": "m_corrected = scale*m_file + offset",
            "applied": True,
            "scale": scale,
            "offset_da": offset,
            "file_calibration": {
                "model": "timebin = a*sqrt(m_file) + b",
                "a": A,
                "b": B,
            },
            "anchors": [
                {
                    "name": "water_cluster",
                    "target_mz": 37.028405,
                    "status": "accepted",
                    "reason": "",
                    "observed_file_mz": (37.028405 - offset) / scale,
                    "corrected_mz": 37.028405,
                    "timebin": A * ((37.028405 - offset) / scale) ** 0.5 + B,
                    "prominence": 100.0,
                    "snr": 100.0,
                    "persistence": {
                        "available": True,
                        "blocks": 8,
                        "accepted_blocks": 8,
                        "fraction": 1.0,
                        "statuses": ["accepted"] * 8,
                    },
                },
                {
                    "name": "iodobenzene",
                    "target_mz": 203.942993,
                    "status": "accepted",
                    "reason": "",
                    "observed_file_mz": (203.942993 - offset) / scale,
                    "corrected_mz": 203.942993,
                    "timebin": A * ((203.942993 - offset) / scale) ** 0.5 + B,
                    "prominence": 100.0,
                    "snr": 100.0,
                    "persistence": {
                        "available": True,
                        "blocks": 8,
                        "accepted_blocks": 8,
                        "fraction": 1.0,
                        "statuses": ["accepted"] * 8,
                    },
                },
            ],
        },
    )


def legacy_axis(scale=SCALE, offset=OFFSET):
    diagnostics = valid_axis(scale=scale, offset=offset).to_dict()
    targets = {"water_cluster": 37.033, "iodobenzene": 204.951}
    for anchor in diagnostics["anchors"]:
        target = targets[anchor["name"]]
        observed = (target - offset) / scale
        anchor.update(
            {
                "target_mz": target,
                "observed_file_mz": observed,
                "corrected_mz": target,
                "timebin": A * observed**0.5 + B,
            }
        )
    return ptrms.MassAxisCalibration(
        A,
        B,
        scale=scale,
        offset=offset,
        diagnostics=diagnostics,
    )


@pytest.mark.parametrize(
    "mutation",
    [
        lambda axis: axis.diagnostics.update(
            model="m_corrected = scale*m_file + offset + forged"
        ),
        lambda axis: axis.diagnostics["file_calibration"].update(
            model="timebin = a*sqrt(m_file) + forged"
        ),
        lambda axis: axis.diagnostics.update(fallback_reason="contradictory"),
        lambda axis: axis.diagnostics["anchors"][0].update(reason="failed"),
        lambda axis: axis.diagnostics["anchors"][0].update(prominence=0.0),
        lambda axis: axis.diagnostics["anchors"][1].update(snr=0.0),
        lambda axis: axis.diagnostics["anchors"][0].update(timebin=0.0),
        lambda axis: axis.diagnostics["anchors"][1].update(
            observed_file_mz=203.942993 + ptrms.INTERNAL_ANCHOR_SEARCH_DA + 0.001
        ),
        lambda axis: axis.diagnostics["anchors"][0].update(target_mz=37.034),
        lambda axis: axis.diagnostics["anchors"].__setitem__(
            1, axis.diagnostics["anchors"][0].copy()
        ),
    ],
    ids=[
        "correction-model",
        "file-model",
        "fallback-reason",
        "anchor-reason",
        "prominence",
        "snr",
        "timebin",
        "search-window",
        "target",
        "duplicate-anchor",
    ],
)
def test_mass_axis_validator_rejects_forged_anchor_evidence(mutation):
    axis = valid_axis()
    mutation(axis)
    with pytest.raises(ptrms.MassCalibrationError):
        ptrms.validate_mass_axis(axis)


@pytest.mark.parametrize(
    ("scale", "offset"),
    [(2.0, OFFSET), (SCALE, -37.028405)],
    ids=["reviewer-scale", "reviewer-offset"],
)
def test_mass_axis_validator_rejects_implausible_affine_correction(scale, offset):
    axis = valid_axis(scale=scale, offset=offset)
    with pytest.raises(ptrms.MassCalibrationError):
        ptrms.validate_mass_axis(axis)


def test_mass_axis_validator_rejects_anchor_near_search_boundary():
    axis = valid_axis(scale=1.0, offset=-0.18)
    with pytest.raises(ptrms.MassCalibrationError):
        ptrms.validate_mass_axis(axis)


def test_mass_axis_validator_rejects_bad_persistence_arithmetic():
    axis = valid_axis()
    persistence = axis.diagnostics["anchors"][0]["persistence"]
    persistence["accepted_blocks"] = 7
    persistence["fraction"] = 1.0
    with pytest.raises(ptrms.MassCalibrationError):
        ptrms.validate_mass_axis(axis)


def test_mass_axis_validator_rejects_nonfinite_and_unphysical_coefficients():
    for coefficient, value in (("a", 0.0), ("a", np.inf), ("b", np.nan)):
        axis = valid_axis()
        setattr(axis, coefficient, value)
        with pytest.raises(ptrms.MassCalibrationError):
            ptrms.validate_mass_axis(axis)


def test_mass_axis_validator_rejects_missing_or_misordered_anchor_identity():
    axis = valid_axis()
    axis.diagnostics["anchors"] = list(reversed(axis.diagnostics["anchors"]))
    # Names make a deliberately reversed list unambiguous and therefore valid.
    ptrms.validate_mass_axis(axis)

    axis = valid_axis()
    axis.diagnostics["anchors"][0]["name"] = "not_water"
    with pytest.raises(ptrms.MassCalibrationError):
        ptrms.validate_mass_axis(axis)


def test_missing_required_anchor_is_a_structured_calibration_error():
    with (
        _file(_spectrum(include_iodobenzene=False)) as handle,
        pytest.raises(ptrms.MassCalibrationError) as caught,
    ):
        ptrms.load_mass_axis(handle)

    error = caught.value
    assert "mass calibration failed" in str(error)
    assert error.diagnostics["anchors"][1]["name"] == "iodobenzene"
    assert error.diagnostics["anchors"][1]["status"] != "accepted"


def test_raw_cycle_persistence_is_required_when_available():
    spectrum = _spectrum()
    with _file(spectrum) as handle:
        raw = handle["SPECdata/Intensities"]
        raw[0:4, :] = 4.0
        with pytest.raises(ptrms.MassCalibrationError) as caught:
            ptrms.load_mass_axis(handle)

    assert "persistence" in caught.value.diagnostics["anchors"][1]
    assert caught.value.diagnostics["anchors"][1]["persistence"]["accepted_blocks"] < 5


def test_old_config_migrates_absolute_masses_and_widths_once():
    axis = valid_axis()
    old = {
        "peaks": [
            {
                "mz": 100.0,
                "window": {"left": 0.1, "right": 0.2},
                "unknown": {"keep": True},
            },
            {"mz": 100.0, "window": 0.4},
        ],
        "analyze": {"primary_mz": 21.022},
        "ranges": [{"label": "sample_01", "start": 1, "end": 2}],
        "extra": "preserved",
    }

    migrated, changed = ptrms.migrate_config_mass_axis(old, axis)
    again, changed_again = ptrms.migrate_config_mass_axis(migrated, axis)

    assert changed is True
    assert changed_again is False
    assert again == migrated
    assert migrated["mass_axis_domain"] == "corrected"
    assert migrated["mass_axis_version"] == 3
    assert migrated["peaks"][0]["mz"] == pytest.approx(100.05)
    assert migrated["peaks"][0]["window"]["left"] == pytest.approx(0.10007)
    assert migrated["peaks"][0]["window"]["right"] == pytest.approx(0.20014)
    assert migrated["peaks"][1]["window"] == pytest.approx(0.40028)
    assert migrated["analyze"]["primary_mz"] == pytest.approx(21.0167146)
    assert migrated["extra"] == "preserved"
    json.dumps(migrated)


def test_migration_refuses_an_uncalibrated_axis():
    axis = ptrms.MassAxisCalibration(
        A, B, diagnostics={"applied": False, "fallback_reason": "missing"}
    )
    with pytest.raises(ptrms.MassCalibrationError, match="not an applied"):
        ptrms.migrate_config_mass_axis({"peaks": [{"mz": 42.0}]}, axis)


def test_direct_cli_config_migration_persists_the_marker(tmp_path):
    axis = valid_axis()
    path = tmp_path / "config.json"
    path.write_text(json.dumps({"peaks": [{"mz": 100.0}]}), encoding="utf-8")
    args = Namespace(config=str(path))

    migrated = analyze._migrate_loaded_config(
        json.loads(path.read_text(encoding="utf-8")), args, axis
    )

    saved = json.loads(path.read_text(encoding="utf-8"))
    assert migrated == saved
    assert saved["mass_axis_domain"] == "corrected"
    assert saved["mass_axis_version"] == 3
    assert saved["peaks"][0]["mz"] == pytest.approx(100.05)


def test_version_one_config_preserves_timebins_on_axis_upgrade():
    old_axis = legacy_axis()
    new_axis = valid_axis(scale=1.0, offset=0.0)
    historical = old_axis.to_dict()
    historical["reference_stability"] = {
        "source": "historical target-relative values",
        "references": [{"median_ppm": 123.0}],
    }
    old = {
        "peaks": [{"mz": old_axis.file_to_corrected(100.0), "window": 0.2}],
        "mass_axis_domain": "corrected",
        "mass_axis_version": 1,
        "mass_axis_calibration": historical,
    }

    migrated, changed = ptrms.migrate_config_mass_axis(old, new_axis)

    assert changed is True
    assert migrated["mass_axis_version"] == 3
    assert migrated["peaks"][0]["mz"] == pytest.approx(100.0)
    assert migrated["peaks"][0]["window"] == pytest.approx(0.2 / old_axis.scale)
    assert migrated["mass_axis_calibration"] == new_axis.to_dict()


def test_version_two_affine_config_preserves_timebins_on_axis_upgrade():
    old_axis = legacy_axis()
    new_axis = valid_axis(scale=1.0, offset=0.0)
    old = {
        "peaks": [{"mz": old_axis.file_to_corrected(100.0), "window": 0.2}],
        "mass_axis_domain": "corrected",
        "mass_axis_version": 2,
        "mass_axis_calibration": old_axis.to_dict(),
    }

    migrated, changed = ptrms.migrate_config_mass_axis(old, new_axis)

    assert changed is True
    assert migrated["mass_axis_version"] == 3
    assert migrated["peaks"][0]["mz"] == pytest.approx(100.0)
    assert migrated["peaks"][0]["window"] == pytest.approx(0.2 / old_axis.scale)


def test_version_two_mapping_config_migrates_from_file_mass_domain():
    new_axis = valid_axis()
    mapping_diagnostics = {
        "model": ptrms.MAPPING_AUTHORITY_MODEL,
        "authority": "CALdata/Mapping",
        "applied": True,
        "scale": 1.0,
        "offset_da": 0.0,
        "file_calibration": {
            "model": ptrms.FILE_MASS_CALIBRATION_MODEL,
            "a": A,
            "b": B,
        },
    }
    old = {
        "peaks": [{"mz": 100.0, "formula": "C2H2O"}],
        "mass_axis_domain": "corrected",
        "mass_axis_version": 2,
        "mass_axis_calibration": mapping_diagnostics,
    }

    migrated, changed = ptrms.migrate_config_mass_axis(old, new_axis)

    assert changed is True
    assert migrated["peaks"][0]["mz"] == pytest.approx(
        new_axis.file_to_corrected(100.0)
    )


def test_cacheless_version_two_config_refuses_ambiguous_axis_upgrade():
    old = {
        "peaks": [{"mz": 100.0, "formula": "C2H2O"}],
        "mass_axis_domain": "corrected",
        "mass_axis_version": 2,
    }

    with pytest.raises(ValueError, match="version-2 corrected config lacks"):
        ptrms.migrate_config_mass_axis(old, valid_axis())


def test_cacheless_version_one_config_refuses_ambiguous_anchor_upgrade():
    axis = valid_axis(scale=1.0, offset=0.0)
    old = {
        "peaks": [{"mz": 42.0, "formula": "C2H2O"}],
        "mass_axis_domain": "corrected",
        "mass_axis_version": 1,
    }

    with pytest.raises(ValueError, match="lacks reconstructable historical"):
        ptrms.migrate_config_mass_axis(old, axis)


def test_identity_migration_only_adds_the_marker():
    axis = valid_axis(scale=1.0, offset=0.0)
    old = {"peaks": [{"mz": 42.0, "window": 0.2}]}
    migrated, changed = ptrms.migrate_config_mass_axis(old, axis)

    assert changed is True
    assert migrated["peaks"] == old["peaks"]
    assert migrated["mass_axis_domain"] == "corrected"
    assert migrated["mass_axis_version"] == 3
