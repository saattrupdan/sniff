"""Tests for the persistent review app: config discovery, the deterministic
bootstrap, the agent hand-off and the routes the start screen drives."""

import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import urllib.error
import urllib.request
from pathlib import Path
from unittest import mock

import h5py
import numpy as np
import pytest
from calibration_helpers import identity_mass_axis

from sniff import app, viz


@pytest.fixture(autouse=True)
def _synthetic_files_have_explicit_calibration():
    with mock.patch.object(
        app.ptrms, "load_mass_axis", return_value=identity_mass_axis()
    ):
        yield


def make_h5(path, *, cycles=4, mz_count=8, signal=True):
    """A minimal IoniTOF-shaped file: enough for the app, not for the science."""
    intensities = np.ones((cycles, mz_count)) if signal else np.zeros((cycles, mz_count))
    with h5py.File(path, "w") as h5:
        h5.create_dataset("SPECdata/Intensities", data=intensities)
        h5.create_dataset("SPECdata/AverageSpec", data=np.ones(mz_count) if signal else np.zeros(mz_count))
        h5.attrs["InstrumentType"] = "test"


def payload_stub(f, peaks, ranges, **kwargs):
    return {
        "file": "stub.h5",
        "peaks": [],
        "ranges": [],
        "meta": {"ncyc": 4},
        "config_base": {},
    }


# --------------------------------------------------------------------------
# where the config lives
# --------------------------------------------------------------------------
def test_default_recent_path_uses_sniff_state_directory(tmp_path):
    with (
        mock.patch.dict(os.environ, {}, clear=True),
        mock.patch.object(Path, "home", return_value=tmp_path),
    ):
        assert app._recent_path() == tmp_path / ".sniff" / "recent.json"
        assert app._active_path() == tmp_path / ".sniff" / "active.json"
        assert app._onboarding_path() == tmp_path / ".sniff" / "onboarding.json"


def test_frozen_windowed_log_works_without_stderr(tmp_path, monkeypatch):
    monkeypatch.setenv("SNIFF_RECENT_PATH", str(tmp_path / "recent.json"))
    monkeypatch.setattr(sys, "frozen", True, raising=False)
    monkeypatch.setattr(sys, "stderr", None)

    app._log("sniff: windowed diagnostic")

    assert (tmp_path / "log.txt").read_text() == "sniff: windowed diagnostic\n"


def test_frozen_background_failure_uses_windowed_log(tmp_path, monkeypatch):
    class ImmediateThread:
        def __init__(self, target, daemon):
            self.target = target

        def start(self):
            self.target()

    monkeypatch.setenv("SNIFF_RECENT_PATH", str(tmp_path / "recent.json"))
    monkeypatch.setattr(sys, "frozen", True, raising=False)
    monkeypatch.setattr(sys, "stderr", None)
    monkeypatch.setattr(app.threading, "Thread", ImmediateThread)

    app._background(mock.Mock(side_effect=RuntimeError("open failed")))

    assert (tmp_path / "log.txt").read_text() == "sniff: RuntimeError: open failed\n"


def test_active_file_survives_shutdown_until_the_reviewer_leaves(tmp_path):
    h5 = tmp_path / "run.h5"
    h5.touch()
    active = tmp_path / "state" / "active.json"
    with mock.patch.object(app, "ACTIVE_PATH", active):
        app.remember_active(str(h5))
        assert app.load_active() == str(h5.resolve())
        app.forget_active()
        assert app.load_active() is None


def test_automatic_tour_state_survives_app_launches(tmp_path):
    onboarding = tmp_path / "state" / "onboarding.json"
    with mock.patch.object(app, "ONBOARDING_PATH", onboarding):
        assert app.auto_tour_pending() is True
        app.remember_tour_seen()
        assert app.auto_tour_pending() is False
    assert onboarding.is_file()


def test_failed_onboarding_directory_does_not_look_completed(tmp_path):
    blocked = tmp_path / "blocked"
    blocked.write_text("not a directory", encoding="utf-8")
    with mock.patch.object(app, "ONBOARDING_PATH", blocked / "onboarding.json"):
        assert app.auto_tour_pending() is True
        with pytest.raises(OSError):
            app.remember_tour_seen()


def test_valid_config_beside_the_file_is_used(tmp_path):
    h5 = tmp_path / "sniff.h5"
    h5.touch()
    beside = tmp_path / "sniff.json"
    beside.write_text(json.dumps({"peaks": [{"mz": 42.0}]}), encoding="utf-8")
    assert app.config_path_for(str(h5)) == beside


def test_unrelated_same_stem_json_is_never_clobbered(tmp_path):
    h5 = tmp_path / "sniff.h5"
    h5.touch()
    (tmp_path / "sniff.json").write_text('{"unrelated": true}', encoding="utf-8")
    target = app.config_path_for(str(h5))
    assert target == tmp_path / "sniff.sniff.json"
    assert json.loads((tmp_path / "sniff.json").read_text(encoding="utf-8")) == {
        "unrelated": True
    }


def test_legacy_config_name_is_honoured(tmp_path):
    h5 = tmp_path / "run.h5"
    h5.touch()
    legacy = tmp_path / "run-analysis-config.json"
    legacy.write_text(json.dumps({"ranges": []}), encoding="utf-8")
    assert app.config_path_for(str(h5)) == legacy


def test_missing_config_targets_the_stem_json(tmp_path):
    h5 = tmp_path / "run.h5"
    h5.touch()
    assert app.config_path_for(str(h5)) == tmp_path / "run.json"


def test_legacy_product_config_is_recognised_without_clobbering_neighbours(tmp_path):
    h5 = tmp_path / "run.h5"
    h5.touch()
    beside = tmp_path / "run.json"
    beside.write_text('{"owned_by": "someone else"}', encoding="utf-8")
    legacy = tmp_path / "run.ptr.json"
    legacy.write_text(json.dumps({"peaks": []}), encoding="utf-8")
    assert app.config_path_for(str(h5)) == legacy
    assert json.loads(beside.read_text(encoding="utf-8")) == {
        "owned_by": "someone else"
    }


def test_open_reuses_calibration_for_the_same_unchanged_h5(tmp_path):
    h5 = tmp_path / "run.h5"
    make_h5(h5)
    session = app.Session()
    stages = []
    assignment_defaults = []

    def observed_payload(f, peaks, ranges, **kwargs):
        stages.append(session.stage)
        assignment_defaults.append(kwargs.get("assign_identity_defaults"))
        return payload_stub(f, peaks, ranges, **kwargs)

    with (
        mock.patch.object(
            app.ptrms, "load_mass_axis", return_value=identity_mass_axis()
        ) as load_axis,
        mock.patch.object(app, "auto_peaks", return_value=[]),
        mock.patch.object(app, "auto_ranges", return_value=[]),
        mock.patch.object(app.viz, "build_viz_data", observed_payload),
    ):
        session.open(str(h5))
        assert load_axis.call_count == 1
        saved = json.loads((tmp_path / "run.json").read_text(encoding="utf-8"))
        assert saved["mass_axis_calibration"]["applied"] is True
        assert saved["mass_axis_h5_fingerprint"] == app._h5_fingerprint(h5)

        load_axis.reset_mock()
        session.open(str(h5))
        load_axis.assert_not_called()

    assert stages == [
        "Computing data for a new review",
        "Loading H5 data for the saved review",
    ]
    assert assignment_defaults == [True, False]


def test_changed_h5_invalidates_saved_calibration(tmp_path):
    h5 = tmp_path / "run.h5"
    make_h5(h5)
    session = app.Session()
    with (
        mock.patch.object(
            app.ptrms, "load_mass_axis", return_value=identity_mass_axis()
        ) as load_axis,
        mock.patch.object(app, "auto_peaks", return_value=[]),
        mock.patch.object(app, "auto_ranges", return_value=[]),
        mock.patch.object(app.viz, "build_viz_data", payload_stub),
    ):
        session.open(str(h5))
        old_fingerprint = app._h5_fingerprint(h5)
        session.close()
        with h5py.File(h5, "r+") as source:
            source.attrs["changed"] = True
        assert app._h5_fingerprint(h5) != old_fingerprint

        load_axis.reset_mock()
        session.open(str(h5))
        assert load_axis.call_count == 1


def test_replaced_h5_during_open_does_not_keep_generated_config(tmp_path):
    h5 = tmp_path / "run.h5"
    replacement = tmp_path / "replacement.h5"
    make_h5(h5, signal=True)
    make_h5(replacement, cycles=5, signal=False)
    config_path = tmp_path / "run.json"

    def replace_during_payload(f, peaks, ranges, **kwargs):
        os.replace(replacement, h5)
        return payload_stub(f, peaks, ranges, **kwargs)

    session = app.Session()
    with (
        mock.patch.object(
            app.ptrms, "load_mass_axis", return_value=identity_mass_axis()
        ),
        mock.patch.object(app, "auto_peaks", return_value=[{"mz": 42.0}]),
        mock.patch.object(app, "auto_ranges", return_value=[]),
        mock.patch.object(app.viz, "build_viz_data", replace_during_payload),
        pytest.raises(RuntimeError, match="H5 file changed"),
    ):
        session.open(str(h5))

    assert not config_path.exists()


def test_replaced_h5_during_config_publication_removes_generated_config(tmp_path):
    h5 = tmp_path / "run.h5"
    replacement = tmp_path / "replacement.h5"
    make_h5(h5, signal=True)
    make_h5(replacement, cycles=5, signal=False)
    config_path = tmp_path / "run.json"
    real_replace = os.replace
    source_replaced = False

    def replace_during_publish(source, destination):
        nonlocal source_replaced
        if Path(destination) == config_path and not source_replaced:
            real_replace(replacement, h5)
            source_replaced = True
        return real_replace(source, destination)

    session = app.Session()
    with (
        mock.patch.object(
            app.ptrms, "load_mass_axis", return_value=identity_mass_axis()
        ),
        mock.patch.object(app, "auto_peaks", return_value=[{"mz": 42.0}]),
        mock.patch.object(app, "auto_ranges", return_value=[]),
        mock.patch.object(app.viz, "build_viz_data", payload_stub),
        mock.patch.object(app.os, "replace", side_effect=replace_during_publish),
        pytest.raises(RuntimeError, match="H5 file changed"),
    ):
        session.open(str(h5))

    assert source_replaced
    assert not config_path.exists()


def test_invalid_saved_calibration_is_recomputed(tmp_path):
    h5 = tmp_path / "run.h5"
    make_h5(h5)
    session = app.Session()
    with (
        mock.patch.object(
            app.ptrms, "load_mass_axis", return_value=identity_mass_axis()
        ) as load_axis,
        mock.patch.object(app, "auto_peaks", return_value=[]),
        mock.patch.object(app, "auto_ranges", return_value=[]),
        mock.patch.object(app.viz, "build_viz_data", payload_stub),
    ):
        session.open(str(h5))
        session.close()
        config_path = tmp_path / "run.json"
        config = json.loads(config_path.read_text(encoding="utf-8"))
        config["mass_axis_calibration"]["scale"] = 1.01
        config_path.write_text(json.dumps(config), encoding="utf-8")

        load_axis.reset_mock()
        session.open(str(h5))
        assert load_axis.call_count == 1

        session.close()
        config = json.loads(config_path.read_text(encoding="utf-8"))
        config["mass_axis_calibration"]["anchors"][0]["persistence"]["fraction"] = (
            10**1000
        )
        config_path.write_text(json.dumps(config), encoding="utf-8")
        load_axis.reset_mock()
        session.open(str(h5))
        assert load_axis.call_count == 1

    refreshed = json.loads(config_path.read_text(encoding="utf-8"))
    assert refreshed["mass_axis_calibration"]["scale"] == 1.0


def test_config_written_on_open_is_reread_on_the_next_open(tmp_path):
    h5 = tmp_path / "run.h5"
    make_h5(h5)
    session = app.Session()
    with (
        mock.patch.object(app, "auto_peaks", return_value=[{"mz": 42.0}]),
        mock.patch.object(app, "auto_ranges", return_value=[{"label": "sample_01", "start": 1, "end": 2}]),
        mock.patch.object(app.viz, "build_viz_data", payload_stub),
    ):
        session.open(str(h5))
    written = json.loads((tmp_path / "run.json").read_text(encoding="utf-8"))
    assert written["peaks"] and written["ranges"]
    assert "checklist" not in written
    assert app.load_active() == str(h5.resolve())
    with mock.patch.object(app, "auto_peaks") as detect_peaks:
        with mock.patch.object(app.viz, "build_viz_data", payload_stub):
            session.open(str(h5))
    detect_peaks.assert_not_called()  # the saved config is reused, not recomputed


def test_legacy_recent_store_is_read_and_copied_on_next_write(tmp_path):
    recent = tmp_path / "new" / "recent.json"
    legacy = tmp_path / "old" / "recent.json"
    legacy.parent.mkdir()
    legacy.write_text(json.dumps(["/data/old.h5"]), encoding="utf-8")
    with (
        mock.patch.object(app, "RECENT_PATH", recent),
        mock.patch.object(app, "LEGACY_RECENT_PATH", legacy),
        mock.patch.object(app, "_recent_path", return_value=recent),
    ):
        assert app.load_recent() == ["/data/old.h5"]
        app.remember_recent("/data/new.h5")
    assert json.loads(recent.read_text(encoding="utf-8")) == [
        str(Path("/data/new.h5").resolve()),
        "/data/old.h5",
    ]
    assert json.loads(legacy.read_text(encoding="utf-8")) == ["/data/old.h5"]


# --------------------------------------------------------------------------
# the deterministic bootstrap
# --------------------------------------------------------------------------
def test_bootstrap_reports_what_it_decided(tmp_path):
    h5 = tmp_path / "run.h5"
    make_h5(h5)
    with (
        mock.patch.object(app, "auto_peaks", return_value=[{"mz": 42.0}]),
        mock.patch.object(app, "auto_ranges", return_value=[{"label": "sample_01", "start": 1, "end": 2}]),
    ):
        config = app.bootstrap_config(str(h5))
    assert config["peaks"] and config["ranges"]
    assert "checklist" not in config
    assert config["diagnostics"]["n_peaks"] == 1
    assert config["diagnostics"]["n_ranges"] == 1


def test_a_blank_file_says_so(tmp_path):
    h5 = tmp_path / "blank.h5"
    make_h5(h5, signal=False)
    with (
        mock.patch.object(app, "auto_peaks", return_value=[]),
        mock.patch.object(app, "auto_ranges", return_value=[]),
    ):
        config = app.bootstrap_config(str(h5))
    assert "checklist" not in config


# --------------------------------------------------------------------------
# the optional agent endpoint
# --------------------------------------------------------------------------
def _agent_answer(config):
    class Response:
        def read(self):
            return json.dumps(config).encode("utf-8")

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    return Response()


def test_existing_config_scrubs_retired_review_fields_and_keeps_unknowns(tmp_path):
    h5 = tmp_path / "run.h5"
    make_h5(h5)
    config_path = tmp_path / "run.json"
    config_path.write_text(
        json.dumps(
            {
                "peaks": [{"mz": 42.0}],
                "ranges": [],
                "checklist": ["old top-level note"],
                "review": {
                    "checklist": ["old nested note"],
                    "future": {"keep": True},
                },
                "future_field": {"keep": True},
            }
        ),
        encoding="utf-8",
    )
    session = app.Session()
    with mock.patch.object(app.viz, "build_viz_data", payload_stub):
        session.open(str(h5))
    saved = json.loads(config_path.read_text(encoding="utf-8"))
    assert "checklist" not in saved
    assert saved["review"] == {"future": {"keep": True}}
    assert saved["future_field"] == {"keep": True}

    saved["checklist"] = ["retired field added to cached config"]
    saved["review"]["checklist"] = ["retired nested field"]
    config_path.write_text(json.dumps(saved), encoding="utf-8")
    with (
        mock.patch.object(app.ptrms, "load_mass_axis") as load_axis,
        mock.patch.object(app.viz, "build_viz_data", payload_stub),
    ):
        session.open(str(h5))
    load_axis.assert_not_called()
    cached = json.loads(config_path.read_text(encoding="utf-8"))
    assert "checklist" not in cached
    assert cached["review"] == {"future": {"keep": True}}
    assert cached["future_field"] == {"keep": True}

    cleaned, changed = app._scrub_retired_review_fields(
        {"peaks": [], "review": {"checklist": ["only retired field"]}}
    )
    assert changed and "review" not in cleaned


def test_agent_answer_replaces_the_deterministic_config(tmp_path):
    h5 = tmp_path / "run.h5"
    make_h5(h5)
    curated = {"peaks": [{"mz": 7.0, "label": "N2"}], "ranges": []}
    session = app.Session()
    with (
        mock.patch.object(app, "auto_peaks", return_value=[{"mz": 42.0}]),
        mock.patch.object(app, "auto_ranges", return_value=[]),
        mock.patch.object(app.viz, "build_viz_data", payload_stub),
        mock.patch.object(app.urllib.request, "urlopen", return_value=_agent_answer(curated)),
    ):
        session.open(str(h5), agent_url="http://agent.invalid/curate")
    assert session.config["peaks"] == curated["peaks"]
    assert session.agent_status == "Agent review applied."
    on_disk = json.loads((tmp_path / "run.json").read_text(encoding="utf-8"))
    assert on_disk["peaks"] == curated["peaks"]
    assert on_disk["mass_axis_calibration"]["applied"] is True
    assert on_disk["mass_axis_h5_fingerprint"] == app._h5_fingerprint(h5)


def test_unreachable_agent_keeps_the_deterministic_config(tmp_path):
    h5 = tmp_path / "run.h5"
    make_h5(h5)
    session = app.Session()
    with (
        mock.patch.object(app, "auto_peaks", return_value=[{"mz": 42.0}]),
        mock.patch.object(app, "auto_ranges", return_value=[]),
        mock.patch.object(app.viz, "build_viz_data", payload_stub),
        mock.patch.object(
            app.urllib.request, "urlopen", side_effect=urllib.error.URLError("offline")
        ),
    ):
        session.open(str(h5), agent_url="http://agent.invalid/curate")
    assert session.status == "ready"
    assert session.config["peaks"] == [{"mz": 42.0}]
    assert session.agent_status.startswith("Agent review failed")
    on_disk = json.loads((tmp_path / "run.json").read_text(encoding="utf-8"))
    assert on_disk["peaks"] == [{"mz": 42.0}]


def test_agent_reply_without_peaks_is_refused(tmp_path):
    h5 = tmp_path / "run.h5"
    make_h5(h5)
    session = app.Session()
    with (
        mock.patch.object(app, "auto_peaks", return_value=[{"mz": 42.0}]),
        mock.patch.object(app, "auto_ranges", return_value=[]),
        mock.patch.object(app.viz, "build_viz_data", payload_stub),
        mock.patch.object(
            app.urllib.request, "urlopen", return_value=_agent_answer({"note": "hi"})
        ),
    ):
        session.open(str(h5), agent_url="http://agent.invalid/curate")
    assert session.config["peaks"] == [{"mz": 42.0}]
    assert "no peaks or ranges" in session.agent_status


# --------------------------------------------------------------------------
# recents
# --------------------------------------------------------------------------
def test_recents_are_ordered_capped_and_tolerant(tmp_path):
    recent = tmp_path / "recent.json"
    with mock.patch.object(app, "RECENT_PATH", recent):
        for i in range(25):
            app.remember_recent(str(tmp_path / f"{i}.h5"))
        values = app.load_recent()
        assert len(values) == 20
        assert values[0].endswith("24.h5")
        recent.write_text("{", encoding="utf-8")
        assert app.load_recent() == []


def test_remembering_a_file_twice_does_not_duplicate_it(tmp_path):
    recent = tmp_path / "recent.json"
    with mock.patch.object(app, "RECENT_PATH", recent):
        app.remember_recent(str(tmp_path / "a.h5"))
        app.remember_recent(str(tmp_path / "b.h5"))
        app.remember_recent(str(tmp_path / "a.h5"))
        assert app.load_recent() == [str(tmp_path / "a.h5"), str(tmp_path / "b.h5")]


# --------------------------------------------------------------------------
# exporting
# --------------------------------------------------------------------------
def test_export_writes_beside_the_file_and_stays_ready(tmp_path):
    h5 = tmp_path / "run.h5"
    make_h5(h5)
    session = app.Session()
    session.path = str(h5)
    session.config = {"peaks": [], "ranges": []}
    session.status = "ready"

    def fake_analysis(path, config, out):
        Path(out).write_text("compound,cycle\n", encoding="utf-8")
        return {"n_rows": 1}

    with mock.patch.object(app, "analyze_config_to_csv", side_effect=fake_analysis):
        result = session.export()
    assert Path(result["out"]) == tmp_path / "run.csv"
    assert session.status == "ready"
    assert session.status_payload() == {"status": "done", "out": str(tmp_path / "run.csv")}


def test_export_setup_failure_releases_the_busy_reservation(tmp_path, monkeypatch):
    h5 = tmp_path / "run.h5"
    session = app.Session()
    session.path = str(h5)
    session.config = {"peaks": [{"mz": 42.0}], "ranges": []}
    session.status = "ready"
    monkeypatch.setattr(app.tempfile, "mkstemp", mock.Mock(side_effect=OSError("disk full")))

    with pytest.raises(OSError, match="disk full"):
        session.export()

    assert session.busy is False
    assert session.status == "ready"
    assert session.export_error == "disk full"


def test_export_uses_the_admitted_request_snapshot(server, tmp_path, monkeypatch):
    api, session = server
    h5 = tmp_path / "run.h5"
    h5.touch()
    session.path = str(h5)
    session.config_path = tmp_path / "run.json"
    session.config = {}
    session.status = "ready"
    clicked = {
        "peaks": [{"mz": 42.0}],
        "ranges": [],
        "mass_axis_domain": "corrected",
        "mass_axis_version": 2,
    }
    later = {**clicked, "peaks": [{"mz": 99.0}]}
    jobs = []
    monkeypatch.setattr(app, "_background", lambda *args, **kwargs: jobs.append((args, kwargs)))

    assert api.post("/export?version=10", clicked)[0] == 202
    assert api.post("/save?version=11", later)[0] == 200
    assert session.config == later
    assert jobs[0][1]["config"] == clicked
    analysed = []
    monkeypatch.setattr(
        app,
        "analyze_config_to_csv",
        lambda path, config, out: analysed.append((path, config, out)) or {},
    )
    function = jobs[0][0][0]
    function(**jobs[0][1])

    assert analysed[0][0] == str(h5)
    assert analysed[0][1] == clicked
    assert session.busy is False


def test_a_second_export_cannot_report_the_previous_one(tmp_path):
    h5 = tmp_path / "run.h5"
    make_h5(h5)
    session = app.Session()
    session.path = str(h5)
    session.config = {"peaks": [], "ranges": []}
    session.status = "ready"
    session.export_result = {"out": str(tmp_path / "run.csv")}
    with mock.patch.object(app, "analyze_config_to_csv", side_effect=RuntimeError("boom")):
        with pytest.raises(RuntimeError):
            session.export()
    assert session.status_payload()["status"] == "error"
    assert session.status == "ready"  # the file is still open and still reviewable


# --------------------------------------------------------------------------
# the routes
# --------------------------------------------------------------------------
def _review_token(page):
    match = re.search(rb'const PAGE_TOKEN = "([^"]+)"', page)
    assert match is not None
    return match.group(1).decode("ascii")


class _Server:
    def __init__(self, base):
        self.base = base

    def get(self, path, headers=None):
        request = urllib.request.Request(self.base + path, headers=headers or {})
        with urllib.request.urlopen(request, timeout=10) as response:
            return response.status, response.read()

    def post(self, path, body=None):
        data = json.dumps(body or {}).encode("utf-8")
        req = urllib.request.Request(
            self.base + path, data, {"Content-Type": "application/json"}, method="POST"
        )
        try:
            with urllib.request.urlopen(req, timeout=10) as r:
                return r.status, json.loads(r.read().decode("utf-8"))
        except urllib.error.HTTPError as exc:
            return exc.code, json.loads(exc.read().decode("utf-8"))


@pytest.fixture
def server(tmp_path):
    httpd, session, url = app.make_server(port=0)  # the OS picks a free port
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    try:
        yield _Server(url), session
    finally:
        # An open runs on its own thread; letting it outlive the test would have it
        # write to a path this test no longer owns.
        for _ in range(200):
            if not session.busy:
                break
            threading.Event().wait(0.02)
        httpd.shutdown()
        thread.join(timeout=5)
        httpd.server_close()
        session.close()


def _wait_ready(api, timeout=20):
    spent = 0.0
    while spent < timeout:
        _, state = api.get("/api/state")
        if json.loads(state)["status"] in ("ready", "error"):
            return json.loads(state)
        threading.Event().wait(0.05)
        spent += 0.05
    raise AssertionError("the file never finished opening")


def test_the_app_says_how_it_is_showing_itself(server):
    """/api/state has to be able to tell a window from a tab.

    A bundle that asked for a window and silently got a tab leaves the user with no
    evidence: a windowed macOS build writes nothing to a console, so the only honest
    report is the one the running app gives. It is read live, not snapshotted, because
    the window either opened or it did not and the page must be able to say which.
    """
    api, session = server
    _, body = api.get("/api/state")
    assert json.loads(body)["surface"] == "browser"

    previous = app.surface()
    app._surface = "window"
    try:
        _, body = api.get("/api/state")
        assert json.loads(body)["surface"] == "window"
    finally:
        app._surface = previous


def test_update_route_checks_once_and_opens_the_verified_installer(server, monkeypatch):
    api, _ = server
    available = app.updates.Update(
        version="1.2.0",
        asset_name="sniff-review-macos-arm64.pkg",
        download_url=(
            "https://github.com/saattrupdan/sniff/releases/download/"
            "v1.2.0/sniff-review-macos-arm64.pkg"
        ),
        digest="sha256:" + "a" * 64,
        release_url="https://github.com/saattrupdan/sniff/releases/tag/v1.2.0",
    )
    checks = []
    opened = []
    monkeypatch.setattr(
        app.updates,
        "find_update",
        lambda version: checks.append(version) or available,
    )
    monkeypatch.setattr(
        app.updates,
        "install_update",
        lambda candidate: opened.append(candidate),
    )

    status, body = api.get("/api/update")
    assert status == 200
    assert json.loads(body) == {
        "available": True,
        "version": "1.2.0",
        "release_url": "https://github.com/saattrupdan/sniff/releases/tag/v1.2.0",
    }
    assert api.get("/api/update")[0] == 200
    status, body = api.post("/update")
    assert status == 200 and body == {"ok": True, "version": "1.2.0"}
    assert checks == [app.__version__]
    assert opened == [available]
    assert json.loads(api.get("/api/update")[1]) == {"available": False}


def test_concurrent_update_requests_cannot_share_a_download(server, monkeypatch):
    api, _ = server
    available = app.updates.Update(
        version="1.2.0",
        asset_name="sniff-review-macos-arm64.pkg",
        download_url=(
            "https://github.com/saattrupdan/sniff/releases/download/"
            "v1.2.0/sniff-review-macos-arm64.pkg"
        ),
        digest="sha256:" + "a" * 64,
        release_url="https://github.com/saattrupdan/sniff/releases/tag/v1.2.0",
    )
    started = threading.Event()
    release = threading.Event()
    first_response = []
    monkeypatch.setattr(app.updates, "find_update", lambda _version: available)

    def install(_candidate):
        started.set()
        assert release.wait(5)

    monkeypatch.setattr(app.updates, "install_update", install)
    thread = threading.Thread(
        target=lambda: first_response.append(api.post("/update")), daemon=True
    )
    thread.start()
    assert started.wait(5)

    second_status, second_body = api.post("/update")
    release.set()
    thread.join(timeout=5)

    assert second_status == 409
    assert second_body == {"error": "an update is already being opened"}
    assert first_response == [(200, {"ok": True, "version": "1.2.0"})]


def test_install_later_defers_the_update_for_this_app_launch(server, monkeypatch):
    api, _ = server
    available = mock.Mock(version="1.2.0", release_url="https://example.invalid")
    monkeypatch.setattr(app.updates, "find_update", lambda _version: available)

    assert json.loads(api.get("/api/update")[1])["available"] is True
    assert api.post("/update/later") == (200, {"ok": True})
    assert json.loads(api.get("/api/update")[1]) == {"available": False}


@pytest.mark.parametrize(
    "origin",
    ["https://example.com", "http://127.0.0.1:bad", "http://[broken"],
)
def test_cross_origin_websites_cannot_open_the_installer(server, origin):
    api, _ = server
    request = urllib.request.Request(
        api.base + "update",
        data=b"{}",
        headers={"Content-Type": "application/json", "Origin": origin},
        method="POST",
    )

    with pytest.raises(urllib.error.HTTPError) as caught:
        urllib.request.urlopen(request, timeout=10)

    assert caught.value.code == 403
    assert json.loads(caught.value.read()) == {
        "error": "cross-origin update request refused"
    }


def test_update_check_failure_is_quiet_for_the_start_screen(server, monkeypatch):
    api, _ = server

    def fail(_version):
        raise OSError("offline")

    monkeypatch.setattr(app.updates, "find_update", fail)

    status, body = api.get("/api/update")

    assert status == 200
    assert json.loads(body) == {"available": False}


def test_start_screen_offers_an_update_or_a_later_choice():
    html = app._START_HTML

    assert 'id="install-update"' in html and "Install update" in html
    assert 'id="install-later"' in html and "Install later" in html
    assert "fetch('/api/update')" in html
    assert "fetch('/update',{method:'POST'})" in html
    assert "fetch('/update/later',{method:'POST'})" in html


def test_start_screen_opens_files_without_rendering_recents(server, tmp_path, monkeypatch):
    api, session = server
    h5 = tmp_path / "run.h5"
    make_h5(h5)
    status, body = api.get("/")
    assert status == 200 and b"Sniff" in body and b"PTR-MS review" in body
    assert b'<html lang="en">' in body
    assert b"--line2:" not in body and b"--ok:" not in body and b"--cream:" not in body
    assert b"/api/recent" not in body and b'id="recent"' not in body
    assert b"Find the story in your spectrum." in body

    with (
        mock.patch.object(app, "auto_peaks", return_value=[{"mz": 42.0}]),
        mock.patch.object(app, "auto_ranges", return_value=[]),
        mock.patch.object(app.viz, "build_viz_data", payload_stub),
    ):
        code, _ = api.post("/open", {"path": str(h5)})
        assert code == 202
        state = _wait_ready(api)
    assert state["status"] == "ready"
    assert state["config"] == str(tmp_path / "run.json")

    _, recent = api.get("/api/recent")
    assert json.loads(recent)[0]["path"] == str(h5)


def test_open_validates_and_stores_compounds_of_interest(server, tmp_path):
    api, _ = server
    h5 = tmp_path / "run.h5"
    make_h5(h5)
    with (
        mock.patch.object(app, "auto_peaks", return_value=[]) as detect_peaks,
        mock.patch.object(app, "auto_ranges", return_value=[]),
        mock.patch.object(app.viz, "build_viz_data", payload_stub),
    ):
        code, _ = api.post(
            "/open",
            {"path": str(h5), "compounds_of_interest": ["Acetone", "propanal"]},
        )
        assert code == 202
        state = _wait_ready(api)

    assert state["status"] == "ready"
    assert detect_peaks.call_args.kwargs["compounds_of_interest"] == [
        {"name": "acetone", "formula": "C3H6O", "mz": 59.0491},
        {"name": "propanal", "formula": "C3H6O", "mz": 59.0491},
    ]
    config = json.loads((tmp_path / "run.json").read_text(encoding="utf-8"))
    assert config["compounds_of_interest"] == [
        {"name": "acetone", "formula": "C3H6O", "mz": 59.0491},
        {"name": "propanal", "formula": "C3H6O", "mz": 59.0491},
    ]


def test_open_rejects_an_unknown_compound_before_starting(server, tmp_path):
    api, session = server
    h5 = tmp_path / "run.h5"
    make_h5(h5)

    code, body = api.post(
        "/open", {"path": str(h5), "compounds_of_interest": ["not a compound"]}
    )

    assert code == 400
    assert body["error"] == "unrecognised compound: not a compound"
    assert session.state()["status"] == "empty"
    assert not (tmp_path / "run.json").exists()


def test_save_canonicalises_compound_priors_and_rejects_unknowns(tmp_path):
    session = app.Session()
    session.config_path = tmp_path / "run.json"
    session.config = {"peaks": []}
    page = session.begin_review_page()

    assert session.save_config(
        {
            "peaks": [],
            "compounds_of_interest": [
                {"name": "Acetone", "formula": "forged", "mz": -1}
            ],
        },
        page=page,
    )
    saved = json.loads(session.config_path.read_text(encoding="utf-8"))
    assert saved["compounds_of_interest"] == [
        {"name": "acetone", "formula": "C3H6O", "mz": 59.0491}
    ]

    with pytest.raises(ValueError, match="unrecognised compound"):
        session.save_config(
            {"peaks": [], "compounds_of_interest": ["invented compound"]},
            page=page,
        )
    assert json.loads(session.config_path.read_text(encoding="utf-8")) == saved


def test_review_page_404s_until_a_file_is_open(server):
    api, _ = server
    with pytest.raises(urllib.error.HTTPError) as exc:
        api.get("/review")
    assert exc.value.code == 404


def test_onboarding_write_failure_is_reported_without_consuming_state(
    server, tmp_path, monkeypatch
):
    api, _ = server
    blocked = tmp_path / "blocked"
    blocked.write_text("not a directory", encoding="utf-8")
    monkeypatch.setattr(app, "ONBOARDING_PATH", blocked / "onboarding.json")

    code, body = api.post("/onboarding", {})
    assert code == 500
    assert "could not save onboarding state" in body["error"]
    assert app.auto_tour_pending() is True


def test_review_page_is_rendered_in_app_mode(server, tmp_path, monkeypatch):
    api, _ = server
    monkeypatch.setattr(app, "ONBOARDING_PATH", tmp_path / "onboarding.json")
    h5 = tmp_path / "run.h5"
    make_h5(h5)
    with (
        mock.patch.object(app, "auto_peaks", return_value=[]),
        mock.patch.object(app, "auto_ranges", return_value=[]),
        mock.patch.object(app.viz, "build_viz_data", payload_stub),
    ):
        api.post("/open", {"path": str(h5)})
        _wait_ready(api)
    status, html = api.get("/review")
    assert status == 200
    assert b"const APPMODE = true" in html  # Export, not Done
    assert b"const AUTO_TOUR = true" in html
    first_token = _review_token(html)
    assert first_token.startswith("1-") and len(first_token) > 20
    assert b"Open another file" in html
    # Merely rendering the page does not consume a tour that has not appeared yet.
    _, refreshed = api.get("/review")
    assert b"const AUTO_TOUR = true" in refreshed
    second_token = _review_token(refreshed)
    assert second_token.startswith("2-") and second_token != first_token

    code, _ = api.post("/onboarding", {})
    assert code == 200
    _, acknowledged = api.get("/review")
    assert b"const AUTO_TOUR = false" in acknowledged
    third_token = _review_token(acknowledged)
    assert third_token.startswith("3-") and third_token not in {first_token, second_token}


def test_review_render_cannot_mix_payload_and_target_across_an_open(
    server, tmp_path, monkeypatch
):
    api, session = server
    first = tmp_path / "first.h5"
    second = tmp_path / "second.h5"
    make_h5(first)
    make_h5(second)

    def named_payload(file, peaks, ranges, **kwargs):
        return {
            "file": Path(file.filename).name,
            "peaks": peaks,
            "ranges": ranges,
            "meta": {"ncyc": 4},
            "config_base": kwargs["config_base"],
        }

    with (
        mock.patch.object(app, "auto_peaks", return_value=[]),
        mock.patch.object(app, "auto_ranges", return_value=[]),
        mock.patch.object(app.viz, "build_viz_data", side_effect=named_payload),
    ):
        api.post("/open", {"path": str(first)})
        _wait_ready(api)
        first_page = api.get("/review")[1]
        first_token = _review_token(first_page)
        assert first_token.startswith("1-")
        session.finish_close_save(first_token)
        payload_captured = threading.Event()
        release_render = threading.Event()
        real_review_payload = session.review_payload

        def blocked_review_payload():
            payload = real_review_payload()
            payload_captured.set()
            assert release_render.wait(5)
            return payload

        monkeypatch.setattr(session, "review_payload", blocked_review_payload)
        response = {}

        def refresh_first():
            try:
                response["value"] = api.get("/review")
            except urllib.error.HTTPError as exc:
                response["error"] = exc.code

        refresh = threading.Thread(target=refresh_first)
        refresh.start()
        assert payload_captured.wait(5)
        api.post("/open", {"path": str(second)})
        _wait_ready(api)
        release_render.set()
        refresh.join(5)
        assert response == {"error": 409}
        monkeypatch.setattr(session, "review_payload", real_review_payload)
        _, current = api.get("/review")

    assert b'second.h5' in current
    assert str(tmp_path / "second.json").encode("utf-8") in current
    assert _review_token(current).startswith("2-")


def test_save_reload_and_close_preserve_the_latest_config(server, tmp_path):
    api, _ = server
    h5 = tmp_path / "run.h5"
    make_h5(h5)

    def config_payload(_file, peaks, ranges, **kwargs):
        return {
            "file": "stub.h5",
            "peaks": peaks,
            "ranges": ranges,
            "meta": {"ncyc": 4},
            "config_base": kwargs["config_base"],
        }

    with (
        mock.patch.object(app, "auto_peaks", return_value=[]),
        mock.patch.object(app, "auto_ranges", return_value=[]),
        mock.patch.object(app.viz, "build_viz_data", side_effect=config_payload),
    ):
        api.post("/open", {"path": str(h5)})
        _wait_ready(api)
        _, first = api.get("/review")
        first_token = _review_token(first)
        assert first_token.startswith("1-")
        edited = {
            "peaks": [{"mz": 42.0, "label": "reloaded edit"}],
            "ranges": [],
            "mass_axis_domain": "corrected",
            "mass_axis_version": 2,
        }
        assert api.post(
            f"/save?version=10&page={first_token}&closing=1", edited
        )[0] == 200
        _, refreshed = api.get("/review")

    second_token = _review_token(refreshed)
    assert second_token.startswith("2-") and second_token != first_token
    assert b"reloaded edit" in refreshed
    stale = {**edited, "peaks": [{"mz": 1.0, "label": "stale"}]}
    assert api.post(
        f"/save?version=99&page={first_token}&closing=1", stale
    )[0] == 409
    assert api.post(
        f"/save?version=11&page={second_token}&closing=1", edited
    )[0] == 200
    saved = json.loads((tmp_path / "run.json").read_text(encoding="utf-8"))
    assert {key: saved[key] for key in edited} == edited
    assert saved["mass_axis_calibration"]["applied"] is True
    assert saved["mass_axis_h5_fingerprint"] == app._h5_fingerprint(h5)


def test_stale_page_cannot_write_or_export_a_newly_opened_file(
    server, tmp_path, monkeypatch
):
    api, _ = server
    first = tmp_path / "first.h5"
    second = tmp_path / "second.h5"
    make_h5(first)
    make_h5(second)
    with (
        mock.patch.object(app, "auto_peaks", return_value=[]),
        mock.patch.object(app, "auto_ranges", return_value=[]),
        mock.patch.object(app.viz, "build_viz_data", payload_stub),
    ):
        api.post("/open", {"path": str(first)})
        _wait_ready(api)
        _, first_page = api.get("/review")
        first_token = _review_token(first_page)
        assert first_token.startswith("1-")
        api.post("/open", {"path": str(second)})
        _wait_ready(api)
        _, second_page = api.get("/review")
        second_token = _review_token(second_page)
        assert second_token.startswith("2-") and second_token != first_token

    second_config = tmp_path / "second.json"
    before = json.loads(second_config.read_text(encoding="utf-8"))
    stale = {
        "peaks": [{"mz": 99.0}],
        "ranges": [],
        "mass_axis_domain": "corrected",
        "mass_axis_version": 2,
    }
    analysed = mock.Mock()
    monkeypatch.setattr(app, "analyze_config_to_csv", analysed)
    assert api.post(
        f"/save?version=999&page={first_token}&closing=1", stale
    )[0] == 409
    assert api.post(f"/export?version=1000&page={first_token}", stale)[0] == 409
    assert json.loads(second_config.read_text(encoding="utf-8")) == before
    analysed.assert_not_called()


def test_save_rejects_a_body_that_is_not_a_config(server, tmp_path):
    api, _ = server
    h5 = tmp_path / "run.h5"
    make_h5(h5)
    with (
        mock.patch.object(app, "auto_peaks", return_value=[]),
        mock.patch.object(app, "auto_ranges", return_value=[]),
        mock.patch.object(app.viz, "build_viz_data", payload_stub),
    ):
        api.post("/open", {"path": str(h5)})
        _wait_ready(api)
    code, _ = api.post("/save", {"nothing": True})
    assert code == 400
    code, _ = api.post("/save", {"peaks": [{"mz": 1.0}]})
    assert code == 400
    valid = {
        "peaks": [{"mz": 1.0}],
        "ranges": [],
        "mass_axis_domain": "corrected",
        "mass_axis_version": 2,
    }
    code, _ = api.post("/save", valid)
    assert code == 200
    on_disk = json.loads((tmp_path / "run.json").read_text(encoding="utf-8"))
    assert on_disk["peaks"] == valid["peaks"]


def test_open_reports_a_missing_file_without_touching_the_session(server):
    api, session = server
    code, body = api.post("/open", {"path": "/nope/nothing.h5"})
    assert code == 404
    assert session.status == "empty"


def test_export_keeps_the_server_and_the_file_open(server, tmp_path):
    api, session = server
    h5 = tmp_path / "run.h5"
    make_h5(h5)
    with (
        mock.patch.object(app, "auto_peaks", return_value=[]),
        mock.patch.object(app, "auto_ranges", return_value=[]),
        mock.patch.object(app.viz, "build_viz_data", payload_stub),
    ):
        api.post("/open", {"path": str(h5)})
        _wait_ready(api)

    def fake_analysis(path, config, out):
        Path(out).write_text("compound\n", encoding="utf-8")
        return {"n_rows": 0}

    with mock.patch.object(app, "analyze_config_to_csv", side_effect=fake_analysis):
        code, _ = api.post("/export")
        assert code == 202
        for _ in range(200):
            _, st = api.get("/status")
            if json.loads(st)["status"] in ("done", "error"):
                break
            threading.Event().wait(0.02)
    assert json.loads(st)["status"] == "done"
    assert session.path == str(h5)  # still open, still reviewable
    _, page = api.get("/")
    assert b"Sniff" in page and b"PTR-MS review" in page  # and the app is still serving


def test_close_waits_for_ready_to_be_published_atomically(server, tmp_path):
    api, session = server
    h5 = tmp_path / "run.h5"
    make_h5(h5)
    publishing = threading.Event()
    release = threading.Event()
    close_waiting = threading.Event()

    class PausingLock:
        def __init__(self, lock):
            self.lock = lock

        def __enter__(self):
            if publishing.is_set() and not release.is_set():
                close_waiting.set()
            self.lock.acquire()
            return self

        def __exit__(self, *_args):
            if session.status == "ready" and not session._opening:
                publishing.set()
                assert release.wait(5)
            self.lock.release()

    session._lock = PausingLock(session._lock)
    with (
        mock.patch.object(app, "auto_peaks", return_value=[]),
        mock.patch.object(app, "auto_ranges", return_value=[]),
        mock.patch.object(app.viz, "build_viz_data", payload_stub),
    ):
        assert api.post("/open", {"path": str(h5)})[0] == 202
        assert publishing.wait(5)
        response = {}

        def close_review():
            response["value"] = api.post("/close")

        close_thread = threading.Thread(target=close_review, name="close-request")
        close_thread.start()
        assert close_waiting.wait(5)
        assert close_thread.is_alive()
        release.set()
        close_thread.join(5)

    assert response["value"][0] == 200
    assert session.status == "empty"
    assert session.path is None
    assert app.load_active() is None


def test_routes_reserve_background_work_before_returning_accepted(
    server, tmp_path, monkeypatch
):
    api, session = server
    first = tmp_path / "first.h5"
    second = tmp_path / "second.h5"
    first.touch()
    second.touch()
    jobs = []
    monkeypatch.setattr(app, "_background", lambda *args, **kwargs: jobs.append((args, kwargs)))

    assert api.post("/open", {"path": str(first)})[0] == 202
    assert api.post("/open", {"path": str(second)})[0] == 409
    assert api.post("/close")[0] == 409
    assert len(jobs) == 1
    session._opening = False

    session.path = str(first)
    session.config_path = tmp_path / "first.json"
    session.config = {
        "peaks": [{"mz": 42.0}],
        "ranges": [],
        "mass_axis_domain": "corrected",
        "mass_axis_version": 2,
    }
    session.status = "ready"
    assert api.post("/export", session.config)[0] == 202
    assert api.post("/export", session.config)[0] == 409
    assert api.post("/open", {"path": str(second)})[0] == 409
    assert api.post("/close")[0] == 409
    assert len(jobs) == 2
    session._exporting = False
    session.status = "ready"


def test_a_client_that_acts_on_ready_is_never_told_busy(server, tmp_path):
    """The page offers its buttons as soon as the state says ready, so ready has to
    mean the work is really finished — not finished except for the busy flag."""
    api, session = server
    h5 = tmp_path / "run.h5"
    make_h5(h5)
    with (
        mock.patch.object(app, "auto_peaks", return_value=[]),
        mock.patch.object(app, "auto_ranges", return_value=[]),
        mock.patch.object(app.viz, "build_viz_data", payload_stub),
    ):
        api.post("/open", {"path": str(h5)})
        _wait_ready(api)
    assert session.busy is False
    assert api.post("/close")[0] == 200


def test_closing_a_file_leaves_the_server_up(server, tmp_path):
    api, session = server
    h5 = tmp_path / "run.h5"
    make_h5(h5)
    with (
        mock.patch.object(app, "auto_peaks", return_value=[]),
        mock.patch.object(app, "auto_ranges", return_value=[]),
        mock.patch.object(app.viz, "build_viz_data", payload_stub),
    ):
        api.post("/open", {"path": str(h5)})
        _wait_ready(api)
    assert app.load_active() == str(h5.resolve())
    code, _ = api.post("/close")
    assert code == 200 and session.status == "empty"
    assert app.load_active() is None
    with pytest.raises(urllib.error.HTTPError) as exc:
        api.get("/review")
    assert exc.value.code == 404


def test_opening_another_file_closes_the_first_one(tmp_path):
    a, b = tmp_path / "a.h5", tmp_path / "b.h5"
    make_h5(a)
    make_h5(b)
    session = app.Session()
    handles = []
    real_open = h5py.File

    class TrackedFile:
        def __init__(self, path):
            self._f = real_open(path, "r")
            self.closed = False
            handles.append(self)

        def __getattr__(self, name):
            return getattr(self._f, name)

        def __getitem__(self, key):
            return self._f[key]

        def close(self):
            self.closed = True
            self._f.close()

    with (
        mock.patch.object(app, "auto_peaks", return_value=[]),
        mock.patch.object(app, "auto_ranges", return_value=[]),
        mock.patch.object(app.viz, "build_viz_data", payload_stub),
        mock.patch.object(app.h5py, "File", lambda p, *a, **k: TrackedFile(p)),
    ):
        session.open(str(a))
        first = handles[0]
        session.open(str(b))
    assert first.closed
    assert session.path == str(b)


# --------------------------------------------------------------------------
# what the review page posts must be what the CSV describes
# --------------------------------------------------------------------------
def test_export_uses_the_config_the_page_posted(server, tmp_path):
    api, _ = server
    h5 = tmp_path / "run.h5"
    make_h5(h5)
    with (
        mock.patch.object(app, "auto_peaks", return_value=[]),
        mock.patch.object(app, "auto_ranges", return_value=[]),
        mock.patch.object(app.viz, "build_viz_data", payload_stub),
    ):
        api.post("/open", {"path": str(h5)})
        _wait_ready(api)

    seen = {}

    def fake_analysis(path, config, out):
        seen["peaks"] = [p["mz"] for p in config["peaks"]]
        Path(out).write_text("compound\n", encoding="utf-8")
        return {"n_rows": 0}

    edited = {
        "peaks": [{"mz": 99.0, "label": "edited by the reviewer"}],
        "ranges": [],
        "mass_axis_domain": "corrected",
        "mass_axis_version": 2,
    }
    with mock.patch.object(app, "analyze_config_to_csv", side_effect=fake_analysis):
        code, _ = api.post("/export", edited)
        assert code == 202
        for _ in range(200):
            _, st = api.get("/status")
            if json.loads(st)["status"] in ("done", "error"):
                break
            threading.Event().wait(0.02)
    assert seen["peaks"] == [99.0]  # not whatever the last autosave left on disk
    on_disk = json.loads((tmp_path / "run.json").read_text(encoding="utf-8"))
    assert on_disk["peaks"][0]["mz"] == 99.0


def test_export_never_overwrites_someone_elses_table(tmp_path):
    h5 = tmp_path / "run.h5"
    make_h5(h5)
    foreign = tmp_path / "run.csv"
    foreign.write_text("Reviewers own table, NOT a sniff output\n1,2,3\n", encoding="utf-8")
    session = app.Session()
    session.path = str(h5)
    session.config = {"peaks": [], "ranges": []}
    session.status = "ready"

    def fake_analysis(path, config, out):
        Path(out).write_text("File;Variable;Range\n", encoding="utf-8")
        return {}

    with mock.patch.object(app, "analyze_config_to_csv", side_effect=fake_analysis):
        result = session.export()
    assert Path(result["out"]) == tmp_path / "run-sniff.csv"
    assert foreign.read_text(encoding="utf-8").startswith("Reviewers own table")


def test_export_replaces_its_own_previous_summary(tmp_path):
    h5 = tmp_path / "run.h5"
    make_h5(h5)
    ours = tmp_path / "run.csv"
    ours.write_text("File;Variable;Range;Average(Corrected)\nold,row\n", encoding="utf-8")
    session = app.Session()
    session.path = str(h5)
    session.config = {"peaks": [], "ranges": []}
    session.status = "ready"

    def fake_analysis(path, config, out):
        Path(out).write_text(
            "File;Variable;Range;Average(Corrected)\nnew,row\n", encoding="utf-8"
        )
        return {}

    with mock.patch.object(app, "analyze_config_to_csv", side_effect=fake_analysis):
        result = session.export()
    assert Path(result["out"]) == ours
    assert "new,row" in ours.read_text(encoding="utf-8")
    assert not (tmp_path / "run-sniff.csv").exists()


def test_export_leaves_no_temp_files_behind(tmp_path):
    h5 = tmp_path / "run.h5"
    make_h5(h5)
    session = app.Session()
    session.path = str(h5)
    session.config = {"peaks": [], "ranges": []}
    session.status = "ready"
    with mock.patch.object(app, "analyze_config_to_csv", side_effect=RuntimeError("boom")):
        with pytest.raises(RuntimeError):
            session.export()
    assert list(tmp_path.glob("*.tmp")) == []


def test_export_is_refused_while_one_is_running(server, tmp_path):
    api, session = server
    h5 = tmp_path / "run.h5"
    make_h5(h5)
    with (
        mock.patch.object(app, "auto_peaks", return_value=[]),
        mock.patch.object(app, "auto_ranges", return_value=[]),
        mock.patch.object(app.viz, "build_viz_data", payload_stub),
    ):
        api.post("/open", {"path": str(h5)})
        _wait_ready(api)
    session.status = "exporting"
    session._exporting = True
    code, _ = api.post(
        "/export",
        {
            "peaks": [{"mz": 1.0}],
            "ranges": [],
            "mass_axis_domain": "corrected",
            "mass_axis_version": 2,
        },
    )
    assert code == 409
    session._exporting = False
    session.status = "ready"


# --------------------------------------------------------------------------
# durability of the saved config
# --------------------------------------------------------------------------
def test_older_delayed_save_cannot_replace_the_close_time_snapshot(
    server, tmp_path, monkeypatch
):
    api, session = server
    session.path = str(tmp_path / "run.h5")
    session.config_path = tmp_path / "run.json"
    session.config = {}
    session.status = "ready"
    base = {
        "ranges": [],
        "mass_axis_domain": "corrected",
        "mass_axis_version": 2,
    }
    newer = {**base, "peaks": [{"mz": 2.0}]}
    older = {**base, "peaks": [{"mz": 1.0}]}
    newer_writing = threading.Event()
    release_newer = threading.Event()
    older_arrived = threading.Event()
    real_write = app._write_json
    real_save = session.save_config

    def blocked_write(path, value):
        if value.get("peaks") == newer["peaks"]:
            newer_writing.set()
            assert release_newer.wait(5)
        real_write(path, value)

    def observed_save(config, version=None, page=None):
        if version == 100:
            older_arrived.set()
        return real_save(config, version=version, page=page)

    monkeypatch.setattr(app, "_write_json", blocked_write)
    session.save_config = observed_save
    responses = {}
    new_thread = threading.Thread(
        target=lambda: responses.setdefault("new", api.post("/save?version=200", newer))
    )
    old_thread = threading.Thread(
        target=lambda: responses.setdefault("old", api.post("/save?version=100", older))
    )
    new_thread.start()
    assert newer_writing.wait(5)
    old_thread.start()
    assert older_arrived.wait(5)
    release_newer.set()
    new_thread.join(5)
    old_thread.join(5)

    assert responses["new"][0] == 200 and responses["new"][1]["saved"] is True
    assert responses["old"][0] == 200 and responses["old"][1]["saved"] is False
    saved = json.loads(session.config_path.read_text(encoding="utf-8"))
    assert saved["peaks"] == [{"mz": 2.0}]


def test_one_shot_done_rejects_a_snapshot_older_than_autosave(tmp_path):
    config_path = tmp_path / "run.json"
    analysed = []
    result = {}
    ready = threading.Event()
    urls = []

    def publish_url(url):
        urls.append(url)
        ready.set()

    def run_server():
        result["value"] = viz.serve(
            "<html></html>",
            str(config_path),
            port=0,
            timeout=5,
            open_browser=False,
            run_analysis=lambda config: analysed.append(config) or {},
            ready_callback=publish_url,
        )

    thread = threading.Thread(target=run_server, daemon=True)
    thread.start()
    assert ready.wait(5)
    api = _Server(urls[0].rstrip("/"))

    newer = {"peaks": [{"mz": 2.0}], "ranges": []}
    older = {"peaks": [{"mz": 1.0}], "ranges": []}
    assert api.post("/save?version=200", newer)[0] == 200
    assert api.post("/done?version=100", older)[0] == 409
    assert analysed == []
    assert json.loads(config_path.read_text(encoding="utf-8")) == newer
    assert api.post("/done?version=300", newer)[0] == 200
    assert api.post("/ack")[0] == 200
    thread.join(timeout=5)
    assert not thread.is_alive()
    assert analysed == [newer]
    assert result["value"][0] == newer


def test_two_tabs_saving_at_once_do_not_publish_each_other(tmp_path):
    target = tmp_path / "cfg.json"
    errors = []

    def save(worker):
        for _ in range(80):
            try:
                app._write_json(target, {"peaks": [{"mz": 1.0}], "who": worker})
            except OSError as exc:
                errors.append(exc)

    threads = [threading.Thread(target=save, args=(w,)) for w in ("a", "b")]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert errors == []
    assert not list(tmp_path.glob("*.tmp"))
    assert "who" in json.loads(target.read_text(encoding="utf-8"))


def test_a_broken_recents_file_cannot_lose_an_opened_file(tmp_path, monkeypatch):
    h5 = tmp_path / "run.h5"
    make_h5(h5)
    session = app.Session()

    def explode(path):
        raise OSError("read-only home")

    monkeypatch.setattr(app, "remember_recent", explode)
    with (
        mock.patch.object(app, "auto_peaks", return_value=[]),
        mock.patch.object(app, "auto_ranges", return_value=[]),
        mock.patch.object(app.viz, "build_viz_data", payload_stub),
    ):
        session.open(str(h5))
    assert session.status == "ready"
    assert session.path == str(h5)
    assert app.load_active() == str(h5.resolve())


def test_close_forgets_the_last_export(tmp_path):
    h5 = tmp_path / "run.h5"
    make_h5(h5)
    session = app.Session()
    session.path = str(h5)
    session.config = {"peaks": [], "ranges": []}
    session.status = "ready"
    session.export_result = {"out": str(tmp_path / "run.csv")}
    session.close()
    assert session.status_payload() == {"status": "idle"}


def test_recents_ignore_entries_that_are_not_paths(tmp_path):
    (tmp_path / "recent.json").write_text('[123, null, "/x/run.h5"]', encoding="utf-8")
    with mock.patch.object(app, "RECENT_PATH", tmp_path / "recent.json"):
        assert app.load_recent() == ["/x/run.h5"]


def test_an_agent_answer_of_empty_lists_is_refused(tmp_path):
    h5 = tmp_path / "run.h5"
    make_h5(h5)
    session = app.Session()
    with (
        mock.patch.object(app, "auto_peaks", return_value=[{"mz": 42.0}]),
        mock.patch.object(app, "auto_ranges", return_value=[]),
        mock.patch.object(app.viz, "build_viz_data", payload_stub),
        mock.patch.object(
            app.urllib.request, "urlopen", return_value=_agent_answer({"peaks": []})
        ),
    ):
        session.open(str(h5), agent_url="http://agent.invalid/curate")
    assert session.config["peaks"] == [{"mz": 42.0}]
    assert "no peaks or ranges" in session.agent_status


def _review_js(mode="app"):
    data = {"file": "x.h5", "peaks": [], "ranges": [], "meta": {}}
    page = viz.render_html(data, mode=mode)
    return "\n".join(re.findall(r"<script>(.*?)</script>", page, re.S))


def test_review_header_handles_windows_paths_and_omits_interval_note():
    data = {
        "peaks": [],
        "ranges": [],
        "merge_note": "joined 2 gaps, 1 fell to baseline",
        "meta": {"file": r"C:\runs\sample.h5"},
    }
    page = viz.render_html(data, mode="app")

    assert 'id="intnote"' not in page
    assert r'String(M.file||"").split(/[\\/]/).pop()' in page


def test_the_export_dialog_hands_the_review_back():
    """Export used to end on a disabled button reading "Opened - you can close this
    tab", which did nothing and left the modal as the only screen. The dialog has to
    get out of the way once the CSV is revealed, and must be escapable without
    revealing at all, in the failure case too."""
    js = _review_js()
    assert 'id="keepreviewing"' in js, "no way back to the review from the dialog"
    assert ".remove()" in js and 'getElementById("doneov")' in js, (
        "the dialog is never dismissed"
    )
    # revealing the CSV must return to the review before anything is said about the
    # tab, or the user is left holding a dead button
    start = js.index('fetch(isExport?"/reveal":"/open"')
    end = js.index("else { ob.disabled=false", start)
    handler = js[start:end]
    assert handler.index("closeOverlay()") < handler.index("Opened"), (
        "app mode still ends on a relabelled button instead of closing the dialog"
    )
    assert "you can close this tab" in handler, "the one-shot flow lost its wording"
    failure = js[js.index('st.status==="error"'):js.index('st.status==="error"') + 1200]
    assert "keepreviewing" in failure, "a failed export leaves no way out of the modal"


@pytest.mark.skipif(shutil.which("node") is None, reason="node is not installed")
def test_the_embedded_scripts_are_valid_javascript():
    """The page is one very long Python string: a missing ``+`` between two literals
    parses as Python and only fails in the browser, where it takes the whole UI down.
    Parse it here instead."""
    data = {"file": "x.h5", "peaks": [], "ranges": [], "meta": {}}
    pages = [
        viz.render_html(data),
        viz.render_html(data, mode="app"),
        app._START_HTML,
    ]
    for page in pages:
        scripts = re.findall(r"<script>(.*?)</script>", page, re.S)
        assert scripts, "a page with no script did not render"
        for script in scripts:
            fd, name = tempfile.mkstemp(suffix=".mjs")
            with os.fdopen(fd, "w", encoding="utf-8") as handle:
                handle.write(script)
            try:
                proc = subprocess.run(
                    ["node", "--check", name], capture_output=True, timeout=60
                )
                assert proc.returncode == 0, proc.stderr.decode("utf-8", "replace")[:400]
            finally:
                os.unlink(name)


def test_render_html_mode_flag_is_not_confused_with_the_review_page():
    data = {"file": "x.h5", "peaks": [], "ranges": [], "meta": {}}
    app_page = viz.render_html(data, mode="app", auto_tour=False)
    assert "const APPMODE = true" in app_page
    assert "const AUTO_TOUR = false" in app_page
    assert "const APPMODE = false" in viz.render_html(data)


def test_tour_migrates_browser_state_and_tolerates_storage_errors():
    js = _review_js()
    assert "let _onboarded=false" in js
    assert "if(_onboarded){ rememberTour();" in js
    assert 'fetch("/onboarding",{method:"POST"})' in js
    assert "setTimeout(startAutoTour,650)" in js
    assert "checkpanel" not in js
    assert "ptrms-cl" not in js


def test_review_flushes_the_latest_edit_when_the_page_closes():
    data = {"file": "x.h5", "peaks": [], "ranges": [], "meta": {}}
    html = viz.render_html(data, config_path="/tmp/x.json", mode="app")
    assert 'window.addEventListener("pagehide",flushSave)' in html
    assert 'window.addEventListener("beforeunload",flushSave)' in html
    assert "navigator.sendBeacon(saveUrl(\"/save\",version,true),blob)" in html
    assert "if(!r.ok) throw new Error(\"save rejected\")" in html
    assert 'if(version===pendingSaveVersion) setStat("save failed")' in html


def test_no_template_marker_survives_rendering():
    """A marker left in the page parses as an identifier and only fails in the
    browser, where it kills the whole script: catch it here instead."""
    data = {"file": "x.h5", "peaks": [], "ranges": [], "meta": {}}
    for mode in ("review", "app"):
        html = viz.render_html(data, config_path="/tmp/x.json", mode=mode)
        assert "/*__" not in html, f"unreplaced template marker in {mode} mode"


HOOK = Path(__file__).resolve().parents[1] / "packaging" / "runtime_hook.py"


def _run_hook(argv, frozen):
    """Exec the hook the way PyInstaller does, and report what argv became."""
    saved = sys.argv
    try:
        sys.argv = list(argv)
        if frozen:
            sys.frozen = True
        exec(compile(HOOK.read_text(), str(HOOK), "exec"), {"__name__": "runtime_hook"})
        return list(sys.argv)
    finally:
        sys.argv = saved
        if frozen:
            del sys.frozen


def test_a_bare_double_click_becomes_app_mode():
    # Finder starts Contents/MacOS/sniff with no arguments at all; the plain CLI would
    # answer that with usage text and exit 2, which in a windowed bundle is invisible.
    assert _run_hook(["sniff"], frozen=True) == ["sniff", "app"]


def test_the_hook_leaves_named_commands_and_a_plain_cli_alone():
    assert _run_hook(["sniff", "app", "--no-browser"], frozen=True) == ["sniff", "app", "--no-browser"]
    assert _run_hook(["sniff"], frozen=False) == ["sniff"]


def test_the_spec_still_installs_the_hook():
    """A spec that stopped listing the hook would break double-click silently."""
    spec = HOOK.with_name("sniff-app.spec").read_text()
    assert "runtime_hook.py" in spec


def test_the_windows_bundle_separates_desktop_and_cli_launchers():
    """The desktop stays quiet without swallowing terminal command output."""
    spec = HOOK.with_name("sniff-app.spec").read_text()
    assert 'make_executable(name="sniff", console=False)' in spec
    assert 'make_executable(name="sniff-cli", console=True)' in spec
    assert "console=not IS_MAC" not in spec


def test_the_recent_api_flags_the_open_file(server, tmp_path, monkeypatch):
    """The backend keeps the open-file flag for API clients, independently of the UI."""
    api, session = server
    # This test is the only writer of its own recents file; the session-wide one would
    # otherwise see these paths too.
    monkeypatch.setattr(app, "RECENT_PATH", tmp_path / "recent.json")
    first = tmp_path / "run.h5"
    other = tmp_path / "other.h5"
    make_h5(first)
    make_h5(other)
    app.remember_recent(other)
    app.remember_recent(first)

    with (
        mock.patch.object(app, "auto_peaks", return_value=[]),
        mock.patch.object(app, "auto_ranges", return_value=[]),
        mock.patch.object(app.viz, "build_viz_data", payload_stub),
    ):
        api.post("/open", {"path": str(first)})
        _wait_ready(api)

    _, recent = api.get("/api/recent")
    entries = json.loads(recent)
    # Keyed by full path: the session-wide recents file also holds run.h5 from other
    # tests, and a name lookup would land on one of those.
    by_path = {e["path"]: e for e in entries}
    assert by_path[str(first.resolve())]["is_open"] is True
    assert by_path[str(other.resolve())]["is_open"] is False

    session.close()
    _, recent = api.get("/api/recent")
    assert all(not e.get("is_open") for e in json.loads(recent))


def test_browse_says_so_when_the_system_has_no_file_dialog(server, monkeypatch):
    """A missing native dialog is reported instead of offering manual path entry."""
    api, _ = server
    monkeypatch.setattr(app.sys, "platform", "linux")
    monkeypatch.setattr(app.os, "name", "posix")
    monkeypatch.setattr(app.shutil, "which", lambda _cmd: None)
    code, body = api.post("/browse", {})
    assert code == 501
    assert body["error"] == "Native file browsing is unavailable."


def test_kde_browse_uses_kdialogs_open_file_command(monkeypatch):
    """A KDE-only installation must still be able to choose an H5 file."""
    monkeypatch.setattr(app.sys, "platform", "linux")
    monkeypatch.setattr(app.os, "name", "posix")
    monkeypatch.setattr(
        app.shutil,
        "which",
        lambda command: "/usr/bin/kdialog" if command == "kdialog" else None,
    )
    run = mock.Mock(return_value=mock.Mock(stdout="/data/run.h5\n"))
    monkeypatch.setattr(app.subprocess, "run", run)

    assert app._pick_file() == "/data/run.h5"
    run.assert_called_once_with(
        [
            "kdialog",
            "--getopenfilename",
            app.os.path.expanduser("~"),
            "*.h5|IONICON runs (*.h5)",
        ],
        capture_output=True,
        text=True,
        timeout=300,
    )


def test_the_start_screen_uses_browse_for_files_and_compound_input_only():
    """The intro must not expose a manual file-path control."""
    html = app._START_HTML
    js = re.findall(r"<script>(.*?)</script>", html, re.DOTALL)[0]

    assert '<input id="compound-input"' in html
    assert 'id="path"' not in html
    assert 'id="go"' not in html
    assert "Open run" not in html
    assert "$('#path')" not in js
    assert "$('#go')" not in js
    assert "fetch('/browse'" in js
    assert "showInterest(body.path)" in js


def test_the_start_screen_contains_the_compound_prior_modal():
    html = app._START_HTML
    js = re.findall(r"<script>(.*?)</script>", html, re.DOTALL)[0]

    assert 'id="interest" role="dialog" aria-modal="true"' in html
    assert 'id="compound-input"' in html
    assert 'role="combobox"' in html and 'role="listbox"' in html
    assert 'id="interest-skip"' in html and 'id="interest-continue"' in html
    assert "splitCompoundBlock" in js and "addCompoundBlock" in js
    assert "Not recognised by the bundled PTR Library." in js
    assert "body.compounds_of_interest=compounds.map(c=>c.name)" in js
    assert "1,2-butadiene" in js, "library isomers were not exposed to autocomplete"


def test_the_start_screen_shows_browse_unavailable_error():
    """The browse failure is concise and does not point users at a removed field."""
    js = re.findall(r"<script>(.*?)</script>", app._START_HTML, re.DOTALL)[0]

    assert "Native file browsing is unavailable." in js
    assert "type the path" not in js
    assert "$('#path').focus()" not in js


# --------------------------------------------------------------------------
# an open the reviewer can walk out of
# --------------------------------------------------------------------------
def stalled_payload(f, peaks, ranges, *, progress=None, should_stop=None, **kwargs):
    """A build that runs until it is asked to stop, reporting a fraction meanwhile.

    Cancellation is only real if the work in flight honours ``should_stop``, so the
    stub has to be the thing that raises, exactly like the streaming pass does.
    """
    for _ in range(300):
        if should_stop is not None and should_stop():
            raise app.ptrms.AnalysisCancelled("the analysis was cancelled")
        if progress is not None:
            progress(0.5)
        threading.Event().wait(0.01)
    return payload_stub(f, peaks, ranges)


def _wait_state(api, done, timeout=20):
    spent = 0.0
    while spent < timeout:
        _, state = api.get("/api/state")
        state = json.loads(state)
        if done(state):
            return state
        threading.Event().wait(0.02)
        spent += 0.02
    raise AssertionError("the app never reached the state the test was waiting for")


def _wait_status(api, wanted, timeout=20):
    return _wait_state(api, lambda s: s["status"] in wanted, timeout)


def test_cancelling_an_open_returns_the_session_to_empty(server, tmp_path):
    """A cancel is a user changing their mind, not a failure: the session comes back
    to exactly the state it was in before the open, with nothing to apologise for."""
    api, session = server
    h5 = tmp_path / "run.h5"
    make_h5(h5)
    with (
        mock.patch.object(app, "auto_peaks", return_value=[]),
        mock.patch.object(app, "auto_ranges", return_value=[]),
        mock.patch.object(app.viz, "build_viz_data", stalled_payload),
    ):
        assert api.post("/open", {"path": str(h5)})[0] == 202
        opening = _wait_status(api, ("loading",))
        assert opening["cancellable"] is True
        assert api.post("/cancel", {})[0] == 200
        after = _wait_status(api, ("empty", "error"))
    assert after["status"] == "empty", after
    assert after["stage"] == "" and after["error"] is None
    assert session.busy is False


def test_a_cancelled_session_opens_the_file_afterwards(server, tmp_path):
    api, _ = server
    h5 = tmp_path / "run.h5"
    make_h5(h5)
    with (
        mock.patch.object(app, "auto_peaks", return_value=[]),
        mock.patch.object(app, "auto_ranges", return_value=[]),
        mock.patch.object(app.viz, "build_viz_data", stalled_payload),
    ):
        api.post("/open", {"path": str(h5)})
        _wait_status(api, ("loading",))
        api.post("/cancel", {})
        _wait_status(api, ("empty",))
    with (
        mock.patch.object(app, "auto_peaks", return_value=[]),
        mock.patch.object(app, "auto_ranges", return_value=[]),
        mock.patch.object(app.viz, "build_viz_data", payload_stub),
    ):
        assert api.post("/open", {"path": str(h5)})[0] == 202
        assert _wait_status(api, ("ready", "error"))["status"] == "ready"


def test_cancelling_nothing_is_opening_answers_ok_anyway(server):
    """The Cancel button and the poll can cross, so the answer cannot depend on
    which of them arrived first."""
    api, session = server
    code, body = api.post("/cancel", {})
    assert code == 200 and body == {"ok": True}
    assert session.status == "empty"


def test_an_open_reports_progress_and_says_whether_it_can_be_left(server, tmp_path):
    api, _ = server
    h5 = tmp_path / "run.h5"
    make_h5(h5)
    _, before = api.get("/api/state")
    before = json.loads(before)
    assert before["progress"] is None and before["cancellable"] is False

    with (
        mock.patch.object(app, "auto_peaks", return_value=[]),
        mock.patch.object(app, "auto_ranges", return_value=[]),
        mock.patch.object(app.viz, "build_viz_data", stalled_payload),
    ):
        api.post("/open", {"path": str(h5)})
        opening = _wait_state(api, lambda s: (s["progress"] or 0) > 0.4)
        assert opening["status"] == "loading"
        assert opening["progress"] == pytest.approx(0.5)
        assert opening["cancellable"] is True
        api.post("/cancel", {})
        _wait_status(api, ("empty",))
    _, after = api.get("/api/state")
    after = json.loads(after)
    # Nothing outside an open has a fraction to report, and a float the page cannot
    # put on a bar is worse than the null it would have replaced.
    assert after["progress"] is None and after["cancellable"] is False


def test_progress_never_runs_backwards_across_the_phases(server, tmp_path):
    """One phase reporting a smaller fraction than the one before it would make the
    bar jump backwards on screen, which no amount of smoothing explains away."""
    api, _ = server
    h5 = tmp_path / "run.h5"
    make_h5(h5)
    look, go = threading.Event(), threading.Event()

    def out_of_order(f, peaks, ranges, *, progress=None, should_stop=None, **kwargs):
        for frac in (0.9, 0.2):  # the second believes it is further back than one
            if progress is not None:
                progress(frac)
            look.set()
            go.wait(10)
            go.clear()
            look.clear()
        return payload_stub(f, peaks, ranges)

    with (
        mock.patch.object(app, "auto_peaks", return_value=[]),
        mock.patch.object(app, "auto_ranges", return_value=[]),
        mock.patch.object(app.viz, "build_viz_data", out_of_order),
    ):
        api.post("/open", {"path": str(h5)})
        seen = []
        for _ in (0, 1):
            assert look.wait(10), "the build never reported a fraction"
            seen.append(json.loads(api.get("/api/state")[1])["progress"])
            go.set()
        status = _wait_status(api, ("ready", "error"))
    assert seen[0] == pytest.approx(0.9)
    assert all(v is None or v >= seen[0] for v in seen), "the bar moved backwards"
    assert status["status"] == "ready"


def test_a_cancel_during_detection_writes_no_half_made_config(tmp_path):
    """Detection is the one phase whose output is a file. Leaving a config behind
    would tell the next open that a human had already been here."""
    h5 = tmp_path / "run.h5"
    make_h5(h5)
    session = app.Session()

    def detect_and_cancel(_f):
        session._cancel.set()
        return []

    with (
        mock.patch.object(app, "auto_peaks", return_value=[]),
        mock.patch.object(app, "auto_ranges", side_effect=detect_and_cancel),
    ):
        assert session.open(str(h5)) is None
    assert not (tmp_path / "run.json").exists()
    assert session.status == "empty" and session.stage == "" and session.error is None


def test_the_deterministic_pipeline_says_where_detection_is(tmp_path):
    """auto_peaks is 0.1 s of the roughly 1 s band and auto_ranges the other 0.8 s,
    so the bar has to move at the boundary between them rather than at the end."""
    h5 = tmp_path / "run.h5"
    make_h5(h5)
    seen, order = [], []
    with (
        mock.patch.object(
            app, "auto_peaks", side_effect=lambda _f: order.append("peaks") or []
        ),
        mock.patch.object(
            app, "auto_ranges", side_effect=lambda _f: order.append("ranges") or []
        ),
    ):
        app.bootstrap_config(str(h5), progress=seen.append)
    assert order == ["peaks", "ranges"]
    assert seen == [0.0, pytest.approx(0.1), 1.0]


def test_the_start_screen_holds_an_open_behind_a_modal_sheet():
    """An open is the one thing on this page that takes long enough to be worth
    leaving, so it gets a real bar and a Cancel button instead of a line of text."""
    html = app._START_HTML
    assert 'role="dialog"' in html and 'aria-modal="true"' in html
    assert "#ov{position:fixed;inset:0" in html, "the sheet is not full screen"
    assert 'role="progressbar"' in html, "the bar is not announced as a bar"
    assert 'id="cancel"' in html and "fetch('/cancel'" in html
    assert "s.progress" in html and "s.cancellable" in html
    assert ".ovstage{margin:10px 0 12px" in html
    assert 'points="8,105 74,105 99,105 108,34 117,105 251,105"' in html
    for retired in ("nose", "mascot", "lab-coat", "magnifier"):
        assert retired not in html.lower()
    assert 'aria-hidden="true"' in html and 'class="ovart"' in html
    assert "iondrift" in html and "breathe" in html
    assert 'left(rest)' in html, "no ETA is derived from the fraction"
    assert "minutes left" in html and "s left" in html, "the ETA is a raw number"


def test_a_finished_open_navigates_and_offers_no_button_about_it():
    """The sheet must not stay standing on a page the user has already left, and
    Back must not return to it."""
    js = re.findall(r"<script>(.*?)</script>", app._START_HTML, re.S)[0]
    ready = js[js.index("}else if(watching){") :]
    ready = ready[: ready.index("if(s.status==='error')")]
    assert "handoff()" in ready, "a finished open does not use the completion handoff"
    assert "return;" in ready, "the poll keeps running into a page that is going away"
    assert "function handoff()" in js
    assert "animationend" in js and "setTimeout(settle,760)" in js
    assert "sessionStorage.setItem('sniff-review-entrance','1')" in js
    assert "prefers-reduced-motion: reduce" in js
    assert "location.replace('/review')" in js
    assert "Open the review" not in ready


def test_the_start_screen_no_longer_stops_the_app_from_its_main_panel():
    """A browser tab has no window to close, so the quit route survives as one quiet
    footer link on that surface only; it is no longer a button on the start screen."""
    html = app._START_HTML
    assert '>Stop the app</button>' not in html
    assert '<a class="link" id="quit" href="#" hidden>Stop the app</a>' in html
    assert "s.surface!=='browser'" in html, "the link is not tied to the surface"
    assert "'/shutdown'" in html, "the route lost its only page caller"


def test_new_configs_enable_smart_models_and_adapt_explicit_tables(tmp_path):
    h5 = tmp_path / "adapted.h5"
    make_h5(h5)
    template = [{"mz": 59.0, "label": "acetone", "formula": "C3H6O"}]
    with (
        mock.patch.object(app, "auto_peaks", return_value=[{"mz": 59.004}, {"mz": 73.0}]),
        mock.patch.object(app, "auto_ranges", return_value=[]),
    ):
        config = app.bootstrap_config(str(h5), template_peaks=template)

    assert config["analysis_schema_version"] == 3
    assert config["analyze"]["peak_fit"] == "empirical-v1"
    assert config["analyze"]["isotope_mode"] == "formula-envelope-v2"
    assert config["peaks"][0]["formula"] == "C3H6O"
    assert [peak["mz"] for peak in config["peaks"]] == [59.004, 73.0]
    adaptation = config["diagnostics"]["peak_table_adaptation"]
    assert adaptation["n_matched"] == 1
    assert adaptation["n_new"] == 1


def test_adaptation_refuses_to_replace_an_existing_target_review(tmp_path):
    h5 = tmp_path / "reviewed.h5"
    make_h5(h5)
    h5.with_suffix(".json").write_text(
        json.dumps({"peaks": [{"mz": 59.0}], "ranges": []}),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="already has a saved review"):
        app.Session().open(str(h5), template_peaks=[{"mz": 59.0}])


def test_review_offers_explicit_table_adaptation_and_peak_preview():
    page = viz.render_html(
        {
            "meta": {},
            "config_base": {},
            "per_cycle": {},
            "spectrum": [],
            "peaks": [],
            "ranges": [],
            "transmission": {"masses": [], "factors": []},
            "rate_constants": [],
        },
        mode="app",
    )

    assert "Use table on another file" in page
    assert "adapt_peaks:cfg.peaks" in page
    assert "'/peak-preview?lo='" in page
    assert '"/nist-webbook?mz="' not in page
    assert "Bundled catalogue proposals" in page
    assert "names/structures, not PTR-MS proof" in page
