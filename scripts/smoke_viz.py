#!/usr/bin/env python3
"""Browser regression test for the generated PTR-MS review page (Sniff).

This deliberately uses the installed ``agent-browser`` CLI rather than a browser
framework.  The page is generated from deterministic synthetic data, served from
an in-process HTTP server, and exercised in headless Chrome through real DOM
interactions.  Keep this test small: it protects the identification contract
without making the package grow a frontend test dependency.
"""

from __future__ import annotations

import copy
import http.server
import json
import shutil
import socketserver
import subprocess
import sys
import tempfile
import threading
from pathlib import Path
from typing import Any, ClassVar
from urllib.parse import urlparse

from sniff import app, viz

SESSION = "sniff-ms-viz-regression"


def _candidate(formula: str, name: str, probability: float) -> dict[str, Any]:
    """Return the smallest candidate record accepted by the review UI."""
    return {
        "formula": formula,
        "name": name,
        "ion_mz": 100.0,
        "delta_mDa": 0.0,
        "dbe": 1.0,
        "k": None,
        "k_estimated": True,
        "flags": [],
        "iso_pred": [0.01, 0.01],
        "iso_obs": None,
        "iso_used": False,
        "probability": probability,
    }


def _synthetic_data() -> dict[str, Any]:
    """Build a stable, standalone review payload without an HDF5 fixture."""
    return {
        "meta": {
            "file": r"C:\measurements\synthetic-review.h5",
            "ncyc": 4,
            "dur": 1.0,
            "x_axis_unit": "cycle",
            "x_axis": {
                "cycle": [1, 2, 3, 4],
                "relative": [0.0, 2.0, 9.0, 15.0],
                "absolute": [
                    1700000000.0004,
                    1700000000.0008,
                    1700000000.0012,
                    1700000000.0016,
                ],
                "absolute_available": True,
            },
            "a": 1000.0,
            "b": 0.0,
            "R": 1500.0,
            "R_phys": 3100.0,
            "primary_mz": 19.022,
            "proton": 1.0073,
            "k_anchor": 1.7,
            "kinetic": True,
            "humidity_correct": True,
            "humidity_p": 0.6,
            "humidity_ref": 1.3,
            "whole_run_windows": False,
            "sources": {
                "R": "config.analyze",
                "R_phys": "config.analyze",
                "primary_mz": "config.analyze",
                "whole_run_windows": "config.analyze",
                "K": "config.analyze",
                "molar_volume": "config.analyze",
                "kinetic": "config.analyze",
                "k_anchor": "config.analyze",
                "humidity_correct": "config.analyze",
                "humidity_p": "config.analyze",
                "humidity_ref": "config.analyze",
            },
            "K_source": "config.analyze",
            "K_file": 0.8,
            "K_file_source": "file acquisition calibration",
            "K_default": 1.0,
            "molar_volume": 24.5,
            "molar_volume_source": "config.analyze",
            "molar_volume_file": 24.0,
            "molar_volume_file_source": "file drift temperature",
            "primary_available": True,
            "humidity_ref_default": 1.0,
            "concentration_available": True,
        },
        "transmission": {"masses": [50.0, 200.0], "factors": [1.0, 1.0]},
        "per_cycle": {
            # not flat: the concentration conversion divides cycle by cycle, so a
            # drifting primary current separates mean(K/I_p) from K/mean(I_p)
            "primary": [100.0, 104.0, 98.0, 102.0],
            "humidity": [1.0, 1.0, 1.0, 1.0],
            "discriminator": [1.0, 1.0, 1.0, 1.0],
        },
        # The values are only for painting the deterministic canvas.  The review
        # assertions below concern the identification card and its interactions.
        "spectrum": [10] * 13000,
        "peaks": [
            {
                "_config_original": {
                    "mz": 100.0,
                    "label": "Unassigned sole candidate",
                    "apex": 99.999,
                    "center": 100.001,
                    "mass": 100.002,
                    "win_l": 0.04,
                    "win_r": 0.04,
                    "future": {"keep": [1, 2, 3]},
                    "ion_role": {"future_role": {"keep": True}},
                },
                "id": 0,
                "mz": 100.0,
                "apex": 100.0,
                "label": "Unassigned sole candidate",
                "formula": "",
                "k": None,
                "k_estimated": False,
                "flags": [],
                "ion_role": {
                    "model": "ion-role-v1",
                    "kind": "isotope",
                    "status": "supported",
                    "label": "possible isotope",
                    "evidence": ["test evidence"],
                    "exclude_from_analyte_assignment": True,
                },
                "clustered": False,
                "win_l": 0.04,
                "win_r": 0.04,
                "win_manual": True,
                "trace": [10.0, 10.0, 20.0, 20.0],
                "abundance": 15.0,
                "candidates": [
                    _candidate("C2H6O", "Ethanol candidate", 1.0),
                ],
                "id_confidence": 1.0,
                "id_ambiguous": False,
                "overlap": None,
            },
            {
                "id": 1,
                "mz": 110.0,
                "apex": 110.0,
                "label": "Ambiguous mix",
                "formula": "",
                "k": None,
                "k_estimated": False,
                "flags": [],
                "clustered": False,
                "win_l": 0.04,
                "win_r": 0.04,
                "win_manual": True,
                "trace": [20.0, 20.0, 20.0, 20.0],
                "abundance": 20.0,
                "candidates": [
                    _candidate("C4H10O", "Candidate one", 0.72),
                    _candidate("C5H12", "Candidate two", 0.28),
                ],
                "id_confidence": 0.72,
                "id_ambiguous": True,
                "overlap": None,
            },
            {
                "id": 2,
                "mz": 120.0,
                "apex": 120.0,
                "label": "Curated solvent",
                "formula": "C3H8O",
                "k": None,
                "k_estimated": False,
                "flags": [],
                "clustered": False,
                "win_l": 0.04,
                "win_r": 0.04,
                "win_manual": True,
                "trace": [30.0, 30.0, 30.0, 30.0],
                "abundance": 30.0,
                "candidates": [
                    _candidate("C3H8O", "Library solvent", 0.65),
                    _candidate("C4H10", "Other generated candidate", 0.35),
                ],
                "id_confidence": 0.65,
                "id_ambiguous": True,
                "overlap": None,
            },
            {
                "id": 3,
                "mz": 130.0,
                "apex": 130.0,
                "label": "Cluster component A",
                "formula": "",
                "k": None,
                "k_estimated": False,
                "flags": [],
                "clustered": True,
                "win_l": 130.0 / (2.0 * 1200.0),
                "win_r": 130.0 / (2.0 * 1200.0),
                "win_manual": False,
                "trace": [40.0, 40.0, 40.0, 40.0],
                "abundance": 40.0,
                "candidates": [],
                "id_confidence": None,
                "id_ambiguous": False,
                "overlap": {
                    "neighbor": 130.05,
                    "sep_mDa": 50.0,
                    "level": "deconvolved",
                },
            },
            {
                "id": 4,
                "mz": 130.05,
                "apex": 130.05,
                "label": "Cluster component B",
                "formula": "",
                "k": None,
                "k_estimated": False,
                "flags": [],
                "clustered": True,
                "win_l": 130.05 / (2.0 * 1200.0),
                "win_r": 130.05 / (2.0 * 1200.0),
                "win_manual": False,
                "trace": [35.0, 35.0, 35.0, 35.0],
                "abundance": 35.0,
                "candidates": [],
                "id_confidence": None,
                "id_ambiguous": False,
                "overlap": {"neighbor": 130.0, "sep_mDa": 50.0, "level": "deconvolved"},
            },
            {
                "id": 5,
                "mz": 140.0,
                "apex": 140.0,
                "label": "Isolated control",
                "formula": "",
                "k": None,
                "k_estimated": False,
                "flags": [],
                "clustered": False,
                "win_l": 140.0 / (2.0 * 1200.0),
                "win_r": 140.0 / (2.0 * 1200.0),
                "win_manual": False,
                "trace": [25.0, 25.0, 25.0, 25.0],
                "abundance": 25.0,
                "candidates": [],
                "id_confidence": None,
                "id_ambiguous": False,
                "overlap": None,
            },
        ],
        "ranges": [
            {
                "_config_original": {
                    "label": "sample_01",
                    "start": 1,
                    "end": 2,
                    "unit": "cycle",
                    "merged_gaps": [{"start": 2, "end": 2}],
                    "future": {"keep": True},
                },
                "label": "sample_01",
                "start": 1,
                "end": 2,
                "class": "sample",
            },
            {"label": "sample_02", "start": 3, "end": 4, "class": "sample"},
        ],
        "merge_note": "joined 2 gaps, 1 fell to baseline",
        "config_base": {
            "unknown_top_level": {"keep": True},
            "viz": {"unknown_setting": "keep"},
            "analyze": {"unknown_setting": "keep"},
            "checklist": ["retired top-level note"],
            "review": {"checklist": ["retired nested note"]},
        },
        # two library entries far from the fixture peaks: enough for name/formula
        # consistency and hand-drawn naming, without moving any existing mDa column
        "rate_constants": [
            {
                "name": "toluene",
                "formula": "C7H8",
                "mz": 93.0699,
                "k": 2.2,
                "k_estimated": False,
                "flags": [],
            },
            {
                "name": "acetone",
                "formula": "C3H6O",
                "mz": 59.0491,
                "k": 3.9,
                "k_estimated": False,
                "flags": [],
            },
        ],
    }


class _ReviewHandler(http.server.BaseHTTPRequestHandler):
    """Serve the generated page and deterministic whole/interval spectra."""

    html = ""
    spectrum = b"[]"
    interval_spectrum = b"[]"
    posts: ClassVar[list[tuple[str, dict[str, Any]]]] = []
    reject_next_done = False

    def log_message(self, *_args: Any) -> None:
        pass

    def do_POST(self) -> None:
        length = int(self.headers.get("Content-Length", "0"))
        body = json.loads(self.rfile.read(length) or b"{}")
        path = urlparse(self.path).path
        self.posts.append((path, body))
        rejected = path == "/done" and self.reject_next_done
        if rejected:
            type(self).reject_next_done = False
        payload = b'{"error":"stale"}' if rejected else b'{"ok":true}'
        self.send_response(409 if rejected else 200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(payload)))
        self.end_headers()
        self.wfile.write(payload)

    def do_GET(self) -> None:
        if self.path in ("/", "/index.html"):
            body, content_type = self.html.encode("utf-8"), "text/html; charset=utf-8"
        elif self.path.startswith("/spectrum"):
            body, content_type = self.interval_spectrum, "application/json"
        else:
            self.send_error(404)
            return
        self.send_response(200)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)


def _browser(session: str, *args: str, stdin: str | None = None) -> str:
    """Run one agent-browser command and return its stdout."""
    completed = subprocess.run(
        ["agent-browser", "--allow-file-access", "--session", session, *args],
        input=stdin,
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode:
        detail = (completed.stderr or completed.stdout).strip()
        raise AssertionError(
            "agent-browser {} failed: {}".format(" ".join(args), detail)
        )
    return completed.stdout.strip()


def _eval_json(session: str, expression: str) -> Any:
    """Evaluate ``JSON.stringify(expression)`` in the page and return the value."""
    raw = _browser(
        session, "eval", "--stdin", stdin="JSON.stringify(" + expression + ")"
    )
    value: Any = json.loads(raw)
    if isinstance(value, str):
        value = json.loads(value)
    return value


def _eval(session: str, expression: str) -> dict[str, Any]:
    """Evaluate JSON.stringify(expression) in the actual browser page."""
    value: Any = _eval_json(session, expression)
    if not isinstance(value, dict):
        raise TypeError("browser expression did not return an object")
    return value


def _freeze_animations(session: str) -> None:
    """Disable CSS transitions/animations for deterministic layout measurements.

    The sidebar animates ``grid-template-columns`` (~0.28 s), and headless Chrome
    only advances a transition when it is handed rendering frames, so geometry read
    straight after a click can report a mid-flight - or frozen - column width. The
    layout assertions here are about the settled geometry the app asks for, so the
    cosmetics are switched off instead of raced.
    """
    _browser(
        session,
        "eval",
        "(() => { const s=document.createElement('style'); "
        "s.textContent='*,*::before,*::after{transition:none!important;"
        "animation:none!important}'; document.head.appendChild(s); "
        "void document.body.offsetHeight; })()",
    )


def _open(session: str, url: str) -> None:
    """Open a page; callers freeze animations once the page has settled."""
    _browser(session, "open", url)


def _assert(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def _assert_eq(actual: Any, expected: Any, message: str) -> None:
    """Assert equality, reporting the actual value so layout drift is diagnosable."""
    if actual != expected:
        raise AssertionError(
            "{}: got {}".format(message, json.dumps(actual, sort_keys=True))
        )


def _assert_complete_posts(expected: dict[str, Any], path: str, start: int) -> int:
    """Require every new request body to equal the browser's full config snapshot."""
    recent = [
        body
        for posted_path, body in _ReviewHandler.posts[start:]
        if posted_path == path
    ]
    _assert(recent, f"browser did not POST {path} after the edit")
    _assert(
        all(body == expected for body in recent),
        f"{path} body differs from the browser's complete buildConfig()",
    )
    return len(_ReviewHandler.posts)


def _assert_config_round_trip(config: dict[str, Any]) -> None:
    """Check fields whose loss would make a saved review non-reproducible."""
    _assert(config["unknown_top_level"]["keep"], "unknown top-level field was dropped")
    _assert(
        config["viz"]["x_axis_unit"] in {"cycle", "relative", "absolute"},
        "x-axis unit was not saved",
    )
    _assert(config["viz"]["unknown_setting"] == "keep", "unknown viz field was dropped")
    _assert(
        config["analyze"]["unknown_setting"] == "keep",
        "unknown analyze field was dropped",
    )
    _assert(
        [peak["mz"] for peak in config["peaks"]] == [100, 110, 120, 130, 130.05, 140],
        "peaks did not round-trip",
    )
    _assert(
        [(item["label"], item["start"], item["end"]) for item in config["ranges"]]
        == [("sample_01", 1, 2), ("sample_02", 3, 4)],
        "ranges did not round-trip",
    )
    _assert("checklist" not in config, "retired checklist field was saved")
    _assert("review" not in config, "empty retired review object was saved")


def _standalone_browser_pass(data: dict[str, Any]) -> None:
    """Exercise the offline app through a real ``file:`` browser URL."""
    session = f"{SESSION}-file-{threading.get_ident()}"
    with tempfile.TemporaryDirectory(prefix="sniff-ms-viz-") as directory:
        html_path = Path(directory) / "review.html"
        html_path.write_text(viz.render_html(data), encoding="utf-8")
        download_path = Path(directory) / "downloaded-config.json"
        try:
            _open(session, html_path.as_uri())
            _browser(session, "wait", "--load", "networkidle")
            _browser(session, "eval", "localStorage.removeItem('ptrms-onboarded')")
            _browser(session, "reload")
            _browser(session, "wait", "--load", "networkidle")
            _browser(session, "wait", "800")
            _assert(
                _eval_json(
                    session,
                    "document.querySelector('#tourblock')?.hidden===false",
                ),
                "first analysis did not open the guided tour",
            )
            _browser(session, "eval", "document.querySelector('#tskip').click()")
            _assert(
                _eval_json(
                    session, "localStorage.getItem('ptrms-onboarded')==='1'"
                ),
                "skipping the guided tour was not remembered",
            )
            _browser(session, "reload")
            _freeze_animations(session)  # the reload dropped the injected override
            _browser(session, "wait", "--load", "networkidle")
            _browser(session, "wait", "800")
            _assert(
                _eval_json(
                    session,
                    "!document.querySelector('#tourblock') || "
                    "document.querySelector('#tourblock').hidden",
                ),
                "the guided tour reopened after it was skipped",
            )
            _assert(
                _eval(
                    session,
                    "({display:document.querySelector('#xaxiswrap').style.display, tab:tab})",
                )
                == {"display": "", "tab": "trace"},
                "standalone x-axis selector is not visible on Signal over time",
            )
            _browser(
                session,
                "eval",
                "document.querySelector('#maintabs button[data-tab=spec]').click()",
            )
            _assert(
                _eval(
                    session,
                    "({display:document.querySelector('#xaxiswrap').style.display, tab:tab})",
                )
                == {"display": "none", "tab": "spec"},
                "standalone x-axis selector remains visible on Mass spectrum",
            )
            _browser(
                session,
                "eval",
                "document.querySelector('#maintabs button[data-tab=trace]').click()",
            )
            _browser(
                session,
                "eval",
                "document.querySelector('#xaxisunit').value='relative'; document.querySelector('#xaxisunit').dispatchEvent(new Event('change',{bubbles:true}))",
            )
            relative_axis = _eval(
                session,
                "({unit:xAxisUnit, text:document.querySelector('#rngtbl').innerText, config:buildConfig()})",
            )
            _assert(
                relative_axis["unit"] == "relative"
                and "0.0 s" in relative_axis["text"],
                "standalone relative-time labels did not update",
            )
            _browser(
                session,
                "eval",
                "document.querySelector('#xaxisunit').value='absolute'; document.querySelector('#xaxisunit').dispatchEvent(new Event('change',{bubbles:true}))",
            )
            absolute_axis = _eval(
                session,
                "({unit:xAxisUnit, text:document.querySelector('#rngtbl').innerText, config:buildConfig()})",
            )
            _assert(
                absolute_axis["unit"] == "absolute"
                and "UTC" not in absolute_axis["text"]
                and "22:13:20" in absolute_axis["text"],
                "standalone absolute interval labels did not update",
            )
            _browser(
                session,
                "eval",
                "document.querySelector('#xaxisunit').value='cycle'; document.querySelector('#xaxisunit').dispatchEvent(new Event('change',{bubbles:true}))",
            )
            _browser(session, "eval", "document.querySelector('#methodBtn').click()")
            initial = _eval(
                session,
                "({body:document.body.innerText, text:document.querySelector('#methodlive').innerText, "
                "panelText:document.querySelector('#methodpanel').innerText, "
                "buttons:Array.from(document.querySelectorAll('button')).map(b=>b.innerText), "
                "served:SERVED, protocol:location.protocol})",
            )
            _assert(
                not initial["served"] and initial["protocol"] == "file:",
                "standalone smoke did not open through file:",
            )
            _assert(
                "Download config" in initial["body"] and "Done" not in initial["body"],
                "standalone controls still use served-mode wording",
            )
            _assert(
                "Download config" in initial["text"] and "Done" not in initial["text"],
                "standalone Methods still use served-mode wording",
            )
            _assert(
                "Mass calibration & drift" in initial["panelText"],
                "Methods omits the detailed scientific sections",
            )
            _assert(
                "Download config" in initial["buttons"]
                and "Done" not in initial["buttons"],
                "standalone export control is wrong",
            )

            _browser(
                session,
                "eval",
                "(() => { const set=(id,v)=>{const e=document.querySelector('#'+id);"
                "e.value=v; e.dispatchEvent(new Event('change',{bubbles:true}));}; "
                "set('R','1700'); set('Rphys','3300'); set('primarymz','20.022'); "
                "set('K','2.5'); set('Vm','25.5'); set('kanchor','2.2'); "
                "set('hump','0.8'); set('href','1.7'); "
                "document.querySelector('#kinetic').click(); document.querySelector('#humid').click(); "
                "document.querySelector('#wholewindows').click(); })()",
            )
            _browser(session, "wait", "800")
            edited = _eval(
                session,
                "({text:document.querySelector('#methodlive').innerText, "
                "banner:document.querySelector('#stalebanner').innerText, "
                "config:buildConfig(), stale:staleSettings()})",
            )
            _assert(
                "R = 1700 (browser edit)" in edited["text"],
                "standalone Methods did not mark R as edited",
            )
            _assert(
                "K: 2.50 (browser edit)" in edited["text"]
                and "molar volume: 25.50 L/mol (browser edit)" in edited["text"],
                "standalone Methods did not mark calibration edits",
            )
            _assert(
                "Kinetic correction: off (browser edit)" in edited["text"],
                "standalone Methods did not mark kinetic edits",
            )
            _assert(
                "Humidity correction: off (browser edit)" in edited["text"],
                "standalone Methods did not mark humidity edits",
            )
            _assert(
                "k_anchor = 2.2" in edited["text"]
                and "(browser edit)" in edited["text"],
                "standalone Methods did not mark k_anchor as edited",
            )
            _assert(
                "p = 0.8 (browser edit)" in edited["text"],
                "standalone Methods did not mark humidity p as edited",
            )
            _assert(
                "reference = 1.7 (browser edit)" in edited["text"],
                "standalone Methods did not mark humidity reference as edited",
            )
            _assert(
                "export and re-extraction" in edited["text"]
                and "Download config" in edited["text"]
                and "Done" not in edited["text"],
                "standalone stale Methods wording is inaccurate",
            )
            _assert(
                "export and re-extraction" in edited["banner"]
                and "Download config" in edited["banner"]
                and "Done" not in edited["banner"],
                "standalone stale banner wording is inaccurate",
            )
            _assert_config_round_trip(edited["config"])

            _browser(
                session,
                "eval",
                "(() => { const set=(id,v)=>{const e=document.querySelector('#'+id);"
                "e.value=v; e.dispatchEvent(new Event('change',{bubbles:true}));}; "
                "set('R','1500'); set('Rphys','3100'); set('primarymz','19.022'); "
                "set('K','1.0'); set('Vm','24.5'); set('kanchor','1.7'); set('hump','0.6'); "
                "set('href','1.3'); document.querySelector('#kinetic').click(); "
                "document.querySelector('#humid').click(); document.querySelector('#wholewindows').click(); })()",
            )
            _browser(session, "wait", "800")
            reverted = _eval(
                session,
                "({text:document.querySelector('#methodlive').innerText, "
                "banner:document.querySelector('#stalebanner').innerText, stale:staleSettings()})",
            )
            _assert(
                reverted["stale"]
                == {"primary": False, "Rphys": False, "windows": False},
                "standalone revert left stale state behind",
            )
            _assert(
                not reverted["banner"] and "PREVIEW STALE" not in reverted["text"],
                "standalone revert left a stale warning behind",
            )
            _assert(
                "browser edit" not in reverted["text"],
                "standalone revert did not restore initial provenance",
            )

            _browser(
                session,
                "eval",
                "document.querySelector('#R').value='1750'; "
                "document.querySelector('#R').dispatchEvent(new Event('change',{bubbles:true}))",
            )
            _browser(session, "wait", "700")
            payload = _eval(session, "({config:buildConfig()})")["config"]
            _assert_config_round_trip(payload)
            # the panel overlays the top bar, so close it before clicking the export row
            _browser(session, "eval", "document.querySelector('#methodClose').click()")
            _browser(
                session,
                "eval",
                "window.__downloadSnapshot=null; "
                "document.querySelector('#exportrow button').addEventListener("
                "'click',()=>window.__downloadSnapshot=buildConfig(),{once:true})",
            )
            try:
                _browser(session, "download", "#exportrow button", str(download_path))
            except AssertionError as error:
                # Chromium's download event is unavailable for file: pages in some
                # agent-browser daemon versions. Verify the real button handler still
                # built the complete payload before accepting that limitation.
                detail = str(error).lower()
                _assert(
                    "resource temporarily unavailable" in detail
                    or "allow-file-access ignored" in detail,
                    f"standalone download failed unexpectedly: {error}",
                )
                _browser(
                    session,
                    "eval",
                    "document.querySelector('#exportrow button').click()",
                )
                _assert(
                    _eval(session, "({config:window.__downloadSnapshot})")["config"]
                    == payload,
                    "standalone download handler did not use buildConfig()",
                )
            else:
                _assert(
                    download_path.exists(),
                    "standalone Download config did not create a file",
                )
                _assert(
                    json.loads(download_path.read_text(encoding="utf-8")) == payload,
                    "downloaded config differs from buildConfig()",
                )
        finally:
            _browser(session, "close")


def _start_screen_browser_pass() -> None:
    """Check the current-file card at the narrow widths users commonly have."""
    session = f"{SESSION}-start-{threading.get_ident()}"
    with tempfile.TemporaryDirectory(prefix="sniff-ms-start-") as directory:
        html_path = Path(directory) / "start.html"
        html_path.write_text(app._START_HTML, encoding="utf-8")
        try:
            _open(session, html_path.as_uri())
            _browser(session, "wait", "--load", "networkidle")
            _freeze_animations(session)
            _browser(
                session,
                "eval",
                "showInterest('/runs/context.h5'); const i=document.querySelector('#compound-input'); "
                "i.value='acet'; i.dispatchEvent(new Event('input',{bubbles:true}));",
            )
            partial = _eval(
                session,
                "({open:!document.querySelector('#interest').hidden, "
                "invalid:document.querySelector('#compound-input').getAttribute('aria-invalid'), "
                "expanded:document.querySelector('#compound-input').getAttribute('aria-expanded'), "
                "active:document.querySelector('#compound-input').getAttribute('aria-activedescendant'), "
                "suggestions:Array.from(document.querySelectorAll('.suggestion')).map(e=>e.innerText), "
                "dialogScroll:(()=>{const e=document.querySelector('.priorcard'); "
                "return e.scrollHeight>e.clientHeight})(), suggestionScroll:(()=>{const e="
                "document.querySelector('#compound-suggestions'); return e.scrollHeight>e.clientHeight})()})",
            )
            _assert(
                partial["open"]
                and partial["invalid"] == "true"
                and partial["expanded"] == "true"
                and partial["active"] == "compound-option-0"
                and any("acetone" in name for name in partial["suggestions"])
                and not partial["dialogScroll"]
                and partial["suggestionScroll"],
                "compound autocomplete is missing or makes the dialog itself scroll",
            )
            _browser(
                session,
                "eval",
                "const input=document.querySelector('#compound-input'); input.value='not in the library'; "
                "input.dispatchEvent(new Event('input',{bubbles:true})); "
                "window.__invalidTabAllowed=input.dispatchEvent(new KeyboardEvent('keydown',"
                "{key:'Tab',bubbles:true,cancelable:true})); "
                "window.__shiftTabAllowed=input.dispatchEvent(new KeyboardEvent('keydown',"
                "{key:'Tab',shiftKey:true,bubbles:true,cancelable:true})); "
                "input.value='acetone'; input.dispatchEvent(new Event('input',{bubbles:true})); "
                "input.dispatchEvent(new KeyboardEvent('keydown',"
                "{key:'Enter',bubbles:true,cancelable:true})); "
                "const transfer=new DataTransfer(); "
                "transfer.setData('text/plain','1,2-butadiene, methanol\\nnot in the library'); "
                "input.dispatchEvent(new ClipboardEvent('paste',"
                "{clipboardData:transfer,bubbles:true,cancelable:true}));",
            )
            compounds = _eval(
                session,
                "({chips:Array.from(document.querySelectorAll('.chip>span')).map(e=>e.textContent), "
                "error:document.querySelector('#compound-status').innerText, "
                "enabled:!document.querySelector('#interest-continue').disabled, "
                "tabAllowed:window.__invalidTabAllowed,shiftTabAllowed:window.__shiftTabAllowed, "
                "chipRoles:Array.from(document.querySelectorAll('.chip')).map(e=>e.getAttribute('role'))})",
            )
            _assert_eq(
                compounds["chips"],
                ["acetone", "1,2-butadiene", "methanol"],
                "compound entry did not preserve autocomplete and comma-containing names",
            )
            _assert(
                "Rejected: not in the library" in compounds["error"]
                and compounds["enabled"]
                and compounds["tabAllowed"]
                and compounds["shiftTabAllowed"]
                and set(compounds["chipRoles"]) == {"listitem"},
                "bulk entry, invalid Tab handling, or chip semantics regressed",
            )
            _browser(
                session,
                "eval",
                "openFile=(path,compounds)=>{window.__opened={path,"
                "names:compounds&&compounds.map(c=>c.name),hasCompounds:compounds!==undefined}}; "
                "document.querySelector('#interest-continue').click();",
            )
            continued = _eval(session, "window.__opened")
            _assert_eq(
                continued,
                {
                    "path": "/runs/context.h5",
                    "names": ["acetone", "1,2-butadiene", "methanol"],
                    "hasCompounds": True,
                },
                "the compound modal did not submit its validated selection",
            )
            _browser(
                session,
                "eval",
                "showInterest('/runs/skipped.h5'); document.querySelector('#interest-skip').click();",
            )
            skipped = _eval(session, "window.__opened")
            _assert_eq(
                skipped,
                {"path": "/runs/skipped.h5", "names": [], "hasCompounds": True},
                "skipping the compound modal did not clear the contextual prior",
            )
            _browser(
                session,
                "eval",
                "current({file:'/runs/2026/September/very-long-ionicon-run-name.h5'})",
            )
            for width in (320, 375):
                _browser(session, "set", "viewport", str(width), "700")
                layout = _eval(
                    session,
                    "(() => { const card=document.querySelector('#now .now'); "
                    "const txt=card.querySelector('.txt'), name=txt.querySelector('b'), "
                    "sub=txt.querySelector('.sub'); "
                    "const buttons=[...card.querySelectorAll('button')]; "
                    "const cr=card.getBoundingClientRect(); "
                    "return {width:window.innerWidth, name:name.textContent, "
                    "sub:sub.textContent, nameHeight:name.getBoundingClientRect().height, "
                    "subHeight:sub.getBoundingClientRect().height, "
                    "nameFits:name.scrollWidth<=name.clientWidth, "
                    "subFits:sub.scrollWidth<=sub.clientWidth, "
                    "actionsVisible:buttons.every(e=>{const r=e.getBoundingClientRect(); "
                    "return r.width>0&&r.height>0&&r.left>=cr.left&&r.right<=cr.right}), "
                    "noPageOverflow:document.documentElement.scrollWidth<=window.innerWidth}; })()",
                )
                _assert_eq(
                    layout["width"], width, "start screen viewport did not apply"
                )
                _assert(
                    layout["name"] == "very-long-ionicon-run-name.h5"
                    and layout["sub"].endswith("/September"),
                    "current file identity is not rendered at narrow width",
                )
                _assert(
                    layout["nameHeight"] > 0
                    and layout["subHeight"] > 0
                    and layout["nameFits"]
                    and layout["subFits"]
                    and layout["actionsVisible"]
                    and layout["noPageOverflow"],
                    f"current-file card does not fit at {width}px: {layout}",
                )
        finally:
            _browser(session, "close")


def _review_round_browser_pass(session: str) -> None:
    """Regress interval edits, unit-aware values, sample ticks and naming.

    Every mutation made here is undone before the pass returns, so the served
    page keeps handing its original state to the Done checks that follow.
    """
    # the ticks read their scope from the 'average over' control, so this pass owns it
    # from here on and starts from the whole-run view like any fresh review
    _browser(session, "eval", "setSpecRange('all')")
    # --- interval edits: the card must show what the plot now shows, in time order ---
    _browser(
        session,
        "eval",
        "document.querySelector('#maintabs button[data-tab=trace]').click()",
    )
    _browser(
        session,
        "eval",
        "selRange=ranges.find(r=>r.label==='sample_01')._id; renderRanges(); drawMain();",
    )
    # a synthetic drag along the trace plot, dispatched on the canvas so the app's own
    # offset math and edge hit-testing run exactly as they do for a real drag
    drag_js = (
        "(() => { const rect=plotC.getBoundingClientRect(); "
        "const send=(t,x,tgt)=>tgt.dispatchEvent(new MouseEvent(t,{bubbles:true,cancelable:true,"
        "view:window,clientX:rect.left+x,clientY:rect.top+20})); "
        "const xAt=c=>traceX(plotC.clientWidth)(axisAtCycle(c)); "
        "const drag=(x0,x1)=>{ send('mousedown',x0,plotC); send('mousemove',(x0+x1)/2,plotC); "
        "send('mousemove',x1,plotC); send('mouseup',x1,window); }; return {drag,xAt}; })()"
    )
    _browser(session, "eval", "window.__drag=" + drag_js + ";")
    # pull sample_01's start past its own end (the app swaps them), then start the
    # interval later than sample_02 so the card rows have to reorder
    _browser(
        session,
        "eval",
        "window.__drag.drag(window.__drag.xAt(1), window.__drag.xAt(3.9)); "
        "window.__drag.drag(window.__drag.xAt(3.95), window.__drag.xAt(3.99));",
    )
    reorder = _eval(
        session,
        "({dom:Array.from(document.querySelectorAll('#rngtbl tbody tr')).map(tr=>"
        "[tr.querySelector('.lbl').value, tr.querySelector('td.mini').textContent]), "
        "order:ranges.map(r=>r.label), sorted:ranges.every((r,i,a)=>!i||a[i-1].start<=r.start), "
        "cells:ranges.map(r=>formatRange(r))})",
    )
    _assert(
        reorder["sorted"] and [row[0] for row in reorder["dom"]] == reorder["order"],
        "interval rows are not kept in chronological order after a resize: "
        + str(reorder["dom"]),
    )
    _assert(
        [row[1] for row in reorder["dom"]] == reorder["cells"],
        "the intervals table still shows stale ranges after resizing: "
        + str(reorder["dom"]),
    )
    _browser(
        session,
        "eval",
        "ranges.forEach(r=>{ if(r.label==='sample_01'){ r.start=1; r.end=2; } "
        "else { r.start=3; r.end=4; } }); sortRanges(); renderRanges(); redraw();",
    )
    restored = _eval(
        session,
        "({rows:Array.from(document.querySelectorAll('#rngtbl tbody tr td.mini'))"
        ".map(td=>td.textContent), expect:ranges.map(r=>formatRange(r))})",
    )
    _assert(
        restored["rows"] == restored["expect"],
        "interval restore did not return the table to its original ranges",
    )

    # --- the 'average over' names follow interval renames ---
    _browser(
        session,
        "eval",
        "(() => { const i=document.querySelector('#rngtbl tbody tr .lbl'); "
        "i.value='breath_01'; i.dispatchEvent(new Event('change',{bubbles:true})); })()",
    )
    renamed = _eval(
        session,
        "({options:Array.from(document.querySelectorAll('#scoperange option'))"
        ".map(o=>o.textContent), ticks:buildConfig().peaks.map(p=>p.samples||null)})",
    )
    _assert(
        "breath_01" in renamed["options"],
        "the interval list kept a stale interval name: " + str(renamed["options"]),
    )
    _assert(
        all(t is None for t in renamed["ticks"]) and len(renamed["ticks"]) == 6,
        "renaming an interval changed which samples include each compound: "
        + str(renamed["ticks"]),
    )
    _browser(
        session,
        "eval",
        "(() => { const i=document.querySelector('#rngtbl tbody tr .lbl'); "
        "i.value='sample_01'; i.dispatchEvent(new Event('change',{bubbles:true})); })()",
    )
    back = _eval(
        session,
        "({options:Array.from(document.querySelectorAll('#scoperange option'))"
        ".map(o=>o.textContent), config:buildConfig()})",
    )
    _assert(
        "sample_01" in back["options"],
        "interval rename did not round-trip into the interval list",
    )
    _assert_config_round_trip(back["config"])

    # --- an interval edit must never drop a curated compound from the analysis ---
    interval_edit = _eval(
        session,
        "(() => { const snap=ranges.map(r=>({...r})); const before=buildConfig().peaks.length; "
        "ranges.forEach(r=>{ r.class='background'; }); renderRanges(); redraw(); "
        "const allBackground=buildConfig().peaks.length; "
        "ranges.length=0; ranges.push({...snap[0]}); renderRanges(); redraw(); "
        "const oneDeleted=buildConfig().peaks.length; "
        "ranges.length=0; snap.forEach(r=>ranges.push({...r})); sortRanges(); "
        "renderRanges(); syncSpecRange(); redraw(); "
        "return {before, allBackground, oneDeleted, restored:buildConfig().peaks.length}; })()",
    )
    for key, label in (
        ("allBackground", "every interval became a background"),
        ("oneDeleted", "an interval was deleted"),
    ):
        _assert(
            interval_edit[key] == interval_edit["before"],
            "the curated peak list changed because {}: {} peaks before, {} after".format(
                label, interval_edit["before"], interval_edit[key]
            ),
        )
    _assert(
        interval_edit["restored"] == interval_edit["before"],
        "the interval edit round trip lost compounds",
    )

    # --- a new sample interval keeps 'every sample' compounds in every sample ---
    new_interval = _eval(
        session,
        "(() => { const all=peaks[0], part=peaks[1]; setSel(part,[sampleLabels()[0]]); "
        "const was={all:selState(all), part:selState(part)}; "
        "const wasAll=peaks.map(p=>selState(p)==='all'); "
        "ranges.push({label:'sample_03', class:'sample', start:3, end:4, _id:nextRangeId++}); "
        "sortRanges(); adoptSampleKey('sample_03', wasAll); renderRanges(); redraw(); "
        "const now={all:selState(all), part:selState(part)}; "
        "const i=ranges.findIndex(r=>r.label==='sample_03'); if(i>=0) ranges.splice(i,1); "
        "dropSampleKey('sample_03'); peaks.forEach(p=>selAll(p)); sortRanges(); "
        "renderRanges(); syncSpecRange(); renderPeaks(); redraw(); "
        "return {was, now, n:ranges.length, back:peaks.every(p=>selState(p)==='all')}; })()",
    )
    _assert(
        new_interval["was"] == {"all": "all", "part": "some"}
        and new_interval["now"] == {"all": "all", "part": "some"},
        "adding a sample interval changed which compounds were in every sample: "
        + str(new_interval),
    )
    _assert(
        new_interval["n"] == 2 and new_interval["back"],
        "the new-interval check did not restore the review state: " + str(new_interval),
    )

    # --- one sample interval at a time: the tick box speaks about the interval in
    # scope, so dropping a compound out of one sample leaves the others alone ---
    _browser(session, "eval", "document.querySelector('#pkdetails').click()")
    scoped = _eval(
        session,
        "(() => { const p=peaks.find(q=>q.label==='Curated solvent'); "
        "const si=sampleIntervals(); "
        "const row=()=>Array.from(document.querySelectorAll('#peaksbody li'))"
        ".find(e=>e.querySelector('.lbl').value==='Curated solvent'); "
        "const box=()=>row().querySelector('[data-a=use]'); "
        "const setScope=r=>{ const sel=document.getElementById('scoperange'); "
        "  sel.value=r==null?'all':'i'+r._id; "
        "  sel.dispatchEvent(new Event('change',{bubbles:true})); }; "
        "setScope(si[0]); "
        "const onTick=[box().checked,box().indeterminate]; "
        "box().click();                                    /* out of this sample only */ "
        "const written=buildConfig().peaks.find(q=>q.label==='Curated solvent'); "
        "const out={scoped:scopeRange()===si[0], onTick, "
        "  off:{state:selState(p), box:[box().checked,box().indeterminate], "
        "  samples:written.samples, aggregate:selState(p,null)}}; "
        "box().click(); out.back=selState(p,null); "
        "setScope(si[1]); out.other=box().checked; "
        "setScope(null); out.aggregate=selState(p); return out; })()",
    )
    _assert(
        scoped["scoped"] and scoped["onTick"] == [True, False],
        "with one sample in scope the tick box is not that sample's own tick: "
        + str(scoped),
    )
    _assert(
        scoped["off"]["samples"] == ["sample_02"]
        and scoped["off"]["state"] == "none"
        and scoped["off"]["aggregate"] == "some"
        and scoped["off"]["box"] == [False, False],
        "unticking one sample interval changed the other samples too: " + str(scoped),
    )
    _assert(
        scoped["back"] == "all"
        and scoped["other"] is True
        and scoped["aggregate"] == "all",
        "ticking the sample back did not restore the compound everywhere: "
        + str(scoped),
    )

    # ticking one sample and unticking another between renders must not cancel the
    # compound out of the analysis
    two_toggles = _eval(
        session,
        "(() => { const p=peaks.find(q=>q.label==='Curated solvent'); "
        "setSel(p,[sampleLabels()[0]]); renderPeaks(); "
        "const row=()=>Array.from(document.querySelectorAll('#peaksbody li'))"
        ".find(e=>e.querySelector('.lbl').value==='Curated solvent'); "
        "const box=()=>row().querySelector('[data-a=use]'); "
        "const setScope=r=>{ const sel=document.getElementById('scoperange'); "
        "  sel.value=r==null?'all':'i'+r._id; "
        "  sel.dispatchEvent(new Event('change',{bubbles:true})); }; "
        "setScope(sampleIntervals()[1]); box().click();    /* in with the second */ "
        "setScope(sampleIntervals()[0]); box().click();    /* out with the first */ "
        "const written=buildConfig().peaks.find(q=>q.label==='Curated solvent'); "
        "const out={state:selState(p,null), samples:written?written.samples:null, "
        "  listed:!!written}; "
        "selAll(p); setScope('all'); renderPeaks(); redraw(); return out; })()",
    )
    _assert(
        two_toggles["state"] == "some"
        and two_toggles["samples"] == ["sample_02"]
        and two_toggles["listed"],
        "ticking one sample and unticking another lost the compound: "
        + str(two_toggles),
    )

    # --- the interval in scope is one control, and it follows the interval the
    # review is looking at; reclassifying happens in the table, by name ---
    scope = _eval(
        session,
        "(() => { const p=peaks.find(q=>q.label==='Curated solvent'); "
        "p.samples=sampleLabels().slice(); renderPeaks();  /* in both, recorded by hand */ "
        "const before=JSON.stringify(buildConfig()); "
        "const sel=document.getElementById('scoperange'); "
        "const r=sampleIntervals()[sampleIntervals().length-1], first=r.label; "
        "sel.value='i'+r._id; sel.dispatchEvent(new Event('change')); "
        "const follows=SPECWIN.lo===r.start && SPECWIN.hi===r.end && scopeRange()===r; "
        "const row=()=>document.querySelector('#rngtbl tbody tr.sel'); "
        "const cls=()=>row().querySelector('select'); "
        "cls().value='background'; cls().dispatchEvent(new Event('change',{bubbles:true})); "
        "const off={cls:r.class, label:r.label, "
        "  saved:buildConfig().ranges.some(q=>q.label===r.label)}; "
        "cls().value='sample'; cls().dispatchEvent(new Event('change',{bubbles:true})); "
        "const w=buildConfig().peaks.find(q=>q.label==='Curated solvent'); "
        "return {follows, off, first, plotSelected:selRange===r._id, "
        "  clsBack:r.class, labelBack:r.label, "
        "  samples:(w && 'samples' in w) ? w.samples : null, "
        "  untouched:JSON.stringify(buildConfig())===before}; })()",
    )
    _assert(
        scope["follows"] and scope["plotSelected"],
        "the interval selector does not drive the spectrum and the plot together: "
        + str(scope),
    )
    # the class is carried by the interval name, which is what the analysis blanks on,
    # so a class switch that did not rename would be lost on save
    _assert(
        scope["off"]["cls"] == "background"
        and str(scope["off"]["label"]).startswith("background")
        and scope["off"]["saved"],
        "reclassifying did not rename the interval for the analysis: "
        + str(scope["off"]),
    )
    _assert(
        scope["clsBack"] == "sample"
        and scope["labelBack"] == scope["first"]
        and scope["untouched"]
        and scope["samples"] is None,
        "classing an interval background and back changed the config: " + str(scope),
    )

    # --- whole run: partial -> tick -> empty -> tick, never a dead end ---
    cycle = _eval(
        session,
        "(() => { const p=peaks.find(q=>q.label==='Curated solvent'); "
        "setSpecRange('all'); setSel(p,[sampleLabels()[0]]); renderPeaks(); "
        "const state=()=>{ const li=Array.from(document.querySelectorAll('#peaksbody li'))"
        ".find(e=>e.querySelector('.lbl').value==='Curated solvent'); "
        "const b=li.querySelector('[data-a=use]'); "
        "return [selState(p), b.checked, b.indeterminate]; }; "
        "const click=()=>{ Array.from(document.querySelectorAll('#peaksbody li'))"
        ".find(e=>e.querySelector('.lbl').value==='Curated solvent')"
        ".querySelector('[data-a=use]').click(); }; "
        "const was=state(); click(); const all=state(); click(); const none=state(); "
        "click(); const again=state(); setSel(p,sampleLabels()); renderPeaks(); redraw(); "
        "return {was, all, none, again}; })()",
    )
    _assert(
        cycle["was"] == ["some", True, True]
        and cycle["all"] == ["all", True, False]
        and cycle["none"] == ["none", False, False]
        and cycle["again"] == ["all", True, False],
        "the whole-run box does not flick between tick and empty: " + str(cycle),
    )

    # --- the faint composite VOC curve names itself and can be switched off ---
    disc = _eval(
        session,
        "(() => { const b=document.getElementById('legDisc'); "
        "const label=b.textContent.trim(); b.click(); "
        "const off={state:showDisc, cfg:buildConfig().viz.show_disc, "
        "aria:b.getAttribute('aria-pressed')}; "
        "b.click(); return {label, off, on:showDisc, "
        "cfgOn:buildConfig().viz.show_disc}; })()",
    )
    _assert(
        "own scale" in disc["label"],
        "the composite VOC legend does not say the curve is on its own scale: "
        + str(disc),
    )
    _assert(
        disc["off"]["state"] is False
        and disc["off"]["cfg"] is False
        and disc["off"]["aria"] == "false"
        and disc["on"] is True
        and disc["cfgOn"] is True,
        "the composite VOC legend is not a working switch: " + str(disc),
    )

    # --- the context card shrinks as far as it is asked to, on a tall window ---
    was_size = _eval(session, "({w:window.innerWidth, h:window.innerHeight})")
    _browser(session, "set", "viewport", "1280", "1400")
    split = _eval(
        session,
        "(() => { const keep=plotHeight; plotHeight=1e9; relayout(); "
        "const ctx=document.getElementById('intcard'); "
        "const head=ctx.querySelector('h2').offsetHeight; "
        "const out={plot:+plotC.dataset.h, body:ctx.offsetHeight-head, min:MIN_CONTEXT_HEIGHT, "
        "max:Number(plotResize.getAttribute('aria-valuemax'))}; "
        "plotHeight=keep; relayout(); return out; })()",
    )
    _browser(session, "set", "viewport", str(was_size["w"]), str(was_size["h"]))
    _assert(
        split["plot"] > 560
        and split["body"] <= split["min"] + 4
        and abs(split["max"] - split["plot"]) <= 1,
        "the plot stops growing before the card below is full: " + str(split),
    )

    # --- editing one label then clicking another must preserve the new focus ---
    original_ui = _eval(session, "({tab,selId})")
    for active_tab in ("trace", "spec"):
        focus_target = _eval(
            session,
            "(() => { const order=orderedPeaks(); "
            f"setTab({json.dumps(active_tab)}); selId=order[0].id; renderPeaks(); "
            "const source=document.querySelector('#peaksbody li.sel .lbl'); "
            "const target=document.querySelector('#peaksbody li:not(.sel) .lbl'); "
            "source.dataset.editSource='true'; target.dataset.focusTest='true'; "
            "return {sourceId:order[0].id, sourceLabel:order[0].label, "
            "targetId:order[1].id, targetLabel:order[1].label}; })()",
        )
        temporary_label = f"Temporary {active_tab} label"
        _browser(session, "fill", "[data-edit-source]", temporary_label)
        _browser(session, "click", "[data-focus-test]")
        label_focus = _eval(
            session,
            "(() => { const i=document.querySelector('[data-focus-test]'); "
            f"const source=peaks.find(p=>p.id==={focus_target['sourceId']}); "
            "return {exists:!!i, focused:document.activeElement===i, "
            "readonly:i?i.readOnly:null, selectedId:selId, sourceLabel:source.label, "
            "trace:document.getElementById('traceof').textContent}; })()",
        )
        _assert(
            label_focus["exists"]
            and label_focus["focused"]
            and label_focus["readonly"] is False
            and label_focus["selectedId"] == focus_target["targetId"]
            and label_focus["sourceLabel"] == temporary_label
            and label_focus["trace"].startswith(focus_target["targetLabel"]),
            f"editing then switching labels failed on {active_tab}: "
            + str({"target": focus_target, "state": label_focus}),
        )
        _browser(
            session,
            "eval",
            "document.activeElement.blur(); "
            f"peaks.find(p=>p.id==={focus_target['sourceId']}).label="
            f"{json.dumps(focus_target['sourceLabel'])}; renderPeaks(); drawMain(); scheduleSave()",
        )
    _browser(
        session,
        "eval",
        f"selId={json.dumps(original_ui['selId'])}; "
        f"setTab({json.dumps(original_ui['tab'])}); renderPeaks(); drawMain()",
    )

    # --- arrow keys walk the list the user is looking at, not the stored order ---
    nav = _eval(
        session,
        "(() => { const set=o => { peakOrder=o; "
        "document.getElementById('pkorder').value=o; renderPeaks(); "
        "return [...document.querySelectorAll('#peaksbody li .lbl')].map(i=>i.value); }; "
        "const walk=o => { set(o); selId=orderedPeaks()[0].id; const seen=[]; "
        "for(let i=0;i<3;i++){ seen.push((peaks.find(p=>p.id===selId)||{}).label); "
        "document.dispatchEvent(new KeyboardEvent('keydown',{key:'ArrowDown',bubbles:true})); } "
        "return seen; }; "
        "const lists={}, walks={}; "
        "for(const o of ['mz','abundance','label']){ lists[o]=set(o); walks[o]=walk(o); } "
        "peakOrder='mz'; document.getElementById('pkorder').value='mz'; renderPeaks(); "
        "return {lists, walks}; })()",
    )
    for order in ("mz", "abundance", "label"):
        _assert(
            nav["walks"][order] == nav["lists"][order][:3],
            "arrow-down does not follow the {} order: {} down a list of {}".format(
                order,
                nav["walks"][order],
                nav["lists"][order],
            ),
        )
    _assert(
        nav["lists"]["abundance"] != nav["lists"]["mz"]
        and nav["lists"]["label"] != nav["lists"]["mz"],
        "the fixture orders are indistinguishable, so the arrow test proves nothing",
    )

    # --- the unit selector keeps its slot across the two tabs, and the axis selector
    # that belongs to the trace tab sits hard right. The spectrum header carries no
    # right-hand control now that the interval lives in the sidebar, so what is worth
    # pinning down is that the two tab groups refuse to shrink: a flexible sibling
    # used to spread its deficit over every item on the line, which showed up as the
    # unit buttons sliding by 11-22 px when the tab changed. Details are closed first:
    # an expanded row is deliberately wide, and a wide row decides how much of the
    # window the sidebar takes, which is a separate question from this header. ---
    _browser(
        session,
        "eval",
        "(() => { const l=document.querySelector('#peaksbody .plist'); "
        "if(l && l.classList.contains('det')) document.getElementById('pkdetails').click(); })()",
    )
    header = _eval(
        session,
        "(() => { const q=document.querySelector('#qtabs'); const headEl=q.closest('h2'); "
        "const head=headEl.getBoundingClientRect(); "
        "const go=t => { document.querySelector('#maintabs button[data-tab='+t+']')"
        ".click(); return Math.round(q.getBoundingClientRect().left); }; "
        "const onTrace=go('trace'); "
        "const xa=document.getElementById('xaxiswrap').getBoundingClientRect(); "
        "const afterUnit=xa.left>q.getBoundingClientRect().right; "
        "const onSpec=go('spec'); go('trace'); "
        "const kids=[...headEl.children].filter(e=>getComputedStyle(e).display!=='none'); "
        "const rightmost=kids.reduce((a,e)=>Math.max(a,e.getBoundingClientRect().right),0); "
        "return {onTrace, onSpec, afterUnit, "
        "pinned:[...document.querySelectorAll('#maintabs,#qtabs')]"
        ".every(el=>getComputedStyle(el).flexShrink==='0'), "
        "isRightmost:Math.abs(rightmost-xa.right)<1, gap:head.right-xa.right}; })()",
    )
    _assert(
        abs(header["onTrace"] - header["onSpec"]) <= 1,
        "the Raw/Corrected/Conc selector shifts between the two tabs: "
        "{} px vs {} px".format(header["onTrace"], header["onSpec"]),
    )
    _assert(
        header["pinned"],
        "the tab groups are free to shrink, so their slot depends on what else is "
        "on the header line: " + str(header),
    )
    _assert(
        header["isRightmost"] and header["afterUnit"] and 0 < header["gap"] < 24,
        "the x-axis selector is not the right-hand control on the trace tab: "
        + str(header),
    )
    interval_controls = _eval(
        session,
        "({n:document.querySelectorAll('#scoperange,#specrange').length,"
        "head:document.querySelectorAll('#specrangewrap').length})",
    )
    _assert(
        interval_controls == {"n": 1, "head": 0},
        "there is not exactly one interval selector: " + str(interval_controls),
    )
    _browser(session, "eval", "document.getElementById('pkdetails').click()")

    # --- the selected unit drives the sidebar values and the spectrum alike ---
    # Details stays open: the abundance cells and the pills only render there
    _browser(
        session,
        "eval",
        "document.querySelector('#maintabs button[data-tab=spec]').click()",
    )
    units = _eval(
        session,
        "(() => { const p=peaks.find(x=>x.label==='Curated solvent'); "
        "cfg.humid=true; /* never depend on an earlier click */ "
        "setSpecRange('all'); /* the whole-run trace, not an interval's */ "
        "const shown=()=>Array.from(document.querySelectorAll('#peaksbody li'))"
        ".find(li=>li.querySelector('.lbl').value==='Curated solvent')"
        ".querySelector('.abundance').textContent; "
        "const traceMean=a=>Array.from(a).filter(v=>isFinite(v)).reduce((x,y)=>x+y,0)"
        "/a.length; "
        "document.querySelector('#qtabs button[data-q=raw]').click(); "
        "const raw=peakAbundance(p), rawText=shown(); "
        "document.querySelector('#qtabs button[data-q=cor]').click(); "
        "const cor=peakAbundance(p), corText=shown(); "
        "document.querySelector('#qtabs button[data-q=con]').click(); "
        "const con=peakAbundance(p), conText=shown(); "
        # the factor is the mean of K/I_p, not K over the mean I_p
        "const pinv=Array.from(PC.primary).filter(v=>v>0).map(v=>1/v); "
        "const meanRecip=pinv.reduce((a,b)=>a+b,0)/pinv.length; "
        # the humidity correction belongs to near-thermoneutral compounds only, so a
        # flagged and an unflagged compound must each match their own trace
        "const flags=p.flags; p.flags=['humid']; "
        "const conHum=peakAbundance(p), conHumTrace=traceMean(computeTraces(p).con); "
        "p.flags=['not-humid']; "
        "const conDry=peakAbundance(p), conDryTrace=traceMean(computeTraces(p).con); "
        "p.flags=flags; "
        "return {raw,cor,con,rawText,corText,conText,changed:rawText!==conText, "
        "corExpect:raw/interpT(p.mz), conExpect:raw/interpT(p.mz)*cfg.K*meanRecip, "
        "humidDiffers:Math.abs(conHum-conDry)>1e-12, conHum, conHumTrace, conDry, conDryTrace, "
        "qtabsVisible:getComputedStyle(document.querySelector('#qtabs')).display}; })()",
    )
    _assert(
        units["qtabsVisible"] != "none",
        "the Raw/Corrected/Conc selector is missing from the Mass spectrum tab",
    )
    _assert(
        units["changed"],
        "the sidebar abundance does not follow the selected unit: "
        + str([units["rawText"], units["corText"], units["conText"]]),
    )  # transmission is flat in this fixture, so only Conc must differ

    # Every setting the concentration conversion reads must move the sidebar, not
    # only the plots: a cached value from the previous setting is a wrong number.
    for control, value, unit in (
        ("K", "2.5", "con"),
        ("Vm", "49.0", "ug"),  # molar volume only enters through µg/m³
        ("kanchor", "3.4", "con"),
    ):
        edited = _eval(
            session,
            "(() => { const p=peaks.find(x=>x.label==='Curated solvent'); "
            "document.querySelector('#qtabs button[data-q=" + unit + "]').click(); "
            "const el=document.querySelector('#" + control + "'), keep=el.value; "
            "el.value='" + value + "'; "
            "el.dispatchEvent(new Event('change',{bubbles:true})); "
            "const cell=Array.from(document.querySelectorAll('#peaksbody li'))"
            ".find(li=>li.querySelector('.lbl').value==='Curated solvent')"
            ".querySelector('.abundance'); "
            "const dom=cell?cell.textContent:null, shown=peakAbundance(p); "
            "const tr=computeTraces(p)." + unit + "; "
            "const trace=Array.prototype.reduce.call(tr,(a,b)=>a+b,0)/tr.length; "
            "el.value=keep; el.dispatchEvent(new Event('change',{bubbles:true})); "
            "const back=computeTraces(p)." + unit + "; "
            "return {dom, shown, trace, restored:peakAbundance(p), "
            "traceBack:Array.prototype.reduce.call(back,(a,b)=>a+b,0)/back.length}; })()",
        )
        _assert(
            abs(edited["shown"] - edited["trace"]) <= abs(edited["trace"]) * 1e-9,
            "after editing {} the sidebar value is not the mean of the compound's own "
            "trace: {} vs {}".format(control, edited["shown"], edited["trace"]),
        )
        if edited["dom"] is not None:
            _assert(
                abs(float(edited["dom"]) - edited["shown"])
                <= abs(edited["shown"]) * 0.02,
                "the sidebar row kept an old value after {} was edited: DOM says {}, "
                "the trace says {}".format(control, edited["dom"], edited["shown"]),
            )
        _assert(
            abs(edited["restored"] - edited["traceBack"])
            <= abs(edited["traceBack"]) * 1e-9,
            "restoring {} did not restore the sidebar value: {} vs {}".format(
                control, edited["restored"], edited["traceBack"]
            ),
        )

    kinetic = _eval(
        session,
        "(() => { const p=peaks.find(x=>x.label==='Curated solvent'); "
        "const el=document.querySelector('#kinetic'), keep=el.checked; "
        "el.checked=!keep; el.dispatchEvent(new Event('change',{bubbles:true})); "
        "const shown=peakAbundance(p), tr=computeTraces(p).con; "
        "const trace=Array.prototype.reduce.call(tr,(a,b)=>a+b,0)/tr.length; "
        "el.checked=keep; el.dispatchEvent(new Event('change',{bubbles:true})); "
        "return {shown, trace, restored:peakAbundance(p)}; })()",
    )
    for field in ("shown", "restored"):
        _assert(
            abs(kinetic[field] - kinetic["trace"]) <= abs(kinetic["trace"]) * 1e-9,
            "{} the kinetic correction left a stale sidebar value: {} vs {}".format(
                "toggling" if field == "shown" else "restoring",
                kinetic[field],
                kinetic["trace"],
            ),
        )

    for key, expect in (("cor", "corExpect"), ("con", "conExpect")):
        _assert(
            abs(units[key] - units[expect]) <= abs(units[expect]) * 1e-9,
            "sidebar {} does not use the trace conversion: {} vs {}".format(
                key, units[key], units[expect]
            ),
        )
    _assert(
        units["humidDiffers"],
        "the humidity correction is not applied to a near-thermoneutral compound",
    )
    for sidebar, trace, kind in (
        ("conHum", "conHumTrace", "near-thermoneutral"),
        ("conDry", "conDryTrace", "ordinary"),
    ):
        _assert(
            abs(units[sidebar] - units[trace]) <= abs(units[trace]) * 1e-9,
            "the sidebar Conc for a {} compound disagrees with its own trace: "
            "{} vs {}".format(kind, units[sidebar], units[trace]),
        )

    # --- a generated name must not contradict the assigned formula ---
    naming = _eval(
        session,
        "(() => { const p=peaks.find(x=>x.label==='Curated solvent'); "
        "const row=()=>Array.from(document.querySelectorAll('#peaksbody li'))"
        ".find(li=>li.querySelector('.lbl').value===p.label); "
        "const set=v=>{ const i=row().querySelector('.lbl'); i.value=v; "
        "i.dispatchEvent(new Event('change',{bubbles:true})); }; "
        "const pills=()=>Array.from(row().querySelectorAll('.pill')).map(e=>e.textContent); "
        "set('toluene'); const wrongName=pills(); const note=document.querySelector('#idpanel').innerText; "
        "set('unknown m/z 120.000'); const unknownWithFormula=pills(); "
        "set('Curated solvent'); return {wrongName, unknownWithFormula, note, "
        "restored:pills(), label:p.label, formula:p.formula}; })()",
    )
    _assert(
        any("name" in t and "formula" in t for t in naming["wrongName"])
        and "disagree" in naming["note"],
        "a label that belongs to another formula is not flagged: "
        + str(naming["wrongName"]),
    )
    _assert(
        any("name" in t and "formula" in t for t in naming["unknownWithFormula"]),
        "'unknown' next to an assigned formula is presented as if it were identified: "
        + str(naming["unknownWithFormula"]),
    )
    _assert(
        not any("name" in t and "formula" in t for t in naming["restored"]),
        "the corrected label still shows a name/formula warning",
    )

    # --- a hand-drawn peak stays unknown until a formula candidate is accepted ---
    hand = _eval(
        session,
        "(() => { setTab('spec'); const rect=plotC.getBoundingClientRect(); "
        "const send=(t,x)=>plotC.dispatchEvent(new MouseEvent(t,{bubbles:true,cancelable:true,"
        "view:window,ctrlKey:true,clientX:rect.left+x,clientY:rect.top+120})); "
        "const add=(m,span)=>{ vSpec={lo:m-0.07, hi:m+0.07}; clampView(); drawSpec(); "
        "const x0=specXAtMz(m-span/2), x1=specXAtMz(m+span/2); "
        "send('mousedown',x0); send('mousemove',(x0+x1)/2); send('mousemove',x1); "
        "window.dispatchEvent(new MouseEvent('mouseup',{bubbles:true,"
        "clientX:rect.left+x1,clientY:rect.top+120})); "
        "const p=peaks[peaks.length-1]; return {apex:p.apex, label:p.label, formula:p.formula}; }; "
        "const named=add(93.0699, 0.04); const unnamed=add(115.55, 0.04); "
        "peaks.splice(peaks.length-2, 2); "
        "if(!peaks.some(p=>p.id===selId)) selId=peaks.length?peaks[0].id:null; "
        "renderPeaks(); initSpecView(); clampView(); drawSpec(); "
        "return {named, unnamed}; })()",
    )
    _assert(
        hand["named"]["label"].startswith("unknown m/z 93.")
        and hand["named"]["formula"] == "",
        "a hand-drawn library mass was assigned without candidate review: "
        + str(hand["named"]),
    )
    _assert(
        hand["unnamed"]["label"].startswith("unknown m/z 115.")
        and hand["unnamed"]["formula"] == "",
        "a hand-drawn peak off the library carries a name or formula it cannot justify: "
        + str(hand["unnamed"]),
    )
    _browser(session, "eval", "document.querySelector('#pkdetails').click()")


def _provenance_browser_pass() -> None:
    """Regress effective file sources and explicit reset provenance."""
    session = f"{SESSION}-provenance-{threading.get_ident()}"
    with tempfile.TemporaryDirectory(prefix="sniff-ms-viz-provenance-") as directory:
        directory_path = Path(directory)
        omitted = _synthetic_data()
        omitted_meta = omitted["meta"]
        omitted_meta["sources"]["K"] = "legacy default"
        omitted_meta["sources"]["molar_volume"] = "legacy default"
        omitted_meta["K_default"] = omitted_meta["K_file"]
        omitted_meta["K_source"] = omitted_meta["K_file_source"]
        omitted_meta["molar_volume"] = 24.465
        omitted_meta["molar_volume_file"] = 24.465
        omitted_meta["molar_volume_source"] = (
            "25 °C fallback (drift metadata unavailable)"
        )
        omitted_meta["molar_volume_file_source"] = omitted_meta["molar_volume_source"]
        omitted_path = directory_path / "omitted-calibration.html"
        omitted_path.write_text(viz.render_html(omitted), encoding="utf-8")

        equal = _synthetic_data()
        equal_meta = equal["meta"]
        equal_meta["K_file"] = equal_meta["K_default"]
        equal_meta["molar_volume_file"] = equal_meta["molar_volume"]
        equal_path = directory_path / "equal-calibration.html"
        equal_path.write_text(viz.render_html(equal), encoding="utf-8")

        try:
            _open(session, omitted_path.as_uri())
            _browser(session, "wait", "--load", "networkidle")
            _browser(session, "eval", "localStorage.setItem('ptrms-onboarded', '1')")
            _browser(session, "reload")
            _freeze_animations(session)  # the reload dropped the injected override
            _browser(session, "wait", "--load", "networkidle")
            _browser(session, "eval", "document.querySelector('#methodBtn').click()")
            omitted_state = _eval(
                session,
                "({text:document.querySelector('#methodlive').innerText, "
                "config:buildConfig()})",
            )
            _assert(
                "K: 0.800 (file acquisition calibration)" in omitted_state["text"],
                "omitted K reported the resolver default instead of the file source",
            )
            _assert(
                "molar volume: 24.46 L/mol (25 °C fallback "
                "(drift metadata unavailable))" in omitted_state["text"],
                "omitted molar volume reported the resolver default instead of its file source",
            )
            _assert_config_round_trip(omitted_state["config"])

            _open(session, equal_path.as_uri())
            _browser(session, "wait", "--load", "networkidle")
            _browser(session, "eval", "document.querySelector('#methodBtn').click()")
            initial = _eval(
                session,
                "({text:document.querySelector('#methodlive').innerText, "
                "config:buildConfig()})",
            )
            _assert(
                "K: 1.00 (curated config)" in initial["text"]
                and "molar volume: 24.50 L/mol (curated config)" in initial["text"],
                "equal configured/file values lost their initial config provenance",
            )
            _assert_config_round_trip(initial["config"])
            _browser(
                session,
                "eval",
                "(() => { const set=(id,v)=>{const e=document.querySelector('#'+id);"
                "e.value=v; e.dispatchEvent(new Event('change',{bubbles:true}));}; "
                "set('K','2.5'); set('Vm','25.5'); })()",
            )
            _browser(session, "wait", "700")
            edited = _eval(
                session, "({text:document.querySelector('#methodlive').innerText})"
            )
            _assert(
                "K: 2.50 (browser edit)" in edited["text"]
                and "molar volume: 25.50 L/mol (browser edit)" in edited["text"],
                "equal configured/file values did not retain browser edit provenance",
            )
            _browser(session, "eval", "document.querySelector('#resetK').click()")
            _browser(session, "wait", "700")
            reset = _eval(
                session,
                "({text:document.querySelector('#methodlive').innerText, "
                "config:buildConfig()})",
            )
            _assert(
                "K: 1.00 (file acquisition calibration)" in reset["text"]
                and "molar volume: 24.50 L/mol (file drift temperature)"
                in reset["text"],
                "reset-to-file provenance was hidden by equal initial values",
            )
            _assert_config_round_trip(reset["config"])
        finally:
            _browser(session, "close")


def _app_onboarding_browser_pass(data: dict[str, Any]) -> None:
    """Exercise durable app onboarding branches in a real browser."""
    onboarding = copy.deepcopy(data)
    session = f"{SESSION}-onboarding-{threading.get_ident()}"
    _ReviewHandler.posts = []
    _ReviewHandler.html = viz.render_html(
        onboarding, mode="app", auto_tour=False
    )
    server = socketserver.ThreadingTCPServer(("127.0.0.1", 0), _ReviewHandler)
    server.daemon_threads = True
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    try:
        port = server.server_address[1]
        _open(session, f"http://127.0.0.1:{port}/")
        _browser(session, "wait", "--load", "networkidle")
        _browser(
            session,
            "eval",
            "localStorage.setItem('ptrms-onboarded','1');sessionStorage.clear()",
        )

        _ReviewHandler.posts = []
        _ReviewHandler.html = viz.render_html(
            onboarding, mode="app", auto_tour=True
        )
        _browser(session, "reload")
        _browser(session, "wait", "--load", "networkidle")
        _browser(session, "wait", "800")
        migrated = _eval(
            session,
            "({tour:document.querySelector('#tourblock')?.hidden!==false})",
        )
        _assert(migrated["tour"], "legacy browser state reopened the guided tour")
        _assert(
            any(path == "/onboarding" for path, _body in _ReviewHandler.posts),
            "legacy browser onboarding state was not migrated to the app",
        )

        _browser(session, "eval", "document.querySelector('#tourBtn').click()")
        replayed = _eval(
            session,
            "({visible:document.querySelector('#tourblock').hidden===false})",
        )
        _assert(replayed["visible"], "the guided tour could not be replayed manually")
        _browser(session, "eval", "document.querySelector('#tskip').click()")

        _browser(session, "eval", "localStorage.removeItem('ptrms-onboarded')")
        _ReviewHandler.posts = []
        storage_failure = viz.render_html(onboarding, mode="app", auto_tour=True)
        storage_failure = storage_failure.replace(
            "<script>",
            "<script>Storage.prototype.getItem=()=>{throw new Error('blocked')};"
            "Storage.prototype.setItem=()=>{throw new Error('blocked')};</script><script>",
            1,
        )
        _ReviewHandler.html = storage_failure
        _browser(session, "reload")
        _browser(session, "wait", "--load", "networkidle")
        _browser(session, "wait", "800")
        fallback = _eval(
            session,
            "({visible:document.querySelector('#tourblock')?.hidden===false})",
        )
        _assert(
            fallback["visible"],
            "unavailable browser storage suppressed the first guided tour",
        )
        _assert(
            any(path == "/onboarding" for path, _body in _ReviewHandler.posts),
            "starting the automatic tour did not acknowledge app onboarding",
        )
    finally:
        _browser(session, "close")
        server.shutdown()
        thread.join(timeout=5)
        server.server_close()


def _interval_spectrum() -> list[int]:
    """Return an interval spectrum whose shared maximum exposes re-centring bugs."""
    spectrum = [1] * 13000
    # m/z 130.0284 is inside both clustered search neighbourhoods and dominates
    # them. The isolated control has its own shifted maximum at m/z 140.0199.
    spectrum[11403] = 100
    spectrum[11833] = 80
    return spectrum


def main() -> int:
    if shutil.which("agent-browser") is None:
        raise SystemExit(
            "agent-browser is required; install it with: npm i -g agent-browser"
        )

    data = _synthetic_data()
    _app_onboarding_browser_pass(data)
    _ReviewHandler.html = viz.render_html(data)
    _ReviewHandler.spectrum = json.dumps(data["spectrum"]).encode("ascii")
    _ReviewHandler.interval_spectrum = json.dumps(_interval_spectrum()).encode("ascii")
    _ReviewHandler.posts = []
    server = socketserver.ThreadingTCPServer(("127.0.0.1", 0), _ReviewHandler)
    server.daemon_threads = True
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    session = f"{SESSION}-{threading.get_ident()}"

    try:
        port = server.server_address[1]
        _open(session, f"http://127.0.0.1:{port}/")
        # The first-visit tour is useful to people but would make this regression
        # nondeterministic; mark it complete before reloading the generated page.
        _browser(session, "eval", "localStorage.setItem('ptrms-onboarded', '1')")
        _browser(session, "reload")
        _freeze_animations(session)  # the reload dropped the injected override
        _browser(session, "wait", "--load", "networkidle")
        # Discard any delayed request from the previous page/session before the
        # first controlled edit; every request below has a matching snapshot.
        _ReviewHandler.posts = []
        post_cursor = 0

        header = _eval(
            session,
            "({file:document.querySelector('#file').textContent, "
            "intervalNote:!!document.querySelector('#intnote')})",
        )
        _assert(
            header == {"file": "synthetic-review.h5", "intervalNote": False},
            "review header exposed a path or the removed interval merge note",
        )

        # The sidebar keeps the current m/z order by default, then supports a
        # descending mean integrated-signal order without changing config peak order.
        peak_order = _eval(
            session,
            "({value:document.querySelector('#pkorder').value, "
            "order:document.querySelector('.pkorder').innerText, "
            "labels:Array.from(document.querySelectorAll('#peaksbody input.lbl')).map(e=>e.value), "
            "values:Array.from(document.querySelectorAll('#peaksbody .mini')).map(e=>e.className), "
            "config:buildConfig()})",
        )
        _assert(peak_order["value"] == "mz", "peak order default is not m/z")
        _assert("order by" in peak_order["order"], "peak order label is wrong")
        _assert(
            set(peak_order["values"]) == {"mini mz"},
            "compact m/z order does not show only m/z",
        )
        _assert(
            peak_order["labels"]
            == [
                "Unassigned sole candidate",
                "Ambiguous mix",
                "Curated solvent",
                "Cluster component A",
                "Cluster component B",
                "Isolated control",
            ],
            "default peak order is not m/z",
        )
        preserved = peak_order["config"]
        first_peak = preserved["peaks"][0]
        first_range = preserved["ranges"][0]
        _assert(
            first_peak["apex"] == 99.999
            and first_peak["center"] == 100.001
            and first_peak["mass"] == 100.002
            and first_peak["future"] == {"keep": [1, 2, 3]}
            and first_peak["ion_role"]["model"] == "ion-role-v1"
            and first_peak["ion_role"]["kind"] == "isotope"
            and first_peak["ion_role"]["future_role"] == {"keep": True}
            and first_range["merged_gaps"] == [{"start": 2, "end": 2}]
            and first_range["future"] == {"keep": True},
            "nested authored peak/range fields were not preserved on save",
        )
        # Serve the saved snapshot as the next page's authored config, then perform
        # a real browser reload. This catches loss in both buildConfig and the page
        # initialiser, rather than merely checking a string in generated HTML.
        reload_data = copy.deepcopy(data)
        reload_data["config_base"] = preserved
        reload_data["peaks"][0]["_config_original"] = preserved["peaks"][0]
        reload_data["ranges"][0]["_config_original"] = preserved["ranges"][0]
        _ReviewHandler.html = viz.render_html(reload_data)
        _browser(session, "reload")
        _freeze_animations(session)
        _browser(session, "wait", "--load", "networkidle")
        reloaded = _eval(
            session,
            "(() => { const c=buildConfig(); const p=c.peaks[0], r=c.ranges[0]; "
            "return {apex:p.apex,center:p.center,mass:p.mass,future:p.future,"
            "role:p.ion_role,merged:r.merged_gaps,rangeFuture:r.future}; })()",
        )
        _assert(
            reloaded
            == {
                "apex": 99.999,
                "center": 100.001,
                "mass": 100.002,
                "future": {"keep": [1, 2, 3]},
                "role": {
                    "model": "ion-role-v1",
                    "kind": "isotope",
                    "status": "supported",
                    "label": "possible isotope",
                    "evidence": ["test evidence"],
                    "exclude_from_analyte_assignment": True,
                    "future_role": {"keep": True},
                },
                "merged": [{"start": 2, "end": 2}],
                "rangeFuture": {"keep": True},
            },
            "saved authored fields did not survive browser reload",
        )
        header_layout = _eval(
            session,
            "(() => { const head=document.querySelector('.pkhead'); const h=head.getBoundingClientRect(); "
            "const cs=getComputedStyle(head); "
            "const left=h.left+(parseFloat(cs.paddingLeft)||0), right=h.right-(parseFloat(cs.paddingRight)||0); "
            "const t=document.querySelector('.pktitle').getBoundingClientRect(); "
            "const c=document.querySelector('.pkcontrols').getBoundingClientRect(); "
            "const o=document.querySelector('.pkorder'); "
            "return {display:cs.display, "
            "direction:cs.flexDirection, "
            "titleLeft:Math.abs(t.left-left)<1, controlsRight:Math.abs(c.right-right)<1, "
            "twoRows:t.bottom<=c.top, orderDirection:getComputedStyle(o).flexDirection, "
            "controlOrder:Array.from(document.querySelector('.pkcontrols').children).map(e=>e.id||e.className)}; })()",
        )
        _assert(
            header_layout
            == {
                "display": "flex",
                "direction": "column",
                "titleLeft": True,
                "controlsRight": True,
                "twoRows": True,
                "orderDirection": "row",
                "controlOrder": ["pkorder", "pkcheckall", "pkdetails"],
            },
            "peak header controls are not laid out as requested",
        )
        splitter = _eval(
            session,
            "(() => { const r=document.querySelector('#plotresize'); "
            "return {role:r.getAttribute('role'), orientation:r.getAttribute('aria-orientation'), "
            "min:+r.getAttribute('aria-valuemin'), max:+r.getAttribute('aria-valuemax'), "
            "height:+document.querySelector('#plot').dataset.h, handlers:typeof r.onpointerdown==='function', "
            "bottomAligned:Math.abs(document.querySelector('.sidebar .card').getBoundingClientRect().bottom-"
            "document.querySelector('#intcard').getBoundingClientRect().bottom)<1}; })()",
        )
        _assert(
            splitter["role"] == "separator"
            and splitter["orientation"] == "horizontal"
            and splitter["min"] >= 60
            and splitter["max"] > splitter["min"]
            and splitter["height"] >= splitter["min"]
            and splitter["handlers"]
            and splitter["bottomAligned"],
            "plot/card resize splitter is missing or unbounded",
        )
        abundance_format = _eval(
            session,
            "({grouped:fmtAbundance(12345.6), ordinary:fmtAbundance(999.9)})",
        )
        _assert(
            abundance_format == {"grouped": "12,346", "ordinary": "999.90"},
            "abundance values are not formatted with thousands separators",
        )
        peak_toggle = _eval(
            session,
            "({label:document.querySelector('#pkcheckall').textContent, "
            "disabled:document.querySelector('#pkcheckall').disabled, "
            "checked:Array.from(document.querySelectorAll('#peaksbody input[data-a=use]'))"
            ".filter(e=>e.checked).length})",
        )
        _assert(
            peak_toggle == {"label": "Uncheck all", "disabled": False, "checked": 6},
            "peak toggle did not start in the all-checked state",
        )
        _browser(session, "eval", "document.querySelector('#pkcheckall').click()")
        peak_toggle = _eval(
            session,
            "({label:document.querySelector('#pkcheckall').textContent, "
            "checked:Array.from(document.querySelectorAll('#peaksbody input[data-a=use]'))"
            ".filter(e=>e.checked).length})",
        )
        _assert(
            peak_toggle == {"label": "Check all", "checked": 0},
            "peak toggle did not uncheck every peak",
        )
        _browser(
            session,
            "eval",
            "document.querySelector('#peaksbody input[data-a=use]').click()",
        )
        peak_toggle = _eval(
            session,
            "({label:document.querySelector('#pkcheckall').textContent, "
            "checked:Array.from(document.querySelectorAll('#peaksbody input[data-a=use]'))"
            ".filter(e=>e.checked).length})",
        )
        _assert(
            peak_toggle == {"label": "Check all", "checked": 1},
            "peak toggle did not represent a mixed selection",
        )
        _browser(session, "eval", "document.querySelector('#pkcheckall').click()")
        peak_toggle = _eval(
            session,
            "({label:document.querySelector('#pkcheckall').textContent, "
            "checked:Array.from(document.querySelectorAll('#peaksbody input[data-a=use]'))"
            ".filter(e=>e.checked).length})",
        )
        _assert(
            peak_toggle == {"label": "Uncheck all", "checked": 6},
            "peak toggle did not check every peak from a mixed selection",
        )
        _browser(
            session,
            "eval",
            "document.querySelector('#pkorder').value='abundance'; "
            "document.querySelector('#pkorder').dispatchEvent(new Event('change',{bubbles:true}))",
        )
        _browser(session, "wait", "100")
        peak_order = _eval(
            session,
            "({value:document.querySelector('#pkorder').value, "
            "labels:Array.from(document.querySelectorAll('#peaksbody input.lbl')).map(e=>e.value), "
            "values:Array.from(document.querySelectorAll('#peaksbody .mini')).map(e=>e.className), "
            "config:buildConfig()})",
        )
        _assert(peak_order["value"] == "abundance", "abundance order was not selected")
        _assert(
            set(peak_order["values"]) == {"mini abundance"},
            "compact abundance order does not show only abundance",
        )
        _assert(
            peak_order["labels"]
            == [
                "Cluster component A",
                "Cluster component B",
                "Curated solvent",
                "Isolated control",
                "Ambiguous mix",
                "Unassigned sole candidate",
            ],
            "abundance peak order is wrong",
        )
        _assert(
            peak_order["config"]["viz"]["peak_order"] == "abundance",
            "peak order was not saved in viz config",
        )
        _browser(
            session,
            "eval",
            "document.querySelector('#pkorder').value='label'; "
            "document.querySelector('#pkorder').dispatchEvent(new Event('change',{bubbles:true}))",
        )
        label_order = _eval(
            session,
            "({value:document.querySelector('#pkorder').value, "
            "labels:Array.from(document.querySelectorAll('#peaksbody input.lbl')).map(e=>e.value), "
            "config:buildConfig()})",
        )
        _assert(label_order["value"] == "label", "label order was not selected")
        _assert(
            label_order["labels"]
            == [
                "Ambiguous mix",
                "Cluster component A",
                "Cluster component B",
                "Curated solvent",
                "Isolated control",
                "Unassigned sole candidate",
            ],
            "label peak order is wrong",
        )
        _assert(
            label_order["config"]["viz"]["peak_order"] == "label",
            "label peak order was not saved in viz config",
        )
        _browser(
            session,
            "eval",
            "document.querySelector('#pkorder').value='mz'; "
            "document.querySelector('#pkorder').dispatchEvent(new Event('change',{bubbles:true})); "
            "document.querySelector('#pkdetails').click()",
        )
        details = _eval(
            session,
            "({mz:document.querySelectorAll('#peaksbody .mini.mz').length, "
            "abundance:document.querySelectorAll('#peaksbody .mini.abundance').length, "
            "wide:!getComputedStyle(document.querySelector('#app')).gridTemplateColumns"
            ".startsWith('360px'), tags:Array.from(document.querySelectorAll('#peaksbody .dc.pills'))"
            ".every(e=>e.scrollWidth<=e.clientWidth+1), rows:Array.from"
            "(document.querySelectorAll('#peaksbody .plist li')).every(e=>e.scrollWidth<=e.clientWidth+1), "
            "deletion:getComputedStyle(document.querySelector('#peaksbody .dc.del')).width})",
        )
        _assert_eq(
            details,
            {
                "mz": 6,
                "abundance": 6,
                "wide": True,
                "tags": True,
                "rows": True,
                "deletion": "28px",
            },
            "details view does not show both m/z and abundance",
        )
        _browser(session, "eval", "document.querySelector('#pkdetails').click()")

        # Check all display units, including irregular timestamp conversion. Ranges
        # remain integer cycles in the saved config while the card follows the axis.
        axis = _eval(
            session,
            "({unit:xAxisUnit, options:Array.from(document.querySelectorAll('#xaxisunit option')).map(o=>({value:o.value,disabled:o.disabled}))})",
        )
        _assert(
            axis["options"]
            == [
                {"value": "cycle", "disabled": False},
                {"value": "relative", "disabled": False},
                {"value": "absolute", "disabled": False},
            ],
            "x-axis selector options are wrong",
        )
        _assert(
            _eval(
                session,
                "({display:document.querySelector('#xaxiswrap').style.display, tab:tab})",
            )
            == {"display": "", "tab": "trace"},
            "x-axis selector is not visible on Signal over time",
        )
        _browser(
            session,
            "eval",
            "document.querySelector('#maintabs button[data-tab=spec]').click()",
        )
        _assert(
            _eval(
                session,
                "({display:document.querySelector('#xaxiswrap').style.display, tab:tab})",
            )
            == {"display": "none", "tab": "spec"},
            "x-axis selector remains visible on Mass spectrum",
        )
        spectrum_layout = _eval(
            session,
            "(() => { const h=document.querySelector('.main>.card h2').getBoundingClientRect(); "
            "const c=document.querySelector('.main>.card').getBoundingClientRect(); "
            "return {headerFlush:Math.abs(h.right-c.right)<=1, "
            "bottomAligned:Math.abs(document.querySelector('.sidebar .card').getBoundingClientRect().bottom-"
            "document.querySelector('#idcard').getBoundingClientRect().bottom)<1}; })()",
        )
        _assert(
            spectrum_layout == {"headerFlush": True, "bottomAligned": True},
            "Mass spectrum header or card is not aligned to the viewport",
        )
        _browser(
            session,
            "eval",
            "document.querySelector('#maintabs button[data-tab=trace]').click()",
        )
        _assert(
            _eval(
                session,
                "({display:document.querySelector('#xaxiswrap').style.display, tab:tab})",
            )
            == {"display": "", "tab": "trace"},
            "x-axis selector did not return on Signal over time",
        )
        _browser(
            session,
            "eval",
            "document.querySelector('#xaxisunit').value='relative'; document.querySelector('#xaxisunit').dispatchEvent(new Event('change',{bubbles:true}))",
        )
        converted = _eval(session, "({cycle:axisAtCycle(3), back:cycleAtAxis(9)})")
        _assert(
            converted["cycle"] == 9 and converted["back"] == 3,
            "irregular relative timestamps did not map cycles",
        )
        relative_axis = _eval(
            session,
            "({unit:xAxisUnit, text:document.querySelector('#rngtbl').innerText, config:buildConfig()})",
        )
        _assert(
            relative_axis["unit"] == "relative" and "0.0 s" in relative_axis["text"],
            "relative-time interval labels did not update",
        )
        _browser(
            session,
            "eval",
            "document.querySelector('#xaxisunit').value='absolute'; document.querySelector('#xaxisunit').dispatchEvent(new Event('change',{bubbles:true}))",
        )
        absolute_axis = _eval(
            session,
            "({unit:xAxisUnit, text:document.querySelector('#rngtbl').innerText, config:buildConfig()})",
        )
        _assert(
            absolute_axis["unit"] == "absolute"
            and "UTC" not in absolute_axis["text"]
            and "22:13:20" in absolute_axis["text"],
            "absolute interval labels did not update",
        )
        absolute_precision = _eval(
            session,
            "(() => { const vals=AXIS.absolute; "
            "const labels=vals.map(v=>formatAxis(v,true)); "
            "const mid=(vals[1]+vals[2])/2; "
            "return {distinct:vals[0]!==vals[1], mapping:cycleAtAxis(mid), labels, "
            "valid:labels.every(x=>/^\\d{2}:\\d{2}:\\d{2}$/.test(x))}; })()",
        )
        _assert(
            absolute_precision["distinct"]
            and 2.49 < absolute_precision["mapping"] < 2.51
            and absolute_precision["valid"],
            "sub-millisecond absolute axis mapping or labels are invalid",
        )
        _browser(
            session,
            "eval",
            "document.querySelector('#xaxisunit').value='cycle'; document.querySelector('#xaxisunit').dispatchEvent(new Event('change',{bubbles:true}))",
        )
        # Axis changes are already covered above; discard their debounced save before
        # checking that each subsequent edit saves one complete current snapshot.
        post_cursor = len(_ReviewHandler.posts)

        # Methods combines live provenance with detailed scientific explanations:
        # inspect curated non-default settings, edit the controls, and verify both save
        # and Done payloads.
        _browser(session, "eval", "document.querySelector('#methodBtn').click()")
        initial_methods = _eval(
            session,
            "({text:document.querySelector('#methodlive').innerText, "
            "panelText:document.querySelector('#methodpanel').innerText, kinetic:cfg.kinetic, "
            "R:cfg.R, Rphys:cfg.Rphys, primary:cfg.primarymz, humid:cfg.humid, "
            "humidChecked:document.querySelector('#humid').checked})",
        )
        _assert(
            "R integration windows" in initial_methods["text"],
            "Methods omits R: {!r}".format(initial_methods["text"][:200]),
        )
        _assert(
            "Mass calibration & drift" in initial_methods["panelText"],
            "Methods omits the detailed scientific sections",
        )
        _assert("Rphys" in initial_methods["text"], "Methods omits physical resolution")
        _assert(
            "Kinetic correction: on" in initial_methods["text"],
            "kinetic state is stale",
        )
        _assert(
            "curated config" in initial_methods["text"], "configured source is missing"
        )
        _assert(
            initial_methods["humid"] and initial_methods["humidChecked"],
            "humidity checkbox did not initialise from the effective config",
        )
        _assert(
            initial_methods["R"] == 1500 and initial_methods["primary"] == 19.022,
            "curated Methods settings did not initialise the browser",
        )
        _browser(
            session,
            "eval",
            "(() => { const set=(id,v)=>{const e=document.querySelector('#'+id);"
            "e.value=v; e.dispatchEvent(new Event('change',{bubbles:true}));}; "
            "set('R','1700'); set('Rphys','3300'); set('primarymz','20.022'); "
            "set('K','2.5'); set('Vm','25.5'); set('kanchor','2.2'); set('hump','0.8'); "
            "set('href','1.7'); document.querySelector('#kinetic').click(); "
            "document.querySelector('#humid').click(); document.querySelector('#wholewindows').click(); })()",
        )
        _browser(session, "wait", "1000")
        edited = _eval(
            session,
            "({text:document.querySelector('#methodlive').innerText, "
            "banner:document.querySelector('#stalebanner').innerText, "
            "config:buildConfig(), humidChecked:document.querySelector('#humid').checked, "
            "stale:staleSettings(), served:SERVED, protocol:location.protocol})",
        )
        _assert(edited["served"], "browser did not recognise the HTTP review server")
        _assert(
            "Kinetic correction: off" in edited["text"],
            "Methods did not update kinetic state",
        )
        _assert(
            "browser edit" in edited["text"],
            "Methods did not update calibration provenance",
        )
        _assert(
            "PREVIEW STALE" in edited["text"] and "PREVIEW STALE" in edited["banner"],
            "re-extraction edits were not prominently marked stale",
        )
        _assert(
            "Done-only" in edited["text"]
            and "preview 19.022" in edited["text"]
            and "final 20.022" in edited["text"],
            "stale Methods wording did not distinguish preview and final values",
        )
        _assert(
            edited["stale"] == {"primary": True, "Rphys": True, "windows": True},
            "stale settings did not track the edited re-extraction values",
        )
        _assert(
            "m/z 37 / m/z 20.022" in edited["text"],
            "Methods did not update humidity denominator",
        )
        _assert(not edited["humidChecked"], "humidity checkbox did not update")
        _assert(
            edited["config"]["analyze"]["R_phys"] == 3300,
            "R_phys control was not exported",
        )
        _assert(
            edited["config"]["analyze"]["primary_mz"] == 20.022,
            "primary m/z control was not exported",
        )
        _assert(
            edited["config"]["analyze"]["whole_run_windows"],
            "window mode control was not exported",
        )
        _assert_config_round_trip(edited["config"])
        post_cursor = _assert_complete_posts(edited["config"], "/save", post_cursor)
        # Reverting exactly to the immutable preview settings must clear every stale
        # marker and restore the original provenance, not leave a sticky warning.
        post_cursor = len(_ReviewHandler.posts)
        _browser(
            session,
            "eval",
            "(() => { const set=(id,v)=>{const e=document.querySelector('#'+id);"
            "e.value=v; e.dispatchEvent(new Event('change',{bubbles:true}));}; "
            "set('Rphys','3100'); set('primarymz','19.022'); "
            "document.querySelector('#wholewindows').click(); })()",
        )
        _browser(session, "wait", "1000")
        reverted = _eval(
            session,
            "({text:document.querySelector('#methodlive').innerText, "
            "banner:document.querySelector('#stalebanner').innerText, "
            "stale:staleSettings(), config:buildConfig()})",
        )
        _assert(
            reverted["stale"] == {"primary": False, "Rphys": False, "windows": False},
            "reverting to preview settings left stale state behind",
        )
        _assert(
            reverted["banner"] == "" and "PREVIEW STALE" not in reverted["text"],
            "reverting to preview settings left a stale warning behind",
        )
        _assert(
            "primary m/z: 19.022 (curated config)" in reverted["text"]
            and "Rphys" in reverted["text"],
            "reverting did not restore initial setting provenance",
        )
        post_cursor = _assert_complete_posts(reverted["config"], "/save", post_cursor)
        # Leave the final Done payload edited, so the smoke covers the actual
        # authoritative rerun configuration after a stale->fresh transition.
        post_cursor = len(_ReviewHandler.posts)
        _browser(
            session,
            "eval",
            "(() => { const set=(id,v)=>{const e=document.querySelector('#'+id);"
            "e.value=v; e.dispatchEvent(new Event('change',{bubbles:true}));}; "
            "set('Rphys','3300'); set('primarymz','20.022'); "
            "document.querySelector('#wholewindows').click(); })()",
        )
        _browser(session, "wait", "1000")
        final_preview = _eval(
            session, "({config:buildConfig(), stale:staleSettings()})"
        )
        _assert(
            final_preview["stale"] == {"primary": True, "Rphys": True, "windows": True},
            "final edited settings did not become stale again",
        )
        _assert_config_round_trip(final_preview["config"])
        post_cursor = _assert_complete_posts(
            final_preview["config"], "/save", post_cursor
        )
        post_cursor = len(_ReviewHandler.posts)
        _browser(
            session,
            "eval",
            "document.querySelector('#K').value=''; document.querySelector('#K').dispatchEvent(new Event('change',{bubbles:true}))",
        )
        _browser(session, "wait", "700")
        unavailable = _eval(
            session,
            "({methods:document.querySelector('#methodlive').innerText, "
            "note:document.querySelector('#calnote').innerText, config:buildConfig()})",
        )
        _assert(
            "Concentration is unavailable" in unavailable["methods"]
            and "unavailable" in unavailable["note"],
            "concentration availability did not update when K was cleared",
        )
        _assert_config_round_trip(unavailable["config"])
        post_cursor = _assert_complete_posts(
            unavailable["config"], "/save", post_cursor
        )
        post_cursor = len(_ReviewHandler.posts)
        _browser(
            session,
            "eval",
            "document.querySelector('#K').value='2.5'; document.querySelector('#K').dispatchEvent(new Event('change',{bubbles:true}))",
        )
        _browser(session, "wait", "700")
        filled = _eval(session, "({config:buildConfig()})")
        post_cursor = _assert_complete_posts(filled["config"], "/save", post_cursor)
        post_cursor = len(_ReviewHandler.posts)
        _browser(session, "eval", "document.querySelector('#resetK').click()")
        _browser(session, "wait", "700")
        reset = _eval(
            session,
            "({K:cfg.K,Vm:cfg.Vm,text:document.querySelector('#methodlive').innerText,"
            "note:document.querySelector('#calnote').innerText,config:buildConfig()})",
        )
        _assert(
            reset["K"] == 0.8 and reset["Vm"] == 24.0,
            "reset did not restore file-derived calibration",
        )
        _assert_config_round_trip(reset["config"])
        post_cursor = _assert_complete_posts(reset["config"], "/save", post_cursor)
        _assert(
            "file acquisition calibration" in reset["text"]
            and "file drift temperature" in reset["text"],
            "reset provenance is not file-derived",
        )
        _browser(session, "eval", "document.querySelector('#methodClose').click()")
        # Switch from the initial intervals view to the actual identification card.
        _browser(
            session,
            "find",
            "role",
            "button",
            "click",
            "--name",
            "Mass spectrum",
            "--exact",
        )

        # Select the clustered component and load an interval whose shared maximum
        # sits between both model centres. The browser must not turn that maximum
        # into two displayed apexes, while the isolated control still follows it.
        _browser(session, "find", "nth", "3", ".plist li", "click")
        _browser(
            session,
            "eval",
            "(() => { const s=document.querySelector('#scoperange'); "
            "s.value=[...s.options].find(o=>o.value!=='all').value; "
            "s.dispatchEvent(new Event('change',{bubbles:true})); })()",
        )
        _browser(session, "wait", "1000")
        clustered = _eval(
            session,
            "(() => { const ps=[peaks[3],peaks[4]]; "
            "const ws=ps.map(p => { const ax=dispApex(p), "
            "tb=windowTB(ax,p.winL,p.winR); "
            "return {centre:ax,left:tb2m(tb[0]),right:tb2m(tb[1])}; }); "
            "const intersection=Math.max(0,Math.min(ws[0].right,ws[1].right)-"
            "Math.max(ws[0].left,ws[1].left)); "
            "const union=Math.max(ws[0].right,ws[1].right)-"
            "Math.min(ws[0].left,ws[1].left); "
            "const isolatedRow=Array.from(document.querySelectorAll('#peaksbody li')).find(li=>"
            "li.querySelector('.lbl').value==='Isolated control'); "
            "return {windows:ws, overlap:intersection/union, isolated:dispApex(peaks[5]), "
            "isolatedValue:isolatedRow.querySelector('.mini').textContent, "
            "note:document.querySelector('#idpanel').innerText}; })()",
        )
        _assert(
            [item["centre"] for item in clustered["windows"]] == [130, 130.05],
            "clustered displayed centres moved onto the shared interval maximum",
        )
        _assert(
            clustered["overlap"] < 0.6,
            "clustered displayed windows became duplicate-like after interval loading",
        )
        _assert(
            clustered["isolated"] > 140.01,
            "isolated control did not re-centre on its interval maximum",
        )
        _assert(
            clustered["isolatedValue"] == "140.020",
            "sidebar m/z did not follow the selected spectrum",
        )
        # The sidebar number tracks the spectrum actually on screen. Measured in Raw
        # so the value is the interval integral itself, not a unit conversion.
        abundance_raw = _eval(
            session,
            "(() => { const q0=quant; "
            "document.querySelector('#qtabs button[data-q=raw]').click(); "
            "const v=peakAbundance(peaks[5]); "
            "document.querySelector('#qtabs button[data-q='+q0+']').click(); "
            "return {v}; })()",
        )["v"]
        _assert(
            abundance_raw > 25,
            "sidebar abundance did not follow the selected spectrum",
        )
        _assert(
            "deconvolved component" in clustered["note"]
            and "Unresolved fits are withheld" in clustered["note"],
            "clustered-peak wording is missing from the identification card",
        )
        unavailable = _eval(
            session,
            "(() => { const p=peaks[5],old=p.fit; p.fit={status:'unresolved'}; "
            "SHOWSPEC=SPEC; SPECWIN={lo:1,hi:NCYC}; const value=peakAbundance(p), "
            "note=peakAbundanceNote(p); p.fit=old; return {value,note}; })()",
        )
        _assert(
            unavailable["value"] is None and "Unavailable" in unavailable["note"],
            "an unresolved component was displayed as a combined spectrum integral",
        )
        preview = _eval(
            session, "({whole:M.whole_run_windows, first:rawTrace(peaks[0])})"
        )
        _assert(
            preview["whole"] is False and preview["first"]["2"] > preview["first"]["0"],
            "per-interval preview did not preserve interval-specific trace values",
        )
        _browser(session, "find", "nth", "0", ".plist li", "click")

        sole = _eval(
            session,
            "{conf:document.querySelector('#idconf').innerText, "
            "panel:document.querySelector('#idpanel').innerText, "
            "candidateLabel:document.querySelector('#idpanel .p').innerText, "
            "provenance:document.querySelector('#idpanel').innerText}",
        )
        _assert(
            "not assigned" in sole["conf"].lower(),
            "sole candidate is not visibly unassigned",
        )
        _assert(
            "only candidate" in sole["candidateLabel"],
            "sole candidate lost its explicit state",
        )
        _assert(
            "100%" not in sole["conf"], "sole candidate is presented as 100% confidence"
        )
        _assert(
            "formula ranking cannot determine structural isomers" in sole["provenance"],
            "identification card omits the provenance/isomer limitation",
        )

        # Select the unassigned multi-candidate peak, then inspect relative shares.
        _browser(session, "find", "nth", "1", ".plist li", "click")
        multiple = _eval(
            session,
            "{conf:document.querySelector('#idconf').innerText, "
            "shares:Array.from(document.querySelectorAll('#idpanel .p')).map(x=>x.innerText), "
            "provenance:document.querySelector('#idpanel').innerText}",
        )
        _assert(
            "not assigned" in multiple["conf"].lower(),
            "multi-candidate peak is not visibly unassigned",
        )
        _assert(
            multiple["shares"] == ["72% share", "28% share"],
            "candidate percentages are not relative shares",
        )
        _assert(
            "100%" not in multiple["conf"],
            "multi-candidate card shows confidence wording",
        )
        _assert(
            "formula ranking cannot determine structural isomers"
            in multiple["provenance"],
            "multi-candidate card omits provenance text",
        )

        # This is a real click on the second rendered candidate row.  It must
        # update the assignment rather than merely changing source/data strings.
        post_cursor = len(_ReviewHandler.posts)
        _browser(
            session,
            "eval",
            "document.querySelectorAll('#idpanel .cand')[1].click()",
        )
        _browser(session, "wait", "700")
        assigned = _eval(
            session,
            "{conf:document.querySelector('#idconf').innerText, config:buildConfig()}",
        )
        _assert_config_round_trip(assigned["config"])
        post_cursor = _assert_complete_posts(assigned["config"], "/save", post_cursor)
        _assert(
            "assigned" in assigned["conf"].lower(),
            "candidate click did not assign the formula",
        )
        _assert(
            "C5H12" in assigned["conf"], "clicked candidate formula was not retained"
        )

        # Finally, verify a curated assignment remains authoritative while its
        # label and formula are both visible in their respective UI locations.
        _browser(session, "find", "nth", "2", ".plist li", "click")
        curated = _eval(
            session,
            "{label:document.querySelector('.plist li.sel input.lbl').value, "
            "conf:document.querySelector('#idconf').innerText, "
            "chosen:document.querySelector('#idpanel .cand.chosen').innerText}",
        )
        _assert(curated["label"] == "Curated solvent", "curated label was overwritten")
        _assert(
            "assigned" in curated["conf"].lower() and "C3H8O" in curated["conf"],
            "curated formula is not authoritative",
        )
        _assert(
            "C3H8O" in curated["chosen"],
            "assigned formula is not marked in the candidate card",
        )

        # Interval edits, unit-aware values, sample-specific ticks and compound
        # naming are exercised against the same served page.
        _review_round_browser_pass(session)
        _ReviewHandler.reject_next_done = True
        _browser(session, "find", "role", "button", "click", "--name", "Done")
        _browser(session, "wait", "500")
        refused = _eval(
            session,
            "({text:document.querySelector('#doneov').innerText, "
            "retry:!!document.querySelector('#keepreviewing')})",
        )
        _assert(
            "Review was not finished" in refused["text"] and refused["retry"],
            "a rejected one-shot Done left the review spinning or without recovery",
        )
        _browser(session, "find", "role", "button", "click", "--name", "Keep reviewing")
        post_cursor = len(_ReviewHandler.posts)
        _browser(session, "find", "role", "button", "click", "--name", "Done")
        _browser(session, "wait", "800")
        posted = _eval(session, "({config:buildConfig()})")
        _assert_config_round_trip(posted["config"])
        save_posts = [body for path, body in _ReviewHandler.posts if path == "/save"]
        _assert(
            save_posts and save_posts[-1] == posted["config"],
            "latest save body differs from the browser's complete buildConfig()",
        )
        done_posts = [
            body for path, body in _ReviewHandler.posts[post_cursor:] if path == "/done"
        ]
        _assert(
            done_posts and all(body == posted["config"] for body in done_posts),
            "Done body differs from the browser's complete buildConfig()",
        )
        _standalone_browser_pass(data)
        _provenance_browser_pass()
        _start_screen_browser_pass()
        print("viz browser identification/configuration regression: OK")
        return 0
    finally:
        _browser(session, "close")
        server.shutdown()
        server.server_close()


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (AssertionError, json.JSONDecodeError) as exc:
        print(f"viz browser identification regression: FAIL: {exc}", file=sys.stderr)
        raise SystemExit(1)
