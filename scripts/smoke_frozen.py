#!/usr/bin/env python3
"""Smoke-test a packaged bundle: does the frozen app actually serve a review?

``sniff --help`` only proves the bootloader found argparse, so this drives the thing
that matters — the app opening a real ``.h5`` file and serving the review page — with
nothing but a synthetic file and a few seconds of patience.

Usage:  uv run python scripts/smoke_frozen.py dist/sniff/sniff.exe
        uv run python scripts/smoke_frozen.py "dist/Sniff.app/Contents/MacOS/sniff"
        uv run python scripts/smoke_frozen.py ".../sniff" --single-launch
        uv run python scripts/smoke_frozen.py dist/sniff/sniff --expect-browser
"""

from __future__ import annotations

import json
import os
import socket
import subprocess
import sys
import tempfile
import time
import typing as t
import urllib.error
import urllib.request
from pathlib import Path

import h5py
import numpy as np

LAUNCH_TIMEOUT = 30  # seconds for the process to start serving
OPEN_TIMEOUT = 120  # seconds for the app to build the review payload


def make_h5(path: Path) -> Path:
    """A tiny but honest IoniTOF-shaped file, and the config that belongs to it.

    The config is written beside it on purpose: reopening a file that already has a
    config is the path a reviewer hits every day, and it keeps this check out of the
    detection pipeline, which unit tests cover.
    """
    ncyc, nmz = 24, 15000
    anchor_masses = np.array([37.028405, 203.942993])
    a, b = 1000.0, 0.0
    trace = np.ones((ncyc, nmz), dtype=np.float64)
    bins = np.arange(nmz, dtype=np.float64)

    def add_peak(mass, height, cycles=slice(None)):
        centre = a * np.sqrt(mass) + b
        trace[cycles, :] += height * np.exp(-0.5 * ((bins - centre) / 1.5) ** 2)

    add_peak(21.022, 1e6)  # the reagent-ion isotope
    add_peak(31.0, 5e3)
    add_peak(31.0, 4e4, slice(8, None))  # an analyte that arrives halfway through
    add_peak(37.028405, 2e5)  # persistent water-cluster calibration anchor
    add_peak(203.942993, 2e5)  # persistent iodobenzene calibration anchor
    with h5py.File(path, "w") as h5:
        h5.create_dataset("SPECdata/Intensities", data=trace)
        h5.create_dataset("SPECdata/AverageSpec", data=trace.mean(axis=0))
        h5.create_dataset(
            "SPECdata/PCTime", data=np.arange(ncyc, dtype=np.float64)[:, None]
        )
        h5.attrs["InstrumentType"] = "IoniTof"
        h5.attrs["Single Spec Duration (ms)"] = [1000.0]
        h5.attrs["UTC_Offset"] = [0.0]
        h5.attrs["MassAxisType"] = "mz"
        h5.create_group("InstrumentSection")
        cal = h5.create_group("CALdata")
        cal.create_dataset("Mass_Use", data=np.ones(anchor_masses.size, dtype=bool))
        cal.create_dataset("Mass_MZ", data=anchor_masses)
        cal.create_dataset(
            "Mapping",
            data=np.column_stack([anchor_masses, a * np.sqrt(anchor_masses) + b]),
        )
    config = path.with_suffix(".json")
    config.write_text(
        json.dumps(
            {
                "peaks": [{"mz": 31.0, "label": "test analyte"}],
                "ranges": [
                    {"label": "sample_01", "start": 9, "end": 24, "unit": "cycle"}
                ],
                "viz": {"x_axis_unit": "cycle"},
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    return path


def get(url: str, timeout: float = 10.0):
    with urllib.request.urlopen(url, timeout=timeout) as response:
        return response.status, response.read()


def post(url: str, body=None) -> int:
    data = json.dumps(body or {}).encode("utf-8")
    request = urllib.request.Request(
        url, data, {"Content-Type": "application/json"}, method="POST"
    )
    try:
        with urllib.request.urlopen(request, timeout=10) as response:
            return response.status
    except urllib.error.HTTPError as exc:
        return exc.code


def free_port() -> int:
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        return sock.getsockname()[1]


def contents_directory(exe: Path) -> t.Optional[Path]:
    """Return the non-executable PyInstaller contents directory."""
    if exe.parent.name == "MacOS":
        # BUNDLE maps the one-dir payload into the standard app Resources folder.
        candidates = [exe.parent.parent / "Resources"]
    else:
        candidates = [exe.parent / "_internal"]
    return next((path for path in candidates if path.is_dir()), None)


def await_url(proc, port, timeout):
    """Block until the app serves its state endpoint. Returns (url, diagnostics)."""
    url = f"http://127.0.0.1:{port}/"
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if proc.poll() is not None:
            out, err = proc.communicate(timeout=5)
            return None, [text for text in (out, err) if text]
        try:
            status, _ = get(url + "api/state", timeout=0.5)
        except (OSError, urllib.error.URLError):
            time.sleep(0.05)
            continue
        if status == 200:
            return url, []
    return None, []


def await_logged_url(proc, log_path, offset, timeout):
    """Discover a console-free launch from its appended log, then poll its URL."""
    deadline = time.monotonic() + timeout
    prefix = "sniff: app running at "
    while time.monotonic() < deadline:
        if proc.poll() is not None:
            out, err = proc.communicate(timeout=5)
            return None, [text for text in (out, err) if text]
        if log_path.is_file():
            with log_path.open(encoding="utf-8") as handle:
                handle.seek(offset)
                lines = handle.read().splitlines()
            urls = [
                line.removeprefix(prefix) for line in lines if line.startswith(prefix)
            ]
            if urls:
                url = urls[-1]
                try:
                    status, _ = get(url + "api/state", timeout=0.5)
                except (OSError, urllib.error.URLError):
                    pass
                else:
                    if status == 200:
                        return url, []
        time.sleep(0.05)
    return None, []


def main(argv) -> int:
    options = set(argv[2:])
    if len(argv) < 2 or options - {"--single-launch", "--expect-browser"}:
        print(__doc__, file=sys.stderr)
        return 2
    exe = Path(argv[1]).resolve()
    single_launch = "--single-launch" in options
    expect_browser = "--expect-browser" in options
    if not exe.is_file():
        print(
            f"frozen app smoke: FAIL — no such executable file: {exe}", file=sys.stderr
        )
        return 1
    contents = contents_directory(exe)
    if contents is None:
        print(
            "frozen app smoke: FAIL — no PyInstaller contents directory for "
            f"the executable: {exe}",
            file=sys.stderr,
        )
        return 1
    if exe.parent.name != "MacOS":
        cli = exe.with_name(f"sniff-cli{exe.suffix}")
        if not cli.is_file():
            print(
                f"frozen app smoke: FAIL — no terminal launcher beside {exe}",
                file=sys.stderr,
            )
            return 1
        result = subprocess.run(
            [str(cli), "--help"], capture_output=True, text=True, timeout=LAUNCH_TIMEOUT
        )
        if result.returncode != 0 or "usage: sniff" not in result.stdout:
            print(
                f"frozen app smoke: FAIL — {cli.name} did not return CLI help",
                file=sys.stderr,
            )
            print(result.stdout, result.stderr, file=sys.stderr)
            return 1

    work = Path(tempfile.mkdtemp(prefix="sniff-frozen-smoke-"))
    h5 = make_h5(work / "run.h5")
    port = free_port()
    cmd = [str(exe), "app", str(h5), "--no-browser", "--port", str(port)]
    env = dict(os.environ, SNIFF_RECENT_PATH=str(work / "recent.json"))
    print(f"frozen app smoke: {' '.join(cmd)}", file=sys.stderr)

    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        errors="replace",
        env=env,
    )
    try:
        url, lines = await_url(proc, port, LAUNCH_TIMEOUT)
        if url is None:
            print(
                "frozen app smoke: FAIL — the app did not serve within the timeout",
                file=sys.stderr,
            )
            print("\n".join(lines), file=sys.stderr)
            return 1
        log_path = work / "log.txt"
        if (
            not log_path.is_file()
            or f"sniff: app running at {url}" not in log_path.read_text()
        ):
            print(
                "frozen app smoke: FAIL — the frozen app did not write its startup log",
                file=sys.stderr,
            )
            return 1
        base = url.rstrip("/")
        deadline = time.monotonic() + OPEN_TIMEOUT
        state = {}
        while time.monotonic() < deadline:
            state = json.loads(get(base + "/api/state")[1])
            if state.get("status") in ("ready", "error"):
                break
            time.sleep(0.5)
        if state.get("status") != "ready":
            print(
                f"frozen app smoke: FAIL — the file never opened: {state}",
                file=sys.stderr,
            )
            print("\n".join(lines), file=sys.stderr)
            return 1

        status, page = get(base + "/review")
        if status != 200 or b"const APPMODE = true" not in page:
            print(
                f"frozen app smoke: FAIL — /review returned {status} without the "
                "app-mode page",
                file=sys.stderr,
            )
            return 1

        recent = json.loads(get(base + "/api/recent")[1])
        if not recent or not Path(recent[0]["path"]).samefile(h5):
            print(
                f"frozen app smoke: FAIL — recents did not list the file: {recent}",
                file=sys.stderr,
            )
            return 1
        if post(base + "/close") != 200:
            print("frozen app smoke: FAIL — could not close the file", file=sys.stderr)
            return 1
        if single_launch:
            # Hosted macOS runners stall every second invocation of a windowed frozen
            # executable before Python starts. Verify the start screen in the process
            # which already proved that the frozen app can load and serve a review.
            status, start = get(base + "/")
            if status != 200 or b"Open an IONICON run" not in start:
                print(
                    f"frozen app smoke: FAIL — closing the file served {status} "
                    "without returning to the start screen",
                    file=sys.stderr,
                )
                return 1
        if post(base + "/shutdown") != 200:
            print(
                "frozen app smoke: FAIL — could not stop the first app", file=sys.stderr
            )
            return 1
        try:
            proc.wait(timeout=LAUNCH_TIMEOUT)
        except subprocess.TimeoutExpired:
            print(
                "frozen app smoke: FAIL — the first app did not stop",
                file=sys.stderr,
            )
            return 1

        if single_launch:
            print(
                f"frozen app smoke: OK  ({os.path.getsize(exe) // 1024} KiB launcher, "
                f"served {url} and returned to the start screen in one process)"
            )
            return 0

        # Second phase: use the same executable with no arguments at all, which is how
        # Finder and the Start Menu shortcut launch it. The command line would answer
        # that with usage text and exit 2, and a windowed bundle does it invisibly.
        log_offset = log_path.stat().st_size
        bare = subprocess.Popen(
            [str(exe)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
            errors="replace",
            env=dict(env, BROWSER=f"{sys.executable} -c pass"),
        )
        try:
            bare_url, bare_lines = await_logged_url(
                bare, log_path, log_offset, LAUNCH_TIMEOUT
            )
            if bare_url is None:
                print(
                    "frozen app smoke: FAIL — the start-screen launch did not start "
                    "an app (a double-click would do nothing)",
                    file=sys.stderr,
                )
                print("\n".join(bare_lines), file=sys.stderr)
                return 1
            status, start = get(bare_url.rstrip("/") + "/")
            # The start screen itself, not a control on it: what a file is opened
            # from is the thing this launch has to have got right.
            if status != 200 or b"Open an IONICON run" not in start:
                print(
                    f"frozen app smoke: FAIL — the bare launch served {status} "
                    "without the start screen",
                    file=sys.stderr,
                )
                return 1
            if expect_browser:
                surface = json.loads(get(bare_url.rstrip("/") + "/api/state")[1]).get(
                    "surface"
                )
                if surface != "browser":
                    print(
                        "frozen app smoke: FAIL — the bundle did not choose the "
                        f"browser surface: {surface!r}",
                        file=sys.stderr,
                    )
                    return 1
            post(bare_url.rstrip("/") + "/shutdown")
        finally:
            bare.kill()
            bare.wait(timeout=10)

        # Third phase: --window. A runner has no window server, which is exactly the
        # case that must degrade rather than die: the app logs why and serves anyway.
        # On a desktop machine the same command opens a real window and serves too.
        win_port = free_port()
        win = subprocess.Popen(
            [str(exe), "app", str(h5), "--window", "--port", str(win_port)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
            errors="replace",
            env=dict(env, BROWSER=f"{sys.executable} -c pass"),
        )
        windowed = "fell back to a browser tab"
        try:
            win_url, win_lines = await_url(win, win_port, LAUNCH_TIMEOUT)
            if win_url is None:
                print(
                    "frozen app smoke: FAIL — --window neither opened a window nor "
                    "started serving; a desktop without a window server would exit",
                    file=sys.stderr,
                )
                print("\n".join(win_lines), file=sys.stderr)
                return 1
            win_base = win_url.rstrip("/")
            # /review only exists once a file is open, so ask the state first
            deadline = time.monotonic() + OPEN_TIMEOUT
            state = {}
            while time.monotonic() < deadline:
                state = json.loads(get(win_base + "/api/state")[1])
                if state.get("status") in ("ready", "error"):
                    break
                time.sleep(0.5)
            if state.get("status") != "ready":
                print(
                    f"frozen app smoke: FAIL — the windowed app never opened the "
                    f"file: {state}",
                    file=sys.stderr,
                )
                return 1
            status, page = get(win_base + "/review")
            if status != 200:
                print(
                    f"frozen app smoke: FAIL — the windowed app served {status}",
                    file=sys.stderr,
                )
                return 1
            # Ask the app, not the log. A windowed desktop bundle can have no console,
            # so grepping captured output for the word "window" used to report success
            # precisely when nothing had been written — and the failure it was meant to
            # catch contained that word too.
            surface = json.loads(get(win_base + "/api/state")[1]).get("surface")
            if surface not in ("window", "browser"):
                print(
                    "frozen app smoke: FAIL — /api/state did not say how the app is "
                    f"being shown: {surface!r}",
                    file=sys.stderr,
                )
                return 1
            windowed = (
                "opened its own window"
                if surface == "window"
                else "fell back to a browser tab (a runner with no window server "
                "would; a desktop machine should not)"
            )
            post(win_base + "/shutdown")
        finally:
            win.kill()
            win.wait(timeout=10)

        print(
            f"frozen app smoke: OK  ({os.path.getsize(exe) // 1024} KiB launcher, "
            f"served {url}, a bare launch opened the start screen, and --window "
            f"{windowed})"
        )
        return 0
    finally:
        if proc.poll() is None:
            proc.kill()
        proc.wait(timeout=10)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
