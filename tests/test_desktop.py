"""Tests for the optional desktop window.

pywebview is never really started here: the module is replaced by a stand-in that
records the window it was asked for and releases the loop on request. Nothing in this
file may open a window — CI has no display, and a real one would block forever.
"""

import json
import subprocess
import sys
import threading
import types
import urllib.error
import urllib.request

import pytest

from sniff import analyze, app, desktop

TIMEOUT = 10  # seconds for anything that has to hand over to a thread


# --------------------------------------------------------------------------
# the stand-in for pywebview
# --------------------------------------------------------------------------
class _Event:
    """The ``+=`` half of ``webview.event.Event``, plus a way to fire it."""

    def __init__(self):
        self._handlers = []

    def __iadd__(self, handler):
        self._handlers.append(handler)
        return self

    def fire(self):
        for handler in list(self._handlers):
            handler()


class _FakeWindow:
    def __init__(self, webview, title, url=None, width=0, height=0, maximized=False):
        self.webview = webview
        self.title = title
        self.url = url
        self.width = width
        self.height = height
        self.maximized = maximized
        self.events = types.SimpleNamespace(closed=_Event())
        self.dialogs = []
        self.dialog_result = None
        self.destroy_calls = 0

    def create_file_dialog(self, dialog_type=10, file_types=(), **_kwargs):
        self.dialogs.append(
            {"dialog_type": dialog_type, "file_types": tuple(file_types)}
        )
        return self.dialog_result

    def destroy(self):
        self.destroy_calls += 1
        self.webview.close(self)


class _FakeWebview:
    """Enough of the pywebview module to drive the three entry points."""

    FileDialog = types.SimpleNamespace(OPEN=10)

    def __init__(self, fire_closed=True):
        self.windows = []
        self.loops = 0
        self.start_args = ()
        self.fire_closed = fire_closed
        self._released = threading.Event()

    def create_window(self, title, url=None, width=0, height=0, maximized=False, **_kwargs):
        window = _FakeWindow(
            self,
            title,
            url=url,
            width=width,
            height=height,
            maximized=maximized,
        )
        self.windows.append(window)
        return window

    def start(self, *args, **_kwargs):
        self.loops += 1
        self.start_args = args
        if args and callable(args[0]):
            args[0]()
        if not self._released.wait(TIMEOUT):
            raise AssertionError("the fake GUI loop was never released")

    def close(self, window):
        # A toolkit reports a close once: as an event, or only as the loop ending. Both
        # have to be survivable, so a test can pick either.
        if self.fire_closed:
            window.events.closed.fire()
        self._released.set()


class _Server:
    def __init__(self, base):
        self.base = base

    def get(self, path):
        with urllib.request.urlopen(self.base + path, timeout=TIMEOUT) as r:
            return r.status, json.loads(r.read().decode("utf-8"))

    def post(self, path, body=None):
        data = json.dumps(body or {}).encode("utf-8")
        req = urllib.request.Request(
            self.base + path, data, {"Content-Type": "application/json"}, method="POST"
        )
        try:
            with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
                return r.status, json.loads(r.read().decode("utf-8"))
        except urllib.error.HTTPError as exc:
            return exc.code, json.loads(exc.read().decode("utf-8"))


@pytest.fixture(autouse=True)
def _no_window_left_open(monkeypatch):
    """Window state from one test must never be handed to the next one."""
    monkeypatch.setattr(desktop, "_active", None)
    monkeypatch.setattr(app, "_surface", "browser")


@pytest.fixture
def server():
    """The app's own routes on a free port, with no window and no GUI."""
    httpd, session, url = app.make_server(port=0)
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    try:
        yield _Server(url), session
    finally:
        httpd.shutdown()
        thread.join(timeout=TIMEOUT)
        httpd.server_close()
        session.close()


def _wait_for(predicate, timeout=TIMEOUT):
    spare = threading.Event()
    spent = 0.0
    while spent < timeout:
        if predicate():
            return True
        spare.wait(0.02)
        spent += 0.02
    return False


def _no_webview(monkeypatch):
    """Make ``import webview`` fail, whichever way the machine is set up."""
    monkeypatch.setitem(sys.modules, "webview", None)


def _fake_gui(monkeypatch, fire_closed=True):
    """Replace pywebview with the stand-in, and hand it back."""
    fake = _FakeWebview(fire_closed=fire_closed)
    monkeypatch.setattr(desktop, "_import_webview", lambda: fake)
    return fake


def _install_window(monkeypatch, fire_closed=True):
    """A live window with no GUI behind it, as :func:`run_window` would leave it."""
    fake = _fake_gui(monkeypatch, fire_closed=fire_closed)
    window = fake.create_window(desktop.DEFAULT_TITLE)
    monkeypatch.setattr(desktop, "_active", window)
    return fake, window


def _recording_browser(calls):
    """A browser that remembers being opened, since an assert in the served thread
    would be lost with it."""

    def open_url(url=None, *_args):
        calls.append(url)
        return True

    return open_url


def _run_window_in_thread(**kwargs):
    worker = threading.Thread(target=desktop.run_window, kwargs=kwargs)
    worker.start()
    return worker


# --------------------------------------------------------------------------
# availability
# --------------------------------------------------------------------------
def test_availability_is_false_when_webview_cannot_be_imported(monkeypatch):
    _no_webview(monkeypatch)
    assert desktop.available() is False


def test_the_window_reports_a_missing_extra_rather_than_an_import_error(monkeypatch):
    _no_webview(monkeypatch)
    with pytest.raises(desktop.DesktopUnavailable) as excinfo:
        desktop.run_window("http://127.0.0.1:1/")
    assert not isinstance(excinfo.value, ImportError)
    assert "desktop extra" in str(excinfo.value)


def test_an_unimportable_webview_blames_the_packaging_not_the_user(monkeypatch):
    """The installed bundle's failure, and the message that hid it.

    ``webview`` is present, so ``find_spec`` finds it, but importing it raises for a
    reason that has nothing to do with the extra — in the shipped bundle it was
    ``bottle``, a dependency pywebview imports at module scope that PyInstaller never
    saw because this package reaches webview through a string. Telling someone to
    install the extra they already have is worse than saying nothing.
    """
    monkeypatch.setattr(
        desktop.importlib.util, "find_spec",
        lambda name, *a, **k: object() if name == "webview" else None,
    )

    def explode(name):
        raise ModuleNotFoundError("No module named 'bottle'", name=name)

    monkeypatch.setattr(desktop.importlib, "import_module", explode)
    with pytest.raises(desktop.DesktopUnavailable) as excinfo:
        desktop.run_window("http://127.0.0.1:1/")
    message = str(excinfo.value)
    assert "bottle" in message, message
    assert "pip install" not in message, "the fix is in the spec, not on the user's machine"


def test_the_package_imports_without_the_extra():
    """The core install is numpy + h5py + stdlib, so a clean interpreter must not
    have pywebview in it just because the app module was imported."""
    probe = (
        "import sys; import sniff.app, sniff.analyze;"
        "print('webview' in sys.modules)"
    )
    done = subprocess.run(
        [sys.executable, "-c", probe], capture_output=True, text=True, timeout=TIMEOUT
    )
    assert done.returncode == 0, done.stderr
    assert done.stdout.strip() == "False"


# --------------------------------------------------------------------------
# run_window
# --------------------------------------------------------------------------
def test_the_window_is_asked_for_once_with_the_url_and_no_chrome(monkeypatch):
    fake = _fake_gui(monkeypatch)
    worker = _run_window_in_thread(url="http://127.0.0.1:8/", width=1024, height=700)
    assert _wait_for(lambda: fake.loops == 1)
    window = fake.windows[0]
    assert len(fake.windows) == 1
    assert window.url == "http://127.0.0.1:8/"
    assert window.title == desktop.DEFAULT_TITLE
    assert (window.width, window.height) == (1024, 700)
    assert window.maximized is True
    assert fake.start_args == (desktop._set_macos_app_icon,)
    assert desktop.close_window() is True
    worker.join(TIMEOUT)
    assert not worker.is_alive()


def test_the_frozen_macos_window_sets_its_bundled_icon(monkeypatch):
    calls = {}
    image = object()
    loader = types.SimpleNamespace(
        alloc=lambda: types.SimpleNamespace(
            initWithContentsOfFile_=lambda path: calls.setdefault("path", path) and image
        )
    )
    application = types.SimpleNamespace(
        setApplicationIconImage_=lambda value: calls.setdefault("image", value)
    )
    modules = {
        "AppKit": types.SimpleNamespace(
            NSImage=loader,
            NSApplication=types.SimpleNamespace(
                sharedApplication=lambda: application
            ),
        ),
        "Foundation": types.SimpleNamespace(
            NSBundle=types.SimpleNamespace(
                mainBundle=lambda: types.SimpleNamespace(
                    pathForResource_ofType_=lambda name, kind: f"/{name}.{kind}"
                )
            )
        ),
    }
    monkeypatch.setattr(desktop.sys, "platform", "darwin")
    monkeypatch.setattr(desktop.sys, "frozen", True, raising=False)
    monkeypatch.setattr(desktop.importlib, "import_module", modules.__getitem__)

    assert desktop._set_macos_app_icon() is True
    assert calls == {"path": "/sniff.icns", "image": image}


def test_the_close_handler_fires_once_when_the_toolkit_reports_the_close(monkeypatch):
    _fake_gui(monkeypatch)
    closes = []
    worker = _run_window_in_thread(
        url="http://127.0.0.1:8/", on_close=lambda: closes.append(1)
    )
    assert _wait_for(lambda: desktop.current_window() is not None)
    assert desktop.close_window() is True
    worker.join(TIMEOUT)
    assert not worker.is_alive()
    assert closes == [1]
    assert desktop.current_window() is None


def test_the_close_handler_fires_once_when_only_the_loop_reports_it(monkeypatch):
    """A toolkit may end the loop without ever firing ``closed``. The session still
    has to be stopped exactly once, or the server waits on a window that is gone."""
    fake = _FakeWebview(fire_closed=False)
    monkeypatch.setattr(desktop, "_import_webview", lambda: fake)
    closes = []
    worker = _run_window_in_thread(
        url="http://127.0.0.1:8/", on_close=lambda: closes.append(1)
    )
    assert _wait_for(lambda: fake.loops == 1)
    fake.close(fake.windows[0])  # the loop ends, the event never fires
    worker.join(TIMEOUT)
    assert not worker.is_alive()
    assert closes == [1]


def test_closing_a_window_that_is_not_there_is_not_an_error():
    assert desktop.close_window() is False


def test_a_window_that_cannot_be_created_is_reported_as_unavailable(monkeypatch):
    fake = _fake_gui(monkeypatch)
    monkeypatch.setattr(fake, "create_window", lambda *_a, **_k: None)
    with pytest.raises(desktop.DesktopUnavailable):
        desktop.run_window("http://127.0.0.1:8/")
    assert desktop.current_window() is None


# --------------------------------------------------------------------------
# pick_file
# --------------------------------------------------------------------------
def test_pick_file_returns_the_choice_or_none(monkeypatch):
    _fake, window = _install_window(monkeypatch)

    window.dialog_result = ("/data/sniff.h5",)
    assert desktop.pick_file(window) == "/data/sniff.h5"
    assert window.dialogs[0]["file_types"] == desktop.H5_FILE_TYPES
    assert any("*.h5" in f for f in window.dialogs[0]["file_types"])

    window.dialog_result = None
    assert desktop.pick_file(window) is None


def test_a_dialog_that_will_not_open_is_reported_as_unavailable(monkeypatch):
    _fake, window = _install_window(monkeypatch)

    def boom(*_args, **_kwargs):
        raise RuntimeError("no dialog on this desktop")

    monkeypatch.setattr(window, "create_file_dialog", boom)
    with pytest.raises(desktop.DesktopUnavailable):
        desktop.pick_file(window)


# --------------------------------------------------------------------------
# serve_app, which owns the window-vs-browser decision
# --------------------------------------------------------------------------
def _serve_in_thread(monkeypatch, captured, **kwargs):
    """``serve_app()`` on its own thread, with the session it made left reachable."""
    real = app.make_server

    def make_server(**server_kwargs):
        httpd, session, url = real(**server_kwargs)
        captured.update(httpd=httpd, session=session, url=url)
        return httpd, session, url

    monkeypatch.setattr(app, "make_server", make_server)
    thread = threading.Thread(target=lambda: app.serve_app(**kwargs), daemon=True)
    thread.start()
    return thread


def test_the_app_starts_at_opening_screen_after_a_restart(tmp_path, monkeypatch):
    h5 = tmp_path / "run.h5"
    h5.touch()
    app.remember_active(str(h5))
    opened = []

    monkeypatch.setattr(
        app.Session, "open", lambda _session, path, **_kwargs: opened.append(path)
    )
    monkeypatch.setattr(app, "_install_quit_handlers", lambda stop: stop())
    app.serve_app(port=0, open_browser=False)

    assert opened == []
    assert app.load_active() is None


def test_window_shutdown_waits_for_a_delayed_close_time_save(tmp_path, monkeypatch):
    fake = _fake_gui(monkeypatch)
    captured = {}
    thread = _serve_in_thread(
        monkeypatch, captured, port=0, window=True, open_browser=False
    )
    assert _wait_for(lambda: len(fake.windows) == 1 and "session" in captured)
    session = captured["session"]
    session.path = str(tmp_path / "run.h5")
    session.config_path = tmp_path / "run.json"
    session.config = {}
    session.status = "ready"
    config = {
        "peaks": [{"mz": 42.0}],
        "ranges": [],
        "mass_axis_domain": "corrected",
        "mass_axis_version": 3,
    }

    old_page = session.begin_review_page()
    current_page = session.begin_review_page()
    save_entered = threading.Event()
    release_save = threading.Event()
    real_save = session.save_config

    def blocked_save(body, version=None, page=None):
        save_entered.set()
        assert release_save.wait(TIMEOUT)
        return real_save(body, version=version, page=page)

    session.save_config = blocked_save
    assert desktop.close_window() is True
    session.finish_close_save(old_page)  # an earlier refresh finishes during shutdown
    assert thread.is_alive()  # the old page generation must not release the wait
    api = _Server(captured["url"])
    response = {}
    request = threading.Thread(
        target=lambda: response.setdefault(
            "value",
            api.post(
                f"/save?version=1&page={current_page}&closing=1", config
            ),
        )
    )
    request.start()
    assert save_entered.wait(TIMEOUT)
    assert thread.is_alive()  # shutdown is draining the current page's save
    release_save.set()
    request.join(TIMEOUT)
    thread.join(TIMEOUT)

    assert response["value"][0] == 200
    assert not thread.is_alive()
    assert json.loads((tmp_path / "run.json").read_text(encoding="utf-8")) == config


def test_a_closed_window_ends_the_app_without_opening_a_browser(monkeypatch):
    fake = _fake_gui(monkeypatch)
    opened = []
    monkeypatch.setattr(app.webbrowser, "open", _recording_browser(opened))
    captured = {}
    thread = _serve_in_thread(monkeypatch, captured, port=0, window=True)
    assert _wait_for(lambda: len(fake.windows) == 1)
    assert desktop.close_window() is True
    thread.join(TIMEOUT)
    assert not thread.is_alive()
    assert opened == []  # a window that works must not also open a tab
    assert captured["session"].stop.is_set()


def test_a_window_that_cannot_start_loses_nothing(monkeypatch, capsys):
    """One line on stderr, then the browser route, then the same server as always."""
    _no_webview(monkeypatch)
    served = {}
    captured = {}

    def open_a_browser(url):
        with urllib.request.urlopen(url + "api/state", timeout=TIMEOUT) as r:
            served["state"] = r.status, json.loads(r.read().decode("utf-8"))
        captured["session"].stop.set()
        return True

    monkeypatch.setattr(app.webbrowser, "open", open_a_browser)
    thread = _serve_in_thread(monkeypatch, captured, port=0, window=True)
    thread.join(TIMEOUT)
    assert not thread.is_alive()

    lines = [
        line
        for line in capsys.readouterr().err.splitlines()
        if "window" in line.lower() or "browser" in line.lower()
    ]
    assert len(lines) == 1, lines
    assert served["state"][0] == 200
    assert served["state"][1]["status"] == "empty"


def test_a_plain_serve_app_still_just_opens_a_browser(monkeypatch, capsys):
    """The window must not change what a browser run does: no extra line, no window."""
    opened = []
    captured = {}

    def open_a_browser(url):
        opened.append(url)
        captured["session"].stop.set()
        return True

    monkeypatch.setattr(app.webbrowser, "open", open_a_browser)
    thread = _serve_in_thread(monkeypatch, captured, port=0)
    thread.join(TIMEOUT)
    assert not thread.is_alive()
    assert opened == [captured["url"]]
    assert "window" not in capsys.readouterr().err.lower()


# --------------------------------------------------------------------------
# the two directions between the page, the session and the window
# --------------------------------------------------------------------------
def test_shutdown_closes_the_window_exactly_once(server, monkeypatch):
    """The page stops the session and the session closes the window, while the close
    handler only ever stops the session: neither direction may restart the other."""
    api, session = server
    _fake, window = _install_window(monkeypatch)
    closes = []

    def on_close():
        closes.append("closed")
        session.stop.set()

    window.events.closed += on_close  # what serve_app wires when it opens the window

    assert api.post("/shutdown", {})[0] == 200
    assert window.destroy_calls == 1
    assert session.stop.is_set()

    assert api.post("/shutdown", {})[0] == 200  # a stopped app stops quietly
    assert window.destroy_calls == 1
    assert closes == ["closed"]  # the window was closed once, and reported once


def test_closing_the_window_stops_the_session_without_touching_the_shutdown_route(
    monkeypatch,
):
    """The other direction: a close from the window itself must reach nothing but the
    session's stop flag — no second close, and no browser to mop up."""
    fake = _fake_gui(monkeypatch)
    opened = []
    monkeypatch.setattr(app.webbrowser, "open", _recording_browser(opened))
    captured = {}
    thread = _serve_in_thread(monkeypatch, captured, port=0, window=True)
    assert _wait_for(lambda: len(fake.windows) == 1)

    window = fake.windows[0]
    fake.close(window)  # the reviewer used the window's own close box
    thread.join(TIMEOUT)

    assert not thread.is_alive()
    assert captured["session"].stop.is_set()
    assert opened == []
    assert window.destroy_calls == 0  # nothing reached back for another close


# --------------------------------------------------------------------------
# browse
# --------------------------------------------------------------------------
def test_browse_uses_the_window_dialog_when_there_is_a_window(server, monkeypatch):
    api, _session = server
    _fake, window = _install_window(monkeypatch)
    window.dialog_result = ("/data/sniff.h5",)

    def no_subprocess_dialog():
        raise AssertionError("a windowed app owns its own dialog")

    monkeypatch.setattr(app, "_pick_file", no_subprocess_dialog)
    assert api.post("/browse", {}) == (200, {"path": "/data/sniff.h5"})

    window.dialog_result = None
    assert api.post("/browse", {}) == (200, {"cancelled": True})
    assert len(window.dialogs) == 2


def test_browse_asks_the_machine_when_there_is_no_window(server, monkeypatch):
    api, _session = server
    called = []

    def machine_dialog():
        called.append(1)
        return "/data/other.h5"

    monkeypatch.setattr(app, "_pick_file", machine_dialog)
    assert api.post("/browse", {}) == (200, {"path": "/data/other.h5"})
    assert called == [1]
    assert desktop.current_window() is None


def test_a_broken_window_dialog_falls_back_to_the_machine(server, monkeypatch):
    """A window that cannot show a dialog must not cost the user the dialog."""
    api, _session = server
    _fake, window = _install_window(monkeypatch)

    def boom(*_args, **_kwargs):
        raise RuntimeError("no dialog in this window")

    monkeypatch.setattr(window, "create_file_dialog", boom)
    monkeypatch.setattr(app, "_pick_file", lambda: "/data/fallback.h5")
    assert api.post("/browse", {}) == (200, {"path": "/data/fallback.h5"})


# --------------------------------------------------------------------------
# the command line, which is where the window is asked for
# --------------------------------------------------------------------------
def _app_args(**overrides):
    """What argparse would have made of ``sniff app``, with the defaults it uses."""
    argv = {
        "h5": None,
        "port": 8765,
        "agent": None,
        "agent_timeout": 300.0,
        "no_browser": False,
        "window": False,
    }
    argv.update(overrides)
    return types.SimpleNamespace(**argv)


@pytest.mark.parametrize(
    ("platform", "frozen", "overrides", "expected"),
    [
        ("darwin", False, {}, False),  # a checkout keeps the browser it always had
        ("linux", False, {"no_browser": True}, False),
        ("linux", False, {"window": True}, True),
        ("darwin", True, {}, True),
        ("win32", True, {}, True),
        ("linux", True, {}, False),  # the portable bundle opens a browser directly
        ("linux", True, {"no_browser": True}, False),
        ("linux", True, {"window": True, "no_browser": True}, True),
    ],
)
def test_the_window_is_a_bundle_default_and_an_opt_in_elsewhere(
    monkeypatch, platform, frozen, overrides, expected
):
    monkeypatch.setattr(sys, "platform", platform)
    if frozen:
        monkeypatch.setattr(sys, "frozen", True, raising=False)
    else:
        monkeypatch.delattr(sys, "frozen", raising=False)
    assert analyze._wants_window(_app_args(**overrides)) is expected


def test_the_app_command_hands_the_choice_to_the_server(monkeypatch):
    seen = {}

    def fake_serve(**kwargs):
        seen.update(kwargs)

    monkeypatch.setattr(sys, "frozen", True, raising=False)
    monkeypatch.setattr(app, "serve_app", fake_serve)
    assert analyze.cmd_app(_app_args(window=True)) == 0
    assert seen["window"] is True
    assert seen["open_browser"] is True  # a bundle falls back to it if pywebview is out
