"""Persistent local review application.

The CLI flow is agent-driven: one command per file, and the server exits when the
expert clicks Done. This module is the other way round — one long-lived localhost
server that the user lives in, opening one file after another. Opening a file either
loads the config saved beside it or, the first time, runs the deterministic pipeline
to make one, optionally asking an agent endpoint to post-process it. Either way the
result is a config file on disk, which is what the UI then edits.

Measurement files stay local: they are read and written in place. Network access is
limited to the optional agent endpoint a user supplied and a metadata-only GitHub
Release check on startup. Compound-catalogue queries use the bundled offline database.
"""

from __future__ import annotations

import hashlib
import json
import os
import secrets
import shutil
import signal
import subprocess
import sys
import tempfile
import threading
import time
import urllib.error
import urllib.request
import webbrowser
from http import server
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import h5py

from . import __version__, brand, desktop, formula_id, panel, ptrms, viz
from . import update as updates
from .analyze import (
    analyze_config_to_csv,
    auto_peaks,
    auto_ranges,
    auto_ranges_note,
    interval_spectrum,
    resolve_analysis_settings,
    resolve_x_axis_unit,
)


def _recent_path() -> Path:
    """Where recent files are remembered. ``SNIFF_RECENT_PATH`` overrides it so a
    packaged build can be exercised (or a home folder kept clean) without patching
    Python in a frozen bundle."""
    override = os.environ.get("SNIFF_RECENT_PATH")
    return (
        Path(override).expanduser()
        if override
        else Path.home() / ".sniff" / "recent.json"
    )


def _active_path() -> Path:
    """Where the file to resume on the next launch is remembered."""
    override = os.environ.get("SNIFF_ACTIVE_PATH")
    return (
        Path(override).expanduser()
        if override
        else _recent_path().with_name("active.json")
    )


def _onboarding_path() -> Path:
    """Where app-wide onboarding completion is remembered."""
    override = os.environ.get("SNIFF_ONBOARDING_PATH")
    return (
        Path(override).expanduser()
        if override
        else _recent_path().with_name("onboarding.json")
    )


RECENT_PATH = _recent_path()
ACTIVE_PATH = _active_path()
ONBOARDING_PATH = _onboarding_path()
LEGACY_RECENT_PATH = Path.home() / ".ptr-ms" / "recent.json"
RECENT_LIMIT = 20
MASS_AXIS_FINGERPRINT_VERSION = 1
_MASS_AXIS_CACHE_FIELDS = ("mass_axis_calibration", "mass_axis_h5_fingerprint")

# Where an open's phases sit on its progress axis, measured on the real 2 GB /
# 20,725-cycle fixture rather than guessed: opening the h5 file and reading its
# header costs no measurable time, the deterministic pipeline about a second
# (auto_peaks 0.1 s, auto_ranges 0.8 s), and everything that is not the extraction
# pass in viz.build_viz_data another three and a half. The extraction is the rest,
# and it is 89 % of the 33 s an open takes — which is why it gets the bar's whole
# remaining span and reports cycles read, not anything smoother. It reads the run
# twice on a curated file (14.6 s streaming, 13.9 s re-centring the intervals), and
# both belong on the axis.
P_META = 0.03
P_CAL = 0.06
P_DETECT = 0.08
P_BUILD = viz.PREP_FRACTION  # 0.11: where the streaming starts, here and in viz


def _valid_config(value, *, require_mass_axis=False) -> bool:
    """True if a mapping is a current Sniff config.

    Loading still accepts an unmarked legacy config so it can be migrated after the
    file has passed calibration. Writes, however, must carry the corrected mass-axis
    marker; otherwise the next open could interpret the saved m/z values in the wrong
    domain.
    """
    if not isinstance(value, dict) or not ("peaks" in value or "ranges" in value):
        return False
    if not require_mass_axis:
        return True
    return (
        value.get("mass_axis_domain") == ptrms.MASS_AXIS_CONFIG_DOMAIN
        and value.get("mass_axis_version") == ptrms.MASS_AXIS_CONFIG_VERSION
    )


def _read_json(path: Path):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError):
        return None


def _config_key(config) -> str:
    """Return a stable identity for a JSON-compatible review config."""
    return json.dumps(config, sort_keys=True, separators=(",", ":"))


def _h5_fingerprint(path: Path):
    """Return a cheap identity for deciding whether saved calibration is reusable.

    Stable filesystem identity and timestamps catch normal edits and replacements. A
    bounded content sample keeps the check useful on filesystems that report no inode,
    without hashing a multi-gigabyte measurement on every open.
    """
    stat = path.stat()
    if not path.is_file():
        return None
    sample_size = 64 * 1024
    offsets = sorted(
        {
            0,
            max(0, stat.st_size // 2 - sample_size // 2),
            max(0, stat.st_size - sample_size),
        }
    )
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for offset in offsets:
            handle.seek(offset)
            digest.update(handle.read(sample_size))
    return {
        "version": MASS_AXIS_FINGERPRINT_VERSION,
        "path": os.path.normcase(str(path.resolve())),
        "device": int(stat.st_dev),
        "inode": int(stat.st_ino),
        "size": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
        "ctime_ns": int(stat.st_ctime_ns),
        "sample_sha256": digest.hexdigest(),
    }


def _h5_still_matches(path: Path, fingerprint) -> bool:
    """Check that an H5 still has its complete opening fingerprint."""
    if fingerprint is None:
        return False
    try:
        return _h5_fingerprint(path) == fingerprint
    except OSError:
        return False


def _cached_mass_axis(config, fingerprint):
    """Return a fully validated saved axis only for the unchanged source H5."""
    if (
        fingerprint is None
        or not _valid_config(config, require_mass_axis=True)
        or config.get("mass_axis_h5_fingerprint") != fingerprint
    ):
        return None
    diagnostics = config.get("mass_axis_calibration")
    tolerance = (
        diagnostics.get("formula_assignment_tolerance")
        if isinstance(diagnostics, dict)
        else None
    )
    if not isinstance(tolerance, dict) or tolerance.get("model") != (
        ptrms.FORMULA_TOLERANCE_MODEL
    ):
        return None
    try:
        return ptrms.mass_axis_from_dict(diagnostics)
    except ptrms.MassCalibrationError:
        return None


def _with_mass_axis_cache(config, mass_axis, fingerprint):
    """Attach server-owned calibration evidence and its source identity."""
    result = dict(config)
    result["mass_axis_calibration"] = mass_axis.to_dict()
    if fingerprint is None:
        result.pop("mass_axis_h5_fingerprint", None)
    else:
        result["mass_axis_h5_fingerprint"] = fingerprint
    return result


def _scrub_retired_review_fields(config):
    """Remove the retired checklist fields without discarding unknown config data."""
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


def _write_json(path: Path, value) -> None:
    """Durably replace a JSON file without exposing a partial write."""
    tmp = _stage_json(path, value)
    try:
        os.replace(tmp, path)
    finally:
        _remove_staged_json(tmp)


def _write_json_for_source(path: Path, value, source: Path, fingerprint) -> bool:
    """Publish JSON only while its source H5 retains the opening fingerprint."""
    previous = path.read_bytes() if path.exists() else None
    tmp = _stage_json(path, value)
    try:
        if not _h5_still_matches(source, fingerprint):
            return False
        os.replace(tmp, path)
        if _h5_still_matches(source, fingerprint):
            return True
        if previous is None:
            try:
                path.unlink()
            except FileNotFoundError:
                pass
        else:
            _write_bytes(path, previous)
        return False
    finally:
        _remove_staged_json(tmp)


def _stage_json(path: Path, value) -> Path:
    """Write and flush JSON to a private file ready for atomic publication."""
    return _stage_bytes(path, json.dumps(value, indent=2).encode("utf-8"))


def _write_bytes(path: Path, value: bytes) -> None:
    """Durably replace a file with already serialised bytes."""
    tmp = _stage_bytes(path, value)
    try:
        os.replace(tmp, path)
    finally:
        _remove_staged_json(tmp)


def _stage_bytes(path: Path, value: bytes) -> Path:
    """Write and flush bytes to a private file beside their destination."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, raw_tmp = tempfile.mkstemp(
        dir=str(path.parent), prefix=path.name + ".", suffix=".tmp"
    )
    tmp = Path(raw_tmp)
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(value)
            handle.flush()
            os.fsync(handle.fileno())
    except BaseException:
        _remove_staged_json(tmp)
        raise
    return tmp


def _remove_staged_json(path: Path) -> None:
    """Remove a private staged file if it has not already been published."""
    try:
        path.unlink()
    except FileNotFoundError:
        pass


def _replace_from(tmp, target: Path) -> None:
    """Publish a file written elsewhere (a temp name) onto its real path."""
    try:
        os.replace(tmp, target)
    finally:
        if os.path.exists(tmp):
            try:
                os.unlink(tmp)
            except OSError:
                pass


# A sniff summary starts with this header; anything else that turns up under the name
# we would like to write is somebody else's table and stays untouched.
_CSV_MARKERS = ("Variable", "Average(Corrected)")


def _csv_target(h5_path: str) -> Path:
    """Where an export of ``h5_path`` goes: ``<stem>.csv`` beside it, unless a file
    that is not a sniff summary already lives there — a Viewer or Excel export often
    does, and the review may be comparing against it."""
    target = Path(h5_path).with_suffix(".csv")
    if not target.exists():
        return target
    try:
        with target.open("r", encoding="utf-8-sig", errors="replace") as handle:
            head = handle.readline()
    except OSError:
        return target.parent / (target.stem + "-sniff.csv")
    if all(marker in head for marker in _CSV_MARKERS):
        return target
    return target.parent / (target.stem + "-sniff.csv")


def config_path_for(h5_path: str) -> Path:
    """The config file belonging to an h5 file: same name, same folder.

    ``~/d/sniff.h5`` -> ``~/d/sniff.json``. A ``<stem>-analysis-config.json`` written by
    the older CLI flow is honoured when no ``<stem>.json`` exists yet, so opening a
    previously reviewed file does not start a fresh pipeline run. A same-stem JSON
    that is not one of our configs is never overwritten: a ``<stem>.sniff.json`` is
    used instead.
    """
    p = Path(h5_path).expanduser()
    beside = p.with_suffix(".json")
    legacy_cli = p.parent / f"{p.stem}-analysis-config.json"
    legacy_product = p.with_suffix(".ptr.json")
    if beside.exists():
        if _valid_config(_read_json(beside)):
            return beside
        if legacy_product.exists() and _valid_config(_read_json(legacy_product)):
            return legacy_product
        return p.with_suffix(".sniff.json")
    if legacy_cli.exists() and _valid_config(_read_json(legacy_cli)):
        return legacy_cli
    if legacy_product.exists() and _valid_config(_read_json(legacy_product)):
        return legacy_product
    return beside


def _instrument(f) -> str:
    value = f.attrs.get("InstrumentType", "unknown")
    return value.decode(errors="replace") if isinstance(value, bytes) else str(value)


def bootstrap_config(
    h5_path: str,
    f=None,
    *,
    mass_axis=None,
    progress=None,
    should_stop=None,
    template_peaks=None,
    compounds_of_interest=None,
) -> dict:
    """Build a config from the file alone, with no agent and no judgement calls.

    Peaks and intervals come from the deterministic pipeline. ``f`` may be an
    already-open file, since opening a 2 GB run costs tens of seconds.
    ``progress`` and ``should_stop`` are the same pair :func:`ptrms.extract_traces`
    takes: detection is one call each for peaks and intervals, so it reports at the
    boundaries between them (auto_peaks 0.1 s, auto_ranges 0.8 s on the 2 GB run).
    """
    compounds_of_interest = formula_id.normalise_compounds_of_interest(
        compounds_of_interest, strict=False
    ) or []

    def _say(frac):
        if progress is not None:
            progress(max(0.0, min(1.0, float(frac))))

    def _halt():
        if should_stop is not None and should_stop():
            raise ptrms.AnalysisCancelled("the analysis was cancelled")

    own = f is None
    supplied_axis = mass_axis is not None
    source = h5py.File(h5_path, "r") if own else f

    def _detect(function, **kwargs):
        call_kwargs = dict(kwargs)
        if supplied_axis:
            call_kwargs["mass_axis"] = mass_axis
        try:
            return function(source, **call_kwargs)
        except TypeError as exc:
            # Preserve compatibility with callers that replace a detector with a
            # one-argument test double; production detectors accept these keywords.
            unsupported = [name for name in call_kwargs if name in str(exc)]
            if not unsupported:
                raise
            return function(source)

    try:
        _say(0.0)
        if mass_axis is None:
            mass_axis = ptrms.load_mass_axis(
                source, progress=progress, should_stop=should_stop
            )
        else:
            ptrms.validate_mass_axis(mass_axis)
        peak_options = {"assign_all_library": True}
        if compounds_of_interest:
            peak_options["compounds_of_interest"] = compounds_of_interest
        peaks = _detect(auto_peaks, **peak_options)
        adaptation = None
        if template_peaks is not None:
            peaks, adaptation = panel.adapt_peak_table(template_peaks, peaks)
        _say(0.1)
        _halt()  # detection is the only cancellable gap before the review data
        ranges = _detect(auto_ranges)
        _say(1.0)
        ncyc = int(source["SPECdata/Intensities"].shape[0])
        instrument = _instrument(source)
    finally:
        if own:
            source.close()

    settings = resolve_analysis_settings({})
    settings.update(
        {
            "peak_fit": "empirical-v1",
            "isotope_mode": "formula-envelope-v2",
            "isotope_abundance_basis": "unknown",
        }
    )
    config = {
        "analysis_schema_version": 3,
        "peaks": peaks,
        "ranges": ranges,
        "analyze": {k: v for k, v in settings.items() if k != "sources"},
        "viz": {"x_axis_unit": "cycle"},
        "mass_axis_domain": ptrms.MASS_AXIS_CONFIG_DOMAIN,
        "mass_axis_version": ptrms.MASS_AXIS_CONFIG_VERSION,
        "mass_axis_calibration": mass_axis.to_dict(),
        **(
            {"compounds_of_interest": compounds_of_interest}
            if compounds_of_interest
            else {}
        ),
        "diagnostics": {
            "n_peaks": len(peaks),
            "n_ranges": len(ranges),
            "ncyc": ncyc,
            "instrument": instrument,
            **({"peak_table_adaptation": adaptation} if adaptation is not None else {}),
        },
    }
    # Preserve automatic join provenance for CLI/config compatibility without showing
    # it in the review interface. Nothing merged, nothing to record.
    note = auto_ranges_note(ranges)
    if note:
        config["merge_note"] = note
    return config


def load_recent() -> list:
    value = _read_json(RECENT_PATH)
    if value is None and RECENT_PATH == _recent_path():
        # Read the old store without deleting it. The next write publishes the same
        # entries under ~/.sniff, while an interrupted migration leaves the old file
        # available for the previous release.
        value = _read_json(LEGACY_RECENT_PATH)
    if not isinstance(value, list):
        return []
    return [p for p in value if isinstance(p, str)]


def remember_recent(path) -> list:
    path = str(Path(path).expanduser().resolve())
    values = [p for p in load_recent() if p != path]
    values.insert(0, path)
    values = values[:RECENT_LIMIT]
    _write_json(RECENT_PATH, values)
    return values


def load_active():
    """Return the file whose review was active when Sniff last stopped."""
    value = _read_json(ACTIVE_PATH)
    if not isinstance(value, dict):
        return None
    path = value.get("file")
    return path if isinstance(path, str) and path else None


def remember_active(path: str) -> None:
    """Remember a successfully opened file for the next app launch."""
    _write_json(ACTIVE_PATH, {"file": str(Path(path).expanduser().resolve())})


def forget_active() -> None:
    """Forget the resume target after the reviewer explicitly leaves the file."""
    try:
        ACTIVE_PATH.unlink()
    except FileNotFoundError:
        pass


def auto_tour_pending() -> bool:
    """Whether the app has yet to show its one automatic guided tour."""
    try:
        return not ONBOARDING_PATH.is_file()
    except OSError:
        return True


def remember_tour_seen() -> None:
    """Persist that the automatic guided tour has been shown."""
    ONBOARDING_PATH.parent.mkdir(parents=True, exist_ok=True)
    ONBOARDING_PATH.touch(exist_ok=True)


class Session:
    """The one file currently open, plus the work in flight on it.

    Only one file is open at a time: a 2 GB run holds traces for tens of thousands of
    cycles, so a second open would double the footprint. Opening one closes the other.
    """

    def __init__(self):
        self.path = None
        self.config_path = None
        self.config = None
        self.payload = None
        self.mass_axis = None
        self._payload_config = None
        self._session_generation = 0
        self.status = "empty"  # empty | loading | ready | error
        self.stage = ""
        self.error = None
        self.error_details = None
        self.progress = None  # 0..1 while an open runs, None at every other time
        self.started_at = None  # monotonic clock, for the page's ETA
        self.export_result = None
        self.export_error = None
        self.agent_status = None
        self._file = None
        self._lock = threading.Lock()
        self._save_lock = threading.Lock()
        self._save_version = None
        self._page_condition = threading.Condition()
        self._page_counter = 0
        self._current_page = None
        self._closed_pages = set()
        self._opening = False
        self._exporting = False
        self._closing = False
        self._cancel = threading.Event()
        # A double-clicked bundle has no terminal to press Ctrl-C in, so stopping the
        # server is something the page has to be able to ask for.
        self.stop = threading.Event()

    # ---- opening -------------------------------------------------------------
    def _say(self, value):
        """Publish an open's progress. The axis may only ever move forward: one
        phase running after another must not make the bar go backwards."""
        value = max(0.0, min(1.0, float(value)))
        if self.progress is None or value > self.progress:
            self.progress = value

    def _halt(self):
        if self._cancel.is_set():
            raise ptrms.AnalysisCancelled("the open was cancelled")

    def _band(self, lo, hi):
        """A sink that puts one phase's own 0..1 fraction onto the open's axis."""

        def report(frac):
            self._say(lo + max(0.0, min(1.0, float(frac))) * (hi - lo))

        return report

    def _build_sink(self, prep_start):
        """A sink for build_viz_data's fractions, which are already an axis of their
        own: its phases take its first ``viz.PREP_FRACTION`` and the streaming pass
        the rest. From that fraction on the two axes agree — the streaming pass is
        89 % of the work in both — so only what precedes it is re-scaled, into
        whatever the open has not already spent."""

        def report(frac):
            frac = max(0.0, min(1.0, float(frac)))
            if frac <= viz.PREP_FRACTION:
                self._say(
                    prep_start
                    + (frac / viz.PREP_FRACTION) * (P_BUILD - prep_start)
                )
            else:
                self._say(frac)

        return report

    def cancel(self):
        """Ask an in-flight open to stop, and report whether there was one to stop.

        The flag is read by the ``should_stop`` callbacks the phases poll, so the
        open ends on its own thread at the next block boundary and leaves the
        session exactly as an open that never started.
        """
        with self._lock:
            opening = self._opening
        if opening:
            self._cancel.set()
        return opening

    def reserve_open(self) -> bool:
        """Atomically reserve the session for an asynchronous open."""
        with self._lock:
            if self._opening or self._exporting or self._closing:
                return False
            self._opening = True
            self._cancel.clear()
            self.status, self.stage, self.error = "loading", "Opening the file", None
            self.error_details = None
            return True

    def open(
        self,
        path,
        agent_url=None,
        agent_timeout=300.0,
        reserved=False,
        compounds_of_interest=None,
        template_peaks=None,
    ):
        """Load ``path``, making a config first if the file has never been reviewed.

        The h5 file is opened once and reused for detection and for the review data:
        reopening a large file costs the user another 30-90 s for nothing.

        Every phase is given the same progress sink and the same cancel flag, so the
        page can show a bar that reflects the work and get out of the way of a user
        who changed their mind. A cancelled open is not a failure: it leaves the
        session empty, with no error, and ready to open the same file again.
        """
        if not reserved and not self.reserve_open():
            raise RuntimeError("the app is busy with another operation")
        self.progress, self.started_at = 0.0, time.monotonic()
        try:
            self.close(reset_status=False)
            self.status, self.stage, self.error = "loading", "Opening the file", None
            self.error_details = None
            self.agent_status = None
            self.export_result = None
            path = str(Path(path).expanduser().resolve())
            source_path = Path(path)
            config_path = config_path_for(path)
            config = _read_json(config_path) if config_path.exists() else None
            if config is not None and not _valid_config(config):
                raise ValueError(f"{config_path} is not a sniff config")
            if config is not None and template_peaks is not None:
                raise ValueError(
                    "the target already has a saved review; open it normally instead"
                )
            config_prior_changed = False
            if config is not None:
                config, config_prior_changed = _canonicalise_config_compounds(config)
            saved_review = config is not None
            fingerprint = _h5_fingerprint(source_path)
            self._file = h5py.File(path, "r")
            self._say(P_META)

            # A current config may reuse the complete, validated anchor evidence only
            # while it remains tied to the same unchanged H5. Legacy configs still
            # require fresh calibration before any stored masses can be migrated.
            config_needs_write = config_prior_changed
            if config is not None:
                config, scrubbed = _scrub_retired_review_fields(config)
                config_needs_write = config_needs_write or scrubbed
            mass_axis = _cached_mass_axis(config, fingerprint)
            if mass_axis is None:
                self.stage = "Calibrating the mass axis"
                mass_axis = ptrms.load_mass_axis(
                    self._file,
                    progress=self._band(P_META, P_CAL),
                    should_stop=self._cancel.is_set,
                )
                if not _h5_still_matches(source_path, fingerprint):
                    raise RuntimeError("the H5 file changed while it was being opened")
                if config is not None:
                    config, _ = ptrms.migrate_config_mass_axis(config, mass_axis)
                    config = _with_mass_axis_cache(config, mass_axis, fingerprint)
                    config_needs_write = True
            self._say(P_CAL)

            prep_start = P_CAL
            if config is None:
                self.stage = "Detecting peaks and intervals"
                config = bootstrap_config(
                    path,
                    f=self._file,
                    mass_axis=mass_axis,
                    progress=self._band(P_CAL, P_DETECT),
                    should_stop=self._cancel.is_set,
                    template_peaks=template_peaks,
                    compounds_of_interest=compounds_of_interest,
                )
                self._halt()  # a cancel must not leave a half-made config on disk
                config = _with_mass_axis_cache(config, mass_axis, fingerprint)
                config_needs_write = True
                if compounds_of_interest:
                    config["compounds_of_interest"] = compounds_of_interest
                if agent_url:
                    config = self._ask_agent(config, path, agent_url, agent_timeout)
                    self._halt()
                prep_start = P_DETECT
            if compounds_of_interest is not None:
                if compounds_of_interest:
                    config["compounds_of_interest"] = compounds_of_interest
                else:
                    config.pop("compounds_of_interest", None)
                config_needs_write = True
            self.stage = (
                "Loading H5 data for the saved review"
                if saved_review
                else "Computing data for a new review"
            )
            payload = self._payload(
                path,
                config,
                mass_axis=mass_axis,
                progress=self._build_sink(prep_start),
                should_stop=self._cancel.is_set,
                assign_identity_defaults=not saved_review,
            )
            if not _h5_still_matches(source_path, fingerprint):
                raise RuntimeError("the H5 file changed while it was being opened")
            if config_needs_write and not _write_json_for_source(
                config_path, config, source_path, fingerprint
            ):
                raise RuntimeError("the H5 file changed while it was being opened")
            with self._save_lock:
                self.path, self.config_path, self.config = path, config_path, config
                self.payload = payload
                self.mass_axis = mass_axis
                self._payload_config = _config_key(config)
            # The resume pointer belongs to the successful open. Publish it before
            # "ready" lets /close proceed, otherwise an explicit close can delete the
            # pointer just before this thread recreates it.
            try:
                remember_active(path)
            except OSError:
                pass
            try:
                remember_recent(path)
            except OSError:
                # Either piece of bookkeeping may fail independently. In particular, a
                # broken recents store must not leave the previous file as the resume.
                pass
            # Publish readiness and the busy flag under the same lock. The page offers
            # its buttons the moment it sees "ready", and /close must never observe an
            # idle opening before the ready session itself has been published.
            with self._lock:
                self.status, self.stage = "ready", "Ready"
                self._opening = False
            return self.payload
        except ptrms.AnalysisCancelled:
            # Nothing was decided and nothing was half-written: the file closes and
            # the session is as though the open had never been asked for.
            self.close()
            return None
        except Exception as exc:
            self.close()
            self.status, self.error, self.stage = "error", str(exc), "Failed to open"
            self.error_details = (
                exc.diagnostics if isinstance(exc, ptrms.MassCalibrationError) else None
            )
            raise
        finally:
            with self._lock:
                self._opening = False
            self.progress = None

    def _ask_agent(self, config, path, agent_url, agent_timeout):
        """Let an attached agent post-process the automatic config. Its answer is a
        convenience: any failure at all leaves the deterministic config in place,
        because the alternative is losing the pipeline's work over one bad request."""
        self.stage = "Asking the agent to review it"
        body = json.dumps(
            {
                "file": path,
                "config": config,
                "diagnostics": config.get("diagnostics") or {},
            }
        ).encode("utf-8")
        try:
            req = urllib.request.Request(
                agent_url,
                body,
                {"Content-Type": "application/json"},
                method="POST",
            )
            with urllib.request.urlopen(req, timeout=agent_timeout) as response:
                answer = json.loads(response.read().decode("utf-8", "replace"))
            if isinstance(answer, dict):
                answer = answer.get("config", answer)
            if not isinstance(answer, dict) or not (
                answer.get("peaks") or answer.get("ranges")
            ):
                # Key presence is not enough here: an answer of {"peaks": []} would
                # replace real detected work with an empty panel.
                raise ValueError("the agent reply contained no peaks or ranges")
        except Exception as exc:  # any endpoint failure must not cost the user a file
            self.agent_status = (
                f"Agent review failed — using the automatic config "
                f"({type(exc).__name__}: {exc})"
            )[:300]
            return config
        # Agent responses are edits to the already corrected automatic config, not
        # pre-calibration exports. Keep the domain marker authoritative even when an
        # otherwise valid agent omits unknown top-level fields.
        answer, _ = _scrub_retired_review_fields(answer)
        answer.setdefault("mass_axis_domain", config["mass_axis_domain"])
        answer.setdefault("mass_axis_version", config["mass_axis_version"])
        if (
            answer["mass_axis_domain"] != ptrms.MASS_AXIS_CONFIG_DOMAIN
            or answer["mass_axis_version"] != ptrms.MASS_AXIS_CONFIG_VERSION
        ):
            self.agent_status = "Agent review failed — unsupported mass axis marker."
            return config
        for key in _MASS_AXIS_CACHE_FIELDS:
            if key in config:
                answer[key] = config[key]
        self.agent_status = "Agent review applied."
        return answer

    def _payload(
        self,
        path,
        config,
        *,
        mass_axis=None,
        progress=None,
        should_stop=None,
        assign_identity_defaults=False,
    ):
        settings = resolve_analysis_settings(config)
        return viz.build_viz_data(
            self._file,
            config.get("peaks", []),
            config.get("ranges", []),
            mass_axis=mass_axis,
            analysis_settings=settings,
            config_base=config,
            x_axis_unit=resolve_x_axis_unit(config),
            merge_note=config.get("merge_note") or "",
            progress=progress,
            should_stop=should_stop,
            assign_identity_defaults=assign_identity_defaults,
        )

    def close(self, reset_status=True):
        # Waiting for a save already accepted by the server keeps shutdown from clearing
        # the target out from underneath its request thread.
        with self._save_lock:
            if self._file is not None:
                try:
                    self._file.close()
                except OSError:
                    pass
            self._file = None
            self.path = None
            self.config_path = None
            self.config = None
            self.payload = None
            self.mass_axis = None
            self._payload_config = None
            self.progress = None
            self._save_version = None
            self._session_generation += 1
            with self._page_condition:
                self._current_page = None
                self._closed_pages.clear()
            # A closed file has no last export: a tab left open must not be told a run
            # finished when it belongs to a file that is no longer loaded.
            self.export_result = None
            self.export_error = None
            self.agent_status = None
            if reset_status:
                self.status, self.stage = "empty", ""

    def save_config(self, config, version=None, page=None):
        """Persist a current-page edit unless a newer edit arrived first.

        Returns ``False`` for an obsolete request. Browser teardown can race a pending
        debounced request, so the edit-time version rather than thread completion order
        decides which snapshot remains on disk.
        """
        with self._save_lock:
            return self._save_config_locked(config, version=version, page=page)

    def _save_config_locked(self, config, version=None, page=None):
        """Save while the caller holds ``_save_lock``."""
        if not self.accepts_review_page(page):
            raise RuntimeError("this review page is no longer current")
        if not self.config_path:
            raise RuntimeError("no file is open")
        if (
            version is not None
            and self._save_version is not None
            and version < self._save_version
        ):
            return False
        config, _ = _scrub_retired_review_fields(config)
        config, _ = _canonicalise_config_compounds(config)
        for key in _MASS_AXIS_CACHE_FIELDS:
            if key in self.config:
                config[key] = self.config[key]
        _write_json(self.config_path, config)
        self.config = config
        if version is not None:
            self._save_version = version
        return True

    def review_payload(self):
        """Rebuild review data when the saved config changed since the last render."""
        with self._save_lock:
            if not self.path or self.config is None or self.payload is None:
                raise RuntimeError("no file is open")
            config = self.config
            config_key = _config_key(config)
            if config_key == self._payload_config:
                return self.payload
            path = self.path
            mass_axis = self.mass_axis
            generation = self._session_generation
        payload = self._payload(path, config, mass_axis=mass_axis)
        with self._save_lock:
            if (
                generation != self._session_generation
                or path != self.path
                or config_key != _config_key(self.config)
            ):
                raise RuntimeError("the open file changed while refreshing the review")
            self.payload = payload
            self._payload_config = config_key
            return payload

    def peak_preview(self, lo, hi):
        """Return a snapped apex and candidates for a hand-drawn spectrum region."""
        with self._save_lock:
            if not self.path or self.config is None or self.mass_axis is None:
                raise RuntimeError("no file is open")
            path = self.path
            mass_axis = self.mass_axis
            settings = resolve_analysis_settings(self.config)
            compounds = self.config.get("compounds_of_interest")
        with h5py.File(path, "r") as source:
            return viz.preview_peak(
                source,
                lo,
                hi,
                R=settings["R"],
                mass_axis=mass_axis,
                compounds_of_interest=compounds,
            )

    def prepare_review_page(self):
        """Return one generation-consistent payload, target path and page token."""
        with self._save_lock:
            generation = self._session_generation
        payload = self.review_payload()
        with self._save_lock:
            if (
                generation != self._session_generation
                or payload is not self.payload
                or self._payload_config != _config_key(self.config)
            ):
                raise RuntimeError("the open file changed while rendering the review")
            token = self._begin_review_page_locked()
            return payload, str(self.config_path or ""), token

    def begin_review_page(self) -> str:
        """Issue a generation token for the review page being rendered now."""
        with self._save_lock:
            return self._begin_review_page_locked()

    def _begin_review_page_locked(self) -> str:
        """Issue a page token while the caller holds ``_save_lock``."""
        # Versions order requests within one page. A refreshed page starts a new
        # clock, while its token keeps every request from the prior page out.
        self._save_version = None
        with self._page_condition:
            self._page_counter += 1
            token = f"{self._page_counter}-{secrets.token_urlsafe(24)}"
            self._current_page = token
            self._closed_pages.clear()
            return token

    def current_review_page(self):
        """Return the generation token for the most recently rendered review page."""
        with self._page_condition:
            return self._current_page

    def accepts_review_page(self, page) -> bool:
        """Whether a page token belongs to the review currently being shown."""
        if page is None:
            return True  # programmatic API clients predate browser page tokens
        with self._page_condition:
            return page == self._current_page

    def finish_close_save(self, page) -> None:
        """Report that one page generation finished its teardown save attempt."""
        if page is None:
            return
        with self._page_condition:
            self._closed_pages.add(page)
            self._page_condition.notify_all()

    def wait_for_page_close_save(self, page, timeout=1.0) -> bool:
        """Wait briefly for one page generation's teardown snapshot."""
        if page is None:
            return True
        with self._page_condition:
            return self._page_condition.wait_for(
                lambda: page in self._closed_pages, timeout=timeout
            )

    def wait_for_close_save(self, timeout=1.0) -> bool:
        """Wait briefly for the currently rendered page's teardown snapshot."""
        if not self.config_path:
            return True
        return self.wait_for_page_close_save(
            page=self.current_review_page(), timeout=timeout
        )

    def reserve_export(self) -> bool:
        """Atomically reserve the session for an asynchronous export."""
        with self._lock:
            if (
                self._opening
                or self._exporting
                or self._closing
                or not self.path
                or self.config is None
            ):
                return False
            self._exporting = True
            self.status, self.stage, self.error = (
                "exporting",
                "Running the analysis",
                None,
            )
            self.export_result = None
            self.export_error = None
            return True

    def prepare_export(self, config=None, version=None, page=None):
        """Save and capture one immutable export request under the config lock."""
        with self._save_lock:
            if config is not None:
                saved = self._save_config_locked(config, version=version, page=page)
                if not saved:
                    raise RuntimeError("a newer review edit is already saved")
            elif not self.accepts_review_page(page):
                raise RuntimeError("this review page is no longer current")
            if not self.path or self.config is None:
                raise RuntimeError("no file is open")
            return self.path, self.config

    def cancel_export_reservation(self) -> None:
        """Return to a ready review when an admitted export cannot be started."""
        with self._lock:
            self._exporting = False
            self.status, self.stage = "ready", "Ready"

    def reserve_close(self) -> bool:
        """Atomically keep new work out while the current file closes."""
        with self._lock:
            if self._opening or self._exporting or self._closing:
                return False
            self._closing = True
            return True

    def finish_close(self) -> None:
        """Release a close reservation after success or failure."""
        with self._lock:
            self._closing = False

    @property
    def busy(self) -> bool:
        """True while an open or an export is in flight."""
        with self._lock:
            return self._opening or self._exporting or self._closing

    # ---- exporting -----------------------------------------------------------
    def export(self, path=None, config=None, reserved=False):
        """Run the full-precision analysis to ``<stem>.csv`` beside the file.

        Unlike the CLI's Done, nothing shuts down afterwards: the reviewer keeps
        working and exports again. Route callers pass the request snapshot captured
        during admission; direct callers capture it here after reserving the session.
        """
        if not reserved and not self.reserve_export():
            raise RuntimeError("the app is busy or no file is open")
        tmp = None
        try:
            if path is None or config is None:
                path, config = self.prepare_export()
            target = _csv_target(path)
            fd, tmp = tempfile.mkstemp(
                dir=str(target.parent), prefix=target.name + ".", suffix=".tmp"
            )
            os.close(fd)
            result = analyze_config_to_csv(path, config, tmp)
            _replace_from(tmp, target)
        except Exception as exc:
            if tmp is not None:
                try:
                    os.unlink(tmp)
                except OSError:
                    pass
            # The file is still open and still worth reviewing, so a failed export
            # reports the error without dropping the session into an error state.
            with self._lock:
                self.status, self.error, self.stage = "ready", str(exc), "Export failed"
                self._exporting = False
            self.export_error = str(exc)
            raise
        result["out"] = str(target)
        self.export_result = result
        with self._lock:
            self.status, self.stage = "ready", "Ready"
            self._exporting = False
        return result

    # ---- state ---------------------------------------------------------------
    def status_payload(self) -> dict:
        """Export progress in the shape the review page's overlay already polls, so
        the app and the one-shot CLI share one piece of UI code."""
        if self.status == "exporting":
            return {"status": "running"}
        if self.export_error:
            return {"status": "error", "error": self.export_error}
        if self.export_result:
            return {"status": "done", "out": self.export_result.get("out")}
        return {"status": "idle"}

    def state(self) -> dict:
        with self._lock:
            status = self.status
            opening = self._opening
        loading = status == "loading"
        return {
            "status": status,
            "stage": self.stage,
            "error": self.error,
            "error_details": self.error_details,
            "file": self.path,
            "config": str(self.config_path) if self.config_path else None,
            "agent_status": self.agent_status,
            "export": self.export_result,
            # The open's own bar: how far it has got — a float only while it runs,
            # since nothing else on this page has a fraction to report — and whether a
            # Cancel button would do anything at all right now.
            "progress": (
                float(self.progress)
                if loading and self.progress is not None
                else None
            ),
            "cancellable": bool(loading and opening),
            # "window" or "browser": how the user is looking at this app right
            # now. A bundle that meant to open a window and did not has to be able to say
            # so — otherwise the only evidence is a tab the user has to notice.
            "surface": surface(),
        }


def _compound_catalogue():
    """Return every compound name the landing page may accept."""
    return formula_id.compound_catalogue()


def _normalise_compounds_of_interest(raw):
    """Validate user-selected compounds against the bundled PTR Library."""
    return formula_id.normalise_compounds_of_interest(raw)


def _canonicalise_config_compounds(config):
    """Canonicalise a config's optional compound prior or reject it."""
    if "compounds_of_interest" not in config:
        return config, False
    compounds = _normalise_compounds_of_interest(config["compounds_of_interest"])
    result = dict(config)
    if compounds:
        result["compounds_of_interest"] = compounds
    else:
        result.pop("compounds_of_interest", None)
    return result, result != config


_START_TEMPLATE = """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>__PAGE_TITLE__</title>
<style>
:root{
  --bg:#fbf3e8;--card:#fffdf9;--sunk:#fff8ef;--fg:#173c3b;--mut:#5b706d;
  --line:#d8e5df;--acc:#1f6f6b;--accc:#fffdf9;
  --peach:#ffd9a8;--err:#a33b32;--errbg:#fff0e9;
  --ring:rgba(31,111,107,.34);--scrim:rgba(251,243,232,.88);
}
@media(prefers-color-scheme:dark){:root{
  --bg:#102322;--card:#173331;--sunk:#132b29;--fg:#effaf3;--mut:#a9c0b9;
  --line:#31514d;--acc:#71c3ad;--accc:#102322;
  --peach:#ffd9a8;--err:#ff9c8f;--errbg:#3b211e;
  --ring:rgba(113,195,173,.45);--scrim:rgba(16,35,34,.9);
}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);-webkit-font-smoothing:antialiased;
  font:14px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",system-ui,sans-serif}
main{max-width:920px;margin:0 auto;padding:clamp(28px,7vw,72px) 24px 44px}
.head{display:flex;gap:12px;align-items:center;margin-bottom:28px}
svg.brand{flex:none;width:42px;height:42px;border-radius:11px;
  box-shadow:0 1px 4px rgba(0,0,0,.25)}
.tag{font-weight:400;color:var(--mut);font-size:15px;letter-spacing:0}
h1{margin:0;font-size:20px;font-weight:650;letter-spacing:-.015em}
.lede{margin:0;color:var(--mut);max-width:38em}
.card{background:var(--card);border:1px solid var(--line);border-radius:18px;
  box-shadow:0 1px 1px rgba(16,24,40,.04),0 14px 35px -24px rgba(16,24,40,.42)}
.hero{display:grid;grid-template-columns:minmax(0,1fr) 260px;gap:28px;align-items:center;
  margin-bottom:30px}
.eyebrow{display:flex;align-items:center;gap:8px;margin:0 0 13px;color:var(--acc);
  font-size:11px;font-weight:700;letter-spacing:.13em;text-transform:uppercase}
.eyebrow i{width:8px;height:8px;border-radius:50%;background:var(--peach);
  box-shadow:0 0 0 4px rgba(255,217,168,.35)}
.hero h2{max-width:11em;margin:0 0 10px;font-size:clamp(30px,5vw,48px);line-height:1.02;
  letter-spacing:-.045em;font-weight:700}
.spectrum{position:relative;min-height:170px;padding:14px;border-radius:28px;
  background:var(--acc);overflow:hidden;box-shadow:0 16px 35px -22px rgba(31,111,107,.75)}
.spectrum::before,.spectrum::after{content:"";position:absolute;border-radius:50%;
  background:var(--peach);opacity:.9}
.spectrum::before{width:125px;height:125px;right:-32px;top:-44px}
.spectrum::after{width:70px;height:70px;left:-23px;bottom:-28px;background:#f6b89c}
.spectrum svg{position:relative;z-index:1;width:100%;height:140px}
.spectrum .trace{stroke-dasharray:420;stroke-dashoffset:420;animation:trace 1.8s ease-out forwards}
@keyframes trace{to{stroke-dashoffset:0}}
.now{display:flex;gap:14px;align-items:center;padding:14px 16px;margin-bottom:20px;
  border-color:var(--acc)}
.now .txt{min-width:0;flex:1}
.now b{display:block;font-weight:600;overflow:hidden;text-overflow:ellipsis;
  white-space:nowrap}
.now .sub{display:block;color:var(--mut);font-size:12px;overflow:hidden;
  text-overflow:ellipsis;white-space:nowrap}
.now .meta{min-width:0;overflow:hidden}
.now .meta em{overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
.pick{padding:28px 30px;text-align:left}
.pick h2{margin:0 2px 4px;font-size:18px;font-weight:650;letter-spacing:-.02em}
.pick p{margin:0 0 18px;color:var(--mut);font-size:13px}
.quick{display:flex;flex-wrap:wrap;gap:8px 16px;margin-top:18px;color:var(--mut);font-size:11px}
.quick span{display:inline-flex;align-items:center;gap:5px}
.quick i{width:6px;height:6px;border-radius:50%;background:var(--peach)}
.btn{appearance:none;border:0;border-radius:9px;background:var(--acc);color:var(--accc);
  font:inherit;font-weight:550;padding:9px 15px;cursor:pointer}
.btn:hover{filter:brightness(1.07)}
.btn:disabled{opacity:.55;cursor:default;filter:none}
.btn.sec{background:transparent;color:var(--fg);border:1px solid var(--line);font-weight:500}
.btn.sec:hover{background:var(--sunk)}
.btn:focus-visible,.link:focus-visible{outline:2px solid var(--ring);
  outline-offset:2px}
.meta{flex:none;text-align:right;font-size:12px;color:var(--mut)}
.meta em{display:block;font-style:normal}
.note{margin-top:18px;padding:11px 14px;border-radius:10px;background:var(--sunk);
  color:var(--mut);font-size:13px}
.note[hidden]{display:none}
.note.err{background:var(--errbg);color:var(--err)}
.bar{height:2px;margin-top:9px;border-radius:2px;background:var(--line);overflow:hidden}
.bar i{display:block;height:100%;width:35%;background:var(--acc);
  animation:slide 1.5s ease-in-out infinite}
@keyframes slide{from{transform:translateX(-100%)}to{transform:translateX(380%)}}
@media(prefers-reduced-motion:reduce){.bar i{animation:none;width:100%;opacity:.5}
  .spectrum .trace{animation:none;stroke-dashoffset:0}}
footer{display:flex;gap:12px;align-items:center;justify-content:space-between;
  margin-top:32px;color:var(--mut);font-size:12px}
.version{margin-left:auto;opacity:.55}
.link{background:none;border:0;padding:0;color:var(--mut);font:inherit;
  text-decoration:underline;cursor:pointer}
.link:hover{color:var(--fg)}
.link[hidden]{display:none}
#update{position:fixed;inset:0;z-index:10;display:grid;place-items:center;padding:24px;
  background:var(--scrim);overscroll-behavior:contain;
  -webkit-backdrop-filter:blur(3px);backdrop-filter:blur(3px)}
#update[hidden]{display:none}
.updatecard{width:min(470px,100%);background:var(--card);border:1px solid var(--line);
  border-radius:16px;padding:24px;box-shadow:0 24px 60px -28px rgba(16,24,40,.45)}
.updatecard h2{margin:0 0 7px;font-size:21px;letter-spacing:-.025em}
.updatecard p{margin:0;color:var(--mut)}
.updateactions{display:flex;justify-content:flex-end;gap:8px;margin-top:22px}
.updateerror{margin-top:12px;color:var(--err);font-size:13px}
#interest{position:fixed;inset:0;z-index:8;display:grid;place-items:center;padding:24px;
  background:var(--scrim);overscroll-behavior:contain;
  -webkit-backdrop-filter:blur(3px);backdrop-filter:blur(3px)}
#interest[hidden]{display:none}
.priorcard{width:min(570px,100%);max-height:min(720px,calc(100vh - 48px));overflow:auto;
  background:var(--card);border:1px solid var(--line);border-radius:16px;padding:24px;
  box-shadow:0 24px 60px -28px rgba(16,24,40,.45)}
.priorcard h2{margin:0 0 4px;font-size:20px;letter-spacing:-.02em}
.priorcard .intro{margin:0 0 5px;color:var(--mut)}
.priorfile{margin:0 0 18px;color:var(--mut);font-size:12px;overflow-wrap:anywhere}
.combobox{position:relative}
#compound-input{width:100%;border:1px solid var(--line);border-radius:9px;
  background:var(--sunk);color:var(--fg);font:inherit;padding:10px 12px}
#compound-input:focus{outline:2px solid var(--ring);outline-offset:1px;border-color:var(--acc)}
#compound-input[aria-invalid="true"]{border-color:var(--err)}
.suggestions{position:absolute;z-index:2;top:calc(100% + 4px);left:0;right:0;
  max-height:176px;overflow:auto;overscroll-behavior:contain;background:var(--card);
  border:1px solid var(--line);
  border-radius:10px;box-shadow:0 14px 30px -18px rgba(16,24,40,.55)}
.suggestions[hidden]{display:none}
.suggestion{display:flex;width:100%;justify-content:space-between;gap:12px;border:0;
  border-bottom:1px solid var(--line);background:transparent;color:var(--fg);
  padding:9px 11px;text-align:left;font:inherit;cursor:pointer}
.suggestion:last-child{border-bottom:0}.suggestion:hover,.suggestion.active{background:var(--sunk)}
.suggestion small{color:var(--mut);white-space:nowrap}
.inputstatus{min-height:20px;margin-top:5px;color:var(--mut);font-size:12px}
.inputstatus.err{color:var(--err)}
.chips{display:flex;flex-wrap:wrap;gap:7px;min-height:34px;margin:10px 0 4px}
.chip{display:inline-flex;align-items:center;gap:7px;max-width:100%;padding:5px 8px 5px 10px;
  border:1px solid var(--line);border-radius:999px;background:var(--sunk)}
.chip span{overflow:hidden;text-overflow:ellipsis;white-space:nowrap}.chip small{color:var(--mut)}
.chip button{border:0;background:none;color:var(--mut);font:18px/1 sans-serif;padding:0;cursor:pointer}
.priorhint{margin:8px 0 0;color:var(--mut);font-size:12px}
.prioractions{display:flex;align-items:center;justify-content:space-between;gap:10px;margin-top:22px}
.prioractions .right{display:flex;gap:8px;margin-left:auto}
/* An open blocks the whole screen. It is the one thing on this page that takes
   long enough to be worth leaving, so the page has to make leaving possible
   rather than let a second click start a second open behind the first. */
html.lock,html.lock body{overflow:hidden}
#ov{position:fixed;inset:0;z-index:9;display:grid;place-items:center;padding:24px;
  background:var(--scrim);overscroll-behavior:contain;
  -webkit-backdrop-filter:blur(3px);backdrop-filter:blur(3px)}
#ov[hidden]{display:none}
.ovcard{width:min(460px,100%);background:var(--card);border:1px solid var(--line);
  border-radius:14px;padding:20px;box-shadow:0 1px 1px rgba(16,24,40,.04),
  0 24px 60px -28px rgba(16,24,40,.45)}
.ovhead{display:flex;gap:12px;align-items:flex-start}
.ovhead .txt{min-width:0;flex:1}
.ovhead b{display:block;font-weight:600;font-size:15px;overflow:hidden;
  text-overflow:ellipsis;white-space:nowrap}
.ovhead .sub{display:block;color:var(--mut);font-size:12px;overflow:hidden;
  text-overflow:ellipsis;white-space:nowrap}
.ovmark{flex:none;width:16px;height:16px;margin-top:3px;border-radius:50%;
  border:2px solid var(--line);border-top-color:var(--acc);animation:spin .9s linear infinite}
.ovmark.stop{border-color:var(--err);border-top-color:transparent;animation:none;
  border-radius:0;background:none}
@keyframes spin{to{transform:rotate(360deg)}}
.ovart{height:42px;margin:14px 0 4px;border-radius:10px;background:var(--sunk);overflow:hidden}
.ovart svg{display:block;width:100%;height:100%}
.ovart .base{stroke:var(--line);stroke-width:2}
.ovart .peak{stroke:var(--acc);stroke-width:2.5;stroke-linecap:round;stroke-linejoin:round;
  fill:none;stroke-dasharray:90;stroke-dashoffset:90;animation:drawpeak 2.8s ease-in-out infinite}
.ovart .ion{fill:var(--peach);opacity:.7;animation:iondrift 3.2s ease-in-out infinite}
.ovart circle.ion:nth-of-type(2){animation-delay:-1.1s}
.ovart circle.ion:nth-of-type(3){animation-delay:-2.2s}
.ovart .breath{stroke:var(--peach);stroke-width:1.7;fill:none;stroke-linecap:round;
  opacity:.72;animation:breathe 3.4s ease-in-out infinite}
@keyframes drawpeak{0%,100%{stroke-dashoffset:90;opacity:.5}45%,70%{stroke-dashoffset:0;opacity:1}}
@keyframes iondrift{0%,100%{transform:translate(0,3px);opacity:.25}50%{transform:translate(9px,-3px);opacity:.85}}
@keyframes breathe{0%,100%{transform:translateX(-3px);opacity:.2}50%{transform:translateX(4px);opacity:.8}}
.ovstage{margin:10px 0 12px;color:var(--mut);font-size:13px;min-height:20px}
.pbar{height:6px;border-radius:4px;background:var(--line);overflow:hidden}
.pbar i{display:block;height:100%;width:0;background:var(--acc);border-radius:4px;
  transition:width .35s ease}
.ovmeta{display:flex;gap:10px;justify-content:space-between;margin-top:7px;
  color:var(--mut);font-size:12px;font-variant-numeric:tabular-nums}
.ovrow{display:flex;justify-content:flex-end;margin-top:16px}
#overr .msg{margin:10px 0 0;color:var(--err);font-size:13px;white-space:pre-wrap}
#ov.handoff .ovcard{animation:handoffcard .62s cubic-bezier(.22,.75,.25,1) both}
#ov.handoff .ovart{animation:handoffart .62s ease both}
@keyframes handoffcard{to{opacity:0;transform:scale(1.035) translateY(-6px)}}
@keyframes handoffart{to{opacity:0;transform:scale(1.08)}}
@media(prefers-reduced-motion:reduce){.ovmark,.ovart .peak,.ovart .ion,.ovart .breath{animation:none}
  .ovart .peak{stroke-dashoffset:0}.pbar i{transition:none}}
@media(max-width:620px){main{padding-top:28px}.hero{grid-template-columns:1fr;gap:20px}
  .spectrum{min-height:125px}.spectrum svg{height:100px}.pick{padding:23px 20px}
  .quick{margin-top:16px}
  .now{align-items:flex-start;row-gap:10px;flex-wrap:wrap}
  .now .txt,.now .meta{flex:1 1 100%}
  .now b,.now .sub{white-space:normal;overflow-wrap:anywhere;text-overflow:clip}
  #interest{padding:12px}.priorcard{padding:20px;max-height:calc(100vh - 24px)}
  .prioractions{align-items:stretch;flex-direction:column-reverse}
  .prioractions .right{width:100%;margin:0}.prioractions .right .btn{flex:1}
}
</style></head><body><main>
  <div class="head">__MARK__<h1>__APP_NAME__ <span class="tag">__TAGLINE__</span></h1></div>

  <div class="hero">
    <div>
      <h2>Find the story in your spectrum.</h2>
      <p class="lede">Open an IONICON run to review its peaks and intervals. Saved configs
        return exactly as you left them; new runs get a clear starting point.</p>
    </div>
    <div class="spectrum" aria-label="A stylised mass-spectrum peak" role="img">
      <svg viewBox="0 0 260 140" aria-hidden="true" focusable="false">
        <polyline class="trace" points="8,105 74,105 99,105 108,34 117,105 251,105" fill="none"
          stroke="#eafaf6" stroke-width="4" stroke-linecap="round" stroke-linejoin="round"/>
      </svg>
    </div>
  </div>

  <div id="now"></div>

  <div class="card pick">
    <h2>Open an IONICON run</h2>
    <button class="btn" id="browse" type="button">Browse this computer&hellip;</button>
    <div class="quick" id="quick-help"><span><i aria-hidden="true"></i>Runs locally</span>
      <span><i aria-hidden="true"></i>No measurement uploads</span>
      <span><i aria-hidden="true"></i>Measurement files stay on this computer</span>
    </div>
  </div>

  <div class="note" id="state" role="status" aria-live="polite" hidden></div>

  <footer>
    <a class="link" id="quit" href="#" hidden>Stop the app</a>
    <span class="version">v__VERSION__</span>
  </footer>
</main>

<div id="update" role="dialog" aria-modal="true" aria-labelledby="update-title" hidden>
  <div class="updatecard">
    <h2 id="update-title">A new Sniff is ready</h2>
    <p id="update-copy"></p>
    <p class="updateerror" id="update-error" role="alert" hidden></p>
    <div class="updateactions">
      <button class="btn sec" id="install-later" type="button">Install later</button>
      <button class="btn" id="install-update" type="button">Install update</button>
    </div>
  </div>
</div>

<div id="interest" role="dialog" aria-modal="true" aria-labelledby="interest-title" hidden>
  <div class="priorcard">
    <h2 id="interest-title">Any compounds of particular interest?</h2>
    <p class="intro">Add compounds that are especially plausible for this sampling
      context. Sniff will use them to prioritise candidates, not as proof of identity.</p>
    <p class="priorfile" id="interest-file"></p>
    <div class="combobox">
      <input id="compound-input" type="text" autocomplete="off" role="combobox"
        aria-autocomplete="list" aria-controls="compound-suggestions"
        aria-expanded="false" placeholder="Start typing a compound name">
      <div class="suggestions" id="compound-suggestions" role="listbox" hidden></div>
    </div>
    <div class="inputstatus" id="compound-status" role="status" aria-live="polite"></div>
    <div class="chips" id="compound-chips" role="list"
      aria-label="Selected compounds"></div>
    <p class="priorhint">Press Tab, Enter, or comma to add one recognised name. You can
      also paste names separated by commas or new lines.</p>
    <div class="prioractions">
      <button class="btn sec" id="interest-back" type="button">Choose another file</button>
      <div class="right">
        <button class="btn sec" id="interest-skip" type="button">Skip</button>
        <button class="btn" id="interest-continue" type="button" disabled>Continue</button>
      </div>
    </div>
  </div>
</div>

<div id="ov" role="dialog" aria-modal="true" aria-labelledby="ovname" hidden>
  <div class="ovcard">
    <div id="ovload">
      <div class="ovhead">
        <span class="ovmark" aria-hidden="true"></span>
        <div class="txt"><b id="ovname">Opening a file</b>
          <span class="sub" id="ovdir"></span></div>
      </div>
      <div class="ovart" aria-hidden="true">
        <svg viewBox="0 0 320 42" focusable="false">
          <path class="base" d="M12 28H308" fill="none"/>
          <path class="peak" d="M12 28H137L151 8L165 28H308"/>
          <path class="breath" d="M184 18c12-9 24-9 35 0"/>
          <circle class="ion" cx="205" cy="28" r="2.5"/>
          <circle class="ion" cx="248" cy="28" r="2.5"/>
          <circle class="ion" cx="278" cy="28" r="2.5"/>
        </svg>
      </div>
      <p class="ovstage" id="ovstage" role="status" aria-live="polite"></p>
      <div class="pbar" role="progressbar" id="ovbar"
           aria-valuemin="0" aria-valuemax="100" aria-valuenow="0"><i id="ovfill"></i></div>
      <div class="ovmeta"><span id="ovpct">0%</span><span id="oveta"></span></div>
      <div class="ovrow"><button class="btn sec" id="cancel" type="button">Cancel</button>
      </div>
    </div>
    <div id="overr" hidden>
      <div class="ovhead">
        <span class="ovmark stop" aria-hidden="true"></span>
        <div class="txt"><b>That did not open</b></div>
      </div>
      <p class="msg" id="overrmsg"></p>
      <div class="ovrow">
        <button class="btn" id="ovback" type="button">Back to the start screen</button>
      </div>
    </div>
  </div>
</div>
<script>
const $=s=>document.querySelector(s);
const COMPOUNDS=__COMPOUNDS__;
const SEP=String.fromCharCode(92);          // Windows separators, without a literal
function parts(p){const s=String(p).split(SEP).join('/'),i=s.lastIndexOf('/');
  if(i<0)return{name:s,dir:''};
  return{name:s.slice(i+1),dir:i===0?'/':s.slice(0,i)}}
function el(tag,cls,text){const n=document.createElement(tag);
  if(cls)n.className=cls; if(text!=null)n.textContent=text; return n;}

async function checkUpdate(){
  let response;
  try{response=await fetch('/api/update');}catch(e){return;}
  if(!response.ok)return;
  const available=await response.json().catch(()=>null);
  if(!available||!available.available)return;
  $('#update-copy').textContent='Version '+available.version+' is available. Sniff will '
    +'download the verified installer and open it for you.';
  $('#update').hidden=false; lock(true); $('#install-update').focus();
}
async function closeUpdate(){
  try{await fetch('/update/later',{method:'POST'});}catch(e){}
  $('#update').hidden=true;lock(false);
}

const compoundKey=s=>String(s||'').trim().replace(/\\s+/g,' ').toLowerCase();
const COMPOUND_BY_NAME=new Map(COMPOUNDS.map(c=>[compoundKey(c.name),c]));
let pendingPath=null,selectedCompounds=[],shownCompounds=[],activeCompound=0;
function compoundMatch(value){return COMPOUND_BY_NAME.get(compoundKey(value));}
function setCompoundStatus(text,isErr=false){const s=$('#compound-status');
  s.textContent=text||''; s.className='inputstatus'+(isErr?' err':'');}
function renderCompoundChips(){const box=$('#compound-chips'); box.innerHTML='';
  selectedCompounds.forEach((c,index)=>{const chip=el('span','chip');
    chip.setAttribute('role','listitem');
    chip.append(el('span',null,c.name),el('small',null,c.formula||''));
    const remove=el('button',null,'×'); remove.type='button';
    remove.setAttribute('aria-label','Remove '+c.name); remove.onclick=()=>{
      selectedCompounds.splice(index,1); renderCompoundChips(); $('#compound-input').focus();};
    chip.append(remove); box.append(chip); });
  $('#interest-continue').disabled=!selectedCompounds.length;
  $('#interest-continue').textContent=selectedCompounds.length
    ? 'Continue with '+selectedCompounds.length : 'Continue';
}
function renderCompoundSuggestions(){const input=$('#compound-input'),raw=input.value,key=compoundKey(raw);
  const exact=compoundMatch(raw); shownCompounds=key?COMPOUNDS.filter(c=>compoundKey(c.name).includes(key))
    .sort((a,b)=>(compoundKey(a.name).startsWith(key)?0:1)-(compoundKey(b.name).startsWith(key)?0:1)
      ||a.name.localeCompare(b.name)).slice(0,8):[];
  activeCompound=0; const list=$('#compound-suggestions'); list.innerHTML='';
  shownCompounds.forEach((c,index)=>{const option=el('button','suggestion'+(index===0?' active':''));
    option.type='button'; option.id='compound-option-'+index; option.setAttribute('role','option');
    option.setAttribute('aria-selected',index===0?'true':'false');
    option.append(el('span',null,c.name),el('small',null,(c.formula||'')+' · m/z '+c.mz));
    option.onmousedown=event=>{event.preventDefault(); addCompound(c);}; list.append(option); });
  list.hidden=!shownCompounds.length; input.setAttribute('aria-expanded',String(!!shownCompounds.length));
  if(shownCompounds.length)input.setAttribute('aria-activedescendant','compound-option-0');
  else input.removeAttribute('aria-activedescendant');
  input.setAttribute('aria-invalid',String(!!raw&&!exact));
  if(!raw)setCompoundStatus('');
  else if(exact)setCompoundStatus('Recognised: '+exact.name+' ('+exact.formula+').');
  else if(shownCompounds.length)setCompoundStatus('Choose a recognised suggestion.');
  else setCompoundStatus('Not recognised by the bundled PTR Library.',true);
}
function addCompound(compound){if(!compound)return false;
  if(!selectedCompounds.some(c=>compoundKey(c.name)===compoundKey(compound.name))){
    selectedCompounds.push(compound); renderCompoundChips();
  }
  const input=$('#compound-input'); input.value=''; renderCompoundSuggestions(); input.focus(); return true;
}
function addTypedCompound(useSuggestion=false){const input=$('#compound-input');
  const match=compoundMatch(input.value)||(useSuggestion?shownCompounds[activeCompound]:null);
  if(match)return addCompound(match);
  if(input.value.trim())setCompoundStatus('Not recognised by the bundled PTR Library.',true);
  input.setAttribute('aria-invalid','true'); return false;
}
function splitCompoundBlock(text){const good=[],bad=[];
  for(const row of String(text).split(/\\r?\\n/)){const line=row.trim(); if(!line)continue;
    const whole=compoundMatch(line); if(whole){good.push(whole);continue;}
    const pieces=line.split(','); let i=0;
    while(i<pieces.length){let found=null,end=i+1;
      for(let j=pieces.length;j>i;j--){const candidate=pieces.slice(i,j).join(',').trim();
        const match=compoundMatch(candidate); if(match){found=match;end=j;break;}}
      if(found)good.push(found); else if(pieces[i].trim())bad.push(pieces[i].trim());
      i=end;
    }
  }
  return {good,bad};
}
function addCompoundBlock(text){const parsed=splitCompoundBlock(text); parsed.good.forEach(addCompound);
  if(parsed.bad.length)setCompoundStatus('Rejected: '+parsed.bad.join(', ')+'. Not recognised by the bundled PTR Library.',true);
  return parsed;
}
function showInterest(path){pendingPath=path; selectedCompounds=[]; renderCompoundChips();
  const q=parts(path); $('#interest-file').textContent=q.name+(q.dir?' — '+q.dir:'');
  $('#compound-input').value=''; renderCompoundSuggestions(); $('#interest').hidden=false;
  lock(true); $('#compound-input').focus();}
function closeInterest(){pendingPath=null; $('#interest').hidden=true; $('#compound-suggestions').hidden=true;
  lock(false); setCompoundStatus('');}
function beginPendingOpen(compounds){const path=pendingPath;if(!path)return;
  $('#interest').hidden=true; lock(false); pendingPath=null; openFile(path,compounds);}

async function openFile(path,compounds){
  let r; const body={path}; if(compounds!==undefined)body.compounds_of_interest=compounds.map(c=>c.name);
  try{r=await fetch('/open',{method:'POST',headers:{'Content-Type':'application/json'},
    body:JSON.stringify(body)});}
  catch(e){return note('The app is no longer running.',true,false,true);}
  if(!r.ok){
    const b=await r.json().catch(()=>({}));
    return note(b.error||'Could not open that file.',true,false,true);
  }
  // Ask and start watching in the same breath: the server sets "loading" on its own
  // thread, so a poll that arrives first would otherwise show the old screen and
  // wait two and a half seconds before the sheet went up.
  ask=path; watching=true; resetEta(); showSheet(); tick();
}

// ---- the sheet an open runs behind -----------------------------------------
// Reading a 2 GB run takes about half a minute, which is long enough to be worth
// leaving, so the open is shown full screen with a real bar and a Cancel button
// rather than as a line of text the user has to trust.
let watching=false;                    // an open is in flight, from this page
let failed=false;                      // the sheet is showing an error, not a bar
let ask=null;                          // the path this page asked to open
let t0=0;                              // when this page started watching
let navigating=false;                   // ready can be observed by two poll turns
const reducedMotion=window.matchMedia('(prefers-reduced-motion: reduce)').matches;

function resetEta(){t0=0;}
function handoff(){
  if(navigating)return;
  navigating=true;
  watching=false;
  setBar(1);
  $('#ovstage').textContent='Ready';
  $('#cancel').disabled=true; $('#cancel').hidden=true;
  try{sessionStorage.setItem('sniff-review-entrance','1');}catch(e){}
  const go=()=>{
    if(!navigating)return;
    navigating=false;
    location.replace('/review');
  };
  if(reducedMotion)return go();
  const sheet=$('#ov'), card=sheet.querySelector('.ovcard');
  let settled=false;
  const settle=e=>{
    if(settled|| (e && e.animationName!=='handoffcard'))return;
    settled=true; card.removeEventListener('animationend',settle); go();
  };
  card.addEventListener('animationend',settle);
  sheet.classList.add('handoff');
  window.setTimeout(settle,760);
}
function lock(on){
  document.documentElement.classList.toggle('lock',on);
  const m=document.querySelector('main'); if(m)m.inert=on;
}
function showSheet(){
  $('#ovload').hidden=false; $('#overr').hidden=true;
  $('#cancel').disabled=false; $('#cancel').hidden=false;
  setBar(0); $('#oveta').textContent='';
  $('#ov').hidden=false; lock(true); $('#cancel').focus();
}
function closeSheet(){
  $('#ov').hidden=true; $('#ovload').hidden=false; $('#overr').hidden=true;
  lock(false); resetEta(); ask=null; watching=false; failed=false;
}
function failSheet(msg){
  failed=true;
  note(msg,true,false,true);           // so it is still there once the sheet is gone
  $('#ovload').hidden=true; $('#overr').hidden=false;
  $('#overrmsg').textContent=msg;
  $('#ovback').focus();
}
function setBar(p){
  const pct=Math.round(p*100);
  $('#ovfill').style.width=pct+'%';
  $('#ovbar').setAttribute('aria-valuenow',String(pct));
  $('#ovpct').textContent=pct+'%';
}
function left(sec){
  if(sec<10)return'a few seconds left';
  if(sec<55)return'~'+Math.round(sec/5)*5+' s left';
  const m=Math.round(sec/60);
  return m<=1?'about a minute left':'about '+m+' minutes left';
}
function paint(s){
  const q=parts(s.file||ask||'');
  $('#ovname').textContent=q.name||'Opening a file';
  $('#ovname').title=s.file||ask||'';
  $('#ovdir').textContent=q.dir;
  $('#ovstage').textContent=(s.stage||'Working')+((s.agent_status||'')?' — '+s.agent_status:'');
  const now=Date.now(), p=(typeof s.progress==='number')?s.progress:0;
  if(!t0)t0=now;
  // Elapsed time and the fraction it bought, which is the only estimate that needs
  // no guess about how the phases are weighted. It is rough by nature: the first
  // 11 % of the bar is the quick phases, so say "~" and leave it at that.
  const secs=(now-t0)/1000;
  let eta='';
  if(p>=0.04&&secs>=2){
    const rest=secs*(1-p)/p;
    if(isFinite(rest)&&rest>=0&&rest<1800)eta=left(rest);
  }
  $('#oveta').textContent=eta;
  setBar(p);
  $('#cancel').hidden=(s.cancellable===false);   // nothing to cancel, no button
}

let sticky=0;                                // until when the poller must leave this alone
function note(text,isErr,progress,hold){
  if(hold)sticky=Date.now()+12000; else if(!text)sticky=0;
  const box=$('#state'); box.innerHTML='';
  box.className='note'+(isErr?' err':'');
  box.hidden=!text;
  if(!text)return;
  box.append(document.createTextNode(text));
  if(progress){const bar=el('div','bar'); bar.append(el('i')); box.append(bar);}
}

function current(s){
  const box=$('#now'); box.innerHTML='';
  if(!s.file)return;
  const card=el('div','card now'),txt=el('div','txt'),q=parts(s.file);
  txt.append(el('b',null,q.name),el('span','sub',q.dir||q.name));
  txt.title=s.file;
  card.append(txt);
  if(s.export&&s.export.out){
    const out=parts(s.export.out),m=el('div','meta');
    m.append(el('div',null,'exported'),el('em',null,out.name));
    m.title=s.export.out;
    card.append(m);
  }
  const open=el('button','btn','Open the review'),close=el('button','btn sec','Close');
  open.onclick=()=>location='/review';
  close.onclick=async()=>{
    const r=await fetch('/close',{method:'POST'}).catch(()=>null);
    if(!r||!r.ok)return note('Could not close the file.',true,false,true);
    note(''); box.innerHTML='';
  };
  card.append(open,close); box.append(card);
}

async function tick(){
  let s=null;
  try{s=await (await fetch('/api/state')).json();}catch(e){return;}
  $('#quit').hidden=(s.surface!=='browser');   // a tab has no window to close
  if(s.status==='loading'){
    if(!watching){watching=true; resetEta(); showSheet();}
    paint(s);
    note('');
  }else if(watching){
    watching=false;
    if(s.status==='ready'){
      // The handoff lets the user see that the work completed; its guard also
      // handles a ready response arriving twice before navigation finishes.
      handoff();
      return;
    }
    if(s.status==='error')failSheet(s.error||'Could not open that file.');
    else{closeSheet(); note('Opening cancelled.',false,false,true);}
  }
  if(!watching&&!failed){
    if(s.status==='exporting'){
      note((s.stage||'Working')+((s.agent_status||'')?' — '+s.agent_status:''),false,true);
    }else if(s.status==='error'){
      note(s.error||'Could not open that file.',true,false,true);
    }else{
      if(Date.now()>sticky)note('');
      if(s.status==='ready') current(s); else $('#now').innerHTML='';
    }
  }
  setTimeout(tick, watching||s.status==='exporting'?900:2500);
}

$('#compound-input').oninput=renderCompoundSuggestions;
$('#compound-input').onkeydown=event=>{
  if(event.key==='ArrowDown'||event.key==='ArrowUp'){
    if(!shownCompounds.length)return; event.preventDefault();
    activeCompound=(activeCompound+(event.key==='ArrowDown'?1:-1)+shownCompounds.length)%shownCompounds.length;
    document.querySelectorAll('.suggestion').forEach((e,i)=>{
      e.classList.toggle('active',i===activeCompound); e.setAttribute('aria-selected',i===activeCompound?'true':'false');});
    $('#compound-input').setAttribute('aria-activedescendant','compound-option-'+activeCompound);
    return;
  }
  if(event.key==='Enter'){
    if(!$('#compound-input').value.trim())return; event.preventDefault(); addTypedCompound(true); return;
  }
  if(event.key==='Tab'&&!event.shiftKey&&$('#compound-input').value.trim()){
    const match=compoundMatch($('#compound-input').value)||shownCompounds[activeCompound];
    if(match){event.preventDefault();addCompound(match);} return;
  }
  if(event.key===','){
    const input=$('#compound-input'),exact=compoundMatch(input.value);
    const commaInsideName=COMPOUNDS.some(c=>compoundKey(c.name).startsWith(compoundKey(input.value)+','));
    if(exact){event.preventDefault();addCompound(exact);}
    else if(!commaInsideName){event.preventDefault();addTypedCompound(false);}
  }
};
$('#compound-input').onpaste=event=>{const text=event.clipboardData&&event.clipboardData.getData('text');
  if(text&&/[\\r\\n,]/.test(text)){event.preventDefault();addCompoundBlock(text);}};
$('#compound-input').onfocus=renderCompoundSuggestions;
$('#compound-input').onblur=()=>setTimeout(()=>{$('#compound-suggestions').hidden=true;
  $('#compound-input').setAttribute('aria-expanded','false');
  $('#compound-input').removeAttribute('aria-activedescendant');},120);
$('#interest-back').onclick=closeInterest;
$('#interest-skip').onclick=()=>beginPendingOpen([]);
$('#interest-continue').onclick=()=>{
  if($('#compound-input').value.trim()&&!addTypedCompound(false))return;
  beginPendingOpen(selectedCompounds.slice());
};

$('#install-later').onclick=closeUpdate;
$('#install-update').onclick=async()=>{
  const install=$('#install-update'),later=$('#install-later'),error=$('#update-error');
  install.disabled=true; later.disabled=true; error.hidden=true;
  install.textContent='Downloading…';
  let response=null;
  try{response=await fetch('/update',{method:'POST'});}catch(e){}
  const body=response?await response.json().catch(()=>({})):{};
  if(!response||!response.ok){
    error.textContent=body.error||'Could not download the update. Try again later.';
    error.hidden=false; install.disabled=false; later.disabled=false;
    install.textContent='Try again'; return;
  }
  $('#update-title').textContent='Installer opened';
  $('#update-copy').textContent='Close Sniff, finish the installation, then reopen Sniff '
    +'to use version '+body.version+'.';
  install.hidden=true; later.disabled=false; later.textContent='Close'; later.focus();
};

$('#browse').onclick=async()=>{
  const btn=$('#browse'); btn.disabled=true;
  note('Choose a file in the dialog that just opened on this computer.');
  let r=null;
  try{r=await fetch('/browse',{method:'POST'});}catch(e){}
  btn.disabled=false;
  const body=r?await r.json().catch(()=>({})):{};
  if(!r||!r.ok){
    return note((body&&body.error)||'Native file browsing is unavailable.',true,false,true);
  }
  if(body.cancelled)return note('');
  showInterest(body.path);
};
$('#cancel').onclick=async()=>{
  // Disable it here rather than wait for the poll to say the open is over: a second
  // click has nothing to cancel, and a button that answers twice looks like it lied.
  $('#cancel').disabled=true;
  $('#ovstage').textContent='Cancelling';
  try{await fetch('/cancel',{method:'POST'});}catch(e){}
};
$('#ovback').onclick=()=>{closeSheet(); note('');};
$('#quit').onclick=async ev=>{
  ev.preventDefault();
  const ok=await fetch('/shutdown',{method:'POST'}).then(r=>r.ok).catch(()=>false);
  note(ok?'The app has stopped. You can close this tab.':'Could not stop the app.',!ok,false,
       !ok?true:false);
};
checkUpdate();
tick();
</script></body></html>"""

# The brand is spelled once, in brand.py; the page is a template rather than an
# f-string because its CSS is full of braces.
_START_HTML = (
    _START_TEMPLATE.replace('__MARK__', brand.MARK_SVG)
    .replace(
        '__COMPOUNDS__',
        json.dumps(_compound_catalogue(), ensure_ascii=True).replace('<', '\\u003c'),
    )
    .replace('__APP_NAME__', brand.APP_NAME)
    .replace('__TAGLINE__', brand.TAGLINE)
    .replace('__PAGE_TITLE__', brand.PAGE_TITLE)
    .replace('__VERSION__', __version__)
)


def _reveal(path) -> bool:
    """Show a file in the user's own file manager. Best effort: the path is always
    printed in the UI as well, so a missing helper costs nothing."""
    target = str(Path(path))
    if sys.platform == "darwin":
        cmd = ["open", "-R", target]
    elif os.name == "nt":
        cmd = ["explorer", "/select," + target]
    else:
        cmd = ["xdg-open", str(Path(target).parent)]
    try:
        return subprocess.run(cmd, check=False, timeout=15).returncode == 0
    except (OSError, subprocess.SubprocessError):
        return False


def _recent_entries(open_path=None):
    """Build the recent-file API records, including the current-file flag."""
    open_resolved = None
    if open_path:
        try:
            open_resolved = Path(open_path).resolve()
        except OSError:
            open_resolved = None
    entries = []
    for raw in load_recent():
        p = Path(raw)
        cfg = config_path_for(raw)
        try:
            st = p.stat()
            entry = {
                "path": str(p),
                "exists": True,
                "size": st.st_size,
                "mtime": st.st_mtime,
                "config_exists": cfg.exists(),
            }
        except OSError:
            entry = {
                "path": str(p),
                "exists": False,
                "size": 0,
                "mtime": 0,
                "config_exists": cfg.exists(),
            }
        if open_resolved is not None:
            try:
                entry["is_open"] = p.resolve() == open_resolved
            except OSError:
                entry["is_open"] = str(p) == str(open_path)
        entries.append(entry)
    return entries


def _pick_file():
    """Ask the desktop for a path and return it, or ``None`` if cancelled.

    A browser hands over a file's name but never its location, so the only way to give
    the start screen a real file dialog is to ask the machine the server runs on.
    """
    if sys.platform == "darwin":
        cmd = [
            "osascript",
            "-e",
            'POSIX path of (choose file with prompt "Choose an IONICON run")',
        ]
    elif os.name == "nt":
        cmd = [
            "powershell",
            "-NoProfile",
            "-Command",
            "Add-Type -AssemblyName System.Windows.Forms;"
            " $d = New-Object System.Windows.Forms.OpenFileDialog;"
            " $d.Filter = 'IONICON runs (*.h5)|*.h5|All files (*.*)|*.*';"
            " if ($d.ShowDialog() -eq 'OK') { [Console]::Out.Write($d.FileName) }",
        ]
    else:
        for tool, extra in (
            ("zenity", ["--file-selection"]),
            (
                "kdialog",
                [
                    "--getopenfilename",
                    os.path.expanduser("~"),
                    "*.h5|IONICON runs (*.h5)",
                ],
            ),
        ):
            if shutil.which(tool):
                cmd = [tool] + extra
                break
        else:
            raise RuntimeError("Native file browsing is unavailable.")
    done = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    return done.stdout.strip() or None


def _browse():
    """Ask for a path using whichever dialog this machine offers, or ``None`` if
    cancelled.

    A desktop window owns a real dialog, and that is the one the reviewer is looking
    at. A browser tab owns none, so the machine is asked instead — its dialog can land
    behind the window, which is normal for a browser-based start screen.
    """
    window = desktop.current_window()
    if window is not None:
        try:
            return desktop.pick_file(window)
        except desktop.DesktopUnavailable:
            pass  # no dialog in the window after all; ask the machine itself
    return _pick_file()


def _background(fn, *args, **kwargs):
    """Run a long job off the request thread. The session is where the page reads the
    result or failure, so the thread only echoes it to the available diagnostic log."""

    def run():
        try:
            fn(*args, **kwargs)
        except Exception as exc:
            _log(f"sniff: {type(exc).__name__}: {exc}")

    threading.Thread(target=run, daemon=True).start()


def make_server(port=8765, agent_url=None, agent_timeout=300.0):
    """Build the app server on the first free localhost port.

    Returns ``(server, session, url)``. Kept separate from :func:`serve_app` so tests
    can drive the routes without blocking.
    """
    session = Session()
    update_lock = threading.Lock()
    update_install_lock = threading.Lock()
    update_checked = False
    update_deferred = False
    available_update = None

    def latest_update():
        nonlocal update_checked, available_update
        with update_lock:
            if update_deferred:
                return None
            if not update_checked:
                try:
                    available_update = updates.find_update(__version__)
                except (OSError, TypeError, ValueError):
                    available_update = None
                update_checked = True
            return available_update

    def defer_update():
        nonlocal update_deferred
        with update_lock:
            update_deferred = True

    class Handler(server.BaseHTTPRequestHandler):
        def log_message(self, *args):
            pass

        def _send(self, code, body=b"", ctype="application/json"):
            if not isinstance(body, bytes):
                body = json.dumps(body).encode("utf-8")
            self.send_response(code)
            self.send_header("Content-Type", ctype)
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            if body:
                self.wfile.write(body)

        def _same_origin(self):
            expected_host = f"127.0.0.1:{self.server.server_address[1]}"
            if self.headers.get("Host") != expected_host:
                return False
            if self.headers.get("Sec-Fetch-Site") == "cross-site":
                return False
            origin = self.headers.get("Origin")
            if not origin:
                return True
            try:
                parsed = urlparse(origin)
                port = parsed.port
            except ValueError:
                return False
            return (
                parsed.scheme == "http"
                and parsed.hostname == "127.0.0.1"
                and port == self.server.server_address[1]
            )

        def do_GET(self):
            route = urlparse(self.path)
            if route.path in ("/", "/index.html"):
                self._send(200, _START_HTML.encode("utf-8"), "text/html; charset=utf-8")
            elif route.path == "/review":
                if not session.payload:
                    self._send(404, {"error": "no file is open"})
                    return
                previous_page = session.current_review_page()
                session.wait_for_page_close_save(previous_page)
                try:
                    payload, config_path, page_token = session.prepare_review_page()
                except (OSError, RuntimeError, TypeError, ValueError) as exc:
                    self._send(409, {"error": f"could not refresh the review: {exc}"})
                    return
                html = viz.render_html(
                    payload,
                    config_path=config_path,
                    mode="app",
                    page_token=page_token,
                    auto_tour=auto_tour_pending(),
                )
                self._send(200, html.encode("utf-8"), "text/html; charset=utf-8")
            elif route.path == "/api/update":
                available = latest_update()
                if available is None:
                    self._send(200, {"available": False})
                else:
                    self._send(
                        200,
                        {
                            "available": True,
                            "version": available.version,
                            "release_url": available.release_url,
                        },
                    )
            elif route.path == "/api/state":
                self._send(200, session.state())
            elif route.path == "/status":
                self._send(200, session.status_payload())
            elif route.path == "/api/recent":
                self._send(200, _recent_entries(session.path))
            elif route.path == "/peak-preview":
                if not session.path:
                    self._send(404, {"error": "no file is open"})
                    return
                query = parse_qs(route.query)
                try:
                    lo = float(query.get("lo", [""])[0])
                    hi = float(query.get("hi", [""])[0])
                    self._send(200, session.peak_preview(lo, hi))
                except (OSError, RuntimeError, TypeError, ValueError) as exc:
                    self._send(500, {"error": str(exc)})
            elif route.path == "/spectrum":
                if not session.path:
                    self._send(404, {"error": "no file is open"})
                    return
                q = parse_qs(route.query)
                try:
                    lo = int(q.get("lo", ["1"])[0])
                    hi = int(q.get("hi", ["1"])[0])
                    self._send(200, interval_spectrum(session.path, lo, hi))
                except (OSError, RuntimeError, TypeError, ValueError) as exc:
                    self._send(500, {"error": str(exc)})
            else:
                self._send(404, {"error": "not found"})

        def do_POST(self):
            route = urlparse(self.path)
            if route.path in ("/update", "/update/later") and not self._same_origin():
                self._send(403, {"error": "cross-origin update request refused"})
                return
            n = int(self.headers.get("Content-Length", 0) or 0)
            try:
                body = json.loads(self.rfile.read(n) if n else b"{}")
            except (UnicodeError, json.JSONDecodeError, ValueError):
                self._send(400, {"error": "invalid JSON"})
                return
            if route.path == "/update":
                if not update_install_lock.acquire(blocking=False):
                    self._send(409, {"error": "an update is already being opened"})
                    return
                try:
                    available = latest_update()
                    if available is None:
                        self._send(409, {"error": "no compatible update is available"})
                        return
                    try:
                        updates.install_update(available)
                    except (OSError, TypeError, ValueError) as exc:
                        self._send(503, {"error": f"could not open the update: {exc}"})
                        return
                    defer_update()
                    self._send(200, {"ok": True, "version": available.version})
                finally:
                    update_install_lock.release()
            elif route.path == "/update/later":
                defer_update()
                self._send(200, {"ok": True})
            elif route.path == "/open":
                target = str(body.get("path") or "").strip()
                if not target:
                    self._send(400, {"error": "no path given"})
                    return
                if not Path(target).expanduser().is_file():
                    self._send(404, {"error": f"no such file: {target}"})
                    return
                if body.get("adapt_peaks") is not None and config_path_for(target).exists():
                    self._send(
                        409,
                        {
                            "error": "the target already has a saved review; "
                            "open it normally instead"
                        },
                    )
                    return
                try:
                    compounds_of_interest = (
                        _normalise_compounds_of_interest(body["compounds_of_interest"])
                        if "compounds_of_interest" in body
                        else None
                    )
                except ValueError as exc:
                    self._send(400, {"error": str(exc)})
                    return
                template_peaks = body.get("adapt_peaks")
                if template_peaks is not None and (
                    not isinstance(template_peaks, list)
                    or any(
                        not isinstance(peak, dict) or "mz" not in peak
                        for peak in template_peaks
                    )
                ):
                    self._send(400, {"error": "adapt_peaks must be a peak list"})
                    return
                if not session.reserve_open():
                    self._send(409, {"error": "the app is busy with the current file"})
                    return
                _background(
                    session.open,
                    target,
                    agent_url=agent_url,
                    agent_timeout=agent_timeout,
                    reserved=True,
                    compounds_of_interest=compounds_of_interest,
                    template_peaks=template_peaks,
                )
                self._send(202, {"ok": True})
            elif route.path == "/save":
                query = parse_qs(route.query)
                closing = query.get("closing", ["0"])[0] == "1"
                page = query.get("page", [None])[0]
                try:
                    if not session.accepts_review_page(page):
                        self._send(409, {"error": "this review page is no longer current"})
                        return
                    if not _valid_config(body, require_mass_axis=True):
                        self._send(
                            400,
                            {
                                "error": "config must contain peaks or ranges and the "
                                "current mass-axis marker"
                            },
                        )
                        return
                    raw_version = query.get("version", [None])[0]
                    try:
                        version = int(raw_version) if raw_version is not None else None
                    except ValueError:
                        self._send(400, {"error": "save version must be an integer"})
                        return
                    try:
                        saved = session.save_config(body, version=version, page=page)
                    except RuntimeError as exc:
                        self._send(409, {"error": str(exc)})
                        return
                    except ValueError as exc:
                        self._send(400, {"error": str(exc)})
                        return
                    except OSError as exc:
                        # Say so rather than let the request thread die: the page needs a
                        # reason, and the config on disk is still the last good one.
                        self._send(500, {"error": f"could not write the config: {exc}"})
                        return
                    self._send(200, {"ok": True, "saved": saved})
                finally:
                    if closing:
                        session.finish_close_save(page)
            elif route.path == "/export":
                query = parse_qs(route.query)
                page = query.get("page", [None])[0]
                if not session.accepts_review_page(page):
                    self._send(409, {"error": "this review page is no longer current"})
                    return
                if not session.path:
                    self._send(409, {"error": "no file is open"})
                    return
                # The page posts the config it is showing. Autosave is debounced, so
                # exporting whatever happens to be on disk can describe an earlier
                # state than the one the reviewer just looked at. An omitted body is
                # retained for API callers that simply request the last saved export.
                if body and not _valid_config(body, require_mass_axis=True):
                    self._send(400, {"error": "export config has an invalid mass-axis marker"})
                    return
                raw_version = query.get("version", [None])[0]
                try:
                    version = int(raw_version) if raw_version is not None else None
                except ValueError:
                    self._send(400, {"error": "save version must be an integer"})
                    return
                if not session.reserve_export():
                    self._send(409, {"error": "an export is already running"})
                    return
                try:
                    path, config = session.prepare_export(
                        body or None, version=version, page=page
                    )
                except RuntimeError as exc:
                    session.cancel_export_reservation()
                    self._send(409, {"error": str(exc)})
                    return
                except ValueError as exc:
                    session.cancel_export_reservation()
                    self._send(400, {"error": str(exc)})
                    return
                except OSError as exc:
                    session.cancel_export_reservation()
                    self._send(500, {"error": f"could not write the config: {exc}"})
                    return
                _background(
                    session.export,
                    path=path,
                    config=config,
                    reserved=True,
                )
                self._send(202, {"ok": True})
            elif route.path == "/cancel":
                # Cancelling an open that is not running is not an event: the page's
                # Cancel button and a poll can cross, and the answer must not depend
                # on which one arrived first.
                session.cancel()
                self._send(200, {"ok": True})
            elif route.path == "/close":
                if not session.reserve_close():
                    self._send(409, {"error": "wait for the current work to finish"})
                    return
                try:
                    forget_active()
                    session.close()
                except OSError as exc:
                    self._send(500, {"error": f"could not clear the saved session: {exc}"})
                    return
                finally:
                    session.finish_close()
                self._send(200, {"ok": True})
            elif route.path == "/reveal":
                last = (session.export_result or {}).get("out")
                if not last:
                    self._send(409, {"error": "nothing has been exported yet"})
                    return
                self._send(200, {"ok": _reveal(last), "path": last})
            elif route.path == "/ack":
                self._send(200, {"ok": True})
            elif route.path == "/onboarding":
                try:
                    remember_tour_seen()
                except OSError as exc:
                    self._send(
                        500, {"error": f"could not save onboarding state: {exc}"}
                    )
                    return
                self._send(200, {"ok": True})
            elif route.path == "/browse":
                try:
                    picked = _browse()
                except (
                    OSError,
                    RuntimeError,
                    subprocess.SubprocessError,
                    desktop.DesktopUnavailable,
                ) as exc:
                    self._send(501, {"error": str(exc) or "no file dialog here"})
                    return
                self._send(200, {"path": picked} if picked else {"cancelled": True})
            elif route.path == "/shutdown":
                # One direction each, so the two can never chase one another: the page
                # stops the session and closes the window, while closing the window
                # only ever stops the session. Nothing here waits for the other.
                stop_the_app(session)
                self._send(200, {"ok": True})
            else:
                self._send(404, {"error": "not found"})

    httpd = None
    # port=0 asks the OS for a free port, which is what a test wants: probing upward
    # from a guessed number can land on a server a previous test has not released.
    candidates = [0] if port == 0 else range(port, port + 20)
    for candidate in candidates:
        try:
            # ThreadingHTTPServer, not a bare TCPServer: it sets allow_reuse_address,
            # so a restart lands back on the same port instead of drifting.
            httpd = server.ThreadingHTTPServer(("127.0.0.1", candidate), Handler)
            break
        except OSError:
            continue
    if httpd is None:
        raise OSError("no free port found for the app server")
    httpd.daemon_threads = True
    actual = httpd.server_address[1]
    return httpd, session, f"http://127.0.0.1:{actual}/"


# How the running app reached the user: a window of its own, or a browser tab. Set by
# serve_app, read by /api/state, and the difference between "the app opened" and "the
# app opened a window", which no log line a windowed bundle can write would ever show.
_surface = "browser"


def surface() -> str:
    return _surface

def _log(text):
    """Report progress on stderr, and to a log file too when there is no console.

    A windowed frozen bundle has nowhere to print, so its URL and its tracebacks would
    otherwise be lost. PyInstaller sets stderr to ``None`` in Windows windowed mode.
    """
    if sys.stderr is not None:
        print(text, file=sys.stderr, flush=True)
    if not getattr(sys, "frozen", False):
        return
    try:
        path = _recent_path().parent / "log.txt"
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("a", encoding="utf-8") as handle:
            handle.write(text + "\n")
    except OSError:
        pass  # a missing log file must not stop the app


def _install_quit_handlers(on_quit):
    """Make Ctrl-C, ``pkill`` and a bundle's quit all stop the server the same way."""
    for name in ("SIGINT", "SIGTERM"):
        sig = getattr(signal, name, None)
        if sig is None:
            continue
        try:
            signal.signal(sig, lambda *_: on_quit())
        except (OSError, ValueError):
            pass  # not the main thread, or this platform has no such signal


def stop_the_app(session):
    """Stop the server, and take the desktop window with it if there is one.

    This is the one quit path, used by Ctrl-C, ``pkill``, a bundle's quit and the page's
    Stop button. In window mode the native loop owns the main thread, so a quit that
    only set the flag would leave a window standing over a dead server. With no window
    there is nothing to close, and :func:`desktop.close_window` says so by returning
    False; closing it is never a second way into the session, because the window's own
    close handler only ever sets the same flag.
    """
    session.stop.set()
    try:
        desktop.close_window()
    except desktop.DesktopUnavailable as exc:
        # The server is stopping either way; a window left standing is worth a line on
        # stderr, not a failed request or a swallowed quit.
        _log(f"sniff: {exc}; close that window yourself to get rid of it")


def serve_app(
    port=8765,
    open_browser=True,
    agent_url=None,
    agent_timeout=300.0,
    initial=None,
    window=False,
):
    """Serve the app until interrupted. Nothing here closes on its own: an export, a
    closed tab or a closed file all leave the server up.

    ``window=True`` asks for a desktop window instead of a browser tab: the native loop
    takes the main thread, which is why the server runs on a daemon thread behind it.
    A window that cannot start is never fatal — one line on stderr, then the browser
    route, because losing the session over losing the window is the worse trade.
    """
    httpd, session, url = make_server(
        port=port, agent_url=agent_url, agent_timeout=agent_timeout
    )
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    _install_quit_handlers(lambda: stop_the_app(session))

    _log(f"sniff: app running at {url}")
    _log("sniff: a large run takes 30-90 s to open; the app stays up between files.")
    # A fresh app launch always starts at the opening screen. Selecting a file still
    # discovers and reuses its same-basename config, but closing the app is not a request
    # to reopen the last run automatically.
    if initial is None:
        try:
            forget_active()
        except OSError:
            pass
    resume = initial
    if resume:
        if Path(resume).expanduser().is_file():
            if session.reserve_open():
                _background(
                    session.open,
                    resume,
                    agent_url=agent_url,
                    agent_timeout=agent_timeout,
                    reserved=True,
                )
        else:
            _log(f"sniff: cannot resume missing file: {resume}")
    global _surface
    # Decided fresh on every serve: a second call in the same process — a test, or a
    # script that restarts the app — must not inherit the previous run's surface.
    _surface = "browser"
    in_window = False
    if window:
        try:
            # Said before the call, because run_window does not return for the life of
            # the window: recorded afterwards, /api/state would answer "browser" at
            # every moment a window was actually on screen. run_window raises before it
            # blocks when there is no window to make, so this is not a promise it cannot
            # keep. Closing the window is the one thing it may do to the session, and it
            # is the reverse of the /shutdown route, which stops the session first.
            _surface = "window"
            desktop.run_window(url, on_close=session.stop.set)
            in_window = True
        except desktop.DesktopUnavailable as exc:
            _surface = "browser"
            _log(f"sniff: no desktop window ({exc}); opening a browser instead")
    if not in_window and open_browser:
        try:
            if not webbrowser.open(url):
                raise OSError("no browser answered")
        except (OSError, webbrowser.Error) as exc:
            # In a double-clicked bundle this line is the only way the user learns the
            # server is up, so it has to carry the address.
            _log(f"sniff: could not open a browser ({exc}); open {url} in one yourself")
    session.stop.wait()
    # A native window's close event can beat the browser engine's sendBeacon onto the
    # request thread. Keep the localhost endpoint alive briefly for that final snapshot;
    # an ordinary browser tab leaves the server running and never pays this wait.
    if in_window:
        session.wait_for_close_save()
    session.close()
    httpd.shutdown()
    httpd.server_close()
