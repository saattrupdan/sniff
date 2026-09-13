"""Maintainer-only NIST Chemistry WebBook crawling and parsing.

Sniff never imports this module at runtime. The generated package catalogue contains
normalised metadata only: no WebBook HTML or spectra and no PTR-MS evidence.
"""

from __future__ import annotations

import hashlib
import html
import json
import logging
import re
import sqlite3
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from html.parser import HTMLParser
from pathlib import Path

from . import __version__, formula_id, isotopes

logger = logging.getLogger(__name__)

BASE_URL = "https://webbook.nist.gov/cgi/cbook.cgi"
CACHE_PATH = Path.home() / ".sniff" / "nist-webbook.sqlite3"
SCHEDULE_PATH = Path.home() / ".sniff" / "nist-webbook-schedule.sqlite3"
CACHE_VERSION = 1
PARSER_VERSION = 1
CRAWL_DELAY_S = 5.0
REQUEST_TIMEOUT_S = 5.0
MAX_RESPONSE_BYTES = 2 * 1024 * 1024
SUCCESS_TTL_S = 180 * 24 * 60 * 60
EMPTY_TTL_S = 14 * 24 * 60 * 60
MAX_MASS_DETAILS = 4
_NIST_ID = re.compile(r"^[A-Za-z][A-Za-z0-9]+$")
_FORMULA_TEXT = re.compile(r"^[A-Z][A-Za-z0-9]*$")
_NON_DEFAULT_NAME = re.compile(
    r"\b(?:anion|cation|dianion|dication|radical|excited|transition state|"
    r"triplet|zwitterion)\b",
    re.IGNORECASE,
)
_ISOTOPE_NAME = re.compile(
    r"(?:-d\d*\b|\bD\d|deuter(?:ated|ium)|trit(?:iated|ium)|carbon-1[34]|"
    r"nitrogen-15|oxygen-1[78])",
    re.IGNORECASE,
)
_CONFORMER_LABEL = re.compile(r"^[ct]-[A-Z0-9]")
_FORMULA_LABEL = re.compile(r"(?:[A-Z][a-z]?\d*){2,}")
_POLYMER_NAME = re.compile(
    r"(?:^poly(?:mer)?(?:\b|\()|^poly[a-z]|\b(?:homo|co)?polymer\b)",
    re.IGNORECASE,
)
_STRUCTURAL_NOTATION_NAME = re.compile(r"(?:=|:\s*$|[A-Z][a-z]?\d*\.)")
_TRANSIENT_NAMES = {"methylene", "hcoh (hydroxymethylene)"}


class WebBookError(RuntimeError):
    """A recoverable WebBook request, response, or parsing failure."""


class _SearchParser(HTMLParser):
    def __init__(self, kind):
        super().__init__(convert_charrefs=True)
        self.kind = kind
        self.records = []
        self._ol_depth = 0
        self._record = None
        self._in_anchor = False
        self._after_anchor = False
        self._in_strong = False

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "ol":
            self._ol_depth += 1
        elif tag == "li" and self._ol_depth:
            self._record = {
                "name_parts": [],
                "formula_parts": [],
                "weight_parts": [],
            }
            self._after_anchor = False
        elif tag == "a" and self._record is not None:
            nist_id = _id_from_href(attrs.get("href"))
            if nist_id:
                self._record["nist_id"] = nist_id
                self._in_anchor = True
        elif tag == "strong" and self._record is not None:
            self._in_strong = True
        elif tag == "img" and self._record is not None:
            alternate = "".join(str(attrs.get("alt") or "").split())
            if _FORMULA_TEXT.fullmatch(alternate):
                self._record["image_formula"] = alternate

    def handle_endtag(self, tag):
        if tag == "a" and self._in_anchor:
            self._in_anchor = False
            self._after_anchor = True
        elif tag == "strong":
            self._in_strong = False
        elif tag == "li" and self._record is not None:
            record = self._finish_record(self._record)
            if record is not None:
                self.records.append(record)
            self._record = None
            self._after_anchor = False
        elif tag == "ol" and self._ol_depth:
            self._ol_depth -= 1

    def handle_data(self, data):
        if self._record is None:
            return
        if self._in_anchor:
            self._record["name_parts"].append(data)
        elif self._in_strong:
            self._record["weight_parts"].append(data)
        elif self._after_anchor:
            self._record["formula_parts"].append(data)

    def _finish_record(self, raw):
        nist_id = raw.get("nist_id")
        name = " ".join("".join(raw["name_parts"]).split())
        if not nist_id or not name:
            return None
        record = {
            "nist_id": nist_id,
            "name": name,
            "url": _species_url(nist_id),
        }
        if self.kind == "formula":
            formula = raw.get("image_formula")
            if not formula:
                trailing = "".join(raw["formula_parts"])
                match = re.search(r"\(([^()]*)\)", trailing)
                formula = "".join(match.group(1).split()) if match else ""
            if formula:
                record["formula"] = formula
        else:
            weight_text = "".join(raw["weight_parts"]).replace("\xa0", " ").strip()
            try:
                record["search_mass"] = float(weight_text)
            except ValueError:
                return None
        return record


class _DetailParser(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.name_parts = []
        self.items = []
        self._in_name = False
        self._top_seen = False
        self._metadata_ul = False
        self._ul_depth = 0
        self._item = None

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "h1" and attrs.get("id") == "Top":
            self._in_name = True
        elif tag == "ul" and self._top_seen and not self._metadata_ul:
            self._metadata_ul = True
            self._ul_depth = 1
        elif tag == "ul" and self._metadata_ul:
            self._ul_depth += 1
        elif tag == "li" and self._metadata_ul and self._ul_depth == 1:
            self._item = []

    def handle_endtag(self, tag):
        if tag == "h1" and self._in_name:
            self._in_name = False
            self._top_seen = True
        elif tag == "li" and self._item is not None:
            text = " ".join("".join(self._item).split())
            if text:
                self.items.append(text)
            self._item = None
        elif tag == "ul" and self._metadata_ul:
            self._ul_depth -= 1
            if self._ul_depth == 0:
                self._metadata_ul = False

    def handle_data(self, data):
        if self._in_name:
            self.name_parts.append(data)
        if self._item is not None:
            self._item.append(data)

    def record(self, nist_id):
        name = " ".join("".join(self.name_parts).split())
        fields = {}
        for item in self.items:
            for label, key in (
                ("Formula:", "formula"),
                ("Molecular weight:", "molecular_weight"),
                ("IUPAC Standard InChI:", "inchi"),
                ("IUPAC Standard InChIKey:", "inchi_key"),
                ("CAS Registry Number:", "cas"),
            ):
                if item.startswith(label):
                    value = item[len(label) :].strip()
                    if key == "inchi":
                        value = value.removesuffix(" Copy").strip()
                    elif key == "inchi_key":
                        value = value.removesuffix(" Copy").strip()
                    fields[key] = value
                    break
        if not name or not fields.get("formula"):
            raise WebBookError("the WebBook species page had no name or formula")
        return {
            "nist_id": nist_id,
            "name": name,
            "url": _species_url(nist_id),
            **fields,
        }


class WebBookClient:
    """Look up WebBook proposals without making analysis depend on the network."""

    def __init__(
        self,
        cache_path=None,
        *,
        schedule_path=None,
        opener=None,
        clock=None,
        sleeper=None,
        crawl_delay=CRAWL_DELAY_S,
        timeout=REQUEST_TIMEOUT_S,
        max_mass_details=MAX_MASS_DETAILS,
    ):
        self.cache_path = Path(cache_path) if cache_path is not None else CACHE_PATH
        self.schedule_path = (
            Path(schedule_path) if schedule_path is not None else SCHEDULE_PATH
        )
        self._opener = opener or urllib.request.urlopen
        self._clock = clock or time.time
        self._sleep = sleeper or time.sleep
        self.crawl_delay = float(crawl_delay)
        self.timeout = float(timeout)
        self.max_mass_details = int(max_mass_details)
        self._rate_lock = threading.Lock()
        self._query_lock = threading.Lock()
        self._lookup_lock = threading.Lock()
        self._memory_cache = {}

    def enrich_peak(
        self,
        mz,
        *,
        drift=1.0,
        candidates=None,
        obs_ratios=None,
        compounds_of_interest=None,
        tol_mDa=12.0,
    ):
        """Return external proposals or validated formula candidates for one peak."""
        candidates = [dict(candidate) for candidate in candidates or []]
        with self._query_lock:
            unnamed = next(
                (
                    candidate
                    for candidate in candidates
                    if candidate.get("formula")
                    and not candidate.get("preferred_name")
                    and not candidate.get("name")
                ),
                None,
            )
            if unnamed is not None:
                lookup = self.lookup_formula(unnamed["formula"])
                enriched = _attach_species(candidates, lookup["species"])
                return {
                    "status": lookup["status"],
                    "mode": "formula",
                    "candidates": enriched,
                    "excluded": lookup["excluded"],
                }
            if candidates:
                return {
                    "status": "not-needed",
                    "mode": "none",
                    "candidates": candidates,
                    "excluded": 0,
                }
            return self._rescue_mass(
                mz=float(mz),
                drift=float(drift),
                obs_ratios=obs_ratios,
                compounds_of_interest=compounds_of_interest,
                tol_mDa=float(tol_mDa),
            )

    def review_lookup(
        self,
        mz,
        *,
        formula=None,
        drift=1.0,
        obs_ratios=None,
        compounds_of_interest=None,
        tol_mDa=12.0,
    ):
        """Return browser-ready WebBook enrichment for one selected peak."""
        if formula:
            lookup = self.lookup_formula(formula)
            return {
                "status": lookup["status"],
                "mode": "formula",
                "formula": _canonical_formula(formula),
                "species": lookup["species"],
                "excluded": lookup["excluded"],
            }
        result = self.enrich_peak(
            mz,
            drift=drift,
            candidates=[],
            obs_ratios=obs_ratios,
            compounds_of_interest=compounds_of_interest,
            tol_mDa=tol_mDa,
        )
        for candidate in result["candidates"]:
            candidate["isotope_model"] = isotopes.formula_isotope_model(
                candidate["formula"]
            )
        return result

    def lookup_formula(self, formula):
        """Return ordinary neutral species registered for one exact formula."""
        canonical = _canonical_formula(formula)
        key = f"formula:{canonical}"
        url = _query_url(Formula=canonical, NoIon="on", Units="SI")
        payload, status = self._lookup(
            key=key,
            url=url,
            parser=lambda page: _parse_search(page, kind="formula"),
        )
        species = []
        excluded = 0
        for record in payload:
            try:
                record_formula = _canonical_formula(record.get("formula"))
            except ValueError:
                excluded += 1
                continue
            if record_formula != canonical or not _ordinary_name(record["name"]):
                excluded += 1
                continue
            species.append({**record, "formula": record_formula})
        species.sort(key=lambda item: (item["name"].casefold(), item["nist_id"]))
        return {"status": status, "species": species, "excluded": excluded}

    def _rescue_mass(
        self,
        *,
        mz,
        drift,
        obs_ratios,
        compounds_of_interest,
        tol_mDa,
    ):
        neutral = mz / drift - formula_id.PROTON
        tol_da = tol_mDa / 1000.0
        low = neutral - tol_da
        high = neutral + tol_da
        key = f"mass:{low:.6f},{high:.6f}"
        url = _query_url(Value=f"{low:.6f},{high:.6f}", VType="MW", Units="SI")
        search, search_status = self._lookup(
            key=key,
            url=url,
            parser=lambda page: _parse_search(page, kind="mass"),
        )
        eligible = [record for record in search if _ordinary_name(record["name"])]
        eligible.sort(
            key=lambda item: (
                abs(float(item["search_mass"]) - neutral),
                item["name"].casefold(),
                item["nist_id"],
            )
        )
        details = []
        statuses = [search_status]
        for record in eligible[: self.max_mass_details]:
            try:
                detail, detail_status = self.lookup_species(record["nist_id"])
            except (ValueError, WebBookError):
                continue
            statuses.append(detail_status)
            if (
                not _ordinary_name(detail["name"])
                or not detail.get("cas")
                or not detail.get("inchi")
            ):
                continue
            try:
                canonical = _canonical_formula(detail["formula"])
                counts = isotopes.parse_formula(canonical)
            except ValueError:
                continue
            if not formula_id._plausible(counts):
                continue
            exact_mass = formula_id.formula_mass(counts)
            if abs(exact_mass - neutral) > tol_da:
                continue
            detail["formula"] = canonical
            details.append(detail)
        formulas = sorted({detail["formula"] for detail in details})
        candidates = formula_id.score_peak(
            mz,
            drift,
            obs_ratios=obs_ratios,
            tol_mDa=tol_mDa,
            compounds_of_interest=compounds_of_interest,
            extra_formulas=formulas,
            enumerate_candidates=False,
        )
        candidates = _attach_species(candidates, details)
        excluded = max(0, len(search) - len(details))
        return {
            "status": _combined_status(statuses),
            "mode": "mass",
            "neutral_mass": round(neutral, 6),
            "candidates": candidates,
            "excluded": excluded,
        }

    def lookup_species(self, nist_id):
        """Return validated general metadata for one WebBook species ID."""
        if not isinstance(nist_id, str) or not _NIST_ID.fullmatch(nist_id):
            raise ValueError("invalid NIST WebBook species identifier")
        key = f"species:{nist_id}"
        url = _query_url(ID=nist_id, Units="SI")
        payload, status = self._lookup(
            key=key,
            url=url,
            parser=lambda page: _parse_detail(page, nist_id=nist_id),
        )
        if not isinstance(payload, dict):
            raise WebBookError("the WebBook species details were unavailable")
        return payload, status

    def _lookup(self, *, key, url, parser):
        cached = self._cache_get(key)
        now = time.time()
        if cached is not None:
            age = max(0.0, now - cached["fetched_at"])
            ttl = EMPTY_TTL_S if cached["empty"] else SUCCESS_TTL_S
            if age <= ttl:
                return cached["payload"], "cache"
        try:
            with self._lookup_lock:
                refreshed = self._cache_get(key, refresh=True)
                if _fresh_cache_entry(refreshed, now=time.time()):
                    return refreshed["payload"], "cache"
                self._reserve_request_slot()
                refreshed = self._cache_get(key, refresh=True)
                if _fresh_cache_entry(refreshed, now=time.time()):
                    return refreshed["payload"], "cache"
                page = self._fetch(url)
                payload = parser(page)
                self._cache_put(
                    key=key,
                    payload=payload,
                    empty=not bool(payload),
                    response_sha=hashlib.sha256(page).hexdigest(),
                    fetched_at=now,
                )
            return payload, "live"
        except (OSError, ValueError, WebBookError, urllib.error.URLError) as exc:
            logger.info("NIST WebBook lookup failed for %s: %s", key, exc)
            if cached is not None:
                return cached["payload"], "stale"
            return [], "unavailable"

    def _reserve_request_slot(self):
        with self._rate_lock:
            now = self._clock()
            try:
                with self._connect_schedule() as connection:
                    connection.execute("BEGIN IMMEDIATE")
                    row = connection.execute(
                        "SELECT value FROM webbook_meta WHERE key = 'next_request_at'"
                    ).fetchone()
                    previous = float(row[0]) if row is not None else now
                    reserved = max(now, previous)
                    connection.execute(
                        """
                        INSERT INTO webbook_meta (key, value)
                        VALUES ('next_request_at', ?)
                        ON CONFLICT(key) DO UPDATE SET value = excluded.value
                        """,
                        (reserved + self.crawl_delay,),
                    )
            except (OSError, sqlite3.Error, TypeError, ValueError) as exc:
                raise WebBookError(
                    "the shared WebBook request schedule is unavailable"
                ) from exc
            delay = reserved - now
        if delay > 0:
            self._sleep(delay)

    def _fetch(self, url):
        parsed = urllib.parse.urlparse(url)
        if parsed.scheme != "https" or parsed.netloc != "webbook.nist.gov":
            raise WebBookError("refusing a non-WebBook URL")
        request = urllib.request.Request(
            url,
            headers={
                "User-Agent": (
                    f"Sniff/{__version__} catalogue maintainer crawl "
                    "(+https://github.com/saattrupdan/sniff)"
                )
            },
        )
        with self._opener(request, timeout=self.timeout) as response:
            body = response.read(MAX_RESPONSE_BYTES + 1)
        if len(body) > MAX_RESPONSE_BYTES:
            raise WebBookError("the WebBook response exceeded the size limit")
        return body

    def _connect_schedule(self):
        self.schedule_path.parent.mkdir(parents=True, exist_ok=True)
        connection = sqlite3.connect(str(self.schedule_path), timeout=10.0)
        connection.execute(
            """
            CREATE TABLE IF NOT EXISTS webbook_meta (
                key TEXT PRIMARY KEY,
                value REAL NOT NULL
            )
            """
        )
        return connection

    def _connect(self):
        self.cache_path.parent.mkdir(parents=True, exist_ok=True)
        connection = sqlite3.connect(str(self.cache_path), timeout=2.0)
        connection.execute(
            """
            CREATE TABLE IF NOT EXISTS webbook_cache (
                key TEXT PRIMARY KEY,
                payload TEXT NOT NULL,
                fetched_at REAL NOT NULL,
                empty INTEGER NOT NULL,
                response_sha TEXT NOT NULL,
                cache_version INTEGER NOT NULL,
                parser_version INTEGER NOT NULL
            )
            """
        )
        return connection

    def _cache_get(self, key, refresh=False):
        if not refresh and key in self._memory_cache:
            return self._memory_cache[key]
        try:
            with self._connect() as connection:
                row = connection.execute(
                    """
                    SELECT payload, fetched_at, empty, response_sha
                    FROM webbook_cache
                    WHERE key = ? AND cache_version = ? AND parser_version = ?
                    """,
                    (key, CACHE_VERSION, PARSER_VERSION),
                ).fetchone()
        except (OSError, sqlite3.Error):
            return None
        if row is None:
            return None
        try:
            entry = {
                "payload": json.loads(row[0]),
                "fetched_at": float(row[1]),
                "empty": bool(row[2]),
                "response_sha": row[3],
            }
        except (TypeError, ValueError, json.JSONDecodeError):
            return None
        self._memory_cache[key] = entry
        return entry

    def _cache_put(self, *, key, payload, empty, response_sha, fetched_at):
        entry = {
            "payload": payload,
            "fetched_at": float(fetched_at),
            "empty": bool(empty),
            "response_sha": response_sha,
        }
        self._memory_cache[key] = entry
        try:
            with self._connect() as connection:
                connection.execute(
                    """
                    INSERT INTO webbook_cache (
                        key, payload, fetched_at, empty, response_sha,
                        cache_version, parser_version
                    ) VALUES (?, ?, ?, ?, ?, ?, ?)
                    ON CONFLICT(key) DO UPDATE SET
                        payload = excluded.payload,
                        fetched_at = excluded.fetched_at,
                        empty = excluded.empty,
                        response_sha = excluded.response_sha,
                        cache_version = excluded.cache_version,
                        parser_version = excluded.parser_version
                    """,
                    (
                        key,
                        json.dumps(payload, sort_keys=True, separators=(",", ":")),
                        float(fetched_at),
                        int(bool(empty)),
                        response_sha,
                        CACHE_VERSION,
                        PARSER_VERSION,
                    ),
                )
        except (OSError, sqlite3.Error):
            logger.info("NIST WebBook cache is unavailable; using memory only")


def _parse_search(page, *, kind):
    text = page.decode("utf-8", "replace")
    parser = _SearchParser(kind=kind)
    parser.feed(text)
    if parser.records:
        return parser.records
    if (
        "Chemical Formula Not Found" in text
        or "No matching species were found" in text
        or "0 matching species" in text
    ):
        return []
    raise WebBookError("the WebBook search page had an unrecognised structure")


def _parse_detail(page, *, nist_id):
    parser = _DetailParser()
    parser.feed(page.decode("utf-8", "replace"))
    return parser.record(nist_id=nist_id)


def _query_url(**parameters):
    return BASE_URL + "?" + urllib.parse.urlencode(parameters)


def _id_from_href(href):
    if not href:
        return None
    parsed = urllib.parse.urlparse(html.unescape(href))
    if parsed.scheme and (
        parsed.scheme != "https" or parsed.netloc != "webbook.nist.gov"
    ):
        return None
    if parsed.path != "/cgi/cbook.cgi":
        return None
    nist_id = urllib.parse.parse_qs(parsed.query).get("ID", [None])[0]
    return nist_id if nist_id and _NIST_ID.fullmatch(nist_id) else None


def _species_url(nist_id):
    return _query_url(ID=nist_id, Units="SI")


def _canonical_formula(formula):
    counts = isotopes.parse_formula(str(formula or ""))
    return formula_id.formula_str(counts)


def _ordinary_name(name):
    value = " ".join(str(name or "").split())
    return bool(
        value
        and not _NON_DEFAULT_NAME.search(value)
        and not _ISOTOPE_NAME.search(value)
        and not _CONFORMER_LABEL.search(value)
        and not _FORMULA_LABEL.fullmatch(value)
        and not _POLYMER_NAME.search(value)
        and not _STRUCTURAL_NOTATION_NAME.search(value)
        and value.casefold() not in _TRANSIENT_NAMES
    )


def _attach_species(candidates, species):
    by_formula = {}
    for record in species:
        formula = record.get("formula")
        if formula:
            by_formula.setdefault(formula.upper(), []).append(record)
    enriched = []
    for candidate in candidates:
        candidate = dict(candidate)
        proposals = by_formula.get(str(candidate.get("formula") or "").upper(), [])
        if proposals:
            unique = {}
            for proposal in proposals:
                unique[proposal["nist_id"]] = proposal
            candidate["nist_webbook"] = sorted(
                unique.values(),
                key=lambda item: (item["name"].casefold(), item["nist_id"]),
            )
        enriched.append(candidate)
    return enriched


def _fresh_cache_entry(entry, *, now):
    if entry is None:
        return False
    age = max(0.0, now - entry["fetched_at"])
    ttl = EMPTY_TTL_S if entry["empty"] else SUCCESS_TTL_S
    return age <= ttl


def _combined_status(statuses):
    states = set(statuses)
    if not states:
        return "unavailable"
    if "unavailable" in states:
        return "partial" if len(states) > 1 else "unavailable"
    if "live" in states and len(states) > 1:
        return "mixed"
    if "live" in states:
        return "live"
    if "stale" in states:
        return "stale"
    return "cache"
