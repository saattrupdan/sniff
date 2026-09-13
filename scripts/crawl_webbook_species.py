#!/usr/bin/env python3
"""Cache every canonical NIST WebBook species page for offline catalogue builds.

This is a maintainer operation, never a Sniff runtime dependency. Discovery freezes one
WebBook sitemap generation, and the crawl checkpoints each species page in SQLite plus a
content-addressed gzip cache. Interruptions are safe: rerun the same command to resume.
"""

from __future__ import annotations

import argparse
import datetime as dt
import email.utils
import gzip
import hashlib
import json
import os
import re
import signal
import sqlite3
import sys
import threading
import time
import urllib.error
import urllib.parse
import xml.etree.ElementTree as ET
from pathlib import Path

from sniff import nist_webbook

STATE_PATH = Path.home() / ".sniff" / "nist-webbook-full.sqlite3"
CACHE_DIR = Path.home() / ".sniff" / "nist-webbook-pages"
SITEMAP_INDEX = "https://webbook.nist.gov/sitemap_index.xml"
ROBOTS_URL = "https://webbook.nist.gov/robots.txt"
SITEMAP_NS = {"s": "http://www.sitemaps.org/schemas/sitemap/0.9"}
SPECIES_URL = re.compile(
    r"^https://webbook\.nist\.gov/cgi/cbook\.cgi\?ID=([A-Z][A-Z0-9]*)$"
)
SITEMAP_SENTINELS = {"https://webbook.nist.gov/cgi/cbook.cgi?ID=x"}
SITEMAP_URL = re.compile(r"^https://webbook\.nist\.gov/sitemap_[1-9][0-9]*\.xml\.gz$")
STOP_REQUESTED = False
STOP_EVENT = threading.Event()
LOCK_STALE_S = 60 * 60.0
ROBOTS_REFRESH_S = 24 * 60 * 60
MAX_ATTEMPTS = 8


class LocalPersistenceError(RuntimeError):
    """A local cache/state failure for which another HTTP request would be unsafe."""


class CrawlStopped(RuntimeError):
    """The maintainer requested a graceful stop before another HTTP request."""


def _connect(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(str(path), timeout=30.0)
    connection.row_factory = sqlite3.Row
    connection.executescript(
        """
        PRAGMA journal_mode = WAL;
        PRAGMA synchronous = FULL;
        CREATE TABLE IF NOT EXISTS document (
            url TEXT PRIMARY KEY,
            kind TEXT NOT NULL,
            status TEXT NOT NULL,
            response_sha TEXT,
            body_path TEXT,
            last_modified TEXT,
            error TEXT,
            updated_at REAL NOT NULL
        );
        CREATE TABLE IF NOT EXISTS generation (
            manifest_sha TEXT PRIMARY KEY,
            index_sha TEXT NOT NULL,
            discovered_at REAL NOT NULL,
            species_count INTEGER NOT NULL,
            complete INTEGER NOT NULL DEFAULT 0
        );
        CREATE TABLE IF NOT EXISTS generation_species (
            manifest_sha TEXT NOT NULL REFERENCES generation(manifest_sha),
            nist_id TEXT NOT NULL,
            PRIMARY KEY (manifest_sha, nist_id)
        );
        CREATE TABLE IF NOT EXISTS species_job (
            nist_id TEXT PRIMARY KEY,
            url TEXT NOT NULL,
            status TEXT NOT NULL DEFAULT 'pending',
            attempts INTEGER NOT NULL DEFAULT 0,
            retry_at REAL,
            error TEXT,
            response_sha TEXT,
            body_path TEXT,
            detail_json TEXT,
            claim_owner TEXT,
            updated_at REAL
        );
        CREATE INDEX IF NOT EXISTS species_job_status
            ON species_job(status, retry_at, nist_id);
        CREATE TABLE IF NOT EXISTS crawl_meta (
            key TEXT PRIMARY KEY,
            value TEXT NOT NULL
        );
        CREATE TABLE IF NOT EXISTS crawl_lock (
            singleton INTEGER PRIMARY KEY CHECK (singleton = 1),
            owner TEXT NOT NULL,
            heartbeat REAL NOT NULL
        );
        """
    )
    columns = {row[1] for row in connection.execute("PRAGMA table_info(species_job)")}
    if "claim_owner" not in columns:
        connection.execute("ALTER TABLE species_job ADD COLUMN claim_owner TEXT")
    connection.commit()
    return connection


def _sha(body):
    return hashlib.sha256(body).hexdigest()


def _cache_path(cache_dir, digest):
    return cache_dir / digest[:2] / f"{digest}.gz"


def _write_body(cache_dir, body):
    digest = _sha(body)
    path = _cache_path(cache_dir, digest)
    if path.exists():
        try:
            _read_body(path, digest)
            return digest, path
        except (OSError, RuntimeError):
            pass
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_name(path.name + f".{os.getpid()}.tmp")
        with gzip.open(temporary, "wb", compresslevel=6) as stream:
            stream.write(body)
        temporary.replace(path)
        _read_body(path, digest)
    except (OSError, RuntimeError) as exc:
        raise LocalPersistenceError(
            f"could not persist the response cache: {exc}"
        ) from exc
    return digest, path


def _stored_body_path(cache_dir, stored_path, expected_sha):
    path = Path(stored_path)
    if not path.is_absolute():
        return Path(cache_dir) / path
    if path.exists():
        return path
    return _cache_path(Path(cache_dir), expected_sha)


def _relative_body_path(cache_dir, path):
    return str(Path(path).relative_to(Path(cache_dir)))


def _read_body(path, expected_sha):
    with gzip.open(path, "rb") as stream:
        body = stream.read()
    if _sha(body) != expected_sha:
        raise RuntimeError(f"cached response failed SHA-256 verification: {path}")
    return body


def _fetch_document(connection, client, cache_dir, url, kind, *, force=False):
    row = connection.execute(
        "SELECT response_sha, body_path, status FROM document WHERE url=?", (url,)
    ).fetchone()
    if (
        not force
        and row
        and row["status"] == "done"
        and row["response_sha"]
        and row["body_path"]
    ):
        try:
            return _read_body(
                _stored_body_path(cache_dir, row["body_path"], row["response_sha"]),
                row["response_sha"],
            )
        except (OSError, RuntimeError):
            pass
    now = time.time()
    try:
        body = client.fetch_page(url)
        digest, path = _write_body(cache_dir, body)
        connection.execute(
            """
            INSERT INTO document(url,kind,status,response_sha,body_path,error,updated_at)
            VALUES (?,?,'done',?,?,NULL,?)
            ON CONFLICT(url) DO UPDATE SET kind=excluded.kind,status='done',
                response_sha=excluded.response_sha,body_path=excluded.body_path,
                error=NULL,updated_at=excluded.updated_at
            """,
            (url, kind, digest, _relative_body_path(cache_dir, path), now),
        )
        connection.commit()
        return body
    except LocalPersistenceError:
        raise
    except Exception as exc:
        connection.execute(
            """
            INSERT INTO document(url,kind,status,error,updated_at)
            VALUES (?,?,'error',?,?)
            ON CONFLICT(url) DO UPDATE SET status='error',error=excluded.error,
                updated_at=excluded.updated_at
            """,
            (url, kind, str(exc), now),
        )
        connection.commit()
        raise


def _xml_body(body, url):
    if url.endswith(".gz"):
        try:
            return gzip.decompress(body)
        except OSError as exc:
            raise RuntimeError(f"invalid gzip sitemap: {url}") from exc
    return body


def _locations(body, url, root_tag):
    try:
        root = ET.fromstring(_xml_body(body, url))
    except ET.ParseError as exc:
        raise RuntimeError(f"invalid sitemap XML: {url}") from exc
    expected = "{" + SITEMAP_NS["s"] + "}" + root_tag
    if root.tag != expected:
        raise RuntimeError(f"unexpected sitemap root in {url}: {root.tag}")
    return [
        node.text.strip()
        for node in root.findall(
            "s:" + ("sitemap" if root_tag == "sitemapindex" else "url") + "/s:loc",
            SITEMAP_NS,
        )
        if node.text and node.text.strip()
    ]


def _parse_robots(body, agent="Sniff"):
    groups = []
    agents = []
    rules = []
    delay = None
    have_directives = False
    for raw_line in body.decode("utf-8", "replace").splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line or ":" not in line:
            continue
        field, value = (part.strip() for part in line.split(":", 1))
        field = field.casefold()
        if field == "user-agent":
            if have_directives and agents:
                groups.append((agents, rules, delay))
                agents, rules, delay, have_directives = [], [], None, False
            agents.append(value.casefold())
        elif agents and field in ("allow", "disallow"):
            rules.append((field == "allow", value))
            have_directives = True
        elif agents and field == "crawl-delay":
            delay = float(value)
            have_directives = True
    if agents:
        groups.append((agents, rules, delay))
    needle = agent.casefold()
    ranked = []
    for group_agents, group_rules, group_delay in groups:
        matches = [
            0 if candidate == "*" else len(candidate)
            for candidate in group_agents
            if candidate == "*" or candidate in needle
        ]
        if matches:
            ranked.append((max(matches), group_rules, group_delay))
    if not ranked:
        raise RuntimeError(
            "the WebBook robots policy has no applicable user-agent group"
        )
    specificity = max(row[0] for row in ranked)
    selected = [row for row in ranked if row[0] == specificity]
    selected_rules = [rule for _, group_rules, _ in selected for rule in group_rules]
    delays = [value for _, _, value in selected if value is not None]
    return {"delay": max([5.0, *delays]), "rules": selected_rules}


def _robots_allows(policy, url):
    parsed = urllib.parse.urlparse(url)
    target = parsed.path + (("?" + parsed.query) if parsed.query else "")
    matches = []
    for allow, pattern in policy["rules"]:
        if not pattern:
            continue
        expression = re.escape(pattern).replace(r"\*", ".*")
        if expression.endswith(r"\$"):
            expression = expression[:-2] + "$"
        if re.match(expression, target):
            matches.append((len(pattern.replace("*", "")), bool(allow)))
    if not matches:
        return True
    longest = max(length for length, _allow in matches)
    return any(allow for length, allow in matches if length == longest)


def _refresh_robots(connection, client, cache_dir, *, force=False):
    row = connection.execute(
        "SELECT updated_at FROM document WHERE url=? AND status='done'", (ROBOTS_URL,)
    ).fetchone()
    expired = row is None or time.time() - float(row[0]) >= ROBOTS_REFRESH_S
    body = _fetch_document(
        connection, client, cache_dir, ROBOTS_URL, "robots", force=force or expired
    )
    policy = _parse_robots(body)
    if not _robots_allows(policy, ROBOTS_URL):
        raise RuntimeError("the current WebBook robots policy disallows policy refresh")
    client.crawl_delay = policy["delay"]
    connection.execute(
        "INSERT INTO crawl_meta(key,value) VALUES ('robots_policy',?) "
        "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
        (json.dumps(policy, sort_keys=True),),
    )
    connection.commit()
    updated = connection.execute(
        "SELECT updated_at FROM document WHERE url=?", (ROBOTS_URL,)
    ).fetchone()
    return policy, float(updated[0])


def _species_entry(location):
    if location in SITEMAP_SENTINELS:
        return None
    match = SPECIES_URL.fullmatch(location)
    if match:
        return match.group(1), location
    parsed = urllib.parse.urlparse(location)
    if parsed.scheme == "https" and parsed.netloc == "webbook.nist.gov":
        if parsed.path == "/cgi/cbook.cgi":
            raise RuntimeError(f"malformed species URL in sitemap: {location}")
        if parsed.path.startswith("/cgi/inchi/"):
            raw_value = parsed.path.removeprefix("/cgi/inchi/")
            try:
                value = urllib.parse.unquote_to_bytes(raw_value).decode("utf-8")
            except UnicodeDecodeError:
                value = ""
            canonical = urllib.parse.quote(value, safe="/().-*")
            layers = value.split("/")
            if (
                parsed.query
                or parsed.fragment
                or not raw_value.startswith("InChI%3D1S/")
                or canonical != raw_value
                or len(layers) < 2
                or any(not layer for layer in layers)
            ):
                raise RuntimeError(f"malformed InChI URL in sitemap: {location}")
            return "U" + _sha(location.encode("utf-8")), location
    return None


def _webbook_id(key, url):
    parsed = urllib.parse.urlparse(url)
    identifiers = urllib.parse.parse_qs(parsed.query).get("ID", [])
    if parsed.path == "/cgi/cbook.cgi" and identifiers == [key]:
        return key
    return None


def discover(state_path=STATE_PATH, cache_dir=CACHE_DIR, *, minimum_species=100_000):
    connection = _connect(state_path)
    client = nist_webbook.WebBookClient(timeout=30.0)
    policy, _policy_checked_at = _refresh_robots(connection, client, cache_dir)
    if not _robots_allows(policy, SITEMAP_INDEX):
        raise RuntimeError(
            "the current WebBook robots policy disallows sitemap discovery"
        )
    index = _fetch_document(
        connection, client, cache_dir, SITEMAP_INDEX, "sitemap-index"
    )
    sitemap_urls = _locations(index, SITEMAP_INDEX, "sitemapindex")
    if not sitemap_urls:
        raise RuntimeError("the WebBook sitemap index contained no child sitemaps")
    species = {}
    for number, url in enumerate(sitemap_urls, 1):
        nist_webbook._validate_webbook_url(url)
        if not SITEMAP_URL.fullmatch(url):
            raise RuntimeError(f"unexpected child sitemap URL: {url}")
        if not _robots_allows(policy, url):
            raise RuntimeError(f"the current robots policy disallows sitemap: {url}")
        body = _fetch_document(connection, client, cache_dir, url, "sitemap")
        for location in _locations(body, url, "urlset"):
            entry = _species_entry(location)
            if entry:
                key, species_url = entry
                previous = species.setdefault(key, species_url)
                if previous != species_url:
                    raise RuntimeError(f"species manifest key collision: {key}")
        print(
            f"sitemap {number}/{len(sitemap_urls)}: {len(species)} species", flush=True
        )
    if len(species) < minimum_species:
        raise RuntimeError(
            f"refusing an implausibly small WebBook manifest ({len(species)} species)"
        )
    ordered = sorted(species.items())
    manifest = "".join(f"{key}\t{url}\n" for key, url in ordered)
    manifest_sha = _sha(manifest.encode("utf-8"))
    now = time.time()
    connection.execute(
        """
        INSERT OR IGNORE INTO generation(
            manifest_sha,index_sha,discovered_at,species_count,complete
        ) VALUES (?,?,?,?,0)
        """,
        (manifest_sha, _sha(index), now, len(ordered)),
    )
    connection.executemany(
        "INSERT OR IGNORE INTO generation_species(manifest_sha,nist_id) VALUES (?,?)",
        ((manifest_sha, key) for key, _url in ordered),
    )
    connection.executemany(
        """
        INSERT OR IGNORE INTO species_job(nist_id,url)
        VALUES (?,?)
        """,
        (
            (
                key,
                (
                    f"https://webbook.nist.gov/cgi/cbook.cgi?ID={key}&Units=SI"
                    if SPECIES_URL.fullmatch(url)
                    else url
                ),
            )
            for key, url in ordered
        ),
    )
    connection.execute(
        "UPDATE generation SET complete=1 WHERE manifest_sha=?", (manifest_sha,)
    )
    connection.execute(
        "INSERT INTO crawl_meta(key,value) VALUES ('active_manifest',?) "
        "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
        (manifest_sha,),
    )
    connection.commit()
    connection.close()
    print(f"active manifest {manifest_sha}: {len(ordered)} species", flush=True)
    return manifest_sha, len(ordered)


def _active_manifest(connection):
    row = connection.execute(
        "SELECT value FROM crawl_meta WHERE key='active_manifest'"
    ).fetchone()
    if row is None:
        raise RuntimeError("run discover before crawling")
    generation = connection.execute(
        "SELECT * FROM generation WHERE manifest_sha=? AND complete=1", (row[0],)
    ).fetchone()
    if generation is None:
        raise RuntimeError("the active sitemap generation is incomplete")
    return row[0]


def _lock_owner():
    return (
        f"{os.uname().nodename if hasattr(os, 'uname') else sys.platform}:{os.getpid()}"
    )


def _acquire_lock(connection, owner=None):
    owner = owner or _lock_owner()
    now = time.time()
    connection.execute("BEGIN IMMEDIATE")
    row = connection.execute(
        "SELECT owner,heartbeat FROM crawl_lock WHERE singleton=1"
    ).fetchone()
    if row and now - float(row["heartbeat"]) < LOCK_STALE_S:
        connection.rollback()
        raise RuntimeError(f"another crawler owns this state: {row['owner']}")
    connection.execute(
        "INSERT INTO crawl_lock(singleton,owner,heartbeat) VALUES (1,?,?) "
        "ON CONFLICT(singleton) DO UPDATE SET owner=excluded.owner,heartbeat=excluded.heartbeat",
        (owner, now),
    )
    connection.execute(
        "UPDATE species_job SET status='pending',claim_owner=NULL "
        "WHERE status='fetching' AND (claim_owner IS NULL OR claim_owner != ?)",
        (owner,),
    )
    connection.commit()
    return owner


def _heartbeat(connection, owner):
    cursor = connection.execute(
        "UPDATE crawl_lock SET heartbeat=? WHERE singleton=1 AND owner=?",
        (time.time(), owner),
    )
    connection.commit()
    if cursor.rowcount != 1:
        raise RuntimeError("crawler lease ownership was lost; aborting")


def _release_lock(connection, owner):
    connection.execute("DELETE FROM crawl_lock WHERE singleton=1 AND owner=?", (owner,))
    connection.commit()


def _stop(_signum, _frame):
    global STOP_REQUESTED
    STOP_REQUESTED = True
    STOP_EVENT.set()


def _update_claimed(connection, owner, nist_id, assignments, values):
    cursor = connection.execute(
        f"UPDATE species_job SET {assignments} "
        "WHERE nist_id=? AND status='fetching' AND claim_owner=?",
        (*values, nist_id, owner),
    )
    if cursor.rowcount != 1:
        connection.rollback()
        raise RuntimeError("species job ownership was lost; aborting")


def crawl(state_path=STATE_PATH, cache_dir=CACHE_DIR, max_records=None):
    global STOP_REQUESTED
    STOP_REQUESTED = False
    STOP_EVENT.clear()
    connection = _connect(state_path)
    manifest = _active_manifest(connection)
    owner = _lock_owner()
    previous_handlers = {
        signal.SIGINT: signal.signal(signal.SIGINT, _stop),
        signal.SIGTERM: signal.signal(signal.SIGTERM, _stop),
    }
    processed = 0
    try:
        if STOP_EVENT.is_set():
            raise CrawlStopped("crawl stopped before lock acquisition")
        _acquire_lock(connection, owner)
        if STOP_EVENT.is_set():
            raise CrawlStopped("crawl stopped after lock acquisition")
        client = nist_webbook.WebBookClient(
            timeout=30.0,
            waiter=lambda delay: _scheduled_wait(connection, owner, delay),
        )
        policy, policy_checked_at = _refresh_robots(connection, client, cache_dir)
        consecutive_server_failures = 0
        while not STOP_REQUESTED and (max_records is None or processed < max_records):
            row = connection.execute(
                """
                SELECT j.* FROM species_job j
                JOIN generation_species g ON g.nist_id=j.nist_id
                WHERE g.manifest_sha=?
                  AND j.status IN ('pending','error')
                  AND j.attempts < ?
                  AND (j.retry_at IS NULL OR j.retry_at <= ?)
                ORDER BY CASE j.status WHEN 'pending' THEN 0 ELSE 1 END, j.nist_id
                LIMIT 1
                """,
                (manifest, MAX_ATTEMPTS, time.time()),
            ).fetchone()
            if row is None:
                waiting = connection.execute(
                    """
                    SELECT MIN(j.retry_at) FROM species_job j
                    JOIN generation_species g ON g.nist_id=j.nist_id
                    WHERE g.manifest_sha=? AND j.status='error' AND j.attempts < ?
                    """,
                    (manifest, MAX_ATTEMPTS),
                ).fetchone()[0]
                if waiting is None:
                    break
                _sleep_with_heartbeat(
                    connection,
                    owner,
                    max(1.0, min(60.0, float(waiting) - time.time())),
                )
                continue
            if time.time() - policy_checked_at >= ROBOTS_REFRESH_S:
                policy, policy_checked_at = _refresh_robots(
                    connection, client, cache_dir, force=True
                )
            if not _robots_allows(policy, row["url"]):
                raise RuntimeError(
                    f"the current WebBook robots policy disallows {row['url']}"
                )
            connection.execute("BEGIN IMMEDIATE")
            claimed = connection.execute(
                """
                UPDATE species_job SET status='fetching',claim_owner=?,updated_at=?
                WHERE nist_id=? AND status IN ('pending','error')
                """,
                (owner, time.time(), row["nist_id"]),
            ).rowcount
            connection.commit()
            if claimed != 1:
                continue
            now = time.time()
            digest = path = None
            try:
                body = client.fetch_page(row["url"])
                digest, path = _write_body(cache_dir, body)
                detail = nist_webbook._parse_detail(body, nist_id=row["nist_id"])
                detail["url"] = row["url"]
                detail["webbook_id"] = _webbook_id(row["nist_id"], row["url"])
                detail["manifest_key"] = row["nist_id"]
                _update_claimed(
                    connection,
                    owner,
                    row["nist_id"],
                    "status='done',attempts=attempts+1,retry_at=NULL,error=NULL,"
                    "response_sha=?,body_path=?,detail_json=?,claim_owner=NULL,updated_at=?",
                    (
                        digest,
                        _relative_body_path(cache_dir, path),
                        json.dumps(detail, ensure_ascii=True, sort_keys=True),
                        now,
                    ),
                )
                consecutive_server_failures = 0
            except LocalPersistenceError as exc:
                _update_claimed(
                    connection,
                    owner,
                    row["nist_id"],
                    "status='persistence_error',error=?,claim_owner=NULL,updated_at=?",
                    (str(exc), now),
                )
                connection.commit()
                raise
            except urllib.error.HTTPError as exc:
                if exc.code in (401, 403):
                    _update_claimed(
                        connection,
                        owner,
                        row["nist_id"],
                        "status='policy_error',error=?,claim_owner=NULL,updated_at=?",
                        (f"HTTP {exc.code}", now),
                    )
                    connection.commit()
                    raise RuntimeError(
                        f"WebBook refused crawler access with HTTP {exc.code}"
                    )
                if exc.code in (404, 410):
                    _update_claimed(
                        connection,
                        owner,
                        row["nist_id"],
                        "status='missing',attempts=attempts+1,retry_at=NULL,error=?,"
                        "claim_owner=NULL,updated_at=?",
                        (f"HTTP {exc.code}", now),
                    )
                    consecutive_server_failures = 0
                elif exc.code in (408, 425, 429) or 500 <= exc.code < 600:
                    _record_retry(connection, row, exc, now, owner)
                    connection.commit()
                    consecutive_server_failures += 1
                    retry_after = _retry_after(exc)
                    delay = (
                        retry_after
                        if retry_after is not None
                        else min(
                            3600.0,
                            30.0 * (2 ** min(consecutive_server_failures - 1, 7)),
                        )
                    )
                    client.defer_requests(delay)
                    _sleep_with_heartbeat(connection, owner, delay)
                    if consecutive_server_failures >= 8:
                        raise RuntimeError(
                            "repeated WebBook throttling/server failures"
                        )
                else:
                    _update_claimed(
                        connection,
                        owner,
                        row["nist_id"],
                        "status='http_error',attempts=attempts+1,error=?,"
                        "claim_owner=NULL,updated_at=?",
                        (f"HTTP {exc.code}", now),
                    )
                    consecutive_server_failures = 0
            except Exception as exc:
                # A successfully cached but currently unparsable page is terminal for
                # network crawling and can be retried offline with `reparse`.
                if digest is not None and path is not None:
                    _update_claimed(
                        connection,
                        owner,
                        row["nist_id"],
                        "status='parse_error',attempts=attempts+1,retry_at=NULL,error=?,"
                        "response_sha=?,body_path=?,claim_owner=NULL,updated_at=?",
                        (
                            str(exc),
                            digest,
                            _relative_body_path(cache_dir, path),
                            now,
                        ),
                    )
                else:
                    _record_retry(connection, row, exc, now, owner)
            connection.commit()
            processed += 1
            _heartbeat(connection, owner)
            if processed % 100 == 0:
                summary = status(state_path, emit=False)
                print(
                    f"processed {processed}; done={summary['done']} "
                    f"remaining={summary['remaining']} eta={summary['minimum_eta']}",
                    flush=True,
                )
    finally:
        _release_lock(connection, owner)
        connection.close()
        for signum, handler in previous_handlers.items():
            signal.signal(signum, handler)
    summary = status(state_path, emit=False)
    blocking = sum(
        summary.get(key, 0)
        for key in ("exhausted", "persistence_error", "policy_error")
    )
    if blocking:
        raise RuntimeError(f"the crawl has {blocking} blocking failed jobs")
    return processed


def _retry_after(error):
    headers = getattr(error, "headers", None) or {}
    value = headers.get("Retry-After")
    if not value:
        return None
    try:
        delay = float(value)
    except (TypeError, ValueError):
        try:
            parsed = email.utils.parsedate_to_datetime(value)
            if parsed.tzinfo is None:
                parsed = parsed.replace(tzinfo=dt.timezone.utc)
            delay = parsed.timestamp() - time.time()
        except (TypeError, ValueError, OverflowError):
            return None
    return max(0.0, min(24 * 60 * 60.0, delay))


def _sleep_with_heartbeat(connection, owner, delay):
    deadline = time.time() + max(0.0, delay)
    while not STOP_EVENT.is_set() and time.time() < deadline:
        stopped = STOP_EVENT.wait(min(60.0, deadline - time.time()))
        if not stopped:
            _heartbeat(connection, owner)


def _scheduled_wait(connection, owner, delay):
    _sleep_with_heartbeat(connection, owner, delay)
    if STOP_EVENT.is_set():
        raise CrawlStopped("crawl stopped during the shared request wait")


def _record_retry(connection, row, exc, now, owner):
    attempts = int(row["attempts"]) + 1
    exhausted = attempts >= MAX_ATTEMPTS
    retry_at = (
        None if exhausted else now + min(3600.0, 30.0 * (2 ** min(attempts - 1, 7)))
    )
    _update_claimed(
        connection,
        owner,
        row["nist_id"],
        "status=?,attempts=?,retry_at=?,error=?,claim_owner=NULL,updated_at=?",
        (
            "exhausted" if exhausted else "error",
            attempts,
            retry_at,
            str(exc),
            now,
        ),
    )


def reparse(state_path=STATE_PATH, cache_dir=CACHE_DIR):
    connection = _connect(state_path)
    rows = connection.execute(
        """
        SELECT nist_id,url,response_sha,body_path FROM species_job
        WHERE body_path IS NOT NULL AND response_sha IS NOT NULL
        ORDER BY nist_id
        """
    ).fetchall()
    failed = 0
    for index, row in enumerate(rows, 1):
        try:
            body = _read_body(
                _stored_body_path(cache_dir, row["body_path"], row["response_sha"]),
                row["response_sha"],
            )
            detail = nist_webbook._parse_detail(body, nist_id=row["nist_id"])
            detail["url"] = row["url"]
            detail["webbook_id"] = _webbook_id(row["nist_id"], row["url"])
            detail["manifest_key"] = row["nist_id"]
            connection.execute(
                "UPDATE species_job SET status='done',error=NULL,detail_json=?,updated_at=? "
                "WHERE nist_id=?",
                (
                    json.dumps(detail, ensure_ascii=True, sort_keys=True),
                    time.time(),
                    row["nist_id"],
                ),
            )
        except Exception as exc:
            failed += 1
            connection.execute(
                "UPDATE species_job SET status='parse_error',error=?,updated_at=? WHERE nist_id=?",
                (str(exc), time.time(), row["nist_id"]),
            )
        if index % 1000 == 0:
            connection.commit()
    connection.commit()
    connection.close()
    print(f"reparsed {len(rows)} cached species; {failed} failures", flush=True)


def status(state_path=STATE_PATH, *, emit=True):
    connection = _connect(state_path)
    try:
        manifest = _active_manifest(connection)
    except RuntimeError:
        result = {"discovered": False, "state": str(state_path)}
    else:
        counts = dict(
            connection.execute(
                """
                SELECT j.status,COUNT(*) FROM species_job j
                JOIN generation_species g ON g.nist_id=j.nist_id
                WHERE g.manifest_sha=? GROUP BY j.status
                """,
                (manifest,),
            ).fetchall()
        )
        total = sum(counts.values())
        remaining = sum(counts.get(key, 0) for key in ("pending", "fetching", "error"))
        unresolved = sum(
            counts.get(key, 0)
            for key in (
                "parse_error",
                "http_error",
                "exhausted",
                "persistence_error",
                "policy_error",
            )
        )
        seconds = remaining * 5
        result = {
            "discovered": True,
            "state": str(state_path),
            "manifest_sha": manifest,
            "total": total,
            **counts,
            "remaining": remaining,
            "unresolved": unresolved,
            "network_complete": remaining == 0,
            "ready_to_classify": remaining == 0 and unresolved == 0,
            "minimum_eta_seconds": seconds,
            "minimum_eta": f"{seconds / 86400:.2f} days",
        }
    connection.close()
    if emit:
        print(json.dumps(result, indent=2, sort_keys=True))
    return result


def _parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("discover", "crawl", "status", "reparse"))
    parser.add_argument("--state", type=Path, default=STATE_PATH)
    parser.add_argument("--cache-dir", type=Path, default=CACHE_DIR)
    parser.add_argument("--max-records", type=int)
    return parser


def main(argv=None):
    args = _parser().parse_args(argv)
    if args.command == "discover":
        discover(args.state, args.cache_dir)
    elif args.command == "crawl":
        crawl(args.state, args.cache_dir, max_records=args.max_records)
    elif args.command == "reparse":
        reparse(args.state, args.cache_dir)
    else:
        status(args.state)


if __name__ == "__main__":
    main()
