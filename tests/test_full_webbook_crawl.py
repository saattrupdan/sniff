import datetime as dt
import email.utils
import gzip
import importlib.util
import json
import signal
import threading
import time
import urllib.error
from pathlib import Path

import pytest


def _crawler_module():
    path = Path(__file__).parents[1] / "scripts" / "crawl_webbook_species.py"
    spec = importlib.util.spec_from_file_location("crawl_webbook_species", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_sitemap_crawl_checkpoints_raw_pages_and_reparses_offline(
    tmp_path, monkeypatch
):
    crawler = _crawler_module()
    sitemap = "https://webbook.nist.gov/sitemap_1.xml.gz"
    index = f"""<?xml version='1.0'?>
    <sitemapindex xmlns='http://www.sitemaps.org/schemas/sitemap/0.9'>
      <sitemap><loc>{sitemap}</loc></sitemap>
    </sitemapindex>""".encode()
    inchi_url = "https://webbook.nist.gov/cgi/inchi/InChI%3D1S/C2H4O/c1-2-3/h2H%2C1H3"
    urls = gzip.compress(
        f"""<?xml version='1.0'?>
        <urlset xmlns='http://www.sitemaps.org/schemas/sitemap/0.9'>
          <url><loc>https://webbook.nist.gov/cgi/cbook.cgi?ID=C75070</loc></url>
          <url><loc>https://webbook.nist.gov/cgi/cbook.cgi?ID=U109047</loc></url>
          <url><loc>{inchi_url}</loc></url>
          <url><loc>https://example.com/not-a-species</loc></url>
        </urlset>""".encode()
    )
    detail = b"""<!doctype html><html><body><main>
      <h1 id='Top'>Acetaldehyde</h1><ul>
      <li><strong>Formula:</strong> C<sub>2</sub>H<sub>4</sub>O</li>
      <li><strong>Molecular weight:</strong> 44.05</li>
      <li><strong>IUPAC Standard InChI:</strong>
      InChI=1S/C2H4O/c1-2-3/h2H,1H3</li>
      <li><strong>CAS Registry Number:</strong> 75-07-0</li>
      </ul></main></body></html>"""
    responses = {
        crawler.ROBOTS_URL: b"User-agent: *\nDisallow: /cdn-cgi/\nCrawl-delay: 5\n",
        crawler.SITEMAP_INDEX: index,
        sitemap: urls,
        "https://webbook.nist.gov/cgi/cbook.cgi?ID=C75070&Units=SI": detail,
        "https://webbook.nist.gov/cgi/cbook.cgi?ID=U109047&Units=SI": detail,
        inchi_url: detail,
    }

    class FakeClient:
        def __init__(self, **_kwargs):
            self.crawl_delay = 5

        def fetch_page(self, url):
            return responses[url]

    monkeypatch.setattr(crawler.nist_webbook, "WebBookClient", FakeClient)
    state = tmp_path / "state.sqlite3"
    cache = tmp_path / "pages"

    manifest, count = crawler.discover(state, cache, minimum_species=1)
    assert count == 3
    assert crawler.status(state, emit=False)["remaining"] == 3
    assert crawler.crawl(state, cache) == 3

    summary = crawler.status(state, emit=False)
    assert summary["done"] == 3
    assert summary["remaining"] == 0
    with crawler._connect(state) as connection:
        row = connection.execute(
            "SELECT response_sha,body_path,detail_json FROM species_job"
        ).fetchone()
    assert Path(row["body_path"]).is_file()
    assert json.loads(row["detail_json"])["name"] == "Acetaldehyde"
    with crawler._connect(state) as connection:
        u_detail = json.loads(
            connection.execute(
                "SELECT detail_json FROM species_job WHERE nist_id='U109047'"
            ).fetchone()[0]
        )
    assert u_detail["webbook_id"] == "U109047"

    monkeypatch.setattr(
        crawler.nist_webbook,
        "WebBookClient",
        lambda **_kwargs: (_ for _ in ()).throw(AssertionError("network used")),
    )
    crawler.reparse(state)
    assert crawler.status(state, emit=False)["done"] == 3
    with crawler._connect(state) as connection:
        reparsed = json.loads(
            connection.execute(
                "SELECT detail_json FROM species_job WHERE nist_id='U109047'"
            ).fetchone()[0]
        )
    assert reparsed["webbook_id"] == "U109047"
    assert manifest


def test_cache_write_verifies_and_replaces_corrupt_content(tmp_path):
    crawler = _crawler_module()
    body = b"authoritative response"
    digest = crawler._sha(body)
    path = crawler._cache_path(tmp_path, digest)
    path.parent.mkdir(parents=True)
    with gzip.open(path, "wb") as stream:
        stream.write(b"corrupt response")

    written_digest, written_path = crawler._write_body(tmp_path, body)

    assert written_digest == digest
    assert written_path == path
    assert crawler._read_body(path, digest) == body


def test_cache_write_failure_is_not_a_network_retry(tmp_path):
    crawler = _crawler_module()
    blocked = tmp_path / "not-a-directory"
    blocked.write_text("x", encoding="utf-8")

    with pytest.raises(crawler.LocalPersistenceError):
        crawler._write_body(blocked / "pages", b"response")


def test_lost_crawler_lease_aborts_before_more_work(tmp_path):
    crawler = _crawler_module()
    connection = crawler._connect(tmp_path / "state.sqlite3")
    owner = crawler._acquire_lock(connection)
    connection.execute("DELETE FROM crawl_lock")
    connection.commit()

    with pytest.raises(RuntimeError, match="ownership was lost"):
        crawler._heartbeat(connection, owner)


def test_stale_owner_cannot_finish_a_reclaimed_species_job(tmp_path):
    crawler = _crawler_module()
    connection = crawler._connect(tmp_path / "state.sqlite3")
    connection.execute(
        "INSERT INTO species_job(nist_id,url,status,claim_owner) VALUES ('C1','https://webbook.nist.gov/cgi/cbook.cgi?ID=C1&Units=SI','fetching','new-owner')"
    )
    connection.commit()

    with pytest.raises(RuntimeError, match="ownership was lost"):
        crawler._update_claimed(connection, "old-owner", "C1", "status='done'", ())
    assert (
        connection.execute(
            "SELECT status FROM species_job WHERE nist_id='C1'"
        ).fetchone()[0]
        == "fetching"
    )


def test_setup_failure_releases_crawler_lock(tmp_path, monkeypatch):
    crawler = _crawler_module()
    state = tmp_path / "state.sqlite3"
    connection = crawler._connect(state)
    connection.execute("INSERT INTO generation VALUES ('manifest','index',0,0,1)")
    connection.execute("INSERT INTO crawl_meta VALUES ('active_manifest','manifest')")
    connection.commit()
    connection.close()
    monkeypatch.setattr(
        crawler.nist_webbook,
        "WebBookClient",
        lambda **_kwargs: (_ for _ in ()).throw(RuntimeError("setup failed")),
    )

    with pytest.raises(RuntimeError, match="setup failed"):
        crawler.crawl(state, tmp_path / "pages")

    with crawler._connect(state) as connection:
        assert connection.execute("SELECT COUNT(*) FROM crawl_lock").fetchone()[0] == 0

    monkeypatch.setattr(
        crawler.nist_webbook,
        "WebBookClient",
        lambda **_kwargs: (_ for _ in ()).throw(KeyboardInterrupt()),
    )
    with pytest.raises(KeyboardInterrupt):
        crawler.crawl(state, tmp_path / "pages")
    with crawler._connect(state) as connection:
        assert connection.execute("SELECT COUNT(*) FROM crawl_lock").fetchone()[0] == 0

    real_acquire = crawler._acquire_lock

    def interrupted_acquire(connection, owner):
        assert signal.getsignal(signal.SIGTERM) is crawler._stop
        real_acquire(connection, owner)
        raise KeyboardInterrupt()

    monkeypatch.setattr(crawler, "_acquire_lock", interrupted_acquire)
    with pytest.raises(KeyboardInterrupt):
        crawler.crawl(state, tmp_path / "pages")
    with crawler._connect(state) as connection:
        assert connection.execute("SELECT COUNT(*) FROM crawl_lock").fetchone()[0] == 0


def test_stop_request_interrupts_server_cooldown(tmp_path):
    crawler = _crawler_module()
    connection = crawler._connect(tmp_path / "state.sqlite3")
    owner = crawler._acquire_lock(connection)
    crawler.STOP_REQUESTED = False
    crawler.STOP_EVENT.clear()
    sleeper = threading.Thread(
        target=crawler._sleep_with_heartbeat,
        args=(connection, owner, 24 * 60 * 60),
        daemon=True,
    )
    sleeper.start()
    time.sleep(0.02)

    crawler._stop(None, None)
    sleeper.join(timeout=0.5)

    assert not sleeper.is_alive()
    with pytest.raises(crawler.CrawlStopped):
        crawler._scheduled_wait(connection, owner, 24 * 60 * 60)
    crawler._release_lock(connection, owner)
    crawler.STOP_REQUESTED = False
    crawler.STOP_EVENT.clear()


def test_eighth_network_failure_becomes_explicitly_exhausted(tmp_path):
    crawler = _crawler_module()
    connection = crawler._connect(tmp_path / "state.sqlite3")
    connection.execute(
        "INSERT INTO species_job(nist_id,url,attempts,status,claim_owner) VALUES ('C1','https://webbook.nist.gov/cgi/cbook.cgi?ID=C1&Units=SI',7,'fetching','test-owner')"
    )
    row = connection.execute("SELECT * FROM species_job WHERE nist_id='C1'").fetchone()

    crawler._record_retry(
        connection, row, OSError("offline"), time.time(), "test-owner"
    )
    connection.commit()

    status, attempts, owner = connection.execute(
        "SELECT status,attempts,claim_owner FROM species_job WHERE nist_id='C1'"
    ).fetchone()
    assert (status, attempts, owner) == ("exhausted", 8, None)


def test_retry_after_is_honoured_and_capped():
    crawler = _crawler_module()
    error = urllib.error.HTTPError(
        "https://webbook.nist.gov/", 429, "slow down", {"Retry-After": "7200"}, None
    )
    assert crawler._retry_after(error) == 7200
    zero = urllib.error.HTTPError(
        "https://webbook.nist.gov/", 429, "slow down", {"Retry-After": "0"}, None
    )
    assert crawler._retry_after(zero) == 0
    future = email.utils.format_datetime(
        dt.datetime.now(dt.timezone.utc) + dt.timedelta(seconds=120)
    )
    dated = urllib.error.HTTPError(
        "https://webbook.nist.gov/", 503, "wait", {"Retry-After": future}, None
    )
    assert 118 <= crawler._retry_after(dated) <= 120


def test_cached_robots_policy_keeps_its_original_refresh_deadline(tmp_path):
    crawler = _crawler_module()
    calls = []

    class FakeClient:
        crawl_delay = 5

        def fetch_page(self, url):
            calls.append(url)
            return b"User-agent: *\nAllow: /\nCrawl-delay: 5\n"

    connection = crawler._connect(tmp_path / "state.sqlite3")
    client = FakeClient()
    crawler._refresh_robots(connection, client, tmp_path / "pages")
    fetched_at = time.time() - 23 * 60 * 60
    connection.execute(
        "UPDATE document SET updated_at=? WHERE url=?",
        (fetched_at, crawler.ROBOTS_URL),
    )
    connection.commit()

    _policy, checked_at = crawler._refresh_robots(
        connection, client, tmp_path / "pages"
    )

    assert checked_at == fetched_at
    assert calls == [crawler.ROBOTS_URL]


def test_species_manifest_accepts_inchi_urls_and_only_known_sentinel():
    crawler = _crawler_module()
    inchi = "https://webbook.nist.gov/cgi/inchi/InChI%3D1S/H2O/h1H2"

    key, url = crawler._species_entry(inchi)

    assert key.startswith("U") and url == inchi
    assert crawler._species_entry("https://webbook.nist.gov/cgi/cbook.cgi?ID=x") is None
    with pytest.raises(RuntimeError, match="malformed species URL"):
        crawler._species_entry("https://webbook.nist.gov/cgi/cbook.cgi?ID=bad")
    for malformed in (
        "https://webbook.nist.gov/cgi/inchi/InChI=1S/H2O/h1H2",
        "https://webbook.nist.gov/cgi/inchi/InChI%3d1S/H2O/h1H2",
        "https://webbook.nist.gov/cgi/inchi/InChI%3D1S/",
        "https://webbook.nist.gov/cgi/inchi/InChI%3D1S/H2O/%ZZ",
        "https://webbook.nist.gov/cgi/inchi/InChI%3D1S/H2O/%FF",
        "https://webbook.nist.gov/cgi/inchi/InChI%3D1S%2FH2O/h1H2",
        "https://webbook.nist.gov/cgi/inchi/InChI%3D1S/H2O//h1H2",
        "https://webbook.nist.gov/cgi/inchi/InChI%3D1S/H2O/%28h1H2%29",
    ):
        with pytest.raises(RuntimeError, match="malformed InChI URL"):
            crawler._species_entry(malformed)
    assert (
        crawler._webbook_id(
            "U109047",
            "https://webbook.nist.gov/cgi/cbook.cgi?ID=U109047&Units=SI",
        )
        == "U109047"
    )


def test_specific_robots_group_controls_species_requests():
    crawler = _crawler_module()
    policy = crawler._parse_robots(
        b"""User-agent: *
        Allow: /
        Crawl-delay: 2
        User-agent: Sniff
        Allow: /sitemap_
        Disallow: /cgi/
        Crawl-delay: 7
        """
    )

    assert policy["delay"] == 7
    assert crawler._robots_allows(policy, crawler.SITEMAP_INDEX)
    assert not crawler._robots_allows(
        policy, "https://webbook.nist.gov/cgi/cbook.cgi?ID=C75070&Units=SI"
    )
