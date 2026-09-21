import urllib.error
from collections import deque

import pytest

from sniff import formula_id, nist_webbook

FORMULA_RESULTS = b"""<!doctype html><html><body><main><ol>
<li><a href="/cgi/cbook.cgi?ID=C151564&amp;Units=SI">Ethylenimine</a>
(C<sub>2</sub>H<sub>5</sub>N)</li>
<li><a href="/cgi/cbook.cgi?ID=C27770441&amp;Units=SI">EtN radical</a>
(C<sub>2</sub>H<sub>5</sub>N)</li>
<li><a href="/cgi/cbook.cgi?ID=B1002815&amp;Units=SI">t-CD3CD=ND</a>
(C<sub>2</sub>D<sub>5</sub>N)</li>
</ol></main></body></html>"""

MASS_RESULTS = b"""<!doctype html><html><body><main><ol>
<li><strong>&nbsp; 576.66 </strong>
<a href="/cgi/cbook.cgi?ID=C12345&amp;Units=SI">Hentetracontane</a></li>
<li><strong>&nbsp; 576.66 </strong>
<a href="/cgi/cbook.cgi?ID=B999&amp;Units=SI">exotic radical</a></li>
</ol></main></body></html>"""

NON_COMPOUND_RESULTS = b"""<!doctype html><html><body><main><ol>
<li><a href="/cgi/cbook.cgi?ID=C25013862&amp;Units=SI">Ethene, homopolymer</a>
(CH<sub>2</sub>)</li>
<li><a href="/cgi/cbook.cgi?ID=C9002884&amp;Units=SI">Poly(methylene)</a>
(CH<sub>2</sub>)</li>
<li><a href="/cgi/cbook.cgi?ID=C2465567&amp;Units=SI">Methylene</a>
(CH<sub>2</sub>)</li>
<li><a href="/cgi/cbook.cgi?ID=B1000&amp;Units=SI">CH2</a>
(CH<sub>2</sub>)</li>
</ol></main></body></html>"""

DETAIL = b"""<!doctype html><html><body><main>
<h1 id="Top">Hentetracontane</h1><ul>
<li><strong><a>Formula</a>:</strong> C<sub>41</sub>H<sub>84</sub></li>
<li><strong>Molecular weight:</strong> 577.11</li>
<li><strong>IUPAC Standard InChI:</strong>
<span class="inchi-text">InChI=1S/C41H84/c1-3-5-7-9-11-13-15-17-19-21-23-25-27-29-31-33-35-37-39-41-40-38-36-34-32-30-28-26-24-22-20-18-16-14-12-10-8-6-4-2/h3-41H2,1-2H3</span></li>
<li><strong>IUPAC Standard InChIKey:</strong>
<span class="inchi-text">TESTKEY-UHFFFAOYSA-N</span></li>
<li><strong>CAS Registry Number:</strong> 7098-20-6</li>
</ul></main></body></html>"""


class Response:
    def __init__(self, body):
        self.body = body

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        return False

    def read(self, _size=-1):
        return self.body


class Opener:
    def __init__(self, *bodies):
        self.bodies = deque(bodies)
        self.requests = []

    def __call__(self, request, timeout):
        self.requests.append((request, timeout))
        if not self.bodies:
            raise AssertionError("unexpected WebBook request")
        body = self.bodies.popleft()
        if isinstance(body, Exception):
            raise body
        return Response(body)


def client(tmp_path, opener, **kwargs):
    return nist_webbook.WebBookClient(
        cache_path=tmp_path / "webbook.sqlite3",
        schedule_path=tmp_path / "schedule.sqlite3",
        opener=opener,
        crawl_delay=0,
        **kwargs,
    )


def test_detail_parser_distinguishes_unusable_records_from_parser_failures():
    registry_missing = b"<title>Registry Number Not Found</title>"
    search_results = b"<title>Search Results</title>"
    no_formula = b'<h1 id="Top">Coffee ground</h1><ul></ul>'

    for page in (registry_missing, search_results):
        with pytest.raises(nist_webbook.UnusableSpeciesError):
            nist_webbook._parse_detail(page, nist_id="B3000001")
    with pytest.raises(nist_webbook.UnusableSpeciesError):
        nist_webbook._parse_detail(no_formula, nist_id="C123")
    with pytest.raises(nist_webbook.WebBookError) as error:
        nist_webbook._parse_detail(b"<html>unexpected</html>", nist_id="C123")
    assert type(error.value) is nist_webbook.WebBookError


def test_formula_lookup_filters_radicals_and_isotopologues(tmp_path):
    opener = Opener(FORMULA_RESULTS)
    result = client(tmp_path, opener).lookup_formula("C2H5N")

    assert result["status"] == "live"
    assert [species["name"] for species in result["species"]] == ["Ethylenimine"]
    assert result["species"][0]["formula"] == "C2H5N"
    assert result["excluded"] == 2
    assert "Sniff/" in opener.requests[0][0].get_header("User-agent")


def test_formula_lookup_filters_polymers_transients_and_formula_labels(tmp_path):
    result = client(tmp_path, Opener(NON_COMPOUND_RESULTS)).lookup_formula("CH2")

    assert result["species"] == []
    assert result["excluded"] == 4


def test_fetch_policy_refuses_external_and_robots_disallowed_urls(tmp_path):
    opener = Opener(FORMULA_RESULTS)
    webbook = client(tmp_path, opener)

    for url in (
        "https://example.com/cgi/cbook.cgi?ID=C75070",
        "https://webbook.nist.gov/cdn-cgi/example",
    ):
        try:
            webbook.fetch_page(url)
        except nist_webbook.WebBookError:
            pass
        else:
            raise AssertionError(f"unsafe URL was accepted: {url}")
    assert opener.requests == []
    assert (
        nist_webbook._NoRedirect().redirect_request(
            None, None, 302, "redirect", {}, "https://example.com/"
        )
        is None
    )


def test_name_filter_rejects_polymer_transient_and_structural_notation():
    assert not nist_webbook._ordinary_name("Polyoxymethylene")
    assert not nist_webbook._ordinary_name("HCOH (hydroxymethylene)")
    assert not nist_webbook._ordinary_name("(cyc-N=NC)O")
    assert nist_webbook._ordinary_name("methylene chloride")


def test_formula_lookup_uses_persistent_cache_offline(tmp_path):
    first = client(tmp_path, Opener(FORMULA_RESULTS)).lookup_formula("C2H5N")
    offline = Opener(urllib.error.URLError("offline"))
    second = client(tmp_path, offline).lookup_formula("C2H5N")

    assert first["species"] == second["species"]
    assert second["status"] == "cache"
    assert offline.requests == []


def test_mass_lookup_rescues_valid_formula_outside_enumerator_bounds(tmp_path):
    counts = {"C": 41, "H": 84}
    neutral = formula_id.formula_mass(counts)
    mz = neutral + formula_id.PROTON
    assert formula_id.score_peak(mz, 1.0, elements=["C"]) == []
    opener = Opener(MASS_RESULTS, DETAIL)

    result = client(tmp_path, opener).enrich_peak(mz)

    assert result["status"] == "live"
    assert result["mode"] == "mass"
    assert len(opener.requests) == 2
    assert result["candidates"][0]["formula"] == "C41H84"
    assert result["candidates"][0]["formula_source"] == "compound-catalogue"
    assert result["candidates"][0]["nist_webbook"][0]["cas"] == "7098-20-6"


def test_mass_lookup_rejects_detail_without_real_compound_metadata(tmp_path):
    detail = DETAIL.replace(
        b"<li><strong>CAS Registry Number:</strong> 7098-20-6</li>", b""
    )
    counts = {"C": 41, "H": 84}
    mz = formula_id.formula_mass(counts) + formula_id.PROTON

    result = client(tmp_path, Opener(MASS_RESULTS, detail)).enrich_peak(mz)

    assert result["candidates"] == []


def test_enrichment_adds_proposals_without_changing_ptr_name(tmp_path):
    candidate = formula_id.score_peak(44.0495, 1.0)[0]
    assert candidate["formula"] == "C2H5N"
    candidate.pop("name", None)
    candidate.pop("preferred_name", None)
    candidate.pop("names", None)

    result = client(tmp_path, Opener(FORMULA_RESULTS)).enrich_peak(
        44.0495, candidates=[candidate]
    )

    enriched = result["candidates"][0]
    assert enriched.get("name") is None
    assert [item["name"] for item in enriched["nist_webbook"]] == ["Ethylenimine"]


def test_public_fetch_rejects_external_and_robots_disallowed_urls(tmp_path):
    opener = Opener(FORMULA_RESULTS)
    webbook = client(tmp_path, opener)

    for url in (
        "https://example.com/cgi/cbook.cgi?ID=C75070",
        "https://webbook.nist.gov/cdn-cgi/example",
    ):
        try:
            webbook.fetch_page(url)
        except nist_webbook.WebBookError:
            pass
        else:
            raise AssertionError(f"unsafe URL was accepted: {url}")

    assert opener.requests == []


def test_rate_limiter_is_host_wide_across_different_crawl_states(tmp_path):
    schedule = tmp_path / "host-schedule.sqlite3"
    first = nist_webbook.WebBookClient(
        cache_path=tmp_path / "first-state.sqlite3",
        schedule_path=schedule,
        opener=Opener(FORMULA_RESULTS),
        clock=lambda: 10.0,
        sleeper=lambda _delay: None,
        crawl_delay=5,
    )
    slept = []
    second = nist_webbook.WebBookClient(
        cache_path=tmp_path / "second-state.sqlite3",
        schedule_path=schedule,
        opener=Opener(FORMULA_RESULTS),
        clock=lambda: 12.0,
        sleeper=slept.append,
        crawl_delay=5,
    )

    first.lookup_formula("C2H5N")
    second.lookup_formula("C2H6O")

    assert slept == [3.0]


def test_server_cooldown_extends_shared_host_schedule(tmp_path):
    schedule = tmp_path / "host-schedule.sqlite3"
    first = nist_webbook.WebBookClient(
        cache_path=tmp_path / "first.sqlite3",
        schedule_path=schedule,
        opener=Opener(FORMULA_RESULTS),
        clock=lambda: 10.0,
        sleeper=lambda _delay: None,
        crawl_delay=5,
    )
    first.defer_requests(20)
    slept = []
    second = nist_webbook.WebBookClient(
        cache_path=tmp_path / "second.sqlite3",
        schedule_path=schedule,
        opener=Opener(FORMULA_RESULTS),
        clock=lambda: 12.0,
        sleeper=slept.append,
        crawl_delay=5,
    )

    second.lookup_formula("C2H6O")

    assert slept == [18.0]


def test_rate_limiter_waits_between_request_starts(tmp_path):
    times = iter([10.0, 12.0])
    slept = []
    opener = Opener(FORMULA_RESULTS, FORMULA_RESULTS)
    webbook = nist_webbook.WebBookClient(
        cache_path=tmp_path / "webbook.sqlite3",
        schedule_path=tmp_path / "schedule.sqlite3",
        opener=opener,
        clock=lambda: next(times),
        sleeper=slept.append,
        crawl_delay=5,
    )

    webbook.lookup_formula("C2H5N")
    webbook.lookup_formula("C2H6O")

    assert slept == [3.0]


def test_unavailable_shared_cache_disables_live_requests(tmp_path):
    blocked = tmp_path / "not-a-directory"
    blocked.write_text("x", encoding="utf-8")
    opener = Opener(FORMULA_RESULTS)
    webbook = nist_webbook.WebBookClient(
        cache_path=tmp_path / "webbook.sqlite3",
        schedule_path=blocked / "schedule.sqlite3",
        opener=opener,
        crawl_delay=0,
    )

    result = webbook.lookup_formula("C2H5N")

    assert result["status"] == "unavailable"
    assert opener.requests == []


def test_unrecognised_or_offline_response_fails_open(tmp_path):
    malformed = client(tmp_path, Opener(b"<html>changed</html>"))
    offline = client(tmp_path, Opener(urllib.error.URLError("offline")))
    detail_offline = client(
        tmp_path / "detail-offline",
        Opener(MASS_RESULTS, urllib.error.URLError("offline")),
    )
    counts = {"C": 41, "H": 84}
    mz = formula_id.formula_mass(counts) + formula_id.PROTON

    assert malformed.lookup_formula("C2H5N")["status"] == "unavailable"
    assert offline.enrich_peak(44.0495)["candidates"] == []
    assert detail_offline.enrich_peak(mz)["candidates"] == []
