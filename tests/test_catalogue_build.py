import importlib.util
import json
import sqlite3
from pathlib import Path

import pytest

from sniff.catalogue import CompoundCatalogue


def _builder_module():
    path = Path(__file__).parents[1] / "scripts" / "build_compound_catalogue.py"
    spec = importlib.util.spec_from_file_location("build_compound_catalogue", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_export_strips_crawl_state_and_keeps_normalised_metadata(tmp_path):
    builder = _builder_module()
    state_path = tmp_path / "state.sqlite3"
    output_path = tmp_path / "catalogue.sqlite3"
    connection = builder._connect_state(state_path)
    connection.execute(
        """
        INSERT INTO formula_job(
            formula, seed, status, attempts, response_sha, updated_at
        ) VALUES ('C2H4O', 'ptr-library', 'done', 1, 'page-hash', 1)
        """
    )
    detail = {
        "nist_id": "C75070",
        "name": "Acetaldehyde",
        "formula": "C2H4O",
        "cas": "75-07-0",
        "inchi": "InChI=1S/C2H4O/c1-2-3/h2H,1H3",
        "inchi_key": "IKHGUXGNUITLKF-UHFFFAOYSA-N",
        "url": "https://webbook.nist.gov/cgi/cbook.cgi?ID=C75070",
    }
    connection.execute(
        """
        INSERT INTO species_record(
            nist_id, formula, name, url, detail_status, detail_attempts,
            detail_json, response_sha, detail_requested, updated_at
        ) VALUES (?, ?, ?, ?, 'done', 1, ?, 'detail-hash', 1, 1)
        """,
        (
            detail["nist_id"],
            detail["formula"],
            detail["name"],
            detail["url"],
            json.dumps(detail),
        ),
    )
    connection.commit()
    connection.close()

    builder.export_catalogue(state_path, output_path)

    compounds = CompoundCatalogue(output_path)
    assert compounds.lookup_formula("C2H4O")[0]["inchi_key"] == detail["inchi_key"]
    with sqlite3.connect(output_path) as exported:
        tables = {
            row[0]
            for row in exported.execute(
                "SELECT name FROM sqlite_master WHERE type = 'table'"
            )
        }
        assert tables == {"metadata", "formula", "species"}
        assert exported.execute("PRAGMA integrity_check").fetchone()[0] == "ok"


def test_unavailable_formula_lookup_remains_retryable(tmp_path, monkeypatch):
    builder = _builder_module()

    class UnavailableClient:
        def __init__(self, **_kwargs):
            pass

        def lookup_formula(self, _formula):
            return {"status": "unavailable", "species": [], "excluded": 0}

    monkeypatch.setattr(builder, "_ptr_formulas", lambda: ["C2H4O"])
    monkeypatch.setattr(builder.nist_webbook, "WebBookClient", UnavailableClient)
    state_path = tmp_path / "state.sqlite3"

    builder.crawl(state_path)

    with sqlite3.connect(state_path) as connection:
        status, attempts = connection.execute(
            "SELECT status, attempts FROM formula_job WHERE formula = 'C2H4O'"
        ).fetchone()
    assert status == "error"
    assert attempts == 1


def test_ptr_detail_seed_marks_only_one_representative_per_formula(tmp_path):
    builder = _builder_module()
    connection = builder._connect_state(tmp_path / "state.sqlite3")

    builder._seed_ptr_details(connection)

    rows = connection.execute(
        """
        SELECT formula, COUNT(*)
        FROM species_record
        WHERE detail_requested = 1
        GROUP BY formula
        """
    ).fetchall()
    assert rows
    assert all(count == 1 for _formula, count in rows)
    connection.close()


def _full_crawl_state(path, *, pending=False, stale=False):
    connection = sqlite3.connect(path)
    connection.executescript(
        """
        CREATE TABLE crawl_meta (key TEXT PRIMARY KEY, value TEXT NOT NULL);
        CREATE TABLE generation (
            manifest_sha TEXT PRIMARY KEY,
            index_sha TEXT NOT NULL,
            discovered_at REAL NOT NULL,
            species_count INTEGER NOT NULL,
            complete INTEGER NOT NULL
        );
        CREATE TABLE generation_species (
            manifest_sha TEXT NOT NULL,
            nist_id TEXT NOT NULL,
            PRIMARY KEY (manifest_sha, nist_id)
        );
        CREATE TABLE species_job (
            nist_id TEXT PRIMARY KEY,
            url TEXT NOT NULL,
            status TEXT NOT NULL,
            detail_json TEXT,
            parser_version INTEGER
        );
        INSERT INTO crawl_meta VALUES ('active_manifest', 'manifest');
        INSERT INTO generation VALUES ('manifest', 'index', 1, 7, 1);
        """
    )
    rows = [
        (
            "C75070",
            "https://webbook.nist.gov/cgi/cbook.cgi?ID=C75070&Units=SI",
            "done",
            {
                "name": "Acetaldehyde",
                "formula": "C2H4O",
                "cas": "75-07-0",
                "inchi": "InChI=1S/C2H4O/c1-2-3/h2H,1H3",
                "inchi_key": "IKHGUXGNUITLKF-UHFFFAOYSA-N",
                "webbook_id": "C75070",
            },
        ),
        (
            "Uhash",
            "https://webbook.nist.gov/cgi/inchi/InChI%3D1S/C41H84",
            "done",
            {
                "name": "Hentetracontane",
                "formula": "C41H84",
                "inchi": "InChI=1S/C41H84",
                "inchi_key": "TESTKEY-UHFFFAOYSA-N",
                "webbook_id": None,
            },
        ),
        (
            "Bion",
            "https://webbook.nist.gov/cgi/cbook.cgi?ID=Bion&Units=SI",
            "done",
            {"name": "iron oxide anion", "formula": "FeO-", "webbook_id": "Bion"},
        ),
        (
            "Bradical",
            "https://webbook.nist.gov/cgi/cbook.cgi?ID=Bradical&Units=SI",
            "done",
            {"name": "EtN radical", "formula": "C2H5N", "webbook_id": "Bradical"},
        ),
        (
            "Bimplausible",
            "https://webbook.nist.gov/cgi/cbook.cgi?ID=Bimplausible&Units=SI",
            "done",
            {"name": "Impossible", "formula": "CH100", "webbook_id": "Bimplausible"},
        ),
        (
            "missing",
            "https://webbook.nist.gov/cgi/cbook.cgi?ID=missing&Units=SI",
            "missing",
            None,
        ),
        (
            "unusable",
            "https://webbook.nist.gov/cgi/cbook.cgi?ID=unusable&Units=SI",
            "unusable",
            None,
        ),
    ]
    if pending:
        rows[-1] = (*rows[-1][:2], "pending", None)
    connection.executemany(
        "INSERT INTO generation_species VALUES ('manifest', ?)",
        [(row[0],) for row in rows],
    )
    connection.executemany(
        "INSERT INTO species_job VALUES (?, ?, ?, ?, ?)",
        [
            (
                nist_id,
                url,
                status,
                json.dumps(detail) if detail is not None else None,
                (
                    None
                    if stale and nist_id == "C75070"
                    else 2
                    if status in ("done", "unusable")
                    else None
                ),
            )
            for nist_id, url, status, detail in rows
        ],
    )
    connection.commit()
    connection.close()


def test_full_export_validates_and_classifies_complete_manifest(tmp_path):
    builder = _builder_module()
    state_path = tmp_path / "full.sqlite3"
    output_path = tmp_path / "catalogue.sqlite3"
    _full_crawl_state(state_path)

    builder.export_full_catalogue(state_path, output_path)

    compounds = CompoundCatalogue(output_path)
    metadata = compounds.metadata()
    assert metadata["schema_version"] == "2"
    assert metadata["source_manifest_sha"] == "manifest"
    assert metadata["source_species_count"] == "7"
    assert metadata["source_done_count"] == "5"
    assert metadata["source_missing_count"] == "1"
    assert metadata["source_unusable_count"] == "1"
    assert metadata["excluded_unsupported_formula_count"] == "1"
    assert metadata["excluded_non_ordinary_name_count"] == "1"
    assert metadata["excluded_implausible_formula_count"] == "1"
    assert compounds.search("Acetaldehyde")[0]["nist_id"] == "C75070"
    inchi_only = compounds.search("Hentetracontane")[0]
    assert inchi_only["nist_id"] is None
    assert inchi_only["url"].startswith("https://webbook.nist.gov/cgi/inchi/")
    with sqlite3.connect(output_path) as exported:
        tables = {
            row[0]
            for row in exported.execute(
                "SELECT name FROM sqlite_master WHERE type = 'table'"
            )
        }
        assert tables == {"metadata", "formula", "species"}
        assert exported.execute("PRAGMA integrity_check").fetchone()[0] == "ok"
        assert exported.execute("PRAGMA foreign_key_check").fetchall() == []


def test_full_export_refuses_stale_parser_provenance(tmp_path):
    builder = _builder_module()
    state_path = tmp_path / "full.sqlite3"
    _full_crawl_state(state_path, stale=True)

    with pytest.raises(RuntimeError, match="stale parser"):
        builder.export_full_catalogue(state_path, tmp_path / "catalogue.sqlite3")


def test_full_export_refuses_unresolved_jobs_without_replacing_output(tmp_path):
    builder = _builder_module()
    state_path = tmp_path / "full.sqlite3"
    output_path = tmp_path / "catalogue.sqlite3"
    output_path.write_bytes(b"existing")
    _full_crawl_state(state_path, pending=True)

    with pytest.raises(RuntimeError, match="pending=1"):
        builder.export_full_catalogue(state_path, output_path)

    assert output_path.read_bytes() == b"existing"
