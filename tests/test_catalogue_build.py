import importlib.util
import json
import sqlite3
from pathlib import Path

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
