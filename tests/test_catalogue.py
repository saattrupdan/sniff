import sqlite3

import pytest

from sniff import catalogue, formula_id, isotopes, nist_webbook


def _catalogue(path):
    connection = sqlite3.connect(path)
    connection.executescript(
        """
        CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);
        CREATE TABLE formula (
            id INTEGER PRIMARY KEY,
            formula TEXT NOT NULL UNIQUE,
            exact_mass REAL NOT NULL,
            ptr_seed INTEGER NOT NULL
        );
        CREATE TABLE species (
            id INTEGER PRIMARY KEY,
            formula_id INTEGER NOT NULL REFERENCES formula(id),
            nist_id TEXT NOT NULL UNIQUE,
            name TEXT NOT NULL,
            cas TEXT,
            inchi TEXT,
            inchi_key TEXT,
            url TEXT NOT NULL
        );
        INSERT INTO metadata VALUES ('schema_version', '1');
        INSERT INTO formula VALUES (1, 'C2H4O', 44.02621474849, 1);
        INSERT INTO species VALUES (
            1, 1, 'C75070', 'Acetaldehyde', '75-07-0',
            'InChI=1S/C2H4O/c1-2-3/h2H,1H3', 'IKHGUXGNUITLKF-UHFFFAOYSA-N',
            'https://webbook.nist.gov/cgi/cbook.cgi?ID=C75070'
        );
        """
    )
    connection.commit()
    connection.close()
    return catalogue.CompoundCatalogue(path)


def test_missing_catalogue_fails_open_without_creating_file(tmp_path):
    path = tmp_path / "missing.sqlite3"
    compounds = catalogue.CompoundCatalogue(path)

    assert compounds.metadata() == {}
    assert compounds.lookup_formula("C2H4O") == []
    assert compounds.formulas_in_mass_range(44.0262) == []
    assert not path.exists()


def test_formula_mass_and_metadata_queries_use_read_only_catalogue(tmp_path):
    compounds = _catalogue(tmp_path / "catalogue.sqlite3")

    assert compounds.metadata()["schema_version"] == "1"
    result = compounds.lookup_formula("C2H4O")
    assert result[0]["name"] == "Acetaldehyde"
    assert result[0]["cas"] == "75-07-0"
    assert result[0]["inchi_key"] == "IKHGUXGNUITLKF-UHFFFAOYSA-N"
    assert compounds.search("acetald")[0]["formula"] == "C2H4O"
    assert compounds.search("75-07-0")[0]["nist_id"] == "C75070"
    assert compounds.search("InChI=1S/C2H4O")[0]["name"] == "Acetaldehyde"


def test_mass_lookup_recomputes_formula_mass_locally(tmp_path):
    compounds = _catalogue(tmp_path / "catalogue.sqlite3")
    exact = formula_id.formula_mass(isotopes.parse_formula("C2H4O"))

    assert compounds.formulas_in_mass_range(exact, 0.0001) == [
        {"formula": "C2H4O", "exact_mass": pytest.approx(exact)}
    ]
    assert compounds.formulas_in_mass_range(exact + 0.0121, 0.012) == []


def test_candidate_enrichment_keeps_ptr_identity_and_score(tmp_path):
    compounds = _catalogue(tmp_path / "catalogue.sqlite3")
    candidates = [{"formula": "C2H4O", "name": "acetaldehyde", "score": 0.7}]

    enriched = compounds.enrich_candidates(candidates)

    assert enriched[0]["name"] == "acetaldehyde"
    assert enriched[0]["score"] == 0.7
    assert enriched[0]["catalogue"][0]["name"] == "Acetaldehyde"
    assert "catalogue" not in candidates[0]


def test_bundled_catalogue_has_versioned_ptr_seed_and_searchable_metadata():
    compounds = catalogue.CompoundCatalogue()

    metadata = compounds.metadata()
    assert metadata["schema_version"] == "1"
    assert int(metadata["formula_count"]) >= 269
    assert int(metadata["species_count"]) >= 10_900
    acetaldehyde = compounds.search("IKHGUXGNUITLKF-UHFFFAOYSA-N")[0]
    assert acetaldehyde["name"] == "Acetaldehyde"
    assert acetaldehyde["formula"] == "C2H4O"
    assert acetaldehyde["cas"] == "75-07-0"


def test_bundled_catalogue_contains_only_filtered_ordinary_names():
    compounds = catalogue.CompoundCatalogue()
    with sqlite3.connect(compounds.path) as connection:
        names = [row[0] for row in connection.execute("SELECT name FROM species")]

    assert len(names) >= 10_900
    assert all(nist_webbook._ordinary_name(name) for name in names)


def test_cas_is_derived_only_from_valid_webbook_identifier():
    assert catalogue.cas_from_nist_id("C75070") == "75-07-0"
    assert catalogue.cas_from_nist_id("C7732185") == "7732-18-5"
    assert catalogue.cas_from_nist_id("C75071") is None
    assert catalogue.cas_from_nist_id("InChIKey") is None
