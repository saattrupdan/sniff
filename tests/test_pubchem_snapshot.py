"""Tests for bounded, offline PubChem enrichment snapshots."""

import importlib.util
import json
from pathlib import Path


def _module():
    path = Path(__file__).parents[1] / "scripts" / "fetch_pubchem_enrichment.py"
    spec = importlib.util.spec_from_file_location("fetch_pubchem_enrichment", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_snapshot_is_bounded_provenanced_and_resumable(tmp_path, monkeypatch):
    module = _module()
    calls = []

    def request(url):
        calls.append(url)
        if "/cids/" in url:
            return {"IdentifierList": {"CID": [180, 7858]}}, "cid-sha"
        return {
            "PropertyTable": {
                "Properties": [
                    {
                        "CID": 180,
                        "MolecularFormula": "C3H6O",
                        "Title": "Acetone",
                        "IUPACName": "propan-2-one",
                        "InChI": "InChI=1S/C3H6O/c1-3(2)4/h1-2H3",
                        "InChIKey": "CSCPPACGZOOCGX-UHFFFAOYSA-N",
                        "ConnectivitySMILES": "CC(=O)C",
                        "SMILES": "CC(=O)C",
                    },
                    {
                        "CID": 7858,
                        "MolecularFormula": "C3H6O",
                        "Title": "Allyl alcohol",
                    },
                ]
            }
        }, "property-sha"

    monkeypatch.setattr(module, "_request_json", request)
    monkeypatch.setattr(module.time, "sleep", lambda _seconds: None)
    output = tmp_path / "pubchem.json"
    arguments = [
        "--formula",
        "C3H6O",
        "--output",
        str(output),
        "--max-compounds-per-formula",
        "2",
    ]

    module.main(arguments)
    first = json.loads(output.read_text(encoding="utf-8"))
    module.main(arguments)
    second = json.loads(output.read_text(encoding="utf-8"))
    module.main(
        [
            "--formula",
            "C2H6O",
            "--output",
            str(output),
            "--max-compounds-per-formula",
            "2",
        ]
    )
    module.main(arguments)
    extended = json.loads(output.read_text(encoding="utf-8"))

    assert len(calls) == 4
    assert first == second
    assert first["completed_formula_count"] == 1
    assert first["queries"][0]["cid_response_sha256"] == "cid-sha"
    assert first["queries"][0]["property_response_sha256"] == "property-sha"
    assert first["compounds"][0]["canonical_smiles"] == "CC(=O)C"
    assert first["source"] == "PubChem PUG REST"
    assert extended["formula_count"] == 2
    assert extended["completed_formula_count"] == 2
    assert {item["formula"] for item in extended["queries"]} == {"C2H6O", "C3H6O"}
