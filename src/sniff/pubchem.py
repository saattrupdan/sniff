"""Offline access to bounded PubChem structure proposals.

PubChem records enrich formula-supported candidates only. They are not PTR-MS
identifications and never participate in formula scoring or assignment eligibility.
"""

from __future__ import annotations

import json
from functools import lru_cache
from importlib import resources
from pathlib import Path

RESOURCE = "pubchem_enrichment.json"
SCHEMA_VERSION = 1


def default_snapshot_path():
    """Return the bundled PubChem snapshot path."""
    return Path(resources.files("sniff.reference").joinpath(RESOURCE))


@lru_cache(maxsize=4)
def _load(path):
    source = Path(path)
    if not source.is_file():
        return {"metadata": {}, "compounds": []}
    try:
        data = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError):
        return {"metadata": {}, "compounds": []}
    if data.get("schema_version") != SCHEMA_VERSION:
        return {"metadata": {}, "compounds": []}
    compounds = data.get("compounds")
    if not isinstance(compounds, list):
        return {"metadata": {}, "compounds": []}
    return {
        "metadata": {
            key: value
            for key, value in data.items()
            if key not in {"compounds", "queries"}
        },
        "compounds": compounds,
    }


def lookup_formula(formula, *, path=None, limit=None):
    """Return source-labelled PubChem proposals for an exact formula."""
    snapshot = _load(str(path or default_snapshot_path()))
    matches = []
    for item in snapshot["compounds"]:
        if item.get("formula") != formula:
            continue
        cid = int(item["cid"])
        matches.append(
            {
                "formula": formula,
                "exact_mass": None,
                "nist_id": None,
                "name": item.get("name") or item.get("iupac_name") or f"CID {cid}",
                "cas": None,
                "inchi": item.get("inchi"),
                "inchi_key": item.get("inchi_key"),
                "url": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}",
                "source": "PubChem",
                "pubchem_cid": cid,
                "iupac_name": item.get("iupac_name"),
                "canonical_smiles": item.get("canonical_smiles"),
                "isomeric_smiles": item.get("isomeric_smiles"),
            }
        )
        if limit is not None and len(matches) >= max(0, int(limit)):
            break
    return matches


def metadata(*, path=None):
    """Return provenance for the bundled or supplied offline snapshot."""
    return dict(_load(str(path or default_snapshot_path()))["metadata"])
