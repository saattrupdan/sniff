#!/usr/bin/env python3
"""Build a bounded, reproducible PubChem enrichment snapshot via PUG REST."""

import argparse
import hashlib
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

from sniff.catalogue import canonical_formula

API = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PROPERTIES = (
    "Title,IUPACName,MolecularFormula,InChI,InChIKey,"
    "CanonicalSMILES,IsomericSMILES"
)
NOTICE_URL = "https://www.ncbi.nlm.nih.gov/home/about/policies/"


def main(argv=None):
    """Fetch selected formula records and atomically write an offline snapshot."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--peaks-json", type=Path)
    parser.add_argument("--formula", action="append", default=[])
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("src/sniff/reference/pubchem_enrichment.json"),
    )
    parser.add_argument("--max-compounds-per-formula", type=int, default=10)
    parser.add_argument("--requests-per-second", type=float, default=4.0)
    args = parser.parse_args(argv)
    formulas = set(args.formula)
    if args.peaks_json:
        formulas.update(_peak_formulas(args.peaks_json))
    formulas = sorted({canonical_formula(value) for value in formulas if value})
    if not formulas:
        parser.error("provide --formula or --peaks-json with formula candidates")
    if args.max_compounds_per_formula < 1:
        parser.error("--max-compounds-per-formula must be positive")
    delay = 1.0 / min(max(args.requests_per_second, 0.1), 5.0)
    retrieved_at = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    compounds = []
    queries = []
    if args.output.is_file():
        previous = json.loads(args.output.read_text(encoding="utf-8"))
        if (
            previous.get("schema_version") == 1
            and previous.get("max_compounds_per_formula")
            == args.max_compounds_per_formula
        ):
            compounds = list(previous.get("compounds") or [])
            queries = list(previous.get("queries") or [])
            retrieved_at = previous.get("retrieved_at") or retrieved_at
            formulas = sorted(
                set(formulas) | {item["formula"] for item in queries}
            )
    completed = {item["formula"] for item in queries}
    for index, formula in enumerate(formulas):
        if formula in completed:
            continue
        print(f"PubChem {index + 1}/{len(formulas)}: {formula}", flush=True)
        if queries:
            time.sleep(delay)
        query_limit = args.max_compounds_per_formula + 1
        cid_url = (
            f"{API}/compound/fastformula/{urllib.parse.quote(formula)}/cids/JSON"
            f"?MaxRecords={query_limit}&MaxSeconds=30"
        )
        cid_data, cid_sha = _request_json(cid_url)
        cids = sorted({int(value) for value in cid_data["IdentifierList"]["CID"]})
        selected = cids[: args.max_compounds_per_formula]
        property_sha = None
        property_url = None
        if selected:
            time.sleep(delay)
            joined = ",".join(str(value) for value in selected)
            property_url = f"{API}/compound/cid/{joined}/property/{PROPERTIES}/JSON"
            property_data, property_sha = _request_json(property_url)
            for item in property_data["PropertyTable"]["Properties"]:
                try:
                    item_formula = canonical_formula(item.get("MolecularFormula"))
                except ValueError:
                    continue
                if item_formula != formula:
                    continue
                compounds.append(_normalise_compound(formula, item))
        queries.append(
            {
                "formula": formula,
                "cid_url": cid_url,
                "cid_response_sha256": cid_sha,
                "returned_cids": len(cids),
                "selected_cids": len(selected),
                "truncated": len(cids) > len(selected),
                "property_url": property_url,
                "property_response_sha256": property_sha,
            }
        )
        _write_snapshot(
            args.output,
            formulas=formulas,
            compounds=compounds,
            queries=queries,
            retrieved_at=retrieved_at,
            max_compounds_per_formula=args.max_compounds_per_formula,
        )
    digest = _write_snapshot(
        args.output,
        formulas=formulas,
        compounds=compounds,
        queries=queries,
        retrieved_at=retrieved_at,
        max_compounds_per_formula=args.max_compounds_per_formula,
    )
    print(
        f"Wrote {args.output}: {len(queries)}/{len(formulas)} formulae, "
        f"{len(compounds)} compounds, sha256={digest}",
        flush=True,
    )


def _peak_formulas(path):
    data = json.loads(path.read_text(encoding="utf-8"))
    formulas = set()
    for peak in data.get("peaks", []):
        if peak.get("formula"):
            formulas.add(peak["formula"])
        top = peak.get("top_candidate")
        if isinstance(top, dict) and top.get("formula"):
            formulas.add(top["formula"])
        for key in ("candidates", "ion_candidates"):
            values = peak.get(key) or []
            if values and isinstance(values[0], dict) and values[0].get("formula"):
                formulas.add(values[0]["formula"])
        for interpretation in peak.get("interpretation_candidates") or []:
            values = interpretation.get("formula_candidates") or []
            if values and isinstance(values[0], dict) and values[0].get("formula"):
                formulas.add(values[0]["formula"])
    return formulas


def _write_snapshot(
    output,
    *,
    formulas,
    compounds,
    queries,
    retrieved_at,
    max_compounds_per_formula,
):
    payload = {
        "schema_version": 1,
        "source": "PubChem PUG REST",
        "source_url": API,
        "retrieved_at": retrieved_at,
        "notice_url": NOTICE_URL,
        "scope": "formula-supported PTR-MS proposals selected from the documented input",
        "max_compounds_per_formula": max_compounds_per_formula,
        "formula_count": len(formulas),
        "completed_formula_count": len(queries),
        "compound_count": len(compounds),
        "truncated_formula_count": sum(item["truncated"] for item in queries),
        "queries": sorted(queries, key=lambda item: item["formula"]),
        "compounds": sorted(compounds, key=lambda item: (item["formula"], item["cid"])),
    }
    serialised = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    temporary = output.with_name(output.name + ".tmp")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary.write_text(serialised, encoding="utf-8")
    temporary.replace(output)
    return hashlib.sha256(serialised.encode("utf-8")).hexdigest()


def _request_json(url):
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "Sniff-PubChem-snapshot/1 (offline enrichment)"},
    )
    for attempt in range(6):
        try:
            with urllib.request.urlopen(request, timeout=90) as response:
                raw = response.read()
            return json.loads(raw), hashlib.sha256(raw).hexdigest()
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                raw = b'{"IdentifierList":{"CID":[]}}'
                return json.loads(raw), hashlib.sha256(raw).hexdigest()
            if exc.code not in {429, 500, 502, 503, 504} or attempt == 5:
                raise
            retry_after = exc.headers.get("Retry-After")
            wait = float(retry_after) if retry_after else min(2**attempt, 30)
            time.sleep(wait)
        except urllib.error.URLError:
            if attempt == 5:
                raise
            time.sleep(min(2**attempt, 30))
    raise RuntimeError("unreachable PubChem retry state")


def _normalise_compound(formula, item):
    cid = int(item["CID"])
    return {
        "cid": cid,
        "formula": formula,
        "name": _text(item.get("Title")),
        "iupac_name": _text(item.get("IUPACName")),
        "inchi": _text(item.get("InChI")),
        "inchi_key": _text(item.get("InChIKey")),
        "canonical_smiles": _text(item.get("ConnectivitySMILES")),
        "isomeric_smiles": _text(item.get("SMILES")),
    }


def _text(value):
    text = " ".join(str(value or "").split())
    return text[:4096] or None


if __name__ == "__main__":
    main()
