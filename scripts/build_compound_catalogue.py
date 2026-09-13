#!/usr/bin/env python3
"""Build Sniff's bundled PTR-focused NIST Chemistry WebBook catalogue.

The maintainer crawl is deliberately separate from Sniff's runtime. Its SQLite state is
resumable and retains request hashes/checkpoints; the exported package database contains
only normalised compound metadata and no HTML or spectra.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sqlite3
import time
from pathlib import Path

from sniff import formula_id, isotopes, nist_webbook, ptrms
from sniff.catalogue import (
    CATALOGUE_SCHEMA_VERSION,
    canonical_formula,
    cas_from_nist_id,
)

DEFAULT_STATE = Path.home() / ".sniff" / "compound-catalogue-build-state.sqlite3"
DEFAULT_OUTPUT = Path("src/sniff/reference/compound_catalogue.sqlite3")


def _connect_state(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(str(path), timeout=10.0)
    connection.row_factory = sqlite3.Row
    connection.executescript(
        """
        PRAGMA foreign_keys = ON;
        CREATE TABLE IF NOT EXISTS formula_job (
            formula TEXT PRIMARY KEY,
            seed TEXT NOT NULL,
            status TEXT NOT NULL DEFAULT 'pending',
            attempts INTEGER NOT NULL DEFAULT 0,
            last_error TEXT,
            response_sha TEXT,
            updated_at REAL
        );
        CREATE TABLE IF NOT EXISTS species_record (
            nist_id TEXT PRIMARY KEY,
            formula TEXT NOT NULL,
            name TEXT NOT NULL,
            url TEXT NOT NULL,
            detail_status TEXT NOT NULL DEFAULT 'pending',
            detail_attempts INTEGER NOT NULL DEFAULT 0,
            detail_error TEXT,
            detail_json TEXT,
            response_sha TEXT,
            detail_requested INTEGER NOT NULL DEFAULT 0,
            search_seen INTEGER NOT NULL DEFAULT 0,
            updated_at REAL
        );
        CREATE INDEX IF NOT EXISTS species_record_formula
            ON species_record(formula);
        """
    )
    columns = {
        row[1] for row in connection.execute("PRAGMA table_info(species_record)")
    }
    if "detail_requested" not in columns:
        connection.execute(
            "ALTER TABLE species_record "
            "ADD COLUMN detail_requested INTEGER NOT NULL DEFAULT 0"
        )
    if "search_seen" not in columns:
        connection.execute(
            "ALTER TABLE species_record "
            "ADD COLUMN search_seen INTEGER NOT NULL DEFAULT 0"
        )
    connection.commit()
    return connection


def _ptr_formulas():
    compounds = (ptrms.load_rate_constants() or {}).get("compounds", [])
    return sorted({canonical_formula(row["formula"]) for row in compounds})


def _expanded_formulas(elements, max_neutral_mass):
    """Enumerate the selected plausible formula space in narrow mass slices."""
    found = set()
    step = 0.2
    tolerance = step / 2.0 + 1e-9
    for index in range(int(math.ceil(max_neutral_mass / step)) + 1):
        centre = index * step
        for counts, mass in formula_id.enumerate_formulas(
            centre,
            tolerance,
            elements=[element for element in elements if element != "H"],
        ):
            if 0 < mass <= max_neutral_mass:
                found.add(formula_id.formula_str(counts))
    return sorted(found)


def seed_jobs(connection, *, expand_elements=None, max_neutral_mass=300.0):
    rows = [(formula, "ptr-library") for formula in _ptr_formulas()]
    if expand_elements:
        rows.extend(
            (formula, "expanded-" + expand_elements)
            for formula in _expanded_formulas(expand_elements, max_neutral_mass)
        )
    connection.executemany(
        "INSERT OR IGNORE INTO formula_job(formula, seed) VALUES (?, ?)", rows
    )
    connection.commit()
    return connection.execute("SELECT COUNT(*) FROM formula_job").fetchone()[0]


def _cache_hash(connection, key):
    try:
        row = connection.execute(
            "SELECT response_sha FROM webbook_cache WHERE key = ?", (key,)
        ).fetchone()
    except sqlite3.OperationalError:
        return None
    return row[0] if row else None


def _seed_ptr_details(connection):
    representatives = {}
    for compound in (ptrms.load_rate_constants() or {}).get("compounds", []):
        formula = compound.get("formula")
        cas = compound.get("cas")
        nist_id = "C" + str(cas or "").replace("-", "")
        valid_cas = cas_from_nist_id(nist_id)
        if formula and cas and valid_cas == cas and formula not in representatives:
            representatives[formula] = compound
    for formula, compound in representatives.items():
        nist_id = "C" + compound["cas"].replace("-", "")
        connection.execute(
            """
            INSERT INTO species_record(
                nist_id, formula, name, url, detail_requested, updated_at
            ) VALUES (?, ?, ?, ?, 1, ?)
            ON CONFLICT(nist_id) DO UPDATE SET detail_requested=1
            """,
            (
                nist_id,
                canonical_formula(formula),
                compound["name"],
                f"https://webbook.nist.gov/cgi/cbook.cgi?ID={nist_id}",
                time.time(),
            ),
        )
    connection.commit()


def crawl(
    state_path,
    *,
    expand_elements=None,
    max_neutral_mass=300.0,
    details=False,
    ptr_details=False,
):
    connection = _connect_state(state_path)
    count = seed_jobs(
        connection,
        expand_elements=expand_elements,
        max_neutral_mass=max_neutral_mass,
    )
    print(f"Catalogue queue: {count} formulae", flush=True)
    if ptr_details:
        _seed_ptr_details(connection)
    client = nist_webbook.WebBookClient(cache_path=state_path, timeout=20.0)
    jobs = connection.execute(
        "SELECT formula FROM formula_job WHERE status != 'done' ORDER BY formula"
    ).fetchall()
    for index, job in enumerate(jobs, 1):
        formula = job["formula"]
        now = time.time()
        try:
            result = client.lookup_formula(formula)
            if result.get("status") == "unavailable":
                raise nist_webbook.WebBookError(
                    "the formula lookup was unavailable and remains queued for retry"
                )
            for species in result["species"]:
                connection.execute(
                    """
                    INSERT INTO species_record(
                        nist_id, formula, name, url, search_seen, updated_at
                    ) VALUES (?, ?, ?, ?, 1, ?)
                    ON CONFLICT(nist_id) DO UPDATE SET
                        formula=excluded.formula, name=excluded.name,
                        url=excluded.url, search_seen=1, updated_at=excluded.updated_at
                    """,
                    (
                        species["nist_id"],
                        canonical_formula(species["formula"]),
                        " ".join(species["name"].split()),
                        species["url"],
                        now,
                    ),
                )
            connection.execute(
                """
                UPDATE formula_job
                SET status='done', attempts=attempts+1, last_error=NULL,
                    response_sha=?, updated_at=?
                WHERE formula=?
                """,
                (_cache_hash(connection, f"formula:{formula}"), now, formula),
            )
            connection.commit()
            print(
                f"formula {index}/{len(jobs)} {formula}: "
                f"{len(result['species'])} retained, {result['excluded']} excluded",
                flush=True,
            )
        except Exception as exc:
            connection.execute(
                """
                UPDATE formula_job
                SET status='error', attempts=attempts+1, last_error=?, updated_at=?
                WHERE formula=?
                """,
                (str(exc), now, formula),
            )
            connection.commit()
            print(f"formula {formula}: {exc}", flush=True)

    if details or ptr_details:
        detail_filter = "" if details else "AND detail_requested = 1"
        rows = connection.execute(
            f"""
            SELECT nist_id FROM species_record
            WHERE detail_status NOT IN ('done', 'rejected') {detail_filter}
            ORDER BY nist_id
            """
        ).fetchall()
        for index, row in enumerate(rows, 1):
            nist_id = row["nist_id"]
            now = time.time()
            try:
                detail, _status = client.lookup_species(nist_id)
                formula = canonical_formula(detail["formula"])
                name = " ".join(detail["name"].split())
                if not nist_webbook._ordinary_name(name):
                    raise ValueError("species name failed the ordinary-compound filter")
                connection.execute(
                    """
                    UPDATE species_record
                    SET formula=?, name=?, url=?, detail_status='done',
                        detail_attempts=detail_attempts+1, detail_error=NULL,
                        detail_json=?, response_sha=?, updated_at=?
                    WHERE nist_id=?
                    """,
                    (
                        formula,
                        name,
                        detail["url"],
                        json.dumps(detail, sort_keys=True),
                        _cache_hash(connection, f"species:{nist_id}"),
                        now,
                        nist_id,
                    ),
                )
                connection.commit()
                print(f"detail {index}/{len(rows)} {nist_id}: {name}", flush=True)
            except Exception as exc:
                status = "rejected" if isinstance(exc, ValueError) else "error"
                connection.execute(
                    """
                    UPDATE species_record
                    SET detail_status=?, detail_attempts=detail_attempts+1,
                        detail_error=?, updated_at=?
                    WHERE nist_id=?
                    """,
                    (status, str(exc), now, nist_id),
                )
                connection.commit()
                print(f"detail {nist_id}: {exc}", flush=True)
    connection.close()


def _create_output(path):
    if path.exists():
        path.unlink()
    path.parent.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(str(path))
    connection.executescript(
        """
        PRAGMA foreign_keys = ON;
        CREATE TABLE metadata (
            key TEXT PRIMARY KEY,
            value TEXT NOT NULL
        );
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
        CREATE INDEX formula_exact_mass ON formula(exact_mass);
        CREATE INDEX species_formula_id ON species(formula_id);
        CREATE INDEX species_name ON species(name COLLATE NOCASE);
        CREATE INDEX species_cas ON species(cas);
        CREATE INDEX species_inchi ON species(inchi);
        CREATE INDEX species_inchi_key ON species(inchi_key);
        """
    )
    return connection


def export_catalogue(state_path, output_path):
    state = _connect_state(state_path)
    output = _create_output(output_path)
    ptr_formulas = set(_ptr_formulas())
    rows = state.execute(
        """
        SELECT nist_id, formula, name, url, detail_json
        FROM species_record
        WHERE detail_status != 'rejected'
          AND (search_seen = 1 OR detail_status = 'done')
        ORDER BY formula, name COLLATE NOCASE, nist_id
        """
    ).fetchall()
    species = []
    formulas = {
        canonical_formula(row["formula"])
        for row in state.execute(
            "SELECT formula FROM formula_job WHERE status = 'done' ORDER BY formula"
        )
    }
    for row in rows:
        try:
            formula = canonical_formula(row["formula"])
            counts = isotopes.parse_formula(formula)
        except ValueError:
            continue
        name = " ".join(row["name"].split())
        if not nist_webbook._ordinary_name(name):
            continue
        detail = json.loads(row["detail_json"]) if row["detail_json"] else {}
        detail_formula = detail.get("formula")
        if detail_formula and canonical_formula(detail_formula) != formula:
            continue
        formulas.add(formula)
        species.append(
            {
                "nist_id": row["nist_id"],
                "formula": formula,
                "name": name,
                "cas": detail.get("cas") or cas_from_nist_id(row["nist_id"]),
                "inchi": detail.get("inchi"),
                "inchi_key": detail.get("inchi_key"),
                "url": row["url"],
                "exact_mass": formula_id.formula_mass(counts),
            }
        )
    formula_rows = []
    for formula in sorted(formulas):
        mass = formula_id.formula_mass(isotopes.parse_formula(formula))
        formula_rows.append((formula, mass, int(formula in ptr_formulas)))
    output.executemany(
        "INSERT INTO formula(formula, exact_mass, ptr_seed) VALUES (?, ?, ?)",
        formula_rows,
    )
    formula_ids = {
        formula: formula_id
        for formula_id, formula in output.execute("SELECT id, formula FROM formula")
    }
    output.executemany(
        """
        INSERT INTO species(formula_id, nist_id, name, cas, inchi, inchi_key, url)
        VALUES (?, ?, ?, ?, ?, ?, ?)
        """,
        [
            (
                formula_ids[row["formula"]],
                row["nist_id"],
                row["name"],
                row["cas"],
                row["inchi"],
                row["inchi_key"],
                row["url"],
            )
            for row in species
        ],
    )
    now = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    metadata = {
        "schema_version": str(CATALOGUE_SCHEMA_VERSION),
        "catalogue_version": "1",
        "webbook_parser_version": str(nist_webbook.PARSER_VERSION),
        "source": "NIST Chemistry WebBook, SRD 69",
        "source_url": "https://webbook.nist.gov/",
        "built_at": now,
        "formula_count": str(len(formula_rows)),
        "species_count": str(len(species)),
        "scope": "PTR Library formula families plus maintainer-selected expansions",
    }
    output.executemany(
        "INSERT INTO metadata(key, value) VALUES (?, ?)", sorted(metadata.items())
    )
    output.commit()
    integrity = output.execute("PRAGMA integrity_check").fetchone()[0]
    output.close()
    state.close()
    if integrity != "ok":
        raise RuntimeError(f"catalogue integrity check failed: {integrity}")
    digest = hashlib.sha256(output_path.read_bytes()).hexdigest()
    print(
        f"Wrote {output_path}: {len(formula_rows)} formulae, {len(species)} species, "
        f"sha256={digest}",
        flush=True,
    )


def _parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("crawl", "export", "build"))
    parser.add_argument("--state", type=Path, default=DEFAULT_STATE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--expand-elements",
        choices=("CHNO", "CHNOS"),
        help="also queue the plausible formula space for these elements",
    )
    parser.add_argument("--max-neutral-mass", type=float, default=300.0)
    parser.add_argument(
        "--details",
        action="store_true",
        help="crawl every retained species page for InChI metadata (very slow)",
    )
    parser.add_argument(
        "--ptr-details",
        action="store_true",
        help="crawl one PTR Library representative per formula for InChI metadata",
    )
    return parser


def main(argv=None):
    args = _parser().parse_args(argv)
    if args.command in ("crawl", "build"):
        crawl(
            args.state,
            expand_elements=args.expand_elements,
            max_neutral_mass=args.max_neutral_mass,
            details=args.details,
            ptr_details=args.ptr_details,
        )
    if args.command in ("export", "build"):
        export_catalogue(args.state, args.output)


if __name__ == "__main__":
    main()
