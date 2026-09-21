"""Read-only access to Sniff's bundled compound catalogue.

The catalogue contains metadata proposals, not PTR-MS identification evidence. Formulae
are matched by locally computed exact mass, and names never replace PTR Library names or
become automatic structural assignments.
"""

from __future__ import annotations

import sqlite3
from importlib import resources
from pathlib import Path

from . import formula_id, isotopes

CATALOGUE_SCHEMA_VERSION = 2
CATALOGUE_RESOURCE = "compound_catalogue.sqlite3"


def canonical_formula(value):
    """Return *value* in Hill notation, rejecting malformed formulae."""
    return formula_id.formula_str(isotopes.parse_formula(str(value or "")))


def cas_from_nist_id(nist_id):
    """Recover a valid CAS number from the usual WebBook ``C<digits>`` identifier."""
    text = str(nist_id or "")
    if len(text) < 4 or not text.startswith("C") or not text[1:].isdigit():
        return None
    digits = text[1:]
    if len(digits) < 3:
        return None
    body, check = digits[:-1], int(digits[-1])
    total = sum(int(digit) * weight for weight, digit in enumerate(reversed(body), 1))
    if total % 10 != check:
        return None
    head = body[:-2] or "0"
    return f"{int(head)}-{body[-2:]}-{check}"


def default_catalogue_path():
    """Return the installed catalogue path."""
    return Path(
        str(resources.files("sniff").joinpath("reference").joinpath(CATALOGUE_RESOURCE))
    )


class CompoundCatalogue:
    """Small query facade over the versioned, bundled SQLite catalogue."""

    def __init__(self, path=None):
        self.path = Path(path) if path is not None else default_catalogue_path()

    @property
    def available(self):
        return self.path.is_file()

    def _connect(self):
        if not self.available:
            return None
        uri = self.path.resolve().as_uri() + "?mode=ro"
        connection = sqlite3.connect(uri, uri=True)
        connection.row_factory = sqlite3.Row
        return connection

    def metadata(self):
        """Return catalogue build metadata, or an empty mapping when unavailable."""
        connection = self._connect()
        if connection is None:
            return {}
        try:
            return {
                row["key"]: row["value"]
                for row in connection.execute("SELECT key, value FROM metadata")
            }
        finally:
            connection.close()

    def lookup_formula(self, formula, *, limit=None):
        """Return validated species proposals for one exact molecular formula."""
        canonical = canonical_formula(formula)
        connection = self._connect()
        if connection is None:
            return []
        query = """
            SELECT f.formula, f.exact_mass, s.nist_id, s.name, s.cas,
                   s.inchi, s.inchi_key, s.url
            FROM formula AS f
            JOIN species AS s ON s.formula_id = f.id
            WHERE f.formula = ?
            ORDER BY s.name COLLATE NOCASE, s.nist_id, s.url
        """
        parameters = [canonical]
        if limit is not None:
            query += " LIMIT ?"
            parameters.append(max(0, int(limit)))
        try:
            rows = connection.execute(query, parameters).fetchall()
        finally:
            connection.close()
        return [dict(row) for row in rows]

    def formulas_in_mass_range(self, neutral_mass, tol_da=0.012, *, limit=250):
        """Return catalogue formulae inside a strict neutral exact-mass window."""
        target = float(neutral_mass)
        tolerance = float(tol_da)
        connection = self._connect()
        if connection is None:
            return []
        try:
            rows = connection.execute(
                """
                SELECT formula, exact_mass
                FROM formula
                WHERE exact_mass BETWEEN ? AND ?
                ORDER BY ABS(exact_mass - ?), formula
                LIMIT ?
                """,
                (target - tolerance, target + tolerance, target, int(limit)),
            ).fetchall()
        finally:
            connection.close()
        output = []
        for row in rows:
            formula = canonical_formula(row["formula"])
            exact_mass = formula_id.formula_mass(isotopes.parse_formula(formula))
            if abs(exact_mass - target) <= tolerance:
                output.append({"formula": formula, "exact_mass": exact_mass})
        return output

    def search(self, query, *, limit=50):
        """Search names, formulae, CAS numbers, InChI and InChIKey offline."""
        value = " ".join(str(query or "").split())
        if not value:
            return []
        connection = self._connect()
        if connection is None:
            return []
        pattern = f"%{value}%"
        try:
            rows = connection.execute(
                """
                SELECT f.formula, f.exact_mass, s.nist_id, s.name, s.cas,
                       s.inchi, s.inchi_key, s.url
                FROM species AS s
                JOIN formula AS f ON f.id = s.formula_id
                WHERE s.name LIKE ? COLLATE NOCASE
                   OR f.formula = ? COLLATE NOCASE
                   OR s.cas = ?
                   OR s.inchi LIKE ? COLLATE NOCASE
                   OR s.inchi_key = ? COLLATE NOCASE
                ORDER BY CASE WHEN s.name = ? COLLATE NOCASE THEN 0 ELSE 1 END,
                         s.name COLLATE NOCASE, s.nist_id
                LIMIT ?
                """,
                (pattern, value, value, pattern, value, value, int(limit)),
            ).fetchall()
        finally:
            connection.close()
        return [dict(row) for row in rows]

    def enrich_candidates(self, candidates):
        """Attach catalogue proposals without changing PTR names or candidate scores."""
        enriched = []
        for candidate in candidates or []:
            candidate = dict(candidate)
            formula = candidate.get("formula")
            if formula:
                matches = self.lookup_formula(formula)
                if matches:
                    candidate["catalogue"] = matches
            enriched.append(candidate)
        return enriched

    def score_peak(
        self,
        mz,
        *,
        drift=1.0,
        candidates=None,
        obs_ratios=None,
        compounds_of_interest=None,
        elements=None,
        tol_mDa=12.0,
    ):
        """Supplement and enrich a peak using only the strict local mass window."""
        tolerance = float(tol_mDa) / 1000.0
        neutral_mass = float(mz) / float(drift) - formula_id.PROTON
        extra = [
            item["formula"]
            for item in self.formulas_in_mass_range(neutral_mass, tolerance)
        ]
        # Keep already-computed candidate ordering when the catalogue contributes no
        # new formula. This avoids needless score drift in old saved reviews.
        if not extra and candidates is not None:
            return self.enrich_candidates(candidates)
        scored = formula_id.score_peak(
            float(mz),
            float(drift),
            obs_ratios=obs_ratios,
            compounds_of_interest=compounds_of_interest,
            elements=elements,
            extra_formulas=extra,
            tol_mDa=tol_mDa,
        )
        return self.enrich_candidates(scored)
