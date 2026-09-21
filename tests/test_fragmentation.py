"""Fragmentation-profile provenance and conservative co-variation scoring."""

import h5py
import numpy as np

from sniff import fragmentation, ptrms


def _context(e_n=120.0):
    return {
        "model": fragmentation.MODEL,
        "status": "available",
        "reagent": "H3O+",
        "e_n_td": e_n,
    }


def _rate_table(profile_e_n=120.0):
    return {
        "compounds": [
            {
                "formula": "C4H8O",
                "fragmentation_profiles": [
                    {
                        "profile_id": "test-profile",
                        "name": "test ketone",
                        "reagent": "H3O+",
                        "parent_mz": 73.065,
                        "e_n_td": profile_e_n,
                        "products": [
                            {
                                "mz": 43.018,
                                "pathway_type": "fragment",
                                "yield_percent": 70.0,
                            },
                            {
                                "mz": 91.075,
                                "pathway_type": "adduct",
                                "yield_percent": 10.0,
                            },
                        ],
                        "reference": "synthetic test reference",
                        "doi": "doi:test",
                    }
                ],
            }
        ]
    }


def _varying_trace(scale=1.0):
    base = np.tile(np.r_[np.zeros(20), np.full(20, 100.0)], 5)
    return base * scale


def test_reaction_context_reads_measured_primary_ion_and_en():
    with h5py.File("reaction", "w", driver="core", backing_store=False) as h5:
        info = np.array(
            [
                [b"E/N_Act", b"PrimionIdx"],
                [b"Td", b""],
            ],
            dtype="S32",
        )
        h5.create_dataset("AddTraces/PTR-Reaction/Info", data=info)
        h5.create_dataset(
            "AddTraces/PTR-Reaction/Data",
            data=np.array([[119.0, 0.0], [121.0, 0.0], [120.0, 0.0]]),
        )
        h5.create_dataset(
            "PTR-PrimaryIons/Descriptions",
            data=np.array([b"H3O+", b"NO+"], dtype="S32"),
        )

        context = fragmentation.reaction_context(h5)

    assert context["status"] == "available"
    assert context["reagent"] == "H3O+"
    assert context["primary_ion_index"] == 0
    assert context["e_n_td"] == 120.0
    assert context["e_n_range_td"] == [119.0, 121.0]


def test_reaction_context_withholds_mixed_reagent_or_en_runs():
    with h5py.File("mixed", "w", driver="core", backing_store=False) as h5:
        h5.create_dataset(
            "AddTraces/PTR-Reaction/Info",
            data=np.array([[b"E/N_Act", b"PrimionIdx"], [b"Td", b""]], dtype="S32"),
        )
        h5.create_dataset(
            "AddTraces/PTR-Reaction/Data",
            data=np.array([[100.0, 0.0], [125.0, 1.0]]),
        )
        h5.create_dataset(
            "PTR-PrimaryIons/Descriptions",
            data=np.array([b"H3O+", b"NO+"], dtype="S32"),
        )

        context = fragmentation.reaction_context(h5)

    assert context["status"] == "unavailable"
    assert "primary ion changes" in context["reason"]


def test_fragment_support_reranks_existing_candidate_and_links_channel():
    peaks = [
        {
            "mz": 73.065,
            "candidates": [
                {
                    "formula": "C3H4O2",
                    "assignment_eligible": True,
                    "probability": 0.52,
                    "score": 0.52,
                },
                {
                    "formula": "C4H8O",
                    "assignment_eligible": True,
                    "probability": 0.48,
                    "score": 0.48,
                },
            ],
        },
        {"mz": 43.018, "candidates": []},
    ]
    traces = {
        73.065: _varying_trace(),
        43.018: _varying_trace(scale=0.7),
    }

    fragmentation.apply_fragmentation_evidence(
        peaks, traces, _rate_table(), _context()
    )

    assert peaks[0]["candidates"][0]["formula"] == "C4H8O"
    supported = peaks[0]["candidates"][0]
    assert supported["assignment_eligible"] is True
    assert supported["fragmentation_evidence"]["status"] == "support"
    assert supported["fragmentation_factor"] == fragmentation.SUPPORT_FACTOR
    assert peaks[1]["candidates"] == []
    assert "id_confidence" not in peaks[1]
    assert peaks[1]["fragmentation_links"][0]["parent_formula"] == "C4H8O"
    assert all(
        link["expected_fragment_mz"] != 91.075
        for peak in peaks
        for link in peak.get("fragmentation_links", [])
    )


def test_fragment_support_never_promotes_broad_candidate_to_assignment():
    peaks = [
        {
            "mz": 73.065,
            "candidates": [
                {
                    "formula": "C4H8O",
                    "assignment_eligible": False,
                    "mass_match": "broad-proposal-only",
                    "probability": 1.0,
                    "score": 1.0,
                }
            ],
        },
        {"mz": 43.018, "candidates": []},
    ]

    fragmentation.apply_fragmentation_evidence(
        peaks,
        {73.065: _varying_trace(), 43.018: _varying_trace(scale=0.5)},
        _rate_table(),
        _context(),
    )

    candidate = peaks[0]["candidates"][0]
    assert candidate["fragmentation_evidence"]["status"] == "support"
    assert candidate["assignment_eligible"] is False
    assert candidate["mass_match"] == "broad-proposal-only"


def test_condition_mismatch_and_flat_traces_cannot_support_candidate():
    mismatch_peaks = [
        {
            "mz": 73.065,
            "candidates": [
                {
                    "formula": "C4H8O",
                    "assignment_eligible": True,
                    "probability": 1.0,
                    "score": 1.0,
                }
            ],
        },
        {"mz": 43.018, "candidates": []},
    ]
    fragmentation.apply_fragmentation_evidence(
        mismatch_peaks,
        {73.065: _varying_trace(), 43.018: _varying_trace()},
        _rate_table(profile_e_n=80.0),
        _context(e_n=120.0),
    )
    assert (
        mismatch_peaks[0]["candidates"][0]["fragmentation_evidence"]["status"]
        == "unavailable"
    )

    flat_peaks = [
        {
            "mz": 73.065,
            "candidates": [
                {
                    "formula": "C4H8O",
                    "assignment_eligible": True,
                    "probability": 1.0,
                    "score": 1.0,
                }
            ],
        },
        {"mz": 43.018, "candidates": []},
    ]
    fragmentation.apply_fragmentation_evidence(
        flat_peaks,
        {73.065: np.ones(100), 43.018: np.ones(100)},
        _rate_table(),
        _context(),
    )
    evidence = flat_peaks[0]["candidates"][0]["fragmentation_evidence"]
    assert evidence["status"] == "inconclusive"
    assert flat_peaks[0]["candidates"][0]["fragmentation_factor"] == 1.0


def test_bundled_profiles_preserve_pathway_conditions_and_provenance():
    table = ptrms.load_rate_constants()
    propene = next(item for item in table["compounds"] if item["formula"] == "C3H6")
    profile = next(
        item
        for item in propene["fragmentation_profiles"]
        if item["profile_id"] == "ptrlibrary-row-69"
    )

    assert profile["name"] == "propene"
    assert profile["reagent"] == "H3O+"
    assert profile["e_n_td"] == 106.0
    assert profile["products"] == [
        {"mz": 41.039, "pathway_type": "fragment", "yield_percent": 15.0}
    ]
    assert profile["reference"]
    assert profile["doi"]
