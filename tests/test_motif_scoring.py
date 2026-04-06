"""Tests for motif_scoring.py batch scorer."""

from pynetworkin.motif_scoring import score_sequences


def test_pka_substrate_lrraslg():
    """PKA substrate KEMPTIDE (LRRASLG): PKA_group must be the top KIN scorer.

    LRRASLG is the canonical basophilic PKA substrate (kemptide).  The S at
    position 5 of the 7-mer is surrounded by the R-R-x-S consensus; we embed
    it in a neutral A-rich context so the NN windows (up to 13 aa wide) have
    valid residues on both sides.
    """
    seq = "AAAALRRASLGAAAA"
    protein_id = "KEMPTIDE"
    results = score_sequences({protein_id: seq})

    s_pos = seq.upper().index("S") + 1  # 1-based

    assert protein_id in results, "protein not scored at all"
    assert s_pos in results[protein_id], f"position {s_pos} not scored"

    kin_scores = results[protein_id][s_pos].get("KIN", {})
    assert "PKA_group" in kin_scores, "PKA_group kinase model not found in KIN tree"

    _, _, pka_score = kin_scores["PKA_group"]
    assert (
        pka_score > 0.2
    ), f"PKA_group score {pka_score:.4f} is not > 0.2 for canonical PKA substrate LRRASLG"
