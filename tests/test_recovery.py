import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from recovery import recover_false_negatives, CONTEXT_RECOVERY_THRESHOLD


def test_close_pair_is_recovered():
    """
    A pair with no motif score but a short graph distance must be recovered.
    """
    node_index = {"KIN_A": 0, "SUB_B": 1, "SUB_C": 2}
    dist_matrix = np.array([
        [0.0, 0.3, 5.0],   # KIN_A is close to SUB_B (d=0.3), far from SUB_C
        [0.3, 0.0, 5.0],
        [5.0, 5.0, 0.0],
    ], dtype=np.float32)

    candidates = [("KIN_A", "SUB_B"), ("KIN_A", "SUB_C")]
    motif_scores = {}  # motif step scored nothing

    result = recover_false_negatives(candidates, dist_matrix, node_index, motif_scores)
    recovered_pairs = [(r["kinase_id"], r["substrate_uniprot"]) for r in result]

    assert ("KIN_A", "SUB_B") in recovered_pairs, "Close pair should be recovered"
    assert ("KIN_A", "SUB_C") not in recovered_pairs, "Distant pair should not be recovered"


def test_motif_scored_pair_not_duplicated():
    """
    A pair already in motif_scores must not appear in the recovered list.
    """
    node_index = {"KIN_A": 0, "SUB_B": 1}
    dist_matrix = np.array([[0.0, 0.2], [0.2, 0.0]], dtype=np.float32)
    candidates = [("KIN_A", "SUB_B")]
    motif_scores = {("KIN_A", "SUB_B"): 0.9}  # already scored

    result = recover_false_negatives(candidates, dist_matrix, node_index, motif_scores)
    assert len(result) == 0, "Motif-scored pair must not be re-added by recovery"


def test_infinite_distance_not_recovered():
    """
    A pair with no path in the STRING network (infinite distance) must not be recovered.
    """
    node_index = {"KIN_A": 0, "SUB_B": 1}
    dist_matrix = np.array([[0.0, np.inf], [np.inf, 0.0]], dtype=np.float32)
    candidates = [("KIN_A", "SUB_B")]
    motif_scores = {}

    result = recover_false_negatives(candidates, dist_matrix, node_index, motif_scores)
    assert len(result) == 0, "Disconnected pair must not be recovered"


def test_recovered_fields_complete():
    """
    Every recovered prediction must contain all required output fields.
    """
    node_index = {"KIN_A": 0, "SUB_B": 1}
    dist_matrix = np.array([[0.0, 0.2], [0.2, 0.0]], dtype=np.float32)
    candidates = [("KIN_A", "SUB_B")]
    motif_scores = {}

    result = recover_false_negatives(candidates, dist_matrix, node_index, motif_scores)
    assert len(result) == 1
    r = result[0]
    for field in ["kinase_id", "substrate_uniprot", "motif_score",
                  "context_score", "networkin_score", "recovered", "recovery_method"]:
        assert field in r, f"Missing field: {field}"
    assert r["recovered"] is True
    assert r["motif_score"] == -1.0   # sentinel
