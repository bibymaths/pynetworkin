# recovery.py
#
# False negative recovery for NetworKIN.
#
# The original NetworKIN used the STRING network context score only as a
# false-positive filter on top of motif predictions. This means any
# kinase–substrate pair that the motif step missed entirely was permanently
# lost, regardless of strong network context evidence.
#
# This module recovers such pairs using network proximity (context score)
# as the primary signal, inspired by the LinkPhinder knowledge-graph approach.
#
# A pair is recovered if:
#   1. It was NOT scored by the motif step (absent from motif_scores dict).
#   2. Its context_score (1 / (1 + shortest_path_dist)) >= CONTEXT_RECOVERY_THRESHOLD.
#
# Recovered predictions are tagged with recovery_method = "context_proximity"
# and recovered = True in the output.

import numpy as np

CONTEXT_RECOVERY_THRESHOLD = 0.6  # tune empirically; 0.6 corresponds to path length <= 0.67


def recover_false_negatives(
    candidates: list[tuple[str, str]],
    dist_matrix: np.ndarray,
    node_index: dict[str, int],
    motif_scores: dict[tuple[str, str], float],
) -> list[dict]:
    """
    Recover kinase–substrate pairs missed by the motif step.

    Parameters
    ----------
    candidates   : all (kinase_id, substrate_uniprot) pairs from the STRING network
    dist_matrix  : all-pairs shortest path matrix from floyd_warshall()
    node_index   : maps protein ID -> row/col index in dist_matrix
    motif_scores : dict of {(kinase_id, substrate_id): motif_score} from pynetphorest

    Returns
    -------
    List of prediction dicts for recovered pairs, each with:
        kinase_id, substrate_uniprot, motif_score (-1.0 sentinel),
        context_score, networkin_score, recovered (True), recovery_method
    """
    recovered = []
    for kinase_id, substrate_id in candidates:
        # Skip if motif step already scored this pair
        if (kinase_id, substrate_id) in motif_scores:
            continue

        ki = node_index.get(kinase_id)
        si = node_index.get(substrate_id)
        if ki is None or si is None:
            continue

        d = dist_matrix[ki, si]
        if d == np.inf:
            continue

        c_score = float(1.0 / (1.0 + d))
        if c_score >= CONTEXT_RECOVERY_THRESHOLD:
            recovered.append(
                {
                    "kinase_id": kinase_id,
                    "substrate_uniprot": substrate_id,
                    "motif_score": -1.0,  # sentinel: motif step did not score this pair
                    "context_score": c_score,
                    "networkin_score": c_score,  # context-only score for recovered pairs
                    "recovered": True,
                    "recovery_method": "context_proximity",
                }
            )

    return recovered
