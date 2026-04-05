"""motif_scoring.py – batch kinase-motif scorer.

Replaces the legacy subprocess-based C-binary calls for motif scoring.
All sequences are scored in a single batched call; no per-sequence
subprocess or loop is used.
"""

from __future__ import annotations

from typing import Dict, Optional, Tuple

from pynetphorest import core as _core

# Load the atlas once at module import time; subsequent calls reuse this.
_MODELS = _core.load_atlas(None)


def score_sequences(
    id_seq: Dict[str, str],
    id_pos_res: Optional[Dict[str, Dict[int, str]]] = None,
) -> Dict[str, Dict[int, Dict[str, Dict[str, Tuple[str, str, float]]]]]:
    """Score all sequences in *id_seq* against the bundled kinase motif atlas.

    Parameters
    ----------
    id_seq:
        Mapping of protein-id → full protein sequence.
    id_pos_res:
        Optional mapping of protein-id → {1-based position: residue}.
        When provided, only those positions are scored.
        When ``None`` or empty dict, every S/T/Y in each sequence is scored.

    Returns
    -------
    Nested dict: ``result[protein_id][pos][tree][kinase] = (res, peptide, score)``

    where

    - *pos*    is the 1-based sequence position
    - *tree*   is the top-level classifier type (e.g. ``'KIN'``, ``'SH2'``, ``'1433'``)
    - *kinase* is the specific kinase / group name (e.g. ``'PKA_group'``)
    - *res*    is the phospho-residue character (``'S'``, ``'T'``, or ``'Y'``)
    - *peptide* is the display window from the motif scorer
    - *score*  is the posterior probability in ``[0, 1]``
    """
    if id_pos_res is None:
        id_pos_res = {}

    result: Dict[str, Dict[int, Dict[str, Dict[str, Tuple[str, str, float]]]]] = {}

    for protein_id, seq in id_seq.items():
        if id_pos_res and protein_id not in id_pos_res:
            continue

        seq_upper = seq.upper()
        allowed_positions = id_pos_res.get(protein_id, {})

        for i, aa in enumerate(seq_upper):
            if aa not in ("S", "T", "Y"):
                continue
            pos1 = i + 1
            if allowed_positions and pos1 not in allowed_positions:
                continue

            peptide = _core.get_display_window(seq_upper, i)

            for model in _MODELS:
                if aa not in model["residues"]:
                    continue
                score = _core.get_model_posterior(seq_upper, i, model)
                if score <= 0.0:
                    continue

                meta = model["meta"]
                tree = meta["classifier"]   # e.g. 'KIN', 'SH2', '1433'
                kinase = meta["kinase"]     # e.g. 'PKA_group', 'Abl_group'

                if protein_id not in result:
                    result[protein_id] = {}
                if pos1 not in result[protein_id]:
                    result[protein_id][pos1] = {}
                if tree not in result[protein_id][pos1]:
                    result[protein_id][pos1][tree] = {}
                result[protein_id][pos1][tree][kinase] = (aa, peptide, score)

    return result
