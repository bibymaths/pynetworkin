# graph_scoring.py
#
# STRING network–based graph / context scoring for NetworKIN.
#
# The original NetworKIN pipeline combines two independent signals:
#   1. A motif-based likelihood (from pynetphorest, handled in motif_scoring.py).
#   2. A network context likelihood derived from the STRING protein–protein
#      interaction network (handled here).
#
# The final NetworKIN score is the product of these two likelihoods.
# A high score therefore requires both a good sequence motif match AND
# a short network distance between kinase and substrate in STRING.

from __future__ import annotations

from typing import Any

import pandas as pd

from logger import logger


def compute_networkin_score(
    motif_likelihood: float,
    context_likelihood: float,
) -> float:
    """
    Combine motif and STRING context likelihoods into the final NetworKIN score.

    The NetworKIN score is defined as the product of the motif likelihood
    (from pynetphorest sequence scoring) and the STRING network context
    likelihood (derived from graph-based shortest-path distance).

    Parameters
    ----------
    motif_likelihood : float
        Likelihood ratio from the motif scoring step (>= 0).
    context_likelihood : float
        Likelihood ratio from the STRING network context scoring step (>= 0).

    Returns
    -------
    float
        NetworKIN score = motif_likelihood × context_likelihood.

    Examples
    --------
    >>> compute_networkin_score(2.0, 1.5)
    3.0
    >>> compute_networkin_score(0.0, 5.0)
    0.0
    """
    return motif_likelihood * context_likelihood


def filter_and_rank_predictions(
    predictions: list[dict[str, Any]],
    min_networkin: float = 2.0,
    min_motif: float = 0.05,
    top_k: int = 5,
) -> list[dict[str, Any]]:
    """
    Filter and rank kinase–substrate predictions by NetworKIN score.

    Removes predictions below minimum score thresholds, then keeps only
    the top-*k* kinase candidates per (target protein, phosphosite position)
    pair, sorted by descending NetworKIN score.

    Parameters
    ----------
    predictions : list of dict
        Raw prediction rows as returned by ``compile_predictions`` or
        ``recover_predictions``.  Each dict must contain the keys
        ``"Name"``, ``"Position"``, ``"NetworKIN score"``, and
        ``"Motif probability"``.
    min_networkin : float, optional
        Minimum NetworKIN score to retain a prediction.  Default is ``2.0``.
    min_motif : float, optional
        Minimum motif probability to retain a prediction.  Default is ``0.05``.
    top_k : int, optional
        Maximum number of kinase predictions to keep per (protein, position)
        pair.  Default is ``5``.

    Returns
    -------
    list of dict
        Filtered and ranked prediction rows in the same format as the input.
        The list is sorted by (``Name``, ``Position``, ``NetworKIN score``
        descending).

    Notes
    -----
    Predictions recovered by the false-negative recovery step
    (``recovered = True``) may have ``"Motif probability"`` set to ``-1.0``
    as a sentinel value.  Such rows will be excluded by the ``min_motif``
    filter unless you set ``min_motif`` to a negative value.
    """
    df = pd.DataFrame(predictions)
    if df.empty:
        logger.warning("No predictions available after scoring")
        return predictions

    filtered = df[
        (df["NetworKIN score"] > min_networkin) & (df["Motif probability"] > min_motif)
    ].copy()
    filtered = filtered.sort_values(
        ["Name", "Position", "NetworKIN score"], ascending=[True, True, False]
    )
    filtered = filtered.groupby(["Name", "Position"], as_index=False).head(top_k)

    logger.success("Retained {} predictions after ranking and filtering", len(filtered))
    return filtered.to_dict(orient="records")
