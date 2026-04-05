# inputs/string_network.py
#
# Fetches the STRING functional association network via the STRING REST API.
# API docs: https://string-db.org/cgi/help?sessionId=&subpage=api
#
# Cache: Parquet files in .cache/, refreshed if older than 7 days.

from io import StringIO
from pathlib import Path
from datetime import date, timedelta

import httpx
import pandas as pd

CACHE_DIR = Path(".cache")
CACHE_DIR.mkdir(exist_ok=True)

STRING_API_NETWORK  = "https://string-db.org/api/tsv/network"
STRING_API_MAPPING  = "https://string-db.org/api/tsv/get_string_ids"
CACHE_TTL_DAYS = 7


def _cache_path(name: str) -> Path:
    return CACHE_DIR / f"{name}.parquet"


def _cache_valid(path: Path) -> bool:
    if not path.exists():
        return False
    mod_date = date.fromtimestamp(path.stat().st_mtime)
    return (date.today() - mod_date) < timedelta(days=CACHE_TTL_DAYS)


def fetch_string_network(
    proteins: list[str],
    species: int = 9606,
    min_score: int = 400,
    refresh: bool = False,
) -> pd.DataFrame:
    """
    Fetch the STRING functional association network for a list of UniProt IDs.

    Parameters
    ----------
    proteins    : list of UniProt IDs
    species     : NCBI taxonomy ID (default 9606 = Homo sapiens)
    min_score   : combined score threshold [0, 1000] (default 400 = medium confidence)
    refresh     : ignore local cache and re-fetch

    Returns
    -------
    DataFrame with columns: protein_a, protein_b, combined_score (float, normalised to [0,1])
    """
    cache_key = f"string_{species}_{min_score}"
    cache = _cache_path(cache_key)
    if not refresh and _cache_valid(cache):
        return pd.read_parquet(cache)

    response = httpx.post(
        STRING_API_NETWORK,
        data={
            # STRING API requires identifiers separated by carriage returns (\r)
            # as documented at https://string-db.org/cgi/help?subpage=api%23getting-the-network
            "identifiers": "\r".join(proteins),
            "species": species,
            "required_score": min_score,
            "caller_identity": "networkin_refactor",
        },
        timeout=180,
    )
    response.raise_for_status()

    df = pd.read_csv(StringIO(response.text), sep="\t")
    required_cols = {"score", "preferredName_A", "preferredName_B"}
    missing = required_cols - set(df.columns)
    if missing:
        raise KeyError(
            f"STRING API response is missing expected columns: {missing}. "
            "The API format may have changed."
        )
    # Normalise score to [0, 1]
    df["combined_score"] = df["score"] / 1000.0
    df = df.rename(columns={"preferredName_A": "protein_a", "preferredName_B": "protein_b"})
    df = df[["protein_a", "protein_b", "combined_score"]]

    df.to_parquet(cache, index=False)
    return df
