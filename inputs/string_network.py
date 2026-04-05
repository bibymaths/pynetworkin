# inputs/string_network.py
#
# Fetches the STRING functional association network.  Priority order:
#
#   1. Flat file (primary) — data/string_data/9606.links.v12.0.tsv.gz
#      Override path with NETWORKIN_STRING_FLAT_FILE env var.
#   2. STRING REST API (fallback for small queries) — string-db.org
#   3. Bundled sample TSV — data/fallback/string_sample.tsv
#
# API docs: https://string-db.org/cgi/help?sessionId=&subpage=api
#
# Cache: Parquet files in .cache/, refreshed if older than 7 days.

import os
from datetime import date, timedelta
from io import StringIO
from pathlib import Path

import httpx
import pandas as pd

CACHE_DIR = Path(os.environ.get("NETWORKIN_CACHE_DIR", ".cache"))
CACHE_DIR.mkdir(exist_ok=True)

STRING_API_NETWORK = "https://string-db.org/api/tsv/network"
STRING_API_MAPPING = "https://string-db.org/api/tsv/get_string_ids"
CACHE_TTL_DAYS = 7

_REPO_ROOT = Path(__file__).resolve().parent.parent
_DEFAULT_FLAT_FILE = _REPO_ROOT / "data" / "string_data" / "9606.links.v12.0.tsv.gz"
FALLBACK_TSV = _REPO_ROOT / "data" / "fallback" / "string_sample.tsv"


def _cache_path(name: str) -> Path:
    return CACHE_DIR / f"{name}.parquet"


def _cache_valid(path: Path) -> bool:
    if not path.exists():
        return False
    mod_date = date.fromtimestamp(path.stat().st_mtime)
    return (date.today() - mod_date) < timedelta(days=CACHE_TTL_DAYS)


def _load_fallback() -> pd.DataFrame:
    """Return the bundled STRING sample as a last resort."""
    if FALLBACK_TSV.exists():
        return pd.read_csv(FALLBACK_TSV, sep="\t", low_memory=False)
    return pd.DataFrame(columns=["protein_a", "protein_b", "combined_score"])


def _load_flat_file(
    proteins: list[str] | None = None,
    species: int = 9606,
    min_score: int = 400,
) -> pd.DataFrame:
    """
    Load interactions from a local STRING flat-file (gzip TSV).

    The file path defaults to data/string_data/9606.links.v12.0.tsv.gz and
    can be overridden with the NETWORKIN_STRING_FLAT_FILE environment variable.

    The flat file has the format produced by string_network_filter.py:
        tree  kinase_a  kinase_b  string_id_a  string_id_b  score

    Parameters
    ----------
    proteins  : optional list of gene/protein names to filter by
    species   : NCBI taxonomy ID (used only to build the default file path)
    min_score : minimum combined score [0, 1000]
    """
    flat_file = Path(os.environ.get("NETWORKIN_STRING_FLAT_FILE", str(_DEFAULT_FLAT_FILE)))
    if not flat_file.exists():
        raise FileNotFoundError(f"STRING flat file not found: {flat_file}")

    import gzip as _gzip

    open_fn = _gzip.open if flat_file.suffix == ".gz" else open
    with open_fn(flat_file, "rt") as fh:
        df = pd.read_csv(fh, sep="\t", header=None, low_memory=False)

    # The filter script produces: tree kinase_a kinase_b string_id_a string_id_b score
    if df.shape[1] >= 6:
        df.columns = [
            "tree",
            "protein_a",
            "protein_b",
            "string_id_a",
            "string_id_b",
            "score",
        ] + list(df.columns[6:])
    elif df.shape[1] >= 3:
        # Minimal variant: protein_a, protein_b, score
        df.columns = ["protein_a", "protein_b", "score"] + list(df.columns[3:])
    else:
        raise ValueError(f"Unexpected column count in STRING flat file: {df.shape[1]}")

    df["combined_score"] = pd.to_numeric(df["score"], errors="coerce") / 1000.0
    df = df[df["combined_score"] >= min_score / 1000.0]

    if proteins:
        protein_set = set(proteins)
        mask = df["protein_a"].isin(protein_set) | df["protein_b"].isin(protein_set)
        df = df[mask]

    return df[["protein_a", "protein_b", "combined_score"]].dropna()


def _fetch_rest_api(
    proteins: list[str],
    species: int = 9606,
    min_score: int = 400,
) -> pd.DataFrame:
    """
    Fetch a network subset from the STRING REST API.

    Suitable for small queries (< 2 000 proteins).
    """
    response = httpx.post(
        STRING_API_NETWORK,
        data={
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
    df["combined_score"] = df["score"] / 1000.0
    df = df.rename(columns={"preferredName_A": "protein_a", "preferredName_B": "protein_b"})
    return df[["protein_a", "protein_b", "combined_score"]]


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
    proteins    : list of gene/protein names (used to filter flat file or query REST)
    species     : NCBI taxonomy ID (default 9606 = Homo sapiens)
    min_score   : combined score threshold [0, 1000] (default 400 = medium confidence)
    refresh     : ignore local cache and re-fetch

    Returns
    -------
    DataFrame with columns: protein_a, protein_b, combined_score (float, [0,1])
    """
    cache_key = f"string_{species}_{min_score}"
    cache = _cache_path(cache_key)
    if not refresh and _cache_valid(cache):
        return pd.read_parquet(cache)

    # ── 1. Flat file ──────────────────────────────────────────────────────────
    try:
        df = _load_flat_file(proteins=proteins, species=species, min_score=min_score)
        if not df.empty:
            df.to_parquet(cache, index=False)
            return df
    except Exception:
        pass  # fall through to REST API

    # ── 2. REST API ───────────────────────────────────────────────────────────
    try:
        df = _fetch_rest_api(proteins=proteins, species=species, min_score=min_score)
        if not df.empty:
            df.to_parquet(cache, index=False)
            return df
    except Exception:
        pass  # fall through to bundled fallback

    # ── 3. Bundled fallback ───────────────────────────────────────────────────
    return _load_fallback()
