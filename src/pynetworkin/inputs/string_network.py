# inputs/string_network.py
#
# Fetches the STRING functional association network.  Priority order:
#
#   1. NETWORKIN_STRING_FLAT_FILE env var — explicit user-supplied path.
#   2. Cached downloaded flat file — ~/.cache/string_data/ (or NETWORKIN_CACHE_DIR).
#   3. Runtime download from stringdb-downloads.org.
#   4. STRING REST API (fallback for small queries) — string-db.org.
#   5. Bundled sample TSV — pynetworkin/data/fallback/string_sample.tsv.
#
# API docs: https://string-db.org/cgi/help?sessionId=&subpage=api
#
# Cache: Parquet files in .cache/, refreshed if older than 7 days.

import os
import tempfile
from datetime import date, timedelta
from io import StringIO
from pathlib import Path

import httpx
import pandas as pd

from pynetworkin.logger import logger
from pynetworkin.resources import open_fallback_string

CACHE_DIR = Path(os.environ.get("NETWORKIN_CACHE_DIR", ".cache"))
CACHE_DIR.mkdir(exist_ok=True)

STRING_API_NETWORK = "https://string-db.org/api/tsv/network"
STRING_API_MAPPING = "https://string-db.org/api/tsv/get_string_ids"
CACHE_TTL_DAYS = 7

# Full human STRING v12.0 flat file — downloaded at runtime when needed.
STRING_FULL_DOWNLOAD_URL = (
    "https://stringdb-downloads.org/download/"
    "protein.links.v12.0/9606.protein.links.v12.0.txt.gz"
)
_STRING_CACHE_FILENAME = "9606.protein.links.v12.0.txt.gz"
_STRING_DOWNLOAD_TIMEOUT = 300  # seconds


def _string_download_cache_dir() -> Path:
    """Return the directory used to cache the downloaded STRING flat file."""
    return CACHE_DIR / "string_data"


def _downloaded_flat_file_path() -> Path:
    return _string_download_cache_dir() / _STRING_CACHE_FILENAME


def _cache_path(name: str) -> Path:
    return CACHE_DIR / f"{name}.parquet"


def _cache_valid(path: Path) -> bool:
    if not path.exists():
        return False
    mod_date = date.fromtimestamp(path.stat().st_mtime)
    return (date.today() - mod_date) < timedelta(days=CACHE_TTL_DAYS)


def _load_fallback() -> pd.DataFrame:
    """Return the bundled STRING sample as a last resort."""
    logger.warning("Falling back to bundled STRING sample data")
    try:
        with open_fallback_string() as fh:
            return pd.read_csv(fh, sep="\t", low_memory=False)
    except Exception:
        return pd.DataFrame(columns=["protein_a", "protein_b", "combined_score"])


def _download_string_flat_file() -> Path:
    """Download the full human STRING v12.0 flat file to the cache directory.

    Uses streaming download + atomic temp-file rename to avoid partial files.

    Returns
    -------
    Path
        Path to the downloaded (and now cached) file.

    Raises
    ------
    RuntimeError
        If the download fails.
    """
    dest = _downloaded_flat_file_path()
    dest.parent.mkdir(parents=True, exist_ok=True)

    logger.info("Downloading STRING flat file from {}", STRING_FULL_DOWNLOAD_URL)
    try:
        with httpx.stream(
            "GET",
            STRING_FULL_DOWNLOAD_URL,
            timeout=_STRING_DOWNLOAD_TIMEOUT,
            follow_redirects=True,
        ) as response:
            response.raise_for_status()
            # Write to a temp file in the same directory, then rename atomically.
            tmp_fd, tmp_path = tempfile.mkstemp(dir=dest.parent, suffix=".tmp")
            try:
                with open(tmp_fd, "wb") as tmp_fh:
                    for chunk in response.iter_bytes(chunk_size=1024 * 1024):
                        tmp_fh.write(chunk)
                Path(tmp_path).replace(dest)
            except Exception:
                try:
                    Path(tmp_path).unlink(missing_ok=True)
                except OSError:
                    pass
                raise
    except httpx.HTTPError as exc:
        raise RuntimeError(
            f"Failed to download STRING flat file from {STRING_FULL_DOWNLOAD_URL}: {exc}"
        ) from exc

    logger.success("STRING flat file cached at {}", dest)
    return dest


def _resolve_flat_file() -> Path:
    """Return a Path to the full STRING flat file, downloading if necessary.

    Priority:
      1. NETWORKIN_STRING_FLAT_FILE env var (if set and file exists)
      2. Cached downloaded file (if present and non-empty)
      3. Download from stringdb-downloads.org

    Raises
    ------
    RuntimeError
        If download fails and no local copy is available.
    FileNotFoundError
        If the env-var path is set but the file does not exist.
    """
    # 1. Explicit env-var override
    env_path = os.environ.get("NETWORKIN_STRING_FLAT_FILE", "")
    if env_path:
        p = Path(env_path)
        if not p.exists():
            raise FileNotFoundError(
                f"NETWORKIN_STRING_FLAT_FILE is set but the file does not exist: {p}"
            )
        logger.info("Using STRING flat file from env var: {}", p)
        return p

    # 2. Cached downloaded file
    cached = _downloaded_flat_file_path()
    if cached.exists() and cached.stat().st_size > 0:
        logger.info("Using cached STRING flat file: {}", cached)
        return cached

    # 3. Download
    return _download_string_flat_file()


def _load_flat_file(
    proteins: list[str] | None = None,
    species: int = 9606,
    min_score: int = 400,
) -> pd.DataFrame:
    """Load interactions from the full STRING flat-file (gzip TSV).

    Resolves the file path using :func:`_resolve_flat_file` (env-var →
    cache → download).

    The raw STRING protein.links file has columns:
        protein1  protein2  combined_score

    Parameters
    ----------
    proteins  : optional list of gene/protein names to filter by
    species   : NCBI taxonomy ID (for context only; resolution is handled above)
    min_score : minimum combined score [0, 1000]
    """
    import gzip as _gzip

    flat_file = _resolve_flat_file()

    open_fn = _gzip.open if flat_file.suffix == ".gz" else open
    with open_fn(flat_file, "rt") as fh:
        df = pd.read_csv(fh, sep=" ", low_memory=False)

    # The raw protein.links file has a header: protein1 protein2 combined_score
    if "protein1" in df.columns and "protein2" in df.columns:
        df = df.rename(columns={"protein1": "protein_a", "protein2": "protein_b"})
        df["combined_score"] = pd.to_numeric(df["combined_score"], errors="coerce") / 1000.0
    elif df.shape[1] >= 6:
        # Legacy filtered format: tree kinase_a kinase_b string_id_a string_id_b score
        df.columns = [
            "tree",
            "protein_a",
            "protein_b",
            "string_id_a",
            "string_id_b",
            "score",
        ] + list(df.columns[6:])
        df["combined_score"] = pd.to_numeric(df["score"], errors="coerce") / 1000.0
    elif df.shape[1] >= 3:
        df.columns = ["protein_a", "protein_b", "score"] + list(df.columns[3:])
        df["combined_score"] = pd.to_numeric(df["score"], errors="coerce") / 1000.0
    else:
        raise ValueError(f"Unexpected column count in STRING flat file: {df.shape[1]}")

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
    """Fetch a network subset from the STRING REST API.

    Suitable for small queries (< 2 000 proteins).
    """
    logger.info("Falling back to STRING REST API")
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
    """Fetch the STRING functional association network for a list of protein IDs.

    Parameters
    ----------
    proteins    : list of gene/protein names (used to filter flat file or query REST)
    species     : NCBI taxonomy ID (default 9606 = Homo sapiens)
    min_score   : combined score threshold [0, 1000] (default 400 = medium confidence)
    refresh     : ignore local parquet cache and re-fetch

    Returns
    -------
    DataFrame with columns: protein_a, protein_b, combined_score (float, [0,1])
    """
    cache_key = f"string_{species}_{min_score}"
    cache = _cache_path(cache_key)
    if not refresh and _cache_valid(cache):
        return pd.read_parquet(cache)

    # ── 1 & 2 & 3. Flat file (env-var / cached / download) ───────────────────
    try:
        df = _load_flat_file(proteins=proteins, species=species, min_score=min_score)
        if not df.empty:
            df.to_parquet(cache, index=False)
            return df
    except Exception:
        pass  # fall through to REST API

    # ── 4. REST API ───────────────────────────────────────────────────────────
    try:
        df = _fetch_rest_api(proteins=proteins, species=species, min_score=min_score)
        if not df.empty:
            df.to_parquet(cache, index=False)
            return df
    except Exception:
        pass  # fall through to bundled fallback

    # ── 5. Bundled fallback ───────────────────────────────────────────────────
    return _load_fallback()
