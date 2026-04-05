# inputs/phosphosites.py
#
# Fetches phosphorylation site data from:
#   - PhosphoSitePlus: https://www.phosphosite.org/downloads/
#   - Phospho.ELM:     http://phospho.elm.eu.org/dumps/
#
# Cache: Parquet files in .cache/ at the repo root, refreshed if older than 7 days.
# Trigger refresh with: refresh=True argument or --refresh CLI flag.

import gzip
import io
import tarfile
from datetime import date, timedelta
from pathlib import Path

import httpx
import pandas as pd

CACHE_DIR = Path(".cache")
CACHE_DIR.mkdir(exist_ok=True)

PHOSPHOSITE_URL = "https://www.phosphosite.org/downloads/Phosphorylation_site_dataset.gz"
# Phospho.ELM only exposes its dump archive over plain HTTP; HTTPS is not available
# for the dump endpoint (as of 2024). The download is verified by raise_for_status().
PHOSPHOELM_URL  = "http://phospho.elm.eu.org/dumps/phosphoELM_all_latest.dump.tgz"

CACHE_TTL_DAYS = 7


def _cache_path(name: str) -> Path:
    return CACHE_DIR / f"{name}.parquet"


def _cache_valid(path: Path) -> bool:
    if not path.exists():
        return False
    mod_date = date.fromtimestamp(path.stat().st_mtime)
    return (date.today() - mod_date) < timedelta(days=CACHE_TTL_DAYS)


def fetch_phosphosite(refresh: bool = False) -> pd.DataFrame:
    """
    Download and parse the PhosphoSitePlus phosphorylation site dataset.

    Returns a DataFrame with columns:
        uniprot_id, gene, site_residue, position, sequence_window, organism
    """
    cache = _cache_path("phosphosite")
    if not refresh and _cache_valid(cache):
        return pd.read_parquet(cache)

    response = httpx.get(PHOSPHOSITE_URL, follow_redirects=True, timeout=120)
    response.raise_for_status()
    with gzip.open(io.BytesIO(response.content), "rt") as f:
        # PhosphoSitePlus TSV has a multi-line header — skip first 3 rows
        df = pd.read_csv(f, sep="\t", skiprows=3, low_memory=False)

    # Normalise column names to snake_case
    df.columns = [c.lower().replace(" ", "_").replace("-", "_") for c in df.columns]

    # Keep only human records by default
    if "organism" in df.columns:
        df = df[df["organism"].str.contains("human", case=False, na=False)]

    df.to_parquet(cache, index=False)
    return df


def fetch_phospho_elm(refresh: bool = False) -> pd.DataFrame:
    """
    Download and parse the Phospho.ELM dump.

    Returns a DataFrame with columns:
        uniprot_id, position, residue, kinase, pmid, species
    """
    cache = _cache_path("phosphoelm")
    if not refresh and _cache_valid(cache):
        return pd.read_parquet(cache)

    response = httpx.get(PHOSPHOELM_URL, follow_redirects=True, timeout=120)
    response.raise_for_status()
    # Phospho.ELM uses tab-separated format inside a .tgz — extract in memory
    df = None
    with tarfile.open(fileobj=io.BytesIO(response.content), mode="r:gz") as tar:
        for member in tar.getmembers():
            if member.name.endswith(".dump") or member.name.endswith(".tab"):
                f = tar.extractfile(member)
                df = pd.read_csv(f, sep="\t", low_memory=False)
                break

    if df is None:
        raise RuntimeError(
            "No .dump or .tab file found in the Phospho.ELM archive. "
            "The archive format may have changed."
        )

    df.columns = [c.lower().replace(" ", "_") for c in df.columns]
    df.to_parquet(cache, index=False)
    return df
