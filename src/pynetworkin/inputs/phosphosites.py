# inputs/phosphosites.py
#
# Fetches phosphorylation site data.  Priority order:
#
#   1. OmniPath (primary) — public REST API, no auth required.
#   2. Local PhosphoSitePlus export — download manually from
#      https://www.phosphosite.org/downloads/ and set
#      NETWORKIN_PSP_LOCAL_FILE=/path/to/Phosphorylation_site_dataset.tsv
#   3. Bundled fallback TSV — data/fallback/phosphosites_sample.tsv
#
# NOTE: Phospho.ELM has been removed — the server is no longer maintained.
# NOTE: Fetching PhosphoSitePlus directly requires login; use a local copy instead.
#
# Cache: Parquet files in .cache/ at the repo root, refreshed if older than 7 days.
# Trigger refresh with: refresh=True argument or --refresh CLI flag.

import io
import os
from datetime import date, timedelta
from pathlib import Path

import httpx
import pandas as pd

CACHE_DIR = Path(os.environ.get("NETWORKIN_CACHE_DIR", ".cache"))
CACHE_DIR.mkdir(exist_ok=True)

# OmniPath enzyme–substrate (kinase–substrate) endpoint — no auth needed.
OMNIPATH_ENZSUB_URL = "https://omnipathdb.org/enzsub"

# Bundled offline fallback (always present in the repository).
_REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
FALLBACK_TSV = _REPO_ROOT / "data" / "fallback" / "phosphosites_sample.tsv"

CACHE_TTL_DAYS = 7

# Expected output columns (normalised schema).
OUTPUT_COLUMNS = ["uniprot_id", "gene", "site_residue", "position", "sequence_window", "organism"]


def _cache_path(name: str) -> Path:
    return CACHE_DIR / f"{name}.parquet"


def _cache_valid(path: Path) -> bool:
    if not path.exists():
        return False
    mod_date = date.fromtimestamp(path.stat().st_mtime)
    return (date.today() - mod_date) < timedelta(days=CACHE_TTL_DAYS)


def _load_fallback() -> pd.DataFrame:
    """Return the bundled phosphosite sample as a last resort."""
    if FALLBACK_TSV.exists():
        return pd.read_csv(FALLBACK_TSV, sep="\t", low_memory=False)
    # Absolute minimum: empty frame with correct columns.
    return pd.DataFrame(columns=OUTPUT_COLUMNS)


def _fetch_omnipath() -> pd.DataFrame:
    """
    Fetch enzyme–substrate interactions from OmniPath.

    Uses the /enzsub endpoint which is publicly accessible.  Falls back first
    to KinaseExtra if the main result is empty, then raises on failure.
    """
    resp = httpx.get(
        OMNIPATH_ENZSUB_URL,
        params={
            "organisms": "9606",
            "fields": "sources,references",
            "format": "tsv",
        },
        timeout=60,
        follow_redirects=True,
    )
    resp.raise_for_status()
    df = pd.read_csv(io.StringIO(resp.text), sep="\t", low_memory=False)

    if df.empty:
        # Retry with the KinaseExtra endpoint (subset of enzsub)
        resp2 = httpx.get(
            "https://omnipathdb.org/interactions",
            params={
                "datasets": "kinaseextra",
                "organisms": "9606",
                "format": "tsv",
            },
            timeout=60,
            follow_redirects=True,
        )
        resp2.raise_for_status()
        df = pd.read_csv(io.StringIO(resp2.text), sep="\t", low_memory=False)

    # Normalise column names
    df.columns = [c.lower().replace(" ", "_").replace("-", "_") for c in df.columns]

    # Map OmniPath columns to our schema
    rename = {
        "substrate": "uniprot_id",
        "substrate_genesymbol": "gene",
        "residue_type": "site_residue",
        "residue_offset": "position",
    }
    df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})

    df["organism"] = "human"
    if "sequence_window" not in df.columns:
        df["sequence_window"] = ""

    # Keep only the canonical columns that downstream code expects
    for col in OUTPUT_COLUMNS:
        if col not in df.columns:
            df[col] = ""
    return df[OUTPUT_COLUMNS].dropna(subset=["uniprot_id"])


def _load_local_phosphositeplus() -> pd.DataFrame:
    """
    Load a locally stored PhosphoSitePlus export.

    The path is read from the NETWORKIN_PSP_LOCAL_FILE environment variable.
    PhosphoSitePlus requires a registered account; this function never
    attempts a direct download.
    """
    local_path = os.environ.get("NETWORKIN_PSP_LOCAL_FILE", "")
    if not local_path:
        raise FileNotFoundError(
            "NETWORKIN_PSP_LOCAL_FILE is not set.  "
            "Download Phosphorylation_site_dataset.gz from "
            "https://www.phosphosite.org/downloads/ and point the env var at it."
        )
    p = Path(local_path)
    if not p.exists():
        raise FileNotFoundError(f"PhosphoSitePlus local file not found: {p}")

    import gzip

    open_fn = gzip.open if p.suffix == ".gz" else open
    with open_fn(p, "rt") as fh:
        # PhosphoSitePlus TSV has a 3-row preamble
        df = pd.read_csv(fh, sep="\t", skiprows=3, low_memory=False)

    df.columns = [c.lower().replace(" ", "_").replace("-", "_") for c in df.columns]
    if "organism" in df.columns:
        df = df[df["organism"].str.contains("human", case=False, na=False)]

    # Map PSP columns to our schema
    rename = {
        "acc_id": "uniprot_id",
        "gene": "gene",
        "mod_rsd": "site_residue",
        "domain": "position",
        "site_+/-7_aa": "sequence_window",
    }
    df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})
    df["organism"] = "human"
    for col in OUTPUT_COLUMNS:
        if col not in df.columns:
            df[col] = ""
    return df[OUTPUT_COLUMNS].dropna(subset=["uniprot_id"])


def fetch_phosphosite(refresh: bool = False) -> pd.DataFrame:
    """
    Return phosphorylation site data using the best available source.

    Priority:
      1. OmniPath REST (public, no auth)
      2. Local PhosphoSitePlus file (NETWORKIN_PSP_LOCAL_FILE env var)
      3. Bundled fallback TSV (data/fallback/phosphosites_sample.tsv)

    Parameters
    ----------
    refresh : bool
        When True, bypass any local cache and re-fetch from the live source.

    Returns
    -------
    DataFrame with columns: uniprot_id, gene, site_residue, position,
                            sequence_window, organism
    """
    cache = _cache_path("phosphosite")
    if not refresh and _cache_valid(cache):
        return pd.read_parquet(cache)

    # ── 1. OmniPath ──────────────────────────────────────────────────────────
    try:
        df = _fetch_omnipath()
        if not df.empty:
            df.to_parquet(cache, index=False)
            return df
    except Exception:
        pass  # fall through to next source

    # ── 2. Local PhosphoSitePlus ─────────────────────────────────────────────
    try:
        df = _load_local_phosphositeplus()
        if not df.empty:
            df.to_parquet(cache, index=False)
            return df
    except Exception:
        pass  # fall through to bundled fallback

    # ── 3. Bundled fallback ───────────────────────────────────────────────────
    return _load_fallback()
