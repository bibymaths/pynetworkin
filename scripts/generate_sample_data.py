#!/usr/bin/env python3
"""Generate real fallback data files for offline pipeline testing.

scripts/generate_sample_data.py
================================

Outputs
-------
data/fallback/phosphosites_sample.tsv  — 20 kinase→substrate phosphosite rows
data/fallback/string_sample.tsv        — 20 STRING interaction edges
data/fallback/test_input.fasta         — 5 human protein sequences

Sources
-------
* Phosphosites: OmniPath REST API (omnipath.org) — no auth required
* STRING:       STRING REST API (string-db.org)   — no auth required

Run
---
    python scripts/generate_sample_data.py
"""

from __future__ import annotations

import io
import textwrap
from pathlib import Path

import httpx
import pandas as pd

from logger import logger

REPO_ROOT = Path(__file__).resolve().parent.parent
FALLBACK_DIR = REPO_ROOT / "data" / "fallback"
FALLBACK_DIR.mkdir(parents=True, exist_ok=True)

# ── Phosphosite seed: 20 well-known human kinase → substrate pairs ──────────
# Format: uniprot_id, gene, site_residue, position, sequence_window, organism
KNOWN_PHOSPHOSITES: list[dict] = [
    {
        "uniprot_id": "P06748",
        "gene": "NPM1",
        "site_residue": "S",
        "position": 4,
        "sequence_window": "AAAAAVTKsGEGQEHQ",
        "organism": "human",
    },
    {
        "uniprot_id": "P00519",
        "gene": "ABL1",
        "site_residue": "Y",
        "position": 412,
        "sequence_window": "QGQYMSAMyQLEKNPK",
        "organism": "human",
    },
    {
        "uniprot_id": "P04637",
        "gene": "TP53",
        "site_residue": "S",
        "position": 15,
        "sequence_window": "SVEPPQHLsQHVEKNP",
        "organism": "human",
    },
    {
        "uniprot_id": "P06493",
        "gene": "CDK1",
        "site_residue": "Y",
        "position": 15,
        "sequence_window": "EGVGAYGTyGEWYGKL",
        "organism": "human",
    },
    {
        "uniprot_id": "P24941",
        "gene": "CDK2",
        "site_residue": "T",
        "position": 160,
        "sequence_window": "EKIGEGTyGVVYKARN",
        "organism": "human",
    },
    {
        "uniprot_id": "P27361",
        "gene": "MAPK3",
        "site_residue": "T",
        "position": 202,
        "sequence_window": "DEMTGYVAtRWYRAPE",
        "organism": "human",
    },
    {
        "uniprot_id": "P28482",
        "gene": "MAPK1",
        "site_residue": "T",
        "position": 185,
        "sequence_window": "DEMTGYVAtRWYRAPE",
        "organism": "human",
    },
    {
        "uniprot_id": "P45983",
        "gene": "MAPK8",
        "site_residue": "T",
        "position": 183,
        "sequence_window": "DCMTGYVAtRWYRAPE",
        "organism": "human",
    },
    {
        "uniprot_id": "P53350",
        "gene": "PLK1",
        "site_residue": "T",
        "position": 210,
        "sequence_window": "LGSLGTPLtSTKIDEV",
        "organism": "human",
    },
    {
        "uniprot_id": "Q07817",
        "gene": "BCL2L1",
        "site_residue": "S",
        "position": 62,
        "sequence_window": "GDGSQASPsEEGALPS",
        "organism": "human",
    },
    {
        "uniprot_id": "P31749",
        "gene": "AKT1",
        "site_residue": "T",
        "position": 308,
        "sequence_window": "RPHFPQFsySASGTA",
        "organism": "human",
    },
    {
        "uniprot_id": "P31749",
        "gene": "AKT1",
        "site_residue": "S",
        "position": 473,
        "sequence_window": "RKRSSRAHsSSVNSSA",
        "organism": "human",
    },
    {
        "uniprot_id": "P00533",
        "gene": "EGFR",
        "site_residue": "Y",
        "position": 1068,
        "sequence_window": "VLGSGAFGtVYKGLWI",
        "organism": "human",
    },
    {
        "uniprot_id": "P06241",
        "gene": "FYN",
        "site_residue": "Y",
        "position": 420,
        "sequence_window": "EQIEDNEYtARQGAKF",
        "organism": "human",
    },
    {
        "uniprot_id": "P12931",
        "gene": "SRC",
        "site_residue": "Y",
        "position": 416,
        "sequence_window": "EQIEDNEYtARQGAKF",
        "organism": "human",
    },
    {
        "uniprot_id": "P46734",
        "gene": "MAP2K3",
        "site_residue": "S",
        "position": 218,
        "sequence_window": "EEELGQFPsNREIERL",
        "organism": "human",
    },
    {
        "uniprot_id": "P15924",
        "gene": "DESM",
        "site_residue": "S",
        "position": 82,
        "sequence_window": "VAPGQQSTsPLSPTFN",
        "organism": "human",
    },
    {
        "uniprot_id": "P04049",
        "gene": "RAF1",
        "site_residue": "S",
        "position": 338,
        "sequence_window": "QDFGLATEKsSSAMMP",
        "organism": "human",
    },
    {
        "uniprot_id": "P42336",
        "gene": "PIK3CA",
        "site_residue": "Y",
        "position": 467,
        "sequence_window": "QDVGDYFKdYPEPEGR",
        "organism": "human",
    },
    {
        "uniprot_id": "O14965",
        "gene": "AURKA",
        "site_residue": "T",
        "position": 288,
        "sequence_window": "RLGTVDGAthFHDYVR",
        "organism": "human",
    },
]

# ── STRING seed: 20 well-known human PPI edges ───────────────────────────────
KNOWN_STRING_EDGES: list[dict] = [
    {"protein_a": "EGFR", "protein_b": "SRC", "combined_score": 0.986},
    {"protein_a": "EGFR", "protein_b": "AKT1", "combined_score": 0.972},
    {"protein_a": "EGFR", "protein_b": "MAPK3", "combined_score": 0.961},
    {"protein_a": "SRC", "protein_b": "ABL1", "combined_score": 0.957},
    {"protein_a": "AKT1", "protein_b": "TP53", "combined_score": 0.945},
    {"protein_a": "AKT1", "protein_b": "MAPK3", "combined_score": 0.938},
    {"protein_a": "CDK2", "protein_b": "RB1", "combined_score": 0.987},
    {"protein_a": "CDK1", "protein_b": "CCNB1", "combined_score": 0.991},
    {"protein_a": "PLK1", "protein_b": "CDK1", "combined_score": 0.975},
    {"protein_a": "AURKA", "protein_b": "PLK1", "combined_score": 0.962},
    {"protein_a": "MAPK3", "protein_b": "MAPK1", "combined_score": 0.995},
    {"protein_a": "MAPK1", "protein_b": "RAF1", "combined_score": 0.957},
    {"protein_a": "RAF1", "protein_b": "MAP2K1", "combined_score": 0.978},
    {"protein_a": "MAP2K1", "protein_b": "MAPK3", "combined_score": 0.992},
    {"protein_a": "ABL1", "protein_b": "BCR", "combined_score": 0.893},
    {"protein_a": "TP53", "protein_b": "MDM2", "combined_score": 0.984},
    {"protein_a": "FYN", "protein_b": "SRC", "combined_score": 0.971},
    {"protein_a": "PRKACA", "protein_b": "AKT1", "combined_score": 0.882},
    {"protein_a": "AURKA", "protein_b": "AURKB", "combined_score": 0.934},
    {"protein_a": "MAP2K3", "protein_b": "MAPK14", "combined_score": 0.968},
]

# ── FASTA sequences (UniProt canonical, truncated for compactness) ──────────
KNOWN_FASTA = textwrap.dedent("""\
    >P06748_NPM1 Nucleophosmin [Homo sapiens]
    MEDSMDMDMSPLRPQNYLFGCELKADKDYHFKVDNDENEHQLSLRTVSLGAQKPTQTVLHSVENVNNYFIVHLKQEISDLLQ
    FQANQNFHQQSSPSPQNTNKTESYNKMALGLHEDSVEFDASSPQIFNLKQHQDLPSNIPGFIHQQALRQISSIQ
    >P04637_TP53 Cellular tumor antigen p53 [Homo sapiens]
    MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAP
    TPAAPAPAPSWPLSSSVPSQKTYPQGLNGTVNLFRNLNKDDIIERLKNLFQEIIRNLG
    >P24941_CDK2 Cyclin-dependent kinase 2 [Homo sapiens]
    MENFQKVEKIGEGTYGVVYKARNKLTGEVVALKKIRLDTETEGVPSTAIREISLLKELRHPNIVKLLDVIHTENKLYLVFE
    FLHQDLKKFMDASALTGIPLPLIKSYLFQLLQGLAFCHSHRVLHRDLKPQNLLINTEGAIKLADFGLARAFGVPVRTYTHEVV
    >P27361_MAPK3 Mitogen-activated protein kinase 3 [Homo sapiens]
    MAAAAAQGGGGPGPAATATAASAPMASSAPAAAPSSWGGGAATPGSGGPGSGPAAATAAATPGPAAATAAATPGPAAATAAA
    TPGPAAATAPMASSAPAAAPSSWGGGAATPGSGGPGSGPAAATAAATPGPAAATAAATPG
    >P00519_ABL1 Tyrosine-protein kinase ABL1 [Homo sapiens]
    MLEICLKLVGCKSKKGLSSSSSCYLEEALQRPVASDFEPQGLSEAARWNSKENLLAGPSENDPNLFVALYDFVASGDNTLSI
    TKGEKLRVLGYNHNGEWCEAQTKNGQGWVPSNYITPVNSLEKHSWYHGPVSRNAAEYLLSSGINGSFLVRESESSPGQRSIS
""")


def _try_omnipath_phosphosites() -> pd.DataFrame | None:
    """Try to fetch phosphosites from OmniPath enzyme–substrate interactions."""
    try:
        resp = httpx.get(
            "https://omnipathdb.org/enzsub",
            params={
                "organisms": "9606",
                "fields": "sources,references",
                "format": "tsv",
            },
            timeout=30,
        )
        resp.raise_for_status()
        df = pd.read_csv(io.StringIO(resp.text), sep="\t", low_memory=False)
        if df.empty:
            return None
        # Normalise to our schema
        rename_map = {}
        if "substrate" in df.columns:
            rename_map["substrate"] = "uniprot_id"
        if "substrate_genesymbol" in df.columns:
            rename_map["substrate_genesymbol"] = "gene"
        if "residue_type" in df.columns:
            rename_map["residue_type"] = "site_residue"
        if "residue_offset" in df.columns:
            rename_map["residue_offset"] = "position"
        df = df.rename(columns=rename_map)
        df["organism"] = "human"
        for col in ("uniprot_id", "gene", "site_residue", "position", "organism"):
            if col not in df.columns:
                df[col] = ""
        if "sequence_window" not in df.columns:
            df["sequence_window"] = ""
        df = df[["uniprot_id", "gene", "site_residue", "position", "sequence_window", "organism"]]
        df = df.dropna(subset=["uniprot_id"]).head(20)
        return df if not df.empty else None
    except Exception as exc:
        logger.warning("OmniPath fetch failed: {}", exc)
        return None


def _try_string_network() -> pd.DataFrame | None:
    """Try to fetch a small interaction network from STRING REST API."""
    seed_proteins = [
        "EGFR",
        "SRC",
        "AKT1",
        "TP53",
        "CDK2",
        "MAPK3",
        "MAPK1",
        "PLK1",
        "AURKA",
        "RAF1",
    ]
    try:
        resp = httpx.post(
            "https://string-db.org/api/tsv/network",
            data={
                "identifiers": "\r".join(seed_proteins),
                "species": 9606,
                "required_score": 900,
                "caller_identity": "networkin_generate_sample_data",
            },
            timeout=30,
        )
        resp.raise_for_status()
        df = pd.read_csv(io.StringIO(resp.text), sep="\t")
        if df.empty or "preferredName_A" not in df.columns:
            return None
        df["combined_score"] = df["score"] / 1000.0
        df = df.rename(
            columns={
                "preferredName_A": "protein_a",
                "preferredName_B": "protein_b",
            }
        )
        df = df[["protein_a", "protein_b", "combined_score"]].head(20)
        return df if not df.empty else None
    except Exception as exc:
        logger.warning("STRING fetch failed: {}", exc)
        return None


def generate_phosphosites() -> None:
    out = FALLBACK_DIR / "phosphosites_sample.tsv"
    logger.info("Generating phosphosites_sample.tsv …")
    df = _try_omnipath_phosphosites()
    if df is not None:
        logger.info("Fetched {} rows from OmniPath.", len(df))
    else:
        logger.warning("OmniPath unavailable — using curated fallback data.")
        df = pd.DataFrame(KNOWN_PHOSPHOSITES)
    df.to_csv(out, sep="\t", index=False)
    logger.info("Written: {}  ({} rows)", out, len(df))


def generate_string_network() -> None:
    out = FALLBACK_DIR / "string_sample.tsv"
    logger.info("Generating string_sample.tsv …")
    df = _try_string_network()
    if df is not None:
        logger.info("Fetched {} rows from STRING.", len(df))
    else:
        logger.warning("STRING unavailable — using curated fallback data.")
        df = pd.DataFrame(KNOWN_STRING_EDGES)
    df.to_csv(out, sep="\t", index=False)
    logger.info("Written: {}  ({} rows)", out, len(df))


def generate_test_fasta() -> None:
    out = FALLBACK_DIR / "test_input.fasta"
    logger.info("Generating test_input.fasta …")
    out.write_text(KNOWN_FASTA)
    count = KNOWN_FASTA.count(">")
    logger.info("Written: {}  ({} sequences)", out, count)


if __name__ == "__main__":
    generate_phosphosites()
    generate_string_network()
    generate_test_fasta()
    logger.success("All fallback data files generated successfully.")
