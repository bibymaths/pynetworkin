"""Tests for core pynetworkin pipeline functions.

Covers:
* load_conversion_tables  – Parquet-based likelihood table loader
* detect_site_file_type   – phosphosite file format detector
* run_pipeline            – integration smoke test (import + config only)
"""

from __future__ import annotations

import io
import tempfile
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest

from pynetworkin.networkin import (
    MAX_QUANT_DIRECT_OUTPUT_FILE,
    MS_MCMC_FILE,
    NETWORKIN_SITE_FILE,
    PROTEOME_DISCOVERER_SITE_FILE,
    LEGACY_SITE_FILE,
    AppConfig,
    NetworkinError,
    detect_site_file_type,
    load_conversion_tables,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_CONV_TABLE_TEXT = (
    "Score\tLower bound\tUpper bound\tLikelihood\tTPR\tFPR\tPPV\tFDR\t"
    "No. positives\tNo. negatives\n"
    "2.50000\t2.00000\t3.00000\t5.00000\t0.80\t0.10\t0.75\t0.25\t40\t10\n"
    "1.00000\t0.50000\t1.50000\t2.00000\t0.50\t0.20\t0.60\t0.40\t25\t20\n"
    "0.10000\t0.00000\t0.20000\t0.50000\t0.20\t0.40\t0.30\t0.70\t10\t40\n"
)


def _make_parquet_bytes(species: str = "human") -> bytes:
    """Build an in-memory Parquet file containing one sample conversion table."""
    df = pd.DataFrame(
        [
            {
                "path_mode": "direct",
                "score_kind": "motif",
                "species": species,
                "tree": "KIN",
                "player_name": "PKA_group",
                "raw_data": _CONV_TABLE_TEXT,
            },
            {
                "path_mode": "direct",
                "score_kind": "string",
                "species": species,
                "tree": "KIN",
                "player_name": "PKA_group",
                "raw_data": _CONV_TABLE_TEXT,
            },
        ]
    )
    buf = io.BytesIO()
    df.to_parquet(buf, index=False)
    buf.seek(0)
    return buf.read()


# ---------------------------------------------------------------------------
# load_conversion_tables
# ---------------------------------------------------------------------------


def test_load_conversion_tables_builds_nested_dict(tmp_path: Path) -> None:
    """load_conversion_tables must build tables[species][tree][name][score_kind]."""
    parquet_file = tmp_path / "conversion_direct.parquet"
    parquet_file.write_bytes(_make_parquet_bytes("human"))

    tables = load_conversion_tables(str(parquet_file), species="human")

    assert "human" in tables, "top-level species key missing"
    assert "KIN" in tables["human"], "tree key 'KIN' missing"
    assert "PKA_group" in tables["human"]["KIN"], "player_name 'PKA_group' missing"
    assert "motif" in tables["human"]["KIN"]["PKA_group"], "'motif' score_kind missing"
    assert "string" in tables["human"]["KIN"]["PKA_group"], "'string' score_kind missing"


def test_load_conversion_tables_species_filter(tmp_path: Path) -> None:
    """Only rows matching the requested species should appear in the result."""
    # Build a Parquet with both human and yeast rows
    df = pd.DataFrame(
        [
            {
                "path_mode": "direct",
                "score_kind": "motif",
                "species": "human",
                "tree": "KIN",
                "player_name": "PKA_group",
                "raw_data": _CONV_TABLE_TEXT,
            },
            {
                "path_mode": "direct",
                "score_kind": "motif",
                "species": "yeast",
                "tree": "KIN",
                "player_name": "Cdc28",
                "raw_data": _CONV_TABLE_TEXT,
            },
        ]
    )
    parquet_file = tmp_path / "conversion_direct.parquet"
    df.to_parquet(parquet_file, index=False)

    tables = load_conversion_tables(str(parquet_file), species="human")

    assert "human" in tables
    assert "yeast" not in tables, "yeast rows must not appear when species='human'"


def test_load_conversion_tables_entries_are_lists(tmp_path: Path) -> None:
    """Each leaf value must be a non-empty list of CConvEntry objects."""
    parquet_file = tmp_path / "conversion_direct.parquet"
    parquet_file.write_bytes(_make_parquet_bytes("human"))

    tables = load_conversion_tables(str(parquet_file), species="human")

    motif_tbl = tables["human"]["KIN"]["PKA_group"]["motif"]
    assert isinstance(motif_tbl, list), "conversion table must be a list"
    assert len(motif_tbl) > 0, "conversion table must not be empty"
    # Each entry must have a .score and .L attribute
    entry = motif_tbl[0]
    assert hasattr(entry, "score"), "CConvEntry must have a 'score' attribute"
    assert hasattr(entry, "L"), "CConvEntry must have an 'L' attribute"


def test_load_conversion_tables_sorted_descending(tmp_path: Path) -> None:
    """Conversion table entries must be sorted by score descending."""
    parquet_file = tmp_path / "conversion_direct.parquet"
    parquet_file.write_bytes(_make_parquet_bytes("human"))

    tables = load_conversion_tables(str(parquet_file), species="human")
    motif_tbl = tables["human"]["KIN"]["PKA_group"]["motif"]

    scores = [e.score for e in motif_tbl]
    assert scores == sorted(scores, reverse=True), "entries must be sorted by score descending"


# ---------------------------------------------------------------------------
# detect_site_file_type
# ---------------------------------------------------------------------------


def _write_site_file(content: str, suffix: str = ".tsv") -> str:
    """Write *content* to a temp file and return its path."""
    tf = tempfile.NamedTemporaryFile(
        mode="w", suffix=suffix, delete=False, encoding="utf-8"
    )
    tf.write(content)
    tf.flush()
    tf.close()
    return tf.name


def test_detect_networkin_site_file(tmp_path: Path) -> None:
    """Three-column TSV → NETWORKIN_SITE_FILE."""
    path = _write_site_file("P12345\t42\tS\nP12345\t77\tT\n")
    assert detect_site_file_type(path) == NETWORKIN_SITE_FILE


def test_detect_proteome_discoverer_file(tmp_path: Path) -> None:
    """Two-column TSV → PROTEOME_DISCOVERER_SITE_FILE."""
    path = _write_site_file("P12345\tABCpDEF\n")
    assert detect_site_file_type(path) == PROTEOME_DISCOVERER_SITE_FILE


def test_detect_maxquant_file(tmp_path: Path) -> None:
    """MaxQuant header line → MAX_QUANT_DIRECT_OUTPUT_FILE."""
    header = "\t".join(["Proteins", "col2", "col3", "col4", "Leading", "col6"])
    path = _write_site_file(header + "\n")
    assert detect_site_file_type(path) == MAX_QUANT_DIRECT_OUTPUT_FILE


def test_detect_legacy_site_file(tmp_path: Path) -> None:
    """Space-separated line with >3 tokens where tokens[1]=='phospho' → LEGACY_SITE_FILE."""
    # Must be space-separated (not tab) and have >3 tokens so it doesn't
    # match the len==3 NETWORKIN_SITE_FILE branch first.
    path = _write_site_file("P12345 phospho S42 extra_column\n")
    assert detect_site_file_type(path) == LEGACY_SITE_FILE


def test_detect_ms_mcmc_file(tmp_path: Path) -> None:
    """Filename starting with 'MS' and single-token header → MS_MCMC_FILE."""
    # Single-token line doesn't match the len==2 or len==3 branches, so the
    # filename-based check fires.
    path = _write_site_file("P12345\n", suffix=".tsv")
    ms_path = str(Path(path).parent / "MS_test_output.tsv")
    Path(path).rename(ms_path)
    assert detect_site_file_type(ms_path) == MS_MCMC_FILE


def test_detect_unknown_format_raises(tmp_path: Path) -> None:
    """An unrecognised format must raise NetworkinError."""
    path = _write_site_file("col1\tcol2\tcol3\tcol4\n")
    with pytest.raises(NetworkinError):
        detect_site_file_type(path)


def test_detect_empty_file_raises(tmp_path: Path) -> None:
    """An empty file must raise NetworkinError."""
    path = _write_site_file("")
    with pytest.raises(NetworkinError):
        detect_site_file_type(path)


# ---------------------------------------------------------------------------
# run_pipeline – import + config smoke test
# ---------------------------------------------------------------------------


def test_run_pipeline_importable() -> None:
    """run_pipeline must be importable from pynetworkin.networkin."""
    from pynetworkin.networkin import run_pipeline  # noqa: F401


def test_appconfig_defaults() -> None:
    """AppConfig must accept the minimum required fields and expose properties."""
    cfg = AppConfig(
        organism="9606",
        fasta_path="/tmp/test.fasta",
        sites_path=None,
        datadir="/tmp/data",
        blast_dir="",
    )
    assert cfg.species_name == "human"
    assert cfg.fasta_stem == "test.fasta"
    assert cfg.blast_output_path.endswith(".blast.out")


def test_run_pipeline_raises_on_missing_fasta(tmp_path: Path) -> None:
    """run_pipeline must raise when the FASTA file does not exist."""
    from pynetworkin.networkin import run_pipeline

    cfg = AppConfig(
        organism="9606",
        fasta_path=str(tmp_path / "nonexistent.fasta"),
        sites_path=None,
        datadir=str(tmp_path),
        blast_dir="",
    )
    with pytest.raises(Exception):
        run_pipeline(cfg)
