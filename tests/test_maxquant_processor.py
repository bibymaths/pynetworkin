"""Tests for MaxQuantProcessor – MaxQuant Site Table processing pipeline."""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from pynetworkin.inputs.maxquant_processor import (
    MaxQuantProcessor,
    ProcessingReport,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def processor() -> MaxQuantProcessor:
    return MaxQuantProcessor(rate_limit_delay=0.0)


SAMPLE_FASTA = """>sp|P12345|PROT_HUMAN Protein
MSEQENCELINES
"""

MINIMAL_SITE_TABLE_ROWS = [
    # (Proteins, Positions within proteins, ...)
    "CON__P01044-1;CON__P01045-1\t331;329",
    "ENSP00000449404;CON__P05787\t49;21",
    "P12345\t100",
    "NP_001234.1\t50",
    "REV__Q99999\t1",
]


def _make_site_table_file(tmp_path: Path, rows: list[str] | None = None) -> Path:
    if rows is None:
        rows = MINIMAL_SITE_TABLE_ROWS
    p = tmp_path / "Phospho (STY)Sites.txt"
    header = "Proteins\tPositions within proteins"
    content = "\n".join([header] + rows) + "\n"
    p.write_text(content, encoding="utf-8")
    return p


# ---------------------------------------------------------------------------
# ProteinIDInfo / clean_protein_ids
# ---------------------------------------------------------------------------


class TestCleanProteinIds:
    def test_uniprot_bare(self, processor: MaxQuantProcessor) -> None:
        result = processor.clean_protein_ids("P12345")
        assert len(result) == 1
        r = result[0]
        assert r.cleaned_id == "P12345"
        assert r.id_type == "uniprot"
        assert not r.is_contaminant
        assert not r.is_reverse

    def test_contaminant_stripped(self, processor: MaxQuantProcessor) -> None:
        result = processor.clean_protein_ids("CON__P01044-1")
        assert len(result) == 1
        r = result[0]
        assert r.cleaned_id == "P01044"
        assert r.is_contaminant
        assert r.isoform == "1"

    def test_reverse_stripped(self, processor: MaxQuantProcessor) -> None:
        result = processor.clean_protein_ids("REV__P12345-1")
        assert len(result) == 1
        r = result[0]
        assert r.is_reverse
        assert r.cleaned_id == "P12345"

    def test_swissprot_format(self, processor: MaxQuantProcessor) -> None:
        result = processor.clean_protein_ids("sp|P12345|PROT_HUMAN")
        assert len(result) == 1
        assert result[0].cleaned_id == "P12345"
        assert result[0].id_type == "uniprot"

    def test_trembl_format(self, processor: MaxQuantProcessor) -> None:
        result = processor.clean_protein_ids("tr|Q67890|Q67890_HUMAN")
        assert len(result) == 1
        assert result[0].cleaned_id == "Q67890"
        assert result[0].id_type == "uniprot"

    def test_ensembl_format(self, processor: MaxQuantProcessor) -> None:
        result = processor.clean_protein_ids("ENSP00000123456")
        assert len(result) == 1
        assert result[0].cleaned_id == "ENSP00000123456"
        assert result[0].id_type == "ensembl"

    def test_refseq_version_removed(self, processor: MaxQuantProcessor) -> None:
        result = processor.clean_protein_ids("NP_001234.1")
        assert len(result) == 1
        assert result[0].cleaned_id == "NP_001234"
        assert result[0].id_type == "refseq"

    def test_multiple_semicolon_separated(self, processor: MaxQuantProcessor) -> None:
        result = processor.clean_protein_ids("CON__P01044-1;CON__P01045-1;Q2KJ62")
        assert len(result) == 3
        cleaned = [r.cleaned_id for r in result]
        assert "P01044" in cleaned
        assert "P01045" in cleaned
        assert "Q2KJ62" in cleaned

    def test_mixed_ids_from_sample(self, processor: MaxQuantProcessor) -> None:
        cell = "ENSP00000449404;CON__P05787;ENSP00000293308"
        result = processor.clean_protein_ids(cell)
        assert len(result) == 3
        assert result[0].id_type == "ensembl"
        assert result[1].is_contaminant
        assert result[2].id_type == "ensembl"

    def test_problem_statement_example_1(self, processor: MaxQuantProcessor) -> None:
        """CON__P12345-1;CON__P67890-2 → P12345;P67890."""
        result = processor.clean_protein_ids("CON__P12345-1;CON__P67890-2")
        assert result[0].cleaned_id == "P12345"
        assert result[1].cleaned_id == "P67890"

    def test_problem_statement_example_2(self, processor: MaxQuantProcessor) -> None:
        """sp|P12345|PROT_HUMAN;tr|Q67890|Q67890_HUMAN → P12345;Q67890."""
        result = processor.clean_protein_ids("sp|P12345|PROT_HUMAN;tr|Q67890|Q67890_HUMAN")
        assert result[0].cleaned_id == "P12345"
        assert result[1].cleaned_id == "Q67890"

    def test_problem_statement_example_3(self, processor: MaxQuantProcessor) -> None:
        """ENSP00000123456;ENSP00000789012 remain as-is."""
        result = processor.clean_protein_ids("ENSP00000123456;ENSP00000789012")
        assert all(r.id_type == "ensembl" for r in result)

    def test_problem_statement_example_4(self, processor: MaxQuantProcessor) -> None:
        """NP_001234.1;XP_005678.1 → NP_001234;XP_005678."""
        result = processor.clean_protein_ids("NP_001234.1;XP_005678.1")
        assert result[0].cleaned_id == "NP_001234"
        assert result[1].cleaned_id == "XP_005678"

    def test_empty_string(self, processor: MaxQuantProcessor) -> None:
        assert processor.clean_protein_ids("") == []

    def test_whitespace_only(self, processor: MaxQuantProcessor) -> None:
        assert processor.clean_protein_ids("  ;  ") == []


# ---------------------------------------------------------------------------
# _clean_proteins_column
# ---------------------------------------------------------------------------


class TestCleanProteinsColumn:
    def test_removes_contaminants_and_reverses(self, processor: MaxQuantProcessor) -> None:
        cleaned = processor._clean_proteins_column("CON__P01044-1;P12345;REV__Q99999")
        ids = cleaned.split(";")
        assert "P12345" in ids
        assert not any("CON__" in i for i in ids)
        assert not any("REV__" in i for i in ids)

    def test_deduplicates(self, processor: MaxQuantProcessor) -> None:
        # Both ENSP IDs are different so no dedup expected, but if same acc appears twice
        cleaned = processor._clean_proteins_column("P12345;P12345")
        assert cleaned.count("P12345") == 1

    def test_all_contaminants_returns_empty(self, processor: MaxQuantProcessor) -> None:
        assert processor._clean_proteins_column("CON__P01044;REV__P99999") == ""


# ---------------------------------------------------------------------------
# ProcessingReport
# ---------------------------------------------------------------------------


class TestProcessingReport:
    def test_to_dict_structure(self) -> None:
        r = ProcessingReport(input_file="test.txt", total_rows=10)
        d = r.to_dict()
        assert d["input_file"] == "test.txt"
        assert d["total_rows"] == 10
        assert "id_types" in d
        assert set(d["id_types"]) == {"uniprot", "refseq", "ensembl", "other"}


# ---------------------------------------------------------------------------
# process_site_table – offline tests (no network)
# ---------------------------------------------------------------------------


class TestProcessSiteTable:
    def test_missing_input_file(self, processor: MaxQuantProcessor, tmp_path: Path) -> None:
        result = processor.process_site_table(
            input_file=tmp_path / "nonexistent.txt",
            output_dir=tmp_path / "out",
        )
        report = result["report"]
        assert report.errors, "Expected errors for missing file"

    def test_missing_proteins_column(
        self, processor: MaxQuantProcessor, tmp_path: Path
    ) -> None:
        p = tmp_path / "bad.txt"
        p.write_text("Other\tColumns\nA\tB\n", encoding="utf-8")

        with patch.object(processor, "fetch_sequences_from_uniprot", return_value={}):
            result = processor.process_site_table(input_file=p, output_dir=tmp_path / "out")

        report = result["report"]
        assert any("Proteins" in e for e in report.errors)

    def test_output_files_created(
        self, processor: MaxQuantProcessor, tmp_path: Path
    ) -> None:
        site_file = _make_site_table_file(tmp_path)
        out_dir = tmp_path / "out"

        # Patch network calls so test is offline
        def _fake_fetch(ids: list[str], id_types: dict | None = None) -> dict[str, str]:
            return {pid: f">sp|{pid}|FAKE_HUMAN Fake\nMSEQ\n" for pid in ids}

        with patch.object(processor, "fetch_sequences_from_uniprot", side_effect=_fake_fetch):
            result = processor.process_site_table(input_file=site_file, output_dir=out_dir)

        assert result["fasta_path"].exists()
        assert result["sites_path"].exists()
        assert result["mapping_path"].exists()
        assert result["report_path"].exists()

    def test_report_json_valid(
        self, processor: MaxQuantProcessor, tmp_path: Path
    ) -> None:
        site_file = _make_site_table_file(tmp_path)
        out_dir = tmp_path / "out"

        def _fake_fetch(ids: list[str], id_types: dict | None = None) -> dict[str, str]:
            return {}

        with patch.object(processor, "fetch_sequences_from_uniprot", side_effect=_fake_fetch):
            result = processor.process_site_table(input_file=site_file, output_dir=out_dir)

        with result["report_path"].open() as fh:
            data = json.load(fh)
        assert "total_rows" in data
        assert "id_types" in data

    def test_id_mapping_csv_has_correct_columns(
        self, processor: MaxQuantProcessor, tmp_path: Path
    ) -> None:
        site_file = _make_site_table_file(tmp_path)
        out_dir = tmp_path / "out"

        def _fake_fetch(ids: list[str], id_types: dict | None = None) -> dict[str, str]:
            return {}

        with patch.object(processor, "fetch_sequences_from_uniprot", side_effect=_fake_fetch):
            result = processor.process_site_table(input_file=site_file, output_dir=out_dir)

        df = pd.read_csv(result["mapping_path"])
        assert "original_id" in df.columns
        assert "cleaned_id" in df.columns

    def test_contaminants_counted_in_report(
        self, processor: MaxQuantProcessor, tmp_path: Path
    ) -> None:
        site_file = _make_site_table_file(tmp_path)
        out_dir = tmp_path / "out"

        def _fake_fetch(ids: list[str], id_types: dict | None = None) -> dict[str, str]:
            return {}

        with patch.object(processor, "fetch_sequences_from_uniprot", side_effect=_fake_fetch):
            result = processor.process_site_table(input_file=site_file, output_dir=out_dir)

        report: ProcessingReport = result["report"]
        # Sample rows include CON__ and REV__ entries
        assert report.contaminants_removed > 0
        assert report.reverse_removed > 0

    def test_cleaned_sites_excludes_contaminants(
        self, processor: MaxQuantProcessor, tmp_path: Path
    ) -> None:
        site_file = _make_site_table_file(tmp_path)
        out_dir = tmp_path / "out"

        def _fake_fetch(ids: list[str], id_types: dict | None = None) -> dict[str, str]:
            return {}

        with patch.object(processor, "fetch_sequences_from_uniprot", side_effect=_fake_fetch):
            result = processor.process_site_table(input_file=site_file, output_dir=out_dir)

        cleaned_df = pd.read_csv(result["sites_path"], sep="\t")
        for val in cleaned_df["Proteins"]:
            assert "CON__" not in str(val)
            assert "REV__" not in str(val)

    def test_fasta_content(
        self, processor: MaxQuantProcessor, tmp_path: Path
    ) -> None:
        site_file = _make_site_table_file(tmp_path)
        out_dir = tmp_path / "out"

        def _fake_fetch(ids: list[str], id_types: dict | None = None) -> dict[str, str]:
            return {"P12345": ">sp|P12345|PROT_HUMAN Fake\nMSEQ\n"}

        with patch.object(processor, "fetch_sequences_from_uniprot", side_effect=_fake_fetch):
            result = processor.process_site_table(input_file=site_file, output_dir=out_dir)

        fasta_content = result["fasta_path"].read_text()
        assert ">sp|P12345|PROT_HUMAN" in fasta_content

    def test_real_sample_data(
        self, processor: MaxQuantProcessor, tmp_path: Path
    ) -> None:
        """Smoke-test with the actual sample_MQ.res file if present."""
        sample = Path(__file__).parent.parent / "data_MaxQuant_sample" / "sample_MQ.res"
        if not sample.exists():
            pytest.skip("sample_MQ.res not present")

        def _fake_fetch(ids: list[str], id_types: dict | None = None) -> dict[str, str]:
            return {}

        with patch.object(processor, "fetch_sequences_from_uniprot", side_effect=_fake_fetch):
            result = processor.process_site_table(input_file=sample, output_dir=tmp_path / "out")

        report: ProcessingReport = result["report"]
        assert report.total_rows > 0
        assert report.unique_original_ids > 0


# ---------------------------------------------------------------------------
# fetch_sequences_from_uniprot – unit tests with mocked ProtMapper
# ---------------------------------------------------------------------------


class TestFetchSequences:
    def test_returns_fasta_for_valid_id(self, processor: MaxQuantProcessor) -> None:
        result_df = pd.DataFrame([{
            "From": "P12345",
            "Entry": "P12345",
            "Entry Name": "PROT_HUMAN",
            "Protein names": "Test protein",
            "Sequence": "MSEQENCELINES",
        }])

        with patch.object(processor._mapper, "get", return_value=(result_df, [])):
            result = processor.fetch_sequences_from_uniprot(
                ["P12345"], id_types={"P12345": "uniprot"}
            )

        assert "P12345" in result
        assert result["P12345"].startswith(">sp|P12345|")

    def test_returns_empty_for_failed_id(self, processor: MaxQuantProcessor) -> None:
        empty_df = pd.DataFrame(columns=["From", "Entry", "Entry Name", "Protein names", "Sequence"])

        with patch.object(processor._mapper, "get", return_value=(empty_df, ["INVALID"])):
            result = processor.fetch_sequences_from_uniprot(
                ["INVALID"], id_types={"INVALID": "uniprot"}
            )

        assert result == {}
