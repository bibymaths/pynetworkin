"""pipeline.py — adapter layer between cli.py and NetworKIN.py.

This module provides the run_pipeline(...) function expected by cli.py.
It currently delegates execution to the legacy NetworKIN.py script as a
subprocess, then normalizes/copies its output into the location requested
by the Typer CLI.

Assumptions:
- NetworKIN.py is located in the project root next to this file.
- NetworKIN.py writes its main TSV results to: results/<input_filename>.result.tsv
- The CLI currently passes a FASTA file as input.
- Output formats supported here:
    - tsv: copies the NetworKIN result TSV
    - sif: converts TSV rows to a simple edge list
"""

from __future__ import annotations

import csv
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any


def _coerce_species(species: int | str) -> str:
    """Return species as a string accepted by NetworKIN.py."""
    return str(species)


def _build_command(
    networkin_script: Path,
    input_file: Path,
    species: int | str,
    refresh: bool,
    verbose: bool,
) -> list[str]:
    """Build the subprocess command for NetworKIN.py."""
    cmd = [
        sys.executable,
        str(networkin_script),
        _coerce_species(species),
        str(input_file),
    ]

    if refresh:
        cmd.append("--refresh")
    if verbose:
        cmd.append("--verbose")

    return cmd


def _expected_networkin_output(input_file: Path) -> Path:
    """Return the default output path used by NetworKIN.py."""
    return Path("results") / f"{input_file.name}.result.tsv"


def _count_prediction_rows(tsv_path: Path) -> int:
    """Count actual prediction rows in the NetworKIN TSV output.

    Skips:
    - comment/header lines beginning with '#'
    - diagnostic counter lines like 'c_np = 123'
    - blank lines
    """
    count = 0
    with tsv_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                continue
            if stripped.startswith("c_") and "=" in stripped:
                continue
            count += 1
    return count


def _count_recovered_rows(tsv_path: Path) -> int:
    """Count rows marked as recovered in the TSV, if the column exists."""
    recovered = 0

    with tsv_path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(
            (row for row in handle if row.strip() and not row.startswith("c_")),
            delimiter="\t",
        )

        if not reader.fieldnames:
            return 0

        if "recovered" not in reader.fieldnames:
            return 0

        for row in reader:
            value = str(row.get("recovered", "")).strip().lower()
            if value in {"true", "1", "yes"}:
                recovered += 1

    return recovered


def _copy_tsv_output(src: Path, dst: Path) -> None:
    """Copy TSV output to the requested destination."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, dst)


def _write_sif_from_tsv(src_tsv: Path, dst_sif: Path) -> int:
    """Convert NetworKIN TSV output to a simple SIF-like edge list.

    Format:
        <substrate>\tpredicted_by\t<kinase>

    Returns the number of edges written.
    """
    dst_sif.parent.mkdir(parents=True, exist_ok=True)
    written = 0

    with src_tsv.open("r", encoding="utf-8", errors="replace", newline="") as infile, \
         dst_sif.open("w", encoding="utf-8", newline="") as outfile:
        reader = csv.DictReader(
            (row for row in infile if row.strip() and not row.startswith("c_")),
            delimiter="\t",
        )

        if not reader.fieldnames:
            return 0

        for row in reader:
            substrate = str(row.get("Name", "")).strip()
            kinase = str(
                row.get("Kinase/Phosphatase/Phospho-binding domain", "")
            ).strip()

            if not substrate or not kinase:
                continue

            outfile.write(f"{substrate}\tpredicted_by\t{kinase}\n")
            written += 1

    return written


def run_pipeline(
    input_file: str,
    output_path: str,
    output_format: str = "tsv",
    refresh: bool = False,
    use_kg_embedding: bool = False,
    string_score: int = 400,
    species: int = 9606,
    verbose: bool = False,
) -> dict[str, Any]:
    """Run the legacy NetworKIN pipeline and return a CLI summary.

    Parameters are shaped to match cli.py.

    Notes:
    - use_kg_embedding and string_score are currently accepted for compatibility
      with cli.py, but are not yet consumed by NetworKIN.py.
    - output_format must be 'tsv' or 'sif'.
    """
    del use_kg_embedding
    del string_score

    input_path = Path(input_file).expanduser().resolve()
    final_output = Path(output_path).expanduser().resolve()
    output_format = output_format.lower().strip()

    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    if output_format not in {"tsv", "sif"}:
        raise ValueError(f"Unsupported output format: {output_format}")

    project_root = Path(__file__).resolve().parent
    networkin_script = project_root / "NetworKIN.py"

    if not networkin_script.exists():
        raise FileNotFoundError(
            f"NetworKIN.py not found at expected location: {networkin_script}"
        )

    # NetworKIN.py writes relative paths like results/<input>.result.tsv
    # so run it from the project root.
    cmd = _build_command(
        networkin_script=networkin_script,
        input_file=input_path,
        species=species,
        refresh=refresh,
        verbose=verbose,
    )

    proc = subprocess.run(
        cmd,
        cwd=project_root,
        capture_output=True,
        text=True,
    )

    if proc.returncode != 0:
        stderr = proc.stderr.strip()
        stdout = proc.stdout.strip()
        message = stderr or stdout or "NetworKIN.py failed without error output."
        raise RuntimeError(message)

    raw_tsv = project_root / _expected_networkin_output(input_path)

    if not raw_tsv.exists():
        raise FileNotFoundError(
            "NetworKIN.py completed but expected result file was not found: "
            f"{raw_tsv}"
        )

    total_predictions = _count_prediction_rows(raw_tsv)
    recovered_predictions = _count_recovered_rows(raw_tsv)

    if output_format == "tsv":
        _copy_tsv_output(raw_tsv, final_output)
        motif_scored = total_predictions
    else:
        motif_scored = _write_sif_from_tsv(raw_tsv, final_output)

    return {
        "total": total_predictions,
        "motif_scored": motif_scored,
        "recovered": recovered_predictions,
        "output_file": str(final_output),
        "stdout": proc.stdout,
        "stderr": proc.stderr,
    }