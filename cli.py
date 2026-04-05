"""cli.py – Typer + Rich CLI entry-point for pynetworkin."""

from __future__ import annotations

import time
from pathlib import Path

import typer
from rich import box
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

console = Console()

APP_NAME = "pynetworkin"
VERSION = "0.1.0"

BANNER = f"""
[bold white on dark_blue]  pynetworkin v{VERSION}  [/]
[dim]Kinase–substrate network prediction[/]
[dim]Author: Abhinav Mishra[/]
"""

app = typer.Typer(
    name=APP_NAME,
    help="[bold cyan]pynetworkin[/] — kinase–substrate network prediction.",
    rich_markup_mode="rich",
    no_args_is_help=True,
)


@app.command()
def predict(
    input_file: Path = typer.Argument(..., help="Input FASTA or phosphosite file."),  # noqa: B008
    output: Path | None = typer.Option(None, "--output", "-o", help="Output file path."),  # noqa: B008
    format: str = typer.Option("tsv", "--format", "-f", help="Output format: tsv or sif."),
    refresh: bool = typer.Option(False, "--refresh", "-r", help="Force refresh of cached data."),
    use_kg_embedding: bool = typer.Option(
        False, "--use-kg-embedding", help="Use KG embedding scores."
    ),
    string_score: int = typer.Option(
        400, "--string-score", help="STRING combined score threshold (0-1000)."
    ),
    species: int = typer.Option(
        9606, "--species", help="NCBI taxonomy ID (default: 9606 = human)."
    ),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable verbose logging."),
) -> None:
    """[bold green]Predict[/] kinase–substrate interactions from an input file."""
    console.print(Panel(BANNER, border_style="dark_blue", expand=False))

    if not input_file.exists():
        console.print(f"[bold red]Error:[/] input file not found: {input_file}")
        raise typer.Exit(code=1)

    if output is None:
        output = input_file.with_suffix(f".networkin.{format}")

    t0 = time.time()

    from rich.progress import Progress, SpinnerColumn, TextColumn, TimeElapsedColumn

    results: dict = {}

    with Progress(
        SpinnerColumn(style="cyan"),
        TextColumn("[bold cyan]{task.description}"),
        TimeElapsedColumn(),
        console=console,
        transient=False,
    ) as progress:
        task = progress.add_task("Running pipeline", total=None)

        try:
            import csv
            import shutil

            from NetworKIN import AppConfig
            from NetworKIN import run_pipeline as core_run_pipeline

            del use_kg_embedding
            del string_score

            project_root = Path(__file__).resolve().parent
            output_format = format.lower().strip()
            final_output = output.expanduser().resolve()

            if output_format not in {"tsv", "sif"}:
                raise ValueError(f"Unsupported output format: {output_format}")

            progress.update(task, description="Preparing configuration")

            config = AppConfig(
                organism=str(species),
                fasta_path=str(input_file.expanduser().resolve()),
                sites_path=None,
                datadir=str(project_root / "data"),
                blast_dir=str(project_root / "bin") + "/",
                verbose=verbose,
                refresh=refresh,
                result_dir=str(project_root / "results"),
                temp_dir=str(project_root / "tmp"),
            )

            progress.update(task, description="Executing NetworKIN pipeline")
            core_results = core_run_pipeline(config)

            progress.update(task, description="Collecting output")
            raw_tsv = Path(core_results["output_path"]).expanduser().resolve()

            if not raw_tsv.exists():
                raise FileNotFoundError(
                    f"NetworKIN reported success but output file was not found: {raw_tsv}"
                )

            def _count_prediction_rows(tsv_path: Path) -> int:
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
                recovered_count = 0
                with tsv_path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
                    reader = csv.DictReader(
                        (row for row in handle if row.strip() and not row.startswith("c_")),
                        delimiter="\t",
                    )

                    if not reader.fieldnames or "recovered" not in reader.fieldnames:
                        return 0

                    for row in reader:
                        value = str(row.get("recovered", "")).strip().lower()
                        if value in {"true", "1", "yes"}:
                            recovered_count += 1
                return recovered_count

            def _write_sif_from_tsv(src_tsv: Path, dst_sif: Path) -> int:
                dst_sif.parent.mkdir(parents=True, exist_ok=True)
                written = 0

                with (
                    src_tsv.open("r", encoding="utf-8", errors="replace", newline="") as infile,
                    dst_sif.open("w", encoding="utf-8", newline="") as outfile,
                ):
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

            progress.update(task, description="Summarizing predictions")
            total_predictions = _count_prediction_rows(raw_tsv)
            recovered_predictions = _count_recovered_rows(raw_tsv)

            progress.update(task, description="Writing final output")
            if output_format == "tsv":
                final_output.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(raw_tsv, final_output)
                motif_scored = total_predictions
            else:
                motif_scored = _write_sif_from_tsv(raw_tsv, final_output)

            results = {
                "total": total_predictions,
                "motif_scored": motif_scored,
                "recovered": recovered_predictions,
                "output_file": str(final_output),
            }

            progress.update(task, description="[bold green]Pipeline complete[/]", completed=1)

        except Exception as exc:
            progress.update(task, description="[bold red]Pipeline failed[/]")
            console.print(f"[bold red]Pipeline error:[/] {exc}")
            raise typer.Exit(code=1) from exc

    elapsed = time.time() - t0

    table = Table(title="Prediction Summary", box=box.ROUNDED, border_style="cyan")
    table.add_column("Metric", style="bold cyan")
    table.add_column("Value", style="bold white")
    table.add_row("Total predictions", str(results.get("total", 0)))
    table.add_row("Motif-scored", str(results.get("motif_scored", 0)))
    table.add_row("Recovered (FN)", str(results.get("recovered", 0)))
    table.add_row("Output file", str(results.get("output_file", output)))
    table.add_row("Elapsed time", f"{elapsed:.2f}s")
    console.print(table)


@app.command()
def info() -> None:
    """Show package and runtime information."""
    import os
    import platform
    import sys

    console.print(Panel(BANNER, border_style="dark_blue", expand=False))

    table = Table(title="Runtime Information", box=box.ROUNDED, border_style="cyan")
    table.add_column("Key", style="bold cyan")
    table.add_column("Value", style="white")

    table.add_row("Package", f"pynetworkin v{VERSION}")
    table.add_row("Python", sys.version.split()[0])
    table.add_row("Platform", platform.platform())
    table.add_row(
        "NETWORKIN_DATABASE_URL",
        os.environ.get("NETWORKIN_DATABASE_URL", "[dim]not set[/]"),
    )
    table.add_row(
        "NETWORKIN_CACHE_DIR",
        os.environ.get("NETWORKIN_CACHE_DIR", ".cache (default)"),
    )

    try:
        import typer as _typer

        table.add_row("typer", _typer.__version__)
    except Exception:
        pass
    try:
        import rich as _rich

        table.add_row("rich", _rich.__version__)
    except Exception:
        pass
    try:
        import pandas as _pd

        table.add_row("pandas", _pd.__version__)
    except Exception:
        pass
    try:
        import numpy as _np

        table.add_row("numpy", _np.__version__)
    except Exception:
        pass

    console.print(table)


@app.command()
def cache(
    clear: bool = typer.Option(False, "--clear", help="Clear all cached files."),
) -> None:
    """Show cache contents or clear the cache."""
    import os
    import shutil
    from pathlib import Path

    cache_dir = Path(os.environ.get("NETWORKIN_CACHE_DIR", ".cache"))

    if clear:
        if cache_dir.exists():
            shutil.rmtree(cache_dir)
            console.print(f"[bold green]Cache cleared:[/] {cache_dir}")
        else:
            console.print("[dim]Cache directory does not exist.[/]")
        return

    if not cache_dir.exists() or not any(cache_dir.iterdir()):
        console.print(f"[dim]Cache is empty or does not exist: {cache_dir}[/]")
        return

    table = Table(title=f"Cache Contents ({cache_dir})", box=box.ROUNDED, border_style="cyan")
    table.add_column("File", style="bold cyan")
    table.add_column("Size", style="white")
    table.add_column("Modified", style="dim")

    import datetime

    for f in sorted(cache_dir.iterdir()):
        if f.is_file():
            size = f.stat().st_size
            mtime = datetime.datetime.fromtimestamp(f.stat().st_mtime).strftime("%Y-%m-%d %H:%M")
            size_str = f"{size / 1024:.1f} KB" if size < 1_048_576 else f"{size / 1_048_576:.1f} MB"
            table.add_row(f.name, size_str, mtime)

    console.print(table)


if __name__ == "__main__":
    app()
