"""cli.py – Typer + Rich CLI entry-point for pynetworkin."""

from __future__ import annotations

import time
from pathlib import Path
from typing import Optional

import typer
from rich import box
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

console = Console()

APP_NAME = "pynetworkin"
VERSION = "0.1.0"

BANNER = """
[bold white on dark_blue]  pynetworkin v{version}  [/]
[dim]Kinase–substrate network prediction[/]
[dim]Author: Abhinav Mishra[/]
""".format(version=VERSION)

app = typer.Typer(
    name=APP_NAME,
    help="[bold cyan]pynetworkin[/] — kinase–substrate network prediction.",
    rich_markup_mode="rich",
    no_args_is_help=True,
)


@app.command()
def predict(
    input_file: Path = typer.Argument(..., help="Input FASTA or phosphosite file."),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Output file path."),
    format: str = typer.Option("tsv", "--format", "-f", help="Output format: tsv or sif."),
    refresh: bool = typer.Option(False, "--refresh", "-r", help="Force refresh of cached data."),
    use_kg_embedding: bool = typer.Option(False, "--use-kg-embedding", help="Use KG embedding scores."),
    string_score: int = typer.Option(400, "--string-score", help="STRING combined score threshold (0-1000)."),
    species: int = typer.Option(9606, "--species", help="NCBI taxonomy ID (default: 9606 = human)."),
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

    steps = [
        ("Reading input", None),
        ("Fetching phosphosites", None),
        ("Fetching STRING network", None),
        ("Motif scoring", None),
        ("Context scoring", None),
        ("False-negative recovery", None),
        ("Writing output", None),
    ]

    from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TimeElapsedColumn

    results: dict = {}

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TimeElapsedColumn(),
        console=console,
        transient=True,
    ) as progress:
        for step_name, _ in steps:
            task = progress.add_task(step_name, total=None)
            progress.advance(task)

        try:
            import pipeline
            results = pipeline.run_pipeline(
                input_file=str(input_file),
                output_path=str(output),
                output_format=format,
                refresh=refresh,
                use_kg_embedding=use_kg_embedding,
                string_score=string_score,
                species=species,
                verbose=verbose,
            )
        except ImportError:
            console.print("[yellow]Warning:[/] pipeline module not found — running in demo mode.")
            results = {
                "total": 0,
                "motif_scored": 0,
                "recovered": 0,
                "output_file": str(output),
            }
        except Exception as exc:
            console.print(f"[bold red]Pipeline error:[/] {exc}")
            raise typer.Exit(code=1)

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
    import platform
    import sys
    import os

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
