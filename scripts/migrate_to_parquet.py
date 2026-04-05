#!/usr/bin/env python3
"""Migrate likelihood conversion table text files to Parquet format.

This script reads every conversion-table ``.txt`` file from the two legacy
directories:

* ``data/likelihood_conversion_table_direct/``
* ``data/likelihood_conversion_table_indirect/``

and writes them into two compact, snappy-compressed Parquet files:

* ``data/conversion_direct.parquet``
* ``data/conversion_indirect.parquet``

Each row in the resulting Parquet files represents one original text file and
carries its parsed metadata together with the raw file contents so that no
numerical precision is lost on round-trip.

Columns
-------
path_mode   : str  – ``"direct"`` or ``"indirect"``
score_kind  : str  – ``"string"`` (from STRING) or ``"motif"`` (from NetPhorest)
species     : str  – species extracted from filename (e.g. ``"human"``)
tree        : str  – domain/tree tag (e.g. ``"KIN"``, ``"SH2"``)
player_name : str  – kinase / domain name (e.g. ``"CDK7"``, ``"general"``)
raw_data    : str  – complete file content (header + data rows) as a string

Usage
-----
Run from the repository root::

    python scripts/migrate_to_parquet.py [--datadir data]

The ``--datadir`` argument defaults to ``data/`` relative to the repository
root (the directory that contains the
``likelihood_conversion_table_direct/`` sub-directories).
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

FILENAME_PATTERN = re.compile(
    r"conversion_tbl_([a-z]+)_smooth_([a-z]+)_([A-Z0-9]+)_([a-zA-Z0-9_/-]+)"
)

_DIR_TO_MODE: dict[str, str] = {
    "likelihood_conversion_table_direct": "direct",
    "likelihood_conversion_table_indirect": "indirect",
}

_OUTPUT_NAMES: dict[str, str] = {
    "direct": "conversion_direct.parquet",
    "indirect": "conversion_indirect.parquet",
}

PARQUET_COMPRESSION = "snappy"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _parse_filename(stem: str) -> tuple[str, str, str, str] | None:
    """Return ``(score_src, species, tree, player_name)`` or ``None``."""
    match = FILENAME_PATTERN.findall(stem)
    if not match:
        return None
    score_src, species, tree, player_name = match[0]
    return score_src, species, tree, player_name


def _score_kind(score_src: str) -> str:
    return "string" if score_src == "string" else "motif"


def _collect_records(source_dir: Path, path_mode: str) -> list[dict]:
    """Iterate over ``source_dir`` and return a list of row dicts."""
    records: list[dict] = []

    for filepath in sorted(source_dir.glob("conversion_tbl_*_smooth*")):
        if not filepath.is_file():
            continue

        parsed = _parse_filename(filepath.stem)
        if parsed is None:
            print(f"  [WARN] skipping unmatched filename: {filepath.name}", file=sys.stderr)
            continue

        score_src, species, tree, player_name = parsed

        raw_data = filepath.read_text(encoding="utf-8")

        records.append(
            {
                "path_mode": path_mode,
                "score_kind": _score_kind(score_src),
                "species": species,
                "tree": tree,
                "player_name": player_name,
                "raw_data": raw_data,
            }
        )

    return records


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def migrate(datadir: str | os.PathLike) -> None:
    """Perform the full migration for *datadir*.

    Args:
        datadir: Path to the ``data/`` directory that contains the
            ``likelihood_conversion_table_direct/`` and
            ``likelihood_conversion_table_indirect/`` sub-directories.
    """
    datadir = Path(datadir)

    # Group files by path_mode
    mode_records: dict[str, list[dict]] = {"direct": [], "indirect": []}

    for subdir_name, path_mode in _DIR_TO_MODE.items():
        source_dir = datadir / subdir_name
        if not source_dir.is_dir():
            print(
                f"[WARN] directory not found, skipping: {source_dir}", file=sys.stderr
            )
            continue

        records = _collect_records(source_dir, path_mode)
        mode_records[path_mode].extend(records)
        print(f"  Collected {len(records)} files from {source_dir}")

    # Write one Parquet file per path_mode
    for path_mode, records in mode_records.items():
        if not records:
            print(f"[WARN] no records for path_mode={path_mode!r}, skipping.", file=sys.stderr)
            continue

        df = pd.DataFrame(records)
        out_path = datadir / _OUTPUT_NAMES[path_mode]
        df.to_parquet(out_path, index=False, compression=PARQUET_COMPRESSION)
        print(f"  Wrote {len(df)} rows → {out_path}  ({out_path.stat().st_size:,} bytes)")


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--datadir",
        default="data",
        help="Path to the data directory (default: %(default)s)",
    )
    return parser


if __name__ == "__main__":
    args = _build_arg_parser().parse_args()
    migrate(args.datadir)
