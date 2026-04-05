"""Package-safe resource helpers for pynetworkin.

All bundled data files are accessed via :mod:`importlib.resources` so the
package works correctly both in editable installs and from a built wheel.
"""

from __future__ import annotations

from collections.abc import Generator
from contextlib import contextmanager
from importlib.resources import as_file, files
from pathlib import Path
from typing import IO

_PKG = files("pynetworkin")


@contextmanager
def conversion_parquet_path(mode: str) -> Generator[Path, None, None]:
    """Yield a filesystem path to the bundled conversion-table Parquet file.

    Parameters
    ----------
    mode:
        ``"direct"`` or any other value (selects ``"indirect"``).

    Yields
    ------
    Path
        Guaranteed to exist for the duration of the ``with`` block.
    """
    name = "conversion_direct.parquet" if mode == "direct" else "conversion_indirect.parquet"
    with as_file(_PKG / "data" / name) as p:
        yield p


def open_fallback_phosphosites() -> IO[str]:
    """Return an open text file-like for the bundled phosphosites fallback TSV."""
    return (_PKG / "data" / "fallback" / "phosphosites_sample.tsv").open("rt", encoding="utf-8")


def open_fallback_string() -> IO[str]:
    """Return an open text file-like for the bundled STRING fallback TSV."""
    return (_PKG / "data" / "fallback" / "string_sample.tsv").open("rt", encoding="utf-8")
