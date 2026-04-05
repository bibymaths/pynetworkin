"""Tests for package-safe resource loading and STRING download logic.

Covers:
* conversion_parquet_path – yields a real, readable file
* open_fallback_phosphosites / open_fallback_string – package-relative TSV access
* no dependency on repo-root or cwd for bundled resources
* _resolve_flat_file priority order:
    - env-var path provided and exists
    - cached file already present
    - download triggered when nothing cached
    - download failure falls through to REST API / bundled fallback
"""

from __future__ import annotations

import gzip
import os
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _write_minimal_string_gz(path: Path) -> None:
    """Write a minimal STRING protein.links gzip file at *path*."""
    path.parent.mkdir(parents=True, exist_ok=True)
    content = "protein1 protein2 combined_score\n9606.ENSP001 9606.ENSP002 700\n"
    with gzip.open(path, "wt") as fh:
        fh.write(content)


# ---------------------------------------------------------------------------
# resources.py – bundled resource access
# ---------------------------------------------------------------------------


def test_conversion_parquet_path_direct_exists() -> None:
    """conversion_parquet_path('direct') must yield an existing, non-empty file."""
    from pynetworkin.resources import conversion_parquet_path

    with conversion_parquet_path("direct") as p:
        assert p.exists(), f"Bundled direct parquet not found at {p}"
        assert p.stat().st_size > 0, "Bundled direct parquet is empty"


def test_conversion_parquet_path_indirect_exists() -> None:
    """conversion_parquet_path('indirect') must yield an existing, non-empty file."""
    from pynetworkin.resources import conversion_parquet_path

    with conversion_parquet_path("indirect") as p:
        assert p.exists(), f"Bundled indirect parquet not found at {p}"
        assert p.stat().st_size > 0, "Bundled indirect parquet is empty"


def test_conversion_parquet_path_not_repo_relative(tmp_path: Path) -> None:
    """The parquet path must not depend on the repo root or current working directory."""
    from pynetworkin.resources import conversion_parquet_path

    original_cwd = Path.cwd()
    try:
        os.chdir(tmp_path)
        with conversion_parquet_path("direct") as p:
            assert p.exists(), "Parquet not accessible after cwd change"
    finally:
        os.chdir(original_cwd)


def test_conversion_parquet_readable_as_dataframe() -> None:
    """The bundled parquet must be loadable as a pandas DataFrame."""
    from pynetworkin.resources import conversion_parquet_path

    with conversion_parquet_path("direct") as p:
        df = pd.read_parquet(p)
    assert not df.empty, "Bundled direct parquet produced an empty DataFrame"


def test_fallback_phosphosites_tsv_readable() -> None:
    """open_fallback_phosphosites() must return a readable text stream."""
    from pynetworkin.resources import open_fallback_phosphosites

    with open_fallback_phosphosites() as fh:
        df = pd.read_csv(fh, sep="\t", low_memory=False)
    assert not df.empty, "Bundled phosphosites fallback TSV produced an empty DataFrame"


def test_fallback_string_tsv_readable() -> None:
    """open_fallback_string() must return a readable text stream."""
    from pynetworkin.resources import open_fallback_string

    with open_fallback_string() as fh:
        df = pd.read_csv(fh, sep="\t", low_memory=False)
    assert not df.empty, "Bundled STRING fallback TSV produced an empty DataFrame"


def test_fallback_tsvs_not_repo_relative(tmp_path: Path) -> None:
    """Fallback TSVs must be accessible regardless of cwd."""
    from pynetworkin.resources import open_fallback_phosphosites, open_fallback_string

    original_cwd = Path.cwd()
    try:
        os.chdir(tmp_path)
        with open_fallback_phosphosites() as fh:
            assert fh.read(), "phosphosites fallback is empty after cwd change"
        with open_fallback_string() as fh:
            assert fh.read(), "string fallback is empty after cwd change"
    finally:
        os.chdir(original_cwd)


# ---------------------------------------------------------------------------
# string_network.py – _resolve_flat_file priority logic
# ---------------------------------------------------------------------------


def test_resolve_flat_file_env_var_exists(tmp_path: Path) -> None:
    """When NETWORKIN_STRING_FLAT_FILE is set and the file exists, use it."""
    from pynetworkin.inputs import string_network as sn

    flat = tmp_path / "custom_string.tsv.gz"
    _write_minimal_string_gz(flat)

    with patch.dict(os.environ, {"NETWORKIN_STRING_FLAT_FILE": str(flat)}):
        result = sn._resolve_flat_file()
    assert result == flat


def test_resolve_flat_file_env_var_missing_raises(tmp_path: Path) -> None:
    """When env var points to a non-existent file, raise FileNotFoundError."""
    from pynetworkin.inputs import string_network as sn

    with patch.dict(os.environ, {"NETWORKIN_STRING_FLAT_FILE": str(tmp_path / "nope.gz")}):
        with pytest.raises(FileNotFoundError):
            sn._resolve_flat_file()


def test_resolve_flat_file_uses_cache_when_present(tmp_path: Path) -> None:
    """When the cached file exists and is non-empty, return it without downloading."""
    from pynetworkin.inputs import string_network as sn

    cached = tmp_path / "string_data" / sn._STRING_CACHE_FILENAME
    _write_minimal_string_gz(cached)

    env = {
        "NETWORKIN_CACHE_DIR": str(tmp_path),
        "NETWORKIN_STRING_FLAT_FILE": "",
    }
    with patch.dict(os.environ, env, clear=False):
        # Reload the CACHE_DIR-derived constant used inside _resolve_flat_file
        with patch.object(sn, "_downloaded_flat_file_path", return_value=cached):
            result = sn._resolve_flat_file()
    assert result == cached


def test_resolve_flat_file_downloads_when_no_cache(tmp_path: Path) -> None:
    """When no cache exists and no env var is set, _download_string_flat_file is called."""
    from pynetworkin.inputs import string_network as sn

    dest = tmp_path / "string_data" / sn._STRING_CACHE_FILENAME

    def fake_download() -> Path:
        _write_minimal_string_gz(dest)
        return dest

    env = {"NETWORKIN_STRING_FLAT_FILE": ""}
    with patch.dict(os.environ, env, clear=False):
        with patch.object(sn, "_downloaded_flat_file_path", return_value=tmp_path / "missing.gz"):
            with patch.object(sn, "_download_string_flat_file", side_effect=fake_download) as mock_dl:
                result = sn._resolve_flat_file()
    mock_dl.assert_called_once()
    assert result == dest


def test_fetch_string_network_falls_back_to_bundled_on_all_failures(tmp_path: Path) -> None:
    """When flat-file and REST API both fail, bundled sample TSV is returned."""
    from pynetworkin.inputs import string_network as sn

    env = {
        "NETWORKIN_CACHE_DIR": str(tmp_path),
        "NETWORKIN_STRING_FLAT_FILE": "",
    }

    with patch.dict(os.environ, env, clear=False):
        with patch.object(sn, "_resolve_flat_file", side_effect=RuntimeError("no file")):
            with patch.object(sn, "_fetch_rest_api", side_effect=RuntimeError("no api")):
                df = sn.fetch_string_network(proteins=["PROT1"], refresh=True)

    assert isinstance(df, pd.DataFrame)
    # Should match bundled fallback columns
    assert set(df.columns) >= {"protein_a", "protein_b", "combined_score"} or not df.empty


def test_download_string_flat_file_atomic_rename(tmp_path: Path, monkeypatch) -> None:
    """_download_string_flat_file must use atomic rename (no partial files left)."""
    import httpx
    from pynetworkin.inputs import string_network as sn

    dest = tmp_path / "string_data" / sn._STRING_CACHE_FILENAME

    # Build a minimal gzip payload to stream
    import gzip as gz
    import io

    payload = gz.compress(b"protein1 protein2 combined_score\n9606.A 9606.B 500\n")

    mock_response = MagicMock()
    mock_response.raise_for_status = MagicMock()
    mock_response.iter_bytes = MagicMock(return_value=iter([payload]))
    mock_response.__enter__ = MagicMock(return_value=mock_response)
    mock_response.__exit__ = MagicMock(return_value=False)

    monkeypatch.setattr(sn, "_downloaded_flat_file_path", lambda: dest)

    with patch("httpx.stream", return_value=mock_response):
        result = sn._download_string_flat_file()

    assert result == dest
    assert dest.exists()
    assert dest.stat().st_size > 0
    # No temp files left behind
    remaining_tmp = list(dest.parent.glob("*.tmp"))
    assert not remaining_tmp, f"Temp files not cleaned up: {remaining_tmp}"


def test_download_string_flat_file_raises_on_http_error(tmp_path: Path, monkeypatch) -> None:
    """_download_string_flat_file must raise RuntimeError on HTTP failure."""
    import httpx
    from pynetworkin.inputs import string_network as sn

    dest = tmp_path / "string_data" / sn._STRING_CACHE_FILENAME
    monkeypatch.setattr(sn, "_downloaded_flat_file_path", lambda: dest)

    mock_response = MagicMock()
    mock_response.raise_for_status = MagicMock(
        side_effect=httpx.HTTPStatusError("404", request=MagicMock(), response=MagicMock())
    )
    mock_response.__enter__ = MagicMock(return_value=mock_response)
    mock_response.__exit__ = MagicMock(return_value=False)

    with patch("httpx.stream", return_value=mock_response):
        with pytest.raises(RuntimeError, match="Failed to download STRING"):
            sn._download_string_flat_file()
