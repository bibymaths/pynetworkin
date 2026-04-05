"""conftest.py – pytest configuration for the pynetworkin test suite.

Adds the ``src/`` directory to *sys.path* so that ``import pynetworkin``
resolves correctly without requiring an editable install.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Ensure the src/ layout is importable from any working directory.
SRC_DIR = Path(__file__).resolve().parent / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))
