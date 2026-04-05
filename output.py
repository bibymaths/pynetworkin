# output.py
#
# Output writers for NetworKIN predictions.
#
# Formats:
#   - TSV: tab-separated, one row per kinase–substrate prediction.
#   - Cytoscape SIF: source \t interaction \t target, one row per edge.
#
# The `recovered` column is new in this refactor.
# Existing tools reading the TSV will ignore unknown trailing columns.
#
# Column order in STANDARD_COLUMNS matches the original NetworKIN output exactly.
# Do NOT reorder these columns — downstream parsers depend on positional order.

import pandas as pd
from pathlib import Path

# Column names match the original NetworKIN CSV output exactly
# (as written by printResult in NetworKIN.py).
STANDARD_COLUMNS = [
    "Name",
    "Position",
    "Tree",
    "Motif Group",
    "Kinase/Phosphatase/Phospho-binding domain",
    "NetworKIN score",
    "Motif probability",
    "STRING score",
    "Target STRING ID",
    "Kinase STRING ID",
    "Target Name",
    "Kinase Name",
    "Target description",
    "Kinase description",
    "Peptide sequence window",
    "Intermediate nodes",
    "recovered",          # new column — False for original predictions, True for recovered
    "recovery_method",    # new column — empty string or method name e.g. "context_proximity"
]


def write_tsv(predictions: list, path) -> None:
    """
    Write predictions to a tab-separated file.
    Columns follow STANDARD_COLUMNS order.
    Missing optional fields are filled with None.
    """
    df = pd.DataFrame(predictions)
    for col in STANDARD_COLUMNS:
        if col not in df.columns:
            df[col] = None
    df = df[STANDARD_COLUMNS]
    df.to_csv(path, sep="\t", index=False)


def write_cytoscape(predictions: list, path) -> None:
    """
    Write Cytoscape-compatible SIF format:
        source \t pp \t target
    where pp indicates a phosphorylation interaction.
    Weighted by NetworKIN score.
    """
    rows = []
    for p in predictions:
        rows.append({
            "source":      p.get("Kinase/Phosphatase/Phospho-binding domain", ""),
            "interaction": "pp",
            "target":      p.get("Name", ""),
            "weight":      p.get("NetworKIN score", 0.0),
        })
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def write_output(predictions: list, path, fmt: str = "tsv") -> None:
    """Dispatcher — call write_tsv or write_cytoscape based on fmt."""
    fmt = fmt.lower().strip()
    if fmt in ("tsv", "tab"):
        write_tsv(predictions, path)
    elif fmt in ("cytoscape", "sif"):
        write_cytoscape(predictions, path)
    else:
        raise ValueError(f"Unknown output format: {fmt!r}. Use 'tsv' or 'cytoscape'.")
