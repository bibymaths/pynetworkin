# PyNetworKIN

**PyNetworKIN** is a Bayesian kinase–substrate prediction pipeline for
phosphoproteomics. It integrates sequence-motif scoring (via
[NetPhorest](http://netphorest.info)) with protein-interaction context (via the
[STRING](https://string-db.org) network) to predict which kinases, phosphatases, or
phospho-binding domains are responsible for observed phosphorylation events.

This repository is a modernised Python 3 port of the original NetworKIN 3.0 tool
(Linding, Jensen, Horn & Kim, 2005–2013), extended to support STRING v12 protein
interaction data.

---

## Features

- Predicts kinase/phosphatase/phospho-binding domain substrates from FASTA +
  phosphosite input.
- Supports human (9606) and yeast (4932) proteomes.
- Accepts multiple phosphosite input formats: NetworKIN TSV, ProteomeDiscoverer,
  MaxQuant, and custom formats.
- Integrates sequence motif posterior probabilities with STRING network proximity
  scores using pre-calibrated Bayesian likelihood-ratio tables.
- Outputs per-site predictions as a TSV file in the `results/` directory.

---

## Requirements

| Dependency | Version | Notes |
|---|---|---|
| Python | ≥ 3.10 | |
| NumPy | ≥ 1.26 | |
| Pandas | ≥ 2.2 | |
| pynetphorest | ≥ 0.1.1 | Motif scoring atlas |
| NCBI BLAST+ | ≥ 2.9 | `blastp` must be on `PATH` or supplied via `--blast-dir` |

---

## Installation

```bash
pip install -e .
```

Or with the optional knowledge-graph extras:

```bash
pip install -e ".[kg]"
```

---

## Usage

### CLI

```bash
pynetworkin predict <FASTA-file> [options]
```

| Argument / Option | Default | Description |
|---|---|---|
| `FASTA-file` | (required) | Input FASTA or phosphosite file |
| `--output` / `-o` | `<input>.networkin.tsv` | Output file path |
| `--format` / `-f` | `tsv` | Output format: `tsv` or `sif` |
| `--species` | `9606` | NCBI taxonomy ID (`9606` = human, `4932` = yeast) |
| `--refresh` / `-r` | off | Force re-fetch of cached network data |
| `--verbose` / `-v` | off | Enable verbose logging |

### Example

```bash
pynetworkin predict data_MaxQuant_sample/test.fasta --output results/test.networkin.tsv
```

Results are written to `results/<fasta-filename>.result.tsv`.

### Other commands

```bash
pynetworkin info       # Show runtime/package information
pynetworkin cache      # Show cache contents
pynetworkin cache --clear  # Clear cached network data
```

### Python API

```python
from pynetworkin import AppConfig, run_pipeline

config = AppConfig(
    organism="9606",
    fasta_path="data_MaxQuant_sample/test.fasta",
    sites_path=None,
    datadir="data",
    blast_dir="",
)
results = run_pipeline(config)
print(results["prediction_count"], "predictions written to", results["output_path"])
```

---

## Input formats

### FASTA file

Standard FASTA format. Protein IDs are taken as everything between `>` and the
first `_` on the header line.

### Sites file (auto-detected)

| Format | Detection | Description |
|---|---|---|
| NetworKIN TSV | 3-column TSV | `protein_id \t position \t residue` |
| ProteomeDiscoverer | 2-column | `protein_id \t phosphopeptide` (phosphosites in lowercase) |
| MaxQuant | Column header `Proteins` + `Leading` | Direct MaxQuant phosphosite output |
| Space-separated | column 2 = `phospho` | Space-separated with residue+position in col 2 |

---

## Output format

Results TSV columns:

| Column | Description |
|---|---|
| Name | Target protein ID |
| Position | Phosphosite position in the protein |
| Tree | NetPhorest tree (KIN, SH2, PTP, 1433, …) |
| Motif Group | NetPhorest classifier group |
| Kinase/Phosphatase/Phospho-binding domain | Predicted enzyme |
| NetworKIN score | Integrated Bayesian score (≥ 0.02 reported) |
| Motif probability | Raw NetPhorest posterior |
| STRING score | STRING best-path proximity score |
| Target STRING ID | Ensembl protein ID of the substrate |
| Kinase STRING ID | Ensembl protein ID of the enzyme |
| Target Name | Human-readable substrate name |
| Kinase Name | Human-readable enzyme name |
| Target description | STRING functional description of substrate |
| Kinase description | STRING functional description of enzyme |
| Peptide sequence window | ±7 aa window around the phosphosite |
| Intermediate nodes | Best-path intermediate proteins in STRING |
| recovered | `True` if recovered by the false-negative recovery step |
| recovery_method | Method used for recovery (e.g. `context_proximity`) |

---

## Repository structure

```
src/
  pynetworkin/          # Core pipeline package
    __init__.py         # Public API (AppConfig, run_pipeline)
    networkin.py        # Main pipeline: AppConfig, run_pipeline, detect_site_file_type, …
    motif_scoring.py    # pynetphorest batch scorer wrapper
    graph_scoring.py    # STRING network context scoring & prediction ranking
    likelihood.py       # Bayesian likelihood conversion tables
    logger.py           # Loguru/Rich logging wrapper
    output.py           # TSV / Cytoscape SIF output writers
    recovery.py         # False-negative recovery via network proximity
    cli.py              # Typer CLI entry-point
    inputs/
      phosphosites.py   # OmniPath / PhosphoSitePlus / fallback fetcher
      string_network.py # STRING flat-file / REST API / fallback fetcher
scripts/
  backup.py                  # Legacy NetworKIN 3.0 reference script (Python 3 port)
  cleanup_HGNC_mapping.py    # HGNC symbol–Ensembl ID reconciliation utility
  generate_sample_data.py    # Generate offline fallback data files
  migrate_to_parquet.py      # Migrate legacy .txt conversion tables → Parquet
data/
  conversion_direct.parquet   # Pre-built likelihood tables (direct STRING paths)
  conversion_indirect.parquet # Pre-built likelihood tables (indirect STRING paths)
  fallback/                   # Bundled offline sample data
  string_data/                # STRING interaction flat files
tests/
  conftest.py             # pytest path setup (adds src/ to sys.path)
  test_motif_scoring.py
  test_output.py
  test_recovery.py
  test_networkin.py       # Tests for load_conversion_tables, detect_site_file_type, run_pipeline
```

See [ARCHITECTURE.md](ARCHITECTURE.md) for a detailed description of the execution flow.

---

## Data sources

- **pynetphorest**: kinase-group motif models (Python package).
- **STRING v12**: human protein interactions and sequences.
  Downloaded from [string-db.org](https://string-db.org).
- **OmniPath**: phosphorylation site reference data (fetched live, cached locally).

---

This repository provides a modern reimplementation of the NetworKIN framework.

- Original NetworKIN was described in:
  Linding et al., Cell 2007

- This implementation:
  - Does NOT reuse original NetworKIN source code
  - Replaces NetPhorest with pynetphorest
  - Uses a rewritten likelihood model
  - Implements a new modular pipeline

---

License: MIT
