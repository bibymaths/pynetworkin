# RUN.md — NetworKIN Execution Manual

A complete guide to setting up, configuring, and running the NetworKIN
kinase–substrate prediction pipeline.

---

## Table of Contents

1. [Requirements](#1-requirements)
2. [Installation](#2-installation)
3. [Database setup](#3-database-setup)
4. [Data sources](#4-data-sources)
5. [Minimal offline run](#5-minimal-offline-run)
6. [Full run with live data](#6-full-run-with-live-data)
7. [Output options](#7-output-options)
8. [Testing](#8-testing)
9. [Troubleshooting](#9-troubleshooting)
10. [Directory structure](#10-directory-structure)

---

## 1. Requirements

| Dependency | Minimum version | Notes |
|---|---|---|
| Python | 3.10 | Type-union syntax (`X \| Y`) used throughout |
| BLAST+ | any | `makeblastdb` and `blastp` must be on `$PATH` |
| SQLite **or** PostgreSQL | SQLite 3.x / PG 13+ | SQLite is fine for development |
| (Optional) Numba | 0.57+ | JIT-compiled recovery step; auto-disabled if absent |

---

## 2. Installation

```bash
# Clone the repository
git clone https://github.com/bibymaths/pynetworkin.git
cd pynetworkin

# Create and activate a virtual environment (recommended)
python -m venv .venv
source .venv/bin/activate          # Linux / macOS
# .venv\Scripts\activate           # Windows

# Install the package with development dependencies
pip install -e ".[dev]"

# Copy environment template
cp .env.example .env
# Edit .env and set NETWORKIN_DATABASE_URL (at minimum)
```

---

## 3. Database setup

NetworKIN uses Alembic for schema migrations.

### SQLite (development)

```bash
export NETWORKIN_DATABASE_URL="sqlite:///./dev.db"
alembic upgrade head
```

### PostgreSQL (production)

```bash
export NETWORKIN_DATABASE_URL="postgresql://user:password@localhost:5432/networkin"
alembic upgrade head
```

> **Tip:** Store these in your `.env` file and `source .env` before running.

---

## 4. Data sources

### Phosphorylation sites (priority order)

1. **OmniPath** (primary, public REST API)  
   Fetched automatically from `https://omnipathdb.org/enzsub`.  
   Requires internet access.  Cached for 7 days in `.cache/`.

2. **Local PhosphoSitePlus** (manual download required)  
   PhosphoSitePlus requires a registered account — download is not automated.  
   1. Log in at <https://www.phosphosite.org/downloads/>  
   2. Download `Phosphorylation_site_dataset.gz`  
   3. Set the path: `export NETWORKIN_PSP_LOCAL_FILE=/path/to/file.tsv`

3. **Bundled fallback TSV** (always available, offline)  
   `data/fallback/phosphosites_sample.tsv` — 20 curated human phosphosites.  
   Used automatically when all live sources fail.

### STRING network (priority order)

1. **Flat file** (primary, bundled with repo)  
   `data/string_data/9606.links.v12.0.tsv.gz`  
   Override path: `export NETWORKIN_STRING_FLAT_FILE=/path/to/file.tsv.gz`

2. **STRING REST API** (fallback, requires internet)  
   Used automatically for small queries when the flat file is unavailable.

3. **Bundled sample** (offline fallback)  
   `data/fallback/string_sample.tsv` — 20 curated interaction edges.

### Kinase motifs

Motif scoring uses the **`pynetphorest`** Python package exclusively.  
No external files need to be downloaded — the package bundles all model data.

---

## 5. Minimal offline run

This run uses only bundled data and requires no internet access or database.

```bash
python -m cli predict data/fallback/test_input.fasta \
  --output results.tsv --format tsv
```

The fallback data files in `data/fallback/` are used automatically when live
sources are unavailable.

---

## 6. Full run with live data

```bash
# Force re-download of all live data (OmniPath + STRING)
python -m cli predict data/fallback/test_input.fasta \
  --output results.tsv --refresh
```

### With a custom phosphosite file

```bash
python -m cli predict my_sequences.fasta \
  --sites data/my_phosphosites.tsv \
  --output results.tsv
```

### With ProteomeDiscoverer or MaxQuant input

```bash
# ProteomeDiscoverer format
python -m cli predict proteome.fasta \
  --sites proteome_discoverer_output.txt \
  --format pd \
  --output results.tsv

# MaxQuant format
python -m cli predict proteome.fasta \
  --sites maxquant_phospho_sites.txt \
  --format mq \
  --output results.tsv
```

---

## 7. Output options

### TSV (default)

```bash
python -m cli predict sequences.fasta --output results.tsv --format tsv
```

Output columns:  
`Name`, `Position`, `Tree`, `Motif Group`, `Kinase/Phosphatase/Phospho-binding domain`,  
`NetworKIN score`, `Motif probability`, `STRING score`,  
`Target STRING ID`, `Kinase STRING ID`, `Target Name`, `Kinase Name`,  
`Target description`, `Kinase description`, `Peptide sequence window`,  
`Intermediate nodes`, `recovered`, `recovery_method`

### Cytoscape-compatible edge list

```bash
python -m cli predict sequences.fasta --output network.tsv --format cytoscape
```

Output columns: `source`, `interaction`, `target`, `weight`

### Recovery (enabled by default)

False-negative recovery is enabled by default.  It uses STRING graph distance
to rescue kinase–substrate pairs where the motif score fell below threshold but
the proteins are close in the interaction network.

```bash
# Disable recovery
python -m cli predict sequences.fasta --output results.tsv --no-recovery
```

### KG embedding (optional, requires additional setup)

```bash
python -m cli predict sequences.fasta --output results.tsv --kg-embedding
```

---

## 8. Testing

```bash
# Run all tests with coverage
pytest tests/ -v --cov=.

# Run a specific test file
pytest tests/test_motif_scoring.py -v

# Quick smoke test using fallback data
python -m cli predict data/fallback/test_input.fasta --output /tmp/test.tsv
wc -l /tmp/test.tsv     # should be > 1
```

### Regenerating fallback data

If you have internet access and want to refresh the bundled sample data:

```bash
python scripts/generate_sample_data.py
```

---

## 9. Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `NETWORKIN_DATABASE_URL is not set` | Missing env var | Run `export NETWORKIN_DATABASE_URL=sqlite:///./dev.db` |
| `alembic: command not found` | Alembic not installed | Run `pip install -e ".[dev]"` |
| `OmniPath fetch failed` | No internet / OmniPath down | Pipeline falls back to local PSP file, then bundled TSV automatically |
| `STRING flat file not found` | Flat file missing or wrong path | Check `data/string_data/9606.links.v12.0.tsv.gz` exists; set `NETWORKIN_STRING_FLAT_FILE` if needed |
| Empty results (`0 predictions`) | All kinase scores below threshold, or empty FASTA | Check FASTA has ≥1 S/T/Y residue per protein; lower `--threshold` |
| Slow first run | Numba JIT compilation | Normal; subsequent runs are fast.  Set `NUMBA_DISABLE_JIT=1` to skip JIT |
| `BLAST not found` | `blastp` not on PATH | Install BLAST+: `sudo apt install ncbi-blast+` or `conda install -c bioconda blast` |
| `ModuleNotFoundError: pynetphorest` | Package not installed | Run `pip install -e ".[dev]"` |
| `KeyError: 'motif'` | Old conversion table format | Check `data/likelihood_conversion_table_direct/` contains the expected `.txt` files |

---

## 10. Directory structure

```
pynetworkin/
├── NetworKIN.py                 # Main orchestration script & CLI entry point
├── cli/                         # CLI module (python -m cli predict …)
├── inputs/
│   ├── phosphosites.py          # Phosphosite fetching (OmniPath → PSP → fallback)
│   └── string_network.py        # STRING network (flat file → REST → fallback)
├── data/
│   ├── fallback/                # Bundled offline data (always present)
│   │   ├── phosphosites_sample.tsv
│   │   ├── string_sample.tsv
│   │   └── test_input.fasta
│   ├── string_data/             # STRING v12 flat files (large; not in git LFS)
│   │   └── 9606.links.v12.0.tsv.gz
│   ├── likelihood_conversion_table_direct/   # Per-kinase score→likelihood tables
│   └── 9606.protein.sequences.v12.0.fa      # Human proteome FASTA (BLAST DB)
├── tests/                       # pytest test suite
├── scripts/
│   └── generate_sample_data.py  # Regenerate data/fallback/ from live sources
├── legacy_c/                    # Archived C source & headers (not compiled)
├── legacy_web/                  # Archived Perl scripts & shell batch runners
├── .env.example                 # Environment variable template
├── requirements.txt             # Runtime dependencies
├── NON_CODE_AUDIT.md            # Non-Python file inventory & decisions
└── RUN.md                       # This file
```

---

*Last updated: 2026-04-05*
