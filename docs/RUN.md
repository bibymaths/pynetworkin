# RUN.md — NetworKIN Execution Manual

A complete guide to setting up, configuring, and running the NetworKIN
kinase–substrate prediction pipeline.

---

## Table of Contents

1. [Requirements](#1-requirements)
2. [Installation](#2-installation)
3. [Data sources](#3-data-sources)
4. [Minimal offline run](#4-minimal-offline-run)
5. [Full run with live data](#5-full-run-with-live-data)
6. [Output options](#6-output-options)
7. [Testing](#7-testing)
8. [Troubleshooting](#8-troubleshooting)
9. [Directory structure](#9-directory-structure)

---

## 1. Requirements

| Dependency | Minimum version | Notes |
|---|---|---|
| Python | 3.10 | Type-union syntax (`X \| Y`) used throughout |
| BLAST+ | any | `makeblastdb` and `blastp` must be on `$PATH` |

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
```

---

## 3. Data sources

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

## 4. Minimal offline run

This run uses only bundled data and requires no internet access.

```bash
pynetworkin predict data/fallback/test_input.fasta \
  --output results.tsv --format tsv
```

The fallback data files in `data/fallback/` are used automatically when live
sources are unavailable.

---

## 5. Full run with live data

```bash
# Force re-download of all live data (OmniPath + STRING)
pynetworkin predict data/fallback/test_input.fasta \
  --output results.tsv --refresh
```

### With a custom phosphosite file

```bash
pynetworkin predict my_sequences.fasta \
  --sites data/my_phosphosites.tsv \
  --output results.tsv
```

### With ProteomeDiscoverer or MaxQuant input

```bash
# ProteomeDiscoverer format — pass the FASTA and the peptide list (.txt)
# File format is auto-detected (2-column protein-ID / peptide file)
pynetworkin predict proteome.fasta \
  --sites proteome_discoverer_peptides.txt \
  --output results.tsv

# MaxQuant format — pass the FASTA and the MaxQuant phosphosite table (.res)
# File format is auto-detected (tab-separated with "Proteins" and
# "Positions within proteins" header columns)
pynetworkin predict proteome.fasta \
  --sites maxquant_phospho_sites.res \
  --output results.tsv
```

---

## 6. Output options

### TSV (default)

```bash
pynetworkin predict sequences.fasta --output results.tsv --format tsv
```

Output columns:
`Name`, `Position`, `Tree`, `Motif Group`, `Kinase/Phosphatase/Phospho-binding domain`,
`NetworKIN score`, `Motif probability`, `STRING score`,
`Target STRING ID`, `Kinase STRING ID`, `Target Name`, `Kinase Name`,
`Target description`, `Kinase description`, `Peptide sequence window`,
`Intermediate nodes`, `recovered`, `recovery_method`

### Cytoscape-compatible edge list

```bash
pynetworkin predict sequences.fasta --output network.tsv --format cytoscape
```

Output columns: `source`, `interaction`, `target`, `weight`

### Recovery (enabled by default)

False-negative recovery is enabled by default.  It uses STRING graph distance
to rescue kinase–substrate pairs where the motif score fell below threshold but
the proteins are close in the interaction network.

---

## 7. Testing

```bash
# Run all tests with coverage
pytest tests/ -v --cov=.

# Run a specific test file
pytest tests/test_motif_scoring.py -v

# Quick smoke test using fallback data
pynetworkin predict data/fallback/test_input.fasta --output /tmp/test.tsv
wc -l /tmp/test.tsv     # should be > 1
```

### Regenerating fallback data

If you have internet access and want to refresh the bundled sample data:

```bash
python scripts/generate_sample_data.py
```

---

## 8. Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `OmniPath fetch failed` | No internet / OmniPath down | Pipeline falls back to local PSP file, then bundled TSV automatically |
| `STRING flat file not found` | Flat file missing or wrong path | Check `data/string_data/9606.links.v12.0.tsv.gz` exists; set `NETWORKIN_STRING_FLAT_FILE` if needed |
| Empty results (`0 predictions`) | All kinase scores below threshold, or empty FASTA | Check FASTA has ≥1 S/T/Y residue per protein |
| `BLAST not found` | `blastp` not on PATH | Install BLAST+: `sudo apt install ncbi-blast+` or `conda install -c bioconda blast` |
| `ModuleNotFoundError: pynetphorest` | Package not installed | Run `pip install -e ".[dev]"` |
| `KeyError: 'motif'` | Missing or corrupt Parquet conversion table | Check `data/conversion_direct.parquet` and `data/conversion_indirect.parquet` exist; run `python scripts/migrate_to_parquet.py` to regenerate |

---

## 9. Directory structure

```
pynetworkin/
├── src/pynetworkin/             # Core pipeline package
│   ├── cli.py                   # CLI entry-point (pynetworkin predict …)
│   ├── networkin.py             # Pipeline orchestration, site-file parsers
│   ├── motif_scoring.py         # pynetphorest batch scorer
│   ├── graph_scoring.py         # STRING context scoring and ranking
│   ├── likelihood.py            # Bayesian likelihood conversion tables
│   ├── output.py                # TSV / Cytoscape output writers
│   ├── recovery.py              # False-negative recovery
│   └── inputs/
│       ├── phosphosites.py      # Phosphosite fetching (OmniPath → PSP → fallback)
│       ├── string_network.py    # STRING network (flat file → REST → fallback)
│       └── maxquant_processor.py # MaxQuant site-table pre-processor
├── data/
│   ├── fallback/                # Bundled offline data (always present)
│   │   ├── phosphosites_sample.tsv
│   │   ├── string_sample.tsv
│   │   └── test_input.fasta
│   ├── string_data/             # STRING v12 flat files
│   │   └── 9606.links.v12.0.tsv.gz
│   ├── conversion_direct.parquet    # Likelihood tables (direct paths)
│   ├── conversion_indirect.parquet  # Likelihood tables (indirect paths)
│   └── 9606.protein.sequences.v12.0.fa  # Human proteome FASTA (BLAST DB)
├── tests/                       # pytest test suite
├── scripts/
│   ├── backup.py                # Legacy NetworKIN 3.0 reference (not executed)
│   ├── cleanup_HGNC_mapping.py  # HGNC symbol reconciliation utility
│   ├── generate_sample_data.py  # Regenerate data/fallback/ from live sources
│   └── migrate_to_parquet.py    # Migrate legacy .txt conversion tables → Parquet
├── .env.example                 # Environment variable template
├── requirements.txt             # Runtime dependencies
├── NON_CODE_AUDIT.md            # Non-Python file inventory & decisions
└── RUN.md                       # This file
```

---

*Last updated: 2026-04-05*
