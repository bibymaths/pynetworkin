# PyNetworKIN — Architecture Overview

## 1. Project Summary

**PyNetworKIN** is a Bayesian kinase–substrate prediction pipeline for
phosphoproteomics. It integrates two complementary sources of evidence:

1. **pynetphorest** – a Python package that scores peptide-sequence motifs using
   per-kinase classifiers, producing a posterior probability for each
   (site, kinase-group) pair.
2. **STRING network** – a pre-computed best-path proximity score between proteins
   in the STRING protein–protein interaction network.

The two evidence streams are combined via pre-calibrated Bayesian likelihood-ratio
conversion tables to produce a single **NetworKIN score** per (phosphosite, kinase/
phosphatase/phospho-binding domain) pair.

A **false-negative recovery** step rescues kinase–substrate pairs missed by the
motif step but well-supported by STRING network proximity.

Supported organisms: **human (9606)** and **yeast (4932)**.

---

## 2. Repository Map

```
.
├── src/pynetworkin/              # Core pipeline package
│   ├── __init__.py               # Public API: AppConfig, run_pipeline
│   ├── cli.py                    # Typer CLI entry-point (pynetworkin predict …)
│   ├── networkin.py              # Main pipeline orchestration + site-file parsers
│   ├── motif_scoring.py          # pynetphorest batch scorer wrapper
│   ├── graph_scoring.py          # STRING context scoring and prediction ranking
│   ├── likelihood.py             # Bayesian likelihood-ratio conversion table utilities
│   ├── logger.py                 # Loguru/Rich logging wrapper
│   ├── output.py                 # TSV / Cytoscape SIF output writers
│   ├── recovery.py               # False-negative recovery via network proximity
│   └── inputs/
│       ├── phosphosites.py       # OmniPath / PhosphoSitePlus / fallback fetcher
│       ├── string_network.py     # STRING flat-file / REST API / fallback fetcher
│       └── maxquant_processor.py # MaxQuant site-table pre-processor
│
├── data/
│   ├── conversion_direct.parquet    # Likelihood tables for direct STRING paths
│   ├── conversion_indirect.parquet  # Likelihood tables for indirect STRING paths
│   ├── group_human_protein_name_map.tsv         # NetPhoREST group → protein name (base)
│   ├── group_human_protein_name_map_curated.tsv # Curated group → name map (v12 aliases)
│   ├── 9606.protein.sequences.v12.0.fa          # STRING v12 human proteome (BLAST DB)
│   ├── 9606.protein.sequences.v9.0.fa           # STRING v9 human proteome (legacy)
│   ├── 9606.protein.aliases.v12.0.txt.gz        # STRING v12 protein aliases
│   ├── 9606.protein.info.v12.0.txt.gz           # STRING v12 protein descriptions
│   ├── 9606.alias_best.v9.0.tsv.gz              # STRING v9 best aliases (human)
│   ├── 9606.text_best.v9.0.tsv.gz               # STRING v9 protein descriptions (human)
│   ├── 4932.*                                   # Equivalent files for yeast
│   ├── string_data/
│   │   ├── 9606.links.v12.0.tsv.gz              # Primary STRING v12 flat file (human)
│   │   ├── 9606.bestpath_0340_0950.v9.tsv.gz    # Legacy STRING v9 best-path scores
│   │   └── 4932.string_000_0170_1000.tsv.gz     # Yeast STRING interactions
│   └── fallback/
│       ├── phosphosites_sample.tsv # Bundled offline phosphosite sample
│       ├── string_sample.tsv       # Bundled offline STRING edge sample
│       └── test_input.fasta        # Bundled offline test FASTA
│
├── scripts/
│   ├── backup.py                  # Legacy NetworKIN 3.0 reference (Python 3 port, not executed)
│   ├── cleanup_HGNC_mapping.py    # HGNC symbol–Ensembl ID reconciliation utility
│   ├── generate_sample_data.py    # Regenerate data/fallback/ from live sources
│   └── migrate_to_parquet.py      # Migrate legacy .txt conversion tables → Parquet
│
├── tests/
│   ├── conftest.py                 # Adds src/ to sys.path (no editable install needed)
│   ├── test_motif_scoring.py
│   ├── test_output.py
│   ├── test_recovery.py
│   ├── test_networkin.py
│   └── test_maxquant_processor.py
│
├── docs/                          # MkDocs documentation source
├── pyproject.toml                 # Package build and tool configuration
└── README.md
```

---

## 3. Execution Flow (End-to-End)

### 3.1 Entry Point

```bash
pynetworkin predict <FASTA-file> [--sites <sites-file>] [options]
```

The CLI is implemented with Typer in `cli.py`.  It constructs an `AppConfig`
dataclass and calls `run_pipeline(config)`.

### 3.2 Input Ingestion

| Stage | Function | Description |
|---|---|---|
| FASTA | `read_fasta_file(path)` | Reads protein sequences; IDs parsed from header (UniProt `sp|ID|` or first whitespace-delimited token) |
| Sites | `detect_site_file_type(path)` | Auto-detects one of 5 formats (see below) |
| Sites | `read_networkin_sites` / `read_proteome_discoverer_sites` / `read_max_quant_sites` | Format-specific parsers |

Detected site-file formats:

| Constant | Description |
|---|---|
| `NETWORKIN_SITE_FILE` | 3-column TSV: `protein \t position \t residue` |
| `PROTEOME_DISCOVERER_SITE_FILE` | 2-column: `protein \t phosphopeptide` (lowercase = phosphosite) |
| `MAX_QUANT_DIRECT_OUTPUT_FILE` | MaxQuant phosphosite output (header starts with `Proteins … Leading`) |
| `LEGACY_SITE_FILE` | Space-separated; second token is `phospho` |
| `MS_MCMC_FILE` | Filename starts with `MS` |

### 3.3 Reference Data Loading

1. **Group → Domain map** (`read_group_to_domain_map`): loads
   `data/group_human_protein_name_map_curated.tsv` (if present, else falls back to
   `group_human_protein_name_map.tsv`), which maps NetPhoREST classifier groups to
   HGNC-approved protein names.
2. **Alias / description hashes** (`read_alias_files`): reads STRING v9 text
   descriptions and STRING v12 protein aliases; builds `string_alias` (STRING ID
   → best name), `string_desc` (STRING ID → description), and `name_hash`
   (STRING ID → list of aliases).

### 3.4 Sequence → STRING Mapping (BLAST)

`map_peptides_to_string(config, id_pos_res, id_seq)`:

- Writes query sequences to a temp FASTA.
- Runs `blastp` against the STRING sequence database
  (`data/9606.protein.sequences.v12.0.fa`).
- The best BLAST hit (sorted by bit-score) is kept; identity < 90 % or
  e-value > 1e-40 trigger a warning.
- Builds bidirectional dictionaries `incoming2string` and `string2incoming`.

### 3.5 STRING Network Loading

`load_string_data(config, string2incoming, string_alias)`:

- Reads `data/string_data/9606.links.v12.0.tsv.gz` (primary flat file).
- Falls back to the STRING REST API for small queries.
- Populates `tree_pred_string_data[substrate_STRING_id][kinase_STRING_id]`
  with `{"_name": name, "_score": float, "_path": path}`.

### 3.6 Motif Scoring (pynetphorest)

`score_sequences(id_seq, id_pos_res)` in `motif_scoring.py`:

- Loads the pynetphorest atlas at import time (once per process).
- For each protein and each S/T/Y site, calls `_core.get_model_posterior`.
- Builds `id_pos_tree_pred[protein_id][position][tree][kinase] = (res, peptide, score)`.

### 3.7 Bayesian Score Integration

`compile_predictions(...)`:

- For each (protein, position, tree, group) tuple from pynetphorest:
  1. Look up the STRING mapping for the substrate.
  2. Find matching STRING edge(s) between the substrate and kinase.
  3. Load the appropriate likelihood-ratio conversion table from the Parquet file
     (`data/conversion_direct.parquet` or `data/conversion_indirect.parquet`).
  4. Convert raw pynetphorest and STRING scores to likelihood ratios via
     `ConvertScore2L(score, conv_tbl)` (linear interpolation).
  5. Compute unified likelihood: `L_final = L_motif × L_string`.
  6. Collect results where `NetworKIN score ≥ 0.02`.

### 3.8 False-Negative Recovery

`recover_predictions(...)` in `networkin.py` (calls `recover_false_negatives`
from `recovery.py`):

- Computes all-pairs shortest-path distances on the STRING network
  (Floyd–Warshall via `networkx`).
- For each candidate (kinase, substrate) pair absent from motif scoring,
  computes `context_score = 1 / (1 + path_dist)`.
- Recovers pairs where `context_score ≥ CONTEXT_RECOVERY_THRESHOLD` (default 0.6).
- Recovered predictions are tagged `recovered=True`, `recovery_method="context_proximity"`.

### 3.9 Output Format

TSV columns:

| Column | Description |
|---|---|
| Name | Target protein ID |
| Position | Phosphosite position |
| Tree | NetPhoREST tree (KIN, SH2, PTP, 1433, …) |
| Motif Group | NetPhoREST classifier group |
| Kinase/Phosphatase/Phospho-binding domain | Enzyme name |
| NetworKIN score | Integrated Bayesian score |
| Motif probability | Raw pynetphorest posterior |
| STRING score | Raw STRING proximity score |
| Target STRING ID | Ensembl protein ID of target |
| Kinase STRING ID | Ensembl protein ID of enzyme |
| Target Name | Human-readable target name |
| Kinase Name | Human-readable enzyme name |
| Target description | STRING description of target |
| Kinase description | STRING description of enzyme |
| Peptide sequence window | ±7 aa window around phosphosite |
| Intermediate nodes | Best-path intermediate STRING proteins |
| recovered | `True` if recovered by false-negative recovery |
| recovery_method | Method used (e.g. `context_proximity`) |

---

## 4. Module Reference

### `networkin.py`

Main pipeline module (~970 lines). Responsibilities:

- `AppConfig` dataclass — all pipeline configuration.
- `detect_site_file_type` — site file format detection.
- `read_fasta_file`, `read_networkin_sites`, `read_proteome_discoverer_sites`,
  `read_max_quant_sites` — input parsers.
- `read_group_to_domain_map`, `read_alias_files` — reference data loading.
- `map_peptides_to_string` — BLAST-based sequence → STRING ID mapping.
- `load_string_data` — STRING network loading.
- `compile_predictions` — Bayesian score integration.
- `recover_predictions` — false-negative recovery orchestration.
- `run_pipeline` — top-level pipeline entry point.

### `motif_scoring.py`

- `score_sequences(id_seq, id_pos_res)` — pynetphorest batch scorer.
  Returns `result[protein_id][pos][tree][kinase] = (res, peptide, score)`.

### `graph_scoring.py`

- `compute_networkin_score(motif_likelihood, context_likelihood)` — score product.
- `filter_and_rank_predictions(predictions, ...)` — post-scoring filter/rank.

### `likelihood.py`

Utility module for likelihood-ratio calibration. Provides:

| Function | Description |
|---|---|
| `GenerateLikelihoodConversionTbl` | Bin predictions and compute per-bin likelihood ratios |
| `LocalSmooth` / `SmoothUpPeak` / `SmoothDownPeak` | Smooth the conversion table |
| `WriteConversionTableBin` | Persist a conversion table to disk |
| `ReadConversionTableBin` | Load a conversion table from disk |
| `ReadConversionTableFromMemory` | Parse a conversion table from an in-memory string |
| `ConvertScore2L` | Interpolate a likelihood ratio from a loaded table |

### `output.py`

- `STANDARD_COLUMNS` — authoritative column list for TSV output.
- `write_tsv(predictions, path)` — write tab-separated output.
- `write_cytoscape(predictions, path)` — write Cytoscape SIF edge list.
- `write_output(predictions, path, fmt)` — dispatcher.

### `recovery.py`

- `CONTEXT_RECOVERY_THRESHOLD` — default threshold (0.6).
- `recover_false_negatives(candidates, dist_matrix, node_index, motif_scores)` —
  recovers missed kinase–substrate pairs.

### `inputs/phosphosites.py`

- `fetch_phosphosite(refresh)` — OmniPath → local PSP → bundled fallback.

### `inputs/string_network.py`

- `fetch_string_network(proteins, species, min_score, refresh)` — flat file
  → STRING REST API → bundled fallback.

### `inputs/maxquant_processor.py`

- `MaxQuantProcessor` — cleans MaxQuant phosphosite tables, normalises
  protein IDs, downloads FASTA sequences from UniProt.

---

## 5. Key Dependencies

| Dependency | Role | Notes |
|---|---|---|
| Python ≥ 3.10 | Runtime | |
| pynetphorest | Motif scoring | Python package; bundles all kinase models |
| NumPy / pandas | Data manipulation | |
| networkx | Floyd–Warshall shortest paths | Used in false-negative recovery |
| httpx | HTTP requests | OmniPath and STRING REST API |
| NCBI BLAST+ (`blastp`) | Sequence → STRING ID mapping | Install separately; must be on `PATH` or set via `--blast-dir` |
| typer / rich | CLI and progress output | |
| loguru | Logging | |
| pyarrow | Parquet conversion tables | |
