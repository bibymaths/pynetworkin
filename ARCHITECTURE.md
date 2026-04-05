# PyNetworKIN — Architecture Overview

## 1. Project Summary

**PyNetworKIN** (formerly NetworKIN 3.0) is a Bayesian kinase–substrate prediction
pipeline for phosphoproteomics. It integrates two complementary sources of evidence:

1. **NetPhorest** – an ANSI-C binary that scores peptide-sequence motifs using
   per-kinase position-specific scoring matrices (PSSMs), producing a posterior
   probability for each (site, kinase-group) pair.
2. **STRING network** – a pre-computed best-path proximity score between proteins in
   the STRING protein–protein interaction network.

The two evidence streams are combined via pre-calibrated Bayesian likelihood-ratio
conversion tables to produce a single **NetworKIN score** per (phosphosite, kinase/
phosphatase/phospho-binding domain) pair.

Supported organisms: **human (9606)** and **yeast (4932)**.

---

## 2. Repository Map

```
.
├── NetworKIN.py                   # Main entry point (CLI + pipeline orchestration)
├── likelihood.py                  # Likelihood-ratio conversion table utilities
├── filter_sites.py                # Standalone utility: filter MS phosphosites → .tsv
├── batch_networkin.sh             # SLURM batch script for cluster execution
│
├── netphorest/                    # NetPhorest C source + compiled binary
│   ├── netphorest.c               # Main source
│   ├── netphorest                 # Compiled binary (Linux x86-64)
│   ├── motifs.c / pssm_data.h … # Supporting C source/header files
│   └── hugo_pelm.tsv             # Kinase-group → HUGO name mapping
│
├── data/                          # All reference data
│   ├── 9606.protein.sequences.v12.0.fa          # STRING v12 human proteome (BLAST DB)
│   ├── 9606.protein.sequences.v9.0.fa           # STRING v9 human proteome (legacy)
│   ├── 9606.protein.aliases.v12.0.txt.gz        # STRING v12 protein aliases
│   ├── 9606.protein.info.v12.0.txt.gz           # STRING v12 protein descriptions
│   ├── 9606.alias_best.v9.0.tsv.gz              # STRING v9 best aliases (human)
│   ├── 9606.text_best.v9.0.tsv.gz               # STRING v9 protein descriptions (human)
│   ├── 4932.*                                   # Equivalent files for yeast
│   ├── string_data/
│   │   ├── 9606.bestpath_0340_0950.v9.tsv.gz    # Pre-computed best-path STRING scores (human, v9)
│   │   ├── 4932.string_000_0170_1000.tsv.gz     # Pre-computed STRING scores (yeast)
│   │   ├── 9606.links.v12.0.tsv.gz              # Raw STRING v12 links (human)
│   │   └── 9606.protein.links.v12.0.txt.gz      # Full STRING v12 network
│   ├── likelihood_conversion_table_direct/      # Per-tree, per-species conversion tables
│   │   └── conversion_tbl_{netphorest|string}_smooth_{species}_{tree}_{kinase}.txt
│   ├── likelihood_conversion_table_indirect/    # Same, for indirect-path STRING scores
│   ├── group_human_protein_name_map.tsv         # NetPhorest group → protein name (human)
│   ├── hanno_group_human_protein_name_map.tsv   # Updated group → name map (v12 aliases)
│   └── new_name_map.py                          # Script that generated hanno_ map
│
├── binary_human_public/           # Alternative NetPhorest binary (human-only trees)
├── binary_universal_public/       # Alternative NetPhorest binary (all trees)
├── HGNC_Symbol_mapping/           # HGNC symbol ↔ Ensembl mapping files
│
├── results/                       # Output directory (auto-created at runtime)
├── tmp/                           # Temporary working directory (auto-created)
│
├── test.fas / test.tsv            # Minimal test input (SRC protein, 3 sites)
├── test1.fas / test1.tsv          # Secondary test input
├── phospho.tsv                    # Sample phosphosite input
├── cured_morpho_seqs_v2.fa        # Morpho-proteomics FASTA sequences
├── MS_MCMC.csv                    # Mass-spec MCMC sample data
│
├── README.txt                     # Original NetworKIN 3.0 installation notes
├── README_new.txt                 # Updated notes (STRING v12, Python 3 migration)
└── README.md                      # Project README
```

---

## 3. Execution Flow (End-to-End)

### 3.1 Entry Point

```
python3 NetworKIN.py [options] <organism> <FASTA-file> [sites-file]
```

All configuration is parsed by `optparse.OptionParser` at the bottom of
`NetworKIN.py`; the `Main()` function is then called.

### 3.2 Input Ingestion

| Stage | Function | Description |
|---|---|---|
| FASTA | `readFasta(fastafile)` | Reads protein sequences; IDs taken as the first `_`-delimited field |
| Sites | `CheckInputType(sitesfile)` | Auto-detects one of 5 formats (see below) |
| Sites | `readPhosphoSites` / `readPhosphoSitesProteomeDiscoverer` / `readPhosphoSitesMaxQuant` / `readRunessitesfile` / `readMCMCssitesfile` | Format-specific parsers |

Detected site-file formats:

| Constant | Description |
|---|---|
| `NETWORKIN_SITE_FILE` | 3-column TSV: `protein \t position \t residue` |
| `PROTEOME_DISCOVERER_SITE_FILE` | 2-column: `protein \t phosphopeptide` (lowercase = phosphosite) |
| `MAX_QUANT_DIRECT_OUTPUT_FILE` | MaxQuant phosphosite output |
| `RUNES_SITE_FILE` | Space-separated; residue+position in column 2 |
| `MS_MCMC_FILE` | Filename starts with `MS` |

### 3.3 Reference Data Loading

1. **Group → Domain map** (`ReadGroup2DomainMap`): loads
   `data/hanno_group_human_protein_name_map.tsv` which maps NetPhorest
   classifier groups to HGNC-approved protein names.
2. **Alias / description hashes** (`readAliasFiles`): reads STRING v9 text
   descriptions and STRING v12 protein aliases; builds `alias_hash` (STRING ID
   → best name), `desc_hash` (STRING ID → description), and `name_hash`
   (STRING ID → list of aliases).

### 3.4 Sequence → STRING Mapping (BLAST)

`mapPeptides2STRING(...)`:

- Writes query sequences to a temp FASTA.
- Runs `blastp` against the STRING sequence database (`9606.protein.sequences.v12.0.fa`).
- The best BLAST hit (sorted by bit-score) is kept; identity < 90 % or
  e-value > 1e-40 trigger a warning.
- Builds bidirectional dictionaries `incoming2string` and `string2incoming`.

### 3.5 STRING Network Loading

`loadSTRINGdata(string2incoming, ...)`:

- Reads the pre-computed best-path score file
  `data/string_data/9606.bestpath_0340_0950.v9.tsv.gz`.
- Supports 3-, 6-, 7-, or 8-column formats depending on whether indirect-path
  scores and path annotations are included.
- Populates `tree_pred_string_data[substrate_STRING_id][kinase_STRING_id]`
  with `{"_name": name, "_score": float, "_path": path}`.

### 3.6 Motif Scoring (NetPhorest)

`runNetPhorest_one_instance(id_seq, id_pos_res)` / `runNetPhorest(...)`:

- Writes sequences to a temp FASTA file.
- Invokes the compiled NetPhorest binary (path from `-n` flag or
  `NETPHOREST_PATH` env var).
- Parses tab-delimited output:
  `Name \t Position \t Residue \t Peptide \t Method \t Organism \t Tree \t Classifier \t Posterior`
- Builds `id_pos_tree_pred[protein_id][position][tree][group]` = `(residue, peptide, score)`.

### 3.7 Bayesian Score Integration

`printResult(...)`:

- For each (protein, position, tree, group) tuple from NetPhorest:
  1. Look up the STRING mapping for the substrate.
  2. Find the matching STRING edge(s) between the substrate and kinase.
  3. Load the appropriate likelihood-ratio conversion table from
     `data/likelihood_conversion_table_{direct|indirect}/`.
  4. Convert raw NetPhorest and STRING scores to likelihood ratios via
     `ConvertScore2L(score, conv_tbl)` (linear interpolation on the
     pre-calibrated table).
  5. Compute unified likelihood: `L_final = L_netphorest × L_string`.
  6. Emit results where `NetworKIN score ≥ 0.02`.
- Results are written to `results/<fastafile>.result.csv`.

### 3.8 Output Format

CSV columns (header row written at the start):

| Column | Description |
|---|---|
| Name | Target protein ID |
| Position | Phosphosite position |
| Tree | NetPhorest tree (KIN, SH2, PTP, 1433, …) |
| NetPhorest Group | NetPhorest classifier group |
| Kinase/Phosphatase/Phospho-binding domain | Enzyme name |
| NetworKIN score | Integrated Bayesian score |
| NetPhorest probability | Raw NetPhorest posterior |
| STRING score | Raw STRING proximity score |
| Target STRING ID | Ensembl protein ID of target |
| Kinase STRING ID | Ensembl protein ID of enzyme |
| Target Name | Human-readable target name |
| Kinase Name | Human-readable enzyme name |
| Target description | STRING description of target |
| Kinase description | STRING description of enzyme |
| Peptide sequence window | ±7 aa window around phosphosite |
| Intermediate nodes | Best-path intermediate STRING proteins |

---

## 4. Module Reference

### `NetworKIN.py`

The single-file pipeline script (~1478 lines). Responsibilities:
- CLI argument parsing (`optparse`)
- Input format detection and parsing
- BLAST-based protein mapping
- STRING network loading
- NetPhorest execution and output parsing
- Bayesian score integration
- CSV output

Notable globals set at startup: `organism`, `fastafile`, `sitesfile`,
`blastDir`, `netphorest_bin`, `fn_blast_output`, `fn_netphorest_output`.

### `likelihood.py`

Utility module for likelihood-ratio calibration. Provides:

| Function | Description |
|---|---|
| `GenerateLikelihoodConversionTbl` | Bin predictions and compute per-bin likelihood ratios |
| `LocalSmooth` / `SmoothUpPeak` / `SmoothDownPeak` | Smooth the conversion table (ensure monotonicity) |
| `WriteConversionTableBin` | Persist a conversion table to disk |
| `ReadConversionTableBin` | Load a conversion table from disk |
| `ConvertScore2L` | Interpolate a likelihood ratio from a loaded table |

### `filter_sites.py`

Standalone script (not imported by `NetworKIN.py`). Reads a mass-spec CSV,
filters to valid S/T/Y phosphosites, and outputs a `.tsv` suitable for use
as a `SitesFile`.

---

## 5. Key Dependencies

| Dependency | Role | Notes |
|---|---|---|
| Python ≥ 3.6 | Runtime | Migrated from Python 2.7 |
| NumPy | `np.unique` in result counting | |
| Pandas | Used in `filter_sites.py` | |
| NCBI BLAST+ (`blastp`) | Sequence → STRING ID mapping | Install separately; set `BLAST_PATH` env var or use `-b` |
| NetPhorest binary | Motif scoring | C source in `netphorest/`; compile with `cc -O3 -o netphorest netphorest.c -lm` |

---

## 6. Known Issues / Technical Debt

1. **`readMCMCssitesfile` bug** – `res_pos` is referenced before assignment in
   the filter condition; the condition is also logically always `True`.
2. **`from string import *`** – broad wildcard import in `NetworKIN.py`; safe but
   not idiomatic Python 3.
3. **STRING v12 coverage gap** – fewer STRING-backed predictions are obtained
   with STRING v12 than v9 (see `README_new.txt`); the v9 best-path file is
   still used as default.
