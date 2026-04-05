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
- Outputs per-site predictions as a CSV file in the `results/` directory.

---

## Requirements

| Dependency | Version | Notes |
|---|---|---|
| Python | ≥ 3.6 | |
| NumPy | any | |
| Pandas | any | (used by `filter_sites.py`) |
| NCBI BLAST+ | ≥ 2.9 | `blastp` must be on `PATH` or supplied via `-b` |
| NetPhorest | — | C binary in `netphorest/`; see compilation instructions below |

---

## Installation

### 1. Compile NetPhorest

```bash
cd netphorest
cc -O3 -o netphorest netphorest.c -lm
cd ..
```

> **Note:** Avoid GCC 4.x — it causes silent crashes in the NetPhorest binary.
> GCC 5+ or Clang are recommended.

### 2. Install BLAST+

Download from [NCBI BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).
Ensure `blastp` and `makeblastdb` are on your `PATH`, or pass the directory with `-b`.

### 3. Install Python dependencies

```bash
pip install numpy pandas
```

---

## Usage

```
python3 NetworKIN.py [options] <organism> <FASTA-file> [sites-file]
```

| Argument | Description |
|---|---|
| `organism` | NCBI taxon ID: `9606` (human) or `4932` (yeast) |
| `FASTA-file` | FASTA file with protein sequences; IDs must match the sites file |
| `sites-file` | Phosphosite file (optional; if omitted, all S/T/Y residues are predicted) |

### Key options

| Flag | Default | Description |
|---|---|---|
| `-n` / `--netphorest` | `$NETPHOREST_PATH` | Path to the NetPhorest binary |
| `-b` / `--blast` | `$BLAST_PATH` | Directory containing BLAST+ binaries |
| `-d` / `--data` | `./data` | Directory with reference data files |
| `-p` / `--path` | `direct` | STRING path type: `direct` or `indirect` |
| `-t` / `--threads` | `1` | Number of BLAST threads |
| `-v` / `--verbose` | off | Print detailed progress to stderr |
| `-u` / `--uncovered` | off | Use STRING likelihood for kinases not covered by NetPhorest |

### Example

```bash
python3 NetworKIN.py \
    -n netphorest/netphorest \
    -d data \
    9606 \
    test.fas \
    test.tsv
```

Results are written to `results/<fasta-filename>.result.csv`.

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
| Rune's format | column 2 = `phospho` | Space-separated with residue+position in col 2 |

---

## Output format

Results CSV columns:

| Column | Description |
|---|---|
| Name | Target protein ID |
| Position | Phosphosite position in the protein |
| Tree | NetPhorest tree (KIN, SH2, PTP, 1433, …) |
| NetPhorest Group | NetPhorest classifier group |
| Kinase/Phosphatase/Phospho-binding domain | Predicted enzyme |
| NetworKIN score | Integrated Bayesian score (≥ 0.02 reported) |
| NetPhorest probability | Raw NetPhorest posterior |
| STRING score | STRING best-path proximity score |
| Target STRING ID | Ensembl protein ID of the substrate |
| Kinase STRING ID | Ensembl protein ID of the enzyme |
| Target Name | Human-readable substrate name |
| Kinase Name | Human-readable enzyme name |
| Target description | STRING functional description of substrate |
| Kinase description | STRING functional description of enzyme |
| Peptide sequence window | ±7 aa window around the phosphosite |
| Intermediate nodes | Best-path intermediate proteins in STRING |

---

## Repository structure

See [ARCHITECTURE.md](ARCHITECTURE.md) for a detailed description of the code
structure and execution flow.

---

## Data sources

- **NetPhorest**: kinase-group motif models.
  Source: `netphorest/`; original publication: Miller *et al.*, 2008.
- **STRING v12**: human protein interactions and sequences.
  Downloaded from [string-db.org](https://string-db.org).
- **HGNC symbol mapping**: `HGNC_Symbol_mapping/`.

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