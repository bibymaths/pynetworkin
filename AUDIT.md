# Codebase Audit

## File Discovery

### Python files (`find . -type f -name "*.py" | sort`)
```
./HGNC_Symbol_mapping/cleanup_HGNC_mapping.py
./NetworKIN.py
./data/new_name_map.py
./data/string_data/string_network_filter.py
./filter_sites.py
./likelihood.py
./results/csv_column_rewrite.py
./results/test.py
```

### C files (`find . -type f -name "*.c" | sort`)
```
./binary_human_public/motifs.c
./binary_human_public/netphorest.c
./binary_universal_public/motifs.c
./netphorest/motifs.c
./netphorest/netphorest.c
./netphorest/netphorest_unsupported_aa.c
```

### Perl files (`find . -type f -name "*.pl" | sort`)
```
./binary_human_public/insr2h.pl
./binary_human_public/nn2h.pl
./binary_human_public/pssm2h.pl
./netphorest/insr2h.pl
./netphorest/nn2h.pl
./netphorest/pssm2h.pl
```

### Shell scripts (`find . -type f -name "*.sh" | sort`)
```
./batch_networkin.sh
./binary_human_public/run.sh
./binary_universal_public/run.sh
./results/batch_NK_results.sh
```

### SQL files
*(none found)*

### Config files (`.cfg` / `.ini` / `.conf`)
*(none found)*

### Top Python files by line count (`wc -l … | sort -rn | head -40`)
```
2359 total
1477 ./NetworKIN.py
 367 ./likelihood.py
 229 ./data/new_name_map.py
 132 ./HGNC_Symbol_mapping/cleanup_HGNC_mapping.py
  72 ./data/string_data/string_network_filter.py
  43 ./results/csv_column_rewrite.py
  32 ./filter_sites.py
   7 ./results/test.py
```

*Only 8 Python files exist in the repository; all are audited below.*

---

## Python File Audits

## NetworKIN.py
- **Role**: Main orchestration script and entry point for the NetworKIN kinase–substrate prediction pipeline. Accepts a FASTA file of protein sequences and an optional phospho-sites file, runs NetPhoREST (via subprocess) to score kinase motifs, maps peptides to STRING interaction partners using BLAST, combines motif scores with STRING co-expression/interaction likelihoods into a unified NetworKIN score, and writes the final predictions to stdout/CSV.
- **MySQL calls**: None. No database driver is imported or used.
- **Subprocess calls**:
  - Line 36: `import subprocess`
  - Lines 106–141: `myPopen()` (defined twice; second definition supersedes the first) — wraps `subprocess.Popen` to run shell commands and capture stdout/stderr.
  - Line 426: `myPopen('gzip -cd …alias_best.tsv.gz')`
  - Line 435: `myPopen('gzip -cd …text_best.tsv.gz')`
  - Line 460: `myPopen('gzip -cd …text_best.v9.0.tsv.gz')`
  - Line 472: `myPopen('gzip -cd …protein.aliases.v12.0.txt.gz')`
  - Line 523: `myPopen(netphorest_bin + ' < ' + file_in.name + ' > ' + file_out.name)` — runs the NetPhoREST binary.
  - Line 654: `threading.Thread(target=myPopen, args=arg).start()` — multi-threaded NetPhoREST invocations.
  - Line 679: `netphorest_results = myPopen(command)` — single-instance NetPhoREST call.
  - Line 766: `myPopen(command)` — BLAST `makeblastdb`.
  - Line 775: `blast_out = myPopen(command)` — BLAST `blastp`.
  - Line 866: `mappingFile = myPopen(command)` — BLAST mapping step.
  - Line 908: `data = myPopen(command)` — additional BLAST query.
- **ctypes**: None.
- **NetPhosK / NetPhoREST / C extension imports**: No Python-level `import` of a C extension. NetPhoREST is invoked as an external compiled binary via `subprocess` (path resolved from the `NETPHOREST_PATH` environment variable or the `--netphorest` CLI option).
- **Exports called from other files**:
  - This is the main script; it is not imported by other files in the repository.
  - It imports `ReadConversionTableBin` and `ConvertScore2L` from `likelihood.py` (lines 43–44).

---

## likelihood.py
- **Role**: Statistical likelihood / score-conversion library. Provides functions to build a score-to-likelihood conversion table from a set of positive/negative predictions (binning, smoothing up-peaks and down-peaks), serialise/deserialise that table to/from disk, and look up the likelihood value for a given raw NetPhoREST probability score.
- **MySQL calls**: None.
- **Subprocess calls**: None.
- **ctypes**: None.
- **NetPhosK / NetPhoREST / C extension imports**: None.
- **Exports called from other files**:
  - `ReadConversionTableBin(path_conversion_table)` — imported and used by `NetworKIN.py` (line 43, used extensively in `printResult`).
  - `ConvertScore2L(score, conv_tbl)` — imported and used by `NetworKIN.py` (line 44, called at lines 1239 and similar).
  - Additional functions (`GenerateLikelihoodConversionTbl`, `WriteConversionTableBin`, `WriteConversionTableFDR`, `LocalSmooth`, helper smoothing functions) are defined but not referenced by any other file in this repository; they appear to be used by an external training/calibration workflow.

---

## data/new_name_map.py
- **Role**: Standalone utility / data-preparation script for mapping STRING protein aliases to kinase group/domain names. Reads the STRING v12.0 protein alias file (gzip-compressed), builds a `name_hash` dictionary, then reads a `group_human_protein_name_map.tsv` file and attempts multi-pass fuzzy name normalisation to produce a curated `hanno_group_human_protein_name_map.tsv`. Executed as a script at import time (module-level code at lines 229–230 calls its own functions immediately).
- **MySQL calls**: None.
- **Subprocess calls**:
  - Line 1: `import subprocess`
  - Lines 6–35: `myPopen()` — wraps `subprocess.Popen` identically to the copy in `NetworKIN.py`.
  - Line 14: `myPopen('gzip -cd …protein.aliases.v12.0.txt.gz')` — called from `readAliasFiles`.
  - Line 41: same `myPopen` gzip call inside `readAliasFiles`.
- **ctypes**: None.
- **NetPhosK / NetPhoREST / C extension imports**: None.
- **Exports called from other files**: No file in this repository imports from `data/new_name_map.py`. The module is self-contained and run directly.

---

## HGNC_Symbol_mapping/cleanup_HGNC_mapping.py
- **Role**: One-off data-cleanup script to produce a revised `9606.alias_best.v9.0.17032014_2.tsv` HGNC symbol mapping file. Uses the `Lindinglab.idmapper` library (external/proprietary) to merge Ensembl–Biomart and HGNC ID maps across multiple Ensembl releases, resolves ambiguous mappings against a manual override dictionary, then integrates results with the existing alias file.
- **MySQL calls**: None.
- **Subprocess calls**: None.
- **ctypes**: None.
- **NetPhosK / NetPhoREST / C extension imports**: Line 1: `from Lindinglab.idmapper import *` — imports from the `Lindinglab` package (not present in this repository; an external/private dependency).
- **Exports called from other files**: None — this is a standalone script not imported elsewhere.
- **Notes**: Uses Python 2-style `print` statements (lines 97, 103, 116, 127) and `.has_key()` / `.iteritems()` dict methods (lines 80, 82, 86, 91, 100, 123, 126) — **Python 2 only**; will fail under Python 3.

---

## data/string_data/string_network_filter.py
- **Role**: One-off data-preparation script that filters the STRING v12.0 protein interaction links file to retain only interactions where one of the proteins is a known kinase, attaches kinase group/domain annotations, and writes the result to a gzip-compressed TSV (`9606.links.v12.0.tsv.gz`).
- **MySQL calls**: None.
- **Subprocess calls**: None.
- **ctypes**: None.
- **NetPhosK / NetPhoREST / C extension imports**: None.
- **Exports called from other files**: None — standalone script not imported elsewhere.

---

## results/csv_column_rewrite.py
- **Role**: One-off result-post-processing script. Reads an existing NetworKIN output CSV (`cured_morpho_seqs_v2.fa.result_old.csv`), selects and renames a specific subset of columns to match the KinomeXplorer schema, filters rows with `networkin_score >= 0.00001`, adds an `Iteration` column, and writes a new CSV (`KinomeXplorer_all_predictions_v3.csv`).
- **MySQL calls**: None.
- **Subprocess calls**: None.
- **ctypes**: None.
- **NetPhosK / NetPhoREST / C extension imports**: None.
- **Exports called from other files**: None — standalone script not imported elsewhere.

---

## filter_sites.py
- **Role**: One-off site-filtering script. Reads phospho-site data from a CSV file (`MS_Gaussian_updated_09032023.csv`) and a FASTA file (`cured_morpho_seqs_v2.fa`), cross-references the two datasets, and writes a TSV (`MS_Gaussian_updated_09032023.tsv`) containing only rows where the modified residue is S, T, or Y and the protein is present in the FASTA.
- **MySQL calls**: None.
- **Subprocess calls**: None.
- **ctypes**: None.
- **NetPhosK / NetPhoREST / C extension imports**: None.
- **Exports called from other files**: None — standalone script not imported elsewhere.
- **Notes**: Uses the deprecated `'rU'` file-open mode (line 14), which raises a `ValueError` in Python 3.11+.

---

## results/test.py
- **Role**: Ad-hoc query/inspection script. Loads the NetworKIN output CSV and prints a filtered view for a specific kinase–substrate pair (MAP2K2 / MAPK3).
- **MySQL calls**: None.
- **Subprocess calls**: None.
- **ctypes**: None.
- **NetPhosK / NetPhoREST / C extension imports**: None.
- **Exports called from other files**: None — standalone throwaway script.

---

## Non-Python Files (C, Perl, Shell)

### C files

#### netphorest/netphorest.c and binary_human_public/netphorest.c
- **Role**: Core NetPhoREST prediction engine implementing a random-forest classifier for kinase recognition motifs. Reads peptide sequences from stdin and writes tab-separated predictions (protein ID, position, tree, group, probability) to stdout. Invoked exclusively via subprocess from `NetworKIN.py`.
- **Subprocess / ctypes in Python**: n/a (is itself the subprocess target).

#### netphorest/motifs.c, binary_human_public/motifs.c, binary_universal_public/motifs.c
- **Role**: Sequence motif scoring helper compiled into the `netphorest` binary. Implements position-specific scoring matrix (PSSM) and other motif-scoring functions used by `netphorest.c`.

#### netphorest/netphorest_unsupported_aa.c
- **Role**: Variant of `netphorest.c` that handles non-standard amino acid codes gracefully (used for organisms/sequences where unusual residues may appear).

### Perl scripts

#### netphorest/insr2h.pl, binary_human_public/insr2h.pl
- **Role**: Converts NetPhoREST insulin receptor PSSM data to a binary format consumed by the C motif-scoring code.

#### netphorest/nn2h.pl, binary_human_public/nn2h.pl
- **Role**: Converts neural-network weight files to binary format for the C scoring code.

#### netphorest/pssm2h.pl, binary_human_public/pssm2h.pl
- **Role**: Converts PSSM text files to binary format for the C scoring code.

### Shell scripts

#### batch_networkin.sh
- **Role**: Batch wrapper that iterates over a list of FASTA/sites file pairs and calls `NetworKIN.py` for each, suitable for running large-scale prediction jobs.

#### binary_human_public/run.sh and binary_universal_public/run.sh
- **Role**: Build/compile scripts that invoke the Perl formatters to regenerate binary model files and then compile the C NetPhoREST binary.

#### results/batch_NK_results.sh
- **Role**: Post-processing batch script that collects multiple NetworKIN result files and merges/filters them.

---

## Cross-Cutting Observations

| Concern | Files affected |
|---|---|
| No MySQL / relational DB calls anywhere | — |
| `subprocess.Popen` (shell=True) | `NetworKIN.py` lines 108, 127; `data/new_name_map.py` line 14 |
| `myPopen` duplicated verbatim | `NetworKIN.py` (two definitions, lines 106 and 119) and `data/new_name_map.py` line 6 |
| Python 2-only syntax | `HGNC_Symbol_mapping/cleanup_HGNC_mapping.py` (`.has_key()`, `.iteritems()`, `print` statements) |
| Deprecated `'rU'` open mode (Python ≥ 3.11 raises ValueError) | `filter_sites.py` line 14 |
| External proprietary dependency | `HGNC_Symbol_mapping/cleanup_HGNC_mapping.py` — `Lindinglab.idmapper` |
| NetPhoREST invoked as compiled binary subprocess | `NetworKIN.py` lines 523, 654, 679 |
| No ctypes usage anywhere | — |
| No SQL files | — |
| No config/ini/conf files | — |

---

## Live Data Migration

| Stale data source | Original location | External DB | Pipeline stage | Status |
|---|---|---|---|---|
| `%s/%s.bestpath_0340_0950.v9.tsv.gz` (STRING network) | `NetworKIN.py` `loadSTRINGdata()` | STRING v9/v12 | Context scoring | [x] live data migration complete |
| `phospho.tsv` / user-provided phospho flat files | `NetworKIN.py` `Main()` lines 1076–1090 | PhosphoSitePlus, Phospho.ELM | Input loading | [x] live data migration complete |
| `%s/%s.text_best.v9.0.tsv.gz` (STRING aliases) | `NetworKIN.py` `readAliasFiles()` | STRING v9 | Input loading | static file (alias lookup; not replaced) |
| `%s/%s.protein.aliases.v12.0.txt.gz` (STRING aliases v12) | `NetworKIN.py` `readAliasFiles()` | STRING v12 | Input loading | static file (alias lookup; not replaced) |
| `9606.protein.sequences.v12.0.fa` (BLAST DB) | `NetworKIN.py` `mapPeptides2STRING()` | STRING v12 | BLAST mapping | static file (BLAST DB; not replaced) |

### New modules

- `inputs/phosphosites.py` — live fetch from PhosphoSitePlus and Phospho.ELM; 7-day Parquet cache in `.cache/`
- `inputs/string_network.py` — live fetch via STRING REST API; 7-day Parquet cache in `.cache/`
- `requirements.txt` — added `httpx`, `pandas`, `pyarrow`

### New CLI flag

`--refresh` added to `NetworKIN.py`: forces re-download of live data even when a valid 7-day cache exists.
