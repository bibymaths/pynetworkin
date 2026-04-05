# Codebase Audit

## File Discovery

### Python files (`find . -type f -name "*.py" | sort`)
```
./HGNC_Symbol_mapping/cleanup_HGNC_mapping.py
./conftest.py
./data/new_name_map.py
./data/string_data/string_network_filter.py
./scripts/backup.py
./scripts/cleanup_HGNC_mapping.py
./scripts/generate_sample_data.py
./scripts/migrate_to_parquet.py
./src/pynetworkin/__init__.py
./src/pynetworkin/cli.py
./src/pynetworkin/graph_scoring.py
./src/pynetworkin/inputs/__init__.py
./src/pynetworkin/inputs/maxquant_processor.py
./src/pynetworkin/inputs/phosphosites.py
./src/pynetworkin/inputs/string_network.py
./src/pynetworkin/likelihood.py
./src/pynetworkin/logger.py
./src/pynetworkin/motif_scoring.py
./src/pynetworkin/networkin.py
./src/pynetworkin/output.py
./src/pynetworkin/recovery.py
./tests/test_maxquant_processor.py
./tests/test_motif_scoring.py
./tests/test_networkin.py
./tests/test_output.py
./tests/test_recovery.py
```

---

## Python File Audits

## src/pynetworkin/networkin.py
- **Role**: Main pipeline orchestration module. Defines `AppConfig`, `run_pipeline`,
  `detect_site_file_type`, input parsers, BLAST mapping, STRING network loading,
  Bayesian score integration, and false-negative recovery orchestration.
- **Entry point**: `run_pipeline(config: AppConfig) -> dict`

## src/pynetworkin/cli.py
- **Role**: Typer CLI entry-point. Defines `pynetworkin predict`, `pynetworkin info`,
  `pynetworkin cache`, and `pynetworkin prepare-maxquant` commands.

## src/pynetworkin/motif_scoring.py
- **Role**: pynetphorest batch scorer wrapper. `score_sequences(id_seq, id_pos_res)`
  returns `{protein_id: {pos: {tree: {kinase: (res, peptide, score)}}}}`.

## src/pynetworkin/graph_scoring.py
- **Role**: STRING context scoring utilities. `compute_networkin_score` multiplies
  motif and context likelihoods. `filter_and_rank_predictions` applies thresholds
  and keeps top-k per site.

## src/pynetworkin/likelihood.py
- **Role**: Likelihood-ratio calibration utilities. Provides `ReadConversionTableBin`,
  `ReadConversionTableFromMemory`, `ConvertScore2L`, and table-generation helpers.

## src/pynetworkin/output.py
- **Role**: Output writers. `write_tsv`, `write_cytoscape`, `write_output`.
  `STANDARD_COLUMNS` defines the canonical TSV column order.

## src/pynetworkin/recovery.py
- **Role**: False-negative recovery. `recover_false_negatives` uses Floyd-Warshall
  shortest-path distances to rescue kinase-substrate pairs missed by motif scoring.

## src/pynetworkin/logger.py
- **Role**: Loguru/Rich logging wrapper.

## src/pynetworkin/inputs/phosphosites.py
- **Role**: Phosphosite data fetcher. Priority: OmniPath -> local PSP -> bundled fallback.
  Caches to `.cache/phosphosite.parquet` (7-day TTL).

## src/pynetworkin/inputs/string_network.py
- **Role**: STRING network fetcher. Priority: flat file -> REST API -> bundled fallback.
  Caches to `.cache/string_<species>_<min_score>.parquet` (7-day TTL).

## src/pynetworkin/inputs/maxquant_processor.py
- **Role**: MaxQuant site-table pre-processor. Cleans protein IDs, downloads FASTA
  sequences from UniProt, writes `cleaned_proteins.fasta`, `cleaned_sites.txt`,
  `id_mapping.csv`, `processing_report.json`.

## scripts/backup.py
- **Role**: Legacy Python 3 port of the original NetworKIN 3.0 script. Not executed
  by the package; kept as historical reference.

## scripts/cleanup_HGNC_mapping.py
- **Role**: HGNC symbol-Ensembl ID reconciliation utility. Requires the
  `Lindinglab.idmapper` library (external/not public). Run manually, not imported
  by the package.

## scripts/generate_sample_data.py
- **Role**: Generates `data/fallback/` sample files from live sources. Run manually
  when updating the bundled offline test data.

## scripts/migrate_to_parquet.py
- **Role**: Migrates per-kinase `.txt` conversion tables to the consolidated
  `data/conversion_direct.parquet` and `data/conversion_indirect.parquet` formats.

## data/new_name_map.py
- **Role**: One-off data-preparation script that generated
  `data/group_human_protein_name_map_curated.tsv` from STRING v12 aliases.
  Not imported by the package; run manually when refreshing the name map.

## data/string_data/string_network_filter.py
- **Role**: One-off data-preparation script that filtered the STRING v12 links file
  to retain only interactions involving known kinases. Not imported by the package.

## HGNC_Symbol_mapping/cleanup_HGNC_mapping.py
- **Role**: Legacy Python 2 script (`.has_key()`, `.iteritems()`) for HGNC symbol
  reconciliation. Superseded by `scripts/cleanup_HGNC_mapping.py`. Requires
  `Lindinglab.idmapper` (external, not public).

---

## Cross-Cutting Observations

| Concern | Files affected |
|---|---|
| No SQL/database calls in core pipeline | -- |
| Person-name identifiers removed | `networkin.py` (RUNES_SITE_FILE -> LEGACY_SITE_FILE, hanno_path -> curated_path); `backup.py`; `new_name_map.py`; `data/hanno_group_human_protein_name_map.tsv` renamed |
| Legacy Python 2 syntax | `HGNC_Symbol_mapping/cleanup_HGNC_mapping.py` (`.has_key()`, `.iteritems()`, `print` statements) |
| External proprietary dependency | `HGNC_Symbol_mapping/cleanup_HGNC_mapping.py`, `scripts/cleanup_HGNC_mapping.py` -- `Lindinglab.idmapper` |
| pynetphorest used via Python API | `motif_scoring.py` -- no C binary subprocess needed |
