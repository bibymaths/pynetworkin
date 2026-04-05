# Data Sources

pynetworkin uses the following data sources for phosphosite and network data.

## OmniPath (Primary)

[OmniPath](https://omnipathdb.org) provides enzyme–substrate interactions via a public REST API with no authentication required.

- Endpoint: `https://omnipathdb.org/enzsub`
- Species: Homo sapiens (9606)
- No registration required
- Cached locally for 7 days (`.cache/phosphosite.parquet`)

## PhosphoSitePlus (Optional Local)

[PhosphoSitePlus](https://www.phosphosite.org) requires a registered account. pynetworkin supports loading a **locally downloaded** export.

```bash
export NETWORKIN_PSP_LOCAL_FILE=/path/to/Phosphorylation_site_dataset.gz
pynetworkin predict input.fasta
```

> **Note:** Phospho.ELM is no longer supported — the server is offline/unmaintained.

## STRING v12.0

[STRING](https://string-db.org) provides functional protein association networks.

- Flat file (preferred): `data/string_data/9606.links.v12.0.tsv.gz`
  - Override: `NETWORKIN_STRING_FLAT_FILE=/path/to/file`
- REST API fallback: `https://string-db.org/api/tsv/network`
- Cached locally for 7 days

## pynetphorest Atlas (Bundled)

The NetPhoREST kinase motif atlas is bundled via the `pynetphorest` package. Over 500 kinase classifiers are loaded at import time.

## Bundled Fallback Data

Offline fallback TSVs are provided in `data/fallback/` for testing without network access:

- `phosphosites_sample.tsv`
- `string_sample.tsv`
