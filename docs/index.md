# pynetworkin

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19433562.svg)](https://doi.org/10.5281/zenodo.19433562)
[![PyPI version](https://img.shields.io/pypi/v/pynetworkin-bio.svg)](https://pypi.org/project/pynetworkin-bio/)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/pynetworkin-bio.svg)](https://pypi.org/project/pynetworkin-bio/)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/pynetworkin-bio.svg)](https://pypi.org/project/pynetworkin-bio/)
[![Docker Image Version](https://img.shields.io/docker/v/bibymaths/pynetworkin?label=ghcr.io&logo=docker)](https://ghcr.io/bibymaths/pynetworkin)
[![GitHub Release](https://img.shields.io/github/v/release/bibymaths/pynetworkin)](https://github.com/bibymaths/pynetworkin/releases)
[![GitHub Actions](https://img.shields.io/github/actions/workflow/status/bibymaths/pynetworkin/release.yml?label=release)](https://github.com/bibymaths/pynetworkin/actions/workflows/release.yml)
[![License](https://img.shields.io/github/license/bibymaths/pynetworkin)](https://github.com/bibymaths/pynetworkin/blob/main/LICENSE) 

**Kinase–substrate network prediction** — a modernised implementation of the NetworKIN algorithm.

## Overview

pynetworkin predicts kinase–substrate relationships by combining:

- **Motif scoring** via [pynetphorest](https://github.com/bibymaths/pynetworkin) (NetPhoREST atlas)
- **Network context scoring** using the STRING v12.0 protein interaction network
- **False-negative recovery** to rescue pairs missed by motif scoring alone
- **Live data sources**: OmniPath REST API and STRING API with local caching

## Quick Start

```bash
pip install pynetworkin
pynetworkin predict input.fasta --output results.tsv
```

Or with uv:

```bash
uvx pynetworkin predict input.fasta -o results.tsv --refresh
```

## Key Features

| Feature | Description |
|---------|-------------|
| Motif scoring | PyNetPhoREST atlas with >500 kinase models |
| Network context | STRING v12.0 combined scores |
| FN recovery | Context-proximity based recovery |
| Live data | OmniPath REST, STRING API, local caching |
| CLI | Typer + Rich progress bars and summary tables |
| Output | TSV and Cytoscape SIF formats |

## Architecture

See [ARCHITECTURE.md](https://github.com/bibymaths/pynetworkin/blob/main/docs/ARCHITECTURE.md) for full technical details.
