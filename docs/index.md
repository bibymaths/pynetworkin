# pynetworkin

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

See [ARCHITECTURE.md](https://github.com/bibymaths/pynetworkin/blob/main/ARCHITECTURE.md) for full technical details.
