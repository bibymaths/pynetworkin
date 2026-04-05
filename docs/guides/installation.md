# Installation

## PyPI (recommended)

```bash
pip install pynetworkin
```

## Source install with uv

```bash
git clone https://github.com/bibymaths/pynetworkin
cd pynetworkin
uv sync
uv run pynetworkin --help
```

## Development install

```bash
uv sync --group dev
```

## Docs dependencies

```bash
uv sync --group docs
uv run --group docs mkdocs serve
```

## Verify installation

```bash
pynetworkin info
```

See [RUN.md](../RUN.md) for full run instructions including database setup, data downloads, and example commands.
