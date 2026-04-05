#!/bin/bash
# Open an interactive bash shell inside the running 'networkin' container.
# Make executable first: chmod +x docker_interactive.sh
# The source is mounted at /app and /app/src is on PYTHONPATH automatically.
# To register the pynetworkin CLI entry-point inside the container, run:
#   pip install -e .
# Usage example inside the container:
#   pynetworkin predict data/fallback/test_input.fasta --output /tmp/out.tsv
#   pynetworkin predict proteome.fasta --sites peptides.txt --output results.tsv
set -e
docker compose exec networkin bash
