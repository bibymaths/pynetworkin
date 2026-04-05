#!/bin/bash
# Open an interactive bash shell inside the running 'networkin' container.
# Make executable first: chmod +x docker_interactive.sh
# Usage example inside the container:
#   python -m src.pynetworkin.networkin --help
set -e
docker compose exec networkin bash
