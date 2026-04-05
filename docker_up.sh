#!/usr/bin/env bash
# Build and start the pynetworkin Docker service in detached mode.
set -euo pipefail

docker compose up --build -d
echo "PyNetworKIN service is running in the background."