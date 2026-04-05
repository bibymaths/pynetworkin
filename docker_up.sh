#!/bin/bash
# Build and start the pynetworkin Docker environment in detached mode.
# Make executable first: chmod +x docker_up.sh
set -e
docker compose up --build -d
echo "Container 'networkin' is running in the background."
