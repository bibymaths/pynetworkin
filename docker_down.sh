#!/bin/bash
# Stop and remove the pynetworkin Docker containers.
# Make executable first: chmod +x docker_down.sh
set -e
docker compose down
echo "Container 'networkin' has been stopped and removed."
