#!/usr/bin/env bash
# Stop and remove the pynetworkin Docker services.
set -euo pipefail

docker compose down
echo "PyNetworKIN containers stopped and removed."
