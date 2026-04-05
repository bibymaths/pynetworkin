#!/usr/bin/env bash
# Build and push PyNetworKIN image to Docker Hub.
# Usage:
#   ./docker_push.sh 0.1.0
set -euo pipefail

IMAGE_NAME="bibymaths/pynetworkin"
VERSION="${1:-}"

if [[ -z "${VERSION}" ]]; then
  echo "Usage: $0 <version>"
  exit 1
fi

docker build -t "${IMAGE_NAME}:${VERSION}" -t "${IMAGE_NAME}:latest" .
docker push "${IMAGE_NAME}:${VERSION}"
docker push "${IMAGE_NAME}:latest"

echo "Pushed:"
echo "  ${IMAGE_NAME}:${VERSION}"
echo "  ${IMAGE_NAME}:latest"