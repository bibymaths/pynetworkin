#!/usr/bin/env bash
# Build and push PyNetworKIN image to GHCR.
# Usage: ./docker_push.sh <version>
# Example: ./docker_push.sh 0.2.0
#
# NOTE: Normally handled automatically by the release workflow on git tag.
# Use this only for manual/emergency pushes outside of CI.
set -euo pipefail

IMAGE_NAME="ghcr.io/bibymaths/pynetworkin"
VERSION="${1:-}"

if [[ -z "${VERSION}" ]]; then
    echo "Usage: $0 <version>"
    echo "Example: $0 0.2.0"
    exit 1
fi

# Ensure logged in to GHCR
if ! docker info 2>/dev/null | grep -q "ghcr.io"; then
    echo "Logging in to GHCR..."
    echo "${CR_PAT:?Set CR_PAT env var to a GitHub PAT with write:packages scope}" \
        | docker login ghcr.io -u bibymaths --password-stdin
fi

docker build \
    -t "${IMAGE_NAME}:${VERSION}" \
    -t "${IMAGE_NAME}:latest" \
    .

docker push "${IMAGE_NAME}:${VERSION}"
docker push "${IMAGE_NAME}:latest"

echo ""
echo "Pushed:"
echo "  ${IMAGE_NAME}:${VERSION}"
echo "  ${IMAGE_NAME}:latest"
echo ""
echo "NOTE: The release workflow applies additional tags (major.minor)."
echo "      Prefer tagging via git: git tag v${VERSION} && git push --tags"