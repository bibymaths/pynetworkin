#!/usr/bin/env bash
# Local test harness for pynetworkin Docker setup.
# Tests build, compose-dev (build), compose (image pull), interactive run,
# and down — WITHOUT pushing or pulling from Docker Hub.
# Usage: ./docker_test_local.sh [version_tag]
set -euo pipefail

IMAGE_NAME="ghcr.io/bibymaths/pynetworkin"
VERSION="${1:-local-test}"
FULL_TAG="${IMAGE_NAME}:${VERSION}"

PASS="\033[0;32m✔\033[0m"
PASS="\033[0;32m✔\033[0m"
FAIL="\033[0;31m✘\033[0m"
INFO="\033[0;34m➜\033[0m"

header() { echo -e "\n\033[1;33m=== $* ===\033[0m"; }
pass()   { echo -e "  ${PASS} $*"; }
fail()   { echo -e "  ${FAIL} $*"; exit 1; }
info()   { echo -e "  ${INFO} $*"; }

# ── 1. Build image locally ──────────────────────────────────────────────────
header "1. Building image: ${FULL_TAG}"
if docker build -t "${FULL_TAG}" -t "${IMAGE_NAME}:latest-local" .; then
    pass "Image built: ${FULL_TAG}"
else
    fail "docker build failed"
fi

# ── 2. Smoke-test: entrypoint --help ────────────────────────────────────────
header "2. Smoke test — pynetworkin --help"
if docker run --rm "${FULL_TAG}" --help; then
    pass "ENTRYPOINT responds to --help"
else
    fail "--help returned non-zero"
fi

# ── 3. BLAST+ sanity check inside image ─────────────────────────────────────
header "3. BLAST+ binary check"
if docker run --rm --entrypoint blastp "${FULL_TAG}" -version; then
    pass "blastp is on PATH and working"
else
    fail "blastp not found or errored"
fi

# ── 4. docker-compose.dev.yml — build + up + down ───────────────────────────
header "4. docker-compose.dev.yml (build mode)"
info "Starting dev compose (tail -f /dev/null keeps it alive)..."
docker compose -f docker-compose.dev.yml up -d --build
sleep 3

info "Checking container is running..."
RUNNING=$(docker compose -f docker-compose.dev.yml ps --services --filter status=running)
if echo "${RUNNING}" | grep -q "networkin"; then
    pass "Service 'networkin' is running in dev compose"
else
    fail "Service 'networkin' not running"
fi

info "Running pynetworkin --help inside dev container..."
if docker compose -f docker-compose.dev.yml exec networkin pynetworkin --help; then
    pass "pynetworkin reachable inside dev container"
else
    fail "pynetworkin not reachable inside dev container"
fi

info "Tearing down dev compose..."
docker compose -f docker-compose.dev.yml down
pass "Dev compose down cleanly"

# ── 5. docker-compose.yml — using local image (no Hub pull) ─────────────────
header "5. docker-compose.yml (image mode, local tag)"
# Temporarily patch compose to use our local-test tag so no Hub pull happens
TMP_COMPOSE=$(mktemp /tmp/compose-test-XXXX.yml)
sed "s|bibymaths/pynetworkin:latest|${IMAGE_NAME}:latest-local|g" docker-compose.yml > "${TMP_COMPOSE}"

info "Starting compose with local image..."
docker compose -f "${TMP_COMPOSE}" up -d
sleep 3

RUNNING2=$(docker compose -f "${TMP_COMPOSE}" ps --services --filter status=running)
if echo "${RUNNING2}" | grep -q "networkin"; then
    pass "Service 'networkin' is running in image-mode compose"
else
    fail "Service 'networkin' not running in image-mode compose"
fi

docker compose -f "${TMP_COMPOSE}" down
rm -f "${TMP_COMPOSE}"
pass "Image-mode compose down cleanly"

# ── 6. Interactive shell test ────────────────────────────────────────────────
header "6. Interactive shell sanity (non-TTY)"
info "Running python --version inside container..."
if docker run --rm --entrypoint python "${FULL_TAG}" --version; then
    pass "Python available in image"
else
    fail "python not found"
fi

info "Checking /work mount point exists..."
if docker run --rm -v "$(pwd):/work" --entrypoint ls "${FULL_TAG}" /work > /dev/null 2>&1; then
    pass "/work volume mount works"
else
    fail "/work mount check failed"
fi

# ── 7. docker_down.sh script ────────────────────────────────────────────────
header "7. docker_down.sh (idempotent teardown)"
info "Spinning up a compose stack to test down script..."
docker compose -f docker-compose.dev.yml up -d --build > /dev/null 2>&1
sleep 2
if bash docker_down.sh; then
    pass "docker_down.sh ran successfully"
else
    fail "docker_down.sh failed"
fi

# ── Summary ──────────────────────────────────────────────────────────────────
header "All local tests passed 🎉"
echo -e "  Image built and tagged as: \033[1m${FULL_TAG}\033[0m"
echo -e "  Use \033[1m./docker_push.sh <version>\033[0m only when triggered by a git tag."