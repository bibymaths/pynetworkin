FROM python:3.11-slim

# Install runtime dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    ca-certificates \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# Install NCBI BLAST+ 2.17.0 from official binaries (cross-platform)
ARG TARGETARCH
RUN set -eux; \
    case "$(uname -m)" in \
        x86_64|amd64) BLAST_ARCH="x64" ;; \
        aarch64|arm64) BLAST_ARCH="aarch64" ;; \
        *) echo "Unsupported architecture: $(uname -m)" && exit 1 ;; \
    esac; \
    BLAST_VERSION="2.17.0+"; \
    BLAST_TARBALL="ncbi-blast-${BLAST_VERSION}-${BLAST_ARCH}-linux.tar.gz"; \
    echo "Downloading BLAST+ ${BLAST_VERSION} for ${BLAST_ARCH}..."; \
    curl -fSL "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/${BLAST_TARBALL}" \
        -o /tmp/blast.tar.gz \
        || { echo "ERROR: Failed to download ${BLAST_TARBALL}. Check architecture (${BLAST_ARCH}) and network access."; exit 1; }; \
    mkdir -p /opt/blast; \
    tar -xzf /tmp/blast.tar.gz -C /opt/blast --strip-components=1 \
        || { echo "ERROR: Failed to extract ${BLAST_TARBALL}. Archive may be corrupt."; exit 1; }; \
    rm /tmp/blast.tar.gz; \
    blastn -version || { echo "ERROR: BLAST+ installation verification failed."; exit 1; }

ENV PATH="/opt/blast/bin:${PATH}"
ENV PYTHONPATH="/app/src"

WORKDIR /app

# Install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Source code is mounted at runtime via docker-compose volume
