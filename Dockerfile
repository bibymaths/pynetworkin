FROM python:3.11-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    ca-certificates \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

RUN set -eux; \
    case "$(uname -m)" in \
        x86_64|amd64) BLAST_ARCH="x64" ;; \
        aarch64|arm64) BLAST_ARCH="aarch64" ;; \
        *) echo "Unsupported architecture: $(uname -m)" && exit 1 ;; \
    esac; \
    BLAST_VERSION="2.17.0+"; \
    BLAST_TARBALL="ncbi-blast-${BLAST_VERSION}-${BLAST_ARCH}-linux.tar.gz"; \
    curl -fSL "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/${BLAST_TARBALL}" \
        -o /tmp/blast.tar.gz; \
    mkdir -p /opt/blast; \
    tar -xzf /tmp/blast.tar.gz -C /opt/blast --strip-components=1; \
    rm /tmp/blast.tar.gz; \
    /opt/blast/bin/blastp -version

ENV PATH="/opt/blast/bin:${PATH}"
ENV PYTHONUNBUFFERED=1

WORKDIR /app

COPY pyproject.toml README.md LICENSE ./
COPY src ./src

COPY data/conversion_direct.parquet ./data/
COPY data/conversion_indirect.parquet ./data/
COPY data/fallback ./data/fallback

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir .

ENTRYPOINT ["pynetworkin"]
CMD ["--help"]