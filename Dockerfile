# Build stage for samtools
FROM ubuntu:plucky-20250714 AS samtools-build
ARG SAMTOOLS_VERSION=1.22.1
ARG MAKE_JOBS=4

# Install build dependencies (cached as long as this line and base image don't change)
RUN --mount=type=cache,target=/var/cache/apt \
    apt update && \
    apt install -y \
        build-essential zlib1g-dev libncurses5-dev \
        libbz2-dev liblzma-dev libcurl4-openssl-dev wget && \
    apt clean && \
    rm -rf /var/lib/apt/lists/*

# Only this part will rebuild if SAMTOOLS_VERSION changes
WORKDIR /tmp
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && \
    make -j${MAKE_JOBS} && \
    make install && \
    cd .. && \
    rm -rf samtools-${SAMTOOLS_VERSION}*


FROM ubuntu:plucky-20250714
LABEL authors="Peter Kruczkiewicz" \
      version="2.0.0" \
      description="Docker image for TMAP and TVC the CFIA-NCFAD/nf-ionampliseq pipeline" \
      org.opencontainers.image.source="https://github.com/CFIA-NCFAD/nf-ionampliseq"

# Copy tmap, tvc and other related binaries from peterk87/tvc-tmap-torrent-suite:v1.0.0
COPY --from=peterk87/tvc-tmap-torrent-suite:v1.0.0 /usr/local/bin/tvc /usr/local/bin
COPY --from=peterk87/tvc-tmap-torrent-suite:v1.0.0 /usr/local/bin/tvcutils /usr/local/bin
COPY --from=peterk87/tvc-tmap-torrent-suite:v1.0.0 /usr/local/bin/tvcassembly /usr/local/bin
COPY --from=peterk87/tvc-tmap-torrent-suite:v1.0.0 /usr/local/bin/tmap /usr/local/bin

COPY --from=samtools-build /usr/local/bin/samtools /usr/local/bin/samtools

# Install runtime dependencies for samtools and OpenBLAS for tvc
RUN --mount=type=cache,target=/var/cache/apt \
    apt update && \
    apt install -y \
        libcurl4t64 \
        libncurses6 \
        libopenblas-dev && \
    apt clean && \
    rm -rf /var/lib/apt/lists/*

