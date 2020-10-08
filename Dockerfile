FROM iontorrent/tsbuild:bionic-5.12 AS builder
# Download and extract TorrentSuite v5.12.1sp1
RUN curl --silent -L https://github.com/iontorrent/TS/archive/TorrentSuite_5.12.1.sp1.tar.gz | tar --strip-components=1 -xvzf -
# Build Analysis modules including tmap and tvc
RUN MODULES=Analysis buildTools/build.sh

FROM nfcore/base:1.10.2
LABEL authors="Peter Kruczkiewicz" \
      version="1.0.0" \
      description="Docker image containing all software requirements for the peterk87/ionampliseq pipeline"


# Install OpenBLAS for tvc
RUN apt update && \
    apt install -y libopenblas-dev libopenblas-base && \
    apt clean

COPY --from=builder /src/build/Analysis/tvc /usr/local/bin
COPY --from=builder /src/build/Analysis/tvcutils /usr/local/bin
COPY --from=builder /src/build/Analysis/tvcassembly /usr/local/bin
COPY --from=builder /src/build/Analysis/TMAP/tmap /usr/local/bin

# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-ionampliseq-1.0dev/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-ionampliseq-1.0dev > nf-core-ionampliseq-1.0dev.yml

# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron
