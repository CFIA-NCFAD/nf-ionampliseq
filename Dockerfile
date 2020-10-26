FROM nfcore/base:1.10.2
LABEL authors="Peter Kruczkiewicz" \
      version="1.0.0" \
      description="Docker image containing all software requirements for the peterk87/ionampliseq pipeline"

# Copy tmap, tvc and other related binaries from peterk87/tvc-tmap-torrent-suite:v1.0.0
COPY --from=peterk87/tvc-tmap-torrent-suite:v1.0.0 /usr/local/bin/tvc /usr/local/bin
COPY --from=peterk87/tvc-tmap-torrent-suite:v1.0.0 /usr/local/bin/tvcutils /usr/local/bin
COPY --from=peterk87/tvc-tmap-torrent-suite:v1.0.0 /usr/local/bin/tvcassembly /usr/local/bin
COPY --from=peterk87/tvc-tmap-torrent-suite:v1.0.0 /usr/local/bin/tmap /usr/local/bin

# Install OpenBLAS for tvc
RUN apt update && \
    apt install -y libopenblas-dev libopenblas-base && \
    apt clean

# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-ionampliseq-1.0.0/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-ionampliseq-1.0.0 > nf-core-ionampliseq-1.0.0.yml

# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron
