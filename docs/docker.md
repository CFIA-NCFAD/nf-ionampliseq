# Docker container for CFIA-NCFAD/nf-ionampliseq

- **Container name:** `ghcr.io/cfia-ncfad/nf-ionampliseq`

This container hosted on the GitHub Container Registry at `ghcr.io` and specified in the [Dockerfile](../Dockerfile) provides the Thermo Fisher TMAP and TVC binaries for read mapping and variant calling of Ion Torrent data.

## Container Contents

The Docker container includes the following key components:

### Core Tools

- **TMAP (Torrent Mapping Alignment Program)**: Thermo Fisher's read mapping tool for Ion Torrent data
- **TVC (Torrent Variant Caller)**: Thermo Fisher's variant calling tool for Ion Torrent data
- **TVC Utils**: Additional utilities for TVC processing
- **TVC Assembly**: Assembly tools for TVC

### Supporting Tools

- **Samtools v1.22.1**: For BAM file manipulation and indexing
- **OpenBLAS**: Mathematical library required for TVC operations

### Base Image

- **Ubuntu 22.04 (Plucky)**: Provides the underlying Linux environment

## Usage in the Pipeline

The container is used by the following Nextflow processes:

- **TMAP process**: Read mapping and alignment of Ion Torrent sequencing data
- **TVC process**: Variant calling and post-processing of aligned reads

## Build Docker image

```bash
DOCKER_BUILDKIT=1 docker build -t ghcr.io/cfia-ncfad/nf-ionampliseq:latest .
```

Should see something like the following:

```text
[+] Building 68.6s (16/16) FINISHED                                                                                                                                                            docker:default
 => [internal] load build definition from Dockerfile                                                                                                                                                     0.0s
 => => transferring dockerfile: 1.78kB                                                                                                                                                                   0.0s
 => [internal] load metadata for docker.io/peterk87/tvc-tmap-torrent-suite:v1.0.0                                                                                                                        0.8s
 => [internal] load metadata for docker.io/library/ubuntu:plucky-20250714                                                                                                                                0.0s
 => [internal] load .dockerignore                                                                                                                                                                         0.0s
 => => transferring context: 2B                                                                                                                                                                           0.0s
 => CACHED [samtools-build 1/4] FROM docker.io/library/ubuntu:plucky-20250714                                                                                                                            0.0s
 => FROM docker.io/peterk87/tvc-tmap-torrent-suite:v1.0.0@sha256:4f8384a0cd3cac975c8d80ef67f163023d08ae25711f96733a9bf189a78a1018                                                                        3.7s
 => => resolve docker.io/peterk87/tvc-tmap-torrent-suite:v1.0.0@sha256:4f8384a0cd3cac975c8d80ef67f163023d08ae25711f96733a9bf189a78a1018                                                                  0.0s
 => => sha256:81d180f088e42bbb8d1e17d9c784ccd05504d3bc67aac9a7786851954491bfa2 4.63kB / 4.63kB                                                                                                           0.0s
 => => sha256:171857c49d0f5e2ebf623e6cb36a8bcad585ed0c2aa99c87a055df034c1e5848 26.70MB / 26.70MB                                                                                                           0.9s
 => => sha256:61e52f862619ab016d3bcfbd78e5c7aaaa1989b4c295e6dbcacddd2d7b93e1f5 162B / 162B                                                                                                               0.2s
 => => sha256:4f8384a0cd3cac975c8d80ef67f163023d08ae25711f96733a9bf189a78a1018 2.00kB / 2.00kB                                                                                                           0.0s
 => => sha256:419640447d267f068d2f84a093cb13a56ce77e130877f5b8bdb4294f4a90a84f 852B / 2.00kB                                                                                                               0.4s
 => => sha256:10c6a7e2448635e15f0a64f1097a8d567f20bf2d22217aad02c6ba73f80de294 15.00MB / 15.00MB                                                                                                           0.8s
 => => sha256:7d0145695d8d7ab3d361195ecd15b155a8b7362de4b2323b184cc0753e7ca125 3.49MB / 3.49MB                                                                                                           0.8s
 => => sha256:24eada07801e98b586b5f9123c89ccd7f2aa198b81aaf13027deafa0c337c6f1 2.57MB / 2.57MB                                                                                                           1.1s
 => => sha256:3ed2c085af2450dee706909a0f6b888d6fec33942fe9533aea742fe40d87fab3 1.87MB / 1.87MB                                                                                                           1.0s
 => => extracting sha256:171857c49d0f5e2ebf623e6cb36a8bcad585ed0c2aa99c87a055df034c1e5848                                                                                                                0.8s
 => => sha256:bb00ac787006b3830b434f08203251779a569f771a9b998964618e226807f225 49.53MB / 49.53MB                                                                                                           2.2s
 => => extracting sha256:419640447d267f068d2f84a093cb13a56ce77e130877f5b8bdb4294f4a90a84f                                                                                                                0.0s
 => => extracting sha256:61e52f862619ab016d3bcfbd78e5c7aaaa1989b4c295e6dbcacddd2d7b93e1f5                                                                                                                0.0s
 => => extracting sha256:10c6a7e2448635e15f0a64f1097a8d567f20bf2d22217aad02c6ba73f80de294                                                                                                               0.3s
 => => extracting sha256:7d0145695d8d7ab3d361195ecd15b155a8b7362de4b2323b184cc0753e7ca125                                                                                                               0.1s
 => => extracting sha256:24eada07801e98b586b5f9123c89ccd7f2aa198b81aaf13027deafa0c337c6f1                                                                                                               0.0s
 => => extracting sha256:3ed2c085af2450dee706909a0f6b888d6fec33942fe9533aea742fe40d87fab3                                                                                                               0.0s
 => => extracting sha256:bb00ac787006b3830b434f08203251779a569f771a9b998964618e226807f225                                                                                                               0.7s
 => [samtools-build 2/4] RUN --mount=type=cache,target=/var/cache/apt     apt update &&     apt install -y         build-essential zlib1g-dev libncurses5-dev         libbz2-dev liblzma-dev libcurl4-  30.2s
 => => [stage-1 2/7] COPY --from=peterk87/tvc-tmap-torrent-suite:v1.0.0 /usr/local/bin/tvc /usr/local/bin                                                                                                   0.1s
 => => [stage-1 3/7] COPY --from=peterk87/tvc-tmap-torrent-suite:v1.0.0 /usr/local/bin/tvcutils /usr/local/bin                                                                                              0.1s
 => => [stage-1 4/7] COPY --from=peterk87/tvc-tmap-torrent-suite:v1.0.0 /usr/local/bin/tvcassembly /usr/local/bin                                                                                           0.1s
 => => [stage-1 5/7] COPY --from=peterk87/tvc-tmap-torrent-suite:v1.0.0 /usr/local/bin/tmap /usr/local/bin                                                                                                  0.1s
 => [samtools-build 3/4] WORKDIR /tmp                                                                                                                                                                    0.1s
 => => [samtools-build 4/4] RUN wget https://github.com/samtools/samtools/releases/download/1.22.1/samtools-1.22.1.tar.bz2 &&     tar -xjf samtools-1.22.1.tar.bz2 &&     cd samtools-1.22.1 &&     make   28.7s
 => [stage-1 6/7] COPY --from=samtools-build /usr/local/bin/samtools /usr/local/bin/samtools                                                                                                             0.7s
 => [stage-1 7/7] RUN --mount=type=cache,target=/var/cache/apt     apt update &&     apt install -y libopenblas-dev &&     apt clean &&     rm -rf /var/lib/apt/lists/*                                  7.8s
 => => exporting to image                                                                                                                                                                                 0.7s
 => => exporting to image                                                                                                                                                                                 0.7s
 => => writing image sha256:10eb4d86f6a6aac7e2b57f460892da21536cb4010a8d5e4ce56f5eaa58c8944d                                                                                                             0.0s
 => => naming to docker.io/cfia-ncfad/nf-ionampliseq:2.0.0                                                                                                                                               0.0s
```

## Multi-stage Build

The Dockerfile uses a multi-stage build approach:

1. **samtools-build stage**: Compiles Samtools from source to ensure compatibility
2. **Final stage**: Combines TMAP/TVC binaries with compiled Samtools and system dependencies

This approach optimizes the build process and ensures all components are properly integrated.

## Pushing new image to ghcr.io

1. Setup a GitHub personal access token
2. Docker login to the `ghcr.io`

    ```bash
    export CR_PAT="GH personal access token"  # 
    echo $CR_PAT | docker login ghcr.io -u $GH_USERNAME --password-stdin
    ```

3. Tag image appropriately

    ```bash
    docker tag cfia-ncfad/nf-ionampliseq:latest ghcr.io/cfia-ncfad/nf-ionampliseq:latest
    ```

4. Push the image to GitHub Container Registry

    ```bash
    docker push ghcr.io/cfia-ncfad/nf-ionampliseq:latest
    ```

The new image should be able to be used within the Nextflow pipeline using the new tag. You can test with:

```bash
# Test the container locally
docker run --rm ghcr.io/cfia-ncfad/nf-ionampliseq:latest tmap --version
docker run --rm ghcr.io/cfia-ncfad/nf-ionampliseq:latest tvc --version
docker run --rm ghcr.io/cfia-ncfad/nf-ionampliseq:latest samtools --version
```

## Container Configuration

### Environment Variables

The container doesn't require specific environment variables to function, but you can override the default settings through Nextflow parameters.

### Resource Requirements

- **Memory**: Recommended minimum 8GB RAM for TMAP/TVC processes
- **CPU**: Multi-threading supported for both TMAP and TVC
- **Storage**: Temporary storage needed for intermediate BAM files

### Volume Mounts

When running the container directly (not through Nextflow), you may need to mount:

- Input data directories
- Output directories
- Reference genome files
- BED files for targeted analysis

## Troubleshooting

### Common Issues

1. **Permission Errors**: Ensure the container has appropriate permissions to read/write mounted volumes
2. **Memory Issues**: TMAP and TVC can be memory-intensive; increase container memory limits if needed
3. **Version Compatibility**: Ensure the container version matches the pipeline version requirements

### Debugging

To debug issues within the container:

```bash
# Run interactive shell
docker run -it --rm ghcr.io/cfia-ncfad/nf-ionampliseq:latest /bin/bash

# Check tool versions
tmap --version
tvc --version
samtools --version

# Verify file permissions
ls -la /usr/local/bin/
```

## Security Considerations

- The container runs as root by default (Ubuntu base image)
- Consider using a non-root user for production deployments
- Regularly update the base Ubuntu image for security patches
- The container only contains the necessary tools and libraries for the pipeline

## Alternative Container Engines

While this documentation focuses on Docker, the pipeline also supports:

- **Singularity/Apptainer**: For HPC environments
- **Podman**: Docker-compatible alternative
- **Charliecloud**: For secure container execution

See the pipeline configuration files for specific container engine settings.
