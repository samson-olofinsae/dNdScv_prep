# Dockerfile
FROM mambaorg/micromamba:1.5.8

# Create env with all tools
USER root
RUN micromamba create -y -n dndscv-prep -c conda-forge -c bioconda \
    python=3.10 pandas bwa samtools bcftools htslib tabix pigz \
 && micromamba clean --all --yes

# Activate env by default
SHELL ["/bin/bash", "-lc"]
ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/envs/dndscv-prep/bin:$PATH

# Workdir inside container
WORKDIR /work

# Default entrypoint (still lets you override command)
ENTRYPOINT ["/bin/bash", "-lc"]
CMD ["python3 mutation-caller.py --help"]
