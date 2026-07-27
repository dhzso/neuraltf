FROM python:3.11-slim-bookworm

# -----------------------------------------------------------------------------
# System packages: build toolchain and utilities needed for scientific libraries
# -----------------------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    g++ \
    make \
    git \
    curl \
    wget \
    pkg-config \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# -----------------------------------------------------------------------------
# Python tooling: upgrade pip, setuptools, wheel
# -----------------------------------------------------------------------------
RUN pip install --upgrade pip setuptools wheel

# -----------------------------------------------------------------------------
# JupyterLab: development notebook interface
# -----------------------------------------------------------------------------
RUN pip install "jupyterlab>=4.0,<5.0"

# -----------------------------------------------------------------------------
# Non-root user: bioforge (UID 1000, GID 1000)
# -----------------------------------------------------------------------------
RUN groupadd --gid 1000 bioforge \
    && useradd --uid 1000 --gid 1000 --create-home --home-dir /home/bioforge bioforge

# -----------------------------------------------------------------------------
# Workspace: mounted at runtime via docker-compose
# -----------------------------------------------------------------------------
RUN mkdir -p /workspace && chown bioforge:bioforge /workspace

WORKDIR /workspace

# -----------------------------------------------------------------------------
# BioForge: editable install with bioinformatics extra (Layer 4)
# -----------------------------------------------------------------------------
COPY --chown=bioforge:bioforge . /workspace
RUN pip install -e ".[bio,dev]"

# -----------------------------------------------------------------------------
# Runtime: JupyterLab as non-root user
# -----------------------------------------------------------------------------
EXPOSE 8888

USER bioforge

CMD jupyter lab --ip=0.0.0.0 --port=8888 --no-browser --ServerApp.token='' --ServerApp.password='' --ServerApp.root_dir=/workspace
