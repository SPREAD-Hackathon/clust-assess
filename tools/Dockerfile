# clust-assess image for CLIMB (-profile climb / docker).
#   docker build -t clust-assess:0.1.0 -f tools/Dockerfile tools/
# Pushing to a registry (e.g. GHCR for the team) is a separate, explicit step.
FROM mambaorg/micromamba:1.5.10
COPY --chown=$MAMBA_USER:$MAMBA_USER env-clust-assess-full.yml /tmp/env.yml
RUN micromamba install -y -n base -f /tmp/env.yml && micromamba clean -a -y
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN mob_init && amrfinder -u
ENV PATH=/opt/conda/bin:$PATH
