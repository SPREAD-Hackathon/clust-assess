# clust-assess image (-profile docker). Not needed on CLIMB, which uses conda envs on /shared/team.
#   docker build -t clust-assess:0.1.0 -f tools/Dockerfile tools/
# Two envs: mob_suite 3.1.x and pling 3 have conflicting numpy pins.
# ACCIO is installed editable from a pinned clone: its setup.py does not package data/repetitive.dna.fas.
FROM mambaorg/micromamba:1.5.10
ARG ACCIO_COMMIT=6f6386ada75e0086a2d8d698271913d5b60547ad
COPY --chown=$MAMBA_USER:$MAMBA_USER env-clust-assess.yml env-mobsuite.yml /tmp/
RUN micromamba create -y -n clust-assess -f /tmp/env-clust-assess.yml \
 && micromamba create -y -n mobsuite -f /tmp/env-mobsuite.yml \
 && micromamba clean -a -y
ENV PATH=/opt/conda/envs/clust-assess/bin:/opt/conda/envs/mobsuite/bin:$PATH
RUN git clone --quiet https://github.com/mpgriffith/accio.git /opt/accio \
 && git -C /opt/accio checkout --quiet $ACCIO_COMMIT \
 && python -m pip install --quiet --no-deps -e /opt/accio
RUN mob_init && amrfinder -u
