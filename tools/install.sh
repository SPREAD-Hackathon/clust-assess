#!/usr/bin/env bash
# Install clust-assess software for the `standard` / `climb` profiles. Two conda envs are needed,
# because mob_suite 3.1.x (numpy<1.23.5) and pling 3 (numpy>=1.26) cannot share one:
#   CORE_PREFIX  <- tools/env-clust-assess.yml + ACCIO (editable install of tools/accio at a pinned commit)
#   MOB_PREFIX   <- tools/env-mobsuite.yml      (MOB-suite + AMRFinderPlus, plus their databases)
# ACCIO is installed editable because its setup.py does not package data/repetitive.dna.fas.
# Existing envs are left alone. Pass the prefixes to Nextflow:
#   --clust_assess_env $CORE_PREFIX --mobsuite_env $MOB_PREFIX --amrfinder_env $MOB_PREFIX
# On CLIMB, put both prefixes and this clone under /shared/team (task pods cannot see /home/jovyan).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
MAMBA="${MAMBA:-$(command -v mamba || command -v conda || echo "$HOME/miniforge3/bin/mamba")}"
CORE_PREFIX="${CORE_PREFIX:-$HOME/miniforge3/envs/clust-assess}"
MOB_PREFIX="${MOB_PREFIX:-$HOME/miniforge3/envs/clust-assess-mobsuite}"
ACCIO_REPO="https://github.com/mpgriffith/accio.git"
ACCIO_COMMIT="${ACCIO_COMMIT:-6f6386ada75e0086a2d8d698271913d5b60547ad}"   # 2026-11-04 "read support"

if [ ! -x "$CORE_PREFIX/bin/python" ]; then
  "$MAMBA" env create -y -p "$CORE_PREFIX" -f "$HERE/env-clust-assess.yml"
fi
[ -d "$HERE/accio/.git" ] || git clone --quiet "$ACCIO_REPO" "$HERE/accio"
git -C "$HERE/accio" fetch --quiet origin
git -C "$HERE/accio" checkout --quiet "$ACCIO_COMMIT"
"$CORE_PREFIX/bin/python" -m pip install --quiet --no-deps -e "$HERE/accio"

if [ ! -x "$MOB_PREFIX/bin/mob_typer" ]; then
  "$MAMBA" env create -y -p "$MOB_PREFIX" -f "$HERE/env-mobsuite.yml"
  "$MOB_PREFIX/bin/mob_init"
  "$MOB_PREFIX/bin/amrfinder" -u
fi

echo "accio:     $("$CORE_PREFIX/bin/accio" --help >/dev/null 2>&1 && echo ok || echo MISSING) (commit ${ACCIO_COMMIT:0:7})"
echo "pling:     $("$CORE_PREFIX/bin/pling" --version 2>&1 | tail -1)"
echo "skani:     $("$CORE_PREFIX/bin/skani" --version 2>&1 | tail -1)"
echo "nucmer:    $("$CORE_PREFIX/bin/nucmer" --version 2>&1 | tail -1)"
echo "mob_typer: $("$MOB_PREFIX/bin/mob_typer" --version 2>&1 | tail -1)"
echo "amrfinder: $("$MOB_PREFIX/bin/amrfinder" --version 2>&1 | tail -1)"
echo
echo "nextflow run main.nf ... --clust_assess_env $CORE_PREFIX --mobsuite_env $MOB_PREFIX --amrfinder_env $MOB_PREFIX"
