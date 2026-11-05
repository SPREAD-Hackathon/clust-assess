#!/usr/bin/env bash
# Install clust-assess software for the local `standard` / `conda` profiles.
#   - conda env `clust-assess` from tools/env-clust-assess.yml
#   - ACCIO from github.com/mpgriffith/accio at a pinned commit (editable install into that env)
# MOB-suite and AMRFinderPlus are reused from existing envs (see nextflow.config params).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAMBA="${MAMBA:-$HOME/miniforge3/bin/mamba}"
ENV_NAME="${ENV_NAME:-clust-assess}"
ENV_PREFIX="${ENV_PREFIX:-$HOME/miniforge3/envs/$ENV_NAME}"
ACCIO_REPO="https://github.com/mpgriffith/accio.git"
ACCIO_COMMIT="${ACCIO_COMMIT:-6f6386ada75e0086a2d8d698271913d5b60547ad}"   # 2026-11-04 "read support"

if [ ! -d "$ENV_PREFIX" ]; then
  "$MAMBA" env create -y -n "$ENV_NAME" -f "$HERE/env-clust-assess.yml"
fi

if [ ! -d "$HERE/accio/.git" ]; then
  git clone --quiet "$ACCIO_REPO" "$HERE/accio"
fi
git -C "$HERE/accio" fetch --quiet origin
git -C "$HERE/accio" checkout --quiet "$ACCIO_COMMIT"
"$ENV_PREFIX/bin/python" -m pip install --quiet --no-deps -e "$HERE/accio"

echo "accio:  $("$ENV_PREFIX/bin/accio" --help >/dev/null 2>&1 && echo ok) (commit $ACCIO_COMMIT)"
echo "pling:  $("$ENV_PREFIX/bin/pling" --version 2>&1 | tail -1)"
echo "skani:  $("$ENV_PREFIX/bin/skani" --version 2>&1 | tail -1)"
echo "nucmer: $("$ENV_PREFIX/bin/nucmer" --version 2>&1 | tail -1)"
