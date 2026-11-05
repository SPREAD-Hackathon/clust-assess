#!/usr/bin/env bash
# clust-assess test on a CLIMB-BIG-DATA notebook (v3).
#   - two conda envs on /shared/team (task pods cannot see /home/jovyan):
#       CORE   = ACCIO deps + Pling 3 + plasnet + MUMmer 4, plus ACCIO as an editable install of a pinned
#                clone (ACCIO's setup.py does not package data/repetitive.dna.fas)
#       MOBENV = MOB-suite + AMRFinderPlus + their DBs  (reuses the env built by the first run)
#   - stub run on the notebook (local executor)
#   - test run as Kubernetes pods (CLIMB's executor); local run as a fallback if k8s fails
#   - inventory of the teammates' ACCIO databases and input sets for the real run
# Self-contained: writes its own env file and k8s config, so it works on the current `main`.
set -uo pipefail

SHARED="$(realpath "$HOME/shared-team" 2>/dev/null || echo "$HOME")"
WD="$(realpath -m "${WD:-$SHARED/clust-assess-run}")"
BRANCH="${BRANCH:-main}"
CORE="${CORE:-$WD/env-core}"
MOBENV="${MOBENV:-$WD/env-clust-assess}"
IMAGE="${IMAGE:-${JUPYTER_IMAGE_SPEC:-${JUPYTER_IMAGE:-ubuntu:24.04}}}"
ACCIO_COMMIT=6f6386ada75e0086a2d8d698271913d5b60547ad
REPORT="$WD/climb_report.txt"
mkdir -p "$WD" && cd "$WD" || exit 1
: > "$REPORT"
log()  { echo "== $*" | tee -a "$REPORT"; }
errs() { grep -E "ERROR|Caused by|Command error|exit status|denied|Exception" -A12 "$1" | head -45 >> "$REPORT"; }
CONDA_BIN="$(command -v mamba || command -v conda)"

log "machine"
{ date; echo "nproc: $(nproc)  cpu.max: $(cat /sys/fs/cgroup/cpu.max 2>/dev/null)  memory.max: $(cat /sys/fs/cgroup/memory.max 2>/dev/null)"
  echo "task image: $IMAGE"; } >> "$REPORT" 2>&1

log "clone $BRANCH"
[ -d clust-assess ] || git clone -q -b "$BRANCH" https://github.com/SPREAD-Hackathon/clust-assess.git
cd "$WD/clust-assess" && git pull -q; git log --oneline -1 >> "$REPORT"
echo "CLIMB k8s settings applied without -c: $(nextflow config -profile test 2>/dev/null | grep -cE "cephfs-shared-team|executor = 'k8s'") matching lines" >> "$REPORT"

log "core env ($CORE)"
if [ ! -x "$CORE/bin/python" ]; then
    cat > "$WD/env-core.yml" <<'EOF'
channels: [conda-forge, bioconda]
dependencies:
  - python=3.11
  - pandas
  - biopython
  - scipy
  - numpy
  - intervaltree
  - pysam
  - blast=2.17.0
  - mash=2.3
  - mummer4=4.0.1
  - samtools
  - bwa
  - minimap2
  - skani=0.4.0
  - pling=3.0.3
  - plasnet=0.7.0
  - pyyaml
  - pytest
  - git
  - pip
EOF
    "$CONDA_BIN" env create -y -p "$CORE" -f "$WD/env-core.yml" > "$WD/env-core.log" 2>&1
    echo "env create exit $? (log: $WD/env-core.log)" >> "$REPORT"
fi

log "ACCIO: editable install of a pinned clone"
ACCIO_DIR="$WD/clust-assess/tools/accio"
accio_data_ok() { "$CORE/bin/python" -c 'import importlib.util as u, pathlib, sys; o = u.find_spec("accio").origin; sys.exit(0 if (pathlib.Path(o).parent / "data" / "repetitive.dna.fas").exists() else 1)' 2>/dev/null; }
[ -d "$ACCIO_DIR/.git" ] || git clone -q https://github.com/mpgriffith/accio.git "$ACCIO_DIR"
git -C "$ACCIO_DIR" checkout -q "$ACCIO_COMMIT" && git -C "$ACCIO_DIR" log --oneline -1 >> "$REPORT"
accio_data_ok || "$CORE/bin/python" -m pip install -q --no-deps -e "$ACCIO_DIR" >> "$REPORT" 2>&1
{ accio_data_ok && echo "accio data files: ok" || echo "accio data files: MISSING"
  echo "accio imported from: $("$CORE/bin/python" -c 'import importlib.util as u; print(u.find_spec("accio").origin)')"; } >> "$REPORT" 2>&1

log "MOB-suite + AMRFinder env ($MOBENV)"
if [ ! -x "$MOBENV/bin/mob_typer" ]; then
    printf 'channels: [conda-forge, bioconda]\ndependencies: [python=3.11, mob_suite=3.1.9, ncbi-amrfinderplus=4.2.7]\n' > "$WD/env-mobsuite.yml"
    { "$CONDA_BIN" env create -y -p "$MOBENV" -f "$WD/env-mobsuite.yml" \
      && "$MOBENV/bin/mob_init" && "$MOBENV/bin/amrfinder" -u; } > "$WD/env-mobsuite.log" 2>&1
    echo "mobsuite env exit $? (log: $WD/env-mobsuite.log)" >> "$REPORT"
else
    echo "reusing existing env (databases already initialised)" >> "$REPORT"
fi

ENVARGS=(--clust_assess_env "$CORE" --mobsuite_env "$MOBENV" --amrfinder_env "$MOBENV")

log "unit tests"
"$CORE/bin/python" -m pytest -q tests 2>&1 | tail -1 >> "$REPORT"

log "stub run (notebook, local executor)"
nextflow run main.nf -profile standard,stub -stub-run -ansi-log false --max_cpus 2 "${ENVARGS[@]}" > "$WD/stub.log" 2>&1
echo "stub_exit=$?" >> "$REPORT"; errs "$WD/stub.log"

log "test run on Kubernetes pods (build 33-plasmid DB + assess 6 plasmids)"
cat > "$WD/climb_k8s.config" <<EOF
params.max_cpus = 4
process {
    executor     = 'k8s'
    container    = '$IMAGE'
    beforeScript = 'export PATH=$CORE/bin:$MOBENV/bin:\$PATH'
}
EOF
ETC=(); [ -f /etc/nextflow.config ] && ETC=(-c /etc/nextflow.config)
START=$(date +%s)
timeout 60m nextflow run main.nf -profile test "${ETC[@]}" -c "$WD/climb_k8s.config" -ansi-log false "${ENVARGS[@]}" > "$WD/test_k8s.log" 2>&1
K8S=$?; echo "test_k8s_exit=$K8S  wall=$(( $(date +%s) - START ))s" >> "$REPORT"; errs "$WD/test_k8s.log"

if [ "$K8S" -ne 0 ]; then
    log "fallback: test run on the notebook (local executor)"
    NCPU=$(nproc); [ "$NCPU" -gt 4 ] && NCPU=4
    nextflow run main.nf -profile standard,test -ansi-log false --max_cpus "$NCPU" "${ENVARGS[@]}" > "$WD/test_local.log" 2>&1
    echo "test_local_exit=$?" >> "$REPORT"; errs "$WD/test_local.log"
fi

log "results"
{ echo "--- final_output.csv"; head -8 results_test/final_output.csv
  echo "--- inputs"; cut -f1-4 results_test/inputs.tsv
  echo "--- subcommunities"; cut -f1,3,9,10,15 results_test/subcommunities.tsv
  echo "--- software_versions.yml"; cat results_test/runningstat/software_versions.yml
  echo "--- run stats"; grep WorkflowStats .nextflow.log | tail -1 \
      | grep -oE "succeededCount=[0-9]+|failedCount=[0-9]+|succeedDuration=[^;]+|peakCpus=[0-9]+|peakMemory=[^;]+"
  echo "--- task times"; cut -f4,6,9,10 results_test/runningstat/trace.txt 2>/dev/null | head -20
} >> "$REPORT" 2>&1

log "teammate databases and inputs (for the real run)"
{ for d in "$SHARED"/clust-assess/accio_db_*; do
      echo "--- $d"; ls "$d" | tr '\n' ' '; echo
      "$CORE/bin/python" bin/clust_assess.py verify-db --db "$d" --out /dev/stdout 2>&1 | tail -2; echo
  done
  for l in "$SHARED"/clust-assess/accio_create_*.log; do echo "--- tail $l"; tail -4 "$l"; done
  for s in Danish_plasmids_fasta raw_15_plasmid_set pling_15_plasmid_set; do
      p="$SHARED/clust-assess/$s"; [ -e "$p" ] || continue
      echo "--- $s: $(find "$p" -maxdepth 2 -type f | wc -l) files, $(cat $(find "$p" -maxdepth 2 -type f -name '*.f*a*') 2>/dev/null | grep -c '^>') FASTA records"
      find "$p" -maxdepth 2 -type f | head -4
  done
  echo "--- denmark_plsdb_accession.tsv"; head -3 "$SHARED/clust-assess/denmark_plsdb_accession.tsv"; wc -l < "$SHARED/clust-assess/denmark_plsdb_accession.tsv"
} >> "$REPORT" 2>&1

echo; echo "Done. Paste the contents of: $REPORT"
