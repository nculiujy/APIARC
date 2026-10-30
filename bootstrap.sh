#!/usr/bin/env bash
set -euo pipefail

SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export APIARC_REPO="$SRC"
source "$SRC/scripts/lib/common.sh"
source "$SRC/runtime/conda.sh"

FULL=0

for arg in "$@"; do
    case "$arg" in
        --full)      FULL=1 ;;
        --no-full)   FULL=0 ;;
        -h|--help)
            cat <<EOF
usage: ./bootstrap.sh [--full]

  --full   also build the analysis environment (snakemake + bioconda tools)

environment:
  APIARC_RUNTIME_DIR   runtime directory (default $APIARC_RUNTIME_DIR)
  APIARC_LAUNCH_ENV    launcher prefix  (default $APIARC_LAUNCH_ENV)
  APIARC_ANALYSIS_ENV  analysis prefix  (default $APIARC_ANALYSIS_ENV)
EOF
            exit 0 ;;
        *) die "unknown option '$arg'; see ./bootstrap.sh --help" ;;
    esac
done

cat <<EOF

============================================================
APIARC v$(apiarc_version) bootstrap
  repo     $SRC
  runtime  $APIARC_RUNTIME_DIR
============================================================

EOF

ensure_runtime_dirs
printf '%s\n' "$SRC" > "$APIARC_STATE_DIR/repo"

runtime_install
bash "$SRC/scripts/install-cli.sh"

if [ "$FULL" = "1" ]; then
    log "building analysis environment (snakemake + bioconda tools) ..."
    provision_analysis
fi

cat <<EOF

  1) start the dashboard
        apiarc-app
     first run without a full bootstrap:
        ./bootstrap.sh --full

  2) prepare data in the browser
        workflow/data/RNAseq/<dataset>/
        workflow/data/ChIPseq/<dataset>/
     or download runs from GEO on the same page

  3) run the pipelines from the dashboard, then read the results

  status   apiarc-app-status
  logs     apiarc-app-logs -f
  update   apiarc-app-update all
  remove   apiarc-app-uninstall

============================================================
EOF
