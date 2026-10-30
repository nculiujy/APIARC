#!/usr/bin/env bash
set -euo pipefail

SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export APIARC_REPO="$SRC"
source "$SRC/scripts/lib/common.sh"
source "$SRC/runtime/conda.sh"

mkdir -p "$APIARC_STATE_DIR"
printf '%s\n' "$SRC" > "$APIARC_STATE_DIR/repo"

runtime_install
provision_analysis

banner
ok "analysis environment is ready; start the dashboard with apiarc-app"
