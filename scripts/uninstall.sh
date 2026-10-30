#!/usr/bin/env bash
set -euo pipefail

SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export APIARC_REPO="$SRC"
source "$SRC/scripts/lib/common.sh"

ASSUME_YES=0
KEEP_RUNTIME=0
KEEP_ANALYSIS=0

while [ $# -gt 0 ]; do
    case "$1" in
        -y|--yes)      ASSUME_YES=1; shift ;;
        --keep-runtime) KEEP_RUNTIME=1; shift ;;
        --keep-analysis) KEEP_ANALYSIS=1; shift ;;
        -h|--help)
            cat <<EOF
usage: scripts/uninstall.sh [-y] [--keep-runtime] [--keep-analysis]

  -y               do not ask for confirmation
  --keep-runtime   keep $APIARC_RUNTIME_DIR
  --keep-analysis  keep $APIARC_ANALYSIS_ENV

removes: ~/.local/bin/apiarc-app*, the .bashrc cd banner block,
         $APIARC_RUNTIME_DIR and $APIARC_ANALYSIS_ENV
keeps:   result/, workflow/data/, logs/, config/, workflow/resources/
EOF
            exit 0 ;;
        *) die "unknown option '$1'" ;;
    esac
done

LINK_DIR="${PREFIX:-$HOME/.local}/bin"
RC_FILE="${APIARC_SHELL_RC:-$HOME/.bashrc}"

if [ "$ASSUME_YES" != "1" ]; then
    cat <<EOF
This will remove:
  - $LINK_DIR/apiarc-app* (symlinks)
  - the apiarc cd banner block in $RC_FILE
  - ${APIARC_RUNTIME_DIR}$([ "$KEEP_RUNTIME" = "1" ] && echo " (kept)")
  - ${APIARC_ANALYSIS_ENV}$([ "$KEEP_ANALYSIS" = "1" ] && echo " (kept)")
Data in workflow/data, result, logs and workflow/resources stays untouched.
EOF
    read -r -p "Continue? [y/N] " answer
    case "$answer" in
        y|Y|yes|YES) ;;
        *) warn "aborted"; exit 0 ;;
    esac
fi

if [ -f "$APIARC_STATE_DIR/app.pid" ] && kill -0 "$(app_pid)" 2>/dev/null; then
    source "$SRC/runtime/conda.sh"
    runtime_stop
fi

removed=0
for file in "$LINK_DIR"/apiarc-app*; do
    [ -e "$file" ] || [ -L "$file" ] || continue
    rm -f "$file" && removed=$((removed + 1))
done
ok "removed $removed symlink(s) from $LINK_DIR"

if [ -f "$RC_FILE" ] && [ -w "$RC_FILE" ]; then
    sed -i '/^# >>> apiarc cd banner >>>$/,/^# <<< apiarc cd banner <<<$/d' "$RC_FILE" 2>/dev/null || true
    ok "removed cd banner block from $RC_FILE"
fi

if [ "$KEEP_RUNTIME" = "1" ]; then
    warn "kept $APIARC_RUNTIME_DIR"
else
    rm -rf "$APIARC_RUNTIME_DIR"
    ok "removed $APIARC_RUNTIME_DIR"
fi

if [ "$KEEP_ANALYSIS" = "1" ]; then
    warn "kept $APIARC_ANALYSIS_ENV"
else
    rm -rf "$APIARC_ANALYSIS_ENV"
    ok "removed $APIARC_ANALYSIS_ENV"
fi

ok "uninstall finished"
