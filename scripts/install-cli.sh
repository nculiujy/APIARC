#!/usr/bin/env bash
set -euo pipefail

SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export APIARC_REPO="$SRC"
source "$SRC/scripts/lib/common.sh"

RT_BIN="$APIARC_RUNTIME_DIR/bin"
RT_LIB="$APIARC_RUNTIME_DIR/lib"
RT_RUNTIME="$APIARC_RUNTIME_DIR/runtime"
LINK_DIR="${PREFIX:-$HOME/.local}/bin"

log "installing CLI"
log "  repo        $SRC"
log "  runtime dir $APIARC_RUNTIME_DIR"
log "  link dir    $LINK_DIR"

ensure_runtime_dirs
mkdir -p "$RT_BIN" "$RT_LIB" "$RT_RUNTIME"

install -m 644 "$SRC/scripts/lib/common.sh" "$SRC/scripts/lib/cd-banner.sh" "$RT_LIB/"
install -m 644 "$SRC"/runtime/*.sh "$RT_RUNTIME/"

for source_file in "$SRC"/cli/apiarc-app*; do
    name=$(basename "$source_file")
    install -m 755 "$source_file" "$RT_BIN/$name"
    sed -i "s|@@LIB@@|$RT_LIB/common.sh|g" "$RT_BIN/$name"
    sed -i "s|@@RUNTIME_DIR@@|$RT_RUNTIME|g" "$RT_BIN/$name"
    sed -i "s|@@REPO@@|$SRC|g" "$RT_BIN/$name"
done

mkdir -p "$LINK_DIR" 2>/dev/null || warn "cannot create $LINK_DIR, use $RT_BIN directly"
linked=0
for file in "$RT_BIN"/apiarc-app*; do
    [ -e "$file" ] || continue
    if ln -sfn "$file" "$LINK_DIR/$(basename "$file")" 2>/dev/null; then
        linked=$((linked + 1))
    else
        warn "symlink failed: $LINK_DIR/$(basename "$file")"
    fi
done
ok "$linked symlink(s) installed in $LINK_DIR"

case ":$PATH:" in
    *":$LINK_DIR:"*) ;;
    *) warn "$LINK_DIR is not in PATH; add: export PATH=\"$LINK_DIR:\$PATH\"" ;;
esac

RC_FILE="${APIARC_SHELL_RC:-$HOME/.bashrc}"
BANNER="$RT_LIB/cd-banner.sh"
if [ -f "$RC_FILE" ] && [ -w "$RC_FILE" ]; then
    sed -i '/^# >>> apiarc cd banner >>>$/,/^# <<< apiarc cd banner <<<$/d' "$RC_FILE" 2>/dev/null || true
    cat >> "$RC_FILE" <<EOF
# >>> apiarc cd banner >>>
if [ -n "\${PS1:-}" ] && [ -f "$BANNER" ]; then
    . "$BANNER"
    if declare -p PROMPT_COMMAND 2>/dev/null | grep -q 'declare -a'; then
        case " \${PROMPT_COMMAND[*]:-} " in
            *" _apiarc_cd_banner "*) ;;
            *) PROMPT_COMMAND+=(_apiarc_cd_banner) ;;
        esac
    else
        case "\${PROMPT_COMMAND:-}" in
            *_apiarc_cd_banner*) ;;
            *) PROMPT_COMMAND="_apiarc_cd_banner\${PROMPT_COMMAND:+;\$PROMPT_COMMAND}" ;;
        esac
    fi
fi
# <<< apiarc cd banner <<<
EOF
    ok "cd banner hook written to $RC_FILE"
else
    warn "skipped cd banner: $RC_FILE missing or not writable"
fi

ok "CLI installed in $RT_BIN"
ls -1 "$RT_BIN"/apiarc-app* 2>/dev/null | sed 's|.*/|  |'
