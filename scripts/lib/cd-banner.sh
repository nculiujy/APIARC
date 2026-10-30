#!/usr/bin/env bash
_APIARC_BANNER_LIB="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
_APIARC_BANNER_REPO="$(cat "$_APIARC_BANNER_LIB/../state/repo" 2>/dev/null || true)"
_APIARC_BANNER_LAST=""

_apiarc_cd_banner() {
    [ "${APIARC_NO_CD_BANNER:-0}" = "1" ] && return 0
    [ -n "$_APIARC_BANNER_REPO" ] || return 0
    [ "$PWD" = "$_APIARC_BANNER_REPO" ] || return 0
    [ "$_APIARC_BANNER_LAST" = "$PWD" ] && return 0
    _APIARC_BANNER_LAST="$PWD"
    local version
    version=$(sed -n 's/^APP_VERSION = "\(.*\)"$/\1/p' "$_APIARC_BANNER_REPO/streamlit_app/utils/version.py" 2>/dev/null | head -n1)
    printf '\n\033[36mAPIARC v%s\033[0m\n' "${version:-?}"
    printf '  %-22s %s\n' "apiarc-app" "start the dashboard"
    printf '  %-22s %s\n' "apiarc-app-status" "launcher / analysis / process status"
    printf '  %-22s %s\n' "apiarc-app-logs -f" "follow the server log"
    printf '  %-22s %s\n' "apiarc-app-update all" "rebuild launcher and analysis env"
    printf '  %-22s %s\n' "apiarc-app-shell" "shell with analysis tools on PATH"
    printf '  %-22s %s\n' "apiarc-app-stop" "stop the dashboard"
    printf '  data: workflow/data/<omics>/<dataset>   results: result/<omics>/<dataset>\n\n'
}
