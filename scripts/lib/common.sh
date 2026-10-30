#!/usr/bin/env bash
set -euo pipefail

APIARC_REPO="${APIARC_REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
export APIARC_REPO

export APIARC_RUNTIME_DIR="${APIARC_RUNTIME_DIR:-$APIARC_REPO/.apiarc-runtime}"
export APIARC_LAUNCH_ENV="${APIARC_LAUNCH_ENV:-$APIARC_RUNTIME_DIR/launch}"
export APIARC_ANALYSIS_ENV="${APIARC_ANALYSIS_ENV:-$APIARC_REPO/analysis-env}"
export APIARC_PKGS_DIRS="${APIARC_PKGS_DIRS:-$APIARC_RUNTIME_DIR/pkgs}"
export APIARC_ENVS_DIRS="${APIARC_ENVS_DIRS:-$APIARC_RUNTIME_DIR/envs}"
export APIARC_CACHE_DIR="${APIARC_CACHE_DIR:-$APIARC_RUNTIME_DIR/cache}"
export APIARC_STATE_DIR="${APIARC_STATE_DIR:-$APIARC_RUNTIME_DIR/state}"
export APIARC_SERVER_LOG="${APIARC_SERVER_LOG:-$APIARC_STATE_DIR/streamlit.out}"
export APIARC_CONFIG="${APIARC_CONFIG:-$APIARC_REPO/config/config.yaml}"

export CONDA_PKGS_DIRS="$APIARC_PKGS_DIRS"
export CONDA_ENVS_DIRS="$APIARC_ENVS_DIRS"
export PIP_CACHE_DIR="$APIARC_CACHE_DIR/pip"

if [ -t 1 ]; then
    C_RESET=$'\033[0m'; C_INFO=$'\033[36m'; C_OK=$'\033[32m'
    C_WARN=$'\033[33m'; C_FAIL=$'\033[31m'
else
    C_RESET=""; C_INFO=""; C_OK=""; C_WARN=""; C_FAIL=""
fi

log()  { printf '%s[apiarc]%s %s\n' "$C_INFO" "$C_RESET" "$*"; }
ok()   { printf '%s[  ok  ]%s %s\n' "$C_OK" "$C_RESET" "$*"; }
warn() { printf '%s[ warn ]%s %s\n' "$C_WARN" "$C_RESET" "$*" >&2; }
die()  { printf '%s[fail ]%s %s\n' "$C_FAIL" "$C_RESET" "$*" >&2; exit 1; }

apiarc_version() {
    sed -n 's/^APP_VERSION = "\(.*\)"$/\1/p' "$APIARC_REPO/streamlit_app/utils/version.py" 2>/dev/null | head -n1
}

launcher_python() { echo "$APIARC_LAUNCH_ENV/bin/python"; }
launcher_streamlit() { echo "$APIARC_LAUNCH_ENV/bin/streamlit"; }
analysis_marker() { echo "$APIARC_STATE_DIR/analysis_env"; }
pid_file() { echo "$APIARC_STATE_DIR/app.pid"; }
port_file() { echo "$APIARC_STATE_DIR/app.port"; }

analysis_prefix() {
    local marker
    marker=$(analysis_marker)
    if [ -f "$marker" ]; then
        local value
        value=$(head -n1 "$marker" | tr -d '[:space:]')
        [ -n "$value" ] && { echo "$value"; return; }
    fi
    echo "$APIARC_ANALYSIS_ENV"
}

ensure_runtime_dirs() {
    mkdir -p "$APIARC_RUNTIME_DIR" "$APIARC_STATE_DIR" "$APIARC_PKGS_DIRS" \
             "$APIARC_ENVS_DIRS" "$APIARC_CACHE_DIR/pip" \
             "$APIARC_REPO/logs/runs" "$APIARC_REPO/workflow/data/RNAseq" \
             "$APIARC_REPO/workflow/data/ChIPseq"
}

write_analysis_marker() {
    local prefix="${1:-$APIARC_ANALYSIS_ENV}"
    ensure_runtime_dirs
    printf '%s\n' "$prefix" > "$(analysis_marker)"
}

read_runtime() {
    local file="$APIARC_STATE_DIR/runtime"
    if [ -f "$file" ]; then
        head -n1 "$file" | tr -d '[:space:]'
    else
        echo "conda"
    fi
}

write_runtime() {
    ensure_runtime_dirs
    printf '%s\n' "${1:-conda}" > "$APIARC_STATE_DIR/runtime"
}

conda_load() {
    local base
    if [ -n "${CONDA_EXE:-}" ] && [ -f "$(dirname "$(dirname "$CONDA_EXE")")/etc/profile.d/conda.sh" ]; then
        . "$(dirname "$(dirname "$CONDA_EXE")")/etc/profile.d/conda.sh"
        return 0
    fi
    for base in "$HOME/miniconda3" "$HOME/anaconda3" "$HOME/miniforge3" \
                "$HOME/mambaforge" "/opt/conda" "/usr/local/conda"; do
        if [ -f "$base/etc/profile.d/conda.sh" ]; then
            . "$base/etc/profile.d/conda.sh"
            return 0
        fi
    done
    command -v conda >/dev/null 2>&1 && return 0
    return 1
}

base_python() {
    local candidate
    for candidate in python3.10 python3 python; do
        if command -v "$candidate" >/dev/null 2>&1; then
            command -v "$candidate"
            return 0
        fi
    done
    return 1
}

find_free_port() {
    local start="${1:-8501}"
    local port="$start"
    while [ "$port" -lt $((start + 200)) ]; do
        if ! (exec 3<>"/dev/tcp/127.0.0.1/$port") 2>/dev/null; then
            echo "$port"
            return 0
        fi
        exec 3>&- 2>/dev/null || true
        port=$((port + 1))
    done
    echo "$start"
}

wait_ready() {
    local url="$1"
    local seconds="${2:-90}"
    local waited=0
    while [ "$waited" -lt "$seconds" ]; do
        if curl -fsS -o /dev/null "$url/_stcore/health" 2>/dev/null; then
            return 0
        fi
        sleep 2
        waited=$((waited + 2))
    done
    return 1
}

app_running() {
    local file pid
    file=$(pid_file)
    [ -f "$file" ] || return 1
    pid=$(head -n1 "$file" | tr -d '[:space:]')
    [ -n "$pid" ] || return 1
    kill -0 "$pid" 2>/dev/null
}

app_pid() {
    head -n1 "$(pid_file)" 2>/dev/null | tr -d '[:space:]'
}

banner() {
    cat <<EOF
APIARC v$(apiarc_version)
  repo      $APIARC_REPO
  runtime   $APIARC_RUNTIME_DIR
  launcher  $APIARC_LAUNCH_ENV
  analysis  $(analysis_prefix)
EOF
}
