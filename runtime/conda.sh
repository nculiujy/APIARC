#!/usr/bin/env bash
_APIARC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
for _candidate in "$_APIARC_DIR/../scripts/lib/common.sh" "$_APIARC_DIR/../lib/common.sh"; do
    if [ -f "$_candidate" ]; then
        source "$_candidate"
        break
    fi
done
[ -n "${APIARC_REPO:-}" ] || { echo "[fail] cannot locate scripts/lib/common.sh" >&2; exit 1; }

LAUNCHER_PACKAGES=(
    "streamlit>=1.36"
    "pandas>=2.0"
    "numpy>=1.24"
    "pyyaml>=6.0"
    "altair>=5.0"
    # The signal tab needs client-side zoom, pan and a crosshair; altair only
    # renders a static chart, so that tab is drawn with plotly.
    "plotly>=5.20"
)

LAUNCHER_OPTIONAL=(
    "pyBigWig>=0.3.22"
    "pymupdf>=1.24"
)

launcher_ready() {
    [ -x "$(launcher_streamlit)" ] && [ -x "$(launcher_python)" ]
}

_install_launcher_pip() {
    local python="$1"
    "$python" -m pip install --upgrade pip wheel setuptools >/dev/null 2>&1 || warn "pip self-upgrade failed"
    "$python" -m pip install "${LAUNCHER_PACKAGES[@]}" || return 1
    "$python" -m pip install "${LAUNCHER_OPTIONAL[@]}" || warn "optional packages (pyBigWig/pymupdf) skipped"
}

_create_launcher_conda() {
    conda_load || return 1
    log "creating launcher environment with conda: $APIARC_LAUNCH_ENV"
    conda create -y -p "$APIARC_LAUNCH_ENV" --override-channels -c conda-forge python=3.10 pip >/dev/null
    _install_launcher_pip "$(launcher_python)"
}

_create_launcher_venv() {
    local base
    base=$(base_python) || return 1
    log "creating launcher environment with venv: $APIARC_LAUNCH_ENV"
    "$base" -m venv --clear "$APIARC_LAUNCH_ENV"
    _install_launcher_pip "$(launcher_python)"
}

runtime_install() {
    ensure_runtime_dirs
    write_runtime conda
    if launcher_ready; then
        ok "launcher already present: $APIARC_LAUNCH_ENV"
        return 0
    fi
    if _create_launcher_conda; then
        ok "launcher installed (conda)"
    elif _create_launcher_venv; then
        ok "launcher installed (venv)"
    else
        die "launcher installation failed; install conda or python3-venv and retry"
    fi
    "$(launcher_python)" -c "import streamlit" || die "streamlit import failed in launcher"
    ok "launcher python: $("$(launcher_python)" -V 2>&1)"
}

_conda_env_create() {
    local mode="$1"
    local args=(-y -p "$APIARC_ANALYSIS_ENV" -f "$APIARC_REPO/environment.yml")
    [ "$mode" = "offline" ] && args+=(--offline)
    conda_load || return 1
    log "conda env create ($mode) -> $APIARC_ANALYSIS_ENV"
    conda env create "${args[@]}"
}

provision_analysis() {
    ensure_runtime_dirs
    if [ -x "$APIARC_ANALYSIS_ENV/bin/snakemake" ]; then
        write_analysis_marker "$APIARC_ANALYSIS_ENV"
        repair_prefix_scripts
        ok "analysis environment already present: $APIARC_ANALYSIS_ENV"
    else
        if _conda_env_create offline; then
            ok "analysis environment created (offline cache)"
        elif _conda_env_create online; then
            ok "analysis environment created (network)"
        else
            die "analysis environment creation failed; check environment.yml and network access"
        fi
        write_analysis_marker "$APIARC_ANALYSIS_ENV"
        verify_analysis
    fi
    if ! prebuild_stage_envs; then
        warn "stage environments incomplete; the pipeline creates them on first use"
    fi
}

repair_prefix_scripts() {
    local prefix bin file first magic interp target fixed=0
    prefix=$(analysis_prefix)
    bin="$prefix/bin"
    if [ ! -d "$bin" ]; then
        return 0
    fi
    for file in "$bin"/*; do
        if [ ! -f "$file" ]; then
            continue
        fi
        magic=""
        read -r -n 2 magic 2>/dev/null < "$file" || continue
        if [ "$magic" != "#!" ]; then
            continue
        fi
        first=""
        IFS= read -r first < "$file" || continue
        interp="${first:2}"
        interp="${interp%% *}"
        if [ -x "$interp" ]; then
            continue
        fi
        case "$interp" in
            */bin/*) ;;
            *) continue ;;
        esac
        target="$bin/${interp##*/}"
        if [ ! -x "$target" ]; then
            continue
        fi
        sed -i -E "1s|^#![^ ]*/bin/|#!$bin/|" "$file"
        fixed=$((fixed + 1))
    done
    if [ "$fixed" -gt 0 ]; then
        warn "repaired $fixed scripts in $bin (interpreter was outside $prefix)"
    fi
    return 0
}

prebuild_stage_envs() {
    local prefix python env_root conda_bin="" run_path
    prefix=$(analysis_prefix)
    python="$prefix/bin/python"
    env_root="$APIARC_REPO/workflow/envs"
    if [ ! -x "$python" ]; then
        warn "analysis environment python missing: $python"
        return 1
    fi
    if [ ! -d "$env_root" ]; then
        warn "no stage environment definitions: $env_root"
        return 1
    fi
    if ! conda_load; then
        warn "conda not found; stage environments were not built"
        return 1
    fi
    if [ -n "${CONDA_EXE:-}" ]; then
        conda_bin=$(dirname "$CONDA_EXE")
    fi
    run_path="$prefix/bin:$PATH"
    if [ -n "$conda_bin" ]; then
        run_path="$conda_bin:$run_path"
    fi
    log "building stage environments -> $prefix/envs"
    if PATH="$run_path" "$python" "$APIARC_REPO/scripts/prepare-envs.py" \
        --env-root "$env_root" --prefix "$prefix/envs"; then
        ok "stage environments ready: $prefix/envs"
        return 0
    fi
    return 1
}

verify_analysis() {
    local prefix bin
    prefix=$(analysis_prefix)
    bin="$prefix/bin"
    local missing=()
    for tool in snakemake conda; do
        [ -x "$bin/$tool" ] || command -v "$tool" >/dev/null 2>&1 || missing+=("$tool")
    done
    for tool in prefetch fasterq-dump; do
        [ -x "$bin/$tool" ] || missing+=("$tool")
    done
    if [ "${#missing[@]}" -gt 0 ]; then
        warn "analysis environment is missing: ${missing[*]}"
        warn "rerun scripts/build-analysis-env.sh after fixing environment.yml"
        return 1
    fi
    ok "analysis tools: $("$bin/snakemake" --version 2>/dev/null | head -n1)"
}

runtime_build_all() {
    runtime_install
    provision_analysis
}

runtime_run() {
    local port="${1:-}"
    local bind="${2:-0.0.0.0}"
    ensure_runtime_dirs
    launcher_ready || die "launcher missing; run ./bootstrap.sh first"
    if app_running; then
        warn "APIARC is already running (pid $(app_pid), port $(cat "$(port_file)" 2>/dev/null || echo '?'))"
        return 0
    fi
    [ -n "$port" ] || port=$(find_free_port 8501)
    local prefix
    prefix=$(analysis_prefix)
    : > "$APIARC_SERVER_LOG"
    (
        cd "$APIARC_REPO"
        APIARC_REPO="$APIARC_REPO" \
        APIARC_RUNTIME_DIR="$APIARC_RUNTIME_DIR" \
        APIARC_LAUNCH_ENV="$APIARC_LAUNCH_ENV" \
        APIARC_ANALYSIS_ENV="$prefix" \
        nohup "$(launcher_python)" -m streamlit run "$APIARC_REPO/streamlit_app/app.py" \
            --server.port "$port" \
            --server.address "$bind" \
            --server.headless true \
            --server.enableCORS false \
            --server.enableXsrfProtection false \
            --browser.gatherUsageStats false \
            >> "$APIARC_SERVER_LOG" 2>&1 &
        echo $! > "$(pid_file)"
    )
    printf '%s\n' "$port" > "$(port_file)"
    log "starting streamlit on $bind:$port (pid $(app_pid)) ..."
    if wait_ready "http://127.0.0.1:$port" 90; then
        ok "APIARC ready: http://127.0.0.1:$port"
        if [ "$bind" = "0.0.0.0" ]; then
            local host
            host=$(hostname -I 2>/dev/null | awk '{print $1}')
            [ -n "$host" ] && ok "LAN address: http://$host:$port"
        fi
    else
        warn "streamlit did not become ready; last log lines:"
        tail -n 40 "$APIARC_SERVER_LOG" >&2 || true
        return 1
    fi
}

runtime_stop() {
    if ! app_running; then
        rm -f "$(pid_file)"
        warn "APIARC is not running"
        return 0
    fi
    local pid
    pid=$(app_pid)
    log "stopping APIARC (pid $pid) ..."
    kill -TERM "$pid" 2>/dev/null || true
    local waited=0
    while kill -0 "$pid" 2>/dev/null && [ "$waited" -lt 15 ]; do
        sleep 1
        waited=$((waited + 1))
    done
    if kill -0 "$pid" 2>/dev/null; then
        warn "graceful stop timed out, sending SIGKILL"
        kill -KILL "$pid" 2>/dev/null || true
    fi
    rm -f "$(pid_file)"
    ok "APIARC stopped"
}

runtime_logs() {
    local lines="${1:-80}"
    local follow="${2:-0}"
    [ -f "$APIARC_SERVER_LOG" ] || die "no server log yet: $APIARC_SERVER_LOG"
    if [ "$follow" = "1" ]; then
        tail -n "$lines" -f "$APIARC_SERVER_LOG"
    else
        tail -n "$lines" "$APIARC_SERVER_LOG"
    fi
}

runtime_status() {
    banner
    echo
    if launcher_ready; then
        ok "launcher: ready"
    else
        warn "launcher: missing (run ./bootstrap.sh)"
    fi
    local prefix
    prefix=$(analysis_prefix)
    if [ -x "$prefix/bin/snakemake" ]; then
        ok "analysis: ready -> $prefix"
    else
        warn "analysis: missing (run scripts/build-analysis-env.sh)"
    fi
    if app_running; then
        ok "app: running (pid $(app_pid), port $(cat "$(port_file)" 2>/dev/null || echo '?'))"
    else
        warn "app: stopped"
    fi
    echo
    for tool in snakemake conda prefetch fasterq-dump pigz; do
        local found
        found=$(PATH="$prefix/bin:$PATH" command -v "$tool" 2>/dev/null || true)
        if [ -n "$found" ]; then
            printf '  %-14s %s\n' "$tool" "$found"
        else
            printf '  %-14s %s\n' "$tool" "-"
        fi
    done
}

runtime_shell() {
    local prefix
    prefix=$(analysis_prefix)
    log "entering APIARC shell (analysis prefix: $prefix)"
    (
        cd "$APIARC_REPO"
        export PATH="$prefix/bin:$PATH"
        export APIARC_REPO APIARC_RUNTIME_DIR APIARC_ANALYSIS_ENV="$prefix" APIARC_LAUNCH_ENV
        export CONDA_PKGS_DIRS CONDA_ENVS_DIRS PIP_CACHE_DIR
        exec "${SHELL:-/bin/bash}" -i
    )
}

runtime_update() {
    local scope="${1:-all}"
    case "$scope" in
        all|launcher) runtime_install ;;
    esac
    case "$scope" in
        all|analysis) provision_analysis ;;
    esac
    case "$scope" in
        all|launcher|analysis) ;;
        *) die "unknown scope '$scope'; use launcher|analysis|all" ;;
    esac
    ok "update finished ($scope)"
}
