#!/usr/bin/env bash
set -euo pipefail

SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export APIARC_REPO="${APIARC_REPO:-$SRC}"
source "$SRC/scripts/lib/common.sh"

# A synthetic build (scripts/mock-hg19.sh) points APIARC_ANNO_ROOT at a scratch
# directory so the dashboard never mistakes it for a downloaded reference.
ANNO_ROOT="${APIARC_ANNO_ROOT:-$APIARC_REPO/workflow/resources/anno}"
MINE=""
ASSEMBLY=""
THREADS="8"
FORCE=0
PROGRESS_FILE=""

usage() {
    cat <<EOF
usage: scripts/build-annotation.sh --mine rna|chip --assembly GRCh37|GRCh38|GRCm38 [--threads N] [--force] [--progress-file PATH]

  --mine rna|chip        rna: genome + GTF + gene BED + HISAT2 index; chip: bowtie2 index
  --assembly NAME        GRCh37 (hg19), GRCh38 (hg38) or GRCm38 (mm10)
  --threads N            build threads (default 8)
  --force                delete the existing index and rebuild
  --progress-file PATH   write stage/percent JSON for the dashboard; omit it and
                         a plain CLI run behaves exactly as before

File names, index prefixes and download URLs come from workflow/scripts/common/annotation.py,
the same registry the pipeline and the dashboard read.
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        --mine)          MINE="${2:-}"; shift 2 ;;
        --assembly)      ASSEMBLY="${2:-}"; shift 2 ;;
        --threads)       THREADS="${2:-}"; shift 2 ;;
        --progress-file) PROGRESS_FILE="${2:-}"; shift 2 ;;
        --force)         FORCE=1; shift ;;
        -h|--help)       usage; exit 0 ;;
        *) die "unknown option '$1'; see --help" ;;
    esac
done

if [ -z "$MINE" ] || [ -z "$ASSEMBLY" ]; then
    usage
    die "--mine and --assembly are required"
fi

case "$MINE" in
    rna|chip) ;;
    *) die "unsupported --mine '$MINE'; use rna or chip" ;;
esac

# ---- progress reporting -----------------------------------------------------
# The dashboard polls one JSON file instead of parsing the log: log lines are
# written for humans and change freely, while the step list below is the
# contract the page renders. Weights, not step counts, drive the bar because
# the index build owns most of the wall clock while "prepare" is instant.
case "$MINE" in
    rna)
        STEPS=(prepare genome gtf contigs gene_bed splice index)
        WEIGHTS=(1 22 8 3 3 3 60)
        ;;
    chip)
        STEPS=(prepare genome link index)
        WEIGHTS=(1 12 1 86)
        ;;
esac

STEP_NO=0
CUR_FRAC="0"
CUR_MSG=""
CUR_DETAIL=""
DONE_STATE="false"
DONE_MSG="构建完成"
ERR_TEXT=""
BUILD_PID=0
RUN_START=$(date '+%Y-%m-%d %H:%M:%S')

human_bytes() {
    awk -v n="${1:-0}" 'BEGIN{split("B KB MB GB TB PB",u," ");i=1;while(n>=1024&&i<6){n/=1024;i++}printf "%.1f %s",n,u[i]}'
}

json_escape() {
    printf '%s' "${1:-}" | sed -e 's/\\/\\\\/g' -e 's/"/\\"/g' -e 's/\t/ /g' | tr -d '\r\n'
}

json_list() {
    local sep="" item
    for item in "$@"; do
        printf '%s"%s"' "$sep" "$item"
        sep=","
    done
}

dir_bytes() {
    local out
    out=$(du -sb "$1" 2>/dev/null | cut -f1) || true
    printf '%s' "${out:-0}"
}

# Atomic replace: the page may read while we write, and half a JSON document
# reads as "no progress at all".
emit_progress() {
    [ -n "$PROGRESS_FILE" ] || return 0
    local index=-1 step=""
    if [ "$STEP_NO" -gt 0 ]; then
        index=$((STEP_NO - 1))
        step="${STEPS[$index]:-}"
    fi
    mkdir -p "$(dirname "$PROGRESS_FILE")" 2>/dev/null || true
    local tmp="$PROGRESS_FILE.tmp"
    {
        printf '{\n'
        printf '  "kind": "%s",\n' "$MINE"
        printf '  "assembly": "%s",\n' "$ASSEMBLY"
        printf '  "steps": [%s],\n' "$(json_list "${STEPS[@]}")"
        printf '  "weights": [%s],\n' "$(json_list "${WEIGHTS[@]}")"
        printf '  "step": "%s",\n' "$step"
        printf '  "step_index": %s,\n' "$index"
        printf '  "step_progress": %s,\n' "$CUR_FRAC"
        printf '  "message": "%s",\n' "$(json_escape "$CUR_MSG")"
        printf '  "detail": "%s",\n' "$(json_escape "$CUR_DETAIL")"
        printf '  "started_at": "%s",\n' "$RUN_START"
        printf '  "updated_at": "%s",\n' "$(date '+%Y-%m-%d %H:%M:%S')"
        printf '  "done": %s,\n' "$DONE_STATE"
        printf '  "error": "%s"\n' "$(json_escape "$ERR_TEXT")"
        printf '}\n'
    } >"$tmp" 2>/dev/null || return 0
    mv -f "$tmp" "$PROGRESS_FILE" 2>/dev/null || true
}

enter_step() {  # enter_step NAME MESSAGE
    STEP_NO=$((STEP_NO + 1))
    CUR_FRAC="0"
    CUR_MSG="${2:-}"
    CUR_DETAIL=""
    emit_progress
    log "step ${STEP_NO}/${#STEPS[@]}: ${1}"
}

report_step() {  # report_step FRACTION [DETAIL]
    CUR_FRAC="${1:-0}"
    CUR_DETAIL="${2:-}"
    emit_progress
}

finish_step() {
    CUR_FRAC="1"
    CUR_DETAIL=""
    emit_progress
}

on_exit() {
    local code=$?
    trap - EXIT
    if [ -n "$PROGRESS_FILE" ]; then
        if [ "$code" -eq 0 ]; then
            DONE_STATE="true"
            CUR_FRAC="1"
            CUR_MSG="$DONE_MSG"
            CUR_DETAIL=""
            STEP_NO=${#STEPS[@]}
            emit_progress
        else
            DONE_STATE="false"
            ERR_TEXT="退出码 $code，停在 ${STEPS[0]}"
            if [ "$STEP_NO" -gt 0 ]; then
                ERR_TEXT="退出码 $code，停在 ${STEPS[$((STEP_NO - 1))]}"
            fi
            emit_progress
        fi
    fi
    exit "$code"
}
trap on_exit EXIT

# Bytes on disk are the only honest signal an index build gives: neither
# hisat2-build nor bowtie2-build prints a percentage, but both write parts
# whose total lands near a fixed multiple of the FASTA.
index_expect_bytes() {  # index_expect_bytes FASTA_BYTES PERCENT
    awk -v s="${1:-0}" -v p="${2:-100}" 'BEGIN{printf "%d", s*p/100}'
}

watch_index() {  # watch_index DIR EXPECTED_BYTES
    local dir="$1" expected="$2"
    local got frac
    while kill -0 "$BUILD_PID" 2>/dev/null; do
        got=$(dir_bytes "$dir")
        frac=$(awk -v g="$got" -v e="$expected" 'BEGIN{r=(e>0)?g/e:0; if(r<0)r=0; if(r>0.95)r=0.95; printf "%.3f", r}')
        report_step "$frac" "$(human_bytes "$got") / 约 $(human_bytes "$expected")"
        sleep 3
    done
}

PREFIX_DIR=$(analysis_prefix)
export PATH="$PREFIX_DIR/bin:$PATH"
PYTHON="$PREFIX_DIR/bin/python"
[ -x "$PYTHON" ] || PYTHON=$(command -v python3) || die "no python found"

# One line of metadata per assembly, tab separated so the shell can read it.
annotation_meta() {
    "$PYTHON" - "$ASSEMBLY" <<'PY'
import os
import sys

sys.path.insert(0, os.path.join(os.environ["APIARC_REPO"], "workflow", "scripts"))
from common import annotation

try:
    info = annotation.single(sys.argv[1])
except ValueError as error:
    sys.exit(str(error))
print("\t".join([
    os.path.basename(info["genome"]),
    os.path.basename(info["gtf"]),
    info["index"],
    info["base_url"],
    info["release"],
]))
PY
}

META=$(annotation_meta) || die "unknown assembly '$ASSEMBLY'"
IFS=$'\t' read -r GENOME GTF IDX BASE_URL RELEASE <<<"$META"

RNA_DIR="$ANNO_ROOT/RNAseq_anno/$ASSEMBLY"
CHIP_DIR="$ANNO_ROOT/ChIPseq_anno/$ASSEMBLY"
BED="${GTF%.gtf}.bed"
enter_step prepare "准备目录与工具"
mkdir -p "$RNA_DIR" "$CHIP_DIR" "$APIARC_REPO/logs"
LOG_FILE="$APIARC_REPO/logs/annotation_${MINE}_${ASSEMBLY}.log"
export LOG_FILE

log "mine=$MINE assembly=$ASSEMBLY threads=$THREADS force=$FORCE"
log "release=$RELEASE"
log "tools=$PREFIX_DIR/bin"
log "target=$([ "$MINE" = "rna" ] && echo "$RNA_DIR" || echo "$CHIP_DIR")"
finish_step

download() {  # download URL TARGET LABEL
    local url="$1"
    local target="$2"
    local label="${3:-$(basename "$target")}"
    if [ -s "$target" ]; then
        ok "reuse $target"
        return 0
    fi
    log "download $url"
    # Content-Length is the only way to a real percentage; a server that hides
    # it (or a mirror that re-chunks) falls back to bytes and speed only.
    local total
    # A failed HEAD (offline mirror, redirect loop) must not kill the build:
    # the download itself still knows how to fail properly.
    total=$(curl -sIL --max-time 20 "$url" 2>/dev/null \
        | awk 'BEGIN{IGNORECASE=1} /^content-length:/{v=$2} END{gsub(/\r/,"",v); printf "%d", v+0}') || total=0
    [ -n "$total" ] || total=0
    if command -v wget >/dev/null 2>&1; then
        wget -q -c -O "$target" "$url" &
    else
        curl -fL -sS -C - -o "$target" "$url" &
    fi
    local pid=$!
    local start now elapsed got speed frac eta detail
    start=$(date +%s)
    while kill -0 "$pid" 2>/dev/null; do
        got=$(stat -c%s "$target" 2>/dev/null || echo 0)
        now=$(date +%s)
        elapsed=$((now - start))
        if [ "$elapsed" -le 0 ]; then elapsed=1; fi
        speed=$((got / elapsed))
        if [ "${total:-0}" -gt 0 ]; then
            frac=$(awk -v g="$got" -v t="$total" 'BEGIN{r=g/t; if(r<0)r=0; if(r>0.85)r=0.85; printf "%.3f", r}')
            eta=$(awk -v g="$got" -v t="$total" -v s="$speed" 'BEGIN{print (s>0)? int((t-g)/s) : -1}')
            detail="$label · $(human_bytes "$got") / $(human_bytes "$total") · $(human_bytes "$speed")/s"
            if [ "$eta" -ge 0 ]; then detail="$detail · 剩余 ${eta}s"; fi
        else
            frac="0"
            detail="$label · 已下载 $(human_bytes "$got") · $(human_bytes "$speed")/s"
        fi
        report_step "$frac" "$detail"
        sleep 2
    done
    if ! wait "$pid"; then
        rm -f "$target"
        die "download failed: $url"
    fi
    got=$(stat -c%s "$target" 2>/dev/null || echo 0)
    ok "downloaded $got bytes"
    report_step "0.85" "$label · 已下载 $(human_bytes "$got")"
}

ensure_genome() {
    local fasta="$RNA_DIR/$GENOME"
    enter_step genome "下载并解压基因组"
    if [ -s "$fasta" ]; then
        ok "reuse $fasta"
        finish_step
        return 0
    fi
    download "$BASE_URL/$GENOME.gz" "$fasta.gz" "下载 $GENOME"
    report_step "0.90" "解压 $GENOME"
    gunzip -c "$fasta.gz" > "$fasta"
    rm -f "$fasta.gz"
    ok "genome ready: $fasta ($(stat -c%s "$fasta") bytes)"
    finish_step
}

ensure_gtf() {
    local gtf="$RNA_DIR/$GTF"
    enter_step gtf "下载并解压注释"
    if [ -s "$gtf" ]; then
        ok "reuse $gtf"
        finish_step
        return 0
    fi
    download "$BASE_URL/$GTF.gz" "$gtf.gz" "下载 $GTF"
    report_step "0.90" "解压 $GTF"
    gunzip -c "$gtf.gz" > "$gtf"
    rm -f "$gtf.gz"
    ok "annotation ready: $gtf"
    finish_step
}

ensure_contigs() {
    # bowtie2/HISAT2 count rows by counting reads; a GTF that names chr1 while
    # the FASTA names 1 yields zero counts with no error, so align the naming.
    enter_step contigs "统一染色体命名"
    "$PYTHON" "$APIARC_REPO/workflow/scripts/common/align_gtf_contigs.py" \
        --fasta "$RNA_DIR/$GENOME" --gtf "$RNA_DIR/$GTF" >> "${LOG_FILE:-/dev/null}" 2>&1 \
        || die "GTF and genome contig names cannot be aligned; see the log above"
    finish_step
}

ensure_gene_bed() {
    # Same 7 columns the previous hand-made BED had: chr, start, end, gene_id,
    # score, strand, gene_name. Downstream R code merges DEG tables on gene_id.
    local bed="$RNA_DIR/$BED"
    enter_step gene_bed "从 GTF 生成基因 BED"
    if [ -s "$bed" ] && [ "$FORCE" != "1" ]; then
        ok "reuse $bed"
        finish_step
        return 0
    fi
    log "derive gene BED from $GTF"
    awk -F '\t' 'BEGIN { OFS = "\t" }
        $3 == "gene" {
            if (!match($9, /gene_id "[^"]+"/)) next
            gid = substr($9, RSTART + 9, RLENGTH - 10)
            name = gid
            if (match($9, /gene_name "[^"]+"/)) {
                name = substr($9, RSTART + 11, RLENGTH - 12)
            }
            print $1, $4 - 1, $5, gid, 0, $7, name
        }' "$RNA_DIR/$GTF" > "$bed.tmp"
    mv "$bed.tmp" "$bed"
    ok "gene BED ready: $bed ($(wc -l < "$bed") genes)"
    finish_step
}

build_rna() {
    ensure_genome
    ensure_gtf
    ensure_contigs
    ensure_gene_bed
    local index_dir="$RNA_DIR/$IDX"
    if [ "$FORCE" = "1" ]; then
        rm -rf "$index_dir"
    elif [ -s "$index_dir/genome.1.ht2" ] && [ -s "$index_dir/genome.8.ht2" ]; then
        ok "HISAT2 index already present: $index_dir"
        DONE_MSG="索引已存在，未重建"
        return 0
    fi
    command -v hisat2-build >/dev/null 2>&1 || die "hisat2-build not found in $PREFIX_DIR/bin"
    command -v extract_splice_sites.py >/dev/null 2>&1 || die "extract_splice_sites.py not found in $PREFIX_DIR/bin"
    mkdir -p "$index_dir"
    enter_step splice "提取剪接位点与外显子"
    log "extract splice sites and exons from $GTF"
    extract_splice_sites.py "$RNA_DIR/$GTF" > "$RNA_DIR/splicesites_$IDX.txt"
    extract_exons.py "$RNA_DIR/$GTF" > "$RNA_DIR/exons_$IDX.txt"
    finish_step
    enter_step index "构建 HISAT2 索引"
    log "hisat2-build -p $THREADS $GENOME -> $index_dir/genome"
    hisat2-build -p "$THREADS" \
        --ss "$RNA_DIR/splicesites_$IDX.txt" \
        --exon "$RNA_DIR/exons_$IDX.txt" \
        "$RNA_DIR/$GENOME" "$index_dir/genome" &
    BUILD_PID=$!
    # 190% of the FASTA, measured on the GRCm38 build: 5.2 GB of index for a
    # 2.8 GB genome. An estimate only feeds the bar; the build's own exit code
    # decides whether the step is finished.
    watch_index "$index_dir" "$(index_expect_bytes "$(stat -c%s "$RNA_DIR/$GENOME")" 190)"
    if ! wait "$BUILD_PID"; then
        die "hisat2-build failed; see the log above"
    fi
    ok "HISAT2 index ready: $index_dir/genome"
    finish_step
}

build_chip() {
    ensure_genome
    enter_step link "链接基因组 FASTA"
    local link="$CHIP_DIR/$GENOME"
    if [ ! -e "$link" ]; then
        ln -sf "../../RNAseq_anno/$ASSEMBLY/$GENOME" "$link"
        ok "linked genome: $link"
    fi
    finish_step
    enter_step index "构建 bowtie2 索引"
    if [ "$FORCE" = "1" ]; then
        rm -f "$CHIP_DIR/$IDX".*.bt2
    elif [ -s "$CHIP_DIR/$IDX.1.bt2" ] && [ -s "$CHIP_DIR/$IDX.rev.2.bt2" ]; then
        ok "bowtie2 index already present: $CHIP_DIR/$IDX"
        DONE_MSG="索引已存在，未重建"
        return 0
    fi
    command -v bowtie2-build >/dev/null 2>&1 || die "bowtie2-build not found in $PREFIX_DIR/bin"
    log "bowtie2-build --threads $THREADS $GENOME -> $CHIP_DIR/$IDX"
    bowtie2-build --threads "$THREADS" "$RNA_DIR/$GENOME" "$CHIP_DIR/$IDX" &
    BUILD_PID=$!
    # 136% of the FASTA, measured on the GRCm38 build: 3.8 GB for 2.8 GB.
    watch_index "$CHIP_DIR" "$(index_expect_bytes "$(stat -c%s "$RNA_DIR/$GENOME")" 136)"
    if ! wait "$BUILD_PID"; then
        die "bowtie2-build failed; see the log above"
    fi
    ok "bowtie2 index ready: $CHIP_DIR/$IDX"
    finish_step
}

case "$MINE" in
    rna)  build_rna ;;
    chip) build_chip ;;
esac

ok "annotation build finished: $MINE $ASSEMBLY"
