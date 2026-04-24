#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
WORKFLOW_ANNO_DIR="$SCRIPT_DIR"

RNASEQ_DIR="$WORKFLOW_ANNO_DIR/RNAseq_anno"
CHIPSEQ_DIR="$WORKFLOW_ANNO_DIR/ChIPseq_anno"


TOTAL_THREADS=$(nproc)
THREADS=$(( TOTAL_THREADS * 3 / 10 ))


if [ $THREADS -lt 2 ]; then THREADS=2; fi
if [ $THREADS -gt 8 ]; then THREADS=8; fi

echo "=========================================================="
echo " Stable annotation download and build (CPU ~30%)"
echo "Working directory: $WORKFLOW_ANNO_DIR"
echo "Threads used: $THREADS"
echo "=========================================================="

mkdir -p "$RNASEQ_DIR" "$CHIPSEQ_DIR"


download() {
    local url="$1"
    local file="$2"
    if [ ! -f "$file" ]; then
        if command -v axel &>/dev/null; then
            axel -n 8 -o "$file" "$url"
        else
            wget -c "$url" -O "$file"
        fi
    fi
}

# ============================

# ============================
echo "[1/4] Processing Human GRCh38..."
mkdir -p "$RNASEQ_DIR/GRCh38"
cd "$RNASEQ_DIR/GRCh38"

fa="GRCh38.primary_assembly.genome.fa"
gtf="gencode.v44.annotation.gtf"

if [ ! -f "$fa" ]; then
    download "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz" "$fa.gz"
    gunzip "$fa.gz" &
fi

if [ ! -f "$gtf" ]; then
    download "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz" "$gtf.gz"
    gunzip "$gtf.gz" &
fi
wait


if [ ! -f grch38/genome.1.ht2 ]; then
    extract_splice_sites.py "$gtf" > splicesites_grch38.txt
    extract_exons.py "$gtf" > exons_grch38.txt
    mkdir -p grch38
    hisat2-build -p "$THREADS" --ss splicesites_grch38.txt --exon exons_grch38.txt "$fa" grch38/genome
fi


mkdir -p "$CHIPSEQ_DIR/GRCh38"
cd "$CHIPSEQ_DIR/GRCh38"
if [ ! -f grch38.1.bt2 ]; then
    ln -sf ../../RNAseq_anno/GRCh38/$fa ./
    bowtie2-build --threads "$THREADS" "$fa" grch38
fi

# ============================

# ============================
echo "[2/4] Processing Mouse GRCm38..."
mkdir -p "$RNASEQ_DIR/GRCm38"
cd "$RNASEQ_DIR/GRCm38"

fa="GRCm38.primary_assembly.genome.fa"
gtf="gencode.vM25.annotation.gtf"

if [ ! -f "$fa" ]; then
    download "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz" "$fa.gz"
    gunzip "$fa.gz" &
fi

if [ ! -f "$gtf" ]; then
    download "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz" "$gtf.gz"
    gunzip "$gtf.gz" &
fi
wait


if [ ! -f grcm38/genome.1.ht2 ]; then
    extract_splice_sites.py "$gtf" > splicesites_grcm38.txt
    extract_exons.py "$gtf" > exons_grcm38.txt
    mkdir -p grcm38
    hisat2-build -p "$THREADS" --ss splicesites_grcm38.txt --exon exons_grcm38.txt "$fa" grcm38/genome
fi


mkdir -p "$CHIPSEQ_DIR/GRCm38"
cd "$CHIPSEQ_DIR/GRCm38"
if [ ! -f grcm38.1.bt2 ]; then
    ln -sf ../../RNAseq_anno/GRCm38/$fa ./
    bowtie2-build --threads "$THREADS" "$fa" grcm38
fi

echo "=========================================================="
echo "Annotation download completed"
echo "=========================================================="