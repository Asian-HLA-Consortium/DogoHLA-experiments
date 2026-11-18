#!/usr/bin/env bash
# fetch_resources.sh: Download and verify standard resources for DogoHLA workflow.
# Public references: IPD-IMGT/HLA 3.65.0, GRCh38 decoy reference, T1K v1.0.6-r216, precomputed bait 31-mers.
set -euo pipefail

DEST_DIR="${1:-resources}"
mkdir -p "$DEST_DIR"

echo "=== Staging DogoHLA workflow resources in $DEST_DIR ==="

# Expected Checksums
EXPECTED_IPD_SHA256="9202c3d37a89045a8f01f08dfd3ba21c465835b72ac2e9e7a731ef4dfc2c1d6a"
EXPECTED_BAIT_SHA256="50b3b7324a5f5bca95272e69f83ba6cec592663363a79d549f551a209bb78170"

# 1. Bait 31-mers (888,049 canonical 31-mers)
BAIT_FILE="$DEST_DIR/bait_kmers.txt"
if [ -f "$BAIT_FILE" ]; then
    echo "Found bait_kmers.txt, checking checksum..."
    echo "$EXPECTED_BAIT_SHA256  $BAIT_FILE" | sha256sum -c - || {
        echo "Bait checksum mismatch! Re-downloading..."
        rm -f "$BAIT_FILE"
    }
fi

if [ ! -f "$BAIT_FILE" ]; then
    echo "Fetching bait 31-mer set..."
    # If on IBEX, copy from stage storage
    if [ -f "/encrypted/e3008/hohndor-dogohla-20261020/reads-repull-20261114/work/bait/bait_kmers.txt" ]; then
        cp "/encrypted/e3008/hohndor-dogohla-20261020/reads-repull-20261114/work/bait/bait_kmers.txt" "$BAIT_FILE"
    else
        echo "Please place verified bait_kmers.txt (sha256: $EXPECTED_BAIT_SHA256) into $DEST_DIR"
    fi
    echo "$EXPECTED_BAIT_SHA256  $BAIT_FILE" | sha256sum -c -
fi

# 2. IPD-IMGT/HLA 3.65.0 genomic.fa
IPD_FILE="$DEST_DIR/genomic.fa"
if [ -f "$IPD_FILE" ]; then
    echo "Found genomic.fa, checking checksum..."
    echo "$EXPECTED_IPD_SHA256  $IPD_FILE" | sha256sum -c - || {
        echo "IPD genomic checksum mismatch! Re-downloading..."
        rm -f "$IPD_FILE"
    }
fi

if [ ! -f "$IPD_FILE" ]; then
    echo "Fetching IPD-IMGT/HLA 3.65.0 genomic reference..."
    if [ -f "/ibex/scratch/projects/c2014/rob/dogohla-benchmark/hprc-zenodo-20261004/baselines/assets/references/genomic.fa" ]; then
        cp "/ibex/scratch/projects/c2014/rob/dogohla-benchmark/hprc-zenodo-20261004/baselines/assets/references/genomic.fa" "$IPD_FILE"
    else
        echo "Please place verified genomic.fa (sha256: $EXPECTED_IPD_SHA256) into $DEST_DIR"
    fi
    echo "$EXPECTED_IPD_SHA256  $IPD_FILE" | sha256sum -c -
fi

echo "=== All workflow resources successfully verified ==="
