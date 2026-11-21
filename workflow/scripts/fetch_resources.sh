#!/usr/bin/env bash
# fetch_resources.sh: download and build the public reference resources of the DogoHLA workflow.
#
#   bash scripts/fetch_resources.sh [DEST]          (default DEST = resources)
#
# Run it inside the workflow image so that all tools are there (curl, unzip, perl, python3, samtools,
# jellyfish and T1K):
#   singularity exec dogohla-workflow_latest.sif bash scripts/fetch_resources.sh resources
#
# Result (the paths used in README.md and NIG.md):
#   DEST/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa (+ .fai)   GRCh38 with decoys and HLA contigs
#   DEST/ipd/genomic.fa                                             IPD-IMGT/HLA 3.65.0 genomic alleles in T1K format
#                                                                   (built from the public hla.dat, byte-identical to the pinned file)
#   DEST/bait/bait_kmers.txt                                        31-mer bait set (888,049 k-mers)
#   DEST/immuannot/Data-ipd3590-kir2130                             Immuannot reference data (IPD 3.59.0, KIR 2.13.0)
# Every file is checked against its expected sha256; a mismatch stops the script.
set -euo pipefail

DEST="${1:-resources}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG="$(cd "$HERE/.." && pwd)"
WORK="${TMPDIR:-$HOME/tmp}/dogohla-fetch.$$"
mkdir -p "$DEST"/ref "$DEST"/ipd "$DEST"/bait "$DEST"/immuannot "$WORK"
trap 'rm -rf "$WORK"' EXIT

# public sources
GRCH38_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
# IPD-IMGT/HLA release 3.65.0: the exact commit of the IMGTHLA repository that holds it
GRCH38_URL2="https://ftp-trace.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"   # same file, second mirror
IPD_DAT_URL="https://raw.githubusercontent.com/ANHIG/IMGTHLA/5b915f27f7f620361cf83cb626eeac8a03c0247c/hla.dat.zip"
IPD359_URL="https://raw.githubusercontent.com/ANHIG/IMGTHLA/3590/hla.dat.zip"
KIR_URL="https://raw.githubusercontent.com/ANHIG/IPDKIR/2130/kir.dat"

# expected checksums
EXPECTED_GRCH38_SHA256="3b103f4742abfd54938fb0333e19ad067635c8eb86f1dbf0ce44b165c4292b50"
EXPECTED_IPD_SHA256="9202c3d37a89045a8f01f08dfd3ba21c465835b72ac2e9e7a731ef4dfc2c1d6a"
EXPECTED_IPD_ZIP_SHA256="01aa70fed700cd3dc6e59a800ee2479582d0ac35e5e27df39b5db42ae54a1668"
EXPECTED_IPD_DAT_SHA256="817df2fd8e83036ac461208be4e04640f1a18dcdbb4a1bedebaaec26036bd583"
EXPECTED_PARSER_SHA256="8dc86688ebf55827679cbcccdbbe98849db3b12e3ffdc5f88f7ebc190f002983"   # T1K v1.0.6 ParseDatFile.pl
EXPECTED_GENOME_MODE_SHA256="9c19673640340e554942ca8549bb2bef8e0ccb3ad7bfd0511cc8191fd833f10a"
T1K_DIR="${T1K_DIR:-/opt/t1k}"
EXPECTED_BAIT_SHA256="50b3b7324a5f5bca95272e69f83ba6cec592663363a79d549f551a209bb78170"

check() {  # check FILE SHA256
    echo "$2  $1" | sha256sum -c - >/dev/null 2>&1
}

echo "=== Staging DogoHLA workflow resources in $DEST ==="

# 1. GRCh38 reference (CRAM decoding, bait construction)
REF="$DEST/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa"
if [ -f "$REF" ] && check "$REF" "$EXPECTED_GRCH38_SHA256"; then
    echo "GRCh38: found, checksum ok"
else
    echo "GRCh38: downloading"
    curl -fsSL --retry 5 -o "$REF" "$GRCH38_URL" || curl -fsSL --retry 5 -o "$REF" "$GRCH38_URL2"
    check "$REF" "$EXPECTED_GRCH38_SHA256" || { echo "GRCh38 checksum mismatch"; exit 1; }
fi
[ -f "$REF.fai" ] || samtools faidx "$REF"

# 2. IPD-IMGT/HLA 3.65.0 genomic alleles in T1K format, built from the public hla.dat
#    full genomic entries (T1K ParseDatFile.pl --mode genome) replace the truncated T1K entries (DNA mode)
#    of the same name; the 2,684 alleles without a full genomic sequence stay as they are; the record order
#    is the one of the pinned file (resources/ipd/genomic_order.txt.gz, allele names only).
IPD="$DEST/ipd/genomic.fa"
if [ -f "$IPD" ] && check "$IPD" "$EXPECTED_IPD_SHA256"; then
    echo "IPD genomic.fa: found, checksum ok"
else
    echo "IPD genomic.fa: building from IPD-IMGT/HLA 3.65.0"
    check "$T1K_DIR/ParseDatFile.pl" "$EXPECTED_PARSER_SHA256" || { echo "T1K ParseDatFile.pl (v1.0.6) not found in $T1K_DIR"; exit 1; }
    curl -fsSL --retry 5 -o "$WORK/hla.dat.zip" "$IPD_DAT_URL"
    check "$WORK/hla.dat.zip" "$EXPECTED_IPD_ZIP_SHA256" || { echo "hla.dat.zip checksum mismatch"; exit 1; }
    (cd "$WORK" && unzip -o -q hla.dat.zip)
    check "$WORK/hla.dat" "$EXPECTED_IPD_DAT_SHA256" || { echo "hla.dat checksum mismatch"; exit 1; }
    perl "$T1K_DIR/ParseDatFile.pl" "$WORK/hla.dat" --mode dna > "$WORK/ipd_dna_mode.fa"
    perl "$T1K_DIR/ParseDatFile.pl" "$WORK/hla.dat" --mode genome --gene HLA > "$WORK/ipd_genome_mode.fa" 2> "$WORK/parse.log"
    check "$WORK/ipd_genome_mode.fa" "$EXPECTED_GENOME_MODE_SHA256" || { echo "genome-mode parse differs from the expected one"; exit 1; }
    python3 "$HERE/build_genomic_reference.py" --baseline "$WORK/ipd_dna_mode.fa" --genomic "$WORK/ipd_genome_mode.fa" \
        --order "$PKG/resources/ipd/genomic_order.txt.gz" --out "$IPD"
    check "$IPD" "$EXPECTED_IPD_SHA256" || { echo "genomic.fa checksum mismatch"; exit 1; }
fi

# 3. bait 31-mers: IPD 31-mers that do not occur in GRCh38 outside the MHC
BAIT="$DEST/bait/bait_kmers.txt"
if [ -f "$BAIT" ] && check "$BAIT" "$EXPECTED_BAIT_SHA256"; then
    echo "bait_kmers.txt: found, checksum ok"
else
    echo "bait_kmers.txt: building"
    python3 "$HERE/build_bait.py" --ipd-genomic "$IPD" --grch38-ref "$REF" \
        --out-bait "$BAIT" --out-json "$DEST/bait/BAIT.json" --work-dir "$WORK/bait"
    check "$BAIT" "$EXPECTED_BAIT_SHA256" || { echo "bait checksum mismatch"; exit 1; }
fi

# 4. Immuannot reference data (same commands as containers/Dockerfile)
IMM="$DEST/immuannot/Data-ipd3590-kir2130"
if [ -f "$IMM/gen.fa.gz" ]; then
    echo "Immuannot data: found"
else
    echo "Immuannot data: downloading IPD 3.59.0 and KIR 2.13.0 and building"
    P="$PKG/dogohla/frozen/Immuannot/scripts.pub.v3/prepare-reference"
    curl -fsSL --retry 5 -o "$WORK/hla.dat.zip" "$IPD359_URL"
    (cd "$WORK" && unzip -o -q hla.dat.zip && mv hla.dat hla.3.59.0.dat)
    curl -fsSL --retry 5 -o "$WORK/kir.2.13.0.dat" "$KIR_URL"
    mkdir -p "$IMM"
    bash "$PKG/dogohla/frozen/Immuannot/scripts.pub.v3/predata.sh" "$WORK/hla.3.59.0.dat" "$WORK/kir.2.13.0.dat" "$P/c4a.gb" "$IMM"
    cat "$P/curated/DRB3.ele.csv" | gzip >> "$IMM/alleles.csv.gz"
    cat "$P/curated/DRB3.new.fa" | gzip >> "$IMM/gen.fa.gz"
    cat "$P/curated/DRB3.new.cds.fa" | gzip >> "$IMM/CDSseq/HLA-DRB3.fa.gz"
fi

echo "=== All workflow resources staged in $DEST ==="
