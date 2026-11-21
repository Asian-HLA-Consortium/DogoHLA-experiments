#!/usr/bin/env bash
# run_panel.sh: annotate the fetched spans, extract the 8 genes and build the panel (parts 1 and 2).
#
#   bash scripts/run_panel.sh [PANELDIR] [OUTDIR] [JOBS]
#
#   PANELDIR  folder written by fetch_hprc_spans.sh (resources/panel, default)
#   OUTDIR    result folder (out/panel, default): reference.fa, aliases.json, panel_summary.json
#   JOBS      haplotypes processed at the same time (4 by default; each uses 4 threads and about 4 GB)
#   RESOURCES (environment variable) the folder made by fetch_resources.sh (default ./resources)
#
# Run it on the host (not inside the image): it starts one `cwltool --singularity` process per haplotype, so that
# the processes do not share state (cwltool's own --parallel is not safe with Singularity), then one more for the panel.
# Needs CWL_SINGULARITY_CACHE (see README step 0). 172 haplotypes take about 4.5 minutes each.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG="$(cd "$HERE/.." && pwd)"
PANELDIR="$(cd "${1:-resources/panel}" && pwd)"
OUT="$(mkdir -p "${2:-out/panel}" && cd "${2:-out/panel}" && pwd)"
JOBS="${3:-4}"
RES="$(cd "${RESOURCES:-resources}" && pwd)"
export PKG OUT RES
mkdir -p "$OUT/part1"
run_one() {      # donor haplotype span
    local d="$OUT/part1/$1.$2"
    [ -s "$d/$1_$2_records.json" ] && return 0
    mkdir -p "$d"
    cat > "$d/inputs.yml" <<EOT
assembly: {class: File, path: "$3"}
immuannot_ref: {class: Directory, path: "$RES/immuannot/Data-ipd3590-kir2130"}
sample: "$1"
haplotype: "$2"
EOT
    cwltool --singularity --outdir "$d" "$PKG/part1.cwl" "$d/inputs.yml" > "$d/cwltool.log" 2>&1 || { echo "FAILED $1 $2 (see $d/cwltool.log)"; return 1; }
}
export -f run_one
xargs -d '\n' -n 1 -P "$JOBS" bash -c 'IFS=$'"'"'\t'"'"' read -r d h s <<< "$0"; run_one "$d" "$h" "$s"' < "$PANELDIR/panel_jobs.tsv"
n=$(ls "$OUT"/part1/*/*_records.json | wc -l)
echo "$n record files"
{
  echo "ipd_genomic: {class: File, path: \"$RES/ipd/genomic.fa\"}"
  echo "records:"
  ls "$OUT"/part1/*/*_records.json | while read -r f; do echo "  - {class: File, path: \"$f\"}"; done
  [ -n "${EXCLUDE_SAMPLES:-}" ] && echo "exclude_samples: \"$EXCLUDE_SAMPLES\""
} > "$OUT/build_panel_inputs.yml"
cwltool --singularity --outdir "$OUT" "$PKG/tools/build_panel.cwl" "$OUT/build_panel_inputs.yml"
echo "panel written to $OUT (reference.fa, aliases.json, panel_summary.json)"
