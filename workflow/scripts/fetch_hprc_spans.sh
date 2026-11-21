#!/usr/bin/env bash
# fetch_hprc_spans.sh: fetch the MHC spans of the public HPRC Asian assemblies and write the inputs of panel.cwl.
#
#   bash scripts/fetch_hprc_spans.sh [TABLE] [OUTDIR] [EXCLUDE_SAMPLES]
#
#   TABLE    resources/panel/hprc_asian_haplotypes.tsv (default): 86 donors, 172 haplotypes of the HPRC release 2
#            (donor, population, haplotype, public assembly URL, contig, span start and end)
#   OUTDIR   resources/panel (default); the spans go to OUTDIR/spans, the list of haplotypes to OUTDIR/panel_jobs.tsv
#   EXCLUDE_SAMPLES  comma-separated sample ids to leave out (for example the sample you are going to type)
#
# Run it inside the workflow image (it needs samtools; nothing else is installed on the host):
#   singularity exec $CWL_SINGULARITY_CACHE/ghcr.io_asian-hla-consortium_dogohla-workflow:latest.sif \
#       bash scripts/fetch_hprc_spans.sh
# samtools reads only the span (about 3.4 Mb per haplotype) from each 3 Gb assembly, through the .fai and .gzi
# index files that sit next to the public assembly. The span is the region from 100 kb before the first to
# 100 kb after the last of the 8 genes, from the HPRC Zenodo annotation (record 15285071).
# Next: bash scripts/run_panel.sh (annotation of the spans, gene extraction, panel).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG="$(cd "$HERE/.." && pwd)"
TABLE="${1:-$PKG/resources/panel/hprc_asian_haplotypes.tsv}"
OUT="$(mkdir -p "${2:-$PKG/resources/panel}" && cd "${2:-$PKG/resources/panel}" && pwd)"
EXCL=",${3:-},"
mkdir -p "$OUT/spans"
JOBS="${JOBS:-8}"          # parallel downloads
fetch_one() {              # donor hap url contig start end
    local f="$OUT/spans/$1.$2.fa"
    if [ ! -s "$f" ]; then
        samtools faidx "$3" "$4:$5-$6" > "$f.part" && mv "$f.part" "$f"
    fi
}
export -f fetch_one
export OUT
n=0
while IFS=$'\t' read -r donor pop hap url contig s e; do
    [ "$donor" = donor ] && continue
    case "$EXCL" in *",$donor,"*) continue;; esac
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$donor" "$hap" "$url" "$contig" "$s" "$e"
    n=$((n + 1))
done < "$TABLE" > "$OUT/spans.todo"
tr '\t' '\n' < "$OUT/spans.todo" | xargs -d '\n' -n 6 -P "$JOBS" bash -c 'fetch_one "$@"' _
rm -f "$OUT/spans.todo"
echo "$n spans in $OUT/spans"
python3 - "$TABLE" "$OUT" "$EXCL" <<'EOP'
import sys
from pathlib import Path
table, out, excl = sys.argv[1], Path(sys.argv[2]), sys.argv[3]
rows = [l.rstrip('\n').split('\t') for l in open(table)][1:]
rows = [r for r in rows if f',{r[0]},' not in excl]
with open(out / 'panel_jobs.tsv', 'w') as f:
    for r in rows:
        f.write(f'{r[0]}\t{r[2]}\t{out}/spans/{r[0]}.{r[2]}.fa\n')
print('wrote', out / 'panel_jobs.tsv', f'({len(rows)} haplotypes)')
EOP
