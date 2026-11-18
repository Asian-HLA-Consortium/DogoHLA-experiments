"""Development go/no-go gate (criteria fixed by the supervisor before pass-2 numbers were seen):
 (a) false MOS final calls on exact-IPD copies <= 1% of those copies
 (b) arm C recovers the exact sequence of >= 20% of the synthetic mosaic copies, and more than arm A
 (c) arm E recovers >= 50% of the synthetic mosaic copies
Reads results/score_rows_dev.tsv (score.py --tag dev). 'Recovers' = the score.py 'recovered'
definition (true sequence among the called candidates of the gene). Exact-IPD copies = IPD-present
copies of the development donors (SYN partners/exact cases and D-sim), all arms C.
"""
import collections
import csv
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C

import argparse
ap = argparse.ArgumentParser()
ap.add_argument('--tag', default='dev')
ap.add_argument('--syn', default='SYN')
a_ = ap.parse_args()
MOS = {'mosaic1_partner_source', 'mosaic1_partner_other', 'mosaic2'}
rows = list(csv.DictReader(open(C.S / f'results/score_rows_{a_.tag}.tsv'), delimiter='\t'))
T = lambda x: x == 'True'
exact_c = [r for r in rows if r['arm'] == 'C' and not T(r['ipd_absent'])]
fa = sum(T(r['false_mos']) for r in exact_c)
syn = {a: [r for r in rows if r['readset'] == a_.syn and r['arm'] == a and r['class_full'] in MOS] for a in 'ACE'}
rec = {a: sum(T(r['recovered']) for r in syn[a]) for a in 'ACE'}
n = len(syn['C'])
a_ok = fa <= 0.01 * len(exact_c)
b_ok = n > 0 and rec['C'] >= 0.2 * n and rec['C'] > rec['A']
c_ok = n > 0 and rec['E'] >= 0.5 * n
out = dict(false_mos_on_exact=fa, exact_copies=len(exact_c), syn_mosaic_copies=n, recovered=rec,
           a=a_ok, b=b_ok, c=c_ok, gate='PASS' if (a_ok and b_ok and c_ok) else 'FAIL',
           incomplete=sorted({r['donor'] for r in rows if r['arm_complete'] != 'True'}))
by = collections.Counter()
for r in exact_c:
    by[(r['readset'], r['gene'])] += T(r['false_mos'])
out['false_mos_on_exact_by_readset_gene'] = {f'{k[0]}:{k[1]}': v for k, v in by.items() if v}
per_gene = collections.defaultdict(lambda: [0, 0, 0, 0])
for a, i in (('C', 0), ('E', 1)):
    for r in syn[a]:
        per_gene[r['gene']][i] += T(r['recovered'])
for r in syn['C']:
    per_gene[r['gene']][2] += 1
out['syn_mosaic_recovery_per_gene_C_E_n'] = {g: v[:3] for g, v in per_gene.items()}
out['four_field_named_C_vs_A'] = {rs: {a: [sum(T(r['named4']) for r in rows if r['readset'] == rs and r['arm'] == a and r['truth4'] and not T(r['ipd_absent'])),
                                            sum(1 for r in rows if r['readset'] == rs and r['arm'] == a and r['truth4'] and not T(r['ipd_absent']))] for a in 'ACE'}
                                    for rs in sorted({r['readset'] for r in rows})}
out['dsim_absent_recovery'] = {a: [sum(T(r['recovered']) for r in rows if r['readset'] == 'D-sim' and r['arm'] == a and T(r['ipd_absent'])),
                                   sum(1 for r in rows if r['readset'] == 'D-sim' and r['arm'] == a and T(r['ipd_absent']))] for a in 'ACE'}
out['syn_final_call_recovery'] = {a: sum(T(r['recovered_final']) for r in syn[a]) for a in 'ACE'}
C.save_json(C.S / f'results/DEV_GATE_{a_.tag}.json', out)
print(out)
