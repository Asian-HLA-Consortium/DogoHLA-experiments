"""Summarize $S/tmp/dev_false.json (counts only) into results/false_mosaics_summary.tsv and stdout."""
import collections
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C

D = C.load_json(C.S / 'tmp/dev_false.json')
cols = ['tag', 'readset', 'donor', 'gene', 'truth4_or_kind', 'candidate', 'truth_len', 'cand_len',
        'n_events', 'n_interior_events', 'interior_features', 'end_events_bp',
        'interior_reads_truth', 'interior_reads_cand', 'interior_cand_both_strands',
        'cand_reads_placed', 't1k_pass1', 't1k_pass2']
out = []
for tag in ('false', 'wrong'):
    for x in D[tag]:
        L = x['truth_len']
        inter = [e for e in x['events'] if not e['feature'].startswith('UTR')]
        endbp = sum(e['length'] for e in x['events'] if e['feature'].startswith('UTR'))
        sup = collections.Counter()
        placed = collections.Counter()
        for e in inter:
            for k, v in e['support'].items():
                a, b = k.split(':', 1)
                if a == 'cand_placed':
                    placed[b] += v
                else:
                    sup[(a, b)] += v
        out.append(dict(tag=tag, readset=x.get('readset', 'SYN'), donor=x['donor'], gene=x['gene'],
                        truth4_or_kind=x.get('truth4') or x.get('kind'), candidate=x['candidate'],
                        truth_len=L, cand_len=x['cand_len'], n_events=x['n_events'],
                        n_interior_events=len(inter),
                        interior_features=','.join(sorted({e['feature'] for e in inter})),
                        end_events_bp=endbp,
                        interior_reads_truth=sup[('truth', 'fwd')] + sup[('truth', 'rev')],
                        interior_reads_cand=sup[('cand', 'fwd')] + sup[('cand', 'rev')],
                        interior_cand_both_strands=sup[('cand', 'fwd')] > 0 and sup[('cand', 'rev')] > 0,
                        cand_reads_placed=';'.join(f'{k}:{v}' for k, v in placed.most_common()),
                        t1k_pass1=' | '.join(','.join(r[1].split(',')[:2]) + f' a={r[2]} q={r[3]}' for r in x['t1k_pass1']),
                        t1k_pass2=' | '.join(f'{r[1].split(",")[0]} a={r[2]} q={r[3]} / {r[4].split(",")[0]} a={r[5]} q={r[6]}' for r in x['t1k_pass2'])))
with open(C.S / 'results/false_mosaics_summary.tsv', 'w') as f:
    f.write('\t'.join(cols) + '\n')
    for r in out:
        f.write('\t'.join(str(r[c]) for c in cols) + '\n')
for r in out:
    print('\t'.join(str(r[c]) for c in cols))
print('DSIM_GAP')
for g in D['dsim_gap']:
    print(g)
