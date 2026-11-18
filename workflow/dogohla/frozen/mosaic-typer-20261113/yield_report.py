"""YIELD.md from results/yield_rarefaction.tsv (Y1) and, if present, results/score_rows_test.tsv (Y2).
  python yield_report.py
Descriptive analyses added after the freeze at the user's request (DESIGN_ISSUES 14)."""
import collections
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C
import score as SC

CL = ('exact', 'mosaic', 'near-mosaic', 'private')


def q(v):
    v = np.asarray(v, dtype=float)
    return v.mean(), np.percentile(v, 2.5), np.percentile(v, 97.5)


def fmt(m):
    return f'{100 * m[0]:.1f}% [{100 * m[1]:.1f}-{100 * m[2]:.1f}]'


def y1():
    rows = list(csv.DictReader(open(C.S / 'results/yield_rarefaction.tsv'), delimiter='\t'))
    for r in rows:
        for k in CL + ('copies', 'donors'):
            r[k] = int(r[k])
    L = ['## Y1. Rarefaction of the truth set (panel = IPD + n cohort donors; no reads)', '',
         'Whole FOLDS components are added in random order until >= n donors (20 subsets per n, seed 20261113);',
         'only copies whose component is outside the subset are classified. Values: mean over subsets [2.5-97.5%].',
         'n = 0 is IPD alone (one panel); "all" = every other component for each copy (the atlas full panel).', '']
    curves = {}
    for scope in ('ipd_absent', 'all'):
        for gene in ('ALL', 'DRB1', 'nonDRB1'):
            L += [f'### mode any, scope {scope}, gene {gene}', '',
                  '| n target | donors (mean) | copies held out (mean) | exact (panel alone) | mosaic (reachable only by a proposer) | near-mosaic + private (neither) |',
                  '|---|---:|---:|---|---|---|']
            for nn in ['0', '10', '20', '40', '80', 'all']:
                sel = [r for r in rows if r['mode'] == 'any' and r['n_target'] == nn and r['scope'] == scope
                       and (r['gene'] == gene or (gene == 'nonDRB1' and r['gene'] not in ('ALL', 'DRB1')))]
                if not sel:
                    continue
                if gene == 'nonDRB1':
                    agg = collections.defaultdict(lambda: collections.Counter())
                    for r in sel:
                        agg[r['rep']].update({k: r[k] for k in CL + ('copies',)})
                        agg[r['rep']]['donors'] = r['donors']
                    sel = list(agg.values())
                fr = {k: [r[k] / r['copies'] for r in sel if r['copies']] for k in CL}
                neither = [(r['near-mosaic'] + r['private']) / r['copies'] for r in sel if r['copies']]
                curves[(scope, gene, nn)] = (q(fr['exact']), q(fr['mosaic']), q(neither))
                L.append(f"| {nn} | {np.mean([r['donors'] for r in sel]):.0f} | {np.mean([r['copies'] for r in sel]):.0f} | "
                         f"{fmt(q(fr['exact']))} | {fmt(q(fr['mosaic']))} | {fmt(q(neither))} |")
            L.append('')
    L += ['### Same population group vs other groups (scope ipd_absent, all genes)', '',
          '| query group | n target | same group: exact | same: mosaic | other groups: exact | other: mosaic |',
          '|---|---|---|---|---|---|']
    for g in sorted({r['group'] for r in rows if r['mode'] == 'same'}):
        for nn in ['10', '20', '40', '80']:
            def m(mode, k):
                sel = [r for r in rows if r['mode'] == mode and r['group'] == g and r['n_target'] == nn
                       and r['scope'] == 'ipd_absent' and r['gene'] == 'ALL' and r['copies']]
                return fmt(q([r[k] / r['copies'] for r in sel])) + f' (subsets {len(sel)}, donors {np.mean([r["donors"] for r in sel]):.0f})' if sel else 'n/a'
            L.append(f'| {g} | {nn} | {m("same", "exact")} | {m("same", "mosaic")} | {m("other", "exact")} | {m("other", "mosaic")} |')
    L.append('')
    summ = C.load_json(C.S / 'tmp/y1/summary.json')
    L += [f'Classifier check: Y1 class at the full panel vs atlas full_label_L150 (Y1|atlas: count): '
          f'{summ["full_vs_atlas"]}. Population groups (donors): {summ["groups"]}.', '']
    return L, curves


def y2():
    p = C.S / 'results/score_rows_test.tsv'
    if not p.exists():
        return ['## Y2. Realised yield', '', 'Pending: test scoring (job 55555798) has not finished.', ''], None
    rows = list(csv.DictReader(open(p), delimiter='\t'))
    T = lambda x: x == 'True'
    excl = C.exclusion_map()
    group_of = {d: tuple(sorted(v)) for d, v in excl.items()}
    idx = {(r['readset'], r['donor'], r['gene'], r['haplotype'], r['arm']): r for r in rows}
    L = ['## Y2. Realised yield on the test read sets (frozen arms; never pooled)', '']
    out = {}
    for rs in ('T-sim', 'T-real', 'D-real'):
        sel = [r for r in rows if r['readset'] == rs and T(r['ipd_absent'])]
        donors = len({r['donor'] for r in rows if r['readset'] == rs})
        keys = sorted({(r['donor'], r['gene'], r['haplotype']) for r in sel})
        rec = {a: {k: T(idx[(rs,) + k + (a,)]['recovered']) for k in keys if (rs,) + k + (a,) in idx}
               for a in 'ABCDE'}
        L += [f'### {rs}: {donors} donors, {len(keys)} IPD-absent copies', '',
              '| arm | recovered | per 100 donors |', '|---|---:|---:|']
        for a in 'ABCDE':
            n = sum(rec[a].values())
            L.append(f'| {a} | {n}/{len(rec[a])} | {100 * n / donors:.1f} |')

        def diff(a, b, filt=lambda k: True):
            pairs = [dict(donor=k[0], difference=int(rec[a][k]) - int(rec[b][k])) for k in keys
                     if k in rec[a] and k in rec[b] and filt(k)]
            u = SC.uncertainty(pairs, group_of)
            return u, sum(p['difference'] for p in pairs)
        L += ['', '| contrast | copies (net) | per 100 donors | clustered 95% CI (pp of IPD-absent copies) |', '|---|---:|---:|---|']
        for name, a, b in [('panel alone (B-A)', 'B', 'A'), ('mosaics alone (C-A)', 'C', 'A'),
                           ('both (D-A)', 'D', 'A'), ('mosaic plus on top of the panel (D-B)', 'D', 'B')]:
            u, net = diff(a, b)
            ci = u.get('bootstrap95_pp')
            L.append(f'| {name} | {net} | {100 * net / donors:.1f} | '
                     f'{u["difference_pp"]:+.2f} [{ci[0]:+.2f}, {ci[1]:+.2f}] |' if ci else f'| {name} | {net} | {100 * net / donors:.1f} | n/a |')
        both = sum(1 for k in keys if rec['B'].get(k) and rec['C'].get(k))
        L.append(f'\nRecovered by both B and C (overlap): {both}.')
        L += ['', '| realised fraction of the ceiling | DRB1 | other genes |', '|---|---|---|']
        for arm, cls_field, lab in [('C', 'class_ipd', 'arm C on IPD-only-mosaic copies'),
                                    ('D', 'class_full', 'arm D on full-panel-mosaic copies'),
                                    ('E', 'class_ipd', 'oracle E on IPD-only-mosaic copies')]:
            cells = []
            for gsel in (lambda g: g == 'DRB1', lambda g: g != 'DRB1'):
                ks = [k for k in keys if gsel(k[1]) and idx[(rs,) + k + ('A',)][cls_field] == 'mosaic']
                n = sum(rec[arm].get(k, False) for k in ks)
                cells.append(f'{n}/{len(ks)}')
            L.append(f'| {lab} | {cells[0]} | {cells[1]} |')
        L.append('')
        out[rs] = dict(donors=donors, rec={a: sum(rec[a].values()) for a in 'ABCDE'}, n=len(keys))
    return L, out


def main():
    L = ['# Yield: what adding people and recovering mosaics each deliver', '',
         'Added after the freeze at the user\'s request (DESIGN_ISSUES 14). Y1 uses truth sequences only (no reads,',
         'no typing): it gives ceilings. Y2 uses the frozen test arms: it gives realised recovery.', '']
    a, curves = y1()
    b, y2o = y2()
    L += a
    notes = C.S / 'results/YIELD_NOTES.md'
    if notes.exists():
        L += notes.read_text().splitlines() + ['']
    L += b
    (C.S / 'results/YIELD.md').write_text('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
