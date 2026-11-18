"""Part Y1 (added after the freeze at the user's request): rarefaction of the truth set.

Panel(S) = IPD genomic alleles + truth cores of the cohort donors in S, where S is a union of whole
FOLDS.json relatedness components (never splitting a family). Each held-out copy (its component not in
S) is classified exact / mosaic (L=150) / near-mosaic / private as in the atlas classifier:
  - the per-row reach of a query is computed once (vectorized; checked equal to the atlas
    _precompute_query_panel_matches on a sample), and every panel is a row subset of it;
  - k from the atlas solve_min_switches_with_min_length (imported, unchanged);
  - near-mosaic: not exact/mosaic, but mosaic/exact after masking 'private' columns. Here private
    columns = query positions where no panel row has the query's character, an approximation of the
    atlas's event-based private variants; agreement with the atlas at the full panel is reported.
Modes: any (components from all groups), same (components from the query's population group),
other (components from the other groups). n = target donors (whole components added in random
order until >= n); 20 subsets per (mode, group, n); seed 20261113. 'all' = for each query every
other component (one panel per query; the atlas full panel).

  python yield_y1.py --gene G   -> $S/tmp/y1/<G>.json
  python yield_y1.py --combine  -> results/yield_rarefaction.tsv (+ tmp/y1/summary.json)
"""
from __future__ import annotations

import argparse
import collections
import json
import random
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C

SEED = 20261113
NS = [0, 10, 20, 40, 80]
REPS = 20
EAS = {'CHS', 'CHB', 'JPT', 'KHV', 'CDX', 'CHN'}
SAS = {'PJL', 'GIH', 'ITU', 'STU', 'BEB'}


def group_of(meta):
    r = meta.get('resource')
    if r in ('Arab', 'Saudi'):
        return r
    p = meta.get('population')
    return 'HPRC_EAS' if p in EAS else 'HPRC_SAS' if p in SAS else 'HPRC_other'


def components():
    f = C.load_json(C.S / 'inputs/FOLDS.json')
    excl = C.exclusion_map()
    recs = C.load_records()
    truth_donors = sorted({r['donor'] for r in recs})
    comps = {}
    for d in truth_donors:
        comps.setdefault(tuple(sorted(excl[d])), []).append(d)
    comp_list = sorted(comps.items())
    grp = {d: group_of(f['donors'][d]) for d in truth_donors}
    return comp_list, grp, excl


def subsets():
    """Deterministic subsets: {(mode, group, n, rep): [component index]}."""
    comp_list, grp, _ = components()
    cgroups = [sorted({grp[d] for d in ds}) for _, ds in comp_list]
    groups = sorted(set(grp.values()))
    out = {}
    for n in NS:
        for rep in range(REPS):
            rng = random.Random(f'{SEED}-any-{n}-{rep}')
            out[('any', 'all', n, rep)] = pick(list(range(len(comp_list))), comp_list, n, rng)
            for g in groups:
                same = [i for i, cg in enumerate(cgroups) if cg == [g]]
                other = [i for i, cg in enumerate(cgroups) if g not in cg]
                rng = random.Random(f'{SEED}-same-{g}-{n}-{rep}')
                out[('same', g, n, rep)] = pick(same, comp_list, n, rng)
                rng = random.Random(f'{SEED}-other-{g}-{n}-{rep}')
                out[('other', g, n, rep)] = pick(other, comp_list, n, rng)
    return out


def pick(pool, comp_list, n, rng):
    pool = list(pool)
    rng.shuffle(pool)
    chosen, donors = [], 0
    for i in pool:
        if donors >= n:
            break
        chosen.append(i)
        donors += len(comp_list[i][1])
    return chosen


def reach_all(q_al: str, rows_al):
    """Vectorized equivalent of atlas _precompute_query_panel_matches (reach matrix only)."""
    q = np.frombuffer(q_al.encode(), dtype=np.uint8)
    nz = np.nonzero(q != ord('-'))[0]
    st, en = nz[0], nz[-1]
    bc = nz
    n = len(bc)
    P = np.vstack([np.frombuffer(s.encode(), dtype=np.uint8) for s in rows_al])
    isb = P != ord('-')
    first = np.where(isb.any(1), isb.argmax(1), -1)
    last = np.where(isb.any(1), P.shape[1] - 1 - isb[:, ::-1].argmax(1), -1)
    col = np.arange(P.shape[1])[None, :]
    valid = (first[:, None] != -1) & (col >= first[:, None]) & (col <= last[:, None])
    eq = (P == q[None, :]) & valid
    return bc, n, eq


def reach_from_eq(eq, bc):
    n = len(bc)
    mb = eq[:, bc]
    # step_ok[u]: all columns strictly between bc[u] and bc[u+1] equal
    bad = (~eq).astype(np.int32)
    pre = np.concatenate([np.zeros((eq.shape[0], 1), np.int32), np.cumsum(bad, axis=1)], axis=1)
    lo = bc[:-1] + 1
    hi = bc[1:]
    step_ok = (pre[:, hi] - pre[:, lo]) == 0
    brk = np.ones((eq.shape[0], n), dtype=bool)
    brk[:, :-1] = ~(step_ok & mb[:, 1:])
    idx = np.where(brk, np.arange(n)[None, :], n)
    nxt = np.minimum.accumulate(idx[:, ::-1], axis=1)[:, ::-1]
    return np.where(mb, nxt, -1).astype(np.int32)


def classify(reach, eq_sub, bc, solve):
    r = solve(reach, [''] * reach.shape[0], min_seg_bp=150)
    if r is not None:
        return 'exact' if r[0] == 0 else 'mosaic'
    if reach.shape[0] == 0:
        return 'private'
    priv = ~eq_sub.any(axis=0)
    priv_cols = np.nonzero(priv)[0]
    if len(priv_cols) == 0:
        return 'private'
    eq2 = eq_sub.copy()
    eq2[:, priv_cols] = True
    r2 = solve(reach_from_eq(eq2, bc), [''] * reach.shape[0], min_seg_bp=150)
    if r2 is not None and r2[0] >= 1:
        return 'near-mosaic'
    return 'private'


def per_gene(gene, chunk=None, nchunks=None):
    recs = C.load_records()
    classes = C.load_classes()
    msa = C.load_msas(recs, [gene])[gene]
    sys.path.insert(0, str(C.ATLAS / 'src'))
    mc = C.import_by_path('mosaic_classifier_atlas', C.ATLAS / 'src/mosaic_classifier.py')
    solve = mc.solve_min_switches_with_min_length
    comp_list, grp, excl = components()
    comp_of = {d: i for i, (_, ds) in enumerate(comp_list) for d in ds}
    subs = subsets()
    ipd_rows = [i for i, n in enumerate(msa.ipd_names) if n]
    coh_rows = [i for i, n in enumerate(msa.ipd_names) if not n]
    row_comps = {i: {comp_of[d] for d in msa.cohort_donors[i] if d in comp_of} for i in range(len(msa.uids))}
    out = dict(gene=gene, copies=[], checks=[])
    rng = random.Random(SEED)
    gene_recs = [x for x in recs if x['gene'] == gene]
    if chunk is not None:
        gene_recs = gene_recs[chunk::nchunks]
    for r in gene_recs:
        qrow = msa.row_of_source[r['source_id']]
        qa = msa.aligned[qrow]
        qc = comp_of[r['donor']]
        cand = ipd_rows + [i for i in coh_rows if row_comps[i] - {qc}]
        # rows carried ONLY by the query's own component never enter any panel
        bc, n, eq = reach_all(qa, [msa.aligned[i] for i in cand])
        reach = reach_from_eq(eq, bc)
        pos = {row: k for k, row in enumerate(cand)}
        ipd_idx = np.array([pos[i] for i in ipd_rows])
        rec = dict(source_id=r['source_id'], donor=r['donor'], group=grp[r['donor']], comp=qc,
                   ipd_absent=not r['genomic_labels'], atlas_full=classes[r['source_id']]['full_label_L150'],
                   atlas_ipd=classes[r['source_id']]['ipd_label_L150'], res={})
        c0 = classify(reach[ipd_idx], eq[ipd_idx], bc, solve)
        rec['ipd_only'] = c0
        if rng.random() < 0.05:   # validate the vectorized reach against the atlas
            _, ra, _, _ = mc._precompute_query_panel_matches(qa, ['x'] * len(ipd_rows[:200]),
                                                           [msa.aligned[i] for i in ipd_rows[:200]])
            out['checks'].append(bool(np.array_equal(ra, reach[ipd_idx[:200]])))
        full_idx = np.arange(len(cand))
        rec['full'] = classify(reach, eq, bc, solve) if c0 != 'exact' else 'exact'
        if c0 != 'exact':
            for key, chosen in subs.items():
                mode, g, nn, rep = key
                if qc in chosen or nn == 0:
                    continue
                if mode != 'any' and g != grp[r['donor']]:
                    continue
                cs = set(chosen)
                idx = np.array(list(ipd_idx) + [pos[i] for i in coh_rows if i in pos and row_comps[i] & cs])
                rec['res']['|'.join(map(str, key))] = classify(reach[idx], eq[idx], bc, solve)
        out['copies'].append(rec)
    C.save_json(C.S / 'tmp/y1' / (f'{gene}.json' if chunk is None else f'{gene}.part{chunk}.json'), out)
    print(gene, 'copies', len(out['copies']), 'reach checks equal', sum(out['checks']), '/', len(out['checks']))
    agree = collections.Counter((c['full'], c['atlas_full']) for c in out['copies'])
    print('full vs atlas', dict(agree))
    agree0 = collections.Counter((c['ipd_only'], c['atlas_ipd']) for c in out['copies'])
    print('ipd vs atlas', dict(agree0))


def combine():
    subs = subsets()
    comp_list, grp, _ = components()
    data = []
    for g in C.GENES:
        whole = C.S / 'tmp/y1' / f'{g}.json'
        parts = sorted((C.S / 'tmp/y1').glob(f'{g}.part*.json'))
        if whole.exists():
            data += C.load_json(whole)['copies']
        else:
            assert parts, 'missing Y1 output for ' + g
            for p_ in parts:
                data += C.load_json(p_)['copies']
    assert len(data) == 2351, len(data)
    gene_of = lambda c: c['source_id'].split('#HLA-')[1]
    rows = []
    CL = ('exact', 'mosaic', 'near-mosaic', 'private')

    def cls(c, key, mode, g):
        """Class of copy c under subset key, or None if c is not held out by it."""
        m, gg, nn, rep = key
        chosen = subs[key]
        if c['comp'] in chosen:
            return None
        if m != 'any' and gg != c['group']:
            return None
        if c['ipd_only'] == 'exact':
            return 'exact'
        if nn == 0:
            return c['ipd_only']
        return c['res'].get('|'.join(map(str, key)))
    for key in subs:
        mode, g, nn, rep = key
        donors = sum(len(comp_list[i][1]) for i in subs[key])
        if mode != 'any' and donors < nn and nn > 0:
            pass  # recorded as actual donors; small groups cannot reach n
        for scope in ('all', 'ipd_absent'):
            for gene in ('ALL',) + C.GENES:
                cnt = collections.Counter()
                for c in data:
                    if scope == 'ipd_absent' and not c['ipd_absent']:
                        continue
                    if gene != 'ALL' and gene_of(c) != gene:
                        continue
                    k = cls(c, key, mode, g)
                    if k:
                        cnt[k] += 1
                tot = sum(cnt.values())
                rows.append(dict(mode=mode, group=g, n_target=nn, rep=rep, donors=donors, scope=scope,
                                 gene=gene, copies=tot, **{k: cnt[k] for k in CL}))
    # 'all available' (atlas full panel per copy)
    for scope in ('all', 'ipd_absent'):
        for gene in ('ALL',) + C.GENES:
            cnt = collections.Counter()
            for c in data:
                if scope == 'ipd_absent' and not c['ipd_absent']:
                    continue
                if gene != 'ALL' and gene_of(c) != gene:
                    continue
                cnt[c['full']] += 1
            rows.append(dict(mode='any', group='all', n_target='all', rep=0,
                             donors=len({d for _, ds in comp_list for d in ds}) - 1, scope=scope, gene=gene,
                             copies=sum(cnt.values()), **{k: cnt[k] for k in CL}))
    cols = ['mode', 'group', 'n_target', 'rep', 'donors', 'scope', 'gene', 'copies'] + list(CL)
    with open(C.S / 'results/yield_rarefaction.tsv', 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows:
            f.write('\t'.join(str(r[k]) for k in cols) + '\n')
    agree = collections.Counter((c['full'], c['atlas_full']) for c in data)
    C.save_json(C.S / 'tmp/y1/summary.json', dict(full_vs_atlas={f'{a}|{b}': v for (a, b), v in agree.items()},
                                                 groups=collections.Counter(grp.values())))
    print('full vs atlas', dict(agree))
    print('groups', collections.Counter(grp.values()))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--gene')
    ap.add_argument('--combine', action='store_true')
    ap.add_argument('--chunk', type=int)
    ap.add_argument('--nchunks', type=int)
    a = ap.parse_args()
    C.verify_inputs()
    combine() if a.combine else per_gene(a.gene, a.chunk, a.nchunks)
