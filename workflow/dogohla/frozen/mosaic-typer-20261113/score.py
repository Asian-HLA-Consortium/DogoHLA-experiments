"""Scorer for stage mosaic-typer-20261113 (frozen with FREEZE.md before any test scoring).

  python score.py --tag dev  --readsets SYN,D-sim       --arms A,C,E
  python score.py --tag test --readsets T-sim,T-real,D-real --arms A,B,C,D,E

Per truth copy and arm:
  recovered          exact-sequence recovery: the copy's core sha256 equals the sequence of a novel
                     candidate (MOS_/PG/ORC_) present in ANY called slot of that gene (the
                     score_augmented.py 'called_cores' definition, so B and C/D are scored alike)
  recovered_final    ... and that candidate is the slot's FINAL call (IPD-prior rule: in a slot
                     holding both a named IPD allele and a novel candidate the named allele is the
                     call; rerank.py ipd_prior)
  named4             four-field named accuracy (IPD-present copies with one four-field label);
                     slots are assigned to copies by the best permutation; novel calls are not names
  named3/named2      three-/two-field, strict (novel = no name) and decoded (novel candidate named by
                     CDS / protein identity with an IPD genomic allele)
  final_mos          the copy's assigned slot's final call is a MOS_ candidate; false_mos when its
                     sequence is not the copy's truth sequence
Statistics: paired by copy, FOLDS-component-clustered bootstrap (4000, seed 20261020) and
family sign-flip test, copied from interpretation-20261020/score_ablation.py:uncertainty.
"""
from __future__ import annotations

import argparse
import collections
import itertools
import json
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C
import part0 as P0

NOVEL_TAGS = ('*MOS_', '*PG', '*ORC_')


def is_novel(label):
    return any(t in label for t in NOVEL_TAGS)


# ----------------------------------------------------------------------------- statistics

def uncertainty(pairs, group_of):
    """Verbatim logic of interpretation-20261020/score_ablation.py:uncertainty (pairs carry
    'donor' and 'difference')."""
    clusters = collections.defaultdict(list)
    ungrouped = sorted({r["donor"] for r in pairs if r["donor"] not in group_of})
    for r in pairs:
        clusters[group_of.get(r["donor"], (r["donor"],))].append(r["difference"])
    values = [(sum(v), len(v)) for v in clusters.values()]
    if not values:
        return {"family_clusters": 0, "difference_pp": None, "bootstrap95_pp": None,
                "sign_flip_p": None, "singleton_fallback_donors": ungrouped}
    rng = random.Random(20261020)
    boots = []
    for _ in range(4000):
        draw = [rng.choice(values) for _ in values]
        boots.append(100 * sum(x[0] for x in draw) / sum(x[1] for x in draw))
    boots.sort()
    differences = [v[0] for v in values if v[0]]
    observed = abs(sum(differences))
    if len(differences) <= 18:
        outcomes = [sum(s * d for s, d in zip(signs, differences))
                    for signs in itertools.product((-1, 1), repeat=len(differences))]
        pvalue = sum(abs(x) >= observed for x in outcomes) / len(outcomes)
        pmethod = "exact_family_sign_flip"
    else:
        hits = sum(abs(sum(rng.choice((-1, 1)) * d for d in differences)) >= observed
                   for _ in range(20000))
        pvalue = (hits + 1) / 20001
        pmethod = "monte_carlo_family_sign_flip_20000"
    return {"family_clusters": len(values),
            "singleton_fallback_donors": ungrouped,
            "difference_pp": 100 * sum(x[0] for x in values) / sum(x[1] for x in values),
            "bootstrap95_pp": [boots[100], boots[3899]],
            "sign_flip_p": pvalue, "sign_flip_method": pmethod}


# ----------------------------------------------------------------------------- truth

class Truth:
    def __init__(self, readsets):
        self.records = C.load_records()
        self.classes = C.load_classes()
        self.msas = C.load_msas(self.records)
        self.maps = {g: P0.ipd_name_maps(self.msas[g]) for g in C.GENES}
        excl = C.exclusion_map()
        self.group_of = {d: tuple(sorted(v)) for d, v in excl.items()}
        self.copies = collections.defaultdict(list)   # (readset, donor, gene) -> copies
        runs = C.readset_runs()
        dsim = set(runs['D-sim'])
        members = {'D-sim': sorted(dsim), 'T-sim': sorted(runs['T-sim']),
                   'T-real': sorted(d for d in runs['real'] if d not in dsim),
                   'D-real': sorted(d for d in runs['real'] if d in dsim)}
        by_donor = collections.defaultdict(list)
        for r in self.records:
            by_donor[r['donor']].append(r)
        for rs in readsets:
            if rs in C.SYN_SETS:
                for cid in C.load_json(C.S / 'work' / rs / 'CASES.json')['cases']:
                    case = C.load_json(C.syn_dir(cid) / 'CASE.json')
                    tr = C.load_json(C.syn_dir(cid) / 'TRUTH_SEQ.json')
                    g = case['gene']
                    for h in sorted(tr):
                        core = tr[h]['core']
                        kind = 'exact' if (h == '2' or case['kind'] == 'exact') else case['kind']
                        labels = []
                        if core in self.msas[g].row_of_seq and self.msas[g].ipd_names[self.msas[g].row_of_seq[core]]:
                            labels = self.msas[g].ipd_names[self.msas[g].row_of_seq[core]]
                        self.copies[(rs, cid, g)].append(self._copy(
                            rs, cid, g, h, core, tr[h]['exons'], labels, kind, kind,
                            cohort='SYN'))
                    self.group_of[cid] = (cid,)
                continue
            for d in members[rs]:
                for r in by_donor[d]:
                    cl = self.classes[r['source_id']]
                    self.copies[(rs, d, r['gene'])].append(self._copy(
                        rs, d, r['gene'], r['haplotype'], r['core_sequence'], r['core_exons'],
                        r['genomic_labels'], cl['full_label_L150'], cl['ipd_label_L150'],
                        cohort=r['resource']))

    def _copy(self, rs, d, g, h, core, exons, labels, cfull, cipd, cohort):
        cds_names, prot_names = self.maps[g]
        cds = C.cds_of(core, exons)
        t4 = sorted({C.normalized(a, 4) for a in labels} - {None})
        return dict(readset=rs, donor=d, gene=g, haplotype=h, cohort=cohort,
                    sha=C.sha256_str(core), ipd_absent=not labels,
                    truth4=t4[0] if len(t4) == 1 else None,
                    truth3=set(cds_names.get(cds, ())), truth2=set(prot_names.get(C.translate(cds), ())),
                    class_full=cfull, class_ipd=cipd)


# ----------------------------------------------------------------------------- calls

class NovelIndex:
    """label -> (sha256, cds-3field names, protein-2field names) for novel candidates."""

    def __init__(self, maps):
        self.maps = maps
        self.cache = {}

    def add_fasta(self, path):
        if str(path) in self.cache:
            return
        self.cache[str(path)] = True
        for name, ints, seq in C.read_fasta_headers_and_seqs(path):
            if not is_novel(name):
                continue
            g = C.gene_of_name(name)
            if g not in C.GENES:
                continue
            cds = C.cds_of(seq, C.header_exons(ints))
            cn, pn = self.maps[g]
            self.cache[name] = (C.sha256_str(seq), set(cn.get(cds, ())), set(pn.get(C.translate(cds), ())))

    def get(self, label):
        return self.cache.get(label)


def arm_table(rs, d, arm):
    """(calls table path, [reference/addition fasta paths for novel labels]) or (None, reason)."""
    if arm == 'A':
        if rs in C.SYN_SETS:
            p = C.syn_dir(d) / 'A'
            ok = (p / 'COMPLETE.json').exists()
            return (p / 'calls_genotype.tsv', []) if ok else (None, 'missing')
        base = {'D-sim': C.HPRC_CONTROL / 't1k/runs', 'T-sim': C.REGIONAL / 't1k/runs',
                'T-real': C.ARMS_REAL / 'runs', 'D-real': C.ARMS_REAL / 'runs'}[rs]
        return C.verified_table(base / d / 'genomic'), []
    if arm == 'B':
        base = {'T-sim': C.AUGMENTED / 'runs', 'T-real': C.ARMS_REAL / 'runs',
                'D-real': C.ARMS_REAL / 'runs'}.get(rs)
        if base is None:
            return None, 'no_arm_B_for_readset'
        run = base / d / 'genomic_core_hprc_asian'
        rec = C.load_json(run / 'COMPLETE.json')
        import run_donor as RD
        ref = RD.to_host(RD.parse_cmd(rec['command'])['-f'])
        return C.verified_table(run), [ref]
    w = (C.syn_dir(d) / arm) if rs in C.SYN_SETS else (C.S / 'work' / rs / d / arm)
    if not (w / 'COMPLETE.json').exists():
        return None, 'failed' if (w / 'FAILURE.json').exists() else 'missing'
    rec = C.load_json(w / 'COMPLETE.json')
    t = Path(rec['calls_table'])
    assert C.sha256_file(t) == rec['calls_table_sha256'], 'table changed: ' + str(t)
    fas = [w / 'additions.fa']
    if arm == 'D':
        fas.append(Path(rec['pass1_reference']))
    return t, fas


ACCEPTANCE = True


def apply_acceptance(rs, d, arm, g, calls):
    """Round-3 final call rule: a MOS final call rejected by accept.py reverts to the displaced
    named allele; the MOS label stays in the slot only as a detected sequence (recovery)."""
    w = (C.syn_dir(d) / arm) if rs in C.SYN_SETS else (C.S / 'work' / rs / d / arm)
    a = w / 'ACCEPT.json'
    if not a.exists():
        return calls
    rej = {x['slot']: x for x in C.load_json(a)['genes'].get(g, []) if not x['accepted']}
    if not rej:
        return calls
    calls = dict(calls)
    slots = [list(s) for s in calls.get(g, [[], []])]
    for i, x in rej.items():
        if i < len(slots):
            slots[i] = [x['displaced']] + slots[i]   # named + novel -> IPD-prior makes the named the call
    calls[g] = slots
    return calls


def typed(slot):
    named = [l for l in slot if not is_novel(l)]
    novel = [l for l in slot if is_novel(l)]
    return named if (named and novel) else list(slot)


def slot_eval(slot, novel_idx):
    t = typed(slot)
    out = dict(labels=slot, typed=t, final_novel=bool(t) and all(is_novel(l) for l in t),
               final_mos=bool(t) and all('*MOS_' in l for l in t))
    out['v4'] = P0.slot_value([l for l in t if not is_novel(l)], 4) if not any(is_novel(l) for l in t) else None
    for k in (2, 3):
        out[f'v{k}_strict'] = P0.slot_value(t, k) if not any(is_novel(l) for l in t) else None
        vals = set()
        for l in t:
            if is_novel(l):
                n = novel_idx.get(l)
                names = (n[1] if k == 3 else n[2]) if n else set()
                vals.add(next(iter(names)) if len(names) == 1 else None)
            else:
                vals.add(C.normalized(l, k))
        out[f'v{k}_decoded'] = next(iter(vals)) if len(vals) == 1 and None not in vals else None
    out['final_shas'] = {novel_idx.get(l)[0] for l in t if is_novel(l) and novel_idx.get(l)}
    out['all_shas'] = {novel_idx.get(l)[0] for l in slot if is_novel(l) and novel_idx.get(l)}
    out['unknown_novel'] = [l for l in slot if is_novel(l) and not novel_idx.get(l)]
    return out


def copy_eval(c, s):
    return dict(exact=c['sha'] in s['final_shas'],
                named4=c['truth4'] is not None and s['v4'] == c['truth4'],
                named3_strict=s['v3_strict'] in c['truth3'], named3_decoded=s['v3_decoded'] in c['truth3'],
                named2_strict=s['v2_strict'] in c['truth2'], named2_decoded=s['v2_decoded'] in c['truth2'])


def score_gene(copies, slots, novel_idx):
    ev = [slot_eval(s, novel_idx) for s in slots]
    gene_shas = set().union(*(e['all_shas'] for e in ev)) if ev else set()
    best = None
    for perm in itertools.permutations(range(len(ev)), len(copies)):
        ms = [copy_eval(c, ev[p]) for c, p in zip(copies, perm)]
        sc = sum(1000 * m['exact'] + 100 * m['named4'] + 10 * m['named3_decoded'] + m['named2_decoded']
                 for m in ms)
        if best is None or sc > best[0]:
            best = (sc, perm, ms)
    rows = []
    for c, p, m in zip(copies, best[1], best[2]):
        s = ev[p]
        rows.append(dict(m, recovered=c['sha'] in gene_shas, recovered_final=m['exact'],
                         final_mos=s['final_mos'], false_mos=s['final_mos'] and c['sha'] not in s['final_shas'],
                         detected_mos_any_slot=any('*MOS_' in l for e in ev for l in e['labels']),
                         no_call=not s['labels'], unknown_novel=len(s['unknown_novel']),
                         call='|'.join(s['typed']) if s['typed'] else '.'))
    return rows


def score(readsets, arms, tag):
    T = Truth(readsets)
    nidx = NovelIndex(T.maps)
    rows, availability = [], []
    for (rs, d, g), copies in sorted(T.copies.items()):
        copies = sorted(copies, key=lambda c: c['haplotype'])
        for arm in arms:
            if arm == 'B' and rs in ('D-sim',) + C.SYN_SETS:
                continue
            if arm == 'D' and rs in ('D-sim',) + C.SYN_SETS:
                continue
            key = (rs, d, arm)
            if key not in nidx.cache:
                try:
                    t, fas = arm_table(rs, d, arm)
                except Exception as e:
                    t, fas = None, repr(e)
                if t is None:
                    nidx.cache[key] = None
                    availability.append(dict(readset=rs, donor=d, arm=arm, state=fas))
                else:
                    for fa in fas:
                        if Path(fa).exists():
                            nidx.add_fasta(fa)
                    nidx.cache[key] = C.read_t1k(t)
                    availability.append(dict(readset=rs, donor=d, arm=arm, state='complete'))
            calls = nidx.cache[key]
            complete = calls is not None
            if complete and arm in ('C', 'D') and ACCEPTANCE:
                calls = apply_acceptance(rs, d, arm, g, calls)
            slots = calls.get(g, [[], []]) if complete else [[], []]
            for c, r in zip(copies, score_gene(copies, slots, nidx)):
                rows.append(dict(readset=rs, donor=d, gene=g, haplotype=c['haplotype'], arm=arm,
                                 cohort=c['cohort'], ipd_absent=c['ipd_absent'],
                                 truth4=c['truth4'] or '', has_truth3=bool(c['truth3']),
                                 has_truth2=bool(c['truth2']), class_full=c['class_full'],
                                 class_ipd=c['class_ipd'], arm_complete=complete, **r))
    return T, rows, availability


def summarize(T, rows, arms, readsets):
    idx = {(r['readset'], r['donor'], r['gene'], r['haplotype'], r['arm']): r for r in rows}
    out = collections.OrderedDict()

    def frac(sel, field):
        return dict(n=len(sel), correct=sum(bool(r[field]) for r in sel))

    def paired(rs, arm, base, field, filt):
        pairs = []
        for r in rows:
            if r['readset'] != rs or r['arm'] != arm or not filt(r):
                continue
            b = idx.get((rs, r['donor'], r['gene'], r['haplotype'], base))
            if b is None:
                continue
            pairs.append(dict(donor=r['donor'], difference=int(bool(r[field])) - int(bool(b[field])),
                              gene=r['gene']))
        u = uncertainty(pairs, T.group_of)
        u.update(n=len(pairs), gains=sum(p['difference'] > 0 for p in pairs),
                 losses=sum(p['difference'] < 0 for p in pairs))
        return u

    absent = lambda r: r['ipd_absent']
    reach = lambda r: r['ipd_absent'] and r['class_ipd'] == 'mosaic'
    present4 = lambda r: (not r['ipd_absent']) and bool(r['truth4'])
    for rs in readsets:
        S = out.setdefault(rs, collections.OrderedDict())
        for arm in arms:
            sel = [r for r in rows if r['readset'] == rs and r['arm'] == arm]
            if not sel:
                continue
            S[arm] = dict(
                recovery_all_absent=frac([r for r in sel if absent(r)], 'recovered'),
                recovery_final_all_absent=frac([r for r in sel if absent(r)], 'recovered_final'),
                recovery_ipd_only_mosaic=frac([r for r in sel if reach(r)], 'recovered'),
                named4_present=frac([r for r in sel if present4(r)], 'named4'),
                named3_strict=frac([r for r in sel if r['has_truth3']], 'named3_strict'),
                named3_decoded=frac([r for r in sel if r['has_truth3']], 'named3_decoded'),
                named2_strict=frac([r for r in sel if r['has_truth2']], 'named2_strict'),
                named2_decoded=frac([r for r in sel if r['has_truth2']], 'named2_decoded'),
                final_mos=sum(r['final_mos'] for r in sel),
                false_mos=sum(r['false_mos'] for r in sel),
                false_mos_by_class={k: sum(r['false_mos'] for r in sel if r['class_full'] == k)
                                    for k in ('exact', 'mosaic', 'near-mosaic', 'private',
                                              'mosaic1_partner_source', 'mosaic1_partner_other',
                                              'mosaic2', 'near_mosaic')},
                copies_by_class={k: sum(1 for r in sel if r['class_full'] == k)
                                 for k in sorted({r['class_full'] for r in sel})},
                incomplete_donors=sorted({r['donor'] for r in sel if not r['arm_complete']}),
                per_gene={g: dict(recovery=frac([r for r in sel if absent(r) and r['gene'] == g], 'recovered'),
                                  named4=frac([r for r in sel if present4(r) and r['gene'] == g], 'named4'),
                                  false_mos=sum(r['false_mos'] for r in sel if r['gene'] == g))
                          for g in C.GENES})
        cmp = S.setdefault('contrasts', collections.OrderedDict())
        for arm, base in [('C', 'A'), ('D', 'A'), ('B', 'A'), ('C', 'B'), ('D', 'B'), ('E', 'A'), ('C', 'E')]:
            if arm not in S or base not in S:
                continue
            cmp[f'{arm}-{base}'] = dict(
                recovery_all_absent=paired(rs, arm, base, 'recovered', absent),
                recovery_ipd_only_mosaic=paired(rs, arm, base, 'recovered', reach),
                named4_present=paired(rs, arm, base, 'named4', present4),
                named2_decoded=paired(rs, arm, base, 'named2_decoded', lambda r: r['has_truth2']),
                named3_decoded=paired(rs, arm, base, 'named3_decoded', lambda r: r['has_truth3']))
    return out


def write(tag, T, rows, avail, summ):
    res = C.S / 'results'
    cols = ['readset', 'donor', 'gene', 'haplotype', 'arm', 'cohort', 'ipd_absent', 'truth4',
            'class_full', 'class_ipd', 'arm_complete', 'recovered', 'recovered_final', 'named4',
            'named3_strict', 'named3_decoded', 'named2_strict', 'named2_decoded', 'final_mos',
            'false_mos', 'detected_mos_any_slot', 'no_call', 'unknown_novel', 'call']
    with open(res / f'score_rows_{tag}.tsv', 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows:
            f.write('\t'.join(str(r.get(c, '')) for c in cols) + '\n')
    lost = [r for r in rows if r['arm'] in ('C', 'D')]
    C.save_json(res / f'SCORES_{tag}.json', dict(summary=summ, availability=avail,
                code_sha256={p.name: C.sha256_file(p) for p in sorted((C.S / 'src').glob('*.py'))}))
    print(json.dumps(summ, indent=1, default=str)[:20000])


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--tag', required=True)
    ap.add_argument('--readsets', required=True)
    ap.add_argument('--arms', required=True)
    a = ap.parse_args()
    C.verify_inputs()
    rss, arms = a.readsets.split(','), a.arms.split(',')
    T, rows, avail = score(rss, arms, a.tag)
    summ = summarize(T, rows, arms, rss)
    write(a.tag, T, rows, avail, summ)
