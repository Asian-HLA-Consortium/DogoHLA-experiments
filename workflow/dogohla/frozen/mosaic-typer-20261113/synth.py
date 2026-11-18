"""Synthetic diploid development donors (Part 2 development only).

Built as in atlas Part C (splice two truth haplotypes at columns where they agree, segments
>= 150 bp), but parents and partners are drawn ONLY from the 20 D-sim donors' IPD-present copies,
so no test donor's haplotype is used for development. Each case is one gene locus (core + the
2,000 bp observed flanks of the first/last source), simulated with the pinned Phase 2 simulator
(art_modern, NovaSeq X profile, 150 bp PE, 15x per haplotype).

  python synth.py --build          writes $S/work/SYN/<case>/{CASE.json,TRUTH_SEQ.json},
                                   $S/work/SYN/CASES.json and results/synthetic_cases.tsv
Kinds per gene: 2 x mosaic1 with partner = a source, 2 x mosaic1 with partner = third allele,
2 x mosaic2 (partner third allele), 2 x exact (two real copies), 2 x near-mosaic (mosaic1 + one
private exon-2 SNV, partner third allele).
"""
from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C
import proposer as PR

SEED = 20261113
P2LIB = Path('/encrypted/e3008/hohndor-dogohla-20261020/mosaic-read-mapping-20261104/src/phase2_lib.py')
KINDS = ['mosaic1_partner_source'] * 2 + ['mosaic1_partner_other'] * 2 + ['mosaic2'] * 2 + \
        ['exact'] * 2 + ['near_mosaic'] * 2


def dsim_donors():
    return sorted(p.name for p in (C.HPRC_CONTROL / 't1k/runs').iterdir())


def build(seed=SEED, setname='SYN'):
    records = C.load_records()
    msas = C.load_msas(records)
    excl = C.exclusion_map()
    dsim = set(dsim_donors())
    rows_out = []
    cases = []
    for g in C.GENES:
        msa = msas[g]
        ipd_seqs = {msa.seqs[i] for i, n in enumerate(msa.ipd_names) if n}
        pool = [r for r in records if r['gene'] == g and r['donor'] in dsim and r['genomic_labels']
                and r['core_sequence'] in ipd_seqs]
        for k, kind in enumerate(KINDS):
            rng = random.Random(f'{seed}-{g}-{k}')
            case = None
            for _ in range(500):
                case = try_case(g, kind, pool, msa, excl, ipd_seqs, rng)
                if case:
                    break
            cid = f'{setname}-{g}-{k:02d}'
            if not case:
                rows_out.append(dict(case=cid, gene=g, kind=kind, status='could_not_build'))
                continue
            truth, meta = case
            meta.update(case_id=cid, gene=g, kind=kind,
                        truth_sha={h: C.sha256_str(truth[h]['core']) for h in truth},
                        ipd_present_shas=sorted({C.sha256_str(truth[h]['core']) for h in truth
                                                 if truth[h]['core'] in ipd_seqs}),
                        seeds={'1': rng.randrange(2 ** 31), '2': rng.randrange(2 ** 31)})
            d = C.syn_dir(cid)
            d.mkdir(parents=True, exist_ok=True)
            C.save_json(d / 'CASE.json', meta)
            C.save_json(d / 'TRUTH_SEQ.json', truth)
            cases.append(cid)
            rows_out.append(dict(case=cid, gene=g, kind=kind, status='built',
                                 sources=';'.join(meta['sources']),
                                 switch_columns=';'.join(map(str, meta['switch_columns'])),
                                 partner=meta['partner'], hap1_in_ipd=truth['1']['core'] in ipd_seqs,
                                 hap2_in_ipd=truth['2']['core'] in ipd_seqs,
                                 components=';'.join(meta['component_donors'])))
    C.save_json(C.S / 'work' / setname / 'CASES.json', dict(cases=cases, seed=seed))
    cols = ['case', 'gene', 'kind', 'status', 'sources', 'switch_columns', 'partner', 'hap1_in_ipd',
            'hap2_in_ipd', 'components']
    with open(C.S / f'results/synthetic_cases{"" if setname == "SYN" else "_" + setname}.tsv', 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows_out:
            f.write('\t'.join(str(r.get(c, '')) for c in cols) + '\n')
    print(len(cases), 'cases built of', len(rows_out))


def comp_ok(donors, excl):
    for i, a in enumerate(donors):
        for b in donors[i + 1:]:
            if b in excl[a]:
                return False
    return True


def try_case(g, kind, pool, msa, excl, ipd_seqs, rng):
    if kind == 'exact':
        a, b = rng.sample(pool, 2)
        if a['core_sequence'] == b['core_sequence'] or not comp_ok([a['donor'], b['donor']], excl):
            return None
        truth = {'1': hap(a, a, a['core_sequence'], a['core_exons']),
                 '2': hap(b, b, b['core_sequence'], b['core_exons'])}
        return truth, dict(sources=[a['source_id']], switch_columns=[], partner=b['source_id'],
                           component_donors=sorted({a['donor'], b['donor']}))
    a, b = rng.sample(pool, 2)
    if a['core_sequence'] == b['core_sequence']:
        return None
    ra, rb = msa.row_of_source[a['source_id']], msa.row_of_source[b['source_id']]
    panel = PR.Panel(['a', 'b'], [msa.aligned[ra], msa.aligned[rb]])
    P, Pb, Cn = panel.P, panel.Pbase, panel.C
    diff = (P[0] != P[1])
    agree = [s for s in range(1, Cn) if P[0, s] == P[1, s] and P[0, s] < 4]
    nsw = 2 if kind == 'mosaic2' else 1
    if len(agree) < 10:
        return None
    for _ in range(50):
        sw = sorted(rng.sample(agree, nsw))
        bounds = [0] + sw + [Cn]
        segs = [(i % 2, bounds[i], bounds[i + 1]) for i in range(len(bounds) - 1)]
        if not all(Pb[r, e] - Pb[r, s] >= 150 for r, s, e in segs):
            continue
        # each segment must change the sequence: a difference inside every non-first segment
        if not all(diff[s:e].any() for _, s, e in segs):
            continue
        path = PR.Path(segs)
        core = PR.path_sequence(panel, path)
        if core in ipd_seqs:
            continue
        ex = PR.path_exons(panel, path, {0: a['core_exons'], 1: b['core_exons']})
        if ex is None:
            continue
        break
    else:
        return None
    last = a if segs[-1][0] == 0 else b
    if kind == 'near_mosaic':
        ex2 = ex[1]
        pos = rng.randrange(ex2[0], ex2[1])
        core_cols = [c for (r, s, e) in segs for c in range(s, e) if P[r, c] < 4]
        col = core_cols[pos]
        present = {row[col] for row, n in zip(msa.aligned, msa.ipd_names) if n}
        alts = [x for x in 'ACGT' if x not in present]
        if not alts:
            return None
        core = core[:pos] + rng.choice(alts) + core[pos + 1:]
    if kind == 'mosaic1_partner_source':
        partner = b
    else:
        others = [x for x in pool if x['core_sequence'] not in (a['core_sequence'], b['core_sequence'])]
        if not others:
            return None
        partner = rng.choice(others)
    donors = sorted({a['donor'], b['donor'], partner['donor']})
    if not comp_ok(sorted({a['donor'], b['donor']}), excl) or \
            (partner is not b and not comp_ok(donors, excl)):
        return None
    truth = {'1': hap(a, last, core, ex),
             '2': hap(partner, partner, partner['core_sequence'], partner['core_exons'])}
    return truth, dict(sources=[(a if r == 0 else b)['source_id'] for r, _, _ in segs],
                       switch_columns=sw, partner=partner['source_id'], component_donors=donors)


def hap(first, last, core, exons):
    return dict(core=core, exons=[list(e) for e in exons],
                full=first['sequence'][:C.FLANK] + core + last['sequence'][-C.FLANK:])


def ensure_reads(cid):
    d = C.syn_dir(cid)
    r1, r2 = d / 'sim/sample_R1.fq', d / 'sim/sample_R2.fq'
    if r1.exists() and r2.exists() and (d / 'sim/DONE').exists():
        return r1, r2
    p2 = C.import_by_path('phase2_lib_pinned', P2LIB)
    p2.verify_simulator_contract()
    case = C.load_json(d / 'CASE.json')
    truth = C.load_json(d / 'TRUTH_SEQ.json')
    # template order is neutral: the proposer never sees read names anyway
    p2.simulate_diploid_locus(truth['1']['full'], truth['2']['full'], C.FLANK,
                              len(truth['1']['core']), case['seeds']['1'], case['seeds']['2'],
                              d / 'sim', threads=2)
    for f in (d / 'sim').iterdir():
        if f.name not in ('sample_R1.fq', 'sample_R2.fq'):
            if f.is_dir():
                import shutil
                shutil.rmtree(f)
            else:
                f.unlink()
    (d / 'sim/DONE').write_text(C.sha256_file(P2LIB) + '\n')
    return r1, r2


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--build', action='store_true')
    ap.add_argument('--seed', type=int, default=SEED)
    ap.add_argument('--set', default='SYN', choices=C.SYN_SETS)
    a = ap.parse_args()
    C.verify_inputs()
    if a.build:
        if (C.S / 'work' / a.set / 'CASES.json').exists():
            raise SystemExit('set already built: ' + a.set)
        build(a.seed, a.set)
