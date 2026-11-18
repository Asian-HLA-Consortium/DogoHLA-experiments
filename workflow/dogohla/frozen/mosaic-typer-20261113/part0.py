"""Part 0: what happens to mosaic, near-mosaic, private and (sampled) exact truth copies today.

  python part0.py --gene G      per-gene: truth 2/3-field naming from CDS/protein identity, switch
                                geography for the full panel (from mosaic_classes.tsv) and for the
                                IPD-only panel (recomputed with the atlas classifier), linkage reach,
                                partner identity.  -> $S/tmp/part0/<G>.json
  python part0.py --combine     joins existing arm-A (genomic) T1K calls of every read set and writes
                                results/part0_mosaics_today.tsv, results/part0_switches.tsv,
                                results/PART0.md (tables; the interpretation paragraph is added by hand).
Only existing runs are read. No sequences are written to results.
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

SEED = 20261113
NONEXACT = ('mosaic', 'near-mosaic', 'private')


def features(pos, exons, core_len):
    """Feature of a core position given CDS exons [(s,e)] in core coordinates."""
    if pos < exons[0][0]:
        return 'UTR5'
    if pos >= exons[-1][1]:
        return 'UTR3'
    for k, (s, e) in enumerate(exons):
        if s <= pos < e:
            return f'exon{k + 1}'
        if k + 1 < len(exons) and e <= pos < exons[k + 1][0]:
            return f'intron{k + 1}'
    return 'UTR3'


def abd(feat, gene):
    return feat in (('exon2', 'exon3') if gene in C.CLASS_I else ('exon2',))


def diff_events(a, b, cols):
    """Number of variant events (runs of consecutive differing columns) between rows a and b
    restricted to the given sorted column list."""
    n, prev = 0, None
    for c in cols:
        if a[c] != b[c]:
            if prev is None or c != prev + 1:
                n += 1
            prev = c
    return n


def ipd_name_maps(msa: C.GeneMSA):
    cds_names, prot_names = collections.defaultdict(set), collections.defaultdict(set)
    for i, names in enumerate(msa.ipd_names):
        if not names:
            continue
        ex = msa.ipd_exons[i]
        cds = C.cds_of(msa.seqs[i], ex)
        prot = C.translate(cds)
        for n in names:
            v3, v2 = C.normalized(n, 3), C.normalized(n, 2)
            if v3:
                cds_names[cds].add(v3)
            if v2:
                prot_names[prot].add(v2)
    return cds_names, prot_names


def analyse_switches(gene, rec, partner, msa, segs, switches, panel_tag):
    """segs: list of source row indices per segment; switches: list of dicts with
    switch_base, agree_base_start, agree_base_end (query core coords, inclusive)."""
    qrow = msa.row_of_source[rec['source_id']]
    qa = msa.aligned[qrow]
    bcols = msa.base_cols(qrow)
    L = len(bcols)
    exons = rec['core_exons']
    prow = msa.row_of_source[partner['source_id']] if partner else None
    out = []
    for k, sw in enumerate(switches):
        la, ra = msa.aligned[segs[k]], msa.aligned[segs[k + 1]]
        s0, s1 = sw['agree_base_start'], sw['agree_base_end']
        ldiv, rdiv = s0 - 1, s1 + 1
        feats = sorted({features(p, exons, L) for p in range(max(0, s0), min(L, s1 + 1))},
                       key=lambda f: (f[0] != 'U', f))
        row = dict(gene=gene, source_id=rec['source_id'], panel=panel_tag, switch_index=k,
                   left_source=';'.join(msa.ipd_names[segs[k]][:1] or msa.cohort_sources[segs[k]][:1]),
                   right_source=';'.join(msa.ipd_names[segs[k + 1]][:1] or msa.cohort_sources[segs[k + 1]][:1]),
                   switch_base=sw['switch_base'], agree_base_start=s0, agree_base_end=s1,
                   agree_interval_bp=s1 - s0 + 1, interval_features=','.join(feats),
                   interval_has_noncoding=any(not f.startswith('exon') for f in feats),
                   interval_all_exonic=all(f.startswith('exon') for f in feats),
                   left_div_feature=features(ldiv, exons, L) if ldiv >= 0 else 'start',
                   right_div_feature=features(rdiv, exons, L) if rdiv < L else 'end')
        row['bounding_divergence_in_abd'] = (abd(row['left_div_feature'], gene) or
                                             abd(row['right_div_feature'], gene))
        for w in (150, 500):
            lo, hi = max(0, sw['switch_base'] - w), min(L - 1, sw['switch_base'] + w)
            row[f'n_div_variants_pm{w}'] = diff_events(la, ra, bcols[lo:hi + 1])
            lcols = bcols[max(0, s0 - w):max(0, s0)]
            rcols = bcols[min(L, s1 + 1):min(L, s1 + 1 + w)]
            row[f'n_div_left_within{w}'] = diff_events(la, ra, lcols)
            row[f'n_div_right_within{w}'] = diff_events(la, ra, rcols)
        span = (rdiv - ldiv + 1) if (ldiv >= 0 and rdiv < L) else None
        row['bounding_span_bp'] = span
        row['linkable_by_one_read_150'] = span is not None and span <= 150
        row['linkable_by_pair_le450'] = span is not None and span <= 450
        if prow is not None and ldiv >= 0 and rdiv < L:
            c0, c1 = bcols[ldiv], bcols[rdiv]
            pa = msa.aligned[prow]
            eq_l = pa[c0:c1 + 1] == la[c0:c1 + 1]
            eq_r = pa[c0:c1 + 1] == ra[c0:c1 + 1]
            row['partner_equals'] = 'left' if eq_l and not eq_r else 'right' if eq_r and not eq_l \
                else 'both' if eq_l and eq_r else 'neither'
        else:
            row['partner_equals'] = 'NA'
        out.append(row)
    return out


def per_gene(gene):
    records = C.load_records()
    classes = C.load_classes()
    msa = C.load_msas(records, [gene])[gene]
    atlas_src = C.ATLAS / 'src'
    sys.path.insert(0, str(atlas_src))
    mc = C.import_by_path('mosaic_classifier_atlas', atlas_src / 'mosaic_classifier.py')
    cds_names, prot_names = ipd_name_maps(msa)
    recs = [r for r in records if r['gene'] == gene]
    by_dg = collections.defaultdict(list)
    for r in recs:
        by_dg[r['donor']].append(r)
    ipd_panel = {}
    for i, names in enumerate(msa.ipd_names):
        if names:
            ipd_panel[names[0]] = msa.aligned[i]
    name_to_row = {names[0]: i for i, names in enumerate(msa.ipd_names) if names}
    copies, switch_rows, checks = [], [], []
    for r in recs:
        cl = classes[r['source_id']]
        seq, ex = r['core_sequence'], r['core_exons']
        cds = C.cds_of(seq, ex)
        prot = C.translate(cds)
        t4 = sorted({C.normalized(a, 4) for a in r['genomic_labels']} - {None})
        copy = dict(source_id=r['source_id'], donor=r['donor'], gene=gene,
                    haplotype=r['haplotype'], cohort=r['resource'],
                    ipd_absent=not r['genomic_labels'],
                    full_label=cl['full_label_L150'], full_detailed=cl['full_detailed_L150'],
                    ipd_label=cl['ipd_label_L150'], ipd_detailed=cl['ipd_detailed_L150'],
                    truth4=t4[0] if len(t4) == 1 else None,
                    truth3=sorted(cds_names.get(cds, ())), truth2=sorted(prot_names.get(prot, ())),
                    cds_in_ipd=cds in cds_names, protein_in_ipd=prot in prot_names,
                    cds_len=len(cds))
        copies.append(copy)
        selected = cl['full_label_L150'] in NONEXACT or cl['ipd_label_L150'] in NONEXACT
        if not selected:
            continue
        partners = [x for x in by_dg[r['donor']] if x['source_id'] != r['source_id']]
        partner = partners[0] if len(partners) == 1 else None
        # recompute both panels with the atlas classifier so both carry identical fields
        qa = msa.aligned[msa.row_of_source[r['source_id']]]
        excl = None
        for tag in ('ipd', 'full'):
            if tag == 'ipd':
                if cl['ipd_label_L150'] not in ('mosaic', 'near-mosaic'):
                    continue
                panel = ipd_panel
            else:
                if cl['full_label_L150'] not in ('mosaic', 'near-mosaic'):
                    continue
                if excl is None:
                    excl = C.exclusion_map()[r['donor']]
                panel = dict(ipd_panel)
                for i in range(len(msa.uids)):
                    ok = [sid for sid in msa.cohort_sources[i] if msa.donor_of_source[sid] not in excl]
                    if ok and not msa.ipd_names[i]:
                        panel[ok[0]] = msa.aligned[i]
            res = mc.classify_query_against_panel(gene, r['source_id'], qa, panel, min_seg_bp=150)
            expect_k = cl[f'{tag}_k_L150']
            checks.append(dict(source_id=r['source_id'], panel=tag, label=res.detailed_label,
                               k=res.k_switches, atlas_label=cl[f'{tag}_detailed_L150'],
                               atlas_k=expect_k))
            if res.k_switches and res.switches:
                segs = []
                for s in res.segment_sources:
                    segs.append(name_to_row[s] if s in name_to_row else msa.row_of_source[s])
                sws = [dict(switch_base=x.switch_base_0based, agree_base_start=x.agreeing_base_start,
                            agree_base_end=x.agreeing_base_end) for x in res.switches]
                switch_rows += analyse_switches(gene, r, partner, msa, segs, sws, tag)
                copy[f'{tag}_segment_sources'] = ';'.join(res.segment_sources)
    out = C.S / 'tmp/part0' / f'{gene}.json'
    C.save_json(out, dict(gene=gene, copies=copies, switches=switch_rows, checks=checks))
    bad = [c for c in checks if str(c['k']) != str(c['atlas_k']) and c['atlas_k'] != '']
    print(gene, 'copies', len(copies), 'switch rows', len(switch_rows), 'checks', len(checks),
          'k disagreements with atlas', len(bad))
    for b in bad[:20]:
        print('  disagreement', b)


def slot_value(slot, k):
    vals = {C.normalized(l, k) for l in slot}
    return next(iter(vals)) if len(vals) == 1 and None not in vals else None


def match(copy, slot):
    """Per-field correctness of one called slot against one truth copy."""
    out = {}
    for k in (2, 3, 4):
        if not slot:
            out[k] = (False, False)
            continue
        truth = {copy['truth4']} - {None} if k == 4 else set(copy[f'truth{k}'])
        v = slot_value(slot, k)
        contains = bool(truth & {C.normalized(l, k) for l in slot})
        out[k] = (v is not None and v in truth, contains)
    return out


def assign(copies, slots):
    """Best assignment of two slots to two copies (4-field weighted over 3 over 2)."""
    w = {4: 100, 3: 10, 2: 1}
    best = None
    for perm in itertools.permutations(range(len(slots)), len(copies)):
        ms = [match(c, slots[p]) for c, p in zip(copies, perm)]
        score = sum(w[k] * m[k][0] for m in ms for k in (2, 3, 4))
        if best is None or score > best[0]:
            best = (score, perm, ms)
    return best


def combine():
    genes = {}
    for g in C.GENES:
        genes[g] = C.load_json(C.S / 'tmp/part0' / f'{g}.json')
    copies = [c for g in C.GENES for c in genes[g]['copies']]
    switches = [s for g in C.GENES for s in genes[g]['switches']]
    checks = [s for g in C.GENES for s in genes[g]['checks']]
    rng = random.Random(SEED)
    sel = [c for c in copies if c['full_label'] in NONEXACT or c['ipd_label'] in NONEXACT]
    pool = collections.defaultdict(list)
    for c in copies:
        if c['full_label'] == 'exact' and c['ipd_label'] == 'exact' and not c['ipd_absent']:
            pool[c['gene']].append(c)
    need = collections.Counter(c['gene'] for c in sel)
    sample = []
    for g in C.GENES:
        p = sorted(pool[g], key=lambda c: c['source_id'])
        sample += rng.sample(p, min(need[g], len(p)))
    for c in sample:
        c['exact_sample'] = True
    chosen = {c['source_id'] for c in sel} | {c['source_id'] for c in sample}
    by_dg = collections.defaultdict(list)
    for c in copies:
        by_dg[(c['donor'], c['gene'])].append(c)
    runs = C.readset_runs()
    hprc_dsim = set(runs['D-sim'])
    rows, missing = [], []
    for rs, donors in runs.items():
        for donor, run in donors.items():
            name = rs if rs != 'real' else ('D-real' if donor in hprc_dsim else 'T-real')
            try:
                calls = C.read_t1k(C.verified_table(run))
            except Exception as e:  # recorded, never dropped silently
                missing.append(dict(readset=name, donor=donor, error=repr(e)))
                calls = None
            for g in C.GENES:
                cps = sorted(by_dg[(donor, g)], key=lambda c: c['haplotype'])
                if not cps:
                    continue
                slots = calls.get(g, [[], []]) if calls else [[], []]
                best = assign(cps, slots)
                for c, p, m in zip(cps, best[1], best[2]):
                    if c['source_id'] not in chosen:
                        continue
                    row = dict(readset=name, **{k: v for k, v in c.items() if k not in ('truth2', 'truth3')})
                    row['truth3'] = '|'.join(c['truth3'])
                    row['truth2'] = '|'.join(c['truth2'])
                    row['run_complete'] = calls is not None
                    row['called_slot'] = ','.join(slots[p]) if slots[p] else '.'
                    row['no_call'] = not slots[p]
                    for k in (2, 3, 4):
                        row[f'correct{k}'] = m[k][0]
                        row[f'slot_contains_truth{k}'] = m[k][1]
                        row[f'call_resolved{k}'] = slot_value(slots[p], k) is not None
                    row['exact_sample'] = bool(c.get('exact_sample'))
                    rows.append(row)
    res = C.S / 'results'
    cols = list(rows[0].keys())
    with open(res / 'part0_mosaics_today.tsv', 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows:
            f.write('\t'.join('' if r.get(k) is None else str(r.get(k)) for k in cols) + '\n')
    scols = list(switches[0].keys())
    with open(res / 'part0_switches.tsv', 'w') as f:
        f.write('\t'.join(scols) + '\n')
        for r in switches:
            f.write('\t'.join('' if r.get(k) is None else str(r.get(k)) for k in scols) + '\n')
    C.save_json(C.S / 'tmp/part0/combine_meta.json', dict(missing=missing, checks=checks,
                n_rows=len(rows), n_selected=len(sel), n_sample=len(sample)))
    write_md(rows, switches, checks, missing, copies)


def pct(a, b):
    return f'{a}/{b} ({100 * a / b:.0f}%)' if b else '0/0'


def write_md(rows, switches, checks, missing, copies):
    L = ['# Part 0: what T1K (arm A, IPD genomic) does with mosaic alleles today', '',
         'Generated by `src/part0.py --combine` from existing T1K `genomic` runs only. '
         'Copy classes are the atlas classes (L=150). "full" = IPD + leave-family-out cohort panel; '
         '"ipd" = IPD genomic alleles only (the panel arm C uses). Truth 3-field name of a copy = '
         'the 3-field names of IPD genomic.fa alleles with identical CDS; truth 2-field = the '
         '2-field names of IPD alleles with identical translated protein (to the first stop). '
         'A call is correct at k fields when T1K\'s slot resolves to one k-field value in that set; '
         '"contains" = the truth is among T1K\'s tied labels but the slot is unresolved.', '']
    if missing:
        L += [f'Runs that could not be read: {len(missing)} (listed in tmp/part0/combine_meta.json).', '']
    groups = [('full-panel mosaic', lambda r: r['full_label'] == 'mosaic'),
              ('IPD-only mosaic (reachable by arm C)', lambda r: r['ipd_label'] == 'mosaic'),
              ('IPD-only mosaic, IPD-absent', lambda r: r['ipd_label'] == 'mosaic' and r['ipd_absent']),
              ('near-mosaic (full)', lambda r: r['full_label'] == 'near-mosaic'),
              ('near-mosaic (ipd)', lambda r: r['ipd_label'] == 'near-mosaic'),
              ('private (full)', lambda r: r['full_label'] == 'private'),
              ('private (ipd)', lambda r: r['ipd_label'] == 'private'),
              ('exact sample (IPD-present)', lambda r: r['exact_sample'])]
    L += ['## 1. Truth naming of the selected copies', '',
          '| group | copies | IPD-absent | CDS equals an IPD CDS (3-field name exists) | '
          'protein equals an IPD protein (2-field name exists) | no 3-field name | no 2-field name |',
          '|---|---:|---:|---:|---:|---:|---:|']
    seen = {}
    for r in rows:
        seen[r['source_id']] = r
    for name, f in groups:
        cs = [r for r in seen.values() if f(r)]
        L.append(f"| {name} | {len(cs)} | {sum(c['ipd_absent'] for c in cs)} | "
                 f"{sum(c['cds_in_ipd'] for c in cs)} | {sum(c['protein_in_ipd'] for c in cs)} | "
                 f"{sum(not c['cds_in_ipd'] for c in cs)} | {sum(not c['protein_in_ipd'] for c in cs)} |")
    absent_all = [c for c in copies if c['ipd_absent']]
    L += ['', f"All {len(absent_all)} IPD-absent truth copies: {sum(not c['cds_in_ipd'] for c in absent_all)} "
          f"have no IPD allele with identical CDS (no 3-field name), "
          f"{sum(not c['protein_in_ipd'] for c in absent_all)} have no IPD allele with identical protein "
          '(no 2-field name).', '']
    L += ['## 2. What T1K called (per copy and read set)', '',
          'Correct counts use only copies for which a truth name exists at that resolution '
          '(4-field: IPD-present copies; 3/2-field: copies with a CDS/protein-identical IPD allele).', '',
          '| read set | group | copies | no-call | 4-field correct | 3-field correct | 2-field correct | '
          '2-field: truth among tied labels but unresolved |', '|---|---|---:|---:|---:|---:|---:|---:|']
    for rs in ('D-sim', 'D-real', 'T-sim', 'T-real'):
        for name, f in groups:
            cs = [r for r in rows if r['readset'] == rs and f(r)]
            if not cs:
                continue
            e4 = [c for c in cs if c['truth4']]
            e3 = [c for c in cs if c['truth3']]
            e2 = [c for c in cs if c['truth2']]
            L.append(f"| {rs} | {name} | {len(cs)} | {sum(c['no_call'] for c in cs)} | "
                     f"{pct(sum(c['correct4'] for c in e4), len(e4))} | "
                     f"{pct(sum(c['correct3'] for c in e3), len(e3))} | "
                     f"{pct(sum(c['correct2'] for c in e2), len(e2))} | "
                     f"{sum(c['slot_contains_truth2'] and not c['correct2'] for c in e2)} |")
    L += ['', '## 3. Where the switches lie', '']
    for tag in ('ipd', 'full'):
        sw = [s for s in switches if s['panel'] == tag]
        if not sw:
            continue
        n = len(sw)
        L += [f'### panel = {tag}: {n} switches in {len({s["source_id"] for s in sw})} copies', '',
              f"- agreeing interval contains intron/UTR positions (switch can be non-coding): "
              f"{pct(sum(s['interval_has_noncoding'] for s in sw), n)}; entirely exonic: "
              f"{pct(sum(s['interval_all_exonic'] for s in sw), n)}",
              f"- a bounding divergent site lies in an antigen-binding exon: "
              f"{pct(sum(s['bounding_divergence_in_abd'] for s in sw), n)}",
              f"- distinguishing variants within +-150 bp of the switch: median "
              f"{sorted(s['n_div_variants_pm150'] for s in sw)[n // 2]}; zero in "
              f"{sum(s['n_div_variants_pm150'] == 0 for s in sw)}; within +-500 bp: median "
              f"{sorted(s['n_div_variants_pm500'] for s in sw)[n // 2]}; zero in "
              f"{sum(s['n_div_variants_pm500'] == 0 for s in sw)}",
              f"- both bounding divergent sites within one 150 bp read: "
              f"{pct(sum(s['linkable_by_one_read_150'] for s in sw), n)}; within a 450 bp fragment: "
              f"{pct(sum(s['linkable_by_pair_le450'] for s in sw), n)}",
              f"- partner haplotype across the switch region equals: " +
              ', '.join(f'{k} {v}' for k, v in sorted(collections.Counter(s['partner_equals'] for s in sw).items())),
              '- interval features: ' + ', '.join(
                  f'{k} {v}' for k, v in collections.Counter(s['interval_features'] for s in sw).most_common(12)), '']
    bad = [c for c in checks if str(c['k']) != str(c['atlas_k']) and c['atlas_k'] != '']
    L += ['## 4. Consistency with the atlas', '',
          f'Re-running the atlas classifier for {len(checks)} copy-panel pairs: {len(bad)} differ in k '
          'from mosaic_classes.tsv' + ('' if not bad else ' (see tmp/part0/combine_meta.json)') + '.', '']
    (C.S / 'results/PART0.md').write_text('\n'.join(L) + '\n')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--gene')
    ap.add_argument('--combine', action='store_true')
    a = ap.parse_args()
    C.verify_inputs()
    if a.combine:
        combine()
    else:
        per_gene(a.gene)
