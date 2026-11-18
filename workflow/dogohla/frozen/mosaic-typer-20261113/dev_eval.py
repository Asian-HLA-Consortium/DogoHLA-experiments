"""Development evaluation of proposals (arm P or C) against truth: per gene, did the proposer
propose, does any proposal equal a truth copy, and how many proposals are false. Counts only.
  python dev_eval.py ARM READSET[,READSET]"""
import collections
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C

arm, readsets = sys.argv[1], sys.argv[2].split(',')
records = C.load_records()
classes = C.load_classes()
truth = collections.defaultdict(dict)
for r in records:
    truth[(r['donor'], r['gene'])][r['haplotype']] = (C.sha256_str(r['core_sequence']),
                                                    classes[r['source_id']]['full_label_L150'],
                                                    classes[r['source_id']]['ipd_label_L150'],
                                                    not r['genomic_labels'])
msas = C.load_msas(records)


def syn_span(case):
    """True if every switch's bounding distinguishing sites (between the two source rows) lie within
    450 bp of each other (spanned by one read pair)."""
    msa = msas[case['gene']]
    rows_ = [msa.row_of_source[s] for s in case['sources'][:2]]
    a, b = msa.aligned[rows_[0]], msa.aligned[rows_[1]]
    ok = True
    for sw in case['switch_columns']:
        l = sw - 1
        while l >= 0 and a[l] == b[l]:
            l -= 1
        r = sw
        while r < len(a) and a[r] == b[r]:
            r += 1
        if l < 0 or r >= len(a):
            ok = False
            continue
        bases = sum(1 for c in range(l, r + 1) if a[c] != '-')
        ok = ok and bases <= 450
    return ok


tot = collections.Counter()
strata = collections.Counter()
rows = []
for rs in readsets:
    for d in sorted((C.S / 'work' / rs).iterdir()):
        w = d / arm
        if not (w / 'COMPLETE.json').exists():
            tot['missing_or_failed'] += 1
            continue
        p = C.load_json(w / 'proposals.json')
        for g, rec in p['genes'].items():
            if rs == 'SYN':
                case = C.load_json(d / 'CASE.json')
                if case['gene'] != g:
                    continue
                t = {h: (sha, case['kind'] if h == '1' else 'exact', '', sha not in case['ipd_present_shas'])
                     for h, sha in case['truth_sha'].items()}
            else:
                t = truth[(d.name, g)]
            props = [c for c in rec['candidates'] if 'name' in c]
            hits = [c for c in props if c['sha256'] in {v[0] for v in t.values()}]
            if rs == 'SYN' and case['kind'] in ('mosaic1_partner_source', 'mosaic1_partner_other', 'mosaic2'):
                hit = case['truth_sha']['1'] in {c['sha256'] for c in props}
                key = ('pair_spanned' if syn_span(case) else 'not_spanned',
                       'partner=source' if case['kind'] == 'mosaic1_partner_source' else 'partner=other')
                strata[key + ('n',)] += 1
                strata[key + ('hit',)] += hit
            if rs == 'SYN' and case['kind'] == 'exact' or rs != 'SYN':
                n_exact = sum(1 for v in t.values() if not v[3])
                false_on_exact = len(props) - len(hits) if n_exact == len(t) else 0
                strata[('exact_gene_with_false_proposal', g)] += bool(false_on_exact)
                strata[('exact_gene', g)] += n_exact == len(t)
            classes_here = '/'.join(sorted(v[1] for v in t.values()))
            reach = [h for h, v in t.items() if v[3] and v[2] == 'mosaic'] if rs != 'SYN' else \
                [h for h, v in t.items() if v[1].startswith('mosaic')]
            tot['genes'] += 1
            tot['reachable_copies'] += len(reach)
            tot['proposed_genes'] += bool(props)
            tot['proposals'] += len(props)
            tot['true_proposals'] += len(hits)
            tot['false_proposals'] += len(props) - len(hits)
            rows.append((rs, d.name, g, classes_here, rec['decision'], len(props), len(hits),
                         rec.get('pass1_evidence_cost'), rec.get('zero_switch_evidence_cost'),
                         rec.get('best_evidence_cost'), rec.get('best_switches'), rec.get('best_linkage'),
                         rec.get('diag_zero_unexplained_n')))
for r in rows:
    if r[5] or 'mosaic' in r[3]:
        print(*r, sep='\t')
print(dict(tot))
for k in sorted(strata):
    print('STRATUM', k, strata[k])
