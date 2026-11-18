"""Development: where synthetic mosaics are lost.
 - proposed in arm C (pass 1) but not recovered after pass 2: why (tie / end length / other)
 - proposed in round-1 arm P but suppressed in arm C by the final thresholds (which decision)
 - oracle arm E misses: same classification.  Counts and names only."""
import collections
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import edlib
import common as C

MOS = {'mosaic1_partner_source', 'mosaic1_partner_other', 'mosaic2'}
ipd = {n: s for n, i, s in C.read_fasta_headers_and_seqs(C.GENOMIC_FA) if C.gene_of_name(n) in C.GENES}


def classify(table, gene, truth_seq, novel_name_tag):
    rows = [l.split('\t') for l in Path(table).read_text().splitlines() if l.startswith(f'HLA-{gene}\t')]
    labels = []
    for v in rows:
        for x in (v[2], v[5]):
            labels += [l for l in x.split(',') if l != '.']
    if any(novel_name_tag in l for l in labels):
        return 'novel_candidate_in_output_other_seq'
    best = None
    for l in labels:
        s = ipd.get(l)
        if not s:
            continue
        hw = edlib.align(truth_seq, s, mode='HW', task='distance')['editDistance']
        if best is None or hw < best[0]:
            best = (hw, l, len(s))
    if best and best[0] <= 5 and best[2] > len(truth_seq):
        return 'end_length_truth_shorter_than_near_identical_called_IPD'
    if best and best[0] <= 5:
        return 'near_identical_IPD_called_not_shorter'
    return 'other_IPD_called'


records = C.load_records()
msas = C.load_msas(records)


def syn_span(case):
    msa = msas[case['gene']]
    rows_ = [msa.row_of_source[x] for x in case['sources'][:2]]
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
        ok = ok and sum(1 for c in range(l, r + 1) if a[c] != '-') <= 450
    return ok


out = collections.Counter()
detail = []
for cid in C.load_json(C.S / 'work/SYN/CASES.json')['cases']:
    d = C.S / 'work/SYN' / cid
    case = C.load_json(d / 'CASE.json')
    if case['kind'] not in MOS:
        continue
    g = case['gene']
    tsha = case['truth_sha']['1']
    tseq = C.load_json(d / 'TRUTH_SEQ.json')['1']['core']
    pP = C.load_json(d / 'P/proposals.json')['genes'][g]
    pC = C.load_json(d / 'C/proposals.json')['genes'][g]
    hitP = tsha in {c['sha256'] for c in pP['candidates'] if 'name' in c}
    hitC = tsha in {c['sha256'] for c in pC['candidates'] if 'name' in c}
    rc = C.load_json(d / 'C/COMPLETE.json')
    tabC = rc['calls_table']
    recC = any(f'MOS_{tsha[:12]}' in l for l in Path(tabC).read_text().splitlines() if l.startswith(f'HLA-{g}\t'))
    re_ = C.load_json(d / 'E/COMPLETE.json')
    recE = any(f'ORC_{tsha[:12]}' in l for l in Path(re_['calls_table']).read_text().splitlines()
               if l.startswith(f'HLA-{g}\t'))
    out['mosaics'] += 1
    out['round1_P_hit'] += hitP
    out['C_proposal_hit'] += hitC
    out['C_recovered'] += recC
    out['E_recovered'] += recE
    if hitP and not hitC:
        out['suppressed_by_final_params:' + pC['decision']] += 1
    if hitC and not recC:
        out['lost_in_pass2:' + classify(tabC, g, tseq, 'MOS_')] += 1
    if not recE:
        out['oracle_miss:' + classify(re_['calls_table'], g, tseq, 'ORC_')] += 1
    st = ('spanned' if syn_span(case) else 'not_spanned') + '|' + \
        ('partner=source' if case['kind'] == 'mosaic1_partner_source' else 'partner=other')
    out['stratum_n:' + st] += 1
    out['stratum_C:' + st] += recC
    out['stratum_E:' + st] += recE
    detail.append((cid, case['kind'], st, hitP, hitC, recC, recE, pC['decision']))
for r in detail:
    print(*r, sep='\t')
for k in sorted(out):
    print('COUNT', k, out[k])
