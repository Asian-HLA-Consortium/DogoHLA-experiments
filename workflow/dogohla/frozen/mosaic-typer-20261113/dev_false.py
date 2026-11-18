"""Development diagnosis of false and wrong mosaics (round 3, step 1). Counts/positions/names only.

For every arm-C false MOS final call on an exact-IPD copy (score_rows_dev.tsv) and every synthetic
mosaic whose proposal set lacks the truth but holds another MOS candidate:
  - differences candidate vs true allele: events, positions (truth coordinates), feature (UTR/exon/intron)
  - read support at each event: reads (pass-1 T1K candidate reads) closer to the true allele vs closer
    to the candidate in a +-200 bp window, by strand; for candidate-supporting reads, the gene of
    their bwa primary placement (same gene / partner copy / paralog)
  - T1K genotype rows (pass 1 and pass 2) for the gene: labels, abundance, quality
Also: D-sim IPD-absent copies recovered by arm E but not C -> arm-C proposer decision for that gene.
Writes $S/tmp/dev_false.json and prints a summary.
"""
import collections
import csv
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import edlib
import pysam
import common as C
import part0 as P0

ENV = Path('/encrypted/e3008/hohndor-dogohla-20261020/mosaic-read-mapping-20261104/env/bin')
MOS = {'mosaic1_partner_source', 'mosaic1_partner_other', 'mosaic2'}
W = 200


def diff_events(t, m):
    r = edlib.align(m, t, mode='NW', task='path')
    ev, tp, mp = [], 0, 0
    for n, op in re.findall(r'(\d+)([=XID])', r['cigar']):
        n = int(n)
        if op == '=':
            tp += n; mp += n
        elif op == 'X':
            ev.append(('X', tp, mp, n)); tp += n; mp += n
        elif op == 'I':   # extra in candidate (query)
            ev.append(('I', tp, mp, n)); mp += n
        else:             # missing in candidate
            ev.append(('D', tp, mp, n)); tp += n
    return ev


def read_fq(path):
    out = []
    with open(path) as f:
        while True:
            h = f.readline()
            if not h:
                break
            s = f.readline().strip(); f.readline(); f.readline()
            out.append((h[1:].split()[0], s))
    return out


def t1k_rows(table, g):
    out = []
    for l in Path(table).read_text().splitlines():
        v = l.split('\t')
        if v[0] == f'HLA-{g}':
            out.append(v[1:8])
    return out


def analyse(tag, d, g, truth_seq, truth_exons, partner_seq, cand_name, cand_seq, pass1_dir, pass2_table, work):
    ev = diff_events(truth_seq, cand_seq)
    L = len(truth_seq)
    feats = [P0.features(min(max(e[1], 0), L - 1), truth_exons, L) for e in ev]
    # mapping target: truth, candidate, partner, pass-1 called alleles of all genes
    geno = C.read_t1k_raw(pass1_dir / 'calls_genotype.tsv')
    wanted = set()
    for gg, row in geno.items():
        for sl in row['slots'][:max(row['copies'], 0)]:
            wanted |= set(sl['labels'][:3])
    seqs = {}
    for name, ints, s in C.read_fasta_headers_and_seqs(C.GENOMIC_FA):
        if name in wanted:
            seqs[name] = s
    seqs['TRUTH'] = truth_seq
    seqs['CAND'] = cand_seq
    if partner_seq and partner_seq != truth_seq:
        seqs['PARTNER'] = partner_seq
    work.mkdir(parents=True, exist_ok=True)
    fa = work / 't.fa'
    fa.write_text(''.join(f'>{n}\n{s}\n' for n, s in seqs.items()))
    subprocess.run([str(ENV / 'bwa'), 'index', str(fa)], check=True, capture_output=True)
    sam = work / 'a.sam'
    with open(sam, 'w') as f:
        subprocess.run([str(ENV / 'bwa'), 'mem', '-t', '2', '-v', '1', str(fa),
                        str(pass1_dir / 'calls_candidate_1.fq'), str(pass1_dir / 'calls_candidate_2.fq')],
                       stdout=f, stderr=subprocess.DEVNULL, check=True)
    gene_of = lambda n: g if n in ('TRUTH', 'CAND') else ('PARTNER' if n == 'PARTNER' else C.gene_of_name(n))
    near = collections.defaultdict(list)   # event index -> reads
    with pysam.AlignmentFile(str(sam)) as af:
        for rec in af.fetch(until_eof=True):
            if rec.is_secondary or rec.is_supplementary or rec.is_unmapped:
                continue
            pg = gene_of(rec.reference_name)
            near_gene = pg in (g, 'PARTNER') or C.gene_of_name(rec.reference_name)[:2] == g[:2] or \
                (g in C.CLASS_I and pg in ('A', 'B', 'C', 'E', 'F', 'G', 'H', 'J', 'K', 'L', 'N', 'P', 'S', 'T', 'U', 'V', 'W', 'Y'))
            if not near_gene:
                continue
            for i, e in enumerate(ev):
                near[i].append((rec.query_sequence, rec.is_reverse, pg))
    events = []
    for i, e in enumerate(ev):
        t0, t1 = max(0, e[1] - W), min(L, e[1] + W)
        m0, m1 = max(0, e[2] - W), min(len(cand_seq), e[2] + W)
        tw, mw = truth_seq[t0:t1], cand_seq[m0:m1]
        sup = collections.Counter()
        for q, rev, pg in near[i]:
            dt = edlib.align(q, tw, mode='HW', task='distance', k=8)['editDistance']
            dm = edlib.align(q, mw, mode='HW', task='distance', k=8)['editDistance']
            if dt < 0 and dm < 0:
                continue
            dt = 99 if dt < 0 else dt
            dm = 99 if dm < 0 else dm
            if dt == dm:
                continue
            side = 'truth' if dt < dm else 'cand'
            sup[(side, 'rev' if rev else 'fwd')] += 1
            if side == 'cand':
                sup[('cand_placed', pg)] += 1
                if min(dt, dm) > 2:
                    sup[('cand_read_ed_gt2', '')] += 1
        events.append(dict(type=e[0], truth_pos=e[1], length=e[3], feature=feats[i],
                           support={f'{a}:{b}': v for (a, b), v in sup.items()}))
    for p in work.iterdir():
        p.unlink()
    return dict(tag=tag, donor=d, gene=g, candidate=cand_name, truth_len=L, cand_len=len(cand_seq),
                n_events=len(ev), events=events,
                only_end_or_utr=all(f.startswith('UTR') for f in feats) if feats else None,
                t1k_pass1=t1k_rows(pass1_dir / 'calls_genotype.tsv', g), t1k_pass2=t1k_rows(pass2_table, g))


def main():
    S = C.S
    rows = list(csv.DictReader(open(S / 'results/score_rows_dev.tsv'), delimiter='\t'))
    T = lambda x: x == 'True'
    records = C.load_records()
    rec_by = {(r['donor'], r['gene'], r['haplotype']): r for r in records}
    out = dict(false=[], wrong=[], dsim_gap=[])
    for r in rows:
        if not (r['arm'] == 'C' and T(r['false_mos']) and not T(r['ipd_absent'])):
            continue
        d, g, h = r['donor'], r['gene'], r['haplotype']
        crec = C.load_json(S / 'work' / r['readset'] / d / 'C/COMPLETE.json')
        cands = {n: s for n, i, s in C.read_fasta_headers_and_seqs(S / 'work' / r['readset'] / d / 'C/additions.fa')}
        name = r['call'].split('|')[0]
        if r['readset'] == 'SYN':
            tr = C.load_json(S / 'work/SYN' / d / 'TRUTH_SEQ.json')
            tseq, tex = tr[h]['core'], tr[h]['exons']
            pseq = tr['1' if h == '2' else '2']['core']
            p1 = S / 'work/SYN' / d / 'A'
        else:
            rr = rec_by[(d, g, h)]
            tseq, tex = rr['core_sequence'], rr['core_exons']
            pr = rec_by.get((d, g, '1' if h == '2' else '2'))
            pseq = pr['core_sequence'] if pr else None
            p1 = C.HPRC_CONTROL / 't1k/runs' / d / 'genomic'
        res = analyse('false', d, g, tseq, tex, pseq, name, cands[name], p1, Path(crec['calls_table']),
                      S / 'tmp/devfalse' / d)
        res.update(readset=r['readset'], haplotype=h, truth4=r['truth4'], class_full=r['class_full'])
        out['false'].append(res)
        print('FALSE', res['readset'], d, g, h, r['truth4'], name, 'events', res['n_events'],
              'only_end/UTR', res['only_end_or_utr'], 'lens', res['truth_len'], res['cand_len'], flush=True)
    for cid in C.load_json(S / 'work/SYN/CASES.json')['cases']:
        case = C.load_json(S / 'work/SYN' / cid / 'CASE.json')
        if case['kind'] not in MOS:
            continue
        g = case['gene']
        prop = C.load_json(S / 'work/SYN' / cid / 'C/proposals.json')['genes'][g]
        props = [c for c in prop['candidates'] if 'name' in c]
        if not props or case['truth_sha']['1'] in {c['sha256'] for c in props}:
            continue
        tr = C.load_json(S / 'work/SYN' / cid / 'TRUTH_SEQ.json')
        cands = {n: s for n, i, s in C.read_fasta_headers_and_seqs(S / 'work/SYN' / cid / 'C/additions.fa')}
        crec = C.load_json(S / 'work/SYN' / cid / 'C/COMPLETE.json')
        for c in props:
            res = analyse('wrong', cid, g, tr['1']['core'], tr['1']['exons'], tr['2']['core'], c['name'],
                          cands[c['name']], S / 'work/SYN' / cid / 'A', Path(crec['calls_table']),
                          S / 'tmp/devfalse' / cid)
            res.update(kind=case['kind'], decision=prop['decision'], best_switches=prop.get('best_switches'),
                       best_linkage=prop.get('best_linkage'), best_evidence=prop.get('best_evidence_cost'))
            out['wrong'].append(res)
            print('WRONG', cid, case['kind'], c['name'], 'events', res['n_events'], 'only_end/UTR',
                  res['only_end_or_utr'], flush=True)
    # D-sim gap
    idx = {(r['readset'], r['donor'], r['gene'], r['haplotype'], r['arm']): r for r in rows}
    for r in rows:
        if r['readset'] != 'D-sim' or r['arm'] != 'E' or not T(r['ipd_absent']) or not T(r['recovered']):
            continue
        c = idx[('D-sim', r['donor'], r['gene'], r['haplotype'], 'C')]
        if T(c['recovered']):
            continue
        p = C.load_json(S / 'work/D-sim' / r['donor'] / 'C/proposals.json')['genes'][r['gene']]
        out['dsim_gap'].append(dict(donor=r['donor'], gene=r['gene'], haplotype=r['haplotype'],
                                    class_full=r['class_full'], class_ipd=r['class_ipd'],
                                    decision=p['decision'], n_proposed=p.get('n_proposed', 0),
                                    pass1_cost=p.get('pass1_evidence_cost'), zero=p.get('zero_switch_evidence_cost'),
                                    best=p.get('best_evidence_cost'), switches=p.get('best_switches'),
                                    offtarget=p['evidence_stats'].get('excluded_offtarget_edit_gt_max', 0),
                                    reads=p['evidence_stats'].get('reads_assigned')))
    C.save_json(S / 'tmp/dev_false.json', out)
    print('DSIM_GAP', collections.Counter((x['decision'], x['class_ipd']) for x in out['dsim_gap']))


if __name__ == '__main__':
    main()
