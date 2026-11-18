"""Read-evidence acceptance of MOS final calls after pass 2 (round-3 deviation from the design's
final call rule; set on development data before any test result).

For every gene slot whose final call (IPD-prior rule) is a MOS_ candidate M:
  N = the named pass-1 allele (first label of either pass-1 slot) closest to M  (the displaced allele)
  O = the other slot's final allele (named or novel)
Differences between M and N are found by global alignment; at each difference event, reads (pass-1
T1K candidate reads with a bwa primary placement on M, N, O or another allele of the gene) are
classified by edit distance to +-200 bp windows of M, N and O: a read counts for M (N) only if it is
strictly closer to M (N) than to both others.
M is ACCEPTED iff at some event M-unique reads >= acc_min_strand on each strand and >= acc_min_frac of
the classified reads there, and at no event do N-unique reads exceed acc_max_contra. Otherwise the slot's
final call reverts to N (M stays a detected sequence).
"""
from __future__ import annotations

import collections
import re
import subprocess
from pathlib import Path

import edlib
import pysam

import common as C

ENV = Path('/encrypted/e3008/hohndor-dogohla-20261020/mosaic-read-mapping-20261104/env/bin')
W = 200


def events(n_seq, m_seq):
    r = edlib.align(m_seq, n_seq, mode='NW', task='path')
    ev, npos, mpos = [], 0, 0
    for k, op in re.findall(r'(\d+)([=XID])', r['cigar']):
        k = int(k)
        if op == '=':
            npos += k; mpos += k
            continue
        ev.append((npos, mpos))
        if op == 'X':
            npos += k; mpos += k
        elif op == 'I':
            mpos += k
        else:
            npos += k
    return ev


def window(seq, p):
    return seq[max(0, p - W):min(len(seq), p + W)]


def decide(gene, m_name, m_seq, n_name, n_seq, o_seq, fq1, fq2, other_alleles, work: Path, prm):
    work.mkdir(parents=True, exist_ok=True)
    tgt = dict(other_alleles)
    tgt.update({'ACC_M': m_seq, 'ACC_N': n_seq})
    if o_seq and o_seq not in (m_seq, n_seq):
        tgt['ACC_O'] = o_seq
    fa = work / 't.fa'
    fa.write_text(''.join(f'>{k}\n{v}\n' for k, v in tgt.items()))
    subprocess.run([str(ENV / 'bwa'), 'index', str(fa)], check=True, capture_output=True)
    sam = work / 'a.sam'
    with open(sam, 'w') as f:
        subprocess.run([str(ENV / 'bwa'), 'mem', '-t', '2', '-v', '1', str(fa), str(fq1), str(fq2)],
                       stdout=f, stderr=subprocess.DEVNULL, check=True)
    reads = []
    with pysam.AlignmentFile(str(sam)) as af:
        for rec in af.fetch(until_eof=True):
            if rec.is_secondary or rec.is_supplementary or rec.is_unmapped:
                continue
            n = rec.reference_name
            if n.startswith('ACC_') or C.gene_of_name(n) == gene:
                reads.append((rec.query_sequence, rec.is_reverse))
    for p in work.iterdir():
        p.unlink()
    ev = events(n_seq, m_seq)
    out = []
    for npos, mpos in ev:
        wm, wn = window(m_seq, mpos), window(n_seq, npos)
        wo = window(o_seq, npos) if o_seq else None
        cnt = collections.Counter()
        for q, rev in reads:
            dm = edlib.align(q, wm, mode='HW', task='distance', k=6)['editDistance']
            dn = edlib.align(q, wn, mode='HW', task='distance', k=6)['editDistance']
            if dm < 0 and dn < 0:
                continue
            do = edlib.align(q, wo, mode='HW', task='distance', k=6)['editDistance'] if wo else -1
            dm, dn, do = [99 if x < 0 else x for x in (dm, dn, do)]
            cnt['classified'] += 1
            if dm < dn and dm < do:
                cnt['M_rev' if rev else 'M_fwd'] += 1
            elif dn < dm and dn < do:
                cnt['N'] += 1
        out.append(dict(n_pos=npos, m_pos=mpos, **cnt))
    support = any(e.get('M_fwd', 0) >= prm.acc_min_strand and e.get('M_rev', 0) >= prm.acc_min_strand and
                  (e.get('M_fwd', 0) + e.get('M_rev', 0)) >= prm.acc_min_frac * max(1, e.get('classified', 0))
                  for e in out)
    contra = max([e.get('N', 0) for e in out] or [0])
    accepted = bool(ev) and support and contra <= prm.acc_max_contra
    return dict(mos=m_name, displaced=n_name, n_events=len(ev), events=out[:50], support=support,
                max_named_unique_reads=contra, accepted=accepted)


def is_novel(l):
    return any(t in l for t in ('*MOS_', '*PG', '*ORC_'))


def typed(slot):
    named = [l for l in slot if not is_novel(l)]
    novel = [l for l in slot if is_novel(l)]
    return named if (named and novel) else list(slot)


def run(pass1_dir: Path, pass1_ref: Path, pass2_table: Path, additions: Path, out: Path, prm):
    p1 = C.read_t1k(pass1_dir / 'calls_genotype.tsv')
    p2 = C.read_t1k(pass2_table)
    mos = {n: s for n, i, s in C.read_fasta_headers_and_seqs(additions)}
    need = set()
    for g in C.GENES:
        for sl in p1.get(g, []) + p2.get(g, []):
            need |= {l for l in sl if not is_novel(l)}
    named = {n: s for n, i, s in C.read_fasta_headers_and_seqs(pass1_ref) if n in need or n.split()[0] in need}
    seq_of = lambda l: mos.get(l) or named.get(l)
    res = {}
    for g in C.GENES:
        slots = p2.get(g, [[], []])
        for i, sl in enumerate(slots):
            t = typed(sl)
            if not t or not all('*MOS_' in l for l in t):
                continue
            m = t[0]
            cands = [s[0] for s in p1.get(g, []) if s and not is_novel(s[0])]
            if not cands:
                res.setdefault(g, []).append(dict(slot=i, mos=m, accepted=True, reason='no_named_pass1_allele'))
                continue
            n = min(cands, key=lambda x: edlib.align(seq_of(m), seq_of(x), mode='NW')['editDistance'])
            o_slot = typed(slots[1 - i]) if len(slots) == 2 else []
            o_seq = seq_of(o_slot[0]) if o_slot else None
            others = {k: v for k, v in named.items() if C.gene_of_name(k) == g}
            d = decide(g, m, seq_of(m), n, seq_of(n), o_seq, pass1_dir / 'calls_candidate_1.fq',
                       pass1_dir / 'calls_candidate_2.fq', others, out / 'acc_work', prm)
            d['slot'] = i
            res.setdefault(g, []).append(d)
    C.save_json(out / 'ACCEPT.json', dict(rule='read_evidence_acceptance_v1', params=dict(
        acc_min_strand=prm.acc_min_strand, acc_min_frac=prm.acc_min_frac, acc_max_contra=prm.acc_max_contra),
        genes=res))
    try:
        (out / 'acc_work').rmdir()
    except OSError:
        pass
    return res
