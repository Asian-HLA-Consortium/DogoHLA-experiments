"""Toy unit tests for the proposer (run on the login node; < 1 minute, no real data).

Each toy builds a small aligned panel, a diploid truth, simulates error-free paired reads (plus
optional injected errors), projects them through evidence.build_evidence exactly as real reads are
projected, and runs proposer.propose.
"""
import random
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import evidence as EV
import panel_select as PS
import proposer as PR

L = 2400
RNG_SEED = 7


def rand_seq(n, rng):
    return ''.join(rng.choice('ACGT') for _ in range(n))


def mutate(s, positions, rng, avoid=()):
    s = list(s)
    for p in positions:
        alts = [b for b in 'ACGT' if b != s[p] and b not in avoid]
        s[p] = rng.choice(alts)
    return ''.join(s)


def make_panel():
    """Rows (ungapped, no indels except row E): A base; B differs at clustered + spread sites;
    C, D partners; E has a 6-bp deletion so the MSA has gaps."""
    rng = random.Random(RNG_SEED)
    A = rand_seq(L, rng)
    # B: variants at 300,320 | 1000-1150 nothing | 1210,1230 | 1800, 2000 ...
    B = mutate(A, [300, 320, 600, 1140, 1210, 1280, 1800, 2000, 2200], rng)
    Cp = mutate(A, [150, 450, 700, 1500, 1650, 2100, 2300], rng)
    D = mutate(A, [250, 520, 900, 1400, 1950], rng)
    E = A[:1600] + '------' + A[1606:]
    names = ['A', 'B', 'C', 'D', 'E']
    return names, [A, B, Cp, D, E]


def sim_reads(haps, rng, cov=15, rl=150, frag=(300, 400), errors=()):
    """haps: list of (seq, cols) where cols[i] is the MSA column of base i. errors: list of
    (hap_index, base_position, new_base, n_reads) applied to that many reads covering it."""
    reads = []
    q = 0
    err = list(errors)
    for h, (seq, cols) in enumerate(haps):
        n_frag = int(cov * len(seq) / (2 * rl))
        for _ in range(n_frag):
            fl = rng.randint(*frag)
            st = rng.randint(0, max(0, len(seq) - fl))
            r1 = (st, min(len(seq), st + rl))
            r2 = (max(0, st + fl - rl), min(len(seq), st + fl))
            q += 1
            for mate, (a, b) in ((1, r1), (2, r2)):
                s = list(seq[a:b])
                for e in err:
                    if e[0] == h and a <= e[1] < b and e[3] > 0:
                        s[e[1] - a] = e[2]
                        e_list = list(e)
                        e_list[3] -= 1
                        err[err.index(e)] = tuple(e_list)
                        break
                reads.append(dict(qname=f'q{q}', seq=''.join(s), reverse=(mate == 2),
                                  c0=int(cols[a]), c1=int(cols[b - 1]) + 1, prefer_row=None))
    return reads


def path_hap(aligned_rows, segs):
    """Build an ungapped haplotype and its base->column map from (row, c0, c1) segments."""
    seq, cols = [], []
    for r, a, b in segs:
        for c in range(a, b):
            ch = aligned_rows[r][c]
            if ch != '-':
                seq.append(ch)
                cols.append(c)
    return ''.join(seq), np.array(cols)


def run_case(hap_segs, pass1, errors=(), private=None, prm=None, panel_rows=None, rng_seed=1):
    names, rows = make_panel()
    if panel_rows is not None:
        names = [names[i] for i in panel_rows]
        rows = [rows[i] for i in panel_rows]
    panel = PR.Panel(names, rows)
    haps = []
    for segs in hap_segs:
        seq, cols = path_hap(rows, segs)
        haps.append((seq, cols))
    truth = [h[0] for h in haps]
    if private is not None:   # (hap index, base position, new base)
        h, p, b = private
        s = list(haps[h][0])
        s[p] = b
        haps[h] = (''.join(s), haps[h][1])
        truth[h] = haps[h][0]
    rng = random.Random(rng_seed)
    reads = sim_reads(haps, rng, errors=errors)
    ev = EV.build_evidence(panel, reads)
    prm = prm or PR.Params()
    row_exons = {i: [(0, 100), (200, 300), (1000, 1200), (2300, len(panel.seq[i]) - 5)]
                 for i in range(panel.R)}
    rec = PR.propose(panel, ev, pass1, prm, set(panel.seq), row_exons, 'TOY')
    props = [c['_sequence'] for c in rec['candidates'] if '_sequence' in c]
    return rec, props, truth


A, B, Cc, D, E = range(5)


def test_exact_homozygote_no_proposal():
    rec, props, _ = run_case([[(A, 0, L)], [(A, 0, L)]], [[A], [A]])
    assert props == [] and rec['decision'] != 'proposed', rec['decision']


def test_exact_heterozygote_no_proposal():
    rec, props, _ = run_case([[(A, 0, L)], [(B, 0, L)]], [[A], [B]])
    assert props == [], rec['decision']


def test_exact_heterozygote_wrong_pass1_no_mosaic():
    # pass 1 called the wrong partner; the right answer is a zero-switch pair -> no proposal
    rec, props, _ = run_case([[(A, 0, L)], [(B, 0, L)]], [[A], [D]])
    assert props == [], rec['decision']


def test_mosaic1_linked_switch():
    # switch between B-variants 600 and 1210? use the 300/320 cluster: A[:310] is too short (<150?)
    # switch at 450..: A up to col 450, B after; A/B differ at 600 (B side) and 320,300 (A side):
    # 300..600 spans 300 bp -> linked by a fragment; partner C
    hap1 = [(A, 0, 450), (B, 450, L)]
    rec, props, truth = run_case([hap1, [(Cc, 0, L)]], [[A], [Cc]])
    assert rec['decision'] == 'proposed', rec
    assert truth[0] in props


def test_mosaic1_unlinked_switch_partner_is_source():
    # A/B agree on 601..1139 (539 bp > any fragment); switch at col 900
    hap1 = [(A, 0, 900), (B, 900, L)]
    rec, props, truth = run_case([hap1, [(B, 0, L)]], [[A], [B]])
    assert rec['decision'] == 'proposed', rec
    assert props == [truth[0]], (props, rec['solutions'][:2])


def test_mosaic1_unlinked_switch_partner_differs_parsimony():
    # Partner C differs from both sources: parsimony (1 switch vs 2) picks the true pair uniquely.
    hap1 = [(A, 0, 900), (B, 900, L)]
    rec, props, truth = run_case([hap1, [(Cc, 0, L)]], [[A], [Cc]])
    assert rec['decision'] == 'proposed', rec
    assert truth[0] in props
    assert rec['n_tied_solutions'] == 1


def test_mosaic2_single_variant_swap():
    # A with only B's allele at 1800: segments A[0:1700) B[1700:1900) A[1900:): A/B agree at 1700 and
    # 1900 and the B segment holds only the 1800 variant (B's others at 1280 and 2000 are outside).
    hap1 = [(A, 0, 1700), (B, 1700, 1900), (A, 1900, L)]
    rec, props, truth = run_case([hap1, [(Cc, 0, L)]], [[A], [Cc]])
    assert rec['decision'] == 'proposed', rec
    assert truth[0] in props


def test_private_snv_no_proposal():
    # base at 1000 changed to a character no panel row has there
    names, rows = make_panel()
    col = rows[A][1000]
    used = {r[1000] for r in rows}
    new = next(b for b in 'ACGT' if b not in used)
    rec, props, _ = run_case([[(A, 0, L)], [(Cc, 0, L)]], [[A], [Cc]], private=(0, 1000, new))
    assert props == [], rec['decision']


def test_sequencing_error_column_no_proposal():
    # one read carries B's allele at 1800 on haplotype A: unsupported -> no proposal
    names, rows = make_panel()
    rec, props, _ = run_case([[(A, 0, L)], [(Cc, 0, L)]], [[A], [Cc]],
                             errors=[(0, 1800, rows[B][1800], 1)])
    assert props == [], rec['decision']


def test_segment_shorter_than_150_no_proposal():
    # B's allele at 1210 only: a B segment covering 1210 but neither B variant 1140 nor 1280 holds at
    # most 139 bp < 150: no legal mosaic explains it.
    names, rows = make_panel()
    rec, props, _ = run_case([[(A, 0, L)], [(Cc, 0, L)]], [[A], [Cc]],
                             private=(0, 1210, rows[B][1210]))
    assert props == [], rec['decision']


def test_never_proposes_reference_sequence():
    # truth is exactly D (in the reference) but pass 1 called A: any proposal equal to D is refused
    rec, props, _ = run_case([[(D, 0, L)], [(Cc, 0, L)]], [[A], [Cc]])
    names, rows = make_panel()
    assert rows[D] not in props


def test_viterbi_respects_min_segment():
    names, rows = make_panel()
    panel = PR.Panel(names, rows)
    rng = np.random.default_rng(0)
    T = rng.random((panel.C, PR.NCODE)).astype(np.float32)
    cost, path = PR.viterbi(panel, T, PR.Params(), set())
    assert PR.seg_lengths_ok(panel, path, 150)
    for (r1, s1, e1), (r2, s2, e2) in zip(path.segs, path.segs[1:]):
        assert e1 == s2 and panel.P[r1, s2] == panel.P[r2, s2] < 4


def test_leakage_guard():
    row_is_ipd = [True, True, False, False]
    carriers = [set(), set(), {'D1'}, {'D2', 'X9'}]
    excl = {'D1', 'D1father'}
    PS.check_no_leakage([0, 1, 3], row_is_ipd, carriers, excl, 'D1', set(), 'D')
    with pytest.raises(PS.LeakageError):
        PS.check_no_leakage([0, 1, 2], row_is_ipd, carriers, excl, 'D1', set(), 'D')
    with pytest.raises(PS.LeakageError):
        PS.check_no_leakage([0, 3], row_is_ipd, carriers, excl, 'D1', set(), 'C')
    with pytest.raises(PS.LeakageError):
        PS.check_no_leakage([0, 1], row_is_ipd, carriers, excl, 'D1', {1}, 'C')
    with pytest.raises(PS.LeakageError):
        PS.check_no_leakage([0], row_is_ipd, carriers, {'X'}, 'D1', set(), 'C')


def test_integration_bwa_placement(tmp_path):
    """Real read path on a toy: reads (with flanks beyond the alleles) -> bwa mem against the
    pass-1 alleles -> clip trimming -> windows -> projection -> proposal."""
    import numpy as np
    import run_donor as RD
    names, rows = make_panel()
    rng = random.Random(3)
    flank_l, flank_r = rand_seq(300, rng), rand_seq(300, rng)
    hap1 = path_hap(rows, [(A, 0, 900), (B, 900, L)])[0]
    hap2 = path_hap(rows, [(Cc, 0, L)])[0]
    fq1, fq2 = tmp_path / 'r1.fq', tmp_path / 'r2.fq'
    comp = str.maketrans('ACGT', 'TGCA')
    with open(fq1, 'w') as f1, open(fq2, 'w') as f2:
        q = 0
        for h in (flank_l + hap1 + flank_r, flank_l + hap2 + flank_r):
            for _ in range(int(15 * len(h) / 300)):
                fl = rng.randint(300, 400)
                st = rng.randint(0, len(h) - fl)
                a = h[st:st + 150]
                b = h[st + fl - 150:st + fl].translate(comp)[::-1]
                q += 1
                f1.write(f'@p{q}/1\n{a}\n+\n{"I" * 150}\n')
                f2.write(f'@p{q}/2\n{b}\n+\n{"I" * 150}\n')
    panel = PR.Panel(names, rows)
    target = {'HLA-TOY*A': panel.seq[A], 'HLA-TOY*C': panel.seq[Cc]}
    old = RD.C.GENES
    RD.C.GENES = ('TOY',)
    try:
        reads, st = RD.map_reads(target, fq1, fq2, tmp_path / 'map', 2)
    finally:
        RD.C.GENES = old
    rowof = {'HLA-TOY*A': A, 'HLA-TOY*C': Cc}
    rd = []
    for x in reads['TOY']:
        r = rowof[x['contig']]
        bcols = np.nonzero(panel.isbase[r])[0]
        c0, c1, tl, tr = EV.window_from_alignment(bcols, x['ref_start'], x['ref_end'], x['lclip'],
                                                 x['rclip'], x['contig_len'])
        seq = x['seq'][tl:len(x['seq']) - tr] if tr else x['seq'][tl:]
        if len(seq) >= 30:
            rd.append(dict(qname=x['qname'], seq=seq, reverse=x['reverse'], c0=c0, c1=c1, prefer_row=r))
    ev = EV.build_evidence(panel, rd)
    row_exons = {i: [(0, 100), (200, 300), (1000, 1200), (2300, len(panel.seq[i]) - 5)]
                 for i in range(panel.R)}
    rec = PR.propose(panel, ev, [[A], [Cc]], PR.Params(), set(panel.seq), row_exons, 'TOY')
    props = [c['_sequence'] for c in rec['candidates'] if '_sequence' in c]
    assert ev.stats.get('projected', 0) > 0.9 * len(rd)
    assert rec['decision'] == 'proposed' and props == [hap1], (rec['decision'], ev.stats)
