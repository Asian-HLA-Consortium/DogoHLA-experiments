"""Fewest-switch recombinant proposer (Li-Stephens style copying model with read evidence).

Pure algorithm: no file I/O, no truth. Inputs are
  * panel: aligned rows (MSA strings) of the allowed panel P, with names,
  * Evidence: per-MSA-column observed character counts by strand plus per-fragment projections
    (built by evidence.py from the donor's reads),
  * the pass-1 called alleles as panel row indices.

Model
-----
Column evidence. At column c with depth d_c and count n_c(x) of character x (A,C,G,T,gap,other):
  supported(c,x)   = (n_fwd >= min_strand and n_rev >= min_strand)
                     or (n >= min_reads_frac and n/d >= min_frac)
  contradicted(c,x)= d_c >= min_depth_contra and n_c(x) < max_contra_reads and not supported
A diploid pair of characters (y, x) at c costs
  miss(c) = number of supported characters not in {y, x}  +  contra(y) + contra(x).
A pair of paths through panel rows costs the sum over columns plus lambda per switch; every
segment holds >= min_seg_bp ungapped bases; a switch from row a to row b at column s requires
a[s] == b[s] (a shared base, so the result is a valid graph walk). A tiny per-column prior
(prior_eps) keeps paths on T1K's pass-1 alleles when evidence is tied.

Search. Fix one haplotype, solve the other exactly by a segment Viterbi over rows (with the
minimum-length constraint), alternate until stable; several starts (pass-1 alleles). Phase
alternatives with equal column cost are separated by read/mate linkage; remaining ties are all
emitted (up to n_max). A proposal is made only if (i) the pass-1 pair leaves column evidence
unexplained, (ii) the best solution has >= 1 switch, and (iii) it beats the best zero-switch pair
by >= min_gain evidence units. Sequences already in the reference are never proposed.
"""
from __future__ import annotations

import dataclasses
import hashlib
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np

CODES = {'A': 0, 'C': 1, 'G': 2, 'T': 3, '-': 4}
OTHER = 5
OUTSIDE = 6   # before a row's first / after its last base: no sequence, not a deletion (wildcard)
NCODE = 7
INF = np.float32(1e9)


def encode(s: str) -> np.ndarray:
    lut = np.full(256, OTHER, dtype=np.uint8)
    for k, v in CODES.items():
        lut[ord(k)] = v
        lut[ord(k.lower())] = v
    return lut[np.frombuffer(s.encode(), dtype=np.uint8)]


@dataclasses.dataclass
class Params:
    min_strand: int = 2           # reads per strand for 'supported'
    min_frac: float = 0.2         # alternative: fraction of depth ...
    min_reads_frac: int = 3       # ... with at least this many reads
    min_depth_contra: int = 8     # depth needed before an absent character counts against a path
    max_contra_reads: int = 2     # a character seen in fewer reads than this is 'absent'
    lam: float = 0.1              # cost per switch (evidence units)
    w_link: float = 0.05          # cost per read-fragment phase mismatch (linkage)
    prior_eps: float = 1e-6       # per-column cost of leaving the pass-1 alleles (tie-break only)
    min_seg_bp: int = 150
    min_gain: float = 1.0         # mosaic must beat the best zero-switch pair by this much
    max_iter: int = 6
    n_max: int = 4
    max_read_edit: int = 10 ** 9  # reads farther than this from every panel row are off-target
    max_switches: int = 10 ** 9   # per haplotype
    max_residual: float = 1e9     # mosaic solution may leave at most this much column evidence
    skip_outside_cds_only: bool = False  # never propose a candidate equal to a reference allele
                                         # between first-exon start and last-exon end
    acc_min_strand: int = 2       # acceptance after pass 2: MOS-unique reads per strand at >=1 event
    acc_min_frac: float = 0.2     # ... and as a fraction of the reads classified at that event
    acc_max_contra: int = 1       # at no event more than this many reads uniquely favour the named allele

    def as_dict(self):
        return dataclasses.asdict(self)


@dataclasses.dataclass
class Evidence:
    ncol: int
    counts: np.ndarray            # (ncol, NCODE, 2) int32: [col, code, strand]
    frag_cols: List[np.ndarray]   # per fragment: columns observed (sorted)
    frag_codes: List[np.ndarray]  # per fragment: codes at those columns
    overflow: np.ndarray          # (ncol,) reads with unplaceable insertions (non-panel evidence)
    stats: Dict = dataclasses.field(default_factory=dict)


@dataclasses.dataclass
class Path:
    segs: List[Tuple[int, int, int]]   # (row, col_start, col_end_exclusive)

    @property
    def switches(self) -> int:
        return len(self.segs) - 1

    def codes(self, P: np.ndarray) -> np.ndarray:
        return np.concatenate([P[r, a:b] for r, a, b in self.segs])

    def key(self):
        return tuple(self.segs)


# ----------------------------------------------------------------------------- evidence -> costs

def flat_fragments(ev: Evidence):
    if not hasattr(ev, '_flat'):
        if ev.frag_cols:
            cols = np.concatenate(ev.frag_cols).astype(np.int64)
            codes = np.concatenate(ev.frag_codes).astype(np.int64)
            fid = np.concatenate([np.full(len(c), i) for i, c in enumerate(ev.frag_cols)])
        else:
            cols = codes = fid = np.zeros(0, dtype=np.int64)
        ev._flat = (cols, codes, fid, len(ev.frag_cols))
    return ev._flat


def residual_counts(ev: Evidence, sup, y_codes):
    """Fragments that contradict the fixed haplotype at a supported character must come from the
    other haplotype: return R[c, x] = number of such fragments showing supported x at c."""
    cols, codes, fid, nf = flat_fragments(ev)
    R = np.zeros((sup.shape[0], NCODE), dtype=np.float32)
    if nf == 0:
        return R
    s = sup[cols, codes]
    mism = (y_codes[cols] != codes) & s & (y_codes[cols] != OUTSIDE)
    m1 = np.bincount(fid, weights=mism, minlength=nf)
    keep = (m1[fid] > 0) & s
    np.add.at(R, (cols[keep], codes[keep]), 1)
    return R


def column_tables(ev: Evidence, prm: Params):
    n = ev.counts.sum(axis=2)                      # (ncol, NCODE)
    depth = n.sum(axis=1)
    fwd, rev = ev.counts[:, :, 0], ev.counts[:, :, 1]
    with np.errstate(divide='ignore', invalid='ignore'):
        frac = np.where(depth[:, None] > 0, n / np.maximum(depth[:, None], 1), 0.0)
    sup = ((fwd >= prm.min_strand) & (rev >= prm.min_strand)) | \
          ((n >= prm.min_reads_frac) & (frac >= prm.min_frac))
    sup[:, OTHER] = False   # 'other' characters are never a panel allele; kept as overflow-like
    sup[:, OUTSIDE] = False
    contra = (depth[:, None] >= prm.min_depth_contra) & (n < prm.max_contra_reads) & ~sup
    contra[:, OUTSIDE] = False
    return sup, contra, depth


def pair_cost_table(sup, contra, y_codes: np.ndarray, ev: Evidence = None, w_link: float = 0.0) -> np.ndarray:
    """T[c, x] = cost at column c when the fixed haplotype has y_codes[c] and the other has x;
    plus w_link per residual fragment (one contradicting the fixed haplotype) not showing x."""
    ncol = sup.shape[0]
    T = np.zeros((ncol, NCODE), dtype=np.float32)
    nsup = sup.sum(axis=1).astype(np.float32)
    idx = np.arange(ncol)
    y_sup = sup[idx, y_codes].astype(np.float32)
    y_out = y_codes == OUTSIDE
    for x in range(NCODE):
        x_sup = sup[:, x].astype(np.float32)
        both = np.where(y_codes == x, y_sup, y_sup + x_sup)
        miss = np.where(y_out | (x == OUTSIDE), 0.0, nsup - both)
        T[:, x] = miss + contra[:, x] + contra[idx, y_codes]
    if ev is not None and w_link:
        R = residual_counts(ev, sup, y_codes)
        T += w_link * (R.sum(axis=1, keepdims=True) - R)
    return T


def pair_evidence_cost(sup, contra, a: np.ndarray, b: np.ndarray) -> float:
    nsup = sup.sum(axis=1)
    idx = np.arange(sup.shape[0])
    sa, sb = sup[idx, a], sup[idx, b]
    covered = np.where(a == b, sa, sa.astype(int) + sb.astype(int))
    miss = np.where((a == OUTSIDE) | (b == OUTSIDE), 0, nsup - covered)
    return float(miss.sum() + contra[idx, a].sum() + contra[idx, b].sum())


def unexplained_columns(sup, contra, a, b) -> np.ndarray:
    idx = np.arange(sup.shape[0])
    s = sup.copy()
    s[idx, a] = False
    s[idx, b] = False
    s[(a == OUTSIDE) | (b == OUTSIDE)] = False
    return np.nonzero(s.any(axis=1) | contra[idx, a] | contra[idx, b])[0]


# ----------------------------------------------------------------------------- Viterbi

class Panel:
    def __init__(self, names: Sequence[str], aligned: Sequence[str]):
        assert len(names) == len(aligned) and len(set(map(len, aligned))) == 1
        self.names = list(names)
        self.P = np.stack([encode(s) for s in aligned])      # (R, C) uint8
        self.R, self.C = self.P.shape
        for r in range(self.R):
            nz = np.nonzero(self.P[r] != 4)[0]
            if len(nz):
                self.P[r, :nz[0]] = OUTSIDE
                self.P[r, nz[-1] + 1:] = OUTSIDE
        isbase = self.P < 4
        self.isbase = isbase
        self.Pbase = np.zeros((self.R, self.C + 1), dtype=np.int32)
        np.cumsum(isbase, axis=1, out=self.Pbase[:, 1:])
        self.seq = [''.join('ACGT'[c] if c < 4 else ('N' if c == OTHER else '') for c in row)
                    for row in self.P]
        self._tcache = {}

    def thresholds(self, L: int) -> np.ndarray:
        """t[r, s] = largest s' with Pbase[r, s'] <= Pbase[r, s] - L (or -1)."""
        if L not in self._tcache:
            t = np.empty((self.R, self.C + 1), dtype=np.int32)
            for r in range(self.R):
                t[r] = np.searchsorted(self.Pbase[r], self.Pbase[r] - L, side='right') - 1
            self._tcache[L] = t
        return self._tcache[L]


def viterbi(panel: Panel, T: np.ndarray, prm: Params, prior_rows: Set[int],
            allow_switch: bool = True) -> Tuple[float, Path]:
    """Best single path through panel rows for the column cost table T (C x NCODE)."""
    R, C, P = panel.R, panel.C, panel.P
    K = T[np.arange(C)[None, :], P]                               # (R, C)
    if prm.prior_eps and prior_rows:
        pr = np.full(R, prm.prior_eps, dtype=np.float32)
        pr[list(prior_rows)] = 0
        K = K + pr[:, None]
    Cum = np.zeros((R, C + 1), dtype=np.float64)
    np.cumsum(K, axis=1, out=Cum[:, 1:])
    if not allow_switch:
        tot = Cum[:, C]
        r = int(np.argmin(tot))
        return float(tot[r]), Path([(r, 0, C)])
    t = panel.thresholds(prm.min_seg_bp)
    G = np.full((R, C + 1), np.inf)
    PM = np.full((R, C + 1), np.inf)
    G[:, 0] = 0.0
    PM[:, 0] = 0.0
    argA = np.full((C + 1, 4), -1, dtype=np.int64)
    rows = np.arange(R)
    for s in range(1, C):
        ts = t[:, s]
        ok = ts >= 0
        H = np.full(R, np.inf)
        H[ok] = Cum[ok, s] + PM[rows[ok], ts[ok]]
        ch = P[:, s]
        Fs = np.full(R, np.inf)
        if np.isfinite(H).any():
            for x in range(4):
                m = ch == x
                if not m.any():
                    continue
                hx = np.where(m, H, np.inf)
                a = int(np.argmin(hx))
                if np.isfinite(hx[a]):
                    argA[s, x] = a
                    Fs[m] = hx[a] + prm.lam
        G[:, s] = Fs - Cum[:, s]
        PM[:, s] = np.minimum(PM[:, s - 1], G[:, s])
    PM[:, C] = PM[:, C - 1]
    # close the last segment (length >= L unless the path is a single segment)
    tC = t[:, C]
    best_end = np.full(R, np.inf)
    okC = tC >= 0
    best_end[okC] = Cum[okC, C] + PM[rows[okC], tC[okC]]
    single = Cum[:, C]                          # G[:,0] == 0: path with no switch
    total = np.minimum(best_end, single)
    r = int(np.argmin(total))
    cost = float(total[r])
    # traceback
    segs = []
    end = C
    if single[r] <= best_end[r]:
        return cost, Path([(r, 0, C)])
    cur_r, limit = r, t[r, C]
    while True:
        g = G[cur_r, :limit + 1]
        s0 = int(np.argmin(g))
        segs.append((cur_r, s0, end))
        if s0 == 0:
            break
        x = int(P[cur_r, s0])
        a = int(argA[s0, x])
        end = s0
        cur_r, limit = a, t[a, s0]
    segs.reverse()
    return cost, Path(segs)


# ----------------------------------------------------------------------------- linkage

def linkage_cost(ev: Evidence, a: np.ndarray, b: np.ndarray, sup=None) -> int:
    """Sum over fragments of min(mismatches to a, mismatches to b), counted at columns where
    a != b and the fragment's character is supported (sequencing errors are ignored)."""
    cols, codes, fid, nf = flat_fragments(ev)
    if nf == 0:
        return 0
    m = (a[cols] != b[cols]) & (a[cols] != OUTSIDE) & (b[cols] != OUTSIDE)
    if sup is not None:
        m &= sup[cols, codes]
    ma = np.bincount(fid, weights=m & (a[cols] != codes), minlength=nf)
    mb = np.bincount(fid, weights=m & (b[cols] != codes), minlength=nf)
    return int(np.minimum(ma, mb).sum())


def seg_lengths_ok(panel: Panel, path: Path, L: int) -> bool:
    if path.switches == 0:
        return True
    return all(panel.Pbase[r, e] - panel.Pbase[r, s] >= L for r, s, e in path.segs)


def merge_path(segs):
    out = []
    for r, s, e in segs:
        if s >= e:
            continue
        if out and out[-1][0] == r and out[-1][2] == s:
            out[-1] = (r, out[-1][1], e)
        else:
            out.append((r, s, e))
    return Path(out)


def split_at(path: Path, col: int):
    left, right = [], []
    for r, s, e in path.segs:
        if e <= col:
            left.append((r, s, e))
        elif s >= col:
            right.append((r, s, e))
        else:
            left.append((r, s, col))
            right.append((r, col, e))
    return left, right


def phase_swaps(panel: Panel, h1: Path, h2: Path, prm: Params):
    """Alternatives that exchange tails at a column where both haplotypes share a base and the
    switch is legal; column evidence is identical by construction, linkage decides."""
    a, b = h1.codes(panel.P), h2.codes(panel.P)
    sw_cols = sorted({s for p in (h1, h2) for _, s, _ in p.segs[1:]})
    out = []
    for col in sw_cols:
        if a[col] != b[col] or a[col] >= 4:
            continue
        l1, r1 = split_at(h1, col)
        l2, r2 = split_at(h2, col)
        n1, n2 = merge_path(l1 + r2), merge_path(l2 + r1)
        if seg_lengths_ok(panel, n1, prm.min_seg_bp) and seg_lengths_ok(panel, n2, prm.min_seg_bp):
            out.append((n1, n2))
    return out


# ----------------------------------------------------------------------------- driver

@dataclasses.dataclass
class Solution:
    h1: Path
    h2: Path
    evidence_cost: float
    switches: int
    linkage: int = 0


def solve_pair(panel: Panel, sup, contra, init: Path, prm: Params, prior_rows: Set[int],
               allow_switch=True, ev: Evidence = None) -> Solution:
    P = panel.P
    h1 = init
    h2 = None
    seen = set()
    for _ in range(prm.max_iter):
        T = pair_cost_table(sup, contra, h1.codes(P), ev, prm.w_link)
        _, h2n = viterbi(panel, T, prm, prior_rows, allow_switch)
        T = pair_cost_table(sup, contra, h2n.codes(P), ev, prm.w_link)
        _, h1n = viterbi(panel, T, prm, prior_rows, allow_switch)
        key = (h1n.key(), h2n.key())
        h1, h2 = h1n, h2n
        if key in seen:
            break
        seen.add(key)
    a, b = h1.codes(P), h2.codes(P)
    so = Solution(h1, h2, pair_evidence_cost(sup, contra, a, b), h1.switches + h2.switches)
    if ev is not None:
        so.linkage = linkage_cost(ev, a, b, sup)
    return so


def total(so: Solution, prm: Params) -> float:
    return so.evidence_cost + prm.lam * so.switches + prm.w_link * so.linkage


def path_exons(panel: Panel, path: Path, row_exons: Dict[int, List[Tuple[int, int]]]):
    """Exon intervals (0-based, end-exclusive) of a path's sequence, taken from the annotation of
    the segment row that holds each boundary column. Returns None if inconsistent."""
    P = panel.P
    nexon = None
    colint = {}
    for r, _, _ in path.segs:
        ex = row_exons.get(r)
        if not ex:
            return None
        if nexon is None:
            nexon = len(ex)
        if len(ex) != nexon:
            return None
        bcols = np.nonzero(panel.isbase[r])[0]
        colint[r] = [(int(bcols[s]), int(bcols[e - 1])) for s, e in ex]
    starts, ends = [], []
    for k in range(nexon):
        st = en = None
        for r, s, e in path.segs:
            cs, ce = colint[r][k]
            if st is None and s <= cs < e:
                st = cs
            if en is None and s <= ce < e:
                en = ce
        if st is None or en is None:
            return None
        starts.append(st)
        ends.append(en)
    isb = np.concatenate([panel.isbase[r, s:e] for r, s, e in path.segs])
    pre = np.concatenate([[0], np.cumsum(isb)])
    out = []
    for st, en in zip(starts, ends):
        if not isb[st] or not isb[en]:
            return None
        out.append((int(pre[st]), int(pre[en + 1])))
    for i in range(len(out)):
        if out[i][0] >= out[i][1] or (i and out[i][0] < out[i - 1][1]):
            return None
    return out


_CDS_SPAN = {}


def cds_span_set(panel: Panel, row_exons):
    k = id(panel)
    if k not in _CDS_SPAN:
        out = set()
        for r in range(panel.R):
            ex = row_exons.get(r)
            if ex:
                out.add(panel.seq[r][ex[0][0]:ex[-1][1]])
        _CDS_SPAN[k] = out
    return _CDS_SPAN[k]


def path_sequence(panel: Panel, path: Path) -> str:
    return ''.join(panel.seq[r][panel.Pbase[r, s]:panel.Pbase[r, e]] for r, s, e in path.segs)


def propose(panel: Panel, ev: Evidence, pass1_rows: Sequence[Sequence[int]], prm: Params,
            reference_seqs: Set[str], row_exons: Dict[int, List[Tuple[int, int]]],
            gene: str) -> Dict:
    """pass1_rows: per pass-1 slot, the panel rows of its (tied) labels, best first.
    Returns a record with the decision, candidates and diagnostics."""
    sup, contra, depth = column_tables(ev, prm)
    P = panel.P
    rec = dict(gene=gene, params=prm.as_dict(), n_panel_rows=panel.R, ncol=panel.C,
               mean_depth=float(depth.mean()), supported_columns=int((sup.sum(1) > 0).sum()),
               overflow_supported=int((ev.overflow >= prm.min_reads_frac).sum()),
               candidates=[], solutions=[], decision=None)
    slots = [list(s) for s in pass1_rows if s]
    if not slots:
        rec['decision'] = 'no_pass1_call'
        return rec
    if len(slots) == 1:
        slots = slots * 2
    prior = {r for s in slots for r in s[:3]}
    # pass-1 pair evidence (best over tied labels)
    p1 = min(((pair_evidence_cost(sup, contra, P[a], P[b]), a, b)
              for a in slots[0][:3] for b in slots[1][:3]), key=lambda x: x[0])
    rec['pass1_evidence_cost'] = p1[0]
    rec['pass1_unexplained_columns'] = int(len(unexplained_columns(sup, contra, P[p1[1]], P[p1[2]])))
    inits = []
    for s in slots:
        for r in s[:2]:
            if r not in inits:
                inits.append(r)
    zero = [solve_pair(panel, sup, contra, Path([(r, 0, panel.C)]), prm, prior, False, ev) for r in inits]
    # T1K's own pass-1 pair is always a zero-switch candidate
    p1sol = Solution(Path([(p1[1], 0, panel.C)]), Path([(p1[2], 0, panel.C)]), p1[0], 0)
    p1sol.linkage = linkage_cost(ev, P[p1[1]], P[p1[2]], sup)
    zero.append(p1sol)
    best0 = min(zero, key=lambda z: total(z, prm))
    rec['zero_switch_evidence_cost'] = best0.evidence_cost
    rec['zero_switch_linkage'] = best0.linkage
    rec['zero_switch_total'] = total(best0, prm)
    rec['zero_switch_pair'] = [panel.names[best0.h1.segs[0][0]], panel.names[best0.h2.segs[0][0]]]
    if p1[0] <= 0:
        rec['decision'] = 'pass1_explains_all_evidence'
        return rec
    starts = inits + [r for r in (best0.h1.segs[0][0], best0.h2.segs[0][0]) if r not in inits]
    sols = [solve_pair(panel, sup, contra, Path([(r, 0, panel.C)]), prm, prior, True, ev) for r in starts]
    # zero-switch refinement started from every row the switching solutions use
    extra = []
    for so in sols:
        for h in (so.h1, so.h2):
            for r, _, _ in h.segs:
                if r not in starts and r not in extra:
                    extra.append(r)
    for r in extra[:12]:
        z = solve_pair(panel, sup, contra, Path([(r, 0, panel.C)]), prm, prior, False, ev)
        if total(z, prm) < total(best0, prm):
            best0 = z
    rec['zero_switch_evidence_cost'] = best0.evidence_cost
    rec['zero_switch_linkage'] = best0.linkage
    rec['zero_switch_total'] = total(best0, prm)
    rec['zero_switch_pair'] = [panel.names[best0.h1.segs[0][0]], panel.names[best0.h2.segs[0][0]]]
    sols.append(best0)
    a0, b0 = best0.h1.codes(P), best0.h2.codes(P)
    gapfrac = (P == 4).mean(axis=0)
    diag = []
    for c in unexplained_columns(sup, contra, a0, b0)[:80]:
        diag.append(dict(col=int(c), depth=int(ev.counts[c].sum()),
                         counts=ev.counts[c].sum(axis=1).tolist(), pair=[int(a0[c]), int(b0[c])],
                         panel_gap_frac=round(float(gapfrac[c]), 3),
                         supported=np.nonzero(sup[c])[0].tolist()))
    rec['diag_zero_unexplained'] = diag
    rec['diag_zero_unexplained_n'] = int(len(unexplained_columns(sup, contra, a0, b0)))
    # phase alternatives of every solution
    alts = []
    for so in sols:
        for n1, n2 in phase_swaps(panel, so.h1, so.h2, prm):
            a, b = n1.codes(P), n2.codes(P)
            alts.append(Solution(n1, n2, pair_evidence_cost(sup, contra, a, b),
                                 n1.switches + n2.switches))
    pool = {}
    for so in sols + alts:
        k = frozenset([path_sequence(panel, so.h1), path_sequence(panel, so.h2)])
        if k not in pool:
            so.linkage = linkage_cost(ev, so.h1.codes(P), so.h2.codes(P), sup)
            pool[k] = so
    ranked = sorted(pool.values(), key=lambda z: (total(z, prm), z.switches))
    top = ranked[0]
    tied = [z for z in ranked if abs(total(z, prm) - total(top, prm)) < 1e-9]
    for z in ranked[:6]:
        rec['solutions'].append(dict(
            evidence_cost=z.evidence_cost, switches=z.switches, linkage=z.linkage,
            h1=[(panel.names[r], int(s), int(e)) for r, s, e in z.h1.segs],
            h2=[(panel.names[r], int(s), int(e)) for r, s, e in z.h2.segs]))
    rec['best_evidence_cost'] = top.evidence_cost
    rec['best_switches'] = top.switches
    rec['n_tied_solutions'] = len(tied)
    if top.switches == 0:
        rec['decision'] = 'zero_switch_solution'
        return rec
    rec['best_linkage'] = top.linkage
    rec['best_total'] = total(top, prm)
    if max(top.h1.switches, top.h2.switches) > prm.max_switches:
        rec['decision'] = 'too_many_switches'
        return rec
    if top.evidence_cost > prm.max_residual:
        rec['decision'] = 'residual_evidence_above_max'
        return rec
    if (best0.evidence_cost - top.evidence_cost < prm.min_gain) or not (total(top, prm) < total(best0, prm)):
        rec['decision'] = 'mosaic_gain_below_min_gain'
        return rec
    cands, seen = [], set()
    for z in tied:
        for h in (z.h1, z.h2):
            if h.switches == 0:
                continue
            seq = path_sequence(panel, h)
            if seq in seen:
                continue
            seen.add(seq)
            entry = dict(sources=[panel.names[r] for r, _, _ in h.segs],
                         segments=[(int(s), int(e)) for _, s, e in h.segs],
                         switch_columns=[int(s) for _, s, _ in h.segs[1:]],
                         length=len(seq), sha256=hashlib.sha256(seq.encode()).hexdigest())
            ex0 = path_exons(panel, h, row_exons)
            if seq in reference_seqs:
                entry['skipped'] = 'sequence_already_in_reference'
            elif prm.skip_outside_cds_only and ex0 is not None and \
                    seq[ex0[0][0]:ex0[-1][1]] in cds_span_set(panel, row_exons):
                entry['skipped'] = 'differs_from_a_reference_allele_only_outside_first_to_last_exon'
            else:
                ex = ex0
                if ex is None:
                    entry['skipped'] = 'exon_annotation_inconsistent'
                else:
                    entry['exons'] = ex
                    a, b = z.h1.codes(P), z.h2.codes(P)
                    entry['supporting_reads_at_unexplained_columns'] = int(sum(
                        ev.counts[c].sum() for c in unexplained_columns(
                            sup, contra, P[p1[1]], P[p1[2]])[:200]))
                    entry['_sequence'] = seq
            cands.append(entry)
    usable = [c for c in cands if 'skipped' not in c][:prm.n_max]
    rec['candidates'] = cands
    rec['n_proposed'] = len(usable)
    rec['decision'] = 'proposed' if usable else 'all_candidates_skipped'
    return rec
