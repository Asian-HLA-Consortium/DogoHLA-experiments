"""Read evidence for the proposer: place reads, project them onto MSA columns.

Placement. Reads are aligned (bwa mem, paired) to T1K's pass-1 called alleles of ALL genes in the
pass-1 table (so paralog reads go to their own gene) plus tied alternatives; a read belongs to the
gene of its primary alignment. Unmapped reads are counted, never silently dropped.

Projection. Within the MSA window around the placement, the read (whole read, end to end; soft
clips are kept and therefore count as mismatches -- except bases that overhang the END of the
target allele sequence, which are trimmed and counted as overhang) is compared with every distinct
panel-row substring of the window (edlib, semi-global 'HW'); the best row's alignment projects each
read base onto an MSA column. Read gaps against a row base become '-' observations; read bases
inserted relative to the row fill that row's gap columns left to right, and any surplus is counted
as an 'overflow' (non-panel sequence) at that column.
"""
from __future__ import annotations

import collections
import re
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

import proposer as PR

CIG = re.compile(r'(\d+)([=XIDM])')
QUANT = 16
PAD_BASES = 30


class WindowIndex:
    """Distinct panel-row substrings per quantized MSA window."""

    def __init__(self, panel: PR.Panel):
        self.panel = panel
        self.cache: Dict[Tuple[int, int], Dict[str, int]] = {}

    def subs(self, w0, w1) -> Dict[str, int]:
        k = (w0, w1)
        if k not in self.cache:
            d = {}
            pb = self.panel.Pbase
            for r in range(self.panel.R):
                s = self.panel.seq[r][pb[r, w0]:pb[r, w1]]
                if s and s not in d:
                    d[s] = r
            self.cache[k] = d
        return self.cache[k]


def quantize(c0, c1, ncol):
    return max(0, (c0 // QUANT) * QUANT), min(ncol, -(-c1 // QUANT) * QUANT)


def project(read: str, panel: PR.Panel, r: int, w0: int, w1: int):
    """Project read onto columns via row r restricted to [w0, w1). Returns (cols, codes, overflow
    columns, edit distance) or None."""
    import edlib
    pb = panel.Pbase
    sub = panel.seq[r][pb[r, w0]:pb[r, w1]]
    if not sub:
        return None
    res = edlib.align(read, sub, mode='HW', task='path')
    if res['editDistance'] < 0 or not res['locations']:
        return None
    t0 = res['locations'][0][0]
    rowcols = np.nonzero(panel.isbase[r, w0:w1])[0] + w0
    cols, codes, overflow = [], [], []
    tpos, rpos = t0, 0
    pending: List[str] = []
    started = False

    def flush_gap(prev_t, next_t):
        # columns strictly between the row bases prev_t and next_t
        if prev_t < 0 or next_t >= len(rowcols):
            return
        gcols = range(rowcols[prev_t] + 1, rowcols[next_t])
        ins = list(pending)
        for i, gc in enumerate(gcols):
            cols.append(gc)
            codes.append(PR.CODES.get(ins[i], PR.OTHER) if i < len(ins) else 4)
        if len(ins) > len(gcols):
            overflow.append(int(rowcols[prev_t]))
        pending.clear()

    for n, op in CIG.findall(res['cigar']):
        n = int(n)
        for _ in range(n):
            if op in '=XM':
                if started:
                    flush_gap(tpos - 1, tpos)
                cols.append(int(rowcols[tpos]))
                codes.append(PR.CODES.get(read[rpos], PR.OTHER))
                tpos += 1
                rpos += 1
                started = True
            elif op == 'D':
                if started:
                    flush_gap(tpos - 1, tpos)
                cols.append(int(rowcols[tpos]))
                codes.append(4)
                tpos += 1
                started = True
            elif op == 'I':
                if started:
                    pending.append(read[rpos])
                else:
                    overflow.append(int(rowcols[tpos]))   # leading insertion: unplaceable
                rpos += 1
    if pending:
        overflow.append(int(rowcols[tpos - 1]))
    return np.array(cols, dtype=np.int64), np.array(codes, dtype=np.uint8), overflow, res['editDistance']


def best_row(read: str, widx: WindowIndex, w0: int, w1: int, prefer_row: Optional[int]):
    """Best-matching panel row for the read. The placement allele's own row is used whenever it is
    among the best (so reads covering the same place are projected through the same gap layout);
    otherwise the representative row of the best distinct substring."""
    import edlib
    subs = widx.subs(w0, w1)
    panel = widx.panel
    pd = None
    if prefer_row is not None:
        pb = panel.Pbase
        ps = panel.seq[prefer_row][pb[prefer_row, w0]:pb[prefer_row, w1]]
        if ps:
            if read in ps:
                return prefer_row, 0
            pd = edlib.align(read, ps, mode='HW', task='distance')['editDistance']
    best, bestd = None, None
    for s, r in subs.items():
        if read in s:
            best, bestd = r, 0
            break
    if best is None:
        for s, r in subs.items():
            k = -1 if bestd is None else bestd
            d = edlib.align(read, s, mode='HW', task='distance', k=k)['editDistance']
            if d >= 0 and (bestd is None or d < bestd):
                best, bestd = r, d
    if pd is not None and pd >= 0 and (bestd is None or pd <= bestd):
        return prefer_row, pd
    return best, bestd


def build_evidence(panel: PR.Panel, reads: Sequence[Dict], max_read_edit: int = 10 ** 9) -> PR.Evidence:
    """reads: dicts with qname, seq (reference orientation), reverse, c0, c1 (MSA columns of the
    placement incl. kept clips, end-exclusive), prefer_row (row of the placement allele or None)."""
    C = panel.C
    counts = np.zeros((C, PR.NCODE, 2), dtype=np.int32)
    overflow = np.zeros(C, dtype=np.int32)
    widx = WindowIndex(panel)
    frag = collections.defaultdict(list)
    st = collections.Counter()
    eds = []
    for rd in reads:
        w0, w1 = quantize(rd['c0'], rd['c1'], C)
        r, d = best_row(rd['seq'], widx, w0, w1, rd.get('prefer_row'))
        if r is None:
            st['unprojectable'] += 1
            continue
        if d > max_read_edit:
            st['excluded_offtarget_edit_gt_max'] += 1
            continue
        pr = project(rd['seq'], panel, r, w0, w1)
        if pr is None:
            st['unprojectable'] += 1
            continue
        cols, codes, ovf, ed = pr
        st['projected'] += 1
        st['exact' if ed == 0 else 'inexact'] += 1
        eds.append(ed)
        np.add.at(counts, (cols, codes, np.full(len(cols), 1 if rd['reverse'] else 0)), 1)
        for c in ovf:
            overflow[c] += 1
        frag[rd['qname']].append((cols, codes))
    fc, fk = [], []
    for q, parts in frag.items():
        if len(parts) == 1:
            cols, codes = parts[0]
        else:
            m = {}
            bad = set()
            for cols, codes in parts:
                for c, k in zip(cols.tolist(), codes.tolist()):
                    if c in m and m[c] != k:
                        bad.add(c)
                    m[c] = k
            keys = sorted(c for c in m if c not in bad)
            cols = np.array(keys, dtype=np.int64)
            codes = np.array([m[c] for c in keys], dtype=np.uint8)
        fc.append(cols)
        fk.append(codes)
    st['fragments'] = len(fc)
    st['mean_edit_distance'] = float(np.mean(eds)) if eds else None
    return PR.Evidence(ncol=C, counts=counts, frag_cols=fc, frag_codes=fk, overflow=overflow,
                       stats=dict(st))


def window_from_alignment(panel_bcols: np.ndarray, ref_start: int, ref_end: int, lclip: int,
                          rclip: int, contig_len: int):
    """MSA column window [c0, c1) for an alignment on a contig whose row has base columns
    panel_bcols; clips that overhang the contig ends are trimmed. Returns (c0, c1, trim_left,
    trim_right)."""
    lo = ref_start - lclip
    hi = ref_end + rclip
    tl = max(0, -lo)
    tr = max(0, hi - contig_len)
    lo = max(0, lo - PAD_BASES)
    hi = min(contig_len, hi + PAD_BASES)
    return int(panel_bcols[lo]), int(panel_bcols[hi - 1]) + 1, tl, tr
