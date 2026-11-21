"""Part A: Variant atlas and SV classification from per-gene MSA and minimap2 PAF."""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

from panels import CLASS_I_GENES, CLASS_II_GENES


@dataclass(frozen=True)
class VariantEvent:
    gene: str
    col_start: int
    col_end: int
    ref_seq: str
    alt_seq: str
    var_type: str
    size_bp: int
    size_class: str
    sv_subclass: str
    repeat_period: Optional[int]
    location: str
    location_broad: str
    antigen_binding_exon: bool


def sequence_span_in_msa(aligned_seq: str) -> Tuple[int, int]:
    """Return (first_non_gap_col, last_non_gap_col) inclusive; (-1, -1) if all gap."""
    start = -1
    end = -1
    for i, ch in enumerate(aligned_seq):
        if ch != "-":
            if start == -1:
                start = i
            end = i
    return start, end


def build_reference_backbone(
    aligned_sequences: Sequence[str],
    preferred_ref_aligned: Optional[str] = None,
) -> str:
    """Build fixed MSA column reference state R[0..M-1] by majority over covered spans.

    Terminal gaps outside each sequence's [start, end] span are treated as unobserved.
    Ties are broken by preferred_ref_aligned if provided and covering that column.
    """
    if not aligned_sequences:
        raise ValueError("Empty aligned_sequences")
    m = len(aligned_sequences[0])
    spans = [sequence_span_in_msa(s) for s in aligned_sequences]
    pref_span = sequence_span_in_msa(preferred_ref_aligned) if preferred_ref_aligned else (-1, -1)

    backbone: List[str] = []
    for c in range(m):
        counts: Counter[str] = Counter()
        for s, (st, en) in zip(aligned_sequences, spans):
            if st != -1 and st <= c <= en:
                ch = s[c].upper()
                if ch in ("A", "C", "G", "T", "-"):
                    counts[ch] += 1
        if not counts:
            if preferred_ref_aligned and pref_span[0] <= c <= pref_span[1]:
                backbone.append(preferred_ref_aligned[c].upper())
            else:
                backbone.append("-")
            continue
        max_cnt = max(counts.values())
        best = [k for k, v in counts.items() if v == max_cnt]
        if (
            len(best) > 1
            and preferred_ref_aligned
            and pref_span[0] <= c <= pref_span[1]
            and preferred_ref_aligned[c].upper() in best
        ):
            backbone.append(preferred_ref_aligned[c].upper())
        else:
            # Deterministic order: bases before gap
            best.sort(key=lambda ch: (ch == "-", ch))
            backbone.append(best[0])
    return "".join(backbone)


def map_exons_to_msa_columns(
    aligned_seq: str,
    exons_0based_half_open: Sequence[Tuple[int, int]],
) -> List[Tuple[int, int]]:
    """Map 0-based half-open ungapped exon intervals [(s, e), ...] to MSA column intervals [(cs, ce), ...]."""
    base_to_col: List[int] = [i for i, ch in enumerate(aligned_seq) if ch != "-"]
    n = len(base_to_col)
    col_intervals: List[Tuple[int, int]] = []
    for s, e in exons_0based_half_open:
        if s < 0 or e > n or s >= e:
            raise ValueError(f"Invalid exon interval {(s, e)} for sequence of ungapped length {n}")
        col_intervals.append((base_to_col[s], base_to_col[e - 1] + 1))
    return col_intervals


def classify_column_location(
    gene: str,
    col_start: int,
    col_end: int,
    exon_col_intervals: Sequence[Tuple[int, int]],
) -> Tuple[str, str, bool]:
    """Classify an MSA column interval [col_start, col_end) into exon/intron/UTR and antigen-binding exon."""
    if not exon_col_intervals:
        return "unannotated", "unannotated", False

    # Effective probe interval (for pure insertion of 0 col width, use [col_start, col_start+1))
    a = col_start
    b = max(col_end, col_start + 1)

    overlapping_exons: List[int] = []
    for idx, (es, ee) in enumerate(exon_col_intervals, start=1):
        if a < ee and b > es:
            overlapping_exons.append(idx)

    if overlapping_exons:
        ex_num = overlapping_exons[0]
        loc = f"exon_{ex_num}"
        broad = "exon"
        if gene in CLASS_I_GENES:
            ag = any(k in (2, 3) for k in overlapping_exons)
        elif gene in CLASS_II_GENES:
            ag = any(k == 2 for k in overlapping_exons)
        else:
            ag = False
        return loc, broad, ag

    first_es = exon_col_intervals[0][0]
    last_ee = exon_col_intervals[-1][1]
    if b <= first_es:
        return "5_UTR", "UTR", False
    if a >= last_ee:
        return "3_UTR", "UTR", False

    for idx in range(len(exon_col_intervals) - 1):
        prev_end = exon_col_intervals[idx][1]
        next_start = exon_col_intervals[idx + 1][0]
        if a < next_start and b > prev_end:
            return f"intron_{idx + 1}", "intron", False

    return "intron", "intron", False


def size_class_of(size_bp: int) -> str:
    if size_bp <= 1:
        return "1 bp"
    if size_bp < 50:
        return "2-49 bp"
    return ">=50 bp"


def detect_tandem_repeat_period(
    sv_seq: str,
    left_flank_ref: str = "",
    right_flank_ref: str = "",
    min_identity: float = 0.80,
    max_internal_period: int = 100,
) -> Optional[int]:
    """Check if sv_seq (>=50 bp) is a tandem repeat expansion/contraction and return its period."""
    s = sv_seq.upper()
    n = len(s)
    if n < 2:
        return None

    # 1. Check internal periodicity for period p in 1..min(max_internal_period, n // 2)
    # Exclude p=1 homopolymer if it looks like a mobile-element A/T tail (handled first by mobile element check)
    max_p = min(max_internal_period, n // 2)
    for p in range(1, max_p + 1):
        comp_len = n - p
        if comp_len < 20:
            continue
        matches = sum(1 for i in range(comp_len) if s[i] == s[i + p])
        if matches / comp_len >= min_identity:
            return p

    # 2. Check single-unit or large-unit tandem duplication against immediate left or right reference flank
    if len(left_flank_ref) >= n:
        lf = left_flank_ref[-n:].upper()
        matches = sum(1 for a, b in zip(s, lf) if a == b)
        if matches / n >= min_identity:
            return n
    if len(right_flank_ref) >= n:
        rf = right_flank_ref[:n].upper()
        matches = sum(1 for a, b in zip(s, rf) if a == b)
        if matches / n >= min_identity:
            return n

    return None


def is_heuristic_mobile_element(sv_seq: str) -> bool:
    """Heuristic check for Alu (~300 bp) or L1 (~6000 bp) with poly-A tail (or reverse-strand poly-T head)."""
    s = sv_seq.upper()
    n = len(s)
    in_alu_len = 250 <= n <= 380
    in_l1_len = 4500 <= n <= 7000
    if not (in_alu_len or in_l1_len):
        return False

    tail = s[-35:]
    head = s[:35]
    has_a_tail = ("AAAAAA" in tail) or (sum(1 for c in s[-15:] if c == "A") >= 10)
    has_t_head = ("TTTTTT" in head) or (sum(1 for c in s[:15] if c == "T") >= 10)
    return has_a_tail or has_t_head


def classify_sv_subclass(
    var_type: str,
    ref_seq: str,
    alt_seq: str,
    left_flank_ref: str = "",
    right_flank_ref: str = "",
) -> Tuple[str, Optional[int]]:
    """Sub-classify SVs >= 50 bp into mobile element, tandem-repeat expansion/contraction, other insertion/deletion."""
    l_ref = len(ref_seq)
    l_alt = len(alt_seq)
    if max(l_ref, l_alt) < 50:
        return "none", None

    is_ins = l_alt >= l_ref
    sv_seq = alt_seq if is_ins else ref_seq

    if is_heuristic_mobile_element(sv_seq):
        return "mobile element (heuristic)", None

    period = detect_tandem_repeat_period(sv_seq, left_flank_ref=left_flank_ref, right_flank_ref=right_flank_ref)
    if period is not None:
        if is_ins:
            return "tandem-repeat expansion", period
        return "tandem-repeat contraction", period

    if is_ins:
        return "other insertion", None
    return "other deletion", None


def extract_sequence_variants_from_msa(
    gene: str,
    ref_aligned: str,
    seq_aligned: str,
    exon_col_intervals: Sequence[Tuple[int, int]],
) -> List[VariantEvent]:
    """Extract all variant alleles carried by seq_aligned relative to ref_aligned over seq_aligned's span."""
    st, en = sequence_span_in_msa(seq_aligned)
    if st == -1:
        return []
    r_st, r_en = sequence_span_in_msa(ref_aligned)
    if r_st == -1:
        return []

    eff_start = max(st, r_st)
    eff_end = min(en, r_en)
    if eff_start > eff_end:
        return []

    # Project columns in [eff_start, eff_end] skipping shared gaps ('-' in both)
    proj_cols: List[Tuple[int, str, str, str]] = []
    for c in range(eff_start, eff_end + 1):
        rc = ref_aligned[c].upper()
        sc = seq_aligned[c].upper()
        if rc == "-" and sc == "-":
            continue
        if rc == sc:
            op = "M"
        elif rc != "-" and sc != "-":
            op = "X"
        elif rc == "-":
            op = "I"
        else:
            op = "D"
        proj_cols.append((c, rc, sc, op))

    # Group contiguous non-'M' runs
    raw_runs: List[List[Tuple[int, str, str, str]]] = []
    cur: List[Tuple[int, str, str, str]] = []
    for item in proj_cols:
        if item[3] == "M":
            if cur:
                raw_runs.append(cur)
                cur = []
        else:
            cur.append(item)
    if cur:
        raw_runs.append(cur)

    # Refine runs: if a large pure I or D run (>=50 bp) has tiny flanking X (<=3 bp),
    # split the flanking X so a 300 bp insertion/deletion is not misclassified as complex due to 1 adjacent SNV.
    refined_runs: List[List[Tuple[int, str, str, str]]] = []
    for run in raw_runs:
        ops = [x[3] for x in run]
        n_i = ops.count("I")
        n_d = ops.count("D")
        n_x = ops.count("X")
        if (n_i >= 50 or n_d >= 50) and (0 < n_x <= 3) and (n_i == 0 or n_d == 0):
            # Split contiguous sub-runs by op type ('X' vs 'I'/'D')
            sub: List[Tuple[int, str, str, str]] = []
            for item in run:
                if sub and (sub[-1][3] == "X") != (item[3] == "X"):
                    refined_runs.append(sub)
                    sub = [item]
                else:
                    sub.append(item)
            if sub:
                refined_runs.append(sub)
        else:
            refined_runs.append(run)

    # Precompute reference bases for flank lookup
    events: List[VariantEvent] = []
    for run in refined_runs:
        col_start = run[0][0]
        col_end = run[-1][0] + 1
        ref_seq = "".join(x[1] for x in run if x[1] != "-")
        alt_seq = "".join(x[2] for x in run if x[2] != "-")
        l_ref = len(ref_seq)
        l_alt = len(alt_seq)

        if l_ref == 1 and l_alt == 1:
            var_type = "SNV"
            size_bp = 1
        elif l_ref == l_alt and l_ref > 1:
            var_type = "MNV"
            size_bp = l_ref
        elif l_ref == 0 and l_alt > 0:
            var_type = "insertion"
            size_bp = l_alt
        elif l_ref > 0 and l_alt == 0:
            var_type = "deletion"
            size_bp = l_ref
        else:
            var_type = "complex"
            size_bp = max(l_ref, l_alt)

        sclass = size_class_of(size_bp)
        if sclass == ">=50 bp":
            left_ref = "".join(c for c in ref_aligned[max(0, col_start - 2000) : col_start] if c != "-")
            right_ref = "".join(
                c for c in ref_aligned[col_end : min(len(ref_aligned), col_end + 2000)] if c != "-"
            )
            sv_sub, period = classify_sv_subclass(
                var_type, ref_seq, alt_seq, left_flank_ref=left_ref, right_flank_ref=right_ref
            )
        else:
            sv_sub, period = "none", None

        loc, loc_broad, ag = classify_column_location(gene, col_start, col_end, exon_col_intervals)
        events.append(
            VariantEvent(
                gene=gene,
                col_start=col_start,
                col_end=col_end,
                ref_seq=ref_seq,
                alt_seq=alt_seq,
                var_type=var_type,
                size_bp=size_bp,
                size_class=sclass,
                sv_subclass=sv_sub,
                repeat_period=period,
                location=loc,
                location_broad=loc_broad,
                antigen_binding_exon=ag,
            )
        )
    return events


def parse_minimap2_inversions(paf_text: str, min_block_bp: int = 1) -> List[Dict]:
    """Parse minimap2 PAF output and return any reverse-strand ('-') alignment blocks."""
    inversions: List[Dict] = []
    for line in paf_text.splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 12:
            continue
        qname = fields[0]
        qlen = int(fields[1])
        qstart = int(fields[2])
        qend = int(fields[3])
        strand = fields[4]
        tname = fields[5]
        tlen = int(fields[6])
        tstart = int(fields[7])
        tend = int(fields[8])
        nmatch = int(fields[9])
        blen = int(fields[10])
        mapq = int(fields[11])
        if strand == "-" and (qend - qstart) >= min_block_bp:
            inversions.append(
                {
                    "query": qname,
                    "query_len": qlen,
                    "qstart": qstart,
                    "qend": qend,
                    "strand": strand,
                    "target": tname,
                    "target_len": tlen,
                    "tstart": tstart,
                    "tend": tend,
                    "matches": nmatch,
                    "block_len": blen,
                    "size_bp": qend - qstart,
                    "mapq": mapq,
                }
            )
    return inversions
