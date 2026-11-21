"""Part B: Mosaic allele classification over MSA columns."""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Set, Tuple
import numpy as np

from variant_classifier import (
    VariantEvent,
    extract_sequence_variants_from_msa,
    sequence_span_in_msa,
    size_class_of,
)


@dataclass
class SwitchDetail:
    switch_index: int
    switch_base_0based: int
    switch_msa_col: int
    agreeing_base_start: int
    agreeing_base_end: int
    agreeing_col_start: int
    agreeing_col_end: int
    agreeing_interval_bp: int
    left_source: str
    right_source: str
    divergence_variants_pm150: int
    divergence_cols_pm150: int


@dataclass
class PrivateVariantSummary:
    col_start: int
    col_end: int
    var_type: str
    size_bp: int
    size_class: str
    ref_seq: str
    alt_seq: str


@dataclass
class MosaicClassification:
    query_id: str
    gene: str
    label: str  # "exact", "mosaic", "near-mosaic", "private"
    detailed_label: str  # e.g. "exact", "mosaic-1", "near-mosaic", "private"
    k_switches: Optional[int]
    k_switches_L1: Optional[int]
    exact_full_length: bool
    exact_sources: List[str] = field(default_factory=list)
    segment_sources: List[str] = field(default_factory=list)
    segment_lengths_bp: List[int] = field(default_factory=list)
    segment_col_ranges: List[Tuple[int, int]] = field(default_factory=list)
    switches: List[SwitchDetail] = field(default_factory=list)
    private_variants: List[PrivateVariantSummary] = field(default_factory=list)
    post_mask_k_switches: Optional[int] = None
    post_mask_sublabel: str = ""


_SEQ_ARR_CACHE: Dict[str, Tuple[np.ndarray, Tuple[int, int]]] = {}


def _get_seq_arr_and_span(s: str) -> Tuple[np.ndarray, Tuple[int, int]]:
    res = _SEQ_ARR_CACHE.get(s)
    if res is not None:
        return res
    s_up = s.upper()
    arr = np.frombuffer(s_up.encode("ascii"), dtype=np.uint8)
    non_gap = np.where(arr != ord("-"))[0]
    if len(non_gap) == 0:
        span = (-1, -1)
    else:
        span = (int(non_gap[0]), int(non_gap[-1]))
    res = (arr, span)
    _SEQ_ARR_CACHE[s] = res
    return res


def _precompute_query_panel_matches(
    query_aligned: str,
    panel_names: Sequence[str],
    panel_aligned_seqs: Sequence[str],
    masked_msa_cols: Optional[Set[int]] = None,
) -> Tuple[List[int], np.ndarray, np.ndarray, List[Tuple[int, int]]]:
    """Compute base-level reach matrix for query against panel sequences."""
    q_arr, (q_st, q_en) = _get_seq_arr_and_span(query_aligned)
    if q_st == -1:
        raise ValueError("Query aligned sequence has no non-gap bases")

    base_cols_arr = np.where(q_arr[q_st : q_en + 1] != ord("-"))[0].astype(np.int32) + q_st
    base_cols = base_cols_arr.tolist()
    n = len(base_cols)
    m = len(query_aligned)
    num_p = len(panel_aligned_seqs)

    if num_p == 0:
        return base_cols, np.full((0, n), -1, dtype=np.int32), np.zeros((0, m), dtype=bool), []

    p_items = [_get_seq_arr_and_span(s) for s in panel_aligned_seqs]
    p_mat = np.vstack([it[0] for it in p_items])
    panel_spans = [it[1] for it in p_items]

    col_idx = np.arange(m, dtype=np.int32)[None, :]
    st_arr = np.array([sp[0] for sp in panel_spans], dtype=np.int32)[:, None]
    en_arr = np.array([sp[1] for sp in panel_spans], dtype=np.int32)[:, None]
    valid_span = (st_arr != -1) & (col_idx >= st_arr) & (col_idx <= en_arr)

    eq_full_col = (p_mat == q_arr[None, :]) & valid_span
    if masked_msa_cols:
        mask_idx = np.array(sorted(c for c in masked_msa_cols if 0 <= c < m), dtype=np.int32)
        if len(mask_idx) > 0:
            eq_full_col[:, mask_idx] = valid_span[:, mask_idx]

    base_cols_arr = np.array(base_cols, dtype=np.int32)
    match_base = eq_full_col[:, base_cols_arr]  # shape (P, n)

    # Check intermediate gap columns between base_cols[u] and base_cols[u+1]
    step_ok = np.ones((num_p, max(0, n - 1)), dtype=bool)
    for u in range(n - 1):
        c1 = base_cols[u]
        c2 = base_cols[u + 1]
        if c2 > c1 + 1:
            step_ok[:, u] = np.all(eq_full_col[:, c1 + 1 : c2], axis=1)

    reach_mat = np.full((num_p, n), -1, dtype=np.int32)
    reach_mat[:, n - 1] = np.where(match_base[:, n - 1], n - 1, -1)
    for u in range(n - 2, -1, -1):
        can_extend = match_base[:, u] & step_ok[:, u] & (reach_mat[:, u + 1] != -1)
        reach_mat[:, u] = np.where(
            can_extend,
            reach_mat[:, u + 1],
            np.where(match_base[:, u], u, -1),
        )

    return base_cols, reach_mat, eq_full_col, panel_spans


def greedy_longest_match_L1(max_reach: Sequence[int]) -> Optional[Tuple[int, List[int]]]:
    """Greedy longest-match solver for L=1 minimum segment count.

    Returns (k_switches, switch_bases) where switch_bases has length k_switches,
    or None if no valid mosaic cover exists.
    """
    n = len(max_reach)
    if n == 0:
        return 0, []
    if max_reach[0] < 0:
        return None
    if max_reach[0] >= n - 1:
        return 0, []

    switches: List[int] = []
    cur_start = 0
    cur_reach = max_reach[0]

    while cur_reach < n - 1:
        best_next_reach = -1
        best_sw = -1
        # Any base s in [cur_start + 1, cur_reach] is an agreeing column with the current source
        for s in range(cur_start + 1, cur_reach + 1):
            r = max_reach[s]
            if r > best_next_reach:
                best_next_reach = r
                best_sw = s
        if best_next_reach <= cur_reach or best_sw == -1:
            return None
        switches.append(best_sw)
        cur_start = best_sw
        cur_reach = best_next_reach

    return len(switches), switches


def brute_force_min_switches(max_reach: Sequence[int], min_seg_bp: int = 1) -> Optional[int]:
    """Exhaustive BFS/brute-force solver for toy testing greedy optimality."""
    n = len(max_reach)
    if n < min_seg_bp or max_reach[0] < 0:
        return None
    if max_reach[0] >= n - 1 and n >= min_seg_bp:
        return 0

    # BFS over current segment start base u in 0..n-1
    dist = [-1] * n
    dist[0] = 0
    q = deque([0])
    while q:
        u = q.popleft()
        r = max_reach[u]
        if r == n - 1 and (n - u) >= min_seg_bp:
            return dist[u]
        upper = min(r, n - min_seg_bp)
        for v in range(u + min_seg_bp, upper + 1):
            if dist[v] == -1:
                dist[v] = dist[u] + 1
                q.append(v)
    return None


def solve_min_switches_with_min_length(
    reach_mat: np.ndarray,
    panel_names: Sequence[str],
    min_seg_bp: int,
) -> Optional[Tuple[int, List[int], List[int], List[Tuple[int, int]]]]:
    """Find minimum switches k with each segment >= min_seg_bp and agreeing switch columns.

    When multiple switch placements achieve the same minimum k, maximizes the minimum
    segment length (and balances switch positions near the middle of each agreeing interval).

    Returns:
      (k_switches, switch_bases, source_indices, agreeing_base_intervals)
      or None if impossible.
    """
    num_p, n = reach_mat.shape
    if num_p == 0 or n < min_seg_bp:
        return None

    max_reach = np.max(reach_mat, axis=0)
    if max_reach[0] < 0:
        return None

    # Check k = 0 (exact match over query gene interval)
    if max_reach[0] == n - 1:
        p_idx = int(np.argmax(reach_mat[:, 0] == n - 1))
        return 0, [], [p_idx], []

    # Step 1: Find minimum number of switches min_k via vectorized BFS over states 0..n-1
    dist = np.full(n, -1, dtype=np.int32)
    dist[0] = 0
    min_k: Optional[int] = None

    q = deque([0])
    while q:
        u = q.popleft()
        d = int(dist[u])
        if min_k is not None and d >= min_k:
            continue
        r = int(max_reach[u])
        if r == n - 1 and (n - u) >= min_seg_bp:
            min_k = d
            break
        lo = u + min_seg_bp
        upper = min(r, n - min_seg_bp)
        if upper >= lo:
            unvis = np.where(dist[lo : upper + 1] == -1)[0]
            if len(unvis) > 0:
                v_new = lo + unvis
                dist[v_new] = d + 1
                q.extend(v_new.tolist())

    if min_k is None:
        return None

    # Step 2: Among paths of length min_k, maximize minimum segment length via vectorized DP over layers 0..min_k
    layers: List[List[int]] = [
        np.where(dist == step)[0].tolist() for step in range(min_k + 1)
    ]

    # Filter backward to only states that can reach n in (min_k - layer) steps
    can_finish = np.zeros(n, dtype=bool)
    for u in layers[min_k]:
        if max_reach[u] == n - 1 and (n - u) >= min_seg_bp:
            can_finish[u] = True

    for step in range(min_k - 1, -1, -1):
        next_valid = (dist == step + 1) & can_finish
        for u in layers[step]:
            r = int(max_reach[u])
            lo = u + min_seg_bp
            upper = min(r, n - min_seg_bp)
            if upper >= lo and np.any(next_valid[lo : upper + 1]):
                can_finish[u] = True

    # Forward vectorized DP to maximize min segment length along valid states
    best_Bottleneck = np.full(n, -1, dtype=np.int32)
    parent_state = np.full(n, -1, dtype=np.int32)
    best_Bottleneck[0] = n + 1

    for step in range(min_k):
        next_valid = (dist == step + 1) & can_finish
        for u in layers[step]:
            if not can_finish[u] or best_Bottleneck[u] < 0:
                continue
            r = int(max_reach[u])
            lo = u + min_seg_bp
            upper = min(r, n - min_seg_bp)
            if upper < lo:
                continue
            v_cands = np.where(next_valid[lo : upper + 1])[0] + lo
            if len(v_cands) == 0:
                continue
            seg_lens = v_cands - u
            cands = np.minimum(int(best_Bottleneck[u]), seg_lens)
            improve = cands > best_Bottleneck[v_cands]
            if np.any(improve):
                imp_v = v_cands[improve]
                best_Bottleneck[imp_v] = cands[improve]
                parent_state[imp_v] = u

    best_last_u = -1
    best_total_bottleneck = -1
    for u in layers[min_k]:
        if can_finish[u] and best_Bottleneck[u] >= 0:
            last_len = n - u
            cand = min(int(best_Bottleneck[u]), last_len)
            if cand > best_total_bottleneck:
                best_total_bottleneck = cand
                best_last_u = u

    if best_last_u == -1:
        return None

    # Reconstruct switch bases s_1..s_k
    switch_bases: List[int] = []
    cur = best_last_u
    while cur != 0:
        switch_bases.append(cur)
        cur = int(parent_state[cur])
    switch_bases.reverse()

    # Select source for each segment [s_i, s_{i+1}]
    bounds = [0] + switch_bases + [n]
    source_indices: List[int] = []
    for i in range(len(bounds) - 1):
        u = bounds[i]
        req_end = bounds[i + 1] if i + 1 < len(bounds) - 1 else n - 1
        # Choose panel member covering [u, req_end] with longest total match span around [u, req_end]
        candidates = np.where(reach_mat[:, u] >= req_end)[0]
        best_p = int(candidates[0])
        best_span = int(reach_mat[best_p, u]) - u
        for p in candidates[1:]:
            sp = int(reach_mat[int(p), u]) - u
            if sp > best_span:
                best_span = sp
                best_p = int(p)
        source_indices.append(best_p)

    # Compute agreeing base intervals [agree_start, agree_end] for each switch
    agreeing_intervals: List[Tuple[int, int]] = []
    for i in range(len(switch_bases)):
        sw = switch_bases[i]
        p_left = source_indices[i]
        p_right = source_indices[i + 1]
        # Left source covers from bounds[i] up to reach_mat[p_left, bounds[i]] >= sw
        b_end = int(reach_mat[p_left, bounds[i]])
        # Find earliest base a_start <= sw from which p_right reaches at least sw
        a_start = sw
        while a_start > 0 and reach_mat[p_right, a_start - 1] >= sw:
            a_start -= 1
        agreeing_intervals.append((a_start, b_end))

    return min_k, switch_bases, source_indices, agreeing_intervals


def count_divergence_near_switch(
    left_aligned: str,
    right_aligned: str,
    base_cols: Sequence[int],
    switch_base_idx: int,
    window_bp: int = 150,
) -> Tuple[int, int]:
    """Count variant events and differing MSA columns between adjacent sources within +-150 bp of switch."""
    n = len(base_cols)
    b_start = max(0, switch_base_idx - window_bp)
    b_end = min(n - 1, switch_base_idx + window_bp)
    c_start = base_cols[b_start]
    c_end = base_cols[b_end]

    diff_cols = 0
    variant_events = 0
    in_diff = False

    for c in range(c_start, c_end + 1):
        lc = left_aligned[c].upper()
        rc = right_aligned[c].upper()
        if lc == "-" and rc == "-":
            continue
        if lc != rc:
            diff_cols += 1
            if not in_diff:
                variant_events += 1
                in_diff = True
        else:
            in_diff = False

    return variant_events, diff_cols


def find_private_variants_of_query(
    gene: str,
    query_aligned: str,
    panel_aligned_seqs: Sequence[str],
    ref_backbone: str,
    exon_col_intervals: Sequence[Tuple[int, int]],
    precomputed_seq_events: Optional[Dict[str, Set[Tuple[int, int, str, str]]]] = None,
) -> Tuple[List[PrivateVariantSummary], Set[int]]:
    """Identify variant alleles of query absent from all of P(q), and return their MSA columns."""
    q_arr, (q_st, q_en) = _get_seq_arr_and_span(query_aligned)
    if q_st == -1:
        return [], set()

    m = len(query_aligned)
    p_items = [_get_seq_arr_and_span(s) for s in panel_aligned_seqs]
    panel_spans = [it[1] for it in p_items]

    # 1. Column-level private columns: columns c in [q_st, q_en] where no panel sequence covering c has q[c]
    private_cols: Set[int] = set()
    if len(panel_aligned_seqs) == 0:
        private_cols.update(range(q_st, q_en + 1))
    else:
        p_mat = np.vstack([it[0] for it in p_items])
        col_idx = np.arange(q_st, q_en + 1, dtype=np.int32)[None, :]
        st_arr = np.array([sp[0] for sp in panel_spans], dtype=np.int32)[:, None]
        en_arr = np.array([sp[1] for sp in panel_spans], dtype=np.int32)[:, None]
        valid_span = (st_arr != -1) & (col_idx >= st_arr) & (col_idx <= en_arr)
        any_match = np.any((p_mat[:, q_st : q_en + 1] == q_arr[None, q_st : q_en + 1]) & valid_span, axis=0)
        for offset in np.where(~any_match)[0]:
            private_cols.add(q_st + int(offset))

    # 2. Event-level private variants relative to ref_backbone (captures private deletions/insertions/MNVs
    # whose exact allele is absent from P(q) even if individual columns appear in different indels)
    q_events = extract_sequence_variants_from_msa(gene, ref_backbone, query_aligned, exon_col_intervals)
    panel_event_keys: Set[Tuple[int, int, str, str]] = set()
    for s in set(panel_aligned_seqs):
        if precomputed_seq_events is not None and s in precomputed_seq_events:
            panel_event_keys.update(precomputed_seq_events[s])
        else:
            evs = {
                (ev.col_start, ev.col_end, ev.ref_seq, ev.alt_seq)
                for ev in extract_sequence_variants_from_msa(gene, ref_backbone, s, exon_col_intervals)
            }
            if precomputed_seq_events is not None:
                precomputed_seq_events[s] = evs
            panel_event_keys.update(evs)

    summaries: List[PrivateVariantSummary] = []
    covered_by_event: Set[int] = set()

    for ev in q_events:
        key = (ev.col_start, ev.col_end, ev.ref_seq, ev.alt_seq)
        ev_cols = set(range(ev.col_start, ev.col_end))
        if key not in panel_event_keys or (ev_cols & private_cols):
            summaries.append(
                PrivateVariantSummary(
                    col_start=ev.col_start,
                    col_end=ev.col_end,
                    var_type=ev.var_type,
                    size_bp=ev.size_bp,
                    size_class=ev.size_class,
                    ref_seq=ev.ref_seq,
                    alt_seq=ev.alt_seq,
                )
            )
            covered_by_event.update(ev_cols)
            private_cols.update(ev_cols)

    # Any remaining private_cols not in q_events (e.g. where query matched ref_backbone because query was in ref)
    rem_cols = sorted(private_cols - covered_by_event)
    if rem_cols:
        runs: List[List[int]] = []
        cur: List[int] = []
        for c in rem_cols:
            if cur and c > cur[-1] + 1:
                runs.append(cur)
                cur = [c]
            else:
                cur.append(c)
        if cur:
            runs.append(cur)

        for run in runs:
            cs, ce = run[0], run[-1] + 1
            q_sub = "".join(query_aligned[c].upper() for c in run if query_aligned[c] != "-")
            # Pick a representative panel sequence covering [cs, ce] to define ref_sub
            p_sub = ""
            for s, (st, en) in zip(panel_aligned_seqs, panel_spans):
                if st != -1 and st <= cs and en >= ce - 1:
                    p_sub = "".join(s[c].upper() for c in run if s[c] != "-")
                    break
            lr, la = len(p_sub), len(q_sub)
            if lr == 1 and la == 1:
                vt, sz = "SNV", 1
            elif lr == la and lr > 1:
                vt, sz = "MNV", lr
            elif lr == 0 and la > 0:
                vt, sz = "insertion", la
            elif lr > 0 and la == 0:
                vt, sz = "deletion", lr
            else:
                vt, sz = "complex", max(lr, la, len(run))
            summaries.append(
                PrivateVariantSummary(
                    col_start=cs,
                    col_end=ce,
                    var_type=vt,
                    size_bp=sz,
                    size_class=size_class_of(sz),
                    ref_seq=p_sub,
                    alt_seq=q_sub,
                )
            )

    summaries.sort(key=lambda x: (x.col_start, x.col_end))
    return summaries, private_cols


def classify_query_against_panel(
    gene: str,
    query_id: str,
    query_aligned: str,
    panel_aligned_dict: Dict[str, str],
    ref_backbone: Optional[str] = None,
    exon_col_intervals: Sequence[Tuple[int, int]] = (),
    min_seg_bp: int = 150,
    precomputed_seq_events: Optional[Dict[str, Set[Tuple[int, int, str, str]]]] = None,
) -> MosaicClassification:
    """Classify query sequence q against panel P(q) into exact, mosaic-k, near-mosaic, or private."""
    panel_names = list(panel_aligned_dict.keys())
    panel_seqs = [panel_aligned_dict[k] for k in panel_names]

    if ref_backbone is None:
        ref_backbone = query_aligned

    q_ungapped = "".join(c for c in query_aligned.upper() if c != "-")
    exact_full_sources = [
        name
        for name, s in zip(panel_names, panel_seqs)
        if "".join(c for c in s.upper() if c != "-") == q_ungapped
    ]
    exact_full_length = len(exact_full_sources) > 0

    base_cols, reach_mat, _, _ = _precompute_query_panel_matches(
        query_aligned, panel_names, panel_seqs
    )
    n = len(base_cols)

    if len(panel_seqs) == 0:
        priv_vars, _ = find_private_variants_of_query(
            gene, query_aligned, panel_seqs, ref_backbone, exon_col_intervals, precomputed_seq_events
        )
        return MosaicClassification(
            query_id=query_id,
            gene=gene,
            label="private",
            detailed_label="private",
            k_switches=None,
            k_switches_L1=None,
            exact_full_length=False,
            private_variants=priv_vars,
            post_mask_k_switches=None,
            post_mask_sublabel="empty_panel",
        )

    max_reach = np.max(reach_mat, axis=0)
    res_L1 = greedy_longest_match_L1(max_reach.tolist())
    k_L1 = res_L1[0] if res_L1 is not None else None

    res_L = solve_min_switches_with_min_length(reach_mat, panel_names, min_seg_bp=min_seg_bp)

    if res_L is not None:
        k_sw, sw_bases, src_indices, agree_ivs = res_L
        if k_sw == 0:
            # Find all panel members that match across the entire query gene interval [0, n-1]
            matching_idx = np.where(reach_mat[:, 0] == n - 1)[0]
            ex_sources = [panel_names[int(i)] for i in matching_idx]
            return MosaicClassification(
                query_id=query_id,
                gene=gene,
                label="exact",
                detailed_label="exact",
                k_switches=0,
                k_switches_L1=0,
                exact_full_length=exact_full_length,
                exact_sources=ex_sources,
                segment_sources=[ex_sources[0]],
                segment_lengths_bp=[n],
                segment_col_ranges=[(base_cols[0], base_cols[-1] + 1)],
            )

        # k_sw >= 1 -> mosaic-k
        bounds = [0] + sw_bases + [n]
        seg_sources = [panel_names[idx] for idx in src_indices]
        seg_lengths = [bounds[i + 1] - bounds[i] for i in range(len(bounds) - 1)]
        seg_col_ranges = [
            (base_cols[bounds[i]], base_cols[bounds[i + 1] - 1] + 1)
            for i in range(len(bounds) - 1)
        ]

        switches: List[SwitchDetail] = []
        for i, sw_b in enumerate(sw_bases):
            a_st, b_en = agree_ivs[i]
            p_left_seq = panel_seqs[src_indices[i]]
            p_right_seq = panel_seqs[src_indices[i + 1]]
            div_vars, div_cols = count_divergence_near_switch(
                p_left_seq, p_right_seq, base_cols, sw_b, window_bp=150
            )
            switches.append(
                SwitchDetail(
                    switch_index=i + 1,
                    switch_base_0based=sw_b,
                    switch_msa_col=base_cols[sw_b],
                    agreeing_base_start=a_st,
                    agreeing_base_end=b_en,
                    agreeing_col_start=base_cols[a_st],
                    agreeing_col_end=base_cols[b_en],
                    agreeing_interval_bp=b_en - a_st + 1,
                    left_source=seg_sources[i],
                    right_source=seg_sources[i + 1],
                    divergence_variants_pm150=div_vars,
                    divergence_cols_pm150=div_cols,
                )
            )

        return MosaicClassification(
            query_id=query_id,
            gene=gene,
            label="mosaic",
            detailed_label=f"mosaic-{k_sw}",
            k_switches=k_sw,
            k_switches_L1=k_L1,
            exact_full_length=False,
            segment_sources=seg_sources,
            segment_lengths_bp=seg_lengths,
            segment_col_ranges=seg_col_ranges,
            switches=switches,
        )

    # Not exact and not mosaic-k at min_seg_bp: find private variants and test post-masking
    priv_vars, priv_cols = find_private_variants_of_query(
        gene, query_aligned, panel_seqs, ref_backbone, exon_col_intervals, precomputed_seq_events
    )

    if priv_cols:
        _, masked_reach, _, _ = _precompute_query_panel_matches(
            query_aligned, panel_names, panel_seqs, masked_msa_cols=priv_cols
        )
        res_masked = solve_min_switches_with_min_length(
            masked_reach, panel_names, min_seg_bp=min_seg_bp
        )
        if res_masked is not None:
            mk_sw, m_sw_bases, m_src_idx, m_agree = res_masked
            if mk_sw >= 1:
                bounds = [0] + m_sw_bases + [n]
                seg_sources = [panel_names[idx] for idx in m_src_idx]
                seg_lengths = [bounds[i + 1] - bounds[i] for i in range(len(bounds) - 1)]
                seg_col_ranges = [
                    (base_cols[bounds[i]], base_cols[bounds[i + 1] - 1] + 1)
                    for i in range(len(bounds) - 1)
                ]
                return MosaicClassification(
                    query_id=query_id,
                    gene=gene,
                    label="near-mosaic",
                    detailed_label="near-mosaic",
                    k_switches=None,
                    k_switches_L1=k_L1,
                    exact_full_length=False,
                    segment_sources=seg_sources,
                    segment_lengths_bp=seg_lengths,
                    segment_col_ranges=seg_col_ranges,
                    private_variants=priv_vars,
                    post_mask_k_switches=mk_sw,
                    post_mask_sublabel=f"mosaic-{mk_sw}_after_masking",
                )
            else:
                # mk_sw == 0: becomes exact (single panel allele) after removing private variants
                matching_idx = np.where(masked_reach[:, 0] == n - 1)[0]
                ex_sources = [panel_names[int(i)] for i in matching_idx]
                return MosaicClassification(
                    query_id=query_id,
                    gene=gene,
                    label="private",
                    detailed_label="private",
                    k_switches=None,
                    k_switches_L1=k_L1,
                    exact_full_length=False,
                    segment_sources=[ex_sources[0]] if ex_sources else [],
                    segment_lengths_bp=[n],
                    private_variants=priv_vars,
                    post_mask_k_switches=0,
                    post_mask_sublabel="exact_plus_private",
                )

    sublabel = "mosaic_short_segments_only" if (k_L1 is not None and k_L1 >= 1) else "unbridgeable_or_private"
    return MosaicClassification(
        query_id=query_id,
        gene=gene,
        label="private",
        detailed_label="private",
        k_switches=None,
        k_switches_L1=k_L1,
        exact_full_length=False,
        private_variants=priv_vars,
        post_mask_k_switches=None,
        post_mask_sublabel=sublabel,
    )
