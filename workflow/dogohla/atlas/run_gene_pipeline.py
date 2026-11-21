"""Per-gene pipeline for Phase 1:
- Step 1: Build per-gene MSA (cohort core sequences + IPD alleles) via MAFFT --auto
- Step 2 (Part A): Call variant atlas from MSA + check inversions via minimap2
- Step 3 (Part B): Classify every cohort truth copy (full panel, IPD-only, cohort-only, Arab+Saudi, 5-fold CV)
- Step 4 (Part C): Generate 200 positive + 200 negative synthetic mosaics from HPRC pairs and classify them
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import os
import random
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, FrozenSet, List, Optional, Sequence, Set, Tuple
import numpy as np

from panels import (
    GENES,
    build_exclusion_map_from_folds,
    load_augmented_fold_grouping,
    load_ipd_reference,
    load_source_records,
    select_cohort_panel_for_donor,
    sha256_str,
    verify_no_leakage,
)
from variant_classifier import (
    VariantEvent,
    build_reference_backbone,
    extract_sequence_variants_from_msa,
    map_exons_to_msa_columns,
    parse_minimap2_inversions,
    sequence_span_in_msa,
)
from mosaic_classifier import (
    MosaicClassification,
    classify_query_against_panel,
)


STAGE_ROOT = Path("/encrypted/e3008/hohndor-dogohla-20261020/mosaic-variant-atlas-20261103")
AUGMENTED_ROOT = Path("/encrypted/e3008/hohndor-dogohla-20261020/simulation-augmented-regional-20261024")
GENOMIC_FA = Path(
    "/ibex/scratch/projects/c2014/rob/dogohla-benchmark/hprc-zenodo-20261004/baselines/assets/references/genomic.fa"
)


def read_fasta_dict(fasta_path: Path) -> Dict[str, str]:
    records: Dict[str, str] = {}
    name: Optional[str] = None
    chunks: List[str] = []
    with Path(fasta_path).open() as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(chunks).upper()
                name = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            records[name] = "".join(chunks).upper()
    return records


def build_or_load_gene_msa(
    gene: str,
    cohort_records: List[Dict],
    ipd_records: List[Dict],
    msa_dir: Path,
    mafft_bin: Path,
    threads: int = 16,
    force: bool = False,
) -> Tuple[Dict[str, str], Dict[str, str], Dict[str, str], Dict]:
    """Build deduplicated MAFFT MSA for gene and return (cohort_aligned, ipd_aligned, uniq_aligned, stats)."""
    msa_dir.mkdir(parents=True, exist_ok=True)
    uniq_in_fa = msa_dir / f"{gene}_unique_in.fa"
    uniq_out_fa = msa_dir / f"{gene}_unique_out.fa"
    stats_json = msa_dir / f"{gene}_msa_stats.json"

    # Choose canonical IPD allele to lead the FASTA
    canon_allele = f"HLA-{gene}*01:01:01:01"
    ipd_by_name = {x["allele"]: x for x in ipd_records}
    if canon_allele not in ipd_by_name:
        canon_allele = ipd_records[0]["allele"]

    # Map unique ungapped sequence -> unique ID u0, u1, ...
    seq_to_uid: Dict[str, str] = {}
    uid_to_seq: Dict[str, str] = {}

    def _register_seq(s: str) -> str:
        s_up = s.upper()
        if s_up not in seq_to_uid:
            uid = f"U{len(seq_to_uid):05d}"
            seq_to_uid[s_up] = uid
            uid_to_seq[uid] = s_up
        return seq_to_uid[s_up]

    # Register canonical IPD allele first as U00000
    _register_seq(ipd_by_name[canon_allele]["sequence"])
    for x in ipd_records:
        _register_seq(x["sequence"])
    for r in cohort_records:
        _register_seq(r["core_sequence"])

    t0 = time.time()
    if not uniq_out_fa.exists() or uniq_out_fa.stat().st_size == 0 or force:
        with uniq_in_fa.open("w") as f:
            for uid, s in uid_to_seq.items():
                f.write(f">{uid}\n{s}\n")
        cmd = [str(mafft_bin), "--auto", "--thread", str(threads), str(uniq_in_fa)]
        print(f"[{gene}] Running MAFFT: {' '.join(cmd)} on {len(uid_to_seq)} unique sequences...", flush=True)
        with uniq_out_fa.open("w") as out_f:
            subprocess.run(cmd, stdout=out_f, check=True)
        runtime_sec = round(time.time() - t0, 2)
    else:
        print(f"[{gene}] Reusing existing MSA at {uniq_out_fa}", flush=True)
        runtime_sec = 0.0

    uniq_aligned = read_fasta_dict(uniq_out_fa)
    if len(uniq_aligned) != len(uid_to_seq):
        raise RuntimeError(
            f"MSA sequence count mismatch for {gene}: expected {len(uid_to_seq)}, got {len(uniq_aligned)}"
        )

    cohort_aligned: Dict[str, str] = {}
    for r in cohort_records:
        uid = seq_to_uid[r["core_sequence"].upper()]
        cohort_aligned[r["source_id"]] = uniq_aligned[uid]

    ipd_aligned: Dict[str, str] = {}
    for x in ipd_records:
        uid = seq_to_uid[x["sequence"].upper()]
        ipd_aligned[x["allele"]] = uniq_aligned[uid]

    msa_len = len(next(iter(uniq_aligned.values())))
    stats = {
        "gene": gene,
        "n_cohort_copies": len(cohort_records),
        "n_ipd_alleles": len(ipd_records),
        "n_unique_sequences": len(uid_to_seq),
        "msa_columns": msa_len,
        "canonical_allele": canon_allele,
        "mafft_runtime_sec": runtime_sec,
    }
    stats_json.write_text(json.dumps(stats, indent=2))
    return cohort_aligned, ipd_aligned, uniq_aligned, stats


def run_part_a_for_gene(
    gene: str,
    cohort_records: List[Dict],
    ipd_records: List[Dict],
    cohort_aligned: Dict[str, str],
    ipd_aligned: Dict[str, str],
    uniq_aligned: Dict[str, str],
    canonical_allele: str,
    minimap2_bin: Path,
    per_gene_dir: Path,
    tmp_dir: Path,
) -> Tuple[ str, List[Tuple[int, int]], Dict[str, Set[Tuple[int, int, str, str]]] ]:
    """Run Part A variant atlas and minimap2 inversion check for a single gene."""
    print(f"[{gene}] Starting Part A (Variant Atlas & Inversion check)...", flush=True)
    ipd_by_name = {x["allele"]: x for x in ipd_records}
    canon_entry = ipd_by_name[canonical_allele]
    canon_aligned = ipd_aligned[canonical_allele]

    ref_backbone = build_reference_backbone(
        list(uniq_aligned.values()),
        preferred_ref_aligned=canon_aligned,
    )
    exon_col_intervals = map_exons_to_msa_columns(canon_aligned, canon_entry["exons"])

    # Extract variants once per unique aligned sequence
    events_by_aligned_seq: Dict[str, List[VariantEvent]] = {}
    precomputed_seq_events: Dict[str, Set[Tuple[int, int, str, str]]] = {}
    for aln_seq in set(uniq_aligned.values()):
        evs = extract_sequence_variants_from_msa(gene, ref_backbone, aln_seq, exon_col_intervals)
        events_by_aligned_seq[aln_seq] = evs
        precomputed_seq_events[aln_seq] = {
            (e.col_start, e.col_end, e.ref_seq, e.alt_seq) for e in evs
        }

    # Aggregate frequencies per variant allele (col_start, col_end, ref_seq, alt_seq)
    allele_meta: Dict[Tuple[int, int, str, str], VariantEvent] = {}
    counts_by_cohort: Dict[Tuple[int, int, str, str], Counter[str]] = defaultdict(Counter)
    ipd_counts: Counter[Tuple[int, int, str, str]] = Counter()

    for r in cohort_records:
        aln = cohort_aligned[r["source_id"]]
        res_name = r["resource"]  # "HPRC", "Arab", "Saudi"
        for ev in events_by_aligned_seq[aln]:
            key = (ev.col_start, ev.col_end, ev.ref_seq, ev.alt_seq)
            allele_meta[key] = ev
            counts_by_cohort[key][res_name] += 1

    for x in ipd_records:
        aln = ipd_aligned[x["allele"]]
        for ev in events_by_aligned_seq[aln]:
            key = (ev.col_start, ev.col_end, ev.ref_seq, ev.alt_seq)
            allele_meta[key] = ev
            ipd_counts[key] += 1

    # Inversion check via minimap2: align each cohort sequence to its closest IPD allele
    inv_dir = tmp_dir / "inversions"
    inv_dir.mkdir(parents=True, exist_ok=True)

    # Identify closest IPD allele for each unique cohort aligned sequence using MSA column agreement
    ipd_names_list = [x["allele"] for x in ipd_records]
    ipd_mat = np.vstack(
        [np.frombuffer(ipd_aligned[name].encode("ascii"), dtype=np.uint8) for name in ipd_names_list]
    )
    closest_ipd_for_aln: Dict[str, str] = {}
    for r in cohort_records:
        c_aln = cohort_aligned[r["source_id"]]
        if c_aln not in closest_ipd_for_aln:
            c_arr = np.frombuffer(c_aln.encode("ascii"), dtype=np.uint8)
            matches = np.sum((ipd_mat == c_arr[None, :]) & (c_arr[None, :] != ord("-")), axis=1)
            best_idx = int(np.argmax(matches))
            closest_ipd_for_aln[c_aln] = ipd_names_list[best_idx]

    used_ipd_targets = sorted(set(closest_ipd_for_aln.values()) | {canonical_allele})
    ref_fa = inv_dir / f"{gene}_closest_ipd_refs.fa"
    qry_fa = inv_dir / f"{gene}_cohort_queries.fa"
    with ref_fa.open("w") as f:
        for t_name in used_ipd_targets:
            f.write(f">{t_name}\n{ipd_by_name[t_name]['sequence']}\n")
    with qry_fa.open("w") as f:
        for r in cohort_records:
            f.write(f">{r['source_id']}\n{r['core_sequence']}\n")

    query_to_closest_ipd = {
        r["source_id"]: closest_ipd_for_aln[cohort_aligned[r["source_id"]]] for r in cohort_records
    }
    cmd_mm2 = [str(minimap2_bin), "-x", "asm20", "-c", "-N", "50", "-p", "0.1", str(ref_fa), str(qry_fa)]
    paf_proc = subprocess.run(cmd_mm2, capture_output=True, text=True, check=True)
    raw_invs = parse_minimap2_inversions(paf_proc.stdout, min_block_bp=50)
    inversions = [
        inv
        for inv in raw_invs
        if inv["target"] in (query_to_closest_ipd.get(inv["query"]), canonical_allele)
    ]

    inv_tsv = per_gene_dir / f"{gene}_inversions.tsv"
    with inv_tsv.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "gene", "query", "query_len", "qstart", "qend", "strand",
                "target", "target_len", "tstart", "tend", "size_bp", "mapq",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for inv in inversions:
            row = {"gene": gene, **inv}
            writer.writerow({k: row[k] for k in writer.fieldnames})

    # Write per-gene variant atlas TSV (no raw sequences > 1 bp exported so Hard Rule 7 is strictly respected)
    var_tsv = per_gene_dir / f"{gene}_variants.tsv"
    fieldnames = [
        "gene",
        "site_id",
        "allele_id",
        "col_start",
        "col_end",
        "var_type",
        "ref_len",
        "alt_len",
        "size_bp",
        "size_class",
        "sv_subclass",
        "repeat_period",
        "location",
        "location_broad",
        "antigen_binding_exon",
        "hprc_count",
        "arab_count",
        "saudi_count",
        "cohort_total_count",
        "ipd_count",
        "in_cohort",
        "absent_from_ipd",
    ]
    sorted_keys = sorted(allele_meta.keys(), key=lambda k: (k[0], k[1], len(k[2]), len(k[3]), k[2], k[3]))
    with var_tsv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for key in sorted_keys:
            ev = allele_meta[key]
            c_cnt = counts_by_cohort[key]
            h_c = c_cnt.get("HPRC", 0)
            a_c = c_cnt.get("Arab", 0)
            s_c = c_cnt.get("Saudi", 0)
            tot_c = h_c + a_c + s_c
            i_c = ipd_counts.get(key, 0)
            aid = hashlib.sha256(f"{gene}:{ev.col_start}:{ev.col_end}:{ev.ref_seq}:{ev.alt_seq}".encode()).hexdigest()[:12]
            writer.writerow(
                {
                    "gene": gene,
                    "site_id": f"{gene}:{ev.col_start}-{ev.col_end}",
                    "allele_id": f"{gene}_{aid}",
                    "col_start": ev.col_start,
                    "col_end": ev.col_end,
                    "var_type": ev.var_type,
                    "ref_len": len(ev.ref_seq),
                    "alt_len": len(ev.alt_seq),
                    "size_bp": ev.size_bp,
                    "size_class": ev.size_class,
                    "sv_subclass": ev.sv_subclass,
                    "repeat_period": "" if ev.repeat_period is None else ev.repeat_period,
                    "location": ev.location,
                    "location_broad": ev.location_broad,
                    "antigen_binding_exon": ev.antigen_binding_exon,
                    "hprc_count": h_c,
                    "arab_count": a_c,
                    "saudi_count": s_c,
                    "cohort_total_count": tot_c,
                    "ipd_count": i_c,
                    "in_cohort": tot_c > 0,
                    "absent_from_ipd": i_c == 0,
                }
            )

    print(
        f"[{gene}] Finished Part A: {len(sorted_keys)} total variant alleles "
        f"({sum(1 for k in sorted_keys if sum(counts_by_cohort[k].values()) > 0)} in cohort), "
        f"{len(inversions)} minimap2 inversions.",
        flush=True,
    )
    return ref_backbone, exon_col_intervals, precomputed_seq_events


def _dedup_panel_dict(panel_dict: Dict[str, str]) -> Dict[str, str]:
    """Deduplicate identical aligned sequences in panel_dict while preserving deterministic representative names."""
    seen_seqs: Dict[str, str] = {}
    deduped: Dict[str, str] = {}
    for name, seq in panel_dict.items():
        if seq not in seen_seqs:
            seen_seqs[seq] = name
            deduped[name] = seq
    return deduped


def run_part_b_for_gene(
    gene: str,
    cohort_records: List[Dict],
    ipd_records: List[Dict],
    cohort_aligned: Dict[str, str],
    ipd_aligned: Dict[str, str],
    excl_map: Dict[str, Set[str]],
    cv5_excl_map: Dict[str, Set[str]],
    ref_backbone: str,
    exon_col_intervals: List[Tuple[int, int]],
    precomputed_seq_events: Dict[str, Set[Tuple[int, int, str, str]]],
    per_gene_dir: Path,
) -> List[Dict]:
    """Run Part B mosaic classification for every cohort truth copy of gene."""
    print(f"[{gene}] Starting Part B (Mosaic classification on {len(cohort_records)} copies)...", flush=True)
    ipd_dedup = _dedup_panel_dict(ipd_aligned)
    ipd_dedup_seqs = set(ipd_dedup.values())
    ipd_full_hashes = {x["sequence_sha256"] for x in ipd_records}
    part_b_cache: Dict[Tuple[str, FrozenSet[str], str], MosaicClassification] = {}

    def _cached_classify(
        panel_tag: str,
        excl_set: FrozenSet[str],
        qid: str,
        q_aln: str,
        panel_dict: Dict[str, str],
    ) -> MosaicClassification:
        ckey = (panel_tag, excl_set, q_aln)
        if ckey in part_b_cache:
            return part_b_cache[ckey]
        res = classify_query_against_panel(
            gene, qid, q_aln, panel_dict, ref_backbone, exon_col_intervals, 150, precomputed_seq_events
        )
        part_b_cache[ckey] = res
        return res

    def _l1_labels(res: MosaicClassification) -> Tuple[str, str]:
        if res.k_switches_L1 == 0:
            return "exact", "exact"
        if res.k_switches_L1 is not None and res.k_switches_L1 >= 1:
            return "mosaic", f"mosaic-{res.k_switches_L1}"
        if res.label == "near-mosaic":
            return "near-mosaic", "near-mosaic"
        return "private", "private"

    rows: List[Dict] = []
    for idx, r in enumerate(cohort_records, start=1):
        qid = r["source_id"]
        q_donor = r["donor"]
        q_aln = cohort_aligned[qid]
        q_hash = r["core_sha256"]
        ipd_absent = len(r.get("genomic_labels", [])) == 0

        # 1. Leave-donor-and-family-out cohort panel from FOLDS.json
        excl_fset = frozenset(excl_map[q_donor])
        cohort_panel_recs = select_cohort_panel_for_donor(q_donor, cohort_records, excl_map)
        cohort_panel_dict = _dedup_panel_dict({x["source_id"]: cohort_aligned[x["source_id"]] for x in cohort_panel_recs})
        cohort_full_hashes = {x["core_sha256"] for x in cohort_panel_recs}

        # Arab+Saudi-only leave-donor-and-family-out cohort panel
        as_panel_recs = select_cohort_panel_for_donor(
            q_donor, cohort_records, excl_map, allowed_resources={"Arab", "Saudi"}
        )
        as_panel_dict = _dedup_panel_dict({x["source_id"]: cohort_aligned[x["source_id"]] for x in as_panel_recs})
        as_full_hashes = {x["core_sha256"] for x in as_panel_recs}

        # 5-fold CV leave-fold-out cohort panels (for exact reconciliation with RECOVERABILITY.json)
        cv5_excl_fset = frozenset(cv5_excl_map[q_donor])
        cv5_cohort_recs = select_cohort_panel_for_donor(q_donor, cohort_records, cv5_excl_map)
        cv5_as_recs = select_cohort_panel_for_donor(
            q_donor, cohort_records, cv5_excl_map, allowed_resources={"Arab", "Saudi"}
        )
        cv5_cohort_hashes = {x["core_sha256"] for x in cv5_cohort_recs}
        cv5_as_hashes = {x["core_sha256"] for x in cv5_as_recs}
        cv5_as_dict = _dedup_panel_dict({x["source_id"]: cohort_aligned[x["source_id"]] for x in cv5_as_recs})

        # Combined full panel: IPD first, then cohort_panel_dict
        full_panel_dict = dict(ipd_dedup)
        full_panel_seen = set(ipd_dedup_seqs)
        for k_name, v_seq in cohort_panel_dict.items():
            if v_seq not in full_panel_seen:
                full_panel_dict[k_name] = v_seq
                full_panel_seen.add(v_seq)

        res_full_150 = _cached_classify("full", excl_fset, qid, q_aln, full_panel_dict)
        full_lbl_L1, full_det_L1 = _l1_labels(res_full_150)

        res_ipd_150 = _cached_classify("ipd", frozenset(), qid, q_aln, ipd_dedup)
        ipd_lbl_L1, _ = _l1_labels(res_ipd_150)

        res_coh_150 = _cached_classify("cohort", excl_fset, qid, q_aln, cohort_panel_dict)
        coh_lbl_L1, _ = _l1_labels(res_coh_150)

        res_as_150 = _cached_classify("as", excl_fset, qid, q_aln, as_panel_dict)
        res_cv5_as_150 = _cached_classify("cv5_as", cv5_excl_fset, qid, q_aln, cv5_as_dict)

        priv_type_counts = Counter(pv.var_type for pv in res_full_150.private_variants)
        priv_summary_str = ";".join(f"{k}:{v}" for k, v in sorted(priv_type_counts.items())) if priv_type_counts else "none"

        rows.append(
            {
                "source_id": qid,
                "gene": gene,
                "donor": q_donor,
                "haplotype": r.get("haplotype", ""),
                "cohort": r["resource"],
                "ipd_absent": ipd_absent,
                "core_len_bp": len(r["core_sequence"]),
                # Full panel (IPD + FOLDS.json leave-out cohort)
                "full_label_L150": res_full_150.label,
                "full_detailed_L150": res_full_150.detailed_label,
                "full_k_L150": "" if res_full_150.k_switches is None else res_full_150.k_switches,
                "full_label_L1": full_lbl_L1,
                "full_detailed_L1": full_det_L1,
                "full_k_L1": "" if res_full_150.k_switches_L1 is None else res_full_150.k_switches_L1,
                "full_exact_full_length": (q_hash in ipd_full_hashes) or (q_hash in cohort_full_hashes),
                "full_segment_sources": ";".join(res_full_150.segment_sources),
                "full_segment_lengths_bp": ";".join(str(x) for x in res_full_150.segment_lengths_bp),
                "full_switch_msa_cols": ";".join(str(sw.switch_msa_col) for sw in res_full_150.switches),
                "full_switch_base_positions": ";".join(str(sw.switch_base_0based) for sw in res_full_150.switches),
                "full_switch_agree_intervals_bp": ";".join(str(sw.agreeing_interval_bp) for sw in res_full_150.switches),
                "full_switch_div_variants_pm150": ";".join(str(sw.divergence_variants_pm150) for sw in res_full_150.switches),
                "full_switch_div_cols_pm150": ";".join(str(sw.divergence_cols_pm150) for sw in res_full_150.switches),
                "full_n_private_variants": len(res_full_150.private_variants),
                "full_private_var_types": priv_summary_str,
                "full_post_mask_sublabel": res_full_150.post_mask_sublabel,
                # IPD-only panel
                "ipd_label_L150": res_ipd_150.label,
                "ipd_detailed_L150": res_ipd_150.detailed_label,
                "ipd_k_L150": "" if res_ipd_150.k_switches is None else res_ipd_150.k_switches,
                "ipd_label_L1": ipd_lbl_L1,
                "ipd_k_L1": "" if res_ipd_150.k_switches_L1 is None else res_ipd_150.k_switches_L1,
                "ipd_exact_full_length": q_hash in ipd_full_hashes,
                # Cohort-only panel (FOLDS.json leave-out)
                "cohort_label_L150": res_coh_150.label,
                "cohort_detailed_L150": res_coh_150.detailed_label,
                "cohort_k_L150": "" if res_coh_150.k_switches is None else res_coh_150.k_switches,
                "cohort_label_L1": coh_lbl_L1,
                "cohort_k_L1": "" if res_coh_150.k_switches_L1 is None else res_coh_150.k_switches_L1,
                "cohort_exact_full_length": q_hash in cohort_full_hashes,
                # Regional Arab+Saudi-only panels (FOLDS.json vs 5-fold CV for RECOVERABILITY.json reconciliation)
                "arab_saudi_label_L150": res_as_150.label,
                "arab_saudi_exact_full_length": q_hash in as_full_hashes,
                "cv5_arab_saudi_label_L150": res_cv5_as_150.label,
                "cv5_arab_saudi_exact_full_length": q_hash in cv5_as_hashes,
                "cv5_cohort_exact_full_length": q_hash in cv5_cohort_hashes,
            }
        )

    out_tsv = per_gene_dir / f"{gene}_mosaic.tsv"
    with out_tsv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"[{gene}] Finished Part B: wrote {len(rows)} rows to {out_tsv}", flush=True)
    return rows


def run_part_c_for_gene(
    gene: str,
    cohort_records: List[Dict],
    ipd_records: List[Dict],
    cohort_aligned: Dict[str, str],
    ipd_aligned: Dict[str, str],
    excl_map: Dict[str, Set[str]],
    ref_backbone: str,
    exon_col_intervals: List[Tuple[int, int]],
    precomputed_seq_events: Dict[str, Set[Tuple[int, int, str, str]]],
    per_gene_dir: Path,
    synth_dir: Path,
    n_positives: int = 200,
) -> None:
    """Run Part C: Build 200 positive synthetic mosaics + 200 negative controls (+1 private SNV) per gene."""
    gene_idx = GENES.index(gene)
    seed = 20261103 + gene_idx
    rng = random.Random(seed)
    print(f"[{gene}] Starting Part C (Synthetic Mosaics, seed={seed})...", flush=True)

    synth_dir.mkdir(parents=True, exist_ok=True)
    ipd_dedup = _dedup_panel_dict(ipd_aligned)
    ipd_aligned_set = set(ipd_dedup.values())

    # Identify HPRC records whose aligned sequence is still present in IPD or in other unrelated cohort components
    # so that excluding the two parent donors does not make a parent sequence itself donor-private.
    hprc_records = [r for r in cohort_records if r["resource"] == "HPRC"]

    # Map each unique aligned sequence to the set of FOLDS.json root components that carry it
    seq_to_donors: Dict[str, Set[str]] = defaultdict(set)
    for r in cohort_records:
        seq_to_donors[cohort_aligned[r["source_id"]]].add(r["donor"])

    eligible_hprc: List[Dict] = []
    seen_hprc_aln: Set[str] = set()
    for r in hprc_records:
        aln = cohort_aligned[r["source_id"]]
        if aln in seen_hprc_aln:
            continue
        # Supported outside r's own component if present in IPD or in another donor outside excl_map[r["donor"]]
        other_donors = seq_to_donors[aln] - excl_map[r["donor"]]
        if (aln in ipd_aligned_set) or len(other_donors) >= 2:
            eligible_hprc.append(r)
            seen_hprc_aln.add(aln)

    # Precompute candidate pairs (hA, hB) from unrelated HPRC donors with matching outer MSA span
    candidate_pairs: List[Tuple[Dict, Dict, int, List[int]]] = []
    for i in range(len(eligible_hprc)):
        rA = eligible_hprc[i]
        dA = rA["donor"]
        alnA = cohort_aligned[rA["source_id"]]
        stA, enA = sequence_span_in_msa(alnA)
        for j in range(i + 1, len(eligible_hprc)):
            rB = eligible_hprc[j]
            dB = rB["donor"]
            if dB in excl_map[dA]:
                continue
            alnB = cohort_aligned[rB["source_id"]]
            if alnA == alnB:
                continue
            stB, enB = sequence_span_in_msa(alnB)
            if stA != stB or enA != enB:
                continue
            # Ensure both alnA and alnB remain available when (excl_map[dA] | excl_map[dB]) are excluded
            comb_excl = excl_map[dA] | excl_map[dB]
            if (alnA not in ipd_aligned_set) and not (seq_to_donors[alnA] - comb_excl):
                continue
            if (alnB not in ipd_aligned_set) and not (seq_to_donors[alnB] - comb_excl):
                continue

            # Compute differing MSA columns and variant event count between alnA and alnB
            diff_cols: List[int] = []
            n_var_diffs = 0
            in_diff = False
            for c in range(stA, enA + 1):
                ca = alnA[c]
                cb = alnB[c]
                if ca == "-" and cb == "-":
                    continue
                if ca != cb:
                    diff_cols.append(c)
                    if not in_diff:
                        n_var_diffs += 1
                        in_diff = True
                else:
                    in_diff = False
            if n_var_diffs >= 4:
                candidate_pairs.append((rA, rB, n_var_diffs, diff_cols))

    if not candidate_pairs:
        raise RuntimeError(f"No eligible HPRC haplotype pairs found for {gene}")

    # Stratify candidate pairs into 3 divergence tertiles (low, medium, high)
    candidate_pairs.sort(key=lambda x: x[2])
    n_cp = len(candidate_pairs)
    t1 = max(1, n_cp // 3)
    t2 = max(t1 + 1, (2 * n_cp) // 3)
    div_bins = {
        "low": candidate_pairs[:t1],
        "medium": candidate_pairs[t1:t2],
        "high": candidate_pairs[t2:],
    }

    # Build grid of (k_true, min_seg_target, div_stratum)
    grid: List[Tuple[int, int, str]] = []
    for k_t in (1, 2, 3):
        for l_min in (150, 300, 1000):
            for d_str in ("low", "medium", "high"):
                grid.append((k_t, l_min, d_str))

    def _try_splice_pair(
        rA: Dict,
        rB: Dict,
        diff_cols: List[int],
        k_true: int,
        min_seg_target: int,
        panel_sequences_set: Set[str],
    ) -> Optional[Dict]:
        alnA = cohort_aligned[rA["source_id"]]
        alnB = cohort_aligned[rB["source_id"]]
        st, en = sequence_span_in_msa(alnA)

        # Maximal agreeing column intervals between consecutive differing columns d_r and d_{r+1}
        # where both alnA and alnB agree for at least 1 non-gap base
        n_diffs = len(diff_cols)
        if n_diffs < k_true + 1:
            return None

        # Try up to 40 random combinations of k_true split indices among 0..n_diffs-2
        indices = list(range(n_diffs - 1))
        for _ in range(40):
            chosen = sorted(rng.sample(indices, k_true))
            # Pick switch column midpoint inside each agreeing interval (diff_cols[idx], diff_cols[idx+1])
            sw_cols: List[int] = []
            agree_col_intervals: List[Tuple[int, int]] = []
            valid_choice = True
            for idx_d in chosen:
                c_left = diff_cols[idx_d] + 1
                c_right = diff_cols[idx_d + 1] - 1
                if c_left > c_right:
                    valid_choice = False
                    break
                # Choose an agreeing column with a non-gap base in alnA/alnB
                non_gap_cands = [c for c in range(c_left, c_right + 1) if alnA[c] != "-"]
                if not non_gap_cands:
                    valid_choice = False
                    break
                mid_c = non_gap_cands[len(non_gap_cands) // 2]
                sw_cols.append(mid_c)
                agree_col_intervals.append((non_gap_cands[0], non_gap_cands[-1]))

            if not valid_choice:
                continue

            # Construct spliced aligned sequence alternating A, B, A, B...
            col_bounds = [0] + sw_cols + [len(alnA)]
            spliced_parts: List[str] = []
            for seg_i in range(k_true + 1):
                src_aln = alnA if (seg_i % 2 == 0) else alnB
                spliced_parts.append(src_aln[col_bounds[seg_i] : col_bounds[seg_i + 1]])
            spliced_aln = "".join(spliced_parts)

            # Compute ungapped segment lengths and base-level switch coordinates
            base_cols = [c for c in range(st, en + 1) if spliced_aln[c] != "-"]
            col_to_base = {c: b_idx for b_idx, c in enumerate(base_cols)}
            n_bases = len(base_cols)
            if n_bases < (k_true + 1) * min_seg_target:
                return None

            sw_bases = [col_to_base[c] for c in sw_cols]
            base_bounds = [0] + sw_bases + [n_bases]
            seg_lens = [base_bounds[i + 1] - base_bounds[i] for i in range(k_true + 1)]
            if any(sl < min_seg_target for sl in seg_lens):
                continue

            # Ensure spliced sequence is not already identical to a single existing sequence in the leave-out panel
            if spliced_aln in panel_sequences_set:
                continue

            agree_base_intervals = [
                (col_to_base[ac_s], col_to_base[ac_e]) for ac_s, ac_e in agree_col_intervals
            ]
            return {
                "spliced_aln": spliced_aln,
                "base_cols": base_cols,
                "sw_cols": sw_cols,
                "sw_bases": sw_bases,
                "seg_lens": seg_lens,
                "agree_col_intervals": agree_col_intervals,
                "agree_base_intervals": agree_base_intervals,
            }
        return None

    # Precompute observed bases at each MSA column across all sequences of the gene (any base absent globally
    # is guaranteed absent from every leave-pair-out panel P(q))
    all_gene_seqs = list(ipd_aligned_set | {cohort_aligned[r["source_id"]] for r in cohort_records})
    all_mat = np.vstack([np.frombuffer(s.encode("ascii"), dtype=np.uint8) for s in all_gene_seqs])
    m = len(ref_backbone)
    global_col_bases: Dict[int, Set[str]] = {}
    for b_ch in ("A", "C", "G", "T"):
        has_b = np.any(all_mat == ord(b_ch), axis=0)
        for c in np.where(has_b)[0]:
            global_col_bases.setdefault(int(c), set()).add(b_ch)

    # Cache leave-two-parents-out panel dicts by frozenset(excluded_donors)
    panel_cache: Dict[FrozenSet[str], Tuple[Dict[str, str], Set[str], Dict[int, Set[str]]]] = {}

    def _get_leave_pair_panel(dA: str, dB: str) -> Tuple[Dict[str, str], Set[str], Dict[int, Set[str]]]:
        excl = frozenset(excl_map[dA] | excl_map[dB])
        if excl in panel_cache:
            return panel_cache[excl]
        coh_recs = [
            r for r in cohort_records if r.get("eligible_reference", True) and r["donor"] not in excl
        ]
        verify_no_leakage(dA, set(excl), coh_recs)
        verify_no_leakage(dB, set(excl), coh_recs)
        p_dict = dict(ipd_dedup)
        p_seqs = set(ipd_aligned_set)
        for r in coh_recs:
            s_aln = cohort_aligned[r["source_id"]]
            if s_aln not in p_seqs:
                p_dict[r["source_id"]] = s_aln
                p_seqs.add(s_aln)

        panel_cache[excl] = (p_dict, p_seqs, global_col_bases)
        return p_dict, p_seqs, global_col_bases

    synth_rows: List[Dict] = []
    synth_fasta_lines: List[str] = []

    idx_generated = 0
    grid_idx = 0
    max_attempts = n_positives * 40
    attempts = 0

    while idx_generated < n_positives and attempts < max_attempts:
        attempts += 1
        k_true, l_min_target, d_str = grid[grid_idx % len(grid)]
        grid_idx += 1

        pool = div_bins[d_str] if div_bins[d_str] else candidate_pairs
        # Try up to 15 random pairs from pool; if l_min_target=1000 and k_true=3 is impossible for short class I genes,
        # fall back to largest feasible min_seg_target (e.g. 500 or 300) and record both target and effective min_seg
        eff_l_min = l_min_target
        spliced_info = None
        chosen_pair = None

        for fallback_l in (l_min_target, 600, 300, 150):
            for _ in range(15):
                rA, rB, n_var_diffs, diff_cols = rng.choice(pool)
                if rng.random() < 0.5:
                    rA, rB = rB, rA
                p_dict, p_seqs, col_bases = _get_leave_pair_panel(rA["donor"], rB["donor"])
                info = _try_splice_pair(rA, rB, diff_cols, k_true, fallback_l, p_seqs)
                if info is not None:
                    spliced_info = info
                    chosen_pair = (rA, rB, n_var_diffs, p_dict, col_bases)
                    eff_l_min = fallback_l
                    break
            if spliced_info is not None:
                break

        if spliced_info is None or chosen_pair is None:
            continue

        idx_generated += 1
        rA, rB, n_var_diffs, p_dict, col_bases = chosen_pair
        pos_id = f"{gene}_SYNTH_POS_{idx_generated:03d}"
        neg_id = f"{gene}_SYNTH_NEG_{idx_generated:03d}"

        pos_aln = spliced_info["spliced_aln"]
        pos_ungapped = "".join(c for c in pos_aln if c != "-")
        base_cols = spliced_info["base_cols"]
        sw_bases = spliced_info["sw_bases"]
        sw_cols = spliced_info["sw_cols"]
        agree_base_ivs = spliced_info["agree_base_intervals"]
        agree_col_ivs = spliced_info["agree_col_intervals"]

        # Create negative control: inject 1 random private SNV (absent from all of P(q)) inside a segment interior
        # Pick a base index at least 20 bp away from ends and switch points where an unused ACGT base exists
        avoid_bases: Set[int] = set(range(0, min(20, len(base_cols)))) | set(
            range(max(0, len(base_cols) - 20), len(base_cols))
        )
        for sb in sw_bases:
            avoid_bases.update(range(max(0, sb - 20), min(len(base_cols), sb + 21)))

        cand_mut_bases = [
            b_idx
            for b_idx, c in enumerate(base_cols)
            if b_idx not in avoid_bases and len(set("ACGT") - col_bases.get(c, set())) > 0
        ]
        if not cand_mut_bases:
            cand_mut_bases = [
                b_idx
                for b_idx, c in enumerate(base_cols)
                if len(set("ACGT") - col_bases.get(c, set())) > 0
            ]
        mut_base_idx = rng.choice(cand_mut_bases)
        mut_msa_col = base_cols[mut_base_idx]
        unused_acgt = sorted(set("ACGT") - col_bases.get(mut_msa_col, set()))
        priv_alt_base = rng.choice(unused_acgt)
        priv_ref_base = pos_aln[mut_msa_col]

        neg_aln_list = list(pos_aln)
        neg_aln_list[mut_msa_col] = priv_alt_base
        neg_aln = "".join(neg_aln_list)
        neg_ungapped = "".join(c for c in neg_aln if c != "-")

        synth_fasta_lines.append(f">{pos_id} gene={gene} k_true={k_true} control=positive\n{pos_ungapped}\n")
        synth_fasta_lines.append(f">{neg_id} gene={gene} k_true={k_true} control=negative_1snv\n{neg_ungapped}\n")

        # Run Part B classifier on positive control and negative control against leave-parents-out panel
        res_pos = classify_query_against_panel(
            gene, pos_id, pos_aln, p_dict, ref_backbone, exon_col_intervals, 150, precomputed_seq_events
        )
        res_neg = classify_query_against_panel(
            gene, neg_id, neg_aln, p_dict, ref_backbone, exon_col_intervals, 150, precomputed_seq_events
        )

        # Also run against the 2-parent panel {pA, pB} to verify intrinsic switch recovery without 3rd-allele shortcuts
        pair_only_dict = {
            rA["source_id"]: cohort_aligned[rA["source_id"]],
            rB["source_id"]: cohort_aligned[rB["source_id"]],
        }
        res_pos_2parent = classify_query_against_panel(
            gene, pos_id, pos_aln, pair_only_dict, ref_backbone, exon_col_intervals, 150, precomputed_seq_events
        )

        # Check switch position accuracy (within true agreeing interval)
        def _count_switches_in_agreeing_ivs(
            pred_switches: Sequence, true_agree_ivs: Sequence[Tuple[int, int]]
        ) -> int:
            matched = 0
            used_true: Set[int] = set()
            for sw in pred_switches:
                sb = sw.switch_base_0based
                for t_idx, (a_s, b_e) in enumerate(true_agree_ivs):
                    if t_idx not in used_true and a_s <= sb <= b_e:
                        matched += 1
                        used_true.add(t_idx)
                        break
            return matched

        n_sw_hit_full = _count_switches_in_agreeing_ivs(res_pos.switches, agree_base_ivs)
        n_sw_hit_2p = _count_switches_in_agreeing_ivs(res_pos_2parent.switches, agree_base_ivs)

        # Append positive row
        synth_rows.append(
            {
                "synth_id": pos_id,
                "gene": gene,
                "control_type": "positive",
                "seed": seed,
                "parent_a_id": rA["source_id"],
                "parent_a_donor": rA["donor"],
                "parent_b_id": rB["source_id"],
                "parent_b_donor": rB["donor"],
                "parent_divergence_variants": n_var_diffs,
                "divergence_stratum": d_str,
                "target_min_seg_bp": l_min_target,
                "effective_min_seg_bp": eff_l_min,
                "true_min_seg_bp": min(spliced_info["seg_lens"]),
                "true_k_switches": k_true,
                "true_switch_bases": ";".join(str(x) for x in sw_bases),
                "true_switch_msa_cols": ";".join(str(x) for x in sw_cols),
                "true_agree_base_intervals": ";".join(f"{a}-{b}" for a, b in agree_base_ivs),
                "true_agree_col_intervals": ";".join(f"{a}-{b}" for a, b in agree_col_ivs),
                "injected_private_snv_col": "",
                "injected_private_snv": "",
                "pred_label_full_panel": res_pos.label,
                "pred_detailed_full_panel": res_pos.detailed_label,
                "pred_k_full_panel": "" if res_pos.k_switches is None else res_pos.k_switches,
                "pred_switches_in_true_interval_full": n_sw_hit_full,
                "exact_k_match_full_panel": res_pos.k_switches == k_true,
                "all_switches_in_interval_full": (res_pos.k_switches == k_true) and (n_sw_hit_full == k_true),
                "pred_label_2parent_oracle": res_pos_2parent.label,
                "pred_k_2parent_oracle": "" if res_pos_2parent.k_switches is None else res_pos_2parent.k_switches,
                "pred_switches_in_true_interval_2parent": n_sw_hit_2p,
                "exact_k_match_2parent_oracle": res_pos_2parent.k_switches == k_true,
                "all_switches_in_interval_2parent": (res_pos_2parent.k_switches == k_true) and (n_sw_hit_2p == k_true),
                "pred_n_private_variants": len(res_pos.private_variants),
                "pred_post_mask_k": "" if res_pos.post_mask_k_switches is None else res_pos.post_mask_k_switches,
            }
        )

        # Append negative row
        synth_rows.append(
            {
                "synth_id": neg_id,
                "gene": gene,
                "control_type": "negative_1snv",
                "seed": seed,
                "parent_a_id": rA["source_id"],
                "parent_a_donor": rA["donor"],
                "parent_b_id": rB["source_id"],
                "parent_b_donor": rB["donor"],
                "parent_divergence_variants": n_var_diffs,
                "divergence_stratum": d_str,
                "target_min_seg_bp": l_min_target,
                "effective_min_seg_bp": eff_l_min,
                "true_min_seg_bp": min(spliced_info["seg_lens"]),
                "true_k_switches": k_true,
                "true_switch_bases": ";".join(str(x) for x in sw_bases),
                "true_switch_msa_cols": ";".join(str(x) for x in sw_cols),
                "true_agree_base_intervals": ";".join(f"{a}-{b}" for a, b in agree_base_ivs),
                "true_agree_col_intervals": ";".join(f"{a}-{b}" for a, b in agree_col_ivs),
                "injected_private_snv_col": mut_msa_col,
                "injected_private_snv": f"{priv_ref_base}>{priv_alt_base}",
                "pred_label_full_panel": res_neg.label,
                "pred_detailed_full_panel": res_neg.detailed_label,
                "pred_k_full_panel": "" if res_neg.k_switches is None else res_neg.k_switches,
                "pred_switches_in_true_interval_full": 0,
                "exact_k_match_full_panel": False,
                "all_switches_in_interval_full": False,
                "pred_label_2parent_oracle": "near-mosaic",
                "pred_k_2parent_oracle": "",
                "pred_switches_in_true_interval_2parent": 0,
                "exact_k_match_2parent_oracle": False,
                "all_switches_in_interval_2parent": False,
                "pred_n_private_variants": len(res_neg.private_variants),
                "pred_post_mask_k": "" if res_neg.post_mask_k_switches is None else res_neg.post_mask_k_switches,
            }
        )

    assert idx_generated == n_positives, f"Only generated {idx_generated}/{n_positives} synthetic mosaics for {gene}"

    # Write FASTA on IBEX only ($S/results/synthetic/<gene>_synthetic.fa)
    synth_fa = synth_dir / f"{gene}_synthetic.fa"
    synth_fa.write_text("".join(synth_fasta_lines))

    # Write per-gene synthetic truth & evaluation TSV
    synth_tsv = per_gene_dir / f"{gene}_synthetic_truth.tsv"
    with synth_tsv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(synth_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(synth_rows)

    n_neg_near = sum(
        1 for r in synth_rows if r["control_type"] == "negative_1snv" and r["pred_label_full_panel"] == "near-mosaic"
    )
    n_pos_mosaic = sum(
        1 for r in synth_rows if r["control_type"] == "positive" and r["pred_label_full_panel"] == "mosaic"
    )
    print(
        f"[{gene}] Finished Part C: {n_positives} positive ({n_pos_mosaic}/{n_positives} classified mosaic) "
        f"+ {n_positives} negative ({n_neg_near}/{n_positives} classified near-mosaic).",
        flush=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Run Phase 1 per-gene pipeline (Parts A, B, C).")
    parser.add_argument("--gene", type=str, default=None, help="HLA gene name (A, B, C, DPA1, DPB1, DQA1, DQB1, DRB1)")
    parser.add_argument("--threads", type=int, default=16, help="Threads for MAFFT")
    parser.add_argument("--force-msa", action="store_true", help="Force recomputing MAFFT MSA")
    args = parser.parse_args()

    gene = args.gene
    if not gene:
        task_id = os.environ.get("SLURM_ARRAY_TASK_ID")
        if task_id is not None:
            gene = GENES[int(task_id)]
        else:
            raise ValueError("Must specify --gene or SLURM_ARRAY_TASK_ID")
    if gene not in GENES:
        raise ValueError(f"Unknown gene {gene}; expected one of {GENES}")

    msa_dir = STAGE_ROOT / "tmp" / "msa"
    per_gene_dir = STAGE_ROOT / "tmp" / "per_gene"
    synth_dir = STAGE_ROOT / "results" / "synthetic"
    per_gene_dir.mkdir(parents=True, exist_ok=True)

    mafft_bin = STAGE_ROOT / "env" / "bin" / "mafft"
    minimap2_bin = STAGE_ROOT / "env" / "bin" / "minimap2"

    all_records = load_source_records(STAGE_ROOT / "inputs" / "SOURCE_RECORDS.json")
    truth_donors = {r["donor"] for r in all_records}
    excl_map = build_exclusion_map_from_folds(STAGE_ROOT / "inputs" / "FOLDS.json", truth_donors=truth_donors)

    donor_to_fold, fold_excluded, _ = load_augmented_fold_grouping(AUGMENTED_ROOT)
    cv5_excl_map = {d: fold_excluded[donor_to_fold[d]] for d in truth_donors}

    ipd_by_gene = load_ipd_reference(GENOMIC_FA)
    gene_cohort_records = [r for r in all_records if r["gene"] == gene]
    gene_ipd_records = ipd_by_gene[gene]

    cohort_aligned, ipd_aligned, uniq_aligned, stats = build_or_load_gene_msa(
        gene=gene,
        cohort_records=gene_cohort_records,
        ipd_records=gene_ipd_records,
        msa_dir=msa_dir,
        mafft_bin=mafft_bin,
        threads=args.threads,
        force=args.force_msa,
    )

    ref_backbone, exon_col_intervals, precomputed_seq_events = run_part_a_for_gene(
        gene=gene,
        cohort_records=gene_cohort_records,
        ipd_records=gene_ipd_records,
        cohort_aligned=cohort_aligned,
        ipd_aligned=ipd_aligned,
        uniq_aligned=uniq_aligned,
        canonical_allele=stats["canonical_allele"],
        minimap2_bin=minimap2_bin,
        per_gene_dir=per_gene_dir,
        tmp_dir=STAGE_ROOT / "tmp",
    )

    run_part_b_for_gene(
        gene=gene,
        cohort_records=gene_cohort_records,
        ipd_records=gene_ipd_records,
        cohort_aligned=cohort_aligned,
        ipd_aligned=ipd_aligned,
        excl_map=excl_map,
        cv5_excl_map=cv5_excl_map,
        ref_backbone=ref_backbone,
        exon_col_intervals=exon_col_intervals,
        precomputed_seq_events=precomputed_seq_events,
        per_gene_dir=per_gene_dir,
    )

    run_part_c_for_gene(
        gene=gene,
        cohort_records=gene_cohort_records,
        ipd_records=gene_ipd_records,
        cohort_aligned=cohort_aligned,
        ipd_aligned=ipd_aligned,
        excl_map=excl_map,
        ref_backbone=ref_backbone,
        exon_col_intervals=exon_col_intervals,
        precomputed_seq_events=precomputed_seq_events,
        per_gene_dir=per_gene_dir,
        synth_dir=synth_dir,
        n_positives=200,
    )

    print(f"[{gene}] ALL STEPS (MSA, Part A, Part B, Part C) COMPLETED SUCCESSFULLY.", flush=True)


if __name__ == "__main__":
    main()
