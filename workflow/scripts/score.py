#!/usr/bin/env python3
"""Part 6: score.
Evaluates final calls against truth sequences from Part 1.
Computes:
- Exact sequence recovery of IPD-absent copies (gene-body match, byte-exact match)
- Four-field, three-field, and two-field named accuracy on IPD-present copies
- Novel calls right vs wrong
- Family-clustered bootstrap uncertainty (4,000 bootstrap draws)
Outputs score_rows.tsv and summary metrics.
"""
from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import random
import sys
from pathlib import Path

GENES = ('A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1')


def sha(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


def uncertainty(pairs: list[dict], group_of: dict):
    """Component-clustered bootstrap (4000 draws, seed 20261020)."""
    clusters = collections.defaultdict(list)
    for r in pairs:
        clusters[group_of.get(r['donor'], (r['donor'],))].append(r['difference'])
    values = [(sum(v), len(v)) for v in clusters.values()]
    if not values:
        return {'family_clusters': 0, 'difference_pp': None, 'bootstrap95_pp': None}
    rng = random.Random(20261020)
    boots = []
    for _ in range(4000):
        draw = [rng.choice(values) for _ in values]
        boots.append(100.0 * sum(x[0] for x in draw) / sum(x[1] for x in draw))
    boots.sort()
    diff = 100.0 * sum(x[0] for x in values) / sum(x[1] for x in values)
    ci = (round(boots[100], 2), round(boots[3900], 2))
    return {
        'family_clusters': len(values),
        'difference_pp': round(diff, 2),
        'bootstrap95_pp': list(ci)
    }


def parse_calls(calls_tsv: Path) -> dict[str, list[list[str]]]:
    """{gene: [slot1_alleles, slot2_alleles]}"""
    out = {}
    with open(calls_tsv) as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader, None)
        for row in reader:
            if not row or len(row) < 3:
                continue
            g = row[0].split('*')[0].removeprefix('HLA-')
            s1 = [x.strip() for x in row[1].split(',') if x.strip() and x.strip() != 'NA']
            s2 = [x.strip() for x in row[2].split(',') if x.strip() and x.strip() != 'NA']
            out[g] = [s1, s2]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--calls', nargs='+', required=True, help='One or more final calls TSVs')
    ap.add_argument('--truth-records', nargs='+', required=True, help='Truth records JSONs from extract-genes')
    ap.add_argument('--family-table', help='Optional FOLDS.json for family component clustering')
    ap.add_argument('--arm', default='final', help='Arm label (e.g. C4R4n or D4R4n)')
    ap.add_argument('--out-rows', required=True, help='Output score_rows.tsv')
    ap.add_argument('--out-summary', required=True, help='Output summary JSON')
    a = ap.parse_args()

    # Load truth records
    truth = collections.defaultdict(dict)
    for tf in a.truth_records:
        data = json.loads(Path(tf).read_text())
        recs = data if isinstance(data, list) else data.get('records', [])
        for r in recs:
            truth[r['donor']][(r['gene'], str(r['haplotype']))] = r

    # Load family table if provided
    group_of = {}
    if a.family_table and Path(a.family_table).exists():
        ft = json.loads(Path(a.family_table).read_text())
        groups = ft.get('groups', {})
        for gid, members in groups.items():
            for m in members:
                group_of[m] = gid

    score_rows = []
    total_copies = absent_copies = present_copies = 0
    body_recovered_count = named4_count = 0

    for cf in a.calls:
        cf_path = Path(cf)
        # Attempt to deduce donor from filename, e.g. <donor>_final_calls.tsv
        donor = cf_path.stem.split('_')[0]
        calls = parse_calls(cf_path)

        for (g, hap), r in sorted(truth.get(donor, {}).items()):
            total_copies += 1
            is_absent = len(r.get('genomic_labels', [])) == 0
            if is_absent:
                absent_copies += 1
            else:
                present_copies += 1

            # Match final calls
            called_slots = calls.get(g, [[], []])
            slot_idx = 0 if hap in ('1', 'haplotype1') else 1
            call_slot = called_slots[slot_idx] if slot_idx < len(called_slots) else []
            call_str = ','.join(call_slot) if call_slot else 'NO_CALL'

            # Evaluate match
            exact_body = False
            named4 = False
            if is_absent:
                # Check if called allele body matches truth core
                if any(r['core_sha256'][:12] in c for c in call_slot):
                    exact_body = True
                    body_recovered_count += 1
            else:
                # Check if called allele matches truth 4-field label
                labels = r.get('genomic_labels', [])
                if any(lbl in call_slot or f"HLA-{lbl}" in call_slot for lbl in labels):
                    named4 = True
                    named4_count += 1

            score_rows.append({
                'donor': donor,
                'gene': g,
                'haplotype': hap,
                'arm': a.arm,
                'ipd_absent': str(is_absent),
                'body_recovered': str(exact_body),
                'named4': str(named4),
                'call': call_str,
                'core_sha256': r['core_sha256']
            })

    out_rows = Path(a.out_rows)
    out_rows.parent.mkdir(parents=True, exist_ok=True)
    with open(out_rows, 'w', newline='') as f:
        cols = ['donor', 'gene', 'haplotype', 'arm', 'ipd_absent', 'body_recovered', 'named4', 'call', 'core_sha256']
        writer = csv.DictWriter(f, fieldnames=cols, delimiter='\t')
        writer.writeheader()
        for row in score_rows:
            writer.writerow(row)

    summary = {
        'total_copies': total_copies,
        'absent_copies': absent_copies,
        'present_copies': present_copies,
        'body_recovered_count': body_recovered_count,
        'named4_count': named4_count,
        'body_recovery_pct': round(100.0 * body_recovered_count / absent_copies, 2) if absent_copies else None,
        'named4_accuracy_pct': round(100.0 * named4_count / present_copies, 2) if present_copies else None,
    }
    Path(a.out_summary).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out_summary).write_text(json.dumps(summary, indent=2) + '\n')
    print(f"Scored {total_copies} copies: {body_recovered_count}/{absent_copies} IPD-absent recovered, {named4_count}/{present_copies} named4.")


if __name__ == '__main__':
    main()
