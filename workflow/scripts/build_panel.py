#!/usr/bin/env python3
"""Part 2: build-panel.
Constructs a deduplicated composite T1K reference (arm B) from pinned IPD genomic alleles
and cohort HLA gene records from Part 1.
Excludes:
- Samples on the exclusion list (e.g. sample's own family / pedigree component)
- Alleles with missing coding exons
- Alleles containing ambiguous 'N' bases
- Alleles whose core sequence is already present in IPD
- Duplicate identical alleles across samples
Generates internal PG identifiers and T1K-compatible headers.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sys
from pathlib import Path

GENES = ('A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1')


def sha(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


def read_fasta_headers_and_seqs(path: Path):
    name, ints, parts = None, [], []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if name:
                    yield name, ints, ''.join(parts).upper()
                toks = line[1:].strip().split()
                name = toks[0]
                ints = [int(x) for x in toks[2:]] if len(toks) > 2 else []
                parts = []
            else:
                parts.append(line.strip())
        if name:
            yield name, ints, ''.join(parts).upper()


def t1k_header(name: str, exons: list[list[int]], length: int) -> str:
    prev = -1
    pos = []
    if not exons:
        raise ValueError('missing coding annotation')
    for s, e in exons:
        if not (0 <= s < e <= length) or s < prev:
            raise ValueError('invalid coding interval')
        pos.extend((s, e - 1))
        prev = e
    return '>' + name + ' ' + str(len(exons)) + ' ' + ' '.join(map(str, pos)) + '\n'


def build_additions(records: list[dict], ipd_seqs: set[str]):
    adds, seen = [], {}
    stats = dict(total_input=len(records), skipped_noexon=0, skipped_N=0, skipped_ipd=0, dup=0, added=0)

    for r in sorted(records, key=lambda x: (x['gene'], x['core_sequence'])):
        seq = r['core_sequence'].upper()
        if not r.get('core_exons'):
            stats['skipped_noexon'] += 1
            continue
        if 'N' in seq:
            stats['skipped_N'] += 1
            continue
        if seq in ipd_seqs:
            stats['skipped_ipd'] += 1
            continue
        key = (r['gene'], seq)
        h = sha(seq)
        if key in seen:
            stats['dup'] += 1
            continue
        name = f"HLA-{r['gene']}*PG{h[:12]}"
        ex = [list(e) for e in r['core_exons']]
        hdr = t1k_header(name, ex, len(seq))
        item = dict(
            name=name,
            gene=r['gene'],
            sha=h,
            exons=ex,
            length=len(seq),
            donor=r.get('donor', 'UNKNOWN'),
            fa=hdr + seq + '\n'
        )
        seen[key] = item
        adds.append(item)

    stats['added'] = len(adds)
    return adds, stats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ipd-genomic', required=True, help='Path to IPD genomic.fa')
    ap.add_argument('--records', nargs='+', required=True, help='One or more JSON files from extract-genes')
    ap.add_argument('--exclude-samples', default='', help='Comma-separated list of sample IDs to exclude')
    ap.add_argument('--out-reference', required=True, help='Output reference.fa')
    ap.add_argument('--out-aliases', required=True, help='Output aliases.json')
    ap.add_argument('--out-summary', required=True, help='Output summary JSON')
    a = ap.parse_args()

    exclude_set = set(filter(None, a.exclude_samples.split(',')))

    ipd_seqs = set()
    for name, ints, seq in read_fasta_headers_and_seqs(Path(a.ipd_genomic)):
        ipd_seqs.add(seq.upper())

    all_records = []
    for rf in a.records:
        data = json.loads(Path(rf).read_text())
        if isinstance(data, list):
            for r in data:
                if r.get('donor') not in exclude_set:
                    all_records.append(r)
        elif isinstance(data, dict) and 'records' in data:
            for r in data['records']:
                if r.get('donor') not in exclude_set:
                    all_records.append(r)

    adds, stats = build_additions(all_records, ipd_seqs)

    out_ref = Path(a.out_reference)
    out_ref.parent.mkdir(parents=True, exist_ok=True)
    with open(out_ref, 'wb') as fo:
        with open(a.ipd_genomic, 'rb') as fi:
            shutil.copyfileobj(fi, fo)
        for add in adds:
            fo.write(add['fa'].encode())

    aliases = {
        add['name']: {
            'gene': add['gene'],
            'sequence_sha256': add['sha'],
            'exons': add['exons'],
            'length': add['length'],
            'donor': add['donor'],
            'internal_identifier': True
        }
        for add in adds
    }
    Path(a.out_aliases).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out_aliases).write_text(json.dumps(aliases, indent=2) + '\n')

    stats['reference_sha256'] = sha(out_ref.read_text())
    stats['excluded_samples'] = sorted(exclude_set)
    Path(a.out_summary).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out_summary).write_text(json.dumps(stats, indent=2) + '\n')

    print(f"Panel reference built: {stats['added']} novel alleles added (from {stats['total_input']} input copies)")


if __name__ == '__main__':
    main()
