#!/usr/bin/env python3
"""Part 1: extract-genes.
Extracts 8 HLA genes (A, B, C, DPA1, DPB1, DQA1, DQB1, DRB1) from phased assemblies
using Immuannot annotations, 2,000 bp flank slicing (sequence[2000:-2000]), and CDS interval mapping.
Outputs per-sample FASTA + JSON records identical in schema to cohort SOURCE_RECORDS.json.
"""
from __future__ import annotations

import argparse
import collections
import gzip
import hashlib
import json
import os
import re
import sys
from pathlib import Path

GENES = ('A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1')
FLANK = 2000
COMP = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')


def sha(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


def parse_gtf_genes_and_cds(gtf_path: Path):
    """Parse gene coordinates and oriented 0-based exclusive CDS intervals from immuannot.gtf.gz."""
    genes = {}
    raw_cds = collections.defaultdict(list)
    op = gzip.open if str(gtf_path).endswith('.gz') else open

    with op(gtf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            contig, feature, start_s, end_s, strand, attrs = (
                parts[0], parts[2], int(parts[3]), int(parts[4]), parts[6], parts[8]
            )
            m_gid = re.search(r'gene_id "([^"]+)"', attrs)
            m_gname = re.search(r'gene_name "([^"]+)"', attrs)
            if not m_gid or not m_gname:
                continue
            gid = m_gid.group(1)
            gname = m_gname.group(1).removeprefix('HLA-')
            if gname not in GENES:
                continue

            if feature == 'gene':
                t_match = re.search(r'template_allele "([^"]+)"', attrs)
                d_match = re.search(r'template_distance ([^;]+);', attrs)
                genes[gid] = {
                    'contig': contig,
                    'gene': gname,
                    'start_1based': start_s,
                    'end_inclusive': end_s,
                    'strand': strand,
                    'template_allele': t_match.group(1) if t_match else None,
                    'template_distance': d_match.group(1) if d_match else None,
                }
            elif feature == 'CDS':
                raw_cds[gid].append((start_s, end_s))

    out_records = []
    for gid, gmeta in genes.items():
        gstart = gmeta['start_1based']
        gend = gmeta['end_inclusive']
        strand = gmeta['strand']
        intervals = []
        for c_start, c_end in raw_cds[gid]:
            if strand == '+':
                rel_start = c_start - gstart
                rel_end = c_end - gstart + 1
            else:
                rel_start = gend - c_end
                rel_end = gend - c_start + 1
            intervals.append((rel_start, rel_end))
        intervals.sort()
        gmeta['core_exons'] = intervals
        out_records.append(gmeta)

    return out_records


def read_fasta_dict(path: Path) -> dict[str, str]:
    seqs = {}
    name, parts = None, []
    op = gzip.open if str(path).endswith('.gz') else open
    with op(path, 'rt') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if name:
                    seqs[name] = ''.join(parts).upper()
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
        if name:
            seqs[name] = ''.join(parts).upper()
    return seqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--assembly', required=True, help='Assembly FASTA (diploid or hap-specific)')
    ap.add_argument('--gtf', required=True, help='Immuannot GTF annotation (.gtf or .gtf.gz)')
    ap.add_argument('--sample', required=True, help='Sample identifier')
    ap.add_argument('--haplotype', default='1', help='Haplotype identifier (1 or 2)')
    ap.add_argument('--out-fasta', required=True, help='Output FASTA')
    ap.add_argument('--out-records', required=True, help='Output JSON records')
    a = ap.parse_args()

    contigs = read_fasta_dict(Path(a.assembly))
    gene_metas = parse_gtf_genes_and_cds(Path(a.gtf))

    records = []
    fasta_entries = []

    copies_seen = collections.defaultdict(int)
    for gm in sorted(gene_metas, key=lambda x: (x['gene'], x['contig'], x['start_1based'])):
        g = gm['gene']
        ctg = gm['contig']
        if ctg not in contigs:
            continue
        copies_seen[g] += 1
        hap = a.haplotype if copies_seen[g] <= 1 else f"{a.haplotype}_{copies_seen[g]}"

        seq_full = contigs[ctg]
        L = len(seq_full)
        s1 = gm['start_1based'] - 1  # 0-based start
        e1 = gm['end_inclusive']     # 0-based exclusive end
        strand = gm['strand']

        # Slice locus with 2000 bp flanks
        lo = max(0, s1 - FLANK)
        hi = min(L, e1 + FLANK)
        flanked_raw = seq_full[lo:hi]

        # Orient locus in gene strand
        if strand == '+':
            flanked = flanked_raw
            start_in_flanked = s1 - lo
            core_len = e1 - s1
        else:
            flanked = flanked_raw.translate(COMP)[::-1]
            start_in_flanked = hi - e1
            core_len = e1 - s1

        # Strict 2000 bp flank check
        if start_in_flanked >= FLANK and start_in_flanked + core_len + FLANK <= len(flanked):
            flanked_exact = flanked[start_in_flanked - FLANK : start_in_flanked + core_len + FLANK]
            core = flanked_exact[FLANK:-FLANK]
        else:
            flanked_exact = flanked
            core = flanked[start_in_flanked : start_in_flanked + core_len]

        core_sha = sha(core)
        core_id = f"{g}*PG{core_sha}"
        rec = {
            'donor': a.sample,
            'gene': g,
            'haplotype': hap,
            'contig': ctg,
            'start_1based': gm['start_1based'],
            'end_inclusive': gm['end_inclusive'],
            'strand': strand,
            'sequence': flanked_exact,
            'core_sequence': core,
            'core_sha256': core_sha,
            'core_id': core_id,
            'core_exons': gm['core_exons'],
            'template_allele': gm.get('template_allele'),
            'template_distance': gm.get('template_distance'),
            'genomic_labels': [gm['template_allele']] if gm.get('template_allele') else [],
        }
        records.append(rec)
        fasta_entries.append(f">{a.sample}#{hap}#HLA-{g} {ctg}:{gm['start_1based']}-{gm['end_inclusive']}({strand})\n{flanked_exact}\n")

    Path(a.out_fasta).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out_fasta).write_text(''.join(fasta_entries))

    Path(a.out_records).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out_records).write_text(json.dumps(records, indent=2) + '\n')
    print(f"Extracted {len(records)} copies for {a.sample}")


if __name__ == '__main__':
    main()
