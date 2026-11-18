#!/usr/bin/env python3
"""Part 3 helper: build-bait.
Builds the canonical 31-mer bait set from IPD genomic alleles minus those present
anywhere in GRCh38 outside the MHC.
Outputs bait_kmers.txt (expected sha256: 50b3b7324a5f5bca95272e69f83ba6cec592663363a79d549f551a209bb78170).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

K = 31
COMP = str.maketrans('ACGT', 'TGCA')


def canon(k: str) -> str:
    r = k[::-1].translate(COMP)
    return min(k, r)


def read_fasta_seqs(path: Path):
    parts = []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if parts:
                    yield ''.join(parts).upper()
                    parts = []
            else:
                parts.append(line.strip())
        if parts:
            yield ''.join(parts).upper()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ipd-genomic', required=True, help='Path to IPD genomic.fa')
    ap.add_argument('--grch38-ref', required=True, help='GRCh38 full analysis set fasta')
    ap.add_argument('--out-bait', required=True, help='Output bait_kmers.txt')
    ap.add_argument('--out-json', required=True, help='Output BAIT.json')
    ap.add_argument('--work-dir', default='tmp_bait', help='Temporary work directory')
    a = ap.parse_args()

    work = Path(a.work_dir)
    work.mkdir(parents=True, exist_ok=True)

    print("Extracting canonical 31-mers from IPD genomic alleles...")
    ks = set()
    n_alleles = 0
    for seq in read_fasta_seqs(Path(a.ipd_genomic)):
        n_alleles += 1
        for i in range(len(seq) - K + 1):
            k = seq[i : i + K]
            if 'N' not in k:
                ks.add(canon(k))

    ipd_kmers_fa = work / 'ipd_kmers.fa'
    with open(ipd_kmers_fa, 'w') as f:
        for i, k in enumerate(sorted(ks)):
            f.write(f'>{i}\n{k}\n')

    print(f"Collected {len(ks)} IPD 31-mers from {n_alleles} alleles.")

    masked_ref = work / 'ref_outside_mhc.fa'
    print("Masking MHC region (chr6:28,500,000-33,500,000) and filtering alt contigs from GRCh38...")
    with open(a.grch38_ref) as fi, open(masked_ref, 'w') as fo:
        keep, name, pos = False, None, 0
        for line in fi:
            if line.startswith('>'):
                name = line[1:].split()[0]
                keep = not (name.startswith('HLA-') or name.endswith('_alt'))
                pos = 0
                if keep:
                    fo.write(line)
                continue
            if not keep:
                continue
            s = line.rstrip('\n')
            if name == 'chr6' and pos + len(s) > 28500000 and pos < 33500000:
                s = ''.join('N' if 28500000 <= pos + j < 33500000 else c for j, c in enumerate(s))
            pos += len(s)
            fo.write(s + '\n')

    found_jf = work / 'found.jf'
    found_txt = work / 'found_outside_mhc.txt'

    print("Running jellyfish count against masked GRCh38...")
    n = len(ks)
    cmd1 = ['jellyfish', 'count', '-m', '31', '-C', '-t', '8', '-s', str(n * 2),
            '--if', str(ipd_kmers_fa), '-o', str(found_jf), str(masked_ref)]
    subprocess.run(cmd1, check=True)

    print("Dumping matches...")
    with open(found_txt, 'w') as fo:
        subprocess.run(['jellyfish', 'dump', '-c', str(found_jf)], stdout=fo, check=True)

    found = set()
    with open(found_txt) as f:
        for line in f:
            k, c = line.split()
            if int(c) > 0:
                found.add(canon(k.upper()))

    print(f"Found {len(found)} IPD 31-mers occurring outside MHC.")
    n_keep = 0
    out_bait = Path(a.out_bait)
    out_bait.parent.mkdir(parents=True, exist_ok=True)
    with open(ipd_kmers_fa) as f, open(out_bait, 'w') as fo:
        for line in f:
            if line.startswith('>'):
                continue
            k = line.strip()
            if k not in found:
                fo.write(k + '\n')
                n_keep += 1

    digest = hashlib.sha256(out_bait.read_bytes()).hexdigest()
    info = {
        'ipd_alleles': n_alleles,
        'ipd_kmers': len(ks),
        'k': K,
        'reference': str(a.grch38_ref),
        'masked': 'chr6:28500000-33500000 (0-based half-open)',
        'excluded': 'HLA-* and *_alt contigs',
        'ipd_kmers_found_outside_mhc': len(found),
        'bait_kmers': n_keep,
        'bait_sha256': digest
    }
    Path(a.out_json).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out_json).write_text(json.dumps(info, indent=2) + '\n')

    # Cleanup temporary large files
    if masked_ref.exists():
        masked_ref.unlink()
    if found_jf.exists():
        found_jf.unlink()

    print(f"Bait set written: {n_keep} k-mers, sha256: {digest}")


if __name__ == '__main__':
    main()
