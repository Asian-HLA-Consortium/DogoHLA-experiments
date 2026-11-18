#!/usr/bin/env python3
"""Part 3: repull-reads.
Produces the rule-2 read set (MHC regional pull UNION whole-genome pairs with >= 3 bait 31-mers).
Supports extraction starting from WGS CRAM/BAM or from FASTQ inputs.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

K, STRIDE, MIN_HITS = 31, 4, 3
COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 22), b''):
            h.update(chunk)
    return h.hexdigest()


def fq(path: Path):
    op = gzip.open if str(path).endswith('.gz') else open
    with op(path, 'rt') as f:
        while True:
            h = f.readline()
            if not h:
                return
            s, p, q = f.readline(), f.readline(), f.readline()
            yield h, s, q


def qname(h: str) -> str:
    n = h[1:].split()[0]
    return n[:-2] if n.endswith(('/1', '/2')) else n


def hits(seq: str, bait: set[str]) -> bool:
    s = seq.strip().upper()
    r = s[::-1].translate(COMP)
    L = len(s)
    n = 0
    for i in range(0, L - K + 1, STRIDE):
        a = s[i : i + K]
        b = r[L - i - K : L - i]
        if (a if a < b else b) in bait:
            n += 1
            if n >= MIN_HITS:
                return True
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--donor', required=True, help='Donor/sample ID')
    ap.add_argument('--cram', help='Input WGS CRAM or BAM file')
    ap.add_argument('--reference-genome', help='GRCh38 reference FASTA for CRAM decoding')
    ap.add_argument('--ipd-genomic', help='IPD genomic.fa for fastq-extractor')
    ap.add_argument('--regional-r1', help='Existing regional pull FASTQ R1')
    ap.add_argument('--regional-r2', help='Existing regional pull FASTQ R2')
    ap.add_argument('--candidate-r1', help='Pre-extracted T1K candidate FASTQ R1')
    ap.add_argument('--candidate-r2', help='Pre-extracted T1K candidate FASTQ R2')
    ap.add_argument('--bait', required=True, help='Path to bait_kmers.txt')
    ap.add_argument('--out-r1', required=True, help='Output repulled r1.fq.gz')
    ap.add_argument('--out-r2', required=True, help='Output repulled r2.fq.gz')
    ap.add_argument('--out-pairs', required=True, help='Output PAIRS.json')
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--work-dir', default='tmp_repull', help='Work directory for temporaries')
    a = ap.parse_args()

    work = Path(a.work_dir)
    work.mkdir(parents=True, exist_ok=True)

    # 1. Load bait 31-mers
    print("Loading bait 31-mers...")
    bait = set(l.strip() for l in open(a.bait))
    bait_sha = sha(Path(a.bait))

    reg_r1 = a.regional_r1
    reg_r2 = a.regional_r2
    cand_r1 = a.candidate_r1
    cand_r2 = a.candidate_r2

    # If starting from CRAM:
    if a.cram:
        assert a.reference_genome and a.ipd_genomic, "CRAM extraction requires --reference-genome and --ipd-genomic"
        cram_p = Path(a.cram)
        cram_sha = sha(cram_p)

        # A: Extract regional MHC reads
        if not (reg_r1 and reg_r2):
            print("Extracting MHC regional reads from CRAM...")
            reg_r1 = work / f"{a.donor}_reg_r1.fq.gz"
            reg_r2 = work / f"{a.donor}_reg_r2.fq.gz"
            reg_bam = work / f"{a.donor}_mhc.bam"
            # chr6:28500000-33500000 plus HLA/alt contigs
            mhc_regions = "chr6:28500000-33500000"
            cmd_view = f"samtools view -h -b -F 0x900 --reference {a.reference_genome} {a.cram} {mhc_regions} > {reg_bam}"
            subprocess.run(cmd_view, shell=True, check=True)
            cmd_reg_fq = f"samtools collate -@ 3 -u -O {reg_bam} {work}/collate-reg | samtools fastq -@ 3 -F 0x900 -n -1 {reg_r1} -2 {reg_r2} -0 /dev/null -s /dev/null -"
            subprocess.run(cmd_reg_fq, shell=True, check=True)

        # B: Whole-genome T1K candidate extraction
        if not (cand_r1 and cand_r2):
            print("Extracting whole-genome candidate reads via fastq-extractor...")
            wgs_r1 = work / f"{a.donor}_wgs_r1.fq.gz"
            wgs_r2 = work / f"{a.donor}_wgs_r2.fq.gz"
            cmd_wgs_fq = f"samtools collate -@ 3 -u -O --reference {a.reference_genome} {a.cram} {work}/collate-wgs | samtools fastq -@ 3 -F 0x900 -n -1 {wgs_r1} -2 {wgs_r2} -0 /dev/null -s /dev/null -"
            subprocess.run(cmd_wgs_fq, shell=True, check=True)

            cand_pref = work / f"{a.donor}_calls_candidate"
            cmd_ext = f"fastq-extractor -t {a.threads} -f {a.ipd_genomic} -o {cand_pref} -s 0.97 -1 {wgs_r1} -2 {wgs_r2}"
            subprocess.run(cmd_ext, shell=True, check=True)
            cand_r1 = f"{cand_pref}_1.fq"
            cand_r2 = f"{cand_pref}_2.fq"

    # 2. Perform Rule-2 Union
    print("Performing rule-2 union...")
    out_r1 = Path(a.out_r1)
    out_r2 = Path(a.out_r2)
    out_r1.parent.mkdir(parents=True, exist_ok=True)

    names = set()
    n_old = n_add = n_cand = 0

    with gzip.open(out_r1, 'wt', compresslevel=4) as w1, gzip.open(out_r2, 'wt', compresslevel=4) as w2:
        for (h1, s1, q1), (h2, s2, q2) in zip(fq(Path(reg_r1)), fq(Path(reg_r2))):
            n1 = qname(h1)
            assert n1 == qname(h2)
            names.add(n1)
            w1.write(f'{h1}{s1}+\n{q1}')
            w2.write(f'{h2}{s2}+\n{q2}')
            n_old += 1

        for (h1, s1, q1), (h2, s2, q2) in zip(fq(Path(cand_r1)), fq(Path(cand_r2))):
            n_cand += 1
            n1 = qname(h1)
            assert n1 == qname(h2)
            if n1 in names:
                continue
            if hits(s1, bait) or hits(s2, bait):
                names.add(n1)
                w1.write(f'{h1}{s1}+\n{q1}')
                w2.write(f'{h2}{s2}+\n{q2}')
                n_add += 1

    r1_sha = sha(out_r1)
    r2_sha = sha(out_r2)

    pairs_rec = {
        'status': 'paired_fastq_verified',
        'donor': a.donor,
        'readset': 'rule2',
        'pairs': n_old + n_add,
        'old_regional_pairs': n_old,
        'added_pairs': n_add,
        'wgs_candidate_pairs_scanned': n_cand,
        'bait_sha256': bait_sha,
        'r1': str(out_r1),
        'r2': str(out_r2),
        'r1_sha256': r1_sha,
        'r2_sha256': r2_sha,
    }
    Path(a.out_pairs).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out_pairs).write_text(json.dumps(pairs_rec, indent=2) + '\n')
    print(f"Repull complete for {a.donor}: {n_old} regional + {n_add} bait-rescued = {n_old + n_add} total pairs.")


if __name__ == '__main__':
    main()
