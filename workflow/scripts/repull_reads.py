#!/usr/bin/env python3
"""Part 3: repull-reads.

Produces the rule-2 read set: the regional MHC recruitment (with mate rescue) UNION every
whole-genome read pair that carries at least 3 bait 31-mers.

Start from a CRAM/BAM alone (needs its index, and the .fai of the reference), or hand in the
regional reads and/or the T1K candidate reads that were made earlier.

The method is the accepted one and is not reimplemented here:
  * regional recruitment: `extract()` and `test_rescue()` of the vendored recruitment code
    (dogohla/frozen/confirmatory-test-20261113/recruit_test.py, copied unchanged from the
    recruitment used for the development donors), with the same region rule: chr6:28,500,000-33,500,000
    plus every HLA-* and chr6_* contig of the CRAM header over its full length, indexed mate rescue,
    no MAPQ filter;
  * candidate union: `fq`, `qname` and `hits` of the vendored rule-2 code
    (dogohla/frozen/reads-repull-20261114/select3.py); same loop as select_rule2.py.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import io
import json
import subprocess
import sys
from pathlib import Path


def load_frozen(package_dir: Path):
    """Check the vendored files against MANIFEST.sha256 and put the frozen folders on sys.path."""
    spec = importlib.util.spec_from_file_location('dogohla_verify', str(package_dir / 'verify.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    print('frozen code verified:', mod.verify(), flush=True)
    mod.use()


def sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 22), b''):
            h.update(chunk)
    return h.hexdigest()


class GzWriter:
    """gzip text writer with mtime 0, so that equal content gives equal bytes."""

    def __init__(self, path: Path):
        self.raw = open(path, 'wb')
        self.gz = gzip.GzipFile(filename='', mode='wb', fileobj=self.raw, compresslevel=4, mtime=0)
        self.txt = io.TextIOWrapper(self.gz, encoding='utf-8', newline='')

    def write(self, s):
        self.txt.write(s)

    def close(self):
        self.txt.close()
        self.raw.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--donor', required=True, help='Donor/sample ID')
    ap.add_argument('--cram', help='Input WGS CRAM or BAM file (its index must sit next to it)')
    ap.add_argument('--reference-genome', help='GRCh38 reference FASTA for CRAM decoding (its .fai next to it)')
    ap.add_argument('--ipd-genomic', help='IPD genomic.fa for fastq-extractor')
    ap.add_argument('--regional-r1', help='Existing regional recruitment FASTQ R1')
    ap.add_argument('--regional-r2', help='Existing regional recruitment FASTQ R2')
    ap.add_argument('--candidate-r1', help='Pre-extracted T1K candidate FASTQ R1')
    ap.add_argument('--candidate-r2', help='Pre-extracted T1K candidate FASTQ R2')
    ap.add_argument('--bait', required=True, help='Path to bait_kmers.txt')
    ap.add_argument('--out-r1', required=True, help='Output repulled r1.fq.gz')
    ap.add_argument('--out-r2', required=True, help='Output repulled r2.fq.gz')
    ap.add_argument('--out-pairs', required=True, help='Output PAIRS.json')
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--work-dir', default='tmp_repull', help='Work directory for temporaries')
    ap.add_argument('--package-dir', default=str(Path(__file__).resolve().parent.parent / 'dogohla'),
                    help='The workflow package folder with verify.py and frozen/')
    a = ap.parse_args()

    load_frozen(Path(a.package_dir))
    import pysam                      # noqa: E402
    import recruit_test as RT         # vendored recruitment code (extract, test_rescue)  # noqa: E402
    import common as _common          # noqa: F401,E402  (load the vendored copies before select3 touches sys.path)
    import run_donor as _run_donor    # noqa: F401,E402
    import select3 as S               # vendored rule-2 code (fq, qname, hits)  # noqa: E402

    work = Path(a.work_dir)
    work.mkdir(parents=True, exist_ok=True)

    print('Loading bait 31-mers...', flush=True)
    bait = set(l.strip() for l in open(a.bait))
    bait_sha = sha(Path(a.bait))

    reg_r1, reg_r2 = a.regional_r1, a.regional_r2
    cand_r1, cand_r2 = a.candidate_r1, a.candidate_r2
    recruit_info = {}
    cram_sha = None

    if a.cram:
        assert a.reference_genome, 'CRAM extraction requires --reference-genome'
        cram = Path(a.cram)
        ref = Path(a.reference_genome)
        cram_sha = sha(cram)

        # A: regional MHC recruitment with mate rescue (vendored recruit code)
        if not (reg_r1 and reg_r2):
            print('Regional recruitment from the CRAM (region bed + mate rescue)...', flush=True)
            with open(work / 'recruit.log', 'wb') as log:
                RT.test_rescue('samtools', log, work)
                with pysam.AlignmentFile(str(cram), 'rc', reference_filename=str(ref)) as bam:
                    regions = [('chr6', 28500000, 33500000)]
                    regions += [(n, 0, l) for n, l in zip(bam.references, bam.lengths)
                                if n.startswith('HLA-') or n.startswith('chr6_')]
                bed = work / 'recruitment.bed'
                bed.write_text(''.join(f'{n}\t{s}\t{e}\n' for n, s, e in regions))
                metrics = RT.extract('samtools', cram, bed, work / 'reads', ref, log)
            reg_r1 = work / 'reads' / 'r1.fq.gz'
            reg_r2 = work / 'reads' / 'r2.fq.gz'
            recruit_info = dict(recruitment_bed_sha256=sha(bed), recruitment_regions=len(regions),
                                recruitment_metrics=metrics)
            for x in ('selected.bam', 'collated.bam'):
                (work / 'reads' / x).unlink(missing_ok=True)

        # B: whole-genome T1K candidate extraction
        if not (cand_r1 and cand_r2):
            assert a.ipd_genomic, 'candidate extraction requires --ipd-genomic'
            print('Whole-genome candidate reads via fastq-extractor...', flush=True)
            wgs_r1 = work / f'{a.donor}_wgs_r1.fq.gz'
            wgs_r2 = work / f'{a.donor}_wgs_r2.fq.gz'
            cmd_wgs_fq = (f'samtools collate -@ 3 -u -O --reference {ref} {cram} {work}/collate-wgs | '
                          f'samtools fastq -@ 3 -F 0x900 -n -1 {wgs_r1} -2 {wgs_r2} -0 /dev/null -s /dev/null -')
            subprocess.run(['bash', '-c', 'set -o pipefail; ' + cmd_wgs_fq], check=True)
            cand_pref = work / f'{a.donor}_calls_candidate'
            subprocess.run(['fastq-extractor', '-t', str(a.threads), '-f', a.ipd_genomic, '-o', str(cand_pref),
                            '-s', '0.97', '-1', str(wgs_r1), '-2', str(wgs_r2)], check=True)
            cand_r1 = f'{cand_pref}_1.fq'
            cand_r2 = f'{cand_pref}_2.fq'
    else:
        assert reg_r1 and reg_r2 and cand_r1 and cand_r2, \
            'without a CRAM give --regional-r1/-r2 and --candidate-r1/-r2'

    # Rule-2 union (same loop as select_rule2.py)
    print('Performing rule-2 union...', flush=True)
    out_r1, out_r2 = Path(a.out_r1), Path(a.out_r2)
    out_r1.parent.mkdir(parents=True, exist_ok=True)
    names = set()
    n_old = n_add = n_cand = 0
    w1, w2 = GzWriter(out_r1), GzWriter(out_r2)
    try:
        for (h1, s1, q1), (h2, s2, q2) in zip(S.fq(Path(reg_r1)), S.fq(Path(reg_r2))):
            n1 = S.qname(h1)
            assert n1 == S.qname(h2)
            names.add(n1)
            w1.write(f'{h1}{s1}+\n{q1}')
            w2.write(f'{h2}{s2}+\n{q2}')
            n_old += 1
        for (h1, s1, q1), (h2, s2, q2) in zip(S.fq(Path(cand_r1)), S.fq(Path(cand_r2))):
            n_cand += 1
            n1 = S.qname(h1)
            assert n1 == S.qname(h2)
            if n1 in names:
                continue
            if S.hits(s1, bait) or S.hits(s2, bait):
                names.add(n1)
                w1.write(f'{h1}{s1}+\n{q1}')
                w2.write(f'{h2}{s2}+\n{q2}')
                n_add += 1
    finally:
        w1.close()
        w2.close()

    pairs_rec = {
        'status': 'paired_fastq_verified',
        'donor': a.donor,
        'readset': 'rule2',
        'pairs': n_old + n_add,
        'old_regional_pairs': n_old,
        'added_pairs': n_add,
        'wgs_candidate_pairs_scanned': n_cand,
        'bait_sha256': bait_sha,
        'cram_sha256': cram_sha,
        'r1': str(out_r1),
        'r2': str(out_r2),
        'r1_sha256': sha(out_r1),
        'r2_sha256': sha(out_r2),
        **recruit_info,
    }
    Path(a.out_pairs).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out_pairs).write_text(json.dumps(pairs_rec, indent=2) + '\n')
    print(f'Repull complete for {a.donor}: {n_old} regional + {n_add} bait-rescued = {n_old + n_add} total pairs.')


if __name__ == '__main__':
    main()
