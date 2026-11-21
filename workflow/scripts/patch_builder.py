#!/usr/bin/env python3
"""Part 5: patch-builder (thin wrapper around the frozen builder).

Runs the frozen all-gene builder (dogohla/frozen/mosaic-allgenes-20261114: drb1_diag2.load_shared and
drb1_run.run, per gene, with the confounded-window files) on the pass-1 T1K run of one arm, exactly as the
accepted equivalence wrapper did. No builder logic is written here.

  arm C: panel = IPD alleles, pass 1 = T1K on IPD genomic.fa       (final call = C4R4n)
  arm D: panel = IPD + the cohort panel,  pass 1 = T1K on the panel (final call = D4R4n)

Outputs in --out-dir: per gene proposals.json, additions.fa, OVERRIDES.json (arm C4/D4),
OVERRIDES_R4n.json (arm C4R4n/D4R4n; the abstain rule) and the other abstain variants, COMPLETE.json.
--out-final-calls: the pass-1 calls with the overrides of the requested arm applied
(frozen ag_score.new_post), one row per gene.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import shutil
import time
from pathlib import Path

GENES = ('A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1')
READSET = 'D-real3'     # label of the read set in the frozen code (re-pulled reads, rule 2)


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for c in iter(lambda: f.read(1 << 22), b''):
            h.update(c)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--donor', required=True)
    ap.add_argument('--arm', default='C4R4n', help='C4, C4R4n, D4 or D4R4n (C = IPD panel, D = cohort panel)')
    ap.add_argument('--pass1-calls', required=True, help='T1K genotype table of this arm')
    ap.add_argument('--pass1-candidate-r1', required=True, help='T1K candidate reads R1 of this arm')
    ap.add_argument('--pass1-candidate-r2', required=True, help='T1K candidate reads R2 of this arm')
    ap.add_argument('--reference', required=True, help='T1K reference FASTA of this arm (IPD, or IPD + panel)')
    ap.add_argument('--aliases', help='aliases.json of the cohort panel (arm D)')
    ap.add_argument('--ipd-genomic', required=True, help='IPD genomic.fa (sha256 pinned in the frozen code)')
    ap.add_argument('--builder-inputs', required=True, help='cohort resource folder (SOURCE_RECORDS.json, msa/, ...)')
    ap.add_argument('--params', required=True, help='frozen builder parameter file PARAMS.json')
    ap.add_argument('--confounded-dir', required=True, help='folder with confounded_<gene>.json')
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--out-final-calls', required=True)
    ap.add_argument('--threads', type=int, default=4, help='unused: the frozen builder fixes its own thread count')
    ap.add_argument('--package-dir', default=str(Path(__file__).resolve().parent.parent / 'dogohla'))
    a = ap.parse_args()

    arm = a.arm[0]                      # 'C' or 'D'
    assert arm in ('C', 'D') and a.arm in (arm + '4', arm + '4R4n'), a.arm
    out = Path(a.out_dir).resolve()
    work = Path('builder_work').resolve()
    work.mkdir(parents=True, exist_ok=True)

    spec = importlib.util.spec_from_file_location('dogohla_wrapper_env', str(Path(a.package_dir) / 'wrapper_env.py'))
    WE = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(WE)
    C, RD, D2, R, AS = WE.setup(a.package_dir, a.builder_inputs, a.ipd_genomic, a.params, work)

    # 1. pass-1 run folder in the layout the frozen code reads (COMPLETE.json + table + candidate reads)
    run_dir = work / f'pass1_{arm}'
    run_dir.mkdir(parents=True, exist_ok=True)
    for name, src in (('calls_genotype.tsv', a.pass1_calls), ('calls_candidate_1.fq', a.pass1_candidate_r1),
                      ('calls_candidate_2.fq', a.pass1_candidate_r2)):
        dst = run_dir / name
        if dst.is_symlink() or dst.exists():
            dst.unlink()
        dst.symlink_to(Path(src).resolve())
    ref = Path(a.reference).resolve()
    rec = {'status': 'complete', 'donor': a.donor, 'readset': READSET, 'arm': 'A' if arm == 'C' else 'B',
           'truth_loaded': False, 'reference_sha256': sha256_file(ref),
           'native_table_sha256': sha256_file(run_dir / 'calls_genotype.tsv')}
    (run_dir / 'COMPLETE.json').write_text(json.dumps(rec, indent=2))
    aliases = Path(a.aliases).resolve() if a.aliases else None
    if arm == 'D':
        assert aliases and aliases.exists(), 'arm D needs --aliases (aliases.json of the panel)'

    def pass1_spec(readset, donor, arm_):
        return dict(run=run_dir, ref=ref, r1=None, r2=None, threads=4,
                    aliases=aliases if arm_ == 'D' else None, receipt=dict(status='complete', donor=donor))
    RD.pass1_spec = pass1_spec

    # 2. frozen builder, all genes (same loop as the accepted wrapper)
    shutil.rmtree(out, ignore_errors=True)
    out.mkdir(parents=True, exist_ok=True)
    print(f'Building {a.donor} arm {arm}...', flush=True)
    shared = D2.load_shared(READSET, a.donor, out, arm)
    rec = dict(status='running', readset=READSET, donor=a.donor, arm=arm, started=time.time(), genes={},
               map_stats=shared['mstats'])
    for g in GENES:
        go = out / g
        go.mkdir(parents=True, exist_ok=True)
        t1 = time.time()
        grec = dict(gene=g, started=t1)
        cf = Path(a.confounded_dir) / f'confounded_{g}.json'
        res = R.run(READSET, a.donor, go, False, arm, gene=g, shared=shared,
                    confounded_path=cf if cf.exists() else None)
        grec.update(res, status='complete', seconds=round(time.time() - t1, 2))
        C.save_json(go / 'COMPLETE.json', grec)
        rec['genes'][g] = dict(status='complete', seconds=grec['seconds'])
    rec.update(status='complete', finished=time.time())
    C.save_json(out / 'COMPLETE.json', rec)
    print(f'Builder complete for {a.donor} arm {arm}.', flush=True)

    # 3. final calls of the requested arm: pass-1 table + the overrides (frozen ag_score.new_post)
    AS.task_dir = lambda tag, rs, d, arm_: out
    calls = C.read_t1k(run_dir / 'calls_genotype.tsv')
    final = {}
    for g in GENES:
        c = AS.new_post('dry', READSET, a.donor, a.arm, g, calls)
        final[g] = c.get(g, [[], []])
    with open(a.out_final_calls, 'w') as f:
        f.write('gene\tslot1\tslot2\n')
        for g in GENES:
            s = final[g]
            f.write(f"HLA-{g}\t{','.join(s[0]) or '.'}\t{','.join(s[1]) or '.'}\n")
    print(f'Final calls ({a.arm}) written to {a.out_final_calls}')


if __name__ == '__main__':
    main()
