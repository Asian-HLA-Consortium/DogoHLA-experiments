#!/usr/bin/env python3
"""Part 6: score (thin wrapper around the frozen scorer).

Runs the frozen all-gene scorer (dogohla/frozen/mosaic-allgenes-20261114/ag_score.py on top of the frozen
mosaic-typer-20261113/score.py) for one donor, exactly as the accepted equivalence wrapper did, and writes the
rows it produces: one row per truth copy and arm, all columns of the frozen scorer, for the arms
  A, B      the T1K tables of the IPD run and the panel run,
  C4, C4R4n builder on the IPD run, without / with the R4n abstain rule,
  D4, D4R4n builder on the panel run, without / with the R4n abstain rule.
Arms B, D4 and D4R4n are scored only when the panel run is given.
No scoring logic is written here; this file only lays the inputs out in the folder structure the
frozen code reads (a copy of the builder folders, with the path of the pass-1 run in their per-gene
COMPLETE.json relocated to this job's folder).
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import shutil
import sys
from pathlib import Path

READSET = 'D-real3'


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for c in iter(lambda: f.read(1 << 22), b''):
            h.update(c)
    return h.hexdigest()


def stage_run(run_dir: Path, table: Path, donor: str, arm_letter: str, ref: Path | None):
    """pass-1 T1K run folder in the layout the frozen code reads."""
    run_dir.mkdir(parents=True, exist_ok=True)
    dst = run_dir / 'calls_genotype.tsv'
    if dst.is_symlink() or dst.exists():
        dst.unlink()
    dst.symlink_to(table.resolve())
    rec = {'status': 'complete', 'donor': donor, 'readset': READSET, 'arm': arm_letter, 'truth_loaded': False,
           'native_table_sha256': sha256_file(dst)}
    if ref is not None:
        rec['reference_sha256'] = sha256_file(ref)
    (run_dir / 'COMPLETE.json').write_text(json.dumps(rec, indent=2))


def stage_builder(src: Path, dst: Path, pass1_run: Path, pass1_ref: Path | None):
    """copy of a builder output folder; the per-gene COMPLETE.json point at this job's pass-1 run"""
    shutil.rmtree(dst, ignore_errors=True)
    shutil.copytree(src, dst)
    for gc in dst.glob('*/COMPLETE.json'):
        rec = json.loads(gc.read_text())
        rec['pass1_run'] = str(pass1_run)
        rec['pass1_reference'] = str(pass1_ref) if pass1_ref else rec.get('pass1_reference')
        gc.write_text(json.dumps(rec, indent=2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--donor', required=True)
    ap.add_argument('--calls-a', required=True, help='T1K genotype table, IPD run (arm A)')
    ap.add_argument('--calls-b', help='T1K genotype table, panel run (arm B)')
    ap.add_argument('--builder-dir-c', required=True, help='patch-builder output folder, IPD arm')
    ap.add_argument('--builder-dir-d', help='patch-builder output folder, panel arm')
    ap.add_argument('--panel-reference', help='T1K reference FASTA of the panel run (needed with arm B / D)')
    ap.add_argument('--ipd-genomic', required=True)
    ap.add_argument('--builder-inputs', required=True, help='cohort resource folder (SOURCE_RECORDS.json, msa/, ...)')
    ap.add_argument('--params', required=True, help='frozen builder parameter file PARAMS.json')
    ap.add_argument('--out-rows', required=True)
    ap.add_argument('--out-avail', required=True)
    ap.add_argument('--package-dir', default=str(Path(__file__).resolve().parent.parent / 'dogohla'))
    a = ap.parse_args()

    panel = bool(a.calls_b and a.builder_dir_d and a.panel_reference)
    work = Path('score_work').resolve()
    work.mkdir(parents=True, exist_ok=True)

    spec = importlib.util.spec_from_file_location('dogohla_wrapper_env', str(Path(a.package_dir) / 'wrapper_env.py'))
    WE = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(WE)
    C, RD, D2, R, AS = WE.setup(a.package_dir, a.builder_inputs, a.ipd_genomic, a.params, work)
    import common3 as K      # frozen: K.OLD maps the re-pulled read-set label to the old one

    d = a.donor
    TYPE_DIR = work / 'type'
    panel_ref = Path(a.panel_reference).resolve() if a.panel_reference else None
    stage_run(TYPE_DIR / d / 'A', Path(a.calls_a), d, 'A', None)
    stage_builder(Path(a.builder_dir_c), TYPE_DIR / d / 'C', TYPE_DIR / d / 'A', None)
    if panel:
        stage_run(TYPE_DIR / d / 'B', Path(a.calls_b), d, 'B', panel_ref)
        stage_builder(Path(a.builder_dir_d), TYPE_DIR / d / 'D', TYPE_DIR / d / 'B', panel_ref)

    # the truth donors of the frozen scorer are the donors with an existing run: only this donor
    C.readset_runs = lambda: {'D-sim': {d: TYPE_DIR / d}, 'T-sim': {}, 'real': {d: TYPE_DIR / d}}

    class Truth3(AS.SC.Truth):
        def __init__(self, rsl):
            base = sorted({K.OLD.get(r, r) for r in rsl})
            super().__init__(base)
            for (r0, dd, g), cps in list(self.copies.items()):
                for r in rsl:
                    if K.OLD.get(r) == r0:
                        self.copies[(r, dd, g)] = [dict(c, readset=r) for c in cps]
            for k in [k for k in self.copies if k[0] not in rsl]:
                del self.copies[k]

    def custom_arm_table(rs, dd, arm):
        run = TYPE_DIR / dd / arm
        if arm == 'A':
            return C.verified_table(run), []
        return C.verified_table(run), [panel_ref]

    def custom_task_dir(tag_, rs, dd, arm):
        return TYPE_DIR / dd / arm[0]

    AS.SC.Truth = Truth3
    AS.SC.arm_table = custom_arm_table
    AS.task_dir = custom_task_dir
    AS.TEST = ()
    AS.S2 = work / 'eq'
    (AS.S2 / 'results').mkdir(parents=True, exist_ok=True)

    arms = 'A,B,C4,C4R4n,D4,D4R4n' if panel else 'A,C4,C4R4n'
    sys.argv = ['ag_score.py', 'dry', READSET, arms]
    AS.main()

    shutil.copy(AS.S2 / 'results/score_rows_dry.tsv', a.out_rows)
    shutil.copy(AS.S2 / 'results/score_avail_dry.tsv', a.out_avail)
    print(f'Score rows written to {a.out_rows} (arms {arms})')


if __name__ == '__main__':
    main()
