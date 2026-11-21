"""Set-up shared by the builder and scorer wrappers (scripts/patch_builder.py, scripts/score.py).

It does NOT contain any typing logic. It only
  1. verifies the vendored frozen code against MANIFEST.sha256 (verify.py) and puts it on sys.path,
  2. imports the frozen modules in the same order as the accepted equivalence runs,
  3. points the frozen modules' location constants at the folders given on the command line
     (the cohort resource folder, the IPD genomic.fa, the frozen parameter file, the bwa of the container).
The frozen files themselves are never edited.
"""
from __future__ import annotations

import importlib.util
import shutil
import sys
from pathlib import Path


def setup(package_dir, builder_inputs, ipd_genomic, params, workdir):
    package_dir = Path(package_dir).resolve()
    spec = importlib.util.spec_from_file_location('dogohla_verify', str(package_dir / 'verify.py'))
    verify = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verify)
    print('frozen code verified:', verify.verify(), flush=True)
    verify.use()

    # same import order as the accepted wrapper: common, run_donor first, so that the later
    # frozen modules, which prepend their original stage folder to sys.path, find these in sys.modules
    import common as C
    import run_donor as RD
    import part0, score, panel_select, evidence, proposer   # noqa: F401  (cached before drb1_lib prepends the original stage folder)
    import drb1_diag2 as D2
    import drb1_run as R
    import ag_score as AS

    frozen = (package_dir / 'frozen').resolve()
    for m in (C, RD, D2, R, AS, AS.SC, AS.SC.P0):
        f = Path(m.__file__).resolve()
        assert frozen in f.parents, f'{m.__name__} was not loaded from the vendored folder: {f}'

    workdir = Path(workdir).resolve()
    cohort = workdir / 'cohort'
    cohort.mkdir(parents=True, exist_ok=True)
    link = cohort / 'inputs'
    if link.is_symlink() or link.exists():
        link.unlink()
    link.symlink_to(Path(builder_inputs).resolve())
    C.S = cohort                                   # frozen code reads S/'inputs/...'
    C.GENOMIC_FA = Path(ipd_genomic).resolve()     # sha256-pinned inside the frozen code (GENOMIC_SHA)
    RD.FROZEN_PARAMS = Path(params).resolve()
    bwa = shutil.which('bwa')
    assert bwa, 'bwa not found on PATH'
    RD.ENV = Path(bwa).parent                      # the builder maps reads with RD.ENV/bwa
    return C, RD, D2, R, AS
