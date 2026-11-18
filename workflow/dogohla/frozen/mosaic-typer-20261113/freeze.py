"""Write $S/FREEZE.md and $S/FREEZE.json: frozen parameters, code sha256, tool versions, endpoints.
  python freeze.py            (refuses to overwrite an existing FREEZE.json)
"""
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C
import proposer as PR

ENV = Path('/encrypted/e3008/hohndor-dogohla-20261020/mosaic-read-mapping-20261104/env/bin')


def tool(cmd):
    try:
        r = subprocess.run(cmd, capture_output=True, text=True)
        return (r.stdout + r.stderr).strip().splitlines()[0:3]
    except Exception as e:
        return [repr(e)]


def main():
    out = C.S / 'FREEZE.json'
    if out.exists():
        raise SystemExit('FREEZE.json exists; any change after the freeze is a recorded deviation')
    prm = C.load_json(C.S / 'src/PARAMS.json')
    PR.Params(**prm)  # validates the keys
    code = {p.name: C.sha256_file(p) for p in sorted((C.S / 'src').glob('*'))
            if p.is_file() and p.suffix in ('.py', '.sbatch', '.json', '.sh')}
    import numpy, edlib, pysam
    versions = dict(python=sys.version.split()[0], numpy=numpy.__version__, pysam=pysam.__version__,
                    edlib=getattr(edlib, '__version__', 'unknown'),
                    bwa=[l for l in tool([str(ENV / 'bwa')]) if 'Version' in l],
                    t1k=C.load_json(C.AUGMENTED / 'T1K_RUNTIME.json')['version'])
    freeze = dict(params=prm, code_sha256=code, inputs_sha256=C.verify_inputs(), versions=versions,
                  genomic_fa_sha256=C.GENOMIC_SHA)
    C.save_json(out, freeze)
    md = ['# FREEZE — mosaic-typer-20261113', '',
          'Frozen after the round-3 development gate (results/DEV_GATE_R3.md). Nothing below changes; any later',
          'change is a recorded deviation. FREEZE.json carries the same content machine-readably.', '',
          '## Parameters (src/PARAMS.json)', '', '```', json.dumps(prm, indent=1), '```', '',
          '## Code sha256', ''] + [f'- `{k}` {v}' for k, v in code.items()] + [
          '', '## Inputs sha256', ''] + [f'- `{k}` {v}' for k, v in freeze['inputs_sha256'].items()] + [
          '', '## Versions', '', '```', json.dumps(versions, indent=1), '```', '',
          '## Arms', '',
          '- A genomic (existing runs; reproduced byte-for-byte on HG00438), B genomic_core_hprc_asian (existing runs),',
          '  C proposer IPD-only + pass 2, D proposer IPD + arm-B fold cohort + pass 2, E oracle (IPD + own true cores).',
          '- Final call rule: IPD-prior (ties to the named allele), then read-evidence acceptance of MOS final calls',
          '  (src/accept.py; DESIGN_ISSUES 11 — deviation set on development data before any test result).', '',
          '## Endpoints (fixed; each reported separately for T-sim and T-real, never pooled; D-real reported',
          'separately as development donors with real reads)', '',
          '- PRIMARY: exact-sequence recovery of IPD-absent truth copies (truth core sha256 equals a novel',
          '  candidate present in a called slot of the gene), arm C minus arm A, over all IPD-absent copies of',
          '  the read set (private copies included). Also on the IPD-only-mosaic subset, final-call recovery,',
          '  and C as a fraction of E.',
          '- SAFETY 1: four-field named accuracy on IPD-present copies, C minus A; non-inferiority margin',
          '  -1.0 percentage point on the lower bound of the 95% CI. Every lost call listed.',
          '- SAFETY 2: false mosaics = copies whose final call is a MOS_ candidate not equal to the truth,',
          '  by truth class (exact / mosaic / near-mosaic / private): count and rate.',
          '- SECONDARY: C vs B, D vs B, 2/3-field (strict and decoded) per arm, per gene (DRB1 and DPA1',
          '  separately), proposal counts per donor, runtime.',
          '- Statistics: paired by copy, FOLDS-component-clustered bootstrap (4000, seed 20261020) 95% CI and',
          '  family sign-flip test (score.py:uncertainty = interpretation-20261020 logic).',
          '- Scorer: src/score.py (sha above) -> results/SCORES_test.json, score_rows_test.tsv; src/report.py.']
    (C.S / 'FREEZE.md').write_text('\n'.join(md) + '\n')
    print(json.dumps(freeze, indent=1))


if __name__ == '__main__':
    main()
