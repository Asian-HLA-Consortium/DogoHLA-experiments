"""Per-donor worker for stage mosaic-typer-20261113.

  python run_donor.py --manifest M --index I

Each manifest entry: {readset, donor, arms:[...]} with arms among
  A_repro  rerun the existing arm-A T1K command exactly (determinism check)
  A        (synthetic only) pass-1 T1K on IPD genomic.fa
  C        proposer with panel = IPD only on top of arm-A pass 1, then pass-2 T1K
  D        proposer with panel = IPD + the arm-B fold cohort content on top of arm-B pass 1, pass 2
  E        oracle: IPD + the donor's own true core sequences, single T1K pass
Outputs under $S/work/<readset>/<donor>/<arm>/ ; every arm writes COMPLETE.json or FAILURE.json.
Sequences (additions.fa, pass-2 reference) stay on IBEX; proposals.json carries no sequence.
"""
from __future__ import annotations

import argparse
import collections
import json
import os
import shutil
import subprocess
import sys
import time
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C
import evidence as EV
import panel_select as PS
import proposer as PR

SIF = C.BENCH / 'containers/ubuntu24.sif'
BIND = {'/benchmark': C.ZEN, '/baseline': C.REGIONAL, '/control': C.HPRC_CONTROL,
        '/augmented': C.AUGMENTED}
RUN_T1K = '/benchmark/baselines/assets/t1k/run-t1k'
ENV = Path('/encrypted/e3008/hohndor-dogohla-20261020/mosaic-read-mapping-20261104/env/bin')
FROZEN_PARAMS = C.S / 'src/PARAMS.json'


def to_host(p: str) -> Path:
    for k, v in BIND.items():
        if p == k or p.startswith(k + '/'):
            return C.guard_path(Path(str(v) + p[len(k):]))
    return C.guard_path(Path(p))


def to_container(p: Path) -> str:
    s = str(p)
    for k, v in BIND.items():
        if s.startswith(str(v) + '/'):
            return k + s[len(str(v)):]
    return s


def parse_cmd(cmd):
    out = {}
    for i, t in enumerate(cmd):
        if t in ('-1', '-2', '-f', '-t', '-o'):
            out[t] = cmd[i + 1]
    return out


def verify_runtime():
    rt = C.load_json(C.AUGMENTED / 'T1K_RUNTIME.json')
    for a in rt['artifacts']:
        got = C.sha256_file(to_host(a['path']))
        assert got == a['sha256'], 'T1K runtime artifact changed: ' + a['path']
    return rt['version']


def t1k(r1: str, r2: str, ref: Path, out: Path, threads: int, log: Path):
    binds = ','.join(f'{v}:{k}' for k, v in BIND.items()) + f',{C.S}:{C.S}'
    cmd = ['singularity', 'exec', '--cleanenv', '--no-home', '--bind', binds, str(SIF), RUN_T1K,
           '-1', r1, '-2', r2, '--preset', 'hla-wgs', '-f', to_container(ref), '-t', str(threads),
           '--alleleDigitUnits', '4', '--alleleDelimiter', ':', '-o', to_container(out / 'calls')]
    with open(log, 'w') as f:
        subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, check=True)
    table = out / 'calls_genotype.tsv'
    assert table.exists() and table.stat().st_size > 0, 'empty genotype table'
    return cmd


# ----------------------------------------------------------------------------- pass-1 specs

def pass1_spec(readset, donor, arm):
    """Existing pass-1 run for arm C (= arm A run) or arm D (= arm B run)."""
    if readset in C.SYN_SETS:
        d = C.syn_dir(donor) / 'A'
        rec = C.load_json(d / 'COMPLETE.json')
        return dict(run=d, ref=C.GENOMIC_FA, r1=rec['r1'], r2=rec['r2'], threads=4, aliases=None)
    if arm in ('A', 'C', 'E', 'A_repro'):
        base = {'D-sim': C.HPRC_CONTROL / 't1k/runs', 'T-sim': C.REGIONAL / 't1k/runs',
                'T-real': C.ARMS_REAL / 'runs', 'D-real': C.ARMS_REAL / 'runs'}[readset]
        run = base / donor / 'genomic'
    elif arm == 'D':
        base = {'T-sim': C.AUGMENTED / 'runs', 'T-real': C.ARMS_REAL / 'runs',
                'D-real': C.ARMS_REAL / 'runs'}[readset]
        run = base / donor / 'genomic_core_hprc_asian'
    else:
        raise ValueError(arm)
    rec = C.load_json(run / 'COMPLETE.json')
    assert rec['status'] == 'complete' and rec['donor'] == donor and not rec.get('truth_loaded')
    cmd = parse_cmd(rec['command'])
    ref = to_host(cmd['-f'])
    aliases = None
    if arm == 'D':
        aliases = ref.parent / 'aliases.json'
        assert C.sha256_file(ref) == rec['reference_sha256'], 'arm-B reference changed'
    return dict(run=run, ref=ref, r1=cmd['-1'], r2=cmd['-2'], threads=int(cmd['-t']),
                aliases=aliases, receipt=rec)


# ----------------------------------------------------------------------------- reference I/O

def read_reference(ref: Path, wanted_names, genes):
    """Return (sequences of wanted names, {gene: set(sequences)}) from a T1K reference FASTA."""
    seqs, by_gene = {}, collections.defaultdict(set)
    for name, ints, seq in C.read_fasta_headers_and_seqs(ref):
        g = C.gene_of_name(name)
        if name in wanted_names:
            seqs[name] = seq
        if g in genes:
            by_gene[g].add(seq)
    return seqs, by_gene


def map_reads(target: dict, fq1: Path, fq2: Path, work: Path, threads: int):
    import pysam
    work.mkdir(parents=True, exist_ok=True)
    fa = work / 'target.fa'
    with open(fa, 'w') as f:
        for n, s in target.items():
            f.write(f'>{n}\n{s}\n')
    subprocess.run([str(ENV / 'bwa'), 'index', str(fa)], check=True, capture_output=True)
    sam = work / 'aln.sam'
    with open(sam, 'w') as f:
        subprocess.run([str(ENV / 'bwa'), 'mem', '-t', str(threads), '-v', '1', str(fa), str(fq1),
                        str(fq2)], stdout=f, stderr=subprocess.DEVNULL, check=True)
    reads = collections.defaultdict(list)
    stats = collections.Counter()
    with pysam.AlignmentFile(str(sam), 'r') as af:
        lens = dict(zip(af.references, af.lengths))
        for rec in af.fetch(until_eof=True):
            if rec.is_secondary or rec.is_supplementary:
                continue
            stats['reads'] += 1
            if rec.is_unmapped:
                stats['unmapped'] += 1
                continue
            g = C.gene_of_name(rec.reference_name)
            stats['gene:' + g] += 1
            if g not in C.GENES:
                continue
            cig = rec.cigartuples or []
            lclip = cig[0][1] if cig and cig[0][0] == 4 else 0
            rclip = cig[-1][1] if cig and cig[-1][0] == 4 else 0
            reads[g].append(dict(qname=rec.query_name, seq=rec.query_sequence,
                                 reverse=rec.is_reverse, contig=rec.reference_name,
                                 ref_start=rec.reference_start, ref_end=rec.reference_end,
                                 lclip=lclip, rclip=rclip, contig_len=lens[rec.reference_name]))
    for p in (sam,):
        p.unlink()
    return reads, dict(stats)


# ----------------------------------------------------------------------------- proposer arm

def load_params():
    if FROZEN_PARAMS.exists():
        return PR.Params(**C.load_json(FROZEN_PARAMS)), C.sha256_file(FROZEN_PARAMS)
    return PR.Params(), None


def excluded_for(readset, donor, excl_map):
    if readset in C.SYN_SETS:
        case = C.load_json(C.syn_dir(donor) / 'CASE.json')
        out = set()
        for d in case['component_donors']:
            out |= excl_map[d]
        return out | {donor}
    return excl_map[donor]


def proposer_arm(readset, donor, arm, out: Path, ctx):
    proposal_only = arm == 'P'      # development: arm-C proposer without the pass-2 T1K run
    if proposal_only:
        arm = 'C'
    spec = pass1_spec(readset, donor, 'C' if arm == 'C' else 'D')
    table = C.verified_table(spec['run']) if readset not in C.SYN_SETS else spec['run'] / 'calls_genotype.tsv'
    geno = C.read_t1k_raw(table)
    prm, prm_sha = load_params()
    aliases = C.load_json(spec['aliases']) if spec['aliases'] else {}
    wanted = set()
    slots_by_gene = {}
    for g, row in geno.items():
        labs = []
        for sl in row['slots'][:row['copies']] if row['copies'] else []:
            if sl['abundance'] <= 0:
                continue
            labs.append(sl['labels'][:3])
            wanted |= set(sl['labels'][:3])
        slots_by_gene[g] = labs
    seqs, ref_by_gene = read_reference(spec['ref'], wanted, set(C.GENES))
    target = {n: seqs[n] for n in sorted(wanted) if n in seqs}
    reads, mstats = map_reads(target, spec['run'] / 'calls_candidate_1.fq',
                              spec['run'] / 'calls_candidate_2.fq', out / 'map', 4)
    msas = ctx['msas']
    excluded = excluded_for(readset, donor, ctx['excl'])
    pg_shas = {v['sequence_sha256'] for v in aliases.values()} if arm == 'D' else None
    gene_recs, additions, mos_aliases = {}, [], {}
    for g in C.GENES:
        msa = msas[g]
        shas = [C.sha256_str(s) for s in msa.seqs]
        row_is_ipd = [bool(n) for n in msa.ipd_names]
        rows = PS.select_rows(arm, row_is_ipd, shas, pg_shas)
        comp_absent = {msa.row_of_source[r['source_id']] for r in ctx['records']
                       if r['gene'] == g and r['donor'] in excluded and not r['genomic_labels']}
        PS.check_no_leakage(rows, row_is_ipd, msa.cohort_donors, excluded, donor, comp_absent, arm)
        if readset in C.SYN_SETS:
            case = C.load_json(C.syn_dir(donor) / 'CASE.json')
            if case['gene'] == g:
                novel = {v for v in case['truth_sha'].values() if v not in case['ipd_present_shas']}
                if novel & {shas[r] for r in rows}:
                    raise PS.LeakageError('synthetic truth sequence present in the panel')
        sub_names = [(msa.ipd_names[r][0] if msa.ipd_names[r] else 'COH_' + shas[r][:12]) for r in rows]
        panel = PR.Panel(sub_names, [msa.aligned[r] for r in rows])
        local = {r: i for i, r in enumerate(rows)}
        row_by_sha = {shas[r]: local[r] for r in rows}

        def label_row(lab):
            if lab in msa.row_of_allele and msa.row_of_allele[lab] in local:
                return local[msa.row_of_allele[lab]]
            if lab in aliases and aliases[lab]['sequence_sha256'] in row_by_sha:
                return row_by_sha[aliases[lab]['sequence_sha256']]
            return None
        pass1_rows = []
        for labs in slots_by_gene.get(g, []):
            rr = [label_row(l) for l in labs]
            pass1_rows.append([r for r in rr if r is not None])
        # read placement -> windows
        rd = []
        unplaced = 0
        for x in reads.get(g, []):
            r = label_row(x['contig'])
            if r is None:
                unplaced += 1
                continue
            bcols = EV_bcols(panel, r)
            c0, c1, tl, tr = EV.window_from_alignment(bcols, x['ref_start'], x['ref_end'],
                                                     x['lclip'], x['rclip'], x['contig_len'])
            seq = x['seq']
            if tr:
                seq = seq[:len(seq) - tr]
            if tl:
                seq = seq[tl:]
            if len(seq) < 30:
                unplaced += 1
                continue
            rd.append(dict(qname=x['qname'], seq=seq, reverse=x['reverse'], c0=c0, c1=c1,
                           prefer_row=r))
        t0 = time.time()
        ev = EV.build_evidence(panel, rd, prm.max_read_edit)
        row_exons = {i: msa.exons(r) for i, r in enumerate(rows)}
        rec = PR.propose(panel, ev, pass1_rows, prm, ref_by_gene.get(g, set()), row_exons, g)
        rec['evidence_stats'] = dict(ev.stats, reads_assigned=len(reads.get(g, [])),
                                     reads_unplaced=unplaced)
        rec['runtime_s'] = round(time.time() - t0, 1)
        n_used = 0
        for c in rec['candidates']:
            seq = c.pop('_sequence', None)
            if seq is None or n_used >= prm.n_max:
                continue
            name = f'HLA-{g}*MOS_{c["sha256"][:12]}'
            c['name'] = name
            additions.append(C.t1k_header(name, c['exons'], len(seq)) + seq + '\n')
            mos_aliases[name] = dict(sequence_sha256=c['sha256'], gene=g, sources=c['sources'],
                                     switch_columns=c['switch_columns'], internal_identifier=True,
                                     genomic_labels=[])
            n_used += 1
        gene_recs[g] = rec
    out.mkdir(parents=True, exist_ok=True)
    C.save_json(out / 'proposals.json', dict(donor=donor, readset=readset, arm=arm,
                params_sha256=prm_sha, mapping=mstats, genes=gene_recs))
    (out / 'aliases_mos.json').write_text(json.dumps(mos_aliases, indent=2) + '\n')
    (out / 'additions.fa').write_text(''.join(additions))
    shutil.rmtree(out / 'map', ignore_errors=True)
    res = dict(n_candidates=len(mos_aliases), params_sha256=prm_sha,
               pass1_run=str(spec['run']), pass1_table_sha256=C.sha256_file(table),
               pass1_reference=str(spec['ref']))
    if proposal_only:
        res.update(pass2='not_run_proposal_only')
        return res
    if not additions:
        res.update(pass2='not_run_no_candidates_reference_identical_to_pass1',
                   calls_table=str(table), calls_table_sha256=C.sha256_file(table))
        return res
    ref2 = out / 'pass2' / 'reference.fa'
    ref2.parent.mkdir(parents=True, exist_ok=True)
    with open(ref2, 'wb') as f:
        with open(spec['ref'], 'rb') as g0:
            shutil.copyfileobj(g0, f)
        f.write(''.join(additions).encode())
    ref_sha = C.sha256_file(ref2)
    cmd = t1k(spec['r1'], spec['r2'], ref2, out / 'pass2', spec['threads'], out / 'pass2/run.log')
    ref2.unlink()
    table2 = out / 'pass2/calls_genotype.tsv'
    res.update(pass2='complete', pass2_reference_sha256=ref_sha, command=cmd,
               calls_table=str(table2), calls_table_sha256=C.sha256_file(table2))
    import accept
    acc = accept.run(spec['run'], spec['ref'], table2, out / 'additions.fa', out, prm)
    res.update(accept_json_sha256=C.sha256_file(out / 'ACCEPT.json'),
               mos_final_calls=sum(len(v) for v in acc.values()),
               mos_rejected=sum(1 for v in acc.values() for x in v if not x['accepted']))
    return res


def EV_bcols(panel, r):
    import numpy as np
    return np.nonzero(panel.isbase[r])[0]


# ----------------------------------------------------------------------------- other arms

def oracle_arm(readset, donor, out: Path, ctx):
    spec = pass1_spec(readset, donor, 'E')
    ipd = set()
    for name, ints, seq in C.read_fasta_headers_and_seqs(C.GENOMIC_FA):
        if C.gene_of_name(name) in C.GENES:
            ipd.add(seq)
    adds, names = [], {}
    if readset in C.SYN_SETS:
        case = C.load_json(C.syn_dir(donor) / 'CASE.json')
        truth = C.load_json(C.syn_dir(donor) / 'TRUTH_SEQ.json')
        cores = [(case['gene'], truth[h]['core'], truth[h]['exons']) for h in sorted(truth)]
    else:
        cores = [(r['gene'], r['core_sequence'], r['core_exons']) for r in ctx['records']
                 if r['donor'] == donor]
    for g, core, ex in cores:
        if core in ipd:
            continue
        sha = C.sha256_str(core)
        name = f'HLA-{g}*ORC_{sha[:12]}'
        if name in names:
            continue
        names[name] = dict(sequence_sha256=sha, gene=g, internal_identifier=True, genomic_labels=[])
        adds.append(C.t1k_header(name, [tuple(e) for e in ex], len(core)) + core + '\n')
    out.mkdir(parents=True, exist_ok=True)
    (out / 'aliases_orc.json').write_text(json.dumps(names, indent=2) + '\n')
    (out / 'additions.fa').write_text(''.join(adds))
    if not adds:
        table = spec['run'] / 'calls_genotype.tsv'
        return dict(pass2='not_run_no_absent_cores', calls_table=str(table),
                    calls_table_sha256=C.sha256_file(table), n_added=0)
    ref = out / 'run' / 'reference.fa'
    ref.parent.mkdir(parents=True, exist_ok=True)
    with open(ref, 'wb') as f:
        with open(C.GENOMIC_FA, 'rb') as g0:
            shutil.copyfileobj(g0, f)
        f.write(''.join(adds).encode())
    ref_sha = C.sha256_file(ref)
    cmd = t1k(spec['r1'], spec['r2'], ref, out / 'run', spec['threads'], out / 'run/run.log')
    ref.unlink()
    table = out / 'run/calls_genotype.tsv'
    return dict(pass2='complete', reference_sha256=ref_sha, command=cmd, n_added=len(adds),
                calls_table=str(table), calls_table_sha256=C.sha256_file(table))


def repro_arm(readset, donor, out: Path):
    spec = pass1_spec(readset, donor, 'A_repro')
    out.mkdir(parents=True, exist_ok=True)
    cmd = t1k(spec['r1'], spec['r2'], spec['ref'], out, spec['threads'], out / 'run.log')
    got = C.sha256_file(out / 'calls_genotype.tsv')
    want = spec['receipt']['native_table_sha256']
    return dict(command=cmd, table_sha256=got, original_table_sha256=want, identical=got == want)


def synthetic_pass1(donor, out: Path):
    sys.path.insert(0, str(C.S / 'src'))
    import synth
    r1, r2 = synth.ensure_reads(donor)
    out.mkdir(parents=True, exist_ok=True)
    cmd = t1k(str(r1), str(r2), C.GENOMIC_FA, out, 4, out / 'run.log')
    return dict(command=cmd, r1=str(r1), r2=str(r2),
                calls_table_sha256=C.sha256_file(out / 'calls_genotype.tsv'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--manifest', required=True)
    ap.add_argument('--index', type=int, default=None)
    a = ap.parse_args()
    assert os.environ.get('SLURM_JOB_ID'), 'run under Slurm'
    idx = a.index if a.index is not None else int(os.environ['SLURM_ARRAY_TASK_ID'])
    task = C.load_json(a.manifest)['tasks'][idx]
    readset, donor = task['readset'], task['donor']
    version = verify_runtime()
    ctx = {}
    failed = []
    for arm in task['arms']:
        out = C.S / 'work' / readset / donor / arm
        if (out / 'COMPLETE.json').exists():
            continue
        if out.exists():
            shutil.rmtree(out)
        out.mkdir(parents=True)
        rec = dict(status='running', readset=readset, donor=donor, arm=arm,
                   job_id=os.environ.get('SLURM_JOB_ID'), task=idx, t1k_version=version,
                   code_sha256={p.name: C.sha256_file(p) for p in sorted((C.S / 'src').glob('*.py'))},
                   started=time.time())
        C.save_json(out / 'RUN.json', rec)
        try:
            if arm in ('C', 'D', 'E', 'P') and 'records' not in ctx:
                ctx['inputs'] = C.verify_inputs()
                ctx['records'] = C.load_records()
                ctx['msas'] = C.load_msas(ctx['records'])
                ctx['excl'] = C.exclusion_map()
            if arm in ('C', 'D', 'P'):
                res = proposer_arm(readset, donor, arm, out, ctx)
            elif arm == 'E':
                res = oracle_arm(readset, donor, out, ctx)
            elif arm == 'A_repro':
                res = repro_arm(readset, donor, out)
            elif arm == 'A':
                res = synthetic_pass1(donor, out)
            else:
                raise ValueError(arm)
            rec.update(res, status='complete', finished=time.time())
            C.save_json(out / 'COMPLETE.json', rec)
        except BaseException as e:
            rec.update(status='failed', finished=time.time(), error=repr(e),
                       traceback=traceback.format_exc())
            C.save_json(out / 'FAILURE.json', rec)
            failed.append(arm)
            print(traceback.format_exc(), file=sys.stderr)
            if arm == 'A':
                break
    if failed:
        raise SystemExit('failed arms: ' + ','.join(failed))


if __name__ == '__main__':
    main()
