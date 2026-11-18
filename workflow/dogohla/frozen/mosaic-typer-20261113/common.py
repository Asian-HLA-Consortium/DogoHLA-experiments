"""Shared paths, loaders and helpers for stage mosaic-typer-20261113.

Everything here is read-only with respect to other stages. Inputs copied into $S/inputs are
verified against INPUTS.sha256.tsv. Truth (SOURCE_RECORDS) is loaded only by Part 0, the
synthetic generator, leakage checks and the scorer -- never by the proposer.
"""
from __future__ import annotations

import collections
import hashlib
import importlib.util
import json
import os
import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

S = Path('/encrypted/e3008/hohndor-dogohla-20261020/mosaic-typer-20261113')
E = Path('/encrypted/e3008/hohndor-dogohla-20261020')
BENCH = Path('/ibex/scratch/projects/c2014/rob/dogohla-benchmark')
ZEN = BENCH / 'hprc-zenodo-20261004'
GENOMIC_FA = ZEN / 'baselines/assets/references/genomic.fa'
GENOMIC_SHA = '9202c3d37a89045a8f01f08dfd3ba21c465835b72ac2e9e7a731ef4dfc2c1d6a'
ATLAS = E / 'mosaic-variant-atlas-20261103'
REGIONAL = E / 'simulation-regional-20261022'
HPRC_CONTROL = E / 'simulation-hprc-control-20261025'
AUGMENTED = E / 'simulation-augmented-regional-20261024'
ARMS_REAL = ZEN / 'genomic-reference-arms-20261018'
GENES = ('A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1')
CLASS_I = {'A', 'B', 'C'}
FLANK = 2000
FORBIDDEN = 'confirmatory-test-20261113'


def guard_path(p) -> Path:
    p = Path(p)
    if FORBIDDEN in str(p):
        raise PermissionError('off-limits path: ' + str(p))
    return p


def sha256_file(p) -> str:
    h = hashlib.sha256()
    with open(guard_path(p), 'rb') as f:
        for c in iter(lambda: f.read(8 * 1024 ** 2), b''):
            h.update(c)
    return h.hexdigest()


def sha256_str(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


def load_json(p):
    return json.loads(guard_path(p).read_text())


def save_json(p, x):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    t = p.with_suffix(p.suffix + '.tmp')
    t.write_text(json.dumps(x, indent=2, default=str) + '\n')
    t.replace(p)


def verify_inputs() -> Dict[str, str]:
    """Check every copied input against its recorded sha256; return {name: sha}."""
    out = {}
    for line in (S / 'inputs/INPUTS.sha256.tsv').read_text().splitlines():
        copy_sha, src_sha, name, src = line.split('\t')
        assert copy_sha == src_sha, 'copy differs from source at copy time: ' + name
        got = sha256_file(S / 'inputs' / name)
        assert got == copy_sha, 'input changed since copy: ' + name
        out[name] = got
    return out


def import_by_path(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    import sys
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def panels_module():
    """The atlas stage's panels.py (copied verbatim, sha recorded in INPUTS.sha256.tsv)."""
    return import_by_path('panels_atlas', S / 'inputs/panels_atlas.py')


# ----------------------------------------------------------------------------- truth / folds

def load_records() -> List[Dict]:
    return panels_module().load_source_records(S / 'inputs/SOURCE_RECORDS.json')


def exclusion_map(truth_donors=None) -> Dict[str, Set[str]]:
    return panels_module().build_exclusion_map_from_folds(S / 'inputs/FOLDS.json', truth_donors)


def load_classes() -> Dict[str, Dict[str, str]]:
    rows = {}
    with open(S / 'inputs/mosaic_classes.tsv') as f:
        head = f.readline().rstrip('\n').split('\t')
        for line in f:
            v = line.rstrip('\n').split('\t')
            rows[v[0]] = dict(zip(head, v))
    return rows


# ----------------------------------------------------------------------------- IPD reference

def read_fasta_headers_and_seqs(path) -> List[Tuple[str, List[int], str]]:
    out = []
    name, ints, parts = None, [], []
    with open(guard_path(path)) as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    out.append((name, ints, ''.join(parts).upper()))
                tok = line[1:].split()
                name, ints, parts = tok[0], [int(x) for x in tok[1:]], []
            else:
                parts.append(line.strip())
    if name is not None:
        out.append((name, ints, ''.join(parts).upper()))
    return out


def header_exons(ints: Sequence[int]) -> List[Tuple[int, int]]:
    """T1K header ints -> [(start, end_exclusive)]."""
    if not ints:
        return []
    n = ints[0]
    return [(ints[1 + 2 * i], ints[2 + 2 * i] + 1) for i in range(n)]


def t1k_header(name: str, exons: Sequence[Tuple[int, int]], length: int) -> str:
    """Inverse of header_exons, with the same validity checks as t1k_reference.header."""
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


def gene_of_name(name: str) -> str:
    return name.split('*')[0].removeprefix('HLA-')


# ----------------------------------------------------------------------------- MSA

class GeneMSA:
    """Per-gene MSA from the atlas stage, with the U-id -> sequence mapping rebuilt exactly as
    run_gene_pipeline.build_or_load_gene_msa registered it (canonical IPD allele first, then IPD
    records in genomic.fa order, then cohort cores in SOURCE_RECORDS order)."""

    def __init__(self, gene: str, ipd_records: List[Dict], cohort_records: List[Dict]):
        self.gene = gene
        canon = f'HLA-{gene}*01:01:01:01'
        by_name = {x['allele']: x for x in ipd_records}
        if canon not in by_name:
            canon = ipd_records[0]['allele']
        seq_to_uid: Dict[str, str] = {}

        def reg(s):
            s = s.upper()
            if s not in seq_to_uid:
                seq_to_uid[s] = f'U{len(seq_to_uid):05d}'
            return seq_to_uid[s]

        reg(by_name[canon]['sequence'])
        for x in ipd_records:
            reg(x['sequence'])
        for r in cohort_records:
            reg(r['core_sequence'])
        aligned = {}
        name = None
        parts = []
        with open(S / f'inputs/msa/{gene}_unique_out.fa') as f:
            for line in f:
                if line.startswith('>'):
                    if name:
                        aligned[name] = ''.join(parts).upper()
                    name, parts = line[1:].strip().split()[0], []
                else:
                    parts.append(line.strip())
        if name:
            aligned[name] = ''.join(parts).upper()
        assert len(aligned) == len(seq_to_uid), f'{gene}: MSA row count mismatch'
        self.uids = sorted(seq_to_uid.values())
        self.row_of_uid = {u: i for i, u in enumerate(self.uids)}
        self.aligned = [aligned[u] for u in self.uids]
        self.seqs = [None] * len(self.uids)
        for s, u in seq_to_uid.items():
            i = self.row_of_uid[u]
            assert self.aligned[i].replace('-', '') == s, f'{gene}: MSA row does not spell {u}'
            self.seqs[i] = s
        self.ncol = len(self.aligned[0])
        self.row_of_seq = {s: i for i, s in enumerate(self.seqs)}
        self.ipd_names: List[List[str]] = [[] for _ in self.uids]
        self.ipd_exons: List[Optional[List[Tuple[int, int]]]] = [None] * len(self.uids)
        for x in ipd_records:
            i = self.row_of_seq[x['sequence'].upper()]
            self.ipd_names[i].append(x['allele'])
            if self.ipd_exons[i] is None:
                self.ipd_exons[i] = list(x['exons'])
        self.row_of_allele = {x['allele']: self.row_of_seq[x['sequence'].upper()] for x in ipd_records}
        self.cohort_sources: List[List[str]] = [[] for _ in self.uids]  # source_ids
        self.cohort_donors: List[Set[str]] = [set() for _ in self.uids]
        self.cohort_exons: List[Optional[List[Tuple[int, int]]]] = [None] * len(self.uids)
        for r in cohort_records:
            i = self.row_of_seq[r['core_sequence'].upper()]
            self.cohort_sources[i].append(r['source_id'])
            self.cohort_donors[i].add(r['donor'])
            if self.cohort_exons[i] is None:
                self.cohort_exons[i] = [tuple(e) for e in r['core_exons']]
        self.donor_of_source = {r['source_id']: r['donor'] for r in cohort_records}
        self.row_of_source = {r['source_id']: self.row_of_seq[r['core_sequence'].upper()]
                              for r in cohort_records}

    def exons(self, i) -> Optional[List[Tuple[int, int]]]:
        return self.ipd_exons[i] if self.ipd_exons[i] is not None else self.cohort_exons[i]

    def base_cols(self, i) -> List[int]:
        return [c for c, ch in enumerate(self.aligned[i]) if ch != '-']


_MSA_CACHE: Dict[str, GeneMSA] = {}


def load_msas(records: List[Dict], genes=GENES) -> Dict[str, GeneMSA]:
    pm = panels_module()
    assert sha256_file(GENOMIC_FA) == GENOMIC_SHA, 'genomic.fa pin mismatch'
    ipd = pm.load_ipd_reference(GENOMIC_FA)
    out = {}
    for g in genes:
        if g not in _MSA_CACHE:
            _MSA_CACHE[g] = GeneMSA(g, ipd[g], [r for r in records if r['gene'] == g])
        out[g] = _MSA_CACHE[g]
    return out


# ----------------------------------------------------------------------------- CDS / names

CODON = {}
_b = 'TCAG'
_aa = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
for _i, _a in enumerate(_aa):
    CODON[_b[_i // 16] + _b[(_i // 4) % 4] + _b[_i % 4]] = _a


def cds_of(seq: str, exons: Sequence[Tuple[int, int]]) -> str:
    return ''.join(seq[s:e] for s, e in exons)


def translate(cds: str) -> str:
    out = []
    for i in range(0, len(cds) - 2, 3):
        a = CODON.get(cds[i:i + 3], 'X')
        if a == '*':
            break
        out.append(a)
    return ''.join(out)


NAME_RE = re.compile(r'([A-Z0-9]+)\*(\d+(?::\d+)+)([NLSCAQ]?)')


def normalized(label: str, n: int) -> Optional[str]:
    """Same projection as simulation-regional-20261022/score_simulation.normalized."""
    label = label.removeprefix('HLA-')
    m = re.fullmatch(r'([A-Z0-9]+)\*(\d+(?::\d+)+)[NLSCAQ]?', label)
    return m[1] + '*' + ':'.join(m[2].split(':')[:n]) if m and len(m[2].split(':')) >= n else None


# ----------------------------------------------------------------------------- T1K tables

def split_labels(x: str) -> List[str]:
    return [] if x in ('.', '') else x.split(',')


def read_t1k(path) -> Dict[str, List[List[str]]]:
    """Same slot logic as score_simulation.read_t1k (copies==1 duplicates the slot)."""
    calls = {}
    for line in guard_path(path).read_text().splitlines():
        if not line or line.startswith('#'):
            continue
        v = line.split('\t')
        copies = int(v[1])
        a = v[2] if copies and float(v[4]) > 0 else '.'
        b = a if copies == 1 else v[5] if copies == 2 and float(v[7]) > 0 else '.'
        calls[v[0].removeprefix('HLA-')] = [split_labels(a), split_labels(b)]
    return calls


def read_t1k_raw(path) -> Dict[str, Dict]:
    """Full rows including abundance and quality per slot."""
    out = {}
    for line in guard_path(path).read_text().splitlines():
        if not line or line.startswith('#'):
            continue
        v = line.split('\t')
        out[v[0].removeprefix('HLA-')] = dict(copies=int(v[1]), slots=[
            dict(labels=split_labels(v[2]), abundance=float(v[3]), quality=int(v[4])),
            dict(labels=split_labels(v[5]), abundance=float(v[6]), quality=int(v[7]))])
    return out


# ----------------------------------------------------------------------------- read sets

def readset_runs() -> Dict[str, Dict[str, Path]]:
    """Existing arm-A (genomic) T1K run folders per read set and donor."""
    out = {'D-sim': {}, 'T-sim': {}, 'real': {}}
    for p in sorted((HPRC_CONTROL / 't1k/runs').iterdir()):
        out['D-sim'][p.name] = p / 'genomic'
    for p in sorted((REGIONAL / 't1k/runs').iterdir()):
        out['T-sim'][p.name] = p / 'genomic'
    for p in sorted((ARMS_REAL / 'runs').iterdir()):
        out['real'][p.name] = p / 'genomic'
    for k in out:
        for d, p in out[k].items():
            guard_path(p)
    return out


def verified_table(run: Path) -> Path:
    r = load_json(run / 'COMPLETE.json')
    assert r['status'] == 'complete' and not r.get('truth_loaded'), run
    t = run / 'calls_genotype.tsv'
    assert sha256_file(t) == r['native_table_sha256'], 'genotype table changed: ' + str(run)
    return t


SYN_SETS = ('SYN', 'SYN2')


def syn_dir(cid: str) -> Path:
    """Work folder of a synthetic development case (set derived from the case-id prefix)."""
    return S / 'work' / ('SYN2' if cid.startswith('SYN2-') else 'SYN') / cid
