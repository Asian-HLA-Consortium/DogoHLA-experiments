#!/usr/bin/env python3
"""Part 7: build-cohort (resource folder for new-person typing and scoring).

The frozen builder and scorer read a "cohort resource folder": SOURCE_RECORDS.json, FOLDS.json,
mosaic_classes.tsv, panels_atlas.py, msa/<gene>_unique_out.fa and INPUTS.sha256.tsv. Our own
cohort's copy of that folder is protected and is not shipped. This step builds the same folder
from what a new user has:

  * the explicit panel (reference.fa + aliases.json; the arm-B reference),
  * optionally the user's own truth, i.e. the part-1 records (extract-genes) of the user's
    phased assemblies, used only by the scorer,
  * optionally a family table (donor<TAB>family); donors without a family are singletons.

It does not contain any typing or scoring logic. The multiple alignment (MAFFT --auto on the
unique sequences, same call as the atlas stage), the variant backbone and the class labels
("exact", "mosaic", "private", ...) come from the vendored atlas code in dogohla/atlas
(byte-identical copies, checked against dogohla/atlas/MANIFEST.sha256), called as the atlas stage
called them. Panel sequences enter as cohort records of the donor "PANEL", so the frozen leakage
check treats the whole panel as unrelated to the typed donor: the panel MUST NOT contain the
typed sample or its relatives (build it with --exclude-samples, or choose another panel).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sys
from pathlib import Path

GENES = ('A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1')
FLANK = 2000


def sha256_str(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


def sha256_file(p) -> str:
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for c in iter(lambda: f.read(1 << 22), b''):
            h.update(c)
    return h.hexdigest()


def verify_atlas(atlas: Path) -> int:
    n = 0
    for line in (atlas / 'MANIFEST.sha256').read_text().splitlines():
        if not line.strip():
            continue
        sha, name = line.split(None, 1)
        got = sha256_file(atlas / name.strip())
        if got != sha:
            raise SystemExit(f'vendored atlas file changed: {name.strip()}')
        n += 1
    return n


def read_fasta(path: Path) -> dict:
    out, name, parts = {}, None, []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    out[name] = ''.join(parts).upper()
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip())
    if name is not None:
        out[name] = ''.join(parts).upper()
    return out


def panel_records(panel_reference, panel_aliases):
    """Cohort records of the explicit panel: the PG sequences of reference.fa listed in aliases.json."""
    if not panel_reference or not panel_aliases:
        return []
    aliases = json.loads(Path(panel_aliases).read_text())
    ref = read_fasta(Path(panel_reference))
    recs = []
    for name, v in aliases.items():
        seq = ref[name]
        if sha256_str(seq) != v['sequence_sha256']:
            raise SystemExit(f'panel sequence differs from aliases.json: {name}')
        recs.append(dict(
            source_id=f'PANEL#{name}', donor='PANEL', haplotype='1', gene=v['gene'], resource='PANEL',
            sequence='N' * FLANK + seq + 'N' * FLANK,       # flanks of panel entries are not kept
            exons=[[s + FLANK, e + FLANK] for s, e in v['exons']],
            genomic_labels=[], core_sha256=sha256_str(seq), eligible_reference=True,
            core_boundary_policy='panel sequence; flanks not available'))
    return recs


def truth_records(files, ipd_by_gene, resource='USER'):
    """Truth records from part 1 (extract-genes JSON), in the layout of the cohort SOURCE_RECORDS."""
    names_of = {}
    for g, lst in ipd_by_gene.items():
        d = {}
        for x in lst:
            d.setdefault(x['sequence'].upper(), []).append(x['allele'].removeprefix('HLA-'))   # 'A*24:02:01:01'
        names_of[g] = d
    recs = []
    for fn in files:
        data = json.loads(Path(fn).read_text())
        data = data['records'] if isinstance(data, dict) else data
        for r in data:
            g, core = r['gene'], r['core_sequence'].upper()
            if g not in GENES:
                continue
            if sha256_str(core) != r['core_sha256']:
                raise SystemExit(f"core sha256 mismatch in {fn}: {r['donor']} {g}")
            seq = r['sequence'].upper()
            if len(seq) != len(core) + 2 * FLANK or seq[FLANK:-FLANK] != core:
                seq = 'N' * FLANK + core + 'N' * FLANK       # assembly ended inside the flank
            recs.append(dict(
                source_id=f"{r['donor']}#{r['haplotype']}#HLA-{g}", donor=r['donor'], haplotype=str(r['haplotype']),
                gene=g, resource=resource, sequence=seq,
                exons=[[s + FLANK, e + FLANK] for s, e in r['core_exons']],
                genomic_labels=names_of[g].get(core, []), core_id=f"{g}*PG{r['core_sha256']}",
                core_sha256=r['core_sha256'], eligible_reference=True,
                core_boundary_policy='gene interval plus exactly 2000 bp flanks'))
    return recs


def l1_labels(res):
    if res.k_switches_L1 == 0:
        return 'exact', 'exact'
    if res.k_switches_L1 is not None and res.k_switches_L1 >= 1:
        return 'mosaic', f'mosaic-{res.k_switches_L1}'
    if res.label == 'near-mosaic':
        return 'near-mosaic', 'near-mosaic'
    return 'private', 'private'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--donors', nargs='+', required=True, help='sample ids to be typed (and scored if truth is given)')
    ap.add_argument('--ipd-genomic', required=True)
    ap.add_argument('--panel-reference')
    ap.add_argument('--panel-aliases')
    ap.add_argument('--truth-records', nargs='*', default=[], help='part-1 record JSON files (extract-genes)')
    ap.add_argument('--truth-resource', default='USER', help='label written to the `cohort` column of the score rows for the truth samples')
    ap.add_argument('--families', help='optional table: donor<TAB>family id; same id = relatives')
    ap.add_argument('--mafft-prefix', help='folder with bin/mafft (default: mafft on PATH)')
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--package-dir', default=str(Path(__file__).resolve().parent.parent / 'dogohla'))
    a = ap.parse_args()

    atlas = Path(a.package_dir).resolve() / 'atlas'
    print('vendored atlas files verified:', verify_atlas(atlas), flush=True)
    sys.path.insert(0, str(atlas))
    import panels as P                      # noqa: E402
    import run_gene_pipeline as RG          # noqa: E402
    from mosaic_classifier import classify_query_against_panel   # noqa: E402

    mafft = Path(a.mafft_prefix) / 'bin' / 'mafft' if a.mafft_prefix else shutil.which('mafft')
    if a.mafft_prefix:       # conda-style prefix: the wrapper script finds its binaries through this variable
        os.environ['MAFFT_BINARIES'] = str(Path(a.mafft_prefix) / 'libexec' / 'mafft')
    minimap2 = shutil.which('minimap2')
    if not mafft or not Path(mafft).exists():
        raise SystemExit('mafft not found: use an image that contains it (see containers/Dockerfile) '
                         'or give mafft_dir (a folder with bin/mafft, MAFFT 7.526)')
    assert minimap2, 'minimap2 not found'

    out = Path(a.out_dir).resolve()
    shutil.rmtree(out, ignore_errors=True)
    (out / 'msa').mkdir(parents=True)
    work = Path('cohort_work').resolve()
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)

    ipd_by_gene = P.load_ipd_reference(Path(a.ipd_genomic))
    prec = panel_records(a.panel_reference, a.panel_aliases)
    trec = truth_records(a.truth_records, ipd_by_gene, a.truth_resource)
    truth_donors = sorted({r['donor'] for r in trec})
    truth_ids = {r['source_id'] for r in trec}
    donors = sorted(set(a.donors) | set(truth_donors))
    records = prec + trec
    print(f'panel records {len(prec)}, truth records {len(trec)} (donors {truth_donors}); typed donors {a.donors}', flush=True)

    # SOURCE_RECORDS.json
    (out / 'SOURCE_RECORDS.json').write_text(json.dumps(dict(
        status='complete', note='built by build-cohort from the explicit panel and the user truth records',
        records=records)) + '\n')

    # FOLDS.json: groups from the family table, the rest singletons
    groups = []
    if a.families:
        fam = {}
        for line in Path(a.families).read_text().splitlines():
            v = line.rstrip('\n').split('\t')
            if len(v) >= 2 and not line.startswith('#'):
                fam.setdefault(v[1], []).append(v[0])
        groups = [sorted(m) for m in fam.values() if len(m) > 1]
    reg = {d: dict(population='user', resource='USER') for d in donors}
    reg['PANEL'] = dict(population='panel', resource='PANEL')
    for g in groups:
        for d in g:
            reg.setdefault(d, dict(population='user', resource='USER'))
    (out / 'FOLDS.json').write_text(json.dumps(dict(status='complete', groups=groups, donors=reg,
                                                      rules='built by build-cohort')) + '\n')
    shutil.copy(atlas / 'panels.py', out / 'panels_atlas.py')

    # reload through the atlas loaders (checks the core slices), exclusion map
    all_recs = P.load_source_records(out / 'SOURCE_RECORDS.json')
    excl_map = P.build_exclusion_map_from_folds(out / 'FOLDS.json')
    rows_all = []
    for g in GENES:
        grecs = [r for r in all_recs if r['gene'] == g]
        ipd_recs = ipd_by_gene[g]
        msa_dir = work / 'msa'
        cohort_aligned, ipd_aligned, uniq_aligned, stats = RG.build_or_load_gene_msa(
            g, grecs, ipd_recs, msa_dir, Path(mafft), a.threads, True)
        for n in (f'{g}_unique_out.fa', f'{g}_msa_stats.json'):
            shutil.copy(msa_dir / n, out / 'msa' / n)
        queries = [r for r in grecs if r['source_id'] in truth_ids]
        if not queries:
            continue
        per_gene = work / 'per_gene'
        per_gene.mkdir(exist_ok=True)
        backbone, exon_cols, events = RG.run_part_a_for_gene(
            g, grecs, ipd_recs, cohort_aligned, ipd_aligned, uniq_aligned, stats['canonical_allele'],
            Path(minimap2), per_gene, work)
        ipd_dedup = RG._dedup_panel_dict(ipd_aligned)
        ipd_hashes = {x['sequence_sha256'] for x in ipd_recs}
        for r in queries:
            d, qid, q_aln = r['donor'], r['source_id'], cohort_aligned[r['source_id']]
            sel = P.select_cohort_panel_for_donor(d, grecs, excl_map)
            coh = RG._dedup_panel_dict({x['source_id']: cohort_aligned[x['source_id']] for x in sel})
            full = dict(ipd_dedup)
            seen = set(ipd_dedup.values())
            for k, v in coh.items():
                if v not in seen:
                    full[k] = v
                    seen.add(v)
            rf = classify_query_against_panel(g, qid, q_aln, full, backbone, exon_cols, 150, events)
            ri = classify_query_against_panel(g, qid, q_aln, ipd_dedup, backbone, exon_cols, 150, events)
            fl1, _ = l1_labels(rf)
            il1, _ = l1_labels(ri)
            coh_hashes = {x['core_sha256'] for x in sel}
            rows_all.append(dict(
                source_id=qid, gene=g, donor=d, haplotype=r.get('haplotype', ''), cohort=r['resource'],
                ipd_absent=len(r.get('genomic_labels', [])) == 0, core_len_bp=len(r['core_sequence']),
                full_label_L150=rf.label, full_detailed_L150=rf.detailed_label,
                full_k_L150='' if rf.k_switches is None else rf.k_switches,
                full_label_L1=fl1, full_k_L1='' if rf.k_switches_L1 is None else rf.k_switches_L1,
                full_exact_full_length=(r['core_sha256'] in ipd_hashes) or (r['core_sha256'] in coh_hashes),
                ipd_label_L150=ri.label, ipd_detailed_L150=ri.detailed_label,
                ipd_k_L150='' if ri.k_switches is None else ri.k_switches,
                ipd_label_L1=il1, ipd_k_L1='' if ri.k_switches_L1 is None else ri.k_switches_L1,
                ipd_exact_full_length=r['core_sha256'] in ipd_hashes))
        print(f'[{g}] classified {len(queries)} truth copies', flush=True)

    cols = ['source_id', 'gene', 'donor', 'haplotype', 'cohort', 'ipd_absent', 'core_len_bp', 'full_label_L150',
            'full_detailed_L150', 'full_k_L150', 'full_label_L1', 'full_k_L1', 'full_exact_full_length',
            'ipd_label_L150', 'ipd_detailed_L150', 'ipd_k_L150', 'ipd_label_L1', 'ipd_k_L1', 'ipd_exact_full_length']
    with open(out / 'mosaic_classes.tsv', 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows_all:
            f.write('\t'.join(str(r[c]) for c in cols) + '\n')

    # INPUTS.sha256.tsv: copy_sha, source_sha, name, source (copies are the sources here)
    names = ['SOURCE_RECORDS.json', 'FOLDS.json', 'mosaic_classes.tsv', 'panels_atlas.py'] + \
            [f'msa/{p.name}' for p in sorted((out / 'msa').iterdir())]
    with open(out / 'INPUTS.sha256.tsv', 'w') as f:
        for n in names:
            s = sha256_file(out / n)
            f.write(f'{s}\t{s}\t{n}\t{out / n}\n')
    (out / 'COHORT.json').write_text(json.dumps(dict(
        status='complete', typed_donors=a.donors, truth_donors=truth_donors, panel_records=len(prec),
        truth_records=len(trec), classified=len(rows_all), mafft=str(mafft), threads=a.threads,
        atlas_manifest_sha256=sha256_file(atlas / 'MANIFEST.sha256')), indent=2) + '\n')
    print(f'cohort folder written to {out}: {len(names)} input files, {len(rows_all)} classified truth copies')


if __name__ == '__main__':
    main()
