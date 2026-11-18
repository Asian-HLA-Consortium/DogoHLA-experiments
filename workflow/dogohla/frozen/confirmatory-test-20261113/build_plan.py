"""Build PLAN.json from pinned metadata in ../inputs (run locally, before any test-donor job)."""
import json,csv
from pathlib import Path
L=Path(__file__).resolve().parent.parent;I=L/'inputs'
spec=json.load(open(L/'code/RELATEDNESS_SPEC.json'));c=sorted(spec['candidates'])
ped={r['SampleID']:r for r in csv.DictReader(open(I/'20130606_g1k_3202_samples_ped_population.txt'),delimiter=' ')}
idx={}
for fn,sub in [('1000G_2504_high_coverage.sequence.index','data'),('698.index','additional_698_related/data')]:
  for l in open(I/fn):
    if l.startswith('#'):continue
    t=l.rstrip('\n').split('\t')
    if len(t)>10 and t[0].endswith('.final.cram'):
      run=t[2];idx[t[9]]=dict(url=f'https://s3.amazonaws.com/1000genomes/1000G_2504_high_coverage/{sub}/{run}/{t[9]}.final.cram',md5=t[1],run=run,ftp=t[0])
asm={a['donor']:a for a in json.load(open(I/'ASSEMBLIES.json'))}
pins={}
for l in open(I/'INPUT_PINS.txt'):
  x=l.split()
  if len(x)==2 and len(x[0])==64:pins[x[1]]=x[0]
Z='/ibex/scratch/projects/c2014/rob/dogohla-benchmark/hprc-zenodo-20261004'
donors=[]
for i,d in enumerate(c):
  fam=sorted({s for s,r in ped.items() if r['FamilyID']==ped[d]['FamilyID']}|({ped[d]['FatherID'],ped[d]['MotherID']}-{'0'}))
  donors.append(dict(donor=d,population=ped[d]['Population'],superpopulation=ped[d]['Superpopulation'],family_1000G=ped[d]['FamilyID'],pedigree_component=fam,panel_fold=i%5,cram_url=idx[d]['url'],cram_md5=idx[d]['md5'],cram_run=idx[d]['run'],assembly_url=asm[d]['url'],assembly_md5=asm[d]['md5'],assembly_bytes=asm[d]['bytes'],assembly_md5_manifest=asm[d]['md5_manifest']))
refs={'ipd':dict(reference='refs/ipd.original.fa',reference_sha256=pins[Z+'/baselines/assets/references/original.fa'],source=Z+'/baselines/assets/references/original.fa'),
      'genomic':dict(reference='refs/genomic.fa',reference_sha256=pins[Z+'/baselines/assets/references/genomic.fa'],source=Z+'/baselines/assets/references/genomic.fa')}
for a in ['genomic_core_hprc_asian','genomic_core_nonasian']:
  refs[a]={}
  for f in range(5):
    src=f'{Z}/genomic-reference-arms-20261018/fold{f}/linear/{a}'
    refs[a][str(f)]=dict(reference=f'refs/fold{f}/{a}/reference.fa',reference_sha256=pins[src+'/reference.fa'],aliases=f'refs/fold{f}/{a}/aliases.json',aliases_sha256=pins[src+'/aliases.json'],complete_sha256=pins[src+'/COMPLETE.json'],source=src)
R='/ibex/scratch/projects/c2014/rob/dogohla-benchmark'
plan=dict(stage='confirmatory-test-20261113',created='2026-11-13',status='frozen_before_any_test_donor_typing_or_truth_derivation',
 donors=donors,clusters=[[d['donor']] for d in donors],
 arms=['ipd','genomic','genomic_core_hprc_asian','genomic_core_nonasian'],
 references=refs,
 t1k=dict(version='T1K v1.0.6-r216',run_t1k='/benchmark/baselines/assets/t1k/run-t1k',run_t1k_sha256=pins[Z+'/baselines/assets/t1k/run-t1k'],
   other_binaries_sha256={k.split('/')[-1]:v for k,v in pins.items() if '/assets/t1k/' in k},
   container=R+'/containers/ubuntu24.sif',container_sha256=pins[R+'/containers/ubuntu24.sif'],
   python='/home/leechuck/hla/mm/envs/asian50-spechla/bin/python3.12',
   command='run-t1k -1 r1.fq.gz -2 r2.fq.gz --preset hla-wgs -f REF -t 4 --alleleDigitUnits 4 --alleleDelimiter : -o calls'),
 recruitment=dict(reference='/home/leechuck/hla/ref1kg/local/GRCh38_full_analysis_set_plus_decoy_hla.fa',reference_sha256=pins[R+'/hla/ref1kg/local/GRCh38_full_analysis_set_plus_decoy_hla.fa'],regions='chr6:28500000-33500000 (0-based half-open) + every HLA-* and chr6_* contig, full length',source_script='MHC-recruitment-expansion41-20261018/recruit.py (extract/test_rescue copied unchanged)'),
 truth=dict(annotation='Immuannot scripts.pub.v3 at commit 31362c3701 (full hash in receipts/IMMUANNOT_SETUP.txt); reference data built with its prepare-reference recipe from IPD-IMGT/HLA 3.59.0 (ANHIG/IMGTHLA branch 3590 hla.dat) + IPD-KIR 2.13.0 + c4a.gb + curated DRB3, the versions named by Zenodo 15285071 (development truth); run on the whole R3 Verkko hi-c assembly with defaults (--diff 0.03, overlap 0.9), 8 threads',
   gate='receipts/IMMUANNOT_VALIDATION.json status=pass: identical gene start/end/strand for all 8 loci on 16 development haplotypes (8 development donors) versus the published Zenodo GTF. Truth jobs depend afterok on this gate job.',
   matching='match() copied unchanged from asian-assembly-sheet-20261012/match.py against IPD-IMGT/HLA 3.65.0 *_gen.fasta (manifest commit 5b915f27f7f620361cf83cb626eeac8a03c0247c); genomic_labels = matches; core_id = GENE*PG+sha256(annotated interval)',
   eligibility='donor-gene truth eligible iff exactly two annotated copies and not both on the same Verkko haplotypeN prefix; ineligible donor-genes are excluded from all arms identically and reported'),
 named_contrasts=[['genomic','ipd'],['genomic_core_hprc_asian','genomic'],['genomic_core_hprc_asian','genomic_core_nonasian'],['genomic_core_nonasian','genomic']],
 copy_contrasts=[['genomic_core_hprc_asian','genomic_core_nonasian'],['genomic_core_hprc_asian','genomic']],
 endpoints=dict(
   primary_1=dict(contrast=['genomic','ipd'],measure='named four-field genotype accuracy (endpoints.named_pair), paired donor-gene; contrast() from score_track_a.py',confirmed_if='family bootstrap 95% lower bound > 0 AND family sign-flip p < 0.05',development=dict(diff_pp=40.30,ci=[35.99,44.65],p=5.0e-05,n=541)),
   primary_2a=dict(contrast=['genomic_core_hprc_asian','genomic_core_nonasian'],measure='exact-sequence recovery of truth copies without an IPD genomic match (endpoints.unmatched_cores), copy denominator, copy_contrast()/uncertainty() from interpretation-20261020',confirmed_if='bootstrap95_pp lower bound > 0 AND sign_flip_p < 0.05',development=dict(recovered=[39,4],copies=175,diff_pp=20.00,ci=[13.79,26.63],p=5.0e-05)),
   primary_2b=dict(contrast=['genomic_core_hprc_asian','genomic'],measure='named four-field accuracy, contrast()',non_inferiority_margin_pp=-3.0,confirmed_if='family bootstrap 95% lower bound > -3.0 percentage points',development=dict(diff_pp=-0.74,ci=[-1.83,0.19],p=0.293))),
 alpha='each of the three claims tested separately at two-sided 0.05 (2b: one-sided non-inferiority via the 95% interval); no multiplicity adjustment; all secondary results descriptive',
 relatedness_rule='receipts/RELATEDNESS_SCREEN.json (Somalier formula, flag >= 0.0884, 1000G NYGC chr21+22, every 10th SNV with AF 0.05-0.95): a test donor flagged against any development, panel or previously typed donor is excluded before scoring; flagged test-test pairs are merged into one resampling cluster',
 failures='a planned donor-arm without a completed calls table is scored as a no-call on every truth-eligible gene; failed jobs are fixed and resubmitted, never dropped; a donor whose truth cannot be derived is reported and contributes no denominator to any arm')
(L/'code/PLAN.json').write_text(json.dumps(plan,indent=1)+'\n')
print(len(donors),[(d['donor'],d['panel_fold'],d['population']) for d in donors])
