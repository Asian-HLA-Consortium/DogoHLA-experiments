"""Truth-blind regional MHC read recruitment for the confirmatory test donors.

extract() and test_rescue() are copied unchanged from
MHC-recruitment-expansion41-20261018/recruit.py (itself unchanged from the asian-seven batch).
Same region rule (chr6:28,500,000-33,500,000 0-based half-open plus every HLA-* and chr6_* contig,
full length), same samtools, same reference, same mate rescue, no MAPQ or truth filter.
Source: the 1000G NYGC 30x CRAM (same CRAM product as the development donors), MD5-checked
against the 1000G sequence index, read-group SM checked. The CRAM is deleted after recruitment.
Run inside the same container/env as the development recruitment.
"""
import gzip,hashlib,json,os,subprocess,sys,time,urllib.request
from pathlib import Path
import pysam
B=Path(__file__).resolve().parent
sys.path.insert(0,str(B))
from baseline_annotation_v1 import sha,save,verify_pairs
SAMTOOLS=Path('/home/leechuck/hla/mm/envs/asian50-spechla/bin/samtools')
REF=Path('/home/leechuck/hla/ref1kg/local/GRCh38_full_analysis_set_plus_decoy_hla.fa')
REF_SHA='3b103f4742abfd54938fb0333e19ad067635c8eb86f1dbf0ce44b165c4292b50'
def run(cmd,log):
    subprocess.run([str(x) for x in cmd],stdout=log,stderr=log,check=True)
def extract(tool,source,bed,out,ref,log):
    out.mkdir()
    selected=out/'selected.bam'
    run([tool,'view','--fetch-pairs','--region-file',bed,'-F','0x900','-b','-T',ref,'-o',selected,source],log)
    run([tool,'collate','-@','3','-o',out/'collated.bam',selected],log)
    run([tool,'fastq','-@','3','-F','0x900','-n','-1',out/'r1.fq.gz','-2',out/'r2.fq.gz','-0',out/'other.fq.gz','-s',out/'singletons.fq.gz',out/'collated.bam'],log)
    return verify_pairs(out/'r1.fq.gz',out/'r2.fq.gz')
def test_rescue(tool,log,D):
    t=D/'semantic-test';t.mkdir()
    header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr6','LN':10000}]}
    records=[]
    for name,first,second in [('rescued',100,5000),('both',120,160),('excluded',8000,8500)]:
        for mate,pos,other in [(1,first,second),(2,second,first)]:
            a=pysam.AlignedSegment();a.query_name=name;a.query_sequence=('ACGT'*13)[:50]
            a.flag=65 if mate==1 else 129;a.reference_id=0;a.reference_start=pos;a.mapping_quality=60;a.cigarstring='50M';a.next_reference_id=0;a.next_reference_start=other;a.query_qualities=pysam.qualitystring_to_array('I'*50);records.append(a)
    with pysam.AlignmentFile(str(t/'input.bam'),'wb',header=header) as f:
        for a in sorted(records,key=lambda a:a.reference_start):f.write(a)
    pysam.index(str(t/'input.bam'));(t/'seed.bed').write_text('chr6\t90\t220\n')
    (t/'ref.fa').write_text('>chr6\n'+'A'*10000+'\n');pysam.faidx(str(t/'ref.fa'))
    result=extract(tool,t/'input.bam',t/'seed.bed',t/'output',t/'ref.fa',log)
    assert result['pairs']==2
    with gzip.open(t/'output/r1.fq.gz','rt') as f:names={line.strip()[1:] for i,line in enumerate(f) if i%4==0}
    assert names=={'rescued','both'}
    save(t/'RESULT.json',dict(status='mate_rescue_semantics_passed',outside_region_mate_retained=True,overlapping_seed_pair_not_duplicated=True,unrelated_pair_excluded=True))
def md5(p):
    h=hashlib.md5()
    with open(p,'rb') as f:
        for c in iter(lambda:f.read(8*1024**2),b''):h.update(c)
    return h.hexdigest()
def fetch(url,dest):
    for attempt in range(5):
        r=subprocess.run(['curl','-sfL','--retry','5','-C','-','-o',str(dest),url])
        if r.returncode==0:return
        time.sleep(30)
    raise RuntimeError('download failed '+url)
def main():
    assert os.environ.get('SLURM_JOB_ID')
    plan=json.loads((B/'PLAN.json').read_text());row=plan['donors'][int(os.environ['SLURM_ARRAY_TASK_ID'])];donor=row['donor']
    D=B.parent/'reads'/donor;D.mkdir(parents=True,exist_ok=True)
    if (D/'PAIRS.json').exists():print('already complete',donor);return
    for x in ['semantic-test','reads']:
        if (D/x).exists():subprocess.run(['rm','-rf',str(D/x)],check=True)
    assert sha(REF)==REF_SHA
    cram=D/f'{donor}.final.cram';crai=D/f'{donor}.final.cram.crai'
    # The host-side sbatch downloads CRAM and CRAI (no network tools inside the container).
    got=md5(cram);assert got==row['cram_md5'],('CRAM md5 mismatch',got);assert crai.exists()
    src_sha=sha(cram)
    with (D/'run.log').open('wb') as log:
        test_rescue(SAMTOOLS,log,D)
        with pysam.AlignmentFile(str(cram),'rc',reference_filename=str(REF)) as bam:
            rg={r['SM'] for r in bam.header.to_dict()['RG']};assert rg=={donor},rg
            regions=[('chr6',28500000,33500000)]
            regions += [(n,0,l) for n,l in zip(bam.references,bam.lengths) if n.startswith('HLA-') or n.startswith('chr6_')]
        bed=D/'recruitment.bed';bed.write_text(''.join(f'{n}\t{s}\t{e}\n' for n,s,e in regions))
        metrics=extract(SAMTOOLS,cram,bed,D/'reads',REF,log)
    extras={}
    for name in ['singletons.fq.gz','other.fq.gz']:
        with gzip.open(D/'reads'/name,'rb') as f:lines=sum(1 for _ in f)
        assert lines%4==0;extras[name]=lines//4
    record={'status':'paired_fastq_verified','sample_key':f'wgs.{donor}.mhc_recruited','donor':donor,'assay':'wgs','job_id':os.environ['SLURM_JOB_ID'],'source_url':row['cram_url'],'source_md5':got,'source_sha256':src_sha,'source_bytes':cram.stat().st_size,'read_group_sample':sorted(rg),'recruitment':'chr6 fixed MHC interval plus all chr6 alternative and HLA contigs; indexed mate rescue; no MAPQ or truth filter','regions':regions,'bed_sha256':sha(bed),'semantic_test_sha256':sha(D/'semantic-test/RESULT.json'),'extra_records':extras,'truth_loaded':False,'accuracy_scored':False,**metrics}
    for mate in [1,2]:
        path=D/'reads'/f'r{mate}.fq.gz';record[f'r{mate}']=str(path);record[f'r{mate}_sha256']=sha(path)
    for x in ['selected.bam','collated.bam']:(D/'reads'/x).unlink(missing_ok=True)
    cram.unlink();crai.unlink()
    save(D/'PAIRS.json',record);print(json.dumps({k:record[k] for k in ['status','donor','pairs','extra_records']}))
if __name__=='__main__':main()
