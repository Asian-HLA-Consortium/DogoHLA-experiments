"""Native T1K on the confirmatory donors: four pre-registered arms per donor, run sequentially.
Command, binary, preset and options identical to genomic-reference-arms-20261018/run_core.py and
typing-expansion41/run.py: run-t1k -1 r1 -2 r2 --preset hla-wgs -f REF -t CPUS
--alleleDigitUnits 4 --alleleDelimiter : -o calls. References are byte-identical stage copies
of the frozen development files (sha256 asserted against PLAN.json)."""
import hashlib,json,os,subprocess,time
from pathlib import Path
B=Path(__file__).resolve().parent;S=B.parent
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for c in iter(lambda:f.read(8*1024**2),b''):h.update(c)
 return h.hexdigest()
def main():
 assert os.environ.get('SLURM_JOB_ID');plan=json.loads((B/'PLAN.json').read_text())
 row=plan['donors'][int(os.environ['SLURM_ARRAY_TASK_ID'])];donor=row['donor'];fold=row['panel_fold']
 pp=S/'reads'/donor/'PAIRS.json';pairs=json.loads(pp.read_text());assert pairs['donor']==donor and pairs['status']=='paired_fastq_verified'
 r={m:Path(str(pairs[m]).replace('/stage/',str(S)+'/',1)) if not Path(pairs[m]).exists() else Path(pairs[m]) for m in ['r1','r2']}
 for m in ['r1','r2']:assert sha(r[m])==pairs[m+'_sha256']
 t1k=Path(plan['t1k']['run_t1k']);assert sha(t1k)==plan['t1k']['run_t1k_sha256']
 for arm in plan['arms']:
  out=S/'runs'/donor/arm
  if (out/'COMPLETE.json').exists():continue
  spec=plan['references'][arm] if arm in ['ipd','genomic'] else plan['references'][arm][str(fold)]
  ref=S/spec['reference'];assert sha(ref)==spec['reference_sha256'],(arm,ref)
  if out.exists():subprocess.run(['rm','-rf',str(out)],check=True)
  out.mkdir(parents=True);command=[str(t1k),'-1',str(r['r1']),'-2',str(r['r2']),'--preset','hla-wgs','-f',str(ref),'-t',os.environ['SLURM_CPUS_PER_TASK'],'--alleleDigitUnits','4','--alleleDelimiter',':','-o',str(out/'calls')]
  rec=dict(status='running',donor=donor,panel_fold=fold,arm=arm,job_id=os.environ['SLURM_JOB_ID'],pairs_sha256=sha(pp),reference_sha256=sha(ref),command=command,truth_loaded=False,started=time.time());(out/'RUN.json').write_text(json.dumps(rec,indent=2)+'\n')
  with (out/'run.log').open('w') as f:subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,check=True)
  t=out/'calls_genotype.tsv';assert t.stat().st_size>0
  rec.update(status='complete',finished=time.time(),native_table_sha256=sha(t));(out/'COMPLETE.json').write_text(json.dumps(rec,indent=2)+'\n')
 print('TYPE_OK',donor)
if __name__=='__main__':main()
