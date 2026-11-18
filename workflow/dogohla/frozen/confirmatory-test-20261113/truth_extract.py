"""Confirmatory truth: Immuannot (pinned, reproduction-gated) gene intervals on the donor's own
phased R3 Verkko assembly, then direct matching against the pinned IPD-IMGT/HLA 3.65.0 genomic
FASTA with match() copied unchanged from asian-assembly-sheet-20261012/match.py (the development
truth matcher). core = annotated gene interval in gene orientation; core_id = GENE*PG<sha256(core)>.
Per donor-gene eligibility (pre-registered): exactly two annotated copies of the gene, not both on
contigs of the same Verkko haplotype prefix. Never reads typing calls."""
import gzip,hashlib,json,re,sys
from pathlib import Path
import pysam
GENES=['A','B','C','DPA1','DPB1','DQA1','DQB1','DRB1']
digest=lambda s:hashlib.sha256(s.encode()).hexdigest()
def fasta(path):
 label=None;parts=[]
 for line in path.open():
  if line.startswith('>'):
   if label is not None:yield label,''.join(parts).upper()
   label=line.split()[1];parts=[]
  else:parts.append(line.strip())
 if label is not None:yield label,''.join(parts).upper()
def match(query,span,start,refs):
 if not query or set(query)-set('ACGT'):return dict(status='ambiguous_or_empty_assembly_sequence',matches=[])
 direct=sorted(label for label,seq in refs.items() if seq==query)
 if direct:return dict(status='exact_at_annotation_boundaries',matches=direct,matched_sequence_sha256=digest(query),matched_length=len(query))
 matches=[];lengths=set();hashes=set()
 for label,ref in refs.items():
  if len(ref)<=len(query) or set(ref)-set('ACGT'):continue
  offset=ref.find(query)
  while offset!=-1:
   lo=start-offset;hi=lo+len(ref)
   if lo>=0 and hi<=len(span) and span[lo:hi]==ref:
    matches.append(label);lengths.add(len(ref));hashes.add(digest(ref));break
   offset=ref.find(query,offset+1)
 if matches:return dict(status='exact_after_contiguous_boundary_extension',matches=sorted(set(matches)),matched_lengths=sorted(lengths),matched_sequence_sha256=sorted(hashes))
 return dict(status='no_exact_match_boundary_or_database_completeness_review_needed',matches=[])
def hap(contig):
 m=re.match(r'(haplotype[12])-',contig);return m[1] if m else 'unassigned:'+contig
def main(donor,gtf,assembly,imgt,out):
 imgt=Path(imgt);manifest=json.loads((imgt/'manifest.json').read_text());assert manifest['release']=='3.65.0'
 pins={x['path']:x['sha256'] for x in manifest['files']};refs={}
 for g in GENES:
  p=imgt/'fasta'/(g+'_gen.fasta');assert hashlib.sha256(p.read_bytes()).hexdigest()==pins['fasta/'+p.name];refs[g]=dict(fasta(p))
 copies={g:[] for g in GENES}
 with gzip.open(gtf,'rt') as f:
  for line in f:
   if line.startswith('#'):continue
   x=line.rstrip('\n').split('\t')
   if len(x)<9 or x[2]!='gene':continue
   name=re.search(r'gene_name "([^"]+)"',x[8])[1].removeprefix('HLA-')
   if name in copies:
    t=re.search(r'template_allele "([^"]+)"',x[8]);d=re.search(r'template_distance (\S+);',x[8])
    copies[name].append(dict(contig=x[0],start_1based=int(x[3]),end_inclusive=int(x[4]),strand=x[6],template_allele=t[1] if t else None,template_distance=d[1] if d else None))
 fa=pysam.FastaFile(assembly);records=[];genes={}
 comp=str.maketrans('ACGT','TGCA')
 for g in GENES:
  cs=copies[g];haps=[hap(c['contig']) for c in cs]
  eligible=len(cs)==2 and not (haps[0]==haps[1] and haps[0].startswith('haplotype'))
  reason=None if eligible else ('copies=%d'%len(cs) if len(cs)!=2 else 'both_copies_same_haplotype')
  genes[g]=dict(eligible=eligible,reason=reason,copies=len(cs),haplotypes=haps)
  for c,h in zip(cs,haps):
   L=fa.get_reference_length(c['contig']);a=c['start_1based']-1;z=c['end_inclusive']
   lo=max(0,a-2000);hi=min(L,z+2000);span=fa.fetch(c['contig'],lo,hi).upper();raw=span[a-lo:z-lo]
   assert c['strand'] in '+-'
   oriented=span if c['strand']=='+' else span.translate(comp)[::-1];start=(a-lo) if c['strand']=='+' else len(span)-(z-lo)
   query=oriented[start:start+len(raw)]
   res=match(query,oriented,start,refs[g])
   flank=oriented[start-2000:start+len(raw)+2000] if start>=2000 and start+len(raw)+2000<=len(oriented) else None
   records.append(dict(donor=donor,gene=g,haplotype=h,**c,**res,genomic_labels=res['matches'],annotation_length=len(query),core_sha256=digest(query),core_id=g+'*PG'+digest(query),sequence=query,flanked_sequence=flank,truth_eligible_gene=eligible))
 json.dump(dict(donor=donor,release='3.65.0',database_commit=manifest['commit'],genes=genes,records=records,gtf_sha256=hashlib.sha256(Path(gtf).read_bytes()).hexdigest(),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),open(out,'w'),indent=1)
 print(donor,{g:(v['eligible'],v['copies']) for g,v in genes.items()})
if __name__=='__main__':main(*sys.argv[1:])
