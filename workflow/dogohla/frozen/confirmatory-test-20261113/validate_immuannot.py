"""Reproduction gate (development donors only): our Immuannot build vs the published
Zenodo 15285071 HPRC-frz2 annotation that defined the development truth intervals.
Input: per-haplotype HLA span FASTA (header '<contig>__<start0>') + the published GTF.
Pass rule (pre-set before any test donor is annotated): every one of the 8 target genes on
every validation haplotype has identical contig, 1-based start, end and strand."""
import gzip,json,re,sys
from pathlib import Path
GENES=['A','B','C','DPA1','DPB1','DQA1','DQB1','DRB1']
def genes(gtf,shift=None):
    out={}
    with gzip.open(gtf,'rt') as f:
        for line in f:
            if line.startswith('#'):continue
            x=line.rstrip('\n').split('\t')
            if len(x)<9 or x[2]!='gene':continue
            name=re.search(r'gene_name "([^"]+)"',x[8])[1].removeprefix('HLA-')
            if name not in GENES:continue
            contig,s,e=x[0],int(x[3]),int(x[4])
            if shift:
                contig,off=contig.rsplit('__',1);s+=int(off);e+=int(off)
            out.setdefault(name,[]).append((contig,s,e,x[6]))
    return out
def main(manifest,outp):
    rows=[];ok=True
    for r in json.load(open(manifest)):
        a=genes(r['ours'],shift=True);b=genes(r['published'])
        for g in GENES:
            same=sorted(a.get(g,[]))==sorted(b.get(g,[]))
            ok&=same and len(b.get(g,[]))==1
            rows.append(dict(donor=r['donor'],haplotype=r['haplotype'],gene=g,ours=a.get(g,[]),published=b.get(g,[]),identical=same))
    res=dict(status='pass' if ok else 'fail',copies=len(rows),identical=sum(x['identical'] for x in rows),rows=rows)
    Path(outp).write_text(json.dumps(res,indent=1)+'\n');print(res['status'],res['identical'],'/',res['copies'])
    for x in rows:
        if not x['identical']:print('DIFF',x)
    sys.exit(0 if ok else 3)
if __name__=='__main__':main(*sys.argv[1:])
