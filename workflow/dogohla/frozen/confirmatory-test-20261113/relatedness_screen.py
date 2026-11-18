"""Pre-specified genotype relatedness screen (rule of reference-relatedness-screen-rule-20261008):
r = 2*(shared_hets - 2*ibs0) / max(1, hets_a + hets_b), Somalier formula, flag r >= 0.08838834764831845,
abstain when joint informative sites < 5000. Genotypes: 1000G NYGC 3202 phased panel, chr21+chr22,
biallelic SNVs with panel AF in [0.05,0.95], every 10th such site. HLA is not read here."""
import gzip,json,sys,itertools
import numpy as np
TH=0.08838834764831845
def load(vcfs,keep):
    G=[];names=None
    for v in vcfs:
        k=0
        with gzip.open(v,'rt') as f:
            for line in f:
                if line.startswith('##'):continue
                x=line.rstrip('\n').split('\t')
                if line.startswith('#'):
                    hdr=x[9:];idx=[hdr.index(s) for s in keep if s in hdr]
                    names=[hdr[i] for i in idx] if names is None else names;continue
                if len(x[3])!=1 or len(x[4])!=1:continue
                info=dict(kv.split('=',1) for kv in x[7].split(';') if '=' in kv)
                try:af=float(info.get('AF','nan'))
                except ValueError:continue
                if not 0.05<=af<=0.95:continue
                k+=1
                if k%10:continue
                gt=x[9:]
                G.append([ (g[0]=='1')+(g[2]=='1') if len(g)>=3 and g[0] in '01' and g[2] in '01' else -1 for g in (gt[i] for i in idx)])
    return names,np.array(G,dtype=np.int8)
def main(spec,out,*vcfs):
    s=json.load(open(spec));cands=s['candidates'];comps=s['comparators']
    keep=list(dict.fromkeys(cands+comps))
    names,G=load(vcfs,keep);pos={n:i for i,n in enumerate(names)}
    rows=[];flags=[];absent=[x for x in keep if x not in pos]
    for a in cands:
        if a not in pos:continue
        for b in keep:
            if b==a or b not in pos:continue
            if b in cands and b<a:continue
            ga,gb=G[:,pos[a]],G[:,pos[b]];ok=(ga>=0)&(gb>=0);ga,gb=ga[ok],gb[ok]
            ha,hb=(ga==1),(gb==1);shared=int((ha&hb).sum());ibs0=int(((ga==0)&(gb==2)|(ga==2)&(gb==0)).sum())
            het=int(ha.sum()+hb.sum());n=int(ok.sum());r=2*(shared-2*ibs0)/max(1,het)
            row=dict(a=a,b=b,b_role='test' if b in cands else 'development_or_panel',sites=n,shared_hets=shared,ibs0=ibs0,hets_ab=het,relatedness=r,insufficient=n<5000)
            rows.append(row)
            if r>=TH:flags.append(row)
    res=dict(rule='2*(shared_hets-2*ibs0)/max(1,hets_ab) >= %r'%TH,sites=int(G.shape[0]),samples_in_panel=len(names),not_in_1000G_panel=absent,flags=flags,
             max_test_vs_dev=max([r['relatedness'] for r in rows if r['b_role']!='test'],default=None),
             insufficient=[r for r in rows if r['insufficient']],pairs=len(rows))
    json.dump(dict(res,rows=rows),open(out,'w'),indent=1);print(json.dumps({k:v for k,v in res.items() if k not in ['insufficient']})[:3000])
if __name__=='__main__':main(*sys.argv[1:])
