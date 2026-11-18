"""Confirmatory scorer. Endpoint functions are imported unchanged from the development code
(endpoints.py named_pair/unmatched_cores; score_native.py decode_slots). contrast() is copied
verbatim from track-a-scoring-20261018/score_track_a.py; copy_contrast() and uncertainty() are
copied verbatim from interpretation-20261020 (write_report.py / score_ablation.py).
Pre-registered deviations from the development scorer: (1) truth comes from TRUTH.json files;
(2) a planned donor-arm without a completed calls table is scored as a no-call on every gene
(attempted-donor denominators), instead of being dropped; (3) clusters are test-donor family
components from PLAN.json.
Usage: confirm_score.py STAGE_DIR OUT.json"""
import collections,hashlib,itertools,json,math,random,sys
from pathlib import Path
B=Path(__file__).resolve().parent;sys.path.insert(0,str(B))
from endpoints import named_pair,unmatched_cores
from score_native import decode_slots
GENES=['A','B','C','DPA1','DPB1','DQA1','DQB1','DRB1']
load=lambda p:json.loads(Path(p).read_text())
NOCALL=dict(labels=[],label_unknown=True,cores=[],core_unknown=True)
# ---- verbatim from score_track_a.py ----
def sign_test(a,b):
 n=a+b
 if n==0:return 1.0
 k=min(a,b);p=sum(math.comb(n,i) for i in range(k+1))/2**n
 return min(1.0,2*p)
def make_contrast(rows,donors,group_of):
 idx={(r['donor'],r['gene'],r['arm']):r for r in rows}
 def contrast(a,b,cohort='all'):
  pairs=[(d,g) for d in donors for g in GENES if (d,g,a) in idx and (d,g,b) in idx and idx[d,g,a]['four_field']['eligible'] and (cohort=='all' or idx[d,g,a]['cohort']==cohort)]
  ga=sum(idx[d,g,a]['four_field']['correct'] and not idx[d,g,b]['four_field']['correct'] for d,g in pairs)
  gb=sum(idx[d,g,b]['four_field']['correct'] and not idx[d,g,a]['four_field']['correct'] for d,g in pairs)
  clusters=collections.defaultdict(list)
  for d,g in pairs:clusters[group_of.get(d,(d,))].append(int(idx[d,g,a]['four_field']['correct'])-int(idx[d,g,b]['four_field']['correct']))
  family_sums=[sum(v) for v in clusters.values()];obs=abs(sum(family_sums));nonzero=[x for x in family_sums if x]
  rngp=random.Random(20261020);nperm=20000
  family_p=(1+sum(abs(sum(x*rngp.choice([-1,1]) for x in nonzero))>=obs for _ in range(nperm)))/(nperm+1) if nonzero else 1.0
  keys=sorted(clusters);rng=random.Random(20261018);boots=[]
  for _ in range(4000):
   s=[rng.choice(keys) for _ in keys];vals=[v for k in s for v in clusters[k]];boots.append(100*sum(vals)/len(vals) if vals else 0)
  boots.sort()
  return dict(a=a,b=b,cohort=cohort,pairs=len(pairs),a_only=ga,b_only=gb,diff_pp=100*(ga-gb)/len(pairs) if pairs else None,sign_test_p=sign_test(ga,gb),sign_test_scope='exploratory gene-pair sign test',family_sign_flip_p=family_p,families=len(clusters),boot95=[boots[100],boots[3899]] if pairs else None)
 return contrast
# ---- verbatim from score_ablation.py ----
def uncertainty(pairs, group_of):
    clusters = collections.defaultdict(list)
    ungrouped = sorted({r["donor"] for r in pairs if r["donor"] not in group_of})
    for r in pairs:
        clusters[group_of.get(r["donor"], (r["donor"],))].append(r["difference"])
    values = [(sum(v), len(v)) for v in clusters.values()]
    if not values:
        return {"family_clusters": 0, "difference_pp": None, "bootstrap95_pp": None,
                "sign_flip_p": None, "singleton_fallback_donors": ungrouped}
    rng = random.Random(20261020)
    boots = []
    for _ in range(4000):
        draw = [rng.choice(values) for _ in values]
        boots.append(100 * sum(x[0] for x in draw) / sum(x[1] for x in draw))
    boots.sort()
    differences = [v[0] for v in values if v[0]]
    observed = abs(sum(differences))
    if len(differences) <= 18:
        outcomes = [sum(s * d for s, d in zip(signs, differences))
                    for signs in itertools.product((-1, 1), repeat=len(differences))]
        pvalue = sum(abs(x) >= observed for x in outcomes) / len(outcomes)
        pmethod = "exact_family_sign_flip"
    else:
        hits = sum(abs(sum(rng.choice((-1, 1)) * d for d in differences)) >= observed
                   for _ in range(20000))
        pvalue = (hits + 1) / 20001
        pmethod = "monte_carlo_family_sign_flip_20000"
    return {"family_clusters": len(values),
            "singleton_fallback_donors": ungrouped,
            "difference_pp": 100 * sum(x[0] for x in values) / sum(x[1] for x in values),
            "bootstrap95_pp": [boots[100], boots[3899]],
            "sign_flip_p": pvalue, "sign_flip_method": pmethod}
# ---- verbatim from write_report.py ----
def copy_contrast(rows, a, b, group_of):
    """Family-resampled recovery-rate difference; denominator is truth copies."""
    index = {(r["donor"], r["gene"], r["arm"]): r for r in rows}
    pairs, recovered_a, recovered_b = [], 0, 0
    for row in rows:
        if row["arm"] != a:
            continue
        other = index.get((row["donor"], row["gene"], b))
        if other is None:
            continue
        left, right = row["unmatched"], other["unmatched"]
        n = left["eligible_copies"]
        assert n == right["eligible_copies"]
        delta = left["exact_copies_recovered"] - right["exact_copies_recovered"]
        assert abs(delta) <= n
        recovered_a += left["exact_copies_recovered"]
        recovered_b += right["exact_copies_recovered"]
        differences = [1 if delta > 0 else -1] * abs(delta) + [0] * (n - abs(delta))
        pairs.extend({"donor": row["donor"], "difference": d} for d in differences)
    return {"a": a, "b": b, "eligible_copies": len(pairs),
            "recovered_a": recovered_a, "recovered_b": recovered_b,
            "donors_with_unnamed_truth": len({r["donor"] for r in pairs}),
            **uncertainty(pairs, group_of)}
# ---- confirmatory driver ----
def main(stage,outp):
 S=Path(stage);plan=load(B/'PLAN.json');arms=plan['arms']
 donors=[r['donor'] for r in plan['donors']];fold={r['donor']:r['panel_fold'] for r in plan['donors']};pop={r['donor']:r['population'] for r in plan['donors']}
 # Pre-registered relatedness rule (applied from the genotype screen, never from typing results):
 # a test donor flagged against any development/panel donor is excluded; flagged test-test pairs
 # are merged into one resampling cluster.
 screen=load(S/'receipts'/'RELATEDNESS_SCREEN.json');excluded=sorted({f['a'] for f in screen['flags'] if f['b_role']!='test'})
 parent={d:d for d in donors}
 def find(x):
  while parent[x]!=x:x=parent[x]
  return x
 for g in plan['clusters']:
  for d in g[1:]:parent[find(d)]=find(g[0])
 for f in screen['flags']:
  if f['b_role']=='test':parent[find(f['a'])]=find(f['b'])
 comp=collections.defaultdict(list)
 for d in donors:comp[find(d)].append(d)
 group_of={d:tuple(sorted(v)) for v in comp.values() for d in v}
 donors=[d for d in donors if d not in excluded]
 rows=[];truth_status={};run_status={}
 for d in donors:
  tp=S/'truth'/d/'TRUTH.json'
  if not tp.exists():truth_status[d]='truth_missing';continue
  T=load(tp);truth_status[d]={g:v for g,v in T['genes'].items()}
  for arm in arms:
   table=S/'runs'/d/arm/'calls_genotype.tsv';done=(table.parent/'COMPLETE.json').exists()
   if done:
    rec=load(table.parent/'COMPLETE.json');assert rec['status']=='complete' and rec['native_table_sha256']==hashlib.sha256(table.read_bytes()).hexdigest()
   run_status[d+'/'+arm]='complete' if done else 'failed_or_missing_scored_as_no_call'
   if arm in ['ipd','genomic']:aliases={}
   else:aliases=load(S/plan['references'][arm][str(fold[d])]['aliases'])
   calls={l.split('\t')[0].removeprefix('HLA-'):decode_slots(l,aliases) for l in table.read_text().splitlines() if l and not l.startswith('#')} if done else {}
   for g in GENES:
    if not T['genes'][g]['eligible']:continue
    copies=[{k:r[k] for k in ['haplotype','genomic_labels','core_id','core_sha256']} for r in T['records'] if r['gene']==g];assert len(copies)==2
    pred=calls.get(g,[NOCALL]*2)
    rows.append(dict(donor=d,panel_fold=fold[d],population=pop[d],region='SAS' if pop[d] in ['BEB','GIH','ITU','PJL','STU'] else 'EAS',cohort='confirmatory',arm=arm,gene=g,four_field=named_pair(pred,copies),unmatched=unmatched_cores(pred,copies),call_complete=done))
 contrast=make_contrast(rows,donors,group_of)
 def totals(rs):
  return dict(donors=len({r['donor'] for r in rs}),eligible=sum(r['four_field']['eligible'] for r in rs),correct=sum(r['four_field']['correct'] for r in rs),called_on_eligible=sum(r['four_field']['eligible'] and r['four_field']['called'] for r in rs),unmatched_copies=sum(r['unmatched']['eligible_copies'] for r in rs),unmatched_recovered=sum(r['unmatched']['exact_copies_recovered'] for r in rs))
 summary={a:totals([r for r in rows if r['arm']==a]) for a in arms}
 by_gene={a:{g:totals([r for r in rows if r['arm']==a and r['gene']==g]) for g in GENES} for a in arms}
 by_region={reg:{a:totals([r for r in rows if r['arm']==a and r['region']==reg]) for a in arms} for reg in ['EAS','SAS']}
 named=[contrast(a,b) for a,b in plan['named_contrasts']]
 gene_named={g:[make_contrast([r for r in rows if r['gene']==g],donors,group_of)(a,b) for a,b in plan['named_contrasts']] for g in GENES}
 copies=[copy_contrast(rows,a,b,group_of) for a,b in plan['copy_contrasts']]
 copies_nodrb1=[copy_contrast([r for r in rows if r['gene']!='DRB1'],a,b,group_of) for a,b in plan['copy_contrasts']]
 gene_copies={g:[copy_contrast([r for r in rows if r['gene']==g],a,b,group_of) for a,b in plan['copy_contrasts']] for g in GENES}
 out=dict(status='confirmatory_scored',relatedness_excluded=excluded,clusters=sorted(set(group_of.values())),plan_sha256=hashlib.sha256((B/'PLAN.json').read_bytes()).hexdigest(),scorer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),donors_planned=len(donors),truth_status=truth_status,run_status=run_status,summary=summary,by_gene=by_gene,by_region=by_region,named_contrasts=named,named_contrasts_by_gene=gene_named,copy_contrasts=copies,copy_contrasts_excluding_DRB1=copies_nodrb1,copy_contrasts_by_gene=gene_copies,rows=rows)
 Path(outp).write_text(json.dumps(out,indent=1)+'\n')
 for a,t in summary.items():print(f"{a:26s} {t['correct']}/{t['eligible']} unmatched {t['unmatched_recovered']}/{t['unmatched_copies']}")
 for c in named:print('named',c['a'],'-',c['b'],c['pairs'],c['diff_pp'],c['boot95'],c['family_sign_flip_p'])
 for c in copies:print('copies',c['a'],'-',c['b'],c['recovered_a'],c['recovered_b'],c['eligible_copies'],c['difference_pp'],c['bootstrap95_pp'],c['sign_flip_p'])
if __name__=='__main__':main(*sys.argv[1:])
