"""Render CONFIRMATORY_RESULT.md tables from receipts/CONFIRMATORY_SCORES*.json (decision rules from PLAN.json)."""
import json,sys
from pathlib import Path
L=Path(__file__).resolve().parent.parent
plan=json.loads((L/'code/PLAN.json').read_text())
sets={'P (pre-registered exclusions)':json.loads((L/'receipts/CONFIRMATORY_SCORES.json').read_text()),
      'A (all 33 admitted)':json.loads((L/'receipts/CONFIRMATORY_SCORES_ALL33.json').read_text())}
GENES=['A','B','C','DPA1','DPB1','DQA1','DQB1','DRB1']
def f(x,n=2):return 'NA' if x is None else f'{x:+.{n}f}'
def ci(c):return 'NA' if not c else f'[{c[0]:+.2f}, {c[1]:+.2f}]'
out=[]
verdict={}
for name,r in sets.items():
    nc={(c['a'],c['b']):c for c in r['named_contrasts']};cc={(c['a'],c['b']):c for c in r['copy_contrasts']}
    e1=nc[('genomic','ipd')];e2a=cc[('genomic_core_hprc_asian','genomic_core_nonasian')];e2b=nc[('genomic_core_hprc_asian','genomic')]
    ok1=bool(e1['boot95']) and e1['boot95'][0]>0 and e1['family_sign_flip_p']<0.05
    ok2a=bool(e2a['bootstrap95_pp']) and e2a['bootstrap95_pp'][0]>0 and e2a['sign_flip_p']<0.05
    ok2b=bool(e2b['boot95']) and e2b['boot95'][0]>-3.0
    verdict[name]=(ok1,ok2a,ok2b)
    donors=len({x['donor'] for x in r['rows']})
    out.append(f'### Set {name}: {donors} donors scored; relatedness-excluded: {", ".join(r["relatedness_excluded"]) or "none"}\n')
    out.append('| Arm | Named four-field correct / eligible | IPD-absent copies recovered |\n|---|---:|---:|')
    for a,t in r['summary'].items():out.append(f"| {a} | {t['correct']}/{t['eligible']} ({100*t['correct']/max(1,t['eligible']):.1f}%) | {t['unmatched_recovered']}/{t['unmatched_copies']} |")
    out.append('\n| Endpoint | Test result | 95% CI | p | Development | Criterion met |\n|---|---|---|---|---|---|')
    out.append(f"| 1: genomic − ipd, named | {f(e1['diff_pp'])} pp (n={e1['pairs']}, +{e1['a_only']}/−{e1['b_only']}) | {ci(e1['boot95'])} | {e1['family_sign_flip_p']:.3g} | +40.30 [35.99, 44.65], n=541 | {ok1} |")
    out.append(f"| 2a: core_hprc_asian − core_nonasian, IPD-absent copy recovery | {e2a['recovered_a']} vs {e2a['recovered_b']} of {e2a['eligible_copies']} ({f(e2a['difference_pp'])} pp) | {ci(e2a['bootstrap95_pp'])} | {e2a['sign_flip_p'] if e2a['sign_flip_p'] is None else format(e2a['sign_flip_p'],'.3g')} | 39 vs 4 of 175, +20.00 [13.79, 26.63] | {ok2a} |")
    out.append(f"| 2b: core_hprc_asian − genomic, named (margin −3 pp) | {f(e2b['diff_pp'])} pp (n={e2b['pairs']}, +{e2b['a_only']}/−{e2b['b_only']}) | {ci(e2b['boot95'])} | {e2b['family_sign_flip_p']:.3g} | −0.74 [−1.83, 0.19] | {ok2b} |")
    o=nc[('genomic_core_nonasian','genomic')];out.append(f"| (descr.) core_nonasian − genomic, named | {f(o['diff_pp'])} pp | {ci(o['boot95'])} | {o['family_sign_flip_p']:.3g} | | |")
    nd={(c['a'],c['b']):c for c in r['copy_contrasts_excluding_DRB1']}[('genomic_core_hprc_asian','genomic_core_nonasian')]
    out.append(f"| (descr.) 2a excluding DRB1 | {nd['recovered_a']} vs {nd['recovered_b']} of {nd['eligible_copies']} | {ci(nd['bootstrap95_pp'])} | {nd['sign_flip_p'] if nd['sign_flip_p'] is None else format(nd['sign_flip_p'],'.3g')} | 21 vs 3 of 65 | |")
    out.append('\nPer gene (named correct/eligible per arm; IPD-absent copies recovered core_hprc_asian vs core_nonasian):\n')
    out.append('| Gene | ipd | genomic | core_hprc_asian | core_nonasian | genomic−ipd (pp) | copies rec. asian vs nonasian |\n|---|---|---|---|---|---|---|')
    for g in GENES:
        b=r['by_gene'];c=[x for x in r['named_contrasts_by_gene'][g] if x['a']=='genomic' and x['b']=='ipd'][0]
        k=[x for x in r['copy_contrasts_by_gene'][g] if x['a']=='genomic_core_hprc_asian' and x['b']=='genomic_core_nonasian'][0]
        cell=lambda a:f"{b[a][g]['correct']}/{b[a][g]['eligible']}"
        out.append(f"| {g} | {cell('ipd')} | {cell('genomic')} | {cell('genomic_core_hprc_asian')} | {cell('genomic_core_nonasian')} | {f(c['diff_pp'],1)} | {k['recovered_a']} vs {k['recovered_b']} / {k['eligible_copies']} |")
    out.append('\nRegion strata (named correct/eligible): '+'; '.join(f"{reg}: "+', '.join(f"{a} {t['correct']}/{t['eligible']}" for a,t in v.items()) for reg,v in r['by_region'].items()))
    fails=[k for k,v in r['run_status'].items() if v!='complete'];te={d:v for d,v in r['truth_status'].items()}
    inel=sum(1 for d,v in te.items() if isinstance(v,dict) for g,x in v.items() if not x['eligible'])
    out.append(f"\nRuns not completed (scored as no-calls): {len(fails)} {fails[:10]}. Truth-ineligible donor-genes: {inel}. Donors without truth: {[d for d,v in te.items() if not isinstance(v,dict)]}\n")
print('\n'.join(out))
names=['1 (genomic > ipd)','2a (Asian panel recovers IPD-absent copies)','2b (no named-accuracy cost)']
print('\n### Verdict (confirmed only if met in both sets)\n')
for i,n in enumerate(names):
    both=all(v[i] for v in verdict.values())
    print(f"- Claim {n}: **{'CONFIRMED' if both else 'NOT CONFIRMED'}** (set P: {list(verdict.values())[0][i]}, set A: {list(verdict.values())[1][i]})")
