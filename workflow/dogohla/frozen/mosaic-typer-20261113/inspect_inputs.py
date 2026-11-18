"""Print structure (field names, counts; no sequences) of the copied inputs. Runs under Slurm."""
import json, sys, collections
from pathlib import Path
S = Path('/encrypted/e3008/hohndor-dogohla-20261020/mosaic-typer-20261113')
d = json.loads((S/'inputs/SOURCE_RECORDS.json').read_text())
print('top keys', [k for k in d if k != 'records'])
for k in d:
    if k != 'records':
        s = json.dumps(d[k]); print(k, s[:800])
r = d['records'][0]
print('policy', r['core_boundary_policy'], r['core_id']); print('record keys', {k: (len(v) if isinstance(v, str) and len(v) > 60 else v) for k, v in r.items()})
print('n records', len(d['records']), collections.Counter(x.get('resource') for x in d['records']))
f = json.loads((S/'inputs/FOLDS.json').read_text())
print('FOLDS keys', {k: (type(v).__name__, len(v) if hasattr(v,'__len__') else v) for k, v in f.items()})
print('groups sample', [g for g in f['groups'] if len(g) > 1][:8])
dn = f['donors']; print('donors sample', list(dn.items())[:3] if isinstance(dn, dict) else dn[:3])
for g in ['A', 'DRB1']:
    print(g, (S/f'inputs/msa/{g}_msa_stats.json').read_text())
