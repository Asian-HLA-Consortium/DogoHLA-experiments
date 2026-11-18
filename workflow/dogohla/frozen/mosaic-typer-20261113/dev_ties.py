"""Per SYN case (and D-sim gene with proposals): proposer decision, tied solutions, proposals, hit."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C
for cid in C.load_json(C.S / 'work/SYN/CASES.json')['cases']:
    case = C.load_json(C.syn_dir(cid) / 'CASE.json')
    g = case['gene']
    p = C.load_json(C.syn_dir(cid) / 'C/proposals.json')['genes'][g]
    props = [c for c in p['candidates'] if 'name' in c]
    hit = case['truth_sha']['1'] in {c['sha256'] for c in props}
    sols = p.get('solutions', [])
    gap = (sols[1]['evidence_cost'] + 0.1 * sols[1]['switches'] + 0.05 * sols[1]['linkage']) - \
          (sols[0]['evidence_cost'] + 0.1 * sols[0]['switches'] + 0.05 * sols[0]['linkage']) if len(sols) > 1 else None
    print(cid, case['kind'], p['decision'], 'tied', p.get('n_tied_solutions'), 'nprop', len(props), 'hit', hit,
          'best', p.get('best_evidence_cost'), p.get('best_switches'), p.get('best_linkage'),
          'zero', p.get('zero_switch_evidence_cost'), p.get('zero_switch_linkage'),
          'gap_to_2nd', None if gap is None else round(gap, 2))
