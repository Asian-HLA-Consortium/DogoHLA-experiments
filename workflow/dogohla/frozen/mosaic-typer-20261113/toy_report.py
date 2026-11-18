"""Print the decision and costs of every toy case (diagnostic for test_proposer.py)."""
import test_proposer as T
names, rows = T.make_panel()
used = {r[1000] for r in rows}
new = next(b for b in 'ACGT' if b not in used)
A, B, C, D, L = T.A, T.B, T.Cc, T.D, T.L
cases = {
 'exact_hom': ([[(A,0,L)],[(A,0,L)]], [[A],[A]], {}),
 'exact_het': ([[(A,0,L)],[(B,0,L)]], [[A],[B]], {}),
 'exact_het_wrong_pass1': ([[(A,0,L)],[(B,0,L)]], [[A],[D]], {}),
 'mosaic1_linked': ([[(A,0,450),(B,450,L)],[(C,0,L)]], [[A],[C]], {}),
 'mosaic1_unlinked_partner_src': ([[(A,0,900),(B,900,L)],[(B,0,L)]], [[A],[B]], {}),
 'mosaic1_unlinked_partner_diff': ([[(A,0,900),(B,900,L)],[(C,0,L)]], [[A],[C]], {}),
 'mosaic2_single_variant': ([[(A,0,1700),(B,1700,1900),(A,1900,L)],[(C,0,L)]], [[A],[C]], {}),
 'private_snv': ([[(A,0,L)],[(C,0,L)]], [[A],[C]], dict(private=(0,1000,new))),
 'seq_error': ([[(A,0,L)],[(C,0,L)]], [[A],[C]], dict(errors=[(0,1800,rows[B][1800],1)])),
 'segment_lt150': ([[(A,0,L)],[(C,0,L)]], [[A],[C]], dict(private=(0,1210,rows[B][1210]))),
 'truth_in_reference': ([[(D,0,L)],[(C,0,L)]], [[A],[C]], {}),
}
for k,(h,p,kw) in cases.items():
    rec, props, truth = T.run_case(h, p, **kw)
    print(f"{k:32s} decision={rec['decision']:28s} p1={rec.get('pass1_evidence_cost')} zero={rec.get('zero_switch_evidence_cost')} best={rec.get('best_evidence_cost')} sw={rec.get('best_switches')} tied={rec.get('n_tied_solutions')} nprop={len(props)} truth_hit={[t in props for t in truth]}")
