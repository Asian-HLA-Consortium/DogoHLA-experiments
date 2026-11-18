"""Deviation D1 companion analysis (written 2026-11-13 after the genotype screen and BEFORE any
truth derivation or typing output existed): rerun the frozen scorer unchanged, but without the
relatedness-screen exclusions against development/panel donors. Test-test flags still merge
clusters. Output: receipts/CONFIRMATORY_SCORES_ALL33.json."""
import json,sys
from pathlib import Path
B=Path(__file__).resolve().parent;sys.path.insert(0,str(B))
import confirm_score as cs
orig=cs.load
def patched(p):
    x=orig(p)
    if str(p).endswith('RELATEDNESS_SCREEN.json'):
        x=dict(x,flags=[f for f in x['flags'] if f['b_role']=='test'])
    return x
cs.load=patched
cs.main(sys.argv[1],sys.argv[2])
