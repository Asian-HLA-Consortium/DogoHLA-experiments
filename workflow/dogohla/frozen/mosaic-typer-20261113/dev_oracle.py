"""Development diagnostic for the oracle arm E on SYN cases: what T1K called for the case gene,
whether the ORC_ (true) candidate appears anywhere in its output, and how the true sequence differs
from each called allele (length, edit distance, and whether differences lie only at the ends).
Numbers and names only."""
import re
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import edlib
import common as C

ipd = None
for cid in sys.argv[1:]:
    d = C.S / 'work/SYN' / cid
    case = C.load_json(d / 'CASE.json')
    g = case['gene']
    truth = C.load_json(d / 'TRUTH_SEQ.json')
    e = d / 'E'
    rec = C.load_json(e / 'COMPLETE.json')
    tab = Path(rec['calls_table'])
    print('==', cid, case['kind'], 'E', rec['pass2'], 'n_added', rec.get('n_added'))
    for line in tab.read_text().splitlines():
        if line.split('\t')[0] == f'HLA-{g}':
            print('  genotype:', '\t'.join(line.split('\t')[1:8]))
    al = tab.parent / 'calls_allele.tsv'
    rows = [l for l in al.read_text().splitlines() if l.startswith(f'HLA-{g}*')]
    print('  allele table:', rows)
    orc = C.load_json(e / 'aliases_orc.json')
    print('  ORC names:', list(orc), 'ORC in genotype table:', any(n in tab.read_text() for n in orc))
    if ipd is None:
        ipd = {n: s for n, i, s in C.read_fasta_headers_and_seqs(C.GENOMIC_FA) if C.gene_of_name(n) in C.GENES}
    called = set()
    for line in tab.read_text().splitlines():
        v = line.split('\t')
        if v[0] == f'HLA-{g}':
            for x in (v[2], v[5]):
                called |= set(x.split(',')) - {'.'}
    for h in truth:
        t = truth[h]['core']
        for n in sorted(called):
            s = ipd.get(n)
            if s is None:
                continue
            hw = edlib.align(t, s, mode='HW', task='locations')
            print(f'  truth hap{h} inside {n}: infix ed {hw["editDistance"]} at {hw["locations"][:1]} '
                  f'(allele len {len(s)}, truth len {len(t)})')
            r = edlib.align(t, s, mode='NW', task='path')
            ops = re.findall(r'(\d+)([=XID])', r['cigar'])
            pos, diffs = 0, []
            for k, op in ops:
                k = int(k)
                if op != '=':
                    diffs.append((pos, op, k))
                if op in '=XI':
                    pos += k
            inner = [x for x in diffs if 0 < x[0] < len(t) - 1]
            print(f'  truth hap{h} (len {len(t)}) vs {n} (len {len(s)}): ed {r["editDistance"]}, '
                  f'diff events {len(diffs)}, interior events {len(inner)}, first {diffs[:6]}')
