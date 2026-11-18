"""Development diagnostic (SYN cases only): how do proposed MOS sequences differ from the truth?
Prints lengths, edit distance and the positions (in truth coordinates) of differences. No sequences."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import edlib
import common as C

for cid in sys.argv[1:]:
    d = C.S / 'work/SYN' / cid
    truth = C.load_json(d / 'TRUTH_SEQ.json')
    case = C.load_json(d / 'CASE.json')
    print('==', cid, case['kind'], 'truth lengths', {h: len(truth[h]['core']) for h in truth})
    for name, ints, seq in C.read_fasta_headers_and_seqs(d / 'C/additions.fa'):
        for h in truth:
            t = truth[h]['core']
            r = edlib.align(seq, t, mode='NW', task='path')
            cig = r['cigar']
            import re
            pos, diffs = 0, []
            for n, op in re.findall(r'(\d+)([=XID])', cig):
                n = int(n)
                if op in '=X' or op == 'D':
                    if op != '=':
                        diffs.append((pos, op, n))
                    pos += n
                else:
                    diffs.append((pos, op, n))
            print(f'  {name} vs hap{h}: len {len(seq)} ed {r["editDistance"]} substring_of_truth {seq in t} truth_in_cand {t in seq} diffs(first 15) {diffs[:15]}')
