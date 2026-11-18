"""Write a task manifest: python make_manifest.py NAME READSET:ARM1,ARM2 [READSET:...] [--donors d1,d2]"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C


def donors_of(rs):
    runs = C.readset_runs()
    dsim = sorted(runs['D-sim'])
    if rs == 'D-sim':
        return dsim
    if rs == 'T-sim':
        return sorted(runs['T-sim'])
    if rs == 'T-real':
        return sorted(d for d in runs['real'] if d not in dsim)
    if rs == 'D-real':
        return sorted(d for d in runs['real'] if d in dsim)
    if rs in C.SYN_SETS:
        return C.load_json(C.S / 'work' / rs / 'CASES.json')['cases']
    raise ValueError(rs)


def main():
    name = sys.argv[1]
    only = None
    specs = []
    for a in sys.argv[2:]:
        if a.startswith('--donors='):
            only = set(a.split('=', 1)[1].split(','))
        else:
            specs.append(a)
    tasks = []
    for sp in specs:
        rs, arms = sp.split(':')
        for d in donors_of(rs):
            if only and d not in only:
                continue
            C.guard_path(Path(d))
            tasks.append(dict(readset=rs, donor=d, arms=arms.split(',')))
    out = C.S / 'work' / f'manifest_{name}.json'
    C.save_json(out, dict(name=name, tasks=tasks))
    print(out, len(tasks))


if __name__ == '__main__':
    main()
