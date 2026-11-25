# KIT-190 descriptive B+B comparison (audit "Stage 2 (counted) and the decision run", D1), fixed before decide.R's
# result is seen. Descriptive only: no test; decides nothing. Stdlib Python; run: python3 -I bb-compare.py
# N-w (250 seeds): share of seeds with any call under X+B minus the share under B+B, X = F, D.
# E-w (200 seeds): mean over seeds of FDP(X+B) - FDP(B+B), X = F, D.
# "Far worse": both differences exceed 0.10 in that cell.
import csv, os
here = os.path.dirname(os.path.abspath(__file__))
def rows(cell_dir):
    with open(os.path.join(here, 'bound-s2', cell_dir, 'journeys.tsv')) as f:
        return [r for r in csv.DictReader(f, delimiter='\t')]
def by_arm(rs, cell, arm):
    return {int(float(r['r'])): r for r in rs if r['cell'] == cell and r['arm'] == arm and r['error'] in ('', 'NA')}
out = []
for cell, d, kind in (('N-w', 'nw', 'any_call'), ('E-w', 'ew', 'fdp')):
    rs = rows(d)
    bb = by_arm(rs, cell, 'B+B')
    diffs = {}
    for x in ('F', 'D'):
        xb = by_arm(rs, cell, x + '+B')
        seeds = sorted(set(bb) & set(xb))
        if kind == 'any_call':
            val = lambda r: 1.0 if r['any_call'] == 'TRUE' else 0.0
        else:
            val = lambda r: float(r['fdp'])
        diffs[x] = (sum(val(xb[s]) - val(bb[s]) for s in seeds) / len(seeds), len(seeds))
    far = all(v[0] > 0.10 for v in diffs.values())
    out.append((cell, kind, diffs['F'][0], diffs['D'][0], diffs['F'][1], far))
print('cell\tendpoint\tF+B_minus_B+B\tD+B_minus_B+B\tseeds\tfar_worse')
for c, k, f, dd, n, far in out:
    print(f'{c}\t{k}\t{f:.4f}\t{dd:.4f}\t{n}\t{far}')
