# E4 hand values for the constructed table (exact fractions; the jackknife in floats).
from fractions import Fraction as F
import math
reps = {1: dict(a=8, v=2, w=1, m=F(1), b=3, T=30, o=50, VO=4, TO=200, t=F(2, 1000)),
        2: dict(a=10, v=3, w=2, m=F(-1, 2), b=2, T=34, o=52, VO=5, TO=210, t=F(3, 1000)),
        3: dict(a=9, v=1, w=1, m=F(1, 5), b=4, T=32, o=48, VO=3, TO=190, t=F(25, 10000))}
def atoms(r):
    x = reps[r]
    return dict(n0d=120, p05d=x['a'] + x['b'], Vd=x['v'] + 1, n1d=50, Td=x['T'], n0td=120 * x['t'],
                n0o=1000, p05o=x['o'], Vo=x['VO'], n1o=250, To=x['TO'], n0to=1000 * x['t'],
                m6s0=x['m'], m6n0=100, m6s1=F(2, 5), m6n1=20, away=x['w'] + 1, calls=x['v'] + 1)
def add(ws):
    out = {}
    for r, w in ws.items():
        for k, v in atoms(r).items(): out[k] = out.get(k, 0) + w * v
    return out
def f(p0, pw, t): return p0 * t / (p0 * t + (1 - p0) * pw)
def stats(S):
    a_d = F(S['p05d']) / S['n0d']; a_o = F(S['p05o']) / S['n0o']
    fdr_d = F(S['Vd'], S['Vd'] + S['Td']); fdr_o = F(S['Vo'], S['Vo'] + S['To'])
    cal_d = S['n0td'] / (S['n0td'] + S['Td']); cal_o = S['n0to'] / (S['n0to'] + S['To'])
    p0d = F(S['n0d'], S['n0d'] + S['n1d']); p0o = F(S['n0o'], S['n0o'] + S['n1o'])
    pwd = F(S['Td'], S['n1d']); pwo = F(S['To'], S['n1o']); td = S['n0td'] / S['n0d']; to = S['n0to'] / S['n0o']
    return dict(rho_O=a_d / a_o, m6_na=F(S['m6s0']) / S['m6n0'], m6_ell=F(S['m6s1']) / S['m6n1'], fdr_added=fdr_d,
                fdrcal_added=cal_d, E_cal=(fdr_d - cal_d) - (fdr_o - cal_o), E_comp=cal_d - cal_o,
                comp_power_first=f(p0o, pwd, to) - f(p0o, pwo, to), comp_null_second=f(p0d, pwd, td) - f(p0o, pwd, to),
                comp_null_first=f(p0d, pwo, td) - f(p0o, pwo, to), comp_power_second=f(p0d, pwd, td) - f(p0d, pwo, td),
                D_all=F(S['away'], S['calls']))
setting = stats(add({1: 1, 2: 1, 3: 1}))
k = {'A': 3, 'B': 3, 'C': 3, 'D': 1}
pool = stats(add({1: 10, 2: 10, 3: 10}))
assert all(setting[s] == pool[s] for s in setting)
for s, v in pool.items(): print('%s\t%.17g' % (s, float(v)))
def jack(thetas):
    n = len(thetas); mean = sum(thetas) / n
    return math.sqrt((n - 1) / n * sum((t - mean) ** 2 for t in thetas))
for stat in ('rho_O', 'E_comp'):
    units = []
    for c in 'ABCD':
        for i in (1, 2, 3):
            ws = {r: 10 - (k[c] if r == i else 0) for r in (1, 2, 3)}
            units.append(float(stats(add(ws))[stat]))
    print('jack_se_pool_%s\t%.17g' % (stat, jack(units)))
    units = [float(stats(add({r: (0 if r == i else 1) for r in (1, 2, 3)}))[stat]) for i in (1, 2, 3)]
    print('jack_se_setting_%s\t%.17g' % (stat, jack(units)))
