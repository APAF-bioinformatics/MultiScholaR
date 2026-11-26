# KIT-211 E4: the constructed tables of the scorer's lifecycle check and their hand values (standard library only).
# Usage: python3 -I e4hand-k211.py <e4 dir>
# Writes, under <e4 dir>, one row directory per scenario (the format of sim-k211.R's job directories: <cell>.tsv.gz,
# metrics.tsv, range.tsv), hand.tsv (the hand values of the base scenario, exact fractions written as %.17g) and
# expect.tsv (what each scenario must give). Every number below is chosen here; the hand values are computed from the
# definitions of ticket 6.3 and 18.3 directly on the row groups, not by any code of score-k211.R.
import gzip, itertools, math, os, sys
from fractions import Fraction as F

OUT = sys.argv[1]
VD = ['vd5', 'vd10', 'vd25']
ALLM = VD + ['mcar5', 'mcar10', 'mcar25']
CELL_MASKS = {'A': ALLM, 'B': ALLM, 'C': ALLM, 'D': ALLM, 'Av': VD, 'Ag': VD}
POOL = {'A': VD, 'B': VD, 'C': VD, 'D': ['vd25']}
MAIN_ARMS = ['S', 'A', 'B', 'Bp', 'N', 'C0']
AG_ARMS = ['S', 'A', 'B', 'Bp', 'B1', 'Bp1']
STRATA = ['one_treated', 'one_control', 'both', 'other_only', 'none']
SUMS = ['n', 'n_tested', 'n_p05', 'n_call', 'n_bh', 'sum_err', 'sum_s_err', 'sum_sgn_err']
METRIC_COLUMNS = ['cell', 'mask', 'r', 'arm', 'kept', 'added', 'calls', 't_r', 'mono_ok', 'qmethod', 'threshold', 'viol',
                  'e10_dev', 'e3_pass', 'n_planted', 'n_down_all', 'n1_band_all', 'n_onoff_all', 'calls_both_s',
                  'calls_new_s', 'gained', 'lost', 'seconds']

# group: (class, stratum, ell, band, ell0, in_o, in_s, onoff)
G = {
    'a0': ('null', 'one_treated', 0, 0, 0, False, True, False),
    'aL': ('null', 'one_treated', 1, 1, 1, False, True, False),
    'aLb': ('null', 'both', 1, 1, 1, False, True, False),
    'xL': ('null', 'one_treated', 1, 1, 1, False, False, False),
    'o0': ('null', 'none', 0, 0, 0, True, True, False),
    'o1': ('null', 'none', 0, 1, 0, True, True, False),
    'oL': ('null', 'one_control', 1, 0, 1, True, True, False),
    'pO': ('up', 'none', 0, 0, 0, True, True, False),
    'pOb': ('down', 'none', 0, 1, 0, True, True, False),
    'pA': ('up', 'one_treated', 0, 1, 0, False, True, False),
    'pOn': ('down', 'none', 0, 0, 0, True, True, True),
}
NULLS = [g for g in G if G[g][0] == 'null']


def grp(n, p05=0, call=0, bh=0, err=0, serr=0, sgn=0, tested=None):
    return dict(n=F(n), n_tested=F(n if tested is None else tested), n_p05=F(p05), n_call=F(call), n_bh=F(bh),
                sum_err=F(err), sum_s_err=F(serr), sum_sgn_err=F(sgn))


# Per replicate offset of the numerator counts of the subgroup (a deterministic variation: r = 1, 2, 3 -> 0, +1, -1).
OFF = {1: 0, 2: 1, 3: -1}


def base_arms(v, s):
    """Row groups per arm for one replicate-mask. v: the offset of aL's p05 (the replicate variation, 0 when degenerate);
    s: scenario overrides (a dict)."""
    ph_s = s.get('pL_S', 12)    # S's aL p05
    ph_n = s.get('xL_N', 10)    # N's xL p05
    ph_a = s.get('pL_A', 7)
    ph_b = s.get('pL_B', 6)
    ph_bp = s.get('pL_Bp', 6)
    pow_a = s.get('pA_A', 19)   # A's pA calls (S has 20): the power cost of A is (20 - this) / 600
    m6a = s.get('m6_A', 0)      # F3: A's aL s_err beyond the balanced value
    m6b = s.get('m6_B', 0)      # F5: B's aL s_err beyond the balanced value
    arms = {}
    arms['S'] = {
        'a0': grp(100, 5, 1), 'aL': grp(100, ph_s + v, 2, serr=3),
        'aLb': grp(40, 4, 0, serr=5), 'o0': grp(900, 45, 6), 'o1': grp(100, 5, 1), 'oL': grp(30, 3, 0, serr=F(-3, 5)),
        'pO': grp(250, call=190, err=-10, sgn=-10, bh=200), 'pOb': grp(50, call=30, err=-4, sgn=-4, bh=33),
        'pA': grp(50, call=20, err=-3, sgn=-3, bh=22), 'pOn': grp(25, call=5, err=-2, sgn=-2, bh=6)}
    arms['A'] = {
        'a0': grp(100, 5, 1), 'aL': grp(100, ph_a + v, 2, serr=2 + m6a), 'aLb': grp(40, 2, 0, serr=5),
        'o0': grp(900, 45, 5, tested=895), 'o1': grp(100, 5, 1), 'oL': grp(30, 3, 0, serr=-2),
        'pO': grp(250, call=190, err=-10, sgn=-10, bh=200), 'pOb': grp(50, call=30, err=-4, sgn=-4, bh=33),
        'pA': grp(50, call=pow_a, err=-3, sgn=-3, bh=22), 'pOn': grp(25, call=5, err=-2, sgn=-2, bh=6)}
    arms['B'] = {
        'a0': grp(100, 5, 1), 'aL': grp(100, ph_b + v, 2, serr=F(3, 2) + m6b), 'aLb': grp(40, 2, 0, serr=5),
        'xL': grp(60, 4, 0, serr=F(3, 2)), 'o0': grp(900, 45, 6), 'o1': grp(100, 5, 1), 'oL': grp(30, 3, 0, serr=-3),
        'pO': grp(250, call=191, err=-10, sgn=-10, bh=201), 'pOb': grp(50, call=30, err=-4, sgn=-4, bh=33),
        'pA': grp(50, call=20, err=-3, sgn=-3, bh=22), 'pOn': grp(40, call=10, err=-3, sgn=-3, bh=12)}
    arms['Bp'] = {
        'a0': grp(100, 5, 1), 'aL': grp(100, ph_bp + v, 2, serr=F(3, 2)), 'aLb': grp(40, 2, 0, serr=5),
        'xL': grp(55, 4, 0, serr=F(3, 2)), 'o0': grp(900, 45, 6), 'o1': grp(100, 5, 1), 'oL': grp(30, 3, 0, serr=-3),
        'pO': grp(250, call=191, err=-10, sgn=-10, bh=201), 'pOb': grp(50, call=30, err=-4, sgn=-4, bh=33),
        'pA': grp(50, call=20, err=-3, sgn=-3, bh=22), 'pOn': grp(38, call=9, err=-3, sgn=-3, bh=11)}
    arms['N'] = {
        'a0': grp(100, 5, 1), 'aL': grp(100, ph_s + v, 2, serr=3), 'aLb': grp(40, 4, 0, serr=5),
        'xL': grp(200, ph_n, 1, serr=0), 'o0': grp(900, 45, 6), 'o1': grp(100, 5, 1), 'oL': grp(30, 3, 0, serr=F(-3, 5)),
        'pO': grp(250, call=190, err=-10, sgn=-10), 'pOb': grp(50, call=30, err=-4, sgn=-4),
        'pA': grp(50, call=20, err=-3, sgn=-3), 'pOn': grp(25, call=5, err=-2, sgn=-2)}
    arms['C0'] = {
        'a0': grp(100, 5, 0), 'aL': grp(100, 7, 0), 'aLb': grp(40, 2, 0), 'xL': grp(200, 8, 0), 'o0': grp(900, 45, 0),
        'o1': grp(100, 5, 0), 'oL': grp(30, 2, 0),
        'pO': grp(250, call=200, err=-1, sgn=-1), 'pOb': grp(50, call=35, err=-1, sgn=-1),
        'pA': grp(50, call=25, err=-1, sgn=-1), 'pOn': grp(40, call=12, err=-1, sgn=-1)}
    # Ag: B1 and Bp1 (1 group required) keep more of the on/off proteins and power on them.
    arms['B1'] = {k: dict(x) for k, x in arms['B'].items()}
    arms['B1']['pO'] = grp(250, call=192, err=-10, sgn=-10, bh=204)
    arms['B1']['pOn'] = grp(40, call=14, err=-3, sgn=-3, bh=18)
    arms['Bp1'] = {k: dict(x) for k, x in arms['Bp'].items()}
    arms['Bp1']['pOn'] = grp(40, call=13, err=-3, sgn=-3, bh=16)
    return arms


def s_a0(arms, r, s):
    """F4: S's a0 s_err by replicate: +m, -m, 0 (the pooled value is zero; m is the scenario's m6_na)."""
    m = s.get('m6_na', 0)
    arms['S']['a0']['sum_s_err'] = F(m) * {1: 1, 2: -1, 3: 0}[r]
    if s.get('m6_na_const'):
        arms['S']['a0']['sum_s_err'] = F(s['m6_na_const'])


def rows_of(arm, groups, cell, mask, r, ag):
    """The rows of one replicate-mask of one arm, aggregated as sim-k211.R's aggregateSubset does."""
    out = []

    def agg(subset, member, keys):
        buckets = {}
        for name, x in groups.items():
            if x['n'] == 0 or not member(name):
                continue
            c, st, ell, band, ell0, in_o, in_s, onoff = G[name]
            full = dict(cls=c, stratum=st, ell=ell, band=band, ell0=ell0)
            key = tuple(full[k] for k in keys)
            b = buckets.setdefault(key, {k: F(0) for k in SUMS})
            for k in SUMS:
                b[k] += x[k]
        for key in sorted(buckets, key=lambda t: tuple(str(i) for i in t)):
            row = dict(cell=cell, mask=mask, r=r, arm=arm, subset=subset, **{'class': '.'}, stratum='.', ell=-1, band=-1,
                       ell0=-1)
            for k, val in zip(keys, key):
                row['class' if k == 'cls' else k] = val
            row.update(buckets[key])
            out.append(row)

    flag = lambda name, i: G[name][i]
    if arm == 'C0':
        agg('whole', lambda n: True, ['cls', 'stratum', 'ell', 'band'])
        agg('whole0', lambda n: True, ['cls', 'ell0'])
    else:
        agg('added', lambda n: not flag(n, 5), ['cls', 'stratum', 'ell'])
        agg('kept_both', lambda n: flag(n, 5), ['cls', 'stratum', 'ell'])
        agg('whole', lambda n: True, ['cls', 'stratum', 'ell', 'band'])
        agg('whole0', lambda n: True, ['cls', 'ell0'])
        agg('common', lambda n: (not flag(n, 5)) and flag(n, 6), ['cls', 'ell'])
        if ag:
            agg('onoff', lambda n: flag(n, 7), ['cls'])
    return out


def num(x):
    if isinstance(x, F):
        x = float(x)
    if isinstance(x, float):
        return '%.17g' % x
    return str(x)


ROW_COLUMNS = ['cell', 'mask', 'r', 'arm', 'subset', 'class', 'stratum', 'ell', 'band', 'ell0'] + SUMS


def metrics_of(arm, groups, cell, mask, r, s_groups, viol=0):
    kept = sum(x['n'] for x in groups.values())
    added = sum(x['n'] for g, x in groups.items() if not G[g][5])
    calls = sum(x['n_call'] for x in groups.values())
    calls_s = sum(x['n_call'] for x in s_groups.values())
    both = min(calls, calls_s)
    gained = sum(x['n'] for g, x in groups.items() if not G[g][6])
    return dict(cell=cell, mask=mask, r=r, arm=arm, kept=kept, added=added, calls=calls, t_r=F(1, 500), mono_ok='TRUE',
                qmethod='qvalue', threshold=F(31, 2), viol=viol, e10_dev=0, e3_pass='TRUE', n_planted=600,
                n_down_all=300, n1_band_all=100, n_onoff_all=40, calls_both_s=both, calls_new_s=calls - both,
                gained=gained, lost=0, seconds='NA')


def write_job(dirpath, cell_runs, tasks, ranges):
    """cell_runs: list of (cell, r, mask, arm, groups, metrics-dict)"""
    os.makedirs(dirpath, exist_ok=True)
    cells = sorted({c for c, *_ in cell_runs}, key=lambda c: ['A', 'B', 'C', 'D', 'Av', 'Ag'].index(c))
    for cell in cells:
        with gzip.open(os.path.join(dirpath, cell + '.tsv.gz'), 'wt') as fh:
            fh.write('\t'.join(ROW_COLUMNS) + '\n')
            for c, r, mask, arm, groups, m, ag in cell_runs:
                if c != cell:
                    continue
                for row in rows_of(arm, groups, cell, mask, r, ag):
                    fh.write('\t'.join(num(row[k]) for k in ROW_COLUMNS) + '\n')
    with open(os.path.join(dirpath, 'metrics.tsv'), 'w') as fh:
        fh.write('\t'.join(METRIC_COLUMNS) + '\n')
        for c, r, mask, arm, groups, m, ag in cell_runs:
            fh.write('\t'.join(num(m[k]) for k in METRIC_COLUMNS) + '\n')
    with open(os.path.join(dirpath, 'range.tsv'), 'w') as fh:
        fh.write('cell\tfirst\tlast\texpected\tcomplete\n')
        for cell, first, last, expected, complete in ranges:
            fh.write('%s\t%d\t%d\t%d\t%d\n' % (cell, first, last, expected, complete))


def build(scenario_dir, s, cells=('A', 'B', 'C', 'D', 'Av', 'Ag'), vary=True, join=False, f1=False):
    all_runs = {}
    for cell in cells:
        runs = []
        reps = [1, 2, 3]
        for r in reps:
            for mask in CELL_MASKS[cell]:
                v = OFF[r] if vary else 0
                arms = base_arms(v, s)
                s_a0(arms, r, s)
                if s.get('S_override'):
                    arms['S'].update(s['S_override'](v))
                for a, ov in s.get('arm_override', {}).items():
                    arms[a] = ov(arms[a], v)
                if join and cell == 'A' and mask == 'vd5' and r in (2, 3):
                    for k in ('aL', 'aLb'):
                        arms['S'][k] = grp(0)
                ag = cell == 'Ag'
                names = AG_ARMS if ag else MAIN_ARMS
                for a in names:
                    viol = 1 if (f1 and cell == 'A' and r == 1 and mask == 'vd5' and a == 'B') else 0
                    runs.append((cell, r, mask, a, arms[a], metrics_of(a, arms[a], cell, mask, r, arms['S'], viol), ag))
        all_runs[cell] = runs
    jobs = []
    if join:
        # cell C: one range complete (r 1 to 3) and a second incomplete (r 4 and 5 of 4 to 6), joined by the scorer
        main = [x for c in cells if c != 'C' for x in all_runs[c]] + all_runs['C']
        write_job(os.path.join(scenario_dir, 'job1'), main, None, [(c, 1, 3, 3, 3) for c in cells])
        extra = []
        for r in (4, 5):
            for mask in CELL_MASKS['C']:
                arms = base_arms(0, s)
                s_a0(arms, 3, s)
                for a in MAIN_ARMS:
                    extra.append(('C', r, mask, a, arms[a], metrics_of(a, arms[a], 'C', mask, r, arms['S']), False))
        write_job(os.path.join(scenario_dir, 'job2'), extra, None, [('C', 4, 6, 3, 2)])
    else:
        write_job(os.path.join(scenario_dir, 'job1'), [x for c in cells for x in all_runs[c]], None,
                  [(c, 1, 3, 3, 3) for c in cells])
    return all_runs


# ---- hand values (base scenario) ------------------------------------------------------------------------------------
def setting_sums(all_runs, cell, mask, arm, reps=(1, 2, 3)):
    """Summed group records of one arm over the replicates of one setting -> a list of per-replicate group dicts."""
    return {r: next(g for c, rr, mk, a, g, m, ag in all_runs[cell] if rr == r and mk == mask and a == arm) for r in reps}


def S_of(per_rep, drop=None, scale=1):
    """Atoms over replicates (all except `drop`), times scale: a dict of the quantities the hand formulas use."""
    tot = {}

    def add(key, val):
        tot[key] = tot.get(key, F(0)) + scale * val

    for r, groups in per_rep.items():
        if r == drop:
            continue
        for g, x in groups.items():
            c, st, ell, band, ell0, in_o, in_s, onoff = G[g]
            null = c == 'null'
            planted = not null
            sets = ['whole'] + (['kept_both'] if in_o else ['added']) + (['common'] if (not in_o and in_s) else [])
            for sname in sets:
                if null:
                    add('n0_' + sname, x['n_tested'])
                    add('p05_' + sname, x['n_p05'])
                    add('V_' + sname, x['n_call'])
                    if ell in (0, 1):
                        add('n0_%s_l%d' % (sname, ell), x['n_tested'])
                        add('p05_%s_l%d' % (sname, ell), x['n_p05'])
                    if st in ('one_treated', 'one_control'):
                        add('m6sum_%s_%d' % (sname, ell), x['sum_s_err'])
                        add('m6n_%s_%d' % (sname, ell), x['n_tested'])
                else:
                    add('n1_' + sname, x['n_tested'])
                    add('T_' + sname, x['n_call'])
                    add('err_' + sname, x['sum_err'])
                    add('sgn_' + sname, x['sum_sgn_err'])
            if null:
                add('n0_whole0_l%d' % ell0, x['n_tested'])
                add('p05_whole0_l%d' % ell0, x['n_p05'])
            if band == 1:
                if null:
                    add('n0b', x['n_tested'])
                    add('p05b', x['n_p05'])
                else:
                    add('Tb', x['n_call'])
                    add('n1b', x['n_tested'])
                    add('sgnb', x['sum_sgn_err'])
            if planted:
                d = 'up' if c == 'up' else 'down'
                add('T_' + d, x['n_call'])
                add('bh_' + d, x['n_bh'])
            if onoff:
                add('on_kept', x['n'])
                add('on_call', x['n_call'])
                add('on_bh', x['n_bh'])
            add('kept', x['n'])
            add('tested', x['n_tested'])
            add('calls', x['n_call'])
            if not in_s:
                add('gained', x['n'])
        add('one', F(1))
        add('planted', F(600))
        add('down_all', F(300))
        add('n1b_all', F(100))
        add('onoff_all', F(40))
    return tot


def stats(S, T, arm):
    """The ratios of ticket 6.3 from summed atoms S of `arm`; T: S's own atoms (for calls shares)."""
    g = lambda k: S.get(k, F(0))
    a_o = g('p05_kept_both') / g('n0_kept_both')
    a_l = g('p05_added_l1') / g('n0_added_l1')
    out = dict(rho_l=a_l / a_o, alpha_O=a_o, alpha_L=a_l,
               rho_l_common=(g('p05_common_l1') / g('n0_common_l1') / a_o) if g('n0_common_l1') else None,
               Rc=(g('p05_whole_l1') / g('n0_whole_l1')) / (g('p05_whole_l0') / g('n0_whole_l0')),
               rho_l0=(g('p05_whole0_l1') / g('n0_whole0_l1')) / (g('p05_whole0_l0') / g('n0_whole0_l0')),
               pow=g('T_whole') / g('planted'),
               fdr_whole=g('V_whole') / (g('V_whole') + g('T_whole')),
               m6_l_added=g('m6sum_added_1') / g('m6n_added_1'), m6_na_added=g('m6sum_added_0') / g('m6n_added_0') if g('m6n_added_0') else None,
               m6_l_H=g('m6sum_whole_1') / g('m6n_whole_1'),
               alpha_band=g('p05b') / g('n0b'), pow_band=g('Tb') / g('n1b_all'), A_band=g('sgnb') / g('n1b'),
               untested_per_r=(g('kept') - g('tested')) / g('one'),
               powbh=(g('bh_up') + g('bh_down')) / g('planted'), powbh_down=g('bh_down') / g('down_all'),
               pow_down=g('T_down') / g('down_all'),
               onoff_removed=1 - g('on_kept') / g('onoff_all') if g('onoff_all') else None,
               onoff_powbh=g('on_bh') / g('onoff_all'), onoff_pow=g('on_call') / g('onoff_all'),
               calls_per_r=g('calls') / g('one'))
    return out


def derived(st, arms_present):
    """The cross-arm statistics from the per-arm statistics dict st[arm]."""
    out = {}
    if 'N' in st:
        for x in ['S', 'A', 'B', 'Bp']:
            out['Phi.' + x] = st[x]['rho_l'] / st['N']['rho_l']
            out['PhiH.' + x] = st[x]['Rc'] / st['N']['Rc']
        for x in ['A', 'B', 'Bp']:
            out['share.' + x] = (st['S']['rho_l'] - st[x]['rho_l']) / (st['S']['rho_l'] - st['N']['rho_l'])
        for x in ['S', 'A', 'B', 'Bp', 'N']:
            out['bandratio.' + x] = st[x]['alpha_band'] / st['C0']['alpha_band']
        out['rho_l0.C0'] = st['C0']['rho_l0']
    for x in st:
        if x != 'S':
            out['dpow.' + x] = st['S']['pow'] - st[x]['pow']
            out['dfdr.' + x] = st[x]['fdr_whole'] - st['S']['fdr_whole']
    if 'B1' in st:
        for b in ['B', 'Bp']:
            for s in ['pow', 'powbh_down', 'onoff_powbh', 'fdr_whole', 'onoff_removed']:
                out['d1.%s.%s' % (b, s)] = st[b + '1'][s] - st[b][s]
    return out


def jack_se(thetas):
    n = len(thetas)
    mean = sum(thetas) / n
    return math.sqrt((n - 1) / n * float(sum((t - mean) ** 2 for t in thetas)))


def hand_values(all_runs):
    rows = []
    # a setting (A, vd5): the sums over r = 1..3 of every arm
    def arm_atoms(cell, mask, arm, drop=None, scale=1):
        return S_of(setting_sums(all_runs, cell, mask, arm), drop, scale)

    def n0t(per_rep, kind, drop=None, scale=1):
        tot = F(0)
        for r, groups in per_rep.items():
            if r == drop:
                continue
            for g, x in groups.items():
                c, st, ell, band, ell0, in_o, in_s, onoff = G[g]
                if c == 'null' and ((in_o and kind == 'kept_both') or ((not in_o) and kind == 'added')):
                    tot += x['n_tested'] * F(1, 500)
        return tot * scale

    def per_arm_stats(sums_by_arm):
        st = {a: stats(S, sums_by_arm['S'], a) for a, S in sums_by_arm.items()}
        return st

    def all_stats(sums_by_arm):
        st = per_arm_stats(sums_by_arm)
        d = derived(st, list(sums_by_arm))
        flat = {}
        for a, dct in st.items():
            for k, v in dct.items():
                flat[a + '.' + k] = v
        flat.update(d)
        # R-g: E_comp = cal_d - cal_o, cal = n0t / (n0t + T)
        return flat

    main = ['S', 'A', 'B', 'Bp', 'N', 'C0']
    # setting estimates (A, vd5), pool estimates: all settings alike, so the pool's ratios equal the setting's
    est_setting = all_stats({a: arm_atoms('A', 'vd5', a) for a in main})
    est_pool = all_stats({a: arm_atoms('A', 'vd5', a, scale=10) for a in main})
    for k in est_setting:
        v1, v2 = est_setting[k], est_pool[k]
        assert v1 == v2 or (v1 is None and v2 is None), k
    wanted = ['S.rho_l', 'A.rho_l', 'B.rho_l', 'Bp.rho_l', 'N.rho_l', 'Phi.S', 'Phi.A', 'Phi.B', 'Phi.Bp', 'PhiH.A', 'share.A',
              'share.B', 'share.Bp', 'S.alpha_band', 'A.alpha_band', 'C0.alpha_band', 'bandratio.A', 'bandratio.S', 'S.pow_band',
              'S.A_band', 'rho_l0.C0', 'S.rho_l0', 'dpow.A', 'dpow.B', 'dfdr.A', 'dfdr.B', 'S.m6_l_added', 'S.m6_l_H', 'S.m6_na_added',
              'A.m6_l_H', 'B.m6_l_H', 'Bp.m6_l_H', 'S.fdr_whole', 'S.pow', 'B.rho_l_common', 'S.Rc', 'A.untested_per_r',
              'S.calls_per_r']
    for k in wanted:
        for scope, cell, mask, est in (('pool', 'ten', 'vd', est_pool), ('setting', 'A', 'vd5', est_setting)):
            v = est[k]
            rows.append((scope, cell, mask, k, 'estimate', v))
    # call share of A against S: calls both per r / S calls per r ; metrics-based: calls_both_s = min(calls, calls_S)
    def callshare(arm):
        both = F(0)
        calls_s = F(0)
        for r in (1, 2, 3):
            gs = setting_sums(all_runs, 'A', 'vd5', 'S')[r]
            ga = setting_sums(all_runs, 'A', 'vd5', arm)[r]
            cs = sum(x['n_call'] for x in gs.values())
            ca = sum(x['n_call'] for x in ga.values())
            both += min(ca, cs)
            calls_s += cs
        return both / calls_s
    for scope, cell, mask in (('pool', 'ten', 'vd'), ('setting', 'A', 'vd5')):
        rows.append((scope, cell, mask, 'callshare.A', 'estimate', callshare('A')))
    # jackknife SEs: the setting leaves one replicate out (3 units), the pool one replicate of one cell (12 units)
    for stat in ['Phi.A', 'dpow.A']:
        units = []
        for r in (1, 2, 3):
            units.append(float(all_stats({a: arm_atoms('A', 'vd5', a, drop=r) for a in main})[stat]))
        rows.append(('setting', 'A', 'vd5', stat, 'jack_se', jack_se_float(units)))
        units = []
        masks_in_pool = {'A': 3, 'B': 3, 'C': 3, 'D': 1}
        for c in 'ABCD':
            for r in (1, 2, 3):
                sums = {}
                for a in main:
                    tot_all = arm_atoms('A', 'vd5', a, scale=10)
                    leave = S_of({r: setting_sums(all_runs, 'A', 'vd5', a)[r]}, None, masks_in_pool[c])
                    sums[a] = {k: tot_all.get(k, F(0)) - leave.get(k, F(0)) for k in set(tot_all) | set(leave)}
                units.append(float(all_stats(sums)[stat]))
        rows.append(('pool', 'ten', 'vd', stat, 'jack_se', jack_se_float(units)))
    # Ag (setting Ag, vd5): the six arms
    ag = ['S', 'A', 'B', 'Bp', 'B1', 'Bp1']
    est_ag = all_stats({a: arm_atoms('Ag', 'vd5', a) for a in ag})
    for k in ['B1.pow', 'B.pow', 'd1.B.pow', 'd1.Bp.powbh_down', 'B.onoff_removed', 'S.onoff_removed', 'B.onoff_powbh',
              'S.powbh_down', 'd1.B.onoff_powbh', 'dpow.B1', 'dfdr.B1', 'A.fdr_whole']:
        rows.append(('reported', 'Ag', 'vd5', k, 'estimate', est_ag[k]))
    return rows


def jack_se_float(thetas):
    n = len(thetas)
    mean = sum(thetas) / n
    return math.sqrt((n - 1) / n * sum((t - mean) ** 2 for t in thetas))



def scenario_check(s):
    """The designed outcome of a degenerate scenario, from the rules of ticket 6.4 applied to the exact estimates (every
    replicate alike, so every limit equals its estimate)."""
    arms = base_arms(0, s)
    s_a0(arms, 3, s)
    per_rep = {1: arms}
    st = {}
    for a in MAIN_ARMS:
        st[a] = stats(S_of({1: arms[a]}), None, a)
    d = derived(st, MAIN_ARMS)
    pc = d['Phi.S'] > F(6, 5)
    fires = {'A': st['A']['m6_l_H'] != 0, 'B': st['B']['m6_l_H'] != 0, 'Bp': st['Bp']['m6_l_H'] != 0}
    f4 = (st['S']['m6_na_added'] or 0) != 0
    classes = {}
    for x in ['A', 'B', 'Bp']:
        if fires[x]:
            classes[x] = 'U'
            continue
        if d['Phi.' + x] <= F(6, 5):
            power_ok = d['dpow.' + x] <= F(1, 200)
            guard_ok = d['dfdr.' + x] <= F(1, 400)
            classes[x] = 'E' if (power_ok and guard_ok) else 'C'
        else:
            classes[x] = 'R'
    cl = list(classes.values())
    if f4:
        out = 'Stop'
    elif not pc:
        out = 'N0'
    elif cl.count('E') >= 2:
        out = 'K2'
    elif cl.count('E') == 1:
        out = 'K1'
    elif 'C' in cl:
        out = 'K3'
    elif all(c == 'R' for c in cl):
        out = 'K4'
    else:
        out = 'K5'
    return out, classes, pc, f4, fires, {k: float(v) for k, v in d.items() if k.startswith('Phi.') or k.startswith('dpow') or k.startswith('dfdr')}


def main():
    os.makedirs(OUT, exist_ok=True)
    expect = []
    # base: K2 (all three candidates remove with no cost), every replicate-mask varies in aL's p05
    base = build(os.path.join(OUT, 'base'), {})
    hand = hand_values(base)
    with open(os.path.join(OUT, 'hand.tsv'), 'w') as fh:
        fh.write('scope\tcell\tmask\tstat\tfield\tvalue\n')
        for scope, cell, mask, stat, field, v in hand:
            assert v is not None, stat
            fh.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (scope, cell, mask, stat, field, num(v)))
    expect.append(('base', 0, 'K2', 'E', 'E', 'E', '', 'TRUE', 'FALSE', 10000))
    # degenerate scenarios (cells A to D; every replicate alike): the outcome rows
    deg = [
        ('k2w', {'m6_B': 3}, 'K2', 'E', 'U', 'E', 'B', 'TRUE', 'FALSE', 0),    # B withheld (F5), A and B' E -> K2
        ('k1', {'pL_A': 20, 'pL_Bp': 20}, 'K1', 'R', 'E', 'R', '', 'TRUE', 'FALSE', 0),
        ('k3', {'pA_A': 10, 'pL_B': 20, 'pL_Bp': 20}, 'K3', 'C', 'R', 'R', '', 'TRUE', 'FALSE', 0),
        ('k4', {'pL_A': 20, 'pL_B': 20, 'pL_Bp': 20}, 'K4', 'R', 'R', 'R', '', 'TRUE', 'FALSE', 0),
        ('k5', {'pL_B': 20, 'pL_Bp': 20, 'm6_A': 3}, 'K5', 'U', 'R', 'R', 'A', 'TRUE', 'FALSE', 0),
        ('n0a', {'pL_S': 5, 'xL_N': 20}, 'N0', 'E', 'E', 'E', '', 'FALSE', 'FALSE', 0),
        ('n0b', {'pL_S': 6, 'xL_N': 11}, 'N0', 'E', 'E', 'E', '', 'FALSE', 'FALSE', 0),
        ('stopf4', {'m6_na_const': 2}, 'Stop', 'E', 'E', 'E', '', 'TRUE', 'TRUE', 0),
    ]
    base_check = scenario_check({})
    assert base_check[0] == 'K2' and list(base_check[1].values()) == ['E', 'E', 'E'], base_check
    for name, s, out, ca, cb, cbp, wh, pc, f4, boot in deg:
        got = scenario_check(s)
        assert got[0] == out and [got[1]['A'], got[1]['B'], got[1]['Bp']] == [ca, cb, cbp], (name, got)
        assert got[2] == (pc == 'TRUE') and got[3] == (f4 == 'TRUE'), (name, got)
        build(os.path.join(OUT, name), s, cells=('A', 'B', 'C', 'D'), vary=False)
        expect.append((name, 0, out, ca, cb, cbp, wh, pc, f4, 1000))
    build(os.path.join(OUT, 'stopf1'), {}, cells=('A', 'B', 'C', 'D'), vary=False, f1=True)
    expect.append(('stopf1', 1, 'Stop', '', '', '', '', '', '', 1000))
    build(os.path.join(OUT, 'join'), {}, cells=('A', 'B', 'C', 'D'), vary=True, join=True)
    expect.append(('join', 0, 'K2', 'E', 'E', 'E', '', 'TRUE', 'FALSE', 1000))
    with open(os.path.join(OUT, 'expect.tsv'), 'w') as fh:
        fh.write('scenario\texit\toutcome\tclass_A\tclass_B\tclass_Bp\twithheld\tpc_pass\tf4\tB\n')
        for row in expect:
            fh.write('\t'.join(str(x) for x in row) + '\n')
    # what the printed designed values are (for the record, not an assertion of the scorer)
    for scen in ('base', 'k2w', 'k1', 'k3', 'k4', 'k5', 'n0a', 'n0b'):
        pass


if __name__ == '__main__':
    main()
