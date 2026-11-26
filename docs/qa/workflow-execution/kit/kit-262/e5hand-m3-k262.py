"""KIT-262 M3 E5 hand values (local, standard library only): the expected values that score-m3-k262.R's self-test compares its
output with, worked here independently of R. Writes e5hand-m3-expected.tsv (key, value) next to this file; the self-test reads that
file, so no expected number is retyped in R. Run from the repository root:
    python3 -I docs/qa/workflow-execution/kit/kit-262/e5hand-m3-k262.py
Contents: Tango's score statistic and the inverted one-sided bounds (bounds.R's form), Newcombe's paired interval (method 10), Wilson's
interval, a Student t quantile (regularised incomplete beta by continued fraction, bisection), a paired-t interval, Johnson's
statistic, the exact sign-flip p of a two-pair case (all 2^2 sign vectors), a type-7 quantile, the stage-1 plan R values.
Known constants checked on the way: qt(0.975, 5) = 2.5706, qt(0.975, 10) = 2.2281, qt(0.95, 10) = 1.8125 (to 4 decimals)."""
import itertools
import math
import statistics
import sys
from pathlib import Path
from statistics import NormalDist

nd = NormalDist()
Z95 = nd.inv_cdf(0.975)
Z90 = nd.inv_cdf(0.95)  # one-sided 0.05
out = []


def put(key, value):
    out.append((key, repr(float(value)) if not isinstance(value, str) else value))


# ---- Student t ----------------------------------------------------------------------------------------------------------
def betacf(a, b, x):
    tiny = 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    d = tiny if abs(d) < tiny else d
    d = 1.0 / d
    h = d
    for m in range(1, 400):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        d = tiny if abs(d) < tiny else d
        c = 1.0 + aa / c
        c = tiny if abs(c) < tiny else c
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        d = tiny if abs(d) < tiny else d
        c = 1.0 + aa / c
        c = tiny if abs(c) < tiny else c
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1.0) < 1e-15:
            break
    return h


def betainc(a, b, x):
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    bt = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log(1 - x))
    if x < (a + 1) / (a + b + 2):
        return bt * betacf(a, b, x) / a
    return 1 - bt * betacf(b, a, 1 - x) / b


def t_cdf(t, df):
    x = df / (df + t * t)
    p = 0.5 * betainc(df / 2, 0.5, x)
    return 1 - p if t > 0 else p


def qt(p, df):
    lo, hi = -50.0, 50.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if t_cdf(mid, df) < p:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


assert abs(qt(0.975, 5) - 2.5706) < 5e-4 and abs(qt(0.975, 10) - 2.2281) < 5e-4 and abs(qt(0.95, 10) - 1.8125) < 5e-4

# ---- Tango's score statistic and the inverted bounds (bounds.R's tangoZ, negative margins allowed) ------------------------
def tango_z(b, c, n, m):
    def score(q):
        return b / (q + m) + (c / q if c > 0 else 0.0) - 2 * (n - b - c) / (1 - 2 * q - m)
    lo = max(0.0, -m) + 1e-12
    hi = (1 - m) / 2 - 1e-12
    if score(hi) > 0:
        q = hi
    elif score(lo) < 0:
        q = lo
    else:
        a, bb = lo, hi
        for _ in range(300):
            mid = (a + bb) / 2
            if score(mid) > 0:
                a = mid
            else:
                bb = mid
        q = (a + bb) / 2
    return (b - c - n * m) / math.sqrt(n * (2 * q + m * (1 - m)))


def invert(f, lo, hi):
    flo, fhi = f(lo), f(hi)
    assert flo * fhi < 0, "no sign change"
    for _ in range(300):
        mid = (lo + hi) / 2
        if f(mid) * flo > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


put("tango_z.b0_c0_n100_m0.10", tango_z(0, 0, 100, 0.10))
put("tango_z.b40_c0_n100_m0.10", tango_z(40, 0, 100, 0.10))
b, c, n = 10, 5, 100
est = (b - c) / n
put("bound_upper.b10_c5_n100", invert(lambda m: tango_z(b, c, n, m) + Z90, est + 1e-9, 0.99))
put("bound_lower.b10_c5_n100", invert(lambda m: tango_z(b, c, n, m) - Z90, -0.99, est - 1e-9))

# ---- Wilson and Newcombe's paired interval ----------------------------------------------------------------------------------
def wilson(x, n_, z):
    p = x / n_
    d = 1 + z * z / n_
    centre = (p + z * z / (2 * n_)) / d
    half = z * math.sqrt(p * (1 - p) / n_ + z * z / (4 * n_ * n_)) / d
    return centre - half, centre + half


def newcombe_paired(a, b_, c_, d_):
    """Newcombe (1998) method 10 for p1 - p2 with p1 = (a + b)/n (candidate), p2 = (a + c)/n (reference); a both, b candidate only,
    c reference only, d neither."""
    n_ = a + b_ + c_ + d_
    p1, p2 = (a + b_) / n_, (a + c_) / n_
    l1, u1 = wilson(a + b_, n_, Z95)
    l2, u2 = wilson(a + c_, n_, Z95)
    A = (a + b_) * (c_ + d_) * (a + c_) * (b_ + d_)
    if A == 0:
        phi = 0.0
    else:
        diff = a * d_ - b_ * c_
        if diff > n_ / 2:
            C = diff - n_ / 2
        elif diff >= 0:
            C = 0.0
        else:
            C = diff
        phi = C / math.sqrt(A)
    lo = (p1 - p2) - math.sqrt((p1 - l1) ** 2 + (u2 - p2) ** 2 - 2 * phi * (p1 - l1) * (u2 - p2))
    hi = (p1 - p2) + math.sqrt((u1 - p1) ** 2 + (p2 - l2) ** 2 - 2 * phi * (u1 - p1) * (p2 - l2))
    return lo, hi, phi


lo, hi, phi = newcombe_paired(10, 5, 3, 82)
put("newcombe_lo.a10_b5_c3_d82", lo)
put("newcombe_hi.a10_b5_c3_d82", hi)
put("newcombe_phi.a10_b5_c3_d82", phi)
lo, hi, phi = newcombe_paired(0, 7, 0, 3)  # a zero marginal: phi = 0 by convention
put("newcombe_lo.a0_b7_c0_d3", lo)
put("newcombe_hi.a0_b7_c0_d3", hi)
lo, hi = wilson(5, 40, Z95)
put("wilson_lo.x5_n40", lo)
put("wilson_hi.x5_n40", hi)

# ---- paired-t interval and bounds for a fixed vector ----------------------------------------------------------------------
d = [0.1, 0.2, 0.0, 0.3, 0.1, 0.2]
mean, sd = statistics.mean(d), statistics.stdev(d)
put("t_mean.d6", mean)
put("t_sd.d6", sd)
put("t_half95.d6", qt(0.975, 5) * sd / math.sqrt(6))
put("t_upper_bound.d6", mean + qt(0.95, 5) * sd / math.sqrt(6))
put("t_lower_bound.d6", mean - qt(0.95, 5) * sd / math.sqrt(6))

# ---- Johnson's statistic (decide.R's effectTest) -----------------------------------------------------------------------------
dj = [0.0] * 170 + [0.20] * 30
e = [x - 0.05 for x in dj]
m_e = statistics.mean(e)
s_e = statistics.stdev(e)
mu3 = sum((x - m_e) ** 3 for x in e) / len(e)
g = mu3 / (sum((x - m_e) ** 2 for x in e) / len(e)) ** 1.5
put("johnson_skew.d170_30", g)
put("johnson_t.d170_30", (m_e + mu3 / (6 * s_e ** 2 * len(e)) + mu3 * m_e ** 2 / (3 * s_e ** 4)) / (s_e / math.sqrt(len(e))))

# ---- the exact sign-flip p of a two-pair case: d = (0.03, 0.02), margin 0.05, statistic t of e = d - m ---------------------------
def t_stat(v):
    n_ = len(v)
    m_ = sum(v) / n_
    s_ = math.sqrt(sum((x - m_) ** 2 for x in v) / (n_ - 1))
    return m_ / (s_ / math.sqrt(n_))


ee = [0.03 - 0.05, 0.02 - 0.05]
t_obs = t_stat(ee)
cnt = sum(1 for sg in itertools.product([-1, 1], repeat=2) if t_stat([s_ * x for s_, x in zip(sg, ee)]) <= t_obs + 1e-12)
put("signflip_p.n2_003_002", cnt / 4)

# ---- a type-7 quantile (R's default) of 1..101 --------------------------------------------------------------------------------
def q7(xs, p):
    xs = sorted(xs)
    h = (len(xs) - 1) * p
    lo_ = math.floor(h)
    hi_ = min(lo_ + 1, len(xs) - 1)
    return xs[lo_] + (h - lo_) * (xs[hi_] - xs[lo_])


v = list(range(1, 102))
put("quantile7_lo.1to101", q7(v, 0.025))
put("quantile7_hi.1to101", q7(v, 0.975))
v2 = [0.0, 0.1, 0.1, 0.3, 0.9]
put("quantile7_lo.five", q7(v2, 0.025))
put("quantile7_hi.five", q7(v2, 0.975))

# ---- the stage-1 plan (reestimate-m3-k262.R): smallest multiple of 50 >= floor with P(R) >= 0.80^(1/2), at most 1,000 ------------
t_target = 0.80 ** 0.5


def pick(cell_kind, value, floor):
    for R in range(floor, 1001, 50):
        if cell_kind == "effect":
            p = nd.cdf(0.05 * math.sqrt(R) / value - Z90)
        else:
            p = nd.cdf(0.10 * math.sqrt(R) / math.sqrt(value) - Z90)
        if p >= t_target:
            return R
    return 1000


for sdv in (0.265, 0.35, 0.40, 0.10):
    put("plan_R.effect_sd%.3f" % sdv, pick("effect", sdv, 200))
for share in (0.30, 0.01):
    put("plan_R.null_share%.2f" % share, pick("null", share, 250))

# ---- calibration ratios of a constructed table (ratio of sums, arithmetic only) ------------------------------------------------
tested = [[10, 10, 10], [10, 10, 10], [10, 10, 10], [10, 10, 10]]
rejected_p = [[1, 0, 0], [0, 1, 0], [0, 0, 0], [1, 1, 0]]
rejected_c = [[1, 1, 0], [0, 1, 1], [0, 0, 1], [1, 1, 1]]
for i in range(3):
    put("calib_ratio_P.third%d" % (i + 1), sum(r[i] for r in rejected_p) / sum(t[i] for t in tested))
    put("calib_ratio_C1.third%d" % (i + 1), sum(r[i] for r in rejected_c) / sum(t[i] for t in tested))
    put("calib_diff.third%d" % (i + 1), sum(r[i] for r in rejected_c) / sum(t[i] for t in tested) - sum(r[i] for r in rejected_p) / sum(t[i] for t in tested))

path = Path(__file__).with_name("e5hand-m3-expected.tsv")
path.write_text("key\tvalue\n" + "".join("%s\t%s\n" % kv for kv in out))
print("wrote %d expected values to %s" % (len(out), path))
for k, vv in out:
    print("%-34s %s" % (k, vv))
sys.exit(0)
