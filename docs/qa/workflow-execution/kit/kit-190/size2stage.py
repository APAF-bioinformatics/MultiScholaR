# KIT-190 bound-declaration.md section 8: simulated size of the whole two-stage procedure at the margin (stdlib Python
# only; indicative; no package code). Run: python3 -I size2stage.py  (prints the rows of section 8's table).
# Floors are option C's R per cell type (VAL-185): null cells 100 (N-r, N-s) and 250 (N-w); effect cells 200.
# The re-estimation formula is section 8's; the final effect test includes the skewness check of section 5 (R4-01):
# |g| <= 1 the studentised t, |g| > 1 Johnson's skew-corrected t (the t reference approximates the sign-flip test).
import random, math
from statistics import NormalDist, mean, stdev
N = NormalDist(); z = N.inv_cdf(0.95); T = 0.80 ** (1 / 8)

def r_effect(s, floor, rmax=1000, m=0.05):
    if s <= 0: return floor
    for R in range(floor, rmax + 1, 50):
        if N.cdf(m * math.sqrt(R) / s - z) >= T: return R
    return rmax

def r_null(d, floor, rmax=1000, m=0.10):
    if d <= 0: return floor
    for R in range(floor, rmax + 1, 50):
        if N.cdf(m * math.sqrt(R) / math.sqrt(d) - z) >= T: return R
    return rmax

def qhat(b, c, n, d):
    lo, hi = 1e-12, (1 - d) / 2 - 1e-12
    f = lambda q: b / (q + d) + (c / q if c else 0) - 2 * (n - b - c) / (1 - 2 * q - d)
    for _ in range(100):
        mid = (lo + hi) / 2
        if f(mid) > 0: lo = mid
        else: hi = mid
    return (lo + hi) / 2

def tango(b, c, n, m=0.10):
    q = qhat(b, c, n, m)
    return (b - c - n * m) / math.sqrt(n * (2 * q + m * (1 - m))) < -z

def size_null(delta, floor, sims, seed):
    rnd = random.Random(seed); p10 = (delta + 0.10) / 2; p01 = (delta - 0.10) / 2; hits = 0; Rs = []
    def draw(n):
        b = c = 0
        for _ in range(n):
            u = rnd.random()
            if u < p10: b += 1
            elif u < p10 + p01: c += 1
        return b, c
    for _ in range(sims):
        b1, c1 = draw(50); R = r_null((b1 + c1) / 50, floor); b2, c2 = draw(R - 50); Rs.append(R)
        hits += tango(b1 + b2, c1 + c2, R)
    return hits / sims, mean(Rs)

def decide_effect(e):
    n = len(e); eb = mean(e); s = stdev(e)
    m2 = sum((v - eb) ** 2 for v in e) / n; m3 = sum((v - eb) ** 3 for v in e) / n
    g = m3 / m2 ** 1.5 if m2 > 0 else 0.0
    if abs(g) <= 1:
        t = eb / (s / math.sqrt(n))
    else:  # Johnson (1978)
        t = (eb + m3 / (6 * s ** 2 * n) + m3 * eb ** 2 / (3 * s ** 4)) / (s / math.sqrt(n))
    return t < -z

def size_effect(shape, s, sims, seed, floor=200):
    rnd = random.Random(seed); hits = 0; Rs = []
    # standardised draws with mean 0 and SD 1
    zi_p = 0.2; zi_mean = zi_p; zi_sd = math.sqrt(zi_p * 2 - zi_p ** 2)  # X = 0 w.p. 0.8, else Exp(1)
    def unit():
        if shape == 'normal': return rnd.gauss(0, 1)
        if shape == 'exponential': return rnd.expovariate(1.0) - 1.0
        x = rnd.expovariate(1.0) if rnd.random() < zi_p else 0.0
        return (x - zi_mean) / zi_sd
    draw = lambda n: [0.05 + s * unit() for _ in range(n)]
    for _ in range(sims):
        x1 = draw(50); R = r_effect(stdev(x1), floor); x = x1 + draw(R - 50); Rs.append(R)
        hits += decide_effect([v - 0.05 for v in x])
    return hits / sims, mean(Rs)

if __name__ == '__main__':
    for floor, deltas in ((100, (0.10, 0.25)), (250, (0.12, 0.25, 0.40))):  # at the margin the discordance is at least 0.10
        for d in deltas:
            print('null floor', floor, 'discordance', d, 'seed', 1000 + int(d * 1000) + floor, size_null(d, floor, 6000, 1000 + int(d * 1000) + floor))
    for shape in ('normal', 'exponential', 'zero_inflated'):
        for s in (0.20, 0.40):
            seed = 2000 + int(s * 100) + {'normal': 0, 'exponential': 1, 'zero_inflated': 2}[shape]
            print('effect', shape, 'sd', s, 'seed', seed, size_effect(shape, s, 20000, seed))
