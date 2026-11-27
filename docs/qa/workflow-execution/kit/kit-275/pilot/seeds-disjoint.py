# KIT-275 r2.1: the pilot's seeds (r = 9601 to 9650) and L1-r2's (r = 9701 to 9703) are disjoint from every seed any earlier step used or the revision
# plans to use. Adapted from diag/seeds-disjoint.py (standard library only: python3 -I seeds-disjoint.py; exits 1 when a seed is shared; no random draws).
# Streams (sim-k275.R, ticket 5.3): generator 275000000 + 100000 c + 10000 s + r (s = the scenario's index in N0, N1, E0, E1, N2, E2; tier D, the diagnostic and the
# pilot use N1, s = 2); planting 276000000 + 100000 c + r; bootstrap 277000000 + 100 c + t. All cell indices c = 0 to 11 and all scenarios s = 1 to 6 are included
# for every range, as the conservative case. Everything is pooled as one set of integers (a generator seed and a planting seed are both `set.seed` values).
import sys

def generator(c, s, r): return 275000000 + 100000 * c + 10000 * s + r
def planting(c, r): return 276000000 + 100000 * c + r
def bootstrap(c, t): return 277000000 + 100 * c + t

def seeds_of(r_values, cells=range(0, 12), scenarios=range(1, 7), label=""):
    out = {}
    for c in cells:
        for r in r_values:
            for s in scenarios:
                out.setdefault(generator(c, s, r), []).append("%s generator c=%d s=%d r=%d" % (label, c, s, r))
            out.setdefault(planting(c, r), []).append("%s planting c=%d r=%d" % (label, c, r))
    return out
def pool(*dicts):
    out = {}
    for d in dicts:
        for k, v in d.items(): out.setdefault(k, []).extend(v)
    return out

# earlier steps (r1): tier D 1 to 200, other tiers 1 to 100, ID 1 to 20, G1 1 to 20 (c = 5, 6, 7), L1 and the smokes 9101 to 9103, the proving send 9201, 9202,
# the diagnostic 9501 to 9550; the bootstrap band (t = 0 to 199)
earlier = pool(seeds_of(range(1, 201), label="r1"), seeds_of([9101, 9102, 9103, 9201, 9202], label="r1 L1/smoke/prove"), seeds_of(range(9501, 9551), label="diagnostic"))
for c in range(0, 12):
    for t in range(0, 200):
        earlier.setdefault(bootstrap(c, t), []).append("bootstrap c=%d t=%d" % (c, t))
# the revision's counted ranges (r2): D 2001 to 2200, other tiers 2001 to 2100, ID 2001 to 2020 (so 2001 to 2200 covers all)
counted = seeds_of(range(2001, 2201), label="r2 counted")
pilot = seeds_of(range(9601, 9651), cells=(1, 2, 3), scenarios=(2,), label="pilot")
l1r2 = seeds_of(range(9701, 9704), label="L1-r2")

bad = 0
def report(name, mine, others):
    global bad
    shared = sorted(set(mine) & set(others))
    print("%s: %d seeds; shared with the earlier, counted and other ranges: %d" % (name, len(mine), len(shared)))
    for seed in shared[:10]:
        print("  %d: %s; other %s" % (seed, mine[seed], others[seed]))
    bad += len(shared)
assert all(len(v) == 1 for v in pilot.values()), "the pilot's own seeds collide"
assert not (set(range(9601, 9651)) & set(range(1, 201)) | set(range(9601, 9651)) & {9101, 9102, 9103, 9201, 9202} | set(range(9601, 9651)) & set(range(9501, 9551))
            | set(range(9601, 9651)) & set(range(2001, 2201)) | set(range(9601, 9651)) & {9701, 9702, 9703}), "r ranges overlap"
report("pilot (r 9601 to 9650; c = 1, 2, 3; s = 2; generator and planting)", pilot, pool(earlier, counted, l1r2))
report("L1-r2 (r 9701 to 9703; all cells and scenarios)", l1r2, pool(earlier, counted, seeds_of(range(9601, 9651), label="pilot (all cells)")))
print("earlier steps: %d seeds; counted r2 ranges: %d seeds" % (len(earlier), len(counted)))
sys.exit(1 if bad else 0)
