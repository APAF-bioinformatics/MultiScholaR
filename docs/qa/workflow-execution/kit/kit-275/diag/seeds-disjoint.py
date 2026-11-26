# KIT-275 diagnostic: the diagnostic's seeds are disjoint from every seed any earlier step of the ticket used or could have used.
# Standard library only: python3 -I seeds-disjoint.py. Exits 1 when a seed is shared. No random draws.
# Streams (sim-k275.R, section 5.3 of the ticket): generator set.seed(275000000 + 100000 c + 10000 s + r), s = the scenario's index in N0, N1, E0, E1, N2, E2
# (tier D and the diagnostic use N1, s = 2); planting set.seed(276000000 + 100000 c + r); bootstrap 277000000 + 100 c + t (a separate band, listed for completeness).
# c = the cell's seed index: D3 1, D4 2, D6 3, D4M 2 (= D4), X3 5, X4 6, X6 7; the script's other cells (0 to 11) are included below as well.
import sys

def generator(c, s, r): return 275000000 + 100000 * c + 10000 * s + r
def planting(c, r): return 276000000 + 100000 * c + r
def bootstrap(c, t): return 277000000 + 100 * c + t

# every data-set index any step used: tier D 1 to 200, the other counted tiers 1 to 100, ID 1 to 20, G1 1 to 20, L1 and the smokes 9101 to 9103,
# the proving send 9201 and 9202 (runs.tsv, declarations/); taken as 1 to 200 plus those, for every cell index 0 to 11 and every scenario s = 1 to 6
USED_R = set(range(1, 201)) | {9101, 9102, 9103, 9201, 9202}
used = {}
for c in range(0, 12):
    for r in USED_R:
        for s in range(1, 7):
            used.setdefault(generator(c, s, r), []).append("generator c=%d s=%d r=%d" % (c, s, r))
        used.setdefault(planting(c, r), []).append("planting c=%d r=%d" % (c, r))
    for t in range(0, 200):
        used.setdefault(bootstrap(c, t), []).append("bootstrap c=%d t=%d" % (c, t))

# the diagnostic: r = 9500 + t, t = 1 to 50; c = 1, 2, 3 (D3, D4, D6); generator s = 2; planting
mine = {}
for c in (1, 2, 3):
    for t in range(1, 51):
        r = 9500 + t
        mine.setdefault(generator(c, 2, r), []).append("generator c=%d s=2 r=%d" % (c, r))
        mine.setdefault(planting(c, r), []).append("planting c=%d r=%d" % (c, r))
shared = sorted(set(mine) & set(used))
assert all(len(v) == 1 for v in mine.values()), "the diagnostic's own seeds collide"
assert not (set(range(9501, 9551)) & USED_R), "the diagnostic's r is in a used range"
print("diagnostic seeds: %d distinct (generator %d to %d, planting %d to %d)" % (len(mine), min(generator(c, 2, 9501) for c in (1, 2, 3)),
      max(generator(c, 2, 9550) for c in (1, 2, 3)), min(planting(c, 9501) for c in (1, 2, 3)), max(planting(c, 9550) for c in (1, 2, 3))))
print("seeds of every earlier step (all cells 0 to 11, all scenarios, r in 1 to 200, 9101 to 9103, 9201, 9202; bootstrap t 0 to 199): %d" % len(used))
print("shared: %d" % len(shared))
for seed in shared[:10]:
    print("  %d: diagnostic %s; earlier %s" % (seed, mine[seed], used[seed]))
sys.exit(1 if shared else 0)
