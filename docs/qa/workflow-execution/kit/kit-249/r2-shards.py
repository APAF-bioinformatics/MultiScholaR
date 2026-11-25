#!/usr/bin/env python3
"""KIT-249 R2 shard count (R1-01): project each shard's job minutes from the per-file
seconds of the last completed full suite, run 45573760967 (KIT-220 R3, 482 files,
head 1c2d8fdd and base dd769671, ubuntu-22.04, reference BLAS); the three files that
merge train 3 added (test-kit174-ruv-grouping.R, test-kit212-da-current.R,
test-kit215-not-tested-lists.R) are not in it and are read from train 3's run
45573711299 (head 8bce32aa).

Usage: python3 -I r2-shards.py <file list> <setup min> <tests.tsv> [<tests.tsv> ...]
Per file: the largest seconds over the given tables (head and base). Files are sorted heaviest first and
dealt round-robin over the shards (run-tests deals the listed order round-robin).
A shard's projected job minutes = setup + the sum of its files' minutes. Prints, for
k = 3..10 shards, the largest projected shard, and the smallest k with every shard
at or below 50 minutes. stdlib only."""
import csv, sys

def secs(path):
    with open(path, newline="") as fh:
        return {r["file"]: float(r["seconds"] or 0) for r in csv.DictReader(fh, delimiter="\t")}

files = [l.strip() for l in open(sys.argv[1]) if l.strip()]
setup = float(sys.argv[2])
tables = [secs(t) for t in sys.argv[3:]]
missing = [f for f in files if not any(f in t for t in tables)]
if missing:
    sys.exit("in no table: " + ", ".join(missing))
mins = {f: max(t.get(f, 0.0) for t in tables) / 60 for f in files}
order = sorted(files, key=lambda f: -mins[f])
print(f"files {len(files)}; test minutes {sum(mins.values()):.1f}; setup {setup:.1f} min per job")
print("heaviest:", "; ".join(f"{f} {mins[f]:.1f}" for f in order[:8]))
chosen = None
for k in range(3, 11):
    shard = [setup] * k
    for i, f in enumerate(order):
        shard[i % k] += mins[f]
    print(f"k={k}: max {max(shard):.1f}  min {min(shard):.1f}  total {sum(shard):.1f} job-min per side;"
          f" shards {' '.join(f'{s:.1f}' for s in shard)}")
    if chosen is None and max(shard) <= 50:
        chosen = (k, shard)
print(f"chosen k = {chosen[0]} (largest projected shard {max(chosen[1]):.1f} min <= 50)" if chosen else "no k <= 10")
