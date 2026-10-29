#!/usr/bin/env python3
import sys, os, re, collections
b, a = sys.argv[1], sys.argv[2]
norm = lambda p: re.sub(r"[0-9a-f]{10,}", "ID", re.sub(r"([A-Za-z0-9_.-]*?-)[0-9a-f]{6,}(/|$)", r"\1ID\2", p))
def load(path):
    d = collections.defaultdict(list)
    for line in open(path):
        h, _, rel = line.rstrip("\n").partition("  ")
        d[norm(rel)].append(h)
    return d
tot = collections.Counter()
for f in sorted(os.listdir(b)):
    if not f.endswith(".hashes") or not os.path.exists(os.path.join(a, f)): continue
    B, A = load(os.path.join(b, f)), load(os.path.join(a, f))
    keys = set(B) | set(A)
    diff = sorted(k for k in keys if k in B and k in A and sorted(B[k]) != sorted(A[k]))
    ob = sorted(k for k in keys if k not in A); oa = sorted(k for k in keys if k not in B)
    same = len(keys) - len(diff) - len(ob) - len(oa)
    tot.update(same=same, diff=len(diff), only_before=len(ob), only_after=len(oa))
    print(f"{f[:-7]}: same {same} differing {len(diff)} only-before {len(ob)} only-after {len(oa)}")
    for k in diff[:8]: print("   DIFF", k)
    for k in ob[:4]: print("   ONLY-BEFORE", k)
    for k in oa[:4]: print("   ONLY-AFTER", k)
print("TOTAL", dict(tot))
