#!/usr/bin/env python3
# KIT-105 run 8: per-file skip comparison against KIT-104, one load mode. Usage: skipcmp.py <kit-104 tsv> <kit-105 tsv>
# A file is reported when, against KIT-104: it is new or gone; it gains a file-level skip; hidden blocks rise; blocks run
# fall; passed falls; failed or error blocks rise; skip results rise; or a skip category appears. Each report carries the
# Q-ids its file names in a quarantine skip, so the reader can check it is a recorded quarantine.
import csv, re, sys, os
TESTS = "/home/doktersmol/Projects/APAF/MultiScholaR/tests/testthat"
def load(p):
    with open(p) as f:
        return {r["file"]: r for r in csv.DictReader(f, delimiter="\t")}
def num(r, k):
    v = r.get(k, "0")
    try: return float(v)
    except (TypeError, ValueError): return 0.0
cats = ["skip_file_guard", "skip_binding_guard", "skip_environment", "skip_quarantine", "skip_fixture", "skip_other"]
old, new = load(sys.argv[1]), load(sys.argv[2])
def qids(f):
    p = os.path.join(TESTS, f)
    if not os.path.exists(p): return ""
    return ",".join(sorted(set(re.findall(r"Q-\d+", open(p, errors="replace").read()))))
tot = {"old": [0, 0, 0], "new": [0, 0, 0]}
for tag, d in (("old", old), ("new", new)):
    for r in d.values():
        tot[tag][0] += num(r, "passed"); tot[tag][1] += num(r, "failed"); tot[tag][2] += num(r, "error_blocks")
print(f"files 104 {len(old)} 105 {len(new)}; passed {tot['old'][0]:.0f} -> {tot['new'][0]:.0f}; failed {tot['old'][1]:.0f} -> {tot['new'][1]:.0f}; error blocks {tot['old'][2]:.0f} -> {tot['new'][2]:.0f}")
print("file\tchange\tdetail\tQ-ids in file")
for f in sorted(set(old) | set(new)):
    if f not in old: print(f"{f}\tnew file\tpassed {num(new[f],'passed'):.0f} failed {num(new[f],'failed'):.0f} errors {num(new[f],'error_blocks'):.0f} skips {num(new[f],'skip_results'):.0f}\t{qids(f)}"); continue
    if f not in new: print(f"{f}\tgone\t\t"); continue
    a, b = old[f], new[f]; ch = []
    if b["file_skip"] == "TRUE" and a["file_skip"] != "TRUE": ch.append("new file-level skip")
    for k, sign, label in (("hidden_blocks", 1, "hidden blocks"), ("blocks_run", -1, "blocks run"), ("passed", -1, "passed"),
                           ("failed", 1, "failed"), ("error_blocks", 1, "error blocks"), ("skip_results", 1, "skips")):
        d = num(b, k) - num(a, k)
        if d * sign > 0: ch.append(f"{label} {num(a,k):.0f}->{num(b,k):.0f}")
    for c in cats:
        if num(a, c) == 0 and num(b, c) > 0: ch.append(f"new category {c}")
    better = [f"{l} {num(a,k):.0f}->{num(b,k):.0f}" for k, s, l in (("failed", -1, "failed"), ("error_blocks", -1, "error blocks")) if (num(b,k)-num(a,k))*s > 0]
    if ch: print(f"{f}\t{'; '.join(ch)}\t{'(also fewer: ' + ', '.join(better) + ')' if better else ''}\t{qids(f)}")
    elif better: print(f"{f}\timproved\t{', '.join(better)}\t")
