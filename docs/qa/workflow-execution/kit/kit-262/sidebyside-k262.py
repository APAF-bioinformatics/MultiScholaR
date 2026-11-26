"""KIT-262 re-measurement: both versions side by side for the abundance-dependent cells (stdlib only; reads committed
tables only). Version 1 = M2 (floor values before the mask), results/m2-cells.tsv; version 2 = the re-measurement
(KIT-172's order), results-rm/m2-cells.tsv. Writes results-rm/side-by-side.tsv and results-rm/side-by-side.md, and the
exact-binomial counts for the re-measurement. Run from this directory: python3 -I sidebyside-k262.py"""
import csv, math
from math import comb
v1 = {int(r["cell"]): r for r in csv.DictReader(open("results/m2-cells.tsv"), delimiter="\t")}
v2 = {int(r["cell"]): r for r in csv.DictReader(open("results-rm/m2-cells.tsv"), delimiter="\t")}
cols = [("mean_complete_features", "complete_feat", "%.1f"), ("p_no_ruv", "P_noRUV", "%d"), ("c1_no_ruv", "C1_noRUV", "%d"), ("c2_no_ruv", "C2_noRUV", "%d"),
        ("c1_not_p", "C1notP", "%d"), ("c2_not_p", "C2notP", "%d"), ("c2_not_c1", "C2notC1", "%d"), ("no_ruv_completeness_only", "compl_only", "%d"),
        ("c1_differs_among_p", "C1_changed", "%d"), ("p_ruv", "P_withRUV", "%d")]
head = ["cell", "n", "g", "F", "missing"] + [c[1] + s for c in cols for s in ("_v1", "_v2")]
rows = []
for c in sorted(v2):
    a, b = v1[c], v2[c]
    row = [c, a["n"], a["g"], a["f"], a["missing_id"]]
    for k, _, fmt in cols:
        row += [fmt % float(a[k]), fmt % float(b[k])]
    rows.append(row)
with open("results-rm/side-by-side.tsv", "w") as f:
    f.write("\t".join(head) + "\n")
    for r in rows: f.write("\t".join(map(str, r)) + "\n")
with open("results-rm/side-by-side.md", "w") as f:
    f.write("# AD cells, both versions (v1 = M2, floor values before the mask; v2 = KIT-172's order). Counts of 40 data sets per cell.\n\n")
    f.write("| " + " | ".join(head) + " |\n|" + "---|" * len(head) + "\n")
    for r in rows: f.write("| " + " | ".join(map(str, r)) + " |\n")
alpha = 0.05 / 180
tail = lambda x, p=0.05, n=40: sum(comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(x, n + 1))
x = next(i for i in range(41) if tail(i) <= alpha)
gain2 = list(csv.DictReader(open("results-rm/m2-gain.tsv"), delimiter="\t"))
gain1 = {int(r["cell"]): r for r in csv.DictReader(open("results/m2-gain.tsv"), delimiter="\t")}
out = ["exact-binomial threshold: %d of 40" % x]
for k, c in (("a", "c1_not_p"), ("b", "c2_not_p"), ("e", "c2_not_c1")):
    w2 = [r["cell"] for r in gain2 if r[k] == "TRUE"]
    e2 = [r["cell"] for r in gain2 if int(r[c]) >= x]
    w1 = [str(c_) for c_ in v2 if c_ in gain1 and gain1[c_][k] == "TRUE"]
    e1 = [str(c_) for c_ in v2 if c_ in gain1 and int(gain1[c_][c]) >= x]
    out.append("%s: v1 Wilson %d cells (%s); v1 exact %d; v2 Wilson %d cells (%s); v2 exact %d (%s)" % (k, len(w1), ",".join(w1), len(e1), len(w2), ",".join(w2), len(e2), ",".join(e2)))
open("results-rm/gain-both.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
