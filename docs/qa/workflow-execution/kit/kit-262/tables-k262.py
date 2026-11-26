"""KIT-262 record tables (standard library only; reads only committed files; no run). Run from this directory:
python3 -I tables-k262.py
Writes results/m2-pooled.md (pooled view; mixes F = 300 and 1,000, see its header), results/m2-k-dist.tsv,
results/m2-percentage.tsv, results/m2-missing.tsv, results/m2-speed.tsv, results/m2-exact.md."""
import csv, collections, math, glob
from math import comb

R = "results/"
cells = list(csv.DictReader(open(R + "m2-cells.tsv"), delimiter="\t"))
rows = []
for p in sorted(glob.glob("m2/*/summary.tsv")):
    rows += list(csv.DictReader(open(p), delimiter="\t"))
assert len(rows) == 7920

# ---- m2-pooled.md (the tables of the first outcome; pooled over g and F) ----
out = ["# Pooled view of M2 (hand-tallied from results/m2-cells.tsv by tables-k262.py; NOT an output of score-k262.R)", "",
       "**Caveat (audit M02).** These rows pool F = 300 with F = 1,000 (and 2 with 3 groups). The two feature counts give opposite outcomes at the same",
       "complete-feature share (n = 96, MCAR 1 %: P without RUV in 66 of 80 at F = 300 and 0 of 80 at F = 1,000, complete share 0.380 and 0.382). The number of",
       "complete features (share x F) places a data set; the share alone does not. Use `m2-cells.tsv` per cell. **Rows for abundance-dependent (AD)",
       "missingness carry the M01 caveat** (the 2 % floor values are set before the mask; see outcome.md).", ""]
ag = collections.OrderedDict()
fields = ["n_valid", "p_no_ruv", "c1_no_ruv", "c2_no_ruv", "c1_not_p", "no_ruv_rescued_by_c1", "no_ruv_completeness_only", "no_ruv_no_r5",
          "r5_rows", "r5_cap_side", "r5_completeness", "p_ruv", "c1_differs_among_p", "c2_differs_among_p"]
for r in cells:
    a = ag.setdefault((int(r["n"]), r["missing_id"]), collections.Counter())
    for f in fields:
        a[f] += int(float(r[f]))
out.append("| n | missingness | data sets | P no RUV | C1 no RUV | C2 no RUV | no-RUV sets: completeness only | no-RUV sets: rescued by C1 | no-RUV sets: no R5 | R5 rows | R5 rows cap side | R5 rows completeness | C1 differs among P-with-RUV |")
out.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for (n, m), a in ag.items():
    N = a["n_valid"]
    pct = lambda x: "%d (%.0f%%)" % (x, 100 * x / N)
    out.append("| %d | %s | %d | %s | %s | %s | %d | %d | %d | %d | %d | %d | %d of %d |" % (n, m, N, pct(a["p_no_ruv"]), pct(a["c1_no_ruv"]),
               pct(a["c2_no_ruv"]), a["no_ruv_completeness_only"], a["no_ruv_rescued_by_c1"], a["no_ruv_no_r5"], a["r5_rows"], a["r5_cap_side"],
               a["r5_completeness"], a["c1_differs_among_p"], a["p_ruv"]))
d = list(csv.DictReader(open(R + "m2-resid-df.tsv"), delimiter="\t"))
ag = collections.OrderedDict()
for r in d:
    a = ag.setdefault((int(r["n"]), int(r["g"]), r["rule"]), collections.Counter())
    for f in ["with_ruv", "k_at_cap", "resid_df_le1", "resid_df_0"]:
        a[f] += int(float(r[f]))
out += ["", "| n | g | rule | data sets with RUV | k used = cap | n-g-k <= 1 | n-g-k = 0 |", "|---|---|---|---|---|---|---|"]
for (n, g, ru), a in ag.items():
    out.append("| %d | %d | %s | %d | %d | %d | %d |" % (n, g, ru, a["with_ruv"], a["k_at_cap"], a["resid_df_le1"], a["resid_df_0"]))
open(R + "m2-pooled.md", "w").write("\n".join(out) + "\n")

# ---- M05: k distribution by (n, g, rule); chosen percentage; realised missing share per cell ----
def T(x): return x == "TRUE"
kd = collections.defaultdict(collections.Counter)
pc = collections.defaultdict(collections.Counter)
for r in rows:
    for rule in ("p", "c1", "c2"):
        if T(r[rule + "_ruv"]):
            kd[(int(r["n"]), int(r["g"]), rule.upper())][int(float(r[rule + "_k_used"]))] += 1
            pc[(int(r["n"]), int(r["g"]), rule.upper())][int(float(r[rule + "_percentage"]))] += 1
with open(R + "m2-k-dist.tsv", "w") as f:
    f.write("n\tg\trule\tk_used\tdata_sets\n")
    for key in sorted(kd):
        for k in sorted(kd[key]): f.write("%d\t%d\t%s\t%d\t%d\n" % (key + (k, kd[key][k])))
with open(R + "m2-percentage.tsv", "w") as f:
    f.write("n\tg\trule\tpercentage\tdata_sets\n")
    for key in sorted(pc):
        for k in sorted(pc[key]): f.write("%d\t%d\t%s\t%d\t%d\n" % (key + (k, pc[key][k])))
# P's winning sweep k (before the cap), by n
ks = collections.defaultdict(list)
for r in rows:
    if T(r["p_ruv"]): ks[int(r["n"])].append(int(float(r["p_k_sweep"])))
with open(R + "m2-k-sweep.tsv", "w") as f:
    f.write("n\tdata_sets_with_RUV\tmin\tmedian\tmax\n")
    for n in sorted(ks):
        v = sorted(ks[n]); f.write("%d\t%d\t%d\t%d\t%d\n" % (n, len(v), v[0], v[len(v) // 2], v[-1]))
ms = collections.defaultdict(list)
for r in rows: ms[int(r["cell"])].append((float(r["realised_missing_share"]), float(r["expected_missing_share"]), float(r["complete_features"])))
head = {int(r["cell"]): r for r in rows}
with open(R + "m2-missing.tsv", "w") as f:
    f.write("cell\tn\tg\tf\tmissing_id\trealised_mean\trealised_min\trealised_max\texpected_mean\tcomplete_features_mean\n")
    for c in sorted(ms):
        v = ms[c]; h = head[c]
        f.write("%d\t%s\t%s\t%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.2f\n" % (c, h["n"], h["g"], h["f"], h["missing_id"], sum(x[0] for x in v) / 40,
                min(x[0] for x in v), max(x[0] for x in v), sum(x[1] for x in v) / 40, sum(x[2] for x in v) / 40))

# ---- speed per job relative to L2 (seconds_sweep by n and F) ----
l2 = {}
for r in csv.DictReader(open("l2/s1/timing.tsv"), delimiter="\t"): l2[(int(r["n"]), int(r["f"]))] = float(r["seconds_sweep"])
with open(R + "m2-speed.tsv", "w") as f:
    f.write("job\tmin_ratio_to_L2\tmax_ratio_to_L2\n")
    for j in ("j1", "j2", "j3", "j4", "cont"):
        acc = collections.defaultdict(list)
        for r in csv.DictReader(open("m2/%s/summary.tsv" % j), delimiter="\t"): acc[(int(r["n"]), int(r["f"]))].append(float(r["seconds_sweep"]))
        ratios = [sum(v) / len(v) / l2[k] for k, v in acc.items()]
        f.write("%s\t%.2f\t%.2f\n" % (j, min(ratios), max(ratios)))

# ---- M03: the exact-binomial counterpart at the same nominal level, and the rule's operating characteristic ----
alpha = 0.05 / 180
def tail(x, p, n=40): return sum(comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(x, n + 1))
x_exact = next(x for x in range(41) if tail(x, 0.05) <= alpha)
gain = list(csv.DictReader(open(R + "m2-gain.tsv"), delimiter="\t"))
cnt = {k: sum(1 for r in gain if int(r[c]) >= x_exact) for k, c in (("a", "c1_not_p"), ("b", "c2_not_p"), ("e", "c2_not_c1"))}
wil = {k: sum(1 for r in gain if r[k] == "TRUE") for k in "abe"}
size = tail(7, 0.05)
lines = ["# Exact-binomial counterpart of the declared gain rule (description only; the declared rule and its verdict stand)", "",
         "Nominal level per cell: 0.05 / 180 = %.6f. Smallest count of 40 whose upper tail P(X >= x | p = 0.05) is at most that: **%d** (P = %.6f)." % (alpha, x_exact, tail(x_exact, 0.05)),
         "Declared Wilson rule (>= 7 of 40): a %d, b %d, e %d cells. Exact test (>= %d of 40): a %d, b %d, e %d cells." % (wil["a"], wil["b"], wil["e"], x_exact, cnt["a"], cnt["b"], cnt["e"]), "",
         "Operating characteristic of the declared rule (score test at R = 40): per-cell size at a true share of 0.05 is P(X >= 7) = %.4f against a nominal %.6f;" % (size, alpha),
         "family-wise over 180 cells at a true share of 0.05 everywhere: 1 - (1 - %.4f)^180 = %.2f (the ticket's own figure 0.46); at 0.03: %.3f; at 0.02: %.4f." % (size, 1 - (1 - size) ** 180, 1 - (1 - tail(7, 0.03)) ** 180, 1 - (1 - tail(7, 0.02)) ** 180),
         "The rule therefore controls error near a true share of 0.03, not 0.05. The extra gain of C2 over C1 (e) rests on one cell (82) at 8 of 40.", "",
         "Cells with a under the exact test: " + ", ".join(r["cell"] for r in gain if int(r["c1_not_p"]) >= x_exact)]
open(R + "m2-exact.md", "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
