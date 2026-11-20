"""KIT-151 SCI-01 re-score (no compute): paired continuous comparisons of the existing follow-up rows.

Usage: python3 -I k151-rescore.py <followup-rows.tsv> <followup-composition.tsv> <out dir>
Reads the committed per-seed rows of run 44654405980; writes rescore-pairs.tsv, rescore-arms.tsv,
rescore-composition.tsv and rescore.md. Standard library only. Bootstrap: 20,000 seed resamples, fixed RNG seed.
"""
import csv
import math
import random
import statistics
import sys
from pathlib import Path

rows_path, comp_path, out_dir = sys.argv[1], sys.argv[2], Path(sys.argv[3])
out_dir.mkdir(parents=True, exist_ok=True)
ARMS = ["none", "scale_views", "cap250", "cap500", "cap1000", "cap2000", "cap5000"]
LABEL = {"none": "no selection (scale_views off)", "scale_views": "no selection (scale_views on)",
         "cap250": "250", "cap500": "500", "cap1000": "1,000", "cap2000": "2,000", "cap5000": "5,000"}
SCEN = ["n20strong", "n20weak"]
MODES = ["slow", "medium", "fast"]
FACT = ["F1", "F2", "F3"]
B = 20000

r = {}
iters = {}
with open(rows_path, newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        if row["scenario"] not in SCEN or row["arm"] not in ARMS:
            continue
        if row["status"] != "ok":
            raise SystemExit("failed row present: " + str(row))
        key = (row["scenario"], row["mode"], row["arm"], int(row["seed_index"]))
        r.setdefault(key, {})[row["planted"]] = float(row["abs_r"])
        iters[key] = int(row["iterations"])


def vec(scen, mode, arm, f):
    return [r[(scen, mode, arm, s)][f] if f != "min23" else min(r[(scen, mode, arm, s)]["F2"],
            r[(scen, mode, arm, s)]["F3"]) for s in range(1, 11)]


def quantile(xs, q):
    xs = sorted(xs)
    pos = (len(xs) - 1) * q
    lo, hi = math.floor(pos), math.ceil(pos)
    return xs[lo] + (xs[hi] - xs[lo]) * (pos - lo)


def boot_ci(d, rng):
    n = len(d)
    means = sorted(statistics.fmean(rng.choices(d, k=n)) for _ in range(B))
    return means[int(0.025 * B)], means[int(0.975 * B) - 1]


def sign_p(wins, losses):
    n = wins + losses
    if n == 0:
        return 1.0
    k = min(wins, losses)
    p = sum(math.comb(n, i) for i in range(0, k + 1)) / 2 ** n * 2
    return min(1.0, p)


rng = random.Random(151)
arm_rows = []
pair_rows = []
for scen in SCEN:
    for mode in MODES:
        for arm in ARMS:
            for f in FACT + ["min23"]:
                v = vec(scen, mode, arm, f)
                it = [iters[(scen, mode, arm, s)] for s in range(1, 11)]
                arm_rows.append(dict(scenario=scen, mode=mode, arm=arm, factor=f, median=statistics.median(v),
                                     mean=statistics.fmean(v), p10=quantile(v, 0.1), worst=min(v),
                                     worst_seed=v.index(min(v)) + 1, n_ge_0_9=sum(x >= 0.9 for x in v),
                                     iter_median=statistics.median(it), at_maxiter=sum(i >= 1000 for i in it)))
        caps = ARMS
        for i, a in enumerate(caps):
            for b in caps[i + 1:]:
                for f in FACT + ["min23"]:
                    va, vb = vec(scen, mode, a, f), vec(scen, mode, b, f)
                    d = [y - x for x, y in zip(va, vb)]  # b minus a
                    lo, hi = boot_ci(d, rng)
                    wins = sum(x > 1e-12 for x in d)
                    losses = sum(x < -1e-12 for x in d)
                    pair_rows.append(dict(scenario=scen, mode=mode, a=a, b=b, factor=f,
                                          mean_diff=statistics.fmean(d), median_diff=statistics.median(d),
                                          ci_lo=lo, ci_hi=hi, b_better=wins, a_better=losses,
                                          sign_p=sign_p(wins, losses), worst_diff=min(d)))


def write(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        for row in rows:
            w.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in row.items()})


write(out_dir / "rescore-arms.tsv", arm_rows)
write(out_dir / "rescore-pairs.tsv", pair_rows)

# Composition: share of D1 / D2 signal-carrying features kept (recall) and share of kept that carry signal (purity).
comp = {}
with open(comp_path, newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        if row["scenario"] not in SCEN:
            continue
        key = (row["scenario"], row["arm"], row["view"], int(row["seed_index"]))
        comp.setdefault(key, {})[row["carries"]] = (int(row["total"]), int(row["kept"]))
comp_rows = []
for scen in SCEN:
    for arm in ["cap250", "cap500", "cap1000", "cap2000", "cap5000"]:
        for view in ["D1", "D2"]:
            rec, pur, f3rec = [], [], []
            for s in range(1, 11):
                c = comp[(scen, arm, view, s)]
                total_all = sum(c[k][0] for k in ("none",)) + 0
                none_t, none_k = c["none"]
                d_total = {"D1": 20000, "D2": 2000}[view]
                kept_all = min(d_total, int(arm[3:]))
                sig_t = d_total - none_t
                sig_k = kept_all - none_k
                rec.append(sig_k / sig_t)
                pur.append(sig_k / kept_all)
                if view == "D1":
                    f3rec.append(c["F3"][1] / c["F3"][0])
            comp_rows.append(dict(scenario=scen, arm=arm, view=view, signal_kept_share_median=statistics.median(rec),
                                  signal_kept_share_min=min(rec), kept_signal_purity_median=statistics.median(pur),
                                  f3_kept_share_median=statistics.median(f3rec) if f3rec else float("nan")))
write(out_dir / "rescore-composition.tsv", comp_rows)

# Markdown summary: slow mode, both scenarios.
lines = ["# SCI-01 re-score of run 44654405980 (N = 20)", ""]
for scen in SCEN:
    for mode in ["slow"]:
        lines += [f"## {scen}, {mode}: matched |r| per arm, median / 10th percentile / worst seed (seeds >= 0.9)", "",
                  "| arm | F1 | F2 | F3 | min(F2,F3) per seed | iterations median; at limit |", "|---|---|---|---|---|---|"]
        for arm in ARMS:
            cells = []
            for f in FACT + ["min23"]:
                a = next(x for x in arm_rows if x["scenario"] == scen and x["mode"] == mode and x["arm"] == arm
                         and x["factor"] == f)
                cells.append(f"{a['median']:.3f} / {a['p10']:.3f} / {a['worst']:.3f} ({a['n_ge_0_9']})")
            lines.append(f"| {LABEL[arm]} | " + " | ".join(cells) + f" | {a['iter_median']:.0f}; {a['at_maxiter']} |")
        lines += ["", f"### {scen}, {mode}: paired differences (column arm minus row arm), mean [95% seed-bootstrap]; "
                  "seeds better/worse; sign-test p", ""]
        for f in ["F2", "F3", "min23"]:
            lines += [f"Factor {f}:", "", "| comparison | mean diff [95% CI] | median diff | better / worse | sign p | "
                      "worst seed diff |", "|---|---|---|---|---|---|"]
            for a, b in [("cap250", "cap500"), ("cap250", "cap1000"), ("cap500", "cap1000"), ("cap1000", "cap2000"),
                         ("cap1000", "cap5000"), ("scale_views", "cap1000"), ("scale_views", "cap250"),
                         ("none", "scale_views")]:
                p = next(x for x in pair_rows if x["scenario"] == scen and x["mode"] == mode and x["a"] == a
                         and x["b"] == b and x["factor"] == f)
                lines.append(f"| {LABEL[b]} vs {LABEL[a]} | {p['mean_diff']:+.3f} [{p['ci_lo']:+.3f}, {p['ci_hi']:+.3f}]"
                             f" | {p['median_diff']:+.3f} | {p['b_better']} / {p['a_better']} | {p['sign_p']:.3f} | "
                             f"{p['worst_diff']:+.3f} |")
            lines.append("")
lines += ["## Signal-carrying features kept (median over seeds)", "",
          "| scenario | cap | D1 signal kept (share of loaded) | D1 kept that are loaded | D1 F3-loaded kept | "
          "D2 signal kept | D2 kept that are loaded |", "|---|---|---|---|---|---|---|"]
for scen in SCEN:
    for arm in ["cap250", "cap500", "cap1000", "cap2000", "cap5000"]:
        d1 = next(x for x in comp_rows if x["scenario"] == scen and x["arm"] == arm and x["view"] == "D1")
        d2 = next(x for x in comp_rows if x["scenario"] == scen and x["arm"] == arm and x["view"] == "D2")
        lines.append(f"| {scen} | {LABEL[arm]} | {d1['signal_kept_share_median']:.3f} | "
                     f"{d1['kept_signal_purity_median']:.3f} | {d1['f3_kept_share_median']:.3f} | "
                     f"{d2['signal_kept_share_median']:.3f} | {d2['kept_signal_purity_median']:.3f} |")
(out_dir / "rescore.md").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
