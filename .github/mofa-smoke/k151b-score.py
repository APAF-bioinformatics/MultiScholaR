"""KIT-151 selection evidence, run b: assemble the jobs' tables and compute the declared paired comparisons.

Usage: python3 -I k151b-score.py <jobs dir> <k151b-params.json> <follow-up followup-rows.tsv> <out dir>
Declaration: docs/qa/workflow-execution/kit/kit-151/selection-b-declaration.md. Standard library only.
Analysis (as declared): per scenario, arm and planted factor (F1, F2, F3 and "worst" = the smallest of the three in
that seed), the matched |r| over seeds; paired differences between declared arm pairs over the same seeds, with the
mean difference and a 95% seed-bootstrap interval (20,000 resamples, random.Random(151)); no pass/fail per seed.
Reading of one comparison: interval inside [-0.02, +0.02] -> "no practical difference"; else lower bound > 0 ->
"better"; upper bound < 0 -> "worse"; otherwise "inconclusive". A failed or missing cell drops that seed from the
comparisons it is in and is counted.
"""
import csv
import json
import math
import pathlib
import random
import statistics
import sys

jobs_dir, params_path, followup_path, out_dir = map(pathlib.Path, sys.argv[1:5])
out_dir.mkdir(parents=True, exist_ok=True)
params = json.loads(params_path.read_text())
SEEDS = params["seeds"]
FACTORS = ["F1", "F2", "F3", "worst"]
VIEWS = ["D1", "D2", "D3", "D4"]
DELTA = 0.02
B = 20000
BASE = ["none", "scale_views", "raw250", "raw1000", "raw2000", "random1000", "trend1000", "trendratio1000",
        "raw1000_z"]
OBS = [a + "_obs50" for a in ["raw250", "raw1000", "raw2000", "random1000", "trend1000", "trendratio1000",
                              "raw1000_z"]]
SCEN_ARMS = {"prev": BASE, "trend": BASE, "full": BASE + OBS}
if params.get("cells", "all") != "all":
    SCEN_ARMS = {}
    for c in params["cells"].split(","):
        sc, a = c.split(":")
        SCEN_ARMS.setdefault(sc, []).append(a)
COMPARISONS = [
    ("selecting at all", "scale_views", "raw250"), ("selecting at all", "scale_views", "raw1000"),
    ("selecting at all", "scale_views", "raw2000"), ("selecting at all", "scale_views", "random1000"),
    ("selecting at all", "scale_views", "trend1000"), ("selecting at all", "scale_views", "trendratio1000"),
    ("statistic", "raw1000", "random1000"), ("statistic", "raw1000", "trend1000"),
    ("statistic", "raw1000", "trendratio1000"), ("statistic", "random1000", "trend1000"),
    ("statistic", "random1000", "trendratio1000"), ("statistic", "trend1000", "trendratio1000"),
    ("number", "raw1000", "raw250"), ("number", "raw1000", "raw2000"),
    ("z-score", "raw1000", "raw1000_z"), ("z-score", "raw1000_obs50", "raw1000_z_obs50"),
    ("scale_views", "none", "scale_views"),
] + [("eligibility", a[:-6], a) for a in OBS]

rows, comps, jobs = [], [], []
for d in sorted(p for p in jobs_dir.iterdir() if p.is_dir()):
    for name, sink in (("rows.tsv", rows), ("composition.tsv", comps)):
        if (d / name).exists():
            with open(d / name, newline="") as f:
                sink.extend(csv.DictReader(f, delimiter="\t"))
    if (d / "job.json").exists():
        j = json.loads((d / "job.json").read_text())
        jobs.append({"job": d.name, "seed_index": j.get("seed_index"), "status": j.get("status"),
                     "failed_cells": j.get("failed_cells"), "seconds": j.get("seconds_total"),
                     "nproc": j.get("nproc"), "cpu_model": j.get("cpu_model"), "maxiter": j.get("maxiter"),
                     "missing_share": json.dumps(j.get("missing_share"), sort_keys=True),
                     "versions": json.dumps(j.get("versions"), sort_keys=True), "error": j.get("error", "")})


def num(x):
    try:
        v = float(x)
        return v if v == v else None
    except (TypeError, ValueError):
        return None


ok = {}
cellinfo = {}
failed_cells = []
for r in rows:
    key = (r["scenario"], int(r["seed_index"]), r["arm"])
    if r["status"] != "ok":
        if r["planted"] == "F1":
            failed_cells.append({"scenario": r["scenario"], "seed_index": r["seed_index"], "arm": r["arm"],
                                 "error": r["error"]})
        continue
    ok.setdefault(key, {})[r["planted"]] = num(r["abs_r"])
    cellinfo[key] = {"iterations": int(r["iterations"]), "stopped_by": r["stopped_by"],
                     "seconds": num(r["seconds"])}
    ok[key].setdefault("_r2", {})[r["planted"]] = {v: num(r["r2_" + v]) for v in VIEWS}
for key, v in ok.items():
    if all(v.get(f) is not None for f in ("F1", "F2", "F3")):
        v["worst"] = min(v["F1"], v["F2"], v["F3"])
expected = [(sc, s, a) for sc, arms in SCEN_ARMS.items() for s in SEEDS for a in arms]
missing = [k for k in expected if k not in ok]


def q(xs, p):
    xs = sorted(xs)
    pos = (len(xs) - 1) * p
    lo, hi = math.floor(pos), math.ceil(pos)
    return xs[lo] + (xs[hi] - xs[lo]) * (pos - lo)


def r4(x):
    return None if x is None else round(x, 4)


summary = []
for sc, arms in SCEN_ARMS.items():
    for a in arms:
        info = [cellinfo[(sc, s, a)] for s in SEEDS if (sc, s, a) in cellinfo]
        for f in FACTORS:
            vals = [(s, ok[(sc, s, a)][f]) for s in SEEDS if (sc, s, a) in ok and ok[(sc, s, a)].get(f) is not None]
            if not vals:
                continue
            xs = [v for _, v in vals]
            worst_seed = min(vals, key=lambda t: t[1])[0]
            row = {"scenario": sc, "arm": a, "factor": f, "seeds": len(xs), "median": r4(statistics.median(xs)),
                   "mean": r4(statistics.fmean(xs)), "p10": r4(q(xs, 0.1)), "worst": r4(min(xs)),
                   "worst_seed": worst_seed, "seeds_ge_0_9_descriptive": sum(x >= 0.9 for x in xs),
                   "iterations_median": statistics.median([i["iterations"] for i in info]) if info else None,
                   "iterations_max": max([i["iterations"] for i in info]) if info else None,
                   "at_maxiter": sum(i["stopped_by"] == "maxiter" for i in info),
                   "seconds_median": r4(statistics.median([i["seconds"] for i in info])) if info else None}
            if f != "worst":
                for v in VIEWS:
                    r2 = [ok[(sc, s, a)]["_r2"][f][v] for s in SEEDS if (sc, s, a) in ok
                          and ok[(sc, s, a)]["_r2"][f][v] is not None]
                    row["r2_" + v + "_median"] = r4(statistics.median(r2)) if r2 else None
            summary.append(row)

rng = random.Random(151)
pairs = []
for group, a, b in COMPARISONS:
    for sc, arms in SCEN_ARMS.items():
        if a not in arms or b not in arms:
            continue
        for f in FACTORS:
            both = [s for s in SEEDS if (sc, s, a) in ok and (sc, s, b) in ok
                    and ok[(sc, s, a)].get(f) is not None and ok[(sc, s, b)].get(f) is not None]
            if len(both) < 2:
                continue
            d = [ok[(sc, s, b)][f] - ok[(sc, s, a)][f] for s in both]
            n = len(d)
            means = sorted(statistics.fmean(rng.choices(d, k=n)) for _ in range(B))
            lo, hi = means[int(0.025 * B)], means[int(0.975 * B) - 1]
            if -DELTA <= lo and hi <= DELTA:
                reading = "no practical difference"
            elif lo > 0:
                reading = "better"
            elif hi < 0:
                reading = "worse"
            else:
                reading = "inconclusive"
            pairs.append({"group": group, "scenario": sc, "a": a, "b": b, "factor": f, "seeds": n,
                          "dropped_seeds": len(SEEDS) - n, "mean_diff": r4(statistics.fmean(d)),
                          "ci_lo": r4(lo), "ci_hi": r4(hi), "median_diff": r4(statistics.median(d)),
                          "b_better": sum(x > 1e-12 for x in d), "b_worse": sum(x < -1e-12 for x in d),
                          "worst_seed_diff": r4(min(d)), "reading_b_vs_a": reading})

# Composition: share of signal-carrying features kept, share of kept that carry signal, kept noise sd, kept observed.
cidx = {}
for c in comps:
    cidx[(c["scenario"], int(c["seed_index"]), c["arm"], c["view"], c["carries"])] = c
comp_summary = []
for sc, arms in SCEN_ARMS.items():
    for a in arms:
        for v in VIEWS:
            rec, pur, f3, nsd, obs = [], [], [], [], []
            for s in SEEDS:
                al = cidx.get((sc, s, a, v, "all"))
                no = cidx.get((sc, s, a, v, "none"))
                if al is None or no is None:
                    continue
                sig_t = int(al["total"]) - int(no["total"])
                sig_k = int(al["kept"]) - int(no["kept"])
                if sig_t > 0:
                    rec.append(sig_k / sig_t)
                if int(al["kept"]) > 0:
                    pur.append(sig_k / int(al["kept"]))
                c3 = cidx.get((sc, s, a, v, "F3"))
                if c3 is not None and int(c3["total"]) > 0:
                    f3.append(int(c3["kept"]) / int(c3["total"]))
                if num(al["kept_noise_sd_median"]) is not None:
                    nsd.append(num(al["kept_noise_sd_median"]))
                if num(al["kept_observed_median"]) is not None:
                    obs.append(num(al["kept_observed_median"]))
            med = lambda xs: r4(statistics.median(xs)) if xs else None
            comp_summary.append({"scenario": sc, "arm": a, "view": v, "seeds": len(rec),
                                 "signal_kept_share": med(rec), "kept_signal_share": med(pur),
                                 "f3_kept_share": med(f3), "kept_noise_sd_median": med(nsd),
                                 "kept_observed_fraction_median": med(obs)})

# Tie to the follow-up run: scenario prev, seeds 1..10, against n20strong slow (maxiter 1000 there).
tie = []
ARM_MAP = {"none": "none", "scale_views": "scale_views", "raw250": "cap250", "raw1000": "cap1000",
           "raw2000": "cap2000"}
if followup_path.exists():
    old = {}
    with open(followup_path, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["scenario"] == "n20strong" and r["mode"] == "slow" and r["status"] == "ok":
                old[(int(r["seed_index"]), r["arm"], r["planted"])] = (num(r["abs_r"]), int(r["iterations"]),
                                                                        r["stopped_by"])
    for s in [s for s in SEEDS if s <= 10]:
        for a, oa in ARM_MAP.items():
            for f in ("F1", "F2", "F3"):
                if ("prev", s, a) in ok and (s, oa, f) in old:
                    o = old[(s, oa, f)]
                    tie.append({"seed_index": s, "arm": a, "planted": f, "old_abs_r": r4(o[0]),
                                "new_abs_r": r4(ok[("prev", s, a)][f]),
                                "abs_diff": abs(ok[("prev", s, a)][f] - o[0]), "old_iterations": o[1],
                                "new_iterations": cellinfo[("prev", s, a)]["iterations"], "old_stopped_by": o[2]})


def write_tsv(name, data):
    if not data:
        (out_dir / name).write_text("")
        return
    keys = list(data[0].keys())
    for d in data:
        for k in d:
            if k not in keys:
                keys.append(k)
    with open(out_dir / name, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(data)


def raw_tsv(name, data):
    if not data:
        (out_dir / name).write_text("")
        return
    with open(out_dir / name, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(data[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(data)


raw_tsv("selection-b-rows.tsv", rows)
raw_tsv("selection-b-composition.tsv", comps)
write_tsv("selection-b-jobs.tsv", jobs)
write_tsv("selection-b-summary.tsv", summary)
write_tsv("selection-b-pairs.tsv", pairs)
write_tsv("selection-b-composition-summary.tsv", comp_summary)
write_tsv("selection-b-tie.tsv", tie)
write_tsv("selection-b-failed.tsv", failed_cells)

maxiter_hits = sum(1 for k, i in cellinfo.items() if i["stopped_by"] == "maxiter")
cpu = {}
for j in jobs:
    cpu[j["cpu_model"]] = cpu.get(j["cpu_model"], 0) + 1
seconds = sum(j["seconds"] or 0 for j in jobs)
score = {"run": params.get("run"), "seeds": SEEDS, "expected_cells": len(expected), "ok_cells": len(ok),
         "failed_cells": len(failed_cells), "missing_cells": len(missing), "jobs": len(jobs),
         "cells_at_maxiter": maxiter_hits, "maxiter": params.get("maxiter"), "cpu_models": cpu,
         "training_seconds_total": round(seconds, 1)}
if tie:
    conv = [t for t in tie if t["old_stopped_by"] == "convergence"]
    score["tie"] = {"compared": len(tie), "identical_1e-9": sum(t["abs_diff"] <= 1e-9 for t in tie),
                    "within_1e-6": sum(t["abs_diff"] <= 1e-6 for t in tie),
                    "old_stopped_by_convergence": len(conv),
                    "of_those_within_1e-6": sum(t["abs_diff"] <= 1e-6 for t in conv),
                    "max_abs_diff_old_converged": max((t["abs_diff"] for t in conv), default=None)}
(out_dir / "selection-b-score.json").write_text(json.dumps(score, indent=2, sort_keys=True) + "\n")

md = ["# KIT-151 selection evidence, run b: score", "",
      f"Run {params.get('run')}; seeds {SEEDS[0]}..{SEEDS[-1]} ({len(SEEDS)}); cells ok {len(ok)} of {len(expected)};"
      f" failed {len(failed_cells)}; missing {len(missing)}; cells stopped at the round limit ({params.get('maxiter')}):"
      f" {maxiter_hits}. CPU models: {cpu}.", ""]
for sc, arms in SCEN_ARMS.items():
    md += [f"## {sc}: matched |r|, median / 10th percentile / worst seed", "",
           "| arm | F1 | F2 | F3 | worst of three | iterations median (max); at limit |", "|---|---|---|---|---|---|"]
    for a in arms:
        cells = []
        it = ""
        for f in FACTORS:
            s = next((x for x in summary if x["scenario"] == sc and x["arm"] == a and x["factor"] == f), None)
            cells.append("-" if s is None else f"{s['median']:.3f} / {s['p10']:.3f} / {s['worst']:.3f}")
            if s is not None:
                it = f"{s['iterations_median']:.0f} ({s['iterations_max']}); {s['at_maxiter']}"
        md.append(f"| {a} | " + " | ".join(cells) + f" | {it} |")
    md += ["", f"### {sc}: paired differences (b minus a), worst-of-three per seed: mean [95% CI] reading", "",
           "| group | a | b | mean [CI] | b better / worse | worst seed diff | reading |", "|---|---|---|---|---|---|---|"]
    for p in pairs:
        if p["scenario"] == sc and p["factor"] == "worst":
            md.append(f"| {p['group']} | {p['a']} | {p['b']} | {p['mean_diff']:+.3f} [{p['ci_lo']:+.3f}, "
                      f"{p['ci_hi']:+.3f}] | {p['b_better']} / {p['b_worse']} | {p['worst_seed_diff']:+.3f} | "
                      f"{p['reading_b_vs_a']} |")
    md.append("")
if tie:
    md += ["## Tie to the follow-up (prev, seeds 1..10, against n20strong slow)", "", "```",
           json.dumps(score["tie"], indent=2), "```", ""]
(out_dir / "selection-b-score.md").write_text("\n".join(md) + "\n")
print("\n".join(md))
