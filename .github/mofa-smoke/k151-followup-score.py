"""KIT-151 cap evidence, follow-up: assemble the jobs' rows and apply the pre-declared readings (a) to (d).

Usage: python3 -I k151-followup-score.py <jobs dir> <params.json> <first-run results-rows.tsv> <out dir>
Readings (docs/qa/workflow-execution/kit/kit-151/cap-followup-declaration.md):
 recovered = greedy-matched |r| >= 0.9; marks need >= 9 of 10 seeds; a failed or missing cell is not recovered.
 (a) early stopping: on n20weak under "slow", arms none and scale_views each reach F1 >= 9 and F3 >= 9.
 (b) convergence default: the fastest mode (fast, then medium) such that for every capped arm (cap250 .. cap5000),
     every scenario in which both it and "slow" ran, and every planted factor, seeds(mode) >= seeds(slow) - 1;
     otherwise "slow".
 (c) cap: under mode (b), per scenario the caps meeting F2 >= 9 and F3 >= 9; a scenario where some cap meets it
     "discriminates"; the recommended cap is the largest cap meeting the rule in every discriminating scenario;
     none in common: "no single cap"; no discriminating scenario: "no cap meets the rule anywhere". If n20weak does
     not discriminate but another scenario does, the limit is detectability at N = 20 with the weak F2 (stated as a
     limitation) and the cap comes from the discriminating scenarios.
 (d) weight_views: under mode (b), the same two marks per scenario; "meets" only if it meets both in every
     discriminating scenario.
"""
import csv
import json
import pathlib
import statistics
import sys

jobs_dir, params_path, first_path, out_dir = map(pathlib.Path, sys.argv[1:5])
out_dir.mkdir(parents=True, exist_ok=True)
params = json.loads(params_path.read_text())
THRESHOLD, NEED = 0.9, 9
MODES = ["fast", "medium", "slow"]
CAPS = [250, 500, 1000, 2000, 5000]
MAIN_ARMS = ["none", "scale_views", "cap250", "cap500", "cap1000", "cap2000", "cap5000", "weight_views", "zscore"]
NOISY_ARMS = ["noisy_cap1000", "noisy_cap2000"]
VIEWS = ["D1", "D2", "D3", "D4"]
scenarios = params["scenarios"]
seeds = list(range(1, params["seeds"] + 1))


def modes_for(sc):
    return params["modes_n60"] if sc.startswith("n60") else params["modes_n20"]


def arms_for(sc):
    if params.get("cells", "all") != "all":
        return sorted({c.split("@")[0] for c in params["cells"].split(",")})
    return MAIN_ARMS + (NOISY_ARMS if sc.startswith("n20") else [])


rows, comps, jobs = [], [], []
for d in sorted(p for p in jobs_dir.iterdir() if p.is_dir()):
    for name, sink in (("rows.tsv", rows), ("composition.tsv", comps)):
        if (d / name).exists():
            with open(d / name, newline="") as f:
                sink.extend(csv.DictReader(f, delimiter="\t"))
    if (d / "job.json").exists():
        j = json.loads((d / "job.json").read_text())
        jobs.append({"job": d.name, "scenario": j.get("scenario"), "seed_index": j.get("seed_index"),
                     "status": j.get("status"), "failed_cells": j.get("failed_cells"),
                     "seconds": j.get("seconds_total"), "nproc": j.get("nproc"), "cpu_model": j.get("cpu_model"),
                     "versions": json.dumps(j.get("versions"), sort_keys=True), "error": j.get("error", "")})
index = {(r["scenario"], int(r["seed_index"]), r["arm"], r["mode"], r["planted"]): r for r in rows}


def num(x):
    try:
        v = float(x)
        return v if v == v else None
    except (TypeError, ValueError):
        return None


def med(values):
    values = [v for v in values if v is not None]
    return round(statistics.median(values), 4) if values else None


summary, hits = [], {}
for sc in scenarios:
    for mode in modes_for(sc):
        for arm in arms_for(sc):
            for p in ["F1", "F2", "F3"]:
                got = [index.get((sc, s, arm, mode, p)) for s in seeds]
                ok = [g for g in got if g is not None and g["status"] == "ok"]
                absr = [num(g["abs_r"]) for g in ok]
                h = sum(1 for v in absr if v is not None and v >= THRESHOLD)
                hits[(sc, mode, arm, p)] = h
                row = {"scenario": sc, "mode": mode, "arm": arm, "planted": p, "seeds_ok": len(ok),
                       "seeds_abs_r_ge_0.9": h, "median_abs_r": med(absr),
                       "min_abs_r": round(min(absr), 4) if absr else None}
                for v in VIEWS:
                    row["median_r2_" + v] = med([num(g["r2_" + v]) for g in ok])
                for v in VIEWS:
                    row["median_kept_" + v] = med([num(g["kept_" + v]) for g in ok])
                row["median_iterations"] = med([num(g["iterations"]) for g in ok])
                row["seeds_stopped_by_maxiter"] = sum(1 for g in ok if g["stopped_by"] == "maxiter")
                row["median_elbo_final"] = med([num(g["elbo_final"]) for g in ok])
                row["seeds_with_dropped_factors"] = sum(1 for g in ok if (num(g["factors_dropped"]) or 0) > 0)
                summary.append(row)


def H(sc, mode, arm, p):
    return hits.get((sc, mode, arm, p), 0)


readings = {}
# (a)
a = {arm: {"F1": H("n20weak", "slow", arm, "F1"), "F3": H("n20weak", "slow", arm, "F3")} for arm in
     ["none", "scale_views"]}
for arm in a:
    a[arm]["reaches"] = a[arm]["F1"] >= NEED and a[arm]["F3"] >= NEED
readings["a_early_stopping"] = {"per_arm": a, "verdict": all(v["reaches"] for v in a.values())
                                if "n20weak" in scenarios and "slow" in modes_for("n20weak") else None}
# (b)
capped = [f"cap{c}" for c in CAPS]
b = {}
for mode in ["fast", "medium"]:
    checked, failures = 0, []
    for sc in scenarios:
        if mode not in modes_for(sc) or "slow" not in modes_for(sc):
            continue
        for arm in capped:
            for p in ["F1", "F2", "F3"]:
                checked += 1
                if H(sc, mode, arm, p) < H(sc, "slow", arm, p) - 1:
                    failures.append(f"{sc} {arm} {p}: {H(sc, mode, arm, p)} vs slow {H(sc, 'slow', arm, p)}")
    b[mode] = {"comparisons": checked, "passes": checked > 0 and not failures, "failures": failures}
recommended_mode = next((m for m in ["fast", "medium"] if b[m]["passes"]), "slow")
readings["b_convergence"] = {"per_mode": b, "recommended_mode": recommended_mode}
# (c)
c = {}
for sc in scenarios:
    mode = recommended_mode if recommended_mode in modes_for(sc) else "slow"
    meeting = [cap for cap in CAPS if H(sc, mode, f"cap{cap}", "F2") >= NEED and H(sc, mode, f"cap{cap}", "F3") >= NEED]
    c[sc] = {"mode_used": mode, "meeting_caps": meeting, "largest_meeting": max(meeting) if meeting else None,
             "per_cap": {cap: {"F2": H(sc, mode, f"cap{cap}", "F2"), "F3": H(sc, mode, f"cap{cap}", "F3")}
                         for cap in CAPS}}
discriminating = [sc for sc in scenarios if c[sc]["meeting_caps"]]
if not discriminating:
    cap_verdict = "no cap meets the rule anywhere"
else:
    common = set(CAPS)
    for sc in discriminating:
        common &= set(c[sc]["meeting_caps"])
    cap_verdict = max(common) if common else "no single cap"
readings["c_cap"] = {"per_scenario": c, "discriminating_scenarios": discriminating, "recommended_cap": cap_verdict,
                     "detectability_limit_at_n20_weak": ("n20weak" not in discriminating) and bool(discriminating)}
# (d)
d = {}
for sc in scenarios:
    mode = recommended_mode if recommended_mode in modes_for(sc) else "slow"
    f2, f3 = H(sc, mode, "weight_views", "F2"), H(sc, mode, "weight_views", "F3")
    d[sc] = {"mode_used": mode, "F2": f2, "F3": f3, "meets": f2 >= NEED and f3 >= NEED}
readings["d_weight_views"] = {"per_scenario": d, "meets_in_every_discriminating_scenario":
                              bool(discriminating) and all(d[sc]["meets"] for sc in discriminating)}

# Reproduction check against the first run (fast mode, weak F2: n20weak = main, n60weak = n60), not a gate.
first = {}
if first_path.exists():
    with open(first_path, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            first[({"main": "n20weak", "n60": "n60weak"}.get(r["scenario"]), int(r["seed_index"]), r["arm"],
                   r["planted"])] = num(r["abs_r"])
diffs = []
for (sc, s, arm, mode, p), r in index.items():
    if mode == "fast" and (sc, s, arm, p) in first and r["status"] == "ok":
        diffs.append(abs(num(r["abs_r"]) - first[(sc, s, arm, p)]))
readings["reproduction_of_first_run"] = {"compared": len(diffs), "max_abs_difference_abs_r": max(diffs) if diffs
                                         else None, "identical": sum(1 for x in diffs if x == 0)}

# Composition summary (median over seeds).
comp_summary = {}
for r in comps:
    key = (r["scenario"], r["arm"], r["view"], r["carries"])
    comp_summary.setdefault(key, []).append((num(r["total"]), num(r["kept"]), num(r["dropped"])))
comp_rows = [{"scenario": k[0], "arm": k[1], "view": k[2], "carries": k[3],
              "median_total": med([x[0] for x in v]), "median_kept": med([x[1] for x in v]),
              "median_dropped": med([x[2] for x in v]), "seeds": len(v)} for k, v in sorted(comp_summary.items())]

expected = sum(len(modes_for(sc)) * len(arms_for(sc)) * len(seeds) * 3 for sc in scenarios)


def write_tsv(name, data):
    if not data:
        return
    with open(out_dir / name, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(data[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(data)


write_tsv("followup-rows.tsv", sorted(rows, key=lambda r: (r["scenario"], r["mode"], r["arm"], int(r["seed_index"]),
                                                           r["planted"])))
write_tsv("followup-summary.tsv", summary)
write_tsv("followup-composition.tsv", comps)
write_tsv("followup-composition-summary.tsv", comp_rows)
write_tsv("followup-jobs.tsv", jobs)
score = {"params": params, "expected_rows": expected, "rows": len(rows),
         "failed_rows": sum(1 for r in rows if r["status"] != "ok"), "jobs": len(jobs), "readings": readings}
(out_dir / "followup-score.json").write_text(json.dumps(score, indent=2))

md = ["# KIT-151 cap follow-up: scored by k151-followup-score.py", "",
      f"Rows {len(rows)} of {expected} expected; failed rows {score['failed_rows']}; jobs {len(jobs)}.", ""]
for sc in scenarios:
    for mode in modes_for(sc):
        md += [f"## {sc}, {mode}", "", "| arm | F1 | F2 | F3 | median abs r F1/F2/F3 | median iterations | maxiter seeds |",
               "|---|---|---|---|---|---|---|"]
        for arm in arms_for(sc):
            s = {p: next(x for x in summary if x["scenario"] == sc and x["mode"] == mode and x["arm"] == arm
                         and x["planted"] == p) for p in ["F1", "F2", "F3"]}
            md.append(f"| {arm} | {s['F1']['seeds_abs_r_ge_0.9']} | {s['F2']['seeds_abs_r_ge_0.9']} | "
                      f"{s['F3']['seeds_abs_r_ge_0.9']} | {s['F1']['median_abs_r']}/{s['F2']['median_abs_r']}/"
                      f"{s['F3']['median_abs_r']} | {s['F1']['median_iterations']} | "
                      f"{s['F1']['seeds_stopped_by_maxiter']} |")
        md.append("")
md += ["## Readings", "", "```", json.dumps(readings, indent=1), "```"]
(out_dir / "followup-score.md").write_text("\n".join(md) + "\n")
print("\n".join(md))
