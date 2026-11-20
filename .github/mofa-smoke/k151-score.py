"""KIT-151 cap evidence: assemble the jobs' rows and score them by the pre-declared rule.

Usage: python3 -I k151-score.py <jobs dir (one sub-folder per job artifact)> <params.json> <out dir>
Rule (docs/qa/workflow-execution/kit/kit-151/cap-evidence-declaration.md): on the main scenario, a factor is
recovered in a seed when its greedy-matched |r| >= 0.9; an arm "recovers F2" when F2 is recovered in >= 9 of 10
seeds and "keeps F3" when F3 is recovered in >= 9 of 10 seeds (losing F3 = |r| < 0.9 in more than 1 of 10 seeds).
A failed or missing arm-seed counts as not recovered. The recommendation is the largest cap in {1000, 2000, 5000,
10000} that recovers F2 and keeps F3; if none, "no cap meets the rule" and the closest cap (most F2 seeds, then most
F3 seeds, then the larger cap). weight_views is judged by the same two marks. Sensitivity scenarios are scored the
same way and reported, not used for the recommendation.
"""
import csv
import hashlib
import json
import pathlib
import statistics
import sys

jobs_dir, params_path, out_dir = map(pathlib.Path, sys.argv[1:4])
out_dir.mkdir(parents=True, exist_ok=True)
params = json.loads(params_path.read_text())
THRESHOLD, NEED = 0.9, 9
CAPS = [1000, 2000, 5000, 10000]
ALL_ARMS = ["none", "scale_views", "cap1000", "cap2000", "cap5000", "cap10000", "weight_views", "zscore"]
arms = ALL_ARMS if params["arms"] == "all" else params["arms"].split(",")
scenarios = params["scenarios"]
seeds = list(range(1, params["seeds"] + 1))
VIEWS = ["D1", "D2", "D3", "D4"]

rows, jobs, excerpts, checks = [], [], {}, {}
for d in sorted(p for p in jobs_dir.iterdir() if p.is_dir()):
    if (d / "rows.tsv").exists():
        with open(d / "rows.tsv", newline="") as f:
            rows.extend(csv.DictReader(f, delimiter="\t"))
    if (d / "job.json").exists():
        j = json.loads((d / "job.json").read_text())
        jobs.append({"job": d.name, "scenario": j.get("scenario"), "seed_index": j.get("seed_index"),
                     "status": j.get("status"), "failed_arms": j.get("failed_arms"), "seconds": j.get("seconds_total"),
                     "nproc": j.get("nproc"), "cpu_model": j.get("cpu_model"),
                     "versions": json.dumps(j.get("versions"), sort_keys=True), "error": j.get("error", "")})
    if (d / "source-excerpts.txt").exists():
        text = (d / "source-excerpts.txt").read_text()
        excerpts.setdefault(hashlib.sha256(text.encode()).hexdigest(), (d.name, text))
    if (d / "save-data-check.json").exists():
        checks[d.name] = json.loads((d / "save-data-check.json").read_text())

index = {(r["scenario"], int(r["seed_index"]), r["arm"], r["planted"]): r for r in rows}


def num(x):
    try:
        v = float(x)
        return v if v == v else None
    except (TypeError, ValueError):
        return None


def median(values):
    values = [v for v in values if v is not None]
    return round(statistics.median(values), 4) if values else None


summary = []
marks = {}
for sc in scenarios:
    for arm in arms:
        for planted in ["F1", "F2", "F3"]:
            got = [index.get((sc, s, arm, planted)) for s in seeds]
            ok = [g for g in got if g is not None and g["status"] == "ok"]
            absr = [num(g["abs_r"]) for g in ok]
            hits = sum(1 for v in absr if v is not None and v >= THRESHOLD)
            row = {"scenario": sc, "arm": arm, "planted": planted, "seeds": len(seeds), "seeds_ok": len(ok),
                   "seeds_failed_or_missing": len(seeds) - len(ok), "seeds_abs_r_ge_0.9": hits,
                   "median_abs_r": median(absr), "min_abs_r": round(min(absr), 4) if absr else None}
            for v in VIEWS:
                row["median_r2_" + v] = median([num(g["r2_" + v]) for g in ok])
            for v in VIEWS:
                row["median_kept_" + v] = median([num(g["kept_" + v]) for g in ok])
            row["median_iterations"] = median([num(g["iterations"]) for g in ok])
            row["seeds_with_dropped_factors"] = sum(1 for g in ok if (num(g["factors_dropped"]) or 0) > 0)
            summary.append(row)
            marks[(sc, arm, planted)] = hits


def judge(sc, arm):
    f2, f3 = marks.get((sc, arm, "F2"), 0), marks.get((sc, arm, "F3"), 0)
    return {"F2_seeds": f2, "F3_seeds": f3, "recovers_F2": f2 >= NEED, "keeps_F3": f3 >= NEED,
            "meets_rule": f2 >= NEED and f3 >= NEED}


verdicts = {sc: {arm: judge(sc, arm) for arm in arms} for sc in scenarios}
decision = {}
for sc in scenarios:
    caps = [c for c in CAPS if f"cap{c}" in arms]
    meeting = [c for c in caps if verdicts[sc][f"cap{c}"]["meets_rule"]]
    if meeting:
        decision[sc] = {"recommendation": max(meeting), "meeting_caps": meeting}
    elif caps:
        closest = max(caps, key=lambda c: (verdicts[sc][f"cap{c}"]["F2_seeds"], verdicts[sc][f"cap{c}"]["F3_seeds"], c))
        decision[sc] = {"recommendation": "no cap meets the rule", "closest_cap": closest,
                        "closest": verdicts[sc][f"cap{closest}"]}

fields = list(rows[0].keys()) if rows else []
with open(out_dir / "results-rows.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(sorted(rows, key=lambda r: (r["scenario"], r["arm"], int(r["seed_index"]), r["planted"])))
with open(out_dir / "results-summary.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(summary[0].keys()), delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(summary)
with open(out_dir / "jobs.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["job", "scenario", "seed_index", "status", "failed_arms", "seconds", "nproc",
                                      "cpu_model", "versions", "error"], delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(jobs)
for i, (sha, (name, text)) in enumerate(sorted(excerpts.items(), key=lambda kv: kv[1][0])):
    (out_dir / f"source-excerpts-{i + 1}.txt").write_text(f"# from job {name}; sha256 {sha}\n" + text)
first_check = next(iter(sorted(checks.items())), (None, None))
(out_dir / "save-data-check.json").write_text(json.dumps({"from_job": first_check[0], "check": first_check[1],
    "jobs_agreeing_on_ok": sorted({(c.get("save_data_false") or {}).get("variance_explained_ok") is True
                                   for c in checks.values()})}, indent=2, default=str))
expected_rows = len(scenarios) * len(seeds) * len(arms) * 3
score = {"params": params, "expected_rows": expected_rows, "rows": len(rows),
         "failed_rows": sum(1 for r in rows if r["status"] != "ok"), "jobs": len(jobs),
         "distinct_source_excerpts": len(excerpts), "verdicts": verdicts, "decision": decision}
(out_dir / "score.json").write_text(json.dumps(score, indent=2))

md = ["# KIT-151 cap evidence: scored by k151-score.py", "",
      f"Rows {len(rows)} of {expected_rows} expected; failed rows {score['failed_rows']}; jobs {len(jobs)}; "
      f"distinct source excerpts {len(excerpts)}.", ""]
for sc in scenarios:
    md += [f"## Scenario {sc}" + (" (decides)" if sc == "main" else " (sensitivity)"), "",
           "| arm | F1 >=0.9 | F2 >=0.9 | F3 >=0.9 | median abs r F1/F2/F3 | kept D1/D2/D3/D4 | meets rule |",
           "|---|---|---|---|---|---|---|"]
    for arm in arms:
        s = {p: next(x for x in summary if x["scenario"] == sc and x["arm"] == arm and x["planted"] == p)
             for p in ["F1", "F2", "F3"]}
        kept = "/".join(str(s["F1"]["median_kept_" + v]) for v in VIEWS)
        md.append(f"| {arm} | {s['F1']['seeds_abs_r_ge_0.9']}/{len(seeds)} | {s['F2']['seeds_abs_r_ge_0.9']}/{len(seeds)} | "
                  f"{s['F3']['seeds_abs_r_ge_0.9']}/{len(seeds)} | {s['F1']['median_abs_r']}/{s['F2']['median_abs_r']}/"
                  f"{s['F3']['median_abs_r']} | {kept} | {verdicts[sc][arm]['meets_rule']} |")
    md += ["", f"Decision by the rule: {json.dumps(decision.get(sc))}", ""]
(out_dir / "score.md").write_text("\n".join(md) + "\n")
print("\n".join(md))
