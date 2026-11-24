# KIT-209 (branch run-kit209-memory): the deliverable's tables from the jobs' journeys.tsv and steps.tsv (ticket 4.1,
# 7.1 to 7.3). Memory discovery (VAL-151): figures with their spread; the listing aid only orders what a reader sees
# first. Standard library only.
#
#   python3 -I k209-cells.py <output folder> <job folder> [<job folder> ...]
# Writes cells.tsv (unit, step, metric: each build's runs, median, range, n; the difference in MiB and %; listed or not),
# units.tsv (the same for the journey maximum, the retained-heap maximum and the journey wall; the order contrast),
# nonc.tsv (every non-completion with its class) and findings.tsv (the listed rows, largest first).
import csv
import os
import statistics
import sys

METRICS = {"vmhwm_mib": "step process peak (MiB)", "entry_heap_mib": "retained R heap at step entry (MiB)",
           "seconds": "step wall (s)"}
JOURNEY_METRICS = {"e1_vmhwm_mib": "journey peak (MiB)", "e2_heap_mib": "retained-heap maximum (MiB)",
                   "journey_wall_s": "journey wall (s)"}
MEMORY = {"vmhwm_mib", "entry_heap_mib", "e1_vmhwm_mib", "e2_heap_mib"}


def read(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def num(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def stats(values):
    values = [v for v in values if v is not None]
    if not values:
        return None
    return {"n": len(values), "median": round(statistics.median(values), 1), "min": round(min(values), 1),
            "max": round(max(values), 1), "range": round(max(values) - min(values), 1),
            "runs": ";".join(str(round(v, 1)) for v in values)}


def listed(metric, base, cand):
    # The aid applies to the memory figures (ticket 7.3); time is context and is never listed.
    if base is None or cand is None or metric not in MEMORY:
        return False, "" if metric in MEMORY else "time: context, not listed"
    diff = cand["median"] - base["median"]
    floor = max(256.0, 0.10 * base["median"])
    spread = 2 * max(base["range"], cand["range"])
    return abs(diff) > floor and abs(diff) > spread, f"floor {floor:.1f}, twice the larger range {spread:.1f}"


def row(unit, step, key, label, base, cand):
    out = {"unit": unit, "step": step, "metric": label}
    metric = key
    for label, s in (("base", base), ("candidate", cand)):
        for key in ("n", "median", "min", "max", "range", "runs"):
            out[f"{label}_{key}"] = "" if s is None else s[key]
    if base is not None and cand is not None:
        diff = round(cand["median"] - base["median"], 1)
        out["difference"] = diff
        out["difference_pct"] = round(100 * diff / base["median"], 1) if base["median"] else ""
        out["listed"], out["listing_basis"] = listed(metric, base, cand)
    else:
        out.update(difference="", difference_pct="", listed="", listing_basis="one build only")
    return out


def main():
    output, folders = sys.argv[1], sys.argv[2:]
    os.makedirs(output, exist_ok=True)
    journeys, steps = [], []
    for folder in folders:
        journeys += read(os.path.join(folder, "journeys.tsv"))
        steps += read(os.path.join(folder, "steps.tsv"))
    counted = [j for j in journeys if j["id"].split("_")[0][:1] in ("U", "S") and j["build"] in ("base", "C")]
    complete = {j["id"] for j in counted if j["class"] == "completed"}
    cells, units, nonc = [], [], []
    by_unit = {}
    for j in counted:
        by_unit.setdefault(j["unit"], []).append(j)
        if j["class"] != "completed":
            nonc.append({k: j[k] for k in ("id", "unit", "build", "rep", "class", "outcome", "error_step", "error_message")})
    for unit, js in sorted(by_unit.items()):
        for metric, label in JOURNEY_METRICS.items():
            b = stats([num(j[metric]) for j in js if j["build"] == "base" and j["id"] in complete])
            c = stats([num(j[metric]) for j in js if j["build"] == "C" and j["id"] in complete])
            r = row(unit, "journey", metric, label, b, c)
            # Order contrast: the candidate-minus-base difference in pairs run B C against pairs run C B.
            pairs = {}
            for j in js:
                if j["id"] in complete:
                    pairs.setdefault((j["rep"], j["order"]), {})[j["build"]] = num(j[metric])
            bc = [p["C"] - p["base"] for (rep, order), p in pairs.items() if order == "BC" and "C" in p and "base" in p]
            cb = [p["C"] - p["base"] for (rep, order), p in pairs.items() if order == "CB" and "C" in p and "base" in p]
            r["order_contrast"] = (f"B C mean {statistics.mean(bc):.1f} (n {len(bc)}); C B mean {statistics.mean(cb):.1f} (n {len(cb)})"
                                   if bc and cb else "")
            units.append(r)
    step_runs = {}
    for s in steps:
        if s["id"].split("_")[0][:1] not in ("U", "S") or s["build"] not in ("base", "C"):
            continue
        if s["status"] != "completed":
            continue
        key = (s["unit"], int(s["index"]), s["step"])
        for metric in METRICS:
            step_runs.setdefault((key, metric), {"base": [], "C": []})[s["build"]].append(num(s[metric]))
    states = {}
    for s in steps:
        states.setdefault((s["unit"], s["step"], s["build"]), set()).add(s["state_as_expected"])
    for ((unit, index, step), metric), runs in sorted(step_runs.items(), key=lambda x: (x[0][0][0], x[0][0][1], x[0][1])):
        r = row(unit, step, metric, METRICS[metric], stats(runs["base"]), stats(runs["C"]))
        r["index"] = index
        r["states_committed_base"] = ",".join(sorted(states.get((unit, step, "base"), set())))
        r["states_committed_candidate"] = ",".join(sorted(states.get((unit, step, "C"), set())))
        cells.append(r)
    fields = ["unit", "index", "step", "metric"] + [f"{b}_{k}" for b in ("base", "candidate") for k in ("n", "median", "min", "max", "range", "runs")] + \
        ["difference", "difference_pct", "listed", "listing_basis", "states_committed_base", "states_committed_candidate", "order_contrast"]
    for name, rows in (("cells.tsv", cells), ("units.tsv", units)):
        with open(os.path.join(output, name), "w", newline="") as handle:
            w = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
    with open(os.path.join(output, "nonc.tsv"), "w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=["id", "unit", "build", "rep", "class", "outcome", "error_step", "error_message"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(nonc)
    found = [r for r in cells + units if r.get("listed") is True]
    found.sort(key=lambda r: -abs(r["difference"]))
    with open(os.path.join(output, "findings.tsv"), "w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(found)
    print(f"cells {len(cells)}, unit rows {len(units)}, non-completions {len(nonc)}, listed {len(found)}")


if __name__ == "__main__":
    main()
