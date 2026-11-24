# KIT-209 (branch run-kit209-memory): the journey and step rows, the outcome classes (declaration, "Why a run did not
# complete") and, for G0a, its checks. Memory discovery (VAL-151): figures are reported, nothing is judged here except
# G0a's lifecycle checks. Standard library only.
#
#   python3 -I k209-summarise.py job <result folder>
#       journeys.tsv, steps.tsv and (plan g0a) g0a.tsv in the folder.
#   python3 -I k209-summarise.py final <folder of job artifacts> <output folder>
#       journeys.tsv, steps.tsv, g0a.tsv of every job and summary.md.
import csv
import glob
import json
import os
import sys

JOURNEY_FIELDS = ["job", "id", "unit", "build", "rep", "order", "inject", "route", "corpus", "outcome", "class",
                  "after_baseline", "steps_done", "steps_total", "e1_vmhwm_mib", "e1_step", "e2_heap_mib", "e2_boundary",
                  "final_heap_mib", "final_vmrss_mib", "baseline_heap_mib", "baseline_vmrss_mib", "load_s",
                  "journey_wall_s", "process_wall_s", "cgroup_peak_mib", "oom_kill", "memory_limit_bytes",
                  "time_limit_s", "cache_dropped", "e3_all_committed", "max_cgroup_procs", "error_step", "error_class",
                  "error_message", "deepest_owner"]
STEP_FIELDS = ["job", "id", "unit", "build", "rep", "index", "step", "status", "seconds", "vmhwm_mib",
               "vmhwm_after_reset_mib", "entry_vmrss_mib", "entry_heap_mib", "heap_max_used_mib", "exit_vmrss_mib",
               "cgroup_procs", "state", "expected_state", "state_as_expected", "rows", "store_bytes"]


def load(path):
    try:
        with open(path) as handle:
            return json.load(handle)
    except (OSError, ValueError):
        return None


def first(value):
    return value[0] if isinstance(value, list) and value else value


def num(value):
    value = first(value)
    try:
        return None if value is None else float(value)
    except (TypeError, ValueError):
        return None


def deepest_owner(frames):
    owners = [first(f.get("owner")) for f in frames or [] if first(f.get("owner")) in ("lane", "package")]
    return owners[-1] if owners else ""


def classify(driver, record, baseline):
    outcome = driver.get("outcome", "")
    if not driver:
        return "runner lost", ""
    if outcome == "completed":
        return "completed", ""
    if outcome.startswith("not run"):
        return "not run", ""
    if (driver.get("oom_kill") or 0) > 0:
        return "memory", ""
    if outcome == "stopped: time limit":
        return "time", ""
    if outcome == "error" and record:
        error = record.get("error") or {}
        kind = first(error.get("class")) or ""
        owner = deepest_owner(error.get("frames"))
        if kind == "k209_state_not_committed":
            # A tab caught the error and only notified: the signalled error's frames classify it, else E3.
            if owner == "package":
                return "package error", owner
            if owner == "lane":
                return "instrument", owner
            return "state not committed (E3)", owner
        if owner == "package":
            return "package error", owner
        if owner == "lane":
            return "instrument", owner
        return "unclear (to the coordinator)", owner
    if outcome.startswith("failed before the baseline"):
        return "instrument", ""
    return "unclear (to the coordinator)", ""


def journey_rows(folder, job):
    journeys, steps = [], []
    for driver_path in sorted(glob.glob(os.path.join(folder, "journeys", "*.driver.json"))):
        driver = load(driver_path) or {}
        base = driver_path[: -len(".driver.json")]
        record = load(base + ".json")
        baseline = load(base + ".baseline.json")
        marks_file = load(base + ".marks.json") or {}
        source = record or baseline or marks_file or {}
        marks = (record or {}).get("marks") or marks_file.get("marks") or []
        cls, owner = classify(driver, record, baseline)
        row = {k: "" for k in JOURNEY_FIELDS}
        row.update(job=job, id=driver.get("id", os.path.basename(base)), outcome=driver.get("outcome", ""), **{"class": cls})
        for key in ("unit", "build", "rep", "order", "inject", "route", "corpus"):
            row[key] = first(source.get(key)) if source.get(key) is not None else ""
        row["after_baseline"] = baseline is not None
        b = (baseline or {}).get("baseline") or (record or {}).get("baseline") or {}
        row["baseline_heap_mib"] = first(b.get("heap_mib"))
        row["baseline_vmrss_mib"] = first(b.get("vmrss_mib"))
        row["load_s"] = first(b.get("load_s"))
        total = (baseline or {}).get("steps") or []
        row["steps_total"] = len(total)
        row["steps_done"] = sum(1 for m in marks if first(m.get("status")) == "completed")
        peaks = [(num(m.get("vmhwm_mib")), first(m.get("step"))) for m in marks if num(m.get("vmhwm_mib")) is not None]
        if peaks:
            row["e1_vmhwm_mib"], row["e1_step"] = max(peaks)
        final = (record or {}).get("final") or {}
        heaps = [(num(m.get("entry_heap_mib")), "before " + str(first(m.get("step")))) for m in marks
                 if num(m.get("entry_heap_mib")) is not None]
        if num(final.get("heap_mib")) is not None:
            heaps.append((num(final.get("heap_mib")), "after the last step"))
        if heaps:
            row["e2_heap_mib"], row["e2_boundary"] = max(heaps)
        row["final_heap_mib"] = first(final.get("heap_mib"))
        row["final_vmrss_mib"] = first(final.get("vmrss_mib"))
        row["journey_wall_s"] = first(final.get("wall_s"))
        for key in ("process_wall_s", "cgroup_peak_mib", "oom_kill", "memory_limit_bytes", "time_limit_s", "cache_dropped"):
            row[key] = driver.get(key, "")
        judged = [first(m.get("state_as_expected")) for m in marks if first(m.get("state_as_expected")) is not None]
        row["e3_all_committed"] = all(judged) if judged else ""
        procs = [num(m.get("cgroup_procs")) for m in marks if num(m.get("cgroup_procs")) is not None]
        row["max_cgroup_procs"] = int(max(procs)) if procs else ""
        error = (record or {}).get("error") or {}
        row["error_step"] = first(error.get("step")) or ""
        row["error_class"] = first(error.get("class")) or ""
        row["error_message"] = (first(error.get("message")) or "").replace("\t", " ").replace("\n", " ")[:300]
        row["deepest_owner"] = owner
        journeys.append(row)
        for m in marks:
            s = {k: "" for k in STEP_FIELDS}
            s.update(job=job, id=row["id"], unit=row["unit"], build=row["build"], rep=row["rep"])
            for key in STEP_FIELDS[5:]:
                if key in m:
                    s[key] = first(m.get(key))
            steps.append(s)
    return journeys, steps


def g0a_checks(folder, journeys, steps):
    by_id = {r["id"]: r for r in journeys}
    checks = []

    def add(item, passed, evidence):
        checks.append({"item": item, "passed": passed, "evidence": evidence})

    life = by_id.get("U1_C_r1", {})
    add("1 lifecycle completed, every mark through the last step",
        life.get("class") == "completed" and life.get("steps_done") == life.get("steps_total") and life.get("steps_total", 0) > 0,
        f"{life.get('class')}; {life.get('steps_done')} of {life.get('steps_total')} steps marked")
    limits = load(os.path.join(folder, "injected-limits.json")) or {}
    mem = by_id.get("U1_C_r1_fail_memory", {})
    add("1 injected memory limit (0.7 x journey VmHWM): killed after the baseline, classed memory",
        mem.get("class") == "memory" and mem.get("after_baseline") is True,
        f"limit {limits.get('injected_memory_max_bytes')} B; outcome {mem.get('outcome')}; after baseline {mem.get('after_baseline')}; "
        f"steps marked {mem.get('steps_done')}")
    tim = by_id.get("U1_C_r1_fail_time", {})
    add("1 injected time limit (0.5 x process wall): stopped, classed time", tim.get("class") == "time",
        f"limit {limits.get('injected_time_limit_s')} s; outcome {tim.get('outcome')}")
    err = by_id.get("U1_C_r1_package_error", {})
    add("1 injected package error (invalid bound setting): raised in a package frame, classed package error",
        err.get("class") == "package error",
        f"outcome {err.get('outcome')}; step {err.get('error_step')}; deepest owner {err.get('deepest_owner')}; "
        f"{err.get('error_message', '')[:160]}")
    cal = {s["step"]: s for s in steps if s["id"] == "CAL_C_r1"}
    c, n = cal.get("calibrate", {}), cal.get("noop", {})
    try:
        c_ok = float(c["vmhwm_mib"]) >= float(c["entry_vmrss_mib"]) + 480
        n_ok = float(n["vmhwm_mib"]) <= float(n["entry_vmrss_mib"]) + 32
    except (KeyError, TypeError, ValueError):
        c_ok = n_ok = False
    add("2 reset control: calibration VmHWM >= entry VmRSS + 480 MiB; no-op VmHWM <= entry VmRSS + 32 MiB", c_ok and n_ok,
        f"calibrate {c.get('vmhwm_mib')} vs entry {c.get('entry_vmrss_mib')}; no-op {n.get('vmhwm_mib')} vs entry {n.get('entry_vmrss_mib')}")
    for unit in ("U1", "U2", "U3", "U4", "U5", "S1", "S2"):
        base, cand = by_id.get(f"{unit}_base_r1", {}), by_id.get(f"{unit}_C_r2", {})
        add(f"3 gate B3 {unit}: base completes; candidate completes",
            base.get("class") == "completed" and cand.get("class") == "completed",
            f"base {base.get('class')} ({base.get('error_step')} {base.get('error_message', '')[:120]}); "
            f"candidate {cand.get('class')} ({cand.get('error_step')} {cand.get('error_message', '')[:120]})")
    return checks


def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    mode = sys.argv[1]
    if mode == "job":
        folder = sys.argv[2]
        job = os.environ.get("K209_JOB", "")
        journeys, steps = journey_rows(folder, job)
        write_tsv(os.path.join(folder, "journeys.tsv"), journeys, JOURNEY_FIELDS)
        write_tsv(os.path.join(folder, "steps.tsv"), steps, STEP_FIELDS)
        if os.environ.get("K209_PLAN") in ("g0a", "g0a-error"):
            checks = g0a_checks(folder, journeys, steps)
            write_tsv(os.path.join(folder, "g0a.tsv"), checks, ["item", "passed", "evidence"])
            for c in checks:
                print(("PASS " if c["passed"] else "FAIL ") + c["item"] + " | " + c["evidence"])
        for r in journeys:
            print(f"{r['id']}: {r['class']} E1 {r['e1_vmhwm_mib']} MiB at {r['e1_step']}; E2 {r['e2_heap_mib']} MiB; "
                  f"{r['journey_wall_s']} s")
    elif mode == "final":
        source, output = sys.argv[2], sys.argv[3]
        os.makedirs(output, exist_ok=True)
        lines = ["# KIT-209 run summary", ""]
        for name, fields in (("journeys.tsv", JOURNEY_FIELDS), ("steps.tsv", STEP_FIELDS), ("g0a.tsv", ["item", "passed", "evidence"])):
            rows = []
            for path in sorted(glob.glob(os.path.join(source, "*", name))):
                with open(path) as handle:
                    rows.extend(csv.DictReader(handle, delimiter="\t"))
            if rows:
                write_tsv(os.path.join(output, name), rows, list(rows[0].keys()))
            if name == "journeys.tsv":
                lines += ["| journey | class | E1 MiB (step) | E2 MiB | wall s |", "|---|---|---|---|---|"]
                lines += [f"| {r['id']} | {r['class']} | {r['e1_vmhwm_mib']} ({r['e1_step']}) | {r['e2_heap_mib']} | {r['journey_wall_s']} |"
                          for r in rows]
            if name == "g0a.tsv" and rows:
                lines += ["", "| G0a item | passed | evidence |", "|---|---|---|"]
                lines += [f"| {r['item']} | {r['passed']} | {r['evidence'][:200]} |" for r in rows]
        with open(os.path.join(output, "summary.md"), "w") as handle:
            handle.write("\n".join(lines) + "\n")
    else:
        sys.exit("mode: job | final")


if __name__ == "__main__":
    main()
