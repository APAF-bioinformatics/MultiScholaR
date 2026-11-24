# KIT-150 (branch run-mofa-measure-k150): the rows of the stages, and the comparisons and claims the declaration
# fixes (docs/qa/workflow-execution/kit/kit-150/declaration.md). Standard library only.
#
#   python3 -I k150-summarise.py job <result folder>
#       rows.tsv in the folder: one row per stage (the stage's record, the driver's, validity), and the lifecycle check
#       (lifecycle.txt) when the job ran the lifecycle plan.
#   python3 -I k150-summarise.py final <folder of job artifacts> <output folder>
#       rows.tsv (every job), summary.md (tables, comparisons per mechanism, claims, app-path rows, Python against the
#       dense floor, model sanity) and comparisons.tsv, claims.tsv.
import csv
import glob
import json
import os
import statistics
import sys

BUILDS = ["base", "A", "B1", "B3", "C"]
STAGES = ["write", "import", "align", "train", "results"]
SIZE_LABEL = {"S": "S", "M": "M", "L": "L", "Lp": "L'", "P1": "P1", "P2": "P2", "P3": "P3"}
FIELDS = ["job", "id", "size", "build", "rep", "stage", "case", "outcome", "valid", "heap_peak_mib", "vmhwm_mib",
          "baseline_heap_mib", "baseline_vmrss_mib", "hwm_reset_mib", "wall_s", "load_s", "python_peak_mib",
          "cgroup_peak_mib", "largest_selected_mib", "dense_floor_mib", "selected", "iterations", "factors_sha256",
          "tables_mib", "row_group_mib", "inputs_mib", "marks", "exit_status", "oom_kill", "time_limit_s",
          "memory_limit_mib", "error"]


def load(path):
    try:
        with open(path) as handle:
            return json.load(handle)
    except (OSError, ValueError):
        return None


def num(value):
    if value is None or value == "" or value == "NA":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def first(value):
    return value[0] if isinstance(value, list) and value else value


def rows_of(folder, job):
    rows = []
    for driver_path in sorted(glob.glob(os.path.join(folder, "stages", "*.driver.json"))):
        driver = load(driver_path) or {}
        base = driver_path[: -len(".driver.json")]
        record = load(base + ".json")
        baseline = load(base + ".baseline.json")
        ident = driver.get("id", os.path.basename(base))
        parts = ident.split("_")
        row = {k: "" for k in FIELDS}
        row.update(job=job, id=ident, size=parts[0], build=parts[1], rep=parts[2].lstrip("r"),
                   stage="_".join(parts[3:]), outcome=driver.get("outcome", ""))
        row["exit_status"] = driver.get("exit_status", "")
        row["oom_kill"] = driver.get("oom_kill", "")
        row["cgroup_peak_mib"] = driver.get("cgroup_peak_mib", "")
        row["time_limit_s"] = driver.get("time_limit_s", "")
        row["memory_limit_mib"] = driver.get("memory_limit_mib", "")
        source = record or baseline or {}
        row["case"] = first(source.get("case")) or ""
        b = source.get("baseline") or {}
        row["baseline_heap_mib"] = first(b.get("heap_mib"))
        row["baseline_vmrss_mib"] = first(b.get("vmrss_mib"))
        row["hwm_reset_mib"] = first(b.get("vmhwm_after_reset_mib"))
        row["load_s"] = first(b.get("load_s"))
        if record:
            f = record.get("figures") or {}
            for key in ["heap_peak_mib", "vmhwm_mib", "wall_s", "python_peak_mib"]:
                row[key] = first(f.get(key))
            selected = record.get("selected")
            if isinstance(selected, list) and selected:
                sizes = [r["rows"] * r["samples"] * 8 / 2 ** 20 for r in selected]
                row["largest_selected_mib"] = round(max(sizes), 1)
                row["dense_floor_mib"] = round(sum(sizes), 1)
                row["selected"] = ";".join(f'{r["view"]}={r["rows"]}x{r["samples"]}' for r in selected)
            model = record.get("model") or {}
            row["iterations"] = first(model.get("iterations"))
            row["factors_sha256"] = (first(model.get("factors_sha256")) or "")[:16]
            for key in ["tables_mib", "row_group_mib", "inputs_mib"]:
                row[key] = first(record.get(key))
            marks = record.get("marks") or []
            row["marks"] = " ".join(f'{first(m.get("name"))}@{first(m.get("at_s"))}s/{first(m.get("vmhwm_mib"))}MiB'
                                    for m in marks if isinstance(m, dict))
            if record.get("status") != "completed":
                row["error"] = (first(record.get("error_message")) or "")[:200]
        row["valid"] = validity(row, record, baseline, driver)
        rows.append(row)
    return rows


def validity(row, record, baseline, driver):
    """yes / no (why): a completed stage needs every figure; a stopped one its stop reason (declaration, Terms)."""
    outcome = row["outcome"]
    if outcome.startswith("not run"):
        return "yes"
    if outcome == "completed":
        needed = ["heap_peak_mib", "vmhwm_mib", "wall_s", "baseline_heap_mib", "baseline_vmrss_mib"]
        missing = [k for k in needed if num(row[k]) is None]
        if row["stage"] in ("train", "floor", "align") and not row["selected"]:
            missing.append("selected")
        if row["stage"] in ("train", "floor") and (num(row["python_peak_mib"]) is None or num(row["python_peak_mib"]) < 0):
            missing.append("python_peak_mib (absent or negative)")
        if missing:
            return "no (missing " + ",".join(missing) + ")"
        if num(row["hwm_reset_mib"]) is None or num(row["hwm_reset_mib"]) > num(row["baseline_vmrss_mib"]) + 16:
            return "no (VmHWM not reset at the baseline)"
        return "yes"
    if outcome == "error":
        return "yes" if row["error"] else "no (error without message)"
    if outcome == "killed: memory limit":
        return "yes" if num(row["oom_kill"]) else "no (no oom_kill count)"
    if outcome.startswith("stopped") or outcome.startswith("failed before"):
        return "yes"
    return "no (no outcome)"


def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: ("" if row.get(k) is None else " ".join(row.get(k).split()) if isinstance(row.get(k), str)
                                 else row.get(k)) for k in fields})


def lifecycle_check(rows):
    lines = []
    ok = True
    by = {r["id"]: r for r in rows}
    for b in BUILDS:
        stages = STAGES + (["floor", "app_import", "app_preview"] if b in ("base", "C") else [])
        for s in stages:
            r = by.get(f"S_{b}_r1_{s}")
            good = r is not None and r["outcome"] == "completed" and r["valid"] == "yes"
            ok &= good
            lines.append(f"{'ok ' if good else 'BAD'} S {b} {s}: " + (f'{r["outcome"]}, valid {r["valid"]}' if r else "absent"))
    expected = {"fail_memory": "killed: memory limit", "fail_time": "stopped: time limit", "fail_error": "error"}
    for s, want in expected.items():
        r = by.get(f"S_C_r1_{s}")
        good = (r is not None and r["outcome"] == want and r["valid"] == "yes"
                and num(r["baseline_vmrss_mib"]) is not None)
        ok &= good
        lines.append(f"{'ok ' if good else 'BAD'} induced {s}: expected '{want}', "
                     + (f'got \'{r["outcome"]}\' (baseline {"seen" if num(r["baseline_vmrss_mib"]) is not None else "absent"})' if r else "absent"))
    lines.insert(0, f"lifecycle: {'PASS' if ok else 'FAIL'} (every stage observed through completion; each induced failure observed after its baseline with its stop reason)")
    return "\n".join(lines) + "\n"


# ------------------------------------------------------------------------------------------------------------------
# Final: comparisons and claims
# ------------------------------------------------------------------------------------------------------------------

def values(rows, size, build, stage, metric):
    return [num(r[metric]) for r in rows if r["size"] == size and r["build"] == build and r["stage"] == stage
            and r["outcome"] == "completed" and r["valid"] == "yes" and num(r[metric]) is not None]


def status(rows, size, build, stage):
    """completed / did not complete (reason) / not measured."""
    found = [r for r in rows if r["size"] == size and r["build"] == build and r["stage"] == stage]
    if not found:
        return "not measured"
    if any(r["outcome"] == "completed" and r["valid"] == "yes" for r in found):
        return "completed"
    if any(r["outcome"] == "completed" for r in found):
        return "instrumentation invalid"
    return "did not complete (" + "; ".join(sorted({r["outcome"] for r in found})) + ")"


def median(xs):
    return statistics.median(xs) if xs else None


def spread(xs):
    return (max(xs) - min(xs)) if len(xs) >= 3 else None


def noise(rows, size, stage, metric, builds=("base", "C")):
    """The larger range of 3 repeats of the given builds at this size; at a size without repeats, the nearest smaller
    size's, labelled (declaration, Comparisons)."""
    order = ["S", "M", "Lp", "L"] if size in ("S", "M", "Lp", "L") else [size]
    start = order.index(size)
    for s in order[start::-1]:
        ranges = [spread(values(rows, s, b, stage, metric)) for b in builds]
        ranges = [x for x in ranges if x is not None]
        if ranges:
            return max(ranges), ("" if s == size else f"noise of {SIZE_LABEL[s]}")
    return None, "no noise floor"


def tolerance(before, noise_value, metric):
    floor = 64.0 if metric != "wall_s" else 1.0
    return max(floor, 0.05 * before, 2 * (noise_value or 0.0))


def verdict(before, after, tol):
    diff = after - before
    if diff < -tol:
        return "reduced"
    if diff > tol:
        return "increased"
    return "no measurable change"


MECHANISMS = [
    ("A: Parquet inputs, streamed import", "base", "A", ["write", "import"]),
    ("B1: save_data FALSE, one model load, fewer copies", "A", "B1", ["train"]),
    ("B3: product v2 (training write, Results read)", "B1", "B3", ["train", "results"]),
    ("C: partitioned alignment", "B3", "C", ["align", "train"]),
    ("Final (C) against base", "base", "C", STAGES),
]
METRICS = [("heap_peak_mib", "R heap peak (MiB)"), ("vmhwm_mib", "process peak (MiB)"), ("wall_s", "wall (s)")]


def fmt(x, digits=1):
    return "" if x is None else f"{x:.{digits}f}"


def comparisons(rows, sizes):
    out = []
    for size in sizes:
        for title, before_b, after_b, stages in MECHANISMS:
            for stage in stages:
                sb, sa = status(rows, size, before_b, stage), status(rows, size, after_b, stage)
                for metric, label in METRICS:
                    row = dict(size=SIZE_LABEL[size], mechanism=title, stage=stage, metric=label, before_build=before_b,
                               after_build=after_b)
                    if sb != "completed":
                        row["verdict"] = (("base" if before_b == "base" else "before") + " did not complete: " + sb
                                          if sb.startswith("did not") else f"before ({before_b}): {sb}")
                        out.append(row)
                        continue
                    if sa != "completed":
                        row["verdict"] = ("after did not complete: " + sa if sa.startswith("did not")
                                          else f"after ({after_b}): {sa}")
                        out.append(row)
                        continue
                    vb, va = values(rows, size, before_b, stage, metric), values(rows, size, after_b, stage, metric)
                    b, a = median(vb), median(va)
                    n, note = noise(rows, size, stage, metric)
                    tol = tolerance(b, n, metric)
                    row.update(before=fmt(b), after=fmt(a), change=fmt(a - b), noise=fmt(n), noise_note=note,
                               tolerance=fmt(tol), n_before=len(vb), n_after=len(va), verdict=verdict(b, a, tol))
                    out.append(row)
    return out


def pipeline_peaks(rows, size, build):
    """Per repeat: the largest process peak of the capped pipeline's stages (all five must have completed)."""
    peaks = []
    for rep in sorted({r["rep"] for r in rows if r["size"] == size and r["build"] == build}):
        got = [r for r in rows if r["size"] == size and r["build"] == build and r["rep"] == rep and r["stage"] in STAGES]
        if len(got) == len(STAGES) and all(r["outcome"] == "completed" and r["valid"] == "yes" for r in got):
            peaks.append(max(num(r["vmhwm_mib"]) for r in got))
    return peaks


def claims(rows, sizes):
    out = []
    # Claim 1: training R heap peak of C, capped, at M and L or L'.
    for size in [s for s in sizes if s in ("M", "L", "Lp")]:
        got = [r for r in rows if r["size"] == size and r["build"] == "C" and r["stage"] == "train"]
        st = status(rows, size, "C", "train")
        if st != "completed":
            out.append(dict(claim="1", size=SIZE_LABEL[size], verdict="after did not complete" if st.startswith("did not") else st,
                            detail=""))
            continue
        checks = []
        for r in got:
            if r["outcome"] == "completed" and r["valid"] == "yes":
                bound = 2 * num(r["largest_selected_mib"]) + 256
                checks.append((num(r["heap_peak_mib"]), bound))
        verdict_1 = "pass" if all(h <= bnd for h, bnd in checks) else "fail"
        out.append(dict(claim="1", size=SIZE_LABEL[size], verdict=verdict_1,
                        detail="; ".join(f"heap peak {h:.1f} MiB vs bound {bnd:.1f} MiB" for h, bnd in checks)))
    # Claim 2: total memory, C against base, capped, at every size base completes.
    for size in [s for s in sizes if s in ("S", "M", "L", "Lp")]:
        pb, pa = pipeline_peaks(rows, size, "base"), pipeline_peaks(rows, size, "C")
        if not pb:
            out.append(dict(claim="2", size=SIZE_LABEL[size], verdict="base did not complete",
                            detail="; ".join(f"{s}: {status(rows, size, 'base', s)}" for s in STAGES)))
            continue
        if not pa:
            out.append(dict(claim="2", size=SIZE_LABEL[size], verdict="after did not complete", detail=""))
            continue
        b, a = median(pb), median(pa)
        ranges = [x for x in (spread(pb), spread(pa)) if x is not None]
        n = max(ranges) if ranges else None
        note = ""
        if n is None:
            n, note = noise(rows, size, "train", "vmhwm_mib")
        tol = tolerance(b, n, "vmhwm_mib")
        out.append(dict(claim="2", size=SIZE_LABEL[size], verdict="pass" if a <= b + tol else "fail",
                        detail=f"pipeline process peak base {b:.1f} MiB (n={len(pb)}), C {a:.1f} MiB (n={len(pa)}), "
                               f"tolerance {tol:.1f} MiB (noise {fmt(n)} {note})"))
    # Claim 3: import of C, the import pair.
    pair = {s: values(rows, s, "C", "import", "heap_peak_mib") for s in ("P1", "P2", "P3")}
    if all(pair.values()):
        details = []
        bound_ok = True
        for s in ("P1", "P2", "P3"):
            got = [r for r in rows if r["size"] == s and r["build"] == "C" and r["stage"] == "import"
                   and r["outcome"] == "completed" and r["valid"] == "yes"]
            for r in got:
                bound = 4 * num(r["row_group_mib"]) + num(r["tables_mib"]) + 64
                bound_ok &= num(r["heap_peak_mib"]) <= bound
            details.append(f"{s}: median {median(pair[s]):.1f} MiB, bound {4 * num(got[0]['row_group_mib']) + num(got[0]['tables_mib']) + 64:.1f} MiB")
        n12 = max([x for x in (spread(pair["P1"]), spread(pair["P2"])) if x is not None] or [0.0])
        n13 = max([x for x in (spread(pair["P1"]), spread(pair["P3"])) if x is not None] or [0.0])
        p1, p2, p3 = median(pair["P1"]), median(pair["P2"]), median(pair["P3"])
        tol_f = tolerance(p1, n12, "heap_peak_mib")
        tol_r = tolerance(p3, n13, "heap_peak_mib")
        features_ok = p2 - p1 <= tol_f
        runs_ok = p1 - p3 <= p3 + tol_r
        details.append(f"features doubled: +{p2 - p1:.1f} MiB (tolerance {tol_f:.1f})")
        details.append(f"runs doubled: +{p1 - p3:.1f} MiB (allowed {p3 + tol_r:.1f})")
        out.append(dict(claim="3", size="P1/P2/P3",
                        verdict="pass" if (bound_ok and features_ok and runs_ok) else "fail",
                        detail="; ".join(details) + f"; bound {'met' if bound_ok else 'NOT met'}, "
                               f"features {'met' if features_ok else 'NOT met'}, runs {'met' if runs_ok else 'NOT met'}"))
    else:
        out.append(dict(claim="3", size="P1/P2/P3", verdict="after did not complete",
                        detail="; ".join(f"{s}: {status(rows, s, 'C', 'import')}" for s in ("P1", "P2", "P3"))))
    return out


def stage_table(rows, sizes):
    lines = ["| Size | Build | Stage | Outcome | R heap peak MiB (range) | Process peak MiB (range) | Wall s (range) | n |",
             "|---|---|---|---|---|---|---|---|"]
    for size in sizes:
        for build in BUILDS:
            for stage in STAGES + ["floor", "app_import", "app_preview"]:
                found = [r for r in rows if r["size"] == size and r["build"] == build and r["stage"] == stage]
                if not found:
                    continue
                st = status(rows, size, build, stage)
                cells = []
                for metric in ("heap_peak_mib", "vmhwm_mib", "wall_s"):
                    xs = values(rows, size, build, stage, metric)
                    if not xs or (metric == "wall_s" and stage == "floor"):
                        cells.append("n.r." if xs else "")
                        continue
                    rng = spread(xs)
                    cells.append(f"{median(xs):.1f}" + (f" ({rng:.1f})" if rng is not None else ""))
                n = len(values(rows, size, build, stage, "vmhwm_mib"))
                lines.append(f"| {SIZE_LABEL[size]} | {build} | {stage} | {st} | {cells[0]} | {cells[1]} | {cells[2]} | {n} |")
    return "\n".join(lines)


def python_table(rows, sizes):
    lines = ["| Size | Build | Case | Python peak MiB | Dense floor MiB | Largest selected MiB | Process peak MiB | Iterations |",
             "|---|---|---|---|---|---|---|---|"]
    for size in sizes:
        for build in BUILDS:
            for stage, case in (("train", "capped"), ("floor", "uncapped, 10 rounds")):
                got = [r for r in rows if r["size"] == size and r["build"] == build and r["stage"] == stage
                       and r["outcome"] == "completed" and r["valid"] == "yes"]
                if not got:
                    continue
                r = got[0]
                lines.append(f"| {SIZE_LABEL[size]} | {build} | {case} | {r['python_peak_mib']} | {r['dense_floor_mib']} | "
                             f"{r['largest_selected_mib']} | {r['vmhwm_mib']} | {r['iterations']} |")
    return "\n".join(lines)


def sanity(rows, sizes):
    lines = []
    for size in sizes:
        for stage in ("train", "floor"):
            got = {}
            for r in rows:
                if r["size"] == size and r["stage"] == stage and r["outcome"] == "completed":
                    got.setdefault(r["build"], set()).add((r["factors_sha256"], r["iterations"]))
            if got:
                distinct = {v for s in got.values() for v in s}
                lines.append(f"- {SIZE_LABEL[size]} {stage}: " + ("the same factors (rounded to 1e-8) and iterations in every build"
                             if len(distinct) == 1 else "DIFFERENT: " + "; ".join(f"{b} {sorted(s)}" for b, s in sorted(got.items()))))
    return "\n".join(lines)


def main():
    mode = sys.argv[1]
    if mode == "job":
        folder = sys.argv[2]
        rows = rows_of(folder, os.environ.get("K150_JOB", "job"))
        write_tsv(os.path.join(folder, "rows.tsv"), rows, FIELDS)
        for r in rows:
            print(f'{r["id"]}\t{r["outcome"]}\tvalid {r["valid"]}\theap {r["heap_peak_mib"]}\thwm {r["vmhwm_mib"]}\twall {r["wall_s"]}')
        if os.environ.get("K150_PLAN") == "g0i":
            by = {r["id"]: r for r in rows}
            lines = []
            ok = True
            for b in ("RC", "C"):
                for s in STAGES:
                    r = by.get(f"S_{b}_r1_{s}")
                    good = r is not None and r["outcome"] == "completed" and r["valid"] == "yes"
                    ok &= good
                    lines.append(f"{'ok ' if good else 'BAD'} S {b} {s}: " + (f'{r["outcome"]}, valid {r["valid"]}, hwm {r["vmhwm_mib"]}, wall {r["wall_s"]}' if r else "absent"))
            r = by.get("S_RC_r1_fail_memory")
            good = r is not None and r["outcome"] == "killed: memory limit" and num(r["baseline_vmrss_mib"]) is not None
            ok &= good
            lines.append(f"{'ok ' if good else 'BAD'} induced fail_memory on RC: " + (f'{r["outcome"]}, baseline {r["baseline_vmrss_mib"]}' if r else "absent"))
            text = f"g0i: {'held' if ok else 'did not hold'}\n" + "\n".join(lines) + "\n"
            with open(os.path.join(folder, "g0i.txt"), "w") as handle:
                handle.write(text)
            print(text)
        if os.environ.get("K150_PLAN") == "lifecycle":
            text = lifecycle_check(rows)
            with open(os.path.join(folder, "lifecycle.txt"), "w") as handle:
                handle.write(text)
            print(text)
        return
    jobs, out = sys.argv[2], sys.argv[3]
    os.makedirs(out, exist_ok=True)
    rows = []
    for folder in sorted(glob.glob(os.path.join(jobs, "*"))):
        if os.path.isdir(os.path.join(folder, "stages")):
            rows.extend(rows_of(folder, os.path.basename(folder)))
    write_tsv(os.path.join(out, "rows.tsv"), rows, FIELDS)
    order = ["S", "M", "Lp", "L", "P1", "P2", "P3"]
    sizes = [s for s in order if any(r["size"] == s for r in rows)]
    main_sizes = [s for s in sizes if s in ("S", "M", "Lp", "L")]
    comp = comparisons(rows, main_sizes)
    comp_fields = ["size", "mechanism", "stage", "metric", "before_build", "after_build", "before", "after", "change",
                   "noise", "noise_note", "tolerance", "n_before", "n_after", "verdict"]
    write_tsv(os.path.join(out, "comparisons.tsv"), comp, comp_fields)
    cl = claims(rows, main_sizes)
    write_tsv(os.path.join(out, "claims.tsv"), cl, ["claim", "size", "verdict", "detail"])
    invalid = [r for r in rows if r["valid"] != "yes"]
    md = ["# KIT-150 measurement summary", "",
          f"Rows: {len(rows)}; instrumentation invalid: {len(invalid)}" + (" (" + ", ".join(r["id"] for r in invalid) + ")" if invalid else ""),
          "", "## Stages (median of repeats; range in brackets when 3 repeats)", "", stage_table(rows, sizes), "",
          "## Per mechanism (median; tolerance max(64 MiB or 1 s, 5 %, 2 x noise floor))", "",
          "| Size | Mechanism | Stage | Metric | Before | After | Change | Noise | Tolerance | Verdict |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for c in comp:
        md.append(f"| {c['size']} | {c['mechanism']} | {c['stage']} | {c['metric']} | {c.get('before', '')} | "
                  f"{c.get('after', '')} | {c.get('change', '')} | {c.get('noise', '')} {c.get('noise_note', '')} | "
                  f"{c.get('tolerance', '')} | {c['verdict']} |")
    md += ["", "## Claims", "", "| Claim | Size | Verdict | Detail |", "|---|---|---|---|"]
    md += [f"| {c['claim']} | {c['size']} | {c['verdict']} | {c['detail']} |" for c in cl]
    md += ["", "## Claim 4: Python against the dense floor (reported, not claimed)", "", python_table(rows, sizes), "",
           "## Model sanity", "", sanity(rows, main_sizes), ""]
    with open(os.path.join(out, "summary.md"), "w") as handle:
        handle.write("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
