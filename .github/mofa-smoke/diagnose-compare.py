"""KIT-062 diagnosis: compare the diagnostic jobs' result.json files.

Usage: python3 -I diagnose-compare.py <dir holding diag-*/result.json> <stored reference dir> <output dir>
Per combination (threads x precision): every pair of jobs, bit-identical or not, same CPU model or not, the largest
absolute and relative differences of the factors, weights, variance explained and last ELBO, and the iterations.
Across combinations: the same. Path independence: variance explained per view and factor and the planted-factor
correlations rounded to 3 decimals, and the best absolute correlation between matched factors of two jobs. The
stored result (tests/testdata/integration/mofa_reference, 17 significant digits) is compared with each job.
Writes diagnose.md and diagnose.json.
"""
import csv
import itertools
import json
import math
import pathlib
import sys

root, reference, out = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2]), pathlib.Path(sys.argv[3])
out.mkdir(parents=True, exist_ok=True)
jobs = {}
for path in sorted(root.glob("*/result.json")):
    r = json.loads(path.read_text())
    if r.get("trained"):
        jobs[path.parent.name] = r


def combo(r):
    threads = "single" if r["environment"].get("OPENBLAS_NUM_THREADS") == "1" else "default"
    return f"{threads}-{r['precision']}"


def block(r, name):
    m = r["model"]
    if name == "factors":
        return [float.fromhex(v) for v in m["factors"]["values_hex"]]
    if name == "weights":
        return [float.fromhex(v) for view in m["weights"] for v in m["weights"][view]["values_hex"]]
    if name == "r2":
        return [float.fromhex(v) for v in m["r2"]["values_hex"]]
    if name == "elbo":
        return [float.fromhex(m["elbo_last_hex"])]


def matrix(r):
    m = r["model"]["factors"]
    n, k = m["dim"]
    v = [float.fromhex(x) for x in m["values_hex"]]
    return [[v[j * n + i] for j in range(k)] for i in range(n)]


def corr(a, b):
    ma, mb = sum(a) / len(a), sum(b) / len(b)
    sa = math.sqrt(sum((x - ma) ** 2 for x in a))
    sb = math.sqrt(sum((x - mb) ** 2 for x in b))
    return sum((x - ma) * (y - mb) for x, y in zip(a, b)) / (sa * sb)


def matched(a, b):
    """For each factor of a, the largest |r| with any factor of b; the smallest of those."""
    fa, fb = list(zip(*matrix(a))), list(zip(*matrix(b)))
    return min(max(abs(corr(x, y)) for y in fb) for x in fa)


def diff(a, b):
    if len(a) != len(b):
        return None
    d = max(abs(x - y) for x, y in zip(a, b))
    scale = max(max(abs(x) for x in a), 1e-300)
    return d, d / scale


def pair(a, b):
    row = {"same_cpu": jobs[a]["cpu_model"] == jobs[b]["cpu_model"],
           "iterations": f"{jobs[a]['model']['iterations']}/{jobs[b]['model']['iterations']}"}
    exact = True
    for name in ("factors", "weights", "r2", "elbo"):
        x, y = block(jobs[a], name), block(jobs[b], name)
        same = x == y
        exact = exact and same
        dd = diff(x, y)
        row[name] = "identical" if same else ("lengths differ" if dd is None else f"{dd[0]:.3g} ({dd[1]:.3g} rel)")
    row["bit_identical"] = exact
    row["r2_equal_3dp"] = [round(v, 3) for v in block(jobs[a], "r2")] == [round(v, 3) for v in block(jobs[b], "r2")]
    pa, pb = jobs[a]["model"]["planted_abs_r"], jobs[b]["model"]["planted_abs_r"]
    row["planted_equal_3dp"] = all(round(pa[z], 3) == round(pb[z], 3) for z in ("Z1", "Z2", "Z3"))
    row["min_matched_factor_abs_r"] = round(matched(jobs[a], jobs[b]), 6)
    return row


summary = {"jobs": {}, "within": {}, "across": {}, "reference": {}}
lines = ["### KIT-062 MOFA2 diagnosis: jobs", "",
         "| job | combination | CPU | numpy BLAS threads | float32 option (default) | iterations | train s | ELBO last | Z1 | Z2 | Z3 | factors md5 |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|"]
for name, r in jobs.items():
    t = ";".join(f"{d.get('internal_api')}:{d.get('num_threads')}" for d in r.get("python_threads_after") or [])
    p = r["model"]["planted_abs_r"]
    summary["jobs"][name] = {"combination": combo(r), "cpu": r["cpu_model"], "threads": t,
                             "iterations": r["model"]["iterations"], "seconds_train": r["seconds_train"],
                             "elbo_last": r["model"]["elbo_last"], "planted": p, "float32_option": r.get("float32_option"),
                             "factors_md5": r["model"]["factors"]["md5"]}
    lines.append(f"| {name} | {combo(r)} | {r['cpu_model']} | {t} | {r.get('float32_option')} ({r.get('float32_default')}) | "
                 f"{r['model']['iterations']} | {r['seconds_train']} | {r['model']['elbo_last']:.6f} | "
                 f"{p['Z1']:.4f} | {p['Z2']:.4f} | {p['Z3']:.4f} | {r['model']['factors']['md5'][:12]} |")

header = ["| pair | same CPU | bit-identical | iterations | factors | weights | variance explained | ELBO | r2 equal (3 dp) | planted equal (3 dp) | min matched factor abs r |",
          "|---|---|---|---|---|---|---|---|---|---|---|"]


def emit(title, pairs, key):
    lines.extend(["", title, ""] + header)
    for a, b in pairs:
        row = pair(a, b)
        summary[key][f"{a} vs {b}"] = row
        lines.append(f"| {a} vs {b} | {row['same_cpu']} | {row['bit_identical']} | {row['iterations']} | {row['factors']} | "
                     f"{row['weights']} | {row['r2']} | {row['elbo']} | {row['r2_equal_3dp']} | {row['planted_equal_3dp']} | "
                     f"{row['min_matched_factor_abs_r']} |")


combos = sorted({combo(r) for r in jobs.values()})
for c in combos:
    members = [n for n, r in jobs.items() if combo(r) == c]
    emit(f"### Within {c}", list(itertools.combinations(members, 2)), "within")
first = {c: next(n for n, r in jobs.items() if combo(r) == c) for c in combos}
emit("### Across combinations (first job of each)", list(itertools.combinations(first.values(), 2)), "across")

# The stored result (written by the package lane's recording run).
stored = {}
for name in ("factors", "variance_explained"):
    path = reference / f"{name}.tsv"
    if path.exists():
        with path.open() as f:
            rows = list(csv.reader(f, delimiter="\t"))
        stored[name] = rows
lines.extend(["", "### Against the stored result (package lane run 44459109264)", "",
              "| job | factors identical | factors max abs diff |", "|---|---|---|"])
if "factors" in stored:
    rows = stored["factors"][1:]
    k = len(stored["factors"][0]) - 1
    by_sample = {row[0]: [float(v) for v in row[1:]] for row in rows}
    for name, r in jobs.items():
        mat = matrix(r)
        samples = r["model"]["factors"]["samples"]
        if len(mat[0]) != k:
            lines.append(f"| {name} | different factor count | |")
            continue
        d = max(abs(mat[i][j] - by_sample[s][j]) for i, s in enumerate(samples) for j in range(k))
        summary["reference"][name] = d
        lines.append(f"| {name} | {d == 0} | {d:.3g} |")

(out / "diagnose.md").write_text("\n".join(lines) + "\n")
(out / "diagnose.json").write_text(json.dumps(summary, indent=2))
print("\n".join(lines))
