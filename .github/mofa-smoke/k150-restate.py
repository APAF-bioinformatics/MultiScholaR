# KIT-150 restatement of claims 1 and 3a (user decision 2026-11-21, ledger INT-27), made AFTER the results were seen,
# from the counted run's rows (docs/qa/workflow-execution/kit/kit-150/rows.tsv, run 44795728762); no new runs.
# Fixed overhead of a stage = that stage's R heap peak of the same build (C) at size S, median of its 3 repeats.
# Restated test: (R heap peak - overhead) <= original bound + 2 x noise, noise = the larger range of the 3 repeats at S
# and at the size tested. Usage: python3 -I k150-restate.py <rows.tsv>
import csv
import statistics
import sys

rows = [r for r in csv.DictReader(open(sys.argv[1]), delimiter="\t")
        if r["outcome"] == "completed" and r["valid"] == "yes" and r["build"] == "C"]


def pick(size, stage):
    return [r for r in rows if r["size"] == size and r["stage"] == stage]


def heaps(size, stage):
    return [float(r["heap_peak_mib"]) for r in pick(size, stage)]


def rng(xs):
    return max(xs) - min(xs) if len(xs) >= 3 else 0.0


out = []
for stage in ("train", "import"):
    s = heaps("S", stage)
    print(f"overhead {stage}: C at S, median {statistics.median(s):.1f} MiB, range {rng(s):.1f}, n={len(s)}")
o_train, n_train = statistics.median(heaps("S", "train")), rng(heaps("S", "train"))
for size in ("M", "L"):
    got = pick(size, "train")
    noise = max(n_train, rng(heaps(size, "train")))
    checks = [(float(r["heap_peak_mib"]) - o_train, 2 * float(r["largest_selected_mib"]) + 256) for r in got]
    ok = all(g <= b + 2 * noise for g, b in checks)
    out.append(("1 (restated after results)", size, "pass" if ok else "fail",
                f"overhead {o_train:.1f} MiB (C train at S); growth " + ", ".join(f"{g:.1f}" for g, _ in checks)
                + f" MiB vs bound {checks[0][1]:.1f} MiB + 2 x noise {noise:.1f}"))
o_imp, n_imp = statistics.median(heaps("S", "import")), rng(heaps("S", "import"))
details, ok = [], True
for size in ("P1", "P2", "P3"):
    got = pick(size, "import")
    noise = max(n_imp, rng(heaps(size, "import")))
    for r in got:
        bound = 4 * float(r["row_group_mib"]) + float(r["tables_mib"]) + 64
        ok &= float(r["heap_peak_mib"]) - o_imp <= bound + 2 * noise
    details.append(f"{size}: growth {statistics.median(heaps(size, 'import')) - o_imp:.1f} MiB vs bound {bound:.1f} MiB (noise {noise:.1f})")
out.append(("3a (restated after results)", "P1/P2/P3", "pass" if ok else "fail",
            f"overhead {o_imp:.1f} MiB (C import at S, 4 data sets); " + "; ".join(details)))
for row in out:
    print("\t".join(row))
