# KIT-275 diagnostic: the scorer's self-test. Standard library only; run as `python3 -I diag-selftest.py make <work dir>` (writes constructed row
# files under <work>/cases/<case>/rows/ and the predictions <work>/expected.tsv: case, design, norm, field, value) and
# `python3 -I diag-selftest.py check <work dir>` (compares what score-diag.R wrote to <work>/out/<case>/ with the predictions; exit 1 when any differs).
# No simulation and no random-number generator: every value is an integer number of millionths given by arithmetic below, so a prediction follows
# by counting. The exact Clopper-Pearson limits are solved by bisection on the binomial distribution; Tukey's hinges are written out from their
# definition (the rules of score-diag.R are written out again here, not read from it).
import csv, math, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(os.path.dirname(HERE), "diag-k275.R")
DESIGNS = ["D3", "D4", "D6"]
NORMS = ["cyclicloess", "quantile", "none"]
PER = {"D3": 3, "D4": 4, "D6": 6}
EXPECTED_N = 50

def columns():
    text = open(SCRIPT).read()
    body = text[text.index("DIAG_COLUMNS <- c(") + len("DIAG_COLUMNS <- c("):text.index("diag_state <- new.env()")]
    cols = re.findall(r'"([A-Za-z0-9_]+)"', body)
    assert len(cols) == len(set(cols)) and "cy_planted" in cols and "norm" in cols, cols
    return cols
COLS = columns()

def cdf(x, n, p):
    return sum(math.comb(n, k) * p ** k * (1 - p) ** (n - k) for k in range(x + 1))
def cp_upper(x, n, alpha):   # the p at which P(X <= x) = alpha
    if x >= n: return 1.0
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if cdf(x, n, mid) > alpha: lo = mid
        else: hi = mid
    return (lo + hi) / 2
def cp_lower(x, n, alpha):   # the p at which P(X >= x) = alpha
    if x <= 0: return 0.0
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if 1 - cdf(x - 1, n, mid) < alpha: lo = mid
        else: hi = mid
    return (lo + hi) / 2
# closed forms at the ends, as a check of the bisection
assert abs(cp_upper(0, 50, 0.025) - (1 - 0.025 ** (1 / 50))) < 1e-9
assert abs(cp_lower(50, 50, 0.025) - 0.025 ** (1 / 50)) < 1e-9

def hinges(values):   # Tukey's five numbers, from the definition
    x = sorted(values)
    n = len(x)
    n4 = math.floor((n + 3) / 2) / 2
    d = [1, n4, (n + 1) / 2, n + 1 - n4, n]
    return [0.5 * (x[math.floor(v) - 1] + x[math.ceil(v) - 1]) for v in d]

def blank(**kw):
    row = {c: "NA" for c in COLS}
    row.update({"seconds": "1.0", "warnings": "0", "step": "da", "e5_y": "TRUE", "da_calls": "1", "plant_ok": "TRUE", "n_fail_y": "1", "n_removed_y": "1",
                "norm_acted_y": "TRUE", "planted_run": "R000", "genuine_run": "R000"})
    row.update({k: str(v) for k, v in kw.items()})
    return row
def micro(v): return "%.6f" % (v / 1e6)

def planted_micro(k, i):   # cell k (0 to 8), data set index i (0 to 49)
    base = 50000 + 40000 * k
    step = min(4000 + 1000 * k, (990000 - base) // 49)
    return base + step * i
def genuine_micro(i): return 700000 + 5000 * ((i * 7) % 50)

def cell_rows(design, norm, k, spec=None):
    spec = spec or {}
    rows = []
    for i in range(spec.get("n", EXPECTED_N)):
        t = i + 1
        kw = dict(design=design, t=t, r=9500 + t, norm=norm, status="ok", **{"class": "stopped_before_da"}, n_runs=2 * PER[design])
        kw["cy_planted"] = micro(planted_micro(k, i))
        kw["cy_genuine"] = micro(genuine_micro(i))
        kw["planted_below_05"] = 1 if planted_micro(k, i) < 500000 else 0
        if spec.get("exact_boundary") and i == 25:
            kw["cy_planted"] = "0.500000"   # exactly 0.5 is not below 0.5
            kw["planted_below_05"] = 0
        if "all" in spec:
            kw["cy_planted"] = micro(spec["all"] + 1000 * (i % 7))
            kw["planted_below_05"] = 1 if spec["all"] + 1000 * (i % 7) < 500000 else 0
        for kind, idxs in spec.get("not_ok", {}).items():
            if i in idxs:
                kw.update(status=kind, cy_planted="NA", cy_genuine="NA", planted_below_05="NA", step="owner" if kind == "crashed" else "norm")
        if i in spec.get("ok_but_na", []):
            kw.update(cy_planted="NA", planted_below_05="NA")
        if i in spec.get("genuine_na", []):
            kw["cy_genuine"] = "NA"
        rows.append(blank(**kw))
    return rows

def write_rows(path, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

def predict(case, design, norm, rows, expected):
    n_rows = len(rows)
    ok = [r for r in rows if r["status"] == "ok" and r["cy_planted"] != "NA"]
    n_valid = len(ok)
    p = [float(r["cy_planted"]) for r in ok]
    g = [float(r["cy_genuine"]) for r in ok if r["cy_genuine"] != "NA"]
    below = sum(1 for r in ok if float(r["cy_planted"]) < 0.5)
    put = lambda field, value: expected.append((case, design, norm, field, value))
    put("n_rows", str(n_rows)); put("n_valid", str(n_valid)); put("n_not_ok", str(sum(1 for r in rows if r["status"] != "ok")))
    put("coverage", "complete" if n_valid == EXPECTED_N else "%d of %d" % (n_valid, EXPECTED_N))
    put("planted_below_05", str(below))
    put("planted_below_share", "%.6f" % (below / n_valid) if n_valid else "NA")
    put("cp_lower", "%.6f" % cp_lower(below, n_valid, 0.025) if n_valid else "NA")
    put("cp_upper", "%.6f" % cp_upper(below, n_valid, 0.025) if n_valid else "NA")
    names = ["min", "hinge1", "median", "hinge3", "max"]
    for kind, vals in (("planted", p), ("genuine", g)):
        five = hinges(vals) if vals else [None] * 5
        for nm, v in zip(names, five):
            put("%s_%s" % (kind, nm), "%.6f" % v if v is not None else "NA")
    put("genuine_n", str(len(g)))
    put("genuine_below_05", str(sum(1 for v in g if v < 0.5)))

def make(work):
    expected, index = [], []
    # case full: nine cells of 50 valid rows (one holds exactly 0.5; one genuine value NA in an ok row), the assembled file, plus a stale part file that must be ignored
    rows_all = []
    for d, design in enumerate(DESIGNS):
        for m, norm in enumerate(NORMS):
            k = 3 * d + m
            spec = {"exact_boundary": True} if k == 4 else ({"genuine_na": [3, 17]} if k == 7 else {})
            rows = cell_rows(design, norm, k, spec)
            predict("full", design, norm, rows, expected)
            rows_all += rows
    write_rows(os.path.join(work, "cases", "full", "rows", "diag-rows.tsv"), rows_all)
    write_rows(os.path.join(work, "cases", "full", "rows", "part-DIAG-D3-9501.tsv"), cell_rows("D3", "cyclicloess", 8, {"all": 100000}))
    index.append(("full", "0"))
    # case edges: part files only (two files); all below (50 of 50), none below (0 of 50), not-ok rows, an ok row with no planted value, a cell with no valid row,
    # a short cell, a unit row, an instrument row
    cells = {}
    cells[("D3", "cyclicloess")] = cell_rows("D3", "cyclicloess", 0, {"all": 100000})                  # 50 of 50 below
    cells[("D3", "quantile")] = cell_rows("D3", "quantile", 0, {"all": 600000})                        # 0 of 50 below
    cells[("D3", "none")] = cell_rows("D3", "none", 1, {"not_ok": {"instrument": [2, 3, 4, 5, 6], "crashed": [10, 11]}, "ok_but_na": [20]})   # 42 valid
    cells[("D4", "cyclicloess")] = cell_rows("D4", "cyclicloess", 2, {"n": 30})                        # short: 30 rows
    cells[("D4", "quantile")] = cell_rows("D4", "quantile", 3, {"not_ok": {"instrument": list(range(50))}})   # no valid row
    cells[("D4", "none")] = cell_rows("D4", "none", 4)
    cells[("D6", "cyclicloess")] = cell_rows("D6", "cyclicloess", 5, {"all": 495000})                  # values 0.495 to 0.501: the boundary 0.5 is not below
    cells[("D6", "quantile")] = cell_rows("D6", "quantile", 6)
    cells[("D6", "none")] = cell_rows("D6", "none", 7, {"genuine_na": list(range(50))})                # no genuine value at all
    unit = blank(design="D6", t=7, r=9507, norm="unit", status="instrument", **{"class": "unit_incomplete_no_result_no_trace: 0 of 3 rows"})
    first, second = [], []
    for n, ((design, norm), rows) in enumerate(cells.items()):
        predict("edges", design, norm, rows, expected)
        (first if n < 5 else second).extend(rows)
    second.append(unit)
    write_rows(os.path.join(work, "cases", "edges", "rows", "part-DIAG-A.tsv"), first)
    write_rows(os.path.join(work, "cases", "edges", "rows", "part-DIAG-B.tsv"), second)
    index.append(("edges", "0"))
    # case norows: an empty rows directory: the scorer stops with a message and a non-zero exit
    os.makedirs(os.path.join(work, "cases", "norows", "rows"), exist_ok=True)
    index.append(("norows", "1"))
    with open(os.path.join(work, "expected.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["case", "design", "norm", "field", "value"])
        w.writerows(expected)
    with open(os.path.join(work, "index.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["case", "exit"])
        w.writerows(index)
    print("diag-selftest: %d cases, %d predictions" % (len(index), len(expected)))

def check(work):
    expected = list(csv.DictReader(open(os.path.join(work, "expected.tsv")), delimiter="\t"))
    index = list(csv.DictReader(open(os.path.join(work, "index.tsv")), delimiter="\t"))
    bad = 0
    total = 0
    summaries = {}
    for row in index:
        out = os.path.join(work, "out", row["case"])
        code = open(os.path.join(out, "exit.txt")).read().strip() if os.path.exists(os.path.join(out, "exit.txt")) else None
        total += 1
        if code is None or (code == "0") != (row["exit"] == "0"):
            bad += 1
            print("DIFFERS %s exit: got %s, expected %s" % (row["case"], code, "0" if row["exit"] == "0" else "non-zero"))
        path = os.path.join(out, "diag-summary.tsv")
        if row["exit"] == "0":
            if not os.path.exists(path):
                bad += 1
                print("DIFFERS %s: diag-summary.tsv absent" % row["case"])
            else:
                summaries[row["case"]] = {(r["design"], r["norm"]): r for r in csv.DictReader(open(path), delimiter="\t")}
        elif os.path.exists(path):
            bad += 1
            print("DIFFERS %s: diag-summary.tsv written for an empty input" % row["case"])
    for e in expected:
        total += 1
        got = summaries.get(e["case"], {}).get((e["design"], e["norm"]), {}).get(e["field"])
        want = e["value"]
        same = got == want
        if not same and got is not None:
            try: same = abs(float(got) - float(want)) <= 1.5e-6
            except ValueError: same = False
        if not same:
            bad += 1
            print("DIFFERS %s %s %s %s: got %s, expected %s" % (e["case"], e["design"], e["norm"], e["field"], got, want))
    print("diag-selftest: %d predictions checked, %d differ" % (total, bad))
    return 1 if bad else 0

if __name__ == "__main__":
    mode, work = sys.argv[1], sys.argv[2]
    if mode == "make": make(work)
    elif mode == "check": sys.exit(check(work))
    else: sys.exit("usage: diag-selftest.py make|check <work dir>")
