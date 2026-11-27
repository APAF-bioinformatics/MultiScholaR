# KIT-275 r2.1 pilot: the scorer's self-test. Standard library only; run as `python3 -I pilot-selftest.py make <work dir>` (writes constructed row files
# under <work>/cases/<case>/rows/ and the predictions <work>/expected.tsv: case, anomaly, design, field, value; and decisions.tsv) and
# `python3 -I pilot-selftest.py check <work dir>` (compares what score-pilot.R wrote to <work>/out/<case>/ with the predictions; exit 1 when any differs).
# Adapted from diag/selftest/diag-selftest.py. No simulation and no random-number generator: every value is an integer number of millionths given by
# arithmetic, so a prediction follows by counting. The rule of R2.3 and the incomplete-cell rule (audit P-1) are written out here, not read from the scorer:
#   an anomaly is ACCEPTED iff every design has 50 valid journeys and at least 45 of them below 0.5; REJECTED iff every design has 50 valid and one has fewer than 45;
#   INCOMPLETE iff any design has fewer than 50 valid journeys. Outcome: V1a accepted -> "V1a accepted"; V1a incomplete -> no decision (V1a); V1a rejected and V1b
#   accepted -> "V1b accepted (V1a rejected)"; V1b incomplete -> no decision (V1b); both rejected -> "neither accepted ...".
import csv, math, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(os.path.dirname(HERE), "pilot-k275.R")
DESIGNS = ["D3", "D4", "D6"]
ANOMALIES = ["V1a", "V1b"]
PER = {"D3": 3, "D4": 4, "D6": 6}
EXPECTED_N = 50
ACCEPT_AT = 45

def columns():
    text = open(SCRIPT).read()
    body = text[text.index("DIAG_COLUMNS <- c(") + len("DIAG_COLUMNS <- c("):text.index("diag_state <- new.env()")]
    cols = re.findall(r'"([A-Za-z0-9_]+)"', body)
    assert len(cols) == len(set(cols)) and "cy_planted" in cols and "anomaly" in cols and "dev_sd_ok" in cols, cols
    return cols
COLS = columns()

def cdf(x, n, p):
    return sum(math.comb(n, k) * p ** k * (1 - p) ** (n - k) for k in range(x + 1))
def cp_upper(x, n, alpha):
    if x >= n: return 1.0
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if cdf(x, n, mid) > alpha: lo = mid
        else: hi = mid
    return (lo + hi) / 2
def cp_lower(x, n, alpha):
    if x <= 0: return 0.0
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if 1 - cdf(x - 1, n, mid) < alpha: lo = mid
        else: hi = mid
    return (lo + hi) / 2
assert abs(cp_upper(0, 50, 0.025) - (1 - 0.025 ** (1 / 50))) < 1e-9
assert abs(cp_lower(50, 50, 0.025) - 0.025 ** (1 / 50)) < 1e-9

def hinges(values):
    x = sorted(values)
    n = len(x)
    n4 = math.floor((n + 3) / 2) / 2
    d = [1, n4, (n + 1) / 2, n + 1 - n4, n]
    return [0.5 * (x[math.floor(v) - 1] + x[math.ceil(v) - 1]) for v in d]

def blank(**kw):
    row = {c: "NA" for c in COLS}
    row.update({"seconds": "1.0", "warnings": "0", "step": "da", "e5_y": "TRUE", "da_calls": "1", "plant_ok": "TRUE", "dev_mean_ok": "TRUE", "dev_sd_ok": "TRUE",
                "n_fail_y": "1", "n_removed_y": "1", "norm_acted_y": "TRUE", "planted_run": "R000", "genuine_run": "R000"})
    row.update({k: str(v) for k, v in kw.items()})
    return row
def micro(v): return "%.6f" % (v / 1e6)

def cell_rows(anomaly, design, below, valid=EXPECTED_N, extra_not_ok=0, boundary_at=None):
    """`valid` ok journeys of which `below` have the planted value below 0.5 (values 0.30 + 0.002 i), the rest above (0.60 + 0.002 i); `extra_not_ok` instrument rows."""
    rows = []
    for i in range(valid):
        t = i + 1
        planted = 300000 + 2000 * i if i < below else 600000 + 2000 * i
        if boundary_at is not None and i == boundary_at: planted = 500000   # exactly 0.5 is not below 0.5
        rows.append(blank(design=design, t=t, r=9600 + t, anomaly=anomaly, status="ok", **{"class": "stopped_before_da"}, n_runs=2 * PER[design],
                          cy_planted=micro(planted), planted_below_05=1 if planted < 500000 else 0, cy_genuine=micro(700000 + 5000 * ((i * 7) % 50))))
    for i in range(extra_not_ok):
        t = valid + i + 1
        rows.append(blank(design=design, t=t, r=9600 + t, anomaly=anomaly, status="instrument", **{"class": "check_failed: e5"}, n_runs=2 * PER[design]))
    return rows

def write_rows(path, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

def predict(case, anomaly, design, rows, expected):
    ok = [r for r in rows if r["status"] == "ok" and r["cy_planted"] != "NA"]
    n_valid = len(ok)
    p = [float(r["cy_planted"]) for r in ok]
    g = [float(r["cy_genuine"]) for r in ok if r["cy_genuine"] != "NA"]
    below = sum(1 for v in p if v < 0.5)
    put = lambda field, value: expected.append((case, anomaly, design, field, value))
    put("n_rows", str(len(rows))); put("n_valid", str(n_valid)); put("n_not_ok", str(sum(1 for r in rows if r["status"] != "ok")))
    put("coverage", "complete" if n_valid == EXPECTED_N else "%d of %d" % (n_valid, EXPECTED_N))
    put("planted_below_05", str(below))
    put("planted_below_share", "%.6f" % (below / n_valid) if n_valid else "NA")
    put("cp_lower", "%.6f" % cp_lower(below, n_valid, 0.025) if n_valid else "NA")
    put("cp_upper", "%.6f" % cp_upper(below, n_valid, 0.025) if n_valid else "NA")
    for kind, vals in (("planted", p), ("genuine", g)):
        five = hinges(vals) if vals else [None] * 5
        for nm, v in zip(["min", "hinge1", "median", "hinge3", "max"], five):
            put("%s_%s" % (kind, nm), "%.6f" % v if v is not None else "NA")
    put("genuine_n", str(len(g)))
    put("genuine_below_05", str(sum(1 for v in g if v < 0.5)))
    return n_valid, below

def verdict(cells):   # cells: {design: (n_valid, below)}
    if any(cells[d][0] != EXPECTED_N for d in DESIGNS): return "incomplete"
    return "accepted" if all(cells[d][1] >= ACCEPT_AT for d in DESIGNS) else "rejected"
def rule_line(anomaly, cells):
    v = verdict(cells)
    text = "%s: %s -> accepted iff all >= %d: %s" % (anomaly, ", ".join("%s %d/%d" % (d, cells[d][1], cells[d][0]) for d in DESIGNS), ACCEPT_AT, v)
    return text + ("" if v != "incomplete" else " (a cell has fewer than 50 valid journeys: an instrument defect, not a rejection)")
def outcome(va, vb):
    if va == "accepted": return "V1a accepted"
    if va == "incomplete": return "no decision: V1a incomplete (instrument defect; re-send once under VAL-54)"
    if vb == "accepted": return "V1b accepted (V1a rejected)"
    if vb == "incomplete": return "no decision: V1b incomplete (instrument defect; re-send once under VAL-54)"
    return "neither accepted: back to the owner (R2.3 (3)); no other candidate, no further pilot"

def make(work):
    expected, decisions, index = [], [], []
    def case(name, spec, part_files=False, stale=False):
        """spec: {anomaly: {design: dict(below=, valid=, extra_not_ok=, boundary_at=)}}"""
        rows_all, cells_by = [], {}
        for a in ANOMALIES:
            cells_by[a] = {}
            for d in DESIGNS:
                kw = spec[a][d]
                rows = cell_rows(a, d, **kw)
                cells_by[a][d] = predict(name, a, d, rows, expected)
                rows_all += rows
        folder = os.path.join(work, "cases", name, "rows")
        if part_files:
            write_rows(os.path.join(folder, "part-PILOT-A.tsv"), rows_all[:len(rows_all) // 2])
            unit = blank(design="D6", t=7, r=9607, anomaly="unit", status="instrument", **{"class": "unit_incomplete_no_result_no_trace: 0 of 2 rows"})
            write_rows(os.path.join(folder, "part-PILOT-B.tsv"), rows_all[len(rows_all) // 2:] + [unit])
        else:
            write_rows(os.path.join(folder, "pilot-rows.tsv"), rows_all)
        if stale:   # a part file beside the assembled file must be ignored
            write_rows(os.path.join(folder, "part-PILOT-D3-9601.tsv"), cell_rows("V1a", "D3", below=0))
        va, vb = verdict(cells_by["V1a"]), verdict(cells_by["V1b"])
        decisions.append((name, "rule_V1a", rule_line("V1a", cells_by["V1a"])))
        decisions.append((name, "rule_V1b", rule_line("V1b", cells_by["V1b"])))
        decisions.append((name, "outcome", outcome(va, vb)))
        index.append((name, "0"))
    cell = lambda below, **kw: dict(below=below, **kw)
    accepted = {d: cell(b) for d, b in zip(DESIGNS, (50, 47, 45))}                  # exactly 45 accepts
    rejected = {d: cell(b) for d, b in zip(DESIGNS, (50, 50, 44))}                  # one cell at 44 rejects
    weak = {d: cell(b) for d, b in zip(DESIGNS, (30, 20, 5))}
    case("a-accepted", {"V1a": accepted, "V1b": weak}, stale=True)                  # V1a decides whatever V1b holds; a stale part file is ignored
    case("b-fallback", {"V1a": rejected, "V1b": accepted}, part_files=True)         # V1b accepted, part files and a unit row
    case("c-neither", {"V1a": rejected, "V1b": weak})
    case("d-incomplete-a", {"V1a": {"D3": cell(50), "D4": cell(40, valid=40, extra_not_ok=10), "D6": cell(50)}, "V1b": accepted})   # P-1: not a rejection
    case("e-incomplete-b", {"V1a": rejected, "V1b": {"D3": cell(50), "D4": cell(50), "D6": cell(49, valid=49, extra_not_ok=1)}})
    case("f-accepted-b-incomplete", {"V1a": accepted, "V1b": {"D3": cell(50), "D4": cell(50), "D6": cell(0, valid=0, extra_not_ok=50)}})   # V1a complete and accepted decides
    case("g-boundary", {"V1a": {"D3": cell(47, boundary_at=0), "D4": cell(46, boundary_at=1), "D6": cell(46, boundary_at=2)}, "V1b": weak})   # a value of exactly 0.5 replaces a below value: not below, so 46, 45, 45 (all accepted); a scorer using <= would count 47, 46, 46
    os.makedirs(os.path.join(work, "cases", "norows", "rows"), exist_ok=True)
    index.append(("norows", "1"))
    for name, fname, header in (("expected.tsv", None, ["case", "anomaly", "design", "field", "value"]),):
        with open(os.path.join(work, name), "w", newline="") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(header); w.writerows(expected)
    with open(os.path.join(work, "decisions.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["case", "item", "value"]); w.writerows(decisions)
    with open(os.path.join(work, "index.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["case", "exit"]); w.writerows(index)
    print("pilot-selftest: %d cases, %d predictions, %d decision lines" % (len(index), len(expected), len(decisions)))

def check(work):
    expected = list(csv.DictReader(open(os.path.join(work, "expected.tsv")), delimiter="\t"))
    decisions = list(csv.DictReader(open(os.path.join(work, "decisions.tsv")), delimiter="\t"))
    index = list(csv.DictReader(open(os.path.join(work, "index.tsv")), delimiter="\t"))
    bad = total = 0
    summaries, decided = {}, {}
    for row in index:
        out = os.path.join(work, "out", row["case"])
        code = open(os.path.join(out, "exit.txt")).read().strip() if os.path.exists(os.path.join(out, "exit.txt")) else None
        total += 1
        if code is None or (code == "0") != (row["exit"] == "0"):
            bad += 1
            print("DIFFERS %s exit: got %s, expected %s" % (row["case"], code, "0" if row["exit"] == "0" else "non-zero"))
        s, d = os.path.join(out, "pilot-summary.tsv"), os.path.join(out, "pilot-decision.tsv")
        if row["exit"] == "0":
            if not (os.path.exists(s) and os.path.exists(d)):
                bad += 1
                print("DIFFERS %s: output absent" % row["case"])
            else:
                summaries[row["case"]] = {(r["anomaly"], r["design"]): r for r in csv.DictReader(open(s), delimiter="\t")}
                decided[row["case"]] = {r["item"]: r["value"] for r in csv.DictReader(open(d), delimiter="\t")}
        elif os.path.exists(s) or os.path.exists(d):
            bad += 1
            print("DIFFERS %s: output written for an empty input" % row["case"])
    for e in expected:
        total += 1
        got = summaries.get(e["case"], {}).get((e["anomaly"], e["design"]), {}).get(e["field"])
        same = got == e["value"]
        if not same and got is not None:
            try: same = abs(float(got) - float(e["value"])) <= 1.5e-6
            except ValueError: same = False
        if not same:
            bad += 1
            print("DIFFERS %s %s %s %s: got %s, expected %s" % (e["case"], e["anomaly"], e["design"], e["field"], got, e["value"]))
    for e in decisions:
        total += 1
        got = decided.get(e["case"], {}).get(e["item"])
        if got != e["value"]:
            bad += 1
            print("DIFFERS %s %s: got %s, expected %s" % (e["case"], e["item"], got, e["value"]))
    print("pilot-selftest: %d predictions checked, %d differ" % (total, bad))
    return 1 if bad else 0

if __name__ == "__main__":
    mode, work = sys.argv[1], sys.argv[2]
    if mode == "make": make(work)
    elif mode == "check": sys.exit(check(work))
    else: sys.exit("usage: pilot-selftest.py make|check <work dir>")
