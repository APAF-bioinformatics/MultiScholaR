# KIT-281 scorer self-test: compares what score-k281.R wrote for each constructed case with the values predicted in expected.tsv
# (written by make-selftest.py in Python before the scorer ran), and runs the G3 request checker on its fixture requests.
# Standard library only. Usage: python3 -I check-selftest.py <selftest dir> <scored dir>. Exit 1 on any difference.
# Prints one "SELFTEST: <checks> checks, <differences> differ" line last.
import csv, os, subprocess, sys

KEYS = {"decision.tsv": ["id"], "falsifier.tsv": ["falsifier"], "positive-control.tsv": ["stat"], "settings.tsv": ["scope", "cell", "mask", "stat"],
        "designs.tsv": ["scope", "cell", "mask", "stat"], "pooled.tsv": ["scope", "cell", "mask", "stat"], "labels-settings.tsv": ["cell", "mask"],
        "checks.tsv": ["check"], "qmethod-fallback.tsv": ["cell", "mask", "arm"], "replicates.tsv": ["cell"], "outcome.tsv": []}
here, scored = sys.argv[1], sys.argv[2]


def table(case, name):
    with open(os.path.join(scored, case, name), newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def same(got, want):
    if want.startswith(">"):
        try:
            return float(got) > float(want[1:])
        except ValueError:
            return False
    if got == want:
        return True
    try:
        g, w = float(got), float(want)
        return abs(g - w) <= 1e-12 * max(1.0, abs(w))
    except ValueError:
        return got.upper() == want.upper() and got.upper() in ("TRUE", "FALSE")


bad = total = 0
for c in csv.DictReader(open(os.path.join(here, "cases.tsv")), delimiter="\t"):
    case = c["case"]
    log_path = os.path.join(scored, case, "run.log")
    log = open(log_path).read() if os.path.exists(log_path) else ""
    code_path = os.path.join(scored, case, "exit.txt")
    code = open(code_path).read().strip() if os.path.exists(code_path) else "missing"
    want_stop = c["expect"].startswith("stop:")
    ok = (code == "1" and c["expect"][5:] in log) if want_stop else code == "0"
    total += 1
    print("%-18s %s (exit %s; expected %s)" % (case, "as predicted" if ok else "DIFFERS", code, c["expect"]))
    bad += 0 if ok else 1
for line in csv.DictReader(open(os.path.join(here, "expected.tsv")), delimiter="\t"):
    total += 1
    case, f, key, col, want = line["case"], line["file"], line["key"], line["column"], line["value"]
    try:
        rows = table(case, f)
        parts = key.split("|") if key else []
        hit = [r for r in rows if [r[k] for k in KEYS[f]] == parts]
        got = hit[0][col] if len(hit) == 1 else "(%d rows)" % len(hit)
        ok = len(hit) == 1 and same(got, want)
    except Exception as e:  # a missing file or column is a difference
        got, ok = "error: %s" % e, False
    if not ok:
        bad += 1
        print("DIFFERS %s %s %s %s: got %s, predicted %s" % (case, f, key, col, got, want))
# the G3 request checker on its fixtures
req = os.path.join(here, "requests")
for line in csv.DictReader(open(os.path.join(here, "requests.tsv")), delimiter="\t"):
    total += 1
    done = subprocess.run([sys.executable, "-I", os.path.join(here, "..", "check-request-k281.py"), "--g3", os.path.join(req, line["g3_table"]),
                           os.path.join(req, line["request"])], capture_output=True, text=True)
    ok = str(done.returncode) == line["exit"]
    print("%-18s %s (check-request exit %d; expected %s)" % (line["request"], "as predicted" if ok else "DIFFERS", done.returncode, line["exit"]))
    bad += 0 if ok else 1
print("SELFTEST: %d checks, %d differ" % (total, bad))
sys.exit(1 if bad else 0)
