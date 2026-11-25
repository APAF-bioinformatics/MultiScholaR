# E16: compares what score-k188.R wrote for each constructed case with the values predicted in expected.tsv (written by make-e16.py
# before the scorer ran). Standard library only. Usage: python3 -I check-e16.py <e16 dir> <scored dir>. Exit 1 on any difference.
import csv, os, sys

KEYS = {"conditions.tsv": ["cell", "job", "scenario", "mask", "arm", "estimand"], "phi.tsv": ["cell", "job", "scenario", "mask", "arm"],
        "verdicts.tsv": ["group", "arm"], "statuses.tsv": ["part", "cell", "job", "arm", "arm_scale", "status", "class"]}
here, scored = sys.argv[1], sys.argv[2]
def table(case, name):
    path = os.path.join(scored, case, name)
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))
def same(got, want):
    if got == want: return True
    try: return abs(float(got) - float(want)) <= 1e-9 * max(1.0, abs(float(want)))
    except ValueError: return got.upper() == want.upper() and got.upper() in ("TRUE", "FALSE")
bad = 0; total = 0
cases = list(csv.DictReader(open(os.path.join(here, "cases.tsv")), delimiter="\t"))
for c in cases:
    log = open(os.path.join(scored, c["case"], "run.log")).read() if os.path.exists(os.path.join(scored, c["case"], "run.log")) else ""
    code = open(os.path.join(scored, c["case"], "exit.txt")).read().strip()
    want_stop = c["expect"].startswith("stop:")
    ok = (code != "0" and c["expect"][5:] in log) if want_stop else code == "0"
    total += 1
    print("%-8s %s (exit %s; expected %s)" % (c["case"], "as predicted" if ok else "DIFFERS", code, c["expect"]))
    bad += 0 if ok else 1
for line in csv.DictReader(open(os.path.join(here, "expected.tsv")), delimiter="\t"):
    total += 1
    case, f, key, col, want = line["case"], line["file"], line["key"], line["column"], line["value"]
    try:
        if f == "outcome.txt":
            text = open(os.path.join(scored, case, f)).read().strip()
            ok = text.startswith(want) if col == "startswith" else want in text
            got = text
        elif key.startswith("absent:"):
            got = str(sum(1 for r in table(case, "conditions.tsv") if r["cell"] == key[7:] and r["role"] == "deciding"))
            ok = got == want
        else:
            parts = key.split("|")
            hit = [r for r in table(case, f) if [r[k] for k in KEYS[f]] == parts]
            got = hit[0][col] if len(hit) == 1 else "(%d rows)" % len(hit)
            ok = len(hit) == 1 and same(got, want)
    except Exception as e:  # a missing file or column is a difference
        got, ok = "error: %s" % e, False
    if not ok:
        bad += 1
        print("DIFFERS %s %s %s %s: got %s, predicted %s" % (case, f, key, col, got, want))
print("E16: %d checks, %d differ" % (total, bad))
sys.exit(1 if bad else 0)
