# KIT-281: the refusal of a counted request whose R differs from rule G3's (ticket sections 7 and 8; local, standard library only).
# Usage: python3 -I check-request-k281.py --g3 <g3-rule.tsv> [--count <g3-count.tsv>] <request.json> [<request.json> ...]
# Exit 0: the requests together run, for every cell A to D, exactly r = 401 to 400 + R once (R from the table, which --count re-derives with
# Fractions: N-hat = 200 x the sum of the twelve counts; R = 200 when N-hat >= 276,220, else min(571, max(200, ceiling(200 x 276,220 / N-hat)))),
# in jobs of at most 200 replicates with a timeout of at most 75 minutes, each a sim-k281.R script, and every request CITES the table: the first 16 hex
# digits of the g3-rule.tsv's SHA-256 appear in its `requested_by`.
import csv, hashlib, json, sys
from fractions import Fraction

SIM = "docs/qa/workflow-execution/kit/kit-281/sim-k281.R"
CELLS = ["A", "B", "C", "D"]
args = sys.argv[1:]
g3_path, count_path, requests = None, None, []
while args:
    a = args.pop(0)
    if a == "--g3":
        g3_path = args.pop(0)
    elif a == "--count":
        count_path = args.pop(0)
    else:
        requests.append(a)
problems = []
rows = list(csv.DictReader(open(g3_path), delimiter="\t"))
if sorted(r["cell"] for r in rows) != CELLS:
    problems.append("the G3 table does not hold cells A to D once each")
r_table = {int(r["r_cell"]) for r in rows}
first_table = {int(r["first_r"]) for r in rows}
last_table = {int(r["last_r"]) for r in rows}
if len(r_table) != 1 or first_table != {401} or len(last_table) != 1 or last_table != {400 + next(iter(r_table))}:
    problems.append("the G3 table holds one R for every cell with r = 401 to 400 + R")
R = next(iter(r_table)) if len(r_table) == 1 else None


def rule(n_hat):
    n_hat = Fraction(n_hat)
    if n_hat >= 276220:
        return 200
    if n_hat == 0:
        return 571
    q = Fraction(200 * 276220) / n_hat
    return int(min(571, max(200, -((-q.numerator) // q.denominator))))


if count_path:
    counts = list(csv.DictReader(open(count_path), delimiter="\t"))
    if len(counts) != 12:
        problems.append("the count table does not hold twelve settings")
    n_hat = 200 * sum(Fraction(r["m_s"]) for r in counts)
    if R is not None and rule(n_hat) != R:
        problems.append("the table's R %d differs from rule G3's %d for N-hat %s" % (R, rule(n_hat), n_hat))
seen = {c: [] for c in CELLS}
cite = hashlib.sha256(open(g3_path, "rb").read()).hexdigest()[:16]
for path in requests:
    q = json.load(open(path))
    if cite not in q.get("requested_by", ""):
        problems.append("%s: requested_by does not cite the G3 table (its SHA-256 starts %s)" % (path, cite))
    for key, want in (("base_commit", "none"), ("tests", []), ("blas", "reference"), ("extra_packages", [])):
        if key in q and q[key] != want:
            problems.append("%s: %s must be %r" % (path, key, want))
    scripts = q.get("scripts", [])
    if len(scripts) < 1:
        problems.append("%s: no script" % path)
    for s in scripts:
        if s.get("path") != SIM:
            problems.append("%s: a script is not sim-k281.R" % path)
            continue
        first, last = int(s["args"][1]), int(s["args"][2])
        cells = [c for c in s.get("env", {}).get("K211_CELLS", "A,B,C,D").split(",") if c]
        if last - first + 1 > 200:
            problems.append("%s: a job of %d replicates (at most 200)" % (path, last - first + 1))
        if not (1 <= int(s.get("timeout_minutes", 0)) <= 75):
            problems.append("%s: timeout_minutes missing or above 75" % path)
        for c in cells:
            seen[c].append((first, last))
for c in CELLS:
    covered = sorted(r for first, last in seen[c] for r in range(first, last + 1))
    if R is not None and covered != list(range(401, 401 + R)):
        problems.append("cell %s: the requests run %d replicate(s), %s; the rule's range is r = 401 to %d (R = %d), each once" % (
            c, len(covered), ("r = %d to %d" % (covered[0], covered[-1])) if covered else "none", 400 + R, R))
for p in problems:
    print("REFUSED:", p)
print("check-request:", "refused" if problems else "accepted (every cell r = 401 to %d once)" % (400 + R))
sys.exit(1 if problems else 0)
