"""KIT-262 E0 (local, standard library only): the copied blocks of m1-k262.R and m2-k262.R equal their sources
except the marked lines. Run from the repository root: python3 -I docs/qa/workflow-execution/kit/kit-262/e0-k262.py
Sources: the cancor-once test file (tests/testthat/test-workflow-omics-ruv-cancor-once.R), KIT-190's effect.R
(e2PickCells, e2Values) and KIT-196's sim-k188-69fd01bc.R (solveLod). Prints one line per check and exits 1 on a failure."""
import sys
from pathlib import Path

root = Path(".")
here = root / "docs/qa/workflow-execution/kit/kit-262"
kit190 = root / "docs/qa/workflow-execution/kit/kit-190/effect.R"
test = root / "tests/testthat/test-workflow-omics-ruv-cancor-once.R"
sim188 = root / "docs/qa/workflow-execution/kit/kit-196/sim-k188-69fd01bc.R"


def lines(path):
    return path.read_text().split("\n")


def block(text, start_marker):
    """The lines between a '# ---- copy: ...' marker line containing start_marker and the next end marker."""
    out, on = [], False
    for line in text:
        if line.startswith("# ---- copy:") and start_marker in line:
            on = True
            continue
        if on and line.startswith("# ---- end of the copy ----"):
            return out
        if on:
            out.append(line)
    raise SystemExit("copy block not found: " + start_marker)


def nonblank(seq):
    return [l for l in seq if l.strip()]


def function_span(text, header):
    start = next(i for i, l in enumerate(text) if l.startswith(header))
    end = start
    while text[end] != "}":
        end += 1
    return text[start : end + 1]


failures = 0


def check(name, ok, detail=""):
    global failures
    print(("PASS " if ok else "FAIL ") + name + (" " + detail if detail else ""))
    if not ok:
        failures += 1


m1, m2 = lines(here / "m1-k262.R"), lines(here / "m2-k262.R")
t = lines(test)

# m1: the fixture block equals the test file from `ruvAssay <- function` to the end of `ruvComparable`
a = next(i for i, l in enumerate(t) if l.startswith("ruvAssay <- function"))
b = next(i for i, l in enumerate(t) if l.startswith("ruvComparable <- function"))
while t[b] != "}":
    b += 1
source = nonblank(t[a : b + 1])
copy = nonblank(block(m1, "test-workflow-omics-ruv-cancor-once.R:5-28"))
check("m1 fixture block = test file", source == copy, "%d lines" % len(source))

# m2: ruvComparable
b0 = next(i for i, l in enumerate(t) if l.startswith("ruvComparable <- function"))
b1 = b0
while t[b1] != "}":
    b1 += 1
check("m2 ruvComparable = test file", nonblank(t[b0 : b1 + 1]) == nonblank(block(m2, "test-workflow-omics-ruv-cancor-once.R:43-54")))

# m2: e2PickCells, solveLod
e = lines(kit190)
check("m2 e2PickCells = kit-190 effect.R", nonblank(function_span(e, "e2PickCells <- function")) == nonblank(block(m2, "e2PickCells")))
s = lines(sim188)
i = next(i for i, l in enumerate(s) if l.startswith("solveLod <- function"))
check("m2 solveLod = kit-196 sim-k188", nonblank(s[i : i + 3]) == nonblank(block(m2, "mask solver")))

# m2: m2Values against e2Values: every unmarked line is a line of e2Values (stripped); marked lines carry '# [X]'
src = set(l.strip() for l in function_span(e, "e2Values <- function"))
mine = function_span(m2, "m2Values <- function")
unmarked_new = [l for l in mine if l.strip() and l.strip() not in src and "# [" not in l]
check("m2Values: every line not in e2Values is marked", not unmarked_new, "; ".join(unmarked_new))
kept = [l for l in mine if l.strip() in src]
check("m2Values: %d lines kept verbatim, %d marked" % (len(kept), len([l for l in mine if "# [" in l])), True)

# the order of the kept lines follows e2Values (the stream of draws is the same up to the declared changes)
positions = []
src_lines = [l.strip() for l in function_span(e, "e2Values <- function")]
cursor = 0
ordered = True
for l in kept:
    s_ = l.strip()
    try:
        cursor = src_lines.index(s_, cursor)
    except ValueError:
        # lines that repeat or move are reported
        ordered = False
        print("  order: " + s_)
check("m2Values: kept lines in e2Values' order", ordered)

sys.exit(1 if failures else 0)
