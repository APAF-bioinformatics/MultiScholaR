"""KIT-262 re-measurement E0 (local, stdlib only; run from the repository root): the generator `m2Values` of the
current m2-k262.R differs from the one at ade127d1 (the M2 binding) only in lines marked [G]; and the lines of the old
generator absent from the new one are exactly the old signature and `values[ones] <- 1` (now `<<-` inside the [G] helper).
With floor_order = "before" (the default) the code path is the old one: floor values chosen first, then the mask, in
the same random-number order (e2PickCells, then the mask draws)."""
import re, subprocess, sys
path = "docs/qa/workflow-execution/kit/kit-262/m2-k262.R"
old = subprocess.run(["git", "show", "ade127d110e166af7c7c8b4d19af7943df251d4e:" + path], capture_output=True, text=True, check=True).stdout.split("\n")
new = open(path).read().split("\n")
def span(t):
    a = next(i for i, l in enumerate(t) if l.startswith("m2Values <- function")); b = a
    while t[b] != "}": b += 1
    return t[a:b + 1]
def code(l): return re.sub(r"\s*#.*$", "", l).strip()
o, n = span(old), span(new)
oc = set(code(l) for l in o if code(l)); nc = set(code(l) for l in n if code(l))
bad = [l for l in n if code(l) and code(l) not in oc and "[G]" not in l]
gone = [l.strip() for l in o if code(l) and code(l) not in nc]
print("new lines not in the old generator and not marked [G]:", bad)
print("old lines absent from the new generator:", gone)
expected_gone = {"m2Values <- function(seed, n, g, n_features, mech, rate) {", "values[ones] <- 1"}
ok = not bad and set(code(l) for l in gone) <= expected_gone
# the order of the unchanged lines is unchanged
common = [code(l) for l in n if code(l) in oc and code(l) not in ("}", "{")]
cursor = 0; ordered = True
oc_list = [code(l) for l in o if code(l)]
for c in common:
    try: cursor = oc_list.index(c, cursor)
    except ValueError: ordered = False
print("unchanged lines in the old order:", ordered)
print("E0-order:", "PASS" if ok and ordered else "FAIL")
sys.exit(0 if ok and ordered else 1)
