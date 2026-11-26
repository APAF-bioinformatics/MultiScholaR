# KIT-275 E4: compares what the scorer wrote for each constructed case with the predictions made in make-selftest.py.
# Usage: python3 -I check-selftest.py <work dir>      (the work dir holds expected.tsv, index.tsv and out/<case>/ written by run-selftest.R)
# Prints one line per differing prediction and a summary; exits 1 when any differs. Standard library only.
import csv, os, sys

work = sys.argv[1]

def table(path):
    if not os.path.exists(path): return None
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def lookup(case, key):
    out = os.path.join(work, "out", case)
    parts = key.split("|")
    kind = parts[0]
    if kind == "exit":
        p = os.path.join(out, "exit.txt")
        return open(p).read().strip() if os.path.exists(p) else None
    if kind == "stderr":
        p = os.path.join(out, "stderr.txt")
        return open(p).read() if os.path.exists(p) else None
    if kind == "decision":
        rows = table(os.path.join(out, "decision.tsv"))
        if rows is None: return None
        for r in rows:
            if r["item"] == parts[1]: return r["value"]
        return None
    if kind == "cell":
        rows = table(os.path.join(out, "cells.tsv"))
        if rows is None: return None
        ident, tier, design, version, arm, field = parts[1:7] if len(parts) == 7 else (None,) * 6
        if len(parts) == 7:
            for r in rows:
                if (r["id"], r["tier"], r["design"], r["version"], r["arm"]) == (ident, tier, design, version, arm):
                    return r[field]
        return None
    if kind == "labels":
        rows = table(os.path.join(out, "labels.tsv"))
        if rows is None: return None
        for r in rows:
            if r["id"] == parts[1]: return r[parts[2]]
        return None
    if kind == "scorer":
        rows = table(os.path.join(out, "scorer-%s.tsv" % case))
        if rows is None: return None
        for r in rows:
            if r["item"] == parts[1]: return r[parts[2]]
        return None
    raise ValueError(key)

expected = table(os.path.join(work, "expected.tsv"))
bad = 0
for e in expected:
    got = lookup(e["case"], e["key"])
    kind, want, tol = e["kind"], e["expected"], e["tolerance"]
    if kind == "num":
        ok = False
        if got not in (None, "NA", ""):
            try:
                g, w = float(got), float(want)
                t = float(tol or "1e-8")
                ok = abs(g - w) <= t * max(1.0, abs(w))
            except ValueError:
                ok = False
    elif kind == "between":
        lo_, hi_ = (float(x) for x in want.split(":"))
        try: ok = lo_ - 1e-9 <= float(got) <= hi_ + 1e-9
        except (TypeError, ValueError): ok = False
    elif kind == "str":
        ok = got == want
    elif kind == "substr":
        ok = got is not None and want in got
    elif kind == "exit":
        ok = got is not None and (got != "0" if want == "nonzero" else got == want)
    else:
        raise ValueError(kind)
    if not ok:
        bad += 1
        print("DIFFERS %s %s: got %r, predicted %r" % (e["case"], e["key"], got, want))
print("selftest: %d predictions, %d differ" % (len(expected), bad))
sys.exit(1 if bad else 0)
