# KIT-275 r2.1 (ticket 8.1, 8.3 E0; audit Q1, K2-03, K2-07): the diff of every changed script against r1's audited versions, by hunk, each hunk tagged with the ticket
# section or finding it implements. Writes e0-r2.txt. A hunk is tagged when a line it adds carries `K275-r2: <tag>`, or when it is a pure label substitution
# (V1g -> V1a, the label check of 8.1 / K2-07: the old and new lines differ only by that token). ANYTHING ELSE IS AN UNTAGGED HUNK AND IS NOT ALLOWED: the script exits 1.
# Standard library and git only: `python3 -I e0-r2.py` from this folder (needs git history with 577a6f7b; r1's counted versions equal the audited d6b95a05).
import difflib, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
PIN = "577a6f7b9ee14b263b3e1cacd64d946c67c9130e"   # r1's counted declaration commit: its scripts are the audited ones (d6b95a05)
REL = "docs/qa/workflow-execution/kit/kit-275/"
CHANGED = ["sim-k275.R", "score-k275.R", "selftest/make-selftest.py", "selftest/run-selftest.R"]
UNCHANGED = ["fixture-k275.R", "e0-diff.py", "selftest/check-selftest.py"]
NEW = {"crosscheck-k275.R": "audit Q1 (b), the cross-implementation identity check (the pilot's control IS the counted control)"}

def git_show(path):
    return subprocess.run(["git", "show", "%s:%s%s" % (PIN, REL, path)], check=True, capture_output=True, text=True, cwd=HERE).stdout.split("\n")

def tags_of(lines):
    found = []
    for line in lines:
        i = line.find("K275-r2")
        if i >= 0: found.append(re.sub(r"^[:\s(]+", "", line[i + 7:]).strip()[:70])
    return sorted(set(found))

def label_only(old, new):
    return len(old) == len(new) and all(o.replace("V1g", "V1a") == n for o, n in zip(old, new)) and any("V1g" in o for o in old)

out, bad = [], 0
for path in CHANGED:
    old = git_show(path)
    new = open(os.path.join(HERE, path)).read().split("\n")
    matcher = difflib.SequenceMatcher(None, old, new, autojunk=False)
    hunks = [op for op in matcher.get_opcodes() if op[0] != "equal"]
    out.append("=" * 100)
    out.append("%s: %d hunks against %s (r1's audited version)" % (path, len(hunks), PIN[:8]))
    labels = 0
    for tag, i1, i2, j1, j2 in hunks:
        o, n = old[i1:i2], new[j1:j2]
        if label_only(o, n):
            labels += 1
            out.append("hunk old %d-%d new %d-%d  [8.1 label check (K2-07): V1g -> V1a]" % (i1 + 1, i2, j1 + 1, j2))
            continue
        tags = tags_of(n)
        # a pure deletion carries no new line: it is allowed only when the hunk is immediately followed or preceded by a tagged hunk (none are expected)
        if not tags:
            bad += 1
        out.append("hunk old %d-%d new %d-%d  [%s]" % (i1 + 1, i2, j1 + 1, j2, "; ".join(tags) if tags else "UNTAGGED"))
        out += ["  - " + line for line in o] + ["  + " + line for line in n]
    out.append("%s: %d label-only hunks, %d other hunks" % (path, labels, len(hunks) - labels))
for path in UNCHANGED:
    same = git_show(path) == open(os.path.join(HERE, path)).read().split("\n")
    out.append("%s: %s" % (path, "unchanged against r1" if same else "CHANGED (not allowed here)"))
    bad += 0 if same else 1
for path, why in NEW.items():
    n = len(open(os.path.join(HERE, path)).read().split("\n"))
    out.append("%s: new file, %d lines [%s]" % (path, n, why))
out.append("e0-r2: %d untagged or unexpected hunk(s)" % bad)
open(os.path.join(HERE, "e0-r2.txt"), "w").write("\n".join(out) + "\n")
print("\n".join(l for l in out if l.startswith(("e0-r2", "hunk", "=")) or ": " in l and not l.startswith("  ")))
sys.exit(1 if bad else 0)
