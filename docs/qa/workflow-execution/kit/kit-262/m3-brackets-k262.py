"""KIT-262 M3 local check (standard library only): bracket balance of the R scripts, ignoring comments and string literals.
Run from the repository root: python3 -I docs/qa/workflow-execution/kit/kit-262/m3-brackets-k262.py <file.R> [...]
Prints one line per file and exits 1 on an imbalance. It does not parse R: it catches an unclosed or extra ( [ { only."""
import sys

PAIRS = {")": "(", "]": "[", "}": "{"}


def check(path):
    text = open(path, encoding="utf-8").read()
    stack, i, line = [], 0, 1
    n = len(text)
    while i < n:
        c = text[i]
        if c == "\n":
            line += 1
        elif c == "#":
            while i < n and text[i] != "\n":
                i += 1
            continue
        elif c in "\"'`":
            q = c
            i += 1
            while i < n and text[i] != q:
                if text[i] == "\\":
                    i += 1
                if i < n and text[i] == "\n":
                    line += 1
                i += 1
        elif c in "([{":
            stack.append((c, line))
        elif c in ")]}":
            if not stack or stack[-1][0] != PAIRS[c]:
                return "FAIL %s: unmatched %r at line %d" % (path, c, line)
            stack.pop()
        i += 1
    if stack:
        return "FAIL %s: unclosed %r from line %d" % (path, stack[-1][0], stack[-1][1])
    return "PASS %s: %d lines, brackets balanced" % (path, line)


bad = 0
for p in sys.argv[1:]:
    r = check(p)
    print(r)
    bad += r.startswith("FAIL")
sys.exit(1 if bad else 0)
