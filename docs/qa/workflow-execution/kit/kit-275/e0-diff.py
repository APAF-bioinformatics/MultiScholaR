# KIT-275 E0 (the generator copy): diffs sim-k275.R against `git show 0866982f:docs/qa/workflow-execution/kit/kit-188/sim-k188.R` and lists
# exactly the marked changes. The blocks of KIT-188 that are not carried (the proteomics journeys, KIT-188's estimands, checks, units, plants and
# job mechanics) are removed from the original first, by the line ranges OMITTED below (numbers at 0866982f); what remains of the original must
# equal sim-k275.R's text before its added section except at the marked hunks. Standard library only: python3 -I e0-diff.py [--patch <file>]
# Exits 1 when a hunk of the retained text carries no K275- mark (outside the header comment) or the added section's marker line is missing.
import difflib, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
PIN = "0866982f"
PATH = "docs/qa/workflow-execution/kit/kit-188/sim-k188.R"
OMITTED = [(86, 93, "BLIND, the estimand columns blinded in KIT-188's lifecycle"), (809, 969, "Part P: the protein journeys (pConfig, REQUESTS_P, protMatrixOf, journeyP)"),
           (978, 1079, "KIT-188's estimands: scoreFamily, ASSAY_COLUMNS, ROW_COLUMNS, newRow, armRows"),
           (1105, 1161, "checkE9L (E9 of the LP arm; its direct-call idea is reused in checkE12)"),
           (1180, 2283, "E9P to E4, the units (ARMS_OF, conditionsOf, runCondition*, plants, runUnit), the child and the parent: replaced by the K275 section")]
ADDED_MARK = "# K275: everything below is added (not copied)"
# a pure deletion inside a marked block of the original carries that block's mark (the old line numbers of 0866982f)
MARK_RANGES = [("K275-E3", 428, 453)]

def main():
    original = subprocess.run(["git", "show", "%s:%s" % (PIN, PATH)], check=True, capture_output=True, text=True, cwd=HERE).stdout.split("\n")
    dropped = set()
    for first, last, _ in OMITTED:
        dropped.update(range(first, last + 1))
    kept = [line for number, line in enumerate(original, 1) if number not in dropped]
    new = open(os.path.join(HERE, "sim-k275.R")).read().split("\n")
    while kept and kept[-1] == "":
        kept.pop()
    marker = [i for i, line in enumerate(new) if line.startswith(ADDED_MARK)]
    if len(marker) != 1:
        print("E0 FAILS: the added section's marker line is missing or repeated")
        return 1
    # the added section starts two lines above its marker (a blank line and a rule); the text before it is the retained copy
    cut = marker[0] - 2
    while cut > 0 and new[cut - 1].strip() == "":
        cut -= 1
    retained_new = new[:cut]
    retained_old = kept[:len(kept)]
    # the old text of 1099-1104 and 1162-1179 stays (positional, omicObjectFromMatrices); the old tail beyond 2283 is empty
    matcher = difflib.SequenceMatcher(None, retained_old, retained_new, autojunk=False)
    bad = 0
    out = []
    old_numbers = [n for n in range(1, len(original) + 1) if n not in dropped]
    for tag, i1, i2, j1, j2 in matcher.get_opcodes():
        if tag == "equal":
            continue
        first_old = old_numbers[i1] if i1 < len(old_numbers) else None
        last_old = old_numbers[i2 - 1] if i2 > i1 else None
        text = "\n".join(retained_new[j1:j2])
        marks = sorted({w.strip("();,:") for w in text.replace("#", " ").split() if w.startswith("K275-") or w.startswith("(K275-")})
        header = j2 <= 50 and all((not l.strip()) or l.lstrip().startswith("#") for l in retained_new[j1:j2])
        if not marks and tag == "delete" and first_old is not None:
            for name, lo, hi in MARK_RANGES:
                if lo <= first_old and last_old <= hi:
                    marks = [name]
        if not marks and not header:
            bad += 1
        out.append("%-8s old %s-%s  new %d-%d  %s" % (tag, first_old, last_old, j1 + 1, j2, ",".join(marks) if marks else ("header comment" if header else "UNMARKED")))
    print("\n".join(out))
    print("E0: %d hunks, %d unmarked" % (len(out), bad))
    if "--patch" in sys.argv:
        path = sys.argv[sys.argv.index("--patch") + 1]
        with open(path, "w") as f:
            f.writelines(difflib.unified_diff([l + "\n" for l in retained_old], [l + "\n" for l in retained_new], "sim-k188.R@%s (omitted ranges removed)" % PIN, "sim-k275.R (before the added section)", n=2))
    return 1 if bad else 0

sys.exit(main())
