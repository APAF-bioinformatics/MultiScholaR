# Job D (determinism, E5): the rows of data sets 1 to 5 of L12 and M12 deciding, run again in job D, equal those of the S jobs with the columns `seconds` and `warnings`
# removed. Standard library only; it prints SHA-256 values and a verdict, never a row or a value. Usage: python3 -I e5d-compare.py <first.tsv.gz> <second.tsv.gz> [data-set range lo hi]
import csv, gzip, hashlib, sys
def digest(path, lo, hi):
    rows = []
    with gzip.open(path, "rt", newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if lo <= int(r["rep"]) <= hi and r["arm"] != "unit":
                for c in ("seconds", "warnings"): r.pop(c, None)
                rows.append("\t".join(r[k] for k in sorted(r)))
    rows.sort()
    return len(rows), hashlib.sha256("\n".join(rows).encode()).hexdigest()
lo, hi = (int(sys.argv[3]), int(sys.argv[4])) if len(sys.argv) > 4 else (1, 5)
a, b = digest(sys.argv[1], lo, hi), digest(sys.argv[2], lo, hi)
print("rows %d / %d, sha256 %s / %s: %s" % (a[0], b[0], a[1][:16], b[1][:16], "identical" if a == b else "DIFFERENT"))
sys.exit(0 if a == b else 1)
