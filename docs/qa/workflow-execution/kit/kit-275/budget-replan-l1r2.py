# KIT-275 r2.1 G2: every counted job re-planned from L1-r2's OWN completed-journey means per arm class and design (run 45929680145, `l1r2/out/scripts/4/out/timing-l1r2.tsv`),
# nothing scaled, the pilot's and the diagnostic's timings NOT used (P-3). Rule (ticket 10, r2.1): compute = data sets x seconds per data set / 4 / 60 (4 workers);
# timeout_minutes = ceil(1.5 x compute); ceiling = timeout + 10 job-minutes of set-up; expected = compute + 13; core-hours = job-minutes x 4 / 60; a job whose plan exceeds 150
# minutes is split at the same R, each part's timeout 1.5 x its plan. Standard library only: python3 -I budget-replan-l1r2.py
import csv, math, os
HERE = os.path.dirname(os.path.abspath(__file__))
T = {}
for r in csv.DictReader(open(os.path.join(HERE, "l1r2/out/scripts/4/out/timing-l1r2.tsv")), delimiter="\t"):
    if r["mean_seconds"] != "NA": T.setdefault(r["design"], {})[r["arm_class"]] = float(r["mean_seconds"])
# journeys per data set by arm class (ticket 4; LX has no journey): tier D now has V1w's NO journey (NO = 3)
mix = {"D": dict(NO=3, LM=3, GK=3, LA=3, GA=3), "K": dict(NO=1, LM=8, LG=3, GK=5), "F": dict(NO=1, LM=1, GK=1, LA=1, GA=1), "S": dict(GSC=3, GK=1),
       "V": dict(NO=2, LM=3, GK=3, LA=2, GA=2), "M": dict(NO=1, LM=3, GK=3), "T": dict(NO=1, LM=1, LG=1, GK=1, LA=1, GA=1), "E": dict(NO=2, PEL2=2, PEG2=2),
       "ID": dict(NO=1, LM=1, GK=1, LA=1, GA=1)}
def sec(tiers, d): return sum(n * T[d][c] for t in (["F", "S"] if tiers == "FS" else [tiers]) for c, n in mix[t].items())
rows, ceiling, expected = [], 0, 0.0
def job(name, tiers, d, R, n_jobs):
    global ceiling, expected
    s = sec(tiers, d)
    compute = R * s / 4 / 60
    parts = 1
    while compute / parts > 150: parts += 1   # the split rule (none is expected here)
    part_compute = compute / parts
    timeout = math.ceil(1.5 * part_compute)
    c, e = timeout + 10, part_compute + 13
    ceiling += c * n_jobs * parts; expected += e * n_jobs * parts
    rows.append((name, d, "%d x %d" % (n_jobs * parts, R // parts), round(s, 1), round(part_compute, 1), timeout, c, round(e, 1)))
for d in ("D3", "D4", "D6"):
    job("D", "D", d, 50, 4); job("K", "K", d, 50, 2); job("F+S", "FS", d, 50, 2); job("V", "V", d, 50, 2); job("T", "T", d, 100, 1); job("E", "E", d, 100, 1)
job("M", "M", "D6", 100, 1); job("ID", "ID", "D4", 20, 1)
print("tier design jobs x data sets | seconds per data set | plan (compute) min | timeout | ceiling | expected   (per job)")
for r in rows: print(r)
print("counted jobs: %d, ceiling %d job-min = %.2f core-h; expected %.1f job-min = %.2f core-h; largest ceiling %d" % (sum(int(r[2].split(" x ")[0]) for r in rows), ceiling, ceiling * 4 / 60, expected, expected * 4 / 60, max(r[6] for r in rows)))

# the job list (declaration-r2.md section 3): tier D 12 jobs (4 per design, r = 2001 to 2200), K 6 and F+S 6 and V 6 (2 per design, r = 2001 to 2100), M 1, ID 1, T 3, E 3
jobs = []
for d in ("D3", "D4", "D6"):
    plan = {}
    for name, tiers, R, n in (("D", "D", 50, 4), ("K", "K", 50, 2), ("FS", "F,S", 50, 2), ("V", "V", 50, 2), ("T", "T", 100, 1), ("E", "E", 100, 1)):
        key = "FS" if name == "FS" else tiers
        s_ = sec(key, d)
        timeout = math.ceil(1.5 * R * s_ / 4 / 60)
        for i in range(n):
            first = 2001 + i * R
            jobs.append((name + "-" + d + "-%d-%d" % (first, first + R - 1), "%s %s %d %d" % (tiers, d, first, first + R - 1), timeout, timeout + 10))
for name, tiers, d, R in (("M", "M", "D6", 100), ("ID", "ID", "D4", 20)):
    timeout = math.ceil(1.5 * R * sec(tiers, d) / 4 / 60)
    jobs.append((name + "-" + d + "-2001-%d" % (2000 + R), "%s %s 2001 %d" % (tiers, d, 2000 + R), timeout, timeout + 10))
print("JOBS", len(jobs), sum(j[3] for j in jobs))
for j in jobs: print("JOB|%s|%s|%d|%d" % j)
