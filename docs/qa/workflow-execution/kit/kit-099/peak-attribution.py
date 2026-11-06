# KIT-099: attribute a journey member's anonymous memory to the worker operation running at each moment. The slice trace
# (candidate-1/memory-trace.tsv, kernel uptime) is put on the wall clock with the host's boot time (/proc/stat btime); each
# worker request runs from its request.rds time to its directory's last file time. Per operation: the largest whole-slice
# anonymous and the largest application (controller) anonymous sample inside its window. Usage: python3 peak-attribution.py <proof root> <stage>
import csv, os, sys, glob
root, stage = sys.argv[1], sys.argv[2]
member = os.path.join(root, "kit/runtime", stage, "candidate-1")
btime = int([l for l in open("/proc/stat") if l.startswith("btime")][0].split()[1])
rows = []
with open(os.path.join(member, "memory-trace.tsv")) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        try: rows.append((btime + float(r["uptime"]), float(r["anon"]), float(r["controller_anon"]) if r["controller_anon"] not in ("NA", "") else float("nan")))
        except (ValueError, KeyError): pass
ops = []
requests = [os.path.join(dp, "request.rds") for dp, _, fns in os.walk(os.path.join(member, "project")) if "request.rds" in fns and "/.workflow-export/wfx-" in dp]
for req in requests:
    d = os.path.dirname(req)
    files = [os.path.join(dp, fn) for dp, _, fns in os.walk(d) for fn in fns]
    start, end = os.path.getmtime(req), max(os.path.getmtime(p) for p in files)
    op = os.popen(f"Rscript -e 'cat(as.character(readRDS(\"{req}\")$operation))'").read().strip() or os.path.basename(d)
    ops.append((start, end, op))
ops.sort()
peak = max(rows, key=lambda r: r[1])
print(f"# {stage}: whole-slice anonymous peak {peak[1]/2**30:.2f} GiB at {peak[0]-rows[0][0]:.0f} s into the trace; worker running then:",
      ", ".join(o for s, e, o in ops if s <= peak[0] <= e) or "none (application only)")
print(f"{'operation':24} {'seconds':>8} {'slice anon GiB':>15} {'app anon GiB':>13}")
for s, e, o in ops:
    inside = [r for r in rows if s <= r[0] <= e]
    if inside: print(f"{o:24} {e-s:8.1f} {max(r[1] for r in inside)/2**30:15.2f} {max(r[2] for r in inside)/2**30:13.2f}")
