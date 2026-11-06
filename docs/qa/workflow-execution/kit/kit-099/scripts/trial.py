#!/usr/bin/env python3
# KIT-099 trial.py <output file> <ceiling s> <step limit s> <command...>: run one in-process trial point (KIT-084's timed.py
# plus a step limit). The recorder prints one line per finished step ("  <step> <state> rows <n> <s> s"); a step that runs
# longer than the step limit since the previous line (or the start) ends the point, as the worker's 1,800 s step deadline
# would. Writes the wall time, the largest single-process RSS among the waited children and the stop reason.
# Exit: the command's code; 124 at the ceiling; 125 at the step limit. The process group is stopped by its own id.
import os, re, resource, signal, subprocess, sys, time
out, ceil, limit, cmd = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4:]
step = re.compile(r"^  \S+\s+\S+\s+rows ")
start = last = time.time(); reason = "exit"
with open(out, "w") as f:
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, start_new_session=True, text=True, bufsize=1)
    os.set_blocking(proc.stdout.fileno(), False)
    while True:
        line = proc.stdout.readline()
        if line:
            f.write(line); f.flush()
            if step.match(line): last = time.time()
            continue
        if proc.poll() is not None:
            rest = proc.stdout.read()
            if rest: f.write(rest)
            break
        now = time.time()
        if now - start > ceil: reason = "ceiling"
        elif now - last > limit: reason = "step_limit"
        if reason != "exit":
            os.killpg(proc.pid, signal.SIGTERM); time.sleep(10)
            if proc.poll() is None: os.killpg(proc.pid, signal.SIGKILL)
            proc.wait(); break
        time.sleep(1)
    rc = {"ceiling": 124, "step_limit": 125}.get(reason, proc.returncode)
    rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    f.write(f"\nRSS_KB {rss} ELAPSED {time.time() - start:.0f} STOP {reason} RUNNING_STEP_SECONDS {time.time() - last:.0f}\n")
sys.exit(rc)
