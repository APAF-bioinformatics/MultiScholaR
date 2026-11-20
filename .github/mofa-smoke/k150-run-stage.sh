#!/usr/bin/env bash
# KIT-150 (branch run-mofa-measure-k150): run one stage (k150-stage.R) in a fresh R process, alone in its own memory
# cgroup (cgroup v2, under $K150_CG), with a memory limit (memory.max, no swap) and a time limit (timeout: TERM, then
# KILL 30 s later), and record how it ended. Declaration: docs/qa/workflow-execution/kit/kit-150/declaration.md.
#
# Usage: k150-run-stage.sh <stage arguments (JSON)> <memory limit, bytes> <time limit, seconds>
# Environment: BUNDLE (the launcher), LANE (this lane's checkout), K150_CG, K150_LOGS.
# Writes <out>.driver.json: the exit status, the time limit hit or not, the cgroup's oom_kill count and memory.peak,
# the wall time of the process, and the outcome:
#   completed | error (the stage's own error) | killed: memory limit | stopped: time limit |
#   stopped: exit <status> (after the baseline) | failed before the baseline (exit <status>)
# Exit status 0 only when the stage completed.
set -u

args=$1
memory=$2
limit=$3
id=$(jq -r .id "$args")
out=$(jq -r .out "$args")
cg="$K150_CG/$id"
log="$K150_LOGS/$id.log"

sudo mkdir "$cg"
echo "$memory" | sudo tee "$cg/memory.max" > /dev/null
if [ -e "$cg/memory.swap.max" ]; then echo 0 | sudo tee "$cg/memory.swap.max" > /dev/null; fi

started=$(date +%s.%N)
# The inner shell moves itself into the cgroup, then becomes the launcher (which execs R: one process throughout).
timeout --kill-after=30 "$limit" bash -c 'echo $$ | sudo tee "$1/cgroup.procs" > /dev/null && exec "$2" --vanilla "$3" "$4"' \
    k150 "$cg" "$BUNDLE/bin/Rscript" "$LANE/.github/mofa-smoke/k150-stage.R" "$args" > "$log" 2>&1
status=$?
ended=$(date +%s.%N)

peak=$(cat "$cg/memory.peak" 2> /dev/null || true)
oom=$(awk '$1 == "oom_kill" { print $2 }' "$cg/memory.events" 2> /dev/null || true)
sudo rmdir "$cg" 2> /dev/null || echo "k150: cgroup $cg not removed" >&2

if [ -f "$out.json" ] && [ "$(jq -r .status "$out.json")" = completed ]; then
    outcome="completed"
elif [ -f "$out.json" ]; then
    outcome="error"
elif [ -n "$oom" ] && [ "$oom" != 0 ]; then
    outcome="killed: memory limit"
elif [ "$status" = 124 ]; then
    outcome="stopped: time limit"
elif [ -f "$out.baseline.json" ]; then
    outcome="stopped: exit $status"
else
    outcome="failed before the baseline (exit $status)"
fi

jq -n --arg id "$id" --arg outcome "$outcome" --argjson status "$status" --arg oom "${oom:-}" --arg peak "${peak:-}" \
    --arg started "$started" --arg ended "$ended" --argjson memory "$memory" --argjson limit "$limit" '{
    id: $id, outcome: $outcome, exit_status: $status,
    oom_kill: (if $oom == "" then null else ($oom | tonumber) end),
    cgroup_peak_mib: (if $peak == "" then null else (($peak | tonumber) / 1048576 * 10 | round / 10) end),
    process_wall_s: ((($ended | tonumber) - ($started | tonumber)) * 100 | round / 100),
    memory_limit_mib: ($memory / 1048576 | round), time_limit_s: $limit }' > "$out.driver.json"

echo "k150 stage $id: $outcome (exit $status; cgroup peak ${peak:-?} B; oom_kill ${oom:-?})"
tail -n 3 "$log" | sed 's/^/    /'
[ "$outcome" = completed ]
