#!/usr/bin/env bash
# KIT-209 (branch run-kit209-memory): run one journey (k209-journey.R) in a fresh R process, alone in its own memory
# cgroup (cgroup v2, under $K209_CG), with a memory limit (memory.max, no swap) and a time limit (timeout: TERM, then
# KILL 30 s later), after dropping the page cache; record how it ended. KIT-150's k150-run-stage.sh, adapted.
#
# Usage: k209-run-journey.sh <journey arguments (JSON)> <memory limit, bytes> <time limit, seconds>
# Environment: LANE (this lane's checkout), K209_CG, K209_LOGS, RSTUDIO_PANDOC.
# Writes <out>.driver.json: the exit status, the cgroup's oom_kill count, memory.peak and memory.events, the process's
# wall time, the cache drop, and the outcome:
#   completed | error (the journey's own error record) | killed: memory limit | stopped: time limit |
#   stopped: exit <status> (after the baseline) | failed before the baseline (exit <status>)
# Exit status 0 only when the journey completed.
set -u

args=$1
memory=$2
limit=$3
id=$(jq -r .id "$args")
out=$(jq -r .out "$args")
cg="$K209_CG/$id"
log="$K209_LOGS/$id.log"

sync
if echo 3 | sudo tee /proc/sys/vm/drop_caches > /dev/null; then dropped=true; else dropped=false; fi

sudo mkdir "$cg"
echo "$memory" | sudo tee "$cg/memory.max" > /dev/null
if [ -e "$cg/memory.swap.max" ]; then echo 0 | sudo tee "$cg/memory.swap.max" > /dev/null; fi

started=$(date +%s.%N)
# The inner shell moves itself into the cgroup, then execs R (one process throughout; pandoc runs as its child).
timeout --kill-after=30 "$limit" bash -c 'echo $$ | sudo tee "$1/cgroup.procs" > /dev/null && exec "$2" --vanilla "$3" "$4"' \
    k209 "$cg" /opt/R/4.6.1/bin/Rscript "$LANE/.github/k209/k209-journey.R" "$args" > "$log" 2>&1
status=$?
ended=$(date +%s.%N)

peak=$(cat "$cg/memory.peak" 2> /dev/null || true)
oom=$(awk '$1 == "oom_kill" { print $2 }' "$cg/memory.events" 2> /dev/null || true)
events=$(tr '\n' ';' < "$cg/memory.events" 2> /dev/null || true)
sudo rmdir "$cg" 2> /dev/null || echo "k209: cgroup $cg not removed" >&2

if [ -f "$out.json" ] && [ "$(jq -r .status "$out.json")" = completed ]; then
    outcome="completed"
elif [ -f "$out.json" ]; then
    outcome="error"
elif [ -n "$oom" ] && [ "$oom" != 0 ]; then
    outcome="killed: memory limit"
elif [ "$status" = 124 ] || [ "$status" = 137 ]; then
    outcome="stopped: time limit"
elif [ -f "$out.baseline.json" ]; then
    outcome="stopped: exit $status"
else
    outcome="failed before the baseline (exit $status)"
fi

jq -n --arg id "$id" --arg outcome "$outcome" --argjson status "$status" --arg oom "${oom:-}" --arg peak "${peak:-}" \
    --arg events "${events:-}" --arg started "$started" --arg ended "$ended" --argjson memory "$memory" \
    --argjson limit "$limit" --argjson dropped "$dropped" '{
    id: $id, outcome: $outcome, exit_status: $status,
    oom_kill: (if $oom == "" then null else ($oom | tonumber) end),
    cgroup_peak_mib: (if $peak == "" then null else (($peak | tonumber) / 1048576 * 10 | round / 10) end),
    memory_events: $events, cache_dropped: $dropped,
    process_wall_s: ((($ended | tonumber) - ($started | tonumber)) * 100 | round / 100),
    memory_limit_bytes: $memory, time_limit_s: $limit }' > "$out.driver.json"

echo "k209 journey $id: $outcome (exit $status; cgroup peak ${peak:-?} B; oom_kill ${oom:-?})"
tail -n 4 "$log" | sed 's/^/    /'
[ "$outcome" = completed ]
