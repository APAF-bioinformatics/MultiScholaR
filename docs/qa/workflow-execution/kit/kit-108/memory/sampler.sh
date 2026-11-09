#!/usr/bin/env bash
# KIT-108 W2, instrument I2: a read-only per-process side sampler (D7).
#
# Runs outside the member's slice. Every <interval> s it lists the active user units matching <pattern> (or one named
# unit) through systemctl --user, reads each unit's ControlGroup, the cgroup's memory.stat (anon, file, shmem) and, for
# every PID in the cgroup's cgroup.procs (nested cgroups included), /proc/<pid>/comm, /proc/<pid>/status (VmRSS,
# RssAnon, RssFile, RssShmem, VmHWM) and /proc/<pid>/smaps_rollup (Rss, Pss_Anon, Anonymous). No process matching by
# name, no ptrace, no attach, no write anywhere but <out.tsv>.
#
# Each row carries both clocks (epoch from date, uptime from /proc/uptime, read together), so the driver's epoch
# marks align to the harness's uptime trace without a bracket.
#
# Usage:
#   sampler.sh <out.tsv> [--pattern 'multischolar*'] [--unit NAME] [--interval 1] [--wait 900] [--max 21600] [--grace S]
#   --pattern  glob for `systemctl --user list-units` (default 'multischolar*'); the sampler waits up to --wait s for a
#              match, then samples until no match is active (row event=unit_gone) and exits 0.
#   --unit     one exact unit instead: if it is not active (now or later), the sampler records unit_gone and exits 0.
#   --max      a ceiling in seconds; reaching it records event=max_reached and exits 0.
#   --grace    once sampling has begun, how long no match may last before unit_gone (default 30 s with --pattern, so a
#              unit set that is briefly empty during a launch is not taken for the end; 0 with --unit). A gap that
#              ends inside the grace is recorded as event=unit_gap.
# Exit 0 on unit_gone / max_reached / never_seen (all recorded); 2 on a usage error.
set -uo pipefail

out="${1:-}"; shift || true
[[ -n "$out" ]] || { echo "usage: sampler.sh <out.tsv> [--pattern P] [--unit U] [--interval S] [--wait S] [--max S]" >&2; exit 2; }
pattern='multischolar*'; unit=''; interval=1; wait_s=900; max_s=21600; grace=
while [[ $# -gt 0 ]]; do
    case "$1" in
        --pattern) pattern="$2"; shift 2 ;;
        --unit) unit="$2"; shift 2 ;;
        --interval) interval="$2"; shift 2 ;;
        --wait) wait_s="$2"; shift 2 ;;
        --max) max_s="$2"; shift 2 ;;
        --grace) grace="$2"; shift 2 ;;
        *) echo "unknown option $1" >&2; exit 2 ;;
    esac
done

header='epoch	uptime	event	unit	cgroup	pid	comm	vm_rss_kb	rss_anon_kb	rss_file_kb	rss_shmem_kb	vm_hwm_kb	smaps_rss_kb	smaps_pss_anon_kb	smaps_anonymous_kb	cg_anon	cg_file	cg_shmem	note'
printf '%s\n' "$header" > "$out"

clocks() { local e u; e=$(date +%s.%N); u=$(cut -d' ' -f1 /proc/uptime); printf '%.3f\t%s' "$e" "$u"; }
event() { printf '%s\t%s\t%s\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t%s\n' "$(clocks)" "$1" "${2:-}" "${3:-}" >> "$out"; }

activeUnits() {
    if [[ -n "$unit" ]]; then
        [[ "$(systemctl --user show -p ActiveState --value "$unit" 2>/dev/null)" == "active" ]] && printf '%s\n' "$unit"
    else
        systemctl --user list-units --state=active --no-legend --plain "$pattern" 2>/dev/null | awk 'NF{print $1}'
    fi
}

sampleUnit() {
    local u="$1" cg dir stat_line pid procs
    cg=$(systemctl --user show -p ControlGroup --value "$u" 2>/dev/null)
    [[ -n "$cg" ]] || return 0
    dir="/sys/fs/cgroup${cg}"
    [[ -d "$dir" ]] || return 0
    stat_line=$(awk '$1=="anon"{a=$2} $1=="file"{f=$2} $1=="shmem"{s=$2} END{printf "%s\t%s\t%s", a, f, s}' "$dir/memory.stat" 2>/dev/null)
    printf '%s\tcgroup\t%s\t%s\t\t\t\t\t\t\t\t\t\t\t%s\t\n' "$(clocks)" "$u" "$cg" "$stat_line" >> "$out"
    # A slice's processes are in its services' cgroups, sampled under their own units.
    [[ "$u" == *.slice ]] && return 0
    procs=$(find "$dir" -name cgroup.procs -exec cat {} + 2>/dev/null)
    for pid in $procs; do
        [[ -r "/proc/$pid/status" ]] || continue
        local comm status smaps
        comm=$(cat "/proc/$pid/comm" 2>/dev/null)
        status=$(awk '/^VmRSS:/{r=$2} /^RssAnon:/{a=$2} /^RssFile:/{f=$2} /^RssShmem:/{s=$2} /^VmHWM:/{h=$2} END{printf "%s\t%s\t%s\t%s\t%s", r, a, f, s, h}' "/proc/$pid/status" 2>/dev/null)
        smaps=$(awk '/^Rss:/{r=$2} /^Pss_Anon:/{p=$2} /^Anonymous:/{a=$2} END{printf "%s\t%s\t%s", r, p, a}' "/proc/$pid/smaps_rollup" 2>/dev/null || printf '\t\t')
        [[ -n "$status" ]] || continue
        printf '%s\tprocess\t%s\t%s\t%s\t%s\t%s\t%s\t\t\t\t\n' "$(clocks)" "$u" "$cg" "$pid" "$comm" "$status" "$smaps" >> "$out"
    done
}

[[ -n "$grace" ]] || { [[ -n "$unit" ]] && grace=0 || grace=30; }
event start "" "pattern=${pattern} unit=${unit:-none} interval=${interval} wait=${wait_s} max=${max_s} grace=${grace} pid=$$"
started=$SECONDS; seen=0; empty_since=
while :; do
    units=$(activeUnits)
    if [[ -z "$units" ]]; then
        if [[ "$seen" == 1 || -n "$unit" ]]; then
            [[ -n "$empty_since" ]] || empty_since=$SECONDS
            if (( SECONDS - empty_since >= grace )); then
                event unit_gone "${unit:-$pattern}" "$([[ $seen == 1 ]] && echo "after sampling" || echo "not active at start")"
                exit 0
            fi
            sleep "$interval"; continue
        fi
        if (( SECONDS - started >= wait_s )); then event never_seen "$pattern" "no active unit within ${wait_s}s"; exit 0; fi
        sleep "$interval"; continue
    fi
    if [[ -n "$empty_since" ]]; then event unit_gap "" "no match for $((SECONDS - empty_since))s"; empty_since=''; fi
    if [[ "$seen" == 0 ]]; then seen=1; event unit_seen "" "$(echo $units | tr ' ' ',')"; fi
    for u in $units; do sampleUnit "$u"; done
    if (( SECONDS - started >= max_s )); then event max_reached "" "${max_s}s"; exit 0; fi
    sleep "$interval"
done
