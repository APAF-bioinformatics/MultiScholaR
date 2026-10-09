#!/usr/bin/env bash
# Guarded launcher and 50 ms memory sampler for KIT stages (KIT-001).
# Usage:
#   run_workflow_kit.sh validate <stage> <library>
#   run_workflow_kit.sh run <stage> <library>
#   run_workflow_kit.sh --sample <cgroup> <trace> <deadline_uptime> <slice> <stop_file> <guard_file>
# The proof root comes from docs/qa/workflow-execution/kit/bindings.json (harness.root).
set -euo pipefail
kit_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd -- "$kit_dir/../.." && pwd)"

if [[ "${1:-}" == "--sample" ]]; then
    [[ "$#" == 7 ]] || exit 2
    group="$2"; trace="$3"; deadline="$4"; slice="$5"; stop_file="$6"; guard_file="$7"
    printf 'uptime\tcurrent\tpeak\tanon\tfile\tkernel\tswap_peak\n' > "$trace"
    while [[ ! -e "$stop_file" ]]; do
        if [[ ! -d "$group" ]]; then
            printf '%s\tcgroup_gone\n' "$(cut -d' ' -f1 /proc/uptime)" >> "${trace%.tsv}-gaps.tsv"
            sleep 0.05; continue
        fi
        uptime="$(cut -d' ' -f1 /proc/uptime)"
        current="$(cat "$group/memory.current" 2>/dev/null || echo NA)"
        peak="$(cat "$group/memory.peak" 2>/dev/null || echo NA)"
        swap="$(cat "$group/memory.swap.peak" 2>/dev/null || echo NA)"
        read -r anon file kernel < <(awk '$1=="anon"{a=$2} $1=="file"{f=$2} $1=="kernel"{k=$2} END{print a, f, k}' "$group/memory.stat" 2>/dev/null || echo "NA NA NA")
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$uptime" "$current" "$peak" "$anon" "$file" "$kernel" "$swap" >> "$trace"
        available="$(awk '/^MemAvailable:/ {printf "%.0f", $2 * 1024}' /proc/meminfo)"
        reason=""
        if [[ "$peak" =~ ^[0-9]+$ ]] && (( peak >= 16642998272 )); then reason="memory_guard"; fi
        if [[ "$available" =~ ^[0-9]+$ ]] && (( available < 2147483648 )); then reason="host_guard"; fi
        if awk -v u="$uptime" -v d="$deadline" 'BEGIN{exit !(u >= d)}'; then reason="time_guard"; fi
        if [[ -n "$reason" ]]; then
            printf '%s\n' "$reason" > "$guard_file"
            systemctl --user stop "$slice" || true
            exit 0
        fi
        sleep 0.05
    done
    exit 0
fi

mode="${1:-}"; stage="${2:-}"; library="${3:-}"
[[ ( "$mode" == validate || "$mode" == run ) && "$stage" =~ ^kit-(lifecycle-(complete|fail)|journey)$ && -d "$library" ]] || {
    echo "Usage: run_workflow_kit.sh validate|run <stage> <library>" >&2; exit 2; }
root="$(jq -r '.harness.root' "$repo/docs/qa/workflow-execution/kit/bindings.json")"
# The browser driver runs in the coordinator (outside the slice); chromote needs the executable (as the pinned runner does).
export CHROMOTE_CHROME="${CHROMOTE_CHROME:-$(command -v chromium || command -v google-chrome || command -v brave || true)}"
[[ -x "$CHROMOTE_CHROME" ]] || { echo "refused: no Chromium-based browser executable found" >&2; exit 2; }
[[ "$root" == /* ]] || { echo "bindings carry no proof root" >&2; exit 2; }
mkdir -p -- "$root"
cd -- "$repo"
if [[ "$mode" == run ]]; then
    # Preconditions (plan): clean tree, bindings committed, 12 GiB free, no other managed unit, no live stage.
    [[ -z "$(git status --porcelain)" ]] || { echo "refused: working tree is not clean" >&2; exit 76; }
    available="$(awk '/^MemAvailable:/ {printf "%.0f", $2 * 1024}' /proc/meminfo)"
    (( available >= 12884901888 )) || { echo "refused: host headroom is below 12 GiB" >&2; exit 76; }
    units="$(systemctl --user list-units --state=active --no-legend --plain 'multischolar*' | sed '/^\s*$/d' || true)"
    [[ -z "$units" ]] || { echo "refused: another managed unit is active: $units" >&2; exit 76; }
    [[ ! -e "$root/$stage" ]] || { echo "refused: stage directory exists: $root/$stage" >&2; exit 76; }
    { date; echo "MemAvailable $(awk -v a="$available" 'BEGIN{printf "%.2f GiB", a/1073741824}')"; git rev-parse HEAD; } \
        > "$root/../$stage-start.txt"
    nohup Rscript --vanilla "$kit_dir/workflow_kit_coordinator.R" "$stage" "$library" "$root" \
        > "$root/../$stage.log" 2>&1 &
    echo "launched $stage (pid $!); log $root/../$stage.log"
    exit 0
fi
exec Rscript --vanilla "$kit_dir/workflow_kit_coordinator.R" "$stage" "$library" "$root" --validate-only
