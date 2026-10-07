#!/usr/bin/env bash
set -euo pipefail
repo=$(pwd -P)
proof=/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/perf-108/ordered-groups
base=/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/perf-108/structural/browser-smoke.IsZngp/library
evidence="$repo/docs/qa/workflow-execution/perf-108/ordered-groups"
if [[ "${1:-1}" == compare ]]; then
    library="$proof/browser-smoke.jEmE2F/library"
    for slot in RS0 RS1; do
        [[ ! -e "$evidence/$slot.log" ]]
        printf 'Starting %s\n' "$slot"
        bash tools/profiling/run_workflow_checkpoint_proof.sh --service reader-check "$library" "$proof" "$repo" "perf108:$slot" > "$evidence/$slot.log" 2>&1
    done
    exit 0
fi
attempt=${1:-1}
[[ "$attempt" == 1 || "$attempt" == 2 || "$attempt" == 3 ]]
[[ ! -e "$evidence/S-build-$attempt.log" && ! -e "$evidence/S-check-$attempt.log" ]]
bash tools/profiling/run_workflow_checkpoint_proof.sh --service reader-check "$base" "$proof" "$repo" perf108:M1 > "$evidence/S-build-$attempt.log" 2>&1
unit=$(head -n 1 "$evidence/S-build-$attempt.log" | python3 -c 'import json,sys; print(json.load(sys.stdin)["unit"])')
member=${unit#multischolar-dia003-}
member=${member%.service}
[[ "$member" =~ ^browser-smoke\.[A-Za-z0-9]{6}$ ]]
bash tools/profiling/run_workflow_checkpoint_proof.sh --service reader-check "$proof/$member/library" "$proof" "$repo" perf108:ZS > "$evidence/S-check-$attempt.log" 2>&1
printf 'Ordered reducer build/check complete: %s\n' "$proof/$member/library"
