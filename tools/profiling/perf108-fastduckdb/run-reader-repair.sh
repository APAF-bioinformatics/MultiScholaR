#!/usr/bin/env bash
set -euo pipefail
repo=$(pwd -P)
proof=/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/perf-108/ordered-groups
base=/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/perf-108/structural/browser-smoke.IsZngp/library
evidence="$repo/docs/qa/workflow-execution/perf-108/ordered-groups"
if [[ "${1:-repair}" == compare-and-reducer ]]; then
    reader_library="$proof/browser-smoke.KHVwi8/library"
    for slot in RC0 RC1; do
        [[ ! -e "$evidence/$slot.log" ]]
        printf 'Starting %s\n' "$slot"
        bash tools/profiling/run_workflow_checkpoint_proof.sh --service reader-check "$reader_library" "$proof" "$repo" "perf108:$slot" > "$evidence/$slot.log" 2>&1
    done
    bash tools/profiling/perf108-fastduckdb/run-ordered-build.sh 3
    exit 0
fi
[[ "${1:-repair}" == repair ]]
[[ ! -e "$evidence/C-build-1.log" && ! -e "$evidence/C-check-1.log" ]]
mkdir -p "$proof/bulk-rollup-library"
bash tools/profiling/run_workflow_checkpoint_proof.sh --service reader-check "$base" "$proof" "$repo" perf108:M1 > "$evidence/C-build-1.log" 2>&1
unit=$(head -n 1 "$evidence/C-build-1.log" | python3 -c 'import json,sys; print(json.load(sys.stdin)["unit"])')
member=${unit#multischolar-dia003-}
member=${member%.service}
[[ "$member" =~ ^browser-smoke\.[A-Za-z0-9]{6}$ ]]
bash tools/profiling/run_workflow_checkpoint_proof.sh --service reader-check "$proof/$member/library" "$proof" "$repo" perf108:ZC > "$evidence/C-check-1.log" 2>&1
printf 'Reader repair build/check complete: %s\n' "$proof/$member/library"
