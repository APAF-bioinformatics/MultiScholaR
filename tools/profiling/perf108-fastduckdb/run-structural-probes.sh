#!/usr/bin/env bash
set -euo pipefail
# Fixed diagnostic schedule; abort on any failure, never retry or overwrite.
repo=$(pwd -P)
proof=/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/perf-108/structural
evidence="$repo/docs/qa/workflow-execution/perf-108/structural"
old=/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/perf-108/cache-qualification/browser-smoke.evZ5ZU/library
[[ ! -e "$proof/structural-fixture.rds" ]]
for slot in ZF PA0 PA1 PV0 PV1 PO0 PO1; do
    [[ ! -e "$evidence/$slot.log" ]]
    case "$slot" in
        ZF|PO1) library="$proof/browser-smoke.IsZngp/library" ;;
        PA0|PA1) library="$proof/browser-smoke.DpyLk9/library" ;;
        PV0|PV1) library="$proof/browser-smoke.U3PFHG/library" ;;
        PO0) library="$old" ;;
    esac
    printf 'Starting %s\n' "$slot"
    bash tools/profiling/run_workflow_checkpoint_proof.sh --service reader-check "$library" "$proof" "$repo" "perf108:$slot" > "$evidence/$slot.log" 2>&1
done
printf 'Fixed diagnostic schedule completed.\n'
