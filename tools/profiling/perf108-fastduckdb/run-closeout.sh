#!/usr/bin/env bash
set -euo pipefail
repo=$(pwd -P)
proof=/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/perf-108/closeout
library=/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/perf-108/ordered-groups/browser-smoke.jEmE2F/library
evidence="$repo/docs/qa/workflow-execution/perf-108/closeout"
[[ -f "$evidence/bindings.rds" && ! -e "$evidence/ZQ.log" && ! -e "$proof/runtime/closure-filters" ]]
mkdir -p "$proof/runtime" "$proof/bulk-rollup-library"
for slot in ZQ PC; do
    [[ ! -e "$evidence/$slot.log" ]]
    printf 'Starting %s\n' "$slot"
    bash tools/profiling/run_workflow_checkpoint_proof.sh --service reader-check "$library" "$proof" "$repo" "perf108:$slot" > "$evidence/$slot.log" 2>&1
done
printf 'Starting eight fixed primary/wide filter members\n'
Rscript --vanilla tools/profiling/workflow_perf106_repair.R closure-filters "$library" "$proof/runtime" > "$evidence/filters.log" 2>&1
printf 'Closeout schedule complete; inspect decision before release.\n'
