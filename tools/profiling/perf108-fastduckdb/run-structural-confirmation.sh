#!/usr/bin/env bash
set -euo pipefail
repo=$(pwd -P)
root=/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/perf-108/structural/runtime
library=/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/perf-108/structural/browser-smoke.IsZngp/library
evidence="$repo/docs/qa/workflow-execution/perf-108/structural"
family=structural
if [[ "${1:-structural}" == ordered ]]; then
    family=ordered
    root=/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/perf-108/ordered-groups/runtime
    library=/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/perf-108/ordered-groups/browser-smoke.jEmE2F/library
    evidence="$repo/docs/qa/workflow-execution/perf-108/ordered-groups"
fi
[[ "${1:-structural}" == structural || "${1:-structural}" == ordered ]]
[[ -f "$evidence/bindings.rds" && ! -e "$root/$family-screen" ]]
mkdir -p "$root"
run_stage() {
    local stage="$1"
    [[ ! -e "$root/$stage" && ! -e "$evidence/$stage.log" ]]
    printf 'Starting %s\n' "$stage"
    Rscript --vanilla tools/profiling/workflow_perf106_repair.R "$stage" "$library" "$root" > "$evidence/$stage.log" 2>&1
}
run_stage "$family-screen"
if ! Rscript --vanilla -e 'x<-readRDS(commandArgs(TRUE)[[1]]); quit(status=if(isTRUE(x$science_outputs_and_measurements_valid)) 0L else 1L)' "$root/$family-screen/decision.rds"; then
    printf 'Screen gate did not pass; full confirmation not launched.\n'
    exit 0
fi
run_stage "$family-rollup"
if ! Rscript --vanilla -e 'x<-readRDS(commandArgs(TRUE)[[1]]); quit(status=if(identical(x$status,"rollup_candidate_supported")) 0L else 1L)' "$root/$family-rollup/decision.rds"; then
    printf 'Rollup did not qualify; consumer confirmations not launched.\n'
    exit 0
fi
run_stage "$family-primary"
run_stage "$family-wide"
printf 'Fixed confirmation schedule completed.\n'
