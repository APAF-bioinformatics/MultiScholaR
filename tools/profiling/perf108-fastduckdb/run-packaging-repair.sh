#!/usr/bin/env bash
set -euo pipefail
repo=$(pwd -P)
proof=/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/perf-108/packaging-repair
library=/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/perf-108/ordered-groups/browser-smoke.jEmE2F/library
evidence="$repo/docs/qa/workflow-execution/perf-108/packaging-repair"
[[ ! -e "$evidence/build-1.log" && ! -e "$evidence/check-1.log" ]]
mkdir -p "$proof/bulk-rollup-library"
printf 'Starting clean source installation (600s ceiling).\n'
bash tools/profiling/run_workflow_checkpoint_proof.sh --service reader-check "$library" "$proof" "$repo" perf108:UB > "$evidence/build-1.log" 2>&1
unit=$(head -n 1 "$evidence/build-1.log" | python3 -c 'import json,sys; print(json.load(sys.stdin)["unit"])')
member=${unit#multischolar-dia003-}
member=${member%.service}
[[ "$member" =~ ^browser-smoke\.[A-Za-z0-9]{6}$ ]]
printf 'Starting fresh installed checks (300s ceiling).\n'
bash tools/profiling/run_workflow_checkpoint_proof.sh --service reader-check "$proof/$member/library" "$proof" "$repo" perf108:UT > "$evidence/check-1.log" 2>&1
printf 'Clean installation and installed checks completed: %s\n' "$proof/$member/library"
