#!/usr/bin/env bash
set -euo pipefail
# One authorized, non-retrying correction cycle through the existing runner.
repo=$(pwd -P)
proof=/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/perf-108/structural
evidence="$repo/docs/qa/workflow-execution/perf-108/structural"
base=/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/perf-108/cache-qualification/browser-smoke.VSrEGM/library
[[ -f "$repo/DESCRIPTION" && -d "$proof/bulk-rollup-library" ]]
[[ ! -e "$evidence/O-native-build-2.log" && ! -e "$evidence/O-app-build-2.log" && ! -e "$evidence/O-check-2.log" ]]
started=$SECONDS
export MULTISCHOLAR_PERF108_REMAINING_BUILD_SECONDS=600
bash tools/profiling/run_workflow_checkpoint_proof.sh --service reader-check "$base" "$proof" "$repo" perf108:K1 > "$evidence/O-native-build-2.log" 2>&1
native_unit=$(head -n 1 "$evidence/O-native-build-2.log" | python3 -c 'import sys,json; print(json.load(sys.stdin)["unit"])')
native_run=${native_unit#multischolar-dia003-}
native_run=${native_run%.service}
[[ "$native_run" =~ ^browser-smoke\.[A-Za-z0-9]{6}$ ]]
MULTISCHOLAR_PERF108_REMAINING_BUILD_SECONDS=$((600 - (SECONDS - started)))
(( MULTISCHOLAR_PERF108_REMAINING_BUILD_SECONDS > 0 ))
bash tools/profiling/run_workflow_checkpoint_proof.sh --service reader-check "$proof/$native_run/library" "$proof" "$repo" perf108:M1 > "$evidence/O-app-build-2.log" 2>&1
(( SECONDS - started <= 600 ))
app_unit=$(head -n 1 "$evidence/O-app-build-2.log" | python3 -c 'import sys,json; print(json.load(sys.stdin)["unit"])')
app_run=${app_unit#multischolar-dia003-}
app_run=${app_run%.service}
[[ "$app_run" =~ ^browser-smoke\.[A-Za-z0-9]{6}$ ]]
unset MULTISCHOLAR_PERF108_REMAINING_BUILD_SECONDS
bash tools/profiling/run_workflow_checkpoint_proof.sh --service reader-check "$proof/$app_run/library" "$proof" "$repo" perf108:ZO > "$evidence/O-check-2.log" 2>&1
printf 'Correction complete. Installed library: %s\n' "$proof/$app_run/library"
