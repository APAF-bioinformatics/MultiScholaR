#!/usr/bin/env bash
# KIT-150 (branch run-mofa-measure-k150): the stages one job runs, in order. Declaration:
# docs/qa/workflow-execution/kit/kit-150/declaration.md ("Runs" and "Jobs").
#
# Usage: k150-plan.sh <plan> <size>
#   lifecycle        size S: the five builds' pipelines (floor and app-path rows for base and C), then three induced
#                    failures (memory limit, time limit, an error)
#   probe            the base build's pipeline at <size> (L), with the floor row
#   counted-small    S and M: the five builds once, base and C three times (interleaved); floor rows for base and C;
#                    app-path rows at M; the import pair (P1, P2, P3) on C three times
#   counted-large    <size> (L or Lp): the five builds once, floor and app-path rows for base and C; base and C twice
#                    more only if the job's remaining time holds twice their first pipelines' time
# Environment: WORK, RES, LANE, BUNDLE, LIBROOT (one library per build: $LIBROOT/<build>), DEADLINE (epoch seconds
# after which no stage starts; stage time limits are cut to it), K150_CG, K150_LOGS.
set -u -o pipefail

plan=$1
size_arg=${2:-}
: "${WORK:?}" "${RES:?}" "${LANE:?}" "${BUNDLE:?}" "${LIBROOT:?}" "${DEADLINE:?}" "${K150_CG:?}" "${K150_LOGS:?}"
STAGES="$RES/stages"
mkdir -p "$STAGES" "$K150_LOGS" "$WORK/args" "$WORK/runs"

MEMORY_LIMIT=15569256448 # 14.5 GiB: the stage is stopped inside its cgroup before the runner's 16 GB are exhausted
DISK_CEILING=10000000000 # 10 GB of measurement files per job
CAP=1000
FLOOR_ROUNDS=10
declare -A CEILING=([write]=1800 [import]=3600 [align]=3600 [train]=9000 [results]=1800 [floor]=5400
    [app_import]=3600 [app_preview]=3600 [fail_memory]=600 [fail_time]=240 [fail_error]=600)

# Sizes: samples, features per data set, data sets, companion features (import pair), seed.
declare -A N=([S]=20 [M]=200 [L]=500 [Lp]=300 [P1]=200 [P2]=200 [P3]=100)
declare -A P=([S]=5000 [M]=30000 [L]=60000 [Lp]=40000 [P1]=30000 [P2]=60000 [P3]=30000)
declare -A COMPANION=([S]=null [M]=null [L]=null [Lp]=null [P1]=1000 [P2]=1000 [P3]=1000)
declare -A SEED=([S]=15001 [M]=15002 [L]=15003 [Lp]=15004 [P1]=15011 [P2]=15012 [P3]=15013)
OMICS4='["proteomics","transcriptomics","lipidomics","metabolomics"]'
OMICS2='["proteomics","metabolomics"]'

declare -A COMMIT VARIANT
while IFS=$'\t' read -r label commit variant; do
    [ "$label" = build ] && continue
    COMMIT[$label]=$commit
    VARIANT[$label]=$variant
done < "$LANE/.github/mofa-smoke/k150-builds.tsv"

STOP=""
declare -A PIPELINE_S

skip() { # <id> <reason>
    jq -n --arg id "$1" --arg outcome "$2" '{id: $id, outcome: $outcome}' > "$STAGES/$1.driver.json"
    echo "k150 stage $1: $2"
}

# run_stage <stage> <build> <size> <rep> <case> <dir> [<memory limit> [<ceiling s>]]; 0 when the stage completed.
run_stage() {
    local stage=$1 build=$2 size=$3 rep=$4 case=$5 dir=$6
    local memory=${7:-$MEMORY_LIMIT} ceiling=${8:-${CEILING[$1]}}
    local id="${size}_${build}_r${rep}_${stage}"
    if [ -n "$STOP" ]; then skip "$id" "not run: $STOP"; return 2; fi
    local left=$((DEADLINE - $(date +%s)))
    if [ "$left" -lt 180 ]; then
        STOP="the job's time limit"
        skip "$id" "not run: $STOP"
        return 2
    fi
    local used
    used=$(du -sb "$WORK" | cut -f1)
    if [ "$used" -gt "$DISK_CEILING" ]; then
        STOP="the disk ceiling (${used} bytes in use)"
        skip "$id" "not run: $STOP"
        return 2
    fi
    local limit=$ceiling
    [ $((left - 60)) -lt "$limit" ] && limit=$((left - 60))
    local omics=$OMICS4
    [ "${COMPANION[$size]}" != null ] && omics=$OMICS2
    mkdir -p "$dir/tmp"
    jq -n --arg id "$id" --arg stage "$stage" --arg build "$build" --arg commit "${COMMIT[$build]}" \
        --arg variant "${VARIANT[$build]}" --arg lib "$LIBROOT/$build" --arg size "$size" --argjson rep "$rep" \
        --arg case "$case" --argjson n "${N[$size]}" --argjson p "${P[$size]}" --argjson companion "${COMPANION[$size]}" \
        --argjson seed "${SEED[$size]}" --argjson omics "$omics" --argjson cap "$CAP" --argjson floor_rounds "$FLOOR_ROUNDS" \
        --arg inputs "$dir/inputs" --arg store "$dir/${STORE_NAME:-store}" --arg out "$STAGES/$id" '{
        id: $id, stage: $stage, build: $build, commit: $commit, variant: $variant, lib: $lib, size: $size, rep: $rep,
        case: $case, n: $n, p: $p, companion_p: $companion, seed: $seed, omics: $omics, cap: $cap,
        floor_rounds: $floor_rounds, inputs: $inputs, store: $store, out: $out }' > "$WORK/args/$id.json"
    TMPDIR="$dir/tmp" bash "$LANE/.github/mofa-smoke/k150-run-stage.sh" "$WORK/args/$id.json" "$memory" "$limit"
}

remove_run() { # the run's folder, through the guarded remover
    python3 -I "$LANE/.github/mofa-smoke/k150-remove.py" "$WORK/runs" "$1"
}

# pipeline <size> <build> <rep> <floor 0|1> <app 0|1>: write, import, alignment, training, results read (capped);
# then the floor row and the app-path rows when asked. Each stage needs the ones it reads from to have completed.
pipeline() {
    local size=$1 build=$2 rep=$3 floor=$4 app=$5
    local dir="$WORK/runs/${size}_${build}_r${rep}"
    local t0
    t0=$(date +%s)
    mkdir -p "$dir"
    local wrote=0 imported=0 trained=0
    run_stage write "$build" "$size" "$rep" capped "$dir" && wrote=1
    if [ $wrote = 1 ]; then run_stage import "$build" "$size" "$rep" capped "$dir" && imported=1
    else skip "${size}_${build}_r${rep}_import" "not run: the input write did not complete"; fi
    if [ $imported = 1 ]; then
        run_stage align "$build" "$size" "$rep" capped "$dir"
        run_stage train "$build" "$size" "$rep" capped "$dir" && trained=1
    else
        for s in align train; do skip "${size}_${build}_r${rep}_$s" "not run: the import did not complete"; done
    fi
    if [ $trained = 1 ]; then run_stage results "$build" "$size" "$rep" capped "$dir"
    else skip "${size}_${build}_r${rep}_results" "not run: the training did not complete"; fi
    PIPELINE_S["${size}_${build}_r${rep}"]=$(($(date +%s) - t0))
    if [ "$floor" = 1 ]; then
        if [ $imported = 1 ]; then run_stage floor "$build" "$size" "$rep" floor "$dir"
        else skip "${size}_${build}_r${rep}_floor" "not run: the import did not complete"; fi
    fi
    if [ "$app" = 1 ]; then
        if [ $wrote = 1 ]; then
            if STORE_NAME=app-store run_stage app_import "$build" "$size" "$rep" capped "$dir"; then
                STORE_NAME=app-store run_stage app_preview "$build" "$size" "$rep" capped "$dir"
            else
                skip "${size}_${build}_r${rep}_app_preview" "not run: the app import did not complete"
            fi
        else
            for s in app_import app_preview; do skip "${size}_${build}_r${rep}_$s" "not run: the input write did not complete"; done
        fi
    fi
    remove_run "$dir"
    df -B1 --output=avail "$WORK" | tail -n 1 | awk '{ printf "k150: %.1f GB free after %s\n", $1 / 1e9, "'"${size}_${build}_r${rep}"'" }'
}

# The import pair (claim 3): write and import only, build C.
import_pair() { # <size> <rep>
    local size=$1 rep=$2
    local dir="$WORK/runs/${size}_C_r${rep}"
    mkdir -p "$dir"
    if run_stage write C "$size" "$rep" capped "$dir"; then run_stage import C "$size" "$rep" capped "$dir"
    else skip "${size}_C_r${rep}_import" "not run: the input write did not complete"; fi
    remove_run "$dir"
}

floor_app() { # <build>: base and C carry the floor and app-path rows
    case "$1" in base | C) echo 1 ;; *) echo 0 ;; esac
}

case "$plan" in
    lifecycle)
        for b in base A B1 B3 C; do f=$(floor_app "$b"); pipeline S "$b" 1 "$f" "$f"; done
        dir="$WORK/runs/S_C_r1_failures"
        mkdir -p "$dir"
        run_stage fail_memory C S 1 capped "$dir" 4294967296 600
        run_stage fail_time C S 1 capped "$dir" "$MEMORY_LIMIT" 240
        run_stage fail_error C S 1 capped "$dir"
        remove_run "$dir"
        ;;
    probe)
        pipeline "$size_arg" base 1 1 0
        ;;
    counted-small)
        for b in base A B1 B3 C; do f=$(floor_app "$b"); pipeline S "$b" 1 "$f" 0; done
        for r in 2 3; do pipeline S base "$r" 0 0; pipeline S C "$r" 0 0; done
        for b in base A B1 B3 C; do f=$(floor_app "$b"); pipeline M "$b" 1 "$f" "$f"; done
        for r in 2 3; do pipeline M base "$r" 0 0; pipeline M C "$r" 0 0; done
        for r in 1 2 3; do for s in P1 P2 P3; do import_pair "$s" "$r"; done; done
        ;;
    counted-large)
        for b in base A B1 B3 C; do f=$(floor_app "$b"); pipeline "$size_arg" "$b" 1 "$f" "$f"; done
        need=$((2 * (${PIPELINE_S["${size_arg}_base_r1"]:-0} + ${PIPELINE_S["${size_arg}_C_r1"]:-0})))
        left=$((DEADLINE - $(date +%s)))
        if [ -z "$STOP" ] && [ "$need" -gt 0 ] && [ "$left" -gt $((need + 600)) ]; then
            for r in 2 3; do pipeline "$size_arg" base "$r" 0 0; pipeline "$size_arg" C "$r" 0 0; done
        else
            echo "k150: repeats at $size_arg not run (need ${need} s plus 600 s, ${left} s left; ${STOP:-no stop})" |
                tee "$RES/repeats-not-run.txt"
        fi
        ;;
    g0i-repeat)
        # KIT-209 addendum 1 (A7): one repeat of the S pipeline on RC and on C (the same commits as G0i).
        pipeline S RC 2 0 0
        pipeline S C 2 0 0
        ;;
    m3-M)
        # KIT-209 M3: I-M, 3 per build, order C RC / RC C / C RC (the base of integration is build C).
        pipeline M C 1 0 0; pipeline M RC 1 0 0
        pipeline M RC 2 0 0; pipeline M C 2 0 0
        pipeline M C 3 0 0; pipeline M RC 3 0 0
        ;;
    m3-L)
        # KIT-209 M3: I-L, 2 per build, order C RC / RC C.
        pipeline L C 1 0 0; pipeline L RC 1 0 0
        pipeline L RC 2 0 0; pipeline L C 2 0 0
        ;;
    g0i)
        # KIT-209 G0i: the S pipeline on the head (RC) and on build C, then a memory-limit failure on RC.
        pipeline S RC 1 0 0
        pipeline S C 1 0 0
        dir="$WORK/runs/S_RC_r1_failures"
        mkdir -p "$dir"
        run_stage fail_memory RC S 1 capped "$dir" 4294967296 600
        remove_run "$dir"
        ;;
    *)
        echo "unknown plan $plan" >&2
        exit 2
        ;;
esac
echo "k150: plan $plan done${STOP:+ (stopped: $STOP)}"
