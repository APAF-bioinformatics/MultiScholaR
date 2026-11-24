#!/usr/bin/env bash
# KIT-209 (branch run-kit209-memory): the journeys one job runs, in order. Declaration:
# docs/qa/workflow-execution/kit/kit-209/declaration.md ("G0").
#
# Usage: k209-plan.sh <plan>
#   g0a   corpora (screen), then: the lifecycle journey (lipid_screen, candidate C); three injected failures relative to
#         it (memory.max = 0.7 x its journey VmHWM; time limit = 0.5 x its process wall time; a package error from an
#         invalid bound setting); the reset control (calibration, candidate); gate B: one screen journey per route on
#         the base then the candidate (U1 to U5, S1, S2), with fingerprints.
#   g0p1  the time probes U10 (DIA_LIMPA, fastlimpadpc) and U9 (DIA) on the candidate, one journey each, to the deadline.
#   g0p2  the time probes U7 (metab_large_missing), U11 (lfq_300) and U12 (tmt_300) on the candidate, likewise.
# Environment: WORK, RES, LANE, CHECKOUT (the candidate's checkout), LIBROOT (one library per build: $LIBROOT/<build>),
#   DEADLINE (epoch seconds after which no journey starts; journey time limits are cut to it), K209_CG, K209_LOGS.
set -u -o pipefail

plan=$1
: "${WORK:?}" "${RES:?}" "${LANE:?}" "${CHECKOUT:?}" "${LIBROOT:?}" "${DEADLINE:?}" "${K209_CG:?}" "${K209_LOGS:?}"
JOURNEYS="$RES/journeys"
mkdir -p "$JOURNEYS" "$K209_LOGS" "$WORK/args" "$WORK/runs" "$WORK/corpus"

MEMORY_LIMIT=15569256448 # 14.5 GiB: the journey is stopped inside its cgroup before the runner's 16 GB are exhausted
SCREEN_LIMIT=1800        # a screen journey's time limit in G0a (seconds)

declare -A COMMIT
while IFS=$'\t' read -r label commit; do
    [ "$label" = build ] && continue
    COMMIT[$label]=$commit
done < "$LANE/.github/k209/k209-builds.tsv"

# unit -> route, corpus, format
declare -A ROUTE=([U1]=lipidomics [U2]=metabolomics [U3]=transcriptomics [U4]=DIA [U5]=DIA_LIMPA [S1]=protein [S2]=protein
    [CAL]=calibration [U7]=metabolomics [U9]=DIA [U10]=DIA_LIMPA [U11]=protein [U12]=protein)
declare -A CORPUS=([U1]=lipid_screen [U2]=metab_screen [U3]=trx_screen [U4]=screen [U5]=screen [S1]=lfq_fixture [S2]=tmt_fixture
    [CAL]=none [U7]=metab_large_missing [U9]=full [U10]=full [U11]=lfq_300 [U12]=tmt_300)
declare -A FORMAT=([U4]=diann [U5]=diann [S1]=maxquant [S2]=pd_tmt [U9]=diann [U10]=diann [U11]=maxquant [U12]=pd_tmt)
TO_DEADLINE=86400        # G0p: a probe runs to the job's deadline (cut below)

STOP=""

skip() { # <id> <reason>
    jq -n --arg id "$1" --arg outcome "$2" '{id: $id, outcome: $outcome}' > "$JOURNEYS/$1.driver.json"
    echo "k209 journey $1: $2"
}

corpus() { # <key>: generate once per job
    local key=$1
    [ -d "$WORK/corpus/$key" ] && return 0
    /opt/R/4.6.1/bin/Rscript --vanilla "$LANE/.github/k209/k209-corpus.R" "$key" "$WORK/corpus/$key" "$CHECKOUT" "$RES/corpora.tsv" \
        > "$K209_LOGS/corpus-$key.log" 2>&1 || { echo "k209: corpus $key failed"; tail -n 20 "$K209_LOGS/corpus-$key.log"; return 1; }
}

# journey <unit> <build> <rep> <order> <inject> <fingerprint true|false> [<memory limit> [<time limit s>]]; 0 when completed.
journey() {
    local unit=$1 build=$2 rep=$3 order=$4 inject=$5 fingerprint=$6
    local memory=${7:-$MEMORY_LIMIT} ceiling=${8:-$SCREEN_LIMIT}
    local id="${unit}_${build}_r${rep}${inject:+_$inject}"
    if [ -n "$STOP" ]; then skip "$id" "not run: $STOP"; return 2; fi
    local left=$((DEADLINE - $(date +%s)))
    if [ "$left" -lt 180 ]; then STOP="the job's time limit"; skip "$id" "not run: $STOP"; return 2; fi
    local limit=$ceiling
    [ $((left - 60)) -lt "$limit" ] && limit=$((left - 60))
    local dir="$WORK/runs/$id"
    mkdir -p "$dir/tmp"
    jq -n --arg id "$id" --arg unit "$unit" --arg route "${ROUTE[$unit]}" --arg corpus "${CORPUS[$unit]}" \
        --arg format "${FORMAT[$unit]:-}" --arg corpus_dir "$WORK/corpus/${CORPUS[$unit]}" --arg build "$build" \
        --arg commit "${COMMIT[$build]}" --arg lib "$LIBROOT/$build" --argjson rep "$rep" --arg order "$order" \
        --arg inject "$inject" --argjson fingerprint "$fingerprint" --arg work "$dir/work" \
        --arg harness "$CHECKOUT/tools/profiling/workflow_kit.R" --arg out "$JOURNEYS/$id" '{
        id: $id, unit: $unit, route: $route, corpus: $corpus, format: $format, corpus_dir: $corpus_dir, build: $build,
        commit: $commit, lib: $lib, rep: $rep, order: $order, inject: (if $inject == "" then null else $inject end),
        fingerprint: $fingerprint, work: $work, harness: $harness, out: $out }' > "$WORK/args/$id.json"
    TMPDIR="$dir/tmp" bash "$LANE/.github/k209/k209-run-journey.sh" "$WORK/args/$id.json" "$memory" "$limit"
    local status=$?
    python3 -I "$LANE/.github/k209/k209-remove.py" "$WORK/runs" "$dir"
    return $status
}

case "$plan" in
    g0a)
        for key in lipid_screen metab_screen trx_screen screen lfq_fixture tmt_fixture; do corpus "$key" || STOP="corpus $key"; done
        # 1. The lifecycle journey, then the three injected failures relative to it.
        if journey U1 C 1 lifecycle "" true; then
            hwm=$(jq '[.marks[].vmhwm_mib] | max' "$JOURNEYS/U1_C_r1.json")
            wall=$(jq '.process_wall_s' "$JOURNEYS/U1_C_r1.driver.json")
            mem=$(python3 -I -c "import sys; print(int(0.7 * float(sys.argv[1]) * 1048576))" "$hwm")
            tl=$(python3 -I -c "import sys, math; print(max(1, math.ceil(0.5 * float(sys.argv[1]))))" "$wall")
            jq -n --argjson hwm "$hwm" --argjson wall "$wall" --argjson mem "$mem" --argjson tl "$tl" \
                '{lifecycle_journey_vmhwm_mib: $hwm, lifecycle_process_wall_s: $wall, injected_memory_max_bytes: $mem,
                  injected_time_limit_s: $tl}' | tee "$RES/injected-limits.json"
            journey U1 C 1 lifecycle fail_memory false "$mem"
            journey U1 C 1 lifecycle fail_time false "$MEMORY_LIMIT" "$tl"
        else
            skip "U1_C_r1_fail_memory" "not run: the lifecycle journey did not complete"
            skip "U1_C_r1_fail_time" "not run: the lifecycle journey did not complete"
        fi
        journey U1 C 1 lifecycle package_error false
        # 2. The reset control.
        journey CAL C 1 control "" false
        # 3. Gate B: one screen journey per route, base then candidate.
        for unit in U1 U2 U3 U4 U5 S1 S2; do
            journey "$unit" base 1 BC "" true
            journey "$unit" C 2 BC "" true
        done
        ;;
    g0a-error)
        # G0a item 1's package error alone (VAL-54 instrument run 2 of G0a: the injected setting changed, section 6).
        corpus lipid_screen || STOP="corpus lipid_screen"
        journey U1 C 1 lifecycle package_error false
        ;;
    g0p1)
        corpus full || STOP="corpus full"
        journey U10 C 1 probe "" false "$MEMORY_LIMIT" "$TO_DEADLINE"
        journey U9 C 1 probe "" false "$MEMORY_LIMIT" "$TO_DEADLINE"
        ;;
    g0p2)
        for key in metab_large_missing lfq_300 tmt_300; do corpus "$key" || STOP="corpus $key"; done
        journey U7 C 1 probe "" false "$MEMORY_LIMIT" "$TO_DEADLINE"
        journey U11 C 1 probe "" false "$MEMORY_LIMIT" "$TO_DEADLINE"
        journey U12 C 1 probe "" false "$MEMORY_LIMIT" "$TO_DEADLINE"
        ;;
    *)
        echo "unknown plan $plan" >&2
        exit 2
        ;;
esac
echo "k209: plan $plan done${STOP:+ (stopped: $STOP)}"
