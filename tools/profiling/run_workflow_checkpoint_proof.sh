#!/usr/bin/env bash
set -euo pipefail
# Diagnostic-only budgets. Missing observations fail closed; never infer zero.
checkpointDiagnosticGuard() {
    local available="$1" current="$2"
    if [[ ! "$available" =~ ^[0-9]{1,15}$ || ! "$current" =~ ^[0-9]{1,15}$ ]]; then
        echo "missing_memory_observation"; return 1
    fi
    if (( available < 8589934592 )); then echo "host_headroom"; return 1; fi
    if (( current >= 1610612736 )); then echo "job_early_stop"; return 1; fi
}
checkpointRecordedGuard() {
    local available="$1" current="$2" peak="$3" phase="$4"
    if [[ ! "$available" =~ ^[0-9]{1,15}$ || ! "$current" =~ ^[0-9]{1,15}$ ||
        ! "$peak" =~ ^[0-9]{1,15}$ ]]; then
        echo "missing_memory_observation"; return 1
    fi
    if [[ "$phase" == start ]] && (( available < 10737418240 )); then
        echo "launch_headroom"; return 1
    fi
    if (( available < 2147483648 )); then echo "host_headroom"; return 1; fi
    if (( current >= 8053063680 || peak >= 8053063680 )); then
        echo "job_early_stop"; return 1
    fi
}
# A missing cgroup counter is acceptable only after independently observed exit.
checkpointServiceTerminal() {
    [[ "$1" == *"SubState=exited"* || "$1" == *"SubState=failed"* ]]
}
checkpointMemorySample() {
    local path="$1" stamp="$2" current="$3" peak="$4"
    awk -v stamp="$stamp" -v current="$current" -v peak="$peak" '
        BEGIN { anon=file=kernel=dirty=writeback=-1 }
        $1=="anon" {anon=$2} $1=="file" {file=$2} $1=="kernel" {kernel=$2}
        $1=="file_dirty" {dirty=$2} $1=="file_writeback" {writeback=$2}
        END { if (anon<0 || file<0 || kernel<0 || dirty<0 || writeback<0) exit 2;
            printf "%s\t%s\t%s\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\n",
                stamp,current,peak,anon,file,kernel,dirty,writeback }
    ' "$path"
}
# Individual RSS/PSS are attribution observations, never summed into cgroup peak.
checkpointProcessSample() {
    local group="$1" stamp="$2" output="$3" pid
    while IFS= read -r pid; do
        [[ "$pid" =~ ^[0-9]+$ ]] || continue
        awk -v stamp="$stamp" -v pid="$pid" '
            BEGIN { ppid=rss=anon=file=shmem=pss=private_clean=private_dirty=threads=-1; name="unknown" }
            /^Name:/ {name=$2} /^PPid:/ {ppid=$2} /^VmRSS:/ {rss=$2*1024}
            /^RssAnon:/ {anon=$2*1024} /^RssFile:/ {file=$2*1024}
            /^RssShmem:/ {shmem=$2*1024} /^Pss:/ {pss=$2*1024}
            /^Private_Clean:/ {private_clean=$2*1024} /^Private_Dirty:/ {private_dirty=$2*1024}
            /^Threads:/ {threads=$2}
            END { printf "%s\t%s\t%s\t%s\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%s\n",
                stamp,pid,ppid,name,rss,anon,file,shmem,pss,private_clean,private_dirty,threads }
        ' "/proc/$pid/status" "/proc/$pid/smaps_rollup" >> "$output" 2>> "${output%.tsv}-gaps.log" || true
    done < "$group/cgroup.procs" || printf '%s\tcgroup_process_list_unavailable\n' "$stamp" >> "${output%.tsv}-gaps.log"
}
# Freeze the small PERF-101 execution sources before launching its guarded process.
# Data hashes, package bindings and source-manifest hashes are checked in that process.
checkpointSnapshotPerf101Sources() {
    local source_repo="$1" snapshot="$2" source_file
    [[ ! -e "$snapshot" ]] || return 2
    mkdir -p -- "$snapshot"
    for source_file in tools/profiling/run_workflow_checkpoint_proof.R \
        tools/profiling/run_workflow_checkpoint_proof.sh \
        tools/profiling/workflow_checkpoint_dia_driver.R tools/profiling/workflow_perf101.R \
        tests/testthat/helper-workflow-bulk-rollup.R \
        tests/testdata/workflow-execution/qualification.json \
        tickets/openai/verification/pilot-source-bindings.json \
        tickets/openai/verification/pilot-protocol.json R/func_prot_qc_peptide_methods.R; do
        [[ -f "$source_repo/$source_file" && ! -L "$source_repo/$source_file" ]] || return 2
        mkdir -p -- "$(dirname -- "$snapshot/$source_file")"
        cp -- "$source_repo/$source_file" "$snapshot/$source_file" || return 2
    done
}

proof_script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
if [[ "${1:-}" == "--workflow-lifecycle" ]]; then
    [[ "$#" == 5 && ( "$5" == protocol || "$5" == resources ) ]] || exit 2
    exec Rscript --vanilla "${proof_script_dir}/run_workflow_checkpoint_proof.R" "$@"
fi
if [[ "${1:-}" == "--service" ]]; then
    proof_qvalue=false
    proof_capped_sixteen=false
    proof_ruv_pearson=false
    proof_ruv_small=false
    proof_ruv_stress=false
    proof_stress_samples=""
    proof_native_session=false
    proof_first_qc=false
    proof_reclaim=false
    proof_design_memory=false
    proof_extended120=false
    if [[ "${2:-}" == browser-baseline && "${!#}" == --runtime-120-minutes ]]; then
        proof_extended120=true
        set -- "${@:1:$(($# - 1))}"
    fi
    if [[ "${2:-}" == browser-baseline && "${!#}" == --diagnostic-design-memory ]]; then
        proof_design_memory=true
        set -- "${@:1:$(($# - 1))}"
    fi
    if [[ "${2:-}" == browser-baseline && "${!#}" == --reclaim-memory ]]; then
        proof_reclaim=true
        set -- "${@:1:$(($# - 1))}"
    fi
    if [[ "${2:-}" == browser-baseline && "${!#}" == --first-qc ]]; then
        proof_first_qc=true
        set -- "${@:1:$(($# - 1))}"
    fi
    proof_timed_ruv=false
    proof_recorded=false
    proof_reader=false
    proof_memory_guard=false
    proof_utility_diagnostics=false
    proof_sequence_sampling=false
    proof_validate_only=false
    proof_capacity_treatment=""
    proof_capacity_recipe=""
    proof_capacity_reference=""
    proof_gui_diagnostic=false
    proof_process_sampling=false
    if [[ "${2:-}" == browser-baseline || "${2:-}" == verify-da ]]; then
        if [[ "${!#}" == --native-session ]]; then
            proof_native_session=true
            set -- "${@:1:$(($# - 1))}"
        fi
        if [[ "${!#}" == --timed-ruv-pearson ]]; then
            proof_timed_ruv=true
            set -- "${@:1:$(($# - 1))}" --diagnostic-ruv-pearson
        elif [[ "${!#}" =~ ^--stress-(96|300|500)$ ]]; then
            [[ "$proof_native_session" == true ]] || exit 2
            proof_ruv_stress=true
            proof_stress_samples="${!#}"
            proof_stress_samples="${proof_stress_samples#--stress-}"
            proof_capped_sixteen=true
            proof_timed_ruv=true
            set -- "${@:1:$(($# - 1))}" --diagnostic-ruv-pearson
        fi
    fi
    if [[ "${2:-}" == browser-baseline && "$#" == 9 && "${8:-}" == --capacity-qs2_duckdb && "${9:-}" == --diagnostic-ruv-pearson-small ]]; then
        proof_ruv_pearson=true
        proof_ruv_small=true
        proof_gui_diagnostic=true
        set -- "${@:1:8}"
    elif [[ "${2:-}" == verify-da && "$#" == 6 && "${6:-}" == --diagnostic-ruv-pearson-small ]]; then
        proof_ruv_pearson=true
        proof_ruv_small=true
        proof_memory_guard=true
        set -- "${@:1:5}"
    fi
    if [[ "${2:-}" == verify-da && "$#" == 6 && "${6:-}" == --diagnostic-ruv-pearson ]]; then
        proof_ruv_pearson=true
        proof_memory_guard=true
        set -- "${@:1:5}"
    fi
    if [[ "${2:-}" == browser-baseline && "$#" == 9 && "${8:-}" == --capacity-qs2_duckdb && "${9:-}" == --diagnostic-ruv-pearson ]]; then
        proof_ruv_pearson=true
        proof_gui_diagnostic=true
        set -- "${@:1:8}"
    fi
    if [[ "${2:-}" == browser-baseline && "$#" == 9 && "${8:-}" == --capacity-qs2_duckdb && "${9:-}" == --diagnostic-sixteen-gib ]]; then
        proof_capped_sixteen=true
        proof_gui_diagnostic=true
        set -- "${@:1:8}"
    fi
    if [[ "${2:-}" == browser-baseline && "$#" == 8 && "${8:-}" =~ ^--capacity-(resident|qs2_duckdb)$ ]]; then
        proof_capacity_treatment="${8#--capacity-}"
        proof_memory_guard=true
        proof_sequence_sampling=true
        set -- "${@:1:7}"
    elif [[ "${2:-}" == prepare-companions && "$#" == 5 && "${5:-}" =~ ^capacity:(small|oom_reference|report_1gb|report_2gb|ruv_small|ruv_1gb|ruv_2gb|ruv_300|ruv_500)$ ]]; then
        proof_capacity_recipe="${5#capacity:}"
        proof_memory_guard=true
        set -- "${@:1:4}"
    elif [[ "${2:-}" == accounting-probe* && "$#" == 5 && "${5:-}" == --capacity ]]; then
        proof_memory_guard=true
        proof_sequence_sampling=true
        set -- "${@:1:4}"
    elif [[ "${2:-}" == verify-da && ( "$#" == 6 || "$#" == 7 ) && "${!#}" == --capacity ]]; then
        proof_memory_guard=true
        set -- "${@:1:$(($# - 1))}"
        if [[ "$#" == 6 ]]; then
            proof_capacity_reference="$(realpath -e -- "$6")"
            set -- "${@:1:5}"
        fi
    fi
    if [[ "${2:-}" == utility-member && "${!#}" == --diagnostic-memory ]]; then
        proof_utility_diagnostics=true
        set -- "${@:1:$(($# - 1))}"
    fi
    if [[ "${2:-}" == recorded-* || "${2:-}" == reader-* ]]; then
        proof_recorded=true
        proof_memory_guard=true
        if [[ "${2:-}" == reader-* ]]; then proof_reader=true; fi
        if [[ "${!#}" == --validate-only ]]; then
            proof_validate_only=true
            set -- "${@:1:$(($# - 1))}"
        fi
    fi
    if [[ "${2:-}" == recorded-qc && "${6:-}" == sequence:* ]]; then
        proof_sequence_sampling=true
    fi
    if [[ "${2:-}" == utility-prepare || "${2:-}" == utility-member || "${2:-}" == utility-verify || "${2:-}" == install-corrected ]]; then
        proof_memory_guard=true
    fi
    # Sequential QC extends the existing recorded operation, not a new service.
    if [[ ( "$2" == recorded-qc || "$2" == reader-check ) && "$#" == 6 ]]; then
        [[ ( "$2" == recorded-qc && ( "$6" == sequence-prepare || "$6" =~ ^sequence:(resident|qs2):(diagnostic|timing)$ || "$6" == sequence:qs2_duckdb:timing ) ) || ( "$2" == reader-check && ( "$6" =~ ^(perf108:(B[12]|T[12]|F0|N[12]|I[1-3]|J[12]|K[1-4]|M[12]|Q[CT]|U[BT]|Z[ACOVFSQL]|E[1-8]|R[CS][01]|P[AVO][01]|[YW](S[1-4]|[RPW][1-6])|PC|X(S[1-4]|R[1-6]|F[1-8])|A[SF]1|H[SF]1|C[23]|P[SF][12]|R[SFR][12]|V([1-9]|1[0-6]))|perf101:(D[01]|D0R|V)|perf106:(L[01]|S[1-6]|F[1-8]|R[1-4]|V2-(A|V|I|J|D[12]|Q([1-9]|1[0-9]|2[0-4]))))$ || "$6" == sequence || "$6" == gui-capacity || "$6" == gui-capacity-verifier || "$6" == gui-diagnostic || "$6" == gui-passive || "$6" == gui-lifetime-repair || "$6" == gui-memory-repair || "$6" == gui-memory-verification || "$6" == gui-audit-repair || "$6" == gui-audit-diagnostic || "$6" == gui-audit-member || "$6" == gui-audit-toggle || "$6" == gui-sleep-guard || "$6" == gui-integration-repair || "$6" == gui-worker-member || "$6" == gui-export-member || "$6" == gui-sharing-repair || "$6" == gui-native-session || "$6" == gui-native-outputs || "$6" == gui-reference-first-qc || "$6" == gui-reference-core || "$6" == gui-reference-routing || "$6" == gui-ruv-pearson || "$6" == gui-export-history-member || "$6" == gui-iq-history-member || "$6" == gui-mechanism-preflight || "$6" == gui-rollup-profile || "$6" == gui-raw-replay:small || "$6" == gui-raw-replay:outputs || "$6" == gui-raw-replay:verify || "$6" == gui-cohort-profile || "$6" =~ ^gui-bulk-rollup:(qualification-preflight|qualification-prepare|writer-build|writer-opt-(build|preflight|diagnostic:(primary|contrast))|writer-r-preflight|writer-r-readback:[A-Za-z0-9]{6}|writer-r:(screen|full):(primary|contrast):P1[456]:[12]|io-preflight|io-prepare|io:(duckdb|fastduckdb):[12]|scan-inspection|crossover:(B|P0):(4608|9600|96000|960000):[12]|contrast:(B|P([0-9]|1[0-4]))|diagnostic|diagnostic-preflight|preflight|relational-preflight|materialize-preflight|adaptation-preflight|fast-p3-preflight|writer-preflight|screen:([ABCD]|P([0-9]|1[0-4])):[12]|full:([ABCD]|P([0-9]|1[0-4])))$ || "$6" =~ ^(gui-cohort-size:(256|768|768-threads8)|gui-impute-workers:(preflight|4|8))$ || "$6" =~ ^gui-cohort-threads:(1|2|4|8)$ || "$6" == gui-raw-replay:300 || "$6" == gui-rollup-profile-preflight || "$6" =~ ^gui-rollup-compare:(control|candidate)$ || "$6" =~ ^gui-(rollup|session)-performance:(control|candidate)$ || "$6" =~ ^gui-first-qc-verify:[A-Za-z0-9]{6}:[A-Za-z0-9]{6}$ || "$6" =~ ^gui-history-verify:[A-Za-z0-9]{6}(:[A-Za-z0-9]{6})?$ ) ) ]] || exit 2
        proof_sequence_selector="$6"
        if [[ "$6" == gui-export-history-member || "$6" == gui-iq-history-member ]]; then
            proof_sequence_sampling=true
        fi
        set -- "${@:1:5}"
    else
        proof_sequence_selector=""
    fi
    # The bounded schema diagnostic uses the same arity and guard as preparation.
    if [[ "${2:-}" == reader-inspect && "$#" != 5 ]]; then exit 2; fi
    if [[ "$#" == 6 && "$2" == "exclusion-allocation" && "$6" == "qvalue" ]]; then
        proof_qvalue=true
        set -- "${@:1:5}"
    fi
    if [[ "$#" == 8 && "$2" == "browser-baseline" && "$8" == "--capped-sixteen-gib" ]]; then
        proof_capped_sixteen=true
        set -- "${@:1:7}"
    elif [[ "$#" == 6 && "$2" == "verify-da" && "$6" == "--capped-sixteen-gib" ]]; then
        proof_capped_sixteen=true
        set -- "${@:1:5}"
    fi
    if ! { [[ "$2" == prepare-companions && "$#" == 4 && -n "$proof_capacity_recipe" ]] || [[ "$2" == reader-inspect && "$#" == 5 ]] || [[ "$proof_reader" == true && ( ( "$#" == 5 && ( "$2" == reader-build || "$2" == reader-check || "$2" == reader-prepare || "$2" == reader-verify ) ) || ( "$#" == 6 && ( "$2" == reader-member || "$2" == reader-prepare ) ) ) ]] || [[ "$#" == 5 && ( "$2" == "recorded-prepare" || "$2" == "recorded-verify" || "$2" == "recorded-qc" || "$2" == "verify-normalization" || "$2" == "verify-da" || "$2" == "utility-verify" || "$2" == "exclusion-benchmark" || "$2" == "exclusion-memory" || "$2" == "exclusion-allocation" ) ]] || [[ "$#" == 6 && ( "$2" == "recorded-member" || "$2" == "recorded-timing" ) ]] || [[ "$#" == 4 && ( "$2" == "browser-smoke" || "$2" == "prepare-small" || "$2" == "install-corrected" || "$2" == "accounting-probe" || "$2" == "accounting-probe-failure" || "$2" == "accounting-probe-retained" ) ]] || [[ "$#" == 6 && ( "$2" == "browser-import" || "$2" == "prepare-companions" ) ]] || [[ "$#" == 7 && ( "$2" == "browser-design" || "$2" == "browser-peptides" || "$2" == "browser-proteins" || "$2" == "browser-normalization" || "$2" == "browser-da" || "$2" == "browser-da-accounted" || "$2" == "browser-da-default-qc" || "$2" == "browser-da-profile" || "$2" == "browser-baseline" ) ]] || [[ "$#" == 6 && "$2" == "utility-prepare" ]] || [[ ( "$#" == 9 || "$#" == 10 ) && "$2" == "utility-member" ]]; }; then
        echo "Usage: runner --service browser-smoke LIBRARY OUTPUT_ROOT; or browser-import LIBRARY OUTPUT_ROOT REPORT FASTA" >&2
        exit 2
    fi
    if [[ "$proof_first_qc" == true ]]; then
        [[ "$proof_native_session" == true && "$proof_capacity_treatment" == qs2_duckdb && "$proof_stress_samples" != 500 ]] || exit 2
    fi
    proof_library="$(realpath -e -- "$3")"
    if [[ -n "$proof_sequence_selector" ]]; then set -- "$@" "$proof_sequence_selector"; fi
    proof_root="$(realpath -e -- "$4")"
    proof_memory_max=8589934592
    proof_runtime_max=900
    if [[ "$proof_gui_diagnostic" == true || "$proof_ruv_pearson" == true ]]; then proof_runtime_max=3600; fi
    if [[ "$proof_ruv_small" == true ]]; then proof_runtime_max=900; fi
    if [[ "$proof_extended120" == true ]]; then
        [[ "$proof_stress_samples" == 300 && "$proof_native_session" == true &&
            "$proof_capacity_treatment" == qs2_duckdb && "$proof_reclaim" == true &&
            "$proof_first_qc" == false && "$proof_design_memory" == false ]] || exit 2
        proof_runtime_max=7200
    fi
    if [[ "$proof_recorded" == true ]]; then
        proof_input="$(realpath -e -- "$5")"
        if [[ "$2" != reader-build && "$2" != reader-check ]]; then
            [[ "$proof_input" == "$proof_root/"* ]] || { echo "Recorded input escapes proof root" >&2; exit 2; }
        fi
        proof_recorded_args=("$2" "$proof_input")
        if [[ "$2" == recorded-member || "$2" == recorded-timing || "$2" == reader-member || ( ( "$2" == reader-prepare || "$2" == recorded-qc || "$2" == reader-check ) && "$#" == 6 ) ]]; then proof_recorded_args+=("$6"); fi
        proof_validator=--validate-recorded
        if [[ "$proof_reader" == true ]]; then proof_validator=--validate-reader; fi
        Rscript --vanilla "${proof_script_dir}/run_workflow_checkpoint_proof.R" \
            "$proof_validator" "${proof_recorded_args[@]}"
        if [[ "$proof_validate_only" == true ]]; then
            Rscript --vanilla "${proof_script_dir}/run_workflow_checkpoint_proof.R" --check-service-paths "$proof_library" "$proof_root"
            exit 0
        fi
    fi
    if [[ "$proof_memory_guard" == true ]]; then
        if [[ ! "$proof_sequence_selector" =~ ^(perf108:(B[12]|T[12]|F0|N[12]|I[1-3]|J[12]|K[1-4]|M[12]|Q[CT]|U[BT]|Z[ACOVFSQL]|E[1-8]|R[CS][01]|P[AVO][01]|[YW](S[1-4]|[RPW][1-6])|PC|X(S[1-4]|R[1-6]|F[1-8])|A[SF]1|H[SF]1|C[23]|P[SF][12]|R[SFR][12]|V([1-9]|1[0-6]))|perf101:(D[01]|D0R|V)|perf106:(L[01]|S[1-6]|F[1-8]|R[1-4]|V2-(A|V|I|J|D[12]|Q([1-9]|1[0-9]|2[0-4]))))$ ]]; then checkpointRecordedGuard "$(awk '/^MemAvailable:/ {printf "%.0f", $2 * 1024}' /proc/meminfo)" 0 0 start; fi
        proof_stop_peak=8053063680
        proof_min_host_available=2147483648
    fi
    if [[ "$proof_sequence_selector" =~ ^(perf108:(B[12]|T[12]|F0|N[12]|I[1-3]|J[12]|K[1-4]|M[12]|Q[CT]|U[BT]|Z[ACOVFSQL]|E[1-8]|R[CS][01]|P[AVO][01]|[YW](S[1-4]|[RPW][1-6])|PC|X(S[1-4]|R[1-6]|F[1-8])|A[SF]1|H[SF]1|C[23]|P[SF][12]|R[SFR][12]|V([1-9]|1[0-6]))|perf101:(D[01]|D0R|V)|perf106:(L[01]|S[1-6]|F[1-8]|R[1-4]|V2-(A|V|I|J|D[12]|Q([1-9]|1[0-9]|2[0-4]))))$ ]]; then
        proof_memory_max=17179869184
        proof_stop_peak=16642998272
        proof_min_host_available=2147483648
        proof_runtime_max=600
        if [[ "$proof_sequence_selector" == perf108:K1 || "$proof_sequence_selector" == perf108:M1 ]]; then
            if [[ -n "${MULTISCHOLAR_PERF108_REMAINING_BUILD_SECONDS:-}" ]]; then
                [[ "$MULTISCHOLAR_PERF108_REMAINING_BUILD_SECONDS" =~ ^[1-9][0-9]{0,2}$ ]] || exit 2
                (( MULTISCHOLAR_PERF108_REMAINING_BUILD_SECONDS <= 600 )) || exit 2
                proof_runtime_max=$MULTISCHOLAR_PERF108_REMAINING_BUILD_SECONDS
            fi
        fi
        if [[ "$proof_sequence_selector" == perf108:UT || "$proof_sequence_selector" == perf108:Z[ACOVSQL] || "$proof_sequence_selector" == perf108:E[1-8] ]]; then proof_runtime_max=300; fi
        if [[ "$proof_sequence_selector" == perf108:P[AVO][01] ]]; then proof_runtime_max=300; fi
        if [[ "$proof_sequence_selector" == perf108:R[CS][01] ]]; then proof_runtime_max=300; fi
        if [[ "$proof_sequence_selector" == perf108:[YW]* ]]; then proof_runtime_max=300; fi
        if [[ "$proof_sequence_selector" == perf108:T* || "$proof_sequence_selector" == perf108:F0 || "$proof_sequence_selector" == perf108:N2 || "$proof_sequence_selector" == perf108:I[23] || "$proof_sequence_selector" == perf108:J2 || "$proof_sequence_selector" == perf108:K[24] || "$proof_sequence_selector" == perf108:M2 || "$proof_sequence_selector" == perf108:V* || "$proof_sequence_selector" == perf108:P[SF]* || "$proof_sequence_selector" == perf108:R[SFR]* || "$proof_sequence_selector" == perf108:C[23] ]]; then proof_runtime_max=300; fi
        if [[ "$proof_sequence_selector" == perf106:V2-J ]]; then proof_runtime_max=120; fi
        if [[ "$proof_sequence_selector" == perf106:V2-A || "$proof_sequence_selector" == perf106:V2-V ]]; then proof_runtime_max=300; fi
        if [[ "$proof_sequence_selector" == perf101:D0R || "$proof_sequence_selector" == perf101:V || "$proof_sequence_selector" == perf106:S* ]]; then proof_runtime_max=300; fi
        if [[ "$proof_sequence_selector" == perf106:F* || "$proof_sequence_selector" == perf106:R* ]]; then proof_runtime_max=1800; fi
        proof_sequence_sampling=true
        proof_process_sampling=true
        proof_available="$(awk '/^MemAvailable:/ {printf "%.0f", $2 * 1024}' /proc/meminfo)"
        if [[ ! "$proof_available" =~ ^[0-9]+$ ]] || (( proof_available < 4294967296 )); then
            echo "PERF-101 refused: host headroom is below 4 GiB" >&2
            exit 76
        fi
    fi
    if [[ "$proof_capped_sixteen" == true ]]; then
        proof_memory_max=17179869184
        proof_stop_peak=16642998272
        proof_min_host_available=2147483648
        proof_required_start_available=19327352832
        proof_start_headroom_gib=18
        if [[ "$proof_ruv_stress" == true && ( "$proof_stress_samples" == 96 || "$proof_stress_samples" == 300 ) ]]; then
            # User-approved startup-only exception for the final matched queue.
            proof_required_start_available=17179869184
            proof_start_headroom_gib=16
        fi
        proof_available="$(awk '/^MemAvailable:/ {printf "%.0f", $2 * 1024}' /proc/meminfo)"
        if [[ ! "$proof_available" =~ ^[0-9]+$ ]] ||
            (( proof_available < proof_required_start_available )); then
            echo "16 GiB diagnostic refused: host headroom is below $proof_start_headroom_gib GiB" >&2
            exit 76
        fi
    fi
    if [[ "$proof_sequence_selector" =~ ^gui-bulk-rollup:(qualification-preflight|qualification-prepare|writer-build|writer-opt-(build|preflight|diagnostic:(primary|contrast))|writer-r-preflight|writer-r-readback:[A-Za-z0-9]{6}|writer-r:(screen|full):(primary|contrast):P1[456]:[12]|io-preflight|io-prepare|io:(duckdb|fastduckdb):[12]|scan-inspection|crossover:(B|P0):(4608|9600|96000|960000):[12]|contrast:(B|P([0-9]|1[0-4]))|diagnostic|diagnostic-preflight|preflight|relational-preflight|materialize-preflight|adaptation-preflight|fast-p3-preflight|writer-preflight|screen:([ABCD]|P([0-9]|1[0-4])):[12]|full:([ABCD]|P([0-9]|1[0-4])))$ ]]; then
        proof_memory_max=17179869184
        proof_stop_peak=16642998272
        proof_min_host_available=2147483648
        proof_runtime_max=300
        if [[ "$proof_sequence_selector" == gui-bulk-rollup:full:* || "$proof_sequence_selector" == gui-bulk-rollup:writer-r:full:* ]]; then proof_runtime_max=1800; fi
        if [[ "$proof_sequence_selector" == gui-bulk-rollup:diagnostic || "$proof_sequence_selector" == gui-bulk-rollup:writer-opt-diagnostic:* ]]; then proof_runtime_max=900; fi
        if [[ "$proof_sequence_selector" == gui-bulk-rollup:qualification-prepare || "$proof_sequence_selector" == gui-bulk-rollup:writer-build || "$proof_sequence_selector" == gui-bulk-rollup:writer-opt-build || "$proof_sequence_selector" == gui-bulk-rollup:contrast:* ]]; then proof_runtime_max=1800; fi
        proof_sequence_sampling=true
        proof_process_sampling=true
        proof_available="$(awk '/^MemAvailable:/ {printf "%.0f", $2 * 1024}' /proc/meminfo)"
        if [[ ! "$proof_available" =~ ^[0-9]+$ ]] || (( proof_available < 17179869184 )); then
            echo "bulk rollup refused: host headroom is below 16 GiB" >&2
            exit 76
        fi
    fi
    if [[ "$proof_sequence_selector" == gui-raw-replay:300 || "$proof_sequence_selector" =~ ^(gui-cohort-size:(256|768|768-threads8)|gui-impute-workers:(preflight|4|8))$ ]]; then
        proof_memory_max=17179869184
        proof_stop_peak=16642998272
        proof_min_host_available=2147483648
        if [[ "$proof_sequence_selector" == gui-raw-replay:300 ]]; then proof_runtime_max=5400; fi
        proof_available="$(awk '/^MemAvailable:/ {printf "%.0f", $2 * 1024}' /proc/meminfo)"
        if [[ ! "$proof_available" =~ ^[0-9]+$ ]] || (( proof_available < 17179869184 )); then
            echo "300 replay refused: host headroom is below 16 GiB" >&2
            exit 76
        fi
    fi
    if [[ "$proof_ruv_pearson" == true && "$proof_ruv_small" != true && "$proof_ruv_stress" != true ]]; then
        proof_memory_max=19327352832
        proof_stop_peak=18790481920
        proof_min_host_available=2147483648
        proof_available="$(awk '/^MemAvailable:/ {printf "%.0f", $2 * 1024}' /proc/meminfo)"
        if [[ ! "$proof_available" =~ ^[0-9]+$ ]] || (( proof_available < 21474836480 )); then
            echo "18 GiB diagnostic refused: host headroom is below 20 GiB" >&2
            exit 76
        fi
    fi
    if [[ "$2" == "exclusion-allocation" ]]; then
        proof_memory_max=2147483648
        proof_runtime_max=120
        checkpointDiagnosticGuard "$(awk '/^MemAvailable:/ {printf "%.0f", $2 * 1024}' /proc/meminfo)" 0
    fi
    proof_browser=""
    if [[ "$proof_memory_guard" != true || -n "$proof_capacity_treatment" || "$proof_sequence_selector" == gui-diagnostic || "$proof_sequence_selector" == gui-passive || "$proof_sequence_selector" == gui-reference-first-qc || "$proof_sequence_selector" == gui-reference-routing || "$proof_sequence_selector" == gui-ruv-pearson || "$proof_sequence_selector" == gui-native-session || "$proof_sequence_selector" == gui-lifetime-repair || "$proof_sequence_selector" == gui-memory-repair || "$proof_sequence_selector" == gui-memory-verification ]]; then
        proof_browser="${CHROMOTE_CHROME:-$(command -v chromium || command -v google-chrome || command -v brave)}"
        [[ -x "$proof_browser" ]] || { echo "No browser executable found" >&2; exit 2; }
    fi
    Rscript --vanilla "${proof_script_dir}/run_workflow_checkpoint_proof.R" --check-service-paths "$proof_library" "$proof_root"
    # Keep this shell alive while the service runs so its descriptor owns the lease.
    # Never unlink the lock file: concurrent callers must share the same inode.
    command -v flock >/dev/null || { echo "flock is required for serial proof runs" >&2; exit 2; }
    exec 9>"${proof_root}/.workflow-checkpoint-proof.lock"
    flock --nonblock 9 || { echo "Another proof run owns this output root" >&2; exit 2; }
    proof_run_dir="$(mktemp -d -- "${proof_root}/browser-smoke.XXXXXX")"
    if [[ "$proof_sequence_selector" =~ ^gui-(rollup|session)-performance:(control|candidate)$ ]]; then
        Rscript --vanilla "${proof_script_dir}/run_workflow_checkpoint_proof.R" \
            --condition-mechanism-cache "$proof_root" "$proof_run_dir"
    elif [[ "$proof_timed_ruv" == true && "$2" == browser-baseline ]]; then
        Rscript --vanilla "${proof_script_dir}/run_workflow_checkpoint_proof.R" \
            --condition-gui-cache "$proof_root" "$proof_run_dir" "$5" "$6" "$7"
    fi
    proof_unit="multischolar-dia003-${proof_run_dir##*/}"
    proof_operation=(--internal-browser-smoke "$proof_library" "$proof_run_dir")
    if [[ "$proof_recorded" == true ]]; then
        proof_operation=("--internal-$2" "$proof_library" "$proof_run_dir" "$proof_input")
        if [[ "$2" == recorded-member || "$2" == recorded-timing || "$2" == reader-member || ( ( "$2" == reader-prepare || "$2" == recorded-qc || "$2" == reader-check ) && "$#" == 6 ) ]]; then proof_operation+=("$6"); fi
    elif [[ "$2" == "accounting-probe" || "$2" == "accounting-probe-failure" || "$2" == "accounting-probe-retained" ]]; then
        proof_operation=(--internal-accounting-probe)
        if [[ "$2" == "accounting-probe-failure" ]]; then
            proof_operation=(--internal-accounting-probe-failure)
        elif [[ "$2" == "accounting-probe-retained" ]]; then
            proof_operation=(--internal-accounting-probe-child)
        fi
    elif [[ "$2" == "exclusion-benchmark" || "$2" == "exclusion-memory" || "$2" == "exclusion-allocation" ]]; then
        proof_report="$(realpath -e -- "$5")"
        [[ -f "$proof_report" && "$proof_report" == "$proof_root/"* ]] || exit 2
        proof_operation=(--internal-exclusion-benchmark "$proof_library" "$proof_run_dir" "$proof_report")
        if [[ "$2" == "exclusion-memory" ]]; then
            proof_operation=(--internal-exclusion-memory "$proof_library" "$proof_run_dir" "$proof_report")
        elif [[ "$2" == "exclusion-allocation" ]]; then
            proof_operation=(--internal-exclusion-allocation "$proof_library" "$proof_run_dir" "$proof_report")
        fi
    elif [[ "$2" == "prepare-companions" ]]; then
        if [[ -n "$proof_capacity_recipe" ]]; then
            proof_operation=(--internal-prepare-capacity "$proof_library" "$proof_run_dir" "$proof_capacity_recipe")
        else
        [[ "$5" =~ ^[a-z][a-z0-9_]*$ && "$6" =~ ^[0-9a-f]{64}$ ]] || exit 2
        proof_report="$(realpath -e -- "$proof_root/$5/report.tsv")"
        [[ -f "$proof_report" && "$proof_report" == "$proof_root/"* ]] || exit 2
        proof_operation=(--internal-prepare-companions "$proof_library" "$proof_run_dir" "$5" "$proof_report" "$6")
        fi
    elif [[ "$2" == "install-corrected" ]]; then
        proof_repo="$(realpath -e -- "${proof_script_dir}/../..")"
        mkdir -- "$proof_run_dir/source" "$proof_run_dir/library"
        git -C "$proof_repo" rev-parse HEAD | tee "$proof_run_dir/source-commit.txt"
        git -C "$proof_repo" archive HEAD | tar -x -C "$proof_run_dir/source"
        git -C "$proof_repo" diff --binary HEAD -- R src configure DESCRIPTION NAMESPACE man/readConfigFile.Rd | tee "$proof_run_dir/runtime.patch" >/dev/null
        git -C "$proof_run_dir/source" apply "$proof_run_dir/runtime.patch"
        sha256sum "$proof_run_dir/runtime.patch" | tee "$proof_run_dir/runtime-patch.sha256"
        # Collate already includes these working-tree package files; an installed
        # proof must not silently omit them merely because they are untracked.
        while IFS= read -r -d '' proof_source_file; do
            [[ ( "$proof_source_file" =~ ^R/[A-Za-z0-9_]+\.R$ ||
                "$proof_source_file" =~ ^src/[A-Za-z0-9_]+\.[ch]$ ||
                "$proof_source_file" == src/Makevars || "$proof_source_file" == src/duckdb-LICENSE ||
                "$proof_source_file" == src/build_dia_rsum.R || "$proof_source_file" == configure ) &&
                -f "$proof_repo/$proof_source_file" && ! -L "$proof_repo/$proof_source_file" ]] || exit 2
            mkdir -p -- "$(dirname -- "$proof_run_dir/source/$proof_source_file")"
            cp -- "$proof_repo/$proof_source_file" "$proof_run_dir/source/$proof_source_file"
            sha256sum "$proof_run_dir/source/$proof_source_file" >> "$proof_run_dir/untracked-runtime.sha256"
        done < <(git -C "$proof_repo" ls-files --others --exclude-standard -z -- R src configure)
        if [[ -d "$proof_library/fastduckdb" ]]; then
            cp -a -- "$proof_library/fastduckdb" "$proof_run_dir/library/"
        fi
        # Isolated, pinned experimental dependencies accompany the new installed build.
        for proof_dependency in duckplyr nanoarrow; do
            if [[ -d "$proof_root/bulk-rollup-library/$proof_dependency" ]]; then
                cp -a -- "$proof_root/bulk-rollup-library/$proof_dependency" "$proof_run_dir/library/"
            fi
        done
        Rscript --vanilla -e 'args <- commandArgs(TRUE); jsonlite::write_json(list(baseline_library = args[[1L]], source = normalizePath(file.path(args[[2L]], "source"))), file.path(args[[2L]], "library-origin.json"), auto_unbox = TRUE, pretty = TRUE)' "$proof_library" "$proof_run_dir"
    elif [[ "$2" == "utility-prepare" ]]; then
        proof_contract="$(realpath -e -- "$5")"
        [[ -f "$proof_contract" && "$6" =~ ^[a-z][a-z0-9_]*$ ]] || exit 2
        proof_operation=(--internal-utility-prepare "$proof_library" "$proof_run_dir"
            "$proof_contract" "$6")
    elif [[ "$2" == "utility-member" ]]; then
        proof_contract="$(realpath -e -- "$5")"
        [[ -f "$proof_contract" && "$6" =~ ^[a-z][a-z0-9_]*$ &&
            ( "$7" == "resident" || "$7" == "saved" || "$7" == "qs2" ) &&
            ( "$8" == "independent" || "$8" == "shared_metadata" || "$8" == "retained_alias" ) &&
            "$9" =~ ^[1-9][0-9]*$ ]] || exit 2
        proof_operation=(--internal-utility-member "$proof_library" "$proof_run_dir"
            "$proof_contract" "$6" "$7" "$8" "$9")
        if [[ "$#" == 10 ]]; then
            proof_prepared="$(realpath -e -- "${10}")"
            [[ -f "$proof_prepared/utility-inputs.json" && "$proof_prepared" == "$proof_root/"* ]] || exit 2
            proof_operation+=("$proof_prepared")
        fi
    elif [[ "$2" == "verify-normalization" || "$2" == "verify-da" || "$2" == "utility-verify" ]]; then
        proof_member="$(realpath -e -- "$5")"
        [[ -d "$proof_member" && "$proof_member" == "$proof_root/"* ]] || exit 2
        proof_operation=("--internal-$2" "$proof_library" "$proof_run_dir" "$proof_member")
    elif [[ "$2" == "prepare-small" ]]; then
        proof_operation=(--internal-prepare-small "$proof_library" "$proof_run_dir")
    elif [[ "$2" != "browser-smoke" ]]; then
        proof_report="$(realpath -e -- "$5")"
        proof_fasta="$(realpath -e -- "$6")"
        [[ -f "$proof_report" && -f "$proof_fasta" ]] || exit 2
        proof_operation=(--internal-browser-import "$proof_library" "$proof_run_dir" "$proof_report" "$proof_fasta")
        if [[ "$2" == "browser-design" || "$2" == "browser-peptides" || "$2" == "browser-proteins" || "$2" == "browser-normalization" || "$2" == "browser-da" || "$2" == "browser-da-accounted" || "$2" == "browser-da-default-qc" || "$2" == "browser-da-profile" || "$2" == "browser-baseline" ]]; then
            proof_design="$(realpath -e -- "$7")"
            [[ -f "$proof_design" ]] || exit 2
            proof_operation=("--internal-$2" "$proof_library" "$proof_run_dir" "$proof_report" "$proof_fasta" "$proof_design")
            if [[ "$2" == "browser-da-accounted" || "$2" == "browser-da-default-qc" || "$2" == "browser-da-profile" || "$2" == "browser-baseline" ]]; then
                proof_operation=(--internal-browser-da "$proof_library" "$proof_run_dir" "$proof_report" "$proof_fasta" "$proof_design")
                if [[ "$2" == "browser-da-default-qc" || "$2" == "browser-da-profile" || "$2" == "browser-baseline" ]]; then
                    proof_operation=(--internal-browser-da-default-qc "$proof_library" "$proof_run_dir" "$proof_report" "$proof_fasta" "$proof_design")
                    if [[ "$2" == "browser-da-profile" ]]; then
                        proof_operation=(--internal-browser-da-profile "$proof_library" "$proof_run_dir" "$proof_report" "$proof_fasta" "$proof_design")
                    elif [[ "$2" == "browser-baseline" ]]; then
                        proof_operation=(--internal-browser-baseline "$proof_library" "$proof_run_dir" "$proof_report" "$proof_fasta" "$proof_design")
                    fi
                fi
            fi
        fi
    fi
    if [[ -n "$proof_capacity_treatment" ]]; then
        proof_capacity_argument="capacity:$proof_capacity_treatment"
        if [[ "$proof_ruv_stress" == true ]]; then proof_capacity_argument+=:ruv-pearson-"$proof_stress_samples"
        elif [[ "$proof_ruv_small" == true ]]; then proof_capacity_argument+=:ruv-pearson-small
        elif [[ "$proof_ruv_pearson" == true ]]; then proof_capacity_argument+=:ruv-pearson
        elif [[ "$proof_gui_diagnostic" == true ]]; then proof_capacity_argument+=:diagnostic; fi
        if [[ "$proof_timed_ruv" == true ]]; then proof_capacity_argument+=:timed; fi
        if [[ "$proof_native_session" == true ]]; then proof_capacity_argument+=:native; fi
        if [[ "$proof_first_qc" == true ]]; then
            [[ "$proof_native_session" == true && "$proof_capacity_treatment" == qs2_duckdb && "$proof_stress_samples" != 500 ]] || exit 2
            proof_capacity_argument+=:first-qc
        fi
        if [[ "$proof_reclaim" == true ]]; then
            [[ "$proof_native_session" == true && "$proof_capacity_treatment" == qs2_duckdb ]] || exit 2
            proof_capacity_argument+=:reclaim
        fi
        if [[ "$proof_design_memory" == true ]]; then
            [[ "$proof_stress_samples" == 300 && "$proof_native_session" == true &&
                "$proof_reclaim" == false && "$proof_first_qc" == false ]] || exit 2
            proof_capacity_argument+=:design-memory
        fi
        if [[ "$proof_extended120" == true ]]; then proof_capacity_argument+=:extended120; fi
        proof_operation+=("$proof_capacity_argument")
    fi
    if [[ -n "$proof_capacity_reference" ]]; then
        [[ "$proof_capacity_reference" == "$proof_root/"* ]] || exit 2
        proof_operation+=("$proof_capacity_reference")
    fi
    if [[ "$proof_ruv_pearson" == true && "$2" == verify-da ]]; then
        proof_verify_argument=diagnostic:ruv-pearson
        if [[ "$proof_ruv_stress" == true ]]; then proof_verify_argument+=-"$proof_stress_samples"
        elif [[ "$proof_ruv_small" == true ]]; then proof_verify_argument+=-small; fi
        if [[ "$proof_timed_ruv" == true ]]; then proof_verify_argument+=:timed; fi
        if [[ "$proof_native_session" == true ]]; then proof_verify_argument+=:native; fi
        proof_operation+=("$proof_verify_argument")
    elif [[ "$proof_native_session" == true && "$2" == verify-da ]]; then
        [[ -z "$proof_capacity_reference" ]] || exit 2
        proof_operation+=(native-session)
    fi
    if [[ "$proof_timed_ruv" == true || "$proof_first_qc" == true ]]; then proof_gui_diagnostic=false; fi
    if [[ "$proof_design_memory" == true ]]; then proof_gui_diagnostic=true; fi
    if [[ "$proof_gui_diagnostic" == true || "$proof_reclaim" == true ]]; then proof_process_sampling=true; fi
    proof_command=(Rscript --vanilla "${proof_script_dir}/run_workflow_checkpoint_proof.R" "${proof_operation[@]}")
    if [[ "$proof_sequence_selector" =~ ^(perf108:(B[12]|T[12]|F0|N[12]|I[1-3]|J[12]|K[1-4]|M[12]|Q[CT]|U[BT]|Z[ACOVFSQL]|E[1-8]|R[CS][01]|P[AVO][01]|[YW](S[1-4]|[RPW][1-6])|PC|X(S[1-4]|R[1-6]|F[1-8])|A[SF]1|H[SF]1|C[23]|P[SF][12]|R[SFR][12]|V([1-9]|1[0-6]))|perf101:(D[01]|D0R|V)|perf106:(L[01]|S[1-6]|F[1-8]|R[1-4]|V2-(A|V|I|J|D[12]|Q([1-9]|1[0-9]|2[0-4]))))$ ]]; then
        checkpointSnapshotPerf101Sources "$(realpath -e -- "${proof_script_dir}/../..")" "$proof_run_dir/source"
        if [[ "$proof_sequence_selector" == perf108:* ]]; then
            cp -- "${proof_script_dir}/workflow_perf108_structural.R" "$proof_run_dir/source/tools/profiling/workflow_perf108_structural.R"
            cp -- "${proof_script_dir}/workflow_perf108_ordered.R" "$proof_run_dir/source/tools/profiling/workflow_perf108_ordered.R"
            if [[ "$proof_sequence_selector" == perf108:[YWE]* || "$proof_sequence_selector" == perf108:ZQ ]]; then
                proof_binding_family=structural
                if [[ "$proof_sequence_selector" == perf108:E* || "$proof_sequence_selector" == perf108:ZQ ]]; then proof_binding_family=closeout; fi
                if [[ "$proof_sequence_selector" == perf108:W* ]]; then proof_binding_family=ordered-groups; fi
                mkdir -p "$proof_run_dir/source/docs/qa/workflow-execution/perf-108/$proof_binding_family"
                cp -- "${proof_script_dir}/../../docs/qa/workflow-execution/perf-108/$proof_binding_family/bindings.rds" \
                    "$proof_run_dir/source/docs/qa/workflow-execution/perf-108/$proof_binding_family/bindings.rds"
            fi
            for proof_source_file in R src inst man DESCRIPTION NAMESPACE configure LICENSE COPYING COPYING.LESSER \
                tests/testthat tests/testdata/workflow-checkpoints tools/profiling/perf108-extension tools/profiling/perf108-fastduckdb tools/profiling/workflow_perf106.R tools/profiling/workflow_perf108.R docs/qa/workflow-execution/perf-108/runtime-bindings.rds docs/qa/workflow-execution/perf-108/runtime-repair-bindings.rds docs/qa/workflow-execution/perf-108/runtime-staged-bindings.rds docs/qa/workflow-execution/perf-108/runtime-attribution-bindings.rds docs/qa/workflow-execution/perf-108/runtime-batch-bindings.rds docs/qa/workflow-execution/perf-108/cache-qualification/bindings.rds; do
                if [[ "$proof_source_file" == docs/qa/workflow-execution/perf-108/cache-qualification/bindings.rds && ! -f "${proof_script_dir}/../../$proof_source_file" ]]; then
                    [[ "$proof_sequence_selector" != perf108:X* ]] || exit 2
                    continue
                fi
                if [[ ( "$proof_source_file" == docs/qa/workflow-execution/perf-108/runtime-staged-bindings.rds || "$proof_source_file" == docs/qa/workflow-execution/perf-108/runtime-attribution-bindings.rds || "$proof_source_file" == docs/qa/workflow-execution/perf-108/runtime-batch-bindings.rds ) &&
                    ! -f "${proof_script_dir}/../../$proof_source_file" ]]; then
                    [[ "$proof_sequence_selector" == perf108:M1 || "$proof_sequence_selector" == perf108:M2 ]] || exit 2
                    continue
                fi
                mkdir -p -- "$(dirname -- "$proof_run_dir/source/$proof_source_file")"
                if [[ -d "${proof_script_dir}/../../$proof_source_file" ]]; then
                    mkdir -p -- "$proof_run_dir/source/$proof_source_file"
                    cp -a -- "${proof_script_dir}/../../$proof_source_file/." "$proof_run_dir/source/$proof_source_file/"
                else
                    cp -a -- "${proof_script_dir}/../../$proof_source_file" "$proof_run_dir/source/$proof_source_file"
                fi
            done
            # Compiler products are not source inputs; package installation rebuilds them.
            find "$proof_run_dir/source/src" -maxdepth 1 -type f \( -name '*.o' -o -name '*.so' \) -delete
        fi
        if [[ "$proof_sequence_selector" == perf106:* ]]; then
            for proof_source_file in tools/profiling/workflow_perf106.R \
                tickets/openai/verification/perf106-repair-protocol-v2.json \
                docs/qa/workflow-execution/perf-106/candidate-build.json \
                docs/qa/workflow-execution/perf-106/L0-preserved-primary.rds \
                docs/qa/workflow-execution/perf-106/screen-decision.json \
                docs/qa/workflow-execution/perf-106/confirmation-input-bindings.rds \
                docs/qa/workflow-execution/perf-106/p14-regression-bindings.rds \
                tests/testthat/helper-workflow-checkpoint-fixtures.R \
                tests/testdata/workflow-checkpoints/scenarios.json; do
                mkdir -p -- "$(dirname -- "$proof_run_dir/source/$proof_source_file")"
                cp -- "${proof_script_dir}/../../$proof_source_file" "$proof_run_dir/source/$proof_source_file"
            done
        fi
        proof_command=(Rscript --vanilla "$proof_run_dir/source/tools/profiling/run_workflow_checkpoint_proof.R" "${proof_operation[@]}")
    fi
    if [[ "$proof_utility_diagnostics" == true ]]; then proof_command+=(--diagnostic-memory); fi
    if [[ "$proof_qvalue" == true ]]; then proof_command+=(qvalue); fi
    if [[ "$2" == "install-corrected" ]]; then
        proof_command=(R CMD INSTALL "--library=$proof_run_dir/library" "$proof_run_dir/source")
    fi
    # NUL-separated fields preserve exact argv, including spaces/newlines. This
    # is launch intent, not proof that a service started or remains alive.
    printf '%s\0' dia-controller-request-v1 "${proof_unit}.service" "$2" "$proof_library" \
        "$proof_memory_max" 0 "$proof_runtime_max" "${proof_command[@]}" \
        | tee "$proof_run_dir/controller-request.nul" >/dev/null
    if [[ "$proof_memory_guard" == true || "$2" == "accounting-probe-retained" || "$2" == "browser-da-accounted" || "$2" == "browser-da-default-qc" || "$2" == "browser-da-profile" || "$2" == "browser-baseline" || "$2" == "exclusion-benchmark" || "$2" == "exclusion-memory" || "$2" == "exclusion-allocation" || ( "$proof_capped_sixteen" == true && "$2" == "verify-da" ) ]]; then
        proof_poll_limit=150
        proof_poll_delay=0.1
        if [[ "$2" == "exclusion-benchmark" || "$2" == "exclusion-memory" ]]; then
            proof_poll_limit=930
            proof_poll_delay=1
        fi
        if [[ "$2" == "utility-member" || "$2" == "browser-da-accounted" || "$2" == "browser-da-default-qc" || "$2" == "browser-da-profile" || "$2" == "browser-baseline" ]]; then
            proof_poll_limit=930
            proof_poll_delay=1
        fi
        if [[ "$2" == "exclusion-allocation" || "$proof_capped_sixteen" == true || "$proof_memory_guard" == true ]]; then
            if [[ "$proof_capped_sixteen" == true || "$proof_memory_guard" == true ]]; then
                proof_poll_limit=$((proof_runtime_max * 10 + 300))
            else
                proof_poll_limit=1500
            fi
            proof_poll_delay=0.1
            if [[ "$2" == "exclusion-allocation" ]]; then
                checkpointDiagnosticGuard "$(awk '/^MemAvailable:/ {printf "%.0f", $2 * 1024}' /proc/meminfo)" 0
            fi
            # Stop our service on controller error/interruption as well as guard failure.
            trap 'systemctl --user stop "${proof_unit}.service" >/dev/null 2>&1 || true' EXIT
            trap 'exit 130' INT TERM
        fi
        proof_service_env=(--setenv=NOT_CRAN=true --setenv="CHROMOTE_CHROME=$proof_browser")
        if [[ "$proof_sequence_selector" =~ ^(perf108:(B[12]|T[12]|F0|N[12]|I[1-3]|J[12]|K[1-4]|M[12]|Q[CT]|U[BT]|Z[ACOVFSQL]|E[1-8]|R[CS][01]|P[AVO][01]|[YW](S[1-4]|[RPW][1-6])|PC|X(S[1-4]|R[1-6]|F[1-8])|A[SF]1|H[SF]1|C[23]|P[SF][12]|R[SFR][12]|V([1-9]|1[0-6]))|perf101:(D[01]|D0R|V)|perf106:(L[01]|S[1-6]|F[1-8]|R[1-4]|V2-(A|V|I|J|D[12]|Q([1-9]|1[0-9]|2[0-4]))))$ ]]; then
            proof_service_env+=(--setenv=OMP_NUM_THREADS=1 --setenv=OPENBLAS_NUM_THREADS=1
                --setenv=MKL_NUM_THREADS=1 --setenv=VECLIB_MAXIMUM_THREADS=1)
        fi
        if [[ "$proof_gui_diagnostic" == true || "$proof_ruv_pearson" == true ]]; then
            proof_service_env+=(--setenv="MULTISCHOLAR_DIA_DIAGNOSTIC_SECONDS=$proof_runtime_max")
        fi
        if [[ "$proof_ruv_pearson" == true && "$proof_ruv_small" != true && "$proof_ruv_stress" != true && "$2" == browser-baseline ]]; then
            proof_service_env+=(--setenv=MULTISCHOLAR_DIA_SMALL_GUI_CAP=19327352832)
        elif [[ "$proof_ruv_pearson" == true && "$proof_ruv_small" != true && "$proof_ruv_stress" != true && "$2" == verify-da ]]; then
            proof_service_env+=(--setenv=MULTISCHOLAR_DIA_SMALL_VERIFY_CAP=19327352832)
        elif [[ "$proof_capped_sixteen" == true && "$2" == "browser-baseline" ]]; then
            proof_service_env+=(--setenv=MULTISCHOLAR_DIA_SMALL_GUI_CAP=17179869184)
        elif [[ "$proof_capped_sixteen" == true && "$2" == "verify-da" ]]; then
            proof_service_env+=(--setenv=MULTISCHOLAR_DIA_SMALL_VERIFY_CAP=17179869184)
        fi
        systemd-run --user --unit="$proof_unit" --remain-after-exit --json=short \
            --property=ExitType=cgroup \
            --property="MemoryMax=$proof_memory_max" --property=MemorySwapMax=0 --property="RuntimeMaxSec=$proof_runtime_max" \
            --property=MemoryAccounting=yes --property=IOAccounting=yes \
            "${proof_service_env[@]}" \
            --property="StandardOutput=append:$proof_run_dir/service.log" \
            --property=StandardError=inherit --working-directory="$proof_run_dir" \
            "${proof_command[@]}" | tee "$proof_run_dir/launch.json"
        if [[ "$proof_utility_diagnostics" == true || "$proof_sequence_sampling" == true ]]; then
            printf 'timestamp\tcurrent_bytes\tpeak_bytes\tanon_bytes\tfile_bytes\tkernel_bytes\tdirty_bytes\twriteback_bytes\n' > "$proof_run_dir/memory-samples.tsv"
        fi
        if [[ "$proof_process_sampling" == true ]]; then
            printf 'timestamp\tpid\tppid\tname\trss_bytes\tanon_bytes\tfile_bytes\tshmem_bytes\tpss_bytes\tprivate_clean_bytes\tprivate_dirty_bytes\tthreads\n' > "$proof_run_dir/process-samples.tsv"
        fi
        proof_finished=false
        for ((proof_poll=0; proof_poll<proof_poll_limit; proof_poll++)); do
            proof_state="$(systemctl --user show "${proof_unit}.service" \
                --property=Id,InvocationID,ActiveState,SubState,ExitType,StateChangeTimestampMonotonic,MainPID,TasksCurrent,ControlGroup,Result,ExecMainCode,ExecMainStatus,ExecMainStartTimestampMonotonic,ExecMainExitTimestampMonotonic,MemoryCurrent,MemoryPeak,MemorySwapPeak,CPUUsageNSec,IOReadBytes,IOWriteBytes)"
            if checkpointServiceTerminal "$proof_state"; then
                printf '%s\n' "$proof_state" > "$proof_run_dir/terminal-properties.txt"
                if [[ "$proof_memory_guard" == true ]]; then
                    proof_peak="$(sed -n 's/^MemoryPeak=//p' <<< "$proof_state")"
                    if [[ ! "$proof_peak" =~ ^[0-9]+$ ]] || (( proof_peak >= proof_stop_peak )); then
                        printf 'reason=terminal_peak_guard\n%s\n' "$proof_state" | tee "$proof_run_dir/guard-stop.txt"
                        exit 76
                    fi
                fi
                proof_finished=true
                break
            fi
            if [[ "$proof_capped_sixteen" == true || "$proof_memory_guard" == true ]]; then
                proof_available="$(awk '/^MemAvailable:/ {printf "%.0f", $2 * 1024}' /proc/meminfo)"
                proof_current="$(sed -n 's/^MemoryCurrent=//p' <<< "$proof_state")"
                proof_peak="$(sed -n 's/^MemoryPeak=//p' <<< "$proof_state")"
                if [[ ! "$proof_available" =~ ^[0-9]+$ ||
                    ! "$proof_current" =~ ^[0-9]+$ || ! "$proof_peak" =~ ^[0-9]+$ ]]; then
                    printf 'reason=missing_memory_observation\n%s\n' "$proof_state" \
                        | tee "$proof_run_dir/guard-stop.txt"
                    exit 76
                fi
                if (( proof_available < proof_min_host_available ||
                    proof_current >= proof_stop_peak || proof_peak >= proof_stop_peak )); then
                    printf 'reason=early_memory_stop\ncurrent_bytes=%s\npeak_bytes=%s\nmem_available_bytes=%s\n%s\n' \
                        "$proof_current" "$proof_peak" "$proof_available" "$proof_state" \
                        | tee "$proof_run_dir/guard-stop.txt"
                    exit 76
                fi
            elif [[ "$2" == "exclusion-allocation" ]]; then
                proof_available="$(awk '/^MemAvailable:/ {printf "%.0f", $2 * 1024}' /proc/meminfo)"
                proof_current="$(sed -n 's/^MemoryCurrent=//p' <<< "$proof_state")"
                if ! proof_reason="$(checkpointDiagnosticGuard "$proof_available" "$proof_current")"; then
                    printf '%s\n' "$proof_state" > "$proof_run_dir/guard-properties.txt"
                    printf 'reason=%s\navailable_bytes=%s\ncurrent_bytes=%s\n' "$proof_reason" "$proof_available" "$proof_current" > "$proof_run_dir/guard-stop.txt"
                    systemctl --user stop "${proof_unit}.service"
                    exit 76
                fi
            fi
            if [[ "$proof_utility_diagnostics" == true || "$proof_sequence_sampling" == true ]]; then
                proof_cgroup="$(sed -n 's/^ControlGroup=//p' <<< "$proof_state")"
                [[ "$proof_cgroup" == /* && "$proof_cgroup" != *..* ]] || exit 74
                if [[ "$proof_process_sampling" == true ]] && (( proof_poll % 10 == 0 )); then
                    checkpointProcessSample "/sys/fs/cgroup${proof_cgroup}" "$EPOCHREALTIME" "$proof_run_dir/process-samples.tsv"
                fi
                # Sample existing counters; no extra worker, GC, payload reads or RSS sums.
                if proof_sample="$(checkpointMemorySample "/sys/fs/cgroup${proof_cgroup}/memory.stat" \
                    "$EPOCHREALTIME" "$proof_current" "$proof_peak" 2>> "$proof_run_dir/memory-sample-errors.log")"; then
                    printf '%s\n' "$proof_sample" >> "$proof_run_dir/memory-samples.tsv"
                else
                    # Exit can remove memory.stat after the preceding systemctl poll.
                    # Re-observe exact terminal accounting; never substitute last peak.
                    proof_terminal="$(systemctl --user show "${proof_unit}.service" \
                        --property=Id,InvocationID,ActiveState,SubState,ExitType,StateChangeTimestampMonotonic,MainPID,TasksCurrent,ControlGroup,Result,ExecMainCode,ExecMainStatus,ExecMainStartTimestampMonotonic,ExecMainExitTimestampMonotonic,MemoryCurrent,MemoryPeak,MemorySwapPeak,CPUUsageNSec,IOReadBytes,IOWriteBytes)"
                    if ! checkpointServiceTerminal "$proof_terminal"; then
                        printf 'reason=missing_live_sample\n%s\n' "$proof_terminal" > "$proof_run_dir/guard-stop.txt"
                        exit 74
                    fi
                    printf '%s\tterminal_cgroup_disappeared\n' "$EPOCHREALTIME" >> "$proof_run_dir/memory-sample-gaps.tsv"
                    printf '%s\n' "$proof_terminal" > "$proof_run_dir/terminal-properties.txt"
                    # The normal terminal branch applies the unchanged peak guard.
                    continue
                fi
            fi
            sleep "$proof_poll_delay"
        done
        printf '%s\n' "$proof_state" | tee "$proof_run_dir/terminal-properties.txt"
        if [[ "$proof_finished" != true ]]; then
            if [[ "$proof_capped_sixteen" == true || "$proof_memory_guard" == true ]]; then
                printf 'reason=observation_limit\n%s\n' "$proof_state" \
                    | tee "$proof_run_dir/guard-stop.txt"
                exit 75
            fi
            echo "Probe observation timed out; owned unit is retained for inspection: ${proof_unit}.service" >&2
            exit 75
        fi
        # The main process is terminal. Stop only this request-owned retained
        # probe; never interpret an observation timeout as permission to restart.
        set +e
        Rscript --vanilla "${proof_script_dir}/run_workflow_checkpoint_proof.R" \
            --internal-retained-interval "$proof_run_dir" "${proof_unit}.service"
        proof_interval_status=$?
        systemctl --user stop "${proof_unit}.service"
        proof_stop_status=$?
        if [[ ( "$2" == "exclusion-allocation" || "$proof_capped_sixteen" == true || "$proof_memory_guard" == true ) && "$proof_stop_status" == 0 ]]; then
            trap - EXIT INT TERM
        fi
        set -e
        printf 'unit=%s\ninterval_exit=%s\nstop_exit=%s\n' \
            "$proof_unit" "$proof_interval_status" "$proof_stop_status" \
            | tee "$proof_run_dir/controller-status.txt"
        if [[ "$proof_interval_status" != 0 || "$proof_stop_status" != 0 ]]; then exit 74; fi
        if [[ "$proof_state" != *"Result=success"* || "$proof_state" != *"ExecMainStatus=0"* ]]; then
            exit 1
        fi
        exit 0
    fi
    if [[ "$2" == "accounting-probe" || "$2" == "accounting-probe-failure" ]]; then
        # Probe only: service-owned output avoids an outer tee handling workload
        # logs. Keep exact terminal JSON separate from those logs and exit status.
        set +e
        systemd-run --user --unit="$proof_unit" --wait --json=short \
            --property=MemoryMax=8G --property=MemorySwapMax=0 --property=RuntimeMaxSec=900 \
            --property=MemoryAccounting=yes --property=IOAccounting=yes \
            --property="StandardOutput=append:$proof_run_dir/service.log" \
            --property=StandardError=inherit --working-directory="$proof_run_dir" \
            "${proof_command[@]}" | tee "$proof_run_dir/terminal-accounting.json"
        proof_pipeline_status=("${PIPESTATUS[@]}")
        journalctl --user --unit="${proof_unit}.service" --output=json --no-pager \
            --output-fields=MESSAGE_ID,USER_UNIT,USER_INVOCATION_ID,CPU_USAGE_NSEC,MEMORY_PEAK,MEMORY_SWAP_PEAK,IO_METRIC_READ_BYTES,IO_METRIC_WRITE_BYTES \
            | tee "$proof_run_dir/terminal-journal.jsonl" >/dev/null
        proof_journal_status=("${PIPESTATUS[@]}")
        Rscript --vanilla "${proof_script_dir}/run_workflow_checkpoint_proof.R" \
            --internal-accounting-journal "$proof_run_dir" "${proof_unit}.service"
        proof_metrics_status=$?
        set -e
        printf 'unit=%s\noperation=%s\nservice_exit=%s\nlog_exit=%s\njournal_exit=%s\njournal_log_exit=%s\nmetrics_exit=%s\n' \
            "$proof_unit" "$2" "${proof_pipeline_status[0]}" "${proof_pipeline_status[1]}" \
            "${proof_journal_status[0]}" "${proof_journal_status[1]}" "$proof_metrics_status" \
            | tee "$proof_run_dir/controller-status.txt"
        if [[ "${proof_pipeline_status[0]}" != 0 ]]; then exit "${proof_pipeline_status[0]}"; fi
        if [[ "${proof_pipeline_status[1]}" != 0 ]]; then exit 74; fi
        if [[ "${proof_journal_status[0]}" != 0 || "${proof_journal_status[1]}" != 0 || "$proof_metrics_status" != 0 ]]; then exit 74; fi
        exit 0
    fi
    # Keep both service failure and log-publication failure visible. The parent
    # retains the serial lease until terminal evidence has been written.
    set +e
    systemd-run --user --unit="$proof_unit" --wait --pipe \
        --property=MemoryMax=8G --property=MemorySwapMax=0 --property=RuntimeMaxSec=900 \
        --property=MemoryAccounting=yes --property=IOAccounting=yes \
        --setenv=NOT_CRAN=true --setenv="CHROMOTE_CHROME=$proof_browser" --working-directory="$proof_run_dir" \
        "${proof_command[@]}" 2>&1 | tee "$proof_run_dir/service.log"
    proof_pipeline_status=("${PIPESTATUS[@]}")
    set -e
    printf 'unit=%s\noperation=%s\nservice_exit=%s\nlog_exit=%s\n' \
        "$proof_unit" "$2" "${proof_pipeline_status[0]}" "${proof_pipeline_status[1]}" \
        | tee "$proof_run_dir/controller-status.txt"
    if [[ "${proof_pipeline_status[0]}" != 0 ]]; then
        exit "${proof_pipeline_status[0]}"
    fi
    if [[ "${proof_pipeline_status[1]}" != 0 ]]; then
        exit 74
    fi
    exit 0
fi
exec Rscript --vanilla "${proof_script_dir}/run_workflow_checkpoint_proof.R" "$@"
