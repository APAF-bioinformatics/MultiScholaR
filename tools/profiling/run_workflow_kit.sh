#!/usr/bin/env bash
# Guarded launcher and 50 ms memory sampler for KIT stages (KIT-001).
# Usage:
#   run_workflow_kit.sh validate <stage> <library>
#   run_workflow_kit.sh run <stage> <library>
# Stages: kit-lifecycle-(complete|fail) (journey shape), kit-lifecycle-operation-(complete|fail)
# (operation shape, KIT-005), kit-journey (baseline pair + production default), kit-0NN (operation members).
#   run_workflow_kit.sh --sample <cgroup> <trace> <deadline_uptime> <slice> <stop_file> <guard_file>
# The proof root comes from docs/qa/workflow-execution/kit/bindings.json (harness.root).
set -euo pipefail
kit_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd -- "$kit_dir/../.." && pwd)"

if [[ "${1:-}" == "--sample" ]]; then
    [[ "$#" == 7 ]] || exit 2
    group="$2"; trace="$3"; deadline="$4"; slice="$5"; stop_file="$6"; guard_file="$7"
    printf 'uptime\tcurrent\tpeak\tanon\tfile\tkernel\tswap_peak\n' > "$trace"
    while [[ ! -e "$stop_file" ]]; do
        if [[ ! -d "$group" ]]; then
            printf '%s\tcgroup_gone\n' "$(cut -d' ' -f1 /proc/uptime)" >> "${trace%.tsv}-gaps.tsv"
            sleep 0.05; continue
        fi
        # One awk per sample (KIT-005 correction 6): the earlier six process spawns per sample left gaps
        # above the 250 ms validity rule on a busy host.
        read -r uptime current peak swap anon file kernel available < <(awk '
            FILENAME=="/proc/uptime"{u=$1} FILENAME==ARGV[2]{c=$1} FILENAME==ARGV[3]{p=$1} FILENAME==ARGV[4]{s=$1}
            FILENAME==ARGV[5]{if($1=="anon")a=$2; if($1=="file")f=$2; if($1=="kernel")k=$2}
            FILENAME=="/proc/meminfo" && /^MemAvailable:/{m=$2*1024}
            END{printf "%s %s %s %s %s %s %s %.0f\n", u, (c==""?"NA":c), (p==""?"NA":p), (s==""?"NA":s), (a==""?"NA":a), (f==""?"NA":f), (k==""?"NA":k), m}' \
            /proc/uptime "$group/memory.current" "$group/memory.peak" "$group/memory.swap.peak" "$group/memory.stat" /proc/meminfo 2>/dev/null || echo "NA NA NA NA NA NA NA NA")
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$uptime" "$current" "$peak" "$anon" "$file" "$kernel" "$swap" >> "$trace"
        reason=""
        if [[ "$peak" =~ ^[0-9]+$ ]] && (( peak >= 16642998272 )); then reason="memory_guard"; fi
        if [[ "$available" =~ ^[0-9]+$ ]] && (( available < 2147483648 )); then reason="host_guard"; fi
        if awk -v u="$uptime" -v d="$deadline" 'BEGIN{exit !(u >= d)}'; then reason="time_guard"; fi
        if [[ -n "$reason" ]]; then
            printf '%s\n' "$reason" > "$guard_file"
            systemctl --user stop "$slice" || true
            exit 0
        fi
        sleep 0.05
    done
    exit 0
fi

# KIT-025: the KIT install. The CORE-pinned runner installs the committed tree as before (it is not edited); the
# installed package is then stamped with the commit it was built from (`Config/MultiScholaR/SourceCommit`, in the
# DESCRIPTION file and in Meta/package.rds, which packageDescription reads), so the record's software identity
# names the source commit. Usage: run_workflow_kit.sh install <previous library> <core dir>
if [[ "${1:-}" == "install" ]]; then
    [[ "$#" == 3 && -d "$2" && -d "$3" ]] || { echo "Usage: run_workflow_kit.sh install <previous library> <core dir>" >&2; exit 2; }
    before="$(ls -d "$3"/browser-smoke.* 2>/dev/null | sort || true)"
    bash "$kit_dir/run_workflow_checkpoint_proof.sh" --service install-corrected "$2" "$3"
    run="$(comm -13 <(printf '%s\n' "$before") <(ls -d "$3"/browser-smoke.* | sort) | head -1)"
    [[ -n "$run" && -f "$run/source-commit.txt" && -d "$run/library/MultiScholaR" ]] || { echo "install run directory not found" >&2; exit 1; }
    sha="$(tr -d '[:space:]' < "$run/source-commit.txt")"
    [[ "$sha" =~ ^[0-9a-f]{40}$ ]] || { echo "source commit is not a 40-character sha" >&2; exit 1; }
    Rscript --vanilla -e '
        args <- commandArgs(TRUE); lib <- args[[1L]]; sha <- args[[2L]]
        meta <- file.path(lib, "Meta", "package.rds"); pkg <- readRDS(meta)
        pkg$DESCRIPTION["Config/MultiScholaR/SourceCommit"] <- sha; saveRDS(pkg, meta)
        cat("Config/MultiScholaR/SourceCommit: ", sha, "\n", sep = "", file = file.path(lib, "DESCRIPTION"), append = TRUE)
        stopifnot(identical(unname(read.dcf(file.path(lib, "DESCRIPTION"), fields = "Config/MultiScholaR/SourceCommit")[1L, 1L]), sha))
    ' "$run/library/MultiScholaR" "$sha"
    # KIT-047: the install writes a thin overlay (the package and the exact-sum artifacts); a package the previous
    # kit library carries beyond those (the DPC engines limpa and fastlimpadpc, KIT-047) is carried into the new
    # overlay unchanged, so every later build's library resolves the same engines.
    for extra in "$2"/*/; do
        name="$(basename "$extra")"
        [[ -d "$run/library/$name" ]] || cp -a "$extra" "$run/library/$name"
    done
    echo "INSTALLED $run $sha"
    exit 0
fi

mode="${1:-}"; stage="${2:-}"; library="${3:-}"; rerun="${4:-}"
[[ ( "$mode" == validate || "$mode" == run || ( "$mode" == rerun && "$rerun" =~ ^[a-z0-9-]+(,[a-z0-9-]+)*$ ) ) && "$stage" =~ ^kit-(lifecycle-(operation-)?(complete|fail)|journey(-0[0-9]{2}[a-z]?)?|gui-0[0-9]{2}|0[0-9]{2})$ && -d "$library" ]] || {
    echo "Usage: run_workflow_kit.sh validate|run <stage> <library> | rerun <stage> <library> <member,member>" >&2; exit 2; }
root="$(jq -r '.harness.root' "$repo/docs/qa/workflow-execution/kit/bindings.json")"
# The browser driver runs in the coordinator (outside the slice); chromote needs the executable (as the pinned runner does).
export CHROMOTE_CHROME="${CHROMOTE_CHROME:-$(command -v chromium || command -v google-chrome || command -v brave || true)}"
[[ -x "$CHROMOTE_CHROME" ]] || { echo "refused: no Chromium-based browser executable found" >&2; exit 2; }
[[ "$root" == /* ]] || { echo "bindings carry no proof root" >&2; exit 2; }
mkdir -p -- "$root"
cd -- "$repo"
if [[ "$mode" == run || "$mode" == rerun ]]; then
    # Preconditions (plan): clean tree, bindings committed, 12 GiB free, no other managed unit, no live stage.
    [[ -z "$(git status --porcelain)" ]] || { echo "refused: working tree is not clean" >&2; exit 76; }
    available="$(awk '/^MemAvailable:/ {printf "%.0f", $2 * 1024}' /proc/meminfo)"
    (( available >= 12884901888 )) || { echo "refused: host headroom is below 12 GiB" >&2; exit 76; }
    units="$(systemctl --user list-units --state=active --no-legend --plain 'multischolar*' | sed '/^\s*$/d' || true)"
    [[ -z "$units" ]] || { echo "refused: another managed unit is active: $units" >&2; exit 76; }
    if [[ "$mode" == run ]]; then
        [[ ! -e "$root/$stage" ]] || { echo "refused: stage directory exists: $root/$stage" >&2; exit 76; }
        label="$stage"; extra=()
    else
        # A rerun of named members inside an existing, finished stage (KIT-006): its own start file and
        # log, so the first launch's records are never overwritten.
        [[ -f "$root/$stage/stage.rds" ]] || { echo "refused: no finished stage at $root/$stage" >&2; exit 76; }
        n=1; while [[ -e "$root/../$stage.rerun$n.log" ]]; do n=$((n + 1)); done
        label="$stage.rerun$n"; extra=("--rerun=$rerun")
    fi
    { date; echo "MemAvailable $(awk -v a="$available" 'BEGIN{printf "%.2f GiB", a/1073741824}')"; git rev-parse HEAD; echo "$mode $rerun"; } \
        > "$root/../$label-start.txt"
    nohup Rscript --vanilla "$kit_dir/workflow_kit_coordinator.R" "$stage" "$library" "$root" "${extra[@]}" \
        > "$root/../$label.log" 2>&1 &
    echo "launched $label (pid $!); log $root/../$label.log"
    exit 0
fi
exec Rscript --vanilla "$kit_dir/workflow_kit_coordinator.R" "$stage" "$library" "$root" --validate-only
