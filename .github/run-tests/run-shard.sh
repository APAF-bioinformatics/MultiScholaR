#!/bin/bash
# run-tests: a shard's test files, one R process per file, one after another.
#
# Usage: run-shard.sh <package dir> <out dir> <per-file timeout, minutes> <job deadline, epoch seconds> <file>...
# Per file: <out dir>/logs/<file>.log (the process's console) and the row run-file.R writes; a file that does not
# exist, a process that exits without a row, or one killed by the per-file timeout gets its row here (status crashed
# or timeout). A file that would start after the job deadline is not started and is recorded as timeout. The shard's
# table is <out dir>/tests.tsv. Exit status 0 unless the arguments are wrong: the outcome is in the table.
# Environment: SIDE=base records a file the tree lacks as absent (head: crashed); FILE_TIMEOUTS, a JSON object
# {file: minutes}, overrides the per-file timeout for the files it names.
set -uo pipefail
if [ "$#" -lt 5 ]; then echo "usage: run-shard.sh <package dir> <out dir> <timeout min> <deadline epoch> <file>..." >&2; exit 2; fi
pkg="$(cd "$1" && pwd)"; out="$2"; per_file="$3"; deadline="$4"; shift 4
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$out/rows" "$out/failures" "$out/logs"
out="$(cd "$out" && pwd)"

row() { # file status seconds
    printf 'file\tpassed\tfailed\terrors\tskipped\twarnings\tseconds\tstatus\n%s\t0\t0\t0\t0\t0\t%s\t%s\n' "$1" "$3" "$2" \
        > "$out/rows/${1%.R}.tsv"
}

for file in "$@"; do
    stem="${file%.R}"
    path="$pkg/tests/testthat/$file"
    log="$out/logs/$stem.log"
    if [ ! -f "$path" ] && [ "${SIDE:-head}" = base ]; then
        echo "no such test file in the base tree: tests/testthat/$file" > "$log"
        row "$file" absent 0
        continue
    fi
    if [ ! -f "$path" ]; then
        echo "no such test file in the package tree: tests/testthat/$file" > "$log"
        row "$file" crashed 0
        printf 'file: %s\ntest: <none>\nkind: crashed\nmessage:\n  no such test file in the package tree\n\n' "$file" > "$out/failures/$stem.txt"
        echo "::warning::$file: no such test file (recorded as crashed)"
        continue
    fi
    left=$(( (deadline - $(date +%s)) / 60 ))
    if [ "$left" -lt 2 ]; then
        echo "not started: the job's time budget is spent" > "$log"
        row "$file" timeout 0
        echo "::warning::$file: not started, the job's time budget is spent (recorded as timeout)"
        continue
    fi
    limit="$(jq -r --arg f "$file" --argjson d "$per_file" '.[$f] // $d' <<< "${FILE_TIMEOUTS:-{\}}")"; [ "$left" -lt "$limit" ] && limit="$left"
    start=$(date +%s)
    timeout --kill-after=60 "${limit}m" Rscript --vanilla "$here/run-file.R" "$pkg" "$path" "$out" > "$log" 2>&1
    status=$?
    seconds=$(( $(date +%s) - start ))
    if [ "$status" -eq 124 ] || [ "$status" -eq 137 ]; then
        rm -f "$out/rows/$stem.tsv"
        row "$file" timeout "$seconds"
        echo "::warning::$file: killed after $limit minutes (recorded as timeout)"
    elif [ ! -f "$out/rows/$stem.tsv" ] || [ "$status" -ne 0 ]; then
        rm -f "$out/rows/$stem.tsv"
        row "$file" crashed "$seconds"
        { printf 'file: %s\ntest: <process>\nkind: crashed\nmessage:\n  R exited with status %s and no result row; last lines of its log:\n' "$file" "$status"
          tail -n 15 "$log" | sed 's/^/  /'; echo; } >> "$out/failures/$stem.txt"
        echo "::warning::$file: R exited with status $status and no result row (recorded as crashed)"
    fi
    tail -n 1 "$log"
done

{
    printf 'file\tpassed\tfailed\terrors\tskipped\twarnings\tseconds\tstatus\n'
    for file in "$@"; do tail -n +2 "$out/rows/${file%.R}.tsv"; done
} > "$out/tests.tsv"
cat "$out/tests.tsv"
