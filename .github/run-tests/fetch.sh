#!/bin/bash
# run-tests: wait for a run, download its `results` artifact and print the tables (bash and gh only).
#
# Usage: fetch.sh [<run id> | <request name>] [<out dir>]
#   <run id>        default: the latest run-tests run on the branch run-tests
#   <request name>  the latest run titled "run-tests <name>" (the request file requests/<name>.json); it may take a
#                   minute after the push before the dispatched run exists
#   <out dir>  default: ./run-tests-<run id> (must not exist)
set -euo pipefail
repo="${RUN_TESTS_REPO:-APAF-bioinformatics/MultiScholaR}"
run="${1:-}"
if [ -n "$run" ] && ! [[ "$run" =~ ^[0-9]+$ ]]; then
    name="$run"
    run="$(gh run list -R "$repo" -w run-tests.yml -b run-tests -L 100 --json databaseId,displayTitle \
        -q "[.[] | select(.displayTitle == \"run-tests $name\")][0].databaseId // empty")"
    [ -n "$run" ] || { echo "no run titled 'run-tests $name' (yet)" >&2; exit 1; }
elif [ -z "$run" ]; then
    run="$(gh run list -R "$repo" -w run-tests.yml -b run-tests -L 1 --json databaseId -q '.[0].databaseId')"
    [ -n "$run" ] || { echo "no run-tests run found" >&2; exit 1; }
fi
out="${2:-./run-tests-$run}"
if [ -e "$out" ]; then echo "$out exists: not overwritten" >&2; exit 1; fi
echo "run $run: https://github.com/$repo/actions/runs/$run"
while [ "$(gh run view "$run" -R "$repo" --json status -q .status)" != "completed" ]; do
    gh run view "$run" -R "$repo" --json jobs -q '.jobs[] | "  \(.name): \(.status) \(.conclusion)"'
    sleep 60
done
conclusion="$(gh run view "$run" -R "$repo" --json conclusion -q .conclusion)"
echo "run $run: $conclusion"
gh run view "$run" -R "$repo" --json jobs -q '.jobs[] | "  \(.name): \(.conclusion)"'
mkdir -p "$out"
if ! gh run download "$run" -R "$repo" -n results -D "$out"; then echo "no results artifact" >&2; exit 1; fi
echo; echo "== meta.txt (head)"; sed -n '1,12p' "$out/meta.txt"
echo; echo "== tests.tsv"; column -t -s $'\t' "$out/tests.tsv"
if [ -f "$out/check-summary.txt" ]; then echo; echo "== check-summary.txt (head)"; sed -n '1,40p' "$out/check-summary.txt"; fi
if [ -f "$out/scripts/scripts.tsv" ]; then echo; echo "== scripts"; column -t -s $'\t' "$out/scripts/scripts.tsv"; fi
if [ -f "$out/withheld.txt" ]; then echo; echo "== withheld (not uploaded)"; cat "$out/withheld.txt"; fi
echo; echo "== failures.txt (head; full file: $out/failures.txt)"; sed -n '1,60p' "$out/failures.txt"
