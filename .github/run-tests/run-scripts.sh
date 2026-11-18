#!/bin/bash
# run-tests: the params' scripts, one Rscript process each, from the package tree's root, one after another.
#
# Usage: run-scripts.sh <params.json> <package dir> <out dir>
# Script i (from 1): `Rscript <path> <args...>` with its env, under `timeout` (timeout_minutes, default 120); its
# console in <out dir>/<i>/console.log and its declared output directory (`out`, relative to the package root)
# copied to <out dir>/<i>/out/ (only the files whose path under it matches `out_pattern`, an extended regex, when given). <out dir>/scripts.tsv: index, path, exit, seconds, status (ok, failed, timeout).
set -uo pipefail
params="$1"; pkg="$(cd "$2" && pwd)"; out="$3"
mkdir -p "$out"; out="$(cd "$out" && pwd)"
printf 'index\tpath\texit\tseconds\tstatus\n' > "$out/scripts.tsv"
n="$(jq '.scripts // [] | length' "$params")"
for ((i = 0; i < n; i++)); do
    index=$((i + 1))
    script="$(jq -r ".scripts[$i]" "$params")"
    path="$(jq -r '.path' <<< "$script")"
    dir_out="$(jq -r '.out // ""' <<< "$script")"
    pattern="$(jq -r '.out_pattern // ""' <<< "$script")"
    minutes="$(jq -r '.timeout_minutes // 120' <<< "$script")"
    mapfile -t args < <(jq -r '.args // [] | .[]' <<< "$script")
    mapfile -t envs < <(jq -r '.env // {} | to_entries[] | "\(.key)=\(.value)"' <<< "$script")
    mkdir -p "$out/$index"
    echo "script $index: $path ${args[*]} (env: ${#envs[@]} variables; limit $minutes minutes)"
    start=$(date +%s)
    (cd "$pkg" && env "${envs[@]}" timeout --kill-after=60 "${minutes}m" Rscript "$path" "${args[@]}") > "$out/$index/console.log" 2>&1
    status=$?
    seconds=$(( $(date +%s) - start ))
    case "$status" in 0) verdict=ok ;; 124|137) verdict=timeout ;; *) verdict=failed ;; esac
    tail -n 20 "$out/$index/console.log"
    echo "script $index: exit $status after $seconds s ($verdict)"
    if [ -n "$dir_out" ] && [ -d "$pkg/$dir_out" ]; then
        cp -r "$pkg/$dir_out" "$out/$index/out"
        if [ -n "$pattern" ]; then
            (cd "$out/$index/out" && find . -type f | sed 's#^\./##' | grep -Ev -- "$pattern" | while IFS= read -r f; do rm -f -- "$f"; done)
            echo "kept under out/ (out_pattern $pattern): $(find "$out/$index/out" -type f | wc -l) files"
        fi
    elif [ -n "$dir_out" ]; then
        echo "declared output directory $dir_out was not written" > "$out/$index/out-missing.txt"
    fi
    printf '%s\t%s\t%s\t%s\t%s\n' "$index" "$path" "$status" "$seconds" "$verdict" >> "$out/scripts.tsv"
done
cat "$out/scripts.tsv"
