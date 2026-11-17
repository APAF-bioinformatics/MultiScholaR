#!/bin/bash
# KIT-124 on hosted runners: download a run's per-cell artifacts and lay them out as aggregate.R reads them.
#
# Usage: assemble.sh <run id> <label> <out root>
#
# Downloads the artifacts k124-<label>-<sim>-<cell> of the run (gh run download) and writes, per simulator that has
# any, <out root>/<sim>-<label>/ :
#   the top-level tables (s<i>-runs.tsv and the others), the cells' tables one after another under one header;
#   s<i>-runs.meta.txt and s1-runs.log, the cells' files one after another;
#   cells/     the per-cell restart files (S2 <cell>.rds, S3 s3-<cell>.rds, S4 <cell>.tsv, .notes.tsv, .timing.tsv);
#   runner/    per cell: the library check, the simulator's output, the job's facts.
# aggregate.R reads cells/ where it exists (S2, S3, S4) and the top-level table and log for S1:
#   Rscript docs/qa/workflow-execution/kit/kit-124/aggregate.R <scored dir> S1=<out root>/s1-<label> S2=...
# An existing <out root>/<sim>-<label> stops the script: nothing is overwritten or merged into.
set -euo pipefail
if [ "$#" -ne 3 ]; then echo "usage: assemble.sh <run id> <label> <out root>" >&2; exit 2; fi
run="$1"; label="$2"; root="$3"
repo="${K124_REPO:-APAF-bioinformatics/MultiScholaR}"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
gh run download "$run" -R "$repo" -D "$tmp" -p "k124-$label-s*"
for sim in s1 s2 s3 s4; do
    mapfile -t parts < <(find "$tmp" -maxdepth 1 -type d -name "k124-$label-$sim-*" | sort -V)
    [ "${#parts[@]}" -gt 0 ] || continue
    out="$root/$sim-$label"
    if [ -e "$out" ]; then echo "$out exists: not overwritten" >&2; exit 1; fi
    mkdir -p "$out/runner"
    for part in "${parts[@]}"; do
        cp "$part"/runner/* "$out/runner/"
        if [ -d "$part/cells" ]; then mkdir -p "$out/cells"; cp -n "$part"/cells/* "$out/cells/"; fi
        for file in "$part"/*; do
            [ -f "$file" ] || continue
            name="$(basename "$file")"
            case "$name" in
                *.tsv) if [ -e "$out/$name" ]; then tail -n +2 "$file" >> "$out/$name"; else cp "$file" "$out/$name"; fi ;;
                *) cat "$file" >> "$out/$name" ;;
            esac
        done
    done
    echo "$out: ${#parts[@]} cells ($(for part in "${parts[@]}"; do printf '%s ' "${part##*-}"; done))"
done
