#!/bin/bash
# KIT-124 confirmation round: build the gate's reference tables from the source artifacts (reference-spec.tsv).
#
# Usage: reference-build.sh <artifacts root> <out dir> [<sources list out>]
# <artifacts root> holds each source run's artifacts as <run id>/<artifact name>/ (as `gh run download <run> -D
# <root>/<run>` or actions/download-artifact with run-id writes them). For each line of reference-spec.tsv, every
# artifact matching a source pattern must exist and hold the simulator's runs table; the rows are extracted by
# reference-extract.sh into <out dir>/<file>. Stops with an error naming what is missing. With a third argument,
# writes one line per source table: file, run, artifact, sha256 of the table, rows taken.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
root="$1"; out="$2"; list="${3:-}"
mkdir -p "$out"
[ -n "$list" ] && printf 'file\trun\tartifact\ttable_sha256\trows_taken\n' > "$list"
tail -n +2 "$here/reference-spec.tsv" | while IFS=$'\t' read -r file arms cells last sources; do
  tables=()
  for source in $sources; do
    run="${source%%/*}"; pattern="${source#*/}"
    found=0
    for dir in "$root/$run"/$pattern; do
      [ -d "$dir" ] || continue
      table="$dir/$file"
      if [ ! -f "$table" ]; then echo "::error::reference: artifact $(basename "$dir") of run $run has no $file" >&2; exit 1; fi
      tables+=("$table"); found=1
      if [ -n "$list" ]; then
        taken=$("$here/reference-extract.sh" "$arms" "$cells" "$last" "$table" | tail -n +2 | wc -l)
        printf '%s\t%s\t%s\t%s\t%s\n' "$file" "$run" "$(basename "$dir")" "$(sha256sum "$table" | cut -d' ' -f1)" "$taken" >> "$list"
      fi
    done
    if [ "$found" = 0 ]; then echo "::error::reference: no artifact $pattern of run $run under $root (expired or not downloaded?)" >&2; exit 1; fi
  done
  "$here/reference-extract.sh" "$arms" "$cells" "$last" "${tables[@]}" > "$out/$file"
  echo "$file: $(($(wc -l < "$out/$file") - 1)) rows from ${#tables[@]} tables"
done
