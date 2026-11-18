#!/bin/bash
# KIT-124 confirmation round: the gate's reference rows, extracted from earlier result tables (s<i>-runs.tsv) by a
# fixed text filter, so that the committed copy (scored/confirm-reference/) can be re-made from the source artifacts
# and compared byte for byte (the gate job's `verify` step). No computation: rows are selected and sorted, not changed.
#
# Usage: reference-extract.sh <arms, comma list> <cells, comma list> <last run> <runs.tsv>... > <out.tsv>
# Keeps the header once and the rows whose arm and cell are listed and whose run is at most <last run>; sorts the rows
# by cell, run, arm and subset (LC_ALL=C).
set -euo pipefail
arms="$1"; cells="$2"; last="$3"; shift 3
head -n 1 "$1"
for f in "$@"; do
  if [ "$(head -n 1 "$f")" != "$(head -n 1 "$1")" ]; then echo "reference-extract: header of $f differs" >&2; exit 1; fi
  tail -n +2 "$f"
done | awk -F'\t' -v arms=",$arms," -v cells=",$cells," -v last="$last" \
  'index(arms, "," $6 ",") && index(cells, "," $2 ",") && ($4 + 0) <= last' \
  | LC_ALL=C sort -t "$(printf '\t')" -k2,2 -k4,4n -k6,6 -k7,7
