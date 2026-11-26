#!/bin/bash
# E16 runner (7.4, r1.23): from the package root. Parses both scripts, scores each committed constructed case (cases/<case>/rows.tsv.gz)
# with score-k188.R under the case's K188_SCORE_R and K188_SCORE_CELLS, records each exit status and log, and compares every scored value
# with the prediction committed in expected.tsv. Usage: sh run-e16.sh <e16 dir> <out dir>   (both absolute)
HERE="$1"; OUT="$2"; K188="$HERE/.."
mkdir -p "$OUT"
Rscript -e 'invisible(parse(file.path(commandArgs(TRUE)[[1L]], "sim-k188.R"))); invisible(parse(file.path(commandArgs(TRUE)[[1L]], "score-k188.R"))); cat("parse: both scripts parse\n")' "$K188" || { echo "E16: a script does not parse"; exit 1; }
# the cases are independent: four at a time
score_case() {
    mkdir -p "$OUT/$1"
    K188_SCORE_CELLS="$3" K188_SCORE_R="$4" Rscript "$K188/score-k188.R" "$HERE/cases/$2" "$OUT/$1" > "$OUT/$1/run.log" 2>&1
    echo $? > "$OUT/$1/exit.txt"
}
while IFS="$(printf '\t')" read -r CASE ROWSDIR CELLS RANGES EXPECT; do
    score_case "$CASE" "$ROWSDIR" "$CELLS" "$RANGES" &
    while [ "$(jobs -r | wc -l)" -ge 4 ]; do wait -n; done
done < <(tail -n +2 "$HERE/cases.tsv")
wait
python3 -I "$HERE/check-e16.py" "$HERE" "$OUT"
