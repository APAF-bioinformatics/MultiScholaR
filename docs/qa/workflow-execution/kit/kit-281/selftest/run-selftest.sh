#!/bin/bash
# KIT-281 scorer self-test runner (E4 of L1): from the package root. Parses nothing and simulates nothing: it scores each committed
# constructed case (cases/<case>/) with score-k281.R under the case's environment, records each exit status and log, and compares every
# scored value with the prediction committed in expected.tsv (written by make-selftest.py before the scorer ran).
# Usage: bash run-selftest.sh <selftest dir> <out dir>   (both absolute)
HERE="$1"; OUT="$2"; K281="$HERE/.."
mkdir -p "$OUT"
score_case() {
    mkdir -p "$OUT/$1"
    # a clean environment apart from the case's: no ambient K281_* or K211_* variable reaches the scorer
    env -u K281_RANGE -u K281_G3_TABLE -u K281_GUARD_PROBE -u K211_B_BOOT $(echo "$3" | tr ';' ' ' | sed "s#@HERE@#$HERE#g") \
        Rscript "$K281/score-k281.R" "$HERE/cases/$2" "$OUT/$1" > "$OUT/$1/run.log" 2>&1
    echo $? > "$OUT/$1/exit.txt"
}
while IFS="$(printf '\t')" read -r CASE ROWSDIR ENVSTR EXPECT; do
    score_case "$CASE" "$ROWSDIR" "$ENVSTR" &
    while [ "$(jobs -r | wc -l)" -ge 4 ]; do wait -n; done
done < <(tail -n +2 "$HERE/cases.tsv")
wait
python3 -I "$HERE/check-selftest.py" "$HERE" "$OUT"
