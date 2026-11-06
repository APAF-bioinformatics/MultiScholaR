#!/bin/bash
# KIT-099 record.sh <label> <ticket> <corpus> <stop step|-> <ceiling s>: one in-process record on kit_098 (the large oracle or a
# pair population store; declaration section 5), from a clean worktree at HEAD, one attempt, with the sleep lock. Writes
# kit-099/oracles-<label>_<corpus>.rds and a line in kit-099/oracle-runs.txt.
LABEL=$1; TICKET=$2; CORPUS=$3; STOP=$4; C=$5
R=/home/doktersmol/Projects/APAF/MultiScholaR; P=/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5; S=$P/kit/scripts-099
K=docs/qa/workflow-execution/kit; OUT=$R/$K/kit-099; L=$P/core/browser-smoke.eoJzd6/library
WT=$P/kit/scratch-099-rec; [ -d $WT ] || git -C $R worktree add -q --detach $WT HEAD
git -C $WT checkout -q --detach $(git -C $R rev-parse HEAD); [ -z "$(git -C $WT status --porcelain)" ] || { echo "STOP dirty worktree"; exit 1; }
state() { echo "AC $(cat /sys/class/power_supply/ACAD/online), profile $(cat /sys/firmware/acpi/platform_profile), battery $(cat /sys/class/power_supply/BAT1/capacity)%, load $(cut -d' ' -f1 /proc/loadavg), available $(free -g | awk 'NR==2{print $7}') GiB, swap used $(free -m | awk 'NR==3{print $3}') MiB"; }
systemd-inhibit --what=idle:sleep --who="KIT-099 record" --why="$LABEL" --mode=block bash -c 'while kill -0 '$$' 2>/dev/null; do [ "$(cat /sys/class/power_supply/ACAD/online)" = 0 ] && [ "$(cat /sys/class/power_supply/BAT1/capacity)" -le 15 ] && exit 0; sleep 30; done' &
echo "record $LABEL start $(date +%T) at $(git -C $WT rev-parse --short HEAD): $(state)" | tee -a $OUT/power-state.txt
cd $WT; KIT099_INSTALLED=1 R_LIBS=$L python3 $S/trial.py $S/record-$LABEL.out $C $C Rscript $K/kit-099/oracle-record.R $LABEL $TICKET $OUT/oracles-${LABEL}_$CORPUS.rds $CORPUS $STOP; rc=$?
echo "record $LABEL end $(date +%T): $(state)" | tee -a $OUT/power-state.txt
echo "$LABEL ticket $TICKET corpus $CORPUS stop $STOP library $(basename $(dirname $L)) recorder $(git -C $WT rev-parse --short HEAD) exit $rc | $(grep -a '^recorded' $S/record-$LABEL.out | cut -d'|' -f2-) | $(tail -1 $S/record-$LABEL.out)" | tee -a $OUT/oracle-runs.txt
grep -a "^  " $S/record-$LABEL.out > $OUT/oracle-$LABEL.txt
exit $rc
