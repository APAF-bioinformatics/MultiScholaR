#!/bin/bash
# member.sh <stage> <library>: one member from a clean worktree at HEAD, waited to its decision (KIT-098's member.sh). The wait reads
# the member's log and the unit list: no process matching by name.
STAGE=$1; L=$2; MAX=${3:-4000}
R=/home/doktersmol/Projects/APAF/MultiScholaR; P=/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5; K=docs/qa/workflow-execution/kit; OUT=$R/$K/kit-100
WT=$P/kit/scratch-100-run; [ -d $WT ] || git -C $R worktree add -q --detach $WT HEAD
git -C $WT checkout -q --detach $(git -C $R rev-parse HEAD); cd $WT
state() { echo "AC $(cat /sys/class/power_supply/ACAD/online), profile $(cat /sys/firmware/acpi/platform_profile), battery $(cat /sys/class/power_supply/BAT1/capacity)%, load $(cut -d' ' -f1 /proc/loadavg), available $(free -g | awk 'NR==2{print $7}') GiB"; }
echo "$STAGE start $(date +%T) at $(git rev-parse --short HEAD): $(state)" | tee -a $OUT/power-state.txt
systemd-inhibit --what=idle:sleep --who="KIT-100 member" --why="$STAGE" --mode=block bash -c 'while kill -0 '$$' 2>/dev/null; do [ "$(cat /sys/class/power_supply/ACAD/online)" = 0 ] && [ "$(cat /sys/class/power_supply/BAT1/capacity)" -le 15 ] && exit 0; sleep 30; done' &
bash tools/profiling/run_workflow_kit.sh run $STAGE $L > $P/kit/$STAGE.launch.out 2>&1; rc=$?; echo launch exit $rc; tail -1 $P/kit/$STAGE.launch.out | cut -c1-160
if [ $rc -ne 0 ]; then echo "$STAGE refused at launch $(date +%T): exit $rc, no member ran" | tee -a $OUT/power-state.txt; exit $rc; fi
start=$SECONDS
while [ $((SECONDS-start)) -lt $MAX ]; do
  sleep 20
  if grep -qE "^DECISION|^Error|Execution halted" $P/kit/$STAGE.log 2>/dev/null && ! systemctl --user list-units --state=active --no-legend 'multischolar*' | grep -q .; then break; fi
done
echo "$STAGE end $(date +%T): $(state)" | tee -a $OUT/power-state.txt
grep -E "^DECISION|^Error" $P/kit/$STAGE.log | tail -2 | cut -c1-240
cd $R
KIT_LIMPA_TICKET_DIR=kit-100 R_LIBS=$L Rscript $R/$K/kit-053/stage-summary.R $P/kit/runtime/$STAGE kit100-$STAGE 2>&1 | tail -1 | cut -c1-120
R_LIBS=$L Rscript -e 'args <- commandArgs(TRUE); d <- readRDS(file.path(args[[1]], "decision.rds")); id <- if (dir.exists(file.path(args[[1]], "member"))) "member" else "candidate-1"; j <- readRDS(file.path(args[[1]], id, "journey-result.rds")); `%||%` <- function(a,b) if (is.null(a)) b else a; cat(basename(args[[1]]), d$status, "failing", paste(d$failing_rule %||% "none", collapse = ","), "failed_step", j$failed_step %||% "none", "error", substr(j$error %||% "none", 1, 300), "\n")' $P/kit/runtime/$STAGE 2>&1 | tee -a $OUT/member-results.txt | cut -c1-360
echo MEMBER DONE $STAGE
