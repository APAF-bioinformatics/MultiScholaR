#!/bin/bash
# KIT-099: KIT-084's ops.sh with the kit-099 evidence directory; the decision is written verbatim (kitOperationDecision) as JSON.
# ops.sh <stage> <library>: one operation stage from a clean worktree at HEAD, waited to its decision; the same wait as member.sh with the sleep lock. Reads decision.rds only.
STAGE=$1; L=$2; MAX=${3:-8000}
R=/home/doktersmol/Projects/APAF/MultiScholaR; P=/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5; K=docs/qa/workflow-execution/kit; OUT=$R/$K/kit-099
WT=$P/kit/scratch-099-run; [ -d $WT ] || git -C $R worktree add -q --detach $WT HEAD
git -C $WT checkout -q --detach $(git -C $R rev-parse HEAD); cd $WT
state() { echo "AC $(cat /sys/class/power_supply/ACAD/online), profile $(cat /sys/firmware/acpi/platform_profile), battery $(cat /sys/class/power_supply/BAT1/capacity)%, load $(cut -d' ' -f1 /proc/loadavg), available $(free -g | awk 'NR==2{print $7}') GiB"; }
echo "$STAGE start $(date +%T): $(state)" | tee -a $OUT/power-state.txt
systemd-inhibit --what=idle:sleep --who="KIT-099 member" --why="$STAGE" --mode=block bash -c 'while kill -0 '$$' 2>/dev/null; do [ "$(cat /sys/class/power_supply/ACAD/online)" = 0 ] && [ "$(cat /sys/class/power_supply/BAT1/capacity)" -le 15 ] && exit 0; sleep 30; done' &
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
R_LIBS=$L Rscript -e 'args <- commandArgs(TRUE); d <- readRDS(file.path(args[[1]], "decision.rds")); writeLines(jsonlite::toJSON(d, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = 6), file.path(args[[2]], paste0("operation-decision-", basename(args[[1]]), ".json"))); cat(basename(args[[1]]), d$status, "failing", paste(unlist(d$failing_rule), collapse = ","), "\n"); for (p in d$pairs) cat(" ", p$operation, p$population, "peak", round(p$peak_ratio, 3), "anon", round(p$peak_anon_ratio, 3), "elapsed", round(p$elapsed_ratio, 3), "not_worse", p$not_worse, "reduction", p$reduction, "\n")' $P/kit/runtime/$STAGE $OUT 2>&1 | tee -a $OUT/member-results.txt | cut -c1-300
echo MEMBER DONE $STAGE
