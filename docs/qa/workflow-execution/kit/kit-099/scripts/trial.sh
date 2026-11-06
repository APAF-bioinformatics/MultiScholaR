#!/bin/bash
# KIT-099 trial.sh <samples> [<samples> ...]: the stepped trial (declaration section 2), points in ascending order, one attempt
# each, in process on kit_098 from a clean worktree at HEAD: kit-099/oracle-record.R with ticket 099l's settings on the point's
# corpus (trx_60k_<samples>), under trial.py (ceiling per point, step limit 1,800 s), with the sleep lock. Stops at the first
# point that does not complete inside the limits. Writes kit-099/size-trial.txt and kit-099/power-state.txt.
R=/home/doktersmol/Projects/APAF/MultiScholaR; P=/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5; S=$P/kit/scripts-099
K=docs/qa/workflow-execution/kit; OUT=$R/$K/kit-099; L=$P/core/browser-smoke.eoJzd6/library
WT=$P/kit/scratch-099-rec; [ -d $WT ] || git -C $R worktree add -q --detach $WT HEAD
git -C $WT checkout -q --detach $(git -C $R rev-parse HEAD); [ -z "$(git -C $WT status --porcelain)" ] || { echo "STOP dirty worktree"; exit 1; }
state() { echo "AC $(cat /sys/class/power_supply/ACAD/online), profile $(cat /sys/firmware/acpi/platform_profile), battery $(cat /sys/class/power_supply/BAT1/capacity)%, load $(cut -d' ' -f1 /proc/loadavg), available $(free -g | awk 'NR==2{print $7}') GiB, swap used $(free -m | awk 'NR==3{print $3}') MiB"; }
systemd-inhibit --what=idle:sleep --who="KIT-099 trial" --why="stepped trial" --mode=block bash -c 'while kill -0 '$$' 2>/dev/null; do [ "$(cat /sys/class/power_supply/ACAD/online)" = 0 ] && [ "$(cat /sys/class/power_supply/BAT1/capacity)" -le 15 ] && exit 0; sleep 30; done' &
declare -A CEIL=([24]=900 [100]=1800 [200]=2700 [400]=3600 [600]=5400 [800]=7200)
for N in "$@"; do
  NNN=$(printf %03d $N); LABEL=trx99t$NNN; KEY=trx_60k_$NNN; C=${CEIL[$N]}
  echo "trial $N start $(date +%T) at $(git -C $WT rev-parse --short HEAD): $(state)" | tee -a $OUT/power-state.txt
  cd $WT; KIT099_INSTALLED=1 R_LIBS=$L python3 $S/trial.py $S/trial-$NNN.out $C 1800 Rscript $K/kit-099/oracle-record.R $LABEL 099l $P/kit/runtime/trial-$LABEL.rds $KEY; rc=$?
  echo "trial $N end $(date +%T): $(state)" | tee -a $OUT/power-state.txt
  cd $R
  R_LIBS=$L Rscript -e '
    a <- commandArgs(TRUE); out <- readLines(a[[1]]); tail <- strsplit(out[length(out)], " ")[[1]]; v <- function(k) tail[[match(k, tail) + 1L]]
    steps <- regmatches(out, regexpr("^  \\S+\\s+\\S+\\s+rows .*s$", out, perl = TRUE))
    sec <- as.numeric(sub(".*\\s([0-9]+) s$", "\\1", steps)); names(sec) <- sub("^  (\\S+).*", "\\1", steps)
    rows <- sub(".*rows (\\S+).*", "\\1", steps); names(rows) <- names(sec)
    raw <- as.integer(rows[["design"]] %||% NA); kept <- if ("expression_filter" %in% names(rows)) as.integer(rows[["expression_filter"]]) else NA_integer_
    n <- as.integer(a[[2]]); rss <- as.numeric(v("RSS_KB")) * 1024
    inside <- identical(a[[3]], "0") && all(sec < 1800) && rss < 15.5 * 1024^3 && identical(v("STOP"), "exit")
    cat(sprintf("%s samples %d | raw genes %s retained %s | raw cells %s retained cells %s | exit %s stop %s | elapsed %s s | peak RSS %.2f GiB | steps %s | running step %s s | inside limits %s\n",
      a[[4]], n, raw, kept, format(raw * n, big.mark = ","), format(kept * n, big.mark = ","), a[[3]], v("STOP"), v("ELAPSED"), rss / 1024^3,
      paste(names(sec), sec, sep = "=", collapse = ","), v("RUNNING_STEP_SECONDS"), inside))
    quit(status = if (inside) 0L else 3L)' $S/trial-$NNN.out $N $rc $KEY 2>&1 | grep -v "^Registered\|^  method" | tee -a $OUT/size-trial.txt
  [ ${PIPESTATUS[0]} -eq 0 ] || { echo "TRIAL STOP at $N: not inside the limits"; break; }
done
echo TRIAL DONE
