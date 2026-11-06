#!/bin/bash
# run.sh <stage>: one KIT-100 member on kit_100 (member.sh, KIT-098.s run.sh), then its evidence committed and pushed.
R=/home/doktersmol/Projects/APAF/MultiScholaR; P=/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5; S=$P/kit/scripts-100; K=docs/qa/workflow-execution/kit
L=$(cat $S/install.txt)/library; m=$1
bash $S/member.sh $m $L || { echo "STOP $m"; exit 1; }
[ -f $P/kit/runtime/$m/decision.rds ] || { echo "STOP no-decision-$m"; exit 1; }
cd $R; git add $K/kit-100 && git commit -qm "KIT-100: evidence: $m" -- $K/kit-100 && git push -q origin performance-v2; git log --oneline -1 | cut -c1-80
echo RUN DONE $m
