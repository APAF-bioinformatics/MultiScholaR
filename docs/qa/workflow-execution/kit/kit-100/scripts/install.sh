#!/bin/bash
# install.sh <commit>: the kit_100 install from a pinned worktree at the commit (no uncommitted work from other agents), from kit_098's library.
R=/home/doktersmol/Projects/APAF/MultiScholaR; P=/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5; S=$P/kit/scripts-100; C=$1
WI=$P/kit/scratch-100-install; [ -d $WI ] || git -C $R worktree add -q --detach $WI $C; git -C $WI checkout -q --detach -f $C; git -C $WI status --short | head -3
cd $WI && bash tools/profiling/run_workflow_kit.sh install $P/core/browser-smoke.eoJzd6/library $P/core > $S/install.log 2>&1 || { tail -n 8 $S/install.log | cut -c1-200; echo STOP install; exit 1; }
RUN=$(grep -E "^INSTALLED " $S/install.log | tail -1 | cut -d' ' -f2); echo "$RUN" > $S/install.txt; echo "installed $RUN $(cat $RUN/source-commit.txt)"
