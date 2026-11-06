#!/bin/bash
# bind.sh <message>: the bindings regenerated from the build list in bindings.json plus BIND_EXTRA (scripts-086/bind.sh, log here),
# committed alone and pushed.
MSG=$1
R=/home/doktersmol/Projects/APAF/MultiScholaR; P=/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5; S=$P/kit/scripts-100
cd $R
cp docs/qa/workflow-execution/kit/bindings.json $S/bindings.before.json
ARGS=$(python3 -c "
import json,os
b=json.load(open('docs/qa/workflow-execution/kit/bindings.json'))
print(' '.join(f\"{k}={os.path.dirname(v['library'])}\" for k,v in b['builds'].items() if k.startswith('kit_')))")
ARGS="$ARGS ${BIND_EXTRA:-}"
KIT_065_LIFECYCLE=1 KIT_082_LIFECYCLE=1 R_LIBS=$P/core/browser-smoke.InuLC5/library Rscript docs/qa/workflow-execution/kit/bindings-generator.R $ARGS kit_journey=$P/kit/runtime/kit-journey > $S/bind.log 2>&1 || { echo "BIND FAILED"; tail -5 $S/bind.log; exit 1; }
git diff --quiet -- docs/qa/workflow-execution/kit/bindings.json && { echo "bindings unchanged"; exit 0; }
git add docs/qa/workflow-execution/kit/bindings.json && git commit -qm "$MSG" -- docs/qa/workflow-execution/kit/bindings.json && git push -q origin performance-v2
git log --oneline -1 | cut -c1-110
