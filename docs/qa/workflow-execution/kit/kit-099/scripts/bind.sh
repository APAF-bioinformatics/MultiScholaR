#!/bin/bash
# bind.sh [commit message]: the bindings regenerated from the build list in bindings.json (scripts-098/bind.sh, no new build: kit_098 is bound);
# committed alone and pushed only when a message is given.
MSG=${1:-}
R=/home/doktersmol/Projects/APAF/MultiScholaR; P=/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5; S=$P/kit/scripts-099
cd $R
ARGS=$(python3 -c "
import json,os
b=json.load(open('docs/qa/workflow-execution/kit/bindings.json'))
print(' '.join(f\"{k}={os.path.dirname(v['library'])}\" for k,v in b['builds'].items() if k.startswith('kit_')))")
KIT_065_LIFECYCLE=1 KIT_082_LIFECYCLE=1 R_LIBS=$P/core/browser-smoke.InuLC5/library Rscript docs/qa/workflow-execution/kit/bindings-generator.R $ARGS kit_journey=$P/kit/runtime/kit-journey > $S/bind.log 2>&1 || { echo "BIND FAILED"; tail -5 $S/bind.log; exit 1; }
git diff --stat -- docs/qa/workflow-execution/kit/bindings.json
[ -z "$MSG" ] && exit 0
git add docs/qa/workflow-execution/kit/bindings.json && git commit -qm "$MSG" -- docs/qa/workflow-execution/kit/bindings.json && git push -q origin performance-v2
git log --oneline -1 | cut -c1-110
