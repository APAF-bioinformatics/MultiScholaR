#!/bin/bash
# KIT-099 first step (KIT-098 hand-off): 098a and 098b re-verified offline after the trx *_science products left
# kitOperationIdentityProducts, with kit-099/reverify.R (KIT-092s, passing the member directory). Each stage directory is reflinked into a scratch root
# (the runtime directory is only read), then re-verified and re-decided on its own build's library. No journey runs.
# Usage, from the repository root: reverify.sh <stage> [<stage> ...]
P=/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5
S=$P/kit/scratch-099-reverify
K=docs/qa/workflow-execution/kit
mkdir -p $S/runtime
for STAGE in "$@"; do
  rm -rf "${S:?}/runtime/$STAGE" "$S/repo-$STAGE"
  cp -a --reflink=always "$P/kit/runtime/$STAGE" "$S/runtime/$STAGE" || { echo "REVERIFY $STAGE copy failed"; continue; }
  BUILD=$(python3 -c "import json,sys;b=json.load(open('$K/bindings.json'));print(b['tickets'][sys.argv[1]]['build'])" "${STAGE#kit-journey-}")
  LIB=$(python3 -c "import json,sys;b=json.load(open('$K/bindings.json'));print(b['builds'][sys.argv[1]]['library'])" "$BUILD")
  echo "## $STAGE on $BUILD ($(basename "$(dirname "$LIB")"))"
  R_LIBS=$LIB Rscript $K/kit-099/reverify.R "$STAGE" "$S" 2>&1 | grep -v '^Registered S3\|^  method \|vroom$'
done
