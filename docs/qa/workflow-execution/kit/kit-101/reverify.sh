#!/bin/bash
# KIT-101 (KIT-092's script): re-verifications, offline. Each member's stage directory is reflinked into a scratch root (the runtime
# directory is only read), then re-verified and re-decided by reverify.R on its own build's library. Run from the
# repository root after the bindings are regenerated. Usage: reverify.sh <stage> [<stage> ...]
P=/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5
S=$P/kit/scratch-101-reverify
K=docs/qa/workflow-execution/kit
mkdir -p $S/runtime
for ARG in "$@"; do
  STAGE=${ARG%%:*}; ALIAS=""; [ "$ARG" != "$STAGE" ] && ALIAS=${ARG#*:}
  rm -rf "${S:?}/runtime/$STAGE" "$S/repo-$STAGE"
  cp -a --reflink=always "$P/kit/runtime/$STAGE" "$S/runtime/$STAGE" || { echo "REVERIFY $STAGE copy failed"; continue; }
  BUILD=$(python3 -c "import json,sys;b=json.load(open('$K/bindings.json'));print(b['tickets'][sys.argv[1]]['build'])" "${STAGE#kit-journey-}")
  LIB=$(python3 -c "import json,sys;b=json.load(open('$K/bindings.json'));print(b['builds'][sys.argv[1]]['library'])" "$BUILD")
  echo "## $STAGE on $BUILD ($(basename "$(dirname "$LIB")"))"
  R_LIBS=$LIB Rscript $K/kit-101/reverify.R "$STAGE" "$S" $ALIAS 2>&1 | grep -v '^Registered S3\|^  method \|vroom$'
done
