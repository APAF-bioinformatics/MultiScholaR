#!/bin/bash
# check.sh <commit>: R CMD check --no-tests on the commit, kit_100's library (KIT-099's check-suite.sh, check half).
R=/home/doktersmol/Projects/APAF/MultiScholaR; P=/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5; S=$P/kit/scripts-100
E=$R/docs/qa/workflow-execution/kit/kit-100; L=$(cat $S/install.txt)/library; C=$1
OUT=$P/kit/scratch-100-check; rm -rf "${OUT:?}"; mkdir -p $OUT $E; git -C $R archive --format=tar --prefix=src/ $C | tar x -C $OUT
( cd $OUT; export R_LIBS=$L; R CMD build --no-build-vignettes --no-manual src > build.log 2>&1; T=$(ls MultiScholaR_*.tar.gz | head -1); _R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-tests --no-manual --no-build-vignettes "$T" > check.log 2>&1 )
cp $OUT/MultiScholaR.Rcheck/00check.log $E/check.log; grep -E '^Status' $E/check.log
