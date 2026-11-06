#!/bin/bash
# KIT-099: R CMD check --no-tests on HEAD with the kit_098 library (no install: kit_098 is the bound build, no package
# change), then the full suite at HEAD, one process per file, 8 parallel processes (scripts-095/suite.sh), kit_098 library.
R=/home/doktersmol/Projects/APAF/MultiScholaR; P=/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5; S=$P/kit/scripts-099
E=$R/docs/qa/workflow-execution/kit/kit-099; L=$P/core/browser-smoke.eoJzd6/library; C=$(git -C $R rev-parse HEAD)
OUT=$P/kit/scratch-099-check; rm -rf "${OUT:?}"; mkdir -p $OUT; git -C $R archive --format=tar --prefix=src/ $C | tar x -C $OUT
( cd $OUT; export R_LIBS=$L; R CMD build --no-build-vignettes --no-manual src > build.log 2>&1; T=$(ls MultiScholaR_*.tar.gz | head -1); _R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-tests --no-manual --no-build-vignettes "$T" > check.log 2>&1 )
cp $OUT/MultiScholaR.Rcheck/00check.log $E/check.log; grep -E '^Status' $E/check.log
WT=$P/kit/scratch-084-suite; SO=$E/suite.txt
git -C $WT checkout -q --detach -f $C; git -C $WT clean -qfd tests
cd $WT
{ echo "# KIT-099 full suite checkpoint, one process per file, on the tree at $(git rev-parse --short HEAD), $(R --version | head -1 | cut -d' ' -f1-3), R_LIBS = kit_098 ($(basename $(dirname $L))), 8 parallel processes."; } > $SO
export R_LIBS=$L NOT_CRAN=true
ls tests/testthat/test-*.R | xargs -P 8 -I{} sh -c 'timeout 1800 Rscript '$P'/kit/scripts-078/run1.R . {} 2>&1 | grep -E "^COUNTS|^  BAD" >> '$SO
awk '/^COUNTS/ {f++; p+=$4; x+=$6; e+=$8; s+=$10} END {print "files",f,"pass",p,"fail",x,"error",e,"skip",s}' $SO | tee -a $SO
echo "test files: $(ls tests/testthat/test-*.R | wc -l)" | tee -a $SO
grep -E "^  BAD" $SO | cut -c1-200 | head -30
echo SUITE DONE
