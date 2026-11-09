#!/bin/bash
# KIT-107 run 9 (check mode): the targeted and dependent set (run 1's list), test-baseline.R file check, one process per file, 4 parallel,
# against the installed kit_107 check build; the load_all half is run 1 itself.
export P=/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5; L=$(cat $P/kit/scripts-107/install.txt)/library; CK=$P/kit/scratch-107-check; out=/tmp/claude-1000/-home-doktersmol-Projects-APAF-MultiScholaR/7f1e37b6-a376-45ae-a5c8-0904bea4e8a0/scratchpad/kit-107-lead/skips-check
rm -rf $CK/MultiScholaR.Rcheck/tests $CK/tarx; mkdir -p $CK/tarx "$out/rows" "$out/skips" "$out/logs" && tar xzf $CK/MultiScholaR_*.tar.gz -C $CK/tarx && cp -a $CK/tarx/MultiScholaR/tests $CK/MultiScholaR.Rcheck/tests
cd $CK/MultiScholaR.Rcheck/tests/testthat
export NOT_CRAN=true MULTISCHOLAR_STRICT_SKIPS= TMPDIR=/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/kit/tmp-107 R_LIBS=$CK/MultiScholaR.Rcheck:$L
xargs -a /tmp/claude-1000/-home-doktersmol-Projects-APAF-MultiScholaR/7f1e37b6-a376-45ae-a5c8-0904bea4e8a0/scratchpad/kit-107-lead/targeted-files.txt -P 4 -I{} sh -c 'b=$(basename {} .R); timeout 1800 Rscript $P/kit/scratch-107-suite/tools/testing/test-baseline.R file check '"$CK/MultiScholaR.Rcheck"' '"$CK/MultiScholaR.Rcheck/tests/testthat"'/{} '"$out"' > '"$out"'/logs/$b.log 2>&1; echo "$b rc=$?" >> '"$out"'/times.txt'
echo done
