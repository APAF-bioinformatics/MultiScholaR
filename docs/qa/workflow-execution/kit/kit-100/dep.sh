#!/bin/bash
# dep.sh [library]: KIT-100's dependents (test files naming the defaults, the variants, the log2, filter and logCPM wrappers, and
# every lipid, metab and trx module and declaration test), one process per file, 6 in parallel, on the source tree.
cd /home/doktersmol/Projects/APAF/MultiScholaR
P=/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5; S=$P/kit/scripts-100
export R_LIBS=${1:-$P/core/browser-smoke.eoJzd6/library} NOT_CRAN=true
{ grep -lE "default_variant|intensity_filter|_log2|Log2|partitioned|runWorkflowOmicsReferenceStep|workflowQualificationVariant|trx_expression_filter|trx_tmm_logcpm|TrxReference|trxRunProduct|TrxNorm|TrxQcFilter|kit089WholeIntensity" tests/testthat/test-*.R; ls tests/testthat/test-module-ci-lipid-*.R tests/testthat/test-module-ci-metab-*.R tests/testthat/test-module-ci-trx-*.R tests/testthat/test-workflow-omics-*.R tests/testthat/test-*trx*.R 2>/dev/null; } | sort -u > $S/dep.txt
: > $S/dep.out
cat $S/dep.txt | xargs -P 6 -I{} sh -c "timeout 1800 Rscript $P/kit/scripts-078/run1.R . {} x 2>&1 | grep -E '^COUNTS|^  BAD|^---' >> $S/dep.out"
awk '/^COUNTS/ {f++; p+=$4; x+=$6; e+=$8; s+=$10} END {print "files",f,"pass",p,"fail",x,"error",e,"skip",s}' $S/dep.out
grep -E "^  BAD" $S/dep.out | cut -c1-200 | head -30
for f in $(cat $S/dep.txt); do b=$(basename $f); grep -q "COUNTS $b " $S/dep.out || echo "NO COUNTS $b"; done
rm -rf tests/testthat/_problems tests/testthat/testthat-problems.rds; git checkout -- tests/testthat/results 2>/dev/null; rm -f aa_seq_tbl.tsv; echo DEP DONE
