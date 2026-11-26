# KIT-274 outcome (interim: stopped after R1 on a missed prediction)

R1 run 45723980071 on tests-only 050a0402 (tables: r1-tests.tsv, r1-failures.txt).
Fail by value as predicted: T1 (3 failures), T2 (3), T3 (4), T4 (4), T5 (1), T6 (4), T7 (3), T12 (3). Pass as predicted: T8
(digest of the written file ba839e69b0d2c199d37a7e68ed0207679b0138f7751c5a6897f4e3c233da6108, product trx.tmm_logcpm,
normalised, not corrected), T9 (k 2, corrected, trx.ruviii), T11. X1 base line printed: `ruv_applied=FALSE export=ok entry=listed k=0 text=RUVIII (ruv::RUVIII) was applied ... k = 0 ...`.
MISSED PREDICTION: T10 (guard, predicted pass on base) failed: route A and route B both listed TRUE on the base.
Likely cause (inferred, not shown): k274Run builds a fresh omicNormReactiveState() per run where test-kit264 E1/E2 share one `nd`, so the
second run may not have run; the test has no precondition that the second run made a new trx_ruvs entity. No error text in the log.
Stop rule: a missed prediction stops the work; reported, not corrected.
Not checked: everything after R1 (the fix is drafted in the worktree, uncommitted; R2; R CMD check).
