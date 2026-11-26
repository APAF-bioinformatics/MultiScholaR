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

## R1b (VAL-54 instrument re-send 1 of 3), after the coordinator's ruling (R1 above is kept as it happened)

Ruling quoted: "treated as an instrument defect, re-send 1 of 3 (VAL-54) ... the behaviour T10 guards (a later run supersedes the
listed export) is what KIT-264's E1/E2 assert on the same base, and that file passed 85/0 in your own R1 with one shared `nd`; your T10
differs only in building a fresh omicNormReactiveState() per run, and it has no precondition showing the second run happened. So the run
gave no result on the guarded behaviour." Terms: tests-only; (1) one `nd` per project in k274Run; (2) a precondition in T10 that the
second run committed a new trx_ruvs entity; (3) no assertion changed. Sharing `nd` changes nothing T1 to T9, T11, T12 exercise (only T10
calls k274Run). Predictions for R1b: T10 passes on base with the precondition met; all others as in R1 (T1-T7, T12 fail by value; T8, T9,
T11 pass; X1 base line as in R1). If T10 still fails with the precondition met: a result, stop. If the precondition fails: hand back.
T4 deviation (accepted by the coordinator): the press is made in the summary hook, not through `ruviii = "group"`, so the product set is
read before and after it; it asserts the same observables as the plan register's P04.

R1b result: run 45726224100 on tests-only 60565dd4 (r1b-tests.tsv, r1b-failures.txt): 53 passed, 25 failed, 0 errors (185.5 s). As predicted:
T1 (3), T2 (3), T3 (4), T4 (4), T5 (1), T6 (4), T7 (3), T12 (3) fail by value; T8, T9, T10, T11 pass. T10: both preconditions met
(the second run committed a new trx_ruvs entity), `route A listed: FALSE route B listed: FALSE` on the base: R1's T10 failure was the
instrument (no second run). Core-hours so far: R1 about 1.0 plus R1b about 0.4.
