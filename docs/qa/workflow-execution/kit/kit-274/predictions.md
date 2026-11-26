# KIT-274 R2 predictions (written before the R2 request)

R2: head a0804f80 (fix, merged with staging 7f17547f) against base edb15841 (the same tree with the fix reverted: staging 7f17547f plus the
tests-only file; built on branch kit-274 itself because only that branch may be pushed). 37 files: the ticket's grep gives 36 on the head (34
at 75acd421, plus test-kit213-study-parameters-labels.R from train 23, plus the new file) and test-omics-module-globals.R is added by instruction.
3 shards per side, check both, ceiling 660 job-minutes (44.00 core-hours), expected 403 (26.87).
- test-kit274-trx-no-ruv-export.R: base 53 passed / 25 failed (R1b); head 78 passed / 0 failed; verdict fixed (T1-T7, T12 fixed; T8-T11 same).
  T8's written-file digest equals R1b's ba839e69b0d2c199d37a7e68ed0207679b0138f7751c5a6897f4e3c233da6108 (A5).
- test-kit264-trx-export-current.R: same (85 passed); X1's printed line changes to `export=trx_ruviii_no_ruv/trx_error/...: <message> entry=absent k=NA text=No RUVIII integration export was recorded for this analysis.`
- Every other file, including test-workflow-lipid-declarations.R, test-workflow-metab-declarations.R, test-workflow-analysis-record.R, test-module-ci-trx-summary.R,
  test-kit207-summary-head-guard.R, test-integration-input-transcriptomics.R, test-workflow-trx-declarations.R, test-e2e-*: same.
- R CMD check: identical on both sides (new functions @noRd; no Rd, no NAMESPACE).
Any other verdict stops the work.
