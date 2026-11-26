# KIT-213 R3 (second, on the merge with staging 4d588ba3): predictions, committed before the run

Run: full suite (`tests: all`, 512 test files on the merged tree), 20 shards per side, base = staging 4d588ba3, head =
the commit that adds this file (merge 6cff8354 plus records only); R CMD check on both sides; env
MULTISCHOLAR_TRX_MANIFEST_NOTE=1 and MULTISCHOLAR_KIT105_W3_DUMP=1 (test-prot-norm-core.R is in scope); per-file
timeouts 120 minutes for the four GUI-entry files and test-prot-norm-core.R, 45 otherwise. The ticket's R2 grep on the
merged tree matches 46 files (38 at a68fd69c); all are in this run.

Predictions (a miss stops the ticket and is reported):
1. `test-kit213-study-parameters-labels.R`: only on head, 56 passed, 0 failed.
2. Every other file: verdict `same` (failures, where any, identical on both sides), including
   `test-omic-summary-module.R` and `test-omics-summary-core.R` (on base the three trx `study_parameters` leaves match
   the stored values staging holds; on head the fix's values match the re-frozen ones).
3. Check: the same findings on both sides; nothing only on head.
4. No file times out or crashes on head only.
5. Cost: about 84 core-hours (R3 at e6791275, run 45490237295: 1,219.45 job-minutes = 81.3, plus two check jobs);
   with run 45490237295 the R3 row then holds about 165 of its 186.67 core-hour ceiling.
