# KIT-280 outcome: rule B for the proteomics missing-value filter, the replay marker, the flag per comparison, the corrected texts

Branch `kit-280`; code commit 1eff9c5cf (base e4a6d4aea); tests of the second pass bf6808898. Package version unchanged (0.9.0). Records: `runs.tsv`, `predictions.md` (with its two addenda), `step0.md`, `wording.md`, the run tables `r1-*`, `r1b-*`, `r2-main-*`, `r2-kit247-*`, `t2-*`.

## What was run (observed)

| Run | Id | What | Core-h | Result |
|---|---|---|---|---|
| R1 | 45820775438 | tests-only on the old R/, 4 files | 2.44 | every new test failed by value, no error; T7, T9 passed; test-prot-missing-value-counts.R 64/3 as predicted; test-kit227 3 failures; templates 123/0. Two mistakes of mine in the new tests were found (below) |
| R1b | 45827317090 | new file on the old R/ after the first correction pass | 0.87 | 78 passed / 70 failed, all by value; T5 (b) and (c) pass |
| R2 main | 45830007757 | base e4a6d4ae against head 1eff9c5c, 105 files, 5 shards, R CMD check | 46.50 | stopped on unpredicted verdicts (see "R2's reading") |
| R2 kit247 | 45830003023 | test-kit247-design-resave.R alone | 9.98 | `same`, 188/0 both sides |
| Second pass, head | 45864584883 | five edited files on 1eff9c5c + tests | 3.77 | as predicted: input-current 321/0 (log "KIT-182 F2 before=64 after=58 removed=6"), kit206 116/0, worker-figures 71/0, ruv-factors 312/0, new file 149/0 |
| Second pass, base | 45864580490 | e4a6d4ae + the edited ruv test only | 0.94 | as predicted: 311/1, `expect_length(uncorrected_only, 0L)` fails with length 2 |

Total 64.50 core-h of the ceiling 174.33 (ticket expected 55.37; the scope grew to 106 files and the second pass was added).

## The two correction passes

First (tests only, bb320e7ee; after R1): both were mistakes of mine in the new file, found from the returned values; no expectation weakened.
- Fixture F's fillers: 37 -> 38 (the ticket's F is 40 proteins, n = 240; 37 fillers gave 39 proteins, n = 234 and threshold 13.434). Counts corrected 39 -> 40 kept (T1, T3 "38 FALSE", T3b, T3d 37 -> 38 FALSE, T6, the header).
- T5 (c), a record with no marker at all: 39 -> 38 kept. The earlier comparison also removes P_low2 (one missing value of three is 33.3333 % against the cutoff 33.333 %). The ticket's "39 kept" holds only for the record with just `missing_value_comparison = counts` (case b); the added `expect_false(any(c("P_low1","P_low2") %in% kept(c_)))` states it.

Second (tests only, bf6808898; allowed by the owner, VAL-394 (1)):
- H1: the three "QC change" calls (`test-prot-norm-input-current.R`, `test-kit206-prot-norm-export-guard.R`, `test-prot-norm-worker-figures.R`) run under `MultiScholaR.missing_filter_threshold_rule = "missing"`. Reason: KIT-182's QC change assumed the filter removes proteins from the complete "small" fixture (64 -> 58); under rule B it removes none, which sent two helpers into a never-used fallback (non-automatic row names, no lineage) and left the third's table unchanged. No expectation changed; the fallback is untouched.
- H2, the one declared expected-value change (VAL-394 (1)): `test-prot-ruv-factors-in-model.R` `expect_gt(length(uncorrected_only), 0L)` -> `expect_length(..., 0L)` in the 6-per-group case, `expect_gt` kept for the 3-per-group case; printed on head: 0 and 22; on base: 2 and 24. Stricter, not looser.
- H3: T3c: `keepNA = NA` -> `keepNA = TRUE` (under `keepNA = NA`, `nzchar(NA)` is TRUE and counted the removed protein's NA; the product marks exactly one kept protein) plus `expect_true(is.na(low_intensity_groups[Pmin]))`.
- Predictions for it were committed first (predictions.md).

## R2's reading (explained after the fact, not predicted; diagnosis AUDIT-KIT-280-r2-20261127: no product defect)
- test-prot-norm-input-current.R 39 failures, test-kit206 E2b, test-prot-norm-worker-figures.R 2: the broken KIT-182 premise (H1).
- test-prot-ruv-factors-in-model.R 1: the clean-up after RUV now removes nothing in the 6-per-group case (H2).
- test-prot-norm-core.R 1,555/307 -> 1,543/319: six frozen `workflow` digests (:644, :677, :710, :744, :778, :812), each in W1 (a) and W2: 12 expectations; the clean-up keeps one more protein on the "small" fixture (64, not 63: INFERRED from kit-182/feasibility.md:49). All wait.
- test-gui-entry-proteomics.R 92/101 -> 79/114: the `dynamic_qc_tabs` UI digest changes in 93 steps (prot_dia 29, prot_lfq 21, prot_lfq_fragpipe 21, prot_tmt 21, prot_dia_limpa `import_refusal` 1), not the 73 predicted: the tabset is rendered at every step in all lanes. All wait.
- test-workflow-analysis-record.R 351 -> 353 passed: KIT-177 T5 checks each column of the DA long table, which has two more columns.
- test-omic-qc-module.R 153/7 (3 W5 digests), test-kit227 189/0, check identical both sides (0 errors, 3 warnings, 3 notes each): as predicted.

## Stored values that move (listed, none edited; for the end-of-project re-record)
- tests/testdata/gui-entry prot_dia, prot_lfq, prot_lfq_fragpipe, prot_tmt, prot_dia_limpa: UI digest in 93 steps, plus from the filter on the states digests, the filter text and the DA long-table artifact.
- test-prot-norm-core.R: six frozen `workflow` digests (12 expectations) and the 16 `protNormW1Frozen` values (already failing).
- test-omic-qc-module.R: six HTML digests (three green expectations turn red).
- Edited with authority: test-prot-missing-value-counts.R (option list and three expectations), test-kit227 G1 substitution, the second pass's tests.

## Acceptance criteria
| # | Met? | Evidence |
|---|---|---|
| A1 | met | R1 and R1b: each new test fails by value on the old code (T7, T9 pass); the base sentences printed (G1) |
| A2 | met | T1, T2, T6 pass on head (R2 main: new file 147/1, the one failure T3c's expression, corrected and 149/0 in the second pass) |
| A3 | met | T5, T5b pass: marker row, both declarations, replay by the three-way table, older records reproduce their digests |
| A4 | met | T3, T3b (RUV skip and manual), T3c, T3d, T4, T8, T8b pass |
| A5 | met | T10, T10b pass (draft sentences); NEWS bullet |
| A6 | met | T7: the tab default stays 2 |
| A7 | met | only the declared expectations edited; frozen references untouched |
| A8 | met with explanations | R2 plus the second pass: every scope file `same` or explained above; the explanations of R2's unpredicted verdicts are after the fact; check identical both sides |
| A9 | met | NEWS clause replaced; KIT-197 notes appended; ledger rows; VAL-142 and VAL-149 untouched |
| A10 | met | NEWS, ledger rows `VAL-<merge>` and `VAL-<merge+1>`, this file |
| A11 | met | man/removeRowsWithMissingValuesPercentHelper.Rd edited by hand |
| A12 | met | no file of the omics-globals pattern touched |

## For the change report
- What was wrong: the proteomics missing-value filter (DIA, label-free, TMT) counted a value at or below the intensity threshold (the 1st percentile of all values by default) as missing, while the differential analysis then used it. On `main` the threshold was rounded up to the next whole log2 unit (`main:R/func_prot_qc.R:380`), so the inconsistency was larger there. Proteins holding such a value in one compared group were more often called significant than other proteins (simulated: 1.89 [1.83, 1.96] times as often against the other kept proteins in the pooled ten settings, 1.5 to 3.3 per setting, 1.22 in complete data with no filter; measured at the 1st percentile only).
- Code (base e4a6d4aea -> 1eff9c5cf): `R/func_prot_qc_filtering_helpers.R` kernel `is_missing` (was `:243-247`: `x > threshold`, now `is.na(x)` under the current token); marker row in `R/utils_workflow_analysis_record.R`; flag in `R/utils_workflow_prot_declarations.R` and `R/func_prot_da_results_long_format.R`; both declarations' `fixed`.
- Results: proteins that hold a threshold-level value and failed a group only because of it are now kept (KIT-211: 24.3 more per simulated data set) and marked per comparison. The clean-up after RUV and strict mode count only absent values (not measured; VAL-394 (2): kept as built, stated at the wording pass). Records made before replay with the earlier rule and write no mark.
- Evidence: runs above; KIT-211's results (simulation, one generator, equal group sizes, no normalisation or RUV between filter and test).

## Not checked
- The "; " separator of `low_intensity_groups`: a group label that itself contains "; " would be split (not tested).
- The Glimma, heatmap, enrichment and Excel consumers of the two new columns beyond the scope files that passed.
- iq's MaxLFQ premise of the DIA sentence (that a protein value exists for every run with a positive peptide intensity): not read.
- The 63 -> 64 protein count after RUV on the "small" fixture: INFERRED; no run measured it.
- The `workflow` digest mechanism behind the 12 test-prot-norm-core.R expectations: explained, not re-run.
- The per-step list of the 93 lane steps beyond the digest comparison of R2.
- Real data, unequal groups, other percentiles, Windows.
- The fallback in the KIT-182/206 helpers is defective when reached (non-automatic row names, no lineage): left as it is; a later ticket may replace it by a stop().
- Rule slips: one looping `fetch.sh` run in R1 (stopped); one attempted `git stash list` and one `git checkout -- .` in a single command during the second pass (denied by the permission system before running; nothing changed); one polling-loop Monitor denied; the second-pass base worktree was first created with a wrong argument order (no effect).
