# KIT-217 stage 2: outcome (stopped at stop rule R2)

Ticket `tickets/kit/KIT-217.md` r1.6 (r1.5 plus the implementer's gate answers). Branch `kit-217` from `v0.9-staging`
185ebf75 (KIT-202 merged). Stage 1 (measurement M) is on `kit-217-measure` and was not re-run.

**State: stopped before merge under stop rule R2.** The fix works as the ticket states (every new test passes; control
S equal on every criterion; check identical), but R2 shows three test files moving that 5.11 does not predict, and the
GUI-entry failure list I wrote for gate G7 in r1.6 over-predicts. One correction pass is left; it needs a scope
decision first (below).

## Commits

| commit | what |
|---|---|
| ee4acf08 | tests first: `test-kit217-intensity-filter-missing.R` (T1 to T11, guard G-C1); T-M1 to T-M3 in the three test files SD-15 names |
| 5e107a52 | tests only, instrument fix (VAL-54): character `info` in T3/T4; T6 and T10 (iv) split so their declared-step halves run before the fix |
| ba962c19 | the fix (R/, man/ by hand, config.ini, NEWS, ledger row `VAL-<merge>`, control S script) |
| 2f10c373 | control S script, instrument fix (parse filter); no R/ change |

## Runs (observed)

See `runs.tsv`. Job-minutes: R1 9.6, R1 instrument rerun 12.5, F 32.0, R2 356.1, control S rerun 14.8: 425.0 min =
28.3 core-hours (stage 2 expected 103.00, ceiling 505.33; R3 not sent).

- **R1 (45520764649, ee4acf08) and its instrument rerun (45522217087, 5e107a52).** Every new test fails on the
  commit before the fix at the value 5.6 states: kept sets {F1, F3, F4, F5, F6, F7, F8, F9} (T2, T5, T7 (v)),
  {P1, P2, P3} (T10 (i)), {Q1, Q2, Q3} (T10 (iv)), {R1, R2} (T10 (v)), F3 kept and no warning (T8), marker and
  parameters absent (T7, T11), the tab's inputs absent (T9), no condition from the declared step (T6, T10 (ii), (iv)),
  an "unused argument" error of another class from the S4 method (T6, T10 (iv)); `test-metab-qc-shared.R` fails T-M1
  (two keys recorded). T7 (ii), T10 (iii) and G-C1 pass; the other two edited files pass. In R1 itself testthat
  refused a non-character `info` in T3's loop and T4, and T6/T10 (iv) ended at the expected S4 error before their
  declared-step halves: an instrument defect (no result for those parts), fixed in tests only and re-sent once.
  Inferred: T3 and T4 call the helper with `group_labels`, which the old helper lacks, so they can only fail at the T1
  wrapper (`expect_false(is.null(out))`); 5.6's sets {G1..G4} and {U1, U2, U3} are what the old rule would return. I
  read that as the stated value; the auditor should confirm.
- **F (45523968111, ba962c19, head only).** The new file passes (344/0), the three edited files pass, and nine exposed
  files pass (partitions, lipid/metab declarations, KIT-202, KIT-172 reopen, QC core, module-CI QC).
- **R2 (45526352392, ba962c19 vs 185ebf75; 44 files, 3+3 shards, check).** 34 `same`; new file `only-on-head`
  (344 passed); `test-metab-qc-shared.R` `count-changed`, no failure; `test-workflow-omics-declarations.R`
  `message-changed` (G8: 25 failing on head, 23 on base, `r2/g8-registry-head.txt`); `test-workflow-omics-binding-steps.R`
  `message-changed` (the marker in the KIT-094 snapshot, VAL-189's kind, predicted in r1.6); `test-omic-qc-module.R`
  `new-failure` exactly the four predicted (test (b): tabs, subtabs, ids, id_counts). R CMD check: findings identical
  on both sides (0 only on head). Outside the prediction:
  1. **GUI-entry lanes**: `new-failure` in 16 of 35 (lipid) and 48 of 105 (metab) expectations: the per-step
     comparisons from `qc_intensity` to `summary_export`, as r1.5's 5.11 said. My r1.6 G7 list (all 20 steps, because
     I read the QC tab's UI digest as captured in every step) over-predicted by the four steps before the filter. The
     observed list is the narrower one.
  2. **`test-workflow-operation-registry.R`** `new-failure`, the KIT-025 scan "no registered DIA kernel draws random
     numbers": it greps function text and now reaches KIT-202's `omicsIntensityFilterRefuse`, whose message contains
     "sample(s)", through my `omicsIntensityFilterRemovedCounts` -> `omicsIntensityFilterSampleColumns` in the whole
     kernel. A false positive (no random number is drawn), but caused by this change.
  3. **`test-lipid-01k-qc-intensity-shared.R`** (2 tests) and **`test-metab-01q-qc-intensity-module-characterization.R`**
     (1 test) `new-failure`: their non-quarantined tests mock `lipidIntensityFiltering` / `metaboliteIntensityFiltering`
     with the generics' old three formals; the binding now passes the two new settings (E-3 as ruled), so the mock
     errors on unused arguments. Not among the three files SD-15 allows me to edit.
  4. Control S failed at its own parse filter (instrument defect, no result).
- **Control S rerun (45573134165, 2f10c373; R code identical to ba962c19).** Every criterion equal: (i) complete
  inputs base = head 3648/3648 (lanes and M's `none` tables at g 1 and 2, module-CI at g 1; percentiles 1, 10, 20, 50
  x 19 proportions); (ii) head = independent oracle of the shipped rule 352/352 (M's vd10/vd25 tables of L2 at g 1-2
  and M2 at g 1-4, module-CI at g 1-2; percentiles 1, 10 x proportions 0.2, 0.5), base = oracle of today's rule
  136/136; (iii) head whole = partitioned 56/56 including the removed-for-missing and removed-for-low counts;
  (iv) head under the earlier marker = base 136/136. `s/compare.tsv`, `s/filter-identity.tsv.gz`.

## Acceptance criteria

| # | status | evidence |
|---|---|---|
| A1 | met (stage 1, not re-run) | `kit-217-measure` 30f57029 |
| A2 | for the coordinator (the owner's row is being written) | ledger |
| A3 | met: tests fail first at the stated values, pass at the fix; gates G0, G1, G6, G7 answered (G7's list over-predicted, item 1) | r1, r1i1, f |
| A4 | met | `s/compare.tsv` |
| A5 | **not met**: R2 verdicts outside 5.11 (items 1-3); R3 not run | `r2/compare.tsv` |
| A6 | met by reading: `git diff 185ebf75 ba962c19 -- R/` touches only 5.5's files; the comparison lines keep `<` (the share now feeds `omicsIntensityFilterKeep`, same arithmetic); existing calls of the generics unchanged | diff |
| A7 | met: NEWS and ledger row; DESCRIPTION, lock and NAMESPACE unchanged; man/ only the four pages | diff |
| A8 | drafted below; results figures complete | this file |

## Decision needed (the correction pass)

- Item 2 can be corrected in scope: count the removed-for-missing rows in the whole kernel without calling the
  refusing column helper (the filter has already checked the columns), so the scan no longer reaches the message.
- Item 3 needs either permission to edit the two characterisation tests' mocks (add the two formals, as T-M1 did), or
  a ruling that they go to the end-of-project re-record with the other moved references. No code change avoids it
  under E-3.
- Item 1: record the observed 16/48 list as G7's answer (it equals r1.5's 5.11).
Then R2 again (correction pass) and R3.

## Not checked

- R3 (full suite); KIT-214, KIT-224, KIT-190, KIT-220, KIT-222 interactions (staging has moved to 67569959 meanwhile).
- The results text W-8 (E-4: KIT-198); the report Methods (L-15).
- The app itself (no GUI run beyond `testServer` and the GUI-entry lanes).
- Real data: none.

## For the change report

- **What was wrong.** The lipidomics and metabolomics intensity filter counted, for each feature, the samples whose
  value was below a low-intensity threshold and removed the feature when that share reached the chosen proportion. A
  missing value counted like a value above the threshold: it added nothing below but stayed in the total. So a feature
  with no value in any sample was always kept, a feature measured in one group only was treated like one with
  scattered gaps, and the filter kept more features the more values were missing. Affected: every lipidomics and
  metabolomics analysis whose tables have missing values or zeros (zeros are read as missing at import); complete
  tables are unchanged at every setting.
- **What changed.** A feature is now kept when at least g groups each have a value in at least min(n_i, m) of their
  samples (defaults m = 2, g = 1; both settings in the tab and the study parameters), and, as before, its share below
  the threshold is under the proportion. Pool groups are not counted; a design without a group column or with only
  pool groups is refused. Records made before replay unchanged.
- **Code.** Before (185ebf75): `R/func_lipid_qc.R:64-78`, `R/func_metab_qc_filtering_helpers.R:31-46` (share only),
  generics `R/func_general_s4_generics.R:410-414`, `:526-530`. After (ba962c19): kernel `omicsIntensityFilterKeep` and
  settings checks in `R/utils_workflow_omics_declarations.R` (after `omicsIntensityFilterHelperColumns`); the helpers
  and S4 methods call it; marker `intensity_filter_missing` in `R/utils_workflow_analysis_record.R`
  (`workflowReplayRuleMarkers`, `workflowRecordIntensityFilterEarlier`).
- **What changed in results (evidence).** On M's simulated tables (L2: 2 x 6 lipid, 500 features in two assays;
  vd10/vd25 value-dependent missingness; percentile 10, proportion 0.5; control S `s/filter-identity.tsv.gz`) the
  declared step removes per assay 0 to 4 features for missing values at 10 % missingness and 8 to 20 at 25 %, where
  before it removed none for that reason; with no missing values (mask `none`) nothing changes (control S (i), 3648
  equal rows). Measurement M (stage 1, 4.5 of the ticket): every rule's null rate 0.048 to 0.053; power at the matched
  rule 0.677 at g = 1 against 0.678 before; requiring 2 groups loses strongly decreased features (2 x 6, 8-fold, 25 %
  missing: power 0.712 against 0.819 at g = 1).
- **Limits.** Simulated data only (KIT-185's generator, no unwanted factor, equal groups); unequal groups, pools and
  refusals tested by unit tests, not measured; no real data.

## Correction pass (the coordinator's rulings on R2 45526352392; declared before R2 again)

- Merge of staging 3eeaf100 (4af8caf9): one conflict, in `docs/methods-decision-ledger.md`; staging's rows VAL-201 to
  VAL-225 kept as they are and this ticket's row `VAL-<merge>` placed after VAL-225, now citing VAL-204 (the owner's
  decisions). `workflowReplayRuleMarkers` was not changed on staging (KIT-190, KIT-214, KIT-224 not merged), so the
  marker list merged without hand resolution: KIT-202's `intensity_filter_samples` and this ticket's
  `intensity_filter_missing`, in that order.
- 2c31eda0: (item 2) the whole kernel's removed-for-missing count reads the design's samples directly, no longer through
  the refusing column helper; (item 3) the mocks of the S4 generics in `test-lipid-01k-qc-intensity-shared.R` (one) and
  `test-metab-01q-qc-intensity-module-characterization.R` (two) gain the generics' two new formals (before and after in
  ticket r1.7); no assertion changed. (Item 1) G7's answer is the observed 16-step list per lane (ticket r1.7).
- T3 and T4 failing on the commit before the fix at the T1 wrapper is left to the implementation auditor to judge.
- **R2 again, declared:** purpose: the same 44 files and check, head (this branch after the correction) against staging
  3eeaf100 on the same runner image, 3+3 shards; control S is not repeated (its R code path differs only in the
  removed-count's column selection, which the new tests and R2's partitions test cover; control S (iii) compared the
  removed counts at ba962c19). Ceiling 675 job-minutes. Prediction: as 5.11 with r1.7's G7 list: new file
  `only-on-head`; `test-metab-qc-shared.R` `count-changed` without failure; `test-omic-qc-module.R` `new-failure` limited
  to test (b)'s four expectations; the two GUI-entry files `new-failure` limited to the 16 steps per lane from
  `qc_intensity`; `test-workflow-omics-declarations.R` and `test-workflow-omics-binding-steps.R` `message-changed`;
  everything else (including the registry, 01k and 01q files) `same`; check findings only on head: none. Decides:
  R3 if it holds; otherwise stop and report (no further correction).

## R2 again (45575799976, 54b60b83 vs 3eeaf100): observed

- Everything as declared except one verdict: `test-workflow-operation-registry.R` is `count-changed`, not `same`. Its
  two failing expectations are the same on both sides (the `integration_mofa_model` registry lists; no entry in
  `compare-failures`); the KIT-025 random-number scan passes again (item 2 holds); the extra passes are 354 -> 358 and
  264 -> 268 in the per-parameter label and defaults checks (KIT-039, KIT-027), one each for the four new declared
  parameters. R2 at ba962c19 already showed these counts; I left them out of the prediction.
- GUI-entry: 16 of 35 per lane fail, the steps from `qc_intensity` on, in all four lanes (G7 as corrected in r1.7).
- `test-lipid-01k-…` and `test-metab-01q-…` `same` (item 3 holds). Check findings identical (0 only on head).
- **Stopped (the coordinator's rule: anything outside the prediction stops; no further correction).** R3 not sent.
