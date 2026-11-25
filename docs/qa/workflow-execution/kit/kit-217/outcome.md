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
| ba962c19 | the fix (R/, man/ by hand, config.ini, NEWS, ledger row `VAL-324`, control S script) |
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
  T3 and T4 fail on base at `expect_false(is.null(out))` (the helper before the fix has no `group_labels`), as T1; the
  kept sets {G1, G2, G3, G4} and {U1, U2, U3} are today's rule's, derived (every share is 0), not observable on base
  (implementation audit I09: accepted as met).
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
- **Since when (legacy).** Read on `main` (436cdbe8), the same rule is in four places: `R/func_lipid_qc.R:202-209` and
  `R/func_metab_qc.R:202-209`, and the duplicate helpers in `R/func_lipid_s4_objects.R:4328-4335` and
  `R/func_metab_s4_objects.R:4336-4343`. It is also in the tags: v0.5.0-rc1 (`R/func_lipid_qc.R:209`,
  `R/func_metab_qc_filtering_helpers.R:41`) and v0.9.0 (`R/func_lipid_qc.R:73`, `R/func_metab_qc_filtering_helpers.R:41`).
  For lipids it was first written in ce640728 (2026-01-18). Every released version kept features with no value
  (implementation audit I10, read there).
- **Size (control S, declared step, percentile 10, proportion 0.5, defaults).** Features removed for missing values,
  per assay: 2 groups of 6, 250 features per assay: 0 to 4 at 10 % value-dependent missingness, 8 to 20 at 25 %; 4
  groups of 3, 150 features: 0 to 3 and 9 to 16; complete tables: none.
- **Side effect (owner's decision VAL-279, option 1).** A new proteomics project's `study_parameters.txt` also prints
  the two new metabolomics settings (2 values per group, 1 group), because the proteomics import copies the whole
  package `config.ini` (every omic's sections) into the project and stores it with every state; it already prints two
  older metabolomics settings (legacy). No proteomics number, table, figure or report sentence changes. The five
  proteomics GUI-entry reference files (`tests/testdata/gui-entry/prot_dia.json`, `prot_dia_limpa.json`,
  `prot_tmt.json`, `prot_lfq.json`, `prot_lfq_fragpipe.json`) join the end-of-project re-record; a separate ticket will
  stop proteomics study parameters listing other omics' settings. Reopened records of any omic saved before this change
  report `defaults_changed` TRUE (the shipped `config.ini` changed).
- **Runs.** R1 45520764649 (VAL-54 rerun 45522217087), F 45523968111, R2 45526352392 (control S 45573134165), R2
  again 45575799976, R3 45587576015, probe 45600181863 (first attempt 45597931823, instrument defect); the merge-train
  comparison 45666403751 (13742d43 against 3a6794b0; below).
- **Limits.** Simulated data only (KIT-185's generator, no unwanted factor, equal groups); unequal groups, pools and
  refusals tested by unit tests, not measured; no real data.

## Correction pass (the coordinator's rulings on R2 45526352392; declared before R2 again)

- Merge of staging 3eeaf100 (4af8caf9): one conflict, in `docs/methods-decision-ledger.md`; staging's rows VAL-201 to
  VAL-225 kept as they are and this ticket's row `VAL-324` placed after VAL-225, now citing VAL-204 (the owner's
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
- **Coordinator's ruling (after the fact):** the registry file's `count-changed` is a declaration omission, not a result
  outside the ticket: failures identical on base and head; the extra passes are one per new declared parameter in the
  two per-parameter checks (354 -> 358, 264 -> 268), already visible in the first R2. The R2-again declaration is
  corrected accordingly (registry: `count-changed` with those counts); R3 may be sent.

## R3, declared before sending

- Merge of staging 16112917 (merge train 3: KIT-212, KIT-215, KIT-174 (a), KIT-225) as 0cf6e9e4: one conflict, the
  ledger; staging's rows kept through VAL-229, `VAL-324` after them. No other conflict.
- Does train 3 change a prediction? By reading `git diff 3eeaf100 16112917 -- R/ inst/ tests/`: it touches no function
  of this ticket. `R/utils_workflow_analysis_record.R` gains `assertWorkflowReplayDaSource` (DA replay; T7's replays
  carry no DA) and a report-source entry (`enrichment_lists_not_tested`); `test-workflow-omics-binding-steps.R` strips
  that entry from the KIT-094 comparison: the file already fails on base (VAL-189) and stays `message-changed` here.
  Publication, normalisation-reset, summary and enrichment code (KIT-212, KIT-225, KIT-215) and the RUV grouping text
  (KIT-174 (a)) are off this ticket's route; the four new test files (`test-kit174-ruv-grouping.R`,
  `test-kit212-da-current.R`, `test-kit215-not-tested-lists.R`, `test-kit225-prot-norm-reset-source.R`) are on base and
  head alike and predicted `same`. One uncertainty: if train 3 moved the GUI-entry lanes' base oracles, the two
  GUI-entry files may read `message-changed` instead of `new-failure`; either way head's failures must be the 16 steps
  per lane from `qc_intensity` and nothing else.
- **R3:** full suite, 20+20 shards, check, head = this commit (kit-217 after the merge), base = staging 16112917, same
  runner image; ceiling 2,980 job-minutes. Prediction: new file `only-on-head` (344 passed); `test-metab-qc-shared.R`
  `count-changed` without failure; `test-workflow-operation-registry.R` `count-changed` with the same failures and passes
  354 -> 358 and 264 -> 268; `test-omic-qc-module.R` `new-failure` limited to test (b)'s four expectations; the two
  GUI-entry files limited to the 16 steps per lane from `qc_intensity`; `test-workflow-omics-declarations.R` and
  `test-workflow-omics-binding-steps.R` `message-changed`; every other file `same`; check findings only on head: none.
  Decides A5; then the implementation audit.

## R3 (45587576015, 6b9131dc vs 16112917): observed, stopped

- Every declared verdict held: 476 files `same`; the new file `only-on-head` (344/0); `test-metab-qc-shared.R` and
  `test-workflow-operation-registry.R` `count-changed` (registry 354 -> 358, 264 -> 268, failures identical);
  `test-omic-qc-module.R` the four expectations of test (b); the lipid/metab GUI-entry files the 16 steps per lane from
  `qc_intensity`; `test-workflow-omics-declarations.R` and `test-workflow-omics-binding-steps.R` `message-changed`;
  check findings identical (0 only on head). `test-workflow-omics-binding-steps.R` (`message-changed`) fails one more
  expectation on head (`report_sources`: the two new report-source entries; passed 18 -> 17 in that test), and its
  `records` part moves by the two settings, the marker and the step description.
- **Outside the prediction:**
  1. `test-gui-entry-proteomics.R` `new-failure` (base 189 passed / 4 failed; head 92 / 101). In every proteomics lane
     the per-step comparisons fail from the `design` step on; the first difference is the design state's digest
     (`protein_s4_initial`). Inferred (not shown): the two metabolite keys added to `inst/config/config.ini` (5.3, decided)
     enter the configuration every route reads at the design into its states, so every proteomics state digest moves.
     Nothing in proteomics' filter or route code was changed. The lipid/metab lanes' design steps did not move, so
     their design does not read the shipped `config.ini` the same way; the inference is unproven until a probe compares
     the proteomics design state's `@args` and configuration on base and head.
  2. `test-omic-enrich-module.R` `fixed`: one test fails on base only (a second enrichment Run shows a different text
     digest than the first); proteomics enrichment, off this ticket's route; inferred instability, not shown.
- Per the coordinator's rule no further correction; the ticket goes back with these two verdicts for a ruling.

## Probe P (the coordinator's ruling on R3 item 1): observed

Script `probe-prot-config.R` (here; run from branch `run-k217-probe`, a979ae05, which adds only this script to kit-217's
head), predictions written in it before any value was read. Run 45600181863 (instrument rerun 1 of 3; the first,
45597931823, stopped before any value because the test helpers live outside the namespace). Lanes prot_tmt and prot_lfq,
every step from `import_refusal` to `summary_export`, every committed state, base 16112917 against head; tables in `p/`.

- **Observed.** Every S4 slot other than `args` (the quantitative tables, the design, ids and all other slots) is
  identical on base and head at every step and state of both lanes (1,344 of 1,344 slot comparisons: seven slots at 192 state observations). `args` differs at every state
  from `protein_s4_initial` on (192 of 192), and the whole-state digest the test computes differs exactly there
  (192 of 192). The differing `args` leaves are exactly `metaboliteIntensityFiltering/metabolites_min_values_per_group`
  ("2" on head, absent on base) and `metaboliteIntensityFiltering/metabolites_min_groups` ("1", absent), plus
  `daAnalysisUI/timestamp` at the DA states (the time of the run; it differs between any two runs). The stored
  configuration list of the state record also differs, which I had predicted it would not (130 of 192 observations equal, 62 different): at `protein_s4_initial` by
  the same two keys, at `protein_replicate_filtered` only by `output_file`, a temporary path (differs between any two
  runs). So P1 holds for every data and design component and fails for the stored configuration; P2 holds but for the
  two run-specific leaves; P3 holds.
- **R3's failures read the same way.** In all 101 failing proteomics comparisons on head the only differing part is
  `states` (the state digests); no output (displayed tables, texts, plot data, DA table digests) and no written table's
  sha256 differs, except prot_dia's four `summary_*` steps, whose integration artifacts (`.parquet` against `.tsv`) fail
  on base too.
- **Inferred cause (by reading).** The proteomics import copies the package's `config.ini` into the project and reads it
  into the session's configuration (`R/func_prot_import_detection.R:249-266`; the design import likewise,
  `R/mod_prot_design_import_helpers.R:181-208`); the design stores that configuration in the object's `@args` and the
  state record, so the two metabolite keys this ticket adds to `inst/config/config.ini` (5.3) enter every proteomics
  state. The lipid, metabolite and transcriptomics designs take their configuration from an imported design folder or
  start empty (`R/utils_workflow_omics_design.R:205-221`); their design steps did not move in R2 or R3, and the
  transcriptomics GUI-entry lanes were `same` in R3: they are not hit.
- **Proteomics results.** No proteomics number, table, figure or report sentence differs between base and head on
  these lanes: every data and design component of every state is identical (probe, two lanes), and in all five lanes
  R3 shows no displayed output and no written table moving. What differs is what a new proteomics project stores and
  prints:
  - The state arguments and the design state's stored configuration gain the two metabolite settings.
  - Read in the study-parameters writer, not observed: `study_parameters.txt` then lists `metabolites_min_values_per_group
    = 2` and `metabolites_min_groups = 1` in its `[metaboliteIntensityFiltering]` block. That block already shows the two
    existing metabolite settings, in "Parameters from Final S4 Object" and again in "Additional Configuration
    Parameters".
  - Because the shipped `config.ini` changes, the analysis record of any project, any omic, whose steps were saved
    before this change reports a different installed-defaults digest and `defaults_changed` TRUE when rebuilt.

  Limit: the TMT and LFQ lanes are probed; the DIA, DIA_LIMPA and FragPipe lanes are read from R3's failure messages
  only (state digests only).

## For the owner (item 1; not chosen here)

- **Option 1: keep the fix as it is**; the moved stored references go to the end-of-project re-record. Affected file:
  `tests/testthat/test-gui-entry-proteomics.R`, oracles `tests/testdata/gui-entry/prot_dia.json`,
  `prot_dia_limpa.json`, `prot_tmt.json`, `prot_lfq.json`, `prot_lfq_fragpipe.json`: 97 comparisons newly failing (prot_dia 22
  of its 26 failing steps, prot_dia_limpa 21, prot_tmt 18, prot_lfq 18, prot_lfq_fragpipe 18; the per-step comparisons from
  `design` on) and prot_dia's four `summary_*` comparisons already failing on base, whose message changes. Cost: no code
  change and no further run now; the re-record re-freezes these five oracles. The ledger's 5.3 decision is untouched.
- **Option 2: keep the two keys out of the configuration the proteomics design stores.** Either (a) drop the two keys
  from `inst/config/config.ini` and map them to `no_shipped_default` in the shipped-default map (as the lipid ones are;
  the S4 methods fall back to 2 and 1; the metabolomics workbooks then show no config value for them), which changes
  5.3's "two new `inst/config/config.ini` keys" (the ticket's design text, not one of the owner's SD answers; the
  defaults 2 and 1 themselves are unchanged); or (b) leave `config.ini` and filter the
  `metaboliteIntensityFiltering` section out of the proteomics configuration, a proteomics code change outside 5.5's
  files. Either is a code change needing an extension of the correction pass (the owner's), then R2 and R3 again (about
  20 and 87 core-hours on this ticket's figures). With (a) the proteomics lanes return to `same`; the metabolite
  default provenance in records becomes `no_shipped_default`.
- Transcriptomics: not hit (its design does not read the shipped configuration; R3's transcriptomics lanes `same`).

## Unrelated (R3 item 2), for routing

`test-omic-enrich-module.R`, test "(W3, F107-01) a second Run on the same tab gives the first run's product
parameters, state digest and display" (line 3693): fails on base only in R3 (45587576015): `seen$again$values` against
`seen$run$values`, display lines 2 and 8, second Run `<text> b9f6b656e874200059207826580cc3a8412b32f43527e240901e46cb8eabf978`
against the first Run's `<text> d9625c0c8f0d2e2d05bf97481af72158e56fcbcc0cb47a9ac59b41f286e9dd15`; head passes (1051/19
against base 1050/20, the other 19 failures identical). Shard 4 on both sides (file 264 of 486, dealt round-robin over
20 shards; inferred from the order); base shard ran on an AMD EPYC 9V45, head on an EPYC 9V74 (meta.txt); the file sets
`withr::local_seed(20261109L)` (line 507). Proteomics enrichment, off this ticket's route; possibly an unstable test. Nothing
done about it.

## Records repairs after the implementation audit (AUDIT-KIT-217-impl-20261125) and the owner's VAL-279

- I02, I03, I04, I09, I10: applied above (proteomics paragraph, the 1,344 and 130/62 counts, the KIT-094 snapshot's
  two parts, T3/T4, the change report's legacy, size, side-effect and runs paragraphs); `runs.tsv` corrected.
- I06: NEWS no longer says the workbook gets the rule (KIT-234 withdrew the workbooks), names the all-pool refusal, and
  states the proteomics study-parameters side effect as a known limit until the follow-up ticket.
- I07 (W-5's attribution, 29 % at g = 2 against 4 % at g = 1; the metabolite cell's own figures) and I08 (W-11's note
  is stored in the analysis record's JSON `note` of the two parameters; the record workbook's parameters sheet and the
  reports do not show notes): for KIT-198 and the owner's wording pass; no code change (the owner chose option 1).
- Owner's decision VAL-279 (2026-11-26): option 1; the five proteomics reference files are not edited and join the
  end-of-project re-record.
- Conventions (N6): many added R lines exceed 100 characters, as the edited files' own lines do; errors use
  `rlang::abort(class =)` as KIT-202's (existing-code scope rule).

## Merge of staging 3a6794b0 (bc9ffbc4) and the merge-train comparison run, declared before sending

- **Merge** (the coordinator's authorisation, 2026-11-26). Four conflicts, resolved as ruled:
  - `NEWS.md`: staging's entries kept as they are; this ticket's entry placed directly before KIT-202's (after every
    entry staging added in that section).
  - `docs/methods-decision-ledger.md`: staging's rows through VAL-280 kept; `VAL-324` after them.
  - `R/utils_workflow_analysis_record.R` (the replay-marker list): the roxygen line joined ("KIT-214 adds
    `duplicate_resolution_samples`; KIT-224 adds `ruv_feature_levels`; KIT-201 adds `treat_moderation` ...; KIT-217 adds
    `intensity_filter_missing`"); the rows staging's first (`duplicate_resolution_samples`, `ruv_feature_levels`,
    `treat_moderation`, its closing `))` turned into `),`), then `intensity_filter_missing` closing the list; after the
    list KIT-201's three functions (`workflowTreatModerationCurrent`, `workflowTreatModerationRecord`,
    `workflowDaReplayRuleOptions`), then `workflowRecordIntensityFilterEarlier`.
  - `R/utils_workflow_omics_declarations.R` (adjacency, predicted by the implementation audit, section 6): after
    `omicsIntensityFilterHelperColumns`, staging's KIT-214 block (`omicsDesignSampleColumnsOf`,
    `omicsDesignNumericColumnsOf`, `omicsDuplicateResolutionSampleColumns`, `omicsDuplicateResolutionInput`) first, then
    this ticket's kernel block (`omicsIntensityFilterCountsInForce` to `omicsIntensityFilterRemovedMissing`); no other
    change. Checked: the merged file against staging adds exactly this branch's 220 lines against 185ebf75.
- **Run, declared.** Purpose: the change on the merged tree. Head = `kit-217` after this records commit (bc9ffbc4 plus
  records only), base = staging 3a6794b0 (abfc99ac plus one ledger row); same runner image; R CMD check on both sides;
  4+4 shards; env `MULTISCHOLAR_TRX_MANIFEST_NOTE=1`, `MULTISCHOLAR_KIT105_W3_DUMP=1`; ceiling 900 job-minutes
  (60 core-hours). Scope: 57 files = the ticket's 5.10 grep re-run at staging (41), the 44 files of R2, and the
  implementation audit's additions (proteomics and transcriptomics GUI-entry, `test-omic-enrich-module.R`, the registry
  and analysis-record files, the KIT-201, 214, 220, 222, 224, 234, 235 files), and the files the grep now finds
  (KIT-174, 212, 243). Base counts predicted from staging's own checkpoint-20261126 table.
- **Predictions.** `test-kit217-intensity-filter-missing.R` `only-on-head`, 344 passed, 0 failed; and
  `test-kit214-sample-columns.R` passes in full on head (176/0; the adjacency guard); `test-metab-qc-shared.R`
  `count-changed` without failure; `test-workflow-operation-registry.R` `count-changed`, same failures, the two
  per-parameter checks +4 each; `test-omic-qc-module.R` `new-failure` limited to test (b)'s four expectations;
  `test-gui-entry-lipidomics.R` 16 of 35 and `test-gui-entry-metabolomics.R` 48 of 105 failing on head, the per-step
  comparisons from `qc_intensity` on, `states` only (`new-failure`); `test-gui-entry-proteomics.R`: on head every
  per-step comparison from `design` on in all five lanes fails, where the difference from base is the state digests
  (`states`) only, base's own failing comparisons (25 at the checkpoint) unchanged in every other part (`new-failure`);
  `test-gui-entry-transcriptomics.R` `same`; `test-workflow-omics-declarations.R` and
  `test-workflow-omics-binding-steps.R` `message-changed` (I04's content); `test-omic-enrich-module.R` `same`, or I11's
  known-unstable failure on one side; every other file `same`; check findings only on head: none. Decides: the
  auditor's closure check.

## Merge-train comparison run (45666403751, 13742d43 vs staging 3a6794b0): observed

Every declared prediction held (tables in `mt/`):
- 48 files `same`, among them `test-kit214-sample-columns.R` (176/0 on both sides: the adjacency guard),
  `test-gui-entry-transcriptomics.R` (16/13 on both) and `test-omic-enrich-module.R` (1011/59 on both).
- `test-kit217-intensity-filter-missing.R` `only-on-head`, 344 passed, 0 failed.
- `test-metab-qc-shared.R` `count-changed`, no failure; `test-workflow-operation-registry.R` `count-changed`, the same
  two failures, passes 354 -> 358 and 264 -> 268.
- `test-omic-qc-module.R` `new-failure`, test (b)'s four expectations.
- `test-gui-entry-lipidomics.R` 16 of 35 and `test-gui-entry-metabolomics.R` 48 of 105 failing on head, every lane from
  `qc_intensity` on, `states` only.
- `test-gui-entry-proteomics.R` `new-failure`: base 168 passed / 25 failed, head 92 / 101; every per-step comparison from
  `design` on fails in all five lanes; the only part that differs from base is `states` (every artifact difference in
  head's messages, e.g. the DA long tables of prot_dia, prot_tmt, prot_lfq and prot_lfq_fragpipe at `da`, appears with the
  same digests in base's).
- `test-workflow-omics-declarations.R` and `test-workflow-omics-binding-steps.R` `message-changed`.
- R CMD check: 6 findings, identical on both sides, none only on head.

Budget: stage 2 has used 2,102.0 + 456.9 = 2,558.9 job-minutes = 170.6 core-hours, against 103.00 expected and 505.33
ceiling (the ticket's ceiling, 514.52 with M).
