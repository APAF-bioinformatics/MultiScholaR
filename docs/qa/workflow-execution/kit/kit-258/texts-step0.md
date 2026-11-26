# KIT-258 follow-up "Texts (VAL-306)": step 0 and predictions

Bounded work item (declared before any edit). Outcome: the four in-app texts that still state KIT-214's rule ("highest
average intensity") state the rule as merged (VAL-254, VAL-284, VAL-337): the row measured (neither missing nor 0) in the
most samples; ties by the highest mean of the measured values; then the first row in table order. Owner: ledger VAL-306.
Files: the strings and roxygen comments of the six R files below, two `man/*.Rd` pages by hand (roxygen2 is not run), the
one declared test edit, one new test file, records. No behaviour change. Exit: the new tests pass on head and fail on
staging; the regression scope shows no unpredicted change base against head; records written. Branch `kit-258b` from
origin/v0.9-staging 75acd4219 (train 21 head).

## 1. Every place that states the duplicate rule in words (git grep at 75acd4219, data/benchmark excluded)

In scope (the four texts):

| Text | Site | Today |
|---|---|---|
| (1) Duplicates help | `R/utils_omic_gui_config.R:1263-1264` (`duplicates_help`, both omics through `noun`), shown by `R/mod_omic_qc_tabs.R:234` | "...keeping the feature with the highest average intensity across samples." |
| (2) Result text | `R/mod_omic_qc_tabs.R:331` (`omicQcDuplicatesText`) | "Strategy: Keep feature with highest average intensity" |
| (3) Console | `R/func_omics_qc_shared.R:182`, `:213` (shared kernel) | "...by keeping highest average intensity feature...", "Removed %d lower-intensity duplicate feature row(s)." |
| (3) Console | `R/func_lipid_qc_support_helpers.R:43`, `:74` (lipid kernel) | same two |
| (3) Console | `R/func_lipid_s4_duplicate_helpers.R:294`, `:318` (lipid S4 method, non-ITSD) | "...non-ITSD ID column '%s' by keeping highest average intensity feature...", "Removed %d lower-intensity duplicate non-ITSD feature row(s)." |
| (3) Console | `R/func_metab_s4_qc_methods.R:381`, `:404` (metabolomics S4 method, non-ITSD) | same two |
| (4) Help titles | `R/func_omics_qc_shared.R:143`, `R/func_lipid_qc_support_helpers.R:4`; `man/resolveDuplicateFeaturesByIntensity.Rd:5`, `man/resolveLipidDuplicateFeaturesByIntensity.Rd:5` (`\title`) | "Resolve Duplicate Features by Keeping Highest Average Intensity" |

Already correct (the descriptions, KIT-258): `R/func_omics_qc_shared.R:145-155`, `R/func_lipid_qc_support_helpers.R:6-16`,
the two Rd `\description`/`\value`, NEWS.md:705-708 (states the new rule and what it replaces), `inst/` and README (no
statement of the rule: grep for "average intensity", "highest average", "highest mean", "duplicate" with keep/rule words).

Found, not in the four texts, left untouched (recorded for the coordinator):
- `R/utils_workflow_omics_declarations.R:1853`: the declared step's description "Resolved duplicate %s features by keeping
  highest intensity" is stored in each state's history and in the analysis record
  (`description` in `docs/qa/workflow-execution/kit/kit-235/.../*-record.json`). Changing it is a fifth text and may move
  more stored values; not in VAL-306. Out of scope.
- `R/utils_workflow_lipid_declarations.R:58`: a code comment still says "the row with the highest mean intensity" (not
  shown to a user). Out of scope.
- `tools/refactor/staging/**`: old refactor snapshots, not loaded by the package. Out of scope.
- Test names (not texts): `test-kit214-sample-columns.R:326`, `test-lipid-01-qc-filtering-helpers.R:426`,
  `test-lipid-03-s4-duplicate-helpers.R:107, :216`, `test-metab-01b-...:89`, `test-metab-01m-...:58`. Names only; not edited.

## 2. Every test, stored value, oracle or digest that holds one of the four texts, or a digest of a screen that shows them

Searched: whole `tests/` (git grep, the gui-entry JSON included, `tests/testdata/mci`), by the text's words, the functions
and the module's names.

| Holder | Line(s) | Holds | Fails on staging now? | Effect of the text changes (prediction) |
|---|---|---|---|---|
| `tests/testthat/test-omic-qc-module.R`, test (b), `omicQcFrozen$subtabs` | `:171` (`lipid_duplicates`), `:175` (`metab_duplicates`) | sha256 of `mod_*_qc_duplicates_ui("x")` (shows the help) | the `subtabs` part: YES (KIT-217 changed the intensity sub-tab) | digests of the two duplicates sub-tabs change; the part already fails: message changed only |
| same, `omicQcFrozen$tabs` | `:160` (`lipid_qc`), `:161` (`metab_qc`) | sha256 of the rendered QC tabset (holds each sub-tab's UI, INFERRED) | the `tabs` part: YES (KIT-217) | message changed only |
| same, `ids`, `id_counts` | `:187`, `:194` | element ids | YES (KIT-217) | unchanged (the texts hold no ids) |
| same, W5 tests | `:399-` | proteomics' tabsets only | no | unchanged |
| `tests/testdata/gui-entry/lipid_canonical.json` | 15 values, keys `lipidomics_workflow-qc-duplicates-resolution_results`, lines 773 ... 3789 (steps `qc_duplicates` to the last, 15 of 20) | the result text (2) | the 15 steps: YES, every step from `qc_intensity` on fails (KIT-217: 16 of 35, `states` only) | `outputs` of those steps now differ too: message changed; counts 19/16 unchanged |
| `metab_lc.json` | 15 values, lines 773 ... 3736 | same | YES (48 of 105 over the three lanes = 16 per lane) | same |
| `metab_gc.json` | 15 values, lines 759 ... 3562 | same | YES | same |
| `metab_combined.json` | 15 values, lines 773 ... 3736 | same | YES | same |
| `tests/testthat/test-lipid-01-qc-filtering-helpers.R` | `:445` | the console message of the lipid kernel | no (58 passed, 0 failed, 2 skipped) | DECLARED TEST EDIT (the register and VAL-306): "Removed 1 duplicate feature row" |

The `qc-dynamic_qc_tabs` digest captured in the gui-entry JSON is the same in all four lanes (`cad1d8f9...`), so it cannot
hold the help (the help differs by omic); the captures hold the result text only. The 60 values = 15 steps x 4 files.
No other test, oracle, digest or fixture holds the four texts: `test-module-ci-*-qc.R` and
`test-module-ui-wrapper-direct-shared.R` match button ids and the tab heading only; `test-omic-gui-config.R:588` names
the keys only; the S4 and kernel tests assert row choices, not messages (grep for "Resolving duplicates", "Removed",
"lower-intensity" in `tests/`: only `:445`).

Counts on staging now (newest records: `checkpoint-20261126/tests.tsv` and `merge-train-20261126m-compare.tsv`, train 22
`merge-train-20261126o`): `test-omic-qc-module.R` 156 passed / 4 failed; `test-gui-entry-lipidomics.R` 19/16;
`test-gui-entry-metabolomics.R` 57/48 (same on base and head in train m); `test-lipid-01-qc-filtering-helpers.R` 58/0 (2
skipped); `test-omics-module-globals.R` 14/0.

STORED VALUES ARE NOT EDITED (VAL-306): the five digests above and the 60 captures wait for the end-of-project re-record.

## 3. Regression scope (R2): the grep

```
git grep -l -E 'resolveDuplicateFeatures|resolveLipidDuplicateFeatures|resolveDuplicateFeaturesByIntensity|omicQcDuplicatesText|duplicates_help|duplicates_resolved|qc_duplicates|mod_(lipid|metab|omic)_qc_duplicates|duplicate_resolution|keepDuplicateResolutionRow|omicsDuplicateResolution|Duplicate Resolution|resolution_results' -- tests/testthat
```
gives 39 files: 3 helpers (`helper-gui-entry.R`, `helper-workflow-omics-residency.R`, `helper-workflow-omics-worker.R`; their
readers are in the list or named below) and 36 test files. Six of those hold `testthat::skip("quarantined...` at the top
and are left out: test-lipid-01b, 01c, 01f, metab-01o, 01p, 01s. The R2 list is the other 30 test files, plus the three
named by the task and not found by the grep (`test-gui-entry-lipidomics.R`, `test-gui-entry-metabolomics.R`, which read
the helper's duplicate steps, and `test-omics-module-globals.R`) and the new file: 34 names, in the request.

## 4. Predictions, committed before R1

R1 (tests only: commit holds the new test file, the declared edit of `test-lipid-01-qc-filtering-helpers.R:445`, this
record; nothing under `R/` or `man/`). Head = that commit, base = 75acd4219. One shard per side, no check, files:
`test-kit258b-duplicate-texts.R`, `test-lipid-01-qc-filtering-helpers.R`. Ceiling 40 job-minutes; expected about 26
(2 jobs, install about 10-12 minutes each, files under 3 minutes); about 1.7 core-hours.
- `test-kit258b-duplicate-texts.R`: `only-on-head`; every expectation FAILS, none errors, none passes: 0 passed, 34 failed
  (help 4 x 2, result text 3 x 2, kernel messages 3 x 2, S4 messages 3 x 2, titles 4 x 2).
- `test-lipid-01-qc-filtering-helpers.R`: base 58/0 (2 skipped); head 57/1 (the pin at `:445` fails: no message matches
  "Removed 1 duplicate feature row"); verdict `new-failure`, declared.
- Decides: that each new expectation falsifies the old texts (a pass of any new expectation on staging, or an error, is
  an instrument defect to repair under the brief's four conditions; a fail count other than 34 with all failing is the
  same).

R2 (head = the text commit on top, base = 75acd4219; 34 files, 2 shards per side, check on; ceiling 220 job-minutes;
expected about 190 job-minutes = about 12.7 core-hours; R1 + R2 about 14.4 core-hours, KIT-258 19.62 + 14.4 = 34 of 106):
- `test-kit258b-duplicate-texts.R`: `only-on-head`, all pass (34 passed, 0 failed).
- `test-lipid-01-qc-filtering-helpers.R`: base 58/0 (2 skipped) with the old pin on the old code; head 58/0 (2 skipped)
  with the new pin on the new code: `same`.
- `test-omic-qc-module.R`: 156/4 on both sides; verdict `message-changed`: the `tabs` and `subtabs` expectations (4 of the
  4 failing: tabs, subtabs, ids, id_counts) keep failing; `tabs` and `subtabs` show new head digests, `ids` and `id_counts`
  unchanged. No new failure.
- `test-gui-entry-lipidomics.R` 19/16 and `test-gui-entry-metabolomics.R` 57/48 on both sides; verdict
  `message-changed`: the per-step comparisons of the 15 steps from `qc_duplicates` on differ now in `outputs` as well as
  `states`. No new failure, no fix.
- Every other file in the list: `same` (counts and messages): `test-omics-module-globals.R` 14/0, `test-kit258-duplicate-completeness.R`
  162/0, `test-kit214-sample-columns.R` 176/0, the two declarations files 259/0 each, `test-workflow-omics-declarations.R`
  45/23 (same failures), `test-workflow-kit-harness.R` 734/1 and `test-workflow-operation-registry.R` 1107/2 as at the
  newest staging record (same on both sides; a count differing from my number but equal on both sides is `same`).
- R CMD check: identical findings base and head (the man/ edits change `\title` text only).
- A stop: any new failure, any count change, any message change in a file not named above, or a check finding only on head.

## 5. Addition ruled by the coordinator inside VAL-306's intent: the declared step's description (T02, T03, T06, T07)

Tests first at `f0ca0ef4` (tests only). Predictions on the current head `8496674a` (R/ unchanged), declared before the change
and the run:
- New block (5), 3 expectations x 2 omics (rule's first part, its mean part, none of the earlier words): all 6 FAIL by value,
  the assertion printing "Resolved duplicate lipid features by keeping highest intensity" (and the metabolite text); the
  third expectation fails too because "keeping highest intensity" matches the earlier-words pattern.
- T03, help: "(0|zero)[^.]*not measured" x 2 omics and "keeps its own name and annotation" x 2 omics: all 4 already PASS
  on 8496674a (the applied help says "an intensity of 0 counts as not measured" and "The kept row keeps its own name and
  annotation"); they guard a later rewording, they do not fail first.
- T03, result text: "(0|zero)[^.]*not measured" x 2 prefixes: both already PASS ("a missing value or 0 is not measured").
- So on 8496674a plus these tests: the file reads 40 passed, 6 failed (the six of block 5); on the
  change: 46 passed, 0 failed.

Run (R3 of this follow-up): base = `8496674a`, head = the change commit; five files, 1 shard per side, no check; ceiling 60
job-minutes, expected about 30 (about 2 core-hours):
- `test-kit258b-duplicate-texts.R`: base 34/0, head 46/0 (`count-changed`, no failure).
- `test-workflow-omics-binding-steps.R`: 57/2 on both sides, `message-changed`: the duplicate step's recorded `description`
  differs in both omics (a stored laptop snapshot, `tests/testthat/fixtures/kit094-characterisation.rds`, not edited).
- `test-kit258-duplicate-completeness.R` 162/0, `test-workflow-lipid-declarations.R` 259/0 and
  `test-workflow-metab-declarations.R` 259/0: `same`.
- Nothing else moves. (The register's "the two declarations files" is read as the lipid and the metabolomics declarations
  files; `test-workflow-omics-declarations.R` is not in the run: its stored registry holds no description.)
