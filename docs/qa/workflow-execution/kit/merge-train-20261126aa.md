# Merge train 37 of 2026-11-26: KIT-242, then KIT-260 (full suite, the checkpoint after trains 35 and 36)

**Not stacked.** Branched from `origin/v0.9-staging` at `918174732` (fetched before building; it holds trains 34 to 36: KIT-273, KIT-216, KIT-266 and the ledger commits to VAL-383).
Branch `merge-train-20261126aa`, worktree `/home/doktersmol/Projects/APAF/MultiScholaR-train37`. (The names `merge-train-20261127a/b/c` carry a wrong day by the coordinator's mistake; today is 2026-11-26.)

## Train

| Order | Ticket | Branch head | Merge commit | Ledger rows | Audit |
|---|---|---|---|---|---|
| 1 | KIT-242 | `origin/kit-242` 167df89a | 35f902ea2 | `VAL-<train37-1>` (the ticket's row); the coordinator numbers it at landing (next free: VAL-384) | AUDIT-KIT-242-impl-20261126 (closure check: ready after record fixes) |
| 2 | KIT-260 | `origin/kit-260` d6d927a7 (holds kit-242 167df89a and staging 91817473) | b65f67706 | `VAL-<train37-2>` (the ticket's row) and `VAL-<train37-3>` (the stored-reference row, option (b)) | AUDIT-KIT-260-impl-20261126 (ready after M03, M05 and the owner's by-file answer M01) |

Merge head (code-bearing, bound by the run): `b65f67706`; the docs-only commit after it (`f668c6cbf`, the record fixes) and the later commits change nothing under `R/`, `tests/`, `inst/`, `man/`.

### Conflicts and how they were resolved

- Merge 1 (`git merge --no-ff origin/kit-242`): conflicts in `NEWS.md` and `docs/methods-decision-ledger.md` only; `R/utils_workflow_dia_products.R` and `R/utils_workflow_omics_declarations.R` auto-merged; nothing under `R/`, `tests/`, `inst/`, `man/` conflicted.
  - `NEWS.md`: staging's bullets kept in order (KIT-255's, KIT-269's, KIT-273's), with the one declared edit: KIT-242 removes the sentence at the end of the KIT-255 bullet ("Where a sample's values were missing altogether in one assay ... (see KIT-242; until KIT-242's first part is merged)"); the removal is kept, so the bullet ends "... before the release, KIT-275; VAL-330). (KIT-255)".
  - Ledger: staging's rows through VAL-383 kept, then KIT-242's `VAL-<merge>` row.
- Merge 2 (`git merge --no-ff origin/kit-260`): conflict in the ledger only (one hunk: KIT-260's row added after KIT-242's). `NEWS.md` auto-merged. Both rows kept, KIT-242's first.
- Checks: `git diff origin/kit-260 HEAD -- R tests inst man NAMESPACE DESCRIPTION` is empty (staging has not moved since 91817473: the branch holds it); `git diff origin/kit-260 HEAD -- NEWS.md docs/methods-decision-ledger.md` is empty at the code-bearing merge head `b65f67706` (no line differs from the kit-260 branch's versions). After the record-fix commit the ledger differs from kit-260's only by the placeholder renames, the date cell and the added `VAL-<train37-3>` row below.

### Each member's own change

- KIT-242 (`git diff --name-only <merge base with staging> origin/kit-242`, outside `docs/qa/`): 28 files: `NEWS.md`, `R/func_lipid_norm_ruv_optimization_helpers.R`, `R/func_lipid_s4_normalization_methods.R`, `R/func_lipid_s4_ruv_control_methods.R`, `R/func_lipid_s4_ruv_methods.R`, `R/func_metab_norm.R`, `R/func_metab_s4_design_control_methods.R`, `R/func_metab_s4_ruv_methods.R`, `R/func_metab_s4_transformation_methods.R`, `R/func_omics_norm_shared.R`, `R/func_omics_ruv_factors.R`, `R/func_prot_s4_norm_methods.R`, `R/mod_omic_norm_panels.R`, `R/mod_prot_norm_workflow_helpers.R`, `R/utils_workflow_analysis_record.R`, `R/utils_workflow_dia_products.R`, `R/utils_workflow_omics_actions.R`, `R/utils_workflow_omics_declarations.R`, `R/utils_workflow_prot_declarations.R`, `R/utils_workflow_saved_exports.R`, `docs/methods-decision-ledger.md`, `tests/testthat/test-kit172-downstream.R`, `tests/testthat/test-kit242-normalisation-fallback.R`, `tests/testthat/test-metab-02w-normalise-between-samples-characterization.R`, `tests/testthat/test-module-ci-metab-norm.R`, `tests/testthat/test-workflow-lipid-declarations.R`, `tests/testthat/test-workflow-metab-declarations.R`, `tests/testthat/test-workflow-omics-worker.R`.
- KIT-260 (the train's change against staging, outside `docs/qa/`, minus KIT-242's): 16 files: `R/func_general_filemgmt_study_parameter_writer.R`, `R/func_lipid_import_detection.R`, `R/func_metab_import_detection.R`, `R/mod_omic_import_feature_table.R`, `R/utils_integration_gui.R`, `R/utils_omic_gui_config.R`, `R/utils_workflow_integration_model.R`, `inst/reports/lipidomics/lipidomics_report.rmd`, `inst/reports/metabolomics/metabolomics_report.rmd`, `tests/testthat/test-integration-gui.R`, `tests/testthat/test-kit260-itsd-default-pattern.R`, `tests/testthat/test-metab-00-import-detection-characterization.R`, `tests/testthat/test-metab-import-detection-direct-shared.R`, `tests/testthat/test-omic-norm-module.R`, `tests/testthat/test-omics-norm-core.R`, `tests/testthat/test-workflow-omics-import.R`.
- The train in all against staging 91817473, outside `docs/qa/`: 44 files; plus 52 files under `docs/qa/workflow-execution/kit/kit-242/` and 71 under `kit-260/`.

### Record fixes (one docs-only commit `f668c6cbf`; no test, no R changed)

From AUDIT-KIT-242-impl-20261126 "Closure check, 2026-11-26" and AUDIT-KIT-260-impl-20261126: (1) `kit-242/outcome.md` R1g: "the messages staging gives" replaced by "R2's head's messages, KIT-242 on staging aa119308, not staging's own"
(`kit-242/runs.tsv`'s R1g row already reads "R2 head's"; no change there); (2) `kit-242/outcome.md` "Not checked (final)": no test covers VAL-383's cleared case or the publish step's own write when a sample is removed; (3) the KIT-242 ledger row's date cell reads 2026-11-26;
(4) `kit-260/outcome.md`: the `states` and `final_pca` keys are attributed to KIT-242 Part A's per-assay ITSD arguments, not its RUV outcome; (5) the placeholders: `VAL-<train37-1>` (KIT-242's row), `VAL-<train37-2>` (KIT-260's row), and the stored-reference row `VAL-<train37-3>`
drafted from AUDIT-KIT-260-impl-20261126 MR.3 (5) with option (b) chosen (VAL-376: the sixteen keys left failing, KIT-260's values kept; counts 631/8 to 583/56; the run id is cited by pointing at this record).

## Declaration (committed and pushed before the request)

**Purpose.** The full suite and R CMD check on the merged tree, base against head: the merge check of KIT-242 and KIT-260, and the checkpoint after trains 35 and 36. Base = `918174732` (staging, the commit built from), head = `b65f67706`.

**Scope: the full suite**, `tests` "all": 538 test files on the head (536 on the base; the two new files are KIT-242's and KIT-260's). **Settings** (the form of train 35's request `kit216-r3-full-761b4301.json`): 20 shards per side,
`file_timeout_minutes` 45, `file_timeouts` 120 for the four GUI-entry files and `test-prot-norm-core.R`, 90 for `test-workflow-saved-da.R` and `test-kit266-limpa-combined-qc-figure.R`, 120 for `test-kit247-design-resave.R` (was 100);
`env` `MULTISCHOLAR_TRX_MANIFEST_NOTE=1`; `check` true on both sides; `blas` reference; `os` ubuntu-22.04; no scripts, no extra packages.
Shard planning: train 35's longest shard was 91.5 minutes; with `test-kit266-limpa-combined-qc-figure.R` (about 26 to 44 minutes) and the two new files (a few minutes) the worst case I can build from the committed tables is about 146 job-minutes, under 75 % of 300 (225).

**Ceiling: 2,800 job-minutes = 186.67 core-hours** (core-hours = job-minutes x 4 / 60). **Expected about 1,725 job-minutes = 115 core-hours** (train 35: 1,678.76, plus the new files and the kit266 file).

### Predictions

Base counts come from the newest committed tables, train 36's head table (`merge-train-20261127c-tests-head.tsv`, its 60 files) and otherwise train 35's head column (`merge-train-20261127b-compare.tsv`, the full suite at 761b4301), never from older scoped trains.
The complete list, one row per file, is `merge-train-20261126aa-predictions.tsv`. **A base-side count that differs from the recorded one is REPORTED, not explained away.**

| File | Base | Head | Verdict |
|---|---|---|---|
| test-kit242-normalisation-fallback.R | absent | 571/0 | only-on-head |
| test-kit260-itsd-default-pattern.R | absent | 157/0 | only-on-head |
| test-workflow-omics-worker.R | 289/0 | 265/0, all pass | count-changed |
| test-module-ci-metab-norm.R | 156/0 | 158/0 | count-changed |
| test-metab-02w-normalise-between-samples-characterization.R | 16/0 | 18/0 | count-changed |
| test-omics-norm-core.R | 62/0 (14 skipped) | 63/0 (14 skipped) | count-changed |
| test-omic-norm-module.R | 631/8 | 583/56 | new-failure by declaration: the 48 new failing keys are the owned ones below, and no other |
| test-prot-norm-core.R | 1555/307 | 1555/307 | message-changed: KIT-242's marker rows in the failing digests' messages; no `samples_removed` in any message |
| test-gui-entry-lipidomics.R 19/16; test-gui-entry-metabolomics.R 57/48; test-workflow-omics-declarations.R 43/25; test-workflow-omics-binding-steps.R 57/2 | as stated | as stated | message-changed |
| test-workflow-lipid-declarations.R, test-workflow-metab-declarations.R | 259/0 | 259/0 | same (the closure check corrects the earlier "fixed 256/3") |
| test-omics-module-globals.R 14/0; test-kit172-downstream.R 84/0 | as stated | as stated | same |
| every other file | the newest committed table | same | same (the known one-sided `stringdb_plot` digest flake in `test-omic-enrich-module.R` is recognised as in train 35) |
| R CMD check | 0 ERROR, 3 WARNING, 3 NOTE | identical | nothing only on head |

Totals expected: 538 files; 2 only-on-head, 4 count-changed, 1 new-failure, 5 message-changed, the rest same; no `fixed`; no other new-failure.

**The 48 owned failing keys of `test-omic-norm-module.R`** (kit-260/outcome.md; the 56 failing `omicNormFrozen` keys on the head, 8 failing on the base): the sixteen `log` and `outputs` keys (4 `log`, 12 `outputs`; lipid and metab manual and automatic, `run` and the sessions listed there) are VAL-376's, KIT-260's values kept and left failing until the re-record;
the `provenance` keys are VAL-303 (1)'s; the `states` and `final_pca` keys move with KIT-242 Part A's per-assay ITSD arguments (VAL-352). By scenario: lipid and metab manual: run `log`, `states`, `provenance`, `outputs`, `final_pca`; skip_correlation `states`, `provenance`, `outputs`, `final_pca`; reset `states`, `provenance`, `outputs`; lipid and metab automatic: run `log`, `states`, `provenance`, `outputs`, `final_pca`; `final_plot` and `final_plot_failing` `states`, `provenance`, `outputs`; lipid and metab skip: `provenance` of run, correlation, final_plot, export, reset.
The key-by-key check: the set of failing keys on the head, from `failures.txt`, minus the 8 base failures, equals this list; a key outside it is a result.

**What was never run together.** Test files that staging gained after KIT-242's branch point (from KIT-266, KIT-273 and KIT-216) have never run with KIT-242's code; this run is their first. Any difference outside this list is a RESULT: reported, not repaired.
This run is also the checkpoint after trains 35 and 36.

### Decision and attempts

One run. A shard that times out or an instrument failure with no result: handed back, not re-sent. An HTTP 500/502 at dispatch with no run created is started once by hand.

## Results (run 45801695192, request train37-kit242-kit260-b65f6770, head b65f67706, base 918174732)

The first `gh run watch` ended early on an HTTP 502 from the annotations endpoint (not a run failure); I re-attached once, as allowed. Run success; all 45 jobs ran, every file ran to completion (no timeout, no crash, no withheld file).
Tables beside this file: `merge-train-20261126aa-compare.tsv`, `-check-compare.txt`, `-tests-head.tsv`, `-tests-base.tsv`, `-omic-norm-module-new-keys.txt` (the 48 head-only failing keys, one per line).
Base versus head, 538 files: same 526, only-on-head 2, count-changed 4, message-changed 5, new-failure 1, fixed 0. Tests head: passed 56993, failed 1250, errors 13, skipped 586; base: 56332, 1202, 13, 586.
**Every prediction held: TRAIN GREEN for KIT-242 and KIT-260 on the merged tree, and the checkpoint after trains 35 and 36 shows no unpredicted file.**

| File | Predicted | Observed base | Observed head | Verdict |
|---|---|---|---|---|
| test-kit242-normalisation-fallback.R | absent / 571/0 | absent | 571/0/0/0 | only-on-head: HELD |
| test-kit260-itsd-default-pattern.R | absent / 157/0 | absent | 157/0/0/0 | only-on-head: HELD |
| test-workflow-omics-worker.R | 289/0 / 265/0 | 289/0 | 265/0 | count-changed: HELD |
| test-module-ci-metab-norm.R | 156/0 / 158/0 | 156/0 | 158/0 | count-changed: HELD |
| test-metab-02w-normalise-between-samples-characterization.R | 16/0 / 18/0 | 16/0 | 18/0 | count-changed: HELD |
| test-omics-norm-core.R | 62/0 (14 skipped) / 63/0 (14) | 62/0/0/14 | 63/0/0/14 | count-changed: HELD |
| test-omic-norm-module.R | 631/8 / 583/56 | 631/8 | 583/56 | new-failure by declaration: HELD (key-by-key below) |
| test-prot-norm-core.R | 1555/307 both, message-changed | 1555/307 | 1555/307 | message-changed: HELD; no `samples_removed` in any of its messages (0 occurrences) |
| test-gui-entry-lipidomics.R 19/16; test-gui-entry-metabolomics.R 57/48; test-workflow-omics-declarations.R 43/25; test-workflow-omics-binding-steps.R 57/2 | counts unchanged, message-changed | as stated | as stated | message-changed: HELD |
| the two declarations files 259/0; test-omics-module-globals.R 14/0; test-kit172-downstream.R 84/0 | same | as stated | as stated | same: HELD |
| every other file (526 same) | same | equal to the predictions table | same | HELD; 0 base-side counts differed from the recorded ones (the predictions table, all 536 base files); the `stringdb_plot` flake did not appear (0 mentions in `compare-failures.txt`) |
| R CMD check | identical 0/3/3 | 0 E / 3 W / 3 N | 0 E / 3 W / 3 N | identical, 6 findings on both, only on head none, only on base none: HELD |

**Key-by-key check of `test-omic-norm-module.R`** (Python `-I`; each failure record's last message line names the key, e.g. "lipid skip run provenance"). 56 failing keys on the head, 8 on the base; the 8 are the base's (all `skip ... provenance`, identical keys). The 48 head-only keys:
log 4, outputs 12 (the sixteen of VAL-376), states 12 and final_pca 6 (the eighteen of KIT-242 Part A, VAL-352), provenance 14 (VAL-303 (1); with the 8 on both sides the 22 declared). By scenario: manual lipid and metab 12 each (log, states x3, provenance x3, outputs x3, final_pca x2), automatic 11 each (log, states x3, provenance x3, outputs x3, final_pca x1), skip 5 each (provenance of run, correlation, final_plot, export, reset; 4 already failing on the base, the `run` one new).
That is exactly the declared list, and no key outside it.

**Budget** (jobs API, summed, skipped jobs 0): **1,750.15 job-minutes = 116.68 core-hours** (core-hours = job-minutes x 4 / 60); ceiling 2,800 (186.67), expected about 1,725 (115): 1.5 % above the expectation, 62.5 % of the ceiling.
Tests shards 16.0 to 108.0 minutes (longest `tests head 4`, 48 % of 300, under the 75 % bound); check head 15.32, base 15.58; plan 0.43; results 0.68. Heavy files: `test-kit247-design-resave.R` 4,119.5 s (68.7 min, limit 120), `test-kit266-limpa-combined-qc-figure.R` 1,714.9 s (limit 90), `test-prot-norm-core.R` 1,899.4 s, `test-gui-entry-proteomics.R` 1,114 s.

Decision: KIT-242 and KIT-260 are green on the merged tree; the staging test files from KIT-266, KIT-273 and KIT-216 ran with KIT-242's code for the first time and are `same`.

## What the coordinator resolves at landing

- Land in this order; the ledger rows are KIT-242's `VAL-<train37-1>` (`docs/methods-decision-ledger.md:2556`), KIT-260's `VAL-<train37-2>` (`:2557`) and the stored-reference row `VAL-<train37-3>` (`:2558`); next free number VAL-384 (so 384, 385, 386 if staging has not moved).
  The rows cite VAL-242, VAL-290, VAL-352, VAL-373, VAL-383 (KIT-242) and VAL-241, VAL-275, VAL-303, VAL-352 (4), VAL-376, VAL-181, VAL-189, VAL-259 (KIT-260 and the stored-reference row); fill their status cells where they wait for this merge ("implemented by KIT-242/KIT-260 (VAL-n), merged by train 37").
- The stored-reference row cites this record for the run id (45801695192): replace the pointer by the id if wished.
- Placeholders in the members' records: see the hand-back (grep of `VAL-<` in `kit-242/` and `kit-260/` and `NEWS.md`).
