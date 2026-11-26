# Merge train 29 of 2026-11-26: KIT-247

**This train is NOT stacked.** Branched from `origin/v0.9-staging` at `7bbc344eda21b6bb4a0ab68b81a62e38cbbfdf88` (fetched
before building; it is the newest staging head: no merge train has landed since). It does not contain the merge trains
still out (r, s, t, u and any other agent's), so KIT-255, KIT-271 and KIT-278 are NOT on this base.

## Train

Base: `7bbc344eda21b6bb4a0ab68b81a62e38cbbfdf88` (staging).

| Order | Ticket | Branch head | Merge commit | Ledger row | Audit |
|---|---|---|---|---|---|
| 1 | KIT-247 | `origin/kit-247` a52a808d (contains staging 12701142, not the newer merges) | 3f7826d1 | `VAL-366` (replaces `VAL-<merge>`) | AUDIT-KIT-247-impl-20261126, closure check 2 (ready once the D.4 record fixes land) |

Head (bound by the run): `3f7826d122d0a4c9b82324588c1455f1f4783839` (the merge commit, the last code-bearing commit). The
record-fix commit and this declaration come after it and change `docs/` only.

### Conflicts met and their resolutions

- **`NEWS.md`**: staging's bullets first (the KIT-249 bullet), then KIT-247's one bullet; nothing removed.
- **`docs/methods-decision-ledger.md`**: staging's rows kept in order (to VAL-347); KIT-247's row directly after VAL-347 as
  `VAL-366` (the coordinator numbers it at landing); its status cells still read "<merge date>" and "merged into
  `v0.9-staging` <date> by merge train <n>".
- No conflict under `R/`, `tests/`, `inst/`, `man/`, `DESCRIPTION` or `NAMESPACE`; `R/mod_omic_norm.R`,
  `R/utils_workflow_dia_products.R` and `R/utils_workflow_summary.R` merged automatically.

### The merged code differs from the base by exactly KIT-247's change

`git diff --stat 7bbc344e 3f7826d1 -- R tests inst man DESCRIPTION NAMESPACE` gives 13 files, 1,344 insertions and 33
deletions, the same list and the same numbers as `git diff --stat 12701142 origin/kit-247` over those paths (KIT-247's own
change against its last staging merge); the added and removed lines of every file are identical (compared file by file).
Nothing under `inst/`, `man/`, `DESCRIPTION`, `NAMESPACE` or `tests/testdata` differs.

- `R/`: `func_omics_design_shared.R`, `mod_omic_design.R`, `mod_omic_norm.R`, `mod_omic_qc.R`,
  `mod_prot_design_builder_helpers.R`, `mod_prot_design_import_helpers.R`, `utils_workflow_dia_products.R`,
  `utils_workflow_omics_design.R`, `utils_workflow_summary.R`, `utils_workflow_trx_norm.R`.
- `tests/testthat/`: `test-kit247-design-resave.R`, `test-kit247-import-config.R`, `test-kit247-resave-no-session.R`.
- `omicsDesignCanonicalText` in the merged tree (`git grep`): the definition and its two call sites, all in
  `R/utils_workflow_omics_design.R`:
  - line 337: `design[] <- lapply(design, omicsDesignCanonicalText)`
  - line 349: `omicsDesignCanonicalText <- function(column) {`
  - line 368: `table[] <- lapply(table, omicsDesignCanonicalText)`

Lines added by the merge in `R/` contain no `<<-`, `.GlobalEnv` or `globalenv()`; the `assign(` and `get0(` lines are in
`R/utils_workflow_omics_design.R` (`session$userData` keys), a file `test-omics-module-globals.R` does not scan (it scans
`R/mod_omic*` and `R/func_omics*`; the four of those KIT-247 changes have none).

## Declaration (written before the run is sent)

**Purpose.** The full regression comparison the ticket never ran on its final head (KIT-247's last full comparison,
45698414166, was on the second-correction head; the withdrawal run 45722093901 ran only the three KIT-247 files): merged
tree against staging 7bbc344e, per file, plus R CMD check on both sides.

**Scope recount, on the merged tree.** Each group counted with `git grep` and `ls` on 3f7826d1:

- the ticket's R2 grep over `tests/testthat` (helpers excluded; the register's widened form with `save_workflow_args`,
  `saveWorkflowSummaryParameters`, `runProtDesignImportObserverShell`) plus the ticket's named additions (kit-212, kit-222,
  kit-225, `test-prot-norm-input-current.R`, the four GUI-entry lanes, the three KIT-247 files) and the audit's
  helper-mediated and later additions: the **76 of `kit-247/r2pp-files.txt`** (all present; 69 from the grep and its named
  additions, 7 helper-mediated or added: `test-kit244-limpa-pass-through.R`, `test-kit256-tech-rep-selection-cleared.R`,
  `test-module-ci-lipid-norm.R`, `test-module-ci-metab-norm.R`, `test-module-ci-trx-design.R`, `test-module-ci-trx-qc.R`,
  `test-prot-norm-input-current.R`);
- **grep hits new on the merged tree (5):** `test-kit213-study-parameters-labels.R`, `test-kit254-omics-correlation-skip.R`,
  `test-kit259b-second-correction.R`, `test-kit263-prot-study-parameters-sections.R`, `test-kit264-trx-export-current.R`;
- **partners on shared files (4):** `test-kit257-norm-own-files.R`, `test-kit259-pearson-histogram-figure-failure.R`,
  `test-kit248-correlation-skip-record.R`, `test-kit261-norm-run-button.R`;
- **globals scan (1):** `test-omics-module-globals.R`; **helper-mediated (1):** `test-omics-cyclic-loess-minimum.R`;
- `test-kit255-omics-correlation-assays.R` and `test-kit271-limpa-table-name.R`: NOT in scope (KIT-255 and KIT-271 are not
  on this base: the files do not exist in the tree).

76 + 5 + 4 + 1 + 1 = **87 files**. The register's closure check 2 and its D.4 text for `outcome.md` say 88; its own list
of added files sums to 87 with the 76. No file is missing from the list of 87 that the register names, so the 88 is read as
an addition slip; the count is recorded here as recounted and `kit-247/outcome.md` carries D.4's text, "88 files", unedited
(the coordinator to read it as 87 at landing). `test-prot-norm-core.R` is in scope (it is in the 76), so
`MULTISCHOLAR_KIT105_W3_DUMP` is not set (the last request set no env).

**Run settings.** Copied from KIT-247's last full request, `kit247-final-166ef4f5.json` (run 45698414166): `env` `{}`,
`extra_packages` `[]`, `check` true, `blas` `reference`, `scripts` `[]`, `fetch_refs` `[]`, `tools_ref` `""`,
`file_timeout_minutes` 45, `file_timeouts` `{"test-kit247-design-resave.R": 100, "test-gui-entry-proteomics.R": 60}`.
Changed: the 87 files, **6 shards per side**, base = `7bbc344e`, head = `3f7826d1`. One run.

**Shards and the 75-minute arithmetic.** Each file's time is the longest recorded in the newest committed tables (KIT-247's
final run for the 76, the latest run on disk for the added files; `test-kit213-study-parameters-labels.R` and
`test-kit264-trx-export-current.R` have no passing timing on disk and are taken at 300 s each, `test-kit263` at 15 s). The
files sum to 343.6 minutes on head (the base side lacks the three KIT-247 files, 56.3 + 0.2 minutes: 287.1). Files are dealt
round-robin in the listed order, so the list is ordered to give the partition below (files minutes per shard, head side):
shard 1 = 57.9 (the 56.3-minute design-resave file with 14 files of 6 to 7 seconds), shards 2 to 6 = 57.3, 57.2, 57.1,
57.1, 57.1. A job's setup (library install, native compile) ran 533 to 705 s in train 22, so about 12 minutes: the longest
job is about 57.9 + 12 = **70 minutes**, under 75. The single longest file keeps a 100-minute timeout and the
GUI-entry proteomics file 60; a file killed at its limit is a result, not a shard overrun at the committed timings.

**Ceiling.** 1,000 job-minutes (66.67 core-hours; core-hours = job-minutes x 4 / 60). Expected about 750: the final run's
76 files with 4 + 4 shards took 615.4; the 11 added files add about 83 (41.5 per side), the 4 extra jobs about 45 of setup.
**Stop rule:** any `new-failure` or changed message outside the three KIT-247 files stops the train for this member.

### Predictions

Counts are passed/failed/errors/skipped. "Base count" is the register's figure (closure check 2, D.2 (7); counts at
staging 12701142 from run 45698414166's base side, and the added files from the latest runs on disk), checked against the
newest committed comparison tables. **Where the two differ, the train's own base run is the reference** and the verdict
prediction (`same`, identical failure messages) is what is judged; that holds above all for files not run since their last
change on staging and the files whose newest table is another train's head (stacked on older staging heads). Differences
found: `test-gui-entry-lipidomics.R`, `test-gui-entry-metabolomics.R`, `test-gui-entry-proteomics.R`, `test-integration-input-transcriptomics.R`, `test-kit129-run-path.R`, `test-kit224-ruv-feature-levels.R`, `test-module-ci-trx-qc.R`, `test-module-ci-trx-summary.R`, `test-omics-norm-state-guards.R`, `test-trx-da-without-ruv.R` (the note column gives the newest table's count). The three KIT-247 files are the only
ones predicted to move: absent on base, passing on head. `test-omics-module-globals.R` is predicted at the base's own
count, 13/1 (`shown_idle <<- idle`, `R/mod_omic_norm.R:282`, not yet allowed: KIT-278 is not on this base). R CMD check:
findings identical on both sides, nothing only on head.

| Shard | File | Base count (register; passed/failed/errors/skipped) | Head | Verdict | Note |
|---|---|---|---|---|---|
| 1 | test-kit247-design-resave.R | absent | 188/0/0/0 | only-on-head | fixed/new on head only |
| 1 | test-kit247-import-config.R | absent | 2/0/0/0 | only-on-head | fixed/new on head only |
| 1 | test-kit247-resave-no-session.R | absent | 7/0/0/0 | only-on-head | fixed/new on head only |
| 1 | test-lipid-12-norm-module-contracts.R | 0/0/0/0 | same as base | same |  |
| 1 | test-lipid-14a-summary-module-characterization.R | 0/0/0/0 | same as base | same |  |
| 1 | test-lipid-15-summary-module-seams.R | 0/0/0/0 | same as base | same |  |
| 1 | test-metab-03b-norm-module-helper-characterization.R | 0/0/0/0 | same as base | same |  |
| 1 | test-metab-03c-norm-module-seams.R | 0/0/0/0 | same as base | same |  |
| 1 | test-metab-05-summary-module-characterization.R | 0/0/0/0 | same as base | same |  |
| 1 | test-prot-04a-design-module-characterization.R | 109/0/0/18 | same as base | same |  |
| 1 | test-prot-12-summary-module-contracts.R | 0/0/0/0 | same as base | same |  |
| 1 | test-prot-design-save-preflight.R | 8/0/0/0 | same as base | same |  |
| 1 | test-prot-qc-coordinator-direct-shared.R | 4/0/0/3 | same as base | same |  |
| 1 | test-workflow-memory-mechanism-member.R | 0/0/0/13 | same as base | same |  |
| 1 | test-workflow-neutral-names.R | 39/0/0/0 | same as base | same |  |
| 2 | test-gui-entry-transcriptomics.R | 16/13/0/0 | same as base | same |  |
| 2 | test-integration-input-lipid-metab.R | 208/0/0/0 | same as base | same |  |
| 2 | test-kit244-limpa-pass-through.R | 204/0/0/0 | same as base | same |  |
| 2 | test-kit251-replicate-filter-one-click.R | 9/0/0/0 | same as base | same |  |
| 2 | test-kit264-trx-export-current.R | 85/0/0/0 | same as base | same | no committed table on staging; register's figure from the run named in the register |
| 2 | test-metab-03d-norm-module-characterization.R | 13/0/0/3 | same as base | same |  |
| 2 | test-module-ci-trx-qc.R | 34/1/1/0 | same as base | same | newest committed table 38/0/1/0 (kit-213/r3/compare.tsv); the train's own base run is the reference |
| 2 | test-omics-cyclic-loess-minimum.R | 72/0/0/0 | same as base | same |  |
| 2 | test-omics-module-globals.R | 13/1/0/0 | same as base | same | KIT-278 (14/0) is not on this base; KIT-247 adds no `<<-` or global access in the scanned files |
| 2 | test-prot-05b-norm-module-contracts.R | 440/0/0/22 | same as base | same |  |
| 2 | test-prot-annotation-uniprot-empty.R | 327/0/0/0 | same as base | same |  |
| 2 | test-prot-norm-core.R | 1547/307/0/0 | same as base | same | lead's figure 1547/307; register 1548/306; newest table 1547/307 |
| 2 | test-trx-norm-module.R | 109/0/0/0 | same as base | same |  |
| 2 | test-workflow-checkpoint-dia-state.R | 174/0/0/4 | same as base | same |  |
| 2 | test-workflow-native-session.R | 734/0/0/2 | same as base | same |  |
| 3 | test-gui-entry-metabolomics.R | 105/0/0/0 | same as base | same | newest committed table 57/48/0/0 (merge-train-20261126p/compare.tsv); the train's own base run is the reference |
| 3 | test-integration-input-transcriptomics.R | 67/3/0/0 | same as base | same | newest committed table 79/0/0/0 (merge-train-20261126q-compare.tsv); the train's own base run is the reference |
| 3 | test-kit222-replicate-keys.R | 456/0/0/0 | same as base | same |  |
| 3 | test-kit254-omics-correlation-skip.R | 379/0/0/0 | same as base | same |  |
| 3 | test-kit256-tech-rep-selection-cleared.R | 20/0/0/0 | same as base | same |  |
| 3 | test-kit257-norm-own-files.R | 133/0/0/0 | same as base | same |  |
| 3 | test-lipid-01-qc-filtering-helpers.R | 58/0/0/2 | same as base | same |  |
| 3 | test-module-ci-lipid-design.R | 109/0/0/0 | same as base | same |  |
| 3 | test-module-ci-trx-summary.R | 83/2/0/0 | same as base | same | newest committed table 87/1/0/0 (merge-train-20261126q-compare.tsv); the train's own base run is the reference |
| 3 | test-omic-design-module.R | 160/0/0/0 | same as base | same |  |
| 3 | test-omic-norm-module.R | 639/0/0/0 | same as base | same |  |
| 3 | test-omics-norm-itsd-selection.R | 22/0/0/0 | same as base | same |  |
| 3 | test-prot-norm-input-current.R | 321/0/0/0 | same as base | same |  |
| 3 | test-workflow-figure-paths.R | 181/0/0/0 | same as base | same |  |
| 3 | test-workflow-integration-input.R | 241/0/0/0 | same as base | same |  |
| 4 | test-kit129-run-path.R | 687/6/0/1 | same as base | same | newest committed table 711/0/0/1 (kit-213/r3/compare.tsv); the train's own base run is the reference |
| 4 | test-kit172-downstream.R | 84/0/0/0 | same as base | same |  |
| 4 | test-kit206-prot-norm-export-guard.R | 116/0/0/0 | same as base | same |  |
| 4 | test-kit224-ruv-feature-levels.R | 65/1/0/0 | same as base | same | newest committed table 69/0/0/0 (merge-train-20261126q-compare.tsv); the train's own base run is the reference |
| 4 | test-kit248-correlation-skip-record.R | 206/0/0/0 | same as base | same |  |
| 4 | test-kit259b-second-correction.R | 90/0/0/0 | same as base | same |  |
| 4 | test-kit263-prot-study-parameters-sections.R | 11/0/0/0 | same as base | same |  |
| 4 | test-lipid-12a-norm-module-characterization.R | 5/0/0/5 | same as base | same |  |
| 4 | test-module-ci-metab-qc.R | 63/0/0/1 | same as base | same |  |
| 4 | test-omic-norm-ruv-design.R | 1153/0/0/0 | same as base | same |  |
| 4 | test-omics-summary-core.R | 676/25/0/2 | same as base | same |  |
| 4 | test-prot-design-core.R | 37/0/0/0 | same as base | same |  |
| 4 | test-trx-norm-core.R | 58/0/0/0 | same as base | same |  |
| 4 | test-workflow-omics-application-reads.R | 230/0/0/0 | same as base | same |  |
| 5 | test-gui-entry-proteomics.R | 168/25/0/0 | same as base | same | newest committed table 92/101/0/0 (merge-train-20261126p/compare.tsv); the train's own base run is the reference |
| 5 | test-kit207-summary-head-guard.R | 292/0/0/0 | same as base | same |  |
| 5 | test-kit212-da-current.R | 309/0/0/0 | same as base | same |  |
| 5 | test-kit225-prot-norm-reset-source.R | 56/0/0/0 | same as base | same |  |
| 5 | test-kit261-norm-run-button.R | 105/0/0/0 | same as base | same |  |
| 5 | test-module-ci-lipid-norm.R | 156/0/0/0 | same as base | same |  |
| 5 | test-module-ci-lipid-qc.R | 68/0/0/0 | same as base | same |  |
| 5 | test-module-ci-metab-design.R | 96/0/0/0 | same as base | same |  |
| 5 | test-module-ci-metab-norm.R | 156/0/0/0 | same as base | same |  |
| 5 | test-omics-design-core.R | 64/0/0/10 | same as base | same |  |
| 5 | test-prot-04-design.R | 488/0/0/118 | same as base | same |  |
| 5 | test-trx-da-without-ruv.R | 95/2/0/0 | same as base | same | newest committed table 103/0/0/0 (kit-213/r3/compare.tsv); the train's own base run is the reference |
| 5 | test-workflow-analysis-record.R | 351/0/0/0 | same as base | same |  |
| 5 | test-workflow-checkpoint-harness.R | 628/0/1/8 | same as base | same |  |
| 6 | test-gui-entry-lipidomics.R | 35/0/0/0 | same as base | same | newest committed table 19/16/0/0 (merge-train-20261126p/compare.tsv); the train's own base run is the reference |
| 6 | test-kit213-study-parameters-labels.R | 56/0/0/0 | same as base | same |  |
| 6 | test-kit243-summary-current-outputs.R | 246/0/0/0 | same as base | same |  |
| 6 | test-kit259-pearson-histogram-figure-failure.R | 85/0/0/0 | same as base | same |  |
| 6 | test-module-ci-trx-design.R | 66/0/1/0 | same as base | same |  |
| 6 | test-module-ui-wrapper-direct-shared.R | 45/0/0/0 | same as base | same |  |
| 6 | test-omic-gui-config.R | 1911/0/0/0 | same as base | same |  |
| 6 | test-omic-qc-module.R | 160/0/0/0 | same as base | same |  |
| 6 | test-omic-summary-module.R | 1675/28/0/0 | same as base | same |  |
| 6 | test-omics-design-import-refusal.R | 39/0/0/0 | same as base | same |  |
| 6 | test-omics-norm-core.R | 62/0/0/14 | same as base | same |  |
| 6 | test-omics-norm-state-guards.R | 647/0/0/0 | same as base | same | newest committed table 651/0/0/0 (merge-train-20261126m-compare.tsv); the train's own base run is the reference |
| 6 | test-prot-05d-norm-module-characterization.R | 82/0/0/10 | same as base | same |  |
| 6 | test-prot-design-log2.R | 226/0/0/0 | same as base | same |  |

| R CMD check | | | identical findings, nothing only on head | |

### Decision and attempts

- Every prediction holds: TRAIN GREEN for KIT-247.
- Any other outcome stops the train for this member: recorded with its files and counts, not repaired, not re-sent.
- One run; a re-send only for an instrument defect with no result judged. An HTTP 500/502 at dispatch with no run
  created: started once by hand with the README's command, recorded.

## What the coordinator resolves at landing

- The ledger: `VAL-366` takes its number (next free after the rows landed before it); its status cells ("<merge
  date>", "merged into `v0.9-staging` <date> by merge train <n>", and the train's run id) are filled; VAL-239's, VAL-277's
  and VAL-331's status cells as the register names them.
- `NEWS.md`: KIT-247's bullet in landing order.
- `kit-247/outcome.md`: names `VAL-<merge>` in its text and carries "88 files" (D.4's text; the recount is 87).
- Records owed at landing: plan and DISPATCH rows, I08 (the browser observation) before the release tag, the KIT-216 order
  (the second of KIT-216 and KIT-247 to merge runs the other's new and edited test files).

## Results (run 45734833773, request `merge-train-20261126v-3f7826d1`)

Run: success; base `7bbc344e`, head `3f7826d1`; tables `merge-train-20261126v-compare.tsv`, `-tests.tsv`, `-tests-base.tsv`,
`-check-compare.txt` in this directory. **Verdict: TRAIN GREEN for KIT-247; every verdict prediction held.**

- **Base versus head:** 87 files; same 84, only-on-head 3, new-failure 0, message-changed 0, count-changed 0, fixed 0,
  compare-error 0; `compare-failures.txt`: "no non-identical failures". Head: passed 19,024, failed 558, errors 3, skipped
  236, no timeout, no crash; base: passed 18,827, failed 558, errors 3, skipped 236 (the 197 passed more on head are the
  three KIT-247 files: 188 + 2 + 7).
- **The three KIT-247 files, only on head:** `test-kit247-design-resave.R` 188/0/0/0 ok, `test-kit247-import-config.R`
  2/0/0/0 ok, `test-kit247-resave-no-session.R` 7/0/0/0 ok, as predicted.
- `test-omics-module-globals.R` 13/1/0/0 on both sides (predicted); `test-prot-norm-core.R` 1547/307/0/0 on both sides
  (predicted).
- **R CMD check:** base and head each ERRORs 0, WARNINGs 3, NOTEs 3 (exit 0); findings identical on both, 6; only on head 0.

### Prediction errors of the declaration's absolute base counts (verdicts unaffected)

The verdict `same` held for all 84 files; the stop rule concerns new failures or changed messages against the run's own base
side. The declaration's absolute base counts missed for five files (the declaration had said the train's own base run is the
reference for them); actual counts, the same on base and head:

| File | Declared base | Actual base = head |
|---|---|---|
| test-gui-entry-lipidomics.R | 35/0/0/0 | 19/16/0/0 |
| test-gui-entry-metabolomics.R | 105/0/0/0 | 57/48/0/0 |
| test-gui-entry-proteomics.R | 168/25/0/0 | 92/101/0/0 |
| test-omic-qc-module.R | 160/0/0/0 | 156/4/0/0 |
| test-omics-norm-state-guards.R | 647/0/0/0 | 651/0/0/0 |

The coordinator's note of 2026-11-26 (quoted): "it lists stale base counts for three files: test-gui-entry-lipidomics.R 35/0,
test-gui-entry-metabolomics.R 105/0 and test-omic-qc-module.R 160/0; the committed values on staging are 19/16, 57/48 and
156/4 (these files already fail on stored values since KIT-217). The verdict `same` is unaffected." The run was already out
when it arrived, so the declaration was not changed. Two more files in the declaration's table also carried a register
count that the newest committed table had already moved: `test-gui-entry-proteomics.R` (newest table 92/101, named in the
note column) and `test-omics-norm-state-guards.R` (651/0, from KIT-264's merge-train m table, named in the note column).

### Job-minutes and core-hours

Job minutes (from the jobs' own start and end times): plan 0.4; check base 15.8, check head 16.3; tests head 1 to 6 =
79.0, 58.9, 39.8, 59.9, 70.5, 54.0; tests base 1 to 6 = 11.2, 72.1, 65.8, 45.4, 58.6, 62.7; results 0.5.
**Total 710.7 job-minutes = 47.38 core-hours** (ceiling 1,000 job-minutes; expected about 750).

**Declaration slip: the 75-minute job bound.** At the committed timings no job was to pass 75 minutes (longest estimated
about 70). Head shard 1 took **79.0 minutes**: `test-kit247-design-resave.R` ran 3,999.4 s (66.7 minutes) against the
3,378.7 s (56.3 minutes) of the earlier committed run, on an AMD EPYC 7763 runner (setup 638 s). The file stayed inside its
100-minute timeout; no verdict depends on it. The arithmetic was correct for the committed timings and wrong for the file's
real spread (+18 percent).

### What the coordinator resolves at landing (unchanged from the declaration)

`VAL-366` number and its status cells ("<merge date>", "merged into `v0.9-staging` <date> by merge train <n>", run
id 45734833773), VAL-239's, VAL-277's and VAL-331's status cells as the register names them, `NEWS.md` order,
`kit-247/outcome.md` ("88 files" per D.4 against the recount of 87; `VAL-<merge>` in its text, line 90), plan and DISPATCH rows.

## Merge record (coordinator)

Landed on v0.9-staging on 2026-11-26 by a no-fast-forward merge of this branch at 70ee0e65 onto 8cf4d758 (the train was built on 7bbc344e; since then staging gained merge trains 26, 27, 25, 28 and 30 (KIT-271; KIT-270, KIT-246c, KIT-221; KIT-255; KIT-278, KIT-269; the KIT-258 texts), KIT-211's documents and ledger rows).

- Verdict: green. Run 45734833773 (87 files and R CMD check; 710.7 of 1,000 job-minutes, 47.38 core-hours), the full comparison KIT-247 never ran on its final head: the three KIT-247 files only on head and passing (188/0, 2/0, 7/0); 84 files same; no new failure, no changed message; check identical.
- Errors of the declaration, recorded with the results: five absolute base counts (the three GUI-entry files, test-omic-qc-module.R and test-omics-norm-state-guards.R were declared at counts older than staging's; the counts are the same on both sides, verdict same) and the bound "no job above 75 minutes" (one head shard ran 79.0 minutes; test-kit247-design-resave.R took 3,999 s against 3,379 s committed; inside its 100-minute file timeout). The scope is 87 files; kit-247/outcome.md's "88 files" is the audit's text and reads 87.
- Conflicts at landing: the ledger only. R/mod_omic_norm.R and R/utils_workflow_dia_products.R merged by themselves: in the first, staging's side (KIT-255) adds the correlation notice near lines 598 to 637 and this train changes the Reset near 524 to 571, different functions, read at landing; in the second, staging changed line 1981 and this train lines 983 to 986.
- Ledger at landing: staging's rows kept through VAL-365; KIT-247 = VAL-366. Status cells of VAL-239, VAL-277, VAL-309 and VAL-331 updated. Placeholders in this record and in kit-247/outcome.md replaced.
- Not covered by the train run: the combination with merge trains 25 to 28 and 30 (test-kit255-omics-correlation-assays.R and test-kit271-limpa-table-name.R were not on the train's base). The next checkpoint covers it.
- Owed: the browser observation of the confirmation (audit I08) before the release tag; KIT-216, merging after this, runs KIT-247's three test files and the two files both tickets edit (test-kit206-prot-norm-export-guard.R, test-workflow-analysis-record.R).
