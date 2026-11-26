# Merge train 36 of 2026-11-27: KIT-266 (a combined QC figure on the proteomics DIA limpa route), STACKED on train 35

**Stacked on train 35 (761b4301).** Train 35 (KIT-216, a full-suite run) is out and not landed; if it fails this train is rebuilt (accepted).
Branch `merge-train-20261127c`, worktree `/home/doktersmol/Projects/APAF/MultiScholaR-train36`, base `761b43017` (the merge head of train 35, "Merge KIT-216 into merge train 35"),
member `origin/kit-266` at tip `30810678` taken as it is (it holds two `-s ours` merge commits that change no tree).

## Train

| Order | Ticket | Branch head | Merge commit | Ledger rows | Audit |
|---|---|---|---|---|---|
| 1 | KIT-266 | `origin/kit-266` 30810678 | fa9a7b8c4 | the member's row `VAL-<merge>` (line 2550) and a stored-reference row, text drafted in AUDIT-KIT-266-impl-20261127 section 8 (not added here); the coordinator numbers both at landing | AUDIT-KIT-266-impl-20261127 (ready after record fixes and one ruling) |

### Conflicts

`git merge --no-ff origin/kit-266` conflicted in `docs/methods-decision-ledger.md` only (one hunk at the end of the VAL rows). `NEWS.md` and `R/utils_workflow_dia_products.R` merged
without conflict. Nothing under `R/`, `tests/`, `inst/`, `man/` conflicted. Ledger: every row of both sides kept, the base's first. **The two `VAL-<merge>` rows stay as they are**
(coordinator's ruling): the first (line 2549) is KIT-216's (train 35), the second (line 2550) is KIT-266's; the coordinator tells them apart by ticket at landing.

### The merged tree against the base

`git diff --name-only 761b43017 HEAD` outside `docs/qa/` equals `git diff --name-only 36e3180b origin/kit-266` (36e3180b = merge base with staging) outside `docs/qa/`, file for file, 12 files:
`NEWS.md`, `docs/methods-decision-ledger.md`, `R/mod_prot_norm_workflow_helpers.R`, `R/utils_workflow_dia_norm.R`, `R/utils_workflow_dia_products.R`, `inst/reports/proteomics/DIANN_limpa_report.rmd`,
`tests/testthat/test-kit227-report-ruv-sentences.R`, `test-kit259b-second-correction.R`, `test-kit266-limpa-combined-qc-figure.R` (new), `test-kit272-prot-box-panel-text.R`, `test-prot-norm-core.R`, `test-prot-norm-worker-figures.R`.
Under `R/`, `inst/`, `tests/`: 10 files, 1459 insertions, 71 deletions, the same figures on both sides. Plus 31 files under `docs/qa/workflow-execution/kit/kit-266/`; no other `docs/qa/` file differs.
Files both sides changed: `NEWS.md`, the ledger and `R/utils_workflow_dia_products.R` (KIT-216 and KIT-273 on the base side); changed-line counts of the member's change are equal on the merged tree
(75/75, 16/16, 1/1), so the base's hunks are untouched and the member's are whole.

## Declaration (committed and pushed before the request)

**Purpose.** The merge check of KIT-266 on the stacked tree: base 761b43017 (staging a855e1f0 + train 35) against head (base + KIT-266), R CMD check on both sides. One run.

**Scope: 60 files** (list in `merge-train-20261127c-scope.txt`). Derivation rule: on the merged tree, `git grep -lE` of the ticket's 18 names
(`composite_QC_figure|generateProtNormCompositeQcFigure|resolveProtNormCompositeFigureInputs|writeDIACompositeFigureFiles|saveWorkflowDiaCompositeFigure|publishDIAReferenceRuvCovariates|runDIAReferenceRuvCovariates|runDIAReferenceLimpaNormalisationPassThrough|dispatchProtNormPassThrough|dispatchProtNormRuvCovariates|runProtNormRuvCovariatesObserver|protNormW1LimpaOwner|DIANN_limpa_report|collectDIANativeGuiFiles|ruv_covariates_cancor|removeProtNormOwnFiles|runProtNormNormalizationWorkflow|copyToResultsSummary`)
plus the 10 names of the worker branch
(`normalization-qc-figures|computeWorkflowDiaQcFiguresAction|diaQcFiguresParameters|protNormFiguresAction|protNormW3Managed|protNormAuditManaged|workflowOmicsManagedActions|dispatchProtNormCorrelation|runProtNormApplyCorrelationObserver|runProtNormSkipCorrelationObserver`)
over `tests/testthat` (`:!data/benchmark`): 75 entries, 4 `helper-*` files, 71 test files; minus the 13 files that match only `workflowOmicsManagedActions` (the three that match another worker name stay:
`test-kit252-managed-failure-text.R`, `test-kit259-pearson-histogram-figure-failure.R`, `test-workflow-operation-registry.R`) = 58; plus `test-gui-entry-proteomics.R` and `test-omics-module-globals.R` (neither matches) = **60**.
The audit's register says 59 (74 entries on c91ccfe8); **my count differs by exactly one file, `test-kit216-prot-session-after-enrichment.R`**, which train 35 added on the base side and which matches the grep.

**Settings** (from train 34's request `kit273-train34-9dd61ff8.json`): `file_timeout_minutes` 45; `file_timeouts` `test-kit266-limpa-combined-qc-figure.R` 90, `test-prot-norm-core.R` 60, `test-gui-entry-proteomics.R` 120;
`env` `MULTISCHOLAR_TRX_MANIFEST_NOTE=1`; `extra_packages` none; `check` true (both sides); `blas` reference; `os` ubuntu-22.04; no scripts; `tools_ref` empty.
**12 shards per side.** The planner deals files round-robin by index, so a file cannot be alone in a shard (coordinator's ruling: the three heavy files fall in three different shards with the lightest companions).
Order: five files per shard, balanced longest-first on the largest committed seconds per file (the kit216 file has no committed timing: 300 s assumed); the three heavy files are in shards 1, 2 and 4.

Planned minutes per shard (test minutes from the committed seconds; job minutes = test minutes + about 10 setup):

| Shard | Test min | Job min | Files |
|---|---|---|---|
| 1 | 44.5 | 54.5 | kit266-limpa-combined-qc-figure, e2e-proteomics-dia-limpa, lipid-12-norm-module-contracts, metab-03b-norm-module-helper-characterization, omics-module-globals |
| 2 | 32.5 | 42.5 | prot-norm-core, prot-summary-report-templates, metab-05-summary-module-characterization, module-ci-prot-summary, prot-12-summary-module-contracts |
| 3 | 22.9 | 32.9 | prot-norm-input-current, prot-05b-norm-module-contracts, kit245-github-push-safety, module-ci-prot-enrich, module-ci-prot-norm |
| 4 | 22.1 | 32.1 | gui-entry-proteomics, workflow-saved-qc, kit227-report-ruv-sentences, workflow-gui-memory-lifetimes, general-filemgmt-report-helpers-shared |
| 5 | 21.6 | 31.6 | prot-design-log2, omic-summary-module, kit252-managed-failure-text, omic-norm-ruv-design, kit236-report-sentences |
| 6 | 22.1 | 32.1 | prot-norm-worker-figures, kit216-prot-session-after-enrichment, workflow-integration-input, kit272-prot-box-panel-text, module-ci-prot-da |
| 7 | 21.4 | 31.4 | kit212-da-current, omic-norm-module, enrichment-significance, omics-norm-core, prot-05d-norm-module-characterization |
| 8 | 21.7 | 31.7 | workflow-native-session, workflow-analysis-record, kit205-enrichment-background, go-propagation, module-ci-prot-dia-limpa |
| 9 | 21.4 | 31.4 | workflow-omics-worker, kit273-limpa-ruv-delivery, workflow-omics-figures, kit215-not-tested-lists, workflow-operation-registry |
| 10 | 21.4 | 31.4 | workflow-saved-da, kit206-prot-norm-export-guard, kit259-pearson-histogram-figure-failure, kit225-prot-norm-reset-source, workflow-gui-diagnostic |
| 11 | 21.5 | 31.5 | kit244-limpa-pass-through, kit259b-second-correction, kit271-limpa-table-name, workflow-checkpoint-harness, kit174-ruv-grouping |
| 12 | 21.6 | 31.6 | kit248-correlation-skip-record, kit257-norm-own-files, omics-summary-core, workflow-figure-paths, module-ci-metab-da |

Largest shard: 44.5 test minutes, about 55 job minutes of 300 (18 %); even at the per-file timeouts (90 + 4 light files) no shard reaches 75 % of 300 (225). Sum 294.7 test-minutes per side (the base side is lighter: the new file ran about 572 s on base).

**Ceiling: 1,200 job-minutes = 80.0 core-hours** (core-hours = job-minutes x 4 / 60). **Expected about 865 job-minutes (57.7 core-hours)**: 2 x (295 + 12 x 10) = 830, checks 2 x 16 = 32, plan and results about 3.

### Predictions (passed/failed/errors/skipped; base = the run's base side)

A base-side count that differs from the last recorded count for a file is REPORTED, not explained away. The base side already holds KIT-216 (train 35), whose own effects are on both sides.

| File | Base | Head | Verdict |
|---|---|---|---|
| test-kit266-limpa-combined-qc-figure.R | absent | 210/0/0/0 | only-on-head (E10's Copy call is back to R1's form after the last run: this run is the first evidence for that line; a failure there is a RESULT, reported, nothing repaired) |
| test-prot-norm-core.R | 1547/307 | 1555/307, identical failure messages (+8 passed: the `composite_files` counts) | count-changed; a head-only message is a stop |
| test-gui-entry-proteomics.R | 92/101 | 92/101 | message-changed, only in lane `prot_dia_limpa` (see the masked comparison); classic lanes identical |
| test-kit272-prot-box-panel-text.R | 373/0 | 373/0 | same |
| test-omics-module-globals.R | 14/0 | 14/0 | same |
| test-kit259b-second-correction.R 90/0; test-kit227-report-ruv-sentences.R 186/0; test-prot-norm-worker-figures.R 71/0; test-kit273-limpa-ruv-delivery.R 130/0 | as stated | as stated | same |
| every other file (standing: test-workflow-saved-da.R 755/196, test-omic-summary-module.R 1675/28, test-omics-summary-core.R 676/25, test-omic-norm-module.R 631/8, test-workflow-operation-registry.R 1107/2, test-workflow-gui-diagnostic.R 65/0/3, test-workflow-checkpoint-harness.R 628/0/1; test-kit216-prot-session-after-enrichment.R 152/0/0/1 as in train 35's last table) | newest committed table | same | same |
| R CMD check | 0 ERROR, 3 WARNING, 3 NOTE | identical | nothing only on head |

Totals expected: 60 files; 1 only-on-head, 1 count-changed, 1 message-changed, 57 same; no new-failure, no fixed. (The audit predicted `fixed` for test-prot-norm-core.R at R2 and the tool then said message-changed; count-changed is the likely label here since the passed count moves: the counts and messages are the claim.)

**Masked comparison for the gui-entry lane** (Python `-I`, script outside the downloaded folder, as train 34): on `failures-base.txt` and `failures.txt`, step 1: every failure record of `test-gui-entry-proteomics.R` outside `prot_dia_limpa`
is identical on both sides after the run's path masking, and `test-prot-norm-core.R` has 0 differing records (its 307 failure messages are identical; only passed moves, 1547 to 1555). Step 2, in `prot_dia_limpa` (21 steps, failing on both sides): the head-only lines must be exactly
the R2 set: the composite keys `results/proteomics/protein_qc/proteomics_composite_QC_figure.{pdf,png,rds}` at the 11 steps `norm` to `summary_export` and the two `QC_figures` copies at the three Copy steps; the `pre_norm_*` keys from `norm`;
the eight `*_post_filter_assay1` outputs at `norm` and `correlation` (`<table>` to `<plot>`); the store counts; the nine `ruv_covariates_summary` lines; the `artifacts` length lines (+7, +7, +3 x 6, +5 x 3); and the ordering/chunk-header lines of `names(artifacts)`.
After deleting those, each head record equals its base record; any residue is a missed prediction. KIT-271's and KIT-273's keys are on both sides.

### Decision and attempts

- Every prediction holds: TRAIN GREEN for KIT-266 on the stacked tree. Anything else is a result, recorded and not re-run; no test or expectation is removed, weakened or edited. One run; a shard that times out or an instrument failure with no result: handed back, not re-sent.
  An HTTP 500/502 at dispatch with no run created is started once by hand.
