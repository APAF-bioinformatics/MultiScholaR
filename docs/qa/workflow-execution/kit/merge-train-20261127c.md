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

## Results (run 45766835733, request kit266-train36-fa9a7b8c, head fa9a7b8c4, base 761b43017)

Run success; every file ran to completion (no timeout, no crash, no withheld file). Tables beside this file: `merge-train-20261127c-compare.tsv`, `-check-compare.txt`, `-tests-head.tsv`, `-tests-base.tsv`, `-masked-comparison.txt`.
Base versus head: same 57, only-on-head 1, count-changed 1, message-changed 1, new-failure 0, fixed 0. Tests head: passed 16178, failed 667, errors 4, skipped 90; base: 15960, 667, 4, 90 (the difference 218 = the new file's 210 + 8). **Every prediction held: TRAIN GREEN for KIT-266 (on the stacked tree).**

| File | Predicted | Observed base | Observed head | Verdict |
|---|---|---|---|---|
| test-kit266-limpa-combined-qc-figure.R | absent / 210/0, only-on-head | absent | 210/0/0/0 (1579 s) | only-on-head: HELD. E10 (Copy call restored to R1's form): passed, first evidence for that line is a pass |
| test-prot-norm-core.R | 1547/307 / 1555/307, identical messages, count-changed | 1547/307 | 1555/307 (1949 s, limit 60 min) | count-changed: HELD; its 307 failure records are identical on both sides |
| test-gui-entry-proteomics.R | 92/101 both, message-changed in prot_dia_limpa only | 92/101 | 92/101 | message-changed: HELD (masked comparison below) |
| test-kit272-prot-box-panel-text.R | 373/0 same | 373/0 | 373/0 | same: HELD |
| test-omics-module-globals.R | 14/0 same | 14/0 | 14/0 | same: HELD |
| kit259b 90/0, kit227 186/0, worker-figures 71/0, kit273 130/0, kit216 152/0/0/1 | same | as stated | as stated | same: HELD |
| the standing-failure files (saved-da 755/196, omic-summary-module 1675/28, omics-summary-core 676/25/0/2, omic-norm-module 631/8, operation-registry 1107/2, gui-diagnostic 65/0/3, checkpoint-harness 628/0/1) | same | as stated | same | same: HELD |
| other files | same | no base-side count differed from a prediction | same | same: HELD (I did not compare the unnamed files' base counts with their older tables one by one; the run's own base side is the reference) |
| R CMD check | identical 0/3/3 | 0 E / 3 W / 3 N | 0 E / 3 W / 3 N | identical, 6 findings on both, only on head: none: HELD |

**Masked comparison** (Python `-I`, script kept outside the repository; output `merge-train-20261127c-masked-comparison.txt`). 671 failure records on each side; 660 identical after path masking; 11 differ, all
`test-gui-entry-proteomics.R` lane `prot_dia_limpa` (the 11 steps `norm` to `summary_export`). `test-prot-norm-core.R` and every other file: 0 differing records. At the 11 steps, every differing line is one of the declared kinds:
composite keys, `pre_norm_*` keys, `*_post_filter_assay1` outputs, store counts, the nine `ruv_covariates_summary` lines, `artifacts` length lines and their numbers, `names(artifacts)` ordering and chunk lines. After
deleting those, residue 0, with one note: at the `correlation` step the four keys `results/proteomics/protein_qc/post_norm_{correlation,density,pca,rle}.png` appear as a `+` and a `-` line each (the same key re-sorted
in the `names(artifacts)` chunk, base and head), which the declaration's "ordering lines" mask covers; my first pass of the script did not mask the length numbers or the context line "State saved as: 'ruv_covariates'" and
showed residue only of those lines; the mask was widened to the declared kinds and rerun on the same tables.

**Budget** (jobs API, summed): **805.83 job-minutes = 53.72 core-hours** (core-hours = job-minutes x 4 / 60): ceiling 1,200 (80.0), expected about 865: under both. Longest shard 59.1 job-minutes (`tests head 2`, with test-prot-norm-core.R), 20 % of 300; no shard near 75 %.
Check head 14.2, base 15.57; plan 0.38; results 2.53. A base shard holding the new file (absent on base) took 10.62.

Decision: KIT-266 is green on the stacked tree (stacked on train 35: if train 35 is rebuilt, this train is rebuilt).

## What the coordinator resolves at landing

- Land train 35 first. Ledger: the member's row is the second `VAL-<merge>` (KIT-266; the first is KIT-216's) in `docs/methods-decision-ledger.md` (line 2550); the stored-reference row (draft in AUDIT-KIT-266-impl-20261127 section 8) is added and numbered with it.
  Status cells the member's row cites: the rows the member's row cites: VAL-307, VAL-308, VAL-316, VAL-334, VAL-348, VAL-367 (the coordinator fills those whose status cell waits for this merge).
- Placeholders in the member's records: `docs/qa/workflow-execution/kit/kit-266/outcome.md` (mentions of `VAL-<merge>`); NEWS carries none from this member (to be checked at landing).
