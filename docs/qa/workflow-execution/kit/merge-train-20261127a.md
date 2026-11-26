# Merge train 34 of 2026-11-27: KIT-273

**Not stacked.** Branched from `origin/v0.9-staging` at `1798d5a824661521e51aa28988b2cd8727fcf63e` (fetched before building; it is the newest staging head and the
one the audit named). It holds trains 25 to 33 (VAL-350 to VAL-370); `test-omics-module-globals.R` reads 14/0 on it. Branch `merge-train-20261127a`,
worktree `/home/doktersmol/Projects/APAF/MultiScholaR-train34`.

## Train

Base: `1798d5a824661521e51aa28988b2cd8727fcf63e` (staging).

| Order | Ticket | Branch head | Merge commit | Ledger rows | Audit |
|---|---|---|---|---|---|
| 1 | KIT-273 | `origin/kit-273` e3a0f60a | 9dd61ff8 | `VAL-<train34-1>` (the ticket's row) and `VAL-<train34-2>` (the stored-reference row); the coordinator numbers both at landing | AUDIT-KIT-273-impl-20261126 (verdict: ready after the record fixes and one ruling) |

Head (bound by the run): `9dd61ff8ee1b1b9cf87dd80bc725e14d7179134f`, the merge commit and the last code-bearing commit. The commit after it
(`4698af2ca`, the record fixes) and the later commits are under `docs/` only.

### Conflicts

`git merge --no-ff origin/kit-273` conflicted in two files, both at the insertion points the audit named; nothing under `R/`, `tests/`, `inst/` or `man/` conflicted
(`R/utils_workflow_dia_products.R` auto-merged).
- `NEWS.md`: staging and the branch both insert after KIT-271's bullet. Every bullet kept, staging's first, then KIT-273's.
- `docs/methods-decision-ledger.md`: staging's rows kept through VAL-370; the branch's row follows it directly, id `VAL-<train34-1>`. The stored-reference row
  (register section 3.6) is added as `VAL-<train34-2>` in the record-fix commit.

### The merged tree against the base (`git diff 1798d5a8 9dd61ff8 --name-only`)

Exactly KIT-273's change, 11 files, 840 insertions and 0 deletions: `R/func_general_filemgmt_results_summary_helpers.R` (+12),
`R/utils_workflow_dia_products.R` (+40), `tests/testthat/test-kit273-limpa-ruv-delivery.R` (new, 608 lines), `NEWS.md` (+15),
`docs/methods-decision-ledger.md` (+1 row) and six files under `docs/qa/workflow-execution/kit/kit-273/`. Nothing under `man/`, `NAMESPACE`, `DESCRIPTION` or `inst/`.

### Record fixes (one docs commit, `4698af2ca`; no test, no R changed)

I01 (the coordinator's ruling quoted in `kit-273/outcome.md` and `runs.tsv`; "no expected value changed" corrected; R1b = instrument run 1 of 3 and the spent correction
pass); I02 (the nine returned run tables committed from the artifacts already on disk, no new fetch; the T8 digests and the lane's digest 31bee7fa... in `outcome.md`);
I03 (Not checked; the k = 0 reason; the observed digest; 104 predicted against 130 observed); I04 (tag v0.9.0 lines and the fix's lines); I05 (the ledger row's wording:
the T9 wording, 27.52 core-hours, the limits, REC-082 and LIM-104 open; the status cell now reads "audited; merge train 34"); I09 (the generic display name) and I10 (the
catch-all divergence) in `wording.md` and `outcome.md`; `wording.md` drafts for the owner's pass (the two file names, the console lines, the NEWS "Skip RUV" sentence
that becomes untrue when KIT-276 lands and is then removed).

## Declaration (written before the run is sent)

**Purpose.** The merge check of KIT-273 on the merged tree, base against head, with R CMD check on both sides. The full suite is not repeated.

**Scope: 48 files.** The ticket's grep recounted on the merged tree
(`git grep -lE "copyToResultsSummary|Publication_tables|ruv_adjusted_protein_matrix|ruv_covariates_design|publishDIAReferenceRuvCovariates|runDIAReferenceRuvCovariates|writeDIADesignTable|getDIAActiveManager" -- tests/testthat ':!tests/testthat/helper-*' ':!data/benchmark'`)
gives **46** files: the 44 of R2, `test-kit272-prot-box-panel-text.R` (it joined with train 31) and the new `test-kit273-limpa-ruv-delivery.R` (45 existing files plus the new file; the audit's "45"
is the existing ones). Plus `test-gui-entry-proteomics.R` (the limpa GUI-entry lane; it does not match the grep) and `test-omics-module-globals.R` (continuity; the branch edits no
file of its pattern) = **48**. R2's own 46 differ from this list by exactly `test-kit272-prot-box-panel-text.R` and `test-omics-module-globals.R`.

**Settings.** Copied from R2's request (`kit273-r2-c9caf84d.json`): `file_timeout_minutes` 45; `file_timeouts` `test-gui-entry-proteomics.R` 45; `extra_packages` none; `check` true
on both sides; `blas` reference; `os` ubuntu-22.04; no scripts. Changed against R2: **4 shards per side** (R2 used 3; the audit recommends 4 because R2's largest shard took 72.55 of 80
minutes), and `env` `MULTISCHOLAR_TRX_MANIFEST_NOTE=1` (the standing rule; R2's env was empty). What the variable changes against R2's base counts: only `test-module-ci-trx-summary.R`
(R2 base 83/2, with the variable 87/1, as trains q and z). The variable is read in `tests/testthat/helper-workflow-trx-declarations.R` only; for the other 47 files the committed tables of trains that set it equal R2's base counts, so no other count is predicted to change.
The files are listed in an order that deals them round-robin so that no shard holds `test-prot-norm-core.R` (about 25 minutes) together with `test-gui-entry-proteomics.R` (about 14 to 18)
or `test-prot-norm-input-current.R` (about 15); R2's per-file seconds sum to about 37 minutes per shard on each side.
Base `1798d5a8`, head `9dd61ff8`; one run.

**Ceiling.** 820 job-minutes = 54.67 core-hours (core-hours = job-minutes x 4 / 60): 8 test shards x 80 + 2 checks x 90. Expected about 410 job-minutes (R2: 383.78 with 3 shards per side; the
fourth shard per side adds its setup).

### Predictions (passed/failed/errors/skipped)

Base counts are the newest committed table that holds the file (trains s to z, listed) and otherwise R2's base table (`kit-273/r2-tests-base.tsv`, run 45740193086, committed).
**Where a base count differs from a figure here, the run's own base side is the reference** (staging moved from R2's base 2d46e4a8, which sits on 7b4ae7da, through trains 31 to 33; files
changed in scope: `test-kit236-report-sentences.R`, `test-kit272-prot-box-panel-text.R`, `test-workflow-figure-paths.R`).

| File | Base | Head | Verdict |
|---|---|---|---|
| test-kit273-limpa-ruv-delivery.R | absent | 130/0/0/0 (R2 head, both paths; the 26 conditional expectations run only on head) | only-on-head |
| test-gui-entry-proteomics.R | 92/101 (train v; R2) | 92/101 | message-changed, ONLY in `prot_dia_limpa` at `summary_copy`, `summary_report` and `summary_export`, by the one new key `results_summary/proteomics/Study_report/design_matrix_with_ruv_factors.tab` (digest 31bee7fa...); no other lane, no other step |
| test-prot-norm-core.R | 1547/307 | same | same |
| test-workflow-saved-da.R | 755/196 | same | same |
| test-kit236-report-sentences.R | 321/0 (train y; R2's 289 is older) | same | same |
| test-kit272-prot-box-panel-text.R | 373/0 (trains x, y) | same | same |
| test-module-ci-trx-summary.R | 87/1 (with the variable; 83/2 without it in R2) | same | same |
| test-omics-module-globals.R | 14/0 (train z) | same | same |
| test-omic-summary-module.R | 1675/28 | same | same |
| test-omics-summary-core.R | 676/25/0/2 | same | same |
| test-workflow-saved-exports.R | 271/0/1/4 (R2) | same | same |
| test-workflow-figure-paths.R | 181/0 (train v; R2) | same | same |
| test-module-ci-prot-dia-limpa.R 38/0/0/2; test-workflow-dia-products.R 191/0; test-workflow-saved-qc.R 120/0 | R2's base (no newer table) | same | same |
| the other 36 files | the newest committed table, else R2's base | same | same |
| R CMD check | 0 ERROR, 3 WARNING, 3 NOTE | identical | nothing only on head |

**Masked comparison for the limpa GUI-entry lane** (as train 26, `merge-train-20261126s.md`; script Python `-I`, kept outside the downloaded folder). On the run's `failures-base.txt` and
`failures.txt`: Step 1: every failure record of `test-gui-entry-proteomics.R`, `test-prot-norm-core.R` and `test-workflow-saved-da.R` outside `prot_dia_limpa`'s steps summary_copy, summary_report
and summary_export is identical on both sides after the run's own path masking (and 0 records differ in the other two files). Step 2: at those three steps, delete every line that names the key
`results_summary/proteomics/Study_report/design_matrix_with_ruv_factors.tab` (as a `+` or `-` line or as an "is a character vector" or "is absent" line), the length lines (`length(actual$artifacts)`
and its `vs` forms), the chunk headers (`names(actual$artifacts)[a:b] vs ...`, whose ranges move with the key), the indented quoted path lines of the chunks (context or re-sorted copies) and blank
lines; the rest of each of the three base records must equal the head record. Any residue is a missed prediction. The key's digest in the head record must be 31bee7fa8d6bf07efd85dac297b1a9148868c8cdf07fdad4fdf759bd4593aeb0.

New failures: none. A head that differs from its base, other than as stated, is a result.

### Scope (48 files)

- test-e2e-lipidomics-canonical.R
- test-e2e-metabolomics-lc-gc.R
- test-e2e-proteomics-enrichment.R
- test-e2e-trx-journey.R
- test-general-filemgmt-report-helpers-shared.R
- test-gui-entry-proteomics.R
- test-kit174-ruv-grouping.R
- test-kit201-treat-moderation.R
- test-kit206-prot-norm-export-guard.R
- test-kit212-da-current.R
- test-kit236-report-sentences.R
- test-kit244-limpa-pass-through.R
- test-kit248-correlation-skip-record.R
- test-kit251-replicate-filter-one-click.R
- test-kit257-norm-own-files.R
- test-kit259-pearson-histogram-figure-failure.R
- test-kit259b-second-correction.R
- test-kit271-limpa-table-name.R
- test-kit272-prot-box-panel-text.R
- test-kit273-limpa-ruv-delivery.R
- test-metab-05-summary-module-characterization.R
- test-module-ci-lipid-summary.R
- test-module-ci-metab-da.R
- test-module-ci-metab-summary.R
- test-module-ci-prot-dia-limpa.R
- test-module-ci-prot-enrich.R
- test-module-ci-prot-summary.R
- test-module-ci-trx-summary.R
- test-omic-norm-ruv-design.R
- test-omic-summary-module.R
- test-omics-module-globals.R
- test-omics-summary-core.R
- test-prot-07b-da-handlers-compat.R
- test-prot-12-summary-module-contracts.R
- test-prot-norm-core.R
- test-prot-norm-input-current.R
- test-prot-norm-worker-figures.R
- test-workflow-analysis-record.R
- test-workflow-checkpoint-dia-state.R
- test-workflow-dia-products.R
- test-workflow-figure-paths.R
- test-workflow-gui-memory-lifetimes.R
- test-workflow-integration-input.R
- test-workflow-native-session.R
- test-workflow-omics-figures.R
- test-workflow-saved-da.R
- test-workflow-saved-exports.R
- test-workflow-saved-qc.R

### Decision and attempts

- Every prediction holds: TRAIN GREEN for KIT-273. Anything else is a result, recorded and not re-run; no test or expectation is removed, weakened or edited. One run; a re-send only for an
  instrument defect with no result judged and with the coordinator's word; an HTTP 500/502 at dispatch with no run created is started once by hand.
