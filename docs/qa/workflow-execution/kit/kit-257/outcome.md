# KIT-257 outcome: each normalisation first removes its own earlier files

Ticket `tickets/kit/KIT-257.md` r1.2 (step 0 applied Re-check 1's R1-01 to R1-05). Implements VAL-250. Branch `kit-257`
from `v0.9-staging` 63143aa3. KIT-243 had not merged at the branch point, so the removal core is this ticket's
`removeWorkflowOwnFiles` (`unlink(files, expand = FALSE)`). No R was run on the laptop; every run below is a GitHub
runner (`run-tests`).

## Commits

| Commit | What |
|---|---|
| 1e1b0cee | tests only: `tests/testthat/test-kit257-norm-own-files.R` (nothing under `R/`) |
| 413d739a | the fix (`R/`), NEWS, ledger row `VAL-281` |
| 5ad937e4 | records: R1, wording drafts |
| 0705e897, 18f76b28 | records: R2, outcome, runs |
| 27cc9908 | correction pass, tests first (AUDIT-KIT-257-impl I01, I02, I03, I05): E8b, G2 with a `~` path, E10; the W2 harness of `test-prot-norm-core.R` makes the removal the managed dispatchers make (harness code, no stored value edited) |
| 42fcabd6 | correction pass, fix: I01 (option (a)), I02, I05 (notice), NEWS, ledger row |
| (this commit) | correction pass, records: I04, I06 to I12, the finding for the coordinator |

## What changed (at 413d739a)

- `removeWorkflowOwnFiles(files, message)` (`R/utils_workflow_saved_exports.R:674`): the present, non-directory files
  removed with `unlink(expand = FALSE)`; any still present raises `workflow_figure_unwritten`; one logger line names
  the removed files.
- Transcriptomics, lipidomics, metabolomics: `omicNormApplyRunStart` gains `remove` (`R/mod_omic_norm.R:999`), called
  after `omicNormRunStart`'s refusals and before the previous start is selected, inside the same `tryCatch` (a failure
  is returned as a refusal); one tab-log line when anything was removed. The lists: `omicNormOwnFiles`
  (`R/mod_omic_norm.R:1025`): trx CCA panels (png, pdf) and grid workbook in `publication_graphs_dir`; lipid/metab
  `<assay>_ruv_corrected_{pca,density,rle,correlation}.png` in the omic's QC folder, assay names from the tab, else
  the head's record (R1-03). Callers: `runOmicNormPipeline` (`R/mod_omic_norm.R:456`), `runTrxNormalisation`
  (`R/utils_workflow_trx_norm.R:142`).
- Proteomics: `removeProtNormOwnFiles(paths, stage)` (`R/utils_workflow_dia_norm.R:670`): at a normalisation the five
  RUV figures and both results-table names (TSV and RDS); at a correlation the four table files. Positions:
  in process `runDIAReferenceNormalization` directly after the request (`R/utils_workflow_dia_products.R:748`) and
  `runDIAReferenceCorrelation` after `diaCorrelationSource` (`:1812`); managed normalisation in the actions function
  resolved when the actions start (`R/mod_prot_norm_workflow_helpers.R:638`, R1-01); managed correlation in its first
  action function (`:732`); DIA with limpa in process in the tab's entry before the pass-through (`:510`), managed in
  `dispatchProtNormPassThrough`'s action function (`:270`). The pass-through itself is not changed: its body is pinned
  by digest (`tests/testthat/test-prot-norm-core.R:45`).
- Texts (`omicNormGuardTexts`): `norm_unremoved`, `correlation_unremoved` (drafts, `wording.md`).

## Runs

See `runs.tsv`.

### R1 (45618988736, tests-only 1e1b0cee, 1 job, 5.1 test minutes)

Observed (`r1-tests.tsv`, `r1-failures.txt`, `r1-evidence.txt`): 51 passed, 36 failed, 0 errors, 0 skipped. Every
failure is an assertion. Failing tests: E1 to E9, G3, G5 (proteomics, DIA with limpa, lipid, metab), each at its base
value:

| Row | Base value observed |
|---|---|
| E1 | after the declared run the three CCA files present |
| E2 | at the first dispatch the three planted files present |
| E3 | RUV table present beside the current one; Copy: "Using RUV-normalized file names", `RUV_normalised_results.tsv` and `ruv_normalised_results.RDS` delivered |
| E4 | the five planted RUV figures present |
| E5 | lipid and metab: the four figures present; Copy put all four in `QC_figures/` |
| E6 | lipid and metab: present at the first dispatch |
| E7 | the DA arguments' commit (no state gate) ran at `normalized`; the application's DA accepts `normalized` by reading (`R/mod_prot_da_handlers_init.R:19`, `R/utils_workflow_dia_da.R:47`); the first run's RUV table present; Copy delivered it as `RUV_normalised_results.tsv` |
| E8 | after the second pass-through all four table files present |
| E9 | the correlation's table present when the queued normalisation started |
| G3 | the run proceeded, committed and failed at the publication ("the publication's files could not be copied into the project"); `removeWorkflowOwnFiles` absent |
| G5 | every route proceeded: no refusal; the head moved where the route commits |

G1 (trx, proteomics, lipid) and G2 (trx and proteomics user files; lipid users' files, another omic's folder, the
`Study [pilot]` / `Study p` sibling) passed.

Inferred: the reading of every omic was right (stop rule 1 not met); E7's sequence is reachable in the application by
reading (the proteomics DA accepts `normalized`; AUDIT-KIT-257-impl I07).

### R2 (45622679825, base 1e1b0cee against head 5ad937e4, same runner image, 4+4 shards and check, 486.8 job-minutes)

Declared before sending: the ticket's 37 files (KIT-243's file is not on staging), the new file, and 17 direct
dependents of the functions this fix changed that the ticket's grep does not reach (`test-prot-norm-input-current.R`,
`test-kit206-prot-norm-export-guard.R`, `test-kit225-prot-norm-reset-source.R`, `test-kit174-ruv-grouping.R`,
`test-omic-norm-ruv-design.R`, `test-prot-05c-norm-module-seams.R`, `test-prot-05d-norm-module-characterization.R`,
`test-lipid-12-norm-module-contracts.R`, `test-lipid-12a-norm-module-characterization.R`,
`test-metab-03c-norm-module-seams.R`, `test-metab-03d-norm-module-characterization.R`, `test-workflow-saved-qc.R`,
`test-workflow-gui-memory-lifetimes.R`, `test-module-ci-trx-app.R`, `test-workflow-kit-harness.R`,
`test-workflow-server-shared.R`, `test-workflow-omics-application-reads.R`): 55 files (`r2-files.txt`), balanced by
the longest observed per-file minutes. Not run (single runs, for the merge train): `test-kit129-run-path.R`,
`test-prot-design-log2.R`, `test-workflow-native-session.R`, `test-workflow-saved-da.R`.

Observed (`r2-compare.tsv`, `r2-check-compare.txt`, `r2-compare-failures.txt`, `r2-evidence-head.txt`):
- 53 files `same`, including all four GUI-entry lanes (no per-step file listing moved).
- `test-kit257-norm-own-files.R` `fixed`: 101 passed, 0 failed, 0 errors, 0 skipped on head.
- `test-prot-norm-core.R` `message-changed` (1548 passed, 306 failed on both sides; the failures are the laptop-recorded
  stored references that fail on runners on both sides). One test's message changed, "KIT-111 W1 (a): the DIA
  normalisation and correlation clicks give the recorded states, products, provenance, files and display", in two
  fields only:

  | Case (frozen at) | Field | Frozen | Base actual | Head actual |
  |---|---|---|---|---|
  | `dia_normalization_refused` (`tests/testthat/test-prot-norm-core.R:767`) | `files` | f7b72868... | bbfad929... | 53f7545d... |
  | `dia_normalization_retry` (`tests/testthat/test-prot-norm-core.R:801`) | `files` | f7b72868... | bbfad929... | b51063f4... |

  Every other field of every case, and both sides' `files` digest after the correlation steps (`bbfad929...`), are
  identical base against head.
- R CMD check: base and head 0 ERRORs, 3 WARNINGs, 3 NOTEs; findings identical, none only on head.

Observed (R2 head log, `logs/test-prot-norm-core.log:159`, `:161`, `:162`, the logger line of
`removeWorkflowOwnFiles`; the base has no such line): in W1 (a) the correlation re-apply and Skip each removed the
earlier `ruv_normalised_results_cln_with_replicates.tsv` and wrote it again; the normalisation refused by the kernel
after its start (manual k 50; `runSavedOperation` 1, commit 0) removed `ruv_corrected_cancor.png` and
`ruv_normalised_results_cln_with_replicates.tsv`; the retry removed nothing and its publication wrote
`ruv_corrected_cancor.png` again; no correlation followed, so the table stayed absent. Hence three `files` values
(`protNormCoreFiles`: the names and each TSV's md5): base bbfad929 (both cases), head refused 53f7545d (both files
absent), head retry b51063f4 (the table absent). Every other field of every case is identical base against head. The
head stays `correlation_filtered`: this is R1-02's case (section 3; section 5 "Revert after a re-run"), declared by
file and line under stop rule 3. W2 (W1 (b)) kept bbfad929 for the same two cases because its harness has no removal
(AUDIT-KIT-257-impl I03). The stored values are not edited; they wait for the end-of-project re-record (the
coordinator's row; the owner's by-file decision).

Coordinator's ruling (2026-11-26): stop rule 3 settled; the two digests at `test-prot-norm-core.R:767` and `:801` are
explained by the intended change and nothing else and wait for the end-of-project re-record.

### Correction pass, part 1 (AUDIT-KIT-257-impl-20261126; no run sent)

Done on `kit-257` before merging train 8 (the coordinator's order; R2' follows the merge, part 2):
- **I01 (option (a), the coordinator's choice):** on DIA with limpa the removal runs only when the press designates the
  head (`diaLimpaPassThroughSource(...)$entity_id` against the head's), in process in the tab's entry and managed in
  `dispatchProtNormPassThrough`'s action function; the pinned pass-through body is untouched. Test E8b (a press again
  at `correlation_filtered` and at `ruv_covariates`, in process and managed, keeps the current result's table and the
  covariates' four RUV figures). Note: E8b needs KIT-244's designation rule (train 8) to pass; on this branch alone a
  press at `correlation_filtered` designates the head.
- **I02:** `files <- path.expand(files)` is the first line of `removeWorkflowOwnFiles`. Test: G2 with a `~` path (HOME
  set to a temporary folder inside the test and restored on exit; the row skips, saying so, if this R does not read
  HOME again after start-up, so the real home directory is never touched).
- **I03:** the W2 harness (`protNormW2Normalization`, `protNormW2Correlation`) removes as the managed dispatchers do
  (two lines of harness code; no frozen value edited). Predicted in R2': W2's refused and retry `files` take W1's head
  values (53f7545d, b51063f4), same cause, declared here before the run (INFERRED).
- **I05 (the coordinator's ruling: the user is told):** a proteomics normalisation or correlation that stops after its
  start removed files now carries them (`prot_norm_own_files_removed`) and the tab's error notice gains, on its first
  line, the draft sentence "The earlier RUV figures and normalised table of this analysis were removed when this run
  started; apply or skip the correlation filter again before Copy." (in process and managed; `wording.md`). The
  condition's own message is unchanged (W1 (a)'s `error` field does not move); in W3 the refused run's notice text
  gains the sentence on its first line, which the W3 assertions match by prefix and compare managed against unmanaged
  by first line (INFERRED `same`). NEWS and the ledger row corrected (the table comes back through the correlation
  step; the correlation case added; run ids). Test: E10.
- **I06:** DIA with limpa moved to the change report's appendix.
- Notes recorded below: I07 (E7's row, above), I08, I09, I10, I11, I12.

Notes (AUDIT-KIT-257-impl):
- **I08 (stated limit, INFERRED):** in process the removal follows `diaNormalizationRequest`, which may already have
  marked the pre-normalisation figures as redrawn; if the removal then refuses, the redraw is skipped and the next run
  does not redraw, so the composite's pre-normalisation panels can show the earlier QC table. Needs a locked file and
  a QC change since the last run. Not repaired.
- **I09:** lipid/metab: after RUV then a skip, the tab's "RUV-Corrected" column now shows "Plot not generated yet"
  instead of the earlier run's figures (correct; user-visible; for the owner's wording).
- **I10:** `test-module-ci-trx-da.R`, `test-trx-da-without-ruv.R`, `test-kit196-da-no-value-features.R` reach the
  changed functions and were in neither R2 nor the four declared dependents; each normalises once (predicted `same`);
  added to the scope for the train.
- **I11 (stated limit, INFERRED):** trx/lipid/metab remove at the click; the run rule refuses a busy runner but not a
  pending publication retry, which could write the earlier run's RUV files after the removal. Rare; not repaired.
- **I12:** the logger line in `removeWorkflowOwnFiles` is not wrapped in `logger::skip_formatter` (glue outside the app
  could read `{...}` in a path): not changed (not a line this pass changes). Several added lines exceed 100 characters
  (up to 143), against conventions/R.md's 100; the touched files' idiom runs to 120: recorded divergence. A NULL folder
  is a silent no-op (unreachable from `setupDirectories`).

**Finding for the coordinator (unrelated, no work):** on DIA with limpa the correlation's table is named
`ruv_normalised_results_cln_with_replicates.tsv` before any RUV step (R1 and R2 logs, E8: "the correlation wrote:
ruv_normalised_results_cln_with_replicates.tsv"), because `resolveProtNormCorrelationResultFilenames` treats an absent
RUV result as "RUV applied" (`R/utils_workflow_dia_norm.R:641-659`); Copy then delivers it as
`RUV_normalised_results.tsv`. A misleading name on a route new to 0.9.

## Acceptance criteria

| # | Criterion | Status | Evidence |
|---|---|---|---|
| A1 | Tests-only commit: E1 to E9, G3, G5 fail by assertion at their base values; G1, G2 pass (G4 not applicable) | met | R1 45618988736 |
| A2 | All pass on the fix | met | R2 head: 101 passed, 0 failed |
| A3 | R2 base against head: every file `same` except the new file | not met as stated: one stored digest moved (stop rule 3, above) | R2 `test-prot-norm-core.R` `message-changed` |
| A4 | NEWS, ledger; `outcome.md` with "For the change report" | met (ledger row numbered at merge; F1's route left `<owner/ticket>`) | 413d739a, this file |

## Not checked

- The managed worker paths on a real worker (E2, E6, E9 and the managed G5 cases test the dispatch in the session with
  a recording runner).
- Windows file locks (the case the refusal exists for); the directory stands in for a file that cannot be removed.
- Metabolomics with RUV applied (G1): F172 runs RUV for lipidomics only (KIT-194, KIT-207 precedent); metabolomics
  shares the code path and its removal is shown by E5, E6 and G5.
- G4: KIT-243's file is not on staging at the branch point (not applicable).
- Projects made by earlier builds on real disks.

## For the change report

### Main part (legacy; confirmed with `git show main:<path>`, main 436cdbe8)

**Proteomics: Copy delivered an earlier analysis's normalised table: wrong data in the results package.** The
correlation step writes its results table under a name chosen by the RUV setting
(`main:R/mod_prot_norm.R:2078-2079`, `:2273-2274`: `ruv_normalised_results_cln_with_replicates.tsv` with RUV,
`normalised_results_cln_with_replicates.tsv` without), and nothing removed the other name. Copy prefers the RUV name
whenever that file exists (`main:R/func_general_filemgmt.R:2241-2243`; at the base commit
`R/func_general_filemgmt_results_summary_helpers.R:381-387`). After a normalisation with RUV and the correlation, then
a normalisation without RUV (or with another k, when the differential analysis is then run without a new correlation
step), Copy put the earlier analysis's RUV-normalised abundance table into the results package as
`Publication_tables/RUV_normalised_results.tsv` (and, outside native sessions, its `.RDS`): the package's only
normalised table, while the differential-analysis tables, the report and the study parameters describe the later
analysis. Nothing marked it as earlier. Who was affected: DIA, LFQ and TMT, in every
earlier release. The report's text was not
affected (no proteomics report uses the table's values). How to recognise an affected package: a
`RUV_normalised_results.tsv` beside a report or study parameters saying RUV-III was not applied. Remedy: run the
normalisation (and the correlation step) again, then Copy. The earlier RUV-corrected figures (`ruv_corrected_*.png`,
`ruv_corrected_cancor.png`) also stayed and the native session archive packed them.

**Metabolomics: Copy collected earlier RUV-corrected QC figures.** Copy copies every
`_(pre_norm|post_norm|ruv_corrected)_(pca|density|rle|correlation).png` of the QC folder
(`main:R/func_general_filemgmt.R:2285`; base `R/func_general_filemgmt_results_summary_helpers.R:429`). The RUV
figures are drawn only when RUV runs, so after a run with RUV skipped `QC_figures/` carried the earlier run's
RUV-corrected figures. The composite and the report did not use them.

### Appendix (new to 0.9)

- **Transcriptomics** (no transcriptomics on `main`): the legacy RUV grid writes the CCA panels and the grid workbook;
  a declared design does not (base `R/utils_workflow_trx_declarations.R:829-847`). Copy collected them by name
  (`R/func_general_filemgmt_results_summary_helpers.R:935`, `:939`) and the report showed the panel whenever the file
  existed (`inst/reports/transcriptomics/transcriptomics_report.rmd:131-133`), beside Methods describing the declared
  design.
- **Lipidomics**: the same Copy pattern as metabolomics (base `:527-531`); no such Copy on `main` (git grep).
- **DIA with limpa** (no limpa route on `main`): the same Copy choice when the correlation is applied again after the
  RUV covariates step, whose RUV result decides the table's name.
- **A proteomics run that stops after its start** (for example a manual k larger than the data allow) leaves the
  current analysis without its RUV figures and normalised table until the steps are run again; its error notice says
  so (AUDIT-KIT-257-impl I05). On DIA with limpa a press again on the correlation's or the covariates' result removes
  nothing (I01).

### Both parts

- What changed in results: no number changes. What Copy and the report deliver now matches the analysis run: after
  the fix the earlier files are absent (Copy names a missing file where nothing current replaces it).
- Seams: base 63143aa3 had none; fix 413d739a at the lines listed under "What changed".
- Projects processed before the fix keep their earlier files until the stage is run again. Revert after a re-run
  leaves the removed files absent until the step is run again; the same after a proteomics normalisation or correlation
  that fails after it has started.
- A file that cannot be removed now refuses the run with a message naming it.
- Evidence: R1 45618988736 (base values above); R2 45622679825 (the new tests pass on the fix; 53 regression files
  unchanged; one stored file-list digest of a proteomics sequence moved as the fix intends, declared above). Limits:
  synthetic fixtures (the "small" proteomics recipe, F172, MCI-030), a recording runner for the managed dispatch, a
  directory standing in for a locked file; no Windows run.
- Not in this ticket: F1, the transcriptomics RUVIII export is never cleared, so after a later normalisation without
  RUV the report can still state that RUVIII was applied with an earlier k and show its figure. Routed to its own
  ticket (parked-20261125 row 38 (2)); until it has an owner this report states it.

## Correction pass, part 2 (after trains 8 and 9; declared before R2' is sent)

- Merged `origin/v0.9-staging` 828dd105 into `kit-257` (4a8801db): conflicts only in `NEWS.md` and the ledger (append
  collisions; both kept, staging's first; KIT-257's row `VAL-281` after VAL-272). Code under `R/` and `tests/`
  merged without conflict.
- M1 tests first (a092f2c2): the DA's own earlier figure under a `~` project path (KIT-243 E16's mocked owner; HOME a
  temporary folder, restored on exit). Then M1 (2c05df04): `removeWorkflowDaOwnFigures`'s removal is
  `removeWorkflowOwnFiles(own, <KIT-243's text with %s>)`; the rendered message is KIT-243's unchanged.
- 19d58720: the W2 harness lines (I03) run only where `removeProtNormOwnFiles` exists, so the tests-only base runs W2
  as staging does.
- R2' binding: base `run-kit257-r2p-base` 0cc18658 = staging 828dd105 plus KIT-257's tests only (the new test file and
  the W2 harness lines); head `kit-257` 19d58720. Files: `test-prot-norm-core.R`, `test-kit129-run-path.R`,
  `test-prot-design-log2.R`, `test-kit257-norm-own-files.R`, `test-workflow-native-session.R`,
  `test-kit243-summary-current-outputs.R`, `test-kit244-limpa-pass-through.R`, `test-workflow-saved-da.R`,
  `test-trx-da-without-ruv.R`, `test-kit196-da-no-value-features.R`, `test-module-ci-trx-da.R`; check; env
  `MULTISCHOLAR_TRX_MANIFEST_NOTE=1`, `MULTISCHOLAR_KIT105_W3_DUMP=1` (as merge train 8); 3+3 shards and 2 check
  jobs; ceiling 785 job-minutes (the ticket's correction pass).
- Predictions (staging's current counts, merge-train-20261126a `tests.tsv`):

  | File | Prediction |
  |---|---|
  | `test-kit257-norm-own-files.R` | `fixed`: base fails E1 to E10, E8b, G3, G5, the G2 `~` row and M1 at their base values (G1, G2 pass); head all pass. The two `~` rows skip, saying so, only if this R does not read HOME again (then they observe nothing) |
  | `test-prot-norm-core.R` | 1548/306 both sides; `message-changed` only at W1 (a)'s and W2's `dia_normalization_refused` and `dia_normalization_retry` `files` (head 53f7545d and b51063f4, base bbfad929); everything else, W3 included, identical |
  | `test-kit243-summary-current-outputs.R` | `same`, 246/0 (M1 keeps KIT-243's text and semantics) |
  | `test-kit244-limpa-pass-through.R` | `same`, 204/0 |
  | `test-workflow-saved-da.R` | `same`, 776/175 |
  | `test-kit129-run-path.R` | `same`, 711/0 (1 skipped) |
  | `test-prot-design-log2.R` | `same`, 226/0 |
  | `test-workflow-native-session.R` | `same`, 734/0 (2 skipped) |
  | `test-trx-da-without-ruv.R` | `same`, 103/0 |
  | `test-kit196-da-no-value-features.R` | `same`, 540/0 |
  | `test-module-ci-trx-da.R` | `same`, 52/1 |
  | R CMD check | identical findings |

### R2' (45644271493, base 0cc18658 against head a99e051a, 284.0 job-minutes)

Observed (`r2p-compare.tsv`, `r2p-compare-failures.txt`, `r2p-check-compare.txt`, `r2p-evidence-{base,head}.txt`):
- `test-kit257-norm-own-files.R` `fixed`: base 79 passed, 40 failed; head 133 passed, 0 failed, 0 skipped. Failing on
  base and passing on head: E1 to E7, E9, E8, E10, G3, G5 (proteomics, limpa, lipid, metab), the G2 `~` row and M1. Both
  `~` rows ran (no skip): "path.expand('~') follows HOME: TRUE"; head "own file present: FALSE", "M1 ~ error: none".
  E10 head: the error shown ends "... The earlier RUV figures and normalised table of this analysis were removed when
  this run started; apply or skip the correlation filter again before Copy."
- **E8b passed on both sides, against the prediction** (fails on base). The base is staging plus the tests only, which
  removes nothing at all, so its files stay; E8b's failing-first value against the unconditional removal of the first
  fix (413d739a) was not observed by any run. Head: in process and managed, at `correlation_filtered` and at
  `ruv_covariates`, the table and the four figures present, no error.
  Its value against the first fix is established by reading from observed parts (AUDIT-KIT-257-impl, closure check):
  413d739a removes unconditionally at both limpa sites; E8 shows the same entry removing through this harness; R2'
  shows E8b's presses completing at `correlation_filtered` and `ruv_covariates` with the files present beforehand; so
  on 413d739a E8b fails at the table (both cases) and at the table and figures at `ruv_covariates`. Not run (the
  correction pass is spent); the optional one-job observation is described in the audit's closure check. The
  coordinator's ruling (2026-11-26): no extra run.
- `test-prot-norm-core.R` `message-changed`, 1548/306 both sides, as predicted: W1 (a)'s and W2's
  `dia_normalization_refused` and `dia_normalization_retry` `files` move from bbfad929 (base) to 53f7545d and
  b51063f4 (head); W2 now equals W1 on head (I03 closed). Nothing else in the file differs.
- `same` at staging's counts: `test-kit243-summary-current-outputs.R` 246/0 (M1), `test-kit244-limpa-pass-through.R`
  204/0, `test-workflow-saved-da.R` 776/175, `test-kit129-run-path.R` 711/0 (1 skipped), `test-prot-design-log2.R`
  226/0, `test-workflow-native-session.R` 734/0 (2 skipped), `test-trx-da-without-ruv.R` 103/0,
  `test-kit196-da-no-value-features.R` 540/0, `test-module-ci-trx-da.R` 52/1.
- R CMD check: 0 errors, 3 warnings, 3 notes both sides; findings identical.

Core-hours: R1 1.06 + R2 32.46 + R2' 18.94 = 52.46 of 113.67 (`runs.tsv`'s minutes, 786.6 job-minutes, give 52.44; rounding) (expected 30.0; the correction pass was used).

## Closure check (AUDIT-KIT-257-impl-20261126, 2026-11-26): READY FOR A MERGE TRAIN

- C3: the ledger row's evidence cell gains R2' 45644271493 and M1. Core-hours 52.44 from `runs.tsv` (786.6
  job-minutes).
- **C1 (OPEN, wording defect of the I05 notice):** the notice always says "The earlier RUV figures and normalised
  table of this analysis were removed", but the removed paths vary: E10 itself removed only the table, and a
  correlation start never removes figures. The removed paths are on the condition (`removed_own_files`), so the text
  can follow them. Drafts in `wording.md` for the owner's pass (VAL-41).
- **C2 (OPEN, with C1):** the sentence is joined to the kernel's first line by a single space ("... observed in every
  sample The earlier RUV figures ..."); join with ". " when the first line does not end in punctuation.
- C1 and C2: the small code change (the string choice and the separator) follows the owner's wording pass; this
  ticket's correction pass is spent.
- **C4 (stated limit):** (a) after a partial removal (some files removed, then one that cannot be) the refusal names
  only the file left; the ones already removed are not tagged. (b) I01's condition relies on
  `diaLimpaPassThroughSource` agreeing with the pass-through body's own inline rule. Both read the same states at the
  merged tree (`R/utils_workflow_dia_norm.R:606`, `:630`); the body is digest-pinned, so a change to one shows.

### The train: scope and predictions (the auditor's)

Base staging 828dd105 (or later); head the merge of `kit-257`. NEWS and the ledger conflict; resolve as appends.

| Group | Files | Prediction against staging |
|---|---|---|
| R2's 55 (`r2-files.txt`) | as listed | `same` at staging's counts, except `test-kit257-norm-own-files.R` and `test-prot-norm-core.R` below |
| `test-kit257-norm-own-files.R` | new on head | 133 passed, 0 failed, 0 skipped (R2' head) |
| `test-prot-norm-core.R` | | 1548/306 both sides; `message-changed` at exactly four labels: W1 (a) and W2 `dia_normalization_refused`/`dia_normalization_retry` `files`, bbfad929 to 53f7545d / b51063f4 (the two stored values `:767`, `:801`, declared); nothing else |
| Declared dependents and I10 | `test-kit129-run-path.R`, `test-prot-design-log2.R`, `test-workflow-native-session.R`, `test-workflow-saved-da.R`, `test-trx-da-without-ruv.R`, `test-kit196-da-no-value-features.R`, `test-module-ci-trx-da.R` | `same` (711/0, 226/0, 734/0, 776/175, 103/0, 540/0, 52/1) |
| Partners | `test-kit243-summary-current-outputs.R`, `test-kit244-limpa-pass-through.R` | `same` (246/0, 204/0) |
| I05's other dependent | `test-kit252-managed-failure-text.R` | `same` |
| M1's dependents | `test-enrichment-significance.R`, `test-kit200-prot-da-no-value-proteins.R`, `test-kit201-treat-moderation.R`, `test-kit205-enrichment-background.R`, `test-kit208-gprofiler-identifiers.R`, `test-kit215-not-tested-lists.R`, `test-metab-02ab-da-module-characterization.R`, `test-module-ci-lipid-da.R`, `test-module-ci-metab-da.R`, `test-omic-da-head-fallback.R`, `test-omic-enrich-module.R`, `test-prot-07b-da-handlers-compat.R`, `test-prot-annotation-uniprot-empty.R`, `test-reference-mapping-enrichment.R`, `test-workflow-analysis-record.R` | `same` |
| KIT-259 (if train 10 landed first) | `test-kit259-pearson-histogram-figure-failure.R` | `same` |
| R CMD check | | identical findings |

### Merge resolutions

- **M2 (KIT-254):** the textual conflict after `omicNormApplyRunStart` in `R/mod_omic_norm.R`: keep KIT-254's
  `omicNormCorrelationStart` and `clearOmicNormCorrelation`, then KIT-257's `omicNormOwnFiles`.
- **KIT-259's hunk (`R/mod_prot_norm_workflow_helpers.R:524-527`):** the exact resolved text, replacing the conflict
  block (the two lines after it are already KIT-259's):

```r
      if (!is.null(runner)) dispatchProtNormNormalization(runner, settings, workflowData, normData, experimentPaths,
        omicType, experimentLabel, generatePostNormalizationQcFn, generateRuvCorrectedQcFn, getRuvGroupingVariableFn,
        finish, fail, preQc = generatePreNormalizationQcFn) else
        # KIT-259 (B2): after the completion notice, one warning names the figures that could not be drawn.
        withProtNormFigureNotice(finish(runDIAReferenceNormalization(settings, workflowData, normData, experimentPaths,
          omicType, experimentLabel, generatePostNormalizationQcFn, generateRuvCorrectedQcFn,
          getRuvGroupingVariableFn, pre_qc = generatePreNormalizationQcFn)), showNotificationFn)
```
