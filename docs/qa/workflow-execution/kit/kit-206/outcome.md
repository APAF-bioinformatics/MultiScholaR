# KIT-206 outcome: proteomics session exports refuse a stale table; the native export follows the run

Ticket `tickets/kit/KIT-206.md` revision 1.3 (the implementation audit's correction pass); revision 1.2 applied the plan audit's re-check minors R01 to R06 applied at step 0; the changes
since the audit in its section 19). Authority VAL-134 (1), parity with VAL-105 (2); the coordinator's SD-1 to SD-4 and E-1
to E-5 of 2026-11-24. Branch `kit-206-prot-export-summary` from `v0.9-staging` fe5c9273 (contains KIT-182's merge 3c7fa371).

## Commits

| Commit | What |
|---|---|
| f9f7502f | Tests first: `tests/testthat/test-kit206-prot-norm-export-guard.R` (E1 to E13); nothing under `R/` |
| d268d860 | Instrument repair (VAL-54), tests only: E3 and E12 build their `normalized` head from a run with RUV skipped before the manual run (a run with RUV commits no `normalized` state) |
| 38838e77 | The fix: `R/mod_prot_norm_session_helpers.R` (the guard), `R/mod_prot_norm_workflow_helpers.R` (X1), `R/mod_prot_norm_observer_helpers.R` (X2), `R/utils_workflow_dia_sessions.R` (X6, part (b)), `R/utils_workflow_summary.R` (the Save's named classes); NEWS; ledger row `VAL-<merge>` |
| c3711e2f, fad0311d | Instrument repairs (VAL-54), tests only: E11's setup export and its folder snapshot were lazy arguments read in the wrong order |
| run branch `run-k206-r2-base2` (8f4822a2, 334ff10b) | d268d860 with the same two E11 repairs: the base side of the instrument runs R2-2 and R2-3 (its `R/` equals the tests-only commit's) |

## Runs (`runs.tsv`)

R1 (45355924669) and R1-2 (45358649081, instrument run 1 of 3 on R1's item) on the tests-only side; R2 (45361342951): section 10's
scope, 26 files, fix against tests-only, with R CMD check; R2-2 (45372827072) and R2-3 (45376336123): instrument runs 1 and 2 of 3 on
R2's new file. R1 to R2-3 ran on the pre-merge tree (branch base fe5c9273). R3 (45380286045) ran once on the merge 75360c69 of
this branch with `v0.9-staging` 91e3cdcb (KIT-207 merged, 2ead4fec), against that staging head: the full suite, 20 shards per side,
no check (R2 checked the fix's `R/`; the merge adds only the hand-merged `summaryFailed` condition). Attempt 1 lost head shard 8 in
the library install step (no test ran there); `gh run rerun --failed` repeated it (attempt 2).

## Observed

- **Failing first (A1).** R1-2 at d268d860 fails 28 expectations: the 27 below and E11's set-up check (`added` of length 1; the lazy
  argument repaired in c3711e2f and fad0311d); the base side of R2-3 (334ff10b) fails exactly the 27, 0 errors, all
  in the failing-first list at the values the ticket states: E1, E2, E2b, E3 and E9's second part export successfully (archive
  current state `protein_replicate_filtered`, `protein_intensity_filtered`, `normalized`, `normalized`, `protein_replicate_filtered`;
  E2b's archive carries the first run's correlation result: `TRUE`); E7's S4 export is called with head `protein_replicate_filtered`
  and no notice; E11's refresh writes an archive at `protein_replicate_filtered` and moves the latest pointer; E13 fails with "saved
  table changed during read", leaves `export_pending` set, and the correlation is then refused with "session outputs are pending;
  retry export before dependent actions" (both parts). E12 differs from the ticket's stated base value: the base refuses the export
  for another reason, "native export needs component metadata; explicitly convert the legacy session first" (the native export does
  not support a result-class head); it does not export, so the failing-first expectation still fails at a value, and stop rule 1 does
  not fire. Guards E4, E5, E6, E8, E10 and E9's first two exports pass on the base.
- **On the fix (A1).** R2-3 head fad0311d: 90 passed, 0 failed (verdict `fixed`). Every refusal shows "Error exporting session: " and
  the section 7 text (error, 10 s), writes no archive, summary or latest pointer, leaves nothing pending, and the guard raises
  `prot_norm_correlation_input_stale`. E11's refresh raises the class, writes no archive and leaves the latest pointer; without the
  Revert it writes one at `correlation_filtered`. E12's head is `enrichment_completed`, class `EnrichmentResults`, previous state
  `normalized`, and is refused (the made-from branch). E13: the export with the radio at Automatic after a skipped-RUV run succeeds,
  carries `protein_qc/normalised_results_cln_with_replicates.tsv`, states `ruv_applied` FALSE; the correlation then runs; with that
  file deleted the export is refused with "the file normalised_results_cln_with_replicates.tsv that this export needs is missing from
  <protein_qc folder>; run the step that writes it, then export again." (class `workflow_export_dependency_missing`), nothing pending,
  and the correlation runs.
- **E10 (P02's premise):** on DIA with limpa the `ruv_covariates` record's components differ from `correlation_filtered`'s (`TRUE`),
  so the limpa RUV-state rule is needed; the restored tab exports after the RUV covariates, after the DA's re-commit and after an
  enrichment head (recorder called 3 times) on both commits.
- **Regression scope (A2).** R2: 25 files `same` base against head (including the proteomics GUI-entry lanes, `test-prot-05b`,
  `test-prot-norm-input-current.R` with KIT-182's F1 (h), the native-session, saved-exports and summary-core files); failures present
  on both sides with identical messages: `test-prot-norm-core.R` (306), `test-gui-entry-proteomics.R` (4), `test-kit129-run-path.R`
  (6), `test-omics-summary-core.R` (25), `test-integration-input-transcriptomics.R` (3), `test-workflow-saved-exports.R` (1 error),
  `test-workflow-ruv-pearson-diagnostic.R` (2 errors): stored laptop references, as recorded before (KIT-156). The new file was
  `message-changed` in R2 only for the E11 setup defect; R2-2 and R2-3 repaired it (instrument runs, VAL-54), and R2-3 is `fixed`.
- **Check (A2).** R2: 0 ERROR, 3 WARNING, 3 NOTE on both sides; nothing only on head (`r2-check-compare.txt`).
- **Full suite (A3).** R3 on the merge 75360c69: 464 files `same`, the new file `only-on-head` (90 passed, 0 failed); no
  `new-failure`, no `message-changed`. The files KIT-207's audit names for the second merge are `same`:
  `test-kit207-summary-head-guard.R`, `test-integration-input-lipid-metab.R`, `test-prot-12-summary-module-contracts.R` (quarantined: no expectation runs),
  `test-omics-summary-core.R`, `test-omic-summary-module.R`, `test-gui-entry-proteomics.R` (identical stored-reference failures on
  both). `r3-full-compare.tsv`.
- **A4.** NEWS under the one unreleased heading; ledger row `VAL-<merge>` after VAL-155; no `DESCRIPTION`, `NAMESPACE`, `man/`,
  lock or version change; `helper-workflow-checkpoint-fixtures.R` and `tools/profiling/` untouched.

## Inferred (not run directly)

- The application's Summary "Save Workflow Arguments" shows "The session could not be saved: <text>" on a refusal, since the refresh's
  error reaches `summaryFailed` (read; E11 calls the refresh directly).
- The peptide checkpoint's recorded inputs (search-engine report, IQ output) are now checked before its request; a missing one is
  refused with part (b)'s text instead of failing in staging.
- On the default path part (b) changes no archive member: the run's file equals the radio's, and the lanes (which capture archives as
  present only) are `same`.

## Divergences from the ticket text (recorded)

1. E3 and E12: the ticket's "normalise (RUV `manual`), correlation, Revert twice (head `normalized`)" is impossible on the fixture (a
   run with RUV commits `ruv_corrected` only). The repaired test runs RUV skip, then RUV manual, then the correlation (history ends
   `normalized`, `ruv_corrected`, `correlation_filtered`), and reverts to `normalized`, which is still a table before RUV.
2. E6's enrichment stand-in is checked through the observer with a recorder (the guard's decision only), not a real export.
3. E12's base value is a refusal for another reason, not "success" (above).
4. X2's refusal re-raises with `stop(e)` after the notice (the class is kept).

## Merge with KIT-207 (the shared `summaryFailed`)

This branch changes only the named-classes condition of `summaryFailed` to
`c("workflow_figure_unwritten", "prot_norm_correlation_input_stale", "workflow_export_dependency_missing")`. KIT-207 changes the same
line (adding `integration_input_invalid`) and adds its early-return block for `omic_norm_correlation_input_stale`. The second to merge
writes the union of KIT-207 r1.2 section 3 item 2 (KIT-207's block first, then the four-class condition under "# KIT-052, KIT-206;
KIT-207 (SD-4): a file the save could not write, or refused, is named to the user."). The classes do not overlap. Done in the merge
commit 75360c69 (KIT-207 merged first): KIT-207's refusal block, then the KIT-052 comment, the union comment and the four-class
condition. Since KIT-207, a proteomics Save refused for `integration_input_invalid` also names its reason (KIT-207 I02); no test of this
branch drives a proteomics Save to that refusal.

## Correction pass (implementation audit `AUDIT-KIT-206-impl-20261125.md`; the ticket's one correction)

- **I01.** The guard's result-class block now judges, when the result's `previous_state` names the head itself (a second
  enrichment run), the state before it in the active history (`R/mod_prot_norm_session_helpers.R`). Tests first (run
  45404409911 on 1f69c248 = the branch head plus the tests): E14 fails 5 expectations (the export with a recorder is refused,
  the Summary's refresh raises `prot_norm_correlation_input_stale` and writes no archive, the limpa export is refused; the head's
  `previous_state` is `enrichment_completed`, as the audit read); E12's second enrichment is refused before and after (the
  fallback still judges `normalized`).
- **I02** (coordinator's ruling: fail closed, by parity with VAL-148 SD-5 (b)). `summaryFailed` removes the omic's earlier
  integration input on `prot_norm_correlation_input_stale`. Test first: E11b fails before (after a valid Save the input exists;
  after the Revert the refused Save shows "The session could not be saved: <text>", calls `failed` once, and leaves
  `integration/proteomics/` in place).
- **Run on the merge with staging** (45408209789; head 70dc30fe = the repair commit merged with 8cb536cc; base 8cb536cc; R2's 26
  files, `test-kit207-summary-head-guard.R`, `test-omic-summary-module.R`, with check): as declared. 27 files `same`, the new file
  `only-on-head` with 116 passed and 0 failed (E14 core: `previous_state` `enrichment_completed`, the export proceeds, the refresh
  writes one archive; E14 limpa: the export proceeds; E11b: the valid Save writes `integration_input.json`, the refused Save shows
  "The session could not be saved: <text>", calls `failed` once and removes `integration/proteomics/`; E12's second enrichment is
  refused); R CMD check 0 ERROR, 3 WARNING, 3 NOTE on both sides, nothing only on head (`c2-compare.tsv`, `c2-check-compare.txt`,
  `c2-evidence-head.txt`). Tests-first evidence: `c1-failures.txt`, `c1-evidence.txt`.
- **Module CI on the push of 70dc30fe** (45408115507) failed in "module-ci impacted (transcriptomics)": 10 transcriptomics module-CI
  tests (`test-module-ci-trx-da/design/norm/qc/summary.R`). The same job fails with the same 10 tests on staging 8cb536cc's own push
  (45399571768): pre-existing, not caused by this branch (its `R/` changes are proteomics and the `summaryFailed` classes).

## Hand-offs (not in this ticket; the coordinator opens a ticket; texts from the implementation audit)

- **I04.** The archive's `normalization_method` (and its summary's "Normalization:" line) still come from the tab's selector
  (`R/utils_workflow_dia_sessions.R`, `normalization_method = input$norm_method`), and on reopen `restoreProtNormTabState` sets the
  selector from it. Same class as part (b); the run's method is in the normalisation state's config (`norm_method`). The normalisation
  method is still recorded from the selector.
- **I05.** DIA with limpa: if the RUV covariates step is re-run in another mode after an export, the Summary's refresh passes (the
  limpa rule is by name) and re-uses the earlier export's metadata, so the refreshed archive states the earlier RUV mode. The register's
  five-line repair (re-read `resolveDIARuvRunDescription(manager)` in `refreshDIANativeSummaryArchive`) is not applied: no native
  limpa fixture exists to test it.
- **I09 / F4.** The owned (non-native) and in-memory export routes still record `ruv_mode`, `ruv_applied` and `normalization_method`
  from the tab's inputs (`R/mod_prot_norm_session_helpers.R:136-138`); by reading neither is reached in the shipped application.
- **I10.** At a stale proteomics head the Save still writes `study_parameters.txt` and `analysis_record.json` before the refusal and
  marks the arguments saved (`R/mod_omic_summary.R:263`); "Generate Report" and "Export Session State" are unguarded for proteomics.
  KIT-207 refuses all three in the other omics. A parity decision.
- **D1.** The native export refuses an `enrichment_completed` head with "native export needs component metadata; explicitly convert
  the legacy session first", which is untrue (`R/utils_workflow_dia_sessions.R:1881-1885`; E12's base). The proteomics Summary
  therefore cannot be saved after enrichment unless an export was made earlier in the same application session, and by reading the
  "latest" archive the Summary writes after enrichment cannot be opened by "Load Filtered Session" (adoption reads the head as an
  assay). Owner decision O2 (what a session file's current state may be after enrichment).
- **D2.** Same root cause as I01, pre-existing: after a second enrichment run `protIntegrationInputLineage` stops at the repeated
  name, so the proteomics Save removes the integration input and writes none; the analysis record's `input` for that step names the
  step itself.

## Acceptance criteria

| # | Status | Evidence |
|---|---|---|
| A1 | met | R1-2 (`r1-2-failures.txt`); R2-3 `r2-3-compare.tsv`, `r2-3-failures-base.txt`, `r2-3-evidence-*.txt` (every log line holding `KIT-206` in the run's two logs of the new file, unaltered; regenerated in the correction pass after the first copies were cut) |
| A2 | met | `r2-compare.tsv`, `r2-check-compare.txt`; the new file per R2-3 |
| A3 | met | R3 45380286045 (attempt 2), `r3-full-compare.tsv` |
| A4 | met | the branch diff; this file |

## Not checked

- The running application (no browser): the Summary's notification text, a failed download's appearance.
- A DIA-with-limpa archive reopened into a new owner (E10 restores within one owner).
- The owned (non-native) and in-memory exports still record the RUV mode, `ruv_applied` and the normalisation method from the tab's
  inputs (`R/mod_prot_norm_session_helpers.R:136-138`); by reading they are not reached in the shipped application (the owner attaches
  whenever its required packages are present). Plan candidate F4.
- Directions not executed (implementation audit I08): the S4 download passing at a legitimate head with a reference owner (E7 tests
  the refusal only); the result-class branch on the Summary's refresh at a stale head (S3, DA, enrichment, Save; E12 drives the export
  only, where the base already fails, D1); the limpa exemption's refusal direction (a Revert from `ruv_covariates` to
  `correlation_filtered`; refused by reading); NEWS's `covariates`/`correction`/`none` RUV-mode statement (the limpa owners are not
  native, and the lanes record archives as present only); one Revert to `ruv_corrected` after a single RUV run (E3 tests the two-run
  form).
- That `validateProtDaFilteredSession` accepts `results_file` (read only; the lanes' DA load after export is `same`).
- The S4 export inside a guarded diagnostic runner; Windows; real data.

## User-facing texts for the owner's wording pass (VAL-41; drafts)

- After "Error exporting session: " (export, S4 download) or "The session could not be saved: " (Summary Save): "the current table is
  not the normalisation tab's last result (from the correlation filter, its Skip, or on DIA with limpa the RUV covariates step): a step
  was applied or reverted since that result was made, so the exported session would not describe the current table. Press "Run
  Normalization & RUV" again, apply or skip the correlation filter (on DIA with limpa, then the RUV covariates step), then export."
- Part (b) (draft): "the file <name> that this export needs is missing from <folder>, or is a link or a folder; restore it, or run the
  step that writes it, then export again." Shown after "Error exporting session: " and, in the Summary, after "The session could not be
  saved: ", where "export again" should read "save again"; also shown for the peptide checkpoint's recorded inputs. (The shipped text
  is still "the file <name> that this export needs is missing from <folder>; run the step that writes it, then export again.")
- The Summary Save's pair of notices on any of the four named classes: "The session could not be saved: <text>" (error, persistent),
  then "Study parameters saved with warnings" (warning, 5 s).
- To align with KIT-194's export text (`omicNormGuardTexts()$export`) and KIT-207's Summary texts at the same pass.

## For the change report

- **What was wrong.** In proteomics, after normalisation and the correlation filter, a quality-control step reverted or applied again
  (or a second normalisation without a new correlation step) left the tab's "Export Filtered Data Session" enabled: it wrote a session
  file whose current table was the quality-control (or the second run's) table while carrying the earlier run's correlation result and
  files, and declared itself the normalised data for differential analysis. The Summary's "Save Workflow Arguments" re-wrote the same
  archive from the current table and made it the latest session that "Load Filtered Session" opens. The S4 compatibility download
  (diagnostic runner only) wrote the un-normalised table under the normalised file name. Separately, the native export named its
  results file and its "RUV applied" statement from the current RUV setting, not the run: after a run with RUV skipped and the setting
  changed, the export failed with "saved table changed during read" and left an export pending, which refused every next step until the
  application restarted; where an earlier run's file of that name existed, it was packaged silently. Affected: proteomics users who
  changed quality control after normalising and then exported or saved the Summary, or who changed the RUV setting after a run.
- **Code (base fe5c9273 → fix 38838e77).** Export observer `R/mod_prot_norm_workflow_helpers.R:789` → `:789-793`; S4 download
  `R/mod_prot_norm_observer_helpers.R:219` → `:219-225`; Summary refresh `R/utils_workflow_dia_sessions.R:1685` → `:1724-1725`; the
  native export's metadata `:1865` → `:1902-1911`; its declared results file `:1630-1632` → `:1662-1670`; the dependency check
  before a request `:602` → `:603-605` with `assertDIANativeDependenciesPresent` (`:1615-1643`); the guard
  `R/mod_prot_norm_session_helpers.R:8-55`; the Save's named classes `R/utils_workflow_summary.R:98` → `:98-101`.
- **What changed in results.** No computed value changes. Such exports and Summary saves are now refused with a message; the session
  file's metadata records the run's results file and RUV mode (on DIA with limpa `covariates` or `correction` instead of the setting
  shown); a missing file is refused before anything is pending.
- **Evidence.** R1-2 45358649081 (before), R2-3 45376336123 (before and after, new file), R2 45361342951 (26 files base against head,
  check), R3 45380286045 (full suite on the merge with staging). Synthetic 64-protein fixture only; mechanism, not effect size.

## Budget

Core-hours (job minutes x 4 / 60): R1 0.96, R1-2 0.93, R2 23.25, R2-2 1.91, R2-3 1.75, R3 79.24 (attempt 1 1,148.2 job-minutes with the
lost shard's 8 minutes, attempt 2 40.4). Total 108.04 before the correction pass; the correction pass added the tests-first run 1.17 and the correction run 20.34 (from the
correction row, 222.7): total 129.55 against the ticket's expected 86.8 and ceiling 449.4. Rows: R1 with its instrument
run 1.89 of 4.0; R2 with its two instrument runs 26.91 of 36.0; R3 79.24 of 186.7. Instrument runs used: one of three on R1's item,
two of three on R2's item (VAL-54).
