# KIT-243: outcome

Ticket `tickets/kit/KIT-243.md` r1.2 (Re-check 1's R1-01 to R1-08 applied at step 0). Branch `kit-243` from
`v0.9-staging` 16112917 (merge train 3 landed). No R was run on the laptop; every run below is a GitHub-hosted runner job
(`run-tests`). Nothing under `data/benchmark/` was opened.

## Step 0 (before R1)

**Cited lines at the branch point (SHOWED).** Every line the ticket cites at the train-3 head c719b6c6 was re-read at
16112917 and is unchanged: `R/utils_workflow_dia_da.R` (`workflowDaCandidateFiles` :153-172, `workflowDaMadeFiles`
:183-217, the in-process publication :548-553, `commitDIADaArguments` :605-628, `writeDIAInteractiveVolcanoOutputs`
:1199), `R/utils_workflow_execution.R` (`publishWorkflowDaResult` :2965-3001, its manifest :2982-2985,
`publishWorkflowOmicsResult` :3007-3032), `R/mod_prot_da_handlers_run_analysis.R` (`report_error` :54-61, the accept
:120-124, `experiment_paths` :2), `R/mod_omic_da.R` (:412, :420, :438-441, :470, :497, :509-510),
`R/utils_workflow_omics_actions.R` (:206-224, :585-614, :596, :605-606, :611-612, :620-629),
`R/utils_workflow_summary.R` (:204-210, :281-361, :287-293, :402, :415-418), `R/mod_omic_norm.R` (`omicNormGuardTexts`
:771-814), `R/utils_workflow_saved_exports.R` (:282-294, :325-425, :342-344, :365-384), `R/utils_workflow_analysis_record.R`
(:468-480, :495-529, :591-595), `R/utils_omic_gui_config.R` (:2126, :2127, :2438). Train 3 is the base, so no further
shift.

**Correction to the ticket's mechanism (by reading).** E2's "KIT-194's stub runner `k194Runner`" delivers
`list(display = list())` for a product action (`k194InProcessRunner`, test-kit212-da-current.R:39-40), so a managed DA
through it accepts a manifest with no contrast and the made list would hold no per-contrast heatmap. E2 and G1 (lipid,
managed) use instead a test runner (`k243Runner`) that runs each requested action as the worker would and delivers it:
the contrasts through their declaration on the head, the tab results through `computeWorkflowOmicsDaTabResultsAction`,
the publication through the in-process publication into the project folders. No expected value changes.

**Survey (KIT-212's form over R2's scope, by reading).** R2's scope is KIT-212's 60 files plus the ticket's grep at the
branch point (24 files, of which 14 new) and the named files: 74 files and the new one (`r2-scope.txt`).
- Tests that run two DA publications in one folder in process: `test-workflow-saved-da.R` W1 (a) and W2, whose frozen
  captures hold a hash of the fixture's files (`savedDaW1Files`, :564-569). A later click there that stops writing a file
  an earlier click wrote (the interactive page against its not-available note, or the count barplot) now finds that
  file removed, so its `files` hash would move (stop rule 3: control C, then the owner's word by file). Which clicks do so
  is not known by reading (the comment at :1155-1156 says the q 0.05 clicks publish no page; `dia_changed_arguments`
  sets TREAT 0.5). Declared here as a possible move; R2 decides.
- Copy or Report pressed after a DA following the last Save: the W1 Summary captures press Copy and Report after their
  Save in the same state (the ticket's 6); KIT-212's own file (G7) is in scope.
- After a failed DA publication: the W1 (a) and W2 `dia_failure` cases fail in the kernel, before the accept, so no
  removal or marker happens there.
- Listings of the publication folder: none found that lists hidden files there (`all.files = TRUE` appears only on other
  folders), so the marker, a hidden file deleted at every completed publication, is not seen by any listing.

## State: STOPPED at stop rule 3 (R2)

R2 shows a capture stored by another ticket moving (`test-workflow-saved-da.R`, the `files` hash of three W1 (a) and three
W2 cases; R2 below). As instructed, nothing is re-recorded. The owner's word by file is needed (stop rule 3) before R3
and the merge. The ledger row is not written yet; it waits for that word.

## G0 (base values; runs 45588535137 and, after an instrument repair, 45590598057)

Observed on the base package (R/ = staging 16112917), `g0-base-values.tsv` (G0b) and `g0-run1-values.tsv` (G0):
- **Record after Save (P03):** `analysis_record.json` exists after Save on lipidomics, metabolomics, transcriptomics and
  proteomics. The proteomics Save with the real writers (KIT-206's fixture) completes without failure.
- **G1 and stop rule 4:** right after Save, the record's steps and DA block equal the store's on every omic. Save itself
  changes neither, so stop rule 4 did not fire.
- **G6:** after an identical DA re-run in the same session the record still equals the store (lipid, metabolomics,
  proteomics).
- **E5's change:** a DA at q 0.1 is seen as a difference at `/differential_abundance/parameters/2/value/value`
  (saved 0.05, current 0.1).
- **F172 (lipid, metabolomics):** no heatmap is drawn. The made list names `posmode_groupB-groupA_heatmap.png/.pdf` and
  `all_heatmaps_combined.pdf`. The metabolomics made list also names the `lipids_num_sig_*` files, which never match a
  metabolomics file (as the ticket says).
- **Proteomics, `savedDaW1Owner` at q 0.05:**
  - It writes the not-available note, not the interactive page, and draws the count barplot.
  - At q 1e-300 the barplot's sha256 is unchanged, so it is redrawn identically. E4 is therefore a guard (stop rule 1:
    dropped as failing-first, with this reason; E1 to E3, E4b and E4m still carry part A).
- **KIT-206's fixture** (`k212Prot`) writes the interactive page.

## R1 (run 45590602596, tests-only 0fafcb03)

167 passed, 50 failed. Every E row fails at its base value (`r1-lines.txt`, `r1-tests.tsv`):
- **E1, E2, E3:** the three planted heatmaps are left, and they are among the candidates.
- **E4b and E4m:** the planted page is left beside the note.
- **E5/E6, E7, E8, E9, E10:** Copy and Report proceed (copies=1, render=1).
- **E12:** DA2 fails on base in KIT-052's backup check ("cannot retain previous action output"), yet Copy and Report
  proceed with the earlier files.
- **E13:** the heatmaps are left and Copy proceeds.
- **E14:** all 8 earlier figures are left and Copy proceeds.
- **E4** passes on base (above).

G rows pass, with one exception: the G6 proteomics Save failed. That case had no tab export before the Save, which every
KIT-212 Save path has. This is test mechanics (VAL-54), not a result. It was repaired in 43e106da (the export added; the
Save's errors shown in the info), with no expected value changed. R2's base side is that repaired commit: there the new
file is 169/49 with no G failure.

The fix (4e2f1575) was drafted before R1 was read and committed after R1's verdict.

## R2 (run 45593292336, fix 4e2f1575 against tests-only 43e106da)

75 files and R CMD check (`r2-compare.tsv`, `r2-check-compare.txt`, `r2-head-lines.txt`).
- **73 files: same.**
- **New file: fixed.** Base 169/49, head 236/0.
  - Every E row passes on head (status `changed`, Copy and Report refused with the drafted texts).
  - E12's removal text names the directory.
  - E13, E14 and G14 show the marker refusing and then cleared.
  - G1 passes on every omic in process and on lipid and proteomics managed: status `current`, Copy and Report proceed,
    no marker.
  - G3 and G4 read `none`; G5, G6 (lipid, proteomics) read `current`.
  - G8 to G11 pass.
  - G6b (observation, not a stop) reads `changed`, as 2.4 predicted (one extra Save after a reopen; F6).
- **Check:** identical on both sides (0 errors, 3 warnings, 3 notes), nothing only on head.
- **`test-workflow-saved-da.R`: message-changed (STOP, stop rule 3).** The file fails on both sides already: its frozen
  values were recorded on the laptop (the KIT-156 record).
  - Only the `files` field of `dia_changed_arguments`, `dia_qvalue_only` and `dia_refused` changes, in W1 (a) and in W2
    alike: a06eeaab... on base, 8688a7e5... on head (`r2-saved-da-files.txt`).
  - The other cases' `files` values are identical on both sides.
  - INFERRED: `dia_changed_arguments` is the first click with TREAT 0.5, which publishes the interactive page where the
    earlier q 0.05 clicks wrote the not-available note. On base the note stays beside the page; on head the note is
    removed (E4b's mechanism). The files listing therefore loses one name: the I03 fix acting on a stored capture.
  - The frozen value (f5079a9e..., laptop) is not touched. Re-recording it needs the owner's word by file.

## Not checked

- Proteomics managed through a real worker process under the application (only `savedDaW2`'s in-process worker entry,
  and the handler with its publication step mocked). Also not checked: a browser, Windows file locks on removal, real
  data, and `openDIANativeSession` (not re-read).
- The marker-deletion warning (R1-02) is not exercised by any test. When the publication folder cannot be written at
  all, KIT-131's check names the marker (not tested).
- A repeated identical request: covered by placement (E4b and E4m drive one in process and one managed).
- Two application sessions on one project (R1-08). A crash between publication and the marker's deletion.
- The proteomics managed path without a runner (`runDIAReferenceDaAnalysis(publish = FALSE)` then the publication step,
  `R/mod_prot_da_handlers_run_analysis.R:131-136`). It is not given S2 (b) by the ticket; a worker failure there leaves
  no marker.
- The proteomics native archive does not carry the marker (R1-05, a stated limit).
- R3 (the full suite on the merge) has not run: stopped at stop rule 3.

## For the change report

- **What was wrong.**
  - Copy and Generate Report decide by file name whether a DA figure belongs to the latest analysis. Each DA publication
    wrote over files by name and never removed any. So a heatmap, combined heatmap PDF, proteomics count barplot or
    interactive volcano page/note that the newest analysis did not draw was the earlier analysis's, and it was copied and
    reported as current with no message. This affects every lipidomics and metabolomics user who runs the DA more than
    once in a project, and proteomics users the same way for the barplot and the page/note pair.
  - After a failed publication, the earlier run's tables were copied under the new run's names.
  - Copy and Report also paired the last 'Save Workflow Arguments' (study_parameters.txt, analysis_record.json) with
    newer results when the table or the DA changed after that Save.
- **Code, base 16112917.**
  - Writers: `R/func_lipid_da_export.R:255-283`, `R/func_metab_da_results_output.R:257-285`,
    `R/func_prot_da_results_methods.R:363-396`.
  - Publications: `R/utils_workflow_dia_da.R:548-553`, `R/utils_workflow_execution.R:2965-3001`, `R/mod_omic_da.R:504-513`,
    `R/utils_workflow_omics_actions.R:585-616`; copy by name `R/utils_workflow_saved_exports.R:385-414`.
  - Copy and Report: `R/utils_workflow_summary.R:281-361`, `:402-423`.
  - Legacy on main 436cdbe8: `R/func_lipid_da.R:1797-1800`, `:1870-1873`; `R/func_prot_da.R:4707-4719`;
    `R/func_general_filemgmt.R:2219`, `:2223`, `:2331-2334`. The not-available note and the analysis record are new in
    0.9.
- **Code at the fix, 4e2f1575.**
  - `removeWorkflowDaOwnFigures`, `workflowDaPublicationMarker`, `clearWorkflowDaPublicationPending`
    (`R/utils_workflow_dia_da.R`).
  - The removal points S1 (`runDIAReferenceDaAnalysis`), S2 (a) (`publishWorkflowDaResult`), S2 (b)
    (`R/mod_prot_da_handlers_run_analysis.R`, the accept), S3 (`exportOmicDaInProcess`) and S4
    (`runWorkflowOmicsManagedDa`, new `experiment_paths = NULL`). The marker is cleared in `completeOmicDaRun`.
  - `assertWorkflowDaPublicationComplete`, `workflowSummarySaveIdentity`, `workflowSummarySaveStatus`,
    `assertWorkflowSummarySaveCurrent` and the two press sites (`R/utils_workflow_summary.R`).
  - Two texts in `omicNormGuardTexts` (`R/mod_omic_norm.R`).
- **What changed in results.** No computed value, table or scientific default changed. What changed:
  - which earlier figure files remain in the publication folder;
  - two new refusals at Copy and Report (`workflow_da_publication_pending`, `workflow_summary_save_stale`);
  - one new DA-tab error (`workflow_figure_unwritten`, the removal text).
- **Evidence.** G0 45588535137 and G0b 45590598057 (base values; Save leaves the record equal to the store). R1
  45590602596 (failing first). R2 45593292336 (fix against tests-only: 73 files the same, new file 236/0, check
  identical). Fixtures: F172 (seed 172001), KIT-206's protein fixture, KIT-112's "small" recipe, the trx module-CI
  corpus. All synthetic; no real data.
- **Limits.**
  - Heatmaps saved with Save Heatmap are kept.
  - When the folder cannot be written at all, nothing is removed and Copy is not refused.
  - The proteomics archive does not carry the marker.
  - Part B fails open (logs and proceeds) when no record was saved or the record is unreadable.
  - One extra Save is asked after a proteomics native reopen and an identical re-run (G6b).
  - Proteomics enrichment run after Save asks for Save again until KIT-216.
- **Classification.** Part A: a legacy defect fixed. Part B: a legacy defect fixed with a new refusal (its record half is
  new in 0.9). Part D: a new refusal on a failure path.
