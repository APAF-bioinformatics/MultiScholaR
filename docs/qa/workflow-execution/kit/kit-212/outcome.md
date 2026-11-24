# KIT-212 outcome (ticket r1.4)

Branch `kit-212` from `v0.9-staging` 32766889 (no KIT-216 code in the base: E6's and G4's "base without KIT-216" columns
apply). No R was run on the laptop; every run is listed in `runs.tsv`.

## Step 0 (before R1)

### The survey (KIT-207's form, R2's scope)

Observed by reading at 32766889 (no R run):

- **W1 Summary captures.** The nine cases of section 6 (r1.4) move: `lipid_main`, `lipid_variants`, `lipid_github_fallback`,
  `lipid_global_another`, `metab_main`, `metab_variants`, `metab_github_fallback` (module-CI owners, DA tables on disk, no DA),
  `prot_owner`, `prot_owner_variants` (a session owner, no DA, DA tables on disk). The six proteomics cases on the in-memory
  manager are not checked (S1). `test-omics-summary-core.R` runs lipid, metab and trx cases only.
- **Save, Copy or Report after a table change after a DA, in merged tests.** One move, declared as E14: KIT-207 E11
  (`test-kit207-summary-head-guard.R:604-633`, Copy after a Revert to `lipid_qc_complete`; the DA's source left the history,
  so Copy refuses with `summary_da_stale` before it checks the publication artifacts). Every other KIT-207 row that saves
  after a table change is refused first by KIT-207's own guard (unchanged); its G9 (Reset, then one Revert) keeps the DA's
  source in the history with the same reference, so it stays "current".
- **E15 is not made: KIT-206 E11b does not move.** Its DA is `k206Da` (`test-kit206-prot-norm-export-guard.R:83-86`), which
  only re-commits the DA arguments (`commitDIADaArguments`) and accepts no manifest, so `snapshot$da_results` is NULL, the
  status is "none", and KIT-212's Save guard passes; KIT-206's refusal is reached as before. R2 checks this (the file must be
  `same`). Re-check 3's (m2) assumed a DA manifest.
- **Other files with a DA and a Summary call** (`test-kit172-downstream.R`, `test-workflow-analysis-record.R`,
  `test-workflow-omics-application-reads.R`, `test-workflow-figure-paths.R`, `test-workflow-memory-mechanism-member.R`,
  `test-workflow-checkpoint-harness.R`, `test-workflow-kit-harness.R`, `test-module-ui-wrapper-direct-shared.R`, the GUI-entry
  lanes): the Summary runs at the DA's head, or on a manager that is not a session owner, or without DA files. Not a move by
  reading; R2 judges them.
- **Replay** (the auditor's re-check 3, confirmed): every replaying test record has its DA source among its steps.
- **DA tab observer:** the DA module captures' buses are plain environments, so the observer cannot fire there.
- **F172 threshold:** not needed. E4, E9 to E9c, E10, G11 and G12 do not depend on a sample being removed: a Revert or a
  second run takes the DA's source out of the active history ("stale" by history), and E9c uses the DA on `lipid_log2`
  followed by the normalisation (another table at the head). G11 and G12 print whether the table is the same.

### G0: which DA files the files check judges, per omic

The candidates are the files Copy collects whose names follow the omic's DA publication naming; the made list is computed
from the manifest by the writers' own (now callable) naming. Anything else in the folders (e.g. a fixture's `artifact.txt`)
is "other" and never refused.

| Omic | Folder (Copy) | Candidate naming | Made by (callable) | Kind |
|---|---|---|---|---|
| lipid | `da_output_dir` | `(da\|de)_.+_long_annot.xlsx` | `omicsDaFileNames("lipids", assay, contrast)$table` | per contrast, assay |
| lipid | `Heatmaps` | `<mode>_<contrast>_heatmap.(png\|pdf)`; `all_heatmaps_combined.pdf` | `omicsDaFileNames(...)$plot`; per run | per contrast, assay; per run |
| lipid | `NumSigDeMolecules` | `lipids_num_sig_de_molecules.(tab\|xlsx)`, `lipids_num_sig_barplot.png` | per run | per run |
| metab | `da_output_dir`, `Heatmaps` | as lipid (`metabolites`) | `omicsDaFileNames("metabolites", ...)` | as lipid |
| trx | `da_output_dir` | `de_results_<comparison>.csv` | `trxDaPublicationNames(comparison)$table` | per contrast |
| trx | `Volcano_Plots`, `NumSigDeMolecules` | `volcano_plot_<c>.(png\|pdf)`, `summary_barplot_<c>.(png\|pdf)` | `trxDaPublicationNames` | per contrast |
| prot | `da_output_dir`, `Volcano_Plots`, `NumSigDaMolecules` | tables, `<safe>.(png\|pdf)`, `all_volcano_plots_combined.pdf`, `da_proteins_num_sig_da_molecules.*` | `diaDaPublicationFiles` (the publication's own list, `assertDIADaPublicationFiles`), plus the S4 kernel's `.png`/`.pdf` count plots (`R/func_prot_da_results_methods.R:377-383`) | per contrast; per run |
| prot | `Interactive_Volcano_Plots` | `<friendly>.html`, `<friendly>.unavailable.txt` | `diaInteractiveVolcanoFileName`, `diaInteractiveVolcanoUnavailablePath` | per contrast |
| prot | `Heatmap` (Save Heatmap) | `<prefix>_(heatmap.(png\|pdf)\|heatmap_methods.md\|clusters.csv)` | `omicDaHeatmapSavePrefix("prot_", contrast)` for the manifest's contrasts (any of its three names) | user-saved, by contrast |

Lipid and metab Copy do not collect `Volcano_Plots`; their user-saved heatmaps go to `Heatmap`, which only proteomics' Copy
collects. The lipid and metab names are pinned by G1 and E12 (the real publication's files, none unmade), trx's and
proteomics' by G1.

## Open finding (for the coordinator and the implementation audit): K219-P16, the Report after a later DA

By reading (no R run), r1.4's design does **not** close the route "Save, Copy, run the DA again (other arguments or
contrasts, same table), Report", and no test row of section 10 covers it.
- Report needs only some earlier successful Copy: `values$files_copied` is set by Copy (`R/mod_omic_summary.R:278`) and
  only gates the Report (`:196-205`); nothing resets it when the DA runs again.
- KIT-212's Report guard (`assertWorkflowSummaryDaCurrent(..., files = TRUE)` in `renderWorkflowSummaryReport`,
  `R/utils_workflow_summary.R`, fix 3db806a4) judges the DA status and the **source** folders (`da_output_dir`, the figure
  folders under `publication_graphs_dir`); after the re-run the status is "current" and the source files carry the
  re-run's own names, so it passes.
- The reports read the **publication** folder the earlier Copy filled (`results_summary_dir`:
  `inst/reports/lipidomics/lipidomics_report.rmd:71-74`, `inst/reports/metabolomics/metabolomics_report.rmd:68-71`,
  `inst/reports/transcriptomics/transcriptomics_report.rmd:31-33`, `inst/reports/proteomics/DIANN_report.rmd:68`), so the
  report shows the earlier run's DA tables and figures beside the current Save's parameters.
- Not changed here (coordinator: no design or scope change on this note). Candidate remedies for the ruling: Copy
  records the DA manifest it copied and the Report refuses when the store's manifest differs; or a DA run resets
  `files_copied`.

## What was run, and what was observed

Commits: tests-only e150fd65; fix 3db806a4 (R/, NEWS, ledger, E14); records 668c5920 (control script repair); re-freeze
c854590b; records b004b4db. Runs (ids, minutes, core-hours) in `runs.tsv`.

**Observed.**
- **R1** (45475673660, e150fd65): 176 passed, 95 failed, 0 errors (`r1-tests.tsv`, `r1-failures-by-test.txt`). Every
  failure is in an E row (E1 to E13, E9b, E9c) or the failing-first part of G9 (Q1's managed Save) and G11/G12 (the first
  Save); no precondition failed; every G row passed. Base values printed (`r1-lines.txt`): e.g. E1 lipid `done=1`,
  `source_state=<null>`; E2 `copies=1`; E6 Save "saved with warnings"; E9/E9b/E9c replays complete; E10 the DA ran on the
  loaded `lipid_correlation_filtered`. Stop rule 1 did not fire.
- **R1h** (45475680286, 3db806a4): the new file 303/0; `test-kit207-summary-head-guard.R` (with E14) 292/0.
- **R2** (45479681505, 3db806a4 against e150fd65; `r2-compare.tsv`, `r2-check-compare.txt`): 56 files `same`; the new file
  `fixed`; `test-omic-summary-module.R` and `test-omics-summary-core.R` move, and only in the nine declared cases (the
  module file: the nine; the core file: the eight of the nine it runs, among them `prot_owner` and `prot_owner_variants`, so
  re-check 3's "the core file runs lipid, metab and trx only" was not right: it has a proteomics list,
  `tests/testthat/test-omics-summary-core.R:313`); `test-kit207-summary-head-guard.R` `same` (E14 holds on head), KIT-206's
  file `same` (E15 not needed, as step 0 read); every GUI-entry lane `same`; the DA module captures `same`; check identical
  (0 E, 3 W, 3 N), nothing only on head.
- **Control C** (first run 45489784784: an instrument defect, no case dumped; re-run 45491596893 at 668c5920, `control/`):
  every case dumped in three processes; 152 changed leaves (case/step/field), all at a Copy, Report or later step of the
  nine cases, each equal to the stored value on base (i); the two head dumps identical (ii); (iii) holds. 14 pre-existing
  leaves (the Save steps' `calls`, `files`, `integration` of six cases) left untouched.
- **Re-freeze** c854590b: exactly the 152 changed leaves, each field's text replaced by the head dump's (`git diff
  e150fd65..c854590b -- tests/testthat/test-omic-summary-module.R`; VAL-172 (5) for the seven, VAL-186 for the two).
- **R3** (45494383196, merge 6d23dd52 against staging 9a0bf9f0; `r3-compare.tsv`): 477 files `same`; the new file only on
  head; `test-omic-summary-module.R` `same` (the re-freeze holds); check identical. **One residual, not declared:**
  `test-omics-summary-core.R` 668/33 against 676/25, 8 new failures, each `<case> <report step> notes` of the eight of the
  nine cases the core file runs. The core side drives the Report entry directly (no frame gate,
  `tests/testthat/test-omics-summary-core.R:207-238`) and compares with W1's capture (`omicsSummaryCoreExpected`); on head the
  core's Report shows KIT-212's files refusal, while the re-frozen W1 capture holds the frame's F-11 notice (the tab's Copy
  was refused, so its Report never reaches the core). Not a re-freezable leaf of the core file (it reads the W1 object).
  Stopped for the coordinator (stop rule 2).

**G7 (head, R1h):** status and unmade-file count per W1 fixture: lipid complete none 3, optional_missing none 2, none none 0;
metab complete none 2, optional_missing none 2, none none 0; proteomics in-memory manager none 0 (not checked); proteomics
session owner none 1 (`da_groupB-groupA_long_annot.xlsx`). G15: the in-memory manager copies.

**Inferred.** The refusal paths are reached only where the store's DA left the active history or the table, or where DA
files on disk were not made by the store's DA; every route the ticket lists as allowed (G1 to G15) passed on head, and no
other stored capture moved (R2).

## Acceptance criteria

| # | Evidence | State |
|---|---|---|
| A1 | R1: every E row fails at its base value, G rows pass; R1h/R2: all pass on the fix | met |
| A2 | control C (45491596893): (i) to (iii) hold; the re-freeze commit c854590b edits exactly the 152 changed leaves | met |
| A3 | R2: every file `same` except the new file and the two W1 files at the nine cases; KIT-207's file `same` with E14; no unexpected `new-failure` or `message-changed`; check nothing only on head | met |
| A4 | R3: no difference except A3's, other than the 8 core-side Report notes above | not met: stopped for a ruling |
| A5 | NEWS (Fixes), ledger VAL-<merge>; no version, NAMESPACE or man/ change (all new functions @noRd); this file with "For the change report" | met |

Declared edits of other tickets' files: E14 only (`tests/testthat/test-kit207-summary-head-guard.R:631`, in the fix
commit). E15 was not made (step 0).

## Runs

See `runs.tsv`.

## Not checked

- The managed (worker) DA publication's files on a real worker (G0's made list is the same writer function on both paths;
  the managed lane is the GUI-entry runs, which stayed `same`).
- The lipid G8 export route ran with `k194Input()` as the tab's input; E10's Load Filtered Session was represented by its
  effect on the tab (`tab$current_state_name`), not by an export file and its reopening.
- Trx G8 (the trx DA tab has no Load Filtered Session, `R/utils_omic_gui_config.R:2951`).
- A per-run DA file left by an earlier run when a later run writes none (e.g. `all_heatmaps_combined.pdf` when the later
  DA draws no heatmap): counted as made (G0's per-run rule).
- K219-P16 (above): the Report after a later DA reads the earlier Copy's files; open.
- A browser, Windows, real data.

## For the change report

- **What was wrong.** Every omic's Summary presented the store's last differential analysis whatever table it was run on.
  After the table changed (a Reset and a second normalisation run, a Revert, the correlation filter applied again, a QC
  Revert) without a new DA, 'Save Workflow Arguments' wrote the earlier table's DA block into `study_parameters.txt` and
  `analysis_record.json` (with `source_state` empty or naming a state outside the record's steps), 'Copy to Publication
  Directory' copied the earlier DA's tables and figures, and 'Generate Report' rendered them; Copy and Report also took DA
  files no DA of the session made (a new application session, an opened session file, a later DA with fewer or failed
  contrasts). The lipidomics, metabolomics and transcriptomics DA tab kept showing the earlier results, and a DA run after
  Load Filtered Session ran on the loaded state even after the head had moved. Replay of such a record ran the recorded DA
  on the final replayed table, a different table from the one recorded.
- **Who was affected.** Any user who changed the table after a DA and then saved, copied or reported; who opened a project
  in a new session (or from a session file) with DA files on disk; or who re-ran the DA with fewer contrasts.
- **Code (base 32766889 / fix 3db806a4).** Save `R/utils_workflow_summary.R:124-198` at base (`:157-` at the fix: the DA
  guard `assertWorkflowSummaryDaCurrent`, `:119`, after KIT-207's, and the refusal branch); Copy `:246-316` at base
  (`:281-` at the fix, the guard first); Report `:357-373` at base (`:402-` at the fix); the predicate and the files check
  are new in `R/utils_workflow_dia_da.R` (`workflowDaResultStatus` `:125`, `workflowDaFilesUnmade` `:227` at the fix); the
  DA tab `R/mod_omic_da.R:104-135` at base (the observer and `omicDaFollowHead` `:145` at the fix); replay
  `R/utils_workflow_analysis_record.R:890-917` at base (`assertWorkflowReplayDaSource` `:718`, called at `:918`, at the
  fix); the record's `source_state` rule unchanged (`:504-507` at base).
- **What changed in results.** Nothing computed. Refusals with a message (Save, Report, Copy; Save also removes the omic's
  earlier integration input), a cleared DA display, a DA after a load and a head move taken on the head, and replay
  refusing records whose DA ran on another table (narrowing VAL-99 (3), VAL-139, VAL-154).
- **Evidence.** R1 45475673660 (base values), R1h 45475680286, R2 45479681505, control C 45491596893, R3 (runs.tsv); the
  fixtures are F172 (KIT-172), the module-CI W1 fixtures, the KIT-206 proteomics owner and the MCI-032 trx project.
