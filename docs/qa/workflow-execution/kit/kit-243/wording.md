# KIT-243: user-facing texts (drafts for the owner's wording pass, VAL-41)

Every new text this ticket shows to a user, as written in the fix commit. All four are drafts; the owner sets the final
wording. The ticket's section 7 states which cases each one is true in.

| Key | Where it is shown | Where it is in the code | Draft |
|---|---|---|---|
| `summary_save_stale` (class `workflow_summary_save_stale`) | Copy: after "Copy error: "; Report: after "Report generation failed: " | `omicNormGuardTexts()` (`R/mod_omic_norm.R`) | The saved workflow parameters do not describe the current analysis: a step or the differential analysis was run, reverted or changed after the analysis record (analysis_record.json) was last saved, or it was saved by another version of MultiScholaR. Press 'Save Workflow Arguments' again, then copy the results or generate the report. If this message returns after saving, the analysis record could not be written: see the warning shown when saving. |
| `da_publication_pending` (class `workflow_da_publication_pending`) | Copy and Report, as above | `omicNormGuardTexts()` (`R/mod_omic_norm.R`) | The files of the last differential analysis were not all written to the project folder, so it may still hold files of an earlier analysis. Run the differential analysis again, then copy the results or generate the report. |
| Removal failure (class `workflow_figure_unwritten`) | The DA tab: proteomics "Error in DA analysis: <text>"; lipidomics, metabolomics, transcriptomics in process "Warning: Could not save results to disk: <text>", managed "Analysis error: <text>" | `removeWorkflowDaOwnFigures` (`R/utils_workflow_dia_da.R`) | An earlier file of the differential analysis could not be removed before this analysis's files were written: <file>. Close any program that has it open, then run the differential analysis again. |
| Marker-deletion failure (a warning; R1-02) | Logged as a warning, and a warning notice on a Shiny session | `clearWorkflowDaPublicationPending` (`R/utils_workflow_dia_da.R`) | The differential analysis completed, but <file> could not be removed, so Copy and the Report still refuse. Remove that file, or run the differential analysis again. |

`<file>` is the full path. When the marker itself cannot be written (the publication folder cannot be written to), no new
text is shown: KIT-131's existing figure check (`checkWorkflowFigureFiles`) names the marker file ("required figures were
not written", "missing or empty: .multischolar_da_publication_pending", "in <folder>"), and nothing is removed.

## Refinements proposed by the implementation audit (I10), for the owner's wording pass

The code still shows the four drafts above. The new test asserts them word for word (`k243Texts` in
`tests/testthat/test-kit243-summary-current-outputs.R`), so any rewording edits that list in the same commit.

- **`da_publication_pending`** is also shown while a managed publication is still running, when running the analysis
  again is not the remedy (proteomics refuses a second press). Proposed addition: "If the analysis is still running,
  wait for it to finish." After a crash between a completed publication and the marker's deletion, "were not all
  written" is untrue; the remedy still works.
- **Removal failure.** For a folder under the file's name, or a read-only file on Windows, closing a program is not the
  remedy. Proposed: "... Close any program that has it open, or move it out of the folder, then run the differential
  analysis again."
- **Marker not writable.** KIT-131's text then names an internal file as a figure ("required figures were not written",
  "missing or empty: .multischolar_da_publication_pending", "in <folder>"). Proposed plain text: "The results folder
  <folder> cannot be written to, so the earlier files were not removed and this analysis's files were not written."
- **Hand-off to KIT-236 (report sentences).** With KIT-220 and KIT-243 together, an analysis with nothing significant
  leaves no proteomics count barplot; the reports handle the missing file but print "Figure 2: Number of significant
  ..." unconditionally (`inst/reports/proteomics/DIANN_report.rmd:1187-1189`, likewise LFQ and TMT).
