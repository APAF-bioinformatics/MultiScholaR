# KIT-273 step 0: the ticket (r1.2) re-read against the merged tree

Branch kit-273 from origin/v0.9-staging at 7b4ae7dad (KIT-271 merged by train 26, VAL-348). Read by the implementer; no R, no run.

## Bounded work item

Outcome: after Copy to Results Summary on DIA with limpa, in the correction mode the results package holds the matrix the
differential analysis (DA) fitted, and in the covariates mode the whole design with the RUV-III factors the model used;
in no other state either. Owning ticket KIT-273 (VAL-317 (2), VAL-332). Files: `R/func_general_filemgmt_results_summary_helpers.R`
(one call), `R/utils_workflow_dia_products.R` (one new internal function), `tests/testthat/test-kit273-limpa-ruv-delivery.R`,
NEWS, ledger, this directory. Checks: R1 (tests only, failing first), R2 (45 files and check, base against head).
Exit: ticket section 8.

## What moved since the ticket's lines (15d86bd4 -> 7b4ae7dad; KIT-271, KIT-259, train 26)

| Ticket citation | Now | Note |
|---|---|---|
| Copy proteomics block `func_general_filemgmt_results_summary_helpers.R:373-408` | `:371-427`; the table decision `:381-420`; `copyDIANativeSessionArtifact` call `:424-426` inside `if (isDIANativeSession(...))` | KIT-271 added `limpa_tsv_paths` (`:389`) and the first branch of the table decision |
| writers `utils_workflow_dia_products.R:1798-1802` | `:1802` (design, `ruv_covariates_design.tsv`) and `:1805` (matrix, `ruv_adjusted_protein_matrix.tsv`), inside `publishDIAReferenceRuvCovariates` (`:1782`) | unchanged in content |
| `runDIAReferenceRuvCovariates` `:1661-1673` | `:1664` | |
| DA kernel `func_prot_da_model_methods.R:305-326`, record `:374-376` | `:298-335`, record `:374` | read again: `ruv_mode_declared` `:310`, `w_terms` `:317`, refusal of a double correction `:319-320`, `da_model` with model, w_terms, ruv_mode `:374` |
| archive list `utils_workflow_dia_sessions.R:1684-1688`, removal list `utils_workflow_dia_norm.R:777-785` | not edited by design B; both files moved (KIT-271) | |
| `assertWorkflowSummaryDaCurrent` call `utils_workflow_summary.R:353` | `:353` (also `:262`, `:487`) | |
| `workflowDaResultStatus` `utils_workflow_dia_da.R:125-145`; `accepted$model` `:568`; manifest model `:737-748` | unchanged; the manifest's `model` (model, ruv_mode, w_terms) is `manifest$model`, set by `acceptDIADaManifest` from the first contrast's `model` | the record is `manifest$model`, not a per-contrast field |

The design the ticket decided stands unchanged: one internal function beside `publishDIAReferenceRuvCovariates`, one call in
Copy's proteomics block after the native-session artifact, from the store, by the current DA's own record and source.

## Predictions on top of KIT-271

- The delivered table is `Publication_tables/normalised_results.tsv` in every limpa mode (KIT-271); the new call never touches it.
- The four file-name key lines of the limpa GUI-entry lane (`results/proteomics/protein_qc/normalised_results_cln_with_replicates.tsv`
  and its `ruv_` form, `results_summary/proteomics/Publication_tables/normalised_results.tsv` and `RUV_normalised_results.tsv`)
  are identical on base (tests-only commit) and head of R2, because both contain KIT-271; the ticket's "one new key per Copy
  step" is the only expected difference of the lane: `Study_report/design_matrix_with_ruv_factors.tab` at summary_copy,
  summary_report, summary_export. Compared directly base against head and again with the four key lines masked, as train 26 did.
- Predictions per file: `predictions-r2.md`; R1: `predictions-r1.md`.

## I04 (routed by AUDIT-KIT-271-impl-20261126)

"In correction mode the S4 compatibility download leaves the RUV-corrected object in the protein QC folder under a normalised
file name." Answer: the ticket's scope as written does NOT cover it. The ticket changes Copy's delivery from the store
(`copyToResultsSummary`, one new function in `R/utils_workflow_dia_products.R`) and names the files Copy delivers
(VAL-332); the download is `protNormS4DownloadHandler` (`R/mod_prot_norm_observer_helpers.R:~212-233`), a different
function in a different file, not in the ticket's file list, and no VAL-332 name concerns it. Not widened. Candidate, smallest
repair: in that handler, when the head state is `ruv_corrected`, write the canonical file under the `ruv_normalised_results_cln_with_replicates.RDS`
name (the name that was true before KIT-271) instead of `names$rdsFilename`; one condition in one function, one test.
Limits as the audit states them (the user's download is named `<omic>_<label>_final_s4.RDS`; Copy drops both RDS names in
native sessions; the next normalisation or correlation removes the file).

## Known state respected

"Skip RUV" after a RUV run commits nothing today (KIT-276, VAL-346, not implemented; plan audit P05). Not repaired here.
Where Copy decides by the DA's own record, NEWS says what is true in that state: after a RUV run and "Skip RUV" the next
DA still uses the factors or the corrected matrix, so the delivered file is the one that DA used.
