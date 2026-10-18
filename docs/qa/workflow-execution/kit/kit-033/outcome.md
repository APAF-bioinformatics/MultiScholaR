# KIT-033: existing reports and study_parameters.txt from the record

Status 2026-10-18: **complete.** `R CMD check` at `0120bd3a`: 3 warnings, 3 notes (baseline). Install `kit_033` (`browser-smoke.RV7cne`), installed `source_commit` `0120bd3a…` verified, the overlay present; bindings regenerated (control `kit_036` unchanged). Control `kit_036` (unchanged; no counted member). Host pandoc 3.11.

## What changed (one mechanism: one source of truth)

- **`study_parameters.txt` from the record.** When the owner is the on-disk DIA adapter, `createWorkflowArgsFromConfig` builds the analysis record and:
  - writes every value with the record's own formatter (`workflowRecordValueText`), in the "Parameters from Final S4 Object" blocks and in "Additional Configuration Parameters" (same layout as `formatConfigList`, which is unchanged);
  - adds an "Analysis Record:" section with one `[order. state: operation]` block per step, each declared parameter written as `name [provenance] = value` (and `preceded_by` operations as `operation.name`), plus the DA settings.

  Section layout and parser are unchanged (the reports' `ReadStudyParamsNew` reads `[...]` blocks before key-value lines). Other owners (metabolomics, lipidomics, LFQ/TMT and the in-memory DIA owner) take the unchanged path.
- **Reports read the record.** `workflowRecordReportParameters(record_path, params, study_params)` (`R/utils_workflow_analysis_record.R`):
  - without `analysis_record.json` beside `study_parameters.txt`, it returns the parameters untouched;
  - with it, every report variable the record carries takes the record's value through an explicit map (`workflowRecordReportSources`: report variable, record step and parameter, and file section and key);
  - variables neither the record nor the file carries are listed in one note appended to the methods text: "Parameter not recorded; default shown: …".

  `DIANN_report.rmd` and `DIANN_limpa_report.rmd` call it right after `ExtractStudyParametersNew` in `AnalyzeSamples`.

## Deviations (ticket revision 1.2)

- **Criterion 1 (q-value 0.05):** the package fixes the DIA-NN q-value thresholds at 0.01, as found in KIT-028. The criterion is shown instead with a real value on the real journey: the sample filter's cutoff, 3000 on this run.
- **Templates:** `LFQ_report.rmd` and `TMT_report.rmd` are unchanged. They do not use `ExtractStudyParametersNew`, and their path has no on-disk owner, hence no record, until stage 2 (KIT-020).
- **Defaults note:** printed once in the methods text, listing each default shown with its value, rather than beside each value. The templates use the values numerically in about 80 places.
- **Allowlist:** the overlay function lives with the record code in `R/utils_workflow_analysis_record.R`, outside the allowlist.

## Evidence

| Check | Result |
|---|---|
| Real render on a real journey ([report-check.R](report-check.R) → [report-check-screen.json](report-check-screen.json)): the kit-journey-028 project reflinked; `study_parameters.txt` and `analysis_record.json` written from its store; `DIANN_report.rmd` rendered through `RenderReport` with the pre-ticket template (baseline), this ticket's template with the record, and without it | `study_parameters.txt` 536 lines, 0 truncation marks, the Analysis Record section present. All three renders succeed. **Without the record** the report's text equals the baseline (timestamps masked). **With the record** it also equals the baseline: here the file already carries every value the report reads, so no default is shown and the note is empty. |
| The ticket's defect on the same journey: the file without its `[filterMinNumPeptidesPerSample]` block (as a session whose file lacks it) | Without the record the report states "Samples with fewer than **500** peptides", the template default, a value this run never used. With the record it states **3000**, the recorded cutoff. Report: [DIANN_report-with-record.docx](DIANN_report-with-record.docx). |
| `test-workflow-analysis-record.R` (new: the record-path writer has no truncation marks, writes the 40-entry inclusion list in full and every step parameter as the record carries it; an in-memory owner keeps the current writer, no record section, truncated preview; the overlay takes recorded values, names the defaults shown, and returns the parameters unchanged without a record) | 126 pass |

| Dependents of the changed functions (`createWorkflowArgsFromConfig`, the two DIA-NN templates) and the summaries: `test-e2e-manifest.R`, `test-general-filemgmt-config-direct-shared.R`, `…-config-roundtrip-contracts.R`, `…-report-helpers-shared.R`, `…-s4-contracts.R`, `test-module-ci-prot-{da,import,norm,protein-qc,summary}.R`, `test-prot-02c-…`, `test-prot-05b-…`, `test-prot-12-…`, `test-prot-12a-…`, `test-lipid-14a-…`, `test-metab-05-…` | 108, 33, 24, 84, 30, 245, 90, 100, 68, 110, 651, 562, 446, 62, 129, 160 pass; the three proteomics e2e report tests skip as before |

The rendered reports contain no figures: the fixture's figures are not where the report looks for them. The text is what is compared.

## Claim limits

The DIA on-disk path only. Metabolomics, lipidomics and transcriptomics keep the current writer until their stages declare their steps (KIT-031/032, KIT-040), and LFQ/TMT until KIT-020. The manuscript's reproducibility text should claim the DIA path only. The replay member kind (KIT-028) is not in protocol v1.4; naming it is a user decision.
