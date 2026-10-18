# KIT-029: supplementary rendering, DOCX and Excel from the record

Status 2026-10-18: **implemented; install pending.** Control `kit_036` (unchanged; no counted member). Host pandoc 3.11 (`rmarkdown::pandoc_available()` TRUE).

## What changed (one mechanism: rendering the record for people)

- **One table builder.** `workflowAnalysisRecordTables(record)` turns the record into data frames: software, packages, inputs, design, steps, parameters (one row per step and parameter, `preceded_by` operations included), analyses (the DA settings), external resources. Values are written in full (`workflowRecordValueText`): atomic values comma-joined, lists as JSON, tables as size and hash. The document and the workbook both read it, so they cannot disagree.
- **Template and renderer.** `inst/reports/shared/analysis_record_SI.rmd` (parameter `record_path`) is omic-agnostic, with a `word_document` output like the omic reports. Every section is a table from the builder. Step headings carry the operation, variant, step label and commit time from the record. The template carries no method prose, the external-resources limit from the record, and a citation line (the package citation, else version and source commit). `renderWorkflowAnalysisRecord(record_path, output_dir, format)` renders it with `rmarkdown::render` on the installed template, in a temporary copy. `RenderReport` is not used and is unchanged.
- **Workbook.** `writeWorkflowAnalysisRecordWorkbook(record, path)` (openxlsx, already an import) writes sheets `software`, `inputs`, `design`, `parameters`, `external_resources`, plus `steps`, `analyses`, `packages`.
- **Export.** `exportWorkflowAnalysisRecordBundle(manager, directory)` writes `analysis_record.json`, `.docx` and `.xlsx`. `runProtSummaryReportProgress` calls it after **Generate Report** succeeds, for an on-disk owner, into the project's `study_report_dir`. A failure shows a warning and leaves the report as generated.
- **Harness.** `kitFingerprintFiles` excludes `analysis_record.{json,docx,xlsx,html}`.

Deviation: the ticket names an "Export analysis record" action. A separate button needs the summary UI and server files, which are outside the allowlist, so the bundle is written by the existing Generate Report action, the moment users produce the outputs they attach.

## Evidence

| Check | Result |
|---|---|
| The real journey's record (the KIT-028 origin, [analysis_record-origin-screen.json](../kit-028/analysis_record-origin-screen.json)) rendered on this host: [analysis_record.docx](analysis_record.docx), [analysis_record.xlsx](analysis_record.xlsx) | Rendered in 2.8 s (DOCX 23.6 kB; HTML 705 kB, not kept). Workbook `parameters` has 97 rows = the record's 97 parameters. |
| `test-workflow-analysis-record.R` (new: real DOCX and HTML renders containing every parameter name, every input sha256 and a 300-element vector in full, called with a relative record path; workbook row count equal to the record's parameters, a 300-element value in full, inputs, software and external resources matching; the export hook on success, none on a failed report, and a warning on an export failure) | 108 pass |
| `test-workflow-kit-harness.R` (new: the bundle leaves `kitFingerprintFiles` unchanged) | 341 pass |
| Summary and report regressions: `test-prot-12-summary-module-contracts.R`, `test-prot-12a-…`, `test-general-filemgmt-report-helpers-shared.R` (the `RenderReport` tests), `test-lipid-14a-…`, `test-lipid-15-…`, `test-metab-05-…` | 446, 62, 84 (1 skip), 129, 134, 160 (9 skip) pass |
| Report doubles (`helper-e2e-report-doubles.R`) | used by no test file (its one function, `install_report_template_double`, has no caller), so there is nothing to run |
| Existing omic templates and `RenderReport` | `git diff` empty for `inst/reports/{proteomics,lipidomics,metabolomics}` and `R/func_general_filemgmt_report_render_helpers.R` |
| `R CMD check` at `8ad3bf29` (clean worktree) | 3 warnings, 3 notes (baseline) |

A bug was found by rendering the real record with a relative path: `rmarkdown::render` changes directory before it reads its arguments, so the record was not found. The renderer now normalises both paths, and the render test uses a relative path.

## Claim limits

- Excel limits a cell to 32,767 characters, so a single parameter value longer than that (a very long inclusion list) would be cut in the workbook. The JSON and the DOCX carry it in full.
- The DOCX is checked for content, not layout.
