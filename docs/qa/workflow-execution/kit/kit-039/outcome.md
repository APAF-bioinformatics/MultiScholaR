# KIT-039: The record read by a reviewer

Status 2026-10-19: **completed; implemented at `134aa8f1`** (ticket revision 1.1); real renders on the kit-journey-028 project recorded; check and install below.

## What changed (one mechanism: rendering)

- `workflowRecordReportSources` gains `software` sources (`software_version`, `source_commit` from `sha`, `build`, `r_version`, `blas`) and `outcomes` sources keyed by operation (`ruv_applied`, `best_k`, `best_percentage`, `num_neg_ctrl` from the normalisation step's `outcomes`); `workflowRecordReportParameters` handles the three source kinds explicitly, guards the file branch, writes `"not recorded"` for a software or outcome field a KIT-033-vintage record lacks, and sets `record_present`. **A finding on the way:** on the real journey the report's `study_parameters.txt` parser turns "Best k value", "Best percentage" and "Control genes" into a logical `TRUE`, so the KIT-038 template (and every one before it) printed "using TRUE% of features as empirical negative controls (TRUE proteins) … TRUE unwanted components"; an outcome or software value is therefore taken as the record has it, never coerced to the parser's type (parameters keep KIT-033's coercion).
- `DIANN_report.rmd` and `DIANN_limpa_report.rmd`: the Version Information section prints the record's version, source commit, build, R version and BLAS when a record is present, the git fields otherwise. The DIA-NN report's RUV sentence now receives the realised values.
- `workflowParameterLabels` (`R/utils_workflow_operations.R`): one label per declared parameter (60; 57 without the diagnostic fixture); the registry test checks both directions.
- `workflowAnalysisRecordTables`: `steps` gains `description`; `parameters` gains `label`; `software` gains build, BLAS, LAPACK, threads, OS, locale rows. The SI template prints each step's `description` as its heading sentence (state, operation, variant and label in a footnote line), its outcomes, a labelled parameter table, and a provenance legend.

## Evidence

- Tests: `test-workflow-analysis-record.R` 177 pass (two new tests), `test-workflow-operation-registry.R` 318 pass (one new test); dependents (`test-e2e-manifest.R`, `test-general-filemgmt-report-helpers-shared.R`, `test-module-ci-prot-{summary,da,import}.R`, `test-prot-12-*`, `test-prot-12a-*`, the e2e limpa/tmt-lfq files which skip): 0 failures, 0 errors.
- Real renders ([report-check.R](report-check.R), [report-check-screen.json](report-check-screen.json), [log](report-check-screen.log); pandoc 3.11) on the kit-journey-028 project with the KIT-038 record: with the record the Version Information section reads "MultiScholaR version: 0.5.0, Source commit: bd0b45ad…, Build: 6b06594f…, R: R version 4.6.1, BLAS: not recorded" (that store's environment entry predates KIT-038) and contains no `NA`; the RUV sentence reads "using 10% of features as empirical negative controls (242 proteins) … 2 unwanted components", equal to the record's outcomes (the match was computed after the run from the saved render text, the run's own regex having been malformed; noted in the JSON); with the record removed the document equals the KIT-038 template's render without the record, dates masked (both rendered without the record, so only the template differs). The SI DOCX ([analysis_record_SI-screen.docx](analysis_record_SI-screen.docx)) contains all 16 step sentences, all 40 distinct labels and the legend; the workbook's `parameters` sheet has 44 rows with a `label` column and no empty label. `study_parameters.txt` written twice by this build is identical modulo the Timestamp line (the writer is untouched).

## Deviations and claim limits

- "BLAS: not recorded" is correct for stores written before KIT-038; it reads the value on stores written from `kit_038` on.
- The KIT-033 baseline comparison is against the KIT-038 commit's template (the latest before this ticket), both rendered without a record.

## Check and install

- `R CMD check --no-manual --no-tests` at `134aa8f1` (clean worktree, `_R_CHECK_FORCE_SUGGESTS_=false`, proof library on `R_LIBS`): 3 WARNINGs, 3 NOTEs, the baseline ([log](check-134aa8f1/r-cmd-check.log)).
- Install `kit_039` = `core/browser-smoke.3cm61C` (source commit `134aa8f1` verified on the installed `DESCRIPTION`); bindings regenerated. The control stays `kit_036` (no counted member).

## Follow-up (2026-10-20, user: "fix"): the report parser's bullet defect

Cause: `ReadStudyParamsNew` in the four proteomics templates stripped only the `•` bullet, while `study_parameters.txt` writes `* Best k value: 2` (`buildWorkflowRuvOptimizationLines`), so every key of the "Automatic RUV Optimization Results" section kept a `* ` prefix and the lookups returned nothing; `ruv_applied` then read FALSE on a run where RUV was applied, and the sentence printed "TRUE% of features (TRUE proteins)" through the report's fallbacks. Fix: the bullet regex accepts `*` (`DIANN_report.rmd`, `DIANN_limpa_report.rmd`, `LFQ_report.rmd`, `TMT_report.rmd`), and the trailing `%` the file carries on "Best percentage" is stripped before the sentence appends its own. Evidence ([report-check-parserfix.json](report-check-parserfix.json), [log](report-check-parserfix.log)): on the real journey's `study_parameters.txt`, parsed with the template's own functions, `ruv_applied` TRUE, `best_k` "2", `best_percentage` "10.0", `num_neg_ctrl` "242"; the DIA-NN report rendered **without** a record now reads "using 10.0% of features as empirical negative controls (242 proteins) … 2 unwanted components" and differs from the KIT-038 template's render only in that RUV sentence (the first difference recorded at position 4016); with the record it reads "10% … (242 proteins) … 2" from the record's outcomes. The template-dependent tests are green.
