# KIT-104 findings

Append only. Ids are F104-NN.

| Id | Where | What | Status |
|---|---|---|---|
| F104-01 | `collectDIANativeSessionMetadata` (`R/utils_workflow_dia_sessions.R:1857`) | The native export summary shows "Report Template: NA" on every proteomics route (the DIA_LIMPA e2e lane, and the 101m/l/t journeys). It reads `global$report_template %||% summary$report_template`, and neither has been set since KIT-046 moved the template choice to the report step. The report itself uses the right template. Display only. | queued to KIT-108, which owns the summary and report (fixing it changes every proteomics lane's session-export capture) |
| F104-02 | `test-e2e-proteomics-dia-limpa.R`, the study_parameters content checks (W4) | Stale checks: they expect `DIANN_limpa_report.rmd` in study_parameters and a `proteinMissingValueImputationLimpa` args group, neither of which the native path writes. | **fixed** `ce642947` (test-only; the lane now passes 79/0/0/0, no skip) |
