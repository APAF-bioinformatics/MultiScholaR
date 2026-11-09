# KIT-108 entry-point suite

## D5 test 2: F104-01's lane check fails first

The proteomics lanes' `summary_report` step now carries one more expectation: the export summary's "Report Template:" line (`scripts/proteomics/filtered_session_summary.txt`) must name the template of the report the step rendered (`results_summary/proteomics/<template>_proteomics_<lane>.docx`). On the start tree (`5be34ee5` with only the new helper), `prot_dia` run alone fails exactly that expectation: `actual "NA"`, `expected "DIANN_report"` ([lifecycle/start-prot_dia.txt](lifecycle/start-prot_dia.txt)).

**Instrumentation failure, recorded and rerun once.** The check's first version looked for the report in `Study_report/` and errored ("subscript out of bounds") in every proteomics lane. That run's proteomics file had 1 pass and 5 errors, while lipid, metab and trx passed. [lifecycle/attempt1-prot_dia-error.txt](lifecycle/attempt1-prot_dia-error.txt) records it. `f718c5ea` points the check at the report's real location. The start-tree check, `prot_dia` alone and the proteomics entry file were rerun once.

## Simulation Setup (run 2), on the merged tree

Both lanes were run through `kit-105/lifecycle/lifecycle.R` as committed, with strict skips and `R_LIBS` set to kit_107:

- **trx_screen:** identical to its oracle in 40 s. Traced calls: the frame 1, save 1, the integration export 1, copy 1, render 1, session export 1, the export step 1. These are exactly the ticket's expectations ([lifecycle/merged-trx_screen.txt](lifecycle/merged-trx_screen.txt)).
- **prot_dia:** identical to its oracle in 279 s, with F104-01's lane check passing. Traced calls: the frame 1, save 1, copy 1, render 1, session export 1, `collectDIANativeSessionMetadata` 1, `resolveProtSummaryReportTemplate` 3 (the ticket asks for at least 2: the norm export, the report and the session export) ([lifecycle/merged-prot_dia.txt](lifecycle/merged-prot_dia.txt)).

## Run 3: the four entry files

| File | Pass | Fail | Error | Skip |
|---|---|---|---|---|
| test-gui-entry-transcriptomics.R | 29 | 0 | 0 | 0 |
| test-gui-entry-lipidomics.R | 35 | 0 | 0 | 0 |
| test-gui-entry-metabolomics.R | 105 | 0 | 0 | 0 |
| test-gui-entry-proteomics.R (rerun with the corrected check, 774 s) | **193** | 0 | 0 | 0 |

**Decision: all 10 lanes are identical.** Every lane was recorded and JSON-diffed against its oracle, and each shows 0 step differences: the five proteomics lanes on the rerun, the other five on the first run. Proteomics' 193 is KIT-107's 188 plus F104-01's five lane checks. No oracle changed.
