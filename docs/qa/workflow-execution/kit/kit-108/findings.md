# KIT-108 findings

Append only. Ids are F108-NN. A fix that does not change the methodology is made here, with a test that fails first. Anything else is parked.

| Id | Where | What | Status |
|---|---|---|---|
| F108-01 | `copyToResultsSummary` (`:186`, `:668`, `:693`, `:706`, `:716`) | Failures assigned inside closures are lost (recorded at authoring). | Queued to KIT-109. |
| F108-02 | the non-native proteomics integration branch (scoping risk 6) | W1's probe: the branch is not reachable from the app. The proteomics bus owner is always native (`createWorkflowSessionOwner`; every later replacement attaches a native owner). | Not raised. The branch is kept as configuration (W4's filename and write plugins); evidence in `w1/probes.txt`. |
| F108-03 | trx summary: save before the RUVIII export (`mod_lipid_summary_workflow_helpers.R:119-122`, then `:94`) | It shows the 10 s gate warning and then "Saved Integration S4 Object", but no file is written, and the summary says "not written (no integration export yet)". The success notice is misleading. | **Assigned to W3 (lead decision):** after a refused write there is no success notice. The test fails first, and the one W1 capture (trx_main `save_before_export`) is re-frozen once, with attribution. No result changes. |
| (W1 note) | the summary group | W1 added three report-flow logicals (`render_after_failed_template`, `checking_progress_always`, `record_after_success`) and proteomics' `progress_detect` text beyond D1's table, as today's per-omic behaviour, to keep W1 (b) identical. | Accepted by the lead. |
