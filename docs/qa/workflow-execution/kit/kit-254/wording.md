# KIT-254 user-facing texts (drafts; the owner sets the final wording, VAL-41)

No text on the lipidomics or metabolomics normalisation tab changes: the correlation summary and the notices stay as
they are; after the fix they are true after a Skip pressed after an Apply ("Removed: 0 samples" now describes the
table the analysis uses).

| Where | Code | Draft text |
|---|---|---|
| The tab's log (not a notification), when an Apply or a Skip pressed at a correlation-step head first returns to the run's result | `omicNormCorrelationStart` (`R/mod_omic_norm.R`) | "Correlation filter starts again from '<state>', the table the normalisation produced" |
| Replay of an older record with the filter applied and then skipped (OD-1 (a), VAL-247): in the replay's result (`correlation_statement`) and as a warning of class `workflow_replay_correlation_skip` | `workflowReplayCorrelationStatement` (`R/utils_workflow_analysis_record.R`) | "This record shows the correlation filter applied (threshold <t>) and then skipped: the step 'correlation filtering skipped' was recorded on the filtered table. The replay reproduces the analysis as it was computed, so these samples are left out although the filter was skipped: <assay>: <samples or none>; <assay>: ...." |
| NEWS, `### Fixes` | `NEWS.md` | see the entry (KIT-254); drafted from the ticket's section 12, without the bracket on OD-1 (coordinator's ruling) |

Wording candidates for the owner's pass (ticket F3; not changed here):
- After a Skip the summary still shows the earlier Apply's pair statistics under "=== Correlation Filtering Summary ===" (`R/mod_omic_norm_panels.R:1304-1313`): true of the current table, but it reads as a filter that was applied.
- The exported session's summary text prints "Correlation threshold: <slider>" whether the filter was applied or skipped (`R/utils_workflow_saved_exports.R`; legacy, `main:R/mod_lipid_norm.R:1965`, `:2070`).
- KIT-212's stale text lists "a new normalisation run, a Revert, or the correlation filter applied again", not a Skip (`R/mod_omic_norm.R:807-809`); after this fix a Skip pressed after a DA that followed an Apply also makes the DA stale.
