# KIT-274 wording: drafts for the owner's wording pass

Substance decided by the owner (VAL-318, VAL-332). The texts are drafts; the tests hold each as one constant
(`K274Refusal`, `K274Statement`), the code as one string, so a wording pass is a text edit in these places.

| Text | Where | Draft |
|---|---|---|
| Refusal of "Run RUVIII Export" | `trxAssertRuviiiExportApplied`, R/utils_workflow_trx_declarations.R (class `trx_ruviii_no_ruv`, `trx_error`) | no RUVIII export is made: this analysis applied no correction for unwanted variation (k = 0; {reason}), so RUVIII would return the TMM logCPM unchanged; the integration uses the TMM logCPM. |
| What the tab shows | `registerOmicSummaryExportStep` (unchanged): the step's prefix, the message, the class in brackets | RUVIII export error: no RUVIII export is made: ... (trx_ruviii_no_ruv) |
| Replay statement (and warning of class `workflow_replay_export_skipped`) | `replayWorkflowAnalysisRecord`, R/utils_workflow_analysis_record.R (`replayed$export_statement`) | This record lists a RUVIII export made when no correction for unwanted variation had been applied (k = 0). The replay does not repeat it: such an export is no longer made, and in this state the integration uses the TMM logCPM, which the replay writes as before. |
| Report, an older record that lists such an export | the existing never-exported texts, no new text (template guard, transcriptomics_report.rmd:49) | No RUVIII integration export was recorded for this analysis.  /  Figure not available: ruviii_before_and_after.png |
| NEWS | NEWS.md, Fixes, first entry | as in NEWS |

`{reason}` is the recorded `ruv_not_applied_reason`: "the design none: no correction for unwanted variation", or the k rule's
reason ("the estimated count found no unwanted factor (no component above the permutation null)"); "no reason was recorded"
when the state holds none. The bracketed class in the tab is `trxConditionText`'s, as for every transcriptomics refusal;
changing it is a separate ticket if the owner wants it.

## Points from the implementation audit (I06), drafts for the owner's pass

- Replay statement: it must not say the replay "writes" the TMM logCPM (the replay writes no integration file; it re-creates the TMM
  logCPM product in its store). Draft change of the last clause: "...and in this state the integration uses the TMM logCPM." (the clause
  "which the replay writes as before" removed). The code string in `replayWorkflowAnalysisRecord` and the test constant `K274Statement` both
  carry the current draft and move together at the wording pass.
- NEWS may add: a leftover ruviii_before_and_after.png stays in the project, and Copy delivers it until the normalisation is run again.
