# KIT-254 user-facing texts (drafts; the owner sets the final wording, VAL-41)

No text on the lipidomics or metabolomics normalisation tab changes: the correlation summary and the notices stay as
they are; after the fix they are true after a Skip pressed after an Apply ("Removed: 0 samples" now describes the
table the analysis uses).

| Where | Code | Draft text |
|---|---|---|
| The tab's log (not a notification), when an Apply or a Skip pressed at a correlation-step head first returns to the run's result | `omicNormCorrelationStart` (`R/mod_omic_norm.R`) | "Correlation filter starts again from '<state>', the table the normalisation produced" |
| Replay of a record with the skip step directly after the filter (OD-1 (a), VAL-247): in the replay's result (`correlation_statement`) and as a warning of class `workflow_replay_correlation_skip` | `workflowReplayCorrelationStatement` (`R/utils_workflow_analysis_record.R`) | "This record applies the correlation filter at threshold <t> and then a step named 'correlation filtering skipped' on the filtered table. The replay repeats both as recorded, so its analysis leaves out <samples> (below the threshold in <assays>). To analyse every sample, run the normalisation again and press 'Skip Filtering & Proceed to DE' without applying the filter." (none when no sample left the design) |
| Replay of a record in which the filter was applied more than once (OD-1 (a) as put; coordinator's ruling on AUDIT-KIT-254-impl I08) | the same | "The correlation filter was applied more than once in this record; only the last threshold (<t>) is recorded, so this replay applies it once to the normalised table, and its samples may differ from the analysis that was run." |
| NEWS, `### Fixes` | `NEWS.md` | see the entry (KIT-254); drafted from the ticket's section 12, without the bracket on OD-1 (coordinator's ruling) |

Wording candidates for the owner's pass (ticket F3; not changed here):
- After a Skip the summary still shows the earlier Apply's pair statistics under "=== Correlation Filtering Summary ===" (`R/mod_omic_norm_panels.R:1304-1313`): true of the current table, but it reads as a filter that was applied.
- The exported session's summary text prints "Correlation threshold: <slider>" whether the filter was applied or skipped (`R/utils_workflow_saved_exports.R`; legacy, `main:R/mod_lipid_norm.R:1965`, `:2070`).
- KIT-212's stale text lists "a new normalisation run, a Revert, or the correlation filter applied again", not a Skip (`R/mod_omic_norm.R:807-809`); after this fix a Skip pressed after a DA that followed an Apply also makes the DA stale.
- AUDIT-KIT-254-impl I14 (sharper after the fix): KIT-212's stale text (`R/mod_omic_norm.R:827-829`) names no Skip, yet a
  Skip after a DA that followed an Apply now makes that DA stale; the exported session's summary prints "Correlation
  threshold: <slider>" (`R/utils_workflow_saved_exports.R:224-231`) after a Skip whose data hold every sample.
- NEWS (I07): the multi-assay qualification and "the table the normalisation produced (RUV-corrected when RUV was
  applied)" are drafts under VAL-41 too.
- Closure check C3 (draft for the owner's pass) for the second statement: "The correlation filter was applied more than
  once when this record was made; only the last threshold (<t>) is recorded, so this replay applies it once, to the table
  the normalisation produced (RUV-corrected when RUV was applied), and its samples may differ from the analysis that was
  run." ("in this record" and "the normalised table" are loose; on a record with leftovers of a Reset and a second run,
  plan audit U1, "applied more than once" may be false while "may differ" stays true.)
- Closure check C2: the NEWS replay sentence above replaced "replay as they ran" (false for a filter applied more than
  once); a draft under VAL-41 too.
