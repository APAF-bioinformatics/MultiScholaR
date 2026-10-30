# KIT-071 plan (the planning gate): sessions for lipidomics and metabolomics on the owner

Authored 2026-10-30 from KIT-070's tree (1edcb804). One mechanism: the session files the lipid and metab modules
write and read go through the owner. Kernels untouched. Memory decides, time is reported.

## The per-site list (read from the tree)

| Site | File:line | Today | After |
|---|---|---|---|
| Lipid norm export | `R/mod_lipid_norm_workflow_helpers.R:317-560` (`handleLipidNormExportSession`, inline) | `getState(current)` then `saveRDS` of `lipid_filtered_session_data_<ts>.rds` and `_latest.rds`, the three metadata RDS, the summary text | the session list built by `buildLipidNormExportSessionData` (the same fields), published by `runDIAOwnedSessionExport` (KIT-042/043: the owner's publication request, the journal, `export_complete_<id>.json`, `owned_export` on the context) with `saveLipidNormExportArtifacts` as the artifact writer (`runWorkflowLegacyExport` with the `lipid_` prefix, then the metadata files and the summary); a manager that is not an owner aborts (KIT-048's single path) |
| Metab norm export | `R/mod_metab_norm_session_helpers.R:509-760` (`buildMetabNormExportSessionData`, `saveMetabNormExportSessionRdsFiles`, `saveMetabNormExportMetadataFiles`, `saveMetabNormExportSummaryFile`, `runMetabNormExportSessionWorkflow`) | the same, split | `runMetabNormExportSessionWorkflow` runs `runDIAOwnedSessionExport` with the existing builder as the collector and `saveMetabNormExportArtifacts` (the three writers under the publication request, prefix `metab_`); `saveMetabNormExportSessionRdsFiles` requires the request |
| Shared export helpers | `R/utils_workflow_saved_exports.R` (`runDIAOwnedSessionExport`, `serializeDIAOwnedSessionExport`, `stageWorkflowLegacyExport`, `publishWorkflowLegacyExport`, `runWorkflowLegacyExport`) | file names fixed to `filtered_session_data_*` | a `prefix` argument (default `""`, proteomics unchanged); the journal records `latest_name`; the serializer tolerates a session without `correlation_filtered_s4` |
| Lipid DA Open | `R/mod_lipid_da_session_helpers.R:255-345, 452-540` (`readLipidDaSessionData`, `restoreLipidDaCurrentS4FromSession`, `loadLipidDaSessionFromFile`) | `readRDS` then `saveState` of the object under `r6_current_state_name` | the read resolves the owner's export first (`resolveDIAOwnedExport`; the head's object, `owned_generation = TRUE` on the session list); an owned session reopens the owner (no new commit; the head's state is the DA tab's); a legacy RDS on an owner adopts (`ingestWorkflowOmicsLegacySession`); a test double keeps `saveState` |
| Metab DA Open | `R/mod_metab_da_observer_helpers.R:532-600, 695-770` (`runMetabDaLoadSessionObserverShell`, `restoreMetabDaLoadedSessionState`) | the same | the same |
| Adoption | `R/utils_workflow_dia_sessions.R` (new `ingestWorkflowOmicsLegacySession`, after `ingestDIANativeLegacySession`) | proteomics only | the omic's session (`omic_type` and the declared assay class checked; another shape refused by name, `workflow_session_unsupported`; a descriptor in the object `saved_entity_corrupt`) committed by reference onto a fresh typed owner under `r6_current_state_name`, provenance `legacy` with `adopted = TRUE`, the metadata saved on the owner, the previous owner closed, the new one attached on the module's bus (the session from the reactive domain, as `bindWorkflowOwnerTypeAtDesign`) |
| Lipid summary final state | `R/mod_lipid_summary_state_helpers.R:87-92, 300-320`; the exports at `:222-253` and `R/mod_lipid_summary_workflow_helpers.R:48-51` | ordered scan of `getHistory()` (the in-memory history's names) | the scan over the head's state names (`workflowStateManagerStateNames`: `names(snapshot$states)` on an owner, the history on a double); the `session_state_<date>.RDS` and `<omic>_<label>_final_s4.RDS` writers take the object from that pick (the head's reference) |
| Metab summary final state | `R/mod_metab_summary_observer_helpers.R:26-45, 255-265`; writers at `:173-185, 345-351` | the same | the same |

Not changed, recorded: the `session_state_<date>.RDS` and `_final_s4.RDS` writers remain `saveRDS` of the picked object
(the ticket's "from the head's state": the object is the head's reference, hydrated by the owner; no manifest exists
for these two files today and none is added beyond the summary's existing files); `resolveDIAOwnedExport` is reused
as is (it is omic-neutral); the proteomics Open handler untouched.

## Tests

- New `tests/testthat/test-workflow-omics-sessions.R`: for each omic on a real owner in a temporary project: export
  (files named as today, the completion JSON present, `owned_export` set), Open reopens (same head entity, no new
  commit, component hashes equal); a legacy-shaped RDS (the kit_065 session fields) adopts (fresh owner, provenance
  `adopted`, the state under `r6_current_state_name`, hashes equal to the object's checkpoint); an unsupported shape
  (the other omic's object; a proteomics object) refused by name; the summary's pick equals the head's current state.
- Updated: the lipid/metab norm, DA and summary tests that called the plain writers (they take a publication request
  from a real owner or go through the owner entry); the module-CI export helpers on a real owner.
- Module-CI hash comparison before/after (the KIT-069 tooling); the lipid/metab per-process sweep; the shared files'
  tests `test-workflow-native-session.R`, `test-workflow-saved-recovery.R`, `test-workflow-omics-owner.R`; the KIT-028
  replay (the sessions file changed); R CMD check; the install `kit_071`; the two lifecycle members; the kit_065 legacy
  fixtures (written on the InuLC5 library from the e2e lanes, kept under the proof root, hashes recorded) adopted on
  `kit_071` by a script, recorded.

## Exit

Ticket boxes with evidence; `kit-071/outcome.md`; the decision record; memory.
