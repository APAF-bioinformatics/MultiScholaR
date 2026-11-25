# KIT-256 user-facing texts (DRAFTS; the owner sets the final wording, VAL-41)

Nothing in the interface's text changes in KIT-256: the box simply empties after each successful 'Assign Technical
Replicates', and the notification (`assignment$message`) is unchanged. The texts below are drafts for the owner's
wording pass. The help sentence is not code in this ticket unless the coordinator rules it in.

| Where | Status | Draft text |
|---|---|---|
| Technical-replicate help, appended on both UIs (`R/mod_omic_design_builder.R:129-130`; `R/mod_prot_design_builder_ui.R:142-148`) | DRAFT, not code | "After each assignment the list is emptied. To add runs to a sample you have already marked, select all of its runs again together with the new ones; runs selected on their own become a separate sample." |
| `NEWS.md`, `### Fixes` | in NEWS as a draft (VAL-41) | "Design builder, every omic: the technical-replicate list is emptied after each successful 'Assign Technical Replicates'. Before, in every earlier release with the design builder, it kept the runs just marked; adding the next sample's runs and pressing Assign made two biological samples one, with no warning. To add runs to a sample already marked, select all of its runs again together with the new ones; runs selected on their own become a separate sample. A refused Assign keeps the selection. (KIT-256)" |

The reverse case both drafts state (ticket section 5): after the fix the box no longer remembers the runs just marked.
To add runs to a sample already marked, the user selects all of that sample's runs again together with the new ones,
as the refusal for a partial selection already asks (`R/func_omics_design_shared.R:662-665`). Selecting only the new
runs makes them a further biological sample, with no warning: technical replicates analysed as separate samples
(pseudo-replication; p-values too small). The "Current Technical Replicate Assignments" summary below the button shows
the result either way.

NEWS fact (AUDIT-KIT-256-impl-20261126 I03, the coordinator's correction): "in this and every earlier release" became
"in every earlier release with the design builder". The kept and re-sent selection dates from 33350d50 (2025-07-09);
the tags v0.2 (2025-05-14), v0.73, v0.75, v0.80 (2024) and v0.9 (2025-04-09) have no technical-replicate box. The final words stay the owner's (VAL-41).
