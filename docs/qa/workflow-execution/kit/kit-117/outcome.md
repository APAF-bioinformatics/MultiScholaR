# KIT-117: harness labels after the GUI workover, protocol v1.10 (2026-11-13)

**Outcome:** done. Labels and one comment, no rule, oracle or step changed. Authority: the user, 2026-11-13 (autonomy for the GUI journey).

| Item | Result | Commit |
|---|---|---|
| P-22, P-24 | The trx step table names `mod_omic_<tab> (transcriptomics)` (KIT-109's shared modules) in place of the lipid modules under `trxOmicProfile()` and the deleted trx modules | `9c8c23eb` |
| P-29 | The DIA step map's enrichment row names `mod_omic_enrich.R` and the core (`utils_workflow_enrichment.R`) in place of the deleted `mod_prot_enrich_*` helpers | `9c8c23eb` |
| Driver comment | The volcano-tick comment says the tick is a no-op since R5 (the volcano is on by default); the tick is unchanged | `9c8c23eb` |
| KIT-113 names | No change needed: the generator's KIT-113 block names no function | |
| P-27 (`GeneID`) | Parked: changes an output table's columns (`parked-20261113.md` row 3) | |
| Protocol | v1.10, `protocol.md` section 12 | `9c8c23eb` |
| Bindings | regenerated under v1.10 | `cb9a0dc0` |

## Checks

- Harness tests at KIT-120's exit counts: `test-workflow-kit-harness.R` 735/0, `test-workflow-checkpoint-harness.R` 641/0, `test-e2e-browser-harness.R` 123/0.
- **One member per route on kit_933, under the regenerated bindings:** 117a (lipid), 117b (metab), 117c (trx), 117m (DIA_LIMPA), 117l (LFQ), 117t (TMT) accepted, exact. 117d (DIA, stopped after design) has the DIA journey's known mismatch set; store against store with 933u: 2 states, 15 comparisons, 0 different.

## Recorded for a follow-up (not changed)

- **Stale module names in the step maps** (files the GUI workover moved; the maps are documentation, read by no code): `step-map.json`: `mod_prot_summary_report_helpers`, `mod_prot_summary_server_helpers`; `step-map-limpa.json`: the same two; `step-map-lipid.json`: `mod_lipid_da_analysis_helpers`, `mod_lipid_da_export_helpers`, `mod_lipid_da_session_helpers`, `mod_lipid_design_builder_action_helpers`, `mod_lipid_design_builder_helpers`, `mod_lipid_import_processing_helpers`, `mod_lipid_norm_workflow_helpers`, `mod_lipid_qc_duplicates_server_helpers`, `mod_lipid_qc_intensity_server_helpers`, `mod_lipid_summary_state_helpers`, `mod_lipid_summary_workflow_helpers`; `step-map-metab.json`: `mod_metab_da_observer_helpers`, `mod_metab_da_registration_helpers`, `mod_metab_import_processing_helpers`, `mod_metab_norm_entry_helpers`, `mod_metab_norm_pipeline_helpers`, `mod_metab_norm_render_helpers`, `mod_metab_norm_session_helpers`, `mod_metab_qc_duplicates_server_helpers`, `mod_metab_summary_observer_helpers`. Line numbers in the maps are anchors as of the ticket that wrote them.
- **The harness's `reclaim_memory` launch option is dead** (nothing in `R/` reads it; KIT-121 turned the reclaim on for every app owner).
