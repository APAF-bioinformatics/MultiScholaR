# KIT-106 findings

Append only. Ids are F106-NN. A fix that does not change the methodology is made here, with a test that fails first. Anything else is parked.

| Id | Where | What | Status |
|---|---|---|---|
| F106-01 | the five proteomics entry lanes (`helper-gui-entry.R`) | The DA run updates the contrast selector's choices without selecting one, so every proteomics DA render in the lanes is `<req>`. Proteomics' DA renders have no coverage in the entry lanes, e2e or journeys. | Recorded (ticket). W1's `prot_main` captures, with a contrast selected, are the render proof W5 rests on. |
| F106-02 | `mod_prot_da_handlers_run_analysis.R:60` and `:72` | An in-process proteomics DA run returns before removing its "Running differential abundance analysis..." notice (id `da_working`, duration NULL). Managed sessions remove it (`:49-50`). In-process only. Display only. | **Queued to KIT-111.** Fixing it here re-freezes the five proteomics lanes, which are W5's stop-rule measure. KIT-111 moves proteomics onto the frame and can make one attributed proteomics re-freeze for F106-02, F106-03 and P-25 together. Frozen as today's behaviour in W1 (`prot_main`, `prot_load`). |
| F106-03 | `mod_prot_da_render_heatmap.R:104-105` | In its default state the proteomics cluster summary prints "No clusters defined. ... heatmap.NULL" (`cat()` then `return()` inside `renderPrint`). Display only. | **Queued to KIT-111**, with F106-02, for the same reason. |
| F106-04 | `mod_prot_da_render_table.R:151-160` | Proteomics' "Download All Results" writes only a placeholder line into `DA_results_<date>.zip`. | **Parked (P-26).** A real export is a new feature, not a fix. |
| F106-05 | lipid and metab Run DA | On an owner whose head is `initial`, Run DA raises inside the observer ("state requires explicit legacy conversion") instead of giving the no-data refusal. The shell's tab gating should keep this from being reached. Seen in a W1 scratch probe; not frozen. | **Assigned to W3 (lipid) and W4 (metab):** the frame's run handler refuses with the no-data notification, with a test that fails first. No result changes. |
