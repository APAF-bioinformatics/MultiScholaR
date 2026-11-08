# KIT-106 entry-point suite

## Simulation Setup (run 2): two lanes alone on 9152e613

The lanes ran through `kit-105/lifecycle/lifecycle.R` as committed, with strict skips and `R_LIBS` kit_110. Logs: [lifecycle/](lifecycle/).

- **metab_combined:** identical to its re-frozen oracle, 89 s. Traced calls, exactly as the ticket expected:
  - `mod_omic_da_server` 1;
  - `runMetabReferenceDaAnalysis` 1;
  - `runWorkflowOmicsPublicationInProcess` 1 (the DA publication);
  - `loadWorkflowOmicsDaSession` 2 (the refusal and the load);
  - `workflowOmicsInProcessVolcanoPageFactory` 0 (the box stays off).
- **prot_dia:** identical to its oracle, 327 s. Traced calls:
  - `mod_omic_da_server` 1 (W5's stop rule did not fire);
  - `da_server_load_session_handler` 1 and `da_server_run_analysis_handler` 1 (the plugins);
  - `runDIAReferenceDaAnalysis` 1;
  - `validateProtDaFilteredSession` 1 and `normaliseProtDaContrastsTable` 3 (the moved validators).

Both first attempts ran cleanly, with no instrumentation failure.

## Run 3: the four entry files in parallel on 9152e613

| File | Lanes | Pass | Fail | Error | Skip | Elapsed |
|---|---|---|---|---|---|---|
| test-gui-entry-transcriptomics.R | trx_screen | 29 | 0 | 0 | 0 | 73 s |
| test-gui-entry-lipidomics.R | lipid_canonical | 35 | 0 | 0 | 0 | 134 s |
| test-gui-entry-metabolomics.R | metab_lc, metab_gc, metab_combined | 105 | 0 | 0 | 0 | 300 s |
| test-gui-entry-proteomics.R | prot_dia, prot_dia_limpa, prot_tmt, prot_lfq, prot_lfq_fragpipe | 187 | 0 | 0 | 0 | 1,104 s |

**Decision:** all 10 lanes are identical to the oracles. Four oracles were re-frozen for D3 at W3 and W4 (lipid_canonical, metab_lc, metab_gc, metab_combined), once each, with attribution ([lane-changes.md](lane-changes.md)). The others are unchanged since KIT-110.
