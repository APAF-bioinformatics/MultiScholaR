# KIT-103 run order per omic (D3)

Run order is the order of the runs in a fresh design: a design the builder builds from the imported data when there is no existing design. It is set per omic in the configuration (`omicGuiConfig(<omic>)$design$run_order` and `$run_order_rule`). The field is required and has no shared default. This file and the configuration record the rule. The `design_builder` provenance does not record it (P-08, parked: adding it would change every design state's record and the oracles).

Nothing here changed in KIT-103. Each rule is today's (3c657c14), and it is applied only where it was applied before.

| Omic | `run_order_rule` | `run_order` | Where it is applied | Where it was at 3c657c14 |
|---|---|---|---|---|
| lipidomics | `mixedsort` | `gtools::mixedsort` | `buildOmicsDesignInitialState` (func_omics_design_shared.R), fresh design only | `mod_lipid_design_builder_state_helpers.R:298, :395` |
| metabolomics | `mixedsort` | `gtools::mixedsort` | `buildOmicsDesignInitialState`, fresh design only | `mod_metab_design_builder_server.R:237` |
| transcriptomics | `input` | `identity` (the counts file's column order) | `buildOmicsDesignInitialState`, fresh design only | `trxOmicProfile()$run_order` (`mod_transcriptomics.R:39-41`) |
| proteomics | `mixedsort` | `gtools::mixedsort` | `buildProtDesignInitialState` and `buildDIAReferenceDesignState` (utils_workflow_dia_state.R), both initial-state paths of the `dia_long` adapter | `mod_prot_design_builder_state_helpers.R:348` |

For proteomics, the configuration value records the rule. The DIA core applies `gtools::mixedsort` itself, in `buildProtDesignInitialState` (utils_workflow_dia_state.R), as before W3. `buildDIAReferenceDesignState` calls it for the DIA reference path. The core does not read `config$design$run_order`.

These rules stay as they were:
- **Imported designs:** keep the file's order for every omic. The run-order rule is not applied.
- **Save:** the saved assay columns follow the design (`intersect`).
- **Dropdowns:** the builder's own `mixedsort` of its sample lists (`mod_omic_design_builder.R`) is display only.

## Where each rule is pinned in tests

- **The configuration values**, all four omics: `test-omic-gui-config.R`, "the design group holds today's values per omic". This block checks the `run_order` function, the `run_order_rule` and each one's effect on `c("S10", "S2", "S1")`.
- **Required, no default:** `test-omic-gui-config.R`, "run_order is required with no default; each required design field is refused by class". Each omic's `run_order` is removed and then mistyped, and the validator refuses both.
- **Lipid, metab and trx in the core:** `test-omics-design-core.R`, "C3 the run order per omic: mixedsort for lipid and metab, input order for trx, file order for an imported design". It also checks that an imported design keeps the file's order for each omic, that a missing `run_order` is refused (`omics_design_config_invalid`), and that the configured function is the one applied.
- **Saved columns:** `test-omics-design-core.R`, "C5 the saved assay columns follow the design's run order".
- **Lipid and trx in the shared builder:** `test-omic-design-module.R`, "W4 the builder keeps its reactive `design_matrix`…". It runs the lipid wrapper with no profile (mixedsort) and with `trxOmicProfile()` (input order).
- **Trx in the module CI:** `test-module-ci-trx-design.R`, MCI-027.1, "the builder keeps the counts file's run order and names".
- **Proteomics:** `test-prot-design-core.R`, "run order: a fresh proteomics design is mixedsort on both initial-state paths". It covers the plain table and the DIA reference. The applet block of `test-applets-direct-shared.R` also checks that the fresh design's sample inputs come in mixedsort order.
- **End to end:** the 10 entry lanes' design steps (`tests/testdata/gui-entry/`) are byte-identical to the oracles. These include the design matrix's run order.
