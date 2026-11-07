# KIT-103 findings

Append only. Ids are F103-NN. A fix that does not change the methodology is made here, with a test that fails first. Anything else is parked.

| Id | Where | What | Status |
|---|---|---|---|
| F103-01 | lipid and metab "Import Existing Design" (W2) | In the app, the import fails for both omics. When the import holds vroom's tables (`setWorkflowOmicsBusTables`), the owner rejects vroom's live parser pointer, the `problems` attribute ("payload contains a nonportable resource"), after `config_list`, `column_mapping`, `design_matrix` and `contrasts_tbl` have already been written. The MCI tests inject `read.delim`, so none caught it. | **fixed** `fe1a9742`: `readOmicsDesignTable` drops only that attribute, so values and column types stay vroom's. The "W2 finding" test block fails on today's code. |
| F103-02 | metab design import, error path (W2) | Today's metab import writes `config_list` and `column_mapping` while it is still reading, so a later read failure leaves them replaced. The core reads everything first and writes the bus only then, as lipid already does. Notification texts are identical. This becomes visible when W5 moves metab onto the core. | accepted (error path only; no recorded state) |
| — | W2 | Two kept-as-today differences are now explicit: `design_import$record_design_file` (lipid TRUE, metab FALSE) and the annotation-column pattern per omic in `omicsDesignRules`. | recorded |
| F103-03 | standalone design applet (W6) | Without a DIA owner, the proteomics builder cannot rename or save: the in-memory path has been retired since KIT-048. The applet's `data_tbl` also has no Run column. This is the same before and after W6. | parked as P-15 (the applet acceptance item's "returns a save result" is unmet as written) |
| — | W6 | Kept: proteomics' metadata assignment stores integer replicates (`seq.int`); the shared C4 edit stores doubles. The dia_long adapter keeps proteomics' integer path, so the recorded design states are unchanged. | recorded |
