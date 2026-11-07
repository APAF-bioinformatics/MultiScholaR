# KIT-103 findings

Append only. Ids are F103-NN. A fix that does not change the methodology is made here, with a test that fails first. Anything else is parked.

| Id | Where | What | Status |
|---|---|---|---|
| F103-01 | lipid and metab "Import Existing Design" (W2) | In the app, the import fails for both omics. When the import holds vroom's tables (`setWorkflowOmicsBusTables`), the owner rejects vroom's live parser pointer, the `problems` attribute ("payload contains a nonportable resource"), after `config_list`, `column_mapping`, `design_matrix` and `contrasts_tbl` have already been written. The MCI tests inject `read.delim`, so none caught it. | **fixed** `fe1a9742`: `readOmicsDesignTable` drops only that attribute, so values and column types stay vroom's. The "W2 finding" test block fails on today's code. |
| F103-02 | metab design import, error path (W2) | Today's metab import writes `config_list` and `column_mapping` while it is still reading, so a later read failure leaves them replaced. The core reads everything first and writes the bus only then, as lipid already does. Notification texts are identical. This becomes visible when W5 moves metab onto the core. | accepted (error path only; no recorded state) |
| — | W2 | Two kept-as-today differences are now explicit: `design_import$record_design_file` (lipid TRUE, metab FALSE) and the annotation-column pattern per omic in `omicsDesignRules`. | recorded |
