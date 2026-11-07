# KIT-102 findings

Append only. Ids are F102-NN. A fix that does not change the methodology is made here, with a test that fails first. Anything else is parked in `../gui-workover-parked-decisions.md`.

| Id | Where | What | Status |
|---|---|---|---|
| F102-01 | metab import validation order (W2) | Today's metab import validates before converting the sample columns, so `validateMetabColumnMapping`'s "High proportion of missing values" warning is computed on text: an empty text cell counts as neither NA nor 0. The shared core converts, then checks, then validates, so it counts such a cell as missing. Only the validation warning changes. The tables, artifacts, refusals and status are unchanged, and today the warning appears only in the returned result, never in the GUI. | accepted as the shared order (`d469cea3`; test "F102-01" in `test-workflow-omics-import.R`) |
| — | parked | P-07 (W3): the undefined `clean_acc` breaks the mixed-species organism filter. | parked |
