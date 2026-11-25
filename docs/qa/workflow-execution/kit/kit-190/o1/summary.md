### Tests (head): 1 files; passed 16, failed 0, errors 0, skipped 0; ok 1, failed 0, timeout 0, crashed 0

### Tests (base): 1 files; passed 14, failed 2, errors 0, skipped 0; ok 0, failed 1, timeout 0, crashed 0

### Base versus head: same 0, new-failure 0, message-changed 0, count-changed 0, fixed 1, only-on-head 0, compare-error 0
**Clean under the rule** (no new failure, no changed failure message)

| file | base p/f/e/s | head p/f/e/s | verdict | tests | note |
|---|---|---|---|---|---|
| test-workflow-omics-ruv-cancor-once.R | 14/2/0/0 failed | 16/0/0/0 ok | fixed | fixed: metab sweep is 1:1 with imputing at every percentage, with missing values; fixed: lipid sweep is 1:1 with imputing at every percentage, with missing values; counts: lipid sweep is 1:1 with imputing at every percentage, with missing values (passed 2 -> 3, skipped 0 -> 0); counts: metab sweep is 1:1 with imputing at every percentage, with missing values (passed 2 -> 3, skipped 0 -> 0) |  |

**Run outcome: success** (every head file ran to completion; failing expectations, if any, are in the tables)
