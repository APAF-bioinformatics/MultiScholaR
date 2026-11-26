### Tests (head): 3 files; passed 1605, failed 100, errors 0, skipped 0; ok 0, failed 3, timeout 0, crashed 0

### Tests (base): 3 files; passed 1577, failed 81, errors 0, skipped 0; ok 1, failed 1, timeout 0, crashed 0, absent 1

### Base versus head: same 0, new-failure 2, message-changed 0, count-changed 0, fixed 0, only-on-head 1, compare-error 0
**Not clean under the rule**: see the rows below and compare-failures.txt

| file | base p/f/e/s | head p/f/e/s | verdict | tests | note |
|---|---|---|---|---|---|
| test-kit268-bh-fallback-notice.R | 0/0/0/0 absent | 12/11/0/0 failed | only-on-head |  |  |
| test-kit236-report-sentences.R | 289/0/0/0 ok | 307/6/0/0 failed | new-failure | new-failure (test only on head): KIT-268 T4: the FDR sentence of each report adds the scale sentence after the named Benjamini-Hochberg comparisons; counts: KIT-268 T4: the FDR sentence of each report adds the scale sentence after the named Benjamini-Hochberg comparisons (passed NA -> 18, skipped NA -> 0) |  |
| test-omic-da-module.R | 1288/81/0/0 failed | 1286/83/0/0 failed | new-failure | new-failure: GUI audit A112-03: a managed run whose contrast fell back to BH shows the in-process run's fallback warning; counts: GUI audit A112-03: a managed run whose contrast fell back to BH shows the in-process run's fallback warning (passed 5 -> 3, skipped 0 -> 0) |  |

**Run outcome: success** (every head file ran to completion; failing expectations, if any, are in the tables)
