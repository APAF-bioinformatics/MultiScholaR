### Tests (head): 25 files; passed 4598, failed 13, errors 0, skipped 12; ok 22, failed 3, timeout 0, crashed 0

### Tests (base): 25 files; passed 4477, failed 86, errors 2, skipped 12; ok 22, failed 3, timeout 0, crashed 0

### Base versus head: same 24, new-failure 1, message-changed 0, count-changed 0, fixed 0, only-on-head 0, compare-error 0
**Not clean under the rule**: see the rows below and compare-failures.txt

| file | base p/f/e/s | head p/f/e/s | verdict | tests | note |
|---|---|---|---|---|---|
| test-integration-removed-factors.R | 368/77/2/3 failed | 489/4/0/3 failed | new-failure | fixed: KIT-184 K8: too few samples for the design, a factor that is a combination of the removed factors, an approximate input; fixed: KIT-184 K14: covariates change no count; run_complete and the factor's share by hand; fixed: KIT-184 K15: the test of record is asked with the samples left after one per removed factor; fixed: KIT-184 K11: a result's sentences: the form's label, the note for the mo |  |

### R CMD check, base versus head
```
base: ERRORs 0, WARNINGs 3, NOTEs 3 (check exit 0)
head: ERRORs 0, WARNINGs 3, NOTEs 3 (check exit 0)
findings identical on both: 6; only on head: 0; only on base: 0

== only on head (new)
none

== only on base (gone)
none
```

### Scripts
| # | path | exit | seconds | status |
|---:|---|---:|---:|---|
| 1 | docs/qa/workflow-execution/kit/kit-184/base-store-job.R | 0 | 2096 | ok |

**Run outcome: success** (every head file ran to completion; failing expectations, if any, are in the tables)
