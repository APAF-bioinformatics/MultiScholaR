### Tests (head): 106 files; passed 25656, failed 628, errors 4, skipped 215; ok 87, failed 19, timeout 0, crashed 0

### Tests (base): 106 files; passed 25521, failed 624, errors 4, skipped 215; ok 87, failed 18, timeout 0, crashed 0, absent 1

### Base versus head: same 101, new-failure 1, message-changed 3, count-changed 0, fixed 0, only-on-head 1, compare-error 0
**Not clean under the rule**: see the rows below and compare-failures.txt

| file | base p/f/e/s | head p/f/e/s | verdict | tests | note |
|---|---|---|---|---|---|
| test-prot-norm-core.R | 1548/306/0/0 failed | 1548/306/0/0 failed | message-changed | message-changed: the lipid and metab RUV declaration give the recorded result through the moved owner stripper |  |
| test-omic-enrich-module.R | 1011/59/0/0 failed | 1011/59/0/0 failed | same |  | identical message without the actual digest value: (KIT-113 W1) each W1 (b) session's committed record, results, product, provenance and analysis record step |
| test-ruv-control-anova-grouping.R | 0/0/0/0 absent | 139/0/0/0 ok | only-on-head |  |  |
| test-workflow-omics-binding-steps.R | 58/1/0/0 failed | 58/1/0/0 failed | message-changed | message-changed: KIT-094: the characterisation of the lists the bindings own equals the snapshot taken on the base commit |  |
| test-workflow-omics-ruv-cancor-once.R | 16/0/0/0 ok | 14/2/0/0 failed | new-failure | new-failure: metab sweep is 1:1 with imputing at every percentage, with missing values; new-failure: lipid sweep is 1:1 with imputing at every percentage, with missing values; counts: lipid sweep is 1:1 with imputing at every percentage, with missing values (passed 3 -> 2, skipped 0 -> 0); counts: metab sweep is 1:1 with imputing at every percentage, with missing values (passed 3 -> 2, skipped 0 - |  |
| test-workflow-omics-declarations.R | 45/23/0/0 failed | 43/25/0/0 failed | message-changed | message-changed: KIT-073: the omics declarations and their journeys equal the frozen dump of the two files the core replaced; counts: KIT-073: the omics declarations and their journeys equal the frozen dump of the two files the core replaced (passed 36 -> 34, skipped 0 -> 0) |  |

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
| 1 | docs/qa/workflow-execution/kit/kit-190/norm-dump.R | 0 | 379 | ok |

**Run outcome: success** (every head file ran to completion; failing expectations, if any, are in the tables)
