### Tests (head): 6 files; passed 2716, failed 206, errors 0, skipped 0; ok 2, failed 4, timeout 0, crashed 0

### Tests (base): 6 files; passed 2664, failed 204, errors 0, skipped 0; ok 1, failed 4, timeout 0, crashed 0, absent 1

### Base versus head: same 2, new-failure 1, message-changed 1, count-changed 1, fixed 0, only-on-head 1, compare-error 0
**Not clean under the rule**: see the rows below and compare-failures.txt

| file | base p/f/e/s | head p/f/e/s | verdict | tests | note |
|---|---|---|---|---|---|
| test-kit268-bh-fallback-notice.R | 0/0/0/0 absent | 30/0/0/0 ok | only-on-head |  |  |
| test-omic-da-module.R | 1288/81/0/0 failed | 1286/83/0/0 failed | new-failure | new-failure: GUI audit A112-03: a managed run whose contrast fell back to BH shows the in-process run's fallback warning; counts: GUI audit A112-03: a managed run whose contrast fell back to BH shows the in-process run's fallback warning (passed 5 -> 3, skipped 0 -> 0) |  |
| test-omic-enrich-module.R | 1011/59/0/0 failed | 1011/59/0/0 failed | same |  | identical message without the actual digest value: (KIT-113 W1) each W1 (b) session's committed record, results, product, provenance and analysis record step |
| test-kit236-report-sentences.R | 289/0/0/0 ok | 313/0/0/0 ok | count-changed | counts: KIT-268 T4: the FDR sentence of each report adds the scale sentence after the named Benjamini-Hochberg comparisons (passed NA -> 24, skipped NA -> 0) |  |
| test-gui-entry-lipidomics.R | 19/16/0/0 failed | 19/16/0/0 failed | message-changed | message-changed: KIT-101 entry lane lipid_canonical: two LipidSearch assays |  |

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

**Run outcome: success** (every head file ran to completion; failing expectations, if any, are in the tables)
