### Tests (head): 54 files; passed 10325, failed 232, errors 0, skipped 47; ok 47, failed 7, timeout 0, crashed 0

### Tests (base): 54 files; passed 10309, failed 248, errors 0, skipped 47; ok 46, failed 8, timeout 0, crashed 0

### Base versus head: same 53, new-failure 0, message-changed 0, count-changed 0, fixed 1, only-on-head 0, compare-error 0
**Clean under the rule** (no new failure, no changed failure message)

| file | base p/f/e/s | head p/f/e/s | verdict | tests | note |
|---|---|---|---|---|---|
| test-kit228-integration-input-refused-save.R | 163/16/0/0 failed | 179/0/0/0 ok | fixed | fixed: KIT-228 E2, G4: transcriptomics, the RUVIII export missing at the gate removes the earlier input; with none, nothing is said; fixed: KIT-228 E3, E5, E6, G1, G2, G4: lipid a Save that writes no new input removes the earlier one; fixed: KIT-228 E3, E5, E6, G1, G2, G4: metab a Save that writes no new input removes the earlier one; counts: KIT-228 E2, G4: transcriptomics, the RUVIII export miss |  |

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
