### Tests (head): 135 files; passed 24418, failed 491, errors 3, skipped 287; ok 123, failed 12, timeout 0, crashed 0

### Tests (base): 135 files; passed 24238, failed 671, errors 3, skipped 287; ok 121, failed 14, timeout 0, crashed 0

### Base versus head: same 133, new-failure 0, message-changed 0, count-changed 0, fixed 2, only-on-head 0, compare-error 0
**Clean under the rule** (no new failure, no changed failure message)

| file | base p/f/e/s | head p/f/e/s | verdict | tests | note |
|---|---|---|---|---|---|
| test-kit247-design-resave.R | 170/16/0/0 failed | 186/0/0/0 ok | fixed | fixed: KIT-247 W288a (VAL-331, edited by KIT-277): prot_lfq a formula-only re-save asks, sets the differential analysis back and the next run fits the new formula; fixed: KIT-247 W288a (VAL-331, edited by KIT-277): lipid_canonical a formula-only re-save asks, sets the differential analysis back and the next run fits the new formula; fixed: KIT-247 W288b (VAL-331, edited by KIT-277): a proteomics I |  |
| test-kit277-resave-carried-through.R | 88/164/0/1 failed | 252/0/0/1 ok | fixed | fixed: KIT-277 F1/F2 (S1): prot_lfq a formula-only re-save asks, and Save puts the formula on the DA tab and sets the DA back; fixed: KIT-277 F5, F6 (S1): prot_lfq the same formula saved again asks nothing; the notices say what was done; fixed: KIT-277 F1/F2 (S1): lipid_canonical a formula-only re-save asks, and Save puts the formula on the DA tab and sets the DA back; fixed: KIT-277 F5, F6 (S1):  |  |
| test-omic-enrich-module.R | 1011/59/0/0 failed | 1011/59/0/0 failed | same |  | identical message without the actual digest value: (KIT-113 W1) each W1 (b) session's committed record, results, product, provenance and analysis record step |

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
