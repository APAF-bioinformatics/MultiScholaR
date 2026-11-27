### Tests (head): 19 files; passed 6935, failed 55, errors 0, skipped 3; ok 16, failed 3, timeout 0, crashed 0

### Tests (base): 19 files; passed 6530, failed 55, errors 0, skipped 3; ok 15, failed 3, timeout 0, crashed 0, absent 1

### Base versus head: same 15, new-failure 0, message-changed 0, count-changed 3, fixed 0, only-on-head 1, compare-error 0
**Clean under the rule** (no new failure, no changed failure message)

| file | base p/f/e/s | head p/f/e/s | verdict | tests | note |
|---|---|---|---|---|---|
| test-kit227-report-ruv-sentences.R | 162/0/0/0 ok | 186/0/0/0 ok | count-changed | counts: KIT-227 G1: Automatic on the group column prints the base text byte for byte, without and with a record (passed 42 -> 66, skipped 0 -> 0) |  |
| test-kit235-lipid-metab-report-parser.R | 54/0/0/0 ok | 58/0/0/0 ok | count-changed | counts: KIT-235 T6: '•' and '-' files print the base text, and the same text as the '*' file (M, L1) (passed 20 -> 24, skipped 0 -> 0) |  |
| test-kit236-report-sentences.R | 293/0/0/0 ok | 297/0/0/0 ok | count-changed | counts: KIT-236 G1: the Methods paragraphs outside section 1 print the base text, with and without a record (passed 24 -> 28, skipped 0 -> 0) |  |
| test-kit272-prot-box-panel-text.R | 0/0/0/0 absent | 373/0/0/0 ok | only-on-head |  |  |

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
