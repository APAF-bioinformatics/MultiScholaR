### Tests (head): 53 files; passed 9990, failed 232, errors 0, skipped 47; ok 45, failed 7, timeout 1, crashed 0

### Tests (base): 53 files; passed 9873, failed 232, errors 0, skipped 47; ok 44, failed 7, timeout 1, crashed 0, absent 1

### Base versus head: same 52, new-failure 0, message-changed 0, count-changed 0, fixed 0, only-on-head 1, compare-error 0
**Clean under the rule** (no new failure, no changed failure message)

| file | base p/f/e/s | head p/f/e/s | verdict | tests | note |
|---|---|---|---|---|---|
| test-kit228-integration-input-refused-save.R | 0/0/0/0 absent | 117/0/0/0 ok | only-on-head |  |  |
| test-kit247-design-resave.R | 0/0/0/0 timeout | 0/0/0/0 timeout | same |  | timeout on both |

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

**Run outcome: failure** (a head file crashed or timed out, or a script failed)
