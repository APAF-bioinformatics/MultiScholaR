### Tests (head): 97 files; passed 24754, failed 682, errors 3, skipped 92; ok 79, failed 18, timeout 0, crashed 0

### Tests (base): 97 files; passed 24684, failed 752, errors 3, skipped 92; ok 78, failed 19, timeout 0, crashed 0

### Base versus head: same 95, new-failure 0, message-changed 0, count-changed 0, fixed 2, only-on-head 0, compare-error 0
**Clean under the rule** (no new failure, no changed failure message)

| file | base p/f/e/s | head p/f/e/s | verdict | tests | note |
|---|---|---|---|---|---|
| test-kit276-ruv-superseded.R | 125/70/0/0 failed | 194/1/0/0 failed | fixed | fixed: KIT-276 T1: pass-through, correlation, RUV covariates, correlation again leaves no RUV step; the DA used none; fixed: KIT-276 T3: covariates, correction, covariates leaves the covariates step as the only RUV step; fixed: KIT-276 T4: Skip RUV after covariates takes the RUV result out: head, history, record, DA, texts, Figure 1; fixed: KIT-276 T5: Skip RUV after the correction mode: the DA re |  |
| test-omic-enrich-module.R | 1011/59/0/0 failed | 1011/59/0/0 failed | same |  | identical message without the actual digest value: (KIT-113 W1) each W1 (b) session's committed record, results, product, provenance and analysis record step |
| test-workflow-integration-input.R | 240/1/0/0 failed | 241/0/0/0 ok | fixed | fixed: KIT-061 (K061-1): DIA_LIMPA switched covariates -> correction -> covariates writes the head's adjusted product; counts: KIT-061 (K061-1): DIA_LIMPA switched covariates -> correction -> covariates writes the head's adjusted product (passed 9 -> 10, skipped 0 -> 0) |  |

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
