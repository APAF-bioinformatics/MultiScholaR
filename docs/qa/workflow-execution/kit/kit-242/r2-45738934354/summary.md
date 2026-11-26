### Tests (head): 190 files; passed 36444, failed 1023, errors 3, skipped 272; ok 165, failed 25, timeout 0, crashed 0

### Tests (base): 190 files; passed 36280, failed 1115, errors 7, skipped 272; ok 162, failed 28, timeout 0, crashed 0

### Base versus head: same 184, new-failure 0, message-changed 3, count-changed 0, fixed 3, only-on-head 0, compare-error 0
**Not clean under the rule**: see the rows below and compare-failures.txt

| file | base p/f/e/s | head p/f/e/s | verdict | tests | note |
|---|---|---|---|---|---|
| test-kit242-normalisation-fallback.R | 370/84/4/0 failed | 526/0/0/0 ok | fixed | fixed: KIT-242 E14: lipid legacy RUV-III-C with no control found in A2 keeps A2 uncorrected, with a classed warning; fixed: KIT-242 E15: lipid legacy RUV-III-C with no control in any assay keeps both, each with a classed warning; fixed: KIT-242 E16: lipid legacy RUV-III-C with k = -1 for A2 stops with omics_ruv_failed (invalid_k); fixed: KIT-242 E18: lipid legacy Automatic RUV corrects A2 without  |  |
| test-module-ci-metab-norm.R | 154/4/0/0 failed | 158/0/0/0 ok | fixed | fixed: MCI-015.3 RUV matrix covers skip, manual, automatic, invalid k, controls, and per-assay payloads; counts: MCI-015.3 RUV matrix covers skip, manual, automatic, invalid k, controls, and per-assay payloads (passed 18 -> 22, skipped 0 -> 0) |  |
| test-omic-enrich-module.R | 1011/59/0/0 failed | 1011/59/0/0 failed | same |  | identical message without the actual digest value: (KIT-113 W1) each W1 (b) session's committed record, results, product, provenance and analysis record step |
| test-omic-norm-module.R | 583/56/0/0 failed | 583/56/0/0 failed | message-changed | message-changed: (b) lipid, RUV manual with a manual ITSD selection: run, skip correlation, reset; message-changed: (b) metab, RUV manual with a manual ITSD selection: run, skip correlation, reset; message-changed: (b, c) lipid, RUV automatic: run with a figure missing at the composite, final QC plot, the final PCA failing, skip-correlation with the commit failing; message-changed: (b, c) metab, R |  |
| test-workflow-omics-binding-steps.R | 57/2/0/0 failed | 57/2/0/0 failed | message-changed | message-changed: KIT-094: the characterisation of the lists the bindings own equals the snapshot taken on the base commit |  |
| test-workflow-omics-declarations.R | 43/25/0/0 failed | 43/25/0/0 failed | message-changed | message-changed: KIT-073: the omics declarations and their journeys equal the frozen dump of the two files the core replaced |  |
| test-workflow-omics-worker.R | 261/4/0/0 failed | 265/0/0/0 ok | fixed | fixed: KIT-078: every omics compute declaration runs from a frozen request through the worker's entry point; counts: KIT-078: every omics compute declaration runs from a frozen request through the worker's entry point (passed 92 -> 96, skipped 0 -> 0) |  |

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
