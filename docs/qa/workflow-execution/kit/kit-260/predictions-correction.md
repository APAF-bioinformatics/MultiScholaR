# KIT-260 correction pass: predictions, written before the runs

Head: a516e500 (kit-260 after the correction commits and the merge of staging 42846c0e).

## R1b (the new test file only; base 7af2d51e = the correction's tests commit, before the I02 (a) fix; head a516e500)

- I02 (record, study parameters, report; metab and lipid): **fail on base** (no note: the product names no rows and the
  record had no other source), **pass on head**. The twin cases pass on both.
- I10 (lipid study parameters and report on a run at this version; the integration-tab note; the model record's marked
  rows): pass on both sides (the behaviour exists since 288abf49; these are coverage tests).
- E6's managed arm through the worker's delivery (added at b71b2d62, after the base): only on head; passes.
- I05: on base the test's record step stubs `omicsItsdRowsFromStates`, which does not exist there, so it errors on base.
  On head the result is what the test is for. By reading (INFERRED, not predicted with confidence): the replay harness
  judges every product by sha256 and the replayed ITSD product carries the rows-used fields the older product lacks, so a
  mismatch on `metab_itsd_norm/products/metab_itsd_normalisation_science` (and possibly later states whose products
  chain to it) is the likely outcome; the assays are predicted identical. **If the harness reports any mismatch, the
  coordinator's ruling applies: stop and report (a replay marker would be needed); nothing is added here.**
- Every other row: as R1 (pass on both or as before).

## R2b (130 files, base staging 42846c0e, head a516e500, R CMD check both sides)

- `test-kit260-itsd-default-pattern.R`: only on head; every test passes except, possibly, I05 (above).
- The six files with approved stored-expectation changes: `same` (the base runs staging's own expectations on staging's
  code, the head the new ones on the new code), except `test-omics-norm-core.R`: `count-changed` (one more passing
  expectation on head).
- The eight files that fail on staging (KIT-260's R2 list: `test-omic-summary-module.R`, `test-omics-summary-core.R`,
  `test-prot-norm-core.R`, `test-workflow-omics-declarations.R`, `test-workflow-kit-harness.R`,
  `test-workflow-operation-registry.R`, `test-workflow-omics-binding-steps.R`, `test-workflow-omics-owner.R`, as far as
  staging has not changed them): `same`.
- Every other file: `same`. No `new-failure`, no `message-changed`.
- Check: findings identical on both sides; nothing only on head.
