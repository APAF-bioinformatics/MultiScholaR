# KIT-202 outcome (R2 judged; the unpredicted move decided by VAL-189; R3 next)

Ticket `tickets/kit/KIT-202.md` r1.2 (owner decisions VAL-130). Branch `kit-202-intensity-filter-samples` from
`v0.9-staging` 9453e79b.

## Commits

- f9d7d910: tests only (nothing under `R/` or `man/`): `tests/testthat/test-kit202-intensity-filter-samples.R` (T1 to
  T7), `tests/testthat/test-metab-qc-shared.R` (T8, D7), `kit-202/filter-samples-control.R` (S1, runs in R2).

## R1: run 45326839837 on f9d7d910 (tests only, base none, 1 job, 10.8 job-minutes, ceiling 45)

Observed (`r1/tests.tsv`, `r1/failures.txt`, the two logs):
- `test-kit202-intensity-filter-samples.R`: 123 passed, 32 failed, 0 errors. `test-metab-qc-shared.R`: 18 passed, 2
  failed, 0 errors. Every failure is an assertion; none is an error.
- Values at which the tests fail, each as the ticket states for the base:
  - T1 (a) metabolite seed: threshold 4.9 (head 16410), kept 1 to 8 (head 1, 2, 3, 4, 7, 8). G0's value.
  - T1 (b) lipid: threshold 2.2 (head 100.9), kept {3} (head {2, 3}). The kept set was hand-computed in the ticket.
  - T2, both variants, both omics: kept as T1's base sets.
  - T3 (a): the helper refuses `sample_columns` (result NULL), both omics.
  - T3 (b): no `omics_intensity_filter_samples_guessed` warning; the metabolite case keeps 5 to 8 (the diff prints
    `actual[1:3]` 5 6 7 against `expected[1:7]`), the lipid case keeps {3}. Hand values.
  - T5 (a): no `omics_intensity_filter_samples_invalid` error from either S4 method (absent and text sample) or from the
    partitioned step (text sample).
  - T6 (i): the marker row is absent. T6 (iii) and (v): no marker in the record (in process; managed, lipid and metab).
  - T6 (iv): non-vacuity fails (both rules give the same kept set on base), and record (b) and its replay carry no
    marker. The replay itself has no mismatch on either owner, and record (a) carries no marker, on base.
  - T7: kept 3 to 8 (head 1, 2, 3, 4, 7, 8). Hand value.
  - T8: the mocked helper receives `sample_columns` NULL (head c("S1", "S2")), and the text-sample assay is not refused.
- Passing on base, as the ticket expects: T1 (b)'s guard at 0.5, T4 (six unaffected objects against the oracle at
  percentiles 1 and 10 and proportions 0.2 and 0.5), T5 (b), T6 (ii).

Inferred: R1's stated base values hold, so the ticket goes on to the fix commit after KIT-190 merges (front matter,
RC-03).

## Fix and R2

Commits: 2622fcff (true merge with `v0.9-staging` 9f3b7c49; no 7.3 file changed except the marker list, so no R1b);
6e02de15 (the fix); 5d322d9c (S1 instrument repair, docs only).

**R2, run 45469464049 on 6e02de15**, base e305fab3 (staging head when sent; it differs from the merged 9f3b7c49 only in
`tests/testthat/test-omic-enrich-module.R`, KIT-205's re-freeze, and the ledger; neither in R2's scope, `R/` identical).
363.1 job-minutes (3 + 3 shards, check on both sides, scripts), ceiling 675. Observed (`r2/compare.tsv`):
- 39 files `same`; `test-kit202-intensity-filter-samples.R` `only-on-head`, 163 passed, 0 failed; `test-metab-qc-shared.R`
  `count-changed` (11 -> 12 passed in the edited test), no failure: as 6.2 predicts.
- `test-workflow-omics-declarations.R` `message-changed` (base 46 / 22 failed, head 45 / 23): the head adds one failing
  expectation, `expect_identical(current$registry, frozen$registry)`, whose differences are the intensity filter
  declarations' `fixed` (a list on head, NULL frozen) and `global_statistics` text, both omics (`r2/moved-declarations-registry.txt`).
  Decided by the owner (VAL-181): no exception by file; `fixtures/kit073-omics-declarations-frozen.rds` is deferred to
  the end-of-project re-record; it does not stop the ticket.
- `test-workflow-omics-binding-steps.R` `message-changed` (58 / 1 failed on both sides), **not predicted**: the one
  failing expectation (KIT-094's characterisation against `fixtures/kit094-characterisation.rds`) fails on base for input
  sha256 values only; on head it also shows `steps[[3]]$parameters` of length 3 against 2, both omics: the intensity
  filter step's record carries the marker entry (`r2/moved-binding-steps.txt`). 6.2 predicted no record snapshot moves;
  this one does. By 6.2's rule the ticket stopped before merge. Base message: the four `inputs[[i]]$sha256` differences
  (both omics), then `records`. Head message: the same four, plus `actual$<omic>$steps[[3]]$parameters` length 3 against
  2 and `parameters[[3]]` a list against absent, both omics. Predicted by the ticket: no. Decided: VAL-189 (the owner,
  2026-11-25): the fixture is not edited and is deferred to the end-of-project re-record, as VAL-181; it does not stop
  the ticket.
- Check: findings identical on both sides (0 errors, 3 warnings, 3 notes each); none only on head (`r2/check-compare.txt`).
- S1 died in R2 (instrument defect: the child's loop overwrote its output-directory variable; no rows written, no
  result). Repaired in 5d322d9c (VAL-54 instrument run 1 of 3; no expected value, comparison or stop rule changed).

**S1, run 45482129676 on 5d322d9c** (base argument e305fab3), 16.3 job-minutes, ceiling 75. Every criterion equal
(`s1/compare.tsv`, `s1/filter-identity.tsv`): (i) unaffected base = head 836/836 (11 assays x 76 settings); (ii) head =
design-sample oracle 228/228 and base = numeric-column oracle 228/228 (3 affected objects x 76); (iii) whole =
partitioned 42/42 on head and on base, unaffected base = head 33/33 per variant; (iv) head under the earlier option =
base 228/228.

Stored references deferred to the end-of-project re-record: `tests/testthat/fixtures/kit073-omics-declarations-frozen.rds`
(`test-workflow-omics-declarations.R`, VAL-181) and `tests/testthat/fixtures/kit094-characterisation.rds`
(`test-workflow-omics-binding-steps.R`, VAL-189).

## Acceptance criteria

- A1: met. Base side R1; head side R2: every new test passes (163), T8's file passes.
- A2: met (S1, run 45482129676).
- A3: met with two declared, decided differences: `test-workflow-omics-declarations.R` (predicted by the ticket: no;
  decided: VAL-181) and `test-workflow-omics-binding-steps.R` (predicted by the ticket: no; decided: VAL-189), both
  `message-changed` by the marker and declaration text only. Every other verdict as 6.2. No `compare-error`, no
  `new-failure`.
- A4: met (check findings identical).
- A5: pending (R3 on the merge with staging 393ebb85, d593e351; R1b not needed: the merge changes no 7.3 file and no
  lipid or metabolite filter, import or design-save file).
- A6: met (`git diff 9f3b7c49..6e02de15 -- R/` touches only 7.3's files; the comparison lines unchanged).
- A7: met (NEWS, ledger row `VAL-<merge>`; `DESCRIPTION`, `NAMESPACE`, lock files unchanged; `man/` as 7.3).

## Not checked

- R3 (full suite) and the merge with the current staging.
- Whether other stored record snapshots outside R2's scope carry the intensity filter step's parameters (R3 would show).
- Inferred, not observed: T5 (a)'s partitioned refusal passed on head, so the error class reaches the caller.
