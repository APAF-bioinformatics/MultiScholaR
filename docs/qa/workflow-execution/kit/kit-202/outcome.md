# KIT-202 outcome (interim: R1 done; the fix commit waits for KIT-190's merge)

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

## Acceptance criteria

- A1: base side met (R1). Head side pending (R2).
- A2 to A7: pending (fix commit, R2, R3).

## Not checked

- Any head behaviour: the fix is not written.
- S1 (`filter-samples-control.R`) has not run; it runs first in R2.
- That the error class reaches the caller through the partitioned step's wrappers on head (T5 (a)).
- R1b (only if the merge with a staging that contains KIT-190 changes a file it names).
