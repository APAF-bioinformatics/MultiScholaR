# KIT-267 outcome: the confirmed Enrichment defaults lose the "not yet confirmed" mark; the organism defaults say they are human

Owner decision: ledger VAL-310. Ticket r1.1 and its plan register (Re-check 1: ready, fixes R01 and R02). Branch `kit-267`
from `v0.9-staging` 4d588ba3ed5c6fee76a5a1450c028b409ba419ba. No default, choice, range or computation changed.

## What was run

| Run | Head | Base | Result |
|---|---|---|---|
| R1 45704530515 | 8f334001 (tests only) | 4d588ba3 | T6 (b) 9 failures, T10 5, T3 2: as predicted. Kit240 file otherwise 82 passed, 1 skipped; base 92/0/1 |
| R2 45707020223 | 21af4bc3 (fix) | a6fbb6fc (tests commit plus docs; R code identical to 8f334001) | 22 files: 20 `same`, 2 `fixed`; check identical |

Commits: 8f334001 (tests and R1 predictions), a6fbb6fc (R1 record, R2 predictions; also carries the testers' guide edit,
which is docs only), 21af4bc3 (the change), and the records commit.

## Observed

- R1: failure values as the ticket and predictions state (T6 (b): setequal 1 + eight per-setting; T10: old-mark count 8 vs
  0, organism-mark count 0 vs 3, three label matches; T3: 8 provisional of 13 vs 0, `mark_organism` absent).
- R2 scope recounted on the branch point: 22 files (the ticket's 21 plus the new file). Every file `same` except
  `test-kit240-integration-default.R` (base 82 passed 14 failed; head 96 passed 0 failed, 1 skipped, T9) and
  `test-kit267-enrichment-defaults.R` (base 1/2; head 3/0). KIT-240 T6 (c) passes on both sides.
- R CMD check: 0 ERROR, 3 WARNING, 3 NOTE on both sides; six findings identical; none only on head or only on base.

## Acceptance

| # | Criterion | Evidence | Met |
|---|---|---|---|
| A1 | no "awaiting" in `settled_by`; the eight name VAL-310 and "decided 2026-11-26"; STRING pair unchanged | T1 (T6 (b)) fails 9 on R1, 0 on R2; the STRING pair strings were not edited | yes |
| A2 | `integrationParameterProvisional` FALSE for all 13 entries | T3, R2 | yes |
| A3 | no old mark; exactly three organism marks, each after its label | T2 (T10), R2 | yes |
| A4 | every default and choice equals base | KIT-240 T6 (c) green on base and head | yes |
| A5 | Model, Results and record unchanged | the 20 other files `same` | yes |
| A6 | R CMD check identical | `r2-check-compare.txt` | yes |
| A7 | NEWS replaced, guide sentence removed, ledger row and outcome written | `NEWS.md:17`; guide section 6; ledger row `VAL-<merge>` after VAL-313; this file | yes |

## Deviations and notes

- T1 is written as `expect_true(grepl("decided 2026-11-26", .) && grepl("VAL-310", .))`, one expectation per setting,
  not R01's `expect_match(..., "decided 2026-11-26, VAL-310", fixed = TRUE)`. Re-check 1 arrived after R1 was out; the
  coordinator ruled to leave T1 as sent (same meaning, base count 9). All three new strings do contain the exact
  substring "decided 2026-11-26, VAL-310".
- R1 cost 45.2 job-minutes against the ticket's 25 job-minute step figure (four test jobs, 9 to 14 minutes each, plus
  plan and results; the ticket's 14 is the run's elapsed time). Total ceiling (470) unaffected.
- The base of R2 is a6fbb6fc, not 8f334001: it differs from the tests commit by docs only (R1 record, R2 predictions,
  the guide sentence).
- R02 (Re-check 1): the ticket's dependency note should also name `R/utils_integration_gui.R` (other hunks) among files
  shared with in-flight branches; nothing to do in code.
- `kit-240/default-identity.R` was neither requested nor edited (ruling P01).

## Core-hours

R1 45.2 job-min = 3.01 core-hours; R2 182.6 job-min = 12.17 core-hours; total 227.7 job-min = 15.18 core-hours against
11.27 expected and 31.33 ceiling. No correction pass used. No instrument defect.

## Not checked

- How the label looks in the rendered tab (two adjacent parentheses): the tests match strings, nothing was rendered.
- That no other consumer displays `integrationSettings()$enrichment$settled_by` (grep only).
- Requests on `run-tests` naming `kit-240/default-identity.R` beyond those the ticket lists.
- Whether the testers' guide is bundled into an installer.
- The full test suite (the merge check rides the next train).

## For the change report

Not a change of results. What users are told changed: before, the Enrichment tab said eight defaults were "not yet
confirmed" (a mark that existed only on `v0.9-staging` since KIT-240, in no tag and not on `main`); now it says none,
and three labels (STRING species, KEGG organism, Reactome organism) end "(human; change for other species)". Code
references at base 4d588ba3: `R/utils_integration_enrichment.R:1114-1126` (the strings),
`R/mod_integration_enrichment.R:33-35` (the label); at 21af4bc3 the same functions plus
`INTEGRATION_ENRICHMENT_ORGANISM_SETTINGS` and wording key `mark_organism` in `R/utils_integration_gui.R`. No result,
file, hash or default changed (KIT-240 T6 (c) green). Evidence: runs 45704530515 and 45707020223. KIT-240's
release-checklist item "the eight enrichment defaults confirmed or disclosed (R1-04)" closes here.
