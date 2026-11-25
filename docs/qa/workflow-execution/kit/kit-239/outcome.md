# KIT-239 outcome: the limitations of the MOFA2 integration at the Results tab, in the report and pointed to from the Model tab

Ticket `tickets/kit/KIT-239.md` r1.2 (step 0 applied Re-check 1's R1-01 to R1-05 and recorded the owner's OD-2 answer).
Branch `kit-239` from `v0.9-staging` dd769671 (KIT-234 not merged in: no file of this ticket is changed by KIT-234; merge
order KIT-234, then this). No R was run on the laptop; every computation ran on GitHub-hosted runners.

## Commits

| Commit | What |
|---|---|
| fce22178 | Tests only (nothing under `R/`): L1 to L3 (`test-integration-results.R`), L4, L5 (`test-integration-record.R`), L6 (`test-integration-gui.R`), the guard L8 (`test-integration-record-mofa.R`) |
| 09286c2d | The change, L7, NEWS, the ledger row `VAL-<merge>`, `wording.md`, R1's tables |
| 0813acb7 | Records: R2's and R3's tables, `runs.tsv`, this file |
| 9bb30a1c | Merge of `v0.9-staging` 16112917 (conflict only in the ledger: staging's rows kept, this ticket's after VAL-229) |
| c764d390 | Tests only, for the implementation audit (I01 to I06; nothing under `R/`) |
| c3888b67 | Wording commit: I01, I02 (draft pending the owner), I03, the I06 guard, the ledger row and INT-28's column (I07), `wording.md` |
| (this commit) | Records: R4's tables, `runs.tsv`, this file |

## Runs (`runs.tsv`)

| Run | Commit | Observed |
|---|---|---|
| 45575068720 | fce22178 | Request file pushed blank by mistake; stopped at `plan` (0.1 min), nothing ran. Not a result; the request was sent again under a new name |
| R1 45575122659 | fce22178, base none | `test-integration-record.R` 134 passed, 5 failed; `-results.R` 170 / 6; `-gui.R` 399 / 1. Every failure is a KIT-239 test, at its stated base value (`r1/failures.txt`): L1 testid absent; L2 keys absent and output `""`; L3 `is.function` FALSE (both helpers); L4 `section_limitations` and the GUI keys absent, heading position -1 (Methods heading found at 621785, so the page rendered: not skipped); L5 `is.function(reader)` FALSE; L6 `model_limits_pointer` absent |
| R2 45577927341 | head 09286c2d, base fce22178 | `r2/compare.tsv`: results 170+6 failed on base, 241 passed 0 failed on head (`fixed`); record 134+5 / 155+0 (`fixed`); gui 399+1 / 401+0 (`fixed`); the other 16 scope files `same` (identical counts and status, the six MOFA2 files skipping on both). Check: 0 errors, 3 warnings, 3 notes on both; findings identical (6), none only on head (`r2/check-compare.txt`). Head logs: "KIT-239 L2: 20 samples, 19 in every data set, 1 missing from a data set, 5 data sets averaged"; "KIT-239 L4: Limitations at 626759, after Methods (draft) at 621785, before Data sets combined at 630384; 20 samples, 19 shared" |
| R3 45577930225 | 09286c2d, MOFA lane `run-mofa-smoke-k239` 81882eae | `test-integration-record-mofa.R` 29 passed, 0 failed, not skipped (`r3/tests.tsv`); "KIT-239 L8: likelihoods gaussian, gaussian, gaussian, gaussian, gaussian, gaussian; groups 1" (`r3/console.log`). The Word report rendered from the real MOFA2 model holds the Limitations section after the Methods draft (`r3/report-limitations.txt`) |

| R4 45588000576 | head c3888b67, base c764d390 | `r4/compare.tsv`: on base exactly the eight new assertions fail at their stated values (`r4/failures-base.txt`): L3 count 1 ("samples are missing") and one database ("these external databases"), L7 `limits_averaged_none`'s term, `help_scaling`'s new terms and its two old phrases, L5 unreadable enrichment (`character()`, not NULL). On head 0 failures: results 258 passed, record 167; gui 401, `test-integration-model.R` 263 and `test-integration-removed-factors.R` 526 passed (3 skipped) identical on both sides |

Job minutes 266.1 (R1 17.8 + 0.1, R2 151.1, R3 14.0, R4 83.1) = 17.7 core-hours; expected 21.2, ceiling 96.0. R4's two jobs
went 23.1 minutes above its 60-minute line, taken from the unused instrument line as declared before the request (total
ceiling unchanged). No correction pass, no instrument rerun.

## Implementation audit (AUDIT-KIT-239-impl-20261125) and the wording commit

- I01: `limits_averaged_none` now "No biological sample has more than one run in any data set, so no runs were averaged."
  (the earlier text was false beside the missing-samples sentence when samples are missing and nothing was averaged).
- I02: `help_scaling`'s last sentence is KIT-246 option (a): "Without missing values the difference could not be told from
  none (29 rather than 26 of 30); with 60 samples the pattern was recovered closely in 30 of 30 either way." **Draft; the
  owner's answer on parked row 10 of 2026-11-25 is pending and may change it.**
- I03: agreement-free texts for `limits_missing_samples` and `limits_services_report`.
- I04: L7 holds the must-survive terms. I05: L4 and L5 assert the missing-samples and runs-averaged counts on the report
  path. I06: entries recorded but none readable give the general services sentence (KIT-181's GO-only case is owed to
  KIT-181). I07: ledger row and INT-28's implementation column. I10: ticket corrected (transcriptomics two runs per sample).
- I01 and I03 were corrected after the owner's acceptance and are to be told to the owner; final polish stays the owner's.
- After the merge of staging 16112917 R2 and R3 were not rerun on the merged tree: the merge touched no integration file
  (only the ledger conflicted). R4 ran on the merged tree.

## Acceptance

| A | Criterion | Evidence | Verdict |
|---|---|---|---|
| A1 | Results tab block after "Notes from the model" with the six points and the model's numbers | L1 (position after the advisories, before Variance explained), L2 (heading, lead, every fixed entry, samples 20/19, missing 1, five averaged data sets each with its counts, bias sentence) pass on head (R2) | met |
| A2 | Report "Limitations" directly after "Methods (draft)", computed from the record; services or none; incomplete record gives the fixed sentences and "could not be read", never "None recorded." | L4 (heading by tag between Methods and Data sets combined; fixed entries; services none; record counts), L5 (record facts from JSON; planted KEGG entry gives "KEGG Release 115.0"; incomplete block NULL; NULL text holds `limits_not_computed`, not `none_recorded`) pass (R2); the real-model report in R3 | met |
| A3 | Bias sentence only when runs were averaged; missing-samples sentence only when a sample is missing (both directions) | L2 (present direction on the corpus), L3 (a) to (g) (absent directions) pass (R2) | met |
| A4 | Record, xlsx, `study_parameters.txt`, Methods paragraph, claim text unchanged | Diff (stop rule 3): no line of `integrationClaimLimitText`, `integrationRecordAnalysis`, `integrationRecordTables`, `integrationMethodsText`, `integrationStudyParameterLines` or a record text they use changed (one record text added, `section_limitations`, read only by the new chunk); the existing record tests' expectations unchanged and passing (R2: record base 134 existing passed, head the same plus the new) | met |
| A5 | No stored reference moves | nothing under `tests/testdata/` changed (stop rule 4) | met |
| A6 | Model tab pointer after training; scaling help in numbers, no "weak" | L6 and L7 pass (R2) | met |
| A7 | "Gaussian" and "one group" hold for a real MOFA2 model | L8 passes on the MOFA lane (R3) | met |
| A8 | The texts are the owner's | `wording.md`: the texts are the owner's decision of 2026-11-25 with the three corrections made after it (R1-01, R1-03, R1-04) and told to the owner; final polish the owner's | met, final polish open |

Stop rules: 1 held (R1); 2 to 5 held (R2); 6 held (the six points present; every number as the evidence, checked
against the ticket's section 7 by script); 7 held (R3).

## Observed and inferred

Observed: the run results above; the report text from a real MOFA2 model (`r3/report-limitations.txt`).
Inferred (by reading, not run): the Record tab's docx report carries the section as the html one does (R3 shows the docx
from the Record tab with it); a record built before this ticket renders the section from its own `matching` block.

## Wording observations for the owner's final polish (not changed: the texts are decided)

- `limits_missing_samples` with one sample read "1 of the 20 biological samples are missing" (R3): fixed by I03.
- `limits_averaged_item` reads "the runs of 20 of its 20 biological samples were averaged" when every sample was averaged.
- (Earlier note on KIT-234's point-4 clause removed: KIT-234 1ba403a9 restored it; the six sentences agree word for word.)
- The lead differs: "from the factors." here, "from a model." on the help page (I08, the owner's polish).

## Not checked

- The application in a browser (KIT-238's IJ' is to assert `integration-results-limitations` and the report heading).
- A record with enrichment from a real service rendered end to end (L5 plants one KEGG entry in the record).
- The full suite (shared with the next checkpoint, as the ticket says).
- The merge with KIT-234 (NEWS and ledger rebase by the coordinator); this branch was tested without it.

## For the change report

What changed: users of the Integration tabs and readers of the integration report were not told the limitations of
the factor model (INT-28) anywhere but the integration workbook, which no longer ships (VAL-209). They are now told, on
the Results tab after "Notes from the model" and in the report after "Methods (draft)", six points: how small a pattern
the package's simulations could find at 20 and 60 samples, beside this model's own number of samples and samples in every
data set; how missing values are treated (missing at random), with this model's count of samples missing from a data
set; per data set the runs averaged, and that such averages can be too high for low-abundance features when values are
missing because they were too low to measure (only when runs were averaged); one group, Gaussian, one input per omic,
estimated values treated as measured; what the check against the original analysis script covered and how the defaults
were chosen; the external enrichment services (in the report, those this analysis used with their releases, or none). The
Model tab now points to them after training, and its help for per-feature scaling no longer says "weak factors": it
states the simulated pattern's size, the criterion and the condition in numbers. Who is affected: every user of the
Integration tabs and every integration report; no computed value, record field, Methods sentence or stored reference
changed.

Code (base dd769671 / fix 09286c2d): `help_scaling` `R/utils_integration_gui.R:261-263` / `:262-269`; new entries
`:412-469` at the fix; `integrationResults` return `R/utils_integration_results.R:339-347` / `:339-350`; helpers
`:824`, `:849` (new); `R/mod_integration_results.R:40`, `output$limitations` (new); `R/mod_integration_model.R:177`
(new); `R/utils_integration_record.R:56` and `:1084` (new); `inst/reports/integration/integration_report.rmd:26-35`
(new chunk).

Evidence: R1 45575122659, R2 45577927341, R3 45577930225, R4 45588000576; the numbers in the sentences are KIT-151 c and d's simulations
(one simulated design, 20 and 60 samples, sources in `wording.md`), whose limits the sentences state.
