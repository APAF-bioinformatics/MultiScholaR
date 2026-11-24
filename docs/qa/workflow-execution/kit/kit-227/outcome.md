# KIT-227 outcome (r1.3): report RUV sentences that were false for non-default settings

Branch `kit-227` from `origin/v0.9-staging` cf8af245. Commits: a7600526 (tests only, no package code), f9324d39
(the fix: seven templates, NEWS, `wording.md`). No file under `R/` changed. Owner's decisions of 2026-11-25: OD-1
(drafts W1 to W8 as the basis; final wording later: `wording.md`, marked draft), OD-2 (merge before the release
commit), OD-3 (stated limit, NEWS sentence).

## What was run (observed)

| Run | Commit | Observed |
|---|---|---|
| R1 45506521161 | a7600526 | 66 passed, 100 failed. Every E row (E1, E1b, E2, E3, E3b, E4, E5, E6, E7, E8, E9, E10) failed at its stated base value (e.g. "an automatic optimisation process was run" present 8 times where it must be absent; "NA%" 3; "due to dataset constraints" 12; "auto-determined" 3; "'NA' variable" 4; the covariates caption 2; "choosing the maximum composite score" 1). G1 had no failure: the base fixtures equal the shipped templates. Stop rule 1 not met. `r1-45506521161/` |
| R2 45508506125 | f9324d39 against cf8af245 | New file: 166 passed, 0 failed, 0 skipped. The 26 other files: verdict `same` on every one. R CMD check: 0 errors, 3 warnings, 3 notes on both sides, findings identical (none only on head). `test-module-ci-trx-summary.R` ran on both sides: 82 passed, 3 failed (stored-hash failures, identical on base and head), 0 skipped; the skip message "the KIT-040 references are not shipped" absent. `r2-45508506125/` |

## Acceptance (section 4 rows)

| Row | Evidence | Met |
|---|---|---|
| E1, E1b | R1 fails at "an automatic optimisation process was run" and "NA%"; R2 passes the Manual drafts with "5.0%" (file) and "5%" (record) | yes |
| E2, E6 | R1 fails ("due to dataset constraints", "Following RUV-III batch correction"); R2 passes | yes |
| E3, E3b | R1 fails; R2 passes ("is not recorded" in method and caption) | yes |
| E4, E7 | R1 fails; R2 passes (limpa Skip/not recorded; correction and Skip captions; covariates caption kept) | yes |
| E5, E9, E10 | R1 fails; R2 passes | yes |
| E8 | R1 fails; R2 passes | yes |
| G1 | Method and Results text identical to the base chunks (tests/testdata/kit227/base, verbatim at cf8af245) for LFQ, TMT, DIA-NN, DIA with limpa (file only and with a default record), lipid and metab per-assay (both), transcriptomics Automatic; passes on R1 and R2 | yes |
| Stop rule 2 | No stored reference moved; R2 `same` on all 26 files | not triggered |
| Stop rule 3 | Diff limited to the RUV blocks, the parsers' `ruv_applied` elements (plus `ruv_separation_metric` after `ruv_grouping_variable`), the record-read lines after each overlay call, the Figure 1 caption functions and calls, the transcriptomics `ruv-text` chunk | not triggered |

## Departures from the ticket text (inferred from reading; stated for the audit)

- E1 pins "5.0%", not "5%": the writer writes `sprintf("%.1f%%")`, so a file route prints "5.0%" (the register says so).
- The lipid/metab parser keeps the writer's `*` bullet in its keys, so the writer's file is read only through the
  record. Two consequences handled inside owned lines: the parser's `ruv_applied` element also reads `* Status`
  (otherwise a GUI Skip would print "not recorded"), and the record read also takes `ruv_grouping_variable` (the
  overlay coerces it to NA through the parser's logical NA, so the base prints "'NA' variable" on the default route).
  An extra E5 assertion covers the writer's shape with a default record. The parser itself is not changed.
- The cases R03 requires that the drafts do not word (a percentage, count or k not recorded outside Manual) carry
  implementer wording, listed in `wording.md` for the owner's pass.

## For the change report

- **What was wrong, who was affected.** The report Methods paragraph on RUV described every applied run as an
  automatic optimisation (proteomics LFQ, TMT, DIA-NN and DIA with limpa without the limpa quantification;
  lipid/metab older form; transcriptomics, where a manual run was described as a grid search with the grid's best
  scores); a skipped RUV, and a run with nothing recorded, were said to be "not applied due to dataset constraints";
  a manual run whose plot failed printed "NA%"; the RUV count sentence could follow a Skip; the limpa Results caption
  claimed covariates in the correction mode and after a Skip; lipid/metab printed the grouping as 'NA' (on the default
  GUI route too) and called a non-`group` column "preserving biological signal". Affected: anyone pasting those
  Methods from Manual, Skip, unrecorded or non-`group` runs.
- **Code.** At cf8af245: `inst/reports/proteomics/LFQ_report.rmd:268-273, 470-503, 600-638, 938-939, 1049` (TMT and
  DIA-NN twins), `DIANN_limpa_report.rmd:270, 399-411, 528-566, 949-950, 1012-1017`, `lipidomics_report.rmd:435-445,
  607-664, 716-790, 1014-1015, 1148` (metabolomics twins), `transcriptomics_report.rmd:63-80`. At f9324d39: the same
  blocks, rewritten (see `git diff cf8af245 f9324d39 -- inst/reports`).
- **What changed in results.** No number changed; only the sentences. Automatic runs on the `group` column print the
  earlier text byte for byte (G1). Projects that already generated a report keep their template copy (OD-3).
- **Evidence.** Runs 45506521161 (failing first) and 45508506125 (passing, regressions `same`, check identical).

## Not checked

- No template was rendered end to end with knitr; the chunks were evaluated as extracted (the real parser, overlay,
  record read, method and Results calls).
- The full suite (shared checkpoint, E-2; not sent by this ticket).
- Real writer output for every route; the fixtures follow the writer's code (`support.R:230-317`, `writer.R:586-631`).
- Which figure columns exist after a limpa Skip (the caption is now stage-free).
