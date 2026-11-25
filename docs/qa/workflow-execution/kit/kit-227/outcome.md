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
  (otherwise a Skip read from the writer's section would print "not recorded"), and the record read also takes
  `ruv_grouping_variable` (the overlay coerces the record's value to NA through the parser's logical NA, so the base
  prints "'NA' variable" when the RUV section is in the writer's shape). An extra E5 assertion covers the writer's RUV
  section with a default record; its file keeps the `-` Assay Information lines every lipid/metab case carries,
  because the method text needs a readable 'Number of Assays'. On the file the writer actually produces (`*` on every
  line) that count is not read, and by reading `FormatMethodDetails` stops at `if (num_assays > 1)`, so the report
  prints 'Error running Main function: …' and no Methods, on base and head alike (inferred, not run; audit K227-I01,
  a separate ticket). The writer's per-assay lines are not read either: on the writer's RUV section the paragraph ends
  after the grouping sentence, with no per-assay k or control counts (base and head alike). The parser itself is not
  changed.
- The cases R03 requires that the drafts do not word (a percentage, count or k not recorded outside Manual) carry
  implementer wording, listed in `wording.md` for the owner's pass.

## For the change report

- **What was wrong, who was affected.** The report Methods paragraph on RUV described every applied run as an
  automatic optimisation (proteomics LFQ, TMT, DIA-NN and DIA with limpa without the limpa quantification;
  lipid/metab older form; transcriptomics, where a manual run was described as a grid search with the grid's best
  scores); a skipped RUV, and a run with nothing recorded, were said to be "not applied due to dataset constraints";
  a manual run whose plot failed printed "NA%"; the RUV count sentence could follow a Skip; the limpa Results caption
  claimed covariates in the correction mode and after a Skip; lipid/metab printed the grouping as 'NA' when the RUV
  section was in the writer's shape and a record was present (whether GUI users saw it depends on audit K227-I01: by
  reading, the writer's full file stops the lipid/metab Methods before the RUV paragraph) and called a non-`group` column "preserving biological signal". Affected: anyone pasting those
  Methods from Manual, Skip, unrecorded or non-`group` runs.
- **Code.** At cf8af245: `inst/reports/proteomics/LFQ_report.rmd:268-273, 470-503, 600-638, 938-939, 1049` (TMT and
  DIA-NN twins), `DIANN_limpa_report.rmd:270, 399-411, 528-566, 949-950, 1012-1017`, `lipidomics_report.rmd:435-445,
  607-664, 716-790, 1014-1015, 1148` (metabolomics twins), `transcriptomics_report.rmd:63-80`. At f9324d39:
  `LFQ_report.rmd:270-283, 472, 489-552, 653-683, 1006-1018, 1128` (TMT `:270-283, 472, 489-552, 644-674, 997-1009,
  1119`; DIA-NN `:278-291, 491, 527-590, 700-730, 1052-1064, 1174`), `DIANN_limpa_report.rmd:271-280, 410-474,
  594-629, 1036-1048, 1112-1123`, `lipidomics_report.rmd:435-445, 608-721, 774-827, 1088-1100, 1233` (metabolomics
  `:423-433, 597-708, 761-813, 1062-1074, 1178`), `transcriptomics_report.rmd:66-79`.
- **What changed in results.** No number changed; only the sentences. Automatic runs on the `group` column print the
  earlier text byte for byte (G1). Projects that already generated a report keep their template copy (OD-3).
- **Not changed (stated limits).** `study_parameters.txt` and the normalisation product still give the Skip reason
  'due to dataset constraints' (`R/func_general_filemgmt_study_parameter_support.R:238`,
  `R/utils_workflow_dia_norm.R:111`; ticket E-6). Projects that already generated a report keep their template copy
  (OD-3). The lipid/metab sentences reach GUI users only once the lipid/metab parser reads the writer's `*` lines
  (audit K227-I01).
- **Evidence.** Runs 45506521161 (failing first) and 45508506125 (passing, regressions `same`, check identical).

## Open for the owner (audit register `AUDIT-KIT-227-impl-20261125.md`)

- **K227-I01** (major, outside scope; separate ticket, before the release commit). The lipid/metab parser does not
  read the writer's `*` keys. By reading, the writer's real file stops the Methods at `num_assays`. The per-assay k
  and control counts are lost, and several other fields are wrong. KIT-227's lipid/metab text cannot reach GUI users
  until this is fixed.
- **K227-I04** (owner wording). With a design column other than `group`, the unchanged opening "To remove batch
  effects…" contradicts W5 when the column is a batch column. Proteomics W1 is silent on that column's role as
  RUV-III's replicate sets (parity with W5).
- **K227-I05** (owner wording). The NEWS entry names only "the paragraph". It does not name the Figure 1 caption
  changes or the count sentence after a Skip.
- **K227-I07** (outside scope; a new ticket). The lipid/metab Results captions, applied and Skip ("as today" by W3),
  print fixed findings whatever the data and method. This is the class KIT-055 fixed in proteomics.

## Not checked

- No template was rendered end to end with knitr; the chunks were evaluated as extracted (the real parser, overlay,
  record read, method and Results calls).
- The full suite (shared checkpoint, E-2; not sent by this ticket).
- Real writer output for every route; the fixtures follow the writer's code (`support.R:230-317`, `writer.R:586-631`).
- Which figure columns exist after a limpa Skip (the caption is now stage-free).
- G1 for lipid/metab uses `-` bullets (on the writer's shape the base text was false, so byte identity is not the test
  there).
- The limpa route's covariates caption is byte-identical by reading only.
- K227-I01 is inferred, not run: the lipid/metab parser does not read the writer's `*` keys, so by reading the
  writer's real file stops the Methods at `num_assays`. It needs one runner test of a writer-produced file through the
  lipid and metab templates (separate ticket).
