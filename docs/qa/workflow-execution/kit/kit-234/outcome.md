# KIT-234 outcome: the R Markdown workbooks no longer ship

Ticket `tickets/kit/KIT-234.md` r1.2 (Re-check 1's R1-01 to R1-06 applied, with the coordinator's rulings). Owner's
decision: methods ledger VAL-209. Branch `kit-234` from `v0.9-staging` 67569959.

## Commits

| Commit | What |
|---|---|
| 327e4dfa | tests first: `tests/testthat/test-kit234-workbooks-withdrawn.R` (W1 to W5, W7 to W9) and `w6.R`; nothing under `R/` |
| 652c962b | the change: `git mv inst/workbooks archive/workbooks`, `archive/workbooks/README.md`, `^archive$` in `.Rbuildignore`, `archive/` in `.renvignore`; `R/app_server.R`'s two copy blocks replaced by the folder creation (E-5); the six workbook-reading tests removed or ported; `?integrationTrainModel` Limitations, `?integrationFindInputs`, `?integrationSettings`; README section 2, NEWS, testers' guide; ledger row `VAL-<merge>` |
| a8f0d846 | coordinator's ruling after AUDIT-KIT-239-plan-20261125: the limitation text in `?integrationTrainModel` replaced by that audit's evidence-based draft (section 6), not the workbook's sentences; W9 accepts "synthetic" or "simulated" |
| 450cdece | ledger row reworded to state that draft |
| (this commit) | records |

## Runs (observed)

| Run | Commit | Observed |
|---|---|---|
| R1 45522636588 | 327e4dfa | test file: 11 passed, 18 failed. Failing: W1 (`system.file` = `<pkg>/inst/workbooks`), W2 (five DIA and two metabolomics `.rmd`), W3 (the edited copy's md5 changed; `mine.R` unchanged), W4 (`integration_workflow_starter.rmd`), W7 (workbooks in both recreated folders), W8 (a: "menu() cannot be used non-interactively"; b: the same, no "Using the only ..." message, no pushed file), W9 (all 7 checks absent). W5 passed. W6 exit 1: the installed library's `workbooks/` held the eight files; the tarball held 12 matching paths (the eight files and four folders). Stop rule 1 not met |
| R2 45524570691 | 652c962b against 327e4dfa | `r2/compare.tsv`: `test-enrichment-significance.R` 191 to 181 passed, S11 gone; `test-go-propagation.R` 331 to 321, T10 renamed (the same number of tests; compare lists the old and new titles); `test-integration-removed-factors.R` 526 to 520, same tests; `test-integration-workbook.R` 53 to 43 (KIT-143 static test 7 and KIT-154 39 gone, the new KIT-154 (KIT-234) loop 36); `test-integration-workbook-mofa.R` 3 skipped on both; `test-integration-workbook-api.R` identical; new file `fixed` (29 passed on head); the eight other scope files `same`. `test-workflow-saved-da.R` fails 171 expectations on both sides with identical messages (verdict `same`; the runner-versus-laptop stored references of the KIT-156 record). `r2/check-compare.txt`: 3 WARNINGs and 3 NOTEs on both sides, 0 only on head, 0 only on base. W6 on the fix: exit 0, `system.file("workbooks")` = "", tarball 7,696 paths and none matching |
| L 45524568014 | 652c962b | lane `run-mofa-smoke-k234` (copy of the k184 lane; pandoc 2.9.2.1): `test-integration-workbook-mofa.R` 3 tests, none skipped, 0 failed: the defaults test (8 passed), the equivalence (29 passed; trained model identical to the reference side, the seed-43 control differs) and the ported direct run (12 passed, 153 s; its repeat from the record "equal TRUE, identical TRUE, largest relative 0"; both STRING enrichments repeated with equal ranking and terms). Exit 0 |

Job minutes 152.8 = 10.19 core-hours (4-core runners), against 22.67 expected and 100.0 ceiling. No correction pass
and no instrument rerun used.

## Acceptance criteria

| A | Criterion | Evidence | State |
|---|---|---|---|
| A1 | no workbook in the tarball | W6 (ii) in R2: 0 matching paths of 7,696 | met |
| A2 | none in an installed library | W6 (i) in R2: `""` after `R CMD INSTALL` | met |
| A3 | the application copies none | W2, W4 pass in R2 (fail in R1) | met |
| A4 | existing copies untouched | W3 passes in R2 (md5 identical) | met |
| A5 | `scripts/<omic>` still made where and for whom it is today | W2, W4, W5, W7 pass in R2 | met |
| A6 | kept tests run; per-file counts as stop rule 2 | R2 `compare.tsv` | met |
| A7 | R CMD check identical base against head | R2 `check-compare.txt` | met (on 652c962b; see Not checked) |
| A8 | no stored reference moved | no file under `tests/testdata/` changed (`git diff --stat 67569959 HEAD -- tests/testdata` empty); no stored-reference failure differs | met |
| A9 | the twelve functions exported, documented, present | the KIT-154 (KIT-234) test, 36 passed in R2 | met |
| A10 | the twelve functions end to end with real MOFA2 and pandoc | L: the ported test passed, not skipped | met |
| A11 | "Push to GitHub" as section 2.5 | W8 passes in R2 | met |
| A12 | INT-28's limitations in `?integrationTrainModel`, pointed to from `?integrationFindInputs` | W9 passed in R2 on 652c962b (workbook sentences); the replacing text (a8f0d846) not yet run | met on 652c962b; a8f0d846 to ride in the next run |
| A13 | texts as the owner words them | `wording.md`: every text is a draft | open: owner's wording (OD-1, OD-2) |

## Inferred, not observed

- The end state of a new project equals today's apart from the copied files (the create-path folder creation is
  undone and redone by `setupDirectories(force = TRUE)`); W2 observes the end state, not the intermediate.
- For a non-interactive R with several `.Rmd`, the push fails through `utils::menu`: observed in R1 (W8 base), and the
  same holds in the installed application, which was not run.

## Not checked

- Commits a8f0d846 and 450cdece (the replaced limitation text, W9's relaxed term, the ledger wording) have not run:
  W9 and R CMD check's Rd parse of the new `\section{Limitations}` text are by reading only.
- The installed bundle itself; the application was not started.
- The full suite (shared with the next checkpoint, E-4).
- `project_setup.R` after this change reaches `main` (404 for each workbook; see `wording.md`), not run.

## For the change report

**What changed and who was affected.** Until this change every installation carried eight R Markdown workbooks
(five DIA proteomics, two metabolomics, one integration), and the application copied them into
`<project>/scripts/<omic>/` at every session start, overwriting a user's edited copies in place (KIT-231 G0,
`outcome.md:51` on `origin/kit-231`). Users of the scripted workbook route were affected: by KIT-231 G0 (run
45511537888, observed chunk by chunk), each of the seven omics workbooks stopped with an error at one or more steps
against this version. Now no workbook is installed or copied; copies already in a project stay and are no longer
replaced. Integration workbook users keep the exported `integration*` functions, whose limitations are now on
`?integrationTrainModel`. "Push to GitHub" (Summary tab) now sees only the project's report templates: with none it
says so, with one it uses it (and still needs `<project id>.Rproj`, F7), with several it stops in a non-interactive R
as before.

**Code references.** Base 67569959: `R/app_server.R:327-346` (create path copy) and `:523-549` (watchdog copy);
`inst/workbooks/` (eight files). Fix 652c962b: `R/app_server.R:327-331` and `:508-516` (the folder creation for
proteomics, metabolomics and integration); `archive/workbooks/`; `.Rbuildignore:30`; tests removed or ported:
`test-enrichment-significance.R` (S11), `test-go-propagation.R` (T10's loop), `test-integration-workbook.R` (KIT-143
static test, KIT-154's workbook-code checks), `test-integration-workbook-mofa.R` (render ported to a direct run),
`test-integration-removed-factors.R` (KIT-184's workbook lines).

**Results.** No computed value, record or stored reference changes. Removed with the render test: the only real-pandoc
HTML render of the integration workbook. Kept: the twelve functions end to end with real MOFA2 and the docx record
(L), the module's chain on the lane (`test-integration-record-mofa.R`, `test-integration-enrichment-mofa.R`), and the
wrappers on the corpus with a stand-in trainer (`test-integration-workbook-api.R`).

**Evidence.** R1 45522636588, R2 45524570691, L 45524568014 (tables in `r1/`, `r2/`, `l/`); synthetic corpus only.

## R3, declared before it is sent (2026-11-25)

- **Purpose:** after the owner's answers (r1.3) and the merge of staging dd769671, show on a runner what R2 and L did
  not see: a8f0d846's and 1a40623c's texts (W9 with the final limitation text, R CMD check's parse of the new Rd
  section) and the archiving of `project_setup.R`.
- **Binding:** head = `kit-234` at the declaration's commit (the merge b54e87ca plus this declaration); base = staging
  dd769671. Files: the 14 of R2's scope; 1 shard per side; `check: true`; W6 in the scripts job.
- **Ceiling:** 345 job-minutes (as R2: 2 x 60 tests, 2 x 90 check, 45 scripts). One run; no re-send for text.
- **Expected (base staging against head):** `test-enrichment-significance.R` passed -10, one test fewer;
  `test-go-propagation.R` passed -10 (T10 renamed); `test-integration-removed-factors.R` passed -6;
  `test-integration-workbook.R` 53 -> 43; `test-integration-workbook-mofa.R` 3 skipped on both;
  `test-integration-workbook-api.R` identical; `test-kit234-workbooks-withdrawn.R` only on head, all passing (W9
  included); every other file identical (`test-workflow-saved-da.R` failing identically on both); R CMD check
  findings identical (0 only on head); W6 exit 0.
- **Decides:** A7 and A12 on the final text, and stop rules 2 to 4 on the merged head. Anything else, or a check item
  or test changed by moving `project_setup.R`: stop and report.

## R3, result (run 45573868468)

Observed (`r3/`): every difference is the declared one. `test-enrichment-significance.R` 191 to 181,
`test-go-propagation.R` 331 to 321 (T10 renamed), `test-integration-removed-factors.R` 526 to 520,
`test-integration-workbook.R` 53 to 43, `test-integration-workbook-mofa.R` 3 skipped on both,
`test-kit234-workbooks-withdrawn.R` only on head with 29 passed (W9 on the evidence-based text with the owner's
additions), every other file `same` (`test-workflow-saved-da.R` 171 identical failures on both). R CMD check: 3
WARNINGs and 3 NOTEs on both sides, 0 only on head: moving `project_setup.R` changed no check item. W6 exit 0.
A7 and A12 are met on 51d53171. 111.6 job-minutes (7.44 core-hours); the ticket's total is now 264.4 job-minutes,
17.63 core-hours, against 22.67 expected and 100.0 ceiling.

Not checked by a run: d1fade06 (sent after R3): the two corrections from KIT-239 Re-check 1 (R1-01, R1-03) to the
limitation text in `R/utils_integration_workbook.R` and `man/integrationTrainModel.Rd`. By reading only; W9's terms
are all still present in it, and its Rd parse rides in the merge train's check.
