# KIT-263 outcome: proteomics study parameters stop listing other omics' settings

Branch `kit-263`, from `origin/v0.9-staging` 127011429f5160e6a1ed57fd0552dc16a90d5aa9. Ticket r1.1.
Commits: 8c83cf05 (tests only: T1 to T5, the declared harness edit of `test-omic-summary-module.R`, predictions);
19010df7 (diagnostic D1); 0abd5669 (diagnostic D2); 7be5aa5a (tests only, T1 to T3 fixture rebuilt, expectations
unchanged); 447bdf29 (the fix, the Rd, NEWS). Records: `runs.tsv`, `predictions.md`, `r2-*.tsv/txt`,
`r1-2-failures.txt`. Job-minutes: 313.6 (20.9 core-hours; expected 24.77, ceiling 95.0).

## What was run (observed)

| Run | Commit | Purpose | Result |
|---|---|---|---|
| 45676487966 (R1) | 8c83cf05 | tests only, new file | T1 failed at 3 foreign lines (predicted 6), T2 and T3 failed (predicted pass), T4 errored, T5 failed by value. Stopped and reported (unpredicted R1 verdict; the brief's rule). |
| 45679194881 (D1) | 19010df7 | declared diagnostic, printing only | the T1 fixture's file has the foreign block once (S4 block) and an EMPTY "Additional Configuration Parameters" block: `formatConfigList` errors on the shipped config (below). The module-ci proteomics bus gives the writer no S4 object. |
| 45681727116 (D2) | 0abd5669 | declared diagnostic, the application's owner path | proteomics owner save (236 lines): the foreign block TWICE, 3 + 3 = 6 lines; lipidomics owner save (157 lines): no foreign line. Classification (ii), committed beforehand in `predictions.md`. |
| 45684163296 (R1-2) | 7be5aa5a | tests only, fixture rebuilt on D2's owner construction, expectations unchanged (instrument rule, 1 of 3) | T1 failed at 6 (expected 0), T4 errored (constant absent), T5 failed (`NULL` against "transcriptomics"); T2 and T3 passed (8 passed, 2 failed, 1 error). As predicted. |
| 45686553122 (R2) | 447bdf29 vs 7be5aa5a | 43 files, 3 shards per side, R CMD check both | 42 files `same`; the new file 11 passed, 0 failed on head (base 8/2/1; the table's verdict reads `fixed`, not `only-on-head`, because the tests-only commit already holds the file); `test-omic-summary-module.R` 1675 passed, 28 failed on both sides (identical to merge train 11); `test-gui-entry-proteomics.R` 168/25 on both sides (same); R CMD check: 0 errors, 3 warnings, 3 notes on both sides, 6 findings identical, none only on head or base. |

## Acceptance criteria

| # | Met? | Evidence |
|---|---|---|
| A1 | Met, after one instrument re-send | R1-2: T1 (6 against 0) and T5 fail by value, T4 errors as stated, T2 and T3 pass, the base file's six foreign lines printed (the log of 45684163296). The first R1 (45676487966) did not match because the T1 to T3 fixture did not represent the application (D1, D2). |
| A2 | Met | R2 head: all five tests pass (11 expectations). |
| A3 | Met | R2: 42 files `same` including every pre-existing failure list (`test-omic-summary-module.R` 28, `test-omics-summary-core.R` 25, `test-module-ci-trx-summary.R` 1, `test-gui-entry-proteomics.R` 25); check identical. No stored reference moved. |
| A4 | Met | NEWS, ledger row (`VAL-<merge>`, `VAL-<train18-1>` on train 18), this file. |

## Declared deviations

- The first R1 did not give the predicted values (T1 3 not 6; T2, T3 failed). Two declared, printing-only diagnostics
  followed on the coordinator's ruling (D1, D2). D2 classified the fixture as an instrument defect: the T1 to T3 fixture
  called the writer with a plain object and no state manager, so no analysis record existed and the writer's Additional
  block went through `formatConfigList`, which fails on the shipped config (the finding below). The application's save
  uses an owner (analysis record present) and prints that block with the record's formatter. The fixture was rebuilt on
  the application's owner construction (`omicSummaryW1ProtSetup(owner = TRUE)`, test-omic-summary-module.R:826-841,
  without its tab server), expectations unchanged; R1 re-sent once (instrument rule 1 of 3; the ticket's one correction
  pass was not used).
- Declared harness edit: the KIT-108 W1 recorder in `test-omic-summary-module.R:451-458` gains a final `...`, passed to the
  real writer and not recorded; no stored value edited. R2 base (with the edit) and head give the same 1675/28.
- The filter is written with `which()` and `seq_along()` so a list without names is not emptied.
- The two diagnostic test blocks were removed from the committed test file after use (their output is in runs
  45679194881 and 45681727116).

## Observed counts

- Today's tree (staging 12701142, the base of this branch): a proteomics `study_parameters.txt` written by the
  application (owner path, D2) holds 6 lines that mention metabolomics: in "Parameters from Final S4 Object"
  `[metaboliteIntensityFiltering]`, `metabolites_intensity_cutoff_percentile = 1`,
  `metabolites_proportion_of_samples_below_cutoff = 0.5` (lines 101-103 of the 236-line file); in "Additional
  Configuration Parameters" `metaboliteIntensityFiltering:`, `Metabolites Intensity Cutoff Percentile: 1`,
  `Metabolites Proportion of Samples Below Cutoff: 0.5` (lines 203-205). Two settings, printed twice. No transcriptomics
  line (those sections hold no key).
- With KIT-217 merged (predicted, by reading `origin/kit-217` config.ini:110-113, not run): four settings, printed twice,
  10 lines (2 x (1 header + 4 keys)).
- After this fix: 0 foreign lines in either block (R2 head, the owner construction).
- Lipidomics owner save (D2): 0 foreign lines; headers `globalParameters`, `ITSDNormalization`, `daAnalysisParameters`,
  `ruv_number_k`, `ctrl` only.

## The five proteomics reference files of VAL-279

They do not move again under this ticket: R2 shows every stored reference unchanged (42 files `same`, identical failure
lists), and the change is display only. They move because of KIT-217's keys in the stored arguments, not because of the
writer. This corrects what the owner was told in parked row 31 (`parked-20261125.md:379`, "move once more"): the outcome is
better than stated.

## Separate finding, recorded and not fixed here: the writer's "Additional Configuration Parameters" block is empty whenever `formatConfigList` fails

- **What happens.** `formatConfigList(config_list)` stops with `In index: 17 ... purrr::map_chr(): Result must be length 1, not 2`
  on the shipped `config.ini`. Index 17 is `[plotRle]`, whose value `yaxis_limit=-3,3` (`inst/config/config.ini:135` on staging; `:93` is `main`'s line) the reader
  makes a numeric vector of length 2. The branch for a plain (non-list) value converts with `as.character(value)`, which
  keeps two elements; only character vectors are collapsed first. `purrr::map_chr` needs length one, so the call fails.
- **Where.** `R/func_general_filemgmt_config_io_helpers.R:422` (`formatConfigList`), the `map_chr` at `:431`, the conversion at
  `:488`. The writer calls it inside a `tryCatch` (`R/func_general_filemgmt_study_parameter_writer.R:948`, the call at `:963`,
  the handler at `:973-976`; base 12701142 lines). The handler prints `WORKFLOW ARGS: formatConfigList failed: In index: 17.` to the console and
  assigns `config_lines` inside the handler's own function, so the assignment is lost: the file keeps the header
  "Additional Configuration Parameters:" and its dashes and nothing after them. Nothing tells the user.
- **Log line.** Run 45679194881, the log of `test-kit263-prot-study-parameters-sections.R`, "WORKFLOW ARGS: formatConfigList
  failed: In index: 17."; the same error printed directly by the diagnostic ("formatConfigList error: i In index: 17.
  Caused by error in `purrr::map_chr()`: ... Result must be length 1, not 2."). The file then ends at the header (116 lines,
  lines 115-116).
- **Is it on `main`? Yes, by reading** (`git show main:...`, not run): `main:R/func_general_filemgmt.R:1755` and
  `main:R/func_general_helpers.R:1072` both define `formatConfigList` (the later in the Collate order wins) with the same
  `as.character(value)` branch (`:1078` in helpers, `:79` of the former definition); `main:inst/config/config.ini:93` holds
  `yaxis_limit=-3,3`; the reader converts it by `strsplit(items, split = ",")` and `as.numeric` (`main:R/func_general_filemgmt.R:942`,
  `:1117`); the writer's call is at `main:R/func_general_filemgmt.R:5388` inside a `tryCatch` with the same handler (`:5399`).
- **Who is affected.** Every save that reaches `formatConfigList` for the Additional block, that is, the writer's path
  WITHOUT an analysis record (no state owner with reference states: scripted calls, sessions without a reference owner, or when
  the record is unavailable), with a `config_list` that carries a numeric vector of length greater than one (the shipped
  `plotRle` value, so any project that copied the package configuration). With an owner (the application's path) the writer
  uses the record's formatter and the block prints (D2). Other omics whose `config_list` is read from a config with such a
  value are affected the same way; the shipped proteomics configuration is the one that carries it. Not checked: which
  application routes reach the no-record path today.
- **Why it matters for KIT-263.** On that path the foreign block shows once (S4 block) rather than twice, and the filter on
  `config_list` is exercised only on the record path. The fix filters both sources, so both are covered.

## Not checked

- Windows. The no-record path through the application (which routes). Lipidomics, metabolomics and transcriptomics
  files beyond the lipidomics owner case in D2. KIT-217 merged (its counts are predicted from its config.ini, not run).
  Scientific or performance claims: none made; the change is display only.

## For the change report

- **What was wrong, who was affected.** Every proteomics project copies the package's whole `config.ini`, so the
  protein object's arguments and the project's configuration carry the metabolomics intensity-filter settings
  (percentile 1, proportion 0.5), which no proteomics step uses. In the 0.9 builds (tag v0.9.0, and staging until
  this fix) the application's `study_parameters.txt` printed them twice, in "Parameters from Final S4 Object" and in
  "Additional Configuration Parameters" (observed, D2 run 45681727116: 6 lines; the second print is the analysis
  record's formatter, KIT-033). With KIT-217 merged first, four settings twice (10 lines; by reading). In the
  releases before 0.9 that carry the section (`main` 0.4.1.2; tags v0.3.6.3, v0.5.0-rc1; the section was added to
  `config.ini` on 2025-11-19, `main` 9aa462cc8) they printed once, in the first block
  (`main:R/func_general_filemgmt.R:5180`). The second block (`:5388`) was empty, because `formatConfigList` fails
  on the shipped `plotRle` value and the writer swallows the error (the separate finding above; by reading, not run
  on `main`). A reader of a methods record could take them for settings applied to the proteins.
- **Code references.** Base (staging 12701142): `R/func_general_filemgmt_study_parameter_writer.R` printed `@args` at
  `:802-888` and the merged configuration at `:937-980`; the copy at `R/func_prot_import_detection.R:249-266` and
  `R/mod_prot_design_import_helpers.R:181-209`. Legacy on `main`: `main:R/mod_prot_import.R:908-917`,
  `main:R/func_general_filemgmt.R:5180` (printed); `:5379`, `:5388` (the second block, empty on `main`). Fix (447bdf29): the constant `studyParameterForeignSections` above the
  writer (`:1-6`), the argument `omic_type` and the filter after the `config_list` is resolved (`:101-107`, local copies
  only), the Summary save passing `omic_type = omicType` (`R/utils_workflow_summary.R:266-267`), the Rd.
- **What changed in results.** None. Only the printed file: for proteomics the four sections of config.ini's
  Metabolomics and Transcriptomics blocks are left out of both printed blocks. The object's `@args`, the copied
  `config.ini`, the analysis record, reports and replay are unchanged. The foreign settings still show in the project's
  copied `scripts/proteomics/config.ini` (written back by the design builder and copied into the DIA native export), in
  the protein object's `@args` (integration RDS, session state), and in pre-provenance record entries.
- **Evidence.** Runs 45676487966, 45679194881, 45681727116, 45684163296, 45686553122 (this directory). R2: 42 files
  `same`, check identical; D2's file for today's tree shows the 6 lines; R1-2's base value is 6.
