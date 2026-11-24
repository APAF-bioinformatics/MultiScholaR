# KIT-196 outcome: lipidomics and metabolomics DA on features without values

Ticket `tickets/kit/KIT-196.md` revision 1.2 (step 0 applied the re-check's R01 to R06; see the ticket's sections 15 and
17). Owner's decisions: methods ledger VAL-116 (D1 (a), D2 always, D3, D4 KIT-200, D5 budget, D6, D7 KIT-201) and
VAL-117. Branch `kit-196-da-features-without-values`, from `v0.9-staging` 19997323.

## Commits

| Commit | What |
|---|---|
| 8d095b41 | Tests first (T1 to T9, T10), S1 `da-sparse-control.R`, `sim-k188-69fd01bc.R` (SHA-256 c0fbc4ab34925dbbbfee90169dee3e7af60f24b3d29c1ba94aa83175b52a1c22, byte-identical to KIT-188's at 69fd01bc). Nothing under `R/` |
| bf5ed69f, 648a44f0, 52343cca | T7's three instrument repairs (VAL-54; test mechanics only, no expected value changed) |
| 5b73196b | The fix: `R/func_lipid_da_model.R`, `R/func_metab_da_core_engine.R`, `R/func_metab_da_orchestration.R`, four `man/*.Rd` by hand, NEWS, ledger row `VAL-159` |
| d1f310ed | R1 records; S1 option `--base-records` |
| 1b8fa4f2, a1e46859, 6dcf593c | Merges of `v0.9-staging` (ea93d948; 3d4a585d with KIT-194; fe5c9273). Conflicts only in NEWS.md and the ledger |
| fef003c6, 3225baab | R2 and S1 records |

## Runs (`runs.tsv`)

R 4.6.1, KIT-134 snapshot library, `blas: reference`, ubuntu-22.04 runners, limma 3.68.4 (A7). run-tests runs a
request's scripts one after another in one job, so S1's parts L and M were sent as separate requests to run side by side
as the ticket's "2 scripts" jobs.

- **R1 (gate G0 and failing-first), tests-only commits.** 45316123488 and 45316126666 at 8d095b41: S1 base side, parts L
  and M, base 19997323; the two test files. T7's driver failed on both sides without a result (instrument defect), repaired
  three times: 45325563956 (bf5ed69f), 45327425477 (648a44f0), 45329105038 (52343cca). Three instrument runs, the limit.
- **R2.** First sent on the merge with staging ea93d948 (45330948439, 45330952540), then cancelled by me when KIT-194
  merged and the coordinator asked for S1 on its base; no result. Repeated on the merge with staging 3d4a585d:
  45336286628 (38 test files base 3d4a585d against head a1e46859, and R CMD check), S1 base side on 3d4a585d
  (45336278801 L, 45336282667 M, `--base-only`), then S1 head side at fef003c6 against those base records
  (45346222323 L, 45346224846 M, `--base-records=.../s1/r2-base`). Base and head were split into separate jobs because one
  side alone takes 47 to 68 minutes and both together would pass the 120-minute script ceiling.
- **R3.** Full suite base fe5c9273 against head 6dcf593c, 20 shards, with `check: true` (it also covers R2b: staging
  gained proteomics-only code after R2; 6.4's rule does not repeat S1 for it). Results below.

## Observed

### Gate G0 / S1 base side (R1 at base 19997323; repeated identically at 3d4a585d)

`s1/r1/`, `s1/r2-base/`.

- **A2 (i) holds.** Base's `OR` full stops equal KIT-188's rows per cell, scenario, mask and data set: part L 10 = 10 (L12
  VD 20 % x 9; L64 BA N1 402), part M 27 = 27 (M12 VD 20 % x 17, VD 5 % x 2; M168 VD 20 % x 4, BA x 4). `n_T` and
  `miss_share` of every `OR` journey and assay equal KIT-188's (204 of 204 in L, 102 of 102 in M), compared at 7
  significant digits (KIT-188 writes `signif(x, 7)`; the re-check's "8" corrected in the ticket).
- **A2 (iv), base side, holds.** 120 assays had an error (L 66, of which 26 are one-assay failures in L12 and L64 journeys
  that otherwise completed, which KIT-188 recorded `ok`; M 54, none partial). Every one has at least one feature without
  values in its DA input (L 1 to 4, M 1 to 7). The direct re-run of the contrast half gave "prior.weights contain NA values"
  in all 120, and every one had unequal positive residual df (`base-diagnosis-*.tsv`; call recorded, descriptive only). No
  stop without a feature without values.
- `git diff --name-only 98df468b <base> -- R/` (KIT-188's L1 base against S1's bases) names no lipid or metabolomics filter,
  normalisation, correlation or DA kernel file; it names `R/mod_omic_da.R` (KIT-194's module guard, not on S1's route,
  which calls the reference DA directly) and `R/utils_omic_gui_config.R`. Base results equal on both bases.

### Failing-first tests (A1)

R1e (52343cca; `r1e-tests.tsv`, `r1e-first-failures.txt`, `r1e-test-notes.txt`): every failing-first test fails by an
assertion, at its stated value, both omics; no test errored.

| Test | Base value observed |
|---|---|
| T1 | `failed_assays$pos` = "covariate contains NA or infinite values" |
| T2 | `failed_assays$pos` = "prior.weights contain NA values" |
| T5 | `failed_assays$pos` = "No residual degrees of freedom in linear model fits"; the direct call's condition is `simpleError` |
| T6, T6 (b), T7 | the manifest's errors: `pos` = "prior.weights contain NA values" in each of the three contrasts (T6 (b)'s and T7's preconditions held: design none with k 0 and `ruv_da_method` recorded; the reference's manifest present after the worker's five actions) |
| T8 | seed 196001 (the first; 5 features without values: F0015, F0029, F0040, F0062, F0067); "differential abundance failed for every assay / LCMS_Pos: prior.weights contain NA values" |
| T3, T4, T9 (a), T10 | pass (T10 in R1a and R1c) |
| T9 (b) | measurement: base differs from the reference, largest absolute difference in P.Value 1.05e-05 (both omics). 5.2's trend-off bullet stands and NEWS keeps the D2 sentence |

### The fix (R2 at a1e46859 and fef003c6)

- **A1, head.** test-kit196 passes 540, fails 0 (`r2-tests.tsv`). T1's rounding case was bit-identical (largest relative
  difference 0 in t, P.Value and the q-values; `r2-test-notes.txt`). T8 used seed 196001 on both commits. T9 (b) head equals
  the reference.
- **A3 holds.** Of the 38 files (section 9's scope at the branch point plus test-omics-norm-state-guards.R; listed in `r2-compare.tsv`), 37 are `same` base against head, test-kit196 is `only-on-head` and test-kit128 is `same`. No `compare-error`, no `new-failure`. Files failing identically on both sides (stored references
  recorded on the laptop): test-kit127 (1), test-kit129-run-path (6), test-omic-da-module (53), test-workflow-omics-declarations
  (22), test-workflow-operation-registry (2), test-workflow-saved-da (171).
- **A4 holds.** R CMD check: 0 errors, 3 warnings, 3 notes on both sides, findings identical, none only on head
  (`r2-check-compare.txt`).
- **A2 (ii) to (iv), S1 head side (`s1/r2-head/`).** (ii) every journey that completed on base (158 L, 150 M: 308 journeys) has identical `da` and long table digests on head: 616 of 616 digest comparisons (316 L, 300 M). (iii) every journey that failed on base (L 46: 20 full stops and
  26 one-assay failures; M 54) completes on head; in each of its assays the tested rows' logFC, t, P.Value and q equal the
  contrast half on the DA input without the features without values, identical in all 92 (L) and 54 (M) assays (largest
  relative difference 0), and the rows with `AveExpr` missing are exactly the features without values in every contrast.
  (iv) no head journey fails with either message (204 and 204).
- **KIT-194 on the merged tree.** T6 calls the reference DA wrapper and T7 the DA tab's run (`executeOmicDaRun`), both below
  KIT-194's module guard (`runOmicDaAnalysis`); both reach the DA on a1e46859 and pass. test-omics-norm-state-guards.R is
  `same`.

### D6 (descriptive only; VAL-116 (7); `s1/r2-head/d6-zero-df-*.tsv`, `d6-summary-*.tsv`)

Head side, every journey and assay, per contrast: the reference (the journey's DA input, TRUE, TRUE) against the same
halves without its zero-df features. No computation failed (1,020 rows).

- Zero-df features occur in 130 of 408 L rows (up to 10) and 294 of 612 M rows (up to 12).
- Over the features tested in both fits, the largest absolute change in P.Value is 0.087 (L) and 0.092 (M). In the rows without a zero-df feature (278 of 408 L, 318 of 612 M) the two fits are the same and the change is 0. In the rows with one (130 L, 294 M), the median of the per-row largest change is 0.0086 (L) and 0.0099 (M); it exceeds 0.01 in 58 (L) and 147 (M) rows and 0.05 in 5 and 13.
- q < 0.05 calls made only in the reference: 12 (L), 4 (M); only in the alternative: 16 (L), 2 (M); of 7,420 (L) and 1,515
  (M) reference calls. No zero-df feature was called in the reference.
- (f) also contains the alternative's smaller q-value family; only (e) isolates limma's robust step. Nothing here is
  decided here. **With the owner:** whether zero-df features (one value per observed group) should stay in limma's prior
  estimation is the owner's decision before KIT-188 declares (audit S1; F3), on these figures.

### R3 (full suite; run 45358036892, base fe5c9273, head 6dcf593c)

463 files: 462 `same`, test-kit196-da-no-value-features.R `only-on-head` (ok); no `new-failure`, `message-changed`,
`count-changed`, timeout or crash (`r3-compare.tsv`). R CMD check: findings identical on both sides, none only on head
(`r3-check-compare.txt`). 35 files fail identically on both sides (stored references, as in R2). 42 jobs, 1,205 job-minutes.

### R2b on the merge with staging 241914cb (audit C3), declared before sending

- Purpose: 6.4's merger rule (staging gained code after R2: KIT-203 and others; no 7.3, lipid/metab DA, filter,
  correlation or declaration file, so S1 is not repeated).
- Binds: head = the merge 51b62ca7 plus this records commit; base = staging 241914cb.
- Content: R2's 38 files (`r2-compare.tsv`'s list) and R CMD check, 3 shards per side, `file_timeouts` 120 for the two
  GUI-entry files.
- Ceiling: 6 x 60 + 2 x 90 = 540 job-minutes = 36.0 core-hours.
- Expected: every file `same`, test-kit196 `only-on-head`, check findings identical. A difference is a result: recorded and
  reported, not corrected in the same pass.

**Result (run 45376398677, head 73cd9ee9, base 241914cb):** as expected. 37 files `same`, test-kit196 `only-on-head` (ok),
check findings identical (6 on both, none only on head). 8 jobs, 266 job-minutes = 17.7 core-hours of 36.0
(`r2b-compare.tsv`, `r2b-tests.tsv`, `r2b-check-compare.txt`).

## Core-hours (job minutes x 4 / 60)

| Row | Used | Ceiling |
|---|---|---|
| R1 (incl. three instrument runs) | 12.4 | 13.0 |
| R2 (the first send on the merge with ea93d948, cancelled, no result) | 21.2 | 52.0 |
| R2b under 6.4, tests and check on the merge with KIT-194's staging (45336286628) | 18.6 | 36.0 |
| R2b under 6.4, S1 both sides on KIT-194's base (four runs) | 18.8 | 32.0 (16.0 per part pair as 6.4 allows, both sides) |
| R3 (with check) | 80.3 | 186.7 |
| R2b on the merge with staging 241914cb (audit C3; 45376398677) | 17.7 | 36.0 |
| Total | 169.0 | 490.4 (251.7 without a correction); expected 96.2 |

No correction pass was used.

## Acceptance criteria

| # | Verdict | Evidence |
|---|---|---|
| A1 | met | R1e (base values above); R2 head 540/0 |
| A2 | met | S1 (i) to (iv) above, at both bases |
| A3 | met | `r2-compare.tsv` |
| A4 | met | `r2-check-compare.txt` |
| A5 | met | R3: 462 `same`, the new file only on head; 6.2's verdicts and no other |
| A6 | met | NEWS under the unreleased heading (Fixes); ledger row `VAL-159` after VAL-155; `git diff --name-only` against staging: no DESCRIPTION, packaging, lock or NAMESPACE change; `man/` only the four pages of 7.3 |
| A7 | met | below; limma 3.68.4 on every runner |

**A7, sections 2.3 and 2.5 re-read at 3d4a585d.** The contrast halves call `limma::eBayes(cfit, trend = eBayes_trend,
robust = eBayes_robust)` at `R/func_lipid_da_model.R:115` and `R/func_metab_da_core_engine.R:115` (base); the per-contrast
kernel calls `binding$da_apply(..., TRUE, TRUE)` inside its `tryCatch` (`R/utils_workflow_omics_declarations.R:773-780`) and
the manifest stops only when every contrast failed (`:1603-1605`); `applyOmicsRuvFactorsToDaInput` returns early without W
(`R/func_omics_ruv_factors.R:479`); the intensity filters count with `na.rm = TRUE` (`R/func_lipid_qc.R:66`,
`R/func_metab_qc_filtering_helpers.R:34`); the exported S4 methods still call `runTestsContrasts`
(`R/func_lipid_s4_da_analysis_methods.R:156`, `R/func_metab_s4_da_methods.R:160`), unchanged (KIT-200). Proteomics:
`eBayesFn(..., trend = eBayes_trend, robust = eBayes_robust)` at `R/func_prot_da_model_stats.R:306-307`, trend and robust
fixed TRUE on every DA route including Skip (2.5 corrected in revision 1.2; KIT-200).

## Departures from the ticket's text (recorded)

- **T5's direct call** uses `tryCatch(..., error = function(e) e)` and `expect_s3_class(e, "omics_da_no_observed_features")`
  instead of `expect_error(class =)`, because `expect_error(class =)` lets a different error through as a test error, which A1
  forbids on base. Same claim.
- **T7** drives the DA tab's run (`executeOmicDaRun`, what the module's run button calls) under `shiny::testServer` with the
  managed chain (`omicsTestOwner`, `omicsManagedChainMocks`, `installWorkflowOmicsActions`), as the helper `omicsTabServer`
  does, instead of the module server: the module driver set no state selection (instrument defect). It advances the
  runner's poll and reads the accepted manifest as a state value. The head is `<prefix>_correlation_filtered`, as the ticket
  says.
- **T6 and T6 (b)** also assert the manifest's error messages as text (`unlist(errors)` empty), the same claim as the
  ticket's list assertion, so that the base value is shown.
- **S1**: parts L and M in separate requests (run-tests runs a request's scripts in one job, one after another); in R2 the
  base and head sides in separate jobs (`--base-only`, then `--base-records`), because both sides exceed one 120-minute
  script; the head side reads the base side's records of the same base commit. Child processes: one per unit (R01), native
  code compiled once per side with `pkgbuild::compile_dll(force = TRUE)` after `src/build_dia_rsum.R` as the run-tests jobs
  run it (R02).
- **R3** ran with `check: true` to cover R2b.
- **Ticket text**: A2 (i)'s precision is 7 significant digits (KIT-188's rows), not 8; the scatter's positive residual df are
  7, 8 and 9.

## Inferred, not observed

- That the managed tab path and the in-process path give the same tables (T6 and T7 assert the same values separately; they
  were not compared with each other).
- That replay of an old record whose DA failed for this cause now completes (by reading: replay re-runs the same contrast
  halves; not run).

## Hand-offs

- **KIT-201 (audit I04).** At its next revision, add the guard test: **T-196 (guard; TREAT with features without values;
  passes on base and head).** `kit196Object(omic, scatter = TRUE)` (copied from `test-kit196-da-no-value-features.R`) through
  `runLipidsDA` / `runMetabolitesDA` with `treat_lfc_cutoff = 0.5`: completes; `F1_empty` once per contrast with every
  statistic NA and `AveExpr` NA; the tested rows equal the same call on the object without `F1_empty` (`expect_identical`).
  And replace `KIT-201.md:683`'s "not tested here" with "tested by T-196".
- **KIT-200 (audit I06).** Detection here is by row names (`rownames(cfit$coefficients)[is.na(cfit$Amean)]`); on a fit
  without row names the index is NULL, nothing is removed and limma stops as before. Every app route sets row names, and the
  roxygen requires them. KIT-200, which copies the form, should use a logical index (`empty <- is.na(cfit$Amean)`). Also the
  class and text, the append form, I07 (bit for bit shown under the reference BLAS only) and the corrected NEWS wording.

## Audit notes recorded (I07, I09 to I12)

- I07: "bit for bit" is shown under the reference BLAS of the runs only (`libblas.so.3.10.0`).
- I09: T7 starts below the module's button handler (`runOmicDaAnalysis`), whose checks do not depend on feature values;
  bf5ed69f's added assertion (manifest errors as text) is the same claim.
- I10: D3's "the failure channel adds the assay" holds for the orchestrators; in the app a full stop joins messages without
  assay names (`R/utils_workflow_omics_declarations.R:1604`; F7).
- I11: 122 added lines exceed R.md's 100-character limit, following the files' own idiom (~120); the four `man/*.Rd` are
  hand-edited (the brief forbids `roxygen2::document()`), not regenerated.
- I12 (with the owner): on the app route `features_without_values` reaches neither the manifest nor the record; the user sees
  the rows and the "Not tested" count, which lumps them with features missing a compared group. A reporting candidate.

## What other tickets observed

- KIT-188's lifecycle check 45334020028 on staging ea93d948 recorded one-assay DA failures as `no_da` in its lipid shape L12
  at 20 % missing (arms LA, NO and the oracle), where its earlier lifecycle had recorded them `ok`. Consistent with the stop
  this ticket fixes, but the message was not recorded there: an observation, not evidence here (coordinator's note).

## Not checked

- Windows (no platform code); real data (none used).
- The managed path against the in-process path on the same object (see Inferred).
- Replay of an old failed record (see Inferred).
- The exported S4 `differentialAbundanceAnalysis` on lipid and metabolomics objects (unchanged; KIT-200).
- TREAT (lfc > 0) with features without values: not tested here; KIT-201's guard T-196 tests it.
- An assay of fewer than three features with one without values (the disclosed `B` change): not tested.
- D6's figures decide nothing here.


## Merge record (coordinator)

Merged into `v0.9-staging` on 2026-11-24 from 38c267b2 after the implementation audit and its re-check 1 (register `AUDIT-KIT-196-impl-20261124`: I01 to I12 closed or recorded; mergeable after one run, or now with a recorded deviation from 6.4). The ledger row is VAL-159. **Deviation from ticket section 6.4, recorded by the coordinator:** staging gained KIT-207 (the Summary guard for lipidomics, metabolomics and transcriptomics, merge 2ead4fec) after R2b's base 241914cb, and 6.4 asks for R2b again on that merge. It was not run. The re-check showed that KIT-207 shares no file with this ticket; that its guard reads only state names and lineage and no DA result; that everything the Summary writes comes from the pre-DA head state, so the appended untested rows cannot reach it; that KIT-207's fixtures cannot produce a feature without values (at most 3 of 6 cells removed per feature and group; one assay, one contrast) and that its DA rows ran green on trees without this fix, where such a feature would have stopped the DA; and that this ticket's tests call no Summary function. R2b (run 45376398677, base 241914cb) and R3 (run 45358036892) therefore stand; the validation harness's whole suite, which drives the lipid and metabolite DA, runs next on the staging commit that holds both. Not verified: whether the report bundle and the session state embed DA tables (they could then carry the appended rows, for data whose DA used to fail). Ticket total 169.0 of 490.4 core-hours. With the owner before KIT-188 declares: whether zero-df features stay in limma's prior estimation (D6, pool X151).
