# KIT-218 outcome

The ticket is r1.4. Branch `kit-218` came from `v0.9-staging` f93dfd85. Staging dd769671 (with KIT-200 merged) was merged in at 9c586b58.

Commits (the audit repairs follow in the section "Implementation audit repairs"):
- 23db9b69: tests only (W1 to W4, G1, G3, M1).
- 8fe6057e: R1 records.
- e04d0c02: tests only (W5 to W7, M2; H1).
- 9c586b58: merge of staging dd769671. Nothing of KIT-218's is under `R/`. This is R2's base.
- 5c509642: the fix, NEWS and the ledger row `VAL-233`.

## Step 0

**Against f93dfd85.**
- The ticket was already r1.2. The owner's decisions are ledger row VAL-188, so nothing was duplicated.
- Every cited line holds in:
  - `R/func_prot_da_model_wrapper.R`: `:150-154`, `:218-219`, `:222-233`, `:267-276`, `:304-312`, `:464`, `:499`, `:508-515`, `:522-531`;
  - `R/func_prot_da_results_summary.R`: `:245`, `:254-263`, `:317-322`;
  - `R/func_prot_da_compat.R`: `:115`, `:125`, `:132-133`, `:162-163`, `:597-601`;
  - `R/func_prot_da_model_stats.R:823-824`;
  - `R/func_prot_da_model_methods.R:427`;
  - `NAMESPACE:109`.
- One line moved since 67064baf, cited only in a hand-off: KIT-208 moved the published count table in `R/func_prot_da_results_methods.R` from `:314-340` to `:324-350`.
- Open branches that edit the wrapper file: only `merge-rehearsal`, `sim-k124-06` and `sim-k124-confirm`. All three predate KIT-126.

**At the merge 9c586b58 (R05).**
- `git diff --name-only debed95a 9c586b58 -- R/` lists 16 files: `func_lipid_qc.R`, `func_lipid_s4_normalization_methods.R`, `func_metab_qc_filtering_helpers.R`, `func_metab_s4_qc_methods.R`, `func_prot_da_model_stats.R` (KIT-200's later change: 6 lines added, 1 removed), `mod_integration_enrichment.R`, `utils_integration_enrichment.R`, `utils_integration_gui.R`, `utils_integration_record.R`, `utils_integration_workbook.R`, `utils_workflow_analysis_record.R`, `utils_workflow_integration_input.R`, `utils_workflow_integration_model.R`, `utils_workflow_lipid_declarations.R`, `utils_workflow_metab_declarations.R`, `utils_workflow_omics_declarations.R`.
- The wrapper file and the tail helpers did not change between f93dfd85 and the merge.
- R2's scope is the R2 grep at the merge: 10 files, plus `test-kit196-da-no-value-features.R`. KIT-201's and KIT-204's new files are not in staging, so the second of each pair to merge runs the other's file.

## R1 (run 45491079123, at debed95a)

**How R1 was bound.** R1 did not run on staging. KIT-200 was not in staging then, so R1 ran on the run-only branch `kit218-r1-with-kit200`. That branch is 23db9b69 merged with KIT-200's head 3f8203c2; its ledger conflict was resolved to staging's text. The branch is never merged.

**What R1 showed.**
- W1 to W4 failed at `expect_null`. The caught message ends "Caused by error in `rownames_to_column()`: ! is.data.frame(df) is not TRUE".
- G1 completed with 120, 61, 60 and 60 rows. On W4's input it gave `B_vs_A` only, which is the first-row truncation.
- G2 and G3 passed.
- M1 is recorded in the R2 section.

Files: `r1/`.

## R2 (run 45573395892): fix 5c509642 against base 9c586b58, 11 files, one shard per side, check

### Base side (the tests-only tree, with KIT-200 merged): it reproduces R1, as R05 requires

- W1, W2, W3 and W4 fail at 'is.data.frame(df) is not TRUE'. This is the same as R1.
- G1 completes with 120, 61, 60 and 60 rows. On W4's input it gives `B_vs_A` only. This is the same as R1.
- G3 passes, as in R1.
- **W5:** "Error in `contrasts_tbl[, 1]`: incorrect number of dimensions". This is the stated base value. Testthat counts it as an *error*, not a failure: `expect_error(class =)` re-raises an error of another class.
- **W6:** fails at `expect_null` with the nesting message. Here the table is named `C_vs_A=groupC-groupA`.
- **W7:** fails at `expect_null` with "Caused by error in `fit$coefficients[, coef]`: ! subscript out of bounds". This is the stated base value.
- **M2:** the mean `log2FC` of P11 to P16 (the proteins with a group-B effect) is **0.0764 under the label `B_vs_A` and 1.0322 under `C_vs_A`**. The factor swap is reproduced on base.
- **M1:** with groups 1, 2 and 3 the call stops with "<text>:1:2: unexpected input / 1: 2_". This is limma's `makeContrasts` parsing the name `2_vs_1`, inside `runTestsContrasts` (Re-check 2, RC9). The message is the same on head.

No difference from R1, so stop rules 1 to 3 do not fire.

### Head side

- **The new file:** 64 passed, 0 failed. G1 and M2 are skipped, as declared.
- **Verdicts:** the new file is `fixed`. The other 10 files are `same`. `test-prot-ruv-factors-in-model.R` fails 1 test identically on both sides (M29/M30, the legacy design's DA).
- **Check:** 0 errors, 3 warnings, 3 notes on each side; 6 findings identical; nothing only on head.
- **W2's printed count table (KIT-204 C14, not asserted):** `B_vs_A`, "Not significant", 61. All 61 proteins count as not significant, including P11 to P16. See "Not checked".
- **What the tests show about q-values.** On this fixture the q-values of B_vs_A are missing on base and head alike, because of I01 (the DA step's lost Benjamini-Hochberg fallback, `runTestsContrasts`, handed to its own ticket). So the tests show that the wrapper carries what the DA step returns, not that those q-values are right. Since the audit, `kit218ExpectCarried` also asserts that `log2FC` and `raw_pvalue` are present for every tested protein and prints how many `fdr_qvalue` are missing; `kit218ExpectTruth` checks each label's `log2FC` against the fixture's group means (W1, W4, W6, W7).

Files: `r2/` (`compare.tsv`, `check-compare.txt`, `check-summary.txt`, `tests.tsv`, `tests-base.tsv`, and the new file's base and head logs).

## Acceptance criteria

| # | Verdict | Evidence |
|---|---|---|
| A1 | met | W1 to W4 failed in R1, reproduced in R2 base. W5 to W7 failed at their values in R2 base; W5 shows as an error, not a failure (above). On head W1 to W7, G2 and G3 pass; G1 and M2 are skipped |
| A2 | met | G1 completed on W1 to W4's inputs (R1, R2 base) |
| A3 | met | R2: the new file `fixed`, every other file `same`, nothing only on head in the check |
| A4 | pending | shared full suite after merge |
| A5 | met | NEWS entry; ledger row `VAL-233`; `NAMESPACE`, `man/` and `DESCRIPTION` unchanged (the fix commit changes the wrapper, NEWS and the ledger only) |

Stop rule 4 does not fire: no other difference, and no stored reference moved.

Core-hours: R1 7.0 job-minutes (0.47). R2 tests 21.0 + 21.2 and check 16.1 + 13.2 = 71.5 job-minutes (4.77). Total 5.24 against 13.0 expected and 56.0 ceiling.

## Not checked

- **W2's count table (all 61 'Not significant'): cause established by the implementation audit (`audits/AUDIT-KIT-218-impl-20261125.md`, I01), by reading and reconstruction, not by a run.** On W2's input the comparison B_vs_A has no p-value at or above 0.95 (largest 0.9486; on the complete fixture 0.9475), so qvalue 2.44 stops inside its π0 estimate, and `runTestsContrasts`' error handler assigns the Benjamini-Hochberg fallback to a local copy (`R/func_prot_da_model_stats.R:686`, `:809`): `fdr_qvalue` is missing for every protein of that comparison (the same loss is in R2's head log of `test-prot-07-da-analysis.R`, lines 54-56). The wrapper carries those missing values unchanged; W1 to W7's equality checks on `fdr_qvalue` hold for B_vs_A as missing values. With the fallback working, P11 to P16 would have BH q-values of 6e-6 to 5.4e-3. A defect of the DA step, handed to KIT-249; not of this ticket's change.
- The limpa branch (F6), F7's disagreement, MetaboliteAssayData through the wrapper: not run.
- Whether the wrapper's numbers equal the application's (P10): not added.
- Windows; real data.

## For the change report

**What was wrong, and who was affected.**
- The exported scripted function `daAnalysisWrapperFunction()` gave its per-comparison result tables to the summary step one list level too deep.
- So every call stopped with "is.data.frame(df) is not TRUE" after the differential-abundance step had run, and returned nothing. This was established on complete data in R1 (45491079123) and R2 base (45573395892).
- Only scripted callers were affected. Nothing in the application, replay or `inst/` calls the function.
- History: the same nesting was in `deAnalysisWrapperFunction` at the tags `archive/v0.35` and `v0.3.6.3`, where by reading the workbooks' DA step stopped (not run).

**Code.**
- At base f93dfd85: `R/func_prot_da_model_wrapper.R:273-276` and `:306-309` (the nesting), `:154` (first row only), `:151` (factor passed as is).
- At fix 5c509642: `R/func_prot_da_model_wrapper.R:23-26` (the refusal, class `prot_da_contrasts_invalid`), `:156`, `:159` (read as text, every row), `:272` and `:302` (the list passed once).

**What changed in results.**
- Before the fix nobody had results: every call stopped, so no earlier number changes.
- On the inputs the tests cover, the function now completes, and its tables carry the DA step's `log2FC`, `fdr_qvalue` and `raw_pvalue` unchanged (W1).
- Three inputs would have completed silently and wrongly once the nesting was fixed. They are closed:
  - a data frame without `full_format` (base read its first row only: G1 gives `B_vs_A` only; OD-4 (a), W4);
  - a character vector (now refused, W5);
  - a comparisons column stored as a factor. On base, M2 recorded the swap: the six B-effect proteins average 0.08 under the label `B_vs_A` and 1.03 under `C_vs_A`, so one comparison's statistics stood under the other's name. Now they are read as text (W6, W7). The swap was latent: on every version that holds the nesting, a scripted call stopped with an error before any table was returned, so it was never delivered to a user (audit, question 2). Whether the versions before the nesting (`v0.75`, `v0.80`, `v0.9`) delivered it was not checked.
- In the tests' fixture the comparison B_vs_A carries no q-values, because of the DA step's lost fallback (KIT-249); the wrapper's counts for such a comparison show nothing significant until that ticket merges.

**Still stops or limited.**
- It still stops for bare expressions without `full_format` (F5), for digit-leading group names (F1, M1: limma's `makeContrasts` parse error) and, by reading, for limpa objects (F6).
- Its gene-labelled volcano list is absent (F2).
- Its first count table counts at q < 0.05 whatever threshold is set (F7).

**Evidence and its limits.**
- R1 45491079123 and R2 45573395892. The shared full suite comes after merge.
- One synthetic design: 60 proteins, 12 samples, 3 groups. Limma branch only.

## Implementation audit repairs (`audits/AUDIT-KIT-218-impl-20261125.md`)

- Staging 3eeaf100 merged at 248a75e6. The ledger conflict was resolved by keeping staging's rows through VAL-225, with KIT-218's row `VAL-233` after them.
- **Repairs at c7f105a2 (nothing under `R/`):**
  - I02: `kit218ExpectCarried` asserts that `log2FC` and `raw_pvalue` are present for tested proteins and prints how many `fdr_qvalue` are missing. The texts above and the ledger qualifier come from the register.
  - I03: `kit218ExpectTruth` checks each label against the group-mean truth. It is called in W1, W4, W6 and W7.
  - I04: the NEWS clause was replaced, and NEWS now says the swap was latent.
  - I06: the ledger's line ranges now match this file.
  - I09: recorded as open in `wording.md`.
- **Run A1 (45581173440, c7f105a2, head only, new file):** 80 passed, 0 failed, G1 and M2 skipped. `kit218ExpectTruth` passed everywhere it is called. Printed: `fdr_qvalue` is missing for 60 of 60 tested proteins of B_vs_A and 0 of 60 of C_vs_A. This is I01, seen directly. 10.8 job-minutes, 0.72 core-hours, from the correction allowance.
- **Core-hours in total:** 5.96 against 13.0 expected and 56.0 ceiling.
