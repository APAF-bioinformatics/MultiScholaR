# KIT-221 outcome (PARTIAL: stopped after R1 at a missed prediction)

Base: `origin/v0.9-staging` 176f9d3ff. Tests-only commit 41c69770f (nothing under `R/`). The fix is drafted in the
working tree of the worktree, uncommitted and never run.

## R1 (run 45713440865, request kit221-r1-41c69770, one test job and one script job)

Observed (tables in `r1/`): `test-kit221-da-several-comparisons.R` 43 passed, 38 failed, 0 errors, 0 skipped (predicted
42 / 39); `test-kit218-da-wrapper-tail.R` (G4) 80 passed, 0 failed, 2 skipped (G1 and M2 skip because the fix is
present); `g1.R` exit 1 after 2 s.

| Test | Predicted passed / failed | Observed | Assertion that failed |
|---|---|---|---|
| W1 | 2 / 8 | as predicted | first: `sort(unique(rows$comparison))` is `B_vs_A` only; then 60 not 120 rows (twice), `num_sig_da_molecules` and the first-go table `B_vs_A` only, 1 volcano row not 2, no `log2FC:C_vs_A` column, long table `B_vs_A` only. No error |
| W1b | 0 / 1 | as predicted | `expect_null(run$error)`; M1 = "number of items to replace is not a multiple of replacement length" (the ticket's reading) |
| W2 | 8 / 2 | as predicted | comparisons set; C_vs_A rows (0 against 60); B_vs_A identical |
| W3 | 8 / 5 | as predicted | `length(contrasts_results$results)` 1 not 2; comparisons; 60 not 120 rows; set; C_vs_A rows |
| W4 | 10 / 10 | as predicted | as W1 plus set and C_vs_A rows |
| W5 | 9 / 5 | 10 passed / 4 failed (the `length` assertion passes; see below) | comparisons `B_vs_A` only; `nrow` n not 2n; set; C_vs_A rows |
| W7 | 5 / 5 | as predicted | set; B_vs_A rows; `log2FC`, `raw_pvalue`, `fdr_qvalue` of the rows labelled C_vs_A differ from the one-comparison call; mean `log2FC` of P11..P16 under `C_vs_A` = 1.0322 (B - A's effect: the factor premise holds; stop rule 6 does not apply) |
| W8 | 0 / 1 | as predicted | `expect_null(run$error)` |
| W9 | 0 / 2 | as predicted | both calls stop with "incorrect number of dimensions", not of class `prot_da_contrasts_invalid` |

**Missed prediction (W5, one expectation).** I predicted that `length(multi$value$contrasts_results$results)` is 1
on base. It is 2: `differentialAbundanceAnalysis()` keeps every comparison in `contrasts_results` (ticket 2.1: only
the derived tables are first-only), which the ticket states and I mis-applied to W5 (W3's assertion is right: the
sibling's `contrasts_results` holds one). The row therefore passes on base; the other four W5 expectations fail as
predicted. Per the ticket (a missed prediction: stop and report) no further run was sent and no fix was committed.

**Instrument defect (g1.R).** `ProteinQuantitativeData` is not exported, so `g1.R` (run through the installed package)
stops at the first call: "could not find function". No digest was produced, no value judged. The mechanics fix is
`ns$ProteinQuantitativeData(...)` in `make_object`; not sent (the stop above), so none of the three instrument reruns is used.

## Not checked

Everything after R1: the fix's behaviour, G1 and G3 digests, R2 and R3, the check comparison, NEWS, ledger row.
