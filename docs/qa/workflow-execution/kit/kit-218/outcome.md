# KIT-218 outcome (interim: tests-only stage)

Ticket r1.2. Branch `kit-218` from `v0.9-staging` f93dfd85. The fix commit waits until KIT-200 has merged into staging.

## Step 0 (against f93dfd85)

- Ticket already r1.2: Re-check 1's R01 to R03 applied; the owner's decisions are ledger row VAL-188 on staging (e6791275). Nothing duplicated.
- Every cited line of `R/func_prot_da_model_wrapper.R` holds (`:150-154`, `:218-219`, `:222-233`, `:267-276`, `:304-312`, `:464`, `:499`, `:508-515`, `:522-531`); `R/func_prot_da_results_summary.R:245`, `:254-263`, `:317-322`; `R/func_prot_da_compat.R:115`, `:125`, `:132-133`, `:162-163`, `:597-601`; `R/func_prot_da_model_stats.R:823-824`; `R/func_prot_da_model_methods.R:427`; `NAMESPACE:109`.
- Moved since 67064baf (KIT-208's merge ce59ad23), cited only in hand-offs, not in this ticket's code: `R/func_prot_da_results_methods.R` published count table, `:314-340` is now `:324-350`. `R/utils_workflow_dia_da.R:308-318` holds.
- Open branches editing the wrapper file: only `merge-rehearsal`, `sim-k124-06`, `sim-k124-confirm` (predate KIT-126; not tickets), as the audit found.
- R2 grep at the branch point: the same 8 files as the ticket lists.

## R1 (run 45491079123, at debed95a)

**Binding deviation (stated).** The ticket binds R1 to "the staging commit that holds KIT-200". KIT-200 is not yet in staging, and W2 and G3 need its kernel. R1 therefore ran on run-only branch `kit218-r1-with-kit200` = the tests-only commit 23db9b69 merged with KIT-200's head 3f8203c2. The ledger merge conflict (docs only) was resolved to staging's text. This branch is never merged. R2's base will be the tests-only commit once staging (with KIT-200) is merged into `kit-218`.

**Observed.**
- `test-kit218-da-wrapper-tail.R`: 13 passed, 4 failed. W1, W2, W3 and W4 each failed at `expect_null`, with the caught message ending "Caused by error in `rownames_to_column()`: ! is.data.frame(df) is not TRUE" (purrr's two wrappers above it). This is the stated base value, so stop rules 1 and 2 do not fire.
- G1 completed on every input: W1 120 rows (B_vs_A, C_vs_A), W2 61, W3 60, W4 60 rows with comparison `B_vs_A` only. That is the first-row truncation on base that P02 read. Stop rule 3 does not fire.
- G3 passed (class `omics_da_no_observed_features`).
- G2: `test-prot-07d-da-model-stats-shared.R` 60 passed, 0 failed, 1 skipped.
- M1 (recorded, not judged): with groups `1`, `2`, `3` the caught message was `<text>:1:2: unexpected input / 1: 2_`. This is a parse error of the contrast name `2_vs_1`. The ticket expected a stop "inside `runTestsContrasts`"; where the parse happens was not traced.

**Inferred.** The nesting is the only stop on these four inputs on base; the one-line fix plus OD-4 (a)'s one token is sufficient for them.

Files: `r1/tests.tsv`, `r1/failures.txt` (its message fields are empty: the runner's parser drops multi-line messages; the messages are in `r1/kit218-log.txt`).

## Not checked (so far)

- The fix and R2 (waiting for KIT-200's merge).
- That KIT-200's merged kernel equals its head 3f8203c2, on which R1 ran.
- Where M1's parse error is raised.
