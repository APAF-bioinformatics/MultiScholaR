# KIT-190 R2 after the implementation-audit repairs: predictions (written before the run)

Run `k190-r2-audit-<head>`: the head is the commit after this file is committed (it holds the I01 code ea5c71ad, the
I07 re-dump 09c08b5e and the merge of staging 828dd105, f5ea9cc2). Base: staging 828dd1051d2acab075d6d04f887b7c9e9ccaa48e.
Scope: the 103 files of the final R2 (45619021147) plus `test-workflow-omics-declarations.R`,
`test-gui-entry-lipidomics.R`, `test-gui-entry-metabolomics.R` (I10); the I01 closure files are already among the 103.
106 files, 4 shards, R CMD check, D5 control (`norm-dump.R`, base 828dd105). Ceiling 700 job-minutes. One attempt;
an instrument defect falls under VAL-54. No correction pass is left (outcome.md, "The correction pass").

| File | Prediction | Why |
|---|---|---|
| `test-ruv-control-anova-grouping.R` | only-on-head, all passing (T9 (iv) included) | the new file; I01 code |
| `test-omic-norm-module.R` | same | the 12 provenance values re-dumped (I07); the cap does not bind in its fixtures, so the annotation leaves `normdata` unchanged and no notice is logged or shown |
| `test-workflow-omics-ruv-cancor-once.R` | new-failure at `:73`, exactly as in run 45619021147 | not edited (I06, owner item O1) |
| `test-workflow-omics-binding-steps.R` | message-changed | the KIT-094 snapshot gains the two markers and four outcome fields (`VAL-<merge-b>` draft) |
| `test-prot-norm-core.R` | message-changed | lipid and metab RUV pinned digests, failing on base too (VAL-213) |
| `test-workflow-omics-declarations.R` | message-changed or count-changed (failing on both sides) | on base already moved by KIT-224; on head the registry `fixed` of `lipid_ruv` / `metab_ruv` gains `ruv_control_grouping` and `ruv_k_bound`, and the journeys' Manual RUV states and what follows them move (the control rule). If it passes on base, the result is a new-failure and the prediction failed; it is reported, not repaired |
| `test-gui-entry-lipidomics.R`, `test-gui-entry-metabolomics.R` | same | no RUV value pinned; the cap does not bind |
| `test-module-ci-lipid-norm.R`, `test-module-ci-metab-norm.R`, `test-lipid-12a-norm-module-characterization.R`, `test-metab-03d-norm-module-characterization.R`, `test-workflow-lipid-declarations.R`, `test-workflow-metab-declarations.R` (I01 closure) | same | uncapped summaries keep "Best k: …"; no notice where the cap does not bind |
| every other file of the 103 | same | as in 45619021147 |
| R CMD check | no finding only on head | |
| D5 control | 12 dump lines differ, all provenance digests; the head values equal the 12 stored values of I07 | the markers; nothing else moves |

Budget before this run: 220.39 of 510.4 used (coordinator's figure), plus the I01 failing-first run (1 job).

## The O1 run (VAL-280), declared before it is sent; sent only after run 45651851022 is judged

The owner decided O1 (a) (VAL-280, staging 3a6794b0). Commit 6b465f1e (tests only) replaces `:73` of
`test-workflow-omics-ruv-cancor-once.R` with the register's comment and `expect_true(isTRUE(new$one$success))`, and adds
a printed observation ("KIT-190 observation <omic> <assay>: <percentage>% <status> k=<k> reason=<error_reason>", and
the assay's refusal text where it failed). That observation changes no expectation and reads only the sweep's existing
`optimization_results` and `error`; nothing under `R/` changed. Run 45651851022 (head 687ba6f3) does not contain this
edit, and there `:73` is still predicted to fail as in 45619021147.

Run `k190-o1-<head>`: this file alone, base 687ba6f323eb05692afa8074b8664648bed6b4ef (without the edit) against the head
commit that holds this declaration (the test file is as in 6b465f1e). No check, no scripts, 1 shard, ceiling
40 job-minutes, one attempt (VAL-54 for an instrument defect only).

| Side | Prediction |
|---|---|
| base 687ba6f3 | 14 passed, 2 failed: `:73` for lipid and metab (assay `two` unsuccessful), as in 45619021147 |
| head | 16 passed, 0 failed (3 tests per omic: 1 + 3 + 4 expectations) |
| comparison | `fixed` (failing on base, passing on head) |
| observation | assay `one` accepted at one or more percentages; assay `two` rejected at every percentage, with the reason printed per percentage (recorded as observed; it decides nothing; a separate measurement ticket owns it) |
