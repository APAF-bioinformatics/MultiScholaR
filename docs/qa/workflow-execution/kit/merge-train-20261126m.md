# Merge train 20 of 2026-11-26: KIT-258

**This train is STACKED on merge train 19 (KIT-264), which is stacked on train 18.** It was branched from
`origin/merge-train-20261126l` at its tip `b6b4ed6b76b37e2c96f00ff22f622c9ffa0149b2` (fetched before building; it
carries the KIT-264 merge 9f65a58c and train 19's declaration). Train 18 is landed on staging as 15d86bd4 with trains 15
to 17. **Train 20 lands only if train 19 lands with its code unchanged.** If train 19 is rebuilt or its code changes,
this train is rebuilt and its check repeated. Train 19's own run is not touched by this train.

## Train

Base: `b6b4ed6b76b37e2c96f00ff22f622c9ffa0149b2` (train 19's tip; the code is that of 9f65a58c).

| Order | Ticket | Branch head | Merge commit | Ledger row | Audit |
|---|---|---|---|---|---|
| 1 | KIT-258 | `origin/kit-258` 216ba36b | 77b67929 | `VAL-337` (replaces `VAL-<merge>`) | AUDIT-KIT-258-impl-20261126 and Re-check 1 (verdict ready; I11 = this plan) |

Head (bound by the run): `77b67929b5031dc3b59ddbd5bc13f5b1ab4f112c`. (The request file holds it as well, as
this file is in the declaration commit after it, under `docs/` only.

### Conflicts met and their resolutions

- **`R/utils_workflow_analysis_record.R`** (two hunks), resolved as register I11 says: the base's rows kept in its order
  (`missing_value_comparison`, `enrichment_background`, `intensity_filter_samples`, `duplicate_resolution_samples`,
  `ruv_feature_levels`, `treat_moderation`, `ruv_control_grouping`, `ruv_k_bound`, `intensity_filter_missing`), then
  `duplicate_resolution_rank` last; the list closed once. The base's closing `))` after `intensity_filter_missing`
  became `),`. Header: "; KIT-258 adds `duplicate_resolution_rank`" appended after the KIT-217 clause. The resolved
  tail (braces checked by reading: `list(` of the row opens and closes on its own line pair, the outer `list(` closes
  once, on the last row):

```r
    intensity_filter_missing = list(option = "MultiScholaR.intensity_filter_missing", current = "group_counts",
        earlier = "not_counted"),
    # KIT-258 (VAL-254, VAL-284): among the rows of one duplicated id, the most measured row (a value neither missing nor
    # zero), then the highest mean of its measured values ("measured_then_mean"); records made before ranked by the mean
    # of the non-missing values alone ("mean").
    duplicate_resolution_rank = list(option = "MultiScholaR.duplicate_resolution_rank", current = "measured_then_mean",
        earlier = "mean"))
```

  Header line: `#' replay, `workflowDaReplayRuleOptions`; KIT-217 adds `intensity_filter_missing`; KIT-258 adds `duplicate_resolution_rank`): for`.
- **`docs/methods-decision-ledger.md`**: the base's appended rows (through `VAL-336`) kept in order, then
  KIT-258's row as `VAL-337`. Partly numbered rows inherited from the base are NOT repaired here.
- **`NEWS.md`**: merged without conflict (git auto-merge); not re-read for order, the coordinator orders at landing.
- `kit-258/outcome.md` line 180 names the placeholder; it now reads `VAL-337 (was VAL-<merge>)`.
- No conflict under `R/` other than the file above, nor under `tests/`, `inst/`, `man/`, `DESCRIPTION` or `NAMESPACE`.

### Path lengths

Checked before sending: every path new since the base is at most 72 bytes with the `MultiScholaR/` prefix
(`docs/qa/workflow-execution/kit/kit-258/r2/check-compare.txt`).

## Declaration (written before the run is sent)

**Purpose.** KIT-258's merge check on the merged tree, and the first execution of its test E7.

**Scope recount.** The ticket's grep (`git grep -lE "resolveDuplicateFeaturesByIntensity|resolveLipidDuplicateFeaturesByIntensity|
resolveMetabDuplicateAssayData|resolveDuplicateFeatures[(,]|resolveDuplicateFeaturesForLipidObject|resolve_duplicates|
duplicate_resolution|workflowReplayRuleMarkers" -- tests/testthat`) on the base b6b4ed6b gives 32 files (29 tests and 3
helpers) and on the head 33 (30 tests: the new KIT-258 file). The ticket's 24 are now 29 because train 15 and earlier trains
added `test-kit201-treat-moderation.R`, `test-kit217-intensity-filter-missing.R`, `test-kit224-ruv-feature-levels.R`,
`test-kit252-managed-failure-text.R` and `test-ruv-control-anova-grouping.R`. Added beyond the grep, as the ticket and I11
declare: `test-workflow-omics-binding-steps.R`, `test-workflow-omics-residency.R`, `test-gui-entry-lipidomics.R`,
`test-gui-entry-metabolomics.R`; and, as the new lipid and metabolomics record tests of trains 13 to 15,
`test-kit254-omics-correlation-skip.R` (already in the grep's reach of train 15; `test-kit217-...` and
`test-ruv-control-anova-grouping.R` are in the grep). Scope = **35 test files**, listed in the predictions.

**Run settings.** `MULTISCHOLAR_TRX_MANIFEST_NOTE=1` on both sides; `test-prot-norm-core.R` is not in scope, so
`MULTISCHOLAR_KIT105_W3_DUMP` is not set; extra package `sva` (the declared 139/0 of
`test-ruv-control-anova-grouping.R` holds only with it); no `limpa`, no `mixOmics` (none of the 35 needs them, as KIT-258's R2);
reference BLAS; 3 shards per side; file timeout 60 minutes; check on both sides; no scripts.

**Ceiling.** 600 job-minutes = 40.0 core-hours: KIT-258's R2 row (ticket section 11: 6 x 70 + 2 x 90), its merge-check
ceiling. Expected about 300 job-minutes (about 20 core-hours).

### Predictions

Counts are passed/failed. Where "read from run" is shown the base count was not found in the trains' declarations; the
prediction is the stated verdict and the count is read from the run.

| File | Base | Head | Verdict |
|---|---|---|---|
| test-gui-entry-lipidomics.R | 19/16 (train 15) | same as base | same |
| test-gui-entry-metabolomics.R | 57/48 (train 15) | same as base | same |
| test-kit201-treat-moderation.R | read from run | same as base | same |
| test-kit202-intensity-filter-samples.R | read from run | same as base | same |
| test-kit214-sample-columns.R | 176/0 | 176/0 | same |
| test-kit217-intensity-filter-missing.R | 344/0 | 344/0 | same |
| test-kit224-ruv-feature-levels.R | read from run | same as base | same |
| test-kit252-managed-failure-text.R | read from run | same as base | same |
| test-kit254-omics-correlation-skip.R | 376/3 | 376/3 | same |
| test-kit258-duplicate-completeness.R | absent | 162/0/0/0 | only-on-head |
| test-lipid-01-qc-filtering-helpers.R | read from run | same as base | same |
| test-lipid-01c-qc-duplicate-module-seams.R | read from run | same as base | same |
| test-lipid-01g-qc-duplicates-module-characterization.R | read from run | same as base | same |
| test-lipid-03-s4-duplicate-helpers.R | read from run | same as base | same |
| test-metab-01b-qc-filtering-helpers-characterization.R | read from run | same as base | same |
| test-metab-01m-s4-resolve-duplicate-intensity-characterization.R | read from run | base count + 2 (source-pin test gains two expectations) | count-changed |
| test-metab-01p-qc-duplicates-module-characterization.R | read from run | same as base | same |
| test-metab-01u-qc-duplicates-module-characterization.R | read from run | same as base | same |
| test-metab-02u-resolve-duplicate-features-characterization.R | read from run | same as base (one pin replaced by one) | same |
| test-metab-s4-qc-shared.R | read from run | same as base | same |
| test-module-ci-lipid-qc.R | read from run | same as base | same |
| test-module-ci-metab-qc.R | read from run | same as base | same |
| test-omics-audit-carry-overs.R | read from run | same as base | same |
| test-omics-norm-state-guards.R | read from run | same as base | same |
| test-omics-qc-core.R | read from run | same as base | same |
| test-ruv-control-anova-grouping.R | 139/0 (with sva) | 139/0 | same |
| test-workflow-lipid-declarations.R | read from run | same as base | same |
| test-workflow-metab-declarations.R | read from run | same as base | same |
| test-workflow-omics-application-reads.R | read from run | same as base | same |
| test-workflow-omics-binding-steps.R | 57/2 (train 15) | 57/2 (messages move) | message-changed |
| test-workflow-omics-declarations.R | 43/25 (train 15) | 43/25 | same |
| test-workflow-omics-residency.R | read from run | same as base | same |
| test-workflow-omics-worker.R | read from run | same as base | same |
| test-workflow-operation-registry.R | read from run | same as base | same |
| test-workflow-partitions.R | read from run | same as base | same |
| R CMD check | | | identical findings, nothing only on head |

Reasoning for the three moves. `test-kit258-duplicate-completeness.R`: 142 (R2) + 20 for E7 (seven tests: 12 on the
declared step, 4 on the two kernels, 4 on the S4 methods) = 162/0 (Re-check 1's correction of I11's 149). The 01m file:
its source-pin test pins the old code on the base and the new code on the head, two expectations more on the head.
`test-workflow-omics-binding-steps.R`: the step's recorded parameters gain `duplicate_resolution_rank`; the counts stay at
train 15's 57/2 and only the messages of the already failing expectations move (R2 showed the head section with two
parameters). `test-kit254-omics-correlation-skip.R`: 376/3 on both sides, because this base does not contain train 15's
tests-only edit that is on staging (ledger VAL-325; on staging the file reads 379/0): the three failures are the known
train-15 state, so the verdict is `same` at 376/3.

### Decision and attempts

- Every prediction holds: TRAIN GREEN, landing only after train 19 lands with its code unchanged.
- Any other verdict is a result: recorded with its files and counts, not repaired or re-sent. A failure inside E7 or a parse
  failure of the KIT-258 file: stop and report; no code change.
- One run; one re-send only for an instrument defect with no result judged. An HTTP 500/502 at dispatch with no run
  created: started once by hand with the README's command, recorded.

## What the coordinator resolves at landing

- The ledger: `VAL-337` takes its number; its `<id>` and `<t>` fill; the full stop before "Evidence: R1" (I05,
  "(U1, routed separately). Evidence: R1"); the stored-reference row for the KIT-073 registry and KIT-094 binding-steps
  snapshots (I06, VAL-259 form); the earlier trains' placeholders as left.
- `NEWS.md`: KIT-258's bullet in its place after the landed trains' bullets.
- `kit-258/outcome.md` line 180 placeholder; the ticket's r1.3 line (I07).
- Stacking: lands only with train 19, whose own run is another agent's.

## Result

Run 45712155865 (request `merge-train-20261126m-77b67929`, request commit 880d1ff3), dispatched without error; green
(infrastructure). Head 77b67929 vs base b6b4ed6b. Tables: `merge-train-20261126m-compare.tsv`,
`merge-train-20261126m-check-compare.txt` (same directory). One deviation from the declaration, in the request only:
the two GUI-entry files (lipidomics, metabolomics) carried a 120-minute file timeout (the other files 60), as in
train 15; no file came near either limit.

**Job-minutes 210.9 = 14.06 core-hours** against the ceiling of 600 job-minutes (40.0 core-hours): tests head
12.6/42.2/33.1, tests base 18.2/27.8/48.3, check head 15.2, check base 12.6, plan 0.4, results 0.5.

Totals: 35 files, 32 `same`, 1 `only-on-head`, 1 `count-changed`, 1 `message-changed`, no `new-failure`, no `fixed`.
Tests head 5543 passed / 96 failed; base 5379 / 96.

| File | Predicted | Observed base | Observed head | Verdict |
|---|---|---|---|---|
| test-kit258-duplicate-completeness.R | only-on-head, 162/0 | absent | 162/0/0/0 | only-on-head: HELD (E7's first execution passes) |
| test-metab-01m-s4-resolve-duplicate-intensity-characterization.R | count-changed, base + 2 | 17/0 | 19/0 | count-changed: HELD (source-pin test 8 to 10) |
| test-metab-02u-resolve-duplicate-features-characterization.R | same | 17/0 | 17/0 | same: HELD |
| test-kit214-sample-columns.R | same, 176/0 | 176/0 | 176/0 | same: HELD |
| test-workflow-omics-binding-steps.R | message-changed, 57/2 | 57/2 | 57/2 | message-changed: HELD (KIT-094 snapshot expectation) |
| test-workflow-omics-declarations.R | same, 43/25 | 43/25 | 43/25 | same: HELD |
| test-kit254-omics-correlation-skip.R | same, 376/3 | 376/3 | 376/3 | same: HELD |
| test-gui-entry-lipidomics.R / -metabolomics.R | same (19/16, 57/48) | 19/16, 57/48 | 19/16, 57/48 | same: HELD |
| test-kit217 344/0; test-ruv-control-anova-grouping 139/0 | same | as stated | as stated | same: HELD |
| the other 24 files | same | read from run | same | same: HELD |
| R CMD check | nothing only on head | 0 E / 3 W / 3 N | 0 E / 3 W / 3 N | identical; only on head: none, only on base: none: HELD |

**Decision.** Every prediction holds. TRAIN GREEN for KIT-258, landing only after train 19 lands with its code
unchanged.

## Merge record (coordinator)

Landed on v0.9-staging on 2026-11-26 by a no-fast-forward merge of this branch at 3531caae onto d9aef593 (merge train 19 landed).

- Verdict: green. Run 45712155865 (35 files and R CMD check; 210.9 of 600 job-minutes): the new file 162/0 only on head (test E7 passing on its first train run); the 01m file 17 to 19 count-changed; the binding-steps snapshot message-changed at 57/2; 32 files same; check identical.
- The replay markers list after this landing ends: ruv_control_grouping, ruv_k_bound, intensity_filter_missing, duplicate_resolution_rank.
- Ledger at landing: staging's rows kept through VAL-336; KIT-258 = VAL-337 (the full stop before "Evidence" added, AUDIT-KIT-258-impl I05); the stored-reference row for the KIT-073 registry and the KIT-094 binding-steps snapshot = VAL-338 (I06, in the VAL-259 form).
- Owed after landing by the owner's decision VAL-306: the four in-app texts that still state the earlier rule.
- tests/testthat/test-kit254-omics-correlation-skip.R read 376/3 on both sides on this train's base; staging reads 379/0 (VAL-325).
