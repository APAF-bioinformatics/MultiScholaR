# Merge train 6 of 2026-11-25: KIT-224, KIT-251

This follows the coordinator's instruction after train 5 landed (staging 63143aa3). The procedure is the same as for train 5
(`merge-train-20261125e.md`).

## Train

Base: `63143aa3618eae2729c6864ab10903928bec5071` (`origin/v0.9-staging`, fetched before building).

| Order | Ticket | Branch head | Merge commit | Ledger row | Audit |
|---|---|---|---|---|---|
| 1 | KIT-224 | `origin/kit-224` c2fc484c | cb1e70eb | `VAL-260` | AUDIT-KIT-224-impl-20261125 (sections 5, 6) |
| 2 | KIT-251 | `origin/kit-251` d7a66e3e | e48254be | `VAL-261` | AUDIT-KIT-251-impl-20261125 (section 6) |

The head bound by the run is `e48254be` (full hash in the request). This file and the records commit after the run are
under `docs/qa/` only. The rows are numbered by the coordinator at landing.

### Conflicts met

**KIT-224 conflicted in one code file.** `R/utils_workflow_analysis_record.R` conflicted at the replay-marker list
`workflowReplayRuleMarkers`: staging (KIT-214) and KIT-224 each added a last row. The train stopped and reported it.

The coordinator authorised resolution (a) exactly. Both rows are kept, staging's first; the roxygen line names both. No
other line of the file changes. This is AUDIT-KIT-214-impl I10's predicted resolution. The resolved hunk, against the
base:

```diff
-#' adds `intensity_filter_samples`; KIT-214 adds `duplicate_resolution_samples`): for
+#' adds `intensity_filter_samples`; KIT-214 adds `duplicate_resolution_samples`; KIT-224 adds `ruv_feature_levels`): for
@@ workflowReplayRuleMarkers <- list(
     duplicate_resolution_samples = list(option = "MultiScholaR.duplicate_resolution_samples", current = "design",
-        earlier = "numeric_columns"))
+        earlier = "numeric_columns"),
+    # KIT-224 (decision M8): the lipidomics and metabolomics legacy RUV-III-C correction keeps each feature's level
+    # ("restored"); records made before kept the shift by a constant per feature ("shifted").
+    ruv_feature_levels = list(option = "MultiScholaR.ruv_feature_levels", current = "restored", earlier = "shifted"))
```

- The roxygen line is 120 characters, the width the block already uses.
- The block is `@noRd`. No `man/` page is generated from it: grep of `man/` for `workflowReplayRuleMarkers` and
  "Replay rule markers" finds nothing. So no `man/` edit is needed.

**The ledger conflicted in both merges.**
- In merge 1, KIT-224's side also carried VAL-230 to VAL-233 with their status cells as they stood before train 4 landed
  (its staging merge was at e0747fc1). Staging's rows VAL-230 to VAL-259 were kept as they are on staging; KIT-224's
  stale copies were dropped. KIT-224's row was then added as `VAL-260`.
- In merge 2, KIT-251's row was added as `VAL-261`, with the status form of its register's I04.
- Each `outcome.md` placeholder was replaced likewise.

**No other conflicts.**
- `NEWS.md` merged without conflict, under the one unreleased heading; no entry mentions KIT-222's merged-samples defect.
- `DESCRIPTION`, `NAMESPACE` and `man/` are unchanged.

### Path lengths

Checked by script before sending: 66 paths are new since the base. The longest is 95 bytes with the `MultiScholaR/`
prefix, so none exceeds 100.

## Declaration (written before the run is sent)

- **Purpose:** each member's check on the merged tree, done once on the train.
- **Scope:** 79 files, base and head, with R CMD check on both sides. `MULTISCHOLAR_TRX_MANIFEST_NOTE=1` and `sva` are set
  on both sides (KIT-224's register and its R2).
  - KIT-224: R2's 62 files (`kit-224/r2/compare.tsv`), plus `test-gui-entry-transcriptomics.R`,
    `test-workflow-omics-worker.R` and `test-kit200-prot-da-no-value-proteins.R`.
  - KIT-251: its R2's 12 files (`kit251-r2-26bc37fb.json`).
  - Readers of the marker list beyond those scopes: a grep of `tests/testthat/` for `workflowReplayRuleMarkers` and for
    each marker's option name finds eight files. Three are not in the scopes above and are added:
    `test-kit202-intensity-filter-samples.R`, `test-kit205-enrichment-background.R` and `test-kit214-sample-columns.R`.
    The other five (`test-kit127-…`, `test-kit224-…`, `test-lipid-12c-…`, `test-metab-03a-…`,
    `test-prot-missing-value-counts.R`) are already in.
  - Run settings: 6 shards per side. The per-file limit is 120 minutes for `test-prot-norm-core.R` and
    `test-gui-entry-transcriptomics.R`; these are kill limits only.
- **Ceiling:** 1,340 job-minutes, the sum of the members' merge-check ceilings: KIT-224's R2b, 740 (ticket section 10), and
  KIT-251's R2, 600.
- **Base values:** several register values predate trains 3 to 5. Where they do, the prediction states the change from
  base to head, read on the base side of this run. The base's failure records for the two snapshot tests are train 5's head
  records: `merge-train-20261125e/` run 45598541744. That tree equals 63143aa3 outside the ledger and records. They
  already carry KIT-214's lines.

### Predictions, KIT-224 (register sections 5 and 6, with R2c's values under the trx flag)

**New and unchanged files:**
- `test-kit224-ruv-feature-levels.R`: only on head, 69 passed, 0 failed.
- `test-kit129-run-path.R` and `test-omic-norm-module.R`: `same`.

**Count-changed files** (from R2c, with no new failure):

| File | Passed | Failures |
|---|---|---|
| `test-lipid-12c-ruv-s4-characterization.R` | +8 (R2c 59 to 67) | none on either side |
| `test-metab-03a-ruviii-c-varying-characterization.R` | +8 (32 to 40) | none on either side |
| `test-kit127-omics-ruv-factors.R` | +4 (421 to 425) | its one failure the same on both sides |
| `test-integration-input-transcriptomics.R` | +1 | failures unchanged |

**Message-changed files:**

| File | What changes |
|---|---|
| `test-workflow-trx-declarations.R` | The two `transcriptomics_normalized_data.tsv` SHA-256 leaves fail on base only, so failed −2 and passed +4 (R2c 309/33 to 313/31). `top_de_genes.tsv` and every other failure are identical. |
| `test-module-ci-trx-summary.R` | The export file digest fails on base only: passed +1, failed −1. |
| `test-prot-norm-core.R` | Counts unchanged. Only the `lipid_ruv` and `metab_ruv` result digests move (4ed0ad71 to 9ada5791, 5065333a to cc4276ca). |
| `test-omics-summary-core.R`, `test-omic-summary-module.R` | Counts unchanged. Only the leaf `integration/transcriptomics_normalized_data.tsv` moves, in trx_main and trx_managed. |
| `test-workflow-omics-declarations.R` (deferred snapshot) | Counts unchanged (45/23 on the base). On top of the base's records, which already hold KIT-214's six lines, the head adds only KIT-224's lines. These are the registry's `lipid_ruv` and `metab_ruv` `fixed`, and the states `lipid_`/`metab_ruv_corrected`, `_correlation_filtered`, the correlations and the DA results. Every earlier state's message is identical (R1-01). |
| `test-workflow-omics-binding-steps.R` (deferred snapshot) | Counts unchanged (58/1). On top of the base's records, which hold KIT-214's `steps[[4]]` lines, the head adds only `steps[[9]]` parameters 9 to 10 (the marker), in lipidomics and metabolomics. |
| `test-gui-entry-transcriptomics.R` (VAL-236) | Counts unchanged (16/13). Only the line `actual$artifacts$integration/transcriptomics_normalized_data.tsv` moves, in `summary_save`, `summary_copy`, `summary_report` and `summary_export`. |

**Marker replay unchanged:** `test-kit202-intensity-filter-samples.R`, `test-kit205-enrichment-background.R`,
`test-kit214-sample-columns.R` and `test-prot-missing-value-counts.R` are `same`, passing as on the base.

**Every other file:** `same`.

### Predictions, KIT-251 (register section 6)

- `test-kit251-replicate-filter-one-click.R`: only on head, 9 passed, 0 failed.
- `test-perf107-audit.R`: `same`.
- The other ten files: `same`. Four of them run nothing: `test-prot-02h-…`, `test-module-ci-prot-protein-qc.R`,
  `test-workflow-gui-memory-repairs.R` and `test-workflow-memory-mechanism-member.R`.

### Predictions, R CMD check

- Findings identical on both sides, with nothing only on head.
- The portable-file-names NOTE lists the same number of paths on both sides.

### Decision and attempts

- If every prediction holds, the train is GREEN.
- A member whose prediction fails is named with its files and counts, and is not repaired here.
- One run, with one re-send allowed only for an instrument defect where no result was judged.

## Result

Run 45608483144 (request `merge-train-20261125f-e48254be`, run-tests commit 734367e9) bound head e48254be against base
63143aa3. It was sent once, concluded success, and used 494.8 job-minutes of the 1,340 ceiling (check 28.9). Tables are in
`merge-train-20261125f/`: `compare.tsv`, `tests.tsv`, `tests-base.tsv`, `check-compare.txt`, `check-summary.txt`,
`meta.txt`.

Of the 79 files, 65 are `same` and 14 are not. Every one of the 14 is predicted:

| File | Base p/f/e/s | Head p/f/e/s | Verdict | Prediction |
|---|---|---|---|---|
| test-kit224-ruv-feature-levels.R | absent | 69/0/0/0 | only-on-head | 69/0: held |
| test-lipid-12c-ruv-s4-characterization.R | 59/0/0/0 | 67/0/0/0 | count-changed | +8: held |
| test-metab-03a-ruviii-c-varying-characterization.R | 32/0/0/0 | 40/0/0/0 | count-changed | +8: held |
| test-kit127-omics-ruv-factors.R | 421/1/0/0 | 425/1/0/0 | count-changed | +4, its one failure identical: held |
| test-integration-input-transcriptomics.R | 78/0/0/0 | 79/0/0/0 | count-changed | +1: held |
| test-workflow-trx-declarations.R | 309/33/0/0 | 313/31/0/0 | message-changed | the two normalized_data.tsv leaves base-only: held |
| test-module-ci-trx-summary.R | 86/2/0/0 | 87/1/0/0 | message-changed | export digest base-only: held |
| test-prot-norm-core.R | 1548/306/0/0 | 1548/306/0/0 | message-changed | only the two RUV digests: held |
| test-omics-summary-core.R | 676/25/0/2 | 676/25/0/2 | message-changed | only the one leaf: held |
| test-omic-summary-module.R | 1675/28/0/0 | 1675/28/0/0 | message-changed | only the one leaf: held |
| test-workflow-omics-declarations.R | 45/23/0/0 | 45/23/0/0 | message-changed | KIT-224's lines only: held |
| test-workflow-omics-binding-steps.R | 58/1/0/0 | 58/1/0/0 | message-changed | `steps[[9]]` only: held |
| test-gui-entry-transcriptomics.R | 16/13/0/0 | 16/13/0/0 | message-changed | the four `summary_*` lines: held |
| test-kit251-replicate-filter-one-click.R | absent | 9/0/0/0 | only-on-head | 9/0: held |

### Message-level check

Failure records were compared base against head, with temporary paths masked.

| File | Base-only | Head-only | What moved |
|---|---|---|---|
| `test-prot-norm-core.R` | 2 | 2 | The two RUV result digests: 4ed0ad71… to 9ada5791… and 5065333a… to cc4276ca…. |
| `test-omics-summary-core.R` | 4 | 4 | Only the leaf `integration/transcriptomics_normalized_data.tsv` (trx_main and trx_managed). |
| `test-omic-summary-module.R` | 4 | 4 | Same leaf as above. |
| `test-workflow-trx-declarations.R` | the two `results/transcriptomics/transcriptomics_normalized_data.tsv` SHA-256 failures (actual be38d286…, stored a7045982…) | none | Nothing else. |
| `test-module-ci-trx-summary.R` | the MCI-032.1 export-file digest failure, with the same digests | none | Nothing else. |
| `test-gui-entry-transcriptomics.R` | 4 | 4 | Only the `integration/transcriptomics_normalized_data.tsv` line, base 4645738d… to head 7b1698df…, in the four summary steps. |
| `test-workflow-omics-declarations.R` | 9 | 9 | 14 records identical. See the breakdown below. |
| `test-workflow-omics-binding-steps.R` | (one record each side) | | The head adds only `steps[[9]]$parameters` length 10 against 9, and its `[[10]]` a list against absent, in lipidomics and metabolomics. KIT-214's `steps[[4]]` lines are unchanged on both sides. |

On the two summary files, the register's parenthesis needs a correction: on base the actual equals the stored 7b9d35db…, and
on head it is e0e8d896…. The register gave them the other way round. Its substance holds: only that one leaf moves.

On `test-workflow-omics-declarations.R`, the 14 identical records are every earlier state. The 9 that differ are exactly
KIT-224's:
- The registry, where the head adds only `lipid_ruv` and `metab_ruv` `declaration$fixed` as a list against NULL. KIT-214's
  `duplicate_resolution` lines are present on both sides.
- `lipid_ruv_corrected`, `lipid_correlation_filtered`, `metab_ruv_corrected` and `metab_correlation_filtered`.
- The two correlations records and the two DA-results records.

### Predicted `same`, observed `same`

| File | Passed/failed/skipped on both sides |
|---|---|
| `test-kit129-run-path.R` | 711/0/1 |
| `test-omic-norm-module.R` | 639/0 |
| `test-prot-missing-value-counts.R` | 67/0 |
| `test-kit202-intensity-filter-samples.R` | 163/0 |
| `test-kit205-enrichment-background.R` | 228/0 |
| `test-kit214-sample-columns.R` | 176/0 |
| `test-perf107-audit.R` | 191/0/1 |

The four KIT-251 files that run nothing skipped on both sides: 5, 6, 6 and 13 skips.

So every marker's replay behaviour is unchanged.

### R CMD check

- 0 ERRORs, 3 WARNINGs, 3 NOTEs on both sides.
- All 6 findings identical, none only on head.
- The portable-file-names NOTE lists 238 paths on each side.

### Verdict

Every prediction of KIT-224 and KIT-251 held, and no other file differs. TRAIN GREEN.

Staging was still at 63143aa3 when the result was read. The two rows `VAL-260` and `VAL-261` are numbered by
the coordinator at landing.

### Not checked

- The trees after the run's head (docs only).
- Files outside the 79.
- The members' full suites (KIT-224's R3 on the merged tree; KIT-251's A4 at the next checkpoint).

## Merge record (coordinator)

Landed on `v0.9-staging` on 2026-11-25 by a merge of this branch at 6cfb4823 into staging 63143aa3,
which had not moved since the train was built. The train's rows were numbered in train order:
VAL-260 (KIT-224), VAL-261 (KIT-251); the numbers were written into the ledger, both tickets'
`outcome.md` and this file. The landed tree differs from the train's tip only in those four files.

One code conflict was resolved inside the train on the coordinator's authority, as KIT-214's audit
(I10) had predicted: in `workflowReplayRuleMarkers` both new marker rows are kept,
`duplicate_resolution_samples` (KIT-214, already on staging) first and `ruv_feature_levels`
(KIT-224) after it, and the roxygen line names both tickets; no other line of that file changed.
The run included the three test files that read the marker list beyond the members' scopes
(KIT-202's, KIT-205's and KIT-214's), all unchanged.

Evidence: run 45608483144 (base 63143aa3, head e48254be; 79 test files and R CMD check, the
transcriptomics manifest flag set on both sides): 65 files the same, the 14 differences each as
predicted, down to the failure messages; the check identical on both sides (0 errors, 3 warnings,
3 notes; 238 paths in the portable-file-names note on each side). 494.8 of 1,340 job-minutes.

Stored references that move with KIT-224 and are not edited (they wait for the single
end-of-project re-record): the files named in the owner's decisions VAL-206, VAL-213 and VAL-236.

Not shown by this record: the full suite on the landed tree (next checkpoint).

## Correction note (coordinator, 2026-11-26)

Two items in this record's "Not checked" list are worded wrongly; the record above is left as
written. (1) "KIT-251's A4 at the next checkpoint": KIT-251's A4 is its records (NEWS, ledger,
outcome, hand-off) and is met; the full suite is not one of KIT-251's criteria. (2) "KIT-224's R3
on the merged tree": KIT-224's R3 (run 45573139074) ran at 2871125b against 67569959, not on the
merged tree; the full suite on the merged tree is the staging checkpoint.
