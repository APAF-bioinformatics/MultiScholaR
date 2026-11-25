# Merge train 6 of 2026-11-25: KIT-224, KIT-251

This follows the coordinator's instruction after train 5 landed (staging 63143aa3). The procedure is the same as for train 5
(`merge-train-20261125e.md`).

## Train

Base: `63143aa3618eae2729c6864ab10903928bec5071` (`origin/v0.9-staging`, fetched before building).

| Order | Ticket | Branch head | Merge commit | Ledger row | Audit |
|---|---|---|---|---|---|
| 1 | KIT-224 | `origin/kit-224` c2fc484c | cb1e70eb | `VAL-<train6-1>` | AUDIT-KIT-224-impl-20261125 (sections 5, 6) |
| 2 | KIT-251 | `origin/kit-251` d7a66e3e | e48254be | `VAL-<train6-2>` | AUDIT-KIT-251-impl-20261125 (section 6) |

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
  stale copies were dropped. KIT-224's row was then added as `VAL-<train6-1>`.
- In merge 2, KIT-251's row was added as `VAL-<train6-2>`, with the status form of its register's I04.
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

(filled after the run)
