# Merge train 8 of 2026-11-26: KIT-244, KIT-243, KIT-201

Built on the coordinator's instruction after train 7 landed (staging 42846c0e, ledger rows through VAL-263). The
procedure is the same as for train 7 (`merge-train-20261125g.md`).

## Train

Base: `42846c0e521a346d7755aad9478f3a093961d117` (`origin/v0.9-staging`, fetched before building).

| Order | Ticket | Branch head | Merge commit | Ledger rows | Audit |
|---|---|---|---|---|---|
| 1 | KIT-244 | `origin/kit-244` df6986ac | 092302c0 | `VAL-<train8-1>` | AUDIT-KIT-244-impl-20261125 (section 4 (7), I04) |
| 2 | KIT-243 | `origin/kit-243` 43a3280e | 80175f90 | `VAL-<train8-2>`; stored-value row `VAL-<train8-2b>` (cited, left for the coordinator to write) | AUDIT-KIT-243-impl-20261125, closure check C01-C03 |
| 3 | KIT-201 | `origin/kit-201-treat-follows-record` 84bda25e | aceb7190 | `VAL-<train8-3>`, `VAL-<train8-4>` | AUDIT-KIT-201-impl-20261125, closure check C04 |

The run binds head `aceb7190`; the full hash is in the request. This file and the records commit after the run are
under `docs/qa/` only. All rows are numbered by the coordinator at landing.

KIT-244's merge commit message summarises the change imprecisely. Its NEWS entry is the accurate statement: after the
correlation filter on DIA with limpa, pressing 'Run Normalization & RUV' again keeps the table the first press started
from, and a press before DPC Quantification is refused. The commit message was left as written, with no amend.

### Conflicts met and their resolutions

**`docs/methods-decision-ledger.md`, all three merges.** Staging's rows through VAL-263 were kept, then each member's
rows in train order.
- KIT-244's row carries its I04 text (already on the branch) and a clean status.
- KIT-243's row:
  - I09's correction-run placeholder `<run id>` is filled with 45611244476;
  - its `VAL-<m>` is renamed `VAL-<train8-2b>`;
  - "the train's run" is left to fill after the run.
- KIT-201's two rows: the second (DAS-016) keeps "decided" and adds "merged … with KIT-201".
- Each member's `VAL-<merge>` in its `outcome.md` was replaced; for KIT-201, also in `wording.md`.

**`NEWS.md`, merges 1 and 3.**
- Merge 1: KIT-244's Fixes entry against train 7's entries.
- Merge 3: KIT-201's "Scientific changes" entry against KIT-224's.
- Both entries were kept each time, staging's first, under the one unreleased heading. NEWS is clean in merge 2.

**`R/utils_workflow_analysis_record.R`, merge 3.** This is the authorised hand resolution (KIT-201 closure check C04
item 2; the coordinator's instruction).
- The marker list keeps every row in staging's order (KIT-214's, KIT-224's), then KIT-201's row and its three helper
  functions, unchanged.
- The header names KIT-214, KIT-224 and KIT-201.
- The DA replay block came from the branch already resolved.
- No other line changed. The block is `@noRd`.

Resolved hunk against the base:

```diff
-#' adds `intensity_filter_samples`; KIT-214 adds `duplicate_resolution_samples`; KIT-224 adds `ruv_feature_levels`): for
+#' adds `intensity_filter_samples`; KIT-214 adds `duplicate_resolution_samples`; KIT-224 adds `ruv_feature_levels`;
+#' KIT-201 adds `treat_moderation`, which no declaration fixes: it is read from the DA model record by the DA's
+#' replay, `workflowDaReplayRuleOptions`): for
@@
-    ruv_feature_levels = list(option = "MultiScholaR.ruv_feature_levels", current = "restored", earlier = "shifted"))
+    ruv_feature_levels = list(option = "MultiScholaR.ruv_feature_levels", current = "restored", earlier = "shifted"),
+    # KIT-201 (VAL-131, D3): ... (three comment lines, unchanged)
+    treat_moderation = list(option = "MultiScholaR.treat_moderation_rule", current = "as_main_test", earlier = "constant"))
+
+ workflowTreatModerationCurrent(), workflowTreatModerationRecord(), workflowDaReplayRuleOptions() (unchanged)
```

**No other conflict.** `DESCRIPTION` and `NAMESPACE` are unchanged; the Version line is unchanged. `man/` changes only
by KIT-201's own `man/runTestsContrasts.Rd` (+4 lines, from its branch).

### Owed at landing, outside this branch

- VAL-246's status cell (KIT-244 I04).
- VAL-212's status cell and the stored-value row `VAL-<train8-2b>` (KIT-243 I03, C03).
- `OPEN-DEFECTS.md` O-11 (KIT-201 I05).

### Path lengths

Checked by script before sending: 76 paths are new since the base. The longest is 76 bytes with the `MultiScholaR/`
prefix. None exceeds 100.

## Declaration (written before the run is sent)

**Purpose.** Each member's check on the merged tree, done once on the train. For KIT-243 this is the only comparison
with staging since KIT-220 (C02).

**Scope.** 177 files, deduplicated, base and head, with R CMD check on both sides.
- KIT-244: the 59 files of `kit244-r2-e4a37f2e.json`.
- KIT-243: the 78 files of its correction run, `kit243-corr-47bf368a.json`.
- KIT-201: the 108 files of R2b, `kit201-r2b-50a5023d.json`, plus `test-omic-norm-ruv-design.R`,
  `test-kit214-sample-columns.R`, `test-kit224-ruv-feature-levels.R` and `test-kit251-replicate-filter-one-click.R`.

**Run settings.**
- Environment on both sides: `MULTISCHOLAR_KIT105_W3_DUMP=1` and `MULTISCHOLAR_TRX_MANIFEST_NOTE=1` (KIT-244 section
  4 (7); KIT-243 C01).
- `sva` installed, as KIT-201's R2b.
- 8 shards per side.
- Per-file kill limits of 120 minutes for `test-prot-norm-core.R` and the four GUI-entry files. These are not
  thresholds.

**Ceiling.** 2,280 job-minutes: the sum of the three members' R2 ceilings, which are KIT-244's 660, KIT-243's 1,080 and
KIT-201's 540 (6 × 60 + 2 × 90), each from its ticket's section 10.

**Base values.** They come from the base as it now is: the latest train run on an identical tree for that file, or the
members' runs on staging bases where no train ran it. Where a member's record predates trains 6 and 7, the prediction
states the change.

### Predictions, KIT-244 (register section 4 (7))

- `test-kit244-limpa-pass-through.R` is only on head: 204 passed, 0 failed, 0 skipped. Any skip is an instrument
  defect.
- `test-prot-norm-core.R` is `same`: 1548/306 on both sides, at message level. The body-digest block passes on each
  side against its own pin: base fe852d7a…, head e991f0d2… (VAL-246, re-pin 8fa690d0).
- Its other 57 files are `same`, except the KIT-201 moves listed below.

### Predictions, KIT-243 (closure check C.5)

- `test-kit243-summary-current-outputs.R` is only on head: 246 passed, 0 failed.
- `test-workflow-saved-da.R` is covered by the combined prediction below.
- Every other file is `same` at staging's own counts with the switch, apart from the KIT-201 moves below. For
  example, `test-kit220-da-threshold-consumers.R` is 79/0 on both sides; R2b's 71/2 was the missing switch.

### Predictions, KIT-201 (closure check C04; R2b 45607796537, leaf for leaf)

| File | Base (passed/failed/skipped) | Head (passed/failed/skipped) | Verdict | Moved leaves |
|---|---|---|---|---|
| `test-kit201-treat-moderation.R` | absent | 89/0 | only on head | — |
| `test-prot-ruv-factors-in-model.R` | 306/1 | 312/0 | `fixed` | — |
| `test-gui-entry-proteomics.R` | 189/4 | 168/25 | new-failure | exactly R2b's head-only leaves, including the 12 `prot_lfq_fragpipe` NumSig leaves |
| `test-omic-da-module.R` | 1316/53 | 1289/80 | new-failure | as R2b |
| `test-omic-enrich-module.R` | 1051/19 | 1011/59 | new-failure | as R2b |
| `test-omics-enrich-core.R` | 232/9/1 | 221/20/1 | new-failure | as R2b |

- If I13's flaky leaf fails on base, `test-omic-enrich-module.R` reads 22 to 59. The named flaky leaves are I13's
  three, the "refusals table_removed/{da,open,run_table_removed} values", and plan section 11's `gprofiler_plot`.
- `test-workflow-saved-da.R` is covered by the combined prediction below.
- Every other file is `same`, including every lipid/metab, limpa and transcriptomics file.
- The stop rule is any other moved key. In particular: a NumSig leaf outside `prot_lfq_fragpipe.json`, a lipid/metab
  statistic at cutoff 0, or any limpa or transcriptomics value.

### Combined prediction for `test-workflow-saved-da.R` (KIT-243 with KIT-201)

The prediction derives from both registers and both tickets' run records:
- KIT-243 R2 `r2-saved-da-files.txt`;
- KIT-243 closure C.3 and C.5;
- KIT-201 `r2b/predictions.md` and R2b's `compare-failures.txt`.

**Counts.** 780/171 on the base (staging; R2b's base at 63143aa3, a tree equal in this file's code paths). On head,
776/175, `message-changed`. The four new failures are KIT-201's: `dia_non_head record` and `dia_reopened_rerun record`,
each in W1 (a) and W2. KIT-243 changes no count: equal counts in its R2 and in the correction run.

**Messages.** The head's differences are confined to the union of the two tickets' sets.
- KIT-201's set: in W1 (a) and W2 of `savedDaW1Frozen()`, the already-failing `dependencies`, `display`, `files`,
  `head_manifest`, `manifest`, `products` and `record` of dia_changed_arguments, dia_non_head, dia_one_contrast,
  dia_qvalue_only, dia_refused, dia_reopened_rerun, dia_repeat, dia_two_contrasts and dia_retry, plus the two new
  `record` failures.
- KIT-243's set: the `files` messages of dia_changed_arguments, dia_qvalue_only and dia_refused in W1 (a) and W2.
- KIT-243's set lies inside KIT-201's. So the union equals KIT-201's set.

**What cannot be derived.** The head's actual `files` digests for KIT-243's three cases cannot be derived from the
records, because no run had both changes. They are predicted to differ from R2b's head values for those three cases.
All other moved values are predicted to equal R2b's head values.

**Unchanged.** `fits`, `states` and every limpa case.

### Predictions, R CMD check

- Findings identical on both sides; nothing only on head.
- The portable-file-names NOTE lists the same number of paths on both sides (238 at train 7).

### Decision and attempts

- Every prediction holds: TRAIN GREEN.
- A member whose prediction fails is named with its files and counts. It is not repaired here.
- One run. One re-send is allowed only for an instrument defect where no result was judged.

## Result

(filled after the run)
