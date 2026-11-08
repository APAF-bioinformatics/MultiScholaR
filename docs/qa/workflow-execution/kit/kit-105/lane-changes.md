# KIT-105: entry-lane changes

Every re-frozen entry oracle in KIT-105, the defect or decision it is attributed to, and every capture that changed.

## F105-01: the lipid/trx normalisation grid's image slots and assay labels

- **Defect** (found by KIT-105 W1, `w1/notes-for-phase-b.txt`, `w1/test-first.txt`; display only). The registrars
  `registerLipidNormStaticQcImageOutputs` and `registerLipidNormAssayLabelOutputs`
  (`R/mod_lipid_norm_support_helpers.R`) assign renderers in a `for` loop. The renderer builders'
  returned functions (`buildLipidNormQcImageRenderer`: `function(assaySlot, plotType, stagePrefix)`;
  `buildLipidNormAssayLabelRenderer`: `function(assaySlot)`) left their arguments as lazy promises, and those were first
  read when the outputs rendered, after the loop had ended. So all 24 image slots bound the last spec (assay 2,
  `ruv_corrected`, `correlation`), and all 8 labels bound assay 2. Trx mounts the same registrars
  (`R/mod_trx_norm.R:217-227`), so its slots showed nothing and its labels read "Assay 2: (detecting...)". No result,
  state, file or log changed. Metab's hand-written bindings were already correct.
- **Fix** (`e76e1e14`, lead-assigned under the standing authority of 2026-11-07: a defect fix that does not change the
  scientific method): `force()` the arguments at the top of both returned functions, so each output keeps its own spec.
- **Test first** (`f105-01/test-first.txt`): the W1 block "F105-01: each lipid QC image slot and assay label shows its
  own assay, stage and kind" failed 13 times before the fix (one per captured lipid step) and passes after it. Its
  `omicNormPending()` gate is removed. In `test-omic-norm-module.R`, only the lipid `images` and `labels` digests of
  the 13 (b) steps were re-frozen (26 lines). Every other frozen value is unchanged. The re-freeze script checked
  this against the dump before writing.

### Runs (R 4.6.1, `pkgload::load_all` via `kit/scripts-078/run1.R`, one R process at a time, `MULTISCHOLAR_STRICT_SKIPS=1`)

| Run | File | Pass | Fail | Error | Skip | Elapsed |
|---|---|---|---|---|---|---|
| Verify against the old oracles (fix applied) | test-gui-entry-lipidomics.R | 15 | 20 | 0 | 0 | 74 s |
| Verify against the old oracles (fix applied) | test-gui-entry-transcriptomics.R | 13 | 16 | 0 | 0 | 46 s |
| Re-freeze (`MULTISCHOLAR_GUI_ENTRY_FREEZE=1`), once | test-gui-entry-lipidomics.R | 15 | 0 | 0 | 0 | 71 s |
| Re-freeze (`MULTISCHOLAR_GUI_ENTRY_FREEZE=1`), once | test-gui-entry-transcriptomics.R | 13 | 0 | 0 | 0 | 41 s |
| Verify against the new oracles | test-gui-entry-lipidomics.R | 35 | 0 | 0 | 0 | 73 s |
| Verify against the new oracles | test-gui-entry-transcriptomics.R | 29 | 0 | 0 | 0 | 46 s |

Each file holds only its own lane's `test_that` block, so each freeze run wrote only `lipid_canonical.json` or
`trx_screen.json`. sha256 of every `tests/testdata/gui-entry/*.json` before and after: only those two differ. The other
eight oracles (`metab_{lc,gc,combined}`, `prot_{dia,dia_limpa,tmt,lfq,lfq_fragpipe}`) are byte-identical. The 20 and 16
verify failures are one step comparison per differing step. Each failure is fully covered by the path list below.

### Every changed capture

Path-by-path walk of the old and new oracles (full list: `f105-01/oracle-diff.txt`). Both headers changed
`source_sha` `ddf42372` -> `e76e1e14` and `frozen` 2026-11-07 -> 2026-11-08. All other differences are outputs of the
normalisation tab's grid (`<omic>_workflow-norm-*`). Nothing else changed: heads, histories, `tab_status`, state
digests, artifacts, notifications, modals and all other outputs match the old oracles.

**lipid_canonical** (assays LCMS_Pos, LCMS_Neg; 333 captures over 20 of the lane's steps):

| Capture | Steps | Before | After |
|---|---|---|---|
| `norm-assay1_label_{pca,density,rle,correlation}` (4) | import_refusal, import, design_refusal, design, qc_intensity, qc_duplicates, qc_itsd_refusal, qc_finalize, da_refusal (9) | "Assay 2: (detecting...)" | "Assay 1: (detecting...)" |
| `norm-assay1_label_{pca,density,rle,correlation}` (4) | norm_refusal, norm, correlation, session_export, da_load, da, summary_refusal, summary_save, summary_copy, summary_report, summary_export (11) | "Assay: LCMS_Neg" | "Assay: LCMS_Pos" |
| `norm-{pca,density,rle,correlation}_post_filter_assay{1,2}` (8) | the 11 steps from norm_refusal | the last spec's no-file image (`<table> 10d8284a...`) | the assay's own pre-normalisation figure (`<plot> alt=<kind> - LCMS_Pos` / `LCMS_Neg`) |
| `norm-{pca,density,rle,correlation}_post_norm_assay{1,2}` (8) | norm_refusal | `<table> 10d8284a...` | each slot's own no-file image (a distinct `<table>` digest per slot) |
| `norm-{pca,density,rle,correlation}_post_norm_assay{1,2}` (8) | the 10 steps from norm | `<table> 10d8284a...` | the assay's own post-normalisation figure (`<plot> alt=<kind> - <assay>`) |
| `norm-{pca,density,rle}_ruv_corrected_assay{1,2}`, `norm-correlation_ruv_corrected_assay1` (7) | the 11 steps from norm_refusal | `<table> 10d8284a...` | each slot's own no-file image (distinct `<table>` digests; no RUV-corrected figure exists in this lane) |

The following captures did not change: the 4 assay-2 labels, `norm-correlation_ruv_corrected_assay2` (it was the
last spec, so it was already right), and all image slots in the 9 steps before normalisation. In those steps no
assay is detected, so every slot shows the same "Assay not detected yet" value whichever spec it holds.

**trx_screen** (one assay, `counts`; stage files `raw`/`tmm`/`ruvs`; 256 captures over all 16 of the lane's steps):

| Capture | Steps | Before | After |
|---|---|---|---|
| `norm-assay1_label_{pca,density,rle,correlation}` (4) | all 16 | "Assay 2: (detecting...)" | "Assay: counts" |
| `norm-{pca,rle}_{post_filter,post_norm,ruv_corrected}_assay1` (6) | import_refusal ... norm_refusal (8) | assay 2's "Assay not detected yet" (`<table> ef4895c1...`) | each slot's own no-file image (distinct `<table>` digests) |
| `norm-{pca,rle}_{post_filter,post_norm,ruv_corrected}_assay1` (6) | norm, da, summary_refusal, summary_ruviii, summary_save, summary_copy, summary_report, summary_export (8) | `<table> ef4895c1...` | the stage's own figure (`<plot> alt=pca - counts`, `<plot> alt=rle - counts`) |
| `norm-{density,correlation}_{post_filter,post_norm,ruv_corrected}_assay1` (6) | all 16 | `<table> ef4895c1...` | each slot's own no-file image (distinct `<table>` digests; no such file in any step) |

The 4 assay-2 labels and the 12 assay-2 image slots did not change. Trx has one assay, so these still read
"Assay 2: (detecting...)" / "Assay not detected yet".

### Attribution

All 589 changed captures (333 lipid, 256 trx) and both header changes come from F105-01 alone (`e76e1e14`). They are
the grid outputs that the defect had bound to the wrong spec. No other lane was re-frozen.

### Regressions run after the fix (unchanged expected)

| File | Pass | Fail | Error | Skip | Elapsed |
|---|---|---|---|---|---|
| test-omic-norm-module.R (gate removed, lipid digests re-frozen) | 482 | 0 | 0 | 0 | 120 s |
| test-module-ci-lipid-norm.R | 153 | 0 | 0 | 0 | 47 s |
| test-module-ci-trx-norm.R | 61 | 0 | 0 | 0 | 74 s |
| test-module-ci-metab-norm.R | 150 | 0 | 0 | 0 | 35 s |
| test-workflow-omics-application-reads.R | 230 | 0 | 0 | 0 | 420 s |
| test-workflow-omics-worker.R | 283 | 0 | 0 | 0 | 360 s |
| test-gui-entry-metabolomics.R (3 lanes, identical to their oracles) | 105 | 0 | 0 | 0 | 200 s |

Records: `f105-01/test-first.txt`, `f105-01/targeted.txt`, `f105-01/oracle-diff.txt`.
