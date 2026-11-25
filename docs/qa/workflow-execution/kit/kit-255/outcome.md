# KIT-255 outcome (ticket r2.2; owner's decision VAL-283)

Branch `kit-255b`, branch point staging 988512e6 (KIT-254 on it). Head 146f44a5 at the first hand-back (records only; fix c25e1335; tests commits 03e6dc39, 35eb63b5 and c46afc90); the correction pass is the last section. Scripts-only branch `kit-255-g1` (never merged). No R was run on the laptop; every computation ran on GitHub runners (runs.tsv).

## What was done

- Step 0: ticket text r2.2 (RC2-01 to RC2-05 applied; RC2-06 recorded below, not fixed).
- G1 (observation, base): run 45691779984 (instrument defects in 45685635354 and 45687997124, `runs.tsv`). Observations: `g1-45691779984/observations.tsv`, predictions in `tools/kit255/g1.R` and `g1-*/predictions.tsv`.
- R1 (tests only, 03e6dc39): run 45695453629.
- The fix (c25e1335) and R2: run 45699606601 (50 files, 4 shards per side, check); R2b 45710190574 (the new file alone after one test-mechanics repair).

## G1: observed (base, staging 988512e6 + script)

Held as predicted (RUV applied, with the RUV objects' correlation pairs hand-made: see Deviations): after the press the design lacks B_6, LCMS_Pos holds it, LCMS_Neg lacks it; `args$ruv_uncorrected_assays` keeps B_6 finite and identical before and after (P-a), also on the tab in process and on the worker (managed design = in-process design); lipid DA refuses on both routes ('Assay LCMS_Pos has sample column(s) missing from design: B_6'), metabolomics DA completes (P-b); the integration input has 11 samples, B_6 as a stray feature-table column of the LCMS_Pos view, uncorrected matrices without B_6 (P-c); hand-built (b) state, non-RUV (P-d): kernel and declaration DA complete, LCMS_Neg equal to the without-the-sample state to 1e-10 in every common numeric column, residual df feature by feature equal, tab preflight returns NULL, LCMS_Pos differs, 5 features fall from three or more observed in group B to fewer; RUV-factor kernel on the pipeline's own factors (P-e): both complete, LCMS_Neg equal to 1e-10, W has a B_6 row; K3-b (group B blank in LCMS_Neg, P-f): the table has 520 rows and no finite statistic (limma NA coefficients), LCMS_Pos tested; figures (P-g): PCA and density omit B_6 in LCMS_Neg (11) and show it in LCMS_Pos (12), RLE builds with 12 in both (not predicted), static and Glimma volcano build, the LCMS_Neg heatmap draws B_6 as 50 of 50 zeros, the Combined heatmap draws B_6's LCMS_Neg cells as zeros (27 of 50); integration (P-h): the writer gives 12 samples, 12-column matrices, B_6 numeric all NA in LCMS_Neg, no stray column; the model (lipidomics and metabolomics inputs together) counts B_6 in every view, `drop_incomplete` keeps B_6, whole and partitioned alignments identical (it warns of few samples only); P-i by arithmetic on the metadata reader (12 design rows, 12 per assay, Removed 0); P-k: the LCMS_Neg view's removed factors `not_determined`, 0 features, warning class `integration_input_removed_factors_undetermined`, input written without the factors, LCMS_Pos `exact` (new design and legacy RUV-III-C, both omics); P-j: the base's K2 one-assay-failure project replays exactly on the base (0 mismatches); tarball `tests/testdata/kit-255/base-record.tar.gz` (sha256 818b5980...; written by the base at f4d9c5e7).

Not observed (instrument, corrected by the audit): P-f's RUV arm and P-d's per-sample columns were reported as observations but the script read the wrong names. The script read `ruv_not_observed` and `ruv_not_tested` (`tools/kit255/g1.R:841-842` at 708e5f6a); the kernel returns `ruv_features_not_observed` and `ruv_features_not_tested` (`R/func_lipid_da_model.R:646-648`). P-d looked for a column named `B_6` (`g1.R:309-310`); the long tables carry `intensity.B_6.B` (shown in P-f's column list), so ticket section 6 row 14 is unobserved. No rerun: P-f decides wording only. P-g: RLE keeps B_6 in LCMS_Neg (not predicted).

**Finding (science, for the owner's information; nothing decided):** with real correlations after RUV-III (new design or legacy) a one-sample anomaly is removed by RUV: B_6 (LCMS_Neg values drawn at random) had a best within-group correlation of 0.99 in BOTH assays after RUV (run 45687997124), so the filter, which runs after RUV, failed nothing. The filter fails a sample after RUV only when the anomaly is not removable by k factors.

## R1 (tests only): 45695453629

Per row (`r1/tests.tsv`, `r1/failures.txt`): E1, E2, E3, E5, E6, E7 to E10, E12, E13, E15, E17, E18, E19, E20, E22, E25, G6 failed at their base values by assertion; E4 failed only on its comparison (the test's stand-in object was refused by the class validity, so the stated base value was not reached: an instrument defect); E14 was not reached (E13's marker assertion errored first); E16 had no test (lost when the file was rewritten); E23 lipid was skipped (the base DA refuses) and E23 metab failed; E24 was a test error (matrix width); T-OLD and the G rows that pass on base passed (`r1-predictions.tsv`). E23's lipid defect is shown by G1 P-g and E24's base text by reading of `R/utils_workflow_integration_input.R:141-142` at 988512e6. Test mechanics repaired afterwards, no expected value changed: G4's fixture shape (kit172Values with `n_itsd = 0`, which fails, replaced by `n = 45, n_itsd = 5`), E24's matrix width, E13's marker access, G6 split into its own test.

## R2: 45699606601 (base 988512e6 against head c25e1335)

45 of 50 files `same`; check identical (0 ERROR, 3 WARNING, 3 NOTE both sides); nothing only-on-head except the new file. Differences, all declared:
- `test-metab-s4-qc-shared.R`: count-changed (32 to 35), the declared edit.
- `test-omic-norm-module.R`: new-failure, **the eight provenance digests, NOT edited (VAL-303: they fail by design until the end-of-project re-record)**: lipid `:671`, `:695`, `:712`, `:743` and metab `:896`, `:920`, `:937`, `:968` (the skip scenario's `correlation`, `final_plot`, `export`, `reset` provenance; the run's messages: actual 15e689a4... against expected 204adbf5... for lipid correlation, final_plot and export, c0ed83bb... against baf23fbe... for lipid reset; 9b3649ec... against d3c0ca75... for metab correlation, final_plot, export, 62d62f5a... against 4414a8e8... for metab reset). Predicted before the run.
- `test-workflow-omics-declarations.R`, `test-workflow-omics-binding-steps.R`: message-changed, the KIT-073 and KIT-094 dumps (VAL-181, VAL-189: end-of-project re-record; they fail the same count on both sides, the dump value differs).
- `test-workflow-lipid-declarations.R` and `-metab-` : same; their `fixed` pins (`:47-49`) cover ITSD and the DA only, and this ticket changes the selection step's `fixed` list only, so (VAL-303) they are not re-pinned.
- New file: 288 passed, 2 failed (E4, lipid and metab): the test's without-the-sample object was refused by the class validity (an instrument defect in the test; no result). Repaired; R2b: 294 passed, 0 failed, 0 skipped.
- `test-kit254-omics-correlation-skip.R`: same (its declared E9 and E17 edits pass).
- Largest shard 66.9 minutes (4 shards kept).

## Acceptance (ticket 10.5)

| # | Criterion | Evidence |
|---|---|---|
| A1 | G1 before any code | 45691779984 (and records) |
| A2 | R1: the E rows fail at their base values, G rows and T-OLD pass | 45695453629 (E4's comparison, E14, E16, E23 lipid and E24 did not show their base value: see R1) |
| A3 | R2 base against head, 50 files + check | 45699606601: every file `same` except the new file and the declared moved keys; the eight digests and the two dumps are the declared moves (VAL-303, VAL-181, VAL-189) |
| A4 | the merge train's full suite | not run here |
| A5 | NEWS, ledger row, outcome.md with "For the change report", `man/` (two Rd), NAMESPACE and DESCRIPTION unchanged | in the branch |

## Deviations (every one)

1. The RUV objects' correlation pairs in G1 and in E7 to E10 are hand-made (stop rule (i)): see the finding above. Everything else in those objects is real.
2. G1 P-e and the RUV arm of P-f call the exported DA kernel with the pipeline's own factors and normalised state, with the failing column blank, not the declared RUV DA operation on a blanked state (the declared operation is exercised in P-b on the real state).
3. G1 P-h ran the model on the lipidomics and metabolomics inputs together (it needs two omics); P-i is the arithmetic of the summary on the metadata reader, not the rendered tab.
4. G1 instrument reruns: 2 of 3 used (45687997124, 45691779984); one request (the -3 file) refused at plan (requested_by over 300 characters) and re-sent as -4; one attempt to chain request writing and commit was blocked and redone in separate commands (no rule broken in effect).
5. (Withdrawn in the correction pass.) The notice had no sentence for a sample removed because no assay kept a value for it, on the ground that the kernel's message names the removal; that message never reaches the tab (AUDIT-KIT-255-impl I03). The sentence is now added.
6. `correlation_removed_from_design` is added to the record only when a sample was set missing or removed (so an unchanged press records what it recorded before); the design difference is read from the step's `previous_state`.
7. The exclusion record keeps one scalar `threshold` while every entry shares it, else a list by assay of named thresholds (the ticket left the form open).
8. The shared kernel helper lives in `R/utils_workflow_omics_declarations.R` (a ticket file) beside `omicsPearsonSelectionOperation`, called by both kernels.
9. Four GUI-entry timeouts: only the two GUI-entry files in scope got 120 minutes.
10. The `fixed` pins of the declarations tests were not re-pinned (VAL-303's condition not met).

## Stated limits (in NEWS and the ledger row after the correction pass)

- An RUV-corrected view holding a sample set missing in it has no integration "removed factors" (`not_determined`; reason names the run); the opt-in adjustment is unavailable there.
- The integration model's presence rule still counts such a run as measured (INT-010; a candidate row, not this ticket).
- In the DA heatmaps' Combined view a sample blank in one assay only keeps its column, and its cells there are drawn as zeros (RC2-03): the repair drops columns of the combined matrix with no finite value.
- The per-group minimum is not re-applied after the blanking (parity with every removal today).
- Interim window until KIT-242 Part A merges: where a sample's values were missing altogether in one assay, that assay may not have been normalised between samples; NEWS carries the conditional sentence; no counted run or real-data acceptance run may bind a commit inside the window on "Skip RUV" with ITSD, or its evidence states it.

## RC2-06 (not fixed here; for a separate ticket)

The Combined DA heatmap binds the assays' matrices by position: each assay's matrix takes that assay's own column order (`R/func_lipid_da_plotting.R:427`, `R/func_metab_da_heatmap.R:131` at 988512e6) and `do.call(rbind, ...)` (`:437`; metab `:141`) keeps the first matrix's column names. If two assays list their samples in different orders, one assay's values are drawn under another sample's name and group. A display defect, not caused by this ticket; the KIT-255 repair (dropping columns of the combined matrix with no finite value) does not rely on positional agreement. Smallest repair: order each assay's matrix by the design's sample order before the bind.

## Not checked

At the first hand-back: E16 (not written), E8's equality on the declared route, the managed notice, the two-threshold merge, T-OLD's option-built lanes (E16, E8, the managed notice and the merge are written in the correction pass; the option lanes are not).


The managed worker with a view where RUV applied and the filter failed a sample through real (not hand-made) correlations; Windows and the Electron build; real data (how often a one-assay failure occurs); MOFA2 itself on a blank view (the model is run to its alignment); the report templates rendered with an exclusion (the facts and sentence are by reading and the template chunks; KIT-235's parser read of the new study-parameters block by shape only); the session export's new field and the replay statement's plural forms (tested through KIT-254's E17 and the singular case only); the full suite (the merge train's).

## For the change report

*Fix-commit lines (c25e1335):* the kernel helper `R/utils_workflow_omics_declarations.R:764-857`, the lipid kernel `R/func_lipid_s4_correlation_methods.R:367-371`, the metabolomics kernel `R/func_metab_s4_qc_methods.R:432-436`, the writer `R/utils_workflow_omics_declarations.R:2304`, `:2321` (line numbers at c25e1335, before the correction pass). Instrument use: G1 reruns 2 of 3 and R2b (VAL-54).

**What was wrong.** With two or more assays, the lipidomics and metabolomics sample-correlation filters removed a failing sample's column only from the assay where it failed and the sample from the design if it failed in any assay (`main:R/func_lipid_s4_objects.R:5024`, `:5037`; `main:R/func_metab_s4_objects.R:5031`, `:5044`; at 988512e6 `R/func_lipid_s4_correlation_methods.R:403`, `:416`, metab `R/func_metab_s4_qc_methods.R:468`, `:481`). On `main` both differential analyses then took the design's samples by name in every assay (`main:R/func_lipid_da.R:433`, `main:R/func_metab_da.R:433`), so a sample failing in one assay was analysed in none, silently, and its column stayed in saved objects. New to 0.9 and not on `main`: the lipidomics refusal on every route (0ae6059e, 2026-05-06; for lipidomics a regression from `main`'s result), the metabolomics refusal on the tab in process and completion on a worker, the stray column in the integration input's feature table, and the RUV design and integration input that hold further copies of the values. Affected: lipidomics and metabolomics users with two or more assays (positive and negative mode) in which a sample fails the correlation filter in one assay only. How to recognise affected outputs: a saved state, exported session or integration RDS whose assays hold different sample columns; a `features_<view>.parquet` with a column named after a run. Recompute: apply the filter again (or skip) after the fix, then the differential analysis and Save.

**What changed (code, fix commit c25e1335).** The shared kernel helper `omicsCorrelationFilterWithinAssay` (`R/utils_workflow_omics_declarations.R`), called by both kernels (`R/func_lipid_s4_correlation_methods.R`, `R/func_metab_s4_qc_methods.R`) unless the replay marker holds its earlier token (then the base body runs); the marker `correlation_exclusion` (`R/utils_workflow_analysis_record.R`, `workflowReplayRuleMarkers`) fixed on `omicsPearsonSelectionOperation`; the writer blanks the normalised matrix of an RUV-corrected view (`omicsIntegrationInputSource`); the record outcomes (`workflowRecordCorrelationOutcomes`); the tab's notice and summary line (`R/mod_omic_norm.R`, `R/mod_omic_norm_panels.R`); the replay statement; the study parameters block; the session export field; the report clause; the heatmaps; the removed-factors reason. In results: a sample that fails in one assay stays in the design and in the assays where it passed, set missing where it failed; results in the passing assays include it (standard error of a two-group contrast smaller by 0.89, 0.93, 0.95 at 3, 4, 6 per group; inferred, equal-variance arithmetic) and differ from `main` and from the development builds by that replicate; where the cause of the poor correlation is the specimen itself, its other assays, which passed the same test as every other sample kept, stay in (the owner's stated consequence). A sample empty in every assay with no replicate pair, which the base kept with all-NA columns, is removed (VAL-242 1a). Older records replay as computed (T-OLD, the base-written project, replays exactly on head).

**Evidence.** Runs 45606617982 and 45609166200 (G0, the defect), 45691779984 (G1), 45695453629 (R1), 45699606601 (R2), 45710190574 (R2b); simulation settings: K2 from `kit172Values` seeds 255001 and 255003, replacement seed 255002, 12 samples (two groups of 6), 520 features; limits: synthetic data; the RUV objects' failing pairs are hand-made.

## Hand-offs (as AUDIT-KIT-255-impl-20261126 words them; nothing here is acted on by KIT-255)

- **I08 (for the owner; decided, VAL-330):** after RUV the sample-correlation filter, which judges the RUV-corrected values, will rarely fail a single grossly deviant sample: RUV-III absorbs a one-sample deviation that dominates the within-group residuals and shows in the control features (0.991 after RUV in the anomalous assay against 0.985 in the genuine one). On the new design this is benign (the factor in the model in effect leaves the sample out of that assay's test, at one residual degree of freedom); on the legacy route the sample becomes a pseudo-replicate (standard errors too small by derived factors 0.775, 0.845, 0.905 at 3, 4, 6 per group; derived, not measured). The owner decided (VAL-330): STATED now (NEWS, wording.md) and MEASURED on runners before the release by its own ticket (KIT-275); the method order is not changed.
- **I10 (candidate):** the Combined DA heatmap sets every missing cell to 0 after row scaling (`R/func_lipid_da_plotting.R:462`, `R/func_metab_da_heatmap.R:166` at 146f44a5), so a sample set missing in one assay is drawn and clustered at the feature's average in that block. Candidate: cluster on the filled matrix, draw missing cells in a neutral colour (`na_col`).
- **I11 (candidate):** `excluded[[assay_name]] <- ...` with an unnamed assay list (`R/utils_workflow_omics_declarations.R:808` at 146f44a5) stops with "[[ ]] with missing subscript" when a sample fails, where the base worked; an unnamed pairs list of length 1 is applied to the first assay only (`:787`), where `purrr::map2` recycled it. Candidate: key by position when unnamed.
- **I12 (owner, information; coordinator, schedule):** the integration model's presence rule counts a blank run as measured and `drop_incomplete` keeps it (G1 P-h). MOFA2's likelihood uses only observed entries, so under the default `keep_missing` the blank view is handled correctly (INT-010); what is wrong is the "samples measured" count, the few-samples note and the opt-in's promise. Acceptable on staging as a stated limit; the candidate fix (a per-run finiteness pass before `integrationAlignSamples`) should land before the release commit.
- **RC2-06** (above): the Combined DA heatmap binds assays by position; a candidate row.

## Correction pass (AUDIT-KIT-255-impl I01 option (a), I02, I03, I07; staging 176f9d3f merged)

Tests first (7f7fa56b, predictions `r3-predictions.tsv`), then the fix and records (23fe97d7); one run, head against the tests-first commit: 45718599671 (76.4 job-minutes = 5.09 core-hours; with the earlier 34.92: 40.0 of 111.00; expected 8, ceiling 16). Result tables `r3/`.

**As predicted:** the report rows (kit236: fixed, 293/0); `test-kit235-lipid-metab-report-parser.R` same (54/0); `test-omic-norm-module.R` same, the eight digests failing identically on both sides (VAL-303, unedited; the marker is already in the base commit); R CMD check identical (0 ERROR, 3 WARNING, 3 NOTE both sides); in the new file the rows for I01 (K3 group B and two replicates per group: `not_tested` recorded), I02 (three assays: the LCMS_Neg notice names only LCMS_Pos as where the sample reached the threshold, GCMS as not tested; summary heading), I03 (E25 and G7: `no_value` and the sentence), I07 (E8's equality on the declared route, the managed notice, the two-threshold merge) pass at head and, for the rows predicted to fail, failed at the base.

**Misses (stop and report; no second correction used):**
1. E15, lipid and metab: an ERROR at head (not a failure): the added session assertion called `collectWorkflowOmicsNormSession` with a reactive-values object outside a reactive consumer (`Can't access reactive value 'itsd_selections'`), so the assertions after it (the session field, and the clean-fixture summary) did not run. A test-mechanics defect (it needs `shiny::isolate`); the summary-heading assertion before it passed.
2. E20, lipid and metab: one `expect_true(any(grepl(...)))` FAILED at head. By reading, it is the study-parameters assertion: the test expects "group B in LCMS_Neg" and the writer prints "    - LCMS_Neg: group B" (a test text mismatch); the notice line and the summary assertions precede it and the `not_tested` record assertion passed. Which `any` failed is not shown by the message; confirm from the next run.
Neither is a result about the code under test beyond the rows that ran; both are in the test file. The coordinator rules whether they are instrument defects under VAL-54 (a repair of the two test lines and a rerun of the new file alone, head only).
