# KIT-101 findings

Defects found while writing `design.md` (at `bb87a463`). Append only.

The user, 2026-11-07: "fix all findings". So F-01 to F-05 are fixed in the KIT-101 addendum, not left to their family tickets. Each fix is recorded under its entry, with its commit and its test.

| Id | Where | What | Repro | Was owned by | Status |
|---|---|---|---|---|---|
| F-01 | `R/mod_metab_design.R:228-252` | When the metab design-save preflight fails, `stop("Invalid metabolomics design: ...")` runs outside the `tryCatch` that would show a notification. The result is an observer error, and the "Processing Design Matrix" modal (:206-217) stays open. Lipid shows a notification and returns (`mod_lipid_design_builder_helpers.R:82-88`). Metab also writes the bus (:243-247) before it checks `source_dir` (:252). | `testServer` on `mod_metab_design_server` with a builder result whose design has a run that is not a sample column. Expected: a notification and the modal closed. Today: an observer error. | KIT-103 | **fixed** `0259eb01` |
| F-02 | `R/mod_prot_da_ui.R:346`; `R/mod_prot_da_render_volcano.R:41-42` | The proteomics DA UI defines `heatmap_show_dendro`, which no server code reads. The server reads `volcano_show_labels` and `volcano_label_top_n`, which the UI never defines, so their defaults (TRUE, 10) always apply. | grep for the three ids in `R/mod_prot_da*.R`. | KIT-106 | **fixed** `f8d69b75` |
| F-03 | `R/mod_metab_da.R` (no tree-cut inputs; :358 refers to one); `R/mod_metab_da_registration_helpers.R:90-93,100,165-168` | The metab DA tree-cut defect. The server reads `heatmap_tree_cut_method`, `heatmap_n_clusters`, `heatmap_cut_height` and `heatmap_min_cluster_size`, but the UI defines none of them, so they are always NULL. Lipid defines them at `mod_lipid_da_ui.R:310-331`, with defaults none, 4, 10 and 5. | `grep -c 'ns("heatmap_tree_cut_method")' R/mod_metab_da.R` gives 0. | KIT-106 (parity) | **fixed** `0aee1526` |
| F-04 | `R/mod_lipid_summary_observer_helpers.R:35,64,94,251` | No R code calls the four `runLipidSummary*ObserverShell` functions. MCI-023 (`test-module-ci-lipid-summary.R:17,74`) tests them, not the live handlers, so it passes on dead code. | `grep -lw` for each name in `R/*.R` lists only the file that defines it. | KIT-108 | **fixed** `080c897e` |
| F-05 | `R/mod_prot_norm_workflow_helpers.R:568,572-582,584-585`; `R/mod_prot_qc_peptide_log2.R:61-64` | Dead seams. `runProtNormNormalizationWorkflow` takes 14 parameters it never uses. `runPeptideLog2ApplyStep` takes four unused `Fn` parameters and calls only `runDIAReferencePeptideLog2Matrix`. Tests that inject these parameters never reach production code. | Read the function bodies; grep the tests that pass these arguments. | KIT-105 (norm), KIT-104 (log2) | **fixed** `7792e210` |

F-02, F-03 and F-04 were checked by grep on 2026-11-07. F-01 and F-05 come from reading the code; the addendum's failing-first tests reproduce them.

## Fixes (addendum, 2026-11-07)

Each fix has a test that fails before it. The evidence is in [fixes/test-first.txt](fixes/test-first.txt) and [fixes/targeted.txt](fixes/targeted.txt): 17 targeted files, 0 failures.
- **F-01:** behaviour changes on the failure paths only, where the defect was the behaviour. A preflight failure now shows a notification, closes the modal and returns. `source_dir` is checked before the bus is written, as in lipid. The happy path is unchanged.
- **F-02:** the two volcano inputs were added with the renderer's own defaults (labels shown, 10) in `f8d69b75`. That commit removed `heatmap_show_dendro`, because the heatmap had no dendrogram switch. **Revised at the user's direction** ("f-02 of course id like that connected"): the checkbox is back, checked by default. It drives a new trailing argument, `show_dendrogram = TRUE`, on the exported `generateProtDAHeatmap`, passed through `generateDIAReferenceHeatmap` and `renderProtDAHeatmapProduct` to ComplexHeatmap's `show_row_dend`/`show_column_dend`. It is display only: the cached clustering product and its identity are unchanged, and the default draws exactly as before (test: omitted equals `TRUE`; `FALSE` hides both dendrograms; row order unchanged).
- **F-03:** with `none`/4/10/5, the heatmap is identical to today's (the matrix, the row and column clusters, and the md5 of the drawn PNG). One visible side effect: the saved heatmap's methods markdown gains four rows that were dropped before because the values were NULL.
- **F-04:** quarantined Q-23, MCI-023.3, which tested a dead report implementation. The live `handleLipidSummaryGenerateReport` is now tested by MCI-023.3b.
- **F-05:** one argument line in `R/mod_prot_norm_observer_helpers.R`, the only production caller, also changed; that was a direct consequence of the fix. Q-24 is a record row: four bodies already skipped under Q-22 still inject the removed parameters.

## Findings from the tests work (W2 to W5, at `bb87a463`)

Recorded by the tests work (the entry-point suite, the baseline, the fixtures) and left unfixed; the owner named fixes
each. Repros are the KIT-101 suite's frozen steps (`tests/testdata/gui-entry/<lane>.json`) or a command at `bb87a463`.

| Id | Where | What | Repro | Owner | Status |
|---|---|---|---|---|---|
| F-06 | `tests/testdata/e2e/generate_fixtures.R` | The generator does not reproduce the committed fixtures. It has no generator for `import_tool = "custom"`, so `metab_lc`, `metab_gc` and `metab_combined` fail; and because those lanes no longer draw from the global RNG stream, the three `lipid_canonical` files it writes differ from the committed ones. The proteomics files reproduce byte for byte; the new `trx_screen` uses its own seed. The regenerated lipid files were restored. | `Rscript tests/testdata/e2e/generate_fixtures.R` from the root: three generation errors; `git diff --stat` shows the three LipidSearch files changed. | test infrastructure (the next ticket that regenerates fixtures) | fixing (revision 2.2, user: "fix all") |
| F-07 | `tests/testdata/e2e/prot_dia/report.tsv`, `prot_dia_limpa/report.tsv`; `R/func_prot_qc_peptide_confidence_filters.R:200-205`; `rollUpPrecursorToPeptideHelper` | The committed DIA-NN fixtures lack `Global.Q.Value`, `Global.PG.Q.Value`, `Proteotypic` and `Precursor.Normalised`. Today's Q-value filter (`diann_main_report`) refuses the table and the precursor rollup stops, so the later peptide steps fail or run on an unrolled table. The e2e DIA browser lanes wait for `qvalue_filtered` after this step, so they cannot pass on these fixtures (not run in a browser here; they skip without `CHROMOTE_CHROME`). | Lanes `prot_dia`, `prot_dia_limpa`, steps `qc_peptide_qvalue_filter`, `qc_peptide_rollup` (their notifications). | KIT-104 with the fixture owner | fixing (revision 2.2, user: "fix all") |
| F-08 | `R/utils_workflow_dia_da.R:681-682` (`writeDIAInteractiveVolcanoOutputs`); `R/func_prot_da_volcano_glimma.R:230-233` | On every committed proteomics fixture (5 or 6 proteins) the DA's FDR is all NA (q-values cannot be fitted on so few p-values; compare Q-19), `generateProtDAVolcanoPlotGlimma` returns NULL ("no valid rows after removing Inf/NA"), and the whole DA publication aborts with `saved_entity_outputs_pending` ("interactive volcano output is missing"). The DA load and the publication copy then refuse with "da outputs are pending; retry that action". No proteomics lane reaches a DA result, enrichment or report. The defect question: should an all-NA FDR abort the publication rather than skip the interactive volcano? | Lanes `prot_dia` (frozen), `prot_tmt`, `prot_lfq`, `prot_lfq_fragpipe`: steps `da`, `summary_copy`. | KIT-106 | **fixed** `1d87ee43` |
| F-09 | the DIA_LIMPA protein sub-tabs (`dpc_quant`, `normalisation`, `peptide_count`) | When the peptide stage has not produced its state, each protein sub-tab's apply fails with "object 'parameters' not found" instead of a refusal naming the missing state: an error path refers to a variable that is not yet bound. | Lane `prot_dia_limpa`, steps `qc_protein_dpc_quant`, `qc_protein_normalisation`, `qc_protein_peptide_count`. | KIT-104 | **fixed** `4e2c4c33` |
| F-10 | proteomics normalisation tab | Every action on the tab (run, correlation filter, export) raises the warning "Please complete the Quality Control filtering steps before accessing normalization." even after QC completes (`tab_status$quality_control = "complete"`, head `protein_replicate_filtered`). | Lanes `prot_dia`, `prot_tmt`, `prot_lfq`, `prot_lfq_fragpipe`: steps `norm`, `correlation`, `session_export`. | KIT-105 | **fixed** `1ecd454c` |
| F-11 | the summary tab (every omic); the proteomics import; the trx and protein-route QC tabs | No messaged refusal is reachable through the inputs. The report before the publication copy and the proteomics import with no file are `req()`-gated and say nothing; the trx expression filter and the protein-route QC sub-tabs have no invalid input to refuse. The suite records the first two as silent refusals (head and history unchanged, no notification) and has no QC refusal on those routes. | Every lane's `summary_refusal`; proteomics `import_refusal`. | KIT-102 (shell and refusal parity) | **fixed** `07f72722` |
| F-12 | `tests/testthat/test-general-filemgmt-duplicate-resolution-contracts.R:4`, `test-prot-04-design.R:4` | Under `R CMD check`'s test layout these two path-guarded files do not skip: their first line `repoRoot <- pkgload::pkg_path()` errors outside a source package before the guard runs (13 and 182 blocks become one file-level error each). The other ten path-guarded files skip by construction. | `kit-101/test-baseline.tsv`, mode `check`. | the family ticket that next touches these files | **fixed** `6b0b0e1a` |
| F-13 | the entry-point suite's capture (`tests/testthat/helper-gui-entry.R`, `guiEntryArtifacts`) | prot_tmt, prot_lfq and prot_lfq_fragpipe are not frozen: the count of files in the owner's store after the session export varies by one to three between runs, after the ticket's one correction pass. Every other capture of these lanes was identical over three runs. Their tests fail ("no frozen oracle"). | `entry-suite.md`: runs 3, 4 and the freeze run. | KIT-101 (a decision: drop or narrow the store count, then freeze the three) | fixing (revision 2.2, user: "fix all") |
| F-14 | `kit-099/suite.txt` header; the KIT-101 experiment budget | The ticket gives the baseline "6 in parallel, the KIT-099 checkpoint's settings", but the KIT-099 checkpoint ran 8 in parallel. The baseline ran 6, as the ticket's number says. | `head -1 docs/qa/workflow-execution/kit/kit-099/suite.txt` | none (record) | fixing (revision 2.2, user: "fix all") |

## Fixes, revision 2.2 (code), 2026-11-07

The evidence is in [fixes2/test-first.txt](fixes2/test-first.txt) and [fixes2/targeted.txt](fixes2/targeted.txt). Each fix has a test that fails before it; the targeted files and dependents give 0 failures. No test was quarantined.
- **F-08:** when a contrast has nothing to plot, `writeDIAInteractiveVolcanoOutputs` writes `Interactive_Volcano_Plots/<contrast>.unavailable.txt` holding the reason, and the publication requires that record instead of the page. The DA tab shows the same reason. Every other DA output is unchanged and still required.
- **F-09:** the root cause was that since KIT-048 the three in-process steps passed `parameters`, a name bound only in the observer. Every unmanaged apply failed, whatever the state. Each step now builds its own parameter list. Normalisation and peptide-count name the `protein_dpc_quant` state they need.
- **F-10:** the stale read was `shouldProtNormAutoGeneratePreQc`, which also returns FALSE once pre-QC is generated, and the tab took every FALSE to mean "QC incomplete". The tab now warns only when the state is not a normalisation state and QC is not complete.
- **F-11:** a warning before the unchanged `req()` gates, in the lipid/metab refusal style:
  - the summary report on every omic: no experiment label, or no copy to the publication folder yet;
  - the proteomics import: no file chosen;
  - the trx filter before import.

  The protein-route QC sub-tabs already refused out-of-range inputs with a reason, so they needed no change.
- **F-12:** line 4 of each file uses the idiom of the other ten guarded files. Under check, both now skip at their guard.

| F-15 | proteomics DIA summary, after enrichment completes (lane `prot_dia`) | The analysis record and bundle are not written: `summary_save` and the report warn "analysis record ... not written: component checkpoint is incomplete" on the `enrichment_completed` state. This was hidden while DA failed (F-08). | Lane `prot_dia` on HEAD `07f72722`: steps `summary_save` and `summary_report`. | KIT-108 / KIT-107 | **fixed** `cbd60e56` |

**F-15 fixed (`cbd60e56`).** The record builder was at fault, not the enrichment commit. `buildWorkflowAnalysisRecord` read every state through `readWorkflowStateMetadata`, which accepts only the declared assay classes. `enrichment_completed` is a complete checkpoint of the declared result class `EnrichmentResults` (KIT-003), committed by the same path as every other state.
- The builder now reads each state with `readComponentMetadata`, which accepts every declared class and applies the same checks to assay states.
- The design comes from the last state that is not a result class.

Records without a result-class state, which is every record that could be built before, are unchanged. A test fails first; the record, design and earlier steps are identical to the pre-enrichment record; 24 targeted files give 0 failures. On `prot_dia`, `analysis_record.json` and the report bundle are now written.

| F-16 | `runDIAReferenceEnrichment` (proteomics enrichment) | The enrichment state carries no step provenance. The analysis record lists it as a legacy step: `validateWorkflowAnalysisRecord` reports "legacy step: enrichment_completed" as incomplete, and a replay of a `prot_dia` record that includes enrichment stops there (`workflow_replay_unsupported`). Found by the F-15 fix. | Build and validate the record on lane `prot_dia` after enrichment. | KIT-107 | passed to KIT-107 (user, 2026-11-07: "f16 pass to kit-107") |
| F-17 | the trx QC filter (`mod_trx_qc_filter.R`), lane `trx_screen` | A click on the expression filter before import reaches no observer and is replayed when the design is saved, so the filter runs without a fresh request. Found by the fixtures agent (`26f09136`, which dropped that refusal step from the lane). | Lane `trx_screen`: click the filter before import, then save the design. | KIT-104 / KIT-110 | **fixed** `dfabe0b3` |

## Revision 2.2, phase 1 (fixtures report, 2026-11-07)

The fixtures report found these. F-17 is confirmed: the trx QC sub-servers start only after the design, and the observer is created without `ignoreInit`, so it runs an earlier click when it starts.

| F-18 | proteomics LFQ imports (MaxQuant `proteinGroups.txt`, FragPipe `combined_protein.tsv`; `R/func_prot_import_readers.R`) | Zero intensities appear to be kept as quantities rather than treated as missing. With 2% to 8% zeros in a synthetic fixture, the correlation filter dropped 1 to 4 of 10 samples. MaxQuant and FragPipe write 0 for "not quantified", so a zero kept as a value distorts every downstream step. Scientific; **not yet confirmed by diagnosis**. | Lane `prot_lfq` or `prot_lfq_fragpipe` with zeros in the intensity columns. | KIT-102 (import) | **fixed** `8a2a3581`; `lfq_300` oracles to be re-recorded |
| F-19 | DIA_LIMPA route DA (lane `prot_dia_limpa`, new fixture) | 41 proteins at q < 0.05, of which only 18 are among the 40 planted changes. Every other route has at most 1 false positive on the same synthetic universe. The route ran RUV-III covariates (k = 2) on 8 samples. It is unknown whether the fault is the fixture, limpa/DPC on synthetic peptides, or the RUV step. Scientific; **not yet diagnosed**. | Lane `prot_dia_limpa` on the revision 2.2 fixture. | KIT-105 / KIT-106 | fixing: fixture (both DIA-NN lanes, user); caveats to KIT-114 |
| F-20 | proteomics normalisation tab, DIA_LIMPA route | F-10 residual: after the limpa protein QC completes (head `peptide_count_filtered`), the tab still warns "Please complete the Quality Control filtering steps...". | Lane `prot_dia_limpa`, step `norm`. | KIT-105 | **fixed** `7ff973a3` |
| F-21 | `tests/testthat/helper-e2e-browser.R:675` (`e2e_assert_prot_filtered_session_artifacts`) | The e2e harness looks for the legacy `filtered_session_data_latest.rds`, but the native export writes `filtered_session_*_native-export-*.msdia.tar`. Every proteomics browser lane (DIA, enrichment, TMT, LFQ, FragPipe) reaches the normalisation export and then fails on this check. | Any e2e proteomics lane on the revision 2.2 fixtures. | test harness | **fixed** `cd342970` (e2e evidence `01b60213`) |
| F-22 | entry suite (`guiEntryCanonical`) | Two captures are not yet stable. (a) `prot_dia`'s `enrichment_completed` digest: plotly objects carry random 14-hex widget ids as list names. (b) `prot_tmt`'s store entries: a request id ending in 9 hex characters escapes the 10-hex mask. Without a fix these two lanes cannot be frozen. | Double run of `prot_dia` and `prot_tmt`. | KIT-101 (a third correction; needs approval) | fixing (third correction approved) |

## Diagnosis (`bba662dc`, [diagnosis/](diagnosis/)) and the user's decisions, 2026-11-07

**F-18, confirmed.** No step on the LFQ route turns a 0 into NA.
- **Mechanism:** the readers apply `as.numeric` only. The design checkpoint applies `log2(x + 1)` (`mod_prot_design_state_helpers.R:69-76`), so 0 becomes 0. Cyclic loess then turns those 0s into values from -4.5 to -15.5, which limma treats as measurements.
- **Effect** (synthetic, 8% missing written as 0 vs written empty):
  - MaxQuant: 1 of 40 planted changes found, against 16. The correlation filter drops 3 samples.
  - FragPipe: all five WT samples dropped, and DA fails.
- **Semantics:** MaxQuant and FragPipe write 0 for "not quantified". Perseus and other importers convert 0 to NA. DIA_LIMPA already maps values ≤ 0 to NA.
- **Oracles recorded with zeros kept:** KIT-021 `lfq_300` (8.01% of cells are 0), and very likely KIT-054 `kc_lfq` (real data, not inspected). MCI-004.3 expects a 0 to survive the import.
- **Decision ("Fix + re-record"):**
  - `dplyr::na_if(..., 0)` in `importMaxQuantData` and `importFragPipeData`, with a test that fails first, and MCI-004.3 inverted;
  - re-record the synthetic `lfq_300` oracles under new names, keeping the old ones as history;
  - annotate the KIT-054 real-data oracles as recorded before the fix; real data is re-run at the final real-data gate;
  - put missing values back into the LFQ lane fixtures.

**F-19, mixed cause** (each factor isolated one at a time; the GUI result reproduced exactly):
- **Fixture:** `kitProteinMatrix`'s per-sample loading gives KO − WT = −0.30 log2 at this seed. The generator writes `Precursor.Normalised` without removing it, as real DIA-NN output would. DIA_LIMPA normalises with `none`, so null proteins carry logFC ≈ −0.27, and the unnormalised true matrix alone gives 10 false positives.
  - Centring gives 4 false positives with dpcDE, and 2 on the DIA route (cyclic loess).
  - RUV-III covariates add about 9 to 11 (W uncorrelated with group; 20 controls at n = 8).
  - `qvalue()` is more liberal than BH: 41 calls against 23.
- **Decision ("Both DIA-NN lanes"):** `Precursor.Normalised` is written with the per-run loading removed, for the prot_dia and prot_dia_limpa fixtures. Both lanes are then frozen.

| F-23 | DIA_LIMPA between-sample normalisation (`normaliseBetweenSamples` on the protein object; `dpcDE`) | A normalisation chosen on DIA_LIMPA changes only `protein_quant_table`, but the DA (`dpcDE`) fits `limpa_dpc_quant_results$quantified_elist`. So a chosen normalisation has no effect on the DA results: cyclic loess changed the table by up to 0.51 and the DA calls were unchanged (27/15/12 both ways). Found by the F-19 diagnosis (`diagnosis/f19_norm_check.R`). | `f19_norm_check.R` on the prot_dia_limpa fixture. | KIT-105 | **fixed** `5931d8bf` (default `none` unchanged) |

**Statistical caveats → KIT-114** (user: "Ticket both"): (a) the DA q column uses `qvalue()`, which is more liberal than BH at small n; (b) RUV correction mode gave 88 false positives at n = 8 on this fixture.

**Process note:** the diagnosis agent used a read-only `pgrep -f` once in a wait loop, against the host rule. It stopped and waited on its output files instead.

## Fixes, revision 2.3 ([fixes3/](fixes3/)), 2026-11-07

- **F-17:** the server notes the click value present when it starts, and its observer ignores that value. `ignoreInit = TRUE` was tried first and rejected: it also drops a first click made in the same flush as the start, which broke MCI-028.1/.2/.3 and the F-11 block (recorded). The brief's premise that lipid and metab use `ignoreInit` was wrong; they don't.
- **F-20:** the DIA_LIMPA route never writes `tab_status$quality_control = "complete"`. `mod_prot_qc_server` marks QC complete only at `protein_replicate_filtered`, and the limpa protein stage ends at `peptide_count_filtered`. On DIA_LIMPA, a head at `peptide_count_filtered` or `ruv_covariates` now counts as QC complete for the warning. Nothing else is written.
- **F-21:** the assertion reads the latest-session pointer, checks its schema and status and the archive's name, size and sha256, and checks that the summary names the same archive. Two negative cases each fail exactly one check.
- **The e2e proteomics browser lanes on the revision 2.2 fixtures** (package from source):

| Lane | Result | Time |
|---|---|---|
| DIA | 41 passed, 0 failed | 259 s |
| enrichment (gprofiler2 and clusterProfiler) | 96 passed, 0 failed | 526 s |
| TMT, LFQ and FragPipe | 150 passed, 0 failed | 369 s |
| DIA_LIMPA | 11 passed, 1 skip (no browser QC path yet, by design) | 61 s |

| F-24 | lipid and metab QC apply observers (`mod_lipid_qc_intensity_server_helpers.R:127`, `mod_metab_qc_intensity.R:358`); the trx revert observer | The same class as F-17: a click made before the server starts can be run when it starts. Found while fixing F-17; not reproduced through the public servers yet. | Construct the F-17 test on the lipid and metab QC intensity servers. | KIT-104 | **fixed** `321c2294` (12 sites); 2 display-only sites queued to KIT-104 |
| F-25 | `tests/testthat/test-e2e-proteomics-dia-limpa.R`, the checks after the export | They still expect legacy fields (`use_limpa`, `current_s4_object@args`, `workflow_type` 'DIA'). They are not reachable while the lane skips at QC, and will fail once a limpa browser QC path exists. | Read the lane file after its QC skip. | KIT-104 | queued to KIT-104 (user, 2026-11-07; not reachable until a limpa browser QC path exists) |

## Fixes, F-18 and F-23 ([fixes4/](fixes4/)), 2026-11-07

**F-18 (`8a2a3581`):** `dplyr::na_if(as.numeric(Intensity), 0)` in both readers' pivots. Inputs without zeros give the same result as before.
- **Tests, each failing first:** the reader block (6 fail before); MCI-004.3, inverted because it encoded the defect; and a route block from the owner import through cyclic loess. Before the fix the route block gave log2 0 at the zero cells and a minimum of -11.9 after normalisation. After it, NA sits exactly at those cells and the minimum is above 15.
- **Oracles changed:** KIT-021 `lfq_300_*` (8.01% zeros). Likely also KIT-054 `kc_lfq`/`kc_lfq_dotted` (real data, not read).
- **Unchanged (0 zeros):** every other LFQ/TMT fixture and corpus, and every e2e lane.

**F-23 (`5931d8bf`):** in `normaliseBetweenSamples` (ProteinQuantitativeData), a chosen method is also applied to `limpa_dpc_quant_results$quantified_elist$E`. With `none`, E and the whole DPC result are byte-identical, the owned-args digest is unchanged, and no output or oracle changes.
- **Validity, read from the limpa 1.4.2 source:** `dpcDE` uses E, the `n.observations` imputation flags, and log standard errors as a precision predictor. Only E changes, as RUV-III correction already does.
  - Equivalent for shift-type methods; cyclic loess is approximately one.
  - Approximate for quantile and scale, whose SEs strictly should be rescaled.
  - limpa's documented normalisation point is the precursor matrix before DPC (DIA-NN `Precursor.Normalised`).

| F-26 | `importProteomeDiscovererTMTData` (`R/func_prot_import_tmt.R`) | Latent, as F-18: a literal 0 abundance would be kept as log2 0 (no 0 to NA step; the TMT design branch applies `log2(x + 1)`). Proteome Discoverer normally leaves missing abundances empty, and every synthetic TMT corpus and lane has 0 zeros. Found by the F-18 fix. | Construct a PD table with a 0 abundance. | KIT-102 | **fixed** (see the F-26 fix line below; user: "Fix now") |

**F-26 fixed** ([fixes6/f26.txt](fixes6/f26.txt)): in `importProteomeDiscovererTMTData`, a numeric 0 abundance is read as NA after the pivot, as the LFQ readers now do. A test fails first (2 failures); after the fix, the TMT direct test passes 27/0, MCI import 91/0, and prot-01 15/0. No corpus or lane has a TMT zero, so no oracle changes.
**F-23 caveat** (user: "Keep, note caveat"): quantile and scale are approximate on DIA_LIMPA. Moving its normalisation to the precursor level is added to KIT-114 as item (c).

**F-24 fixed (`321c2294`, [fixes5/](fixes5/)).** Twelve apply, revert and finalise sites were reproduced through the public QC orchestrators and fixed with one helper, `qcClickAfterStart`, which also replaces F-17's inline lines. The sites are lipid and metab intensity apply/revert, duplicates resolve/revert and finalise, and trx revert and finalise. Before the fix the tests failed: lipid 5, metab 5, trx 2. After it: lipid 68/0, metab 64/0, trx 38/0, every dependent 0 failures, and the entry lanes unchanged.
- **Not changed:** lipid and metab *Detect duplicates* and *Analyze internal standards* also run a click made before the sub-modules start, but they commit no state and only fill the display. They are queued to KIT-104, which rebuilds the QC shell.
- **Reachability:** a browser click before the start was not shown. The reproduction was through `testServer`, as for F-17.

| F-27 | the LFQ route's RUV on `lfq_300` after F-18 (`findBestNegCtrlPercentageFn`, the KIT-005 k-validity rule) | The F-18 re-record on `kit_101` stopped at normalisation after 920 s: "No valid percentage found". With zeros read as NA, the synthetic `lfq_300` corpus (8% missing, at random over every protein, 300 samples) leaves no control protein observed in every sample (0.92^300 ≈ 0). The rule then refuses k at every percentage. Before F-18 every protein looked complete, because the zeros counted as values. This is a valid unfavourable result, not an instrumentation fault, and the run was not retried. Two questions follow: (a) the corpus's missingness is not realistic, since real LFQ has a complete high-abundance core and intensity-dependent missingness; (b) the GUI refuses with a bare error rather than saying that no complete control proteins exist. | `scripts-101/record-lfq.sh` on `kit_101` (`record-lfq.out`). | KIT-101 (oracles), KIT-105 / KIT-114 | **fixed**: message `5213a6e2`; corpus `lfq_300_mnar` (`23558bbe`) recorded on kit_101 (`57e249be`, RUV 3%); the control-selection method goes to KIT-114 (d) |

| F-28 | `getUniprotAnnotationsFull`'s error handler (`R/func_prot_annotation_uniprot_full_helpers.R`) | When a live UniProt chunk timed out, the handler returned a "minimal" table without `Entry` or `gene_names`. The DA exports join any non-empty annotation table on `Entry`, so DA failed ("Join columns in `y` must be present"), its outputs stayed pending, and the summary refused, failing the report download. This broke all three TMT/LFQ e2e sub-lanes on the final-build run (`e2e-final.txt`). | Lane `test-e2e-proteomics-tmt-lfq.R` with a UniProt timeout. | KIT-101 | **fixed** `b91a2b13` (returns `createEmptyUniprotTable()`, not cached, so a retry fetches again; annotations are display metadata, no statistic changes); e2e lanes moved off the live network (in progress) |
