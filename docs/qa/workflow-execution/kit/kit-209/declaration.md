# KIT-209 declaration: in-process memory discovery (instrument, settings, units, G0 runs)

Committed alone before any G0 run (ticket section 9.1, "Order of commits", item 1). Ticket revision 1.3
(`tickets/kit/KIT-209.md`, sha256 `bf7cc0d0f21ce6dd6bcf0c7b42c60f8378cd82d2349db616c279666006df9499`). Authority: methods
ledger VAL-145 (an in-process lane on GitHub runners; it cannot show the application/worker split) and VAL-151 (a memory
discovery, not a gate: figures are reported with their spread; nothing is flagged; nothing holds the release).

## 1. Builds

| Label | Commit | Role |
|---|---|---|
| `C` | `c03cb02a1fc6e3c122b39a6d424b4996971e20e1` | the staging head on 2026-11-24 (G0 only; the release candidate is bound later, in the dated section) |
| `base` | `b468a885c28ea7dfe228a0fdaad947825bc4839e` | the base, b468a885 |

G0a runs the base over the candidate's snapshot library (B1 on the candidate's pins). G0a' (the base on its own pins
`fastduckdb` d7851150 and `fast_limpa_dpc` d4f9cc1d) needs the lane's token handling extended by the coordinator
(ticket section 5) and is not run by this declaration.

## 2. The instrument

Run branch `run-kit209-memory`, cut from `run-mofa-measure-k150` (KIT-150's lane). Files:

| File | Role |
|---|---|
| `.github/workflows/k209-measure.yml` | the job: runner facts, the KIT-134 snapshot library at `C` (token only to that action), the bundle's pandoc (distributable.yml's `PANDOC_VERSION`), both builds installed side by side (`R CMD INSTALL` of each commit into its own library), B2's table, the plan, the summary, the upload |
| `.github/k209/k209-run-journey.sh` | one journey: page cache dropped (`sync`; `echo 3 > /proc/sys/vm/drop_caches`), a cgroup v2 child with `memory.max` (14.5 GiB unless injected) and `memory.swap.max` 0, `timeout` (TERM, KILL 30 s later), one fresh `Rscript --vanilla` process; then `memory.peak`, `memory.events`, `oom_kill` and the outcome (k150-run-stage.sh's classes) |
| `.github/k209/k209-journey.R` | the journey: the build's library first on the path (refused if MultiScholaR loads from elsewhere or the library's commit marker differs), the harness's pure functions from the candidate's checkout (`kitJourneySteps`, `kitFindStore`, `kitReadHead`, `kitFingerprintStore`; `tools/profiling/` unchanged since b468a885), the route's steps, the marks |
| `.github/k209/k209-corpus.R` | the corpora, generated on the runner by the candidate's committed generators; each file's SHA-256 compared with `bindings.json`; `kc*` refused by name; a corpus not marked synthetic in the bindings refused |
| `.github/k209/k209-library.R` | B2: each build's Depends, Imports and LinkingTo against the versions its journey's library path finds |
| `.github/k209/k209-plan.sh` | the journeys of one job, in order (plans `g0a`, `g0p1`, `g0p2`) |
| `.github/k209/k209-summarise.py` | journey and step rows, the outcome classes, G0a's checks (standard-library Python) |
| `.github/k209/k209-remove.py` | the guarded remover of a journey's folder |
| `.github/k209/k209-builds.tsv`, `k209-params.json` | the two builds; the run's label, attempt and jobs |

### 2.1 Marks (ticket section 3)

- **Baseline** (after the package, the harness and the arguments are loaded; outside every step): `gc(full = TRUE)` twice
  with `reset = TRUE`; R heap in use (`Ncells + Vcells`, the `(Mb)` column, MiB) and `VmRSS`; `VmHWM` reset by writing
  `5` to `/proc/self/clear_refs`. Written to `<id>.baseline.json`.
- **Each step, at entry:** `gc(full = TRUE)`; heap in use (the retained R heap at that boundary) and `VmRSS`; then
  `gc(reset = TRUE)` and the `VmHWM` reset (the reset value is recorded).
- **Each step, at exit:** `VmHWM` (the step's process peak), `gc()` max used, `VmRSS`, wall seconds, the number of
  processes in the journey's cgroup. Then, outside the step's figures: the head's state (`kitReadHead`), its rows, the
  store's bytes, and whether the state is the harness map's (`kitJourneySteps(route, ruv_mode)`). `<id>.marks.json` is
  rewritten after every step, so a killed journey keeps the steps it finished.
- **Final:** `gc(full = TRUE)`, heap in use and `VmRSS` after the last step.
- **After the figures** (G0a only): `kitFingerprintStore` on the final store, each state's fingerprint reduced to its
  sha256 (reported, base beside candidate).
- Journey maximum = the largest step `VmHWM`; retained-heap maximum = the largest boundary heap (entries and final).

### 2.2 Outcome classes (ticket section 8)

| Class | Rule in `k209-summarise.py` |
|---|---|
| completed | the journey wrote its record with status completed |
| memory | `oom_kill` > 0 in the journey's cgroup (whether the baseline was written is recorded beside it) |
| time | `stopped: time limit` (timeout's status 124, or 137 after the KILL) |
| package error / instrument | the journey wrote an error record: the error's frames at signal time (an R calling handler, so an error a tab catches and only notifies is still seen); the deepest frame owned by the lane or by the package decides: the MultiScholaR namespace = package error, the lane's code = instrument |
| state not committed (E3) | a step completed without error but left a state other than the harness map's, and no error was signalled |
| instrument | failed before the baseline |
| runner lost | no driver record |
| unclear | anything else; reported as it is to the coordinator |

### 2.3 What the lane holds between steps (R07)

Each step body is a function. Between steps the lane holds: `workflow` (the module-equivalent bus; on the omics routes
the module's own `workflow_data`), the session owner inside it, the experiment paths and project directories, the
corpus file paths, the bound settings, the normalisation tab's `norm_data` on the protein routes (as the tab holds it),
and the globals the package itself reads or writes (`project_dirs`, `filtering_progress`). The record lists the lane's
holdings at the end of each journey. On the omics routes the module server's own reactive state is held by
`shiny::testServer` for the whole journey; it is part of what the retained heap measures there.

### 2.4 The routes' step bodies (B3, B5; G0a item 4)

**Lipidomics, metabolomics, transcriptomics: the real module server.** The kit-031, kit-032 and kit-099 recorders call
import and design helpers that exist at neither build (removed by the shared-module rework before b468a885):
`loadLipidImportAssayPreview`, `runLipidImportProcessing`, `prepareMetabImportAssaySelectionState`,
`runMetabImportProcessing`, `trxImportPreview`, `runTrxImportProcessing`, `trxOmicProfile`, `buildLipidDesignInitialState`,
`applyLipidDesignMetadataAssignment`, `applyLipidDesignTechRepAssignment`, `buildLipidDesignContrastState`,
`appendLipidDesignContrast`, `buildLipidDesignSaveResultsPayload`, `runLipidDesignBuilderObserverShell`,
`registerMetabDesignBuilderResultsObserver`, `resolveLipidImportSampleColumns`, `resolveMetabImportSampleColumns`
(`git grep` at HEAD and at b468a885: no definition). The re-check's "every package function they call exists at both
builds" does not hold for these. The lane therefore runs these routes as KIT-075's residency journey does
(`tests/testthat/helper-workflow-omics-residency.R`, at both builds): `shiny::testServer` on `mod_<omic>_server` (the
shared `mod_omic_workflow_server` under `omicGuiConfig(<omic>)`), the import through the tab's processing helper
(`omicImportFeatureTablePreview`, `processOmicImportFeatureTable`, `completeOmicImport`), every later step through the
tabs' own inputs with the kit driver's input ids (`tools/profiling/workflow_kit_driver.R`, `runKitOmicsStep` and
`runKitTrxStep`), each tab's static UI defaults set first (KIT-075's `omicsUiInputDefaults`, copied). This is the routes'
in-process GUI path at both builds; it runs the tabs' observers (no figure is drawn: testServer renders no output that
is not read). Steps and actions:

| Route | Steps (the harness map's ids) and the inputs clicked |
|---|---|
| lipidomics, metabolomics | import (processing helper); design (`design-builder-*`: factors, assign, technical replicates per corpus replicate, contrast second minus first, `~ 0 + group`, save); intensity (`qc-intensity-apply_filter`); duplicates (`qc-duplicates-detect_duplicates`, `resolve_duplicates`); qc_finalize (`qc-itsd-is_pattern` `^IS_`, `qc-itsd-analyze_is`, then `qc-s4_finalize-finalize_qc`); normalisation (`norm-*` settings, `norm-run_normalization`); correlation (`norm-min_pearson_correlation_threshold`, `norm-apply_correlation_filter`); session_export (`norm-export_session`); da (`de-load_filtered_session`, formula and thresholds, `de-run_da_analysis`); report (`summary-experiment_label`, `description`, `save_workflow_args`, `copy_to_publication`, `generate_report`) |
| transcriptomics | import; design; expression_filter (`qc-filter-apply_filter`); qc_finalize (`qc-s4_finalize-finalize_qc`); normalisation (`norm-ruv_mode`, `auto_percentage_min/max/step`, `k_penalty_weight`, `norm-run_normalization`); da (`de-da_q_val_thresh`, `treat_lfc_cutoff`, `de-run_da_analysis`); ruviii_export (`summary-ruviii_m_source`, `summary-run_ruviii_export`); report (as above); session_export (`summary-export_session_state`) |

**DIA_LIMPA: kit-118's body** (`docs/qa/workflow-execution/kit/kit-118/oracle-record-limpa.R`), each step in a function; the
normalisation tab's pass-through (`runDIAReferenceLimpaNormalisationPassThrough`), which kit-118 ran between steps, opens
the correlation step's measurement.

**Protein (LFQ, TMT): kit-021's import, design, normalisation, correlation, session export and DA**; its four protein QC
steps through the QC tabs' in-process step helpers (below) instead of the registry handlers, as the DIA body; the report
through the shared summary core (kit-118's report step; KIT-108 removed the helpers kit-021's report called).

**DIA (IQ): new body (R08).** The DIA tabs' own in-process step helpers with the arguments their observers pass, each
followed by the tab's progress publication (`publishDIAQvalueOutput` on the peptide tabs, `publishDIAProteinProgress` on
the protein tabs), inside the session owner's scope as the observers run them; import and design as kit-118's (route
`DIA`, as the import tab's `dia_route` default):

| Step | Entry point (candidate `C`) | At b468a885 | Arguments (the tab's observer or UI default) |
|---|---|---|---|
| import | `initializeDIANativeImport` `R/utils_workflow_dia_sessions.R:976`, then kit-118's import calls | `:966` | format `diann`, route `DIA` |
| design | `buildProtDesignSaveResultsPayload` `R/utils_workflow_dia_state.R:1229` and kit-118's save calls | `:1225` | groups, technical replicates per corpus replicate (kit-118), one contrast |
| qvalue | `runPeptideQvalueApplyStep` `R/mod_prot_qc_peptide_qvalue.R:78` | `:78` | `.diaNnIdentificationQvalueCutoff` twice, proteotypic TRUE (the observer, `:253-270`) |
| protein_peptide | `runProteinPeptideApplyStep` `R/mod_prot_qc_peptide_protein.R:80` | `:80` | 1, 2 (UI `:41-50`) |
| rollup | `runDIAReferencePrecursorRollup` `R/utils_workflow_dia_qc_rollup.R:237` (the tab's call, `mod_prot_qc_peptide_rollup.R:77`) | `:237` | none |
| intensity | `runPeptideIntensityApplyStep` `R/mod_prot_qc_peptide_intensity.R:201` | `:201` | FALSE, 2, 2, 1 (UI `:42-89`) |
| sample | `runPeptideSampleApplyStep` `R/mod_prot_qc_peptide_sample.R:80` | `:80` | 3000 (the bound journey setting `sample_filter-min_peptides_per_sample`; UI default 500) |
| replicate | `runPeptideReplicateApplyStep` `R/mod_prot_qc_peptide_replicate.R:127` | `:127` | `replicates` (UI `:41-43`) |
| imputation | `runPeptideImputationStep` `R/mod_prot_qc_peptide_impute.R:77` | `:77` | 0.5, no exclusion, not skipped (UI `:41-52`) |
| iq | `runProteinIqRollupApplyStep` `R/mod_prot_qc_protein_rollup.R:63` | `:63` | the experiment paths |
| accession | `runProteinAccessionCleanupStep` `R/mod_prot_qc_protein_cleanup.R:94` | `:94` | `;`, `mean` (UI) |
| protein_intensity | `runProteinIntensityFilterApplyStep` `R/mod_prot_qc_protein_intensity.R:120` | `:120` | FALSE, 2, 2, 1 (UI) |
| duplicates | `runProteinDuplicateRemovalStep` `R/mod_prot_qc_protein_dedup.R:74` | `:74` | `mean` (UI) |
| protein_replicate | `runProteinReplicateFilterApplyStep` `R/mod_prot_qc_protein_replicate.R:80` | `:80` | `group`, 4 (UI) |
| normalisation_ruv | `runDIAReferenceNormalization` `R/utils_workflow_dia_products.R:734` | `:545` | the DIA settings (section 3); no QC figure drawn (kit-021) |
| correlation | `runDIAReferenceCorrelation` `:1586` | `:1014` | 0.5 |
| session_export | `runProtNormExportSessionWorkflow` `R/mod_prot_norm_session_helpers.R:178` | `:178` | |
| da | `runDIAReferenceDaAnalysis` `R/utils_workflow_dia_da.R:357` | `:343` | `~ 0 + group`, 0.05, 0 |
| report | `copyWorkflowSummaryResults` `R/utils_workflow_summary.R:193`, `renderWorkflowSummaryReport` `:304` | `:186`, `:297` | the record export not run (kit-118) |

No GUI-path fingerprint record of the DIA `screen` corpus exists for `C` or b468a885 on the runner image (G0a item 4's
comparison is then "none available"; the base-to-candidate fingerprints are reported).

**Build-conditional branches: none.** Both builds run the same lane code.

## 3. Settings, bound identically on both builds

| Route | Settings |
|---|---|
| lipidomics, metabolomics (`031l`, `032m`, `086p`, `086q`) | ITSD on, `cyclicloess`, RUV `automatic` (1 to 20 %, `max_difference`, k penalty 0.5, adaptive, grouping `group`), correlation 0.5, DA `~ 0 + group` q 0.05 treat 0; the other tab inputs at the UI defaults |
| transcriptomics (`099s`, `099l`) | RUV `automatic`, 5 to 30 % step 5, k penalty 0.05; DA q 0.05 treat 1; RUVIII export M source `group` |
| DIA (journey defaults, KIT-122's full journey) | `cyclicloess`, RUV `automatic` (1 to 20 %, `max_difference`, 0.5, adaptive, `group`), `ruv_percentage` 5, `ruv_k` NA, correlation 0.5; the peptide tabs as section 2.4 |
| DIA_LIMPA (`053`) | between-sample `cyclicloess`, RUV `covariates`, correlation 0.5; DPC engine `fastlimpadpc` bound by name (`RAYON_NUM_THREADS=1`), slope 0.8, chunk 1000; peptide count 1, peptidoforms 2; refused (`workflow_engine_unavailable`) if the engine is not installed |
| protein, U11 and U12 | `cyclicloess`, RUV `manual`, `ruv_percentage` 5, `ruv_k` NA (the tab's manual-mode defaults at both builds: `omicGuiConfig("proteomics")` `ruv_iii` percentage 5, k NULL, which the kernel resolves to 3: `R/utils_workflow_dia_norm.R:82`; both recorded in G0a item 5), correlation 0.5 |
| protein fixtures, S1 and S2 | as U11/U12 but RUV `skip` (the bindings' corpus setting; 5 proteins) |

Environment for both builds: `MULTISCHOLAR_TEST_MODE=true`, `RAYON_NUM_THREADS=1`, other thread settings unset
(recorded); `RSTUDIO_PANDOC` = the bundle's pandoc.

## 4. Corpora

Generated on the runner by the candidate's generators, compared with the bindings' SHA-256: `lipid_screen`,
`metab_screen`, `trx_screen` (`kit_omics_corpus.R`, default recipes), `screen` and `full` (the bound recipe through
`writePerf107SyntheticCohort`, the encoder from `tests/testthat/helper-workflow-checkpoint-fixtures.R`, sourced, not
edited), `lfq_fixture`, `tmt_fixture` (`writeKitFixtureCorpus`: copies of `tests/testdata/e2e/prot_lfq` and `prot_tmt`,
synthetic by `tests/testdata/e2e/README.md`, "All data is synthetic"), `metab_large_missing`, `lfq_300`, `tmt_300` (the bound
recipes' arguments).

## 5. What leaves each job

Only the result folder is uploaded (artifact `k209-<job>`, 30 days): `runner.txt` (host facts, the candidate's commit,
the lock file's sha256, pandoc's version), `builds-installed.tsv`, `library-b2.tsv`, `library-versions.tsv` (package
names and versions of the system library), `corpora.tsv` (file names, sizes, sha256), `script-sha256.txt`, `plan.txt`,
`plan-console.log`, `injected-limits.json` (G0a), `journeys/<id>.{baseline,marks,driver}.json` and `<id>.json` (the
figures, settings, build facts, error messages and frame names, fingerprints' sha256), `logs/<id>.log` (each journey's
console), `journeys.tsv`, `steps.tsv`, `g0a.tsv`. Before the upload the job fails on token-like text, and a console log
or the plan log naming an internal package is withheld (listed in `withheld.txt`). The libraries, the build sources and
logs, the corpora and the stores stay on the runner and are not uploaded or cached.

## 6. G0 runs, declared before each is sent

Instrument failures follow VAL-54: at most three further runs per item, each recorded in `runs.tsv` and `outcome.md`
with what failed and the fix, no expected value or limit changed, inside the G0 re-run row (540 min). A measured time
or peak is a result and is never re-run.

| Run | Purpose | Binds | Jobs, ceiling | What it decides |
|---|---|---|---|---|
| G0a | lane lifecycle; reset control; gate B B1 to B6 at screen size (U1 to U5, S1, S2, base then candidate); DIA body validation; settings | `C`, `base` | 1 job, 150 min | items 1 to 5 below; whether the lane may measure (stop rule 1 and 2) |
| G0p1 | time probes and B7: U10 (`fastlimpadpc`) and U9 on `C`, one journey each, to the job deadline | `C` | 1 job, 240 min | U9 and U10's measured times and peaks: their layouts (9.2) and limits (9.3) |
| G0p2 | time probes and B7: U7, U11, U12 on `C`, one journey each, to the job deadline | `C` | 1 job, 240 min | the same for U7, U11, U12 (R13: whether U7 fits one job) |
| G0i | integration lifecycle on KIT-150's lane with `C` as a row: the S pipeline on `C` and on build C (324e14b4), then a memory-limit failure on `C` | `C`, 324e14b4 | 1 job, 60 min | whether KIT-150's lane runs the head's integration functions; the failure recorded after the baseline |

**G0a items:**
1. The lifecycle journey `U1_C_r1` (lipid_screen, `C`) completes with every mark through the last step; then, relative
   to it in the same job: memory (`memory.max` = 0.7 × its journey `VmHWM`; the kill must come after the baseline mark,
   class memory), time (limit = 0.5 × its process wall time, rounded up; class time), package error (the bound
   `norm_method` set to `k209_invalid_method`; class package error).
2. Reset control `CAL_C_r1`: a calibration step allocates, fills and frees 512 MiB, then a no-op step. Holds if the
   calibration step's `VmHWM` ≥ its entry `VmRSS` + 480 MiB and the no-op step's `VmHWM` ≤ its entry `VmRSS` + 32 MiB.
3. Gate B: B1 (both builds install and load from their own libraries), B2 (`library-b2.tsv`), B3 (each route's screen
   journey completes on the base and on `C`), B4 (the bound settings recorded in every journey record, compared base
   against candidate), B5 (`g0.md`'s table: the lane's steps against each build's GUI route), B6 (no
   `WorkflowExecutionApplication`; the journey's cgroup holds the R process only, pandoc during the report).
4. DIA body: section 2.4's list; the base-to-candidate state fingerprints of U4 reported.
5. Settings: the fixtures' provenance (section 4); both builds' manual-RUV tab defaults (recorded in every protein-route
   journey's build facts).

Screen journeys in G0a run under a 1,800 s limit each (no measured time exists yet; G0a's measured times then set
the screen limits, ticket 9.3). G0p's journeys run to the job's deadline (12 minutes before its ceiling).

## 7. The counted runs (declared in the ticket; not sent by this declaration)

Repeats (ticket 7.1, 7.2): screen 3 runs per build, size 2 per build in the order B C then C B, a third pair (B C) on a
size unit when either build's range of its 2 runs exceeds the base-to-candidate difference of the medians on the
journey maximum or the retained-heap maximum, or when the unit is listed as a finding. Listing aid (7.3): a figure is
listed when the difference of the medians exceeds both max(0.25 GiB, 10 % of the base's median) and twice the larger of
the two builds' ranges; also every step where the committed state differs and every non-completion; an aid, not a
verdict. Re-measure map: ticket section 9.4. Time limits (9.3): 3 × the staging head's G0p time (G0a's for screen
units), cut to the job's remaining time less 12 minutes.
