# KIT-097 audit (audit_only, 2026-11-06)

**Target:** `tickets/kit/KIT-097.md` rev 1.0, status `pending`, `needs_planning: true`, dependency `KIT-096`.
**Mode:** `audit_only`. Nothing was repaired, edited, run or committed apart from this file. No tests, installs, journeys or app launches were run. Nothing under `kit/runtime` and no KC data was read. Two read-only `Rscript` probes were run: `args(Glimma:::glimmaVolcano.DGELRT)` and `requireNamespace()` for edgeR, RUVSeq, EDASeq and ruv.
**Baseline:** branch `performance-v2` at `bbe1fae9`. The tree is dirty only with four untracked `kit-090/*.txt` files, which belong to KIT-090.
**Authorities read:**
- `AGENTS.md`;
- `conventions/R.md`, `R-shiny.md`, `R-golem.md` and `R-bioinformatics.md` (their normative sections);
- `tickets/kit/00-KIT-PLAN.md` v3.8 (rows 231–238, decision 2026-11-06);
- the approved plan `~/.claude/plans/write-a-series-of-atomic-dragon.md` ("KIT-096 and KIT-097: the GUI");
- the notebook `dev/distributable/plans/transcriptomics_workflow.rmd` (lines 160–793);
- neighbours KIT-040 and KIT-093 to KIT-099, plus KIT-060, KIT-061 and KIT-082 for consumer contracts.

**Validator:** `next-openai-ticket.sh KIT-097` passes: valid, 0 errors, 0 warnings, status consistent, planning explicit. The index dependency (`KIT-096`) matches the frontmatter, and KIT-098 depends on KIT-097, so the graph is consistent.

**Framing:** the ticket is planning-gated. Missing detail alone does not block it, because a file-by-file plan is due when it starts. A finding is raised where:
- an acceptance criterion could pass while the scientific claim is false;
- the ticket's own reuse instruction leads into a defect;
- the file allowlist excludes a file that the outcome needs. Under AGENTS.md that is a scope decision, not an implementer's choice.

## Finding register

| ID | Severity | Dimension | Blocked state |
|---|---|---|---|
| KIT097-F01 | blocking | 2 Scientific rigor (producer-consumer) | implementation, promotion |
| KIT097-F02 | blocking | 1 Completeness / 2 Rigor (state order, record) | implementation |
| KIT097-F03 | should-fix | 2 Scientific rigor | implementation |
| KIT097-F04 | should-fix | 2 Scientific rigor (record) | implementation, decision |
| KIT097-F05 | should-fix | 2 Scientific rigor (report wording) | promotion |
| KIT097-F06 | should-fix | 1 Completeness (output map) | implementation |
| KIT097-F07 | should-fix | 1 Completeness (results summary to report) | implementation |
| KIT097-F08 | should-fix | 1 Completeness (MOFA export) | implementation, decision |
| KIT097-F09 | should-fix | 3 Best practice (reuse versus copies) | implementation |
| KIT097-F10 | should-fix | 3 Best practice (scope / allowlist) | implementation |
| KIT097-F11 | should-fix | 3 Verification (no-matrix rule) | execution |
| KIT097-F12 | minor | 3 Test-first | implementation |
| KIT097-F13 | minor | 3 Budget / environment | execution |
| KIT097-F14 | minor | 2 Rigor (volcano semantics) | implementation |
| KIT097-F15 | minor | 2 Rigor (contrast scope) | execution |
| KIT097-F16 | minor | 1 Completeness (harness hand-off) | none (KIT-098 planning) |
| KIT097-F17 | minor | 1 Completeness (technical-replicate M availability) | implementation |
| KIT097-F18 | minor | 4 Reference accuracy | none |

---

### KIT097-F01 (blocking): reusing the omics DA path would make the record and report say limma for an edgeR QL F-test

- **Evidence:**
  - The ticket tells DA to reuse "the omics DA display helpers" (`KIT-097.md` outcome). The shared managed DA is `runWorkflowOmicsManagedDa` (`R/utils_workflow_omics_actions.R:338-371`). It accepts the manifest through `acceptWorkflowOmicsDaManifest` (`R/utils_workflow_omics_declarations.R:1078-1100`), which hard-codes `model = list(model = "limma", ..., eBayes_trend = TRUE, eBayes_robust = TRUE, treat_applied = isTRUE(parameters$lfc > 0))` (lines 1088-1091; the same in `omicsDaOperation`, lines 366-369). Its parameters are `settings$treat_lfc_cutoff` (line 1082).
  - The record copies `model = manifest$model` (`R/utils_workflow_analysis_record.R:308`).
  - The report overlay reads `da_model` (line 806) and the TREAT cutoff as `s4_de_treat_lfc_cutoff` / `ui_de_log_fold_change_cutoff` (lines 755, 788).
  - The notebook fits `estimateDisp(..., robust = TRUE)` and `glmQLFit(..., robust = TRUE)`, then tests with `glmQLFTest` and adjusts with BH through `topTags` (rmd 535-554). It applies no TREAT and no eBayes.
- **Impact:**
  - If the helper is reused unchanged, the analysis record and any report built from it misstate the test that ran.
  - No KIT-097 criterion checks this. AC1 checks only the RUVIII `M`, so the criteria pass while the claim is false.
  - KIT-095 owns `utils_workflow_omics_declarations.R` and could make the model per binding, but neither ticket requires it.
- **Repair:**
  - Add a KIT-097 criterion: the accepted manifest and the record's `analyses$differential_abundance$model` name edgeR, `glmQLFit(robust = TRUE)`, `glmQLFTest`, BH FDR, the realised formula with the `W_` terms, and `treat_applied = FALSE`.
  - The DA tab offers no TREAT/lfc test input. If it keeps an lfc field, label it display-only and record it as such.
  - Record in the planning step that KIT-095's binding supplies the model (a per-binding `model` field). Otherwise KIT-097 needs that file in scope.
- **Closure evidence:** a module CI or in-process test reads `analysis_record.json` after a trx DA and asserts those fields. A grep of the rendered synthetic report finds no `limma`, `eBayes` or `TREAT` in the trx methods text.

### KIT097-F02 (blocking): where RUVIII runs, relative to the DA head, is undecided, and that decides whether AC1 can be met

- **Evidence:**
  - The ticket puts "the RUVIII `M` choice" in the normalisation tab (`KIT-097.md`, Normalisation).
  - The notebook runs RUVIII after DA, on `log_cpm_tmm` with RUVs' controls and `k`, only to export for MOFA (rmd 698-724). DA runs on counts plus W from RUVs (rmd 522-536).
  - The managed DA refuses unless its source is the head (`utils_workflow_omics_actions.R:340-342`).
  - Record steps come from committed states. Publication products such as `*_integration_export` commit nothing (`utils_workflow_omics_declarations.R:1297-1306`), and neither does `omicsManagedProduct` (`utils_workflow_omics_actions.R:179`).
  - KIT-095 lists `trx_ruviii_export` as a worker step with `M` as a parameter, but does not say whether it commits a state. Its store steps are "QC finalise and reset", with no `norm_complete`; the lipid tab commits `lipid_norm_complete` before DA (`R/mod_lipid_norm_workflow_helpers.R:637-639`).
- **Impact:** both outcomes fail:
  - If RUVIII commits a state in the normalisation tab, it becomes the head: DA is then refused, or runs on RUVIII-corrected log values, which is not the notebook's analysis.
  - If RUVIII is a product only, the record has no step carrying `M`, so AC1 ("the record names the RUVIII `M` choice") cannot be met.
- **Repair:**
  - Declare the state sequence in the ticket: `trx_tmm_logcpm`, then `trx_ruv_grid`, then `trx_ruvs`, then (a norm-complete marker or none), then DA on the RUVs state, then `trx_ruviii_export`.
  - Say where the `M` control sits. The summary/export stage matches the notebook's order.
  - Say how `M`, `k` and the control set enter the record: as a committed step's parameters, or as outcomes written by KIT-095.
- **Closure evidence:** an in-process test after a full trx journey asserts:
  - the record's DA `source_state` is the RUVs state;
  - a RUVIII step or outcome carries `M = group` (default), and a second run `M = technical_replicate`;
  - the DA results are unchanged by the RUVIII choice.

### KIT097-F03 (should-fix): the significance column's meaning is unfixed, and the lipid/metab helpers use a Storey q-value

- **Evidence:**
  - Lipid and metab compute `fdr_qvalue` with `qvalue::qvalue`, falling back to BH, and keep BH separately as `fdr_value_bh` (`R/func_lipid_da_model.R:136-162`; the same in `R/func_metab_da_core_engine.R:141-157`).
  - The reusable display filters and plots on `fdr_qvalue` and `raw_pvalue`:
    - the table: `R/mod_lipid_da_results_helpers.R:110-131`;
    - the volcano's `-log10(fdr_qvalue)`, labelled `negLog10FDR`: `R/func_lipid_da_plotting.R:197-205, 309-313`.
  - The notebook's significance is edgeR's BH `FDR` at 0.05 (rmd 554, 571-590, 772).
- **Impact:** reuse forces one of two errors:
  - computing a Storey q-value on edgeR p-values, which changes the DE calls relative to the reference;
  - relabelling BH FDR as a q-value.

  KIT-098 compares final states to the notebook, but the tab's labels and counts are not covered.
- **Repair:**
  - State that the trx significance column is edgeR `topTags` BH `FDR`, identical to the reference, and that tables, volcano and counts label it "FDR (BH)".
  - If helper column names force a mapping, record it and do not compute a q-value.
- **Closure evidence:** on the screen corpus, the DA tab's significant counts at FDR < 0.05 (up and down) equal the reference's, and the reference's `FDR` vector is `identical` to the tab's column.

### KIT097-F04 (should-fix): the RUV `k` and percentage selection mode is neither bounded nor recorded

- **Evidence:**
  - The notebook selects both automatically: percentages 5 to 30 by 5, a `k` cap of 3 when there are fewer than 15 samples and 4 otherwise, penalty 0.05, maximum composite score (rmd 239-302).
  - The lipid tab offers automatic, manual or skip, with a user range and a free `k` (`R/mod_lipid_norm_ui_helpers.R:92-210`).
  - The ticket says "k and percentage selection", but AC1 records only `M`.
  - Today the record's RUV outcomes are lipid/metab/proteomics specific (`R/utils_workflow_analysis_record.R:118-160`). KIT-095's trx outcomes are not yet specified.
- **Impact:**
  - A manual override, or a skip (DA becomes `~ 0 + group` with no W), could run and leave no trace of which choice was made.
  - The report's RUV sentence would then be unverifiable.
  - A manual `k` above the cap departs from the reference without warning.
- **Repair:** add criteria that the record names:
  - the selection mode (automatic or manual);
  - the grid's settings, and the chosen `k` and percentage;
  - the number of controls, and the separation and composite scores;
  - "skip", if offered, which the report must then state as "no RUV covariates".

  Also: the default reproduces the reference, and a manual `k` is bounded by the cap or warned.
- **Closure evidence:** a test of an automatic run and a manual run that asserts the record fields. The automatic run's `k`, percentage and control set are `identical` to the KIT-040 reference.

### KIT097-F05 (should-fix): the template scaffold carries claims this analysis does not make

- **Evidence:**
  - The scaffold candidates are the per-omic templates (1,100 to 1,275 lines). The lipid methods text cites limma (`inst/reports/lipidomics/lipidomics_report.rmd:677`) and RUVIII-C correction of the analysed data (lines 600-656), and parses TREAT and q-value settings (lines 466-477).
  - What the transcriptomics analysis actually does:
    - DA is an edgeR QL F-test with the RUVs W as model covariates (rmd 527-554);
    - the `removeBatchEffect` logCPM is for display only (rmd 420-426);
    - RUVIII is for the MOFA export only (rmd 698-724);
    - the controls are inferred: the least significant genes in a preliminary `~ group` QL GLM, coefficient 2 (rmd 231-236).
  - The notebook's summary hard-codes "Dataset: GSE152840 (Sepsis Transcriptomics Workshop)" (rmd 759).
  - The lipid UI's help text says "most stable lipids used as negative controls" (`mod_lipid_norm_ui_helpers.R:194`).
- **Impact:**
  - A copied scaffold claims that DA ran on RUV-corrected data with limma, TREAT or q-values, and that the controls are "stable" rather than inferred. This conflicts with the series contract: "Inferred control genes are labelled inferred."
  - No criterion checks the report's wording.
- **Repair:** add a wording criterion with required statements:
  - edgeR QL F-test, robust;
  - the realised formula with `W_1..W_k`;
  - BH FDR threshold;
  - RUVs W as covariates;
  - display-only corrected logCPM;
  - RUVIII for integration with its `M`;
  - inferred empirical controls, and how they were chosen.

  Forbidden statements: limma, eBayes, TREAT, q-value, RUVIII-C, "stable/housekeeping controls", any dataset name.
- **Closure evidence:** a grep of the synthetic e2e report for both lists, kept as evidence in `kit-097/`.

### KIT097-F06 (should-fix): no map from notebook outputs to the GUI, results summary or report

- **Evidence:** the notebook's outputs (rmd 160-793):
  - the filtering summary;
  - the grid table and `canonical_correlation_optimization_summary.xlsx`;
  - CCA panels A and B (`canonical_correlation_panels`);
  - `ruv_normalised_log_values.tsv`;
  - `normalization_comparison_grid`: PCA and RLE at raw, TMM and RUVs;
  - `de_results_<contrast>.csv` and `top_de_genes.tsv`;
  - the volcano and the `NumSigDeMolecules` bar plot;
  - the Glimma HTML;
  - the RUVIII TSV and the `ruviii_before_and_after` grid;
  - `ruv_optimization_summary.txt` and `sessionInfo`.

  KIT-095's publication list ("QC figures at three stages, CCA diagnostics, composite, DA publication, final PCA, volcano page, session and integration export") does not cover:
  - the RUVIII before/after grid;
  - the significance bar plot;
  - the top-genes table;
  - the optimisation workbook and summary text.

  The ticket names none of these.
- **Impact:** "every notebook output surfaced" cannot be checked, and the e2e scope is undefined.
- **Repair:** at planning, add a table giving each notebook output's destination: tab display, results-summary path, report section, or "dropped (reason)". Where a needed publication operation is missing, return it to KIT-095.
- **Closure evidence:** the e2e test asserts each mapped output exists, and the map is kept in `kit-097/`.

### KIT097-F07 (should-fix): the results summary has no transcriptomics branch, and the report reads the results summary

- **Evidence:**
  - `copyToResultsSummary` branches only for proteomics, metabolomics and lipidomics (`R/func_general_filemgmt_results_summary_helpers.R:325, 372, 408, 523`). Transcriptomics gets only the common set: contrasts, design, `study_parameters.txt`, the generic `NumSigDaMolecules` tab (line 361-368) and DA workbooks matching `(da|de)_.+_long_annot\.xlsx`.
  - The notebook writes `NumSigDeMolecules/`. Only the lipid branch reads that name (line 585).
  - The lipid template reads its figures and tables from `results_summary_dir` (`lipidomics_report.rmd:71-78, 1192-1261`).
  - `RenderReport` already has the `transcriptomics` template leaf (`R/func_general_filemgmt_report_render_helpers.R:131-142`).
  - The results-summary helper is not in `files_to_modify`.
- **Impact:** AC1's "renders the report" can pass with a report that has no QC, CCA, volcano or table content.
- **Repair:** add `R/func_general_filemgmt_results_summary_helpers.R` to the allowlist with a transcriptomics branch fed by F06's map, or record a decision that the generic path suffices and list what the report then lacks.
- **Closure evidence:** the e2e test asserts the copied file set in `results_summary_dir` and that the rendered report contains those figures.

### KIT097-F08 (should-fix): the MOFA integration export contract is undefined

- **Evidence:** four sources disagree, and KIT-061 is pending and depends on KIT-060:
  - The notebook writes RUVIII `newY` (samples by genes, with a `Run` column) as TSV to `results_dir/transcriptomics_normalized_data.tsv` and to `base_dir/results/Transcriptomics/...` (rmd 714-724).
  - The integration workbook reads `<integration>$data_dir/transcriptomics_normalized_data.tsv` (`inst/workbooks/integration/integration_workflow_starter.rmd:338`).
  - The shared `*_integration_export` writes the head state's S4 object as RDS (`utils_workflow_omics_declarations.R:1297-1306`; lipid caller `R/mod_lipid_summary_workflow_helpers.R:117-120`).
  - KIT-061 owns the declared per-omics integration input.
- **Impact:** the summary's export could write the RUVs state as RDS, not the RUVIII matrix, which silently changes what MOFA would integrate.
- **Repair:** state the exported product (RUVIII `newY`, with `M` recorded), its format and its path, and that KIT-061 will consume it. Alternatively, defer explicitly to KIT-061 and leave the summary without an integration export until then.
- **Closure evidence:** the e2e test asserts the file, and its matrix is `identical` to the reference `newY`.

### KIT097-F09 (should-fix): the reuse is not made concrete, and some reuse conflicts with "other omics' modules unchanged"

- **Evidence:**
  - **Reusable as is:**
    - `installWorkflowOmicsActions`, `omicsManagedCommit` and `omicsManagedProduct`, `workflowOmicsDisplay` (`utils_workflow_omics_actions.R:13, 101, 172, 179`);
    - `workflowOmicsVolcanoPageFactory`, given a `binding$volcano_widget` (line 143);
    - `reportWorkflowOmicsDaExport` (line 373);
    - `exportWorkflowAnalysisRecordBundle` (`utils_workflow_analysis_record.R:722`);
    - `RenderReport`.
  - **Not reusable:**
    - `runWorkflowOmicsManagedNormalisation` hard-codes ITSD, log2, between-sample normalisation, RUV-III and the composite (lines 222-313).
    - `handleLipidSummaryGenerateReport` takes `templateFilename` but hard-codes `system.file("reports", "lipidomics", ...)` and the GitHub `.../reports/lipidomics/` URL (`R/mod_lipid_summary_workflow_helpers.R:401, 412`).
    - The lipid results table keeps only `lipid_id` and `lipid_name` as identifiers (`mod_lipid_da_results_helpers.R:125-130`), so a gene table would lose its id column.
    - The lipid Glimma builder reads `lipid_name` (`func_lipid_da_plotting.R:213`).
  - No omic-neutral summary helper exists. The lipid summary (1,517 lines across 4 files) and the metab summary (1,360 lines across 3 files) are parallel copies.
- **Impact:** an implementer must either edit lipid files (against the backward-compatibility criterion) or copy large blocks (against "reuse, not copies" and AGENTS.md's smallest-sufficient rule). The ticket does not decide which.
- **Repair:**
  - At planning, write a reuse map listing what is reused as is, what is parametrised (with the lipid/metab regression tests that prove them unchanged) and what is copied thin and trx-local.
  - Write the trx normalisation sequence as a small runner on `omicsManagedCommit`, not a copy of the lipid one.
  - Bring any parametrised shared file into the allowlist.
- **Closure evidence:** the reuse map in `kit-097/`; `git diff --stat` against it; the lipid and metab module tests and module CI green.

### KIT097-F10 (should-fix): `files_to_modify` omits files the outcome needs

- **Evidence:** the outcome needs files that the allowlist does not name:
  - `R/mod_transcriptomics.R` holds the tabset (KIT-096 owns it). The tabs must be registered there, as in `R/mod_lipidomics.R:53-117`.
  - `tests/testthat/helper-module-ci-trx-*.R` and `tests/testdata/mci/manifest.json` (KIT-096 has both; lipid has per-tab helpers such as `helper-module-ci-lipid-norm.R`).
  - `tests/testdata/e2e/manifest.json`, a trx fixture directory and `tests/testdata/e2e/report_templates/transcriptomics_report.rmd`. The e2e harness stubs templates from that directory (`tests/testthat/helper-e2e-report-doubles.R:1-43`).
  - `tests/testthat/test-workflow-omics-application-reads.R` (see F11).
  - `R/func_general_filemgmt_results_summary_helpers.R` (see F07).
  - Possibly shared helpers (see F01, F09).
- **Impact:** under AGENTS.md, any of these is a mid-ticket scope expansion that needs a decision.
- **Repair:** extend the allowlist at planning, with a reason for each file.
- **Closure evidence:** the revised ticket frontmatter, and the diff confined to it.

### KIT097-F11 (should-fix): the criterion "No whole matrix held in the application" names no probe

- **Evidence:**
  - KIT-082 established the probe with a positive control: `tests/testthat/test-workflow-omics-application-reads.R` (header lines 1-3; `omicsReadCase` covers lipid and metab only).
  - The lipid normalisation tab hydrates the RUV optimisation result into the app (`utils_workflow_omics_actions.R:268`).
  - The trx grid result carries gene-length `ctl_idx` and `control_genes` vectors plus the `cancor` plot data (rmd 280-290). Neither is a matrix, but each should be bounded and stated.
- **Impact:** the criterion is unfalsifiable as written.
- **Repair:**
  - Extend the KIT-082 probe with a transcriptomics case covering the normalisation displays (CCA, PCA/RLE), the DA tab (table, volcano, download) and the summary, with its positive control.
  - List the non-matrix values the app does hold.
- **Closure evidence:** the probe fails on a deliberately whole-matrix read and passes on the trx tabs.

### KIT097-F12 (minor): the test-first criterion is trivially satisfied

- **Evidence:** "Module CI tests fail before the tabs exist." Any test of a missing function fails.
- **Impact:** the criterion does not falsify behaviour.
- **Repair:** name the behaviour tests (the record's `M`, the DA model fields, the `k`/percentage mode, the results-summary set). Each must fail before the change for the reason it asserts.
- **Closure evidence:** a fail-before log showing the assertion message, not "could not find function".

### KIT097-F13 (minor): the e2e has no ceiling, and its packages are absent

- **Evidence:**
  - The budget is "Module CI, one e2e test; one bounded correction pass", with no runtime ceiling or attempt count (AGENTS.md, Verification).
  - A read-only probe shows RUVSeq and EDASeq are not installed (`requireNamespace` FALSE), while edgeR and ruv are. KIT-093 adds them to `Suggests`.
- **Impact:** a `skip_if_not_installed` could make the module CI and the e2e pass vacuously.
- **Repair:** declare the e2e's ceiling and attempts, and that a skip counts as a failure for acceptance. Install the packages before the counted run.
- **Closure evidence:** the run log shows no skips and the e2e's wall time under its ceiling.

### KIT097-F14 (minor): the volcano's significance and axis are not declared

- **Evidence:**
  - The notebook's static volcano plots `-log10(FDR)` with significance at FDR < 0.05; its ±1 logFC lines are visual only (rmd 596-640).
  - Its Glimma call passes `lfc = 1, fdr = 0.05` (rmd 680-690). These are not formal arguments of Glimma 2.15.2's `glimmaVolcano.DGELRT`, whose `status` defaults to `decideTests` and whose y axis is `negLog10PValue` (read-only `args()` probe), so they fall into `...`.
  - The lipid interactive volcano uses `glimmaXY` with `-log10(fdr_qvalue)` and |logFC| ≥ 1 labels (`func_lipid_da_plotting.R:197-215, 309-321`).
- **Impact:** the interactive and static plots may imply an lfc cutoff that no test applied.
- **Repair:** state the trx volcano's y axis, its status rule (FDR (BH) < threshold, no lfc test) and the label wording.
- **Closure evidence:** a module CI assertion on the widget's status vector against the FDR rule.

### KIT097-F15 (minor): parity covers one contrast, the GUI offers several

- **Evidence:**
  - The notebook evaluates only the first contrast line of `contrast_strings.tab` (rmd 539-545), and its preliminary GLM uses coefficient 2 only (rmd 234).
  - The GUI's "contrasts" go through `omicsDaContrastSteps` per row (`utils_workflow_omics_declarations.R:1065-1071`).
- **Impact:** KIT-098's "identical to the notebook" is defined only for the first contrast.
- **Repair:** state that reference parity is for the first contrast, that the others come from the same `glmQLFit`, and that the record lists them all.
- **Closure evidence:** the record lists each contrast, and the first contrast's table is `identical` to the reference.

### KIT097-F16 (minor): no hand-off of test ids to KIT-098's driver

- **Evidence:** KIT-098 adds a driver step function and the journey steps normalisation, da, session_export and report (plan, KIT-098). Lipid tabs expose `data-testid` (`R/mod_lipidomics.R:65-115`). KIT-097 does not name the ids or buttons its tabs expose.
- **Impact:** KIT-098 may need to edit KIT-097's modules.
- **Repair:** list the stable `data-testid`s and action ids in KIT-097's outcome.
- **Closure evidence:** the ids are listed in `kit-097/outcome.md`, and KIT-098's driver uses them without module edits.

### KIT097-F17 (minor): the rule for the technical-replicate `M` option is unspecified

- **Evidence:** decision 2 makes technical replicates "selectable". KIT-093's class has `technical_replicate_id`. The ticket does not say what the GUI does when the design has no technical replicates; `M` would then be the identity, or the kernel refuses.
- **Impact:** the GUI could offer an option that does nothing or fails without a clear message.
- **Repair:** enable the option only when replicates exist, or show KIT-093's classed refusal; record the attempt.
- **Closure evidence:** a module CI case on a design without technical replicates.

### KIT097-F18 (minor): two of the ticket's phrases name no existing shared code

- **Evidence:**
  - "the omics DA display helpers" names no shared set. The DA display helpers are per omic (`R/mod_lipid_da_display_helpers.R`, `R/mod_metab_da_display_helpers.R`). Only `workflowOmicsDisplay` and `workflowOmicsVolcanoPageFactory` are omic-neutral.
  - "the existing scaffold" could mean the per-omic templates or the shared `inst/reports/shared/analysis_record_SI.rmd` (95 lines).
  - The other citations check out:
    - `omicsManagedCommit` (`utils_workflow_omics_actions.R:172`);
    - the runner-or-wrapper description (lines 1-4);
    - `mod_lipidomics.R:180-196`;
    - `.MODULE_CI_VALID_OMICS` (`tests/testthat/helper-module-ci-manifest.R:1`);
    - `RenderReport`'s transcriptomics leaf;
    - `setupDirectories`' transcriptomics layout (`R/func_general_filemgmt_setupdirectories_helpers.R:78-93`);
    - the check baseline of 3 WARNINGs and 3 NOTEs (`kit-082/check.log`, `kit-086/check.log`).
- **Impact:** the reuse claim cannot be verified as written.
- **Repair:** name the helpers and the scaffold file in the planning revision (see F09).
- **Closure evidence:** the revised ticket text.

---

## Verdict per dimension

1. **Completeness: not ready.** Three gaps:
   - there is no map from notebook outputs to GUI, summary and report (F06);
   - the results summary has no transcriptomics branch, and the report depends on it (F07);
   - the MOFA export product and path are undefined (F08).

   Where RUVIII runs relative to the DA head is undecided, and AC1 depends on it (F02). The e2e scope follows from F06. RenderReport, the directory layout and the record bundle are in place.
2. **Scientific rigor: not ready.**
   - The ticket's reuse direction would make the record and report claim limma for an edgeR QL F-test (F01, blocking).
   - The significance column could become a Storey q-value (F03).
   - The `k`/percentage mode is unrecorded (F04).
   - No criterion checks the report wording against the inferred-controls rule and the RUVs/RUVIII split (F05).

   AC1 records `M` only.
3. **Implementation best practice: should-fix.**
   - Following the existing idiom (managed runner, references in the tabs, golem `mod_` files, module CI plus e2e) is right.
   - The reuse is not made concrete, and some lipid helpers cannot be reused without editing lipid files (F09).
   - The allowlist omits needed files (F10).
   - The no-matrix criterion names no probe (F11).
   - Test-first is trivial (F12), and the e2e budget has no ceiling (F13).
4. **Codebase reference accuracy: mostly accurate.** The cited functions and files exist as described. Two phrases ("omics DA display helpers", "existing scaffold") name no single existing artefact (F18). The undeclared dependencies on `acceptWorkflowOmicsDaManifest`'s limma model and `copyToResultsSummary`'s branches are the material omissions.

**States:** implementation_verified n/a; execution_complete no; decision_supported no (F02, F04 and F08 need decisions); promotion_authorized no. Claim scope: the ticket as written, at rev 1.0.
