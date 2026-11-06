# KIT-097 outcome: the transcriptomics normalisation, DA and Summary & Export tabs (interim GUI), 2026-11-06

Ticket revision 1.1, then **resumed and narrowed** (the user, 2026-11-06: "prove transcriptomics works first"). The tabs are
the **interim GUI**: the lipid norm, DA and summary modules under KIT-096's `profile` hooks. The GUI workover (KIT-101 to
109) replaces them with one shared GUI layer; nothing was refactored or extended beyond the e2e fix below. Plan:
[plan.md](plan.md), written before any code.

**Result: transcriptomics works through the GUI.** The one e2e journey passes on its correction attempt (14 pass, 0 fail,
0 skip, 113 s). `R CMD check` is at the baseline (3 WARNINGs, 3 NOTEs). `kit_097` = `$P/core/browser-smoke.5gh5pn`, from
`35f855bd`.

## Commits

| Commit | What |
|---|---|
| `efb8b605` | the plan |
| `55a53d44` | the tabs, tests and template |
| `8e97b55e` | coordinator review: the rename `trxDaContrasts` -> `trxDaContrastEntries` |
| `fcf1093a` | coordinator: the stop record (kept below) |
| `35f855bd` | the e2e defect: the interactive volcano is built on demand |
| this commit | evidence and this outcome |

## What was built (interim GUI)

- **Normalisation tab** (`R/mod_trx_norm.R`; `mod_lipid_norm_ui/_server(profile =)`): TMM/logCPM (`trx_tmm_logcpm`), the
  RUV grid (`trx_ruv_grid`, `automatic` = the notebook's grid, or `manual` with `k` within the cap; no skip), RUVs
  (`trx_ruvs`, the DA head), then the PCA and RLE at raw, TMM and RUVs, the CCA panels with the grid workbook, and the
  stage grid. It runs on the managed runner (one action each) or in process (the KIT-095 wrappers, plus the same publication
  operations run in process). The lipid tab's RUV QC, PCA, RLE and Log panels and its renderers are reused. The help text
  says "inferred negative controls (the least significant genes in a preliminary group GLM)".
- **DA tab** (`R/mod_trx_da.R`; `mod_lipid_da_ui/_server(profile =)`): the lipid tab's outputs and observers, with the
  model text and the realised formula in place of the formula input. The analysis is edgeR `glmQLFit(robust = TRUE)` and
  `glmQLFTest`, with BH FDR labelled "FDR (BH)". The fold-change field is display-only. There is no TREAT, no heatmap and no
  session load. Every contrast line runs pairwise on one fit. The static volcano is KIT-093's (−log10 FDR (BH)). The
  interactive Glimma volcano is the binding's widget in process, or the worker's page when managed, built on demand.
- **Summary & Export** (`R/mod_trx_summary.R`; `mod_lipid_summary_ui/_server(profile =)`):
  - the RUVIII export after DA (`trx_ruviii_export`, accepted on the head), with `M` = group, or technical replicates
    offered only when the design has them;
  - its before/after figure;
  - the integration export: KIT-095's RUVIII TSV, `transcriptomics_normalized_data.tsv` in `integration_dir`;
  - the record, the results summary and the report from `reports/transcriptomics`.
- **Results summary**: `copyToResultsSummary`'s transcriptomics branch, from the output map, and `DA_results_transcriptomics.xlsx`
  built from the per-contrast `de_results_*.csv`.
- **Report** (`inst/reports/transcriptomics/transcriptomics_report.rmd`): the lipid template's structure (params, paths,
  figures and tables from `results_summary_dir`, the `workflowRecordReportParameters` overlay), with transcriptomics text.
- **Wiring** (`R/mod_transcriptomics.R`): the six tabs, the stepper, and the tab-status chain.

## Defects found and fixed

1. **Name clash (coordinator review, `8e97b55e`).** `trxDaContrasts(contrastsTbl)` in `mod_trx_da.R` collated after KIT-093's
   kernel of the same name and replaced it in the namespace. KIT-093's `trxDaContrast()` call broke, and `R CMD check` showed a
   fourth WARNING. Renamed `trxDaContrastEntries`. A scan of `R/` finds no other function defined twice. After the fix, the
   trx kernel, design, publication, binding and module-CI files pass ([tests.txt](tests.txt)).
2. **The e2e defect (`35f855bd`).** The first counted e2e attempt failed after 1,264 s: `Timed out waiting for
   transcriptomics/differential_analysis to become complete; last status was <NULL>`
   ([e2e-failed-attempt1.txt](e2e-failed-attempt1.txt)). A diagnostic run, not counted, on trx_screen
   ([e2e-diagnosis.txt](e2e-diagnosis.txt)) showed:
   - DA completed about 4 s after the click: `da_status` read "Analysis Complete";
   - the tab status was never lost;
   - `de-volcano_glimma` then arrived, and every later wait for Shiny to be stable failed.

   The cause was the DA finaliser selecting the single assay (`singleAssay = TRUE`) while the "Interactive Plot" checkbox
   defaulted on, so the Glimma volcano of every gene was built and sent on completion.

   **Fix:** under the profile the checkbox is off by default (`da_volcano_interactive = FALSE`), and the profile's render
   requires it, so the static volcano shows and the interactive one is built only when asked for. The lipid HTML and
   behaviour are unchanged.

   **Test first:** `MCI-029.2` (`test-module-ci-trx-app.R`, from the stopped reproduction patch) runs the module in process
   through DA, with the checkbox at the UI's declared default. It asserts that every status is kept and read by the digest,
   and that DA is `complete`. Before the fix it failed on `widgets_on_completion` (1, expected 0); after it, 0 on completion
   and 1 on demand.
   - The e2e's correction attempt then **passed** ([e2e-attempt2.txt](e2e-attempt2.txt)).
   - **Residual risk (not fixed, out of the narrowed scope):** ticking "Interactive Plot" in process on a large gene set
     sends the whole widget to the page; the managed path uses an iframe page. A matter for the GUI workover.

## Evidence

- **Test first** ([test-first.txt](test-first.txt), on a baseline worktree): the e2e helper's transcriptomics id failed on
  its assertion. With the tabs unable to start, the behaviour tests (norm mode and outcomes, the DA model fields, the RUVIII
  M, the results-summary file set) failed on their assertions. The e2e defect's test is first in [tests.txt](tests.txt).
- **Targeted** ([tests.txt](tests.txt), final tree, kit library, one process per file; 0 fail, 0 error, 0 skip each):
  `test-module-ci-trx-norm.R` 61, `-da.R` 50, `-summary.R` 70, `-app.R` 48, `test-workflow-omics-application-reads.R` 230,
  `test-e2e-browser-harness.R` 114. Also KIT-093/095/096's `test-trx-kernels.R` 373, `test-trx-design.R` 70,
  `test-trx-publication.R` 174, `test-workflow-trx-declarations.R` 308, and the trx import, design and QC module CI.
- **Dependents and regressions** ([dependents.txt](dependents.txt)): 56 files, about 9,000 expectations, 0 fail. These
  include:
  - the lipid, metab and proteomics module CI;
  - the lipid and metab declaration tests and the frozen dump (`test-workflow-omics-declarations.R`);
  - the KIT-094 characterisation;
  - every lipid characterisation and seam test that names a parameterised function.

  **No other omic's test file was modified.** The one first-pass error (`test-lipid-12-norm-module-contracts.R`: its seam
  mock takes `buildLipidNormQcTabsetPanel(ns)`) was fixed in the code (the profile is passed only when present). The skips
  are the same counts as at baseline. The lipid DA module CI, characterisation and contracts pass again after the e2e fix.
- **Lipid unchanged:** the rendered HTML of the lipid norm, DA and summary tabs and of `mod_lipidomics_ui` is identical
  to the baseline's.
- **R CMD check** ([check.log](check.log), HEAD `35f855bd`): 3 WARNINGs and 3 NOTEs. The dependency WARNING and the
  code-problems NOTE sections are identical to KIT-096's.
- **Install** ([installed-test.txt](installed-test.txt)): `kit_097` = `$P/core/browser-smoke.5gh5pn` (`35f855bd`); the RUVSeq
  and EDASeq family were carried. `test-module-ci-trx-summary.R` against the installed build: 70 pass, 0 fail, 0 skip. The
  path is in `$P/kit/scripts-097/install.txt`.
- **e2e** ([e2e-attempt2.txt](e2e-attempt2.txt)). The lifecycle check ([e2e-lifecycle.txt](e2e-lifecycle.txt)) showed the app
  process loads MultiScholaR from the kit library, so the e2e runs against the install. The KIT-040 `trx_edge_12` corpus
  (12 runs) was taken from the tile through import, design, QC, normalisation, DA, the RUVIII export (M = group), save,
  copy and report. The record names:
  - `trx_expression_filter`, `trx_qc_finalize` and `trx_ruvs`;
  - DA `source_state = trx_ruvs`, model `edgeR_QL`, `treat_applied` false;
  - the RUVIII outcome `group`.

  The integration TSV, the DA workbook and the RUVIII figure exist; there were no error notifications.

## Acceptance (with evidence)

- **Full GUI journey:** e2e (every state, the RUVIII export, the report). `M = technical_replicate` on a second run and the
  DA's invariance to `M` are shown in module CI `MCI-032.2`, where the DA manifest entity is unchanged by both exports; the e2e
  design has no technical replicates.
- **DA model:** `MCI-031.1` (record: `edgeR_QL`, `glmQLFit(robust = TRUE)`, `QL F-test (glmQLFTest)`, `BH`, the formula
  with `W_1`, `treat_applied = FALSE`). `MCI-031.4` and the e2e cover the tab's text.
- **Significance:** `MCI-031.1` on `trx_edge_12`:
  - the first contrast's topTags product has the KIT-040 `top_genes` sha256, and the tab's rows are `identical` to it;
  - up and down counts at FDR (BH) < 0.05 equal it;
  - the displayed table's `FDR` column (header "FDR (BH)") is `identical` to it.
- **RUV selection:** `MCI-030.1`: automatic records the mode, the grid settings, `k`, percentage, controls ("inferred") and
  scores; the winner's `k`, percentage and control set have the reference's sha256. `MCI-030.2`: manual k 2 at 10 % is
  recorded; k 4 above the cap 3 is refused with `trx_ruv_k_above_cap` and the head kept.
- **Report wording:** [report-grep.txt](report-grep.txt), from the e2e's rendered report: all 16 required statements found,
  all 11 forbidden absent (also asserted in `MCI-032.1`).
- **Output map:** `MCI-032.1` asserts the mapped files in `results_summary_dir` and the report's sections. **Not supplied:**
  the Glimma pages as files under `Interactive_Volcano_Plots/`. KIT-095 keeps them as store assets (the worker's volcano
  page). This goes back to KIT-095 per the ticket.
- **Integration export:** the export's matrix (`newY`) has the reference's sha256 (`MCI-032.1`). The TSV the summary writes
  is byte-identical to KIT-040's stable notebook output. **Declared:** the TSV read back is not bit-identical to `newY`
  (vroom's text precision, the notebook's own writer), so identity is asserted on the product and the file hash.
- **No whole matrix:** the KIT-082 probe's transcriptomics case covers the norm displays, the DA table, both volcanoes, the
  download and the RUVIII export. Managed, nothing reads the count tables and nothing is held. In process the same tabs
  read them (the control). A DA tab variant that holds the state is caught by both the probe and the walk. The application
  holds the grid's score table, the inferred control ids, the CCA plot data and references.
- **Technical-replicate option:** `MCI-032.3`: without technical replicates the option is not offered (`trxRuviiiMChoices`),
  and a request for it shows the `trx_ruviii_replicates_invalid` refusal; no export is made.
- **Hand-off:** below.

## Divergences and notes

- **Outside the file list:** `DESCRIPTION` Collate (the three new files) and the Rd usage of the three exported lipid
  modules (the `profile` argument), as KIT-095/096.
- **Not added:**
  - MCI-030..032 to `tests/testdata/mci/manifest.json`: `test-module-ci-manifest.R`, outside the file list, pins the
    transcriptomics set;
  - a lane to `tests/testdata/e2e/manifest.json`: `test-e2e-manifest.R` pins nine lanes. The e2e reads the hash-checked
    KIT-040 corpus instead.

  Both are a scope decision for the GUI workover.
- **Test corrections while writing the new tests** (no other omic's test was touched):
  - "limma" matched "Glimma" (word boundary added);
  - DT keeps the data's column names and shows the label in the header;
  - the Glimma status is `nonDE`/`upReg`/`downReg`;
  - the `M` choice is server-rendered.
- **KIT-096's `test-module-ci-trx-app.R`** was extended to the six tabs, statuses and ids.
- **The long table** gains `friendly_name` = `comparison` in the tab's copy, so that the lipid selectors filter by contrast.
- **Conventions:** the shinydashboard shell, `wellPanel` styling and managed runner are kept to match the other omics, with
  new functions at about 50 lines or fewer.

## Hand-offs

**KIT-098 (driver ids).** Module `transcriptomics_workflow`, tab input `transcriptomics_tabs`. The tab values are `norm`,
`de` and `summary`, with test ids `trx-tab-norm`, `trx-tab-de` and `trx-tab-summary`.

- **Norm (`norm-`):**
  - inputs: `ruv_mode` (`trx-norm-ruv-mode`), `auto_percentage_min`, `auto_percentage_max`, `auto_percentage_step`,
    `k_penalty_weight`, `ruv_percentage`, `ruv_k`;
  - actions: `run_normalization` (`trx-norm-run`), `reset_normalization` (`trx-norm-reset`);
  - outputs: `ruv_qc_ui`, `cancor_plot_counts`, `ruv_summary_counts`, `ruv_table_counts`, `{pca,rle}_{post_filter,post_norm,ruv_corrected}_assay1`
    (raw, TMM, RUVs), `norm_log`.
- **DA (`de-`):**
  - inputs: `da_q_val_thresh`, `treat_lfc_cutoff` (display only), `volcano_contrast`, `volcano_assay`,
    `volcano_interactive` (off by default), `table_contrast`, `table_assay`, `table_significance`, `table_max_rows`;
  - actions: `run_da_analysis` (`trx-da-run`), `download_da_results` (`trx-da-download`);
  - tabs: `da_results_tabs` (`trx-da-results-tabs`);
  - outputs: `formula_display`, `contrasts_display`, `da_status`, `volcano_static`, `volcano_glimma`, `da_summary_stats`,
    `da_results_table`.
- **Summary (`summary-`):**
  - the RUVIII export: `ruviii_m_source` (rendered in `ruviii_m_ui`), action `run_ruviii_export` (`trx-summary-ruviii-export`),
    outputs `ruviii_status` and `ruviii_before_after`;
  - the inputs `experiment_label` and `description`;
  - the actions `save_workflow_args` (`trx-summary-save-args`), `copy_to_publication` (`trx-summary-copy-publication`),
    `generate_report` (`trx-summary-generate-report`), `download_report` (`trx-summary-download-report`) and
    `export_session_state` (`trx-summary-export-state`).

The e2e helper learned `transcriptomics` (`.E2E_VALID_OMICS`, `e2e_input_id`, `e2e_switch_workflow_tab`).

**KIT-101 to 109 (GUI workover).**
- The reuse map as built is in [plan.md](plan.md) and the list above.
- The hooks are the profile entries in `trxOmicProfile()`/`trxNormDaSummaryProfile()`.
- The interactive-volcano risk is noted above.

**KIT-095.** The Glimma pages as files under `Interactive_Volcano_Plots/`.

Rollback: `git revert` of this ticket's commits; the evidence is kept.

---

## Stop record (coordinator, `fcf1093a`, kept)

The tabs were stopped before closure by the user's decision and later resumed with the narrowed scope above. The user judged
the hook approach (about 67 branches across the lipid DA, normalisation and summary files, a 455-line change to
`mod_lipid_norm_ui_helpers.R`, about 1,070 lines of trx module code) "lazy and incomplete": the lipid modules become an
accidental shared base, metabolomics keeps its own twin copy, and lipid assumptions leak into transcriptomics. The e2e
defect above is one such leak (a lipid default, the interactive volcano on, applied to a single-assay gene table).

**Superseded by** a GUI layer that does for the tabs what KIT-073 did for the declarations: shared omics GUI modules driven
by each omic's binding, with lipid, metab and transcriptomics as configurations (KIT-101 to 109). The hook code stays until
then as a working behavioural reference, not the final structure. The first e2e failure and the unfinished reproduction
(`stopped-wip-test-module-ci-trx-app.patch`, since applied and completed as `MCI-029.2`) are kept here.
