# KIT-095 audit: the transcriptomics binding (audit only, 2026-11-06)

Mode `audit_only`: no ticket, code or document edited; no tests, installs or journeys run. One quick `Rscript` probe of edgeR on a zero-row matrix (F05) and source prints of installed edgeR/ruv functions.

**Inputs audited.**

| Input | Identity |
|---|---|
| Repository HEAD | `bbe1fae9` (branch `performance-v2`). Dirty baseline: four untracked `kit-090/*.txt` files, not touched. |
| `tickets/kit/KIT-095.md` | sha256 `e1d403f3…0e47`, revision 1.0, `pending`, `needs_planning: true`. Validator: valid, 0 errors, 0 warnings, status consistent. |
| Notebook `dev/distributable/plans/transcriptomics_workflow.rmd` | sha256 `44d5f643…652a`, 793 lines. |
| Runbook `transcriptomics_architecture_runbook.md` | sha256 `ad71283c…8e88f01`. |
| Protocol `docs/qa/workflow-execution/kit/protocol.md` | sha256 `1531ad5b…d6bd`, v1.6. |
| Installed packages | R 4.6.1, reference BLAS/LAPACK 3.12.0, edgeR 4.10.1, limma 3.68.4, ruv 0.9.7.1, Glimma 2.15.2. **RUVSeq and EDASeq are not installed** in either library path. |

**Context read.**
- Neighbours: KIT-040 and KIT-093 to KIT-099.
- The approved plan `~/.claude/plans/write-a-series-of-atomic-dragon.md` and `00-KIT-PLAN.md` v3.8.
- The core files named in the brief.
- The lipid binding and its test, and the registry and omics tests.

## Finding register

### KIT095-F01: the DA machinery KIT-095 would reuse only accepts limma

| | |
|---|---|
| Severity | **Blocking** |
| Dimensions | 2 (scientific rigour), 4 (reference accuracy) |
| Blocked state | Implementation (needs a scope decision) |

**Evidence.**
- `newSavedWorkflowDaFit` refuses anything that is not a limma `MArrayLM` (`R/utils_workflow_saved_entities.R:513-528`). The runbook says the same (line 24).
- `omicsDaAssayFit` saves every fit through that envelope (`R/utils_workflow_omics_declarations.R:410-412`). The save is outside the `tryCatch`, so an edgeR `DGEGLM` fit makes the step abort with `saved_entity_invalid`.
- `omicsDaContrastOperation` has limma assumptions built in (`:427`, `:452`):
  - it fixes `eBayes_trend`/`eBayes_robust`;
  - it calls `da_apply(fit, contrast, lfc|NA, TRUE, TRUE)`, which is the treat/eBayes signature.
- `acceptWorkflowOmicsDaManifest` writes `model = list(model = "limma", eBayes_trend = TRUE, eBayes_robust = TRUE, treat_applied = …)` (`:1091-1092`). The analysis record copies `manifest$model` (`R/utils_workflow_analysis_record.R:306`).
- `readWorkflowOmicsDaResults` rebuilds the display list from `fit.eb$coefficients` (`:513`).

**What the notebook does instead.**
- The DA design is derived from the data: `~ 0 + group + W_1..W_k` (notebook :528-533).
- The fit is robust `estimateDisp`/`glmQLFit` followed by `glmQLFTest`, with no treat (:536-554).
- Its table has the columns `logFC, logCPM, F, PValue, FDR`.

**Impact.**
- "DA publication and tab results, reused through adaptors" and `trx.da_results.v2` cannot work unchanged.
- If the manifest were reused, the record would misreport the model as limma/eBayes/treat.
- Replay (`R/utils_workflow_analysis_record.R:604-609`) passes only `formula, qvalue, lfc`. It cannot rebuild a formula that comes from W.
- The fix needs `R/utils_workflow_saved_entities.R`, which is not in `files_to_modify`.

**Smallest repair.** Revise the ticket to declare an explicit edgeR product contract (the runbook asks for one at line 48):
- a DA-fit envelope for a `DGEGLM`, or a declared decision not to save the fit;
- a TRX `da_contrast` declaration whose `fixed` and `model` say `edgeR QL, robust = TRUE`, with no treat;
- the realised formula, including the W columns, recorded as a parameter or outcome;
- a statement on multiple contrasts: the notebook tests only the first line of `contrast_strings.tab`, so more than one is a behaviour extension;
- `utils_workflow_saved_entities.R` added to `files_to_modify` with the user's approval.

**Closure evidence.**
- A ticket revision that names the contract.
- A test showing a TRX contrast saves and re-reads its fit (or omits it, as declared).
- The manifest's `model` reads edgeR.
- The record replays the DA exactly.

---

### KIT095-F02: no ticket owns the publication renderers

| | |
|---|---|
| Severity | **Blocking** |
| Dimensions | 1 (completeness), 3 (best practices) |
| Blocked state | Implementation (needs a scope decision) |

**Evidence.** The shared publication handlers call binding adaptors:
- `qc_plots`, `composite_figure`, `da_export`, `final_pca`, `volcano_widget`, `filtering_metrics` (`R/utils_workflow_omics_declarations.R:1155-1280`);
- the lipid binding fills them with existing lipid functions (`R/utils_workflow_lipid_declarations.R:84-92`).

No ticket supplies the TRX functions behind these adaptors:
- KIT-093's kernels stop at "PCA and RLE inputs; CCA diagnostics" (KIT-093 outcome).
- KIT-095's `files_to_modify` has no `R/func_trx_*.R`.
- KIT-097's files are `mod_trx_*` and `inst/reports/` only.

The notebook's renderers therefore have no owner: `format_cancor_plot`, `format_tuning_curve`, the workbook writer, `generatePcaPlot`/`generateRlePlot` and the 6-panel grid, `generateVolcanoPlot`, `generateBarPlot`, `glimmaVolcano`, the DE CSV and top-genes TSV, and the RUVIII before/after grid.

**Impact.** The outcome "publication operations, reused through adaptors" cannot be met without writing package code outside the allowlist. This is the AGENTS.md rule against materially expanding scope without a decision.

**Smallest repair.** Name one owner and add it to the allowlist:
- either extend KIT-093's `R/func_trx_qc_plots.R` to "renderers", or
- add `R/func_trx_publication.R` to KIT-095.

List the adaptor-to-function map in the ticket.

**Closure evidence.** A ticket revision with the owner and map, plus a test per adaptor on the synthetic corpus.

---

### KIT095-F03: the filter's pre-pass is misstated, and the logCPM pre-pass is incomplete

| | |
|---|---|
| Severity | Should-fix |
| Dimension | 2 (scientific rigour) |
| Blocked state | Implementation |

**Evidence.**
- The notebook runs `calcNormFactors(dge_raw, "TMM")` on the **unfiltered** counts, then `filterByExpr(dge_raw, group = …)` (notebook :166-170).
- Installed edgeR 4.10.1 `filterByExpr.DGEList` uses `lib.size <- y$samples$lib.size * y$samples$norm.factors`.
- `filterByExpr.default` derives three things from the whole data:
  - `MinSampleSize` from the group (`min(tabulate(group)[>0])`, with the `large.n`/`min.prop` adjustment);
  - `MedianLibSize`;
  - `CPM.Cutoff = min.count / median(lib.size) * 1e6`.
- After those, each row's decision depends only on that row.
- The zero-variance filter (`apply(…, 1, var)`, :174-175) is also per row, and the ticket does not mention it.
- For logCPM, `cpm.DGEList` uses `lib.size * norm.factors` of the **filtered** population (:178-180). Given that per-sample vector, `cpm.default` → `.cxx_calculate_cpm_log` is per element, and the scaled prior uses the mean of the lib-size vector. So logCPM is row-independent once the vector is fixed.
- The runbook (line 39) says: "Library sizes and normalization factors belong to the population used at each call."

**Impact.**
- The ticket's "row after the library-size pre-pass" leaves out the initial TMM factors. A pre-pass of `colSums` alone gives a different keep mask whenever any TMM factor ≠ 1.
- The test "partitioned identical to whole" would catch it, but only by spending the single correction pass.
- The TMM pre-pass is itself a whole-matrix `calcNormFactors` read (edgeR calls are not rewritten). The partitioned filter's memory benefit is therefore bounded by that read, which matters for KIT-099's pairs.
- "With the factors as a pre-pass" for logCPM leaves out the filtered library sizes.

**Smallest repair.** Restate the units:
- **`trx_expression_filter`:** rows, given these global statistics:
  - the effective library sizes (`colSums` of the unfiltered counts × the TMM factors on the unfiltered counts, computed whole);
  - their median;
  - `MinSampleSize` from the group;
  - the `filterByExpr` defaults.

  The zero-variance filter follows, per row.
- **`trx_tmm_logcpm`:** TMM whole on the filtered counts. logCPM is per row given `colSums(filtered) × factors` and `prior.count = 2`.

Record these as `global_statistics` in the declarations and the step map.

**Closure evidence.** The step-map rows carry edgeR 4.10.1 source-line evidence, and a test shows the partitioned keep mask and logCPM are `identical` to whole on a corpus whose TMM factors are not 1.

---

### KIT095-F04: the data flow between steps is undefined, and the core cannot carry it

| | |
|---|---|
| Severity | Should-fix |
| Dimensions | 1 (completeness), 3 (best practices) |
| Blocked state | Implementation (the planning step must resolve it) |

**Evidence.** KIT-093's class holds one count assay (`count_data`, "one assay"). The steps need more than that:
- `trx_ruv_grid` needs logCPM and the filtered counts with TMM;
- `trx_ruvs` needs the rounded filtered counts, the controls and k;
- `trx_da_contrast` needs counts and W;
- `trx_ruviii_export` needs logCPM, the controls, k and M;
- the three-stage QC figures need raw logCPM (filtered counts with no norm factors, notebook :490), TMM logCPM and the RUVs-corrected logCPM (:421-426).

The ticket does not say which of these are states and which are products. The core has gaps too:
- `computeWorkflowOmicsOperationAction` builds dependencies only for a `pairs` parameter (`R/utils_workflow_omics_declarations.R:1138`).
- `replayWorkflowOmicsStep` switches only on lipid/metab step names, and the `preceded_by` product replay exists only for `pearson_selection` (`:1330-1351`). A TRX step would hit `workflow_replay_unsupported`.
- The `filtering_metrics` display is keyed on `_intensity_filter`/`_duplicate_resolution` (`:1148`).

KIT-094's list of hard-coded items (KIT-094 outcome) includes none of these.

**Impact.**
- The acceptance criterion "every state identical to KIT-093's kernels" cannot be evaluated until the states are defined.
- The grid → RUVs → DA → RUVIII products would not be dependencies, so their identities would not include them, and replay would fail.

**Smallest repair.** In the planning step, write a state/product table for each step: inputs, outputs, the dependency parameter names, and which output a replay reproduces. Name the core generalisations KIT-095 makes (dependency parameters, the replay switch, the filtering-metrics names) and keep lipid/metab byte-identical through the frozen dump.

**Closure evidence.** The plan's table, plus a test that replays every TRX step, including the product-only ones, exactly.

---

### KIT095-F05: the partitioned machinery cannot run the edgeR kernels

| | |
|---|---|
| Severity | Should-fix |
| Dimension | 4 (reference accuracy) |
| Blocked state | Implementation |

**Evidence.**
- `omicsRunPartitioned` first runs `on_template(template)`, the unchanged kernel on zero-row assays (`R/utils_workflow_omics_declarations.R:63-64`), to get the recorded value and the design.
- Probe: `edgeR::DGEList(matrix(integer(0), 0, 4), group = …)` fails with "'counts' must contain at least one value".

**Impact.** The partitioned filter, and logCPM if partitioned, cannot reuse `omicsRunPartitioned` with an `on_template` that calls the kernel.

**Smallest repair.** Declare a TRX `on_template` that returns the shell without calling edgeR, with the counts recorded from the stored shapes. Alternatively, give the kernels a path that bypasses edgeR on an empty template, tested.

**Closure evidence.** The partitioned test passes, with a test showing the template path never calls edgeR.

---

### KIT095-F06: the integration export does not produce the notebook's MOFA file

| | |
|---|---|
| Severity | Should-fix |
| Dimensions | 1 (completeness), 4 (reference accuracy) |
| Blocked state | Implementation |

**Evidence.**
- `computeWorkflowOmicsIntegrationExportAction` writes `saveRDS(manager$getState(state_name))`, the S4 state (`R/utils_workflow_omics_declarations.R:1297-1305`), and reads no binding.
- The notebook's export is the RUVIII `newY` as a sample × gene TSV with a `Run` column, written to two paths (`results/transcriptomics/…` and `results/Transcriptomics/…`, notebook :715-724).
- `newY` is not a count state. The runbook (line 118) asks for distinct identities for the three expression products.

**Impact.** "Integration export (RUVIII output to MOFA)" reused as-is would export the wrong object.

**Smallest repair.** Add a binding adaptor or TRX handler that writes the saved `newY` product in the notebook's format. State whether both paths are kept (runbook line 47).

**Closure evidence.** The exported TSV is byte-identical to the KIT-040 reference file hash, or the difference is declared before the comparison.

---

### KIT095-F07: the CCA-diagnostics operation does not exist, and the publication list is shared by every binding

| | |
|---|---|
| Severity | Should-fix |
| Dimensions | 1 (completeness), 4 (reference accuracy) |
| Blocked state | Implementation |

**Evidence.**
- The registry registers one fixed list of nine publication operations for every binding (`R/utils_workflow_operations.R:703-714`). It has no CCA operation.
- The ticket lists "CCA diagnostics" as reused.
- KIT-094 makes the steps per binding, but not the publication list.
- `test-workflow-operation-registry.R:94-109` pins the exact set of operation names.

**Impact.** Adding a CCA operation to the shared list registers `lipid_`/`metab_` CCA operations, which changes the lipid/metab registry and breaks the frozen-dump and registry tests. Not adding it leaves the CCA figure and workbook without a worker operation.

**Smallest repair.** Make the publication list a binding field, defaulting to the current nine, and give the TRX binding its own list that includes the CCA diagnostics.

**Closure evidence.** The registry test is updated for the `trx_*` names only, and the frozen dump is unchanged.

---

### KIT095-F08: several notebook outputs map to no declaration

| | |
|---|---|
| Severity | Should-fix |
| Dimension | 1 (completeness) |
| Blocked state | Implementation |

**Evidence: unmapped or implicit outputs.**

| Notebook output | Notebook lines | Ticket mapping |
|---|---|---|
| `canonical_correlation_optimization_summary.xlsx` (two sheets) | :384-393 | none |
| `ruv_normalised_log_values.tsv` (RUVs-corrected logCPM) | :429-432 | none |
| `ruviii_before_and_after.png/.pdf` | :730-745 | none |
| `ruv_optimization_summary.txt` (provenance plus the full grid) | :755-775 | none |
| `cache/ruv_optimization_results.rds`, `cache/de_results.rds` | :778-779 | none. The runbook (line 53) asks for an explicit consumer migration. |
| `top_de_genes.tsv` | :583 | implicit in "DA publication" |
| `summary_barplot_<contrast>` | :662 | implicit in "DA publication" |
| raw logCPM, the first stage of the 6-panel grid | :490 | no product |

Mapped:
- The CCA panels map to "CCA diagnostics" (but see F07).
- The 6-panel grid maps to "QC figures/composite".
- The volcano, the DE CSV and the Glimma page map to "DA publication/volcano page" (but see F01 and F02).

**Impact.** The KIT-040 reference freezes "output file hashes". Without an owner, outputs are either dropped silently or added ad hoc.

**Smallest repair.** Add a notebook-output → operation table to the ticket. Each output is mapped, or explicitly dropped with a reason (for example, the cache RDS is replaced by saved products).

**Closure evidence.** The table appears in the ticket or plan, and each mapped output is checked against its reference hash, or declared non-byte-comparable in advance (PNG/PDF/XLSX carry timestamps).

---

### KIT095-F09: `files_to_modify` misses files the change must touch

| | |
|---|---|
| Severity | Should-fix |
| Dimension | 3 (best practices) |
| Blocked state | Implementation |

**Evidence.**
- `tests/testthat/test-workflow-operation-registry.R:94-109` asserts the exact set of operation names.
- `tests/testthat/test-workflow-omics-declarations.R:60` asserts `names(omicsDeclarationBindings()) == c("lipidomics", "metabolomics")`.
- `tests/testthat/test-workflow-omics-worker.R:44` checks coverage only per lipid/metab.
- `inst/config/config.ini` has only legacy `[transcriptIntensityFiltering]` and `[removeTranscriptsWithOnlyOneReplicate]` sections (:114-120), not the notebook's parameters. Any `workflowOperationDefaults` mapping to a shipped key needs that file.
- The parameter-label and defaults tests (`test-workflow-operation-registry.R:263-285`) demand a label and a default mapping for every declared parameter.

None of these files is listed.

**Impact.** Adding the binding breaks two existing tests outside the allowlist. Shipped defaults are either impossible or silently "no_shipped_default".

**Smallest repair.**
- Add the two test files to `files_to_modify`; the changes are limited to the new names.
- State that the TRX defaults are `fixed`/`no_shipped_default`, or add `config.ini` with the user's approval.
- Do not reuse the legacy transcript sections.

**Closure evidence.** The targeted tests pass, and the lipid/metab assertions are unchanged.

---

### KIT095-F10: the step-map gate cannot be met within this ticket, and it is sequenced inside the declaring ticket

| | |
|---|---|
| Severity | Should-fix |
| Dimensions | 3 (best practices), 2 (scientific rigour) |
| Blocked state | Decision |

**Evidence.**
- Protocol §2, the unit evidence gate (`protocol.md:24`): "a step is declared only when the tracked `kit-001/step-map.json` (hashed into the bindings) records its independence unit … and no open `flag`."
- §5 (`:57`): each route's step map is "hashed in the bindings".
- KIT-095 lists `step-map-trx.json` but neither `docs/qa/workflow-execution/kit/bindings-generator.R` nor `bindings.json`. KIT-098 owns `bindings-generator.R` and protocol v1.7.
- For lipid, the step map was a separate ticket (KIT-068), and its flags were answered by the user before KIT-031 declared anything.

**Impact.**
- The acceptance criterion "the step map passes the protocol unit gate" is not achievable within the ticket's files.
- Writing the map and the declarations together removes the stop point where a flag reaches the user. The filter's TMM-dependent pre-pass (F03) and RUVIII `M` are candidate flags.

**Smallest repair.**
- Order the work: step map first; any `flag` stops the ticket for the user's answer; then the declarations.
- Either add the bindings regeneration to KIT-095, or reword the criterion to "map complete, no open flag; hashed into bindings in KIT-098".

**Closure evidence.** The map is committed before the declarations, with `open_flags: 0` (or the user's answers recorded) and the bindings hash present where the protocol requires it.

---

### KIT095-F11: downstream parameters and record outcomes are not fixed here, though later tickets cannot change them

| | |
|---|---|
| Severity | Should-fix |
| Dimensions | 1 (completeness), 3 (best practices) |
| Blocked state | Implementation |

**Evidence.**
- KIT-097 promises "k and percentage selection" and the "RUVIII `M` choice", and its acceptance criterion says "the record names the RUVIII `M` choice". KIT-097's files exclude the declarations and the record.
- The notebook selects k and the percentage automatically (:295-298). Lipid declares `manual_k`/`manual_percentage` (`omicsRuvParameterNames`, `R/utils_workflow_omics_declarations.R:276-277`).
- The series contract says "Inferred control genes are labelled inferred". The ticket's "record outcomes" names neither that label, nor M, nor the scores.
- The lipid outcomes read a per-assay shape (`R/utils_workflow_analysis_record.R:143-157`).
- The report sources (`:738-790`) would be written before KIT-097's template exists.

**Impact.**
- The GUI tickets would need to reopen KIT-095's files, or ship without the controls.
- The inferred-controls label and the M choice could be missing from the record.

**Smallest repair.** Declare these now:
- the grid/RUVs parameters, including any manual k/percentage override, with the notebook's automatic choice as the default;
- `ruviii_m_source ∈ {group, technical_replicate}`;
- record outcomes: k, percentage, control count, `controls = "inferred"`, separation and composite scores, `M`, filter counts before/after, realised DA formula;
- report-source variable names that reuse the existing `best_k`, `best_percentage`, `num_neg_ctrl`, `ruv_applied` and `s4_de_*` where they fit.

**Closure evidence.** Registry and record tests assert these parameters and outcomes.

---

### KIT095-F12: the verification does not cover the worker path, the products or the publication outputs

| | |
|---|---|
| Severity | Should-fix |
| Dimension | 3 (best practices) |
| Blocked state | Implementation |

**Evidence.**
- The ticket declares every compute step "on the worker", but its acceptance criteria are an in-process journey, partitioned == whole, and the gate. The worker is first exercised in KIT-098's journeys.
- An in-process worker harness exists (`tests/testthat/helper-workflow-omics-worker.R:10-27`, `workflowWorkerMain` under `MULTISCHOLAR_WORKFLOW_WORKER=1`), and the lipid/metab worker test covers every `omics`-role declaration (`test-workflow-omics-worker.R:88-130`).
- Acceptance criterion 1 compares "every state" only. The grid scores, best k/percentage, W, the edgeR table, the corrected logCPM and `newY` are products, not states.
- The approved plan's tests (registry shape; refusal off the owner) are not in the acceptance criteria.

**Impact.** A dispatch, handler or completion defect surfaces only in a counted KIT-098 journey, and product parity is unverified until then.

**Smallest repair.** Add acceptance criteria for:
- each `trx_*` worker step and publication operation run through the in-process worker harness;
- `identical` products against KIT-093 and the KIT-040 reference;
- the registry-shape and refusal-off-owner tests.

All of these stay within the "in-process only" budget.

**Closure evidence.** The new test results.

---

### KIT095-F13: determinism is not stated in the step map

| | |
|---|---|
| Severity | Minor |
| Dimension | 2 (scientific rigour) |
| Blocked state | None |

**Evidence.**
- The path has no RNG:
  - `ruv::ruv_cancorplot` uses `svd`/`cancor` only (printed source);
  - `RUVs` and `RUVIII` use SVD;
  - `estimateDisp`/`glmQLFit`/`glmQLFTest` are deterministic.
- The grid is sequential (`purrr::map`).
- Workers pin BLAS/OMP threads to 1 (`R/utils_workflow_execution.R:1163`). In-process runs are not pinned, but this host uses the reference BLAS.

**Impact.** Low. Worker and in-process runs could differ in the last bits only on a host with multithreaded BLAS.

**Smallest repair.** One `kernel_observations` line in `step-map-trx.json`: no RNG, sequential grid, and the BLAS-thread assumption. No `set.seed` is needed.

**Closure evidence.** The line is present.

---

### KIT095-F14: "as in KIT-085" points at a helper that does not apply

| | |
|---|---|
| Severity | Minor |
| Dimension | 4 (reference accuracy) |
| Blocked state | None |

**Evidence.**
- `prepareOmicsRuvCancorObject` (`R/func_metab_norm.R:288-312`) runs NIPALS imputation for assays with missing values. Counts have none (whole-number policy), so it does not apply.
- The notebook's per-percentage invariants are `t(log_cpm_tmm)` and the preliminary QL ranking (:232-237, :253).
- Precomputing the "All" SVD and passing it as `W1` to `ruv_cancorplot` would change the call. That is an unmeasured optimisation (AGENTS.md).

**Impact.** Ambiguity: someone could apply imputation to counts, or optimise speculatively.

**Smallest repair.** Reword to: "the transpose and the preliminary ranking computed once (the KIT-085 pattern); no imputation; `ruv_cancorplot` called as in the notebook."

**Closure evidence.** The ticket text.

---

### KIT095-F15: terminology and registry defaults are not stated

| | |
|---|---|
| Severity | Minor |
| Dimension | 4 (reference accuracy) |
| Blocked state | None |

**Evidence.**
- The ticket's unit column uses "whole". In the protocol and step maps, `whole_matrix` is the unit and `whole` is the variant (`kit-001/step-map-lipid.json` units).
- The registry's default variant is keyed on the name `<prefix>_intensity_filter` (`R/utils_workflow_operations.R:695`). `trx_expression_filter` will therefore default to `whole`, which is consistent with KIT-099 deciding later, but this is unstated.
- The `TRANSCRIPTOMICS` section map is meaningful only if KIT-093's kernels write `@args` sections (`R/utils_workflow_operations.R:55-67`). Otherwise the worker steps declare `section = NULL`.
- The index title ("partitioned where exact") differs from the file title.

**Impact.** Small.

**Smallest repair.** Use `whole_matrix` for the unit; state "default variant `whole` until KIT-099"; state the section-map content or that it is empty.

**Closure evidence.** The ticket text.

---

### KIT095-F16: environment preconditions and directory roles are not stated

| | |
|---|---|
| Severity | Minor |
| Dimension | 4 (reference accuracy) |
| Blocked state | Execution |

**Evidence.**
- RUVSeq and EDASeq are not installed (`.libPaths()` = user 4.6 library plus `/usr/lib/R/library`).
- R is pinned at 4.6.1 under the series' `IgnorePkg` rule, so installing Bioconductor packages needs care.
- **Verified, positive:** a QC directory role for transcriptomics exists. `gene_qc_dir` is mapped from `feature_qc_dir` (`R/func_general_filemgmt_setupdirectories_helpers.R:216-217`); the transcriptomics config is at `:78-93`. The DA leaf is `da_genes`; the notebook uses `de_genes`.

**Impact.** The in-process journey cannot run on this host until KIT-040/093 install the packages into the library the tests and worker use.

**Smallest repair.** Name `qc_dir_role = "gene_qc_dir"` and the DA directory mapping in the ticket. State the package precondition, inherited from KIT-093, with `skip_if_not_installed` or a hard requirement decided.

**Closure evidence.** The ticket text, plus `packageVersion("RUVSeq")` recorded in the KIT-095 evidence.

## Verdict per dimension

1. **Completeness: not ready.**
   - The six steps and the store steps cover the notebook's computation.
   - The publication side is incomplete. No renderer has an owner (F02). The CCA operation does not exist (F07). The MOFA export would write the wrong object (F06). Six notebook outputs map to nothing (F08).
   - The data flow between steps and its replay are undefined (F04). Downstream parameters and outcomes are not fixed here (F11).
2. **Scientific rigour: mostly sound, with one misstated unit and one engine mismatch.**
   - These claims are right:
     - TMM is whole.
     - logCPM with `prior.count` is row-independent given the per-sample effective library sizes.
     - The edgeR DA is whole.
     - There is no RNG in the grid; the workers pin BLAS threads.
   - The filter's pre-pass leaves out the initial TMM factors (F03). `filterByExpr` is row-independent only given effective library sizes and group/design-derived scalars.
   - The edgeR fit cannot pass through the limma-only DA products, and reusing the manifest would misrecord the model (F01).
3. **Implementation best practices: partly.**
   - Done well: test-first is stated, existing adaptors and the binding pattern are reused, and the budget is bounded.
   - Gaps:
     - `files_to_modify` misses tests and config the change must touch (F09);
     - the step-map gate is sequenced inside the declaring ticket and needs bindings this ticket does not own (F10);
     - the worker path and the products go untested until counted journeys (F12).
4. **Codebase reference accuracy: several reuse claims fail on inspection.**
   - These do not accept TRX as-is: the DA fit envelope (`MArrayLM` only), the omics contrast declaration and manifest (limma), the integration export (S4 RDS), `omicsRunPartitioned` (edgeR refuses its zero-row template), the publication list (no CCA), and the KIT-085 helper (not applicable).
   - These are accurate: `newOmicsOperation`, `omicsRunPartitioned`, the `computeWorkflowOmics*Action` handlers, `workflowOutputDeclarations` and its v2 manifest pattern, `gene_qc_dir`, and the record/report-source tables all exist as described.

**Overall: KIT-095 is not ready to plan against as written.** Two blocking findings need scope decisions before implementation:
- F01, an edgeR DA product contract and `utils_workflow_saved_entities.R`;
- F02, an owner for the renderers.

The ten should-fix findings belong in the ticket revision that the planning gate already requires.
