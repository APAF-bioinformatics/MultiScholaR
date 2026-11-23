# The accepted per-route action scripts (KIT-159 7.4; VAL-41, VAL-26), the data of
# valActionScript. Placed in its own file beside driver.R (7.2 names valActionScript under
# driver.R; the split keeps each file readable and is recorded in outcome.md).

#' One action-script step
#'
#' @param id Step id.
#' @param tab The tabset value the step opens (`NULL`: none).
#' @param declared The state the step must reach: `list(head = <states>)` (any of them),
#'   `list(tab = c(<status key> = <status>))`, or `"diagnostic"` (no state change expected).
#' @param run `function(api)`: the step's actions through the driver's api.
#' @param structural Regular expressions of the input ids the step may set as structural.
#' @param press The buttons the step presses (documentation and the enrichment check).
#' @param module The module id whose rendered outputs are read after each action (`NULL`: all).
#' @return A step (list).
#' @export
valStep <- function(id, tab, declared, run, structural = character(), press = character(),
                    module = NULL) {
    list(id = id, tab = tab, declared = declared, run = run, structural = structural,
         press = press, module = module)
}

#' The ten routes of section 7.4: omic, workflow type and import format (the e2e manifest's)
#' @noRd
valRouteTable <- function() {
    data.frame(
        route = c("prot_dia", "prot_dia_limpa", "prot_lfq", "prot_lfq_fragpipe", "prot_tmt",
                  "trx_screen", "lipid_canonical", "metab_lc", "metab_gc", "metab_combined"),
        omic = c(rep("proteomics", 5L), "transcriptomics", "lipidomics", rep("metabolomics", 3L)),
        workflow = c("DIA", "DIA_LIMPA", "LFQ", "LFQ", "TMT", NA, NA, NA, NA, NA),
        import_tool = c("diann", "diann", "maxquant", "fragpipe", "pd_tmt", "trx_counts",
                        "lipidsearch", "custom", "custom", "custom"),
        stringsAsFactors = FALSE
    )
}

#' A route's row of the route table
#' @noRd
valRoute <- function(route) {
    rt <- valRouteTable()
    row <- rt[rt$route == route, , drop = FALSE]
    if (nrow(row) != 1L) valAbort("no accepted action script for route {route}", "val_spec_invalid")
    as.list(row)
}

#' The state each proteomics QC sub-tab commits (the GUI-entry oracles' heads, tests/testdata/
#' gui-entry; re-read by G1(b))
#' @noRd
PROT_QC_STATES <- c(
    peptide_qvalue_filter = "qvalue_filtered",
    peptide_protein_peptide_filter = "protein_peptide_filtered",
    peptide_rollup = "precursor_rollup",
    peptide_intensity_filter = "intensity_filtered",
    peptide_sample_filter = "sample_filtered",
    peptide_replicate_filter = "replicate_filtered",
    peptide_imputation = "imputed",
    peptide_log2_matrix = "log2_matrix",
    peptide_accession_cleanup = "peptide_accession_cleaned",
    protein_rollup = "protein_s4_created",
    protein_cleanup = "protein_accession_cleaned",
    protein_intensity_filter = "protein_intensity_filtered",
    protein_duplicate_removal = "duplicates_removed",
    protein_replicate_filter = "protein_replicate_filtered",
    protein_dpc_quant = "protein_dpc_quant",
    protein_normalisation = "protein_normalised",
    protein_peptide_count = "peptide_count_filtered"
)

#' The accepted action script of a route (section 7.4)
#'
#' Import (the declared files and the route's format; every other import input at its rendered
#' default), the design builder, the route's QC sub-tabs in the product's spec order, each at its
#' rendered defaults, normalisation (the arm's overrides only), the correlation filter where the
#' tab renders it, DIA_LIMPA's covariates button, DA at rendered defaults; Summary "Save Workflow
#' Arguments" only with `summary`. Never pressed: revert, skip, export, load, report, enrichment,
#' Cancel. The QC sub-tabs are read from `getProtQcPeptideModuleSpecs` and
#' `getProtQcProteinModuleSpecs` of the MultiScholaR namespace.
#'
#' @param route Route (lane) id.
#' @param assays Assay names (lipidomics and metabolomics).
#' @param import_tool The import format (default the route's).
#' @param import_values Further structural import inputs (custom column mapping), by input id.
#' @param contrasts List of group pairs, the first the builder's first contrast.
#' @param summary `list(label, description)` to continue to the Summary step, or `NULL`.
#' @param alternative `NULL` (the accepted script) or one of G1(b)'s alternatives, which no
#'   consumer runs (r1.14): `export_load` ("Export session" after the last normalisation step,
#'   then "Load Filtered Session" on the DA tab before Run; not on transcriptomics, which has no
#'   load button) and `no_covariates` (DIA_LIMPA without "Compute RUV-III Covariates").
#' @return A list: `route`, `omic`, `workflow`, `contrasts`, `alternative`, `steps`.
#' @export
valActionScript <- function(route, assays = NULL, import_tool = NULL, import_values = list(),
                            contrasts = NULL, summary = NULL, alternative = NULL) {
    r <- valRoute(route)
    valCheckAlternative(r, alternative)
    tool <- import_tool %||% r$import_tool
    steps <- switch(r$omic,
                    proteomics = valProtScript(r, tool, summary, alternative),
                    transcriptomics = valTrxScript(summary),
                    valOmicsScript(r, if (r$omic == "lipidomics") "lipid" else "metab", tool,
                                   assays, import_values, summary, alternative))
    list(route = route, omic = r$omic, workflow = r$workflow, contrasts = contrasts,
         alternative = alternative, steps = steps)
}

#' Refuse an alternative the route does not have (G1(b); r1.14)
#' @noRd
valCheckAlternative <- function(r, alternative) {
    if (is.null(alternative)) return(invisible(TRUE))
    ok <- switch(alternative,
                 export_load = !identical(r$omic, "transcriptomics"),
                 no_covariates = identical(r$workflow, "DIA_LIMPA"),
                 FALSE)
    if (!isTRUE(ok)) {
        valAbort("route {r$route} has no alternative {alternative}", "val_spec_invalid")
    }
    invisible(TRUE)
}

#' The export-then-load steps (A17; G1(b)): the tab's export button after the last normalisation
#' step, then the DA tab's "Load Filtered Session" with no file chosen (the latest export); neither
#' changes the head, which stays at the last normalisation state
#' @noRd
valExportLoadSteps <- function(head, norm_tab, norm_module, export_button, da_tab, da_module,
                               da_ui) {
    force(export_button)
    force(da_ui)
    load_button <- paste0(da_module, "-load_filtered_session")
    list(
        valStep("export_session", norm_tab, list(head = head), module = norm_module,
                press = export_button, run = function(api) api$press(export_button)),
        valStep("load_session", da_tab, list(head = head), module = da_module, press = load_button,
                run = function(api) {
                    api$defaults(da_ui, da_module)
                    api$press(load_button)
                })
    )
}

#' The proteomics script
#' @noRd
valProtScript <- function(r, tool, summary, alternative = NULL) {
    ns <- asNamespace("MultiScholaR")
    workflow <- r$workflow
    peptide_level <- workflow %in% c("DIA", "DIA_LIMPA")
    limpa <- identical(workflow, "DIA_LIMPA")
    specs <- c(
        if (peptide_level) lapply(ns$getProtQcPeptideModuleSpecs(workflow), function(s) {
            c(s, list(group = "peptide_qc", tabs = "peptide_filter_tabs"))
        }),
        lapply(ns$getProtQcProteinModuleSpecs(workflow), function(s) {
            c(s, list(group = "protein_qc", tabs = "protein_filter_tabs"))
        })
    )
    qc_tabset <- if (peptide_level) "quality_control-qc_tabs_lfq" else "quality_control-qc_tabs_tmt"
    qc_structural <- paste0("^quality_control-(qc_tabs_lfq|qc_tabs_tmt|",
                            "peptide_qc-peptide_filter_tabs|protein_qc-protein_filter_tabs)$")
    qcStep <- function(spec) {
        force(spec)
        key <- paste0(sub("_qc$", "", spec$group), "_", spec$id)
        state <- PROT_QC_STATES[key]
        if (is.na(state)) {
            valAbort("QC sub-tab {key} has no declared state (flag to G12)", "val_spec_invalid")
        }
        sub_id <- paste0("quality_control-", spec$group, "-", spec$id)
        valStep(paste0("qc_", key), "qc", list(head = unname(state)), module = "quality_control",
                structural = qc_structural, press = paste0(sub_id, "-apply_*"),
                run = function(api) {
                    ui_fn <- get(spec$uiFn, envir = ns)
                    group_label <- if (identical(spec$group, "peptide_qc")) "Peptide QC" else
                        "Protein QC"
                    api$set(stats::setNames(list(group_label), qc_tabset))
                    api$set(stats::setNames(list(spec$title),
                                            paste0("quality_control-", spec$group, "-", spec$tabs)))
                    api$defaults(ui_fn, sub_id)
                    actions <- grep("-apply_", api$helpers$guiEntryActionIds(ui_fn(sub_id)),
                                    value = TRUE)
                    if (!length(actions)) {
                        valAbort("no apply action in {spec$uiFn}", "val_spec_invalid")
                    }
                    api$press(actions[[1L]])
                })
    }
    import <- valStep(
        "import", "setup", list(tab = c(setup_import = "complete")), module = "setup_import",  # nolint
        structural = paste0("^setup_import-(search_results_standard|fasta_file_standard|",
                            "format_override|dia_route)$"),
        press = "setup_import-process_data",  # nolint
        run = function(api) {
            api$defaults(ns$mod_prot_import_ui, "setup_import")
            api$defaults(function(id) {
                ns$buildProtImportFormatSpecificOptionsUi(tool, shiny::NS(id))
            }, "setup_import")
            api$upload("setup_import-search_results_standard", "search_results")
            api$upload("setup_import-fasta_file_standard", "fasta")
            api$set(list(`setup_import-format_override` = tool))
            if (peptide_level) api$set(list(`setup_import-dia_route` = workflow))
            api$press("setup_import-process_data")
            # The organism modal is confirmed at its rendered values (A18); Cancel is never pressed.
            if (api$modal_open()) api$modal_confirm("Confirm Selection")
        })
    design <- valStep(
        "design", "design",  # nolint
        list(head = if (peptide_level) "raw_data_s4" else "protein_s4_initial"),
        module = "design_matrix", structural = "^design_matrix-builder-",
        press = "design_matrix-builder-save_results",
        run = function(api) {
            api$design(prefix = "design_matrix-builder-", tech_rep = FALSE,
                       match = "design_matrix-builder-selected_runs")
        })  # nolint
    norm <- valStep(
        "norm", "normalization",  # nolint
        list(head = if (limpa) c("peptide_count_filtered", "ruv_covariates") else
            c("normalized", "ruv_corrected")),  # nolint
        module = "normalization", press = "normalization-run_normalization",  # nolint
        run = function(api) {
            api$defaults(ns$mod_prot_norm_ui, "normalization")
            api$overrides()
            api$press("normalization-run_normalization")
        })
    correlation <- valStep(
        "correlation", "normalization", list(head = "correlation_filtered"),  # nolint
        module = "normalization", structural = "^normalization-norm_qc_tabs$",
        press = "normalization-apply_correlation_filter",
        run = function(api) {
            api$set(list(`normalization-norm_qc_tabs` = "Correlation Filtering"))
            api$press("normalization-apply_correlation_filter")
        })
    covariates <- if (limpa && !identical(alternative, "no_covariates")) {
        valStep("ruv_covariates", "normalization", list(head = "ruv_covariates"),
                module = "normalization", press = "normalization-apply_ruv_covariates",
                run = function(api) api$press("normalization-apply_ruv_covariates"))
    }
    export_load <- if (identical(alternative, "export_load")) {
        valExportLoadSteps(if (is.null(covariates)) "correlation_filtered" else "ruv_covariates",
                           "normalization", "normalization",
                           "normalization-export_filtered_session", "da",
                           "differential_expression", ns$mod_prot_da_ui)
    }
    da <- valStep(
        "da", "da", list(tab = c(differential_expression = "complete")),  # nolint
        module = "differential_expression", press = "differential_expression-run_da_analysis",
        run = function(api) {
            api$defaults(ns$mod_prot_da_ui, "differential_expression")
            api$press("differential_expression-run_da_analysis")
            # A confirmation the DA raises (the q-value warning) is its primary button.
            if (api$modal_open()) api$modal_confirm(NULL)
        })
    save <- if (!is.null(summary)) {
        valStep("summary_save", "session_summary", "diagnostic", module = "session_summary",
                structural = "^session_summary-(experiment_label|description)$",
                press = "session_summary-save_workflow_args",
                run = function(api) {
                    api$defaults(ns$mod_prot_summary_ui, "session_summary")
                    api$set(list(`session_summary-experiment_label` = summary$label,
                                 `session_summary-description` = summary$description))
                    api$press("session_summary-save_workflow_args")
                })
    }
    Filter(Negate(is.null), c(list(import, design), lapply(specs, qcStep),
                              list(norm, correlation, covariates), export_load, list(da, save)))
}

#' The lipidomics and metabolomics script (`short`: `lipid` or `metab`)
#' @noRd
valOmicsScript <- function(r, short, tool, assays, import_values, summary, alternative = NULL) {
    ns <- asNamespace("MultiScholaR")
    ui <- function(name) get(paste0("mod_", short, "_", name, "_ui"), envir = ns)
    state <- function(x) paste0(short, "_", x)
    import <- valStep(
        "import", "import", list(tab = c(setup_import = "complete")), module = "import",  # nolint
        structural = paste0("^import-(assay[12]_name|vendor_format|assay[12]_file_std|",
                            "metabolite_id_col_custom|annotation_col_custom|",
                            "sample_cols_pattern|is_pattern)$"),
        press = "import-process_import",  # nolint
        run = function(api) {
            api$defaults(ui("import"), "import")
            if (length(assays) >= 1L) api$set(list(`import-assay1_name` = assays[[1L]]))
            api$set(list(`import-vendor_format` = tool))
            api$upload("import-assay1_file_std", "assay1")
            if (api$has_file("assay2") && length(assays) >= 2L) {
                api$set(list(`import-assay2_name` = assays[[2L]]))
                api$upload("import-assay2_file_std", "assay2")
            }
            if (length(import_values)) api$set(import_values)
            api$press("import-process_import")
        })
    design <- valStep("design", "design", list(head = state("raw_data_s4")), module = "design",
                      structural = "^design-builder-", press = "design-builder-save_results",
                      run = function(api) {
                          api$design(prefix = "design-builder-", tech_rep = TRUE,
                                     match = "design-builder-selected_runs")
                      })
    intensity <- valStep("qc_intensity", "qc", list(head = state("intensity_filtered")),
                         module = "qc", press = "qc-intensity-apply_filter",
                         run = function(api) {
                             api$defaults(ui("qc_intensity"), "qc-intensity")
                             api$press("qc-intensity-apply_filter")
                         })
    duplicates <- valStep("qc_duplicates", "qc", list(head = state("duplicates_resolved")),
                          module = "qc",
                          press = c("qc-duplicates-detect_duplicates",
                                    "qc-duplicates-resolve_duplicates"),
                          run = function(api) {
                              api$defaults(ui("qc_duplicates"), "qc-duplicates")
                              api$press("qc-duplicates-detect_duplicates")
                              api$press("qc-duplicates-resolve_duplicates")
                          })
    itsd <- valStep("qc_itsd", "qc", "diagnostic", module = "qc", press = "qc-itsd-analyze_is",
                    run = function(api) {
                        api$defaults(ui("qc_itsd"), "qc-itsd")
                        api$press("qc-itsd-analyze_is")
                    })
    finalize <- valStep("qc_finalize", "qc", list(head = state("qc_complete")), module = "qc",
                        press = "qc-s4_finalize-finalize_qc",
                        run = function(api) {
                            api$defaults(ui("qc_s4"), "qc-s4_finalize")
                            api$press("qc-s4_finalize-finalize_qc")
                        })
    norm <- valStep("norm", "norm", list(head = c(state("normalized"), state("ruv_corrected"))),
                    module = "norm", press = "norm-run_normalization",
                    run = function(api) {
                        api$defaults(ui("norm"), "norm")
                        api$overrides()
                        api$press("norm-run_normalization")
                    })
    correlation <- valStep("correlation", "norm", list(head = state("correlation_filtered")),
                           module = "norm", press = "norm-apply_correlation_filter",
                           run = function(api) api$press("norm-apply_correlation_filter"))
    da <- valStep("da", "de", list(tab = c(differential_analysis = "complete")), module = "de",
                  press = "de-run_da_analysis",
                  run = function(api) {
                      api$defaults(ui("da"), "de")
                      api$press("de-run_da_analysis")
                  })
    save <- if (!is.null(summary)) valSummaryStep(ui("summary"), summary)
    export_load <- if (identical(alternative, "export_load")) {
        valExportLoadSteps(state("correlation_filtered"), "norm", "norm", "norm-export_session",
                           "de", "de", ui("da"))
    }
    Filter(Negate(is.null), c(list(import, design, intensity, duplicates, itsd, finalize, norm,
                                   correlation), export_load, list(da, save)))
}

#' The lipidomics, metabolomics and transcriptomics Summary save step
#' @noRd
valSummaryStep <- function(ui_fn, summary) {
    valStep("summary_save", "summary", "diagnostic", module = "summary",
            structural = "^summary-(experiment_label|description)$",
            press = "summary-save_workflow_args",
            run = function(api) {
                api$defaults(ui_fn, "summary")
                api$set(list(`summary-experiment_label` = summary$label,
                             `summary-description` = summary$description))
                api$press("summary-save_workflow_args")
            })
}

#' The transcriptomics script (the shared modules under the transcriptomics configuration, read
#' when the step runs, so that the cell's RUV design gate is in force)
#' @noRd
valTrxScript <- function(summary) {
    ns <- asNamespace("MultiScholaR")
    ui <- function(fn, group = NULL) {
        force(fn)
        function(id) {
            config <- ns$omicGuiConfig("transcriptomics")
            fn(id, if (is.null(group)) config else config[[group]])
        }
    }
    list_steps <- list(
        valStep("import", "import", list(tab = c(setup_import = "complete")), module = "import",
                structural = "^import-assay1_file_std$", press = "import-process_import",
                run = function(api) {
                    api$defaults(ui(ns$mod_omic_import_ui), "import")
                    api$upload("import-assay1_file_std", "counts")
                    api$press("import-process_import")
                }),
        valStep("design", "design", list(head = "trx_raw_data_s4"), module = "design",
                structural = "^design-builder-", press = "design-builder-save_results",
                run = function(api) {
                          api$design(prefix = "design-builder-", tech_rep = TRUE,  # nolint
                                     match = "design-builder-selected_runs")
                      }),  # nolint
        valStep("qc_filter", "qc", list(head = "trx_expression_filtered"), module = "qc",
                press = "qc-filter-apply_filter",
                run = function(api) {
                    api$defaults(ns$mod_trx_qc_filter_ui, "qc-filter")
                    api$press("qc-filter-apply_filter")
                }),
        valStep("qc_finalize", "qc", list(head = "trx_qc_complete"), module = "qc",
                press = "qc-s4_finalize-finalize_qc",
                run = function(api) {
                    api$defaults(ui(ns$mod_omic_qc_finalize_ui), "qc-s4_finalize")
                    api$press("qc-s4_finalize-finalize_qc")
                }),
        valStep("norm", "norm", list(head = "trx_ruvs"), module = "norm",
                press = "norm-run_normalization",
                run = function(api) {
                    api$defaults(ui(ns$mod_omic_norm_ui, "norm"), "norm")
                    api$overrides()
                    api$press("norm-run_normalization")
                }),
        valStep("da", "de", list(tab = c(differential_analysis = "complete")), module = "de",
                press = "de-run_da_analysis",
                run = function(api) {
                    api$defaults(ui(ns$mod_omic_da_ui, "da"), "de")
                    api$press("de-run_da_analysis")
                }),
        if (!is.null(summary)) valSummaryStep(ui(ns$mod_omic_summary_ui, "summary"), summary)
    )
    Filter(Negate(is.null), list_steps)
}
