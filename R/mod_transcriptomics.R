# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# KIT-096: the transcriptomics workflow module (import, design, QC). The tabs are the lipid import, design-builder and QC
# modules under the transcriptomics profile (`trxOmicProfile()`), which supplies the omic's collaborators at the seams those
# modules already take; nothing is copied. The tabs hold state references only (KIT-075) and run on the managed runner, or
# the binding's wrappers in process (KIT-078). Recorded divergence from conventions/R-shiny.md: the shinydashboard and
# wellPanel shell, the managed runner in place of ExtendedTask and the import's fileInput fallback are kept, to match the
# other omics (the user's instruction) under the existing-code scope rule of conventions/R.md.

# ---------------------------------------------------------------------------------------------------------------
# The profile.
# ---------------------------------------------------------------------------------------------------------------

#' The transcriptomics profile of the lipid import, design and QC modules: labels, test-id prefix, the count table's
#' format, the binding's readers and wrappers, and the collaborators each module takes at its seams.
#' @return A plain list.
#' @noRd
trxOmicProfile <- function() {
    list(omic_type = "transcriptomics", prefix = "trx", testid = "trx", label = "Transcriptomics", label_lower = "transcriptomics",
        class = "TranscriptomicsAssayData", workflow_type = "TRANSCRIPTOMICS", data_format = "trx_counts",
        import_title = "Transcriptomics Data Import", formats = c("Count table (gene ID, gene name, sample counts)" = "trx_counts"),
        assay_name = "counts", id_label = "Gene ID Column", annotation_label = "Gene Name Column",
        import_help = paste("One delimited count table: a gene ID column, a gene name column and one column per sample, named",
            "exactly as in the design. Counts must be whole numbers: a missing, negative or fractional count is refused, with",
            "no rounding. Sample names are kept as they are in the file."),
        load_preview = trxImportPreview, process_import = runTrxImportProcessing,
        design_help = paste("Samples are matched to the count columns by sample name only, in the counts file's order; there",
            "is no positional mapping. Each group needs at least two samples. The order in which you add the factor levels is",
            "the group level order the analysis uses."),
        run_order = identity, group_levels = trxDesignGroupLevels, validate_design = trxDesignCheck, build_object = trxDesignObject,
        record_provenance = recordTrxDesignProvenance, state_name = "trx_raw_data_s4",
        state_description = "Initial TranscriptomicsAssayData S4 object created after design matrix",
        qc_title = "Transcriptomics Quality Control & Filtering", qc_tabs = trxQcTabsUi, qc_init = initializeTrxQcSubmodules,
        read_metadata = function(manager) readWorkflowOmicsMetadata(transcriptomicsDeclarationBinding, manager),
        read_state = readTrxWorkflowState, register_qc_summary = registerTrxQcSummaryOutput, finalize = trxQcFinalizeCollaborators())
}

# ---------------------------------------------------------------------------------------------------------------
# Import: the count table, the count policy, the processing at the lipid import's seams.
# ---------------------------------------------------------------------------------------------------------------

#' The import tab's preview of a count table, in the shape of `loadLipidImportAssayPreview`: read as the notebook reads
#' it (vroom), the id and name columns proposed (`GeneID`, `GeneName` when present), every other column a sample.
#' @noRd
trxImportPreview <- function(assay1File) {
    data <- vroom::vroom(assay1File, show_col_types = FALSE, progress = FALSE)
    # The reader's lazy problems pointer is not a value the session can hold by reference.
    attr(data, "problems") <- NULL
    headers <- names(data)
    id <- if ("GeneID" %in% headers) "GeneID" else headers[[1L]]
    name <- if ("GeneName" %in% headers) "GeneName" else setdiff(headers, id)[[1L]]
    list(headers = headers, detectedFormat = "trx_counts", formatConfidence = 1,
        importResult = list(data = data, sample_columns = setdiff(headers, c(id, name))), assayData = data,
        updates = list(lipidId = list(choices = headers, selected = id), annotation = list(choices = headers, selected = name),
            isPattern = NULL))
}

#' KIT-093's count policy at import on every column but the id and name columns. A refusal keeps KIT-093's class and
#' names it in its message, which the tab shows.
#' @return `assayList`, unchanged, when it passes.
#' @noRd
checkTrxImportCounts <- function(assayList, idCol, annotationCol) {
    table <- assayList[[1L]]
    if (!lipidImportScalarString(annotationCol)) trxAbort("a gene name column is required (trx_counts_invalid)", c("trx_counts_invalid", "trx_error"))
    absent <- setdiff(c(idCol, annotationCol), names(table))
    if (length(absent)) {
        trxAbort(sprintf("the count table lacks the column(s) %s (trx_counts_invalid)", paste(absent, collapse = ", ")),
            c("trx_counts_invalid", "trx_error"))
    }
    withCallingHandlers(validateTrxCounts(as.data.frame(table)[setdiff(names(table), c(idCol, annotationCol))]),
        trx_counts_invalid = function(error) trxAbort(sprintf("%s (%s)", conditionMessage(error), class(error)[[1L]]),
            setdiff(class(error), c("rlang_error", "error", "condition"))))
    invisible(assayList)
}

#' The lipid import's validation under the count policy: the policy first (its class shown), then the mapping checks.
#' @noRd
trxValidateImport <- function(assayList, lipidIdCol, annotationCol, sampleColumns) {
    checkTrxImportCounts(assayList, lipidIdCol, annotationCol)
    validation <- validateLipidImportProcessingState(assayList, lipidIdCol, sampleColumns)
    if (!isTRUE(validation$valid)) {
        trxAbort(sprintf("Invalid transcriptomics import: %s (trx_counts_invalid)", paste(validation$errors, collapse = "; ")),
            c("trx_counts_invalid", "trx_error"))
    }
    validation
}

#' The import tab's processing (the lipid processing with the transcriptomics collaborators): one count table, its samples
#' every column but the id and name columns, no sanitisation, the owner bound to `TRANSCRIPTOMICS`.
#' @param ... The lipid request's other arguments (formats, pattern, sanitisation, second assay) are not used; the
#'   notification and log functions pass through.
#' @noRd
runTrxImportProcessing <- function(workflowData, assay1Name, assay1Data, lipidIdCol, annotationCol, assay1File = NULL,
                                   experimentPaths = NULL, ...) {
    passed <- list(...)
    passed <- passed[intersect(names(passed), c("notify", "removeNotify", "logError"))]
    sampleColumns <- setdiff(names(assay1Data), c(lipidIdCol, annotationCol))
    apply <- function(workflowData, ...) {
        applyLipidImportResultToWorkflow(workflowData, ..., workflowType = "TRANSCRIPTOMICS")
        workflowData$data_type <- "transcriptomics"
    }
    do.call(runLipidImportProcessing, c(list(workflowData = workflowData, assay1Name = assay1Name, assay1Data = assay1Data,
        assay1File = assay1File, assay2File = NULL, assay2Name = "", vendorFormat = "trx_counts", detectedFormat = "trx_counts",
        lipidIdCol = lipidIdCol, annotationCol = annotationCol, sampleColumns = sampleColumns, isPattern = "", sanitizeNames = FALSE,
        experimentPaths = experimentPaths,
        validateImportState = function(assayList, lipidIdCol, sampleColumns) trxValidateImport(assayList, lipidIdCol, annotationCol, sampleColumns),
        writeImportArtifacts = function(...) writeLipidImportSourceArtifacts(..., summaryFileName = "transcriptomics_import_summary.tsv"),
        applyResultToWorkflow = apply), passed))
}

# ---------------------------------------------------------------------------------------------------------------
# Design: the level order, KIT-093's validation at save, the design state.
# ---------------------------------------------------------------------------------------------------------------

#' The group level order the user entered: the assigned groups ordered by the positions of their factor levels in the
#' builder's factor list (the order they were added); a group of several factors by each factor in turn.
#' @noRd
trxDesignGroupLevels <- function(design, factors) {
    assigned <- !is.na(design$group) & design$group != ""
    groups <- unique(design$group[assigned])
    first <- design[match(groups, design$group), , drop = FALSE]
    position <- function(column) if (column %in% names(first)) match(first[[column]], factors) else rep(NA_integer_, length(groups))
    groups[order(position("factor1"), position("factor2"), position("factor3"), seq_along(groups))]
}

#' A refusal's text as the tab shows it: the message and its class.
#' @noRd
trxConditionText <- function(error) {
    text <- conditionMessage(error)
    class <- class(error)[[1L]]
    if (grepl(class, text, fixed = TRUE)) text else sprintf("%s (%s)", text, class)
}

#' KIT-093's design validation at save, by name only, against every imported count column except the runs the user
#' removed (an unassigned run is a count column outside the design), with at least two runs per group (RUVs and RUVIII
#' need replicates). Returns the lipid preflight's shape.
#' @noRd
trxDesignCheck <- function(results, columnMapping) {
    runs <- setdiff(columnMapping$sample_columns, results$excluded_runs)
    tryCatch({
        levels <- results$group_levels %||% trxDesignGroupLevels(results$design_matrix, character())
        checked <- validateTrxDesign(runs, results$design_matrix, levels)
        sizes <- table(checked$group)
        small <- names(sizes)[sizes < 2L]
        if (length(small)) {
            trxAbort(sprintf("group(s) with fewer than two samples: %s", paste(small, collapse = ", ")),
                c("trx_design_group_singleton", "trx_design_invalid", "trx_error"))
        }
        list(valid = TRUE, errors = character())
    }, trx_error = function(error) list(valid = FALSE, errors = trxConditionText(error)))
}

#' The design state's object: the builder's count table and design, the recorded level order, the chosen id and name
#' columns, and the builder's `tech_rep_group` as the technical replicates.
#' @noRd
trxDesignObject <- function(results, columnMapping) {
    createTranscriptomicsAssayData(results$data_cln, results$design_matrix, results$group_levels,
        roles = list(gene_id_column = columnMapping$lipid_id_col, annotation_id_column = columnMapping$annotation_col,
            technical_replicate_id = "tech_rep_group"))
}

# ---------------------------------------------------------------------------------------------------------------
# The module.
# ---------------------------------------------------------------------------------------------------------------

#' @title Transcriptomics Workflow Module
#' @description The transcriptomics workflow's import, design and QC tabs (KIT-096; normalisation, DA and the summary
#'   follow in KIT-097).
#' @param id Module id (`transcriptomics_workflow` in the app).
#' @name mod_transcriptomics
#' @noRd
mod_transcriptomics_ui <- function(id) {
    ns <- shiny::NS(id)
    profile <- trxOmicProfile()
    tab <- function(title, icon, value, ui) {
        shiny::tabPanel(title = shiny::tagList(shiny::icon(icon), title), value = value,
            shiny::div(`data-testid` = paste0("trx-tab-", value), ui))
    }
    shiny::tagList(shiny::fluidRow(shiny::column(12,
        shiny::h2(shiny::icon("chart-line"), " Transcriptomics Workflow", style = "margin-bottom: 20px;"),
        shiny::uiOutput(ns("workflow_progress")),
        testid(shiny::tabsetPanel(id = ns("transcriptomics_tabs"), type = "pills",
            tab(" Setup & Import", "file-import", "import", mod_lipid_import_ui(ns("import"), profile = profile)),
            tab(" Design Matrix", "th", "design", mod_lipid_design_ui(ns("design"), profile = profile)),
            tab(" Quality Control", "check-double", "qc", mod_lipid_qc_ui(ns("qc"), profile = profile))),
        "trx-workflow-tabs"))))
}

#' The bus of the transcriptomics tabs.
#' @noRd
newTrxWorkflowData <- function() {
    shiny::reactiveValues(data_tbl = NULL, data_format = NULL, data_type = "transcriptomics", column_mapping = NULL,
        config_list = list(), state_manager = NULL, contrasts = list(), contrasts_tbl = NULL,
        tab_status = list(setup_import = "pending", design_matrix = "disabled", quality_control = "disabled"),
        state_update_trigger = NULL, processing_log = list())
}

#' @rdname mod_transcriptomics
#' @param project_dirs,omic_type,experiment_label,volumes As `mod_lipidomics_server`.
#' @return The module's bus.
#' @noRd
mod_transcriptomics_server <- function(id, project_dirs, omic_type, experiment_label, volumes = NULL) {
    shiny::moduleServer(id, function(input, output, session) {
        workflow_data <- newTrxWorkflowData()
        registerWorkflowOmicsBus(session, omic_type, workflow_data)
        if (!omic_type %in% names(project_dirs)) {
            logger::log_error(paste("No directory information found for", omic_type))
            return()
        }
        experiment_paths <- project_dirs[[omic_type]]
        # The session's state owner from the start, bound to the route at the import; the managed runner when the
        # application is managed (mod_lipidomics.R's install).
        workflow_data$state_manager <- createWorkflowSessionOwner(workflow_data, experiment_paths, session, omic_type = omic_type)
        installWorkflowOmicsActions(session, workflow_data, omic_type, experiment_label, experiment_paths)
        profile <- trxOmicProfile()
        qc_trigger <- shiny::reactiveVal(NULL)
        mod_lipid_import_server("import", workflow_data, experiment_paths, volumes = NULL, profile = profile)
        mod_lipid_design_server("design", workflow_data, experiment_paths, volumes = NULL, qc_trigger = qc_trigger, profile = profile)
        mod_lipid_qc_server("qc", workflow_data, experiment_paths, omic_type, experiment_label, qc_trigger, profile = profile)
        workflow_data$active_tab <- shiny::reactive(input$transcriptomics_tabs)
        output$workflow_progress <- shiny::renderUI(render_workflow_stepper(list(
            list(name = "Import", key = "setup_import", icon = "file-import"), list(name = "Design", key = "design_matrix", icon = "th"),
            list(name = "QC", key = "quality_control", icon = "check-double")), workflow_data$tab_status))
        shiny::observeEvent(workflow_data$tab_status$design_matrix, {
            status <- workflow_data$tab_status
            if (identical(status$design_matrix, "complete") && identical(status$quality_control, "disabled")) {
                status$quality_control <- "pending"
                workflow_data$tab_status <- status
            }
        }, ignoreNULL = TRUE)
        workflow_data
    })
}
