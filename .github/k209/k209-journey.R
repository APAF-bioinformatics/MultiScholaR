# KIT-209 (branch run-kit209-memory): one in-process journey of one route on one corpus with one build, in a fresh R
# process alone in its own memory cgroup (k209-run-journey.sh). Memory discovery (VAL-151): it measures and reports;
# nothing here judges.
#
# Usage: Rscript --vanilla k209-journey.R <journey arguments (JSON file)>
#
# The step bodies are the in-process recorders' (kit-031 lipid, kit-032 metab, kit-099 trx, kit-118 DIA_LIMPA,
# kit-021 protein), each step in a function scope (R07): between steps the lane holds only the module-equivalent bus
# (`workflow`), the session owner, the experiment paths, the corpus paths and the bound settings (plus the globals the
# package itself reads: `project_dirs`, and `filtering_progress` the package writes). The DIA (IQ) body is new (R08):
# the DIA tabs' own in-process step helpers (mod_prot_qc_peptide_*, mod_prot_qc_protein_*), listed in the declaration.
# The protein QC steps of the DIA and protein routes run through the tabs' in-process step helpers too.
#
# Marks (declaration, "Endpoints"):
#   baseline   after the package, the harness and the arguments are loaded: gc(full = TRUE) twice, then gc(reset = TRUE);
#              R heap in use (Ncells + Vcells, MiB) and VmRSS; VmHWM reset to the current RSS (/proc/self/clear_refs, 5).
#   each step  at entry: gc(full = TRUE), heap in use (MiB) and VmRSS (the retained heap at the boundary), then
#              gc(reset = TRUE) and the VmHWM reset; at exit: VmHWM (MiB; the step's process peak, everything resident at
#              entry included), gc() max used (MiB), VmRSS, seconds, the processes in the journey's cgroup. Then, outside
#              the step's figures: the head's state, its rows and sizes, the store's bytes.
#   final      after the last step: gc(full = TRUE), heap in use and VmRSS.
# Files: <out>.baseline.json at the baseline; <out>.marks.json rewritten after every step (a kill leaves the steps
# done); <out>.json at the end (completed or error, with the frames of the error and their owners).

k209_started <- proc.time()[["elapsed"]]
k209_args <- jsonlite::read_json(commandArgs(trailingOnly = TRUE)[[1L]], simplifyVector = TRUE)
`%||%` <- function(a, b) if (is.null(a)) b else a
k209_rec <- new.env(parent = emptyenv())
k209_rec$marks <- list()

k209_status_mib <- function(field) {
    line <- grep(paste0("^", field, ":"), readLines("/proc/self/status"), value = TRUE)
    if (length(line) != 1L) return(NA_real_)
    round(as.numeric(sub("^[^0-9]*([0-9]+).*$", "\\1", line)) / 1024, 1)
}
k209_write <- function(path, value) {
    tmp <- paste0(path, ".tmp")
    jsonlite::write_json(value, tmp, auto_unbox = TRUE, digits = NA, null = "null", na = "null", pretty = TRUE)
    file.rename(tmp, path)
}
k209_cgroup_procs <- function() {
    line <- grep("^0::", readLines("/proc/self/cgroup"), value = TRUE)
    path <- file.path("/sys/fs/cgroup", sub("^0::/?", "", line), "cgroup.procs")
    if (length(line) != 1L || !file.exists(path)) return(NA_integer_)
    length(readLines(path))
}
k209_common <- function() {
    list(id = k209_args$id, route = k209_args$route, corpus = k209_args$corpus, build = k209_args$build,
        commit = k209_args$commit, unit = k209_args$unit, rep = k209_args$rep, order = k209_args$order,
        inject = k209_args$inject)
}

# --- Load: the build's library first ------------------------------------------------------------------------------
.libPaths(c(k209_args$lib, .libPaths()))
suppressPackageStartupMessages(loadNamespace("MultiScholaR"))
k209_ns <- asNamespace("MultiScholaR")
k209_loaded_from <- normalizePath(find.package("MultiScholaR"), mustWork = TRUE)
if (!identical(dirname(k209_loaded_from), normalizePath(k209_args$lib, mustWork = TRUE)))
    stop(sprintf("MultiScholaR loaded from %s, not from the build's library %s", k209_loaded_from, k209_args$lib))
k209_lib_commit <- trimws(readLines(file.path(k209_args$lib, "k209-build-commit.txt"), n = 1L))
if (!identical(k209_lib_commit, k209_args$commit))
    stop(sprintf("the build's library holds commit %s, the journey names %s", k209_lib_commit, k209_args$commit))
k209_api <- function(name) get(name, k209_ns)
# The harness's pure functions only (kitJourneySteps, kitFindStore, kitReadHead, kitFingerprintStore), from the
# candidate's checkout; tools/profiling/ is unchanged since b468a885.
k209_harness <- new.env(parent = globalenv())
sys.source(k209_args$harness, k209_harness)
k209_h <- new.env(parent = k209_ns)

# --- Settings bound identically on both builds (declaration, "Settings") -----------------------------------------
K209_LABEL <- "dia_baseline_smoke"          # bindings journey$experiment_label (protein routes)
K209_OMICS_LABEL <- "kit_omics_lifecycle"   # the driver's experiment label for the omics routes
K209_DA <- list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0)
K209_SETTINGS <- list(
    lipidomics = list(apply_itsd = TRUE, norm_method = "cyclicloess", ruv_mode = "automatic", correlation_threshold = 0.5,
        intensity_percentile = 10L, intensity_proportion = 0.5, itsd_aggregation = "median", log_offset = 1L),
    transcriptomics = list(ruv_mode = "automatic", auto_percentage_min = 5, auto_percentage_max = 30, auto_percentage_step = 5,
        k_penalty_weight = 0.05, da_q_value = 0.05, da_treat_lfc = 1, ruviii_m_source = "group"),
    DIA = list(norm_method = "cyclicloess", ruv_mode = "automatic", auto_percentage_min = 1, auto_percentage_max = 20,
        separation_metric = "max_difference", k_penalty_weight = 0.5, adaptive_k_penalty = TRUE, ruv_grouping_variable = "group",
        ruv_percentage = 5L, ruv_k = NA, correlation_threshold = 0.5, min_peptides_per_protein = 1, min_peptidoforms_per_protein = 2,
        peptide_intensity = list(strict = FALSE, min_reps = 2, min_groups = 2, percentile = 1), min_peptides_per_sample = 3000,
        replicate_group_column = "replicates", proportion_missing = 0.5),
    DIA_LIMPA = list(norm_method = "cyclicloess", ruv_mode = "covariates", correlation_threshold = 0.5, engine = "fastlimpadpc",
        dpc_slope = 0.8, chunk = 1000, num_peptides_per_protein_thresh = 1, num_peptidoforms_per_protein_thresh = 2,
        ruv_grouping_variable = "group"),
    protein = list(norm_method = "cyclicloess", ruv_mode = "manual", auto_percentage_min = 1, auto_percentage_max = 20,
        separation_metric = "max_difference", k_penalty_weight = 0.5, adaptive_k_penalty = TRUE, ruv_grouping_variable = "group",
        ruv_percentage = 5L, ruv_k = NA, correlation_threshold = 0.5))
K209_SETTINGS$metabolomics <- K209_SETTINGS$lipidomics
# The protein fixtures (S1, S2) bind RUV skip (bindings corpora lfq_fixture, tmt_fixture: 5 proteins).
K209_CORPUS_OVERRIDES <- list(lfq_fixture = list(ruv_mode = "skip"), tmt_fixture = list(ruv_mode = "skip"))

k209_settings <- function() {
    route <- k209_args$route
    if (identical(route, "calibration")) return(list())
    settings <- K209_SETTINGS[[route]]
    override <- K209_CORPUS_OVERRIDES[[k209_args$corpus]]
    if (!is.null(override)) settings <- utils::modifyList(settings, override)
    # G0a's injected package error (declaration, G0a item 1): an invalid bound setting the package must refuse. G0a-2
    # showed the between-sample kernel falls back to unnormalised data on an unknown method (no error), so the setting
    # is the intensity filter's percentile, left undefined (the omics journey clears it after the tab's defaults).
    if (identical(k209_args$inject, "package_error")) settings$intensity_percentile_injected <- "undefined (NULL)"
    settings
}

# --- The routes' step bodies ------------------------------------------------------------------------------------
# Each constructor returns an ordered list of steps (id, fn). `lane` holds the module-equivalent state.
k209_quiet <- function(...) invisible(NULL)

k209_project <- function(lane, omic, label, keyed) {
    base <- file.path(lane$work, "project", if (keyed) label else K209_OMICS_LABEL)
    project_dirs <- setupDirectories(base_dir = base, omic_types = omic, label = if (keyed) label else NULL, force = TRUE,
        reuse_existing = TRUE)
    assign("project_dirs", project_dirs, envir = globalenv())
    lane$project_dirs <- project_dirs
    lane$paths <- project_dirs[[if (keyed) 1L else omic]]
    lane$paths$time_dir <- lane$paths$time_dir %||% file.path(lane$paths$publication_graphs_dir, "filtering_qc")
    dir.create(lane$paths$time_dir, recursive = TRUE, showWarnings = FALSE)
}

k209_fake_session <- function() {
    session <- new.env(parent = emptyenv())
    session$userData <- new.env(parent = emptyenv())
    session$onFlush <- session$onFlushed <- session$onSessionEnded <- function(...) function() invisible(NULL)
    session
}

# The summary tab's copy and report through the shared core (mod_omic_summary.R at both builds), as kit-118's report step.
k209_report <- function(lane, omic, label) {
    summary_config <- omicGuiConfig(omic)$summary
    summary_dirs <- stats::setNames(list(lane$paths), omic)
    copied <- copyWorkflowSummaryResults(summary_config, lane$workflow, summary_dirs, omic, label, "kit-209 journey",
        argsSaved = FALSE, notify = k209_quiet, catFn = k209_quiet)
    if (!is.null(copied$error)) stop("the summary copy failed: ", copied$error)
    report <- renderWorkflowSummaryReport(summary_config, lane$workflow, summary_dirs, omic, label, notify = k209_quiet,
        progress = k209_quiet, exportRecordFn = k209_quiet, catFn = k209_quiet)
    if (!is.character(report$path) || length(report$path) != 1L || !file.exists(report$path)) stop("no report was rendered")
    k209_rec$report <- list(template = report$template, bytes = file.size(report$path))
    invisible(NULL)
}

# Lipidomics, metabolomics and transcriptomics: the real omic module server under shiny::testServer (KIT-075's residency
# journey, tests/testthat/helper-workflow-omics-residency.R), driven with the kit driver's input ids and the bound settings
# (tools/profiling/workflow_kit_driver.R, runKitOmicsStep and runKitTrxStep). The kit-031, 032 and 099 recorders call
# import and design helpers that exist at neither build (the shared-module rework, KIT-102 to KIT-117), so this is the
# routes' in-process GUI path at both builds. Each tab's UI defaults are set before its action (testServer renders no UI).
# The UI-defaults reader is KIT-075's omicsUiInputDefaults, copied.
k209_ui_defaults <- function(ui) {
    out <- list()
    attr_of <- function(tag, name) tag$attribs[[name]]
    visit <- function(tag) {
        if (inherits(tag, "shiny.tag.list") || (is.list(tag) && !inherits(tag, "shiny.tag"))) {
            for (child in tag) visit(child)
            return(invisible())
        }
        if (!inherits(tag, "shiny.tag")) return(invisible())
        id <- attr_of(tag, "id")
        if (identical(tag$name, "input")) {
            type <- attr_of(tag, "type") %||% "text"
            if (identical(type, "checkbox") && !is.null(id)) {
                out[[id]] <<- !is.null(attr_of(tag, "checked"))
            } else if (identical(type, "radio")) {
                if (!is.null(attr_of(tag, "checked"))) out[[attr_of(tag, "name")]] <<- attr_of(tag, "value")
            } else if (!is.null(id) && !is.null(attr_of(tag, "data-from"))) {
                to <- attr_of(tag, "data-to")
                out[[id]] <<- as.numeric(c(attr_of(tag, "data-from"), if (identical(attr_of(tag, "data-type"), "double")) to))
            } else if (!is.null(id) && !is.null(attr_of(tag, "value"))) {
                value <- attr_of(tag, "value")
                out[[id]] <<- if (identical(type, "number")) as.numeric(value) else value
            }
        } else if (identical(tag$name, "select") && !is.null(id)) {
            html <- paste(vapply(tag$children, function(child) paste(as.character(child), collapse = ""), character(1)), collapse = "")
            options <- regmatches(html, gregexpr("<option value=\"[^\"]*\"( selected)?", html))[[1]]
            chosen <- c(grep(" selected$", options, value = TRUE), options)
            if (length(chosen)) out[[id]] <<- sub("^<option value=\"([^\"]*)\".*$", "\\1", chosen[[1]])
        }
        for (child in tag$children) visit(child)
        invisible()
    }
    visit(ui)
    out
}

k209_omics_journey <- function(lane, corpus, settings) {
    omic <- k209_args$route
    trx <- identical(omic, "transcriptomics")
    short <- switch(omic, lipidomics = "lipid", metabolomics = "metab", transcriptomics = "trx")
    module_id <- paste0(omic, "_workflow")
    k209_project(lane, omic, K209_OMICS_LABEL, keyed = FALSE)
    server <- get(paste0("mod_", omic, "_server"), envir = k209_ns)
    # Every input default the omic's UI renders statically, by the module-relative id.
    all_defaults <- k209_ui_defaults(get(paste0("mod_", omic, "_ui"), envir = k209_ns)(module_id))
    names(all_defaults) <- sub(paste0("^", module_id, "-"), "", names(all_defaults))
    lane$defaults_count <- length(all_defaults)
    design_in <- utils::read.delim(corpus$design, stringsAsFactors = FALSE, check.names = FALSE)
    shiny::testServer(server, args = list(id = module_id, project_dirs = lane$project_dirs, omic_type = omic,
        experiment_label = K209_OMICS_LABEL), {
        lane$workflow <- workflow_data
        set <- function(values) do.call(session$setInputs, values)
        tab <- function(value) set(stats::setNames(list(value), paste0(omic, "_tabs")))
        # A tab's defaults: its own UI function rendered under the tab's id (KIT-075's `defaults(ui_name, id)`; the QC
        # sub-tabs render dynamically, so the omic's static UI does not hold them), else the omic's static UI.
        tab_ui <- c(`qc-intensity` = "qc_intensity", `qc-duplicates` = "qc_duplicates", `qc-itsd` = "qc_itsd",
            `qc-s4_finalize` = "qc_s4", norm = "norm", de = "da", summary = "summary", `qc-filter` = "qc_filter")
        defaults <- function(prefix) {
            ui_name <- paste0("mod_", short, "_", tab_ui[[prefix]], "_ui")
            chosen <- if (exists(ui_name, envir = k209_ns, inherits = FALSE)) {
                k209_ui_defaults(get(ui_name, envir = k209_ns)(prefix))
            } else {
                all_defaults[startsWith(names(all_defaults), paste0(prefix, "-"))]
            }
            lane$defaults_set[[prefix]] <- names(chosen)
            if (length(chosen)) set(chosen)
        }
        click <- function(id) set(stats::setNames(list(1L), id))
        builder <- function(name) paste0("design-builder-", name)
        groups <- unique(design_in$group)
        import <- function() {
            import_config <- omicGuiConfig(omic)
            if (trx) {
                preview <- omicImportFeatureTablePreview(corpus$counts, import_config)
                request <- list(assay1_name = "counts", assay1_data = preview$assay_data, assay1_file = corpus$counts,
                    vendor_format = "trx_counts", detected_format = preview$detected_format, id_col = preview$updates$id$selected,
                    annotation_col = preview$updates$annotation$selected, experiment_paths = lane$paths)
            } else {
                preview <- omicImportFeatureTablePreview(corpus$assays[[1L]], import_config)
                exclude <- if (identical(omic, "lipidomics")) OMICS_IMPORT_NORMALISED_COLUMNS
                request <- list(assay1_name = names(corpus$assays)[[1L]], assay1_data = preview$assay_data,
                    assay1_file = corpus$assays[[1L]], assay2_file = corpus$assays[[2L]], assay2_name = names(corpus$assays)[[2L]],
                    vendor_format = corpus$vendor, detected_format = preview$detected_format, id_col = preview$updates$id$selected,
                    annotation_col = preview$updates$annotation$selected,
                    sample_columns = resolveOmicsImportSampleColumns(preview$assay_data, corpus$vendor, "", preview$import_result, exclude),
                    is_pattern = preview$updates$is_pattern %||% "", sanitize_names = FALSE, experiment_paths = lane$paths)
            }
            preview <- NULL
            result <- processOmicImportFeatureTable(workflow_data, request, import_config,
                complete = function(...) completeOmicImport(workflow_data, import_config, notify = FALSE),
                showNotification = k209_quiet, removeNotification = k209_quiet, logInfo = k209_quiet)
            if (!identical(result$status, "success")) stop("the import failed: ", conditionMessage(result$error))
            session$flushReact()
            invisible(NULL)
        }
        design <- function() {
            tab("design")
            for (i in seq_along(groups)) {
                set(stats::setNames(list(groups[[i]]), builder("new_factor")))
                set(stats::setNames(list(i), builder("add_factor")))
            }
            for (i in seq_along(groups)) {
                set(stats::setNames(list(design_in$sample[design_in$group == groups[[i]]], groups[[i]], "", "", 1L),
                    builder(c("selected_runs", "factor1_select", "factor2_select", "factor3_select", "replicate_start"))))
                set(stats::setNames(list(i), builder("assign_metadata")))
            }
            clicks <- 0L
            for (group in groups) for (replicate in unique(design_in$replicates[design_in$group == group])) {
                samples <- design_in$sample[design_in$group == group & design_in$replicates == replicate]
                if (length(samples) < 2L) next
                clicks <- clicks + 1L
                set(stats::setNames(list(samples, "manual", as.integer(replicate)),
                    builder(c("tech_rep_samples", "tech_rep_assignment_mode", "manual_replicate_number"))))
                set(stats::setNames(list(clicks), builder("assign_tech_reps")))
            }
            set(stats::setNames(list(groups[[2L]], groups[[1L]], "~ 0 + group"), builder(c("contrast_group1", "contrast_group2", "formula_string"))))
            click(builder("add_contrast"))
            click(builder("save_results"))
            invisible(NULL)
        }
        report <- function() {
            tab("summary")
            defaults("summary")
            set(list(`summary-experiment_label` = K209_OMICS_LABEL, `summary-description` = "KIT-209 journey"))
            click("summary-save_workflow_args")
            click("summary-copy_to_publication")
            click("summary-generate_report")
            invisible(NULL)
        }
        steps <- if (trx) list(
            list(id = "import", fn = import),
            list(id = "design", fn = design),
            list(id = "expression_filter", fn = function() { tab("qc"); defaults("qc-filter"); click("qc-filter-apply_filter") }),
            list(id = "qc_finalize", fn = function() { defaults("qc-s4_finalize"); click("qc-s4_finalize-finalize_qc") }),
            list(id = "normalisation", fn = function() {
                tab("norm"); defaults("norm")
                set(list(`norm-ruv_mode` = settings$ruv_mode, `norm-auto_percentage_min` = settings$auto_percentage_min,
                    `norm-auto_percentage_max` = settings$auto_percentage_max, `norm-auto_percentage_step` = settings$auto_percentage_step,
                    `norm-k_penalty_weight` = settings$k_penalty_weight))
                click("norm-run_normalization")
            }),
            list(id = "da", fn = function() {
                tab("de"); defaults("de")
                set(list(`de-da_q_val_thresh` = settings$da_q_value, `de-treat_lfc_cutoff` = settings$da_treat_lfc))
                click("de-run_da_analysis")
            }),
            list(id = "ruviii_export", fn = function() {
                tab("summary"); defaults("summary")
                set(list(`summary-ruviii_m_source` = settings$ruviii_m_source))
                click("summary-run_ruviii_export")
            }),
            list(id = "report", fn = report),
            list(id = "session_export", fn = function() click("summary-export_session_state"))) else list(
            list(id = "import", fn = import),
            list(id = "design", fn = design),
            list(id = "intensity", fn = function() {
                tab("qc"); defaults("qc-intensity")
                if (identical(k209_args$inject, "package_error")) set(list(`qc-intensity-intensity_cutoff_percentile` = NULL))
                click("qc-intensity-apply_filter")
            }),
            list(id = "duplicates", fn = function() {
                defaults("qc-duplicates"); click("qc-duplicates-detect_duplicates"); click("qc-duplicates-resolve_duplicates")
            }),
            list(id = "qc_finalize", fn = function() {
                # The driver's internal-standard analysis (commits nothing), then the finalise.
                defaults("qc-itsd"); set(list(`qc-itsd-is_pattern` = "^IS_")); click("qc-itsd-analyze_is")
                defaults("qc-s4_finalize"); click("qc-s4_finalize-finalize_qc")
            }),
            list(id = "normalisation", fn = function() {
                tab("norm"); defaults("norm")
                set(list(`norm-apply_itsd` = isTRUE(settings$apply_itsd), `norm-norm_method` = settings$norm_method,
                    `norm-ruv_mode` = settings$ruv_mode, `norm-auto_percentage_min` = 1, `norm-auto_percentage_max` = 20,
                    `norm-separation_metric` = "max_difference", `norm-k_penalty_weight` = 0.5, `norm-adaptive_k_penalty` = TRUE,
                    `norm-ruv_grouping_variable` = "group"))
                click("norm-run_normalization")
            }),
            list(id = "correlation", fn = function() {
                set(list(`norm-min_pearson_correlation_threshold` = settings$correlation_threshold))
                click("norm-apply_correlation_filter")
            }),
            list(id = "session_export", fn = function() click("norm-export_session")),
            list(id = "da", fn = function() {
                tab("de"); defaults("de")
                click("de-load_filtered_session")
                set(list(`de-formula_string` = K209_DA$formula_string, `de-da_q_val_thresh` = K209_DA$da_q_val_thresh,
                    `de-treat_lfc_cutoff` = K209_DA$treat_lfc_cutoff))
                click("de-run_da_analysis")
            }),
            list(id = "report", fn = report))
        k209_measure(steps, lane)
    })
    invisible(NULL)
}

# The protein-side import and design shared by DIA, DIA_LIMPA and protein (kit-021 and kit-118 bodies).
k209_prot_import <- function(lane, corpus, method) {
    function() {
        k209_project(lane, "proteomics", K209_LABEL, keyed = TRUE)
        workflow <- new.env(parent = emptyenv())
        workflow$state_manager <- WorkflowState$new(audit_enabled = FALSE)
        workflow$tab_status <- list(); workflow$processing_log <- list(); workflow$qc_params <- list()
        lane$workflow <- workflow
        lane$session <- k209_fake_session()
        initializeDIANativeImport(workflow, lane$paths, lane$session, corpus$format, method)
        manager <- workflow$state_manager
        if (!isDIANativeSession(manager)) stop("the import did not open a native session")
        withDIAStateScope(manager, function() {
            inputs <- c(corpus$report, corpus$fasta)
            names(inputs) <- c(paste0(corpus$format, "_report"), "fasta")
            registerDIAInputFiles(manager, inputs)
            imported <- switch(corpus$format, diann = importDIANNData(corpus$report),
                maxquant = importMaxQuantData(corpus$report, use_lfq = TRUE, filter_contaminants = TRUE),
                pd_tmt = importProteomeDiscovererTMTData(corpus$report), stop("format ", corpus$format))
            applyProtImportResultToWorkflow(workflow, imported, corpus$format, corpus$fasta, sanitizeNames = TRUE,
                logInfo = k209_quiet, showNotification = k209_quiet)
            imported <- NULL
            releaseDIAWorkingReference(manager, workflow$data_tbl)
            processProtImportFastaData(workflowData = workflow, fastaPath = corpus$fasta, organismName = "Homo sapiens",
                experimentPaths = lane$paths, logInfo = k209_quiet, logWarn = k209_quiet)
            loadProtImportConfigurationResources(workflowData = workflow, experimentPaths = lane$paths, useShinyFiles = FALSE,
                readConfig = readProtGuiConfig, updateStatus = k209_quiet, messageFn = k209_quiet, showNotification = k209_quiet,
                logInfo = k209_quiet, logError = k209_quiet)
            finalizeProtImportSetupState(workflow, getDIAImportMetadata(manager), corpus$format, basename(corpus$report),
                basename(corpus$fasta), 9606, "Homo sapiens")
            completeDIAInputFiles(manager)
        })
        invisible(NULL)
    }
}

k209_prot_design <- function(lane, corpus) {
    function() {
        workflow <- lane$workflow
        manager <- workflow$state_manager
        withDIAStateScope(manager, function() {
            design_in <- utils::read.delim(corpus$design, stringsAsFactors = FALSE, check.names = FALSE)
            runs <- gtools::mixedsort(as.character(getDIAImportMetadata(manager)$samples$Run))
            raw <- tibble::tibble(Run = runs, group = NA_character_, factor1 = NA_character_, factor2 = NA_character_,
                factor3 = NA_character_, batch = NA_character_, replicates = NA_integer_, tech_reps = NA_integer_)
            for (group in unique(design_in$group)) {
                selected <- design_in$sample[design_in$group == group]
                idx <- which(raw$Run %in% selected)
                if (length(idx) != length(selected)) stop("design samples not among the imported runs")
                raw$factor1[idx] <- group; raw$factor2[idx] <- ""; raw$factor3[idx] <- NA_character_
                raw$replicates[idx] <- seq.int(1L, length.out = length(selected)); raw$group[idx] <- group
            }
            # The driver's technical-replicate actions (kit-118's body; the builder's assign_tech_reps).
            if ("replicates" %in% names(design_in) && !identical(k209_args$route, "protein")) {
                for (group in unique(design_in$group)) for (replicate in unique(design_in$replicates[design_in$group == group])) {
                    samples <- design_in$sample[design_in$group == group & design_in$replicates == replicate]
                    if (length(samples) < 2L) next
                    selected <- which(raw$Run %in% samples)
                    original <- unique(raw$replicates[selected])
                    consolidated <- length(original) - 1
                    highest <- max(original, na.rm = TRUE)
                    adjust <- which(raw$group == group & !is.na(raw$group) & raw$replicates > highest & !raw$Run %in% samples)
                    if (length(adjust)) raw$replicates[adjust] <- raw$replicates[adjust] - consolidated
                    raw$replicates[selected] <- as.integer(replicate)
                    raw$tech_reps[selected] <- stats::setNames(seq_along(samples), samples)[raw$Run[selected]]
                }
                if (any(!is.na(raw$tech_reps))) names(raw$Run) <- raw$Run
            }
            design_groups <- unique(design_in$group)
            contrast_data <- data.frame(contrast_name = paste0(design_groups[[2L]], ".vs.", design_groups[[1L]]),
                numerator = design_groups[[2L]], denominator = design_groups[[1L]], stringsAsFactors = FALSE)
            formula <- workflow$config_list$deAnalysisParameters$formula_string %||% "~ 0 + group"
            results <- buildProtDesignSaveResultsPayload(designMatrix = raw, currentRemovedSamples = character(0),
                dataCln = workflow$data_cln, contrastData = contrast_data, configList = workflow$config_list, formulaString = formula)
            results$data_reference <- saveDIAOwnedValue(results$data_cln, "design-assigned")
            results["data_cln"] <- list(NULL)
            hydrateProtDesignBuilderResults(results, workflow)
            persistProtDesignBuilderArtifacts(results = results, workflowData = workflow, sourceDir = lane$paths$source_dir)
            buildProtDesignStateCheckpoint(workflowData = workflow, workflowType = manager$workflow_type, actionLabel = "Save Design")
        })
        invisible(NULL)
    }
}

# Each tab's QC outputs (the plots the filtering progress panel publishes), as kit-118 publishes them.
k209_mirror_progress <- function(lane) lane$workflow$filtering_progress <- get0("filtering_progress", envir = globalenv(), inherits = FALSE)
k209_publish_peptide <- function(lane, result, step_name) {
    publishDIAQvalueOutput(lane$workflow$state_manager, result, list(data = result$filteringSummary %||% result$filteredS4@peptide_data,
        step_name = step_name, omic_type = "proteomics", experiment_label = K209_LABEL, return_grid = TRUE, overwrite = TRUE,
        project_dirs = lane$project_dirs), updateProteinFiltering)
    k209_mirror_progress(lane)
    invisible(NULL)
}
k209_publish_protein <- function(lane, result, step_name) {
    publishDIAProteinProgress(result, step_name, "proteomics", K209_LABEL, updateProteinFiltering)
    k209_mirror_progress(lane)
    invisible(NULL)
}
k209_scoped <- function(lane, expr_fn) function() withDIAStateScope(lane$workflow$state_manager, expr_fn)

# The protein QC tabs' in-process step helpers (mod_prot_qc_protein_*), each publishing its progress.
k209_protein_qc <- function(lane) {
    w <- function() lane$workflow
    list(
        list(id = "accession", fn = k209_scoped(lane, function() k209_publish_protein(lane,
            runProteinAccessionCleanupStep(w(), ";", "mean"), "10_protein_accession_cleaned"))),
        list(id = "protein_intensity", fn = k209_scoped(lane, function() k209_publish_protein(lane,
            runProteinIntensityFilterApplyStep(w(), FALSE, 2, 2, 1), "11_protein_intensity_filtered"))),
        list(id = "duplicates", fn = k209_scoped(lane, function() k209_publish_protein(lane,
            runProteinDuplicateRemovalStep(w(), "mean"), "12_duplicates_removed"))),
        list(id = "protein_replicate", fn = k209_scoped(lane, function() k209_publish_protein(lane,
            runProteinReplicateFilterApplyStep(w(), lane$paths, "group", 4), "13_protein_replicate_filtered"))))
}

k209_aesthetics <- function() list(color_var = "group", shape_var = "group")

# Normalisation and RUV, correlation, session export, DA and report on the classic routes (kit-021's body).
k209_prot_tail <- function(lane, settings) {
    grouping <- function() getProtNormRuvGroupingVariable(settings$ruv_grouping_variable %||% "group")
    norm_settings <- list(norm_method = settings$norm_method, ruv_mode = settings$ruv_mode,
        auto_percentage_min = settings$auto_percentage_min, auto_percentage_max = settings$auto_percentage_max,
        separation_metric = settings$separation_metric, k_penalty_weight = settings$k_penalty_weight,
        adaptive_k_penalty = settings$adaptive_k_penalty, ruv_percentage = settings$ruv_percentage, ruv_k = settings$ruv_k,
        color_variable = NULL, shape_variable = NULL)
    list(
        list(id = "normalisation_ruv", fn = function() {
            lane$norm_data <- new.env(parent = emptyenv())
            invisible(runDIAReferenceNormalization(norm_settings, lane$workflow, lane$norm_data, lane$paths, "proteomics", K209_LABEL,
                function(reference) NULL, function(reference) NULL, grouping))
        }),
        list(id = "correlation", fn = function() invisible(runDIAReferenceCorrelation(lane$workflow, lane$norm_data, lane$paths,
            settings$correlation_threshold, FALSE, "proteomics", K209_LABEL, grouping, k209_aesthetics))),
        list(id = "session_export", fn = function() invisible(runProtNormExportSessionWorkflow(lane$workflow, lane$norm_data,
            list(norm_method = settings$norm_method, ruv_mode = settings$ruv_mode), lane$paths$source_dir,
            withProgressFn = function(message, value, expr) force(expr), incProgressFn = k209_quiet, messageFn = k209_quiet))),
        k209_prot_da(lane),
        list(id = "report", fn = function() k209_report(lane, "proteomics", K209_LABEL)))
}

k209_prot_da <- function(lane) {
    list(id = "da", fn = k209_scoped(lane, function() {
        workflow <- lane$workflow
        if (!is.data.frame(workflow$contrasts_tbl) || nrow(workflow$contrasts_tbl) == 0L) stop("no contrasts")
        da_data <- new.env(parent = emptyenv())
        da_data$current_s4_object <- getDIAStateReference(workflow$state_manager)
        da_data$contrasts_available <- workflow$contrasts_tbl$contrasts
        invisible(runDIAReferenceDaAnalysis(K209_DA, da_data, workflow, lane$paths, contrasts = workflow$contrasts_tbl,
            annotations = workflow$uniprot_dat_cln, publish = TRUE))
    }))
}

# DIA (IQ): new body (R08). The tabs' own in-process step helpers, with the arguments their observers pass.
k209_route_dia <- function(lane, corpus, settings) {
    w <- function() lane$workflow
    cut <- .diaNnIdentificationQvalueCutoff
    pint <- settings$peptide_intensity
    c(list(
        list(id = "import", fn = k209_prot_import(lane, corpus, "DIA")),
        list(id = "design", fn = k209_prot_design(lane, corpus)),
        list(id = "qvalue", fn = k209_scoped(lane, function() k209_publish_peptide(lane,
            runPeptideQvalueApplyStep(w(), cut, cut, TRUE), "2_qval_filtered"))),
        list(id = "protein_peptide", fn = k209_scoped(lane, function() k209_publish_peptide(lane,
            runProteinPeptideApplyStep(w(), settings$min_peptides_per_protein, settings$min_peptidoforms_per_protein), "3_protein_peptide_filtered"))),
        list(id = "rollup", fn = k209_scoped(lane, function() k209_publish_peptide(lane,
            runDIAReferencePrecursorRollup(w()), "4_precursor_rollup"))),
        list(id = "intensity", fn = k209_scoped(lane, function() k209_publish_peptide(lane,
            runPeptideIntensityApplyStep(w(), pint$strict, pint$min_reps, pint$min_groups, pint$percentile), "5_intensity_filtered"))),
        list(id = "sample", fn = k209_scoped(lane, function() k209_publish_peptide(lane,
            runPeptideSampleApplyStep(w(), settings$min_peptides_per_sample), "6_sample_filtered"))),
        list(id = "replicate", fn = k209_scoped(lane, function() k209_publish_peptide(lane,
            runPeptideReplicateApplyStep(w(), settings$replicate_group_column), "7_replicate_filtered"))),
        list(id = "imputation", fn = k209_scoped(lane, function() k209_publish_peptide(lane,
            runPeptideImputationStep(w(), settings$proportion_missing), "8_imputed"))),
        list(id = "iq", fn = k209_scoped(lane, function() {
            result <- runProteinIqRollupApplyStep(w(), lane$paths, showNotificationFn = k209_quiet)
            publishDIAProteinProgress(result, "9_protein_s4_created", "proteomics", K209_LABEL, updateProteinFiltering, result$proteinObj)
            k209_mirror_progress(lane)
            invisible(NULL)
        }))), k209_protein_qc(lane), k209_prot_tail(lane, settings))
}

# DIA_LIMPA (kit-118 body); the engine bound explicitly (R02), refused when its namespace is missing.
k209_route_limpa <- function(lane, corpus, settings) {
    w <- function() lane$workflow
    grouping <- function() getProtNormRuvGroupingVariable(settings$ruv_grouping_variable)
    ruv_qc <- function(reference) {
        lane$norm_data$qc_plot_paths <- generateProtNormRuvCorrectedQcArtifacts(ruvCorrectedS4 = reference,
            qcDir = lane$norm_data$publication_paths$protein_qc_dir %||% lane$paths$protein_qc_dir, aesthetics = k209_aesthetics(),
            qcPlotPaths = lane$norm_data$qc_plot_paths)
    }
    list(
        list(id = "import", fn = k209_prot_import(lane, corpus, "DIA_LIMPA")),
        list(id = "design", fn = k209_prot_design(lane, corpus)),
        list(id = "qvalue", fn = k209_scoped(lane, function() {
            q <- w()$config_list$srlQvalueProteotypicPeptideClean
            k209_publish_peptide(lane, runDIAReferenceQvalue(w(), as.numeric(q$qvalue_threshold %||% 0.01),
                as.numeric(q$global_qvalue_threshold %||% 0.01), as.logical(as.integer(q$choose_only_proteotypic_peptide %||% 1L)),
                as.numeric(q$global_pg_qvalue_threshold %||% 0.01), q$confidence_provenance_mode %||% "diann_main_report", Sys.time),
                "2_qval_filtered")
        })),
        list(id = "rollup", fn = k209_scoped(lane, function() k209_publish_peptide(lane, runDIAReferencePrecursorRollup(w()), "4_precursor_rollup"))),
        list(id = "peptide_log2", fn = k209_scoped(lane, function() k209_publish_peptide(lane, runDIAReferencePeptideLog2Matrix(w()), "4_log2_matrix"))),
        list(id = "peptide_accession", fn = k209_scoped(lane, function() k209_publish_peptide(lane,
            runDIAReferencePeptideAccessionCleanup(w(), list(delimiter = ";", aggregation_method = "mean")), "7_annotation_cleanup"))),
        list(id = "dpc_quant", fn = k209_scoped(lane, function() k209_publish_protein(lane, suppressMessages(runDIAReferenceProteinDpcQuant(w(),
            list(engine = settings$engine, dpc_slope = settings$dpc_slope, chunk = settings$chunk))), "9_protein_dpc_quant"))),
        list(id = "protein_normalisation_limpa", fn = k209_scoped(lane, function() k209_publish_protein(lane,
            runDIAReferenceProteinNormalisation(w(), list(norm_method = settings$norm_method)), "10_protein_normalised"))),
        list(id = "peptide_count", fn = k209_scoped(lane, function() k209_publish_protein(lane, runDIAReferenceProteinPeptideCountFilter(w(),
            list(num_peptides_per_protein_thresh = settings$num_peptides_per_protein_thresh,
                num_peptidoforms_per_protein_thresh = settings$num_peptidoforms_per_protein_thresh)), "11_peptide_count_filtered"))),
        # The normalisation tab's action on this route is the pass-through (kit-118 ran it between steps; here it opens
        # the correlation step's measurement).
        list(id = "correlation", fn = function() {
            lane$norm_data <- new.env(parent = emptyenv())
            lane$norm_data$qc_plot_paths <- NULL
            invisible(runDIAReferenceLimpaNormalisationPassThrough(list(ruv_mode = settings$ruv_mode), w(), lane$norm_data, k209_quiet))
            invisible(runDIAReferenceCorrelation(w(), lane$norm_data, lane$paths, settings$correlation_threshold, FALSE, "proteomics",
                K209_LABEL, grouping, k209_aesthetics))
            lane$norm_data$correlation_filtering_complete <- TRUE
            invisible(NULL)
        }),
        list(id = "ruv_covariates", fn = function() invisible(runDIAReferenceRuvCovariates(list(ruv_mode = settings$ruv_mode,
            color_variable = NULL, shape_variable = NULL), w(), lane$norm_data, lane$paths, "proteomics", K209_LABEL, ruv_qc, grouping))),
        list(id = "session_export", fn = function() invisible(runProtNormExportSessionWorkflow(w(), lane$norm_data,
            list(norm_method = settings$norm_method, ruv_mode = settings$ruv_mode), lane$paths$source_dir,
            withProgressFn = function(message, value, expr) force(expr), incProgressFn = k209_quiet, messageFn = k209_quiet))),
        k209_prot_da(lane),
        list(id = "report", fn = function() k209_report(lane, "proteomics", K209_LABEL)))
}

# Protein-level LFQ and TMT (kit-021 body; the QC tabs' in-process helpers; manual RUV bound).
k209_route_protein <- function(lane, corpus, settings) {
    c(list(list(id = "import", fn = k209_prot_import(lane, corpus, "DIA")),
        list(id = "design", fn = k209_prot_design(lane, corpus))),
        k209_protein_qc(lane), k209_prot_tail(lane, settings))
}

# The reset control (G0a item 2): allocate and free 512 MiB, then a no-op.
k209_route_calibration <- function(lane, corpus, settings) {
    list(
        list(id = "calibrate", fn = function() {
            held <- numeric(512 * 2^20 / 8)
            held[] <- 1
            held <- NULL
            invisible(gc(full = TRUE))
            invisible(NULL)
        }),
        list(id = "noop", fn = function() invisible(NULL)))
}

# The corpus files by role (generated on the runner by k209-plan.sh's corpus step).
k209_corpus <- function() {
    dir <- k209_args$corpus_dir
    route <- k209_args$route
    if (identical(route, "calibration")) return(list())
    if (route %in% c("lipidomics", "metabolomics")) {
        files <- if (identical(route, "lipidomics")) c(LCMS_Pos = "lipidsearch_lcms_pos.txt", LCMS_Neg = "lipidsearch_lcms_neg.txt") else
            c(LCMS_Pos = "msdial_lcms_pos.txt", GCMS = "msdial_gcms.txt")
        return(list(assays = as.list(stats::setNames(file.path(dir, files), names(files))), design = file.path(dir, "design.tsv"),
            vendor = if (identical(route, "lipidomics")) "lipidsearch" else "msdial"))
    }
    if (identical(route, "transcriptomics")) return(list(counts = file.path(dir, "counts.tsv"), design = file.path(dir, "design.tsv")))
    report <- c(diann = "report.tsv", maxquant = "proteinGroups.txt", pd_tmt = "pd_proteins.tsv")
    format <- k209_args$format
    report_file <- if (identical(k209_args$corpus, "tmt_fixture")) "pd_tmt_peptides.tsv" else report[[format]]
    list(format = format, report = file.path(dir, report_file), fasta = file.path(dir, "proteins.fasta"),
        design = file.path(dir, "design.tsv"))
}

k209_route_steps <- function(lane, corpus, settings) {
    switch(k209_args$route, lipidomics = , metabolomics = , transcriptomics = NULL, DIA = k209_route_dia(lane, corpus, settings),
        DIA_LIMPA = k209_route_limpa(lane, corpus, settings), protein = k209_route_protein(lane, corpus, settings),
        calibration = k209_route_calibration(lane, corpus, settings), stop("unknown route ", k209_args$route))
}

# The owner of a frame's function: the lane, the package, another namespace, or base code.
k209_owner_of <- function(fun) {
    env <- environment(fun)
    if (is.null(env)) return("primitive")
    repeat {
        if (identical(env, k209_h) || identical(env, globalenv())) return("lane")
        if (identical(env, k209_ns)) return("package")
        if (isNamespace(env)) return(paste0("namespace:", getNamespaceName(env)))
        if (identical(env, emptyenv()) || identical(env, baseenv())) return("base")
        env <- parent.env(env)
    }
}

k209_head <- function(lane) {
    project <- file.path(lane$work, "project")
    store <- if (dir.exists(project)) k209_harness$kitFindStore(project) else NULL
    if (is.null(store)) return(list(state = NA_character_, store_bytes = NA_real_))
    head <- k209_harness$kitReadHead(store, k209_api)
    state <- head$current_state
    files <- list.files(store, recursive = TRUE, all.files = TRUE, full.names = TRUE)
    list(state = state, entity_id = head$states[[state]] %||% NA_character_, rows = head$rows[[state]] %||% NA,
        store_bytes = sum(as.numeric(file.size(files)), na.rm = TRUE), history_length = length(head$history), store_dir = store)
}

for (name in c("k209_settings", "k209_quiet", "k209_project", "k209_fake_session", "k209_report", "k209_ui_defaults",
    "k209_omics_journey", "k209_prot_import", "k209_prot_design", "k209_mirror_progress", "k209_publish_peptide",
    "k209_publish_protein", "k209_scoped", "k209_protein_qc", "k209_aesthetics", "k209_prot_tail", "k209_prot_da",
    "k209_route_dia", "k209_route_limpa", "k209_route_protein", "k209_route_calibration", "k209_corpus", "k209_route_steps",
    "k209_head")) {
    fn <- get(name)
    environment(fn) <- k209_h
    assign(name, fn)
}

# --- Prepare (outside the figures) ---------------------------------------------------------------------------------
Sys.setenv(MULTISCHOLAR_TEST_MODE = "true")
k209_lane <- new.env(parent = emptyenv())
k209_lane$work <- k209_args$work
dir.create(k209_lane$work, recursive = TRUE, showWarnings = FALSE)
k209_bound <- k209_settings()
k209_corpus_files <- k209_corpus()
if (identical(k209_args$route, "DIA_LIMPA") && !requireNamespace(k209_bound$engine, quietly = TRUE))
    stop(structure(class = c("workflow_engine_unavailable", "error", "condition"),
        list(message = sprintf("the bound DPC engine %s is not installed", k209_bound$engine), call = NULL)))
k209_steps <- k209_route_steps(k209_lane, k209_corpus_files, k209_bound)
k209_expected <- if (identical(k209_args$route, "calibration")) list() else {
    ruv_mode <- k209_bound$ruv_mode %||% "covariates"
    harness_steps <- k209_harness$kitJourneySteps(k209_args$route, ruv_mode)
    expected <- stats::setNames(lapply(harness_steps, function(step) step$state %||% NA_character_), vapply(harness_steps, `[[`, "", "id"))
    # The harness map's protein route names `ruv_corrected`; a skipped RUV commits `normalized` (kit-021's recorder,
    # "a skipped RUV commits `normalized` where the route names `ruv_corrected`"; diaNormalizationCommit).
    if (k209_args$route %in% c("protein", "DIA") && identical(ruv_mode, "skip") && "normalisation_ruv" %in% names(expected))
        expected$normalisation_ruv <- "normalized"
    expected
}
k209_env <- list(MULTISCHOLAR_TEST_MODE = Sys.getenv("MULTISCHOLAR_TEST_MODE"), RAYON_NUM_THREADS = Sys.getenv("RAYON_NUM_THREADS", NA),
    OMP_NUM_THREADS = Sys.getenv("OMP_NUM_THREADS", NA), OPENBLAS_NUM_THREADS = Sys.getenv("OPENBLAS_NUM_THREADS", NA),
    RSTUDIO_PANDOC = Sys.getenv("RSTUDIO_PANDOC", NA), r = R.version.string, blas = extSoftVersion()[["BLAS"]],
    lapack = La_library(), libpaths = .libPaths())
k209_build_facts <- list(
    dpc_engine_default = tryCatch(resolveLimpaDpcEngine(), error = function(e) paste("error:", conditionMessage(e))),
    fastlimpadpc_installed = requireNamespace("fastlimpadpc", quietly = TRUE),
    ruv_design_panel_option = getOption("MultiScholaR.ruv_design_panel", NA),
    ruv_design_panel_enabled = if (exists("omicGuiRuvDesignEnabled", envir = k209_ns, inherits = FALSE))
        isTRUE(get("omicGuiRuvDesignEnabled", envir = k209_ns)()) else NA,
    prot_manual_ruv_defaults = tryCatch(omicGuiConfig("proteomics")$norm$ruv_iii[c("default", "percentage", "k")],
        error = function(e) paste("error:", conditionMessage(e))),
    managed_runner = tryCatch(is.null(workflowOmicsManagedActions("proteomics")), error = function(e) NA),
    version = as.character(utils::packageVersion("MultiScholaR", lib.loc = k209_args$lib)),
    source_commit = utils::packageDescription("MultiScholaR", lib.loc = k209_args$lib)$`Config/MultiScholaR/SourceCommit` %||% NA)

# Raise-time capture: a tab's observer catches its errors in a tryCatch, which an outer calling handler never sees.
# base::stop, rlang::abort and cli::cli_abort are traced (lane only) so the frames of every raised error are recorded
# at the moment it is raised; the last one before a failed or uncommitted step classifies it.
k209_rec$raised <- 0L
k209_rec$last_raise <- NULL
k209_on_raise <- function() {
    calls <- sys.calls()
    k209_rec$raised <- k209_rec$raised + 1L
    k209_rec$last_raise <- lapply(seq_along(calls), function(j) {
        fun <- tryCatch(sys.function(j), error = function(e) NULL)
        list(call = substr(paste(deparse(calls[[j]][[1L]]), collapse = ""), 1L, 80L),
            owner = if (is.null(fun) || identical(fun, k209_on_raise)) "handler" else k209_owner_of(fun))
    })
}
k209_traced <- vapply(list(c("stop", "base"), c("abort", "rlang"), c("cli_abort", "cli")), function(target) {
    isTRUE(tryCatch({
        suppressMessages(trace(target[[1L]], tracer = quote(k209_on_raise()), print = FALSE, where = asNamespace(target[[2L]])))
        TRUE
    }, error = function(e) FALSE))
}, logical(1))
names(k209_traced) <- c("base::stop", "rlang::abort", "cli::cli_abort")

# --- Measure ----------------------------------------------------------------------------------------------------------
# One step: the entry marks, the resets, the step, the exit marks; then (outside the figures) the head. A step that
# errors, or leaves a state other than the harness map's, ends the journey (E3).
k209_measure <- function(steps, lane) {
    for (i in seq_along(steps)) {
        step <- steps[[i]]
        g_in <- gc(full = TRUE)
        entry <- list(heap_mib = round(sum(g_in[, 2L]), 1), vmrss_mib = k209_status_mib("VmRSS"))
        invisible(gc(reset = TRUE))
        writeLines("5", "/proc/self/clear_refs")
        hwm_reset <- k209_status_mib("VmHWM")
        t <- proc.time()[["elapsed"]]
        # Every error signalled during the step, at signal time (a tab's observer may catch it and only notify): the last
        # one's frames and their owners classify a failed step or a state not committed.
        k209_rec$signal <- NULL
        k209_rec$signals <- 0L
        k209_rec$raised <- 0L
        k209_rec$last_raise <- NULL
        handler <- function(error) {
            if (inherits(error, c("shiny.silent.error", "validation"))) return(invisible(NULL))
            calls <- sys.calls()
            k209_rec$signals <- k209_rec$signals + 1L
            k209_rec$signal <- list(class = class(error)[[1L]], message = substr(conditionMessage(error), 1L, 800L),
                frames = lapply(seq_along(calls), function(j) {
                    fun <- tryCatch(sys.function(j), error = function(e) NULL)
                    list(call = substr(paste(deparse(calls[[j]][[1L]]), collapse = ""), 1L, 80L),
                        owner = if (is.null(fun) || identical(fun, handler)) "handler" else k209_owner_of(fun))
                }))
        }
        failed <- tryCatch({
            withCallingHandlers(step$fn(), error = handler)
            NULL
        }, error = function(error) list(class = class(error)[[1L]], message = substr(conditionMessage(error), 1L, 800L)))
        seconds <- round(proc.time()[["elapsed"]] - t, 2)
        hwm <- k209_status_mib("VmHWM")
        gmax <- round(sum(gc()[, 6L]), 1)
        mark <- list(step = step$id, index = i, entry_heap_mib = entry$heap_mib, entry_vmrss_mib = entry$vmrss_mib,
            vmhwm_after_reset_mib = hwm_reset, vmhwm_mib = hwm, heap_max_used_mib = gmax, exit_vmrss_mib = k209_status_mib("VmRSS"),
            seconds = seconds, cgroup_procs = k209_cgroup_procs(), status = if (is.null(failed)) "completed" else "error")
        state <- tryCatch(k209_head(lane), error = function(e) list(state = NA_character_, head_error = conditionMessage(e)))
        mark$state <- state$state
        mapped <- step$id %in% names(k209_expected)
        mark$expected_state <- if (mapped) k209_expected[[step$id]] else "not in the harness map"
        # A step the harness map has no state for (NA) or does not list (the reset control) is not judged on E3.
        mark$state_as_expected <- if (!mapped) NA else is.na(mark$expected_state) || identical(mark$state, mark$expected_state)
        mark$entity_id <- state$entity_id %||% NA_character_
        mark$rows <- state$rows %||% NA
        mark$store_bytes <- state$store_bytes %||% NA_real_
        mark$head_error <- state$head_error %||% NA_character_
        mark$errors_signalled <- k209_rec$signals
        mark$errors_raised <- k209_rec$raised
        k209_rec$marks[[i]] <- mark
        k209_write(paste0(k209_args$out, ".marks.json"), c(k209_common(), list(marks = k209_rec$marks)))
        cat(sprintf("k209 %-28s %-9s VmHWM %8.1f MiB  heap in %7.1f MiB  %8.1f s  state %s\n", step$id, mark$status, hwm,
            entry$heap_mib, seconds, mark$state %||% NA))
        if (is.null(failed) && isFALSE(mark$state_as_expected))
            failed <- list(class = "k209_state_not_committed", message = sprintf("the step left state %s, the map names %s%s",
                mark$state %||% NA, mark$expected_state, if (is.null(k209_rec$signal)) "" else
                    paste0("; last error signalled: ", k209_rec$signal$message)))
        if (!is.null(failed)) {
            k209_rec$error <- c(failed, list(step = step$id, signalled_class = k209_rec$signal$class %||% NA_character_,
                frames = k209_rec$signal$frames %||% k209_rec$last_raise,
                frames_source = if (!is.null(k209_rec$signal$frames)) "calling handler" else if (!is.null(k209_rec$last_raise))
                    "raise trace (the last error raised in the step)" else "none"))
            break
        }
    }
    invisible(NULL)
}

invisible(gc(full = TRUE))
k209_g0 <- gc(full = TRUE, reset = TRUE)
writeLines("5", "/proc/self/clear_refs")
k209_rec$t0 <- proc.time()[["elapsed"]]
k209_baseline <- list(heap_mib = round(sum(k209_g0[, 2L]), 1), vmrss_mib = k209_status_mib("VmRSS"),
    vmhwm_after_reset_mib = k209_status_mib("VmHWM"), load_s = round(k209_rec$t0 - k209_started, 2),
    loaded_from = k209_loaded_from, cgroup_procs = k209_cgroup_procs())
k209_write(paste0(k209_args$out, ".baseline.json"), c(k209_common(), list(status = "baseline", baseline = k209_baseline,
    settings = k209_bound, environment = k209_env, build_facts = k209_build_facts, traced = as.list(k209_traced),
    steps = if (is.null(k209_steps)) names(k209_expected) else vapply(k209_steps, `[[`, "", "id"))))

# The omics routes run their steps inside the module server's testServer; the others here. A failure outside a step
# (the module server's start) is recorded as the journey's error with no step.
k209_outer <- tryCatch({
    if (is.null(k209_steps)) k209_omics_journey(k209_lane, k209_corpus_files, k209_bound) else k209_measure(k209_steps, k209_lane)
    NULL
}, error = function(error) list(class = class(error)[[1L]], message = substr(conditionMessage(error), 1L, 800L), step = NA_character_))
k209_error <- k209_rec$error %||% k209_outer
k209_gf <- gc(full = TRUE)
k209_final <- list(heap_mib = round(sum(k209_gf[, 2L]), 1), vmrss_mib = k209_status_mib("VmRSS"),
    wall_s = round(proc.time()[["elapsed"]] - k209_rec$t0, 2))
k209_status <- if (is.null(k209_error)) "completed" else "error"

# --- After the figures: fingerprints (G0a only), the lane's holdings --------------------------------------------
k209_fingerprints <- NULL
if (identical(k209_status, "completed") && isTRUE(k209_args$fingerprint)) {
    k209_fingerprints <- tryCatch({
        store <- k209_harness$kitFindStore(file.path(k209_lane$work, "project"))
        prints <- k209_harness$kitFingerprintStore(store, k209_api)
        list(head_sha256 = prints$head_sha256,
            states = lapply(prints$states, function(state) digest::digest(state, algo = "sha256")))
    }, error = function(e) list(error = conditionMessage(e)))
}
k209_holdings <- list(lane = sort(ls(k209_lane)), workflow = if (is.environment(k209_lane$workflow)) sort(ls(k209_lane$workflow)) else NULL)
k209_record <- c(k209_common(), list(status = k209_status, error = k209_error, baseline = k209_baseline, final = k209_final,
    marks = k209_rec$marks, settings = k209_bound, environment = k209_env, build_facts = k209_build_facts,
    report = k209_rec$report, fingerprints = k209_fingerprints, holdings = k209_holdings))
k209_write(paste0(k209_args$out, ".json"), k209_record)
cat(sprintf("k209 journey %s: %s; %.1f s\n", k209_args$id, k209_status, k209_final$wall_s))
if (!is.null(k209_error)) cat("k209 error at", k209_error$step, ":", k209_error$message, "\n")
quit(save = "no", status = if (identical(k209_status, "completed")) 0L else 3L)
