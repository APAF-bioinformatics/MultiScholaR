# KIT-043 render check: the LFQ (MaxQuant) and TMT (Proteome Discoverer) e2e fixtures run on the on-disk owner
# through DA with two non-default values (protein intensity cutoff percentile 2, DA q-value 0.1), the record and
# study_parameters.txt are written by the writer, and each template is rendered three ways:
#   with the record (this ticket's template), without the record (this ticket's template), and without the record
#   with the pre-ticket template (the baseline: the previous commit's template; the control build's template is
#   also rendered, and its difference is the KIT-039 parser fix). Texts are compared with dates masked.
# Usage: R_LIBS=<library> Rscript docs/qa/workflow-execution/kit/kit-043/report-check.R <baseline template dir> <control template dir> <label>
# With KIT043_INSTALLED=1 the installed package on R_LIBS is used instead of load_all, so the record's software block
# carries the install's source commit (the source tree's DESCRIPTION has none, and the record then says "not recorded").
if (identical(Sys.getenv("KIT043_INSTALLED"), "1")) suppressMessages(library(MultiScholaR)) else suppressMessages(devtools::load_all(quiet = TRUE))
args <- commandArgs(TRUE); baseline_dir <- normalizePath(args[[1L]]); control_dir <- normalizePath(args[[2L]]); label <- args[[3L]]
started <- Sys.time()
repo <- normalizePath(".")
ns <- asNamespace("MultiScholaR"); for (f in c("initializeDIANativeImport", "withDIAStateScope", "registerDIAInputFiles", "applyProtImportResultToWorkflow",
    "getDIAImportMetadata", "completeDIAInputFiles", "buildProtDesignStateCheckpoint", "workflowOperationEntry", "getDIAStateReference", "setDIAStepProvenance",
    "completeDIAProteinAction", "runSavedOperation", "diaNormalizationOperation", "readDIAStateMetadata", "commitDIAReferenceState", "diaPearsonPairsOperation",
    "diaPearsonSelectionOperation", "runDIAReferenceDaAnalysis", "writeWorkflowAnalysisRecord", "workflowRecordValueText", "decodeWorkflowRecordValue",
    "closeDIAStateAdapter")) if (!exists(f, inherits = TRUE)) assign(f, get(f, ns))
fixtures <- file.path(repo, "tests/testdata/e2e")
evidence <- file.path(repo, "docs/qa/workflow-execution/kit/kit-043")
work <- file.path(Sys.getenv("HOME"), ".cache", paste0("multischolar-kit043-report-", label))
unlink(work, recursive = TRUE); dir.create(work, recursive = TRUE)
on.exit(unlink(work, recursive = TRUE), add = TRUE)
cases <- list(
    maxquant = list(file = file.path(fixtures, "prot_lfq", "proteinGroups.txt"), fasta = file.path(fixtures, "prot_lfq", "proteins.fasta"),
        design = file.path(fixtures, "prot_lfq", "design.tsv"), read = function(p) importMaxQuantData(p), type = "LFQ", template = "LFQ_report.rmd"),
    pd_tmt = list(file = file.path(fixtures, "prot_tmt", "pd_tmt_peptides.tsv"), fasta = file.path(fixtures, "prot_tmt", "proteins.fasta"),
        design = file.path(fixtures, "prot_tmt", "design.tsv"), read = function(p) importProteomeDiscovererTMTData(p), type = "TMT", template = "TMT_report.rmd"))
fake_session <- function() {
    session <- new.env(parent = emptyenv()); session$userData <- new.env(parent = emptyenv())
    session$onFlush <- session$onFlushed <- session$onSessionEnded <- function(...) function() invisible(NULL)
    session
}
config <- readConfigFile(system.file("config", "config.ini", package = "MultiScholaR"), create_workers = FALSE)
docx_text <- function(path, dir) {
    xml <- utils::unzip(path, files = "word/document.xml", exdir = dir)
    text <- gsub("<[^>]+>", " ", paste(readLines(xml, warn = FALSE), collapse = ""))
    gsub("[[:space:]]+", " ", text)
}
mask <- function(text) gsub("[0-9]{4}-[0-9]{2}-[0-9]{2}( [0-9]{2}:[0-9]{2}:[0-9]{2})?", "<time>", text)
first_difference <- function(a, b) {
    x <- strsplit(a, "")[[1L]]; y <- strsplit(b, "")[[1L]]
    n <- min(length(x), length(y)); first <- which(x[seq_len(n)] != y[seq_len(n)])
    at <- if (length(first)) first[[1L]] else n + 1L
    list(position = at, a = substr(a, max(1L, at - 120L), at + 160L), b = substr(b, max(1L, at - 120L), at + 160L))
}
intensity_percentile <- 2L; da_qvalue <- 0.1
run_case <- function(name, case) {
    out <- list(format = name, type = case$type, template = case$template)
    base <- file.path(work, name); dir.create(base, recursive = TRUE)
    project_dirs <- setupDirectories(base_dir = base, omic_types = "proteomics", label = name, force = TRUE, reuse_existing = TRUE)
    assign("project_dirs", project_dirs, envir = globalenv())
    paths <- project_dirs[[1L]]
    paths$time_dir <- paths$time_dir %||% file.path(paths$publication_graphs_dir, "filtering_qc")
    dir.create(paths$time_dir, recursive = TRUE, showWarnings = FALSE)
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- WorkflowState$new(audit_enabled = FALSE)
    workflow$config_list <- config
    workflow$config_list$globalParameters$workflow_type <- case$type
    workflow$qc_params <- list(protein_qc = list())
    initializeDIANativeImport(workflow, paths, fake_session(), name, "iq")
    manager <- workflow$state_manager
    on.exit(closeDIAStateAdapter(manager), add = TRUE)
    design_in <- utils::read.delim(case$design, stringsAsFactors = FALSE)
    contrasts <- data.frame(contrasts = "groupKO-groupWT", friendly_names = "KO_vs_WT", full_format = "KO_vs_WT=groupKO-groupWT", stringsAsFactors = FALSE)
    withDIAStateScope(manager, function() {
        inputs <- c(case$file, case$fasta); names(inputs) <- c(paste0(name, "_report"), "fasta")
        registerDIAInputFiles(manager, inputs)
        imported <- case$read(case$file)
        applyProtImportResultToWorkflow(workflow, imported, name, case$fasta, logInfo = function(...) NULL, showNotification = function(...) NULL)
        metadata <- getDIAImportMetadata(manager)
        completeDIAInputFiles(manager)
        runs <- as.character(metadata$samples$Run)
        group <- vapply(runs, function(run) design_in$group[[which(endsWith(run, design_in$sample))[[1L]]]], "")
        design <- data.frame(Run = runs, group = unname(group), stringsAsFactors = FALSE)
        design$replicates <- stats::ave(seq_len(nrow(design)), design$group, FUN = seq_along)
        workflow$design_matrix <- design
        workflow$contrasts_tbl <- contrasts
        buildProtDesignStateCheckpoint(workflow, case$type, "Import")
        # The template reads the design matrix beside study_parameters.txt (the summary page writes it in the GUI).
        utils::write.table(design, file.path(paths$source_dir, "design_matrix.tab"), sep = "\t", quote = FALSE, row.names = FALSE)
    })
    steps <- list(
        list(operation = "protein_intensity_filter", parameters = list(strict_mode = FALSE, min_reps_per_group = 2L, min_groups = 2L,
            proteins_intensity_cutoff_percentile = intensity_percentile)),
        list(operation = "protein_duplicate_removal", parameters = list(aggregation_method = "mean")),
        list(operation = "protein_replicate_filter", parameters = list(grouping_variable = "group", parallel_cores = 1L)))
    for (step in steps) {
        entry <- workflowOperationEntry(step$operation)
        variant <- entry$default_variant
        withDIAStateScope(manager, function() {
            directory <- file.path(base, "requests", step$operation); dir.create(directory, recursive = TRUE, showWarnings = FALSE)
            request <- list(operation = step$operation, input_references = list(state = getDIAStateReference(manager)),
                parameters = step$parameters, worker_cap = 12 * 1024^3, directory = directory)
            value <- get(entry$handler, mode = "function")(manager, request, variant)
            value$operation <- step$operation
            setDIAStepProvenance(manager, step$operation, variant, entry$parameters, entry$declaration$fixed, values = step$parameters)
            completeDIAProteinAction(workflow, value)
        })
    }
    norm_parameters <- list(norm_method = "cyclicloess", ruv_mode = "skip", auto_percentage_min = 1, auto_percentage_max = 20,
        separation_metric = "max_difference", k_penalty_weight = 0.5, adaptive_k_penalty = TRUE, ruv_percentage = 5, ruv_k = NA, grouping = "group")
    withDIAStateScope(manager, function() {
        run <- runSavedOperation(manager, diaNormalizationOperation(), getDIAStateReference(manager), norm_parameters)
        saved <- run$value
        final <- readDIAStateMetadata(manager, saved$final)$data
        final$products$normalization_science <- run$reference
        commitDIAReferenceState(manager, "normalized", final, list(norm_method = "cyclicloess", ruv_mode = "skip", ruv_k = saved$best_k,
            percentage_as_neg_ctrl = saved$percentage, ruv_applied = FALSE), "Post-normalization complete: RUV-III skipped by user")
        pairs <- runSavedOperation(manager, diaPearsonPairsOperation(), saved$final, list(group = "group", tech_rep_remove_regex = "pool"))
        selection <- runSavedOperation(manager, diaPearsonSelectionOperation(), saved$final,
            list(pairs = pairs$reference$entity_id, threshold = 0.5, correlations = pairs$value), dependencies = list(pairs$reference))
        spec <- readDIAStateMetadata(manager, selection$value)$data
        commitDIAReferenceState(manager, "correlation_filtered", spec, list(min_pearson_correlation_threshold = 0.5), "Applied final sample correlation filter")
    })
    da_data <- new.env(parent = emptyenv()); da_data$current_s4_object <- getDIAStateReference(manager); da_data$contrasts_available <- contrasts$contrasts
    da_workflow <- new.env(parent = emptyenv()); da_workflow$state_manager <- manager; da_workflow$contrasts_tbl <- contrasts
    runDIAReferenceDaAnalysis(list(formula_string = "~ 0 + group", da_q_val_thresh = da_qvalue, treat_lfc_cutoff = 0), da_data, da_workflow, paths,
        contrasts = contrasts, publish = FALSE)
    # The record and study_parameters.txt, as the summary page writes them (the writer renders from the record on an owner).
    json <- file.path(paths$source_dir, "analysis_record.json")
    writeWorkflowAnalysisRecord(manager, json)
    record <- jsonlite::read_json(json, simplifyVector = FALSE)
    invisible(capture.output(createWorkflowArgsFromConfig(name, source_dir_path = paths$source_dir, final_s4_object = getDIAStateReference(manager),
        contrasts_tbl = contrasts, workflow_data = list(state_manager = manager, config_list = workflow$config_list))))
    study <- readLines(file.path(paths$source_dir, "study_parameters.txt"))
    values_present <- character()
    for (step in Filter(function(s) !identical(s$operation, "legacy") && length(s$parameters), record$steps)) for (parameter in step$parameters) {
        line <- sprintf("  %s [%s] = %s", parameter$name, parameter$provenance, workflowRecordValueText(parameter$value))
        if (!line %in% study) values_present <- c(values_present, paste(step$state, parameter$name))
    }
    da_parameters <- record$analyses$differential_abundance$parameters
    out$study <- list(lines = length(study), truncation_marks = sum(grepl("omitted|\\[list object\\]|\\.\\.\\.", study)),
        record_section = "Analysis Record:" %in% study, values_missing = values_present, record_beside = file.exists(json),
        da_qvalue_recorded = any(vapply(da_parameters, function(p) identical(p$name, "qvalue") && identical(decodeWorkflowRecordValue(p$value), da_qvalue), logical(1))),
        intensity_recorded = any(vapply(Filter(function(s) identical(s$operation, "protein_intensity_filter"), record$steps)[[1L]]$parameters,
            function(p) identical(p$name, "proteins_intensity_cutoff_percentile") && identical(decodeWorkflowRecordValue(p$value), intensity_percentile), logical(1))))
    out$record <- list(workflow = record$workflow, software = record$software[c("version", "sha", "build")],
        operations = vapply(record$steps, function(s) s$operation, ""), da = !is.null(record$analyses$differential_abundance))
    render <- function(template_dir, tag) {
        target <- file.path(base, "scripts", "proteomics", case$template)
        dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
        file.copy(file.path(template_dir, case$template), target, overwrite = TRUE)
        t0 <- Sys.time()
        path <- tryCatch(RenderReport(omic_type = "proteomics", experiment_label = name, rmd_filename = case$template), error = function(e) e)
        if (inherits(path, "error")) return(list(error = conditionMessage(path)))
        keep <- file.path(base, paste0(tag, ".docx")); file.copy(path, keep, overwrite = TRUE)
        list(path = keep, seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")), text = docx_text(keep, file.path(base, paste0(tag, "-xml"))))
    }
    current_dir <- file.path(repo, "inst/reports/proteomics")
    with_record <- render(current_dir, "with-record")
    file.rename(json, file.path(base, "analysis_record.json.removed"))
    without_record <- render(current_dir, "without-record")
    baseline <- render(baseline_dir, "baseline")
    control <- render(control_dir, "control")
    file.rename(file.path(base, "analysis_record.json.removed"), json)
    for (tag in c("with_record", "without_record", "baseline", "control")) {
        r <- get(tag); out$renders[[tag]] <- r[c("error", "seconds")]
        if (!is.null(r$text)) writeLines(r$text, file.path(evidence, sprintf("%s-%s-text.txt", sub("[.]rmd$", "", case$template), gsub("_", "-", tag))))
    }
    if (is.null(with_record$error) && is.null(without_record$error) && is.null(baseline$error)) {
        section <- function(text) regmatches(text, regexpr("Version Information.{0,400}", text, perl = TRUE))
        out$version_with_record <- section(with_record$text)
        out$version_without_record <- section(without_record$text)
        has <- function(pattern, text) is.character(pattern) && length(pattern) == 1L && !is.na(pattern) && nzchar(pattern) &&
            length(text) == 1L && grepl(pattern, text, fixed = TRUE)
        out$with_record_names_software <- has(record$software$version, out$version_with_record) &&
            has(substr(record$software$sha %||% "", 1, 12), out$version_with_record) && length(out$version_with_record) == 1L &&
            !grepl("\\bNA\\b", out$version_with_record)
        # The recorded non-default values in the rendered text: the percentile in the methods sentence, the q-value in the DA text.
        out$with_record_prints_percentile <- grepl(sprintf("below the %dst percentile", intensity_percentile), with_record$text, fixed = TRUE)
        # Without the record the parser reads study_parameters.txt, itself rendered from the record, so the recorded
        # percentile still shows; the template default (1) would show only without both files.
        out$without_record_prints_default_percentile <- grepl("below the 1st percentile", without_record$text, fixed = TRUE)
        out$without_record_prints_recorded_percentile <- grepl(sprintf("below the %dst percentile", intensity_percentile), without_record$text, fixed = TRUE)
        out$with_record_prints_qvalue <- grepl(sprintf("q-values less than %s", da_qvalue), with_record$text, fixed = TRUE)
        out$without_record_prints_default_qvalue <- grepl("q-values less than 0.05", without_record$text, fixed = TRUE)
        out$defaults_note <- regmatches(with_record$text, regexpr("Parameter not recorded; default shown:[^.]*[.]", with_record$text))
        out$without_record_equals_baseline <- identical(mask(baseline$text), mask(without_record$text))
        if (!out$without_record_equals_baseline) out$first_difference_baseline <- first_difference(mask(baseline$text), mask(without_record$text))
        if (is.null(control$error)) {
            out$without_record_equals_control <- identical(mask(control$text), mask(without_record$text))
            if (!out$without_record_equals_control) out$first_difference_control <- first_difference(mask(control$text), mask(without_record$text))
        } else out$control_error <- control$error
        file.copy(with_record$path, file.path(evidence, sprintf("%s-with-record.docx", sub("[.]rmd$", "", case$template))), overwrite = TRUE)
    }
    out
}
results <- list(schema = "kit.report_check.v1", ticket = "KIT-043", label = label, date = format(Sys.Date()), pandoc = as.character(rmarkdown::pandoc_version()),
    non_defaults = list(proteins_intensity_cutoff_percentile = intensity_percentile, da_qvalue = da_qvalue))
for (name in names(cases)) {
    results[[name]] <- tryCatch(run_case(name, cases[[name]]), error = function(e) list(format = name, error = conditionMessage(e)))
    r <- results[[name]]
    if (!is.null(r$error)) { cat(name, "ERROR:", r$error, "\n"); next }
    cat(name, "| record", r$record$workflow, "ops", length(r$record$operations), "DA", r$record$da,
        "| study lines", r$study$lines, "truncation", r$study$truncation_marks, "missing values", length(r$study$values_missing), "record beside", r$study$record_beside,
        "nondefaults recorded", r$study$intensity_recorded, r$study$da_qvalue_recorded,
        "| renders", paste(vapply(names(r$renders), function(t) paste0(t, ":", r$renders[[t]]$error %||% "ok"), ""), collapse = " "),
        "| software", r$with_record_names_software %||% NA, "| percentile", r$with_record_prints_percentile %||% NA, "/", r$without_record_prints_default_percentile %||% NA,
        "| qvalue", r$with_record_prints_qvalue %||% NA, "/", r$without_record_prints_default_qvalue %||% NA,
        "| without = baseline", r$without_record_equals_baseline %||% NA, "| without = control", r$without_record_equals_control %||% NA, "\n")
}
# The DIA-NN templates are untouched by this ticket (byte-identical to the baseline commit's).
results$diann_untouched <- all(vapply(c("DIANN_report.rmd", "DIANN_limpa_report.rmd"), function(f)
    identical(digest::digest(file = file.path(repo, "inst/reports/proteomics", f), algo = "sha256"), digest::digest(file = file.path(baseline_dir, f), algo = "sha256")), logical(1)))
results$library <- list(path = find.package("MultiScholaR"), installed = identical(Sys.getenv("KIT043_INSTALLED"), "1"),
    source_commit = utils::packageDescription("MultiScholaR")$`Config/MultiScholaR/SourceCommit`)
results$head_commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
results$tree_clean <- length(system2("git", c("status", "--porcelain", "--", "R", "tests", "inst", "DESCRIPTION", "NAMESPACE"), stdout = TRUE)) == 0L
results$seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
jsonlite::write_json(results, file.path(evidence, paste0("report-check-", label, ".json")), auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = I(17))
cat("DIANN untouched", results$diann_untouched, "| done", round(results$seconds), "s\n")
