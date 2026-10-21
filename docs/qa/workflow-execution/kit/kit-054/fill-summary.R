# KIT-054 (user request 2026-10-21): fill a lane's results_summary as the Summary tab's "Copy publication outputs" and
# "Generate Report" actions do, on a finished member's project. On the owner path (a store under scripts/proteomics)
# the analysis record and study_parameters.txt are written from the store and the workflow-type template is rendered;
# on the in-memory lane (no store) only the publication copy runs, since the session state is gone with the application.
# Data handling rule: the script prints step names and counts only; its log is not read beyond those markers.
# Usage: R_LIBS=<lane's library> Rscript docs/qa/workflow-execution/kit/kit-054/fill-summary.R <stage dir> <LFQ|TMT> <format>
suppressMessages(library(MultiScholaR))
ns <- asNamespace("MultiScholaR"); api <- function(name) get(name, ns)
`%||%` <- function(a, b) if (is.null(a)) b else a
for (f in ls(ns, all.names = FALSE)) if (!exists(f, envir = globalenv(), inherits = TRUE) && is.function(get(f, ns))) assign(f, get(f, ns), envir = globalenv())
args <- commandArgs(TRUE); stopifnot(length(args) == 3L)
stage <- args[[1L]]; workflow_type <- args[[2L]]; format <- args[[3L]]
label <- "dia_baseline_smoke"
base <- file.path(stage, "candidate-1", "project")
Sys.setenv(MULTISCHOLAR_TEST_MODE = "true")
# The GUI's project layout (<project>/<label>/{scripts,results,results_summary}/proteomics), assembled by hand as the
# KIT-038 record check does; setupDirectories would build the "<omic>_<label>" layout the GUI does not use.
project_base <- file.path(base, label)
results <- file.path(project_base, "results", "proteomics")
summary_base <- file.path(project_base, "results_summary", "proteomics")
timed <- list.dirs(file.path(results, "publication_graphs", "filtering_qc"), recursive = FALSE)
paths <- list(base_dir = project_base, results_dir = results, results_summary_dir = summary_base,
    study_report_dir = file.path(summary_base, "Study_report"), qc_figures_dir = file.path(summary_base, "QC_figures"),
    publication_figures_dir = file.path(summary_base, "Publication_figures"), publication_tables_dir = file.path(summary_base, "Publication_tables"),
    source_dir = file.path(project_base, "scripts", "proteomics"),
    protein_qc_dir = file.path(results, "protein_qc"), qc_dir = file.path(results, "protein_qc"),
    peptide_qc_dir = file.path(results, "peptide_qc"), feature_qc_dir = file.path(results, "protein_qc"), # the GUI points feature_qc_dir at the protein QC on a protein-level route
    da_output_dir = file.path(results, "da_proteins"), publication_graphs_dir = file.path(results, "publication_graphs"),
    pathway_dir = file.path(results, "pathway_enrichment"),
    time_dir = if (length(timed)) timed[[length(timed)]] else file.path(results, "publication_graphs"),
    integration_dir = file.path(project_base, "integration"))
for (d in c("Study_report", "Publication_figures", "Publication_tables", "QC_figures")) dir.create(file.path(summary_base, d), recursive = TRUE, showWarnings = FALSE)
project_dirs <- list(proteomics = paths)
assign("project_dirs", project_dirs, envir = globalenv())
source_dir <- paths$source_dir
store_root <- list.files(source_dir, pattern = "^[.]dia-native-", all.files = TRUE, full.names = TRUE, include.dirs = TRUE)
design <- utils::read.delim(file.path(source_dir, "design_matrix.tab"), stringsAsFactors = FALSE, check.names = FALSE)
cat("project:", basename(stage), "| type", workflow_type, "| design rows", nrow(design), "| store", length(store_root), "\n")
sink_log <- function(expr) invisible(capture.output(suppressMessages(force(expr))))
if (length(store_root) == 1L) {
    identity <- api("readSavedStoreIdentity")(file.path(store_root, ".workflow-entities", "store.json"))
    manager <- api("newDIAStateAdapter")(api("newWorkflowEntityStore")(store_root, identity$project_id, identity$workflow_id), reopen = TRUE)
    on.exit(api("closeDIAStateAdapter")(manager), add = TRUE)
    config <- readProtGuiConfig(file.path(source_dir, "config.ini"))
    api("withDIAStateScope")(manager, function() {
        manifest <- api("resolveDIAStateValue")(manager$.context, manager$.context$snapshot$da_results)
        contrasts <- as.data.frame(manifest$contrasts, stringsAsFactors = FALSE)
        cat("contrasts:", nrow(contrasts), "| head state:", manager$current_state, "\n")
        json <- file.path(source_dir, "analysis_record.json")
        sink_log(writeWorkflowAnalysisRecord(manager, json))
        cat("record written:", file.exists(json), "\n")
        sink_log(createWorkflowArgsFromConfig(format, source_dir_path = source_dir, final_s4_object = getDIAStateReference(manager),
            contrasts_tbl = contrasts, workflow_data = list(state_manager = manager, config_list = config)))
        cat("study_parameters.txt written:", file.exists(file.path(source_dir, "study_parameters.txt")), "\n")
        sink_log(copyToResultsSummary(omic_type = "proteomics", experiment_label = label, contrasts_tbl = contrasts, design_matrix = design, force = TRUE))
        cat("publication copy done\n")
        template <- paste0(workflow_type, "_report.rmd")
        template_source <- system.file("reports", "proteomics", template, package = "MultiScholaR")
        stopifnot(nzchar(template_source))
        file.copy(template_source, file.path(source_dir, template), overwrite = TRUE)
        t0 <- Sys.time()
        rendered <- tryCatch(RenderReport(omic_type = "proteomics", experiment_label = label, rmd_filename = template), error = function(e) e)
        if (inherits(rendered, "error")) cat("report render FAILED:", conditionMessage(rendered), "\n") else
            cat("report rendered:", if (is.character(rendered)) basename(rendered) else paste("<", class(rendered)[[1L]], ">"),
                "in", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s\n")
        study_report <- list.files(paths$study_report_dir, recursive = TRUE)
        cat("Study_report files:", length(study_report), "\n")
    })
} else {
    contrasts_file <- file.path(source_dir, "contrast_strings.tab")
    contrasts <- if (file.exists(contrasts_file)) utils::read.delim(contrasts_file, stringsAsFactors = FALSE, check.names = FALSE) else NULL
    sink_log(copyToResultsSummary(omic_type = "proteomics", experiment_label = label, contrasts_tbl = contrasts, design_matrix = design, force = TRUE))
    cat("publication copy done (no store: no record, no report)\n")
}
summary_dir <- file.path(base, label, "results_summary", "proteomics")
for (d in list.dirs(summary_dir, recursive = FALSE)) cat(sprintf("  %-22s %d files\n", basename(d), length(list.files(d, recursive = TRUE))))
