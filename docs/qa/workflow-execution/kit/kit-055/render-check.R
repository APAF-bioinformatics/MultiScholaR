# KIT-055 render check: each proteomics report template rendered on a copy of a finished kit project (the KIT-021
# fixture lifecycle members for LFQ and TMT, the KIT-036 screen lifecycle candidate for DIA-NN), the record and
# study_parameters.txt written from the copy's store, the publication outputs copied as the Summary tab copies them,
# the template taken from the working tree, the rendered text extracted from the .docx and asserted against the
# ticket's criteria (no canned sentence, no placeholder, no printed warning, the significance figure resolved, the
# enrichment section conditional, the import format named, the software block last). The fixtures are not sensitive;
# the masked texts are kept as evidence.
# Usage: R_LIBS=<library> Rscript docs/qa/workflow-execution/kit/kit-055/render-check.R <label>
suppressMessages(library(MultiScholaR))
ns <- asNamespace("MultiScholaR"); api <- function(name) get(name, ns)
`%||%` <- function(a, b) if (is.null(a)) b else a
for (f in ls(ns, all.names = FALSE)) if (!exists(f, envir = globalenv(), inherits = TRUE) && is.function(get(f, ns))) assign(f, get(f, ns), envir = globalenv())
args <- commandArgs(TRUE); label <- if (length(args)) args[[1L]] else "kit-055"
repo <- normalizePath(".")
proof <- jsonlite::read_json(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))$harness$root
evidence <- file.path(repo, "docs/qa/workflow-execution/kit/kit-055"); dir.create(evidence, showWarnings = FALSE)
work <- file.path(Sys.getenv("HOME"), ".cache", paste0("multischolar-", label)); unlink(work, recursive = TRUE); dir.create(work, recursive = TRUE)
template_dir <- file.path(repo, "inst/reports/proteomics")
Sys.setenv(MULTISCHOLAR_TEST_MODE = "true")
cases <- list(
    LFQ = list(stage = file.path(proof, "kit-journey-021.fixture-lifecycle"), type = "LFQ", format = "maxquant", template = "LFQ_report.rmd"),
    TMT = list(stage = file.path(proof, "kit-journey-021t.fixture-lifecycle"), type = "TMT", format = "pd_tmt", template = "TMT_report.rmd"),
    DIANN = list(stage = file.path(proof, "kit-journey-036-screen-lifecycle"), type = "DIANN", format = "diann", template = "DIANN_report.rmd"))
docx_text <- function(path, dir) {
    xml <- utils::unzip(path, files = "word/document.xml", exdir = dir)
    text <- gsub("<[^>]+>", " ", paste(readLines(xml, warn = FALSE), collapse = ""))
    gsub("[[:space:]]+", " ", text)
}
mask <- function(text) gsub("[0-9]{4}-[0-9]{2}-[0-9]{2}( [0-9]{2}:[0-9]{2}:[0-9]{2})?", "<time>", text)
canned <- c("All samples passed initial QC criteria", "meets industry standards", "batches merged effectively", "Proceed with downstream analysis",
    "demonstrates robust technical quality", "Batch effects were substantially mitigated", "showed notable batch effects")
placeholders <- c("SOP xxx", "Your_fancy_self", "Salesforce", "SOP number(s) and name(s)")
quiet <- function(expr) invisible(capture.output(suppressMessages(suppressWarnings(force(expr)))))
run_case <- function(name, case) {
    out <- list(template = case$template, stage = case$stage)
    source_stage <- file.path(case$stage, "candidate-1", "project")
    if (!dir.exists(source_stage)) return(c(out, list(error = "stage project missing")))
    base <- file.path(work, name); dir.create(base, recursive = TRUE)
    stopifnot(system2("cp", c("-a", "--reflink=auto", file.path(source_stage, "."), base)) == 0L)
    project_base <- file.path(base, "dia_baseline_smoke")
    results <- file.path(project_base, "results", "proteomics"); summary_base <- file.path(project_base, "results_summary", "proteomics")
    timed <- list.dirs(file.path(results, "publication_graphs", "filtering_qc"), recursive = FALSE)
    paths <- list(base_dir = project_base, results_dir = results, results_summary_dir = summary_base,
        study_report_dir = file.path(summary_base, "Study_report"), qc_figures_dir = file.path(summary_base, "QC_figures"),
        publication_figures_dir = file.path(summary_base, "Publication_figures"), publication_tables_dir = file.path(summary_base, "Publication_tables"),
        source_dir = file.path(project_base, "scripts", "proteomics"), protein_qc_dir = file.path(results, "protein_qc"), qc_dir = file.path(results, "protein_qc"),
        peptide_qc_dir = file.path(results, "peptide_qc"), feature_qc_dir = file.path(results, "protein_qc"),
        da_output_dir = file.path(results, "da_proteins"), publication_graphs_dir = file.path(results, "publication_graphs"),
        pathway_dir = file.path(results, "pathway_enrichment"), time_dir = if (length(timed)) timed[[length(timed)]] else file.path(results, "publication_graphs"),
        integration_dir = file.path(project_base, "integration"))
    for (d in c("Study_report", "Publication_figures", "Publication_tables", "QC_figures")) dir.create(file.path(summary_base, d), recursive = TRUE, showWarnings = FALSE)
    assign("project_dirs", list(proteomics = paths), envir = globalenv())
    source_dir <- paths$source_dir
    store_root <- list.files(source_dir, pattern = "^[.]dia-native-", all.files = TRUE, full.names = TRUE, include.dirs = TRUE)
    stopifnot(length(store_root) == 1L)
    identity <- api("readSavedStoreIdentity")(file.path(store_root, ".workflow-entities", "store.json"))
    manager <- api("newDIAStateAdapter")(api("newWorkflowEntityStore")(store_root, identity$project_id, identity$workflow_id), reopen = TRUE)
    on.exit(api("closeDIAStateAdapter")(manager), add = TRUE)
    config <- readProtGuiConfig(file.path(source_dir, "config.ini"))
    design <- utils::read.delim(file.path(source_dir, "design_matrix.tab"), stringsAsFactors = FALSE, check.names = FALSE)
    render <- function(tag) {
        file.copy(file.path(template_dir, case$template), file.path(source_dir, case$template), overwrite = TRUE)
        rendered <- tryCatch(RenderReport(omic_type = "proteomics", experiment_label = "dia_baseline_smoke", rmd_filename = case$template), error = function(e) e)
        if (inherits(rendered, "error") || !is.character(rendered) || !file.exists(rendered)) return(list(error = if (inherits(rendered, "error")) conditionMessage(rendered) else "no output"))
        keep <- file.path(base, paste0(tag, ".docx")); file.copy(rendered, keep, overwrite = TRUE)
        text <- mask(docx_text(keep, file.path(base, paste0(tag, "-xml"))))
        writeLines(text, file.path(evidence, sprintf("%s-%s-text.txt", sub("[.]rmd$", "", case$template), tag)))
        media <- length(grep("^word/media/", utils::unzip(keep, list = TRUE)$Name))
        list(text = text, media = media)
    }
    api("withDIAStateScope")(manager, function() {
        manifest <- api("resolveDIAStateValue")(manager$.context, manager$.context$snapshot$da_results)
        contrasts <- as.data.frame(manifest$contrasts, stringsAsFactors = FALSE)
        quiet(writeWorkflowAnalysisRecord(manager, file.path(source_dir, "analysis_record.json")))
        quiet(createWorkflowArgsFromConfig(case$format, source_dir_path = source_dir, final_s4_object = getDIAStateReference(manager),
            contrasts_tbl = contrasts, workflow_data = list(state_manager = manager, config_list = config)))
        quiet(copyToResultsSummary(omic_type = "proteomics", experiment_label = "dia_baseline_smoke", contrasts_tbl = contrasts, design_matrix = design, force = TRUE))
        without <- render("without-enrichment")
        # The enrichment section with an enrichment output present (a placeholder image in the copied project only).
        dir.create(file.path(paths$publication_figures_dir, "Enrichment_Plots"), showWarnings = FALSE)
        png(file.path(paths$publication_figures_dir, "Enrichment_Plots", "placeholder_enrichment.png"), width = 200, height = 100); par(mar = rep(0, 4)); plot.new(); dev.off()
        with <- render("with-enrichment")
        out$renders <<- list(without_enrichment = without[c("error", "media")], with_enrichment = with[c("error", "media")])
        if (is.null(without$error)) {
            t <- without$text
            has <- function(x) vapply(x, function(p) grepl(p, t, fixed = TRUE), logical(1))
            sections <- regmatches(t, gregexpr("(Sample Details|Method Details|Results|Differential Abundance Analysis|Analyst interpretation|GO Enrichment Analysis|Software and build)", t))[[1L]]
            out$checks <<- list(
                canned_sentences = names(which(has(canned))), placeholders = names(which(has(placeholders))),
                printed_warnings = grepl("## Warning|Warning:", t, fixed = FALSE),
                numsig_de_name = grepl("NumSigDeMolecules", t, fixed = TRUE),
                figure2_caption = grepl("Figure 2: Number of significant", t, fixed = TRUE),
                media_without_enrichment = without$media,
                enrichment_not_run_sentence = grepl("No enrichment analysis was run for this project", t, fixed = TRUE),
                enrichment_methods_absent = !grepl("clusterProfiler R package", t, fixed = TRUE),
                enrichment_methods_present_with_plots = !is.null(with$text) && grepl("clusterProfiler R package|g:Profiler", with$text),
                not_run_sentence_absent_with_plots = !is.null(with$text) && !grepl("No enrichment analysis was run", with$text, fixed = TRUE),
                import_format_named = switch(case$type, LFQ = grepl("MaxQuant proteinGroups export", t, fixed = TRUE), TMT = grepl("Proteome Discoverer", t, fixed = TRUE), DIANN = grepl("DIA-NN", t, fixed = TRUE)),
                da_file_name = grepl("DA_results_proteomics.xlsx", t, fixed = TRUE),
                old_de_file_name = grepl("DE_proteins_results.xlsx", t, fixed = TRUE),
                consensus_clusters = grepl("consensus clusters|pair of time points", t),
                analyst_field = grepl("[To be completed by the analyst", t, fixed = TRUE),
                software_block_last = length(sections) > 0 && identical(sections[[length(sections)]], "Software and build"),
                section_order = sections,
                version_commit_printed = grepl("Source commit:", t, fixed = TRUE),
                qvalue_sentence = grepl("q-values less than", t, fixed = TRUE),
                percentile_sentence = regmatches(t, regexpr("below the [0-9]+(st|nd|rd|th) percentile", t)))
        }
    })
    out
}
results <- list()
for (name in names(cases)) {
    results[[name]] <- tryCatch(run_case(name, cases[[name]]), error = function(e) list(error = conditionMessage(e)))
    r <- results[[name]]
    cat(name, "|", if (!is.null(r$error)) paste("ERROR:", r$error) else paste("renders", paste(names(r$renders), vapply(r$renders, function(x) if (is.null(x$error)) paste0("ok media=", x$media) else paste("ERR", x$error), ""), collapse = "; "),
        "| canned", length(r$checks$canned_sentences), "| placeholders", length(r$checks$placeholders), "| warnings", r$checks$printed_warnings,
        "| fig2", r$checks$figure2_caption, "| enrich not-run", r$checks$enrichment_not_run_sentence, "methods-absent", r$checks$enrichment_methods_absent, "with-plots", r$checks$enrichment_methods_present_with_plots,
        "| format", r$checks$import_format_named, "| DA file", r$checks$da_file_name, "| last", r$checks$software_block_last), "\n")
}
results$template_sha256 <- lapply(stats::setNames(list.files(template_dir, "[.]rmd$"), list.files(template_dir, "[.]rmd$")), function(f) digest::digest(file = file.path(template_dir, f), algo = "sha256"))
results$library <- list(path = find.package("MultiScholaR"), source_commit = utils::packageDescription("MultiScholaR")$`Config/MultiScholaR/SourceCommit` %||% NA_character_)
results$head_commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
jsonlite::write_json(results, file.path(evidence, paste0("render-check-", label, ".json")), auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null")
unlink(work, recursive = TRUE)
cat("written", file.path(evidence, paste0("render-check-", label, ".json")), "\n")
