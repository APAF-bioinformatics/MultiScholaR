# KIT-205 A3 (ii): a record made on the base commit, whose DA has a comparison with annotated untested proteins,
# replays on this commit (head) with the earlier background, byte for byte; a fresh head run on the same input differs
# from it exactly as control C (b) predicts. Evidence code under docs/qa: no package code. Runs on a GitHub-hosted runner
# (run-tests `scripts`) from the package root of the commit under test (head), after head is installed.
#
# Usage: Rscript docs/qa/workflow-execution/kit/kit-205/replay-k205.R <out dir> <base commit>
#
# The fixture: KIT-173's T7 project (helper-go-propagation.R: the committed E. coli fixture of 70 proteins, the annotation
# through getUniprotAnnotations with the download stubbed, the clusterProfiler route on taxon 999999, the real kernel, the
# tab's defaults) with one change: six proteins (the first four `fill` and the first two `p_term` entries of
# t7-fixture-entries.tsv, by name) have no value in group A, so the two-group DA gives them no p-value.
# 1. The base commit's tree is installed into a private library; this script in `--capture` mode runs the project on
#    it and writes the store, the record, the evidence (stored result tables, plot inputs, pathway files' SHA-256), the
#    product and the DA table's untested count.
# 2. On head (this process): (i) the base record is replayed on a copy of that store (helper goPropReplay, the kernel the
#    GUI selects outside test mode) and compared with the base evidence; (ii) the same project is run fresh on head with
#    the kernel's arguments recorded; its input digests must equal the base record's, its results must differ from
#    base's, and they must equal processEnrichments() under the earlier rule on the recorded DA tables without their
#    untested rows (control C (b)'s prediction).
# Writes <out dir>/replay-k205.tsv (check, pass, detail), replay-k205-base.tsv and the logs. Exit status 0.

UNTESTED_ROLES <- list(fill = 4L, p_term = 2L)

#' KIT-173's goPropRunProject with the untested proteins' group-A values removed
k205RunProject <- function(h, root, annotations, process_fn = processEnrichments) {
    roles <- h$goPropReadTsv("t7-fixture-entries.tsv")
    paths <- h$goPropPaths(root)
    for (p in paths) dir.create(p, recursive = TRUE, showWarnings = FALSE)
    saveRDS(annotations, file.path(paths$source_dir, "uniprot_dat_cln.RDS"))
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "go-propagation", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    on.exit(closeDIAStateAdapter(manager), add = TRUE)
    manager$.context$reference_states <- TRUE
    manager$setWorkflowType("DIA")
    entries <- sort(roles$entry)
    untested <- unlist(lapply(names(UNTESTED_ROLES), function(role)
        sort(roles$entry[roles$role == role])[seq_len(UNTESTED_ROLES[[role]])]))
    design <- data.frame(Run = paste0("S", 1:6), group = rep(c("A", "B"), each = 3L),
        batch = rep(c("b1", "b2", "b3"), times = 2L), stringsAsFactors = FALSE)
    quant <- h$goPropQuant(entries, roles$entry[roles$role == "p_term"], roles$entry[roles$role == "d_term"])
    quant[quant$uniprot_acc %in% untested, c("S1", "S2", "S3")] <- NA_real_
    object <- methods::new("ProteinQuantitativeData", protein_quant_table = quant, design_matrix = design,
        sample_id = "Run", group_id = "group", protein_id_column = "uniprot_acc",
        protein_id_table = data.frame(uniprot_acc = entries, gene_names = paste0("GENE", seq_along(entries)),
            stringsAsFactors = FALSE),
        args = list(globalParameters = list(workflow_type = "DIA", use_limpa = FALSE, report_template = "DIA_report"),
            deAnalysisParameters = list(formula_string = "~ 0 + group")))
    setDIAStepProvenance(manager, "design_builder")
    manager$saveState("correlation_filtered", object, list(), "fixture")
    manager$.context$experiment_paths <- paths
    contrasts <- data.frame(contrasts = "groupB-groupA", full_format = "B_vs_A=groupB-groupA",
        friendly_names = "B_vs_A", stringsAsFactors = FALSE)
    workflow <- list2env(list(state_manager = manager, contrasts_tbl = contrasts, tab_status = list(),
        qc_params = list(protein_qc = list()), config_list = list()))
    parameters <- list(selected_contrast = "B_vs_A", organism_taxid = "999999", up_cutoff = 0, down_cutoff = 0,
        q_cutoff = 0.05, correction_method = "gSCS", enable_organism_filter = FALSE)
    withDIAStateScope(manager, function() {
        da_data <- new.env(parent = emptyenv())
        da_data$current_s4_object <- getDIAStateReference(manager)
        da_data$contrasts_available <- contrasts$contrasts
        manifest <- runDIAReferenceDaAnalysis(list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05,
            treat_lfc_cutoff = 0), da_data, workflow, paths, contrasts = contrasts, publish = FALSE)
        produceDIADaPublication(manager, manifest, annotations, 0.05, paths[c("da_output_dir", "publication_graphs_dir")],
            save_assets = FALSE)
        runWorkflowEnrichment(workflow, omicGuiConfig("proteomics")$enrichment, manifest, parameters, paths,
            process_fn = process_fn)
    })
    record <- file.path(root, "analysis_record.json")
    writeWorkflowAnalysisRecord(manager, record)
    state <- withDIAStateScope(manager, function() manager$getState("enrichment_completed"))
    da_files <- list.files(paths$da_output_dir, "_long_annot\\.tsv$", full.names = TRUE)
    da <- lapply(da_files, function(f) utils::read.delim(f, stringsAsFactors = FALSE))
    untested_rows <- sum(vapply(da, function(t) sum(!is.finite(t$raw_pvalue)), integer(1)))
    annotated <- unique(as.character(buildUniprotGoGeneSets(annotations, goPropagationShipped)$term2gene[[2L]]))
    untested_annotated <- sum(vapply(da, function(t) sum(!is.finite(t$raw_pvalue) & t$uniprot_acc %in% annotated),
        integer(1)))
    list(record = record, paths = paths, state = state, product = h$goPropProduct(manager), untested = untested,
        da_files = basename(da_files), untested_rows = untested_rows, untested_annotated = untested_annotated)
}

helpers <- function() {
    h <- new.env(parent = asNamespace("MultiScholaR"))
    sys.source(file.path("tests", "testthat", "helper-go-propagation.R"), envir = h)
    h
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) && identical(args[[1L]], "--capture")) {
    capture <- args[[2L]]
    suppressPackageStartupMessages(library(MultiScholaR))
    environment(k205RunProject) <- asNamespace("MultiScholaR") # the package's internals, as the helper sees them
    Sys.setenv(MULTISCHOLAR_TEST_MODE = "false")
    h <- helpers()
    work <- file.path(tempdir(), "k205-capture")
    dir.create(work, recursive = TRUE)
    annotations <- h$goPropAnnotate(h$goPropT7Raw(), file.path(work, "cache"))
    root <- file.path(work, "project")
    utils::capture.output(result <- suppressMessages(suppressWarnings(k205RunProject(h, root, annotations))))
    saveRDS(h$goPropRunEvidence(result$state, result$paths$pathway_dir), file.path(capture, "evidence.rds"))
    saveRDS(list(parameters = result$product$parameters, service = result$product$service,
        contrasts = result$product$contrasts), file.path(capture, "product.rds"))
    file.copy(result$record, file.path(capture, "record.json"))
    owd <- setwd(root)
    utils::tar(file.path(capture, "project.tar"), files = list.files(".", recursive = TRUE, all.files = TRUE),
        compression = "none", tar = "internal")
    setwd(owd)
    utils::write.table(data.frame(measure = c("package_path", "has_enrichmentTestedRows", "untested_proteins",
        "da_files", "untested_rows", "untested_annotated_rows"), value = c(find.package("MultiScholaR"),
        exists("enrichmentTestedRows", envir = asNamespace("MultiScholaR")), paste(result$untested, collapse = ","),
        paste(result$da_files, collapse = ","), result$untested_rows, result$untested_annotated)),
        file.path(capture, "capture.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    quit(status = 0L)
}

stopifnot(length(args) == 2L)
dir.create(args[[1L]], showWarnings = FALSE, recursive = TRUE)
out <- normalizePath(args[[1L]], mustWork = TRUE)
base <- args[[2L]]
work <- file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "k205-replay")
src <- file.path(work, "base-src")
lib <- file.path(work, "base-lib")
capture <- file.path(work, "capture")
for (d in c(src, lib, capture)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]))
rscript <- file.path(R.home("bin"), "Rscript")
base_check <- list(base_commit = base)
base_check$archive_exit <- system(sprintf("git archive %s | tar -x -C %s", shQuote(base), shQuote(src)))
base_check$install_exit <- system2(file.path(R.home("bin"), "R"), c("CMD", "INSTALL", "--no-test-load",
    paste0("--library=", shQuote(lib)), shQuote(src)), stdout = file.path(out, "replay-base-install.log.txt"),
    stderr = file.path(out, "replay-base-install.log.txt"))
base_check$capture_exit <- system2(rscript, c(shQuote(script), "--capture", shQuote(capture)),
    stdout = file.path(out, "replay-capture.log.txt"), stderr = file.path(out, "replay-capture.log.txt"),
    env = paste0("R_LIBS=", lib))
captured <- if (file.exists(file.path(capture, "capture.tsv"))) utils::read.delim(file.path(capture, "capture.tsv"),
    stringsAsFactors = FALSE) else data.frame(measure = character(), value = character())
utils::write.table(rbind(data.frame(measure = names(base_check), value = vapply(base_check, as.character, "")), captured),
    file.path(out, "replay-k205-base.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

checks <- list()
check <- function(id, pass, detail = "") checks[[length(checks) + 1L]] <<- data.frame(check = id, pass = isTRUE(pass),
    detail = substr(gsub("[\t\r\n]+", " ", paste(detail, collapse = " ")), 1L, 300L), stringsAsFactors = FALSE)
value <- function(measure) captured$value[match(measure, captured$measure)]
result <- tryCatch({
    check("base capture ran the base code", identical(value("has_enrichmentTestedRows"), "FALSE") &&
        startsWith(value("package_path") %||% "", lib), value("package_path"))
    check("the base DA has an annotated untested row", isTRUE(as.integer(value("untested_annotated_rows")) > 0L),
        sprintf("untested rows %s, annotated %s", value("untested_rows"), value("untested_annotated_rows")))
    suppressPackageStartupMessages(library(MultiScholaR))
    environment(k205RunProject) <- asNamespace("MultiScholaR")
    Sys.setenv(MULTISCHOLAR_TEST_MODE = "false")
    h <- helpers()
    record <- jsonlite::read_json(file.path(capture, "record.json"), simplifyVector = FALSE)
    base_evidence <- readRDS(file.path(capture, "evidence.rds"))
    base_step <- Filter(function(s) identical(s$operation, "enrichment"), record$steps)[[1L]]
    check("the base record has no enrichment_background",
        !"enrichment_background" %in% vapply(base_step$parameters, `[[`, "", "name"))
    # (i) the base record replayed on head
    scratch <- file.path(work, "replay")
    project <- file.path(scratch, "origin")
    dir.create(project, recursive = TRUE)
    utils::untar(file.path(capture, "project.tar"), exdir = project, tar = "internal")
    replayed <- NULL
    utils::capture.output(replayed <- suppressMessages(suppressWarnings(h$goPropReplay(project, record,
        file.path(scratch, "run")))))
    check("replay states", identical(replayed$states, c("initial", "correlation_filtered", "enrichment_completed")),
        replayed$states)
    check("replayed provenance omits the marker", is.null(replayed$provenance$fixed$enrichment_background))
    evidence <- h$goPropRunEvidence(replayed$state, replayed$pathway_dir)
    check("replay: result tables identical to base (exact)", identical(evidence$tables, base_evidence$tables))
    check("replay: plot inputs identical", identical(evidence$plot_inputs, base_evidence$plot_inputs))
    check("replay: pathway file names", identical(evidence$files$file, base_evidence$files$file), evidence$files$file)
    check("replay: every TSV's SHA-256 (and each plot file's existence)", identical(evidence$files$sha256,
        base_evidence$files$sha256))
    # (ii) a fresh head run on the same input
    kernel <- new.env()
    recorder <- function(...) {
        kernel$args <- list(...)
        processEnrichments(...)
    }
    annotations <- h$goPropAnnotate(h$goPropT7Raw(), file.path(work, "cache"))
    fresh <- NULL
    utils::capture.output(fresh <- suppressMessages(suppressWarnings(k205RunProject(h, file.path(work, "fresh"),
        annotations, process_fn = recorder))))
    fresh_record <- jsonlite::read_json(fresh$record, simplifyVector = FALSE)
    fresh_step <- Filter(function(s) identical(s$operation, "enrichment"), fresh_record$steps)[[1L]]
    marker <- Filter(function(p) identical(p$name, "enrichment_background"), fresh_step$parameters)
    check("fresh record carries the marker (fixed, tested)", length(marker) == 1L &&
        identical(marker[[1L]]$source, "fixed") && identical(MultiScholaR:::decodeWorkflowRecordValue(marker[[1L]]$value), "tested"))
    check("same input: DA tables' digests equal the base record's", identical(fresh_step$outcomes$da_inputs_sha256,
        base_step$outcomes$da_inputs_sha256))
    fresh_evidence <- h$goPropRunEvidence(fresh$state, fresh$paths$pathway_dir)
    check("fresh head results differ from base's", !identical(fresh_evidence$tables, base_evidence$tables))
    reduce <- function(table) {
        tested <- is.finite(table$fdr_qvalue)
        if ("raw_pvalue" %in% names(table)) tested <- tested & is.finite(table$raw_pvalue)
        if (all(tested)) table else table[tested, , drop = FALSE]
    }
    arguments <- kernel$args
    left_out <- vapply(arguments$da_results@da_data, function(t) nrow(t) - nrow(reduce(t)), integer(1))
    arguments$da_results@da_data <- lapply(arguments$da_results@da_data, reduce)
    arguments$pathway_dir <- file.path(work, "predicted")
    dir.create(arguments$pathway_dir)
    predicted <- NULL
    utils::capture.output(predicted <- suppressMessages(suppressWarnings(withr::with_options(
        list(MultiScholaR.enrichment_background = "all_rows"), do.call(processEnrichments, arguments)))))
    predicted_tables <- lapply(predicted@enrichment_data, function(contrast) lapply(contrast, function(r)
        if (is.null(r)) NULL else r@result))
    check("fresh head = the earlier rule on the tables without their untested rows (control C (b))",
        identical(fresh_evidence$tables, predicted_tables), sprintf("rows left out per comparison: %s",
            paste(names(left_out), left_out, collapse = ", ")))
    differs <- unlist(lapply(names(base_evidence$tables), function(contrast) lapply(c("up", "down"), function(d) {
        b <- base_evidence$tables[[contrast]][[d]]
        f <- fresh_evidence$tables[[contrast]][[d]]
        if (is.null(b) || is.null(f)) return(sprintf("%s/%s: base %s rows, head %s rows", contrast, d,
            if (is.null(b)) "no" else nrow(b), if (is.null(f)) "no" else nrow(f)))
        both <- intersect(b$ID, f$ID)
        sprintf("%s/%s: %d and %d terms, %d shared, largest |log10 p| change %.4g, BgRatio %s -> %s", contrast, d,
            nrow(b), nrow(f), length(both), if (length(both)) max(abs(log10(b$pvalue[match(both, b$ID)]) -
                log10(f$pvalue[match(both, f$ID)]))) else NA_real_, b$BgRatio[[1L]], f$BgRatio[[1L]])
    })))
    check("the change, base record against fresh head (reported)", TRUE, differs)
    "ok"
}, error = function(e) paste("instrument:", conditionMessage(e)))
table <- if (length(checks)) do.call(rbind, checks) else data.frame(check = character(), pass = logical(),
    detail = character())
if (!identical(result, "ok")) table <- rbind(table, data.frame(check = "script", pass = FALSE, detail = result))
utils::write.table(table, file.path(out, "replay-k205.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(table)
