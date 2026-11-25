# KIT-243 G0 (ticket section 10): base values observed on the package at the branch point (v0.9-staging 16112917), before
# any test or fix is written. Evidence code under docs/qa: no package code. Runs on a GitHub-hosted runner (run-tests
# `scripts`) from the package root of the commit under test, whose R/ is the base's.
#
# Observed, per omic: whether `analysis_record.json` exists after "Save Workflow Arguments" with the real writers (P03);
# the record's steps, as (state, result_entity_id), and its DA block against a fresh build from the store right after the
# Save, both through one JSON round trip with the record writer's options (ticket 2.4; G1, stop rule 4), and again after a
# DA re-run with identical settings in the same session (G6); which DA figures each publication writes on the fixtures
# the E rows use (E1 to E4b, E12; F172 draws no heatmap, KIT-212 G1).
#
# The helpers are KIT-212's and KIT-112's (tests/testthat/test-kit212-da-current.R, test-workflow-saved-da.R at the commit
# under test): their top-level assignments only are evaluated; no test of those files runs.
#
# Usage: Rscript docs/qa/workflow-execution/kit/kit-243/g0-base-values.R <out dir>
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
dir.create(args[[1L]], showWarnings = FALSE, recursive = TRUE)
out <- normalizePath(args[[1L]], mustWork = TRUE)
suppressMessages(pkgload::load_all(".", export_all = TRUE, helpers = TRUE, quiet = TRUE))
library(testthat)
g0 <- new.env(parent = pkgload::ns_env("MultiScholaR"))
for (file in c("test-kit212-da-current.R", "test-workflow-saved-da.R")) {
    for (expression in as.list(parse(file.path("tests", "testthat", file), keep.source = FALSE))) {
        if (is.call(expression) && identical(expression[[1L]], as.name("<-")) && is.name(expression[[2L]])) eval(expression, g0)
    }
}
rows <- list()
note <- function(probe, field, value) {
    line <- c(probe = probe, field = field, value = gsub("[\t\n]", " ", paste(format(value), collapse = ",")))
    rows[[length(rows) + 1L]] <<- line
    cat(sprintf("KIT-243 G0\t%s\t%s\t%s\n", line[[1L]], line[[2L]], line[[3L]]))
}

g0$k243Identity <- function(record) list(steps = lapply(record$steps, function(step)
    list(state = step$state, result_entity_id = step$result_entity_id)), differential_abundance = record$analyses$differential_abundance)
g0$k243Trip <- function(x) jsonlite::fromJSON(jsonlite::toJSON(x, auto_unbox = TRUE, digits = I(17), null = "null", na = "null"),
    simplifyVector = FALSE)
g0$k243FirstDiff <- function(a, b, path = "") {
    if (identical(a, b)) return(NA_character_)
    if (is.list(a) && is.list(b)) {
        if (!identical(length(a), length(b)) || !identical(names(a), names(b))) return(sprintf("%s: lengths %d/%d; names %s / %s", path,
            length(a), length(b), paste(names(a), collapse = ","), paste(names(b), collapse = ",")))
        for (i in seq_along(a)) {
            label <- if (!is.null(names(a)) && nzchar(names(a)[[i]])) names(a)[[i]] else as.character(i)
            found <- k243FirstDiff(a[[i]], b[[i]], paste0(path, "/", label))
            if (!is.na(found)) return(found)
        }
        return(sprintf("%s: attributes differ", path))
    }
    substr(sprintf("%s: saved %s; current %s", path, paste(deparse(a), collapse = " "), paste(deparse(b), collapse = " ")), 1L, 600L)
}
#' The on-disk record against the store (2.4): "<no record>", "equal" or the first difference.
g0$k243Compare <- function(manager, source_dir) {
    path <- file.path(source_dir, "analysis_record.json")
    if (!file.exists(path)) return("<no record>")
    saved <- k243Trip(k243Identity(jsonlite::read_json(path, simplifyVector = FALSE)))
    current <- k243Trip(k243Identity(withWorkflowStateScope(manager, function() buildWorkflowAnalysisRecord(manager))))
    found <- k243FirstDiff(saved, current)
    if (is.na(found)) "equal" else found
}
g0$k243Files <- function(dir) {
    if (is.null(dir) || !dir.exists(dir)) return("<none>")
    files <- sort(list.files(dir, recursive = TRUE, all.files = TRUE))
    if (!length(files)) "<none>" else paste(files, collapse = ",")
}
#' The proteomics Save with the real writers (P03): KIT-212's k212ProtSave without its two writer overrides.
g0$k243ProtSave <- function(fixture, dirs) {
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    seen$failed <- list()
    seen$done <- list()
    protNormCoreQuiet(saveWorkflowSummaryParameters(omicGuiConfig("proteomics")$summary, fixture$workflow, dirs,
        "proteomics", "k243", "d", notify = function(message, type = NULL, duration = NULL) {
            seen$notes <- c(seen$notes, as.character(message))
            invisible(NULL)
        }, done = function(saved) seen$done <- c(seen$done, list(saved)),
        failed = function(error) seen$failed <- c(seen$failed, list(error)), catFn = function(...) NULL))
    seen
}

# Instrument repair (run 45588535137, VAL-54): the functions above were closures of the global environment, so
# `k243Compare` could not find `k243Trip`; each is re-enclosed in the helpers' environment.
for (name in c("k243Identity", "k243Trip", "k243FirstDiff", "k243Compare", "k243Files", "k243ProtSave")) {
    f <- get(name, envir = g0)
    environment(f) <- g0
    assign(name, f, envir = g0)
}

#' One probe: its outcome ("ok", the error, or the skip) recorded; a probe's failure does not stop the others.
probe <- function(name, code) {
    outcome <- tryCatch({
        code()
        "ok"
    }, skip = function(condition) paste("skipped:", conditionMessage(condition)),
        error = function(error) paste("error:", conditionMessage(error)))
    note(name, "probe", outcome)
}

eval(quote(for (class in c("lipid", "metab")) probe(paste0(class, "-in-process"), function() {
    root <- withr::local_tempdir()
    seen <- k194Notes()
    omic <- k194Omic(class)
    q <- k207Normalised(file.path(root, "g0"), class, da = TRUE)
    owner <- k194Owner(q$fx)
    dirs <- k212Dirs(file.path(root, "project"), omic)
    k212Publish(q$fx$binding, owner, dirs[[omic]])
    note(class, "publication figures (DA1)", k243Files(dirs[[omic]]$publication_graphs_dir))
    note(class, "made figures (DA1)", basename(withWorkflowStateScope(owner, function()
        workflowDaMadeFiles(owner, dirs[[omic]], omic))))
    rec <- k207Save(q$fx$workflow, omic, dirs)
    note(class, "save done", rec$done)
    note(class, "record exists after Save", file.exists(file.path(dirs[[omic]]$source_dir, "analysis_record.json")))
    note(class, "G1 record vs store after Save", k243Compare(owner, dirs[[omic]]$source_dir))
    k207Da(q$fx)
    note(class, "G6 record vs store after an identical DA re-run", k243Compare(owner, dirs[[omic]]$source_dir))
    kit172Quiet(runWorkflowOmicsReferenceDaAnalysis(q$fx$binding, owner, k194Ref(q$fx),
        list(formula_string = "~ 0 + group", da_q_val_thresh = 0.1, treat_lfc_cutoff = 0), kit172Contrasts()))
    note(class, "E5 record vs store after a DA at q 0.1", k243Compare(owner, dirs[[omic]]$source_dir))
})), g0)

eval(quote(probe("trx-in-process", function() {
    stopifnot(requireNamespace("edgeR", quietly = TRUE))
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    project <- module_ci_trx_project(module_ci_trx_corpus())
    owner <- project$owner
    k194TrxRun(project, omicNormReactiveState(), list(ruv_mode = "manual", ruv_k = 2L, ruv_percentage = 10), seen)
    note("trx", "head after run", owner$current_state)
    kit172Quiet(runTrxReferenceDaAnalysis(owner, getWorkflowStateReference(owner), list(da_q_val_thresh = 0.05, treat_lfc_cutoff = 0),
        trxDaContrastEntries(project$workflow$contrasts_tbl)))
    k212Publish(transcriptomicsDeclarationBinding, owner, project$paths)
    note("trx", "publication figures (DA1)", k243Files(project$paths$publication_graphs_dir))
    rec <- k207Save(project$workflow, "transcriptomics", project$project_dirs)
    note("trx", "save done", rec$done)
    source_dir <- project$project_dirs$transcriptomics$source_dir
    note("trx", "record exists after Save", file.exists(file.path(source_dir, "analysis_record.json")))
    note("trx", "G1 record vs store after Save", k243Compare(owner, source_dir))
})), g0)

eval(quote(probe("prot-k212-in-process", function() {
    stopifnot(workflowOwnerCapable())
    fixture <- k212Prot()
    k212Export(fixture, fixture$norm)
    k212ProtDa(fixture, publish = TRUE)
    note("prot-k212", "publication figures (DA1)", k243Files(fixture$paths$publication_graphs_dir))
    dirs <- k212ProtDirs(fixture, "g0")
    seen <- k243ProtSave(fixture, dirs)
    note("prot-k212", "save done; failed", c(length(seen$done), vapply(seen$failed, conditionMessage, "")))
    note("prot-k212", "save notes", seen$notes)
    note("prot-k212", "record exists after Save", file.exists(file.path(dirs$proteomics$source_dir, "analysis_record.json")))
    note("prot-k212", "G1 record vs store after Save", k243Compare(fixture$manager, dirs$proteomics$source_dir))
    k212ProtDa(fixture)
    note("prot-k212", "G6 record vs store after an identical DA re-run", k243Compare(fixture$manager, dirs$proteomics$source_dir))
})), g0)

eval(quote(probe("prot-saved-da-in-process", function() {
    stopifnot(requireNamespace("Glimma", quietly = TRUE))
    fixture <- savedDaW1Owner("kit243-g0")
    contrasts <- savedDaW1Contrasts(1L)
    fixture$workflow$contrasts_tbl <- contrasts
    fixture$da$contrasts_available <- contrasts$contrasts
    graphs <- fixture$paths$publication_graphs_dir
    barplot <- file.path(graphs, "NumSigDaMolecules", "da_proteins_num_sig_da_molecules.png")
    run <- function(qvalue) suppressWarnings(suppressMessages(utils::capture.output(runDIAReferenceDaAnalysis(
        list(formula_string = "~ 0 + group", da_q_val_thresh = qvalue, treat_lfc_cutoff = 0), fixture$da, fixture$workflow,
        fixture$paths, contrasts, NULL))))
    run(0.05)
    note("prot-saved-da", "publication figures (DA1, q 0.05)", k243Files(graphs))
    before <- if (file.exists(barplot)) digest::digest(file = barplot, algo = "sha256") else "<absent>"
    run(1e-300)
    note("prot-saved-da", "publication figures (DA2, q 1e-300)", k243Files(graphs))
    after <- if (file.exists(barplot)) digest::digest(file = barplot, algo = "sha256") else "<absent>"
    note("prot-saved-da", "count barplot sha256 DA1, DA2", c(before, after))
})), g0)

utils::write.table(do.call(rbind, rows), file.path(out, "g0-base-values.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
