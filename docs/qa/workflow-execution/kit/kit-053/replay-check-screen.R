# KIT-053 replay check (in process) for the DIA_LIMPA route: the record of the route's oracle store on a corpus
# (kit-053/oracle-record.R) is written and read back; the replay runs on a reflinked copy whose store is cut to the
# import and design states (kitPrepareReplayProject, as KIT-028's check does) and drives the route's own step list
# (the record's operations: the peptide stage, the protein stage, the correlation, the RUV covariates, then the DA
# with dpcDE); every state's fingerprint is compared with the origin's.
# Usage: [KIT053_INSTALLED=1] R_LIBS=<library> Rscript docs/qa/workflow-execution/kit/kit-053/replay-check-screen.R <corpus key> [label]
if (identical(Sys.getenv("KIT053_INSTALLED"), "1")) suppressMessages(library(MultiScholaR)) else suppressMessages(devtools::load_all(quiet = TRUE))
started <- Sys.time()
args <- commandArgs(TRUE)
corpus_key <- args[[1L]]
label <- if (length(args) > 1L) args[[2L]] else "tree"
repo <- normalizePath(".")
`%||%` <- function(a, b) if (is.null(a)) b else a
h <- new.env(parent = globalenv())
sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), h)
b <- jsonlite::read_json(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
ns <- asNamespace("MultiScholaR"); api <- function(name) get(name, ns)
# KIT-057: another ticket's oracle record and output directory under `KIT_LIMPA_TICKET_DIR` (default kit-053).
ticket_dir <- Sys.getenv("KIT_LIMPA_TICKET_DIR", "kit-053")
recorded <- readRDS(file.path(repo, "docs/qa/workflow-execution/kit", ticket_dir, paste0("oracles-limpa_", corpus_key, ".rds")))
stopifnot(identical(recorded$route, "DIA_LIMPA"), dir.exists(recorded$store))
origin_root <- dirname(recorded$store)
# The project root is the directory holding the experiment folder (`dia_baseline_smoke/...`), whatever the depth of
# the store below it (the recorder's setupDirectories layout, or a kit member's GUI layout).
origin_project <- origin_root
while (!dir.exists(file.path(origin_project, "dia_baseline_smoke")) && nchar(origin_project) > 1L) origin_project <- dirname(origin_project)
stopifnot(dir.exists(file.path(origin_project, "dia_baseline_smoke")))
store_relative <- substring(origin_root, nchar(origin_project) + 2L)
work <- file.path(b$harness$root, "kit053", paste0("replay-", corpus_key, "-", label)); unlink(work, recursive = TRUE); dir.create(work, recursive = TRUE)
copy_project <- function(name) {
    target <- file.path(work, name); dir.create(target)
    stopifnot(system2("cp", c("-a", "--reflink=auto", file.path(origin_project, "."), target)) == 0L)
    target
}
project_paths <- function(project) {
    # The recorder's project is laid out by setupDirectories (results/proteomics_<label>, scripts/proteomics_<label>);
    # a kit member's by the GUI (results/proteomics, scripts/proteomics): the one that exists is used.
    base <- file.path(project, "dia_baseline_smoke")
    leaf <- if (dir.exists(file.path(base, "scripts", "proteomics_dia_baseline_smoke"))) "proteomics_dia_baseline_smoke" else "proteomics"
    results <- file.path(base, "results", leaf)
    list(base_dir = base, results_dir = results, source_dir = file.path(base, "scripts", leaf),
        protein_qc_dir = file.path(results, "protein_qc"), peptide_qc_dir = file.path(results, "peptide_qc"),
        da_output_dir = file.path(results, "da_proteins"), publication_graphs_dir = file.path(results, "publication_graphs"),
        pathway_dir = file.path(results, "pathway_enrichment"), time_dir = file.path(results, "publication_graphs"),
        integration_dir = file.path(base, "integration"))
}
out <- list(schema = "kit.replay_check.v1", ticket = "KIT-053", route = "DIA_LIMPA", corpus = corpus_key, label = label, date = format(Sys.Date()),
    head_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE), origin = recorded$store, origin_build = recorded$build,
    library = list(path = find.package("MultiScholaR"), source_commit = utils::packageDescription("MultiScholaR")$`Config/MultiScholaR/SourceCommit` %||% NA_character_))
origin_copy <- copy_project("origin")
origin_copy_root <- file.path(origin_copy, store_relative)
identity <- api("readSavedStoreIdentity")(file.path(origin_copy_root, ".workflow-entities", "store.json"))
manager <- api("newDIAStateAdapter")(api("newWorkflowEntityStore")(origin_copy_root, identity$project_id, identity$workflow_id), reopen = TRUE)
json <- file.path(work, "analysis_record.json")
api("writeWorkflowAnalysisRecord")(manager, json)
api("closeDIAStateAdapter")(manager)
analysis <- jsonlite::read_json(json, simplifyVector = FALSE)
states <- vapply(analysis$steps, function(s) s$state, "")
operations <- vapply(analysis$steps, function(s) s$operation, "")
preconditions <- states[seq_len(sum(cumprod(operations %in% c("import", "design_builder"))))]
out$record <- list(workflow = analysis$workflow, operations = operations, states = states, valid = api("validateWorkflowAnalysisRecord")(analysis)$valid,
    da_model = analysis$analyses$differential_abundance$model[c("model", "formula_used", "k", "block_reason", "route")])
corpus <- b$corpora[[corpus_key]]$files
inputs <- c(diann_report = corpus[["report.tsv"]]$path, fasta = corpus[["proteins.fasta"]]$path)
project <- copy_project("replay")
root <- h$kitPrepareReplayProject(file.path(project, store_relative), preconditions, api)
scratch <- file.path(work, "replay-scratch"); dir.create(scratch)
t0 <- Sys.time()
result <- h$kitReplayJourney(analysis, root, origin_copy_root, inputs, project_paths(project), scratch, api)
out$replay <- list(seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")), states = result$states, mismatches = result$mismatches,
    replayed = vapply(result$per_step, function(x) paste(x$operation, x$state), ""))
file.copy(json, file.path(repo, "docs/qa/workflow-execution/kit", ticket_dir, paste0("analysis_record-limpa-", corpus_key, "-", label, ".json")), overwrite = TRUE)
cat("RECORD:", out$record$workflow, "| ops", paste(operations, collapse = ","), "| valid", out$record$valid, "| DA model", paste(unlist(out$record$da_model), collapse = "/"), "\n")
cat("REPLAY:", length(result$states), "states,", length(result$mismatches), "mismatches", if (length(result$mismatches)) paste0("(", paste(result$mismatches, collapse = ", "), ")") else "",
    "| DA replayed", any(grepl("^DA-contrast-output", out$replay$replayed)), "|", round(out$replay$seconds), "s\n")
out$seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
jsonlite::write_json(out, file.path(repo, "docs/qa/workflow-execution/kit", ticket_dir, paste0("replay-check-", corpus_key, "-", label, ".json")), auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = NA)
cat("TOTAL", round(out$seconds), "s\n")
