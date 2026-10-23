# KIT-053: the DIA route's limma DA is unchanged on this build. On a copy of the kit-journey-028 member's project
# (the DIA route, screen corpus, built on kit_027), the DA is re-run in process from the member's own record
# parameters and contrasts, and every contrast's tables are compared with the tables that member accepted.
# Usage: KIT053_INSTALLED=1 R_LIBS=<library> Rscript docs/qa/workflow-execution/kit/kit-053/dia-da-unchanged.R <project copy> <label>
suppressMessages(library(MultiScholaR)); ns <- asNamespace("MultiScholaR")
for (f in ls(ns)) if (!exists(f, envir = globalenv(), inherits = TRUE) && is.function(get(f, ns))) assign(f, get(f, ns), envir = globalenv())
a <- commandArgs(TRUE); project <- a[[1]]; label <- a[[2]]
harness <- new.env(parent = globalenv()); sys.source("tools/profiling/workflow_kit.R", harness)
store <- harness$kitFindStore(project); root <- dirname(store)
id <- readSavedStoreIdentity(file.path(root, ".workflow-entities", "store.json"))
manager <- newDIAStateAdapter(newWorkflowEntityStore(root, id$project_id, id$workflow_id), reopen = TRUE)
sha <- function(df) { df <- as.data.frame(df); rownames(df) <- NULL; digest::digest(list(names(df), lapply(df, function(col) { attributes(col) <- NULL; col })), algo = "sha256") }
tables_of <- function(manifest) withDIAStateScope(manager, function() {
  out <- list(); for (k in names(manifest$results)) { r <- manifest$results[[k]]
    for (t in c("da_proteins_long", "da_proteins_wide", "contrasts_results_table")) {
      v <- r$tables[[t]]; if (is.null(v)) next
      out[[paste(k, t, sep = "/")]] <- sha(if (isWorkflowTableReference(v)) materializeDIATable(manager$.context, v, NULL) else if (isSavedWorkflowReference(v)) hydrateSavedWorkflowEntity(manager$.context$store, v)$value else v) } }
  out })
before <- withDIAStateScope(manager, function() resolveDIAStateValue(manager$.context, manager$.context$snapshot$da_results))
old <- tables_of(before)
# The journey's DA settings (the binding the member ran with) and the contrasts of the member's accepted manifest.
journey_da <- jsonlite::read_json("docs/qa/workflow-execution/kit/bindings.json")$journey$da
settings <- list(formula = journey_da$formula, qvalue = journey_da$q_value, lfc = journey_da$treat_lfc)
contrasts <- as.data.frame(before$contrasts, stringsAsFactors = FALSE)
withDIAStateScope(manager, function() {
  da_data <- new.env(parent = emptyenv()); da_data$current_s4_object <- getDIAStateReference(manager); da_data$contrasts_available <- contrasts$contrasts
  wf <- new.env(parent = emptyenv()); wf$state_manager <- manager; wf$contrasts_tbl <- contrasts
  invisible(runDIAReferenceDaAnalysis(list(formula_string = settings$formula, da_q_val_thresh = settings$qvalue, treat_lfc_cutoff = settings$lfc),
    da_data, wf, list(), contrasts = contrasts, publish = FALSE)) })
after <- withDIAStateScope(manager, function() resolveDIAStateValue(manager$.context, manager$.context$snapshot$da_results))
new <- tables_of(after)
same <- identical(names(old), names(new)) && all(vapply(names(old), function(k) identical(old[[k]], new[[k]]), logical(1)))
out <- list(schema = "kit.dia_da_unchanged.v1", ticket = "KIT-053", label = label, date = format(Sys.Date()),
  head_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE), origin = "kit-journey-028/candidate-1 (copy)",
  parameters = settings, tables = length(old), identical = same, model_before = before$model, model_after = after$model)
jsonlite::write_json(out, file.path("docs/qa/workflow-execution/kit/kit-053", paste0("dia-da-unchanged-", label, ".json")), auto_unbox = TRUE, pretty = TRUE, null = "null")
cat("DIA DA tables:", length(old), "| identical:", same, "| model after:", after$model$model %||% "limma (no kernel model field)", "| route", after$model$route %||% "", "\n")
