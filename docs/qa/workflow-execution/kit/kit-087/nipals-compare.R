# KIT-087: our NIPALS against mixOmics on data already generated, with no chain and no GUI. Stores are read from a reflink copy.
# Usage (R_LIBS = the kit_087 library; run from the repository root):
#   nipals-compare.R extract <label> <project dir> <state or "before_ruv"> <ncomp> <output dir>
#   nipals-compare.R call <input rds> <old|new>
#   nipals-compare.R compare <input rds>
#   nipals-compare.R output <label> <project dir> <prefix> <output file>
suppressMessages(library(MultiScholaR))
args <- commandArgs(TRUE); mode <- args[[1L]]
ns <- asNamespace("MultiScholaR"); api <- function(name) get(name, ns)
for (f in ls(ns)) if (!exists(f, envir = globalenv(), inherits = TRUE) && is.function(get(f, ns))) assign(f, get(f, ns), envir = globalenv())
`%||%` <- function(a, b) if (is.null(a)) b else a
openCopy <- function(project) {
  work <- tempfile("kit087-", tmpdir = Sys.getenv("KIT087_WORK")); dir.create(work)
  stopifnot(system2("cp", c("-a", "--reflink=auto", project, work)) == 0L)
  h <- new.env(parent = globalenv()); sys.source("tools/profiling/workflow_kit.R", h)
  store_dir <- h$kitFindStore(file.path(work, basename(project))); root <- dirname(store_dir)
  identity <- readSavedStoreIdentity(file.path(store_dir, "store.json"))
  list(work = work, manager = newDIAStateAdapter(newWorkflowEntityStore(root, identity$project_id, identity$workflow_id), reopen = TRUE),
    states = names(h$kitReadHead(store_dir, api)$states))
}
assaysOf <- function(object) {
  if (inherits(object, "MetaboliteAssayData")) return(list(tables = object@metabolite_data, id = object@metabolite_id_column))
  if (inherits(object, "LipidomicsAssayData")) return(list(tables = object@lipid_data, id = object@lipid_id_column))
  list(tables = list(protein = object@protein_quant_table), id = object@protein_id_column)
}
if (mode == "extract") {
  label <- args[[2L]]; project <- args[[3L]]; state <- args[[4L]]; ncomp <- as.integer(args[[5L]]); out <- args[[6L]]
  opened <- openCopy(project)
  if (state == "before_ruv") state <- opened$states[[which(opened$states == "ruv_corrected") - 1L]]
  withDIAStateScope(opened$manager, function() {
    object <- opened$manager$getState(state); parts <- assaysOf(object)
    samples <- as.character(object@design_matrix[[object@sample_id]])
    for (assay in names(parts$tables)) {
      table <- as.data.frame(parts$tables[[assay]], check.names = FALSE)
      values <- as.matrix(table[intersect(samples, colnames(table))]); rownames(values) <- as.character(table[[parts$id]])
      values[!is.finite(values)] <- NA_real_
      as_stored <- anyNA(values)
      if (!as_stored) { # The declared pattern, in memory only: the lowest-abundance 30% of features, in 15% of their samples.
        set.seed(20261105L)
        for (feature in order(rowMeans(values))[seq_len(round(nrow(values) * 0.3))]) values[feature, sample.int(ncol(values), round(ncol(values) * 0.15))] <- NA_real_
      }
      empty_rows <- rowSums(!is.na(values)) == 0L
      Y <- t(values[!empty_rows, , drop = FALSE])
      path <- file.path(out, sprintf("%s__%s.rds", label, assay))
      saveRDS(list(Y = Y, ncomp = ncomp, label = label, assay = assay, state = state, missing = if (as_stored) "as stored" else "declared mask"), path, compress = FALSE)
      cat(sprintf("EXTRACT %s %s state %s samples %d features %d missing %.4f %s\n", label, assay, state, nrow(Y), ncol(Y), mean(is.na(Y)), if (as_stored) "stored" else "mask"))
    }
  })
  closeDIAStateAdapter(opened$manager); unlink(opened$work, recursive = TRUE)
} else if (mode == "call") {
  input <- readRDS(args[[2L]]); which_one <- args[[3L]]; gc(reset = TRUE)
  elapsed <- system.time(result <- suppressWarnings(suppressMessages(if (which_one == "old") mixOmics::impute.nipals(input$Y, ncomp = input$ncomp) else imputeNipalsMatrix(input$Y, ncomp = input$ncomp))))[["elapsed"]]
  saveRDS(result, sub("\\.rds$", paste0(".", which_one, ".result"), args[[2L]]), compress = FALSE)
  cat(sprintf("CALL %s %s %s elapsed %.1f s digest %s\n", input$label, input$assay, which_one, elapsed, digest::digest(result)))
} else if (mode == "compare") {
  input <- readRDS(args[[2L]]); old <- readRDS(sub("\\.rds$", ".old.result", args[[2L]])); new <- readRDS(sub("\\.rds$", ".new.result", args[[2L]]))
  difference <- if (identical(old, new)) 0 else max(abs(old - new), na.rm = TRUE)
  cat(sprintf("COMPARE %s %s %s samples %d features %d missing %.4f %s identical %s max_abs_difference %.3g\n", input$label, input$assay, input$state, nrow(input$Y), ncol(input$Y),
    mean(is.na(input$Y)), gsub(" ", "_", input$missing), identical(old, new), difference))
} else if (mode == "output") {
  label <- args[[2L]]; project <- args[[3L]]; prefix <- args[[4L]]; opened <- openCopy(project)
  manager <- opened$manager; omic <- manager$.context$omic_type; binding <- omicsDeclarationBindings()[[omic]]
  provenance <- manager$.context$snapshot$step_provenance[[paste0(prefix, "_ruv_corrected")]]
  parameters <- provenance$parameters %||% provenance$values
  result <- withWorkflowStateScope(manager, function() {
    object <- readWorkflowOmicsState(binding, manager, paste0(prefix, "_normalized"))
    stored <- readWorkflowOmicsState(binding, manager, paste0(prefix, "_ruv_corrected"))
    ruv <- suppressWarnings(suppressMessages(binding$optimise_ruv(object, parameters$ruv_mode, parameters[setdiff(names(parameters), "ruv_mode")])))
    best_k <- binding$extract_best_k(ruv); controls <- binding$extract_controls(ruv)
    corrected <- suppressWarnings(suppressMessages(ruvIII_C_Varying(theObject = object, ruv_grouping_variable = parameters$ruv_grouping_variable, ruv_number_k = best_k, ctrl = controls)))
    slot_name <- if (inherits(object, "MetaboliteAssayData")) "metabolite_data" else "lipid_data"
    a <- methods::slot(corrected, slot_name); b <- methods::slot(stored, slot_name)
    list(identical = mapply(function(x, y) identical(as.data.frame(x), as.data.frame(y)), a[names(b)], b), best_k = unlist(best_k),
      percentage = vapply(ruv, function(x) as.numeric(x$best_percentage %||% NA), numeric(1)))
  })
  line <- sprintf("OUTPUT %s %s corrected_assays_identical_to_stored %s k %s percentage %s", label, prefix, paste(names(result$identical), result$identical, collapse = ", "), paste(result$best_k, collapse = ","), paste(result$percentage, collapse = ","))
  cat(line, "\n"); write(line, args[[5L]], append = TRUE)
  closeDIAStateAdapter(manager); unlink(opened$work, recursive = TRUE)
}
