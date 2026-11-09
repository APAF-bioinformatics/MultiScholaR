# KIT-112 W1, the lead's question (W0): what in the DA click path hydrates the protein matrix or the design in the
# controller, and would D3's request half hydrate anything? A fixture-level measurement on the start commit's code:
#   Rscript docs/qa/workflow-execution/kit/kit-112/w1/request-half-hydration.R <package source> [proteins] [samples]
# A synthetic DIA protein state (normal intensities, two groups, technical-replicate pairs) is accepted at
# `correlation_filtered`; then each statement of `runDIAReferenceDaAnalysis` before the contrast loop (the request half
# of D3) runs alone, then the kernel's projected read (`runProjected` -> `readWorkflowExecutionComponents`), each with the
# R heap's peak over the statement (`gc(reset = TRUE)` before, `gc()`'s "max used" after) and its result's size.
# No large run: the default is 4,000 proteins x 300 samples (a 9.2 MiB matrix).

args <- commandArgs(TRUE)
tree <- normalizePath(args[[1L]])
proteins <- if (length(args) >= 2L) as.integer(args[[2L]]) else 4000L
samples <- if (length(args) >= 3L) as.integer(args[[3L]]) else 300L
setwd(tree)
suppressMessages(pkgload::load_all(tree, quiet = TRUE, export_all = TRUE, helpers = TRUE))

mib <- function(bytes) sprintf("%.2f MiB", bytes / 1024^2)
#' The R heap's peak (Ncells + Vcells, bytes) while `code` runs, above the heap in use before it, and the value's size.
measure <- function(label, code) {
    # Each statement runs twice: the first (cold) call also loads namespaces and caches; the second (warm) is the
    # statement's own cost, as on the DA click after the session's earlier steps.
    expression <- substitute(code)
    frame <- parent.frame()
    peaks <- character()
    for (pass in c("cold", "warm")) {
        invisible(gc(full = TRUE))
        before <- gc(reset = TRUE)
        used_before <- sum(before[, 2L]) * 1024^2
        value <- eval(expression, frame)
        after <- gc()
        peaks[[pass]] <- mib(sum(after[, ncol(after)]) * 1024^2 - used_before)
    }
    cat(sprintf("%-92s heap peak cold +%s, warm +%s; value %s\n", label, peaks[["cold"]], peaks[["warm"]],
        mib(as.numeric(utils::object.size(value)))))
    invisible(value)
}

set.seed(1L)
runs <- sprintf("run_%04d", seq_len(samples))
groups <- rep(c("control", "treated"), each = samples / 2L)
design <- data.frame(Run = runs, group = factor(groups, levels = c("control", "treated")),
    tech_rep_group = paste(groups, rep(seq_len(samples / 4L), each = 2L, length.out = samples), sep = "_"),
    stringsAsFactors = FALSE, row.names = runs)
ids <- sprintf("PG_%06d", seq_len(proteins))
quant <- as.data.frame(c(list(Protein.Group = ids), stats::setNames(lapply(seq_len(samples), function(i)
    stats::rnorm(proteins, 23, 2)), runs)), check.names = FALSE, stringsAsFactors = FALSE)
object <- methods::new("ProteinQuantitativeData", protein_quant_table = quant, protein_id_column = "Protein.Group",
    protein_id_table = data.frame(Protein.Group = ids, gene_names = ids, stringsAsFactors = FALSE), design_matrix = design,
    sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group",
    args = list(globalParameters = list(workflow_type = "DIA")))
cat(sprintf("synthetic state: %d proteins x %d samples; protein table %s in R\n", proteins, samples,
    mib(as.numeric(utils::object.size(quant)))))

root <- tempfile("k112-hydration-")
dir.create(root)
manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "k112-hydration", "DIA"), WorkflowState$new(audit_enabled = FALSE))
manager$.context$reference_states <- TRUE
manager$setWorkflowType("DIA")
manager$saveState("correlation_filtered", object, list(), "fixture")
rm(object, quant)
contrasts <- normaliseProtDaContrastsTable(data.frame(contrasts = "treated.vs.control",
    full_format = "treated.vs.control=grouptreated-groupcontrol", friendly_names = "treated.vs.control"))
settings <- list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0)

withDIAStateScope(manager, function() {
    source <- getDIAStateReference(manager)
    cat("-- the request half (runDIAReferenceDaAnalysis before the contrast loop, utils_workflow_dia_da.R:345-389) --\n")
    measure(":345 assertDIAExportReady(manager, allow_action = 'da')", assertDIAExportReady(manager, allow_action = "da"))
    spec <- measure(":353 readDIAStateMetadata(manager, source)$data (the state record)", readDIAStateMetadata(manager, source)$data)
    prototype <- measure(":354 readDIAExecutionComponents(manager, source, fields = 'design_matrix')",
        readDIAExecutionComponents(manager, source, fields = "design_matrix"))
    cat(sprintf("    the prototype's protein_quant_table: %d rows x %d columns (a zero-row schema)\n",
        nrow(prototype@protein_quant_table), ncol(prototype@protein_quant_table)))
    measure(":355 validateProtDaModelAndContrasts(prototype, formula, contrasts)",
        validateProtDaModelAndContrasts(prototype, settings$formula_string, contrasts))
    prototype <- NULL
    state_record <- measure(":367 readDIAStateMetadata(manager) (the head record)", readDIAStateMetadata(manager))
    checkpoint <- state_record$data
    checkpoint$slots$args$daAnalysisUI <- c(settings[c("formula_string")], list(q_value_threshold = 0.05,
        log_fold_change_cutoff = 0))
    measure(":383 commitDIAReferenceState(manager, current_state, checkpoint, config, description)",
        commitDIAReferenceState(manager, manager$current_state, checkpoint, state_record$config, state_record$description))
    source <- getDIAStateReference(manager)
    spec <- readDIAStateMetadata(manager, source)$data
    cat("-- the compute half's first statements (runSavedOperation, utils_workflow_operations.R:151-174) --\n")
    operation <- diaContrastOperation("iq")
    parameters <- list(formula = settings$formula_string, qvalue = 0.05, lfc = 0, weights = NA,
        contrast = contrasts$full_format[[1L]], annotations = spec$components$protein_id_table$entity_id,
        name = contrasts$contrasts[[1L]], comparison = contrasts$friendly_names[[1L]])
    measure("utils_workflow_operations.R:151 workflowComponentScientificIdentity (the identity)",
        workflowComponentScientificIdentity(manager, source, operation$name, parameters[operation$parameters]))
    measure("utils_workflow_operations.R:174 runProjected's read: readWorkflowExecutionComponents(3 fields)",
        readWorkflowExecutionComponents(manager, source, fields = operation$fields))
    invisible(NULL)
})
closeDIAStateAdapter(manager)
unlink(root, recursive = TRUE)
