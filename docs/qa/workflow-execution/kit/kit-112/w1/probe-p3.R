# KIT-112 W1 probe P3 (scratch, not frozen): an input among the binding's `driver_settings$da` values (formula, q_value,
# treat_lfc) that the controller's validation accepts and the kernel refuses in process, on lfq101_fixture's fixture.
#   Rscript docs/qa/workflow-execution/kit/kit-112/w1/probe-p3.R <package source> <copy of 111l's .dia-native-* store root>
# The store is a copy of `111l`'s (kit_111, the LFQ fixture, synthetic), so its DA source is the journey's own; the copy
# is reopened as a session reopens it and the DA runs in process on it with each candidate.

args <- commandArgs(TRUE)
tree <- normalizePath(args[[1L]])
root <- normalizePath(args[[2L]])
setwd(tree)
suppressMessages(pkgload::load_all(tree, quiet = TRUE, export_all = TRUE, helpers = TRUE))
say <- function(...) cat(sprintf(...), "\n", sep = "")

identity <- jsonlite::read_json(file.path(root, ".workflow-entities", "store.json"))
manager <- newDIAStateAdapter(newWorkflowEntityStore(root, identity$project_id, identity$workflow_id),
    WorkflowState$new(audit_enabled = FALSE), reopen = TRUE)
snapshot <- manager$.context$snapshot
manifest <- resolveWorkflowStateValue(manager$.context, snapshot$da_results)
source <- manifest$source
design <- withDIAStateScope(manager, function() readDIAExecutionComponents(manager, source, fields = "design_matrix")@design_matrix)
say("P3 store %s/%s; head state %s; workflow_type %s; DA source state %s; recorded parameters: %s", identity$project_id,
    identity$workflow_id, manager$current_state, manager$workflow_type %||% "NULL",
    names(snapshot$states)[vapply(snapshot$states, function(r) identical(r$entity_id, source$entity_id), NA)],
    paste(names(manifest$parameters), unlist(lapply(manifest$parameters, format)), sep = "=", collapse = ", "))
say("P3 design: %d samples; columns %s", nrow(design), paste(names(design), collapse = ", "))
print(design)
say("P3 contrasts: %s", paste(manifest$contrasts$full_format, collapse = "; "))

candidates <- list(
    list(label = "the binding's own (control)", formula = "~ 0 + group", q = 0.05, lfc = 0),
    list(label = "q_value 1e-300 (082e's)", formula = "~ 0 + group", q = 1e-300, lfc = 0),
    list(label = "q_value 0", formula = "~ 0 + group", q = 0, lfc = 0),
    list(label = "q_value 1", formula = "~ 0 + group", q = 1, lfc = 0),
    list(label = "q_value -1", formula = "~ 0 + group", q = -1, lfc = 0),
    list(label = "treat_lfc 100", formula = "~ 0 + group", q = 0.05, lfc = 100),
    list(label = "treat_lfc -1", formula = "~ 0 + group", q = 0.05, lfc = -1),
    list(label = "formula with an intercept", formula = "~ group", q = 0.05, lfc = 0))
# Formulas the design can express: one more term per further design column (a saturated model leaves no residual
# degrees of freedom).
for (column in setdiff(names(design), c("group", "Run", "sample", "Sample"))) {
    candidates[[length(candidates) + 1L]] <- list(label = paste("formula + design column", column),
        formula = paste("~ 0 + group +", column), q = 0.05, lfc = 0)
}
for (column in intersect(names(design), c("Run", "sample", "Sample"))) {
    candidates[[length(candidates) + 1L]] <- list(label = paste("formula + the sample column", column),
        formula = paste("~ 0 + group +", column), q = 0.05, lfc = 0)
}

for (candidate in candidates) {
    prototype <- withDIAStateScope(manager, function() readDIAExecutionComponents(manager, source, fields = "design_matrix"))
    validation <- tryCatch({
        validateProtDaModelAndContrasts(prototype, candidate$formula, manifest$contrasts)
        "accepted"
    }, error = function(e) paste("refused:", conditionMessage(e)))
    kernel <- if (startsWith(validation, "accepted")) {
        workflow <- list2env(list(state_manager = manager, contrasts_tbl = manifest$contrasts, tab_status = list()))
        da <- new.env(parent = emptyenv())
        da$current_s4_object <- getDIAStateReference(manager)
        outcome <- NULL
        utils::capture.output(outcome <- tryCatch({
            value <- suppressWarnings(suppressMessages(runDIAReferenceDaAnalysis(list(formula_string = candidate$formula,
                da_q_val_thresh = candidate$q, treat_lfc_cutoff = candidate$lfc), da, workflow, list(), manifest$contrasts,
                NULL, publish = FALSE)))
            sprintf("completed (contrast %s, long table %d rows)", paste(names(value$results), collapse = ","),
                value$results[[1L]]$rows)
        }, error = function(e) paste("FAILED:", conditionMessage(e))))
        outcome
    } else "not run"
    say("P3 %-40s formula '%s', q %s, lfc %s: validation %s; kernel %s", candidate$label, candidate$formula, format(candidate$q),
        format(candidate$lfc), validation, kernel)
}
closeDIAStateAdapter(manager)
