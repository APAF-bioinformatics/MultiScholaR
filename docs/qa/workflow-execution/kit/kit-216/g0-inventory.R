# KIT-216 gate G0 (r1.2, section 3.0): an inventory, no assertion. For every proteomics Summary call made by the
# R2-scope test files below, print the owner's class, whether it is a reference-state `DIAStateAdapter`, its head, its
# history, whether a normalisation state is in the store, whether an export exists in the session, and what the
# ticket's design predicts there (I10 by name, SD-1, SD-2, SD-4; the final object, P06). The base package is not
# changed: the five Summary functions are wrapped in this process only, each wrapper records one row and then runs the
# function's own body. The probe's test counts are not regression evidence (the wrappers drop the functions' srcref).
#
# Usage (from the package root, on a runner): Rscript docs/qa/workflow-execution/kit/kit-216/g0-inventory.R
#   with no argument it runs each file below in its own R process (15 minutes each at most) and writes
#   docs/qa/workflow-execution/kit/kit-216/g0-out/{calls.tsv,files.tsv,<file>.log}.

kit216G0Dir <- "docs/qa/workflow-execution/kit/kit-216"
kit216G0Out <- file.path(kit216G0Dir, "g0-out")
kit216G0Columns <- c("file", "test_line", "call_line", "fn", "omic", "manager_class", "reference_states", "native",
    "head", "history_tail", "normalisation_in_store", "export", "export_policy", "export_kind", "export_guard",
    "final_now", "final_after_p06", "i10_prediction", "save_prediction")
kit216G0Files <- c(
    "test-workflow-integration-input.R", "test-workflow-memory-mechanism-member.R", "test-workflow-native-session.R",
    "test-workflow-analysis-record.R", "test-kit206-prot-norm-export-guard.R", "test-kit207-summary-head-guard.R",
    "test-kit172-downstream.R", "test-kit129-run-path.R", "test-workflow-omics-application-reads.R",
    "test-omics-summary-core.R", "test-omic-summary-module.R", "test-prot-12-summary-module-contracts.R",
    "test-prot-12a-summary-module-characterization.R",
    # Outside the R2 grep, a proteomics Summary module test: recorded to show whether the R2 scope misses a mover.
    "test-module-ci-prot-summary.R")
kit216G0Targets <- c("saveWorkflowSummaryParameters", "renderWorkflowSummaryReport", "writeWorkflowSummarySessionState",
    "protSummaryIntegrationWrite", "resolveProtSummaryFinalS4State")

arguments <- commandArgs(trailingOnly = TRUE)

if (!length(arguments)) {
    dir.create(kit216G0Out, recursive = TRUE, showWarnings = FALSE)
    writeLines(paste(kit216G0Columns, collapse = "\t"), file.path(kit216G0Out, "calls.tsv"))
    self <- file.path(kit216G0Dir, "g0-inventory.R")
    status <- data.frame(file = character(), exit = integer(), seconds = numeric(), stringsAsFactors = FALSE)
    for (file in kit216G0Files) {
        started <- Sys.time()
        log <- file.path(kit216G0Out, paste0(sub("[.]R$", "", file), ".log"))
        exit <- system2("timeout", c("--kill-after=30", "900", "Rscript", "--vanilla", self, file), stdout = log,
            stderr = log)
        status[nrow(status) + 1L, ] <- list(file, as.integer(exit),
            round(as.numeric(difftime(Sys.time(), started, units = "secs")), 1))
        cat(sprintf("KIT-216 G0 %s exit %d after %.0f s\n", file, exit, status$seconds[nrow(status)]))
    }
    utils::write.table(status, file.path(kit216G0Out, "files.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    calls <- utils::read.delim(file.path(kit216G0Out, "calls.tsv"), colClasses = "character")
    cat(sprintf("KIT-216 G0: %d proteomics Summary calls recorded\n", nrow(calls)))
    quit(save = "no", status = 0L)
}

file <- arguments[[1L]]
calls_file <- normalizePath(file.path(kit216G0Out, "calls.tsv"))
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = TRUE))
ns <- asNamespace("MultiScholaR")
kit216G0Four <- c("ruv_covariates", "normalized", "ruv_corrected", "correlation_filtered")
kit216G0Text <- function(x) {
    x <- tryCatch(paste(as.character(unlist(x)), collapse = ","), error = function(e) "<error>")
    if (!length(x) || !nzchar(x)) "<none>" else gsub("[\t\r\n]+", " ", x)
}

#' The test-file locations on the call stack: the outermost (the `test_that` call) and the innermost.
kit216G0Where <- function() {
    locations <- character()
    for (call in sys.calls()) {
        ref <- attr(call, "srcref")
        name <- if (!is.null(ref)) attr(ref, "srcfile")$filename else NULL
        if (is.character(name) && length(name) == 1L && grepl("^test-", basename(name))) {
            locations <- c(locations, sprintf("%s:%d", basename(name), ref[[1L]]))
        }
    }
    if (!length(locations)) return(c(NA_character_, NA_character_))
    c(locations[[1L]], locations[[length(locations)]])
}

kit216G0Probe <- function(fn, frame) {
    tryCatch({
        omic <- if (exists("omicType", envir = frame, inherits = FALSE)) get("omicType", envir = frame) else "proteomics"
        if (!identical(omic, "proteomics")) return(invisible(NULL))
        workflow <- get0("workflowData", envir = frame, inherits = FALSE)
        manager <- tryCatch(shiny::isolate(workflow$state_manager), error = function(e) NULL)
        where <- kit216G0Where()
        reference <- tryCatch(isTRUE(ns$usesDIAReferenceStates(manager)), error = function(e) NA)
        native <- tryCatch(isTRUE(ns$isDIANativeSession(manager)), error = function(e) NA)
        head <- tryCatch(as.character(manager$current_state)[1L], error = function(e) NA_character_)
        history <- tryCatch(ns$diaStateHistory(manager), error = function(e) character())
        context <- tryCatch(if (inherits(manager, "DIAStateAdapter")) manager$.context, error = function(e) NULL)
        stored <- names(context$snapshot$states)
        normalised <- if (is.null(context)) NA else any(c("normalized", "ruv_corrected", "ruv_covariates") %in% stored)
        owned <- context$owned_export
        export <- !is.null(owned$metadata_reference)
        metadata <- if (export) tryCatch(ns$resolveDIAStateValue(context, owned$metadata_reference),
            error = function(e) list(output_policy = "<unreadable>"))
        guard <- if (!export) "<no export>" else if (!identical(fn, "saveWorkflowSummaryParameters")) "<not judged>" else
            if (identical(owned$generation, context$head$entity_id)) "same head (refresh returns)" else
            tryCatch({
                ns$assertProtNormExportHeadCurrent(manager, owned$extra_references$correlation_filtered)
                "pass"
            }, error = function(e) class(e)[[1L]])
        listed <- c("ruv_covariates", "correlation_filtered", "ruv_corrected", "protein_replicate_filtered", "imputed")
        now <- utils::head(listed[listed %in% history], 1L)
        now <- if (length(now)) now else "<none>"
        i10 <- if (!isTRUE(reference)) "pass (not a reference-state owner)" else if (isTRUE(head %in% kit216G0Four))
            "pass" else if (isFALSE(normalised)) "refuse summary_none (SD-1)" else "refuse T3"
        after <- if (isTRUE(reference) && isTRUE(head %in% kit216G0Four)) head else if (startsWith(i10, "refuse"))
            "<refused first>" else now
        save <- if (!identical(fn, "saveWorkflowSummaryParameters") || !startsWith(i10, "pass")) "-" else
            if (!isTRUE(native)) "as today (not native)" else if (!export) "SD-2: no refresh, input written, T2" else
            if (!identical(metadata$output_policy, "native_first_v1")) "SD-2 (not a normalisation export), T2" else
            if (identical(guard, "prot_norm_correlation_input_stale")) "SD-4: no refresh, input written, T2b" else
            "refresh as today"
        row <- c(file, where[[1L]], where[[2L]], fn, omic, kit216G0Text(class(manager)), format(reference),
            format(native), kit216G0Text(head), kit216G0Text(utils::tail(history, 4L)), format(normalised), format(export),
            kit216G0Text(metadata$output_policy), kit216G0Text(metadata$checkpoint_kind), guard, now, after, i10, save)
        cat(paste(row, collapse = "\t"), "\n", file = calls_file, append = TRUE, sep = "")
    }, error = function(e) cat("KIT-216 G0 probe error:", conditionMessage(e), "\n"))
    invisible(NULL)
}

package_env <- tryCatch(as.environment("package:MultiScholaR"), error = function(e) NULL)
for (fn in kit216G0Targets) {
    wrapped <- get(fn, envir = ns)
    body(wrapped) <- bquote({
        kit216G0Probe(.(fn), environment())
        .(body(wrapped))
    })
    attr(wrapped, "srcref") <- NULL
    for (env in list(ns, package_env)) {
        if (is.null(env) || !exists(fn, envir = env, inherits = FALSE)) next
        locked <- bindingIsLocked(fn, env)
        if (locked) unlockBinding(fn, env)
        assign(fn, wrapped, envir = env)
        if (locked) lockBinding(fn, env)
    }
}

reporter <- testthat::ListReporter$new()
outcome <- tryCatch({
    testthat::test_file(file.path("tests", "testthat", file), reporter = reporter, stop_on_failure = FALSE)
    "completed"
}, error = function(e) conditionMessage(e))
results <- as.data.frame(reporter$get_results())
cat(sprintf("KIT-216 G0 file %s: %s; tests %d, failed %d, errors %d, skipped %d\n", file, outcome, nrow(results),
    if (nrow(results)) sum(results$failed) else 0L, if (nrow(results)) sum(results$error) else 0L,
    if (nrow(results)) sum(results$skipped) else 0L))
