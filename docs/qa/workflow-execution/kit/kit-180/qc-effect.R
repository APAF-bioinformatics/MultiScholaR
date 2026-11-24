# KIT-180 script S1: the size of the change on the package's synthetic label-free and TMT fixtures (ticket section 5.3),
# and gate F1 (section 8).
#
# Usage, from the package root on a run-tests runner (the commit under test checked out, full git history):
#   Rscript docs/qa/workflow-execution/kit/kit-180/qc-effect.R <out dir> <arms>
#   <out dir>  default docs/qa/workflow-execution/kit/kit-180/qc-effect
#   <arms>     comma-separated subset of B,Q,B-ruv,Q-ruv (R0: B,B-ruv on the tests-only commit; R3: all four on the fix)
#
# Each arm runs the GUI-entry lanes prot_tmt, prot_lfq and prot_lfq_fragpipe (tests/testthat/helper-gui-entry.R,
# mod_proteomics_server under shiny::testServer) through the DA step; the summary steps after the DA are not run (they
# change no number measured here). No helper file is edited: the helpers are sourced into the global environment and
# four of them are wrapped there (F1):
#   - guiEntrySteps: arm B and B-ruv leave out every step whose id starts with `qc_`; arm Q and Q-ruv leave out only the
#     intensity refusal step; the RUV arms leave the normalisation tab's RUV mode at its rendered default; the
#     normalisation and DA steps are followed by a read of the owner (RUV outcome, DA table, QC provenance).
#   - guiEntryProtQcSettings: an empty list, so the four protein QC steps run at their rendered defaults.
#   - guiEntryStateDigests: the original digests, and for each protein state its distinct protein ids and rows.
# Written to <out dir> (integers, state names, digests and the commit; no data value):
#   qc-effect.tsv      lane, arm, commit, measure, value
#   state-digests.tsv  lane, arm, commit, step, state, digest (the A3 check: arms B and B-ruv base against head)
#   step-texts.tsv     lane, arm, commit, step, output, text (the QC steps' result texts that name their settings)
#   f1.tsv             the F1 checks (a worktree process resolves test_path() in its tree; guiEntryExpectLane() runs
#                      inside test_that() with its record and frozen file redirected)
# A lane that stops is recorded with its condition class; the other lanes still run. The script exits non-zero only
# when it cannot start.

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1L) args[[1L]] else "docs/qa/workflow-execution/kit/kit-180/qc-effect"
arms <- if (length(args) >= 2L) strsplit(args[[2L]], ",", fixed = TRUE)[[1L]] else c("B", "B-ruv")
stopifnot(all(arms %in% c("B", "Q", "B-ruv", "Q-ruv")))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir)
lanes <- c("prot_tmt", "prot_lfq", "prot_lfq_fragpipe")
commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
cat("KIT-180 S1 at", commit, "arms", paste(arms, collapse = ","), "\n")

suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = FALSE))
suppressMessages(testthat::source_test_helpers("tests/testthat", env = globalenv()))

rows <- list()
digest_rows <- list()
text_rows <- list()
put <- function(lane, arm, measure, value) {
    value <- if (is.null(value) || !length(value)) NA_character_ else paste(as.character(value), collapse = ";")
    rows[[length(rows) + 1L]] <<- data.frame(lane = lane, arm = arm, commit = commit, measure = measure, value = value,
        stringsAsFactors = FALSE)
}
writeTables <- function() {
    write <- function(parts, name) {
        if (!length(parts)) return(invisible())
        utils::write.table(do.call(rbind, parts), file.path(out_dir, name), sep = "\t", quote = FALSE, row.names = FALSE)
    }
    write(rows, "qc-effect.tsv")
    write(digest_rows, "state-digests.tsv")
    write(text_rows, "step-texts.tsv")
}

# ------------------------------------------------------------------------------------------------------------------
# The probe: what each run reads from the owner, reset before each lane and arm.
# ------------------------------------------------------------------------------------------------------------------
probe <- new.env(parent = emptyenv())
resetProbe <- function(arm) {
    probe$arm <- arm
    probe$counts <- list()
    probe$norm <- list()
    probe$da <- NULL
    probe$da_id_column <- NA_character_
    probe$id_table <- NULL
    probe$provenance <- list()
}
truthKey <- function(ids) {
    key <- rep(NA_character_, length(ids))
    hit <- grepl("P[0-9]{5}", ids)
    key[hit] <- regmatches(ids, regexpr("P[0-9]{5}", ids))
    key
}
stateCounts <- function(object) {
    if (!isS4(object) || !all(c("protein_quant_table", "protein_id_column") %in% methods::slotNames(object))) return(NULL)
    table <- object@protein_quant_table
    ids <- as.character(table[[object@protein_id_column]])
    list(proteins = length(unique(ids)), rows = nrow(table), ids = unique(ids))
}
readState <- function(manager, name) withWorkflowStateScope(manager, function() manager$getState(name))

original_steps <- guiEntrySteps
original_digests <- guiEntryStateDigests
guiEntryProtQcSettings <- function() list()
guiEntryStateDigests <- function(manager, cache) {
    out <- original_digests(manager)
    for (name in names(out)) if (is.null(probe$counts[[name]])) {
        probe$counts[[name]] <- tryCatch(stateCounts(readState(manager, name)),
            error = function(e) list(error = class(e)[[1L]]))
    }
    out
}
probeNorm <- function(api) {
    manager <- shiny::isolate(api$workflow_data$state_manager)
    head <- manager$current_state
    history <- guiEntryHistory(manager)
    probe$norm$head <- head
    if (!head %in% c("normalized", "ruv_corrected")) return(invisible())
    probe$norm$entering <- if (length(history) >= 2L) history[[length(history) - 1L]] else NA_character_
    probe$norm$read <- tryCatch(withWorkflowStateScope(manager, function() {
        record <- readDIAStateMetadata(manager, getDIAStateReference(manager, head))$data
        saved <- hydrateSavedWorkflowEntity(manager$.context$store, record$products$normalization_science)$value
        normalised <- stateCounts(materializeDIAStateRecord(manager$.context, saved$normalized)$data)
        final <- manager$getState(head)
        cutoffs <- final@args$removeRowsWithMissingValuesPercent
        controls <- saved$controls
        list(k = saved$best_k, percentage = saved$percentage,
            control_count = if (is.logical(controls)) sum(controls, na.rm = TRUE) else length(controls[!is.na(controls)]),
            ruv_applied = saved$ruv_result$ruv_applied, not_applied_reason = saved$ruv_result$ruv_not_applied_reason,
            proteins_normalised = normalised$proteins, proteins_final = stateCounts(final)$proteins,
            cutoff_groupwise = cutoffs$groupwise_percentage_cutoff, cutoff_max_groups = cutoffs$max_groups_percentage_cutoff,
            cutoff_percentile = cutoffs$proteins_intensity_cutoff_percentile)
    }), error = function(e) list(error = class(e)[[1L]]))
    invisible()
}
probeDa <- function(api) {
    manager <- shiny::isolate(api$workflow_data$state_manager)
    files <- list.files(api$inputs$paths$base_dir, pattern = "^da_proteins_.*_long_annot\\.tsv$", recursive = TRUE,
        full.names = TRUE)
    probe$da <- if (length(files)) as.data.frame(readr::read_tsv(files[[1L]], show_col_types = FALSE)) else NULL
    head_state <- tryCatch(readState(manager, manager$current_state), error = function(e) NULL)
    if (isS4(head_state)) probe$da_id_column <- head_state@protein_id_column
    history <- guiEntryHistory(manager)
    if ("protein_accession_cleaned" %in% history) {
        probe$id_table <- tryCatch(readState(manager, "protein_accession_cleaned")@protein_id_table, error = function(e) NULL)
    }
    for (state in intersect(c("protein_accession_cleaned", "protein_intensity_filtered", "duplicates_removed",
        "protein_replicate_filtered"), history)) {
        probe$provenance[[state]] <- manager$.context$snapshot$step_provenance[[state]]$values
    }
    invisible()
}
guiEntrySteps <- function(lane) {
    steps <- original_steps(lane)
    ids <- vapply(steps, `[[`, "", "id")
    steps <- steps[seq_len(which(ids == "da"))]
    ids <- vapply(steps, `[[`, "", "id")
    arm <- probe$arm
    drop <- if (arm %in% c("B", "B-ruv")) startsWith(ids, "qc_") else ids == "qc_protein_intensity_filter_refusal"
    steps <- steps[!drop]
    ids <- ids[!drop]
    norm <- which(ids == "norm")
    norm_run <- steps[[norm]]$run
    if (arm %in% c("B-ruv", "Q-ruv")) {
        # The tab's rendered defaults, RUV mode included (automatic); nothing else is set.
        norm_run <- function(api) {
            api$defaults(mod_prot_norm_ui, "normalization")
            api$set(list(`normalization-run_normalization` = 1L))
        }
    }
    steps[[norm]]$run <- function(api) {
        norm_run(api)
        probeNorm(api)
    }
    da <- which(ids == "da")
    da_run <- steps[[da]]$run
    steps[[da]]$run <- function(api) {
        da_run(api)
        probeDa(api)
    }
    steps
}

# ------------------------------------------------------------------------------------------------------------------
# One lane in one arm.
# ------------------------------------------------------------------------------------------------------------------
da_tables <- list()
runArm <- function(lane, arm) {
    resetProbe(arm)
    started <- Sys.time()
    record <- tryCatch(guiEntryRunLane(lane), error = function(e) e)
    put(lane, arm, "seconds", round(as.numeric(difftime(Sys.time(), started, units = "secs"))))
    if (inherits(record, "error")) {
        put(lane, arm, "status", paste("error", class(record)[[1L]]))
        cat(lane, arm, "stopped:", class(record)[[1L]], conditionMessage(record), "\n")
        return(invisible())
    }
    put(lane, arm, "status", "ok")
    put(lane, arm, "steps", paste(names(record$steps), collapse = ","))
    for (step in names(record$steps)) {
        capture <- record$steps[[step]]
        for (state in names(capture$states)) {
            digest_rows[[length(digest_rows) + 1L]] <<- data.frame(lane = lane, arm = arm, commit = commit, step = step,
                state = state, digest = as.character(capture$states[[state]]), stringsAsFactors = FALSE)
        }
        if (startsWith(step, "qc_protein_")) {
            texts <- Filter(function(value) is.character(value) && grepl("percentile|replicate|group|Proteins remaining",
                value, ignore.case = TRUE) && !startsWith(value, "<"), capture$outputs)
            for (name in names(texts)) {
                text_rows[[length(text_rows) + 1L]] <<- data.frame(lane = lane, arm = arm, commit = commit, step = step,
                    output = name, text = gsub("[\t\r\n]+", " ", texts[[name]]), stringsAsFactors = FALSE)
            }
        }
        if (step %in% c("norm", "correlation", "da")) {
            notes <- vapply(capture$notifications, function(note) paste0(note$type, ": ", substr(note$message, 1L, 160L)), "")
            notes <- grep("^(error|warning)", notes, value = TRUE)
            if (length(notes)) put(lane, arm, paste0("notifications:", step), gsub("[\t\r\n]+", " ", notes))
        }
    }
    put(lane, arm, "head_after_norm", probe$norm$head)
    final_history <- record$steps[[length(record$steps)]]$history
    put(lane, arm, "history", paste(final_history, collapse = ","))
    for (state in names(probe$counts)) {
        counts <- probe$counts[[state]]
        if (is.null(counts)) next
        if (!is.null(counts$error)) { put(lane, arm, paste0("count_error:", state), counts$error); next }
        put(lane, arm, paste0("proteins:", state), counts$proteins)
        put(lane, arm, paste0("rows:", state), counts$rows)
    }
    initial <- probe$counts$protein_s4_initial
    cleaned <- probe$counts$protein_accession_cleaned
    if (!is.null(initial$ids) && !is.null(cleaned$ids)) {
        put(lane, arm, "cleanup_ids_changed", sum(!cleaned$ids %in% initial$ids))
        put(lane, arm, "cleanup_rows_merged", initial$rows - cleaned$rows)
    }
    entering <- probe$norm$entering
    if (!is.null(entering)) {
        put(lane, arm, "norm_source_state", entering)
        put(lane, arm, "proteins_entering_normalisation", probe$counts[[entering]]$proteins)
    }
    read <- probe$norm$read
    if (!is.null(read$error)) {
        put(lane, arm, "norm_read_error", read$error)
        # Without the product's normalised checkpoint: the source's count less the normalisation state's.
        if (!is.null(entering) && probe$norm$head %in% names(probe$counts)) put(lane, arm, "removed_source_to_norm_state",
            probe$counts[[entering]]$proteins - probe$counts[[probe$norm$head]]$proteins)
    }
    if (!is.null(read) && is.null(read$error)) {
        put(lane, arm, "proteins_after_normalisation", read$proteins_normalised)
        put(lane, arm, "proteins_after_normalisation_state", read$proteins_final)
        if (arm %in% c("B-ruv", "Q-ruv")) {
            put(lane, arm, "ruv_k", read$k)
            put(lane, arm, "ruv_control_percentage", read$percentage)
            put(lane, arm, "ruv_control_count", read$control_count)
            put(lane, arm, "ruv_applied", read$ruv_applied)
            put(lane, arm, "ruv_not_applied_reason", read$not_applied_reason)
            put(lane, arm, "ruv_cleanup_removed", read$proteins_normalised - read$proteins_final)
            put(lane, arm, "ruv_cleanup_groupwise_cutoff", read$cutoff_groupwise)
            put(lane, arm, "ruv_cleanup_max_groups_cutoff", read$cutoff_max_groups)
            put(lane, arm, "ruv_cleanup_percentile", read$cutoff_percentile)
        }
    }
    if (arm %in% c("B-ruv", "Q-ruv") && !identical(probe$norm$head, "ruv_corrected")) {
        put(lane, arm, "ruv_outcome", paste("RUV refused on this fixture: head after the normalisation step",
            probe$norm$head %||% "<none>"))
    }
    put(lane, arm, "proteins_after_correlation", probe$counts$correlation_filtered$proteins)
    for (state in names(probe$provenance)) {
        values <- probe$provenance[[state]]
        for (name in names(values)) if (is.atomic(values[[name]]) && length(values[[name]]) <= 5L) {
            put(lane, arm, paste0("provenance:", state, ":", name), values[[name]])
        }
    }
    da <- probe$da
    if (is.null(da)) { put(lane, arm, "da_table", "absent"); return(invisible()) }
    if (!all(c("fdr_qvalue", "log2FC") %in% names(da))) {
        put(lane, arm, "da_columns", names(da))
        return(invisible())
    }
    id_column <- if (probe$da_id_column %in% names(da)) probe$da_id_column else names(da)[[1L]]
    ids <- as.character(da[[id_column]])
    key <- truthKey(ids)
    if (anyNA(key) && is.data.frame(probe$id_table)) {
        list_column <- paste0(probe$da_id_column, "_list")
        if (all(c(probe$da_id_column, list_column) %in% names(probe$id_table))) {
            lookup <- stats::setNames(as.character(probe$id_table[[list_column]]),
                as.character(probe$id_table[[probe$da_id_column]]))
            missing <- is.na(key)
            key[missing] <- truthKey(lookup[ids[missing]])
        }
    }
    put(lane, arm, "da_id_column", id_column)
    put(lane, arm, "da_ids_mapped", sum(is.na(key) | ids != key, na.rm = TRUE))
    put(lane, arm, "da_ids_unmapped", sum(is.na(key)))
    tested <- !is.na(da$fdr_qvalue) & !is.na(key)
    significant <- tested & da$fdr_qvalue < 0.05
    planted <- !is.na(key) & as.integer(substr(key, 2L, 6L)) <= 40L
    put(lane, arm, "da_tested", sum(tested))
    put(lane, arm, "da_significant", sum(significant))
    put(lane, arm, "da_planted_tested", sum(tested & planted))
    put(lane, arm, "da_planted_significant", sum(significant & planted))
    put(lane, arm, "da_null_tested", sum(tested & !planted))
    put(lane, arm, "da_null_significant", sum(significant & !planted))
    da_tables[[paste(lane, arm)]] <<- data.frame(key = key[tested], log2FC = da$log2FC[tested], stringsAsFactors = FALSE)
    invisible()
}

for (arm in arms) for (lane in lanes) {
    runArm(lane, arm)
    writeTables()
}

# The paired comparison: proteins tested in both arms of a pair, the largest absolute difference in log fold change.
for (pair in list(c("B", "Q"), c("B-ruv", "Q-ruv"))) for (lane in lanes) {
    a <- da_tables[[paste(lane, pair[[1L]])]]
    b <- da_tables[[paste(lane, pair[[2L]])]]
    if (is.null(a) || is.null(b)) next
    both <- merge(a, b, by = "key", suffixes = c(".a", ".b"))
    label <- paste(pair, collapse = "_vs_")
    put(lane, label, "tested_in_both", nrow(both))
    put(lane, label, "max_abs_log2fc_difference", if (nrow(both)) signif(max(abs(both$log2FC.a - both$log2FC.b)), 6L) else NA)
}
writeTables()

# ------------------------------------------------------------------------------------------------------------------
# F1: the mechanism S2 needs, on this commit. A git worktree of this commit; a fresh Rscript process whose working
# directory is the worktree's root loads it, resolves test_path() there, and runs guiEntryExpectLane("prot_tmt") inside
# test_that() with its record and frozen file redirected to absolute directories under out (removed afterwards; only
# f1.tsv is kept).
# ------------------------------------------------------------------------------------------------------------------
f1 <- list()
f1_put <- function(check, value) f1[[length(f1) + 1L]] <<- data.frame(check = check, value = as.character(value),
    stringsAsFactors = FALSE)
tree <- file.path(tempdir(), "kit180-f1-tree")
f1_out <- file.path(out_dir, "f1-scratch")
dir.create(file.path(f1_out, "records"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(f1_out, "oracles"), recursive = TRUE, showWarnings = FALSE)
added <- system2("git", c("worktree", "add", "--detach", tree, commit), stdout = TRUE, stderr = TRUE)
f1_put("worktree_added", dir.exists(file.path(tree, "tests", "testthat")))
child <- file.path(tempdir(), "kit180-f1-child.R")
writeLines(c(
    "setwd(commandArgs(trailingOnly = TRUE)[[1L]])",
    "suppressPackageStartupMessages(pkgload::load_all('.', quiet = TRUE, export_all = TRUE, helpers = FALSE))",
    "suppressMessages(testthat::source_test_helpers('tests/testthat', env = globalenv()))",
    "cat('F1 test_path', normalizePath(testthat::test_path('..', 'testdata', 'gui-entry')), '\\n')",
    "ok <- testthat::test_that('KIT-180 F1 prot_tmt inside test_that', { guiEntryExpectLane('prot_tmt') })",
    "cat('F1 test_that', isTRUE(ok), '\\n')",
    "quit(status = if (isTRUE(ok)) 0L else 1L)"), child)
console <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), c(shQuote(child), shQuote(tree)),
    stdout = TRUE, stderr = TRUE, env = c(
        paste0("MULTISCHOLAR_GUI_ENTRY_RECORD_DIR=", shQuote(file.path(f1_out, "records"))),
        paste0("MULTISCHOLAR_GUI_ENTRY_ORACLE_DIR=", shQuote(file.path(f1_out, "oracles"))),
        "MULTISCHOLAR_GUI_ENTRY_FREEZE=1")))
status <- attr(console, "status") %||% 0L
path_line <- grep("^F1 test_path", console, value = TRUE)
f1_put("child_exit_status", status)
f1_put("test_path_inside_worktree", length(path_line) == 1L &&
    startsWith(sub("^F1 test_path ", "", trimws(path_line)), normalizePath(tree)))
f1_put("test_that_returned_true", any(grepl("^F1 test_that TRUE", console)))
f1_put("record_written", file.exists(file.path(f1_out, "records", "prot_tmt.json")))
f1_put("frozen_file_written", file.exists(file.path(f1_out, "oracles", "prot_tmt.json")))
f1_put("checkout_oracle_unchanged", identical(unname(tools::md5sum(file.path(tree, "tests", "testdata", "gui-entry", "prot_tmt.json"))),
    unname(tools::md5sum(file.path("tests", "testdata", "gui-entry", "prot_tmt.json")))))
cat(utils::tail(console, 15L), sep = "\n")
utils::write.table(do.call(rbind, f1), file.path(out_dir, "f1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
unlink(f1_out, recursive = TRUE)
invisible(system2("git", c("worktree", "remove", "--force", tree), stdout = TRUE, stderr = TRUE))
cat("KIT-180 S1 done\n")
