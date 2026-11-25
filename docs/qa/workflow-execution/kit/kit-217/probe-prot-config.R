# KIT-217 probe P (the coordinator's ruling on R3 run 45587576015, item 1): why the proteomics GUI-entry lanes move from
# the `design` step on. Observation only; no package code. Runs as a run-tests script on the head commit (kit-217 after
# the merge of staging 16112917); the base side is a worktree of the base commit.
#
# Usage, from the package root (full git history):
#   Rscript docs/qa/workflow-execution/kit/kit-217/probe-prot-config.R <base commit> [<out dir>]
#
# Per side a fresh Rscript (`--side`) runs `guiEntryRunLane()` (tests/testthat/helper-gui-entry.R, loaded with
# pkgload::load_all(helpers = TRUE)) for the lanes prot_tmt and prot_lfq. The helper's `guiEntryStateDigests` is wrapped
# (the original is still called and returned, so the lane runs exactly as in the test) to write, after every step and
# for every committed state: (a) the state's digest as the test computes it (`guiEntryDigest`); (b) the digest of each
# S4 slot separately (the quantitative tables, the design, `args`, ...), and the state record's stored configuration
# (`readWorkflowStateMetadata(...)$config`) as one digest; (c) every leaf of `args` and of the stored configuration
# (path and canonical value). Reported steps: `design` and every later step (the comparison reads them all).
#
# PREDICTIONS (written before any value is read; from the inference in outcome.md):
#   P1  every slot other than `args`, and the stored configuration list, is identical on base and head at every step
#       and state of both lanes (so no proteomics table, design or result value differs);
#   P2  `args` differs at every state from `protein_s4_initial` on, and the leaves that differ are exactly the two keys
#       added to inst/config/config.ini: metaboliteIntensityFiltering/metabolites_min_values_per_group (head "2", base
#       absent) and metaboliteIntensityFiltering/metabolites_min_groups (head "1", base absent); no other leaf differs;
#   P3  the whole-state digest differs exactly where `args` differs.
# Written: <out>/<side>-slots.tsv, <out>/<side>-leaves.tsv, <out>/compare.tsv (check, lane, step, state, item, equal,
# predicted, holds) and <out>/leaf-diff.tsv. Exit status non-zero when a side fails or a prediction does not hold.

args <- commandArgs(trailingOnly = TRUE)
LANES <- c("prot_tmt", "prot_lfq")

childMain <- function(side, out) {
    # The test helpers are sourced into their own environment by load_all(helpers = TRUE), not into the namespace:
    # take the environment guiEntryRunLane encloses (where it looks up guiEntryStateDigests). (Probe run 45597931823.)
    ns <- environment(get("guiEntryRunLane", mode = "function"))
    sha <- function(x) digest::digest(x, algo = "sha256")
    slots_rows <- list()
    leaf_rows <- list()
    current <- new.env(parent = emptyenv())
    current$lane <- NA_character_
    current$i <- 0L
    leaves <- function(x, path = character()) {
        if (is.list(x) && !is.data.frame(x) && length(x)) {
            nm <- names(x) %||% rep("", length(x))
            nm[!nzchar(nm)] <- paste0("[", which(!nzchar(nm)), "]")
            return(do.call(rbind, lapply(seq_along(x), function(k) leaves(x[[k]], c(path, nm[[k]])))))
        }
        value <- if (is.atomic(x) && length(x) <= 5L && !is.data.frame(x)) paste(format(x), collapse = ",") else
            paste0("<", class(x)[[1L]], " ", sha(get("guiEntryCanonical", ns)(x)), ">")
        data.frame(path = paste(path, collapse = "/"), value = value, stringsAsFactors = FALSE)
    }
    original <- get("guiEntryStateDigests", ns)
    wrapped <- function(manager, cache) {
        current$i <- current$i + 1L
        step <- current$steps[[current$i]]
        history <- get("guiEntryHistory", ns)(manager)
        for (name in unique(history)) {
            res <- tryCatch(withWorkflowStateScope(manager, function() {
                object <- manager$getState(name)
                config <- tryCatch(readWorkflowStateMetadata(manager, getWorkflowStateReference(manager, name))$config,
                    error = function(e) NULL)
                list(object = object, config = config)
            }), error = function(e) NULL)
            if (is.null(res) || !isS4(res$object)) next
            object <- res$object
            base_row <- list(side = side, lane = current$lane, step = step, state = name)
            slots_rows[[length(slots_rows) + 1L]] <<- data.frame(base_row, item = "<state>",
                digest = get("guiEntryDigest", ns)(object), stringsAsFactors = FALSE)
            for (s in methods::slotNames(object)) {
                slots_rows[[length(slots_rows) + 1L]] <<- data.frame(base_row, item = paste0("slot:", s),
                    digest = sha(get("guiEntryCanonical", ns)(methods::slot(object, s))), stringsAsFactors = FALSE)
            }
            slots_rows[[length(slots_rows) + 1L]] <<- data.frame(base_row, item = "config",
                digest = sha(get("guiEntryCanonical", ns)(res$config)), stringsAsFactors = FALSE)
            for (what in c("args", "config")) {
                value <- if (identical(what, "args")) methods::slot(object, "args") else res$config
                if (!is.list(value) || !length(value)) next
                lv <- leaves(value)
                if (is.null(lv) || !nrow(lv)) next
                leaf_rows[[length(leaf_rows) + 1L]] <<- data.frame(base_row, what = what, lv, stringsAsFactors = FALSE)
            }
        }
        original(manager, cache)
    }
    if (bindingIsLocked("guiEntryStateDigests", ns)) unlockBinding("guiEntryStateDigests", ns)
    assign("guiEntryStateDigests", wrapped, envir = ns)
    failed <- FALSE
    for (lane in LANES) {
        cat("KIT-217 P", side, "lane", lane, "\n")
        current$lane <- lane
        current$i <- 0L
        current$steps <- vapply(get("guiEntrySteps", ns)(get("guiEntryLane", ns)(lane)), `[[`, "", "id")
        ok <- tryCatch({ get("guiEntryRunLane", ns)(lane); TRUE }, error = function(e) {
            cat("KIT-217 P", side, "lane", lane, "ERROR", class(e)[[1L]], conditionMessage(e), "\n")
            FALSE
        })
        if (!ok) failed <- TRUE
    }
    utils::write.table(do.call(rbind, slots_rows), file.path(out, paste0(side, "-slots.tsv")), sep = "\t", row.names = FALSE,
        quote = FALSE)
    leaf_table <- do.call(rbind, leaf_rows)
    leaf_table$value <- gsub("[\t\n\r]", " ", leaf_table$value)
    utils::write.table(leaf_table, file.path(out, paste0(side, "-leaves.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
    quit(status = if (failed) 1L else 0L)
}

if (length(args) >= 1L && identical(args[[1L]], "--side")) {
    setwd(args[[2L]])
    suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = TRUE))
    environment(childMain) <- new.env(parent = asNamespace("MultiScholaR"))
    childMain(args[[3L]], args[[4L]])
}

# ---- Parent --------------------------------------------------------------------------------------------------------------
stopifnot(length(args) >= 1L, grepl("^[0-9a-f]{40}$", args[[1L]]))
base_commit <- args[[1L]]
out <- if (length(args) >= 2L) args[[2L]] else "docs/qa/workflow-execution/kit/kit-217/p"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
head_root <- normalizePath(".")
cat("KIT-217 P: head", system2("git", c("rev-parse", "HEAD"), stdout = TRUE), "base", base_commit, "\n")
script <- normalizePath(file.path("docs", "qa", "workflow-execution", "kit", "kit-217", "probe-prot-config.R"))
rscript <- file.path(R.home("bin"), "Rscript")
runSide <- function(tree, side) {
    console <- suppressWarnings(system2(rscript, c(shQuote(script), "--side", shQuote(tree), side, shQuote(out)),
        stdout = TRUE, stderr = TRUE))
    writeLines(console, file.path(out, paste0(side, "-console.txt")))
    cat(grep("^KIT-217 P", console, value = TRUE), sep = "\n")
    cat(utils::tail(console, 10L), sep = "\n")
    attr(console, "status") %||% 0L
}
base_tree <- file.path(tempdir(), "kit217-p-base")
invisible(system2("git", c("worktree", "add", "--detach", shQuote(base_tree), base_commit), stdout = TRUE, stderr = TRUE))
status <- c(base = runSide(base_tree, "base"), head = runSide(head_root, "head"))

read <- function(name) {
    path <- file.path(out, name)
    if (!file.exists(path)) return(NULL)
    utils::read.delim(path, colClasses = "character", quote = "", na.strings = character())
}
sb <- read("base-slots.tsv"); sh <- read("head-slots.tsv")
lb <- read("base-leaves.tsv"); lh <- read("head-leaves.tsv")
rows <- list()
if (!is.null(sb) && !is.null(sh)) {
    key <- function(x) paste(x$lane, x$step, x$state, x$item, sep = "|")
    keys <- union(key(sb), key(sh))
    a <- sb[match(keys, key(sb)), ]; b <- sh[match(keys, key(sh)), ]
    parts <- strsplit(keys, "|", fixed = TRUE)
    item <- vapply(parts, `[`, "", 4L)
    equal <- !is.na(a$digest) & !is.na(b$digest) & a$digest == b$digest
    # args moves wherever a state carries the configuration; the whole-state digest with it (P3)
    args_equal <- stats::setNames(equal[item == "slot:args"], vapply(parts[item == "slot:args"], function(p) paste(p[1:3], collapse = "|"), ""))
    state_key <- vapply(parts, function(p) paste(p[1:3], collapse = "|"), "")
    predicted <- ifelse(item == "slot:args", "differs", ifelse(item == "<state>",
        ifelse(args_equal[state_key] %in% FALSE, "differs", "equal"), "equal"))
    rows[[1L]] <- data.frame(check = ifelse(item == "<state>", "P3", ifelse(item == "slot:args", "P2_args", "P1")),
        lane = vapply(parts, `[`, "", 1L), step = vapply(parts, `[`, "", 2L), state = vapply(parts, `[`, "", 3L), item = item,
        equal = equal, predicted = predicted, holds = (predicted == "equal") == equal, stringsAsFactors = FALSE)
}
leaf_diff <- NULL
if (!is.null(lb) && !is.null(lh)) {
    k <- function(x) paste(x$lane, x$step, x$state, x$what, x$path, sep = "|")
    keys <- union(k(lb), k(lh))
    a <- lb[match(keys, k(lb)), ]; b <- lh[match(keys, k(lh)), ]
    differs <- is.na(a$value) | is.na(b$value) | a$value != b$value
    parts <- strsplit(keys[differs], "|", fixed = TRUE)
    leaf_diff <- data.frame(lane = vapply(parts, `[`, "", 1L), step = vapply(parts, `[`, "", 2L),
        state = vapply(parts, `[`, "", 3L), what = vapply(parts, `[`, "", 4L), path = vapply(parts, `[`, "", 5L),
        base = a$value[differs], head = b$value[differs], stringsAsFactors = FALSE)
    expected_paths <- c("metaboliteIntensityFiltering/metabolites_min_values_per_group",
        "metaboliteIntensityFiltering/metabolites_min_groups")
    ok_leaf <- grepl(paste0("(^|/)(", paste(expected_paths, collapse = "|"), ")$"), leaf_diff$path) & is.na(leaf_diff$base) &
        ((grepl("min_values_per_group$", leaf_diff$path) & leaf_diff$head == "2") |
         (grepl("min_groups$", leaf_diff$path) & leaf_diff$head == "1"))
    leaf_diff$predicted <- ok_leaf
    rows[[2L]] <- data.frame(check = "P2_leaves", lane = "", step = "", state = "", item = paste(nrow(leaf_diff), "differing leaves"),
        equal = NA, predicted = "only the two keys", holds = all(ok_leaf) && nrow(leaf_diff) > 0L, stringsAsFactors = FALSE)
    utils::write.table(leaf_diff, file.path(out, "leaf-diff.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
}
result <- do.call(rbind, rows)
utils::write.table(result, file.path(out, "compare.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
if (!is.null(result)) print(stats::aggregate(holds ~ check, result, function(x) paste0(sum(x), "/", length(x))))
cat("side status:", paste(names(status), status, collapse = ", "), "\n")
quit(status = if (any(status != 0L) || is.null(result) || !all(result$holds)) 1L else 0L)
