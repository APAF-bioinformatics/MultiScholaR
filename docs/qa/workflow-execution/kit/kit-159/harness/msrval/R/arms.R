# Arm kinds' gates and hooks (KIT-159 sections 7.3, 7.5, 7.8).

#' Are two DA tables the same (digest, then values to 1e-12)?
#' @noRd
valSameDaTable <- function(a, b, tolerance = 1e-12) {
    if (identical(digest::digest(a, algo = "sha256"), digest::digest(b, algo = "sha256"))) {
        return(TRUE)
    }
    if (!identical(dim(a), dim(b)) || !identical(names(a), names(b))) return(FALSE)
    for (col in names(a)) {
        x <- a[[col]]
        y <- b[[col]]
        if (is.numeric(x) && is.numeric(y)) {
            if (!identical(is.na(x), is.na(y))) return(FALSE)
            if (any(abs(x - y) > tolerance, na.rm = TRUE)) return(FALSE)
        } else if (!identical(as.character(x), as.character(y))) {
            return(FALSE)
        }
    }
    TRUE
}

#' The kernel identity gate (section 7.5; A21, A30)
#'
#' On units drawn by the declared `identity_seed` (never by hand) plus one unit per refusal step
#' seen in the blind lifecycle, the kernel's DA table equals the workflow arm's (digest, then
#' values to 1e-12) and its refusals occur at the same step with the same class (or message
#' digest).
#'
#' @param units The unit ids.
#' @param kernel,workflow Functions of a unit returning `list(refused, step, class, da)`.
#' @param identity_seed The declared seed.
#' @param n Units drawn (at least three).
#' @param refusal_units Units of the refusal steps seen in the lifecycle.
#' @return A list: `passed`, `units`, `note` (`refusal_path_not_exercised` when no refusal unit).
#' @export
valKernelIdentity <- function(units, kernel, workflow, identity_seed, n = 3L,
                              refusal_units = character()) {
    pool <- setdiff(units, refusal_units)
    drawn <- valWithSeed(identity_seed, pool[sample.int(length(pool), min(n, length(pool)))])
    chosen <- sort(unique(c(drawn, refusal_units)))
    for (u in chosen) {
        w <- workflow(u)
        k <- kernel(u)
        same <- if (isTRUE(w$refused) || isTRUE(k$refused)) {
            isTRUE(w$refused) && isTRUE(k$refused) && identical(w$step, k$step) &&
                identical(w$class %||% w$message_sha256, k$class %||% k$message_sha256)
        } else {
            valSameDaTable(w$da, k$da)
        }
        if (!same) valAbort("the kernel differs from the workflow arm on unit {u}",
                            "val_kernel_identity_failed")
    }
    list(passed = TRUE, units = chosen,
         note = if (length(refusal_units)) "" else "refusal_path_not_exercised")
}

#' Call a post-cell hook outside test mode (section 7.3; revision 1.8 item 2, C22)
#'
#' The positional arguments are as in revision 1.6; `unit`, `arm`, `design_path` and `out_dir` are
#' passed only where the hook's formals name them explicitly (`...` does not receive them).
#'
#' @param hook The ticket's hook function.
#' @param manager,results,summary_dir,tables The positional arguments.
#' @param unit,arm,design_path,out_dir The named arguments.
#' @return The hook's value.
#' @export
valCallHook <- function(hook, manager, results, summary_dir, tables, unit = NULL, arm = NULL,
                        design_path = NULL, out_dir = NULL) {
    named <- list(unit = unit, arm = arm, design_path = design_path, out_dir = out_dir)
    named <- named[names(named) %in% names(formals(hook))]
    args <- c(list(manager, results, summary_dir, tables), named)
    withr::with_options(list(multischolar.test_mode = FALSE), {
        withr::with_envvar(c(MULTISCHOLAR_TEST_MODE = "false"), do.call(hook, args))
    })
}

#' Copy `prepareInputs()$result_tables` into the run folder and check them (revision 1.8 item 3)
#'
#' Each table must be a declared output with `key: "none"`; it is written under `run_dir` and
#' checked as a result table; an undeclared one is `val_publish_forbidden`.
#'
#' @param tables Named list: declared path = data frame.
#' @param outputs The declaration's `outputs`.
#' @param run_dir `runs/<label>/` inside `val-out/`.
#' @param arms,units Ids for placeholder expansion.
#' @param universe Keys of `universe.tsv` for a `not_permitted` lineage, or `NULL`.
#' @return The paths written.
#' @export
valCopyResultTables <- function(tables, outputs, run_dir, arms = character(), units = character(),
                                universe = NULL) {
    declared <- list()
    for (o in outputs) {
        if (!identical(o[["key"]] %||% "none", "none") || isTRUE(o[["handoff"]])) next
        for (p in valExpandOutputPaths(o, arms = arms, units = units)) declared[[p]] <- o
    }
    for (nm in names(tables)) {
        o <- declared[[nm]]
        if (is.null(o)) {
            valAbort("result table {nm} is not a declared output", "val_publish_forbidden")
        }
        path <- file.path(run_dir, nm)
        dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
        utils::write.table(tables[[nm]], path, sep = "\t", quote = FALSE, row.names = FALSE)
        cols <- vapply(o[["columns"]], function(x) x[["type"]], character(1))
        names(cols) <- vapply(o[["columns"]], function(x) x[["name"]], character(1))
        labels <- lapply(o[["labels"]] %||% list(), function(x) as.character(unlist(x)))
        valCheckResultTable(path, nm, cols, labels, o[["max_rows"]] %||% 1000L, universe)
    }
    names(tables)
}
