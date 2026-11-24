# Replay of declared cells (KIT-159 7.3, G5; audit I07, re-check E5): the cell's analysis record,
# a copy of its project store pruned to the record's import and design states by the package's own
# `kitPrepareReplayProject` (tools/profiling/workflow_kit.R, extracted, never edited), the record
# replayed there by `replayWorkflowAnalysisRecord`, and the two stores compared state by state and
# on the head's DA result by the package's own fingerprints (KIT-028's comparison).

#' Replay statuses a cell's provenance may carry (G5 plus the two for cells not replayed)
#' @noRd
REPLAY_STATUSES <- c("identical", "differs", "unsupported", "not_declared", "not_run")

#' The workflow_kit functions a replay uses (extracted with their closure)
#' @noRd
REPLAY_KIT_FUNCTIONS <- c("kitPrepareReplayProject", "kitFingerprintStore",
                          "kitReproduceMismatches", "kitReadHead", "kitOpenStore",
                          "kitFingerprintReference")

#' Is a cell declared for replay? (`replay`: a list of `{arms, units}`; no `units` means all)
#'
#' @param decl The declaration.
#' @param unit,arm The cell.
#' @return `TRUE` or `FALSE`.
#' @export
valReplayDeclared <- function(decl, unit, arm) {
    any(vapply(decl[["replay"]] %||% list(), function(r) {
        arm %in% unlist(r[["arms"]]) && (is.null(r[["units"]]) || unit %in% unlist(r[["units"]]))
    }, logical(1)))
}

#' The replay kit: the package's store preparation and fingerprints, extracted from a checkout's
#' `tools/profiling/workflow_kit.R` (only named function definitions; never top-level code)
#'
#' @param pkg_root The package checkout.
#' @return A list of functions: `prepare(root, preconditions, api)`, `held(root, api)` (the
#'   prepared head's states), `fingerprint(root, api)`, `mismatches(origin, replay)` (state
#'   fingerprints), `da(root, api)` (the head's DA result fingerprint or `NA`); attribute
#'   `helpers` (file, function, sha256_of_source).
#' @export
valReplayKit <- function(pkg_root) {
    env <- valExtractHelpers(file.path(pkg_root, "tools", "profiling", "workflow_kit.R"),
                             REPLAY_KIT_FUNCTIONS, parent = asNamespace("MultiScholaR"))
    entities <- function(root) file.path(root, ".workflow-entities")
    pops <- function(fp) lapply(fp$states, function(x) list(fingerprint = x))
    kit <- list(
        prepare = env$kitPrepareReplayProject,
        held = function(root, api) {
            h <- env$kitReadHead(entities(root), api)
            union(h$history, names(h$references))
        },
        fingerprint = function(root, api) env$kitFingerprintStore(entities(root), api),
        mismatches = function(origin, replay) {
            as.vector(env$kitReproduceMismatches(pops(origin), pops(replay)))
        },
        da = function(root, api) {
            store <- env$kitOpenStore(entities(root), api)
            on.exit(api("closeWorkflowEntityStore")(store), add = TRUE)
            head <- readRDS(file.path(entities(root), "dia-state-head.rds"))
            snapshot <- api("hydrateSavedWorkflowCheckpoint")(store, head)
            if (is.null(snapshot$da_results)) return(NA_character_)
            env$kitFingerprintReference(store, snapshot$da_results, api)$sha256
        }
    )
    attr(kit, "helpers") <- attr(env, "helpers")
    kit
}

#' The record's precondition states: its leading import, design and integration-import steps
#' (as `replayWorkflowAnalysisRecord` reads them)
#' @noRd
valReplayPreconditions <- function(record) {
    ops <- c("import", "design_builder", "integration_import")
    steps <- record$steps %||% list()
    is_pre <- vapply(steps, function(s) as.character(s$operation %||% "") %in% ops, logical(1))
    n <- sum(cumprod(is_pre))
    vapply(steps[seq_len(n)], function(s) as.character(s$state), character(1))
}

#' The record's input files, found by SHA-256 under the cell's folders (the uploaded copies and the
#' design the GUI wrote); a role whose file is not found is left out (the replay then refuses it)
#' @noRd
valReplayInputs <- function(record, roots) {
    files <- unlist(lapply(roots[dir.exists(roots)], list.files, recursive = TRUE,
                           full.names = TRUE, all.files = TRUE))
    files <- files[!grepl("/[.]workflow-entities/", files)]
    out <- character()
    for (inp in record$inputs %||% list()) {
        sha <- inp$sha256
        if (is.null(sha) || is.na(sha) || is.null(inp$role)) next
        cand <- files[basename(files) == as.character(inp$basename %||% "")]
        if (!length(cand)) {
            cand <- files[file.size(files) == suppressWarnings(as.numeric(inp$bytes %||% -1))]
        }
        hit <- Filter(function(f) identical(valSha256(f), sha), cand)
        if (length(hit)) out[[as.character(inp$role)]] <- hit[[1L]]
    }
    out
}

#' Compare an origin store with its replay (the package's fingerprints)
#'
#' @param kit From [valReplayKit()] (or a list with `mismatches`).
#' @param origin,replay `fingerprint()` results of the two stores.
#' @param origin_da,replay_da The head DA result fingerprints.
#' @return A list: `status` (`identical`, `differs`), `identical`, `da_identical`, `differences`.
#' @export
valReplayCompare <- function(kit, origin, replay, origin_da, replay_da) {
    diffs <- kit$mismatches(origin, replay)
    da_same <- identical(origin_da, replay_da)
    if (!da_same) diffs <- c(diffs, "da_results")
    same <- !length(diffs)
    list(status = if (same) "identical" else "differs", identical = same,
         da_identical = da_same, differences = utils::head(as.character(diffs), 50L))
}

#' Replay a record on a pruned copy of the origin store and compare (G5)
#'
#' The origin store is copied, never touched; the copy is pruned to the record's import and design
#' states (`kit$prepare`), the record replayed there, and both stores compared. A refusal or error
#' of the preparation or replay is `unsupported` with its class.
#'
#' @param record The analysis record (as read with `simplifyVector = FALSE`).
#' @param origin_root The origin store root (the folder holding `.workflow-entities`).
#' @param inputs Named character: record input role to file.
#' @param paths The replay's experiment paths.
#' @param scratch A private folder for the copy and the replay's scratch.
#' @param kit From [valReplayKit()].
#' @param api `function(name)` returning a MultiScholaR function.
#' @param replay The replay function (default the package's `replayWorkflowAnalysisRecord`).
#' @return A list: `status`, `class`, `identical`, `da_identical`, `preconditions`,
#'   `prepared_states`, `differences`.
#' @export
valReplayStore <- function(record, origin_root, inputs, paths, scratch, kit,
                           api = function(name) get(name, asNamespace("MultiScholaR")),
                           replay = NULL) {
    replay <- replay %||% api("replayWorkflowAnalysisRecord")
    pre <- valReplayPreconditions(record)
    copy <- file.path(scratch, "store")
    dir.create(copy, recursive = TRUE, showWarnings = FALSE)
    out <- list(status = "unsupported", class = NA_character_, identical = FALSE,
                da_identical = NA, preconditions = pre, prepared_states = character(),
                differences = character())
    got <- tryCatch({
        if (system2("cp", c("-a", file.path(origin_root, "."), copy)) != 0L) {
            valAbort("the origin store could not be copied", "val_replay_copy")
        }
        prepared <- kit$prepare(copy, pre, api)
        held <- kit$held(prepared, api)
        replay(record, prepared, inputs, paths, file.path(scratch, "scratch"))
        list(prepared = prepared, held = held)
    }, error = function(e) e)
    if (inherits(got, "error")) {
        out$class <- class(got)[[1L]]
        return(out)
    }
    out$prepared_states <- got$held
    cmp <- tryCatch(valReplayCompare(kit, kit$fingerprint(origin_root, api),
                                     kit$fingerprint(got$prepared, api),
                                     kit$da(origin_root, api), kit$da(got$prepared, api)),
                    error = function(e) e)
    if (inherits(cmp, "error")) {
        out$class <- class(cmp)[[1L]]
        return(out)
    }
    utils::modifyList(out, cmp)
}

#' Replay a completed workflow cell from its folder (called after its session has closed)
#'
#' @param work_dir The cell's work folder (its project under `project/<route>`, its record
#'   `analysis-record.json`, its uploaded inputs under `inputs/`).
#' @param route,omic The cell's route and omic.
#' @param pkg_root The package checkout (for the replay kit).
#' @param kit From [valReplayKit()] (default from `pkg_root`).
#' @return As [valReplayStore()]; `unsupported` with class `val_replay_no_record` when the cell
#'   wrote no record, with the record's build error class when building it failed.
#' @noRd
valReplayCell <- function(work_dir, route, omic, pkg_root, kit = NULL) {
    rec_path <- file.path(work_dir, "analysis-record.json")
    err_path <- file.path(work_dir, "analysis-record.error")
    none <- function(cls) {
        list(status = "unsupported", class = cls, identical = FALSE, da_identical = NA,
             preconditions = character(), prepared_states = character(),
             differences = character())
    }
    if (!file.exists(rec_path)) {
        cls <- if (file.exists(err_path)) readLines(err_path)[[1L]] else "val_replay_no_record"
        return(none(cls))
    }
    kit <- tryCatch(kit %||% valReplayKit(pkg_root), error = function(e) e)
    if (inherits(kit, "error")) return(none(class(kit)[[1L]]))
    record <- jsonlite::read_json(rec_path, simplifyVector = FALSE)
    project <- file.path(work_dir, "project", route)
    store <- list.files(project, "^store[.]json$", recursive = TRUE, all.files = TRUE,
                        full.names = TRUE)
    if (!length(store)) return(none("val_replay_no_store"))
    origin_root <- dirname(dirname(store[[1L]]))
    inputs <- valReplayInputs(record, c(file.path(work_dir, "inputs"), project))
    scratch <- file.path(work_dir, "replay")
    ns <- asNamespace("MultiScholaR")
    dirs <- suppressMessages(
        ns$setupDirectories(base_dir = file.path(scratch, "project"), omic_types = omic,
                            label = NULL, force = TRUE, reuse_existing = TRUE)
    )
    paths <- dirs[[omic]]
    valReplayStore(record, origin_root, inputs, paths, scratch, kit)
}

#' Build and write a cell's analysis record inside its session (before the owner closes); a build
#' error is kept as its class in `analysis-record.error`
#' @noRd
valWriteCellRecord <- function(manager, work_dir) {
    ns <- asNamespace("MultiScholaR")
    path <- file.path(work_dir, "analysis-record.json")
    got <- tryCatch(ns$withWorkflowStateScope(manager, function() {
        ns$writeWorkflowAnalysisRecord(manager, path)
    }), error = function(e) e)
    if (inherits(got, "error")) {
        unlink(path)
        writeLines(class(got)[[1L]], file.path(work_dir, "analysis-record.error"))
    }
    invisible(path)
}

#' The provenance `replay` block of one cell
#' @noRd
valReplayBlock <- function(replay, declared, completed) {
    if (!declared) return(list(status = "not_declared"))
    if (!completed || is.null(replay)) return(list(status = "not_run"))
    list(status = replay$status, class = replay$class, identical = isTRUE(replay$identical),
         da_identical = replay$da_identical, preconditions = as.list(replay$preconditions),
         prepared_states = as.list(replay$prepared_states),
         differences = as.list(replay$differences))
}
