# Units, chunks and processes (KIT-159 section 7.6; contract rows 3, 4, 6).

UNIT_ROLES <- c("independent", "full", "deletion", "draw", "dependent")
UNIT_FORMS <- c("unit", "replicate", "block")

#' Read and check a unit table (`units-<dataset>.tsv`)
#'
#' @param units A data frame or the path of a unit table.
#' @return The unit table, `run` defaulting to `unit`.
#' @export
valUnits <- function(units) {
    if (is.character(units)) {
        units <- utils::read.delim(units, colClasses = "character", na.strings = character(),
                                   stringsAsFactors = FALSE)
    }
    if (!all(c("unit", "unit_role") %in% names(units))) {
        valAbort("a unit table needs unit and unit_role", "val_spec_invalid")
    }
    if (anyDuplicated(units$unit)) valAbort("a unit id is repeated", "val_spec_invalid")
    if (!all(units$unit_role %in% UNIT_ROLES)) valAbort("an unknown unit_role", "val_spec_invalid")
    if (is.null(units$run)) units$run <- units$unit
    units$run[!nzchar(units$run)] <- units$unit[!nzchar(units$run)]
    if (!is.null(units$form) && !all(units$form %in% c(UNIT_FORMS, ""))) {
        valAbort("an unknown unit form", "val_spec_invalid")
    }
    units
}

#' Delete replicate `r` from every group (form (ii); SCI-08)
#'
#' A group without replicate `r` is left whole, and the deletion is marked `unequal_deletion`
#' (KIT-170 N15). Replicate numbers come from the committed design, never from a run (A41).
#'
#' @param design Data frame with `sample`, `group`, `replicate`.
#' @param r The replicate number to delete.
#' @return A list: `design` (the remaining samples) and `unequal_deletion`.
#' @export
valUnitsDropReplicate <- function(design, r) {
    checkmate::assert_data_frame(design)
    checkmate::assert_names(names(design), must.include = c("sample", "group", "replicate"))
    has_r <- tapply(design$replicate == r, design$group, any)
    list(design = design[design$replicate != r, , drop = FALSE],
         unequal_deletion = !all(has_r))
}

#' Signal a refused step from a cell runner (`val_arm_refused`, recorded with step and class)
#'
#' @param step Step id.
#' @param class The refusal's class (`unclassed` when the package exposes none).
#' @param message_sha256 SHA-256 of the masked message, or `NULL`.
#' @return Never returns.
#' @export
valArmRefused <- function(step, class = "unclassed", message_sha256 = NULL) {
    cnd <- structure(class = c("val_arm_refused", "error", "condition"),
                     list(message = sprintf("refused at %s (%s)", step, class), call = NULL,
                          step = step, refusal_class = class, message_sha256 = message_sha256))
    stop(cnd)
}

#' Peak resident memory of this process (VmHWM), MiB
#' @noRd
valVmHwm <- function() {
    line <- tryCatch(grep("^VmHWM", readLines("/proc/self/status"), value = TRUE),
                     error = function(e) character())
    if (length(line) == 0L) return(NA_real_)
    as.numeric(gsub("[^0-9]", "", line[[1L]])) / 1024
}

#' Run one cell with its own elapsed-time limit, catching its conditions
#' @noRd
valRunOneCell <- function(cell, run_cell, cell_seconds) {
    valWithSeed(as.integer(cell$seed), {
        before <- get(".Random.seed", envir = globalenv())
        t0 <- Sys.time()
        res <- tryCatch({
            setTimeLimit(elapsed = cell_seconds, transient = TRUE)
            value <- run_cell(cell)
            setTimeLimit(elapsed = Inf)
            list(exit = "completed", value = value)
        }, val_arm_refused = function(e) {
            setTimeLimit(elapsed = Inf)
            list(exit = "refused", step = e$step, class = e$refusal_class,
                 message_sha256 = e[["message_sha256"]] %||% NA_character_)
        }, error = function(e) {
            setTimeLimit(elapsed = Inf)
            msg <- conditionMessage(e)
            reason <- if (grepl("elapsed time limit", msg)) "timeout" else "error"
            list(exit = "did_not_complete", stop_reason = reason,
                 message_sha256 = digest::digest(msg, algo = "sha256", serialize = FALSE))
        })
        after <- get(".Random.seed", envir = globalenv())
        res$random_state_changed <- !identical(before, after)
        res$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        res
    })
}

#' Run a chunk of cells in this process (the child side of [valRunCells()])
#'
#' Writes a `<run>.started` marker before each run and `<run>.rds` when it ends, so the parent can
#' tell a killed run from runs never started. `VmHWM` is read at each cell boundary.
#'
#' @param cells_rds RDS of the chunk's cells (one row per run).
#' @param runner R file defining `runCell(cell)`.
#' @param out_dir The chunk's folder.
#' @param cell_seconds Each cell's elapsed-time limit.
#' @return `TRUE`, invisibly.
#' @export
valRunChunk <- function(cells_rds, runner, out_dir, cell_seconds) {
    cells <- readRDS(cells_rds)
    env <- new.env(parent = globalenv())
    sys.source(runner, envir = env)
    writeLines(format(as.numeric(Sys.time()), digits = 15), file.path(out_dir, "chunk.start"))
    for (i in seq_len(nrow(cells))) {
        cell <- as.list(cells[i, ])
        writeLines("started", file.path(out_dir, paste0(cell$run, ".started")))
        res <- valRunOneCell(cell, env$runCell, cell_seconds)
        res$vmhwm_mib <- valVmHwm()
        saveRDS(res, file.path(out_dir, paste0(cell$run, ".rds")))
    }
    writeLines(format(as.numeric(Sys.time()), digits = 15), file.path(out_dir, "chunk.end"))
    invisible(TRUE)
}

#' Start chunk processes at the same time and wait for them
#' @noRd
valStartChunks <- function(groups, runs, runner, work_dir, cell_seconds, tag, chunk_seconds) {
    procs <- lapply(seq_along(groups), function(k) {
        d <- file.path(work_dir, sprintf("chunk-%s%02d", tag, k))
        dir.create(d, recursive = TRUE, showWarnings = FALSE)
        rds <- file.path(d, "cells.rds")
        saveRDS(runs[runs$run %in% groups[[k]], , drop = FALSE], rds)
        code <- sprintf("msrval::valRunChunk('%s', '%s', '%s', %f)", rds, runner, d, cell_seconds)
        list(dir = d, runs = groups[[k]], proc = valRscriptProcess(code, file.path(d, "chunk.log")))
    })
    for (i in seq_along(procs)) {
        procs[[i]]$proc$wait(timeout = if (is.finite(chunk_seconds)) chunk_seconds * 1000 else -1)
        procs[[i]]$killed <- procs[[i]]$proc$is_alive()
        if (procs[[i]]$killed) procs[[i]]$proc$kill()
    }
    procs
}

#' Collect a chunk's run results
#' @noRd
valCollectChunk <- function(p, k, tag, k_total) {
    status <- p$proc$get_exit_status()
    reason <- if (isTRUE(p$killed)) "timeout" else if (identical(status, -9L)) "memory" else "error"
    read_time <- function(f) {
        path <- file.path(p$dir, f)
        if (file.exists(path)) as.numeric(readLines(path, warn = FALSE)[[1L]]) else NA_real_
    }
    rows <- lapply(p$runs, function(r) {
        rds <- file.path(p$dir, paste0(r, ".rds"))
        started <- file.exists(file.path(p$dir, paste0(r, ".started")))
        res <- if (file.exists(rds)) {
            readRDS(rds)
        } else if (started) {
            list(exit = "did_not_complete", stop_reason = reason)
        } else {
            list(exit = "not_started")
        }
        data.frame(run = r, exit = res$exit, stop_reason = res[["stop_reason"]] %||% NA_character_,
                   step = res[["step"]] %||% NA_character_,
                   class = res[["class"]] %||% NA_character_,
                   message_sha256 = res[["message_sha256"]] %||% NA_character_,
                   seconds = res[["seconds"]] %||% NA_real_,
                   random_state_changed = res[["random_state_changed"]] %||% NA,
                   vmhwm_mib = res[["vmhwm_mib"]] %||% NA_real_,
                   chunk = sprintf("%s%02d", tag, k), stringsAsFactors = FALSE)
    })
    chunk <- data.frame(chunk = sprintf("%s%02d", tag, k), started = read_time("chunk.start"),
                        ended = read_time("chunk.end"), exit = status, k = as.integer(k_total),
                        stringsAsFactors = FALSE)
    list(rows = do.call(rbind, rows), chunk = chunk)
}

#' Run cells in concurrent chunks (section 7.6, A37)
#'
#' Rows naming one `run` are executed once and reported for each row (A39). Runs are dealt to `k`
#' chunks that run at the same time, each an `Rscript` with the thread variables at 1. Each cell has
#' its own elapsed-time limit; the chunk's limit is the backstop. A killed chunk's running cell is
#' `did_not_complete` (`memory` for SIGKILL, `timeout` for the chunk limit); its cells never started
#' run in a fresh chunk within the same call (not a retry); any still not started are `not_started`.
#'
#' @param cells Data frame: `cell`, `unit`, `run`, `arm`, `seed` and the runner's own columns.
#' @param runner R file defining `runCell(cell)`.
#' @param k Concurrent chunks.
#' @param cell_minutes Each cell's limit.
#' @param work_dir Working folder (private).
#' @param chunk_minutes The chunk limit (backstop), or `Inf`.
#' @return One row per input cell, with attribute `chunks` (one row per chunk).
#' @export
valRunCells <- function(cells, runner, k = 4L, cell_minutes, work_dir, chunk_minutes = Inf) {
    checkmate::assert_data_frame(cells, min.rows = 1L)
    checkmate::assert_names(names(cells), must.include = c("cell", "unit", "run", "arm", "seed"))
    runs <- cells[!duplicated(cells$run), , drop = FALSE]
    k1 <- min(k, nrow(runs))
    groups <- split(runs$run, rep_len(seq_len(k1), nrow(runs)))
    procs <- valStartChunks(groups, runs, runner, work_dir, cell_minutes * 60, "a",
                            chunk_minutes * 60)
    got <- lapply(seq_along(procs), function(i) valCollectChunk(procs[[i]], i, "a", k1))
    res <- do.call(rbind, lapply(got, `[[`, "rows"))
    chunks <- do.call(rbind, lapply(got, `[[`, "chunk"))
    res$fresh_chunk <- FALSE
    waiting <- res$run[res$exit == "not_started"]
    if (length(waiting) > 0L) {
        p2 <- valStartChunks(list(waiting), runs, runner, work_dir, cell_minutes * 60, "b",
                             chunk_minutes * 60)
        g2 <- valCollectChunk(p2[[1L]], 1L, "b", 1L)
        g2$rows$fresh_chunk <- TRUE
        res <- rbind(res[!res$run %in% waiting, , drop = FALSE], g2$rows)
        chunks <- rbind(chunks, g2$chunk)
    }
    per_cell <- res[match(cells$run, res$run), setdiff(names(res), "run")]
    rownames(per_cell) <- NULL
    base <- cells[, c("cell", "unit", "run", "arm")]
    rownames(base) <- NULL
    out <- cbind(base, per_cell)
    rownames(out) <- NULL
    attr(out, "chunks") <- chunks
    out
}
