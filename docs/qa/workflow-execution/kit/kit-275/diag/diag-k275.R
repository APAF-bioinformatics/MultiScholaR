# KIT-275 diagnostic (VAL-392 (1)): NON-COUNTED. Does the between-sample normalisation explain why the pre-RUV correlation filter at 0.5
# failed the planted sample in only 184, 135 and 45 of 200 data sets (VAL-384)? The arm without RUV (NO) on the planted data sets, at
# D3, D4 and D6, under three between-sample normalisations (the shipped default, quantile, none), 50 data sets per design on fresh
# seeds, the same 50 under the three normalisations. Evidence code under docs/qa: no package code. Nothing is decided by it.
#
# Usage: Rscript diag-k275.R <designs> <first t> <last t> <out dir>      designs: a comma list of D3, D4, D6; t = 1 to 50, r = 9500 + t
#        Rscript diag-k275.R --unit <design> <r> <rds>                    one data set, three journeys (child process)
# Environment: K275_WORKERS child processes at a time (default 4); K275_DIAG_JOB a name for the job's files (default k275-diag);
#   K275_DIAG_DECLARATION the declaration file (its SHA-256 goes to the provenance).
# Re-use, read only: sim-k275.R is SOURCED (`options(k275.library = TRUE)`, as fixture-k275.R does) after its SHA-256 is checked against the
# bound value; generate(), makeVersion(), routeDefaults(), requestK275(), journeyL(), journeyRow(), bestCorrelation(), failsAt() are the
# counted job's own, unchanged. Nothing of sim-k275.R, score-k275.R or fixture-k275.R is copied or edited. The ONE addition is the stub
# below (K275DIAG-stop): journeyL() resolves `runWorkflowOmicsReferenceDaAnalysis` from the global environment, so a global function of
# that name ends each journey at the step the differential analysis would start, after the pre-RUV read and the package's own correlation
# filter. No DA is run or recorded. If the stub did not take (the DA ran), the row is an instrument row (da_calls = 0).

diag_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]))
sim_file <- file.path(dirname(dirname(diag_file)), "sim-k275.R")
SIM_SHA16 <- "2cc87ea569eef7fd"   # the audited and counted sim-k275.R (declaration.md section 1 (b))
sim_sha <- strsplit(system2("sha256sum", shQuote(sim_file), stdout = TRUE), " ")[[1L]][[1L]]
if (!startsWith(sim_sha, SIM_SHA16)) stop("K275DIAG: sim-k275.R is not the bound script (SHA-256 ", sim_sha, ")")
options(k275.library = TRUE)
source(sim_file, local = FALSE)

# ---- constants ---------------------------------------------------------------------------------------------------------
DIAG_NORMS <- c("cyclicloess", "quantile", "none")   # the first is the shipped default (asserted against the route's read default)
DIAG_T <- 1:50
diagSeed <- function(t) 9500L + t   # r = 9501 to 9550; generator 275000000 + 100000 c + 10000 s + r (s = 2, N1), planting 276000000 + 100000 c + r
DIAG_SEEDS <- diagSeed(DIAG_T)
DIAG_DESIGN_ORDER <- c("D6", "D4", "D3")   # within a data set index, the slowest design first
DIAG_COLUMNS <- c("design", "t", "r", "norm", "status", "class", "step", "seconds", "warnings", "n_runs", "planted_run", "genuine_run",
    "cy_planted", "planted_below_05", "cy_genuine", "n_fail_y", "n_removed_y", "e5_y", "da_calls", "plant_ok", "norm_acted_y")
diag_state <- new.env()
diag_state$da_calls <- 0L

# K275DIAG-stop: the stub that ends a journey where the DA would start (see the header). A classed error, caught by journeyL()'s own handler.
runWorkflowOmicsReferenceDaAnalysis <- function(...) {   # K275DIAG-stop
    diag_state$da_calls <- diag_state$da_calls + 1L   # K275DIAG-stop
    stop(structure(class = c("k275diag_stop", "error", "condition"), list(message = "K275DIAG: journey stopped before the DA", call = NULL)))   # K275DIAG-stop
}   # K275DIAG-stop

diagRowFrame <- function(row) {
    stopifnot(setequal(names(row), DIAG_COLUMNS))
    row <- row[DIAG_COLUMNS]
    row[] <- lapply(row, function(x) if (length(x) != 1L) paste(x, collapse = ";") else x)
    as.data.frame(row, stringsAsFactors = FALSE)
}
#' The status of a diagnostic journey: ok only when it ended at the stub as designed (or, rarely, the package filter removed a whole group so
#' that the DA was not reached); every other outcome keeps the journey's own status (instrument, crashed, refused) or is an instrument row.
diagStatus <- function(j, da_calls) {
    status <- as.character(j[["status"]])
    class <- as.character(j[["class"]])
    stopped <- identical(j[["error_class"]], "k275diag_stop") && identical(status, "crashed") && identical(j[["error_step"]], "da")
    if (stopped && da_calls == 1L) return(list(status = "ok", class = "stopped_before_da"))
    if (identical(status, "no_da") && !is.na(class) && startsWith(class, "filter_removed_group:") && da_calls == 0L) return(list(status = "ok", class = class))
    if (stopped || (da_calls != 1L && status %in% c("ok", "no_da"))) return(list(status = "instrument", class = paste0("stub_calls_", da_calls, ": ", class)))
    list(status = status, class = class)
}

#' One data set: V0 and V1g as the counted tier D builds them (mask "none"), then the NO arm's journey under each normalisation.
diagUnit <- function(design, r, rows_path, checks_path) {
    stopifnot(r %in% DIAG_SEEDS)
    p <- CELLS[[design]]
    defaults <- routeDefaults(p[["omic"]])
    stopifnot(identical(defaults[["method"]], DIAG_NORMS[[1L]]))
    per <- p[["per"]]
    n_runs <- p[["groups"]] * per
    a_run <- sprintf("R%03d", n_runs)
    g_run <- sprintf("R%03d", per + 1L)   # the first run of group B: a genuine sample of the planted sample's group
    v0 <- makeVersion(design, r, "V0", "none")
    v1 <- makeVersion(design, r, "V1g", "none", v0)
    sim <- v1[["sim"]]
    y0 <- v0[["sim"]][["y"]]
    y1 <- sim[["y"]]
    in_y <- unname(sim[["assay_of"]][rownames(y1)] == Y_ASSAY)
    planted_cells <- outer(in_y, colnames(y1) == a_run, "&")
    plant_ok <- identical(dim(y0), dim(y1)) && identical(as.vector(y0 != y1), as.vector(planted_cells))
    addCheck("D1_plant", design, r, "V1g against V0", NA, "differs only in the planted sample's Y entries", plant_ok)
    ctx <- list(tier = "DIAG", design = design, cell = design, omic = p[["omic"]], r = r, per = per, smoke = FALSE,
        threshold = defaults[["correlation"]], defaults = defaults, shares = plantedShares(design, r))
    for (norm in DIAG_NORMS) {
        diag_state$da_calls <- 0L
        request <- requestK275("NO", defaults, sim)
        j <- journeyL(sim, v1[["values"]], "NO", request, defaults, tempfile(sprintf("k275diag-%s-", norm)), filter = TRUE, itsd = FALSE,
            method = norm, correlation = TRUE)
        da_calls <- diag_state$da_calls
        row <- journeyRow(ctx, j, sim, "V1g", "NO", FALSE)
        status <- diagStatus(j, da_calls)
        if (identical(status[["status"]], "instrument") && isTRUE(startsWith(status[["class"]], "partial_match_guard:"))) {   # a guard stop inside a journey: leave its trace
            tag <- sub("\\.tsv$", "", sub("^checks-", "", basename(checks_path)))
            try(writeLines(c(status[["class"]], getOption("k275.guard_trace", "")), file.path(dirname(checks_path), sprintf("unit-trace-%s-%s.txt", tag, norm))), silent = TRUE)
        }
        pre <- j[["pairs_pre"]]
        cy <- bestCorrelation(pre[[Y_ASSAY]], sim[["runs"]])
        removed <- j[["corr_removed"]]
        nc <- j[["norm_check"]]
        out <- list(design = design, t = r - 9500L, r = r, norm = norm, status = status[["status"]], class = oneLine(status[["class"]]),
            step = row[["step"]], seconds = row[["seconds"]], warnings = row[["warnings"]], n_runs = row[["n_runs"]], planted_run = a_run, genuine_run = g_run,
            cy_planted = row[["cy_pre_a"]], planted_below_05 = row[["fa_pre_05"]], cy_genuine = unname(cy[[g_run]]),
            n_fail_y = sum(failsAt(cy, THETAS[[1L]])), n_removed_y = if (is.null(removed)) NA_integer_ else length(removed[[Y_ASSAY]]),
            e5_y = row[["e5_y"]], da_calls = da_calls, plant_ok = plant_ok,
            norm_acted_y = if (is.null(nc)) NA else isTRUE(nc[["norm_ok"]][nc[["assay"]] == Y_ASSAY]))
        appendTable(rows_path, diagRowFrame(out))
        addCheck("D2_status", design, r, norm, status[["status"]], "ok", identical(status[["status"]], "ok"), status[["class"]])
        flushChecks(checks_path)
    }
    invisible(NULL)
}

# ---- the child process ---------------------------------------------------------------------------------------------------
diagChildMain <- function() {
    stopifnot(length(args) == 4L)
    rows_path <- Sys.getenv("K275_ROWS_PATH")
    checks_path <- Sys.getenv("K275_CHECKS_PATH")
    last_calls <- NULL
    result <- tryCatch(withCallingHandlers({
        diagUnit(args[[2L]], as.integer(args[[3L]]), rows_path, checks_path)
        list(ok = TRUE)
    }, warning = guardWarning, error = function(e) last_calls <<- paste(vapply(sys.calls(), function(cl) oneLine(deparse(cl)[[1L]], 80L), ""), collapse = " <- ")),
        error = function(e) list(error = oneLine(paste(class(e)[[1L]], conditionMessage(e))), trace = last_calls %||% getOption("k275.guard_trace", "")))
    flushChecks(checks_path)
    saveRDS(result, args[[4L]])
    quit(save = "no", status = 0L)
}

# ---- the job (parent) ------------------------------------------------------------------------------------------------------
diagParentMain <- function() {
    stopifnot(length(args) == 4L)
    designs <- strsplit(args[[1L]], ",", fixed = TRUE)[[1L]]
    first_t <- as.integer(args[[2L]])
    last_t <- as.integer(args[[3L]])
    stopifnot(all(designs %in% c("D3", "D4", "D6")), !is.na(first_t), !is.na(last_t), first_t >= min(DIAG_T), last_t <= max(DIAG_T), first_t <= last_t)
    workers <- as.integer(Sys.getenv("K275_WORKERS", "4"))
    job_name <- Sys.getenv("K275_DIAG_JOB", "k275-diag")
    out_dir <- args[[4L]]
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    out_dir <- normalizePath(out_dir)
    rows_dir <- file.path(out_dir, "rows")
    progress_dir <- file.path(out_dir, "progress")
    dir.create(rows_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)
    units <- do.call(rbind, lapply(seq(first_t, last_t), function(t) do.call(rbind, lapply(intersect(DIAG_DESIGN_ORDER, designs), function(design)
        data.frame(design = design, t = t, r = diagSeed(t), stringsAsFactors = FALSE)))))
    units$tag <- sprintf("DIAG-%s-%d", units$design, units$r)
    cat(sprintf("K275DIAG %s: %d units (%d journeys), %d workers\n", job_name, nrow(units), 3L * nrow(units), workers))
    e3 <- checkE3("lipid")
    addCheck("E3", "lipid", 0L, "defaults", NA, "equal to 5.2", e3$pass, paste(e3$note, "|",
        paste(names(e3$got), vapply(e3$got, function(x) paste(format(x), collapse = "/"), ""), sep = "=", collapse = "; ")))
    flushChecks(file.path(out_dir, sprintf("checks-%s-start.tsv", job_name)))
    if (!e3$pass) stop("K275DIAG: E3 failed: the route defaults at the bound commit differ from section 5.2 (see checks-", job_name, "-start.tsv)")
    declaration <- Sys.getenv("K275_DIAG_DECLARATION")
    library_table <- utils::installed.packages()[, c("Package", "Version")]
    key_packages <- c("MultiScholaR", "RUVIIIC", "limma", "duckdb", "dplyr")
    kit <- dirname(dirname(diag_file))
    shaOrNA <- function(f) if (file.exists(f)) sha(f) else NA
    provenance <- function(final) {
        tables <- if (final) list.files(out_dir, recursive = TRUE, full.names = TRUE, pattern = "\\.(tsv|json)$") else character()
        list(job = job_name, designs = designs, data_sets_t = c(first_t, last_t), seeds_r = range(units$r), norms = DIAG_NORMS, workers = workers,
            units = nrow(units), final = final, counted = FALSE,
            commit = tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE), error = function(e) NA),
            package_tree_vs_counted_tier_d = tryCatch(paste(system2("git", c("diff", "--stat", "577a6f7b9ee14b263b3e1cacd64d946c67c9130e", "HEAD", "--", "R/", "src/", "DESCRIPTION",
                "NAMESPACE", "inst/", "packaging/", ".github/actions/"), stdout = TRUE), collapse = " | "), error = function(e) NA),
            r_version = R.version.string, blas = extSoftVersion()[["BLAS"]], lapack = La_library(), omp_threads = Sys.getenv("OMP_NUM_THREADS"),
            runner_image = paste(Sys.getenv("ImageOS"), Sys.getenv("ImageVersion")),
            processor = tryCatch(sub("^model name\\s*:\\s*", "", grep("^model name", readLines("/proc/cpuinfo"), value = TRUE)[[1L]]), error = function(e) NA),
            cores = parallel::detectCores(),
            seed_rule = "set.seed(275000000 + 100000 c + 10000 s + r, kind = 'Mersenne-Twister', normal.kind = 'Inversion', sample.kind = 'Rejection'), c = 1, 2, 3, s = 2, r = 9500 + t; planting 276000000 + 100000 c + r",
            key_packages = as.list(stats::setNames(library_table[match(key_packages, library_table[, "Package"]), "Version"], key_packages)),
            script_sha256 = sha(diag_file), sim_sha256 = sim_sha, score_k275_sha256 = shaOrNA(file.path(kit, "score-k275.R")),
            fixture_sha256 = shaOrNA(file.path(kit, "fixture-k275.R")), scorer_diag_sha256 = shaOrNA(file.path(dirname(diag_file), "score-diag.R")),
            declaration_file = declaration, declaration_sha256 = if (nzchar(declaration) && file.exists(declaration)) sha(declaration) else NA,
            started = format(started_at, "%Y-%m-%dT%H:%M:%S%z"), finished = if (final) format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z") else NA,
            tables_sha256 = stats::setNames(as.list(vapply(tables, sha, character(1L))), sub(paste0("^", out_dir, "/?"), "", tables)))
    }
    # the provenance is written before any unit runs, and again at the end
    jsonlite::write_json(provenance(FALSE), file.path(out_dir, sprintf("provenance-%s.json", job_name)), auto_unbox = TRUE, pretty = TRUE, digits = NA)
    runChild <- function(i) {
        u <- units[i, ]
        rds <- file.path(progress_dir, sprintf("unit-%s.rds", u$tag))
        rows_path <- file.path(rows_dir, sprintf("part-%s.tsv", u$tag))
        checks_path <- file.path(progress_dir, sprintf("checks-%s.tsv", u$tag))
        env <- c(sprintf("K275_ROWS_PATH=%s", rows_path), sprintf("K275_CHECKS_PATH=%s", checks_path))
        t0 <- Sys.time()
        code <- system2("Rscript", c(shQuote(diag_file), "--unit", u$design, u$r, shQuote(rds)), env = env, stdout = FALSE, stderr = FALSE)
        seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        trace_file <- file.path(progress_dir, sprintf("unit-trace-%s.txt", u$tag))
        res <- if (file.exists(rds)) readRDS(rds) else list(error = paste(sprintf("child process exit %d, no result", code),
            if (file.exists(trace_file)) oneLine(paste(readLines(trace_file, n = 1L), collapse = " "), 200L) else ""))
        unlink(rds)
        if (!is.null(res[["trace"]])) writeLines(c(res[["error"]] %||% "", res[["trace"]]), trace_file)
        one <- data.frame(tag = u$tag, design = u$design, t = u$t, r = u$r, wall_seconds = round(seconds, 1), exit = code,
            error = res[["error"]] %||% "", stringsAsFactors = FALSE)
        # per-unit progress: one line on the console and one small file (no shared file is appended to by two workers)
        utils::write.table(one, file.path(progress_dir, sprintf("done-%s.tsv", u$tag)), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
        cat(sprintf("K275DIAG %s: unit %s done, exit %d, %.0f s%s\n", job_name, u$tag, code, seconds, if (nzchar(one$error)) paste0(", ERROR ", substr(one$error, 1L, 120L)) else ""))
        one
    }
    results <- parallel::mclapply(seq_len(nrow(units)), runChild, mc.cores = workers, mc.preschedule = FALSE)
    timing <- do.call(rbind, lapply(seq_along(results), function(i) if (is.data.frame(results[[i]])) results[[i]] else
        data.frame(tag = units$tag[[i]], design = units$design[[i]], t = units$t[[i]], r = units$r[[i]], wall_seconds = NA, exit = NA,
            error = oneLine(paste("worker failed:", paste(as.character(results[[i]]), collapse = " "))), stringsAsFactors = FALSE)))
    # a unit that errored, or left fewer than three rows, gets a unit row (never a denominator); a unit with no result and no trace is an instrument row
    for (i in seq_len(nrow(units))) {
        path <- file.path(rows_dir, sprintf("part-%s.tsv", units$tag[[i]]))
        have <- if (file.exists(path)) max(0L, length(readLines(path)) - 1L) else 0L
        if (nzchar(timing$error[[i]]) || have < length(DIAG_NORMS)) {
            no_trace <- !nzchar(timing$error[[i]])
            guard <- grepl("partial-match guard", timing$error[[i]], fixed = TRUE)
            row <- stats::setNames(as.list(rep(NA, length(DIAG_COLUMNS))), DIAG_COLUMNS)
            row[c("design", "t", "r", "norm")] <- list(units$design[[i]], units$t[[i]], units$r[[i]], "unit")
            row[["status"]] <- if (guard || no_trace) "instrument" else "crashed"
            row[["class"]] <- if (no_trace) sprintf("unit_incomplete_no_result_no_trace: %d of %d rows", have, length(DIAG_NORMS)) else timing$error[[i]]
            appendTable(path, diagRowFrame(row))
        }
    }
    unit_checks <- lapply(list.files(progress_dir, pattern = "^checks-.*\\.tsv$", full.names = TRUE), readText)
    parent_checks <- if (length(checks_env$rows)) do.call(rbind, lapply(checks_env$rows, function(x) as.data.frame(lapply(x, as.character), stringsAsFactors = FALSE)))
    checks <- do.call(rbind, c(list(parent_checks), unit_checks))
    utils::write.table(checks, file.path(out_dir, sprintf("checks-%s.tsv", job_name)), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    utils::write.table(timing, file.path(out_dir, sprintf("timing-units-%s.tsv", job_name)), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    files <- list.files(rows_dir, pattern = "^part-.*\\.tsv$", full.names = TRUE)
    table <- do.call(rbind, lapply(files, readText))
    table <- table[order(as.integer(table$t), match(table$design, DIAG_DESIGN_ORDER), match(table$norm, c(DIAG_NORMS, "unit"))), , drop = FALSE]
    utils::write.table(table, file.path(rows_dir, "diag-rows.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    unlink(files)
    failed <- checks[checks$pass != "TRUE", , drop = FALSE]
    for (name in unique(failed$check)) {
        x <- failed[failed$check == name, , drop = FALSE]
        cat(sprintf("K275DIAG %s: FAILED check %s (%d); first: %s\n", job_name, name, nrow(x), oneLine(paste(x$condition[[1L]], x$note[[1L]]), 160L)))
    }
    not_ok <- sum(table$status != "ok")
    jsonlite::write_json(provenance(TRUE), file.path(out_dir, sprintf("provenance-%s.json", job_name)), auto_unbox = TRUE, pretty = TRUE, digits = NA)
    cat(sprintf("K275DIAG %s finished: %d units (%d with an error), %d rows, %d not ok, %d failed check rows, wall %.1f min\n", job_name, nrow(units),
        sum(nzchar(timing$error)), nrow(table), not_ok, nrow(failed), as.numeric(difftime(Sys.time(), started_at, units = "mins"))))
    if (not_ok > 0L) quit(save = "no", status = 1L)   # an instrument row or a failed journey: the run's status says so (the scorer step still runs)
    invisible(NULL)
}

if (child) diagChildMain() else diagParentMain()
