# KIT-144: the per-run results of the KIT-124 main scenarios as one compact, committed table.
#
# The simulators' run-level outputs are not committed (they sit in the simulation worktree and are still appended to
# by the extra-hard scenarios). figures.R needs them for the figures that show more than a mean, so this script copies
# what those figures use to kit-144/data/:
#   k124-main-runs.tsv.xz     one row per omic (simulator), main scenario (cell C0..C7), simulated data set (run) and
#                             design (arm), on the design's primary feature subset
#   k124-main-runs.meta.txt   sidecar: script hash, commit, command, each input file's sha256 and size, what was
#                             dropped, the reproduction check, the output's sha256
#
# It reads exactly the files that scored the main scenarios: the inputs named in
# kit-124/scored/s<i>-core/k124-aggregate.meta.txt, each checked against the sha256 recorded there. A table the
# simulator has appended to since (S1 writes every cell to one file) is read up to its scored length, and those bytes
# are checked. Statuses, the rows that are scored and the false discovery proportion of a run are those of
# kit-124/aggregate.R (loadSim, classify, cellTable): the simulator's own `fdp` (false calls / calls, 0 without calls,
# so 1 under the global null when anything is called; one value over both assays in lipidomics/metabolomics), a
# refused run carrying its data set's no-RUV result. Before anything is written the per-cell means are compared with
# scored/s<i>-core/k124-cells.tsv at that table's 6 significant digits, for every design, subset and main scenario.
#
# Usage, from the repository root (one process, under a minute):
#   R_LIBS=<simulation libraries> nice -n 19 Rscript docs/qa/workflow-execution/kit/kit-144/extract-runs.R \
#       <directory that holds s1-confirm .. s4-confirm>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !dir.exists(args[[1L]])) {
    cli::cli_abort(
        "usage: Rscript extract-runs.R <directory that holds s1-confirm .. s4-confirm>",
        class = "kit144_usage"
    )
}
runs_root <- normalizePath(args[[1L]])
scored_dir <- here::here("docs", "qa", "workflow-execution", "kit", "kit-124", "scored")
data_dir <- here::here("docs", "qa", "workflow-execution", "kit", "kit-144", "data")
simulators <- c(S1 = "s1-core", S2 = "s2-core", S3 = "s3-core", S4 = "s4-core")

# As kit-124/aggregate.R.
MEASURES <- c("k", "n_calls", "fdp", "power")
COMPARATOR_ARMS <- c("12", "13", "14", "14E")
FDP_TAIL <- 0.10
KEPT <- c("simulator", "cell", "run", "arm", "subset", "class", MEASURES)

# ---- Inputs: the files the scoring read, verified -----------------------------------------------------------------

# The run tables and message tables among the inputs of one scoring (its sidecar's "inputs" lines).
scoredInputs <- function(folder) {
    path <- file.path(scored_dir, folder, "k124-aggregate.meta.txt")
    lines <- grep("^  [0-9a-f]{64} [0-9]+ /", readLines(path), value = TRUE)
    parts <- strsplit(trimws(lines), " ", fixed = TRUE)
    inputs <- data.frame(
        sha256 = vapply(parts, `[[`, "", 1L),
        bytes = as.numeric(vapply(parts, `[[`, "", 2L)),
        file = sub("^.*/(s[1-4]-confirm/)", "\\1", vapply(parts, \(p) paste(p[-(1:2)], collapse = " "), ""))
    )
    inputs <- inputs[!grepl("[.](meta[.]txt|timing[.]tsv)$", inputs$file), ]
    if (!nrow(inputs)) {
        cli::cli_abort("{.path {path}} names no run table", class = "kit144_inconsistent")
    }
    inputs
}

# A path holding exactly the scored bytes of one input: the file itself, or its first `bytes` bytes in a temporary
# file when the simulator has appended to it since. Stops when those bytes are not the ones that were scored.
scoredBytes <- function(input) {
    path <- file.path(runs_root, input$file)
    if (!file.exists(path)) {
        cli::cli_abort("run-level file not found: {.path {path}}", class = "kit144_missing_file")
    }
    on_disk <- file.size(path)
    read_from <- path
    if (on_disk > input$bytes) {
        read_from <- tempfile()
        con <- file(path, "rb")
        on.exit(close(con))
        writeBin(readBin(con, "raw", input$bytes), read_from)
    }
    if (on_disk < input$bytes || unname(tools::sha256sum(read_from)) != input$sha256) {
        cli::cli_abort(
            "{.path {path}} is not the file that was scored (sha256 of its first {input$bytes} bytes differs)",
            class = "kit144_inconsistent"
        )
    }
    data.frame(input, bytes_on_disk = on_disk, read_from = read_from)
}

# ---- Reading (kit-124/aggregate.R: loadSim) -----------------------------------------------------------------------

readTsv <- function(path, ...) {
    utils::read.delim(path, colClasses = "character", quote = "", comment.char = "", fill = TRUE,
        na.strings = "NA", ...)
}

bindRows <- function(parts) {
    parts <- Filter(\(p) !is.null(p) && nrow(p) > 0L, parts)
    if (!length(parts)) {
        return(NULL)
    }
    do.call(rbind, lapply(parts, \(p) {
        p[] <- lapply(p, as.character)
        p
    }))
}

# The run rows and the message table (cell, run, arm, subset, status, message) of one simulator.
loadSimulator <- function(sim, inputs) {
    paths <- stats::setNames(inputs$read_from, basename(inputs$file))
    pick <- function(pattern) unname(paths[grepl(pattern, names(paths))])
    messages <- NULL
    if (sim == "S1") {
        runs <- readTsv(pick("^s1-runs[.]tsv$"))
        log <- readTsv(pick("^s1-runs[.]log$"), header = FALSE,
            col.names = c("cell", "run", "arm", "status", "message"))
        messages <- data.frame(log[c("cell", "run", "arm")], subset = NA_character_, log[c("status", "message")])
    } else if (sim == "S2") {
        parts <- lapply(pick("[.]rds$"), readRDS)
        runs <- bindRows(lapply(parts, `[[`, "runs"))
        errors <- bindRows(lapply(parts, `[[`, "errors"))
        if (!is.null(errors)) {
            messages <- data.frame(errors[c("cell", "run", "arm")], subset = NA_character_,
                errors[c("status", "message")])
        }
    } else if (sim == "S3") {
        parts <- lapply(pick("[.]rds$"), readRDS)
        runs <- bindRows(lapply(parts, `[[`, "rows"))
        events <- bindRows(lapply(parts, `[[`, "events"))
        if (!is.null(events)) {
            messages <- data.frame(events[c("cell", "run", "arm")], subset = NA_character_,
                events[c("status", "message")])
        }
    } else {
        runs <- bindRows(lapply(pick("^C[0-9]+[.]tsv$"), readTsv))
        notes <- bindRows(lapply(pick("[.]notes[.]tsv$"), readTsv))
        if (!is.null(notes)) {
            messages <- data.frame(notes[c("cell", "run", "arm", "subset", "status")], message = notes$reason)
        }
    }
    if (is.null(messages)) {
        messages <- data.frame(cell = character(), run = character(), arm = character(), subset = character(),
            status = character(), message = character())
    }
    complete_line <- runs$status %in% c("ok", "refused") & !is.na(runs$seconds) | runs$status %in% "crashed"
    runs <- runs[complete_line & !is.na(runs$subset) & nzchar(runs$subset), ]
    for (m in MEASURES) {
        runs[[m]] <- suppressWarnings(as.numeric(runs[[m]]))
    }
    list(runs = runs, messages = messages)
}

# ---- Statuses (kit-124/aggregate.R: classify) ---------------------------------------------------------------------

# ok, designed refusal, method refusal (a comparator's own stop, RUVinv's eigenvalue stop among them) or crash. A
# reclassified RUVinv row takes the no-RUV measures of its data set; S2's design 16 repeats design 7 and is dropped.
classify <- function(runs, messages, sim) {
    if (sim == "S2") {
        runs <- runs[runs$arm != "16", ]
    }
    class <- ifelse(runs$status == "ok", "ok", ifelse(runs$status == "refused",
        ifelse(runs$arm %in% COMPARATOR_ARMS, "refused_method", "refused_designed"), "crashed"))
    crashes <- messages[messages$status == "crashed", ]
    ruvinv <- grepl("fourth-smallest eigenvalue", crashes$message, fixed = TRUE) & crashes$arm == "13"
    keyOf <- function(x) paste(x$cell, x$run, x$arm, if (sim == "S4") x$subset else rep("", nrow(x)), sep = "|")
    all_ruvinv <- tapply(ruvinv, keyOf(crashes), all)
    hit <- class == "crashed" & runs$arm == "13" & keyOf(runs) %in% names(all_ruvinv)[all_ruvinv]
    if (any(hit)) {
        reference <- runs[runs$arm == "1", ]
        at <- match(paste(runs$cell[hit], runs$run[hit], runs$subset[hit]),
            paste(reference$cell, reference$run, reference$subset))
        for (m in setdiff(MEASURES, "k")) {
            runs[[m]][hit] <- reference[[m]][at]
        }
        class[hit] <- "refused_method"
    }
    runs$class <- class
    runs
}

# ---- The check against the scored table (kit-124/aggregate.R: cellTable) ------------------------------------------

meanFinite <- function(v) if (any(is.finite(v))) mean(v[is.finite(v)]) else NA_real_

# One row per cell, design, subset and measure: this script's value beside the scored table's.
compareWithScored <- function(sim, folder, runs) {
    scored_rows <- runs[runs$class != "crashed" & !is.na(runs$fdp), ]
    groups <- split(scored_rows, paste(scored_rows$cell, scored_rows$arm, scored_rows$subset, sep = "\t"))
    mine <- do.call(rbind, lapply(groups, \(r) {
        data.frame(cell = r$cell[[1L]], arm = r$arm[[1L]], subset = r$subset[[1L]], scored = nrow(r),
            mean_fdp = meanFinite(r$fdp), power = meanFinite(r$power), mean_k = mean(r$k, na.rm = TRUE),
            fdp10_share = mean(r$fdp > FDP_TAIL))
    }))
    table <- readTsv(file.path(scored_dir, folder, "k124-cells.tsv"))
    at <- match(paste(table$cell, table$arm, table$subset), paste(mine$cell, mine$arm, mine$subset))
    compared <- c("scored", "mean_fdp", "power", "mean_k", "fdp10_share")
    do.call(rbind, lapply(compared, \(measure) {
        theirs <- suppressWarnings(as.numeric(table[[measure]]))
        ours <- signif(mine[[measure]][at], 6)
        data.frame(simulator = sim, table[c("cell", "arm", "subset")], measure = measure, scored = theirs,
            extracted = ours,
            agree = (is.na(ours) & is.na(theirs)) | (!is.na(ours) & !is.na(theirs) &
                abs(ours - theirs) <= 1e-9 * pmax(abs(theirs), 1e-12)))
    }))
}

# ---- Run ----------------------------------------------------------------------------------------------------------

inputs <- NULL
kept <- NULL
checks <- NULL
dropped <- NULL
for (sim in names(simulators)) {
    sim_inputs <- scoredInputs(simulators[[sim]])
    sim_inputs <- do.call(rbind, lapply(seq_len(nrow(sim_inputs)), \(i) scoredBytes(sim_inputs[i, ])))
    loaded <- loadSimulator(sim, sim_inputs)
    runs <- classify(loaded$runs, loaded$messages, sim)
    checks <- rbind(checks, compareWithScored(sim, simulators[[sim]], runs))
    arms <- readTsv(file.path(scored_dir, simulators[[sim]], "k124-arms.tsv"))
    primary <- paste(runs$arm, runs$subset) %in% paste(arms$arm, arms$subset)[arms$primary == "TRUE"]
    dropped <- rbind(dropped, data.frame(simulator = sim, rows_read = nrow(runs), rows_kept = sum(primary)))
    kept <- rbind(kept, runs[primary, KEPT])
    inputs <- rbind(inputs, data.frame(simulator = sim, sim_inputs))
    cli::cli_inform("{sim}: {nrow(runs)} rows read, {sum(primary)} kept (primary subset)")
}

disagree <- checks[!checks$agree, ]
if (nrow(disagree)) {
    print(disagree, row.names = FALSE)
    cli::cli_abort(
        "{nrow(disagree)} of {nrow(checks)} per-cell values differ from the scored tables; nothing written",
        class = "kit144_inconsistent"
    )
}

arm_order <- c(as.character(1:17), sort(setdiff(unique(kept$arm), as.character(1:17))))
kept$run <- as.integer(kept$run)
kept <- kept[order(kept$simulator, kept$cell, match(kept$arm, arm_order), kept$run), ]
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
table_path <- file.path(data_dir, "k124-main-runs.tsv.xz")
con <- xzfile(table_path, "wb", compression = 9)
utils::write.table(kept, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
close(con)

script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
git <- function(...) {
    tryCatch(
        system2("git", c("-C", shQuote(here::here()), ...), stdout = TRUE, stderr = FALSE),
        error = \(e) NA_character_, warning = \(w) NA_character_
    )
}
cells_checked <- nrow(unique(checks[c("simulator", "cell", "arm", "subset")]))
writeLines(c(
    "KIT-144: per-run results of the KIT-124 main scenarios (C0 to C7), for figures S6 and S7",
    paste("script:", "docs/qa/workflow-execution/kit/kit-144/extract-runs.R"),
    paste("script sha256:", unname(tools::sha256sum(script))),
    paste("commit:", paste(git("rev-parse", "HEAD"), collapse = " ")),
    paste("command: Rscript", script, runs_root),
    paste("R:", R.version.string),
    paste("run at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    paste("output: k124-main-runs.tsv.xz,", nrow(kept), "rows, sha256", unname(tools::sha256sum(table_path))),
    "columns: simulator; cell; run (the data set, 1 to 1000; the seed is 100000 + 10000 x cell index + run);",
    "  arm (the design, KIT-124's numbering); subset (the design's primary feature subset, kit-124's k124-arms.tsv);",
    "  class (ok, refused_designed, refused_method, crashed, as kit-124/aggregate.R classifies);",
    "  k (unwanted factors removed; the mean over the two assays in S4; NA where the design has none);",
    "  n_calls (features called at q < 0.05); fdp (false calls / calls, 0 without calls; one value over both assays",
    "  in S4); power (called true effects / true effects; NA without true effects).",
    "  A refused run carries the no-RUV result of its data set (KIT-124 R-07), as scored.",
    "",
    "dropped:",
    "  extra-hard scenarios (S1 to S14): not read (incomplete for some simulators)",
    sprintf("  rows of a design's non-primary feature subset: %s", paste(sprintf("%s %d of %d",
        dropped$simulator, dropped$rows_read - dropped$rows_kept, dropped$rows_read), collapse = "; ")),
    "  columns: cell_class, seed, status (in class), n_features, n_true, null_p05, null_p01, bias, rmse, seconds",
    "",
    sprintf("reproduction check against kit-124/scored/s<i>-core/k124-cells.tsv, before the subset was dropped: %d",
        cells_checked),
    "  rows (simulator x cell x design x subset), each on scored (count), mean_fdp, power, mean_k and fdp10_share at",
    sprintf("  6 significant digits: %d of %d values agree", sum(checks$agree), nrow(checks)),
    "",
    "inputs (sha256 of the bytes read = the sha256 recorded by the scoring; bytes read; bytes on disk; file):",
    sprintf("  %s %.0f %.0f %s", inputs$sha256, inputs$bytes, inputs$bytes_on_disk, inputs$file)
), file.path(data_dir, "k124-main-runs.meta.txt"))
cli::cli_inform(c(
    "reproduction check: {sum(checks$agree)} of {nrow(checks)} values agree ({cells_checked} scored rows)",
    "written: {.path {table_path}} ({nrow(kept)} rows, {file.size(table_path)} bytes)"
))
