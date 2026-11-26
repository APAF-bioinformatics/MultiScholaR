# KIT-281 lifecycle L1 (ticket section 7, VAL-11): the checks E1', E2 to E10, the scorer self-test (E4), E5, E6, G2, G3 on the scripts as
# declared, blind, once. Evidence code under docs/qa: no package code.
#
# Usage, from the package root of the bound commit on a run-tests runner (full git history):
#   Rscript docs/qa/workflow-execution/kit/kit-281/l1-k281.R <public out dir>
# Blind: the simulated rows and every score computed from them stay in a private work directory under tempdir(), which is not
# uploaded; <public out dir> receives only checks.tsv (one row per check: id, units compared, units passed, pass, mechanics details), timing.tsv
# (seconds per replicate-mask over all arms together, for G2, and E6's wall seconds), g2-plan.tsv, hashes.tsv, provenance-l1.json and progress-l1.tsv
# (written before the first unit and appended after each check), and the two tables of rule G3: g3-count.tsv and g3-rule.tsv. ONE count leaves
# the blind rows (ticket K281-P08; a declared departure from KIT-211's blind L1): per setting, the mean over L1's two replicates of the number of
# tested null proteins that only the 1-group rule adds (subset g_added of arm B1, column n_tested); nothing derived from a p-value, a call or
# a fold change is written. Nothing is printed from the rows or the scores; a failed child prints its check id and its condition's message.
#   E2, E3, E7, E10: sim-k281.R with K211_LIFECYCLE=1 on cells A to D, r = 9201 and r = 9202 as two ranges (rows/r9201 and rows/r9202,
#       each holding every cell), all masks, every arm.
#   E8: K211_DETERMINISM=1 on cell A r = 9201, 9202.
#   E1': S and B alone (K211_ARMS=S,B) on cells A to D, r = 201 to 205, against KIT-211's committed rows (filter-dependent columns exact).
#   E4: the scorer self-test (selftest/run-selftest.sh: constructed tables, predictions by selftest/make-selftest.py in Python, no simulated value).
#   E5: score-k281.R twice on the two ranges' rows (the join): every output file identical.
#   E6: cell A, r = 9201 to 9204 with one worker, killed (SIGKILL to the process tree) once one task is complete; continued from the list of
#       missing tasks; the assembled rows equal an uninterrupted run (four workers).
#   E11: a counted range (r = 401) is refused through the last net before load_all (no G3 table in the tree): exit 1, instrument-stop.txt, no metrics.
#   G0: the three protein filter files equal KIT-211's bound commit's (the replay marker of KIT-280 is asserted absent by every child).
#   G2: the planned job minutes (set-up plus script) of counted jobs of 200, 171 and 77 replicates from this run's seconds per replicate-mask, the
#       timeouts derived from them (g2-plan.tsv); none may pass 75 minutes.
#   G3: the count above, N-hat = 200 x its sum, R = rule G3 (the function g3Rule below, tested against Python's table selftest/g3-cases.tsv).
# E0 (the verbatim blocks) is local: e0-k281.py, l1/e0.txt.

options(warnPartialMatchDollar = TRUE) # script lesson 2; this script loads no package code
# The global handler is only the LAST NET: R runs a global calling handler with every local handler disabled, so its stop() reaches no tryCatch and ends the process at
# top level. The error option below therefore writes the record (the checks so far and the message) before the non-zero exit, so a guard stop is an instrument stop WITH a
# record, never a silent end. (L1 is sequential top-level code: its units are not wrapped in local handlers.)
options(error = function() {
    note <- geterrmessage()
    try({
        done <- if (exists("checks", envir = globalenv()) && is.list(checks) && !is.data.frame(checks) && length(checks)) do.call(rbind, checks) else NULL
        if (!is.null(done)) utils::write.table(done, file.path(public, "checks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
        writeLines(c("instrument stop (L1 ended at top level)", substr(note, 1L, 500L)), file.path(public, "l1-stop.txt"))
    }, silent = TRUE)
    quit(save = "no", status = 1L)
})
globalCallingHandlers(warning = function(w) {
    if (grepl("partial match of", tolower(conditionMessage(w)), fixed = TRUE)) stop(paste("partial-match guard:", conditionMessage(w)), call. = FALSE)
})

# ==== verbatim: l1-k197.R:34-37 at ef63c918 ====
run <- function(script, arguments, env = character(), log) {
    status <- system2(rscript, c(shQuote(script), shQuote(arguments)), env = env, stdout = log, stderr = log)
    status
}
# ==== end verbatim ====

# ==== verbatim: l1-k197.R:38-38 at ef63c918 ====
sha <- function(path) digest::digest(file = path, algo = "sha256")
# ==== end verbatim ====

# ---- KIT-281 L1 (not copied) ---------------------------------------------------------------------------------------------
KIT211_BOUND <- "894b9cc821a70b03baab1e1fa47b3d47440c7ed5" # the commit whose filter files G0 compares against
RANGES <- c(9201L, 9202L)
CELLS_L1 <- "A,B,C,D"
CELLS <- c("A", "B", "C", "D")
VD <- c("vd5", "vd10", "vd25")
G3_THRESHOLD <- 276220 # tested added null proteins needed at R = 200 for P(calibrated | 1.07) >= 0.80 (ticket section 8)
G3_CAP <- 571
G3_FIRST <- 401L
MAX_JOB_MINUTES <- 75
SETUP_MINUTES <- 10
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
public <- args[[1L]]
dir.create(public, recursive = TRUE, showWarnings = FALSE)
here <- "docs/qa/workflow-execution/kit/kit-281"
kit211_rows <- "docs/qa/workflow-execution/kit/kit-211/rows"
sim <- file.path(here, "sim-k281.R")
score <- file.path(here, "score-k281.R")
work <- file.path(tempdir(), "k281-l1")
dir.create(work, showWarnings = FALSE)
rscript <- file.path(R.home("bin"), "Rscript")
checks <- list()
progress_file <- file.path(public, "progress-l1.tsv")
writeLines("check\tunits\tpassed\tpass", progress_file)
bound_commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
# provenance before the first unit (script lesson 3): the scripts' and the declaration's hashes and the bound commit
declaration <- file.path(here, "declaration.md")
jsonlite::write_json(list(ticket = "KIT-281", status = "started", package_commit = bound_commit, kit211_bound_commit = KIT211_BOUND,
    script_sha256 = sha(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]))),
    sim_sha256 = sha(sim), score_sha256 = sha(score), declaration_sha256 = if (file.exists(declaration)) sha(declaration) else NA_character_,
    ranges = RANGES, r_version = R.version.string, started = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
    file.path(public, "provenance-l1.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
addCheck <- function(id, units, passed, detail) {
    checks[[length(checks) + 1L]] <<- data.frame(id = id, units = units, passed = passed, pass = isTRUE(units == passed && units > 0), detail = detail,
        stringsAsFactors = FALSE)
    cat(sprintf("%s\t%s\t%s\t%s\n", id, units, passed, isTRUE(units == passed && units > 0)), file = progress_file, append = TRUE)
}
num <- function(x) ifelse(is.na(x), "NA", sprintf("%.17g", x))
# Tables are read as text and converted by name; a value that is not a number stops L1 (script lesson 1).
readText <- function(path) utils::read.delim(if (grepl("\\.gz$", path)) gzfile(path) else path, colClasses = "character", na.strings = character(0),
    check.names = FALSE, stringsAsFactors = FALSE)
asNum <- function(x, what) {
    if (!all(grepl("^-?[0-9]+(\\.[0-9]+)?([eE][-+]?[0-9]+)?$", x))) stop(sprintf("%s: a value is not a number", what), call. = FALSE)
    as.numeric(x)
}
# Only the condition's class and message of a failed child, never a log tail (the messages carry no estimate).
failureNote <- function(id, log, errors) {
    note <- if (file.exists(errors)) readLines(errors, warn = FALSE)[1L] else grep("^Error", readLines(log, warn = FALSE),
        value = TRUE)[1L]
    cat(sprintf("%s failed: %s\n", id, substr(if (is.na(note)) "no message" else note, 1L, 200L)))
}
cat("KIT-281 L1: start", format(Sys.time(), "%H:%M:%S"), "\n")

# ---- G0: the three filter files equal KIT-211's bound commit's ---------------------------------------------------------------
g0 <- system2("git", c("diff", "--stat", KIT211_BOUND, "HEAD", "--", "R/func_prot_qc_filtering_helpers.R", "R/func_prot_s4_qc_filter_methods.R",
    "R/func_general_parameter_helpers.R"), stdout = TRUE)
addCheck(id = "G0", units = 1L, passed = if (length(g0) == 0L) 1L else 0L, detail = "git diff --stat of the three filter files between 894b9cc8 and the package commit is empty")

# ---- E11 (added after the script audit, AUDIT-KIT-281-scripts RC1-5): a counted range is refused through the last-net writer, before load_all -------------------
# sim-k281.R on r = 401 (a counted seed) with no g3-rule.tsv in the tree must exit 1, leave instrument-stop.txt ("a counted job needs the G3 table"), no metrics.tsv and
# no part file: no draw is made. It proves the G3 refusal and the error option's writer, which no other unit reaches.
e11_dir <- file.path(work, "e11")
status_e11 <- run(script = sim, arguments = c(e11_dir, "401", "401"), env = "K211_CELLS=A", log = file.path(work, "sim-e11.log"))
e11_stop <- file.path(e11_dir, "instrument-stop.txt")
e11_ok <- status_e11 == 1L && file.exists(e11_stop) && any(grepl("a counted job needs the G3 table", readLines(e11_stop, warn = FALSE), fixed = TRUE)) &&
    !file.exists(file.path(e11_dir, "metrics.tsv")) && length(list.files(file.path(e11_dir, "parts"))) == 0L
addCheck(id = "E11", units = 1L, passed = if (e11_ok) 1L else 0L, detail = sprintf("counted range r = 401 with no G3 table: exit %d (1 expected); instrument-stop.txt present %s; no metrics.tsv, no part file",
    status_e11, file.exists(e11_stop)))

# ---- the lifecycle simulation (E2, E3, E7, E10; two ranges per cell, joined by the scorer: E5) -----------------------------
rows_root <- file.path(work, "rows")
range_seconds <- numeric()
status_sim <- vapply(RANGES, \(r) {
    started <- Sys.time()
    status <- run(script = sim, arguments = c(file.path(rows_root, paste0("r", r)), r, r), env = c("K211_LIFECYCLE=1", paste0("K211_CELLS=", CELLS_L1)),
        log = file.path(work, sprintf("sim-%d.log", r)))
    range_seconds[[as.character(r)]] <<- as.numeric(difftime(Sys.time(), started, units = "secs"))
    if (status != 0L) failureNote(id = sprintf("simulation r = %d", r), log = file.path(work, sprintf("sim-%d.log", r)),
        errors = file.path(rows_root, paste0("r", r), "errors.txt"))
    status
}, numeric(1L))
sim_checks <- do.call(rbind, lapply(RANGES, \(r) {
    f <- file.path(rows_root, paste0("r", r), "checks-sim.tsv")
    if (file.exists(f)) utils::read.delim(f, stringsAsFactors = FALSE)
}))
for (id in c("E2", "E3", "E7", "E10")) {
    sub <- if (is.null(sim_checks)) NULL else sim_checks[sim_checks[["check"]] == id, , drop = FALSE]
    ok <- all(status_sim == 0) && !is.null(sub) && nrow(sub) == length(RANGES)
    addCheck(id = id, units = if (ok) sum(sub[["units"]]) else 0L, passed = if (ok) sum(sub[["passed"]]) else 0L,
        detail = if (!ok) sprintf("no result (simulation exits %s)", paste(status_sim, collapse = ", ")) else if (id == "E10")
            sprintf("replicate-masks over both ranges; maximum absolute deviation of logFC %s", num(max(sub[["max_abs_deviation"]]))) else
            "replicate-masks over both ranges")
}

# ---- E8: the determinism mechanics on cell A, r = 9201, 9202 -------------------------------------------------------------
e8_dir <- file.path(work, "e8")
status_e8 <- run(script = sim, arguments = c(e8_dir, "9201", "9202"), env = c("K211_CELLS=A", "K211_DETERMINISM=1"), log = file.path(work, "sim-e8.log"))
det <- if (file.exists(file.path(e8_dir, "determinism.tsv"))) readText(file.path(e8_dir, "determinism.tsv"))
if (status_e8 != 0L) failureNote(id = "E8 simulation", log = file.path(work, "sim-e8.log"), errors = file.path(e8_dir, "errors.txt"))
addCheck(id = "E8", units = if (is.null(det)) 0L else nrow(det), passed = if (is.null(det)) 0L else sum(det[["rows_identical"]] == "TRUE" & det[["metrics_identical"]] == "TRUE"),
    detail = "K211_DETERMINISM=1 on cell A r = 9201, 9202 (rows and metrics, seconds excluded)")

# ---- E1': S and B on KIT-211's r = 201 to 205 against its committed rows ----------------------------------------------------
e1_dir <- file.path(work, "e1")
status_e1 <- run(script = sim, arguments = c(e1_dir, "201", "205"), env = c("K211_CELLS=A,B,C,D", "K211_ARMS=S,B"), log = file.path(work, "sim-e1.log"))
if (status_e1 != 0L) failureNote(id = "E1' simulation", log = file.path(work, "sim-e1.log"), errors = file.path(e1_dir, "errors.txt"))
e1 <- if (status_e1 == 0L) {
    key_columns <- c("cell", "mask", "r", "arm", "subset", "class", "stratum", "ell", "band", "ell0")
    da_columns <- c("n_tested", "n_p05", "n_call", "n_bh", "sum_err", "sum_s_err", "sum_sgn_err")
    all_bad <- NULL
    da_equal <- TRUE
    for (cellname in CELLS) {
        old_dir <- file.path(kit211_rows, paste0(cellname, "-201-300"))
        mine <- readText(file.path(e1_dir, paste0(cellname, ".tsv.gz")))
        old <- readText(file.path(old_dir, paste0(cellname, ".tsv.gz")))
        old <- old[old[["r"]] %in% as.character(201:205) & old[["arm"]] %in% c("S", "B") & old[["mask"]] %in% VD, , drop = FALSE]
        mine <- mine[mine[["arm"]] %in% c("S", "B") & mine[["mask"]] %in% VD, , drop = FALSE]
        for (column in c("n", da_columns)) {
            mine[[column]] <- asNum(mine[[column]], paste("E1' mine", column))
            old[[column]] <- asNum(old[[column]], paste("KIT-211 rows", column))
        }
        key <- function(d) do.call(paste, c(lapply(d[key_columns], as.character), sep = "|"))
        km <- rowsum(mine[c("n", da_columns)], key(mine))
        ko <- rowsum(old[c("n", da_columns)], key(old))
        keys <- union(rownames(km), rownames(ko))
        bad <- is.na(match(keys, rownames(km))) | is.na(match(keys, rownames(ko)))
        both <- !bad
        bad[both] <- km[keys[both], "n"] != ko[keys[both], "n"]
        unit_of <- vapply(strsplit(keys, "|", fixed = TRUE), \(p) paste(p[1:4], collapse = "|"), character(1L))
        metrics_new <- readText(file.path(e1_dir, "metrics.tsv"))
        metrics_old <- readText(file.path(old_dir, "metrics.tsv"))
        metrics_new <- metrics_new[metrics_new[["cell"]] == cellname & metrics_new[["arm"]] %in% c("S", "B"), , drop = FALSE]
        metrics_old <- metrics_old[metrics_old[["cell"]] == cellname & metrics_old[["arm"]] %in% c("S", "B") & metrics_old[["r"]] %in% as.character(201:205) &
            metrics_old[["mask"]] %in% VD, , drop = FALSE]
        m <- merge(metrics_new[c("cell", "mask", "r", "arm", "kept", "added", "threshold")], metrics_old[c("cell", "mask", "r", "arm", "kept", "added", "threshold")],
            by = c("cell", "mask", "r", "arm"), suffixes = c(".new", ".old"), all = TRUE)
        # S04: a unit present on one side only is a recorded failure (not a crash of asNum on NA): compare the rows present on both sides
        both_sides <- !is.na(m[["kept.new"]]) & !is.na(m[["kept.old"]])
        unit_metric_bad <- !both_sides
        unit_metric_bad[both_sides] <- m[["kept.new"]][both_sides] != m[["kept.old"]][both_sides] | m[["added.new"]][both_sides] != m[["added.old"]][both_sides] |
            asNum(m[["threshold.new"]][both_sides], "threshold") != asNum(m[["threshold.old"]][both_sides], "threshold")
        names(unit_metric_bad) <- paste(m[["cell"]], m[["mask"]], m[["r"]], m[["arm"]], sep = "|")
        unit_bad <- c(tapply(bad, unit_of, any), unit_metric_bad)
        all_bad <- c(all_bad, tapply(unit_bad, names(unit_bad), any))
        shared <- intersect(rownames(ko), rownames(km))
        da_equal <- da_equal && isTRUE(all.equal(as.matrix(km[shared, da_columns]), as.matrix(ko[shared, da_columns]), tolerance = 1e-9, check.attributes = FALSE))
    }
    list(units = length(all_bad), passed = sum(!all_bad), da_equal = da_equal)
} else list(units = 0L, passed = 0L, da_equal = NA)
addCheck(id = "E1'", units = e1[["units"]], passed = e1[["passed"]], detail = sprintf(paste("S's and B's filter-dependent columns (membership counts per subset, class, stratum, ell, band and ell0; kept and added",
    "counts; threshold) against KIT-211's committed rows, cells A to D, r = 201 to 205, the three vd masks (120 replicate-masks expected: 4 cells x 3 masks",
    "x 5 replicates x 2 arms); G0 lists DA-file changes, so the DA-dependent columns are written as one flag: all equal %s"), e1[["da_equal"]]))

# ---- E4: the scorer self-test (constructed tables; predictions written by selftest/make-selftest.py in Python before the scorer ran) ----
st_dir <- file.path(here, "selftest")
st_out <- file.path(work, "selftest-out")
status_st <- system2("bash", c(shQuote(file.path(st_dir, "run-selftest.sh")), shQuote(normalizePath(st_dir)), shQuote(st_out)),
    env = paste0("PATH=", R.home("bin"), ":", Sys.getenv("PATH")), stdout = file.path(work, "selftest.log"), stderr = file.path(work, "selftest.log"))
st_log <- if (file.exists(file.path(work, "selftest.log"))) readLines(file.path(work, "selftest.log"), warn = FALSE) else character()
st_line <- grep("^SELFTEST: ", st_log, value = TRUE)
st_counts <- if (length(st_line)) as.integer(regmatches(st_line[[length(st_line)]], gregexpr("[0-9]+", st_line[[length(st_line)]]))[[1L]]) else integer()
addCheck(id = "E4", units = if (length(st_counts) >= 2L) st_counts[[1L]] else 0L, passed = if (length(st_counts) >= 2L) st_counts[[1L]] - st_counts[[2L]] else 0L,
    detail = if (length(st_line)) st_line[[length(st_line)]] else "no SELFTEST summary line (the runner failed before its check)")
if (status_st != 0L) cat("E4: the self-test reports differences (constructed tables only); differing lines:\n",
    paste(grep("DIFFERS", st_log, value = TRUE)[seq_len(min(20L, sum(grepl("DIFFERS", st_log))))], collapse = "\n"), "\n")

# ---- E5: the scoring twice on the lifecycle rows: every output file identical ------------------------------------------------------
score_env <- c(paste0("K281_RANGE=", RANGES[[1L]], ",", RANGES[[2L]]), "K211_B_BOOT=1000")
e5 <- lapply(c("e5a", "e5b"), \(name) run(script = score, arguments = c(rows_root, file.path(work, name)), env = score_env, log = file.path(work, paste0("score-", name, ".log"))))
files_a <- sort(list.files(file.path(work, "e5a")))
files_b <- sort(list.files(file.path(work, "e5b")))
same <- identical(files_a, files_b) && length(files_a) > 0L && all(vapply(files_a, \(f) identical(sha(file.path(work, "e5a", f)),
    sha(file.path(work, "e5b", f))), logical(1L)))
if (any(unlist(e5) != 0L)) failureNote(id = "E5 scoring", log = file.path(work, "score-e5a.log"), errors = "none")
addCheck(id = "E5", units = 2L, passed = if (all(unlist(e5) == 0L) && same) 2L else 0L, detail = sprintf("two scorer runs on the lifecycle rows (mechanics: K281_RANGE, B = 1000) exit %s; %d files each; all identical: %s",
    paste(unlist(e5), collapse = ", "), length(files_a), same))

# ---- E6: cell A, r = 9201 to 9204, killed once one task is complete, continued; equal to an uninterrupted run ----------------------
e6_cut <- file.path(work, "e6-cut")
e6_full <- file.path(work, "e6-full")
dir.create(file.path(e6_cut, "parts"), recursive = TRUE, showWarnings = FALSE)
all_tasks <- file.path(work, "e6-tasks.txt")
writeLines(sprintf("A %d", 9201:9204), all_tasks)
done <- function() sum(file.exists(file.path(e6_cut, "parts", sprintf("A-%04d.metrics.tsv", 9201:9204))))
child <- processx::process$new(rscript, c(sim, e6_cut, "9201", "9204", all_tasks),
    env = c(Sys.getenv(), K211_WORKERS = "1", K211_CELLS = "A"), stdout = file.path(work, "e6-cut.log"), stderr = "2>&1",
    cleanup_tree = TRUE)
deadline <- Sys.time() + 1800
while (child$is_alive() && done() < 1L && Sys.time() < deadline) Sys.sleep(0.2)
killed <- child$is_alive()
invisible(child$kill_tree())
Sys.sleep(5)
complete_after_kill <- done()
missing <- setdiff(9201:9204, 9200L + which(file.exists(file.path(e6_cut, "parts", sprintf("A-%04d.metrics.tsv", 9201:9204)))))
missing_file <- file.path(work, "e6-missing.txt")
writeLines(sprintf("A %d", missing), missing_file)
status_cont <- run(script = sim, arguments = c(e6_cut, "9201", "9204", missing_file), env = "K211_CELLS=A", log = file.path(work, "e6-cont.log"))
started <- Sys.time()
status_full <- run(script = sim, arguments = c(e6_full, "9201", "9204"), env = "K211_CELLS=A", log = file.path(work, "e6-full.log"))
full_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
readGz <- function(dir) if (file.exists(file.path(dir, "A.tsv.gz"))) readLines(file.path(dir, "A.tsv.gz")) else NULL
readMetrics <- function(dir) if (file.exists(file.path(dir, "metrics.tsv"))) {
    m <- readText(file.path(dir, "metrics.tsv"))
    m[setdiff(names(m), "seconds")]
} else NULL
rows_cut <- readGz(e6_cut)
rows_full <- readGz(e6_full)
replicate_counts <- if (is.null(rows_cut)) integer() else {
    cut_rows <- readText(file.path(e6_cut, "A.tsv.gz"))
    table(unique(cut_rows[c("r", "mask")])[["r"]])
}
e6_pass <- killed && complete_after_kill >= 1L && complete_after_kill < 4L && status_cont == 0L && status_full == 0L &&
    !is.null(rows_cut) && identical(rows_cut, rows_full) && identical(readMetrics(e6_cut), readMetrics(e6_full)) &&
    identical(as.integer(replicate_counts), rep(3L, 4L))
addCheck(id = "E6", units = 1L, passed = if (e6_pass) 1L else 0L, detail = sprintf(paste("killed while running: %s; complete tasks after the kill: %d of 4; continued",
    "with %d task(s), exit %d; uninterrupted exit %d; assembled rows identical: %s; metrics identical (seconds",
    "excluded): %s; each replicate once with three masks: %s"), killed, complete_after_kill, length(missing), status_cont,
    status_full, identical(rows_cut, rows_full), identical(readMetrics(e6_cut), readMetrics(e6_full)),
    identical(as.integer(replicate_counts), rep(3L, 4L))))

# ---- G2: the planned minutes of a counted job (200 replicates x 3 masks = 600 replicate-masks) from this run's seconds -------------
seconds_all <- do.call(rbind, lapply(RANGES, \(r) {
    f <- file.path(rows_root, paste0("r", r), "metrics.tsv")
    if (file.exists(f)) readText(f)
}))
s_all <- if (is.null(seconds_all)) NA_real_ else asNum(seconds_all[["seconds"]][seconds_all[["arm"]] == "S"], "seconds")
# S05: G2 plans the counted jobs from this run's seconds: planned script minutes (n replicates x 3 masks x mean seconds / 4 workers), planned job minutes
# (+ 10 minutes of set-up) and timeout_minutes = min(75, 5 x ceiling(1.5 x job / 5)); it passes when every planned job is at most 75 minutes. L1's seconds
# include the lifecycle checks, so the plan errs long (safe). The counted requests use these timeouts and the running ceiling is recomputed from them.
plan_n <- c(200L, 171L, 77L)
plan_script <- plan_n * 3 * mean(s_all) / 4 / 60
plan_job <- plan_script + SETUP_MINUTES
plan_timeout <- pmin(MAX_JOB_MINUTES, 5 * ceiling(1.5 * plan_job / 5))
utils::write.table(data.frame(replicates = plan_n, script_minutes = num(plan_script), job_minutes = num(plan_job), timeout_minutes = plan_timeout,
    stringsAsFactors = FALSE), file.path(public, "g2-plan.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
addCheck(id = "G2", units = length(plan_n), passed = sum(plan_job <= MAX_JOB_MINUTES, na.rm = TRUE), detail = sprintf(paste("planned job minutes (script + %d set-up) at this run's",
    "seconds per replicate-mask for 200, 171 and 77 replicates: %s; timeouts %s; none may pass %d minutes (g2-plan.tsv)"), SETUP_MINUTES,
    paste(sprintf("%.1f", plan_job), collapse = ", "), paste(plan_timeout, collapse = ", "), MAX_JOB_MINUTES))

# ---- G3: the one count, rule G3 as code, its table ------------------------------------------------------------------------------------
# Rule G3 (ticket section 8): N-hat = 200 x (m_1 + ... + m_12); R = 200 when N-hat >= 276,220, else min(571, max(200, ceiling(200 x 276,220 / N-hat))).
g3Rule <- function(n_hat) {
    if (!is.finite(n_hat) || n_hat < 0) stop("rule G3 needs a finite non-negative N-hat", call. = FALSE)
    if (n_hat >= G3_THRESHOLD) 200L else if (n_hat == 0) as.integer(G3_CAP) else as.integer(min(G3_CAP, max(200, ceiling(200 * G3_THRESHOLD / n_hat))))
}
g3_cases <- file.path(here, "selftest", "g3-cases.tsv")
g3_expect <- readText(g3_cases)
g3_got <- vapply(asNum(g3_expect[["n_hat"]], "g3-cases n_hat"), g3Rule, integer(1L))
addCheck(id = "G3rule", units = nrow(g3_expect), passed = sum(g3_got == asNum(g3_expect[["r"]], "g3-cases r")), detail = "rule G3's function against Python's table (Fractions): the examples of the ticket and both edges of every branch")
g3_count <- do.call(rbind, lapply(CELLS, \(cellname) {
    totals <- stats::setNames(numeric(length(VD)), VD)
    for (r in RANGES) {
        d <- readText(file.path(rows_root, paste0("r", r), paste0(cellname, ".tsv.gz")))
        select <- d[["arm"]] == "B1" & d[["subset"]] == "g_added" & d[["class"]] == "null"
        for (mask in VD) totals[[mask]] <- totals[[mask]] + sum(asNum(d[["n_tested"]][select & d[["mask"]] == mask], paste("g3 count", cellname, mask)))
    }
    data.frame(cell = cellname, mask = VD, m_s = unname(totals) / length(RANGES), stringsAsFactors = FALSE)
}))
n_hat <- 200 * sum(g3_count[["m_s"]])
r_cell <- g3Rule(n_hat)
branch <- if (n_hat >= G3_THRESHOLD) "N-hat at or above 276220: R = 200" else if (r_cell == G3_CAP) "R = 571 (the cap of the rule)" else "R = ceiling(200 x 276220 / N-hat)"
g3_rule <- data.frame(cell = CELLS, r_cell = r_cell, first_r = G3_FIRST, last_r = G3_FIRST + r_cell - 1L, n_hat = n_hat, branch = branch, stringsAsFactors = FALSE)
utils::write.table(data.frame(cell = g3_count[["cell"]], mask = g3_count[["mask"]], m_s = num(g3_count[["m_s"]]), stringsAsFactors = FALSE),
    file.path(public, "g3-count.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(data.frame(cell = g3_rule[["cell"]], r_cell = g3_rule[["r_cell"]], first_r = g3_rule[["first_r"]], last_r = g3_rule[["last_r"]],
    n_hat = num(g3_rule[["n_hat"]]), branch = g3_rule[["branch"]], stringsAsFactors = FALSE), file.path(public, "g3-rule.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
addCheck(id = "G3", units = 1L, passed = as.integer(is.finite(n_hat) && r_cell >= 200L && r_cell <= G3_CAP), detail = sprintf("one count over 12 settings written (g3-count.tsv); rule G3's table written (g3-rule.tsv); R and seeds r = %d to %d are in the table, not here",
    G3_FIRST, G3_FIRST + r_cell - 1L))

# ---- public records ------------------------------------------------------------------------------------------------------
checks <- do.call(rbind, checks)
utils::write.table(checks, file.path(public, "checks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
timing <- data.frame(run = c("L1 cells A to D (4 workers; per replicate-mask, all arms together)", "E6 uninterrupted (4 workers, 4 tasks, wall seconds)"),
    replicate_masks = c(length(s_all), NA), mean_seconds = c(mean(s_all), full_seconds), max_seconds = c(max(s_all), NA), stringsAsFactors = FALSE)
utils::write.table(timing, file.path(public, "timing.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
hash_files <- c(file.path(here, c("sim-k281.R", "score-k281.R", "l1-k281.R")), list.files(rows_root, "\\.(tsv|gz)$", recursive = TRUE,
    full.names = TRUE), file.path(e6_full, "A.tsv.gz"))
hashes <- data.frame(file = c("bound commit", "package commit", sub(paste0(work, "/"), "", hash_files, fixed = TRUE)),
    sha256 = c(KIT211_BOUND, bound_commit, vapply(hash_files, \(f) if (file.exists(f)) sha(f) else NA_character_, character(1L))), stringsAsFactors = FALSE)
utils::write.table(hashes, file.path(public, "hashes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("KIT-281 L1: checks\n")
print(checks[c("id", "units", "passed", "pass")], row.names = FALSE)
cat("KIT-281 L1: end", format(Sys.time(), "%H:%M:%S"), "\n")
if (!all(checks[["pass"]])) quit(status = 1L)
