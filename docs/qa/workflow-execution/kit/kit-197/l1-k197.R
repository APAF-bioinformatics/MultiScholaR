# KIT-197 lifecycle L1 (ticket 7.3, VAL-11): the checks E1 to E6 on the scripts as declared, blind. Evidence code under
# docs/qa: no package code.
#
# Usage, from the package root of the bound commit on a run-tests runner (full git history):
#   Rscript docs/qa/workflow-execution/kit/kit-197/l1-k197.R <public out dir>
# Environment: K197_BASE=<40-digit commit> (S2's base, for E1).
# Blind: the simulated rows and every score computed from them stay in a private work directory under tempdir(),
# which is not uploaded; <public out dir> receives only checks.tsv (E1 to E6: pass or fail and mechanics details),
# timing.tsv (seconds per replicate-mask, for G2) and hashes.tsv (SHA-256 of the row files). Nothing is printed from
# the rows or the scores.
#   E1, E2, E3: sim-k197.R with K197_LIFECYCLE=1 on cells A to D, r = 9001, 9002, all six masks.
#   E4: score-k197.R deciding on a constructed table with known sums (r = 1 to 3, every cell and mask alike), against
#       the hand values (exact fractions, e4hand.py in the records) to a relative 1e-12, and the pooled label "bias",
#       outcome O1.
#   E5: score-k197.R deciding twice on the E1 rows: every output file identical.
#   E6: cell A, r = 9001 to 9004 with one worker, killed (SIGKILL to the process tree) once one task is complete;
#       continued from the list of missing tasks; the assembled rows equal an uninterrupted run (four workers).

args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
public <- args[[1L]]
dir.create(public, recursive = TRUE, showWarnings = FALSE)
here <- "docs/qa/workflow-execution/kit/kit-197"
sim <- file.path(here, "sim-k197.R")
score <- file.path(here, "score-k197.R")
work <- file.path(tempdir(), "k197-l1")
dir.create(work, showWarnings = FALSE)
rscript <- file.path(R.home("bin"), "Rscript")
base <- Sys.getenv("K197_BASE")
stopifnot(grepl("^[0-9a-f]{40}$", base))
checks <- list()
addCheck <- function(id, pass, detail) checks[[length(checks) + 1L]] <<- data.frame(id = id, pass = isTRUE(pass),
    detail = detail, stringsAsFactors = FALSE)
run <- function(script, arguments, env = character(), log) {
    status <- system2(rscript, c(shQuote(script), shQuote(arguments)), env = env, stdout = log, stderr = log)
    status
}
sha <- function(path) digest::digest(file = path, algo = "sha256")
cat("KIT-197 L1: start", format(Sys.time(), "%H:%M:%S"), "\n")

# ---- E1, E2, E3 ----------------------------------------------------------------------------------------------------
l1_rows <- file.path(work, "l1")
started <- Sys.time()
status <- run(sim, c(l1_rows, "9001", "9002"), c("K197_LIFECYCLE=1", "K197_DETERMINISM=1", paste0("K197_BASE=", base)),
    file.path(work, "sim-l1.log"))
l1_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
cat("KIT-197 L1: lifecycle simulation exit", status, "after", round(l1_seconds), "s\n")
sim_checks <- if (file.exists(file.path(l1_rows, "checks-sim.tsv"))) utils::read.delim(file.path(l1_rows,
    "checks-sim.tsv"), stringsAsFactors = FALSE) else NULL
for (id in c("E1", "E2a", "E2b", "E3")) {
    row <- if (is.null(sim_checks)) NULL else sim_checks[sim_checks$check == id, , drop = FALSE]
    addCheck(id, status == 0L && !is.null(row) && nrow(row) == 1L && row$pass,
        if (is.null(row) || !nrow(row)) sprintf("no result (simulation exit %d)", status) else
            sprintf("%d of %d units pass; %s", row$passed, row$units, row$detail))
}
if (status != 0L) cat(utils::tail(readLines(file.path(work, "sim-l1.log"), warn = FALSE), 30L), sep = "\n")
det <- if (file.exists(file.path(l1_rows, "determinism.tsv"))) utils::read.delim(file.path(l1_rows, "determinism.tsv"))
addCheck("DET", !is.null(det) && nrow(det) == 2L && all(det$rows_identical) && all(det$metrics_identical),
    if (is.null(det)) "no determinism.tsv" else sprintf("determinism pass (K197_DETERMINISM=1) on cell A r = 9001, 9002: %d of %d identical",
        sum(det$rows_identical & det$metrics_identical), nrow(det)))

# ---- E4 ------------------------------------------------------------------------------------------------------------
e4_rows <- file.path(work, "e4")
dir.create(e4_rows, showWarnings = FALSE)
reps <- data.frame(r = 1:3, a = c(8, 10, 9), v = c(2, 3, 1), w = c(1, 2, 1), m = c(1.0, -0.5, 0.2), b = c(3, 2, 4),
    T = c(30, 34, 32), o = c(50, 52, 48), VO = c(4, 5, 3), TO = c(200, 210, 190), t = c(0.002, 0.003, 0.0025))
row <- function(subset, class, stratum, ell, n, n_tested, n_p05, n_p01, n_call, n_call_away, sum_err, sum_s_err,
                sum_sgn_err, sum_sq_err, sum_resvar, n_resvar, n_s) data.frame(subset = subset, class = class,
    stratum = stratum, ell = ell, n = n, n_tested = n_tested, n_p05 = n_p05, n_p01 = n_p01, n_call = n_call,
    n_call_away = n_call_away, sum_err = sum_err, sum_s_err = sum_s_err, sum_sgn_err = sum_sgn_err,
    sum_sq_err = sum_sq_err, sum_resvar = sum_resvar, n_resvar = n_resvar, n_s = n_s, stringsAsFactors = FALSE)
constructed <- function(x) {
    added0 <- function(s) row(s, "null", "one_treated", 0L, 100, 100, x$a, 1, x$v, x$w, x$m, x$m, 0, 0, 10, 100, 100)
    added1 <- function(s) row(s, "null", "one_treated", 1L, 20, 20, x$b, 0, 1, 1, 0.4, 0.4, 0, 0, 2, 20, 20)
    old_null <- function(s) row(s, "null", "none", 0L, 1000, 1000, x$o, 10, x$VO, 0, 0.3, 0, 0, 2, 80, 1000, 0)
    rbind(added0("added"), added1("added"), row("added", "up", "none", 0L, 50, 50, 40, 30, x$T, 0, -5, 0, -5, 1, 5, 50, 0),
        old_null("kept_both"), row("kept_both", "up", "none", 0L, 250, 250, 230, 220, x$TO, 0, -10, 0, -10, 3, 20, 250, 0),
        old_null("ref"), row("ref", "up", "none", 0L, 250, 250, 230, 220, x$TO, 0, -10, 0, -10, 3, 20, 250, 0),
        added0("whole"), added1("whole"), old_null("whole"),
        row("whole", "up", "none", 0L, 300, 300, 270, 250, x$T + x$TO, 0, -15, 0, -15, 4, 25, 300, 0))
}
masks <- c("vd5", "vd10", "vd25", "mcar5", "mcar10", "mcar25")
metrics <- list()
for (cellname in c("A", "B", "C", "D")) {
    rows <- list()
    for (mask in masks) for (i in seq_len(nrow(reps))) {
        rows[[length(rows) + 1L]] <- cbind(data.frame(cell = cellname, mask = mask, r = reps$r[[i]],
            stringsAsFactors = FALSE), constructed(reps[i, ]))
        metrics[[length(metrics) + 1L]] <- data.frame(cell = cellname, mask = mask, r = reps$r[[i]], kept_h = 0L,
            kept_o = 0L, added = 0L, calls = 0L, power = 0, fdp = 0, added_share_called = NA, added_fdp = NA,
            qmethod = "qvalue", f1_o_not_h = 0L, t_r = reps$t[[i]], mono_ok = TRUE, threshold = 0, seconds = 0,
            stringsAsFactors = FALSE)
    }
    connection <- gzfile(file.path(e4_rows, paste0(cellname, ".tsv.gz")), "w")
    utils::write.table(do.call(rbind, rows), connection, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    close(connection)
}
utils::write.table(do.call(rbind, metrics), file.path(e4_rows, "s2metrics.tsv"), sep = "\t", quote = FALSE,
    row.names = FALSE, na = "NA")
e4_out <- file.path(work, "e4-out")
status <- run(score, c("deciding", e4_rows, e4_out), log = file.path(work, "score-e4.log"))
hand <- c(rho_O = 2, m6_na = 0.0023333333333333335, m6_ell = 0.02, fdr_added = 0.085714285714285715,
    fdrcal_added = 0.0092879256965944269, E_cal = 0.069164195892782068, E_comp = -0.0030577533157512518,
    comp_power_first = 0.0030389363722697054, comp_null_second = -0.0060966896880209577,
    comp_null_first = -0.0049015102778543641, comp_power_second = 0.001843756962103112, D_all = 0.77777777777777779)
hand_jack <- c(rho_O = 0.057560026460742396, E_comp = 0.00027017648151611437)
hand_jack_setting <- c(rho_O = 0.12443840371712694, E_comp = 0.00059750294731672731)
close_enough <- function(x, y) isTRUE(all(abs(x - y) <= 1e-12 * pmax(abs(y), 1e-12)))
e4 <- if (status == 0L) {
    pooled <- utils::read.delim(file.path(e4_out, "pooled.tsv"), stringsAsFactors = FALSE)
    settings <- utils::read.delim(file.path(e4_out, "settings.tsv"), stringsAsFactors = FALSE)
    labels <- utils::read.delim(file.path(e4_out, "labels.tsv"), stringsAsFactors = FALSE)
    outcome <- utils::read.delim(file.path(e4_out, "outcome.tsv"), stringsAsFactors = FALSE)
    pv <- pooled[pooled$mask == "vd", ]
    est_pool <- pv$estimate[match(names(hand), pv$stat)]
    est_settings <- vapply(names(hand), \(s) max(abs(settings$estimate[settings$stat == s] - hand[[s]])), numeric(1L))
    jack_pool <- pv$jack_se[match(names(hand_jack), pv$stat)]
    sa <- settings[settings$cell == "A" & settings$mask == "vd5", ]
    jack_setting <- sa$jack_se[match(names(hand_jack_setting), sa$stat)]
    results <- c(estimates_pool = close_enough(est_pool, hand), estimates_settings = all(est_settings <=
        1e-12 * pmax(abs(hand), 1e-12)), jack_pool = close_enough(jack_pool, hand_jack),
        jack_setting = close_enough(jack_setting, hand_jack_setting),
        pooled_label = identical(labels$label[labels$scope == "pool" & labels$mask == "vd"], "bias"),
        setting_labels = all(labels$label[labels$scope == "setting"] == "bias"),
        attribution = identical(labels$attribution_B[labels$scope == "pool" & labels$mask == "vd"],
            "not specific to value dependence"),
        outcome = identical(outcome$outcome, "O1"))
    list(pass = all(results), detail = paste(names(results), ifelse(results, "ok", "FAIL"), sep = " ", collapse = "; "))
} else list(pass = FALSE, detail = sprintf("score-k197.R exit %d on the constructed table", status))
addCheck("E4", e4$pass, e4$detail)
if (status != 0L) cat(utils::tail(readLines(file.path(work, "score-e4.log"), warn = FALSE), 30L), sep = "\n")

# ---- E5 ------------------------------------------------------------------------------------------------------------
e5 <- lapply(c("e5a", "e5b"), \(name) run(score, c("deciding", l1_rows, file.path(work, name)),
    log = file.path(work, paste0("score-", name, ".log"))))
files_a <- sort(list.files(file.path(work, "e5a")))
files_b <- sort(list.files(file.path(work, "e5b")))
same <- identical(files_a, files_b) && length(files_a) > 0L && all(vapply(files_a, \(f)
    identical(sha(file.path(work, "e5a", f)), sha(file.path(work, "e5b", f))), logical(1L)))
addCheck("E5", all(unlist(e5) == 0L) && same, sprintf("two runs exit %s; %d files each; all identical: %s",
    paste(unlist(e5), collapse = ", "), length(files_a), same))
for (name in c("e5a", "e5b")) if (file.exists(file.path(work, paste0("score-", name, ".log"))) && any(unlist(e5) != 0L))
    cat(grep("Error|error", readLines(file.path(work, paste0("score-", name, ".log")), warn = FALSE), value = TRUE),
        sep = "\n")

# ---- G1's mechanics: the descriptive mode on the lifecycle rows must stop at G1 (they are not S2's draws) -----------
g1_out <- file.path(work, "g1m")
status_g1 <- run(score, c("descriptive", l1_rows, "docs/qa/workflow-execution/kit/kit-187/s2/cells.tsv", g1_out),
    log = file.path(work, "score-g1m.log"))
g1_log <- readLines(file.path(work, "score-g1m.log"), warn = FALSE)
g1_rows <- if (file.exists(file.path(g1_out, "g1.tsv"))) nrow(utils::read.delim(file.path(g1_out, "g1.tsv"))) else 0L
addCheck("G1m", status_g1 != 0L && any(grepl("gate G1 failed", g1_log, fixed = TRUE)) && g1_rows == 96L &&
    file.exists(file.path(g1_out, "g1-mismatch.tsv")) && !file.exists(file.path(g1_out, "descriptive.tsv")),
    sprintf("descriptive mode on r = 9001, 9002: exit %d, stopped at G1: %s, %d comparisons written, no estimand written: %s",
        status_g1, any(grepl("gate G1 failed", g1_log, fixed = TRUE)), g1_rows,
        !file.exists(file.path(g1_out, "descriptive.tsv"))))

# ---- E6 ------------------------------------------------------------------------------------------------------------
e6_cut <- file.path(work, "e6-cut")
e6_full <- file.path(work, "e6-full")
dir.create(file.path(e6_cut, "parts"), recursive = TRUE, showWarnings = FALSE)
all_tasks <- file.path(work, "e6-tasks.txt")
writeLines(sprintf("A %d", 9001:9004), all_tasks)
done <- function() sum(file.exists(file.path(e6_cut, "parts", sprintf("A-%04d.metrics.tsv", 9001:9004))))
child <- processx::process$new(rscript, c(sim, e6_cut, "9001", "9004", all_tasks),
    env = c(Sys.getenv(), K197_WORKERS = "1", K197_CELLS = "A"), stdout = file.path(work, "e6-cut.log"), stderr = "2>&1",
    cleanup_tree = TRUE)
deadline <- Sys.time() + 900
while (child$is_alive() && done() < 1L && Sys.time() < deadline) Sys.sleep(0.2)
killed <- child$is_alive()
child$kill_tree()
Sys.sleep(5)
complete_after_kill <- done()
missing <- setdiff(9001:9004, 9000L + which(file.exists(file.path(e6_cut, "parts",
    sprintf("A-%04d.metrics.tsv", 9001:9004)))))
missing_file <- file.path(work, "e6-missing.txt")
writeLines(sprintf("A %d", missing), missing_file)
status_cont <- run(sim, c(e6_cut, "9001", "9004", missing_file), "K197_CELLS=A", file.path(work, "e6-cont.log"))
started <- Sys.time()
status_full <- run(sim, c(e6_full, "9001", "9004"), "K197_CELLS=A", file.path(work, "e6-full.log"))
full_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
readGz <- function(dir) if (file.exists(file.path(dir, "A.tsv.gz"))) readLines(file.path(dir, "A.tsv.gz")) else NULL
readMetrics <- function(dir) if (file.exists(file.path(dir, "s2metrics.tsv"))) {
    m <- utils::read.delim(file.path(dir, "s2metrics.tsv"), colClasses = "character")
    m[setdiff(names(m), "seconds")]
} else NULL
rows_cut <- readGz(e6_cut)
rows_full <- readGz(e6_full)
replicate_counts <- if (is.null(rows_cut)) integer() else {
    cut_rows <- utils::read.delim(file.path(e6_cut, "A.tsv.gz"), colClasses = "character")
    table(unique(cut_rows[c("r", "mask")])$r)
}
e6_pass <- killed && complete_after_kill >= 1L && complete_after_kill < 4L && status_cont == 0L && status_full == 0L &&
    !is.null(rows_cut) && identical(rows_cut, rows_full) && identical(readMetrics(e6_cut), readMetrics(e6_full)) &&
    identical(as.integer(replicate_counts), rep(6L, 4L))
addCheck("E6", e6_pass, sprintf(paste("killed while running: %s; complete tasks after the kill: %d of 4; continued",
    "with %d task(s), exit %d; uninterrupted exit %d; assembled rows identical: %s; metrics identical (seconds",
    "excluded): %s; each replicate once with six masks: %s"), killed, complete_after_kill, length(missing), status_cont,
    status_full, identical(rows_cut, rows_full), identical(readMetrics(e6_cut), readMetrics(e6_full)),
    identical(as.integer(replicate_counts), rep(6L, 4L))))

# ---- public records ------------------------------------------------------------------------------------------------
checks <- do.call(rbind, checks)
utils::write.table(checks, file.path(public, "checks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
timing <- function(dir, label) if (file.exists(file.path(dir, "s2metrics.tsv"))) {
    m <- utils::read.delim(file.path(dir, "s2metrics.tsv"))
    data.frame(run = label, replicate_masks = nrow(m), mean_seconds = mean(m$seconds), max_seconds = max(m$seconds),
        stringsAsFactors = FALSE)
} else NULL
utils::write.table(rbind(timing(l1_rows, "L1 lifecycle (4 workers, 8 tasks)"),
    timing(e6_full, sprintf("E6 uninterrupted (4 workers, 4 tasks; wall %.0f s)", full_seconds))),
    file.path(public, "timing.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
row_files <- c(file.path(l1_rows, c(paste0(c("A", "B", "C", "D"), ".tsv.gz"), "s2metrics.tsv")),
    file.path(e6_full, "A.tsv.gz"))
utils::write.table(data.frame(file = sub(paste0(work, "/"), "", row_files, fixed = TRUE),
    sha256 = vapply(row_files, \(f) if (file.exists(f)) sha(f) else NA_character_, character(1L))),
    file.path(public, "hashes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("KIT-197 L1: checks\n")
print(checks[c("id", "pass")], row.names = FALSE)
cat("KIT-197 L1: end", format(Sys.time(), "%H:%M:%S"), "\n")
if (!all(checks$pass)) quit(status = 1L)
