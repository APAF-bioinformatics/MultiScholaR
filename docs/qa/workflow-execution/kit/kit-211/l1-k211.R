# KIT-211 lifecycle L1 (ticket 6.7, 7.3, VAL-11): the checks E1 to E10 and G2 on the scripts as declared, blind. Evidence code
# under docs/qa: no package code.
#
# Usage, from the package root of the bound commit on a run-tests runner (full git history):
#   Rscript docs/qa/workflow-execution/kit/kit-211/l1-k211.R <public out dir>
# Blind: the simulated rows and every score computed from them stay in a private work directory under tempdir(), which is
# not uploaded; <public out dir> receives only checks.tsv (one row per check: id, units compared, units passed, pass,
# mechanics details; the maximum absolute deviation of an identity), timing.tsv (seconds per replicate-mask over all arms
# together, for G2, and E6's wall seconds) and hashes.tsv (SHA-256 of the scripts and row files, the bound commit's 40
# digits). Nothing is printed from the rows or the scores; a failed child prints its check id and its condition's message.
#   E2, E3, E7, E9, E10: sim-k211.R with K211_LIFECYCLE=1 on cells A to D, Av and Ag, r = 9101 and r = 9102 as two ranges
#       (rows/r9101 and rows/r9102, each holding every cell), all masks, every arm and reference.
#   E8: K211_DETERMINISM=1 on cell A r = 9101, 9102.
#   E1: S alone (K211_ARMS=S) on cells A to D, r = 101 to 105, against kit-197/sb/rows.
#   E4: score-k211.R on the constructed tables of e4/ (e4hand-k211.py) against the hand values and the designed outcomes.
#   E5: score-k211.R twice on the two ranges' rows (the join of 6.5): every output file identical.
#   E6: cell A, r = 9101 to 9104 with one worker, killed (SIGKILL to the process tree) once one task is complete;
#       continued from the list of missing tasks; the assembled rows equal an uninterrupted run (four workers).
#   G2: the planned minutes of a counted job from this run's seconds per replicate-mask.
# E0 (the verbatim blocks) is local: e0-k211.py, l1/e0.txt.

# ==== verbatim: l1-k197.R:34-37 at ef63c918 ====
run <- function(script, arguments, env = character(), log) {
    status <- system2(rscript, c(shQuote(script), shQuote(arguments)), env = env, stdout = log, stderr = log)
    status
}
# ==== end verbatim ====

# ==== verbatim: l1-k197.R:38-38 at ef63c918 ====
sha <- function(path) digest::digest(file = path, algo = "sha256")
# ==== end verbatim ====

# ---- KIT-211 L1 (not copied) ---------------------------------------------------------------------------------------------
BOUND <- "894b9cc821a70b03baab1e1fa47b3d47440c7ed5"
RANGES <- c(9101L, 9102L)
CELLS_L1 <- "A,B,C,D,Av,Ag"
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
public <- args[[1L]]
dir.create(public, recursive = TRUE, showWarnings = FALSE)
here <- "docs/qa/workflow-execution/kit/kit-211"
old_rows <- "docs/qa/workflow-execution/kit/kit-197/sb/rows"
sim <- file.path(here, "sim-k211.R")
score <- file.path(here, "score-k211.R")
work <- file.path(tempdir(), "k211-l1")
dir.create(work, showWarnings = FALSE)
rscript <- file.path(R.home("bin"), "Rscript")
checks <- list()
addCheck <- function(id, units, passed, detail) checks[[length(checks) + 1L]] <<- data.frame(id = id, units = units,
    passed = passed, pass = isTRUE(units == passed && units > 0), detail = detail, stringsAsFactors = FALSE)
num <- function(x) ifelse(is.na(x), "NA", sprintf("%.17g", x))
# Only the condition's class and message of a failed child, never a log tail (the messages carry no estimate).
failureNote <- function(id, log, errors) {
    note <- if (file.exists(errors)) readLines(errors, warn = FALSE)[1L] else grep("^Error", readLines(log, warn = FALSE),
        value = TRUE)[1L]
    cat(sprintf("%s failed: %s\n", id, substr(if (is.na(note)) "no message" else note, 1L, 200L)))
}
cat("KIT-211 L1: start", format(Sys.time(), "%H:%M:%S"), "\n")
l1_started <- Sys.time()

# ---- the lifecycle simulation (E2, E3, E7, E9, E10; two ranges per cell, joined by the scorer: E5) -----------------------
rows_root <- file.path(work, "rows")
range_seconds <- numeric()
status_sim <- vapply(RANGES, \(r) {
    started <- Sys.time()
    status <- run(sim, c(file.path(rows_root, paste0("r", r)), r, r), c("K211_LIFECYCLE=1", paste0("K211_CELLS=", CELLS_L1)),
        file.path(work, sprintf("sim-%d.log", r)))
    range_seconds[[as.character(r)]] <<- as.numeric(difftime(Sys.time(), started, units = "secs"))
    if (status != 0L) failureNote(sprintf("simulation r = %d", r), file.path(work, sprintf("sim-%d.log", r)),
        file.path(rows_root, paste0("r", r), "errors.txt"))
    status
}, numeric(1L))
sim_checks <- do.call(rbind, lapply(RANGES, \(r) {
    f <- file.path(rows_root, paste0("r", r), "checks-sim.tsv")
    if (file.exists(f)) utils::read.delim(f, stringsAsFactors = FALSE)
}))
warnings_counted <- sum(vapply(RANGES, \(r) {
    f <- file.path(rows_root, paste0("r", r), "e9-warnings.tsv")
    if (file.exists(f)) utils::read.delim(f)[[1L]] else NA_real_
}, numeric(1L)))
for (id in c("E2", "E3", "E7", "E9", "E10")) {
    sub <- if (is.null(sim_checks)) NULL else sim_checks[sim_checks$check == id, , drop = FALSE]
    ok <- all(status_sim == 0) && !is.null(sub) && nrow(sub) == length(RANGES)
    addCheck(id, if (ok) sum(sub$units) else 0L, if (ok) sum(sub$passed) else 0L,
        if (!ok) sprintf("no result (simulation exits %s)", paste(status_sim, collapse = ", ")) else if (id == "E10")
            sprintf("replicate-masks over both ranges; maximum absolute deviation of logFC %s", num(max(sub$max_abs_deviation))) else if (id == "E9")
            sprintf("replicate-masks over both ranges; classed warnings muffled and counted: %s", warnings_counted) else
            "replicate-masks over both ranges")
}

# ---- E8: the determinism mechanics on cell A, r = 9101, 9102 -------------------------------------------------------------
e8_dir <- file.path(work, "e8")
status_e8 <- run(sim, c(e8_dir, "9101", "9102"), c("K211_CELLS=A", "K211_DETERMINISM=1"), file.path(work, "sim-e8.log"))
det <- if (file.exists(file.path(e8_dir, "determinism.tsv"))) utils::read.delim(file.path(e8_dir, "determinism.tsv"))
if (status_e8 != 0L) failureNote("E8 simulation", file.path(work, "sim-e8.log"), file.path(e8_dir, "errors.txt"))
addCheck("E8", if (is.null(det)) 0L else nrow(det), if (is.null(det)) 0L else sum(det$rows_identical & det$metrics_identical),
    "K211_DETERMINISM=1 on cell A r = 9101, 9102 (rows and metrics, seconds excluded)")

# ---- E1: S on KIT-197's r = 101 to 105 against its committed rows -------------------------------------------------------
e1_dir <- file.path(work, "e1")
status_e1 <- run(sim, c(e1_dir, "101", "105"), c("K211_CELLS=A,B,C,D", "K211_ARMS=S"), file.path(work, "sim-e1.log"))
if (status_e1 != 0L) failureNote("E1 simulation", file.path(work, "sim-e1.log"), file.path(e1_dir, "errors.txt"))
e1 <- if (status_e1 == 0L) {
    classes <- c(cell = "character", mask = "character", subset = "character", class = "character", stratum = "character")
    all_bad <- NULL
    da_equal <- TRUE
    da_columns <- c("n_tested", "n_p05", "n_call", "sum_err", "sum_s_err", "sum_sgn_err")
    for (cellname in c("A", "B", "C", "D")) {
        mine <- utils::read.delim(file.path(e1_dir, paste0(cellname, ".tsv.gz")), colClasses = c(classes, arm = "character"))
        old <- utils::read.delim(file.path(old_rows, paste0(cellname, ".tsv.gz")), colClasses = classes)
        old <- old[old$r %in% 101:105 & old$subset != "ref", , drop = FALSE]
        mine <- mine[mine$arm == "S" & mine$subset %in% c("added", "kept_both", "whole"), , drop = FALSE]
        key <- function(d) paste(paste(d$cell, d$mask, d$r, sep = "|"), d$subset, d$class, d$stratum, d$ell, sep = "\t")
        km <- rowsum(mine[c("n", da_columns)], key(mine))
        ko <- rowsum(old[c("n", da_columns)], key(old))
        keys <- union(rownames(km), rownames(ko))
        bad <- is.na(km[keys, "n"]) | is.na(ko[keys, "n"]) | km[keys, "n"] != ko[keys, "n"]
        metrics_new <- utils::read.delim(file.path(e1_dir, "metrics.tsv"), colClasses = c(cell = "character",
            mask = "character", arm = "character"))
        metrics_new <- metrics_new[metrics_new$cell == cellname & metrics_new$arm == "S", , drop = FALSE]
        metrics_old <- utils::read.delim(file.path(old_rows, "s2metrics.tsv"), colClasses = c(cell = "character",
            mask = "character"))
        metrics_old <- metrics_old[metrics_old$cell == cellname & metrics_old$r %in% 101:105, , drop = FALSE]
        m <- merge(metrics_new, metrics_old, by = c("cell", "mask", "r"), suffixes = c(".new", ".old"), all = TRUE)
        unit_metric_bad <- !(m$kept == m$kept_h & m$added.new == m$added.old & m$threshold.new == m$threshold.old) |
            is.na(m$kept) | is.na(m$kept_h)
        names(unit_metric_bad) <- paste(m$cell, m$mask, m$r, sep = "|")
        unit_bad <- c(tapply(bad, sub("\t.*$", "", keys), any), unit_metric_bad)
        all_bad <- c(all_bad, tapply(unit_bad, names(unit_bad), any))
        da_equal <- da_equal && isTRUE(all.equal(as.matrix(km[intersect(rownames(ko), rownames(km)), da_columns]),
            as.matrix(ko[intersect(rownames(ko), rownames(km)), da_columns]), tolerance = 1e-9, check.attributes = FALSE))
    }
    list(units = length(all_bad), passed = sum(!all_bad), da_equal = da_equal)
} else list(units = 0L, passed = 0L, da_equal = NA)
addCheck("E1", e1$units, e1$passed, sprintf(paste("S's filter-dependent columns (membership counts per subset, class, stratum and ell; kept and added",
    "counts; threshold) against KIT-197's committed rows, cells A to D, r = 101 to 105, six masks (120 replicate-masks expected);",
    "G0 lists DA-file changes, so the DA-dependent columns are written as one flag: all equal %s"), e1$da_equal))

# ---- E4: the scorer on the constructed tables (e4hand-k211.py wrote them and the hand values) --------------------------
e4_root <- file.path(here, "e4")
expect <- utils::read.delim(file.path(e4_root, "expect.tsv"), colClasses = c(scenario = "character", outcome = "character",
    class_A = "character", class_B = "character", class_Bp = "character", withheld = "character", pc_pass = "character",
    f4 = "character"))
hand <- utils::read.delim(file.path(e4_root, "hand.tsv"), colClasses = c(scope = "character", cell = "character",
    mask = "character", stat = "character", field = "character"))
e4_results <- logical()
e4Note <- function(name, ok) e4_results[[name]] <<- isTRUE(ok)
close_enough <- function(x, y) isTRUE(all(is.finite(x)) && all(abs(x - y) <= 1e-12 * pmax(abs(y), 1e-12)))
scored <- list()
for (i in seq_len(nrow(expect))) {
    name <- expect$scenario[[i]]
    boot_i <- as.integer(expect$B[[i]])
    out_dir <- file.path(work, "e4", name)
    status <- run(score, c(file.path(e4_root, name), out_dir), if (boot_i != 10000L) paste0("K211_B_BOOT=", boot_i) else
        character(), file.path(work, paste0("score-e4-", name, ".log")))
    scored[[name]] <- list(status = status, dir = out_dir)
    outcome_file <- file.path(out_dir, "outcome.tsv")
    got <- if (file.exists(outcome_file)) utils::read.delim(outcome_file, colClasses = c(outcome = "character"), stringsAsFactors = FALSE)
    e4Note(paste0(name, ":exit"), identical(as.integer(status != 0L), as.integer(expect$exit[[i]])))
    e4Note(paste0(name, ":outcome"), !is.null(got) && identical(got$outcome[[1L]], expect$outcome[[i]]))
    if (expect$exit[[i]] == 0L && !is.null(got)) {
        e4Note(paste0(name, ":classes"), identical(c(got$class_A, got$class_B, got$class_Bp), c(expect$class_A[[i]],
            expect$class_B[[i]], expect$class_Bp[[i]])))
        e4Note(paste0(name, ":withheld"), identical(as.character(got$withheld), expect$withheld[[i]]) ||
            (is.na(got$withheld) && !nzchar(expect$withheld[[i]])))
        e4Note(paste0(name, ":pc_f4"), identical(as.character(got$pc_pass), expect$pc_pass[[i]]) &&
            identical(as.character(got$f4), expect$f4[[i]]))
    }
    if (status != 0L && expect$exit[[i]] == 0L) failureNote(paste("E4 scenario", name), file.path(work, paste0("score-e4-", name, ".log")), "none")
}
base_files <- file.path(scored$base$dir, c("pooled.tsv", "settings.tsv"))
if (all(file.exists(base_files))) {
    pooled <- utils::read.delim(base_files[[1L]], colClasses = c(scope = "character", cell = "character",
        mask = "character", stat = "character"))
    settings <- utils::read.delim(base_files[[2L]], colClasses = c(scope = "character", cell = "character",
        mask = "character", stat = "character"))
    mine <- vapply(seq_len(nrow(hand)), \(i) {
        t <- if (hand$scope[[i]] == "pool") pooled[pooled$scope == "pool" & pooled$cell == hand$cell[[i]] & pooled$mask == hand$mask[[i]], ] else
            settings[settings$scope == hand$scope[[i]] & settings$cell == hand$cell[[i]] & settings$mask == hand$mask[[i]], ]
        v <- t[[if (hand$field[[i]] == "estimate") "estimate" else "jack_se"]][t$stat == hand$stat[[i]]]
        if (length(v) == 1L) as.numeric(v) else NA_real_
    }, numeric(1L))
    e4Note("hand_estimates", close_enough(mine[hand$field == "estimate"], hand$value[hand$field == "estimate"]))
    e4Note("hand_jackknife", close_enough(mine[hand$field == "jack_se"], hand$value[hand$field == "jack_se"]))
} else {
    e4Note("hand_estimates", FALSE)
    e4Note("hand_jackknife", FALSE)
}
join_dir <- scored$join$dir
if (file.exists(file.path(join_dir, "replicates.tsv"))) {
    reps <- utils::read.delim(file.path(join_dir, "replicates.tsv"), colClasses = c(cell = "character"))
    rg <- utils::read.delim(file.path(join_dir, "ranges.tsv"), colClasses = c(cell = "character"))
    st <- utils::read.delim(file.path(join_dir, "settings.tsv"), colClasses = c(scope = "character", cell = "character",
        mask = "character", stat = "character"))
    row <- st[st$scope == "setting" & st$cell == "A" & st$mask == "vd5" & st$stat == "S.rho_l", ]
    pl <- utils::read.delim(file.path(join_dir, "pooled.tsv"), colClasses = c(scope = "character", cell = "character",
        mask = "character", stat = "character"))
    e4Note("join:R", identical(setNames(reps$R, reps$cell), c(A = 3L, B = 3L, C = 5L, D = 3L)) && sum(rg$incomplete) == 1L &&
        pl$R[pl$scope == "pool" & pl$cell == "ten" & pl$mask == "vd" & pl$stat == "S.rho_l"] == 14L)
    e4Note("join:empty_resamples", nrow(row) == 1L && row$boot_nonfinite > 0 && is.finite(row$estimate) && is.finite(row$lo95) &&
        is.finite(row$hi95))
} else e4Note("join:R", FALSE)
addCheck("E4", length(e4_results), sum(e4_results), paste(names(e4_results), ifelse(e4_results, "ok", "FAIL"), sep = " ",
    collapse = "; "))

# ---- E5: the deciding scoring twice on the lifecycle rows: every output file identical ------------------------------------
e5 <- lapply(c("e5a", "e5b"), \(name) run(score, c(rows_root, file.path(work, name)), log = file.path(work, paste0("score-", name, ".log"))))
files_a <- sort(list.files(file.path(work, "e5a")))
files_b <- sort(list.files(file.path(work, "e5b")))
same <- identical(files_a, files_b) && length(files_a) > 0L && all(vapply(files_a, \(f) identical(sha(file.path(work, "e5a", f)),
    sha(file.path(work, "e5b", f))), logical(1L)))
if (any(unlist(e5) != 0L)) failureNote("E5 scoring", file.path(work, "score-e5a.log"), "none")
addCheck("E5", 2L, if (all(unlist(e5) == 0L) && same) 2L else 0L, sprintf("two runs exit %s; %d files each; all identical: %s",
    paste(unlist(e5), collapse = ", "), length(files_a), same))

# ---- E6: cell A, r = 9101 to 9104, killed once one task is complete, continued; equal to an uninterrupted run ----------
e6_cut <- file.path(work, "e6-cut")
e6_full <- file.path(work, "e6-full")
dir.create(file.path(e6_cut, "parts"), recursive = TRUE, showWarnings = FALSE)
all_tasks <- file.path(work, "e6-tasks.txt")
writeLines(sprintf("A %d", 9101:9104), all_tasks)
done <- function() sum(file.exists(file.path(e6_cut, "parts", sprintf("A-%04d.metrics.tsv", 9101:9104))))
child <- processx::process$new(rscript, c(sim, e6_cut, "9101", "9104", all_tasks),
    env = c(Sys.getenv(), K211_WORKERS = "1", K211_CELLS = "A"), stdout = file.path(work, "e6-cut.log"), stderr = "2>&1",
    cleanup_tree = TRUE)
deadline <- Sys.time() + 1800
while (child$is_alive() && done() < 1L && Sys.time() < deadline) Sys.sleep(0.2)
killed <- child$is_alive()
invisible(child$kill_tree())
Sys.sleep(5)
complete_after_kill <- done()
missing <- setdiff(9101:9104, 9100L + which(file.exists(file.path(e6_cut, "parts", sprintf("A-%04d.metrics.tsv", 9101:9104)))))
missing_file <- file.path(work, "e6-missing.txt")
writeLines(sprintf("A %d", missing), missing_file)
status_cont <- run(sim, c(e6_cut, "9101", "9104", missing_file), "K211_CELLS=A", file.path(work, "e6-cont.log"))
started <- Sys.time()
status_full <- run(sim, c(e6_full, "9101", "9104"), "K211_CELLS=A", file.path(work, "e6-full.log"))
full_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
readGz <- function(dir) if (file.exists(file.path(dir, "A.tsv.gz"))) readLines(file.path(dir, "A.tsv.gz")) else NULL
readMetrics <- function(dir) if (file.exists(file.path(dir, "metrics.tsv"))) {
    m <- utils::read.delim(file.path(dir, "metrics.tsv"), colClasses = "character")
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
addCheck("E6", 1L, if (e6_pass) 1L else 0L, sprintf(paste("killed while running: %s; complete tasks after the kill: %d of 4; continued",
    "with %d task(s), exit %d; uninterrupted exit %d; assembled rows identical: %s; metrics identical (seconds",
    "excluded): %s; each replicate once with six masks: %s"), killed, complete_after_kill, length(missing), status_cont,
    status_full, identical(rows_cut, rows_full), identical(readMetrics(e6_cut), readMetrics(e6_full)),
    identical(as.integer(replicate_counts), rep(6L, 4L))))

# ---- G2: the planned minutes of a counted job from this run's seconds per replicate-mask --------------------------------
seconds_of <- function(cells) {
    m <- do.call(rbind, lapply(RANGES, \(r) {
        f <- file.path(rows_root, paste0("r", r), "metrics.tsv")
        if (file.exists(f)) utils::read.delim(f, colClasses = c(cell = "character", mask = "character", arm = "character"))
    }))
    if (is.null(m)) return(NA_real_)
    as.numeric(m$seconds[m$arm == "S" & m$cell %in% cells])
}
s_deciding <- seconds_of(c("A", "B", "C", "D", "Av"))
s_ag <- seconds_of("Ag")
plan_deciding <- 600 * mean(s_deciding) / 4 / 60
plan_ag <- 600 * mean(s_ag) / 4 / 60
addCheck("G2", 2L, sum(c(plan_deciding, plan_ag) <= 150, na.rm = TRUE), sprintf(paste("planned script minutes of a counted job at this run's seconds per",
    "replicate-mask (600 replicate-masks, 4 workers; the limit is 150): deciding cells %.1f, cell Ag %.1f; timeouts 120 and 139 minutes"),
    plan_deciding, plan_ag))

# ---- public records ------------------------------------------------------------------------------------------------------
checks <- do.call(rbind, checks)
utils::write.table(checks, file.path(public, "checks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
timing <- data.frame(run = c("L1 deciding cells A to D and Av (4 workers; per replicate-mask, all arms together)",
    "L1 cell Ag (4 workers; per replicate-mask, all arms together)", "E6 uninterrupted (4 workers, 4 tasks, wall seconds)"),
    replicate_masks = c(length(s_deciding), length(s_ag), NA), mean_seconds = c(mean(s_deciding), mean(s_ag), full_seconds),
    max_seconds = c(max(s_deciding), max(s_ag), NA), stringsAsFactors = FALSE)
utils::write.table(timing, file.path(public, "timing.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
hash_files <- c(file.path(here, c("sim-k211.R", "score-k211.R", "l1-k211.R")), list.files(rows_root, "\\.(tsv|gz)$", recursive = TRUE,
    full.names = TRUE), file.path(e6_full, "A.tsv.gz"))
hashes <- data.frame(file = c("bound commit", "package commit", sub(paste0(work, "/"), "", hash_files, fixed = TRUE)),
    sha256 = c(BOUND, system2("git", c("rev-parse", "HEAD"), stdout = TRUE), vapply(hash_files, \(f) if (file.exists(f)) sha(f) else
        NA_character_, character(1L))), stringsAsFactors = FALSE)
utils::write.table(hashes, file.path(public, "hashes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("KIT-211 L1: checks\n")
print(checks[c("id", "units", "passed", "pass")], row.names = FALSE)
cat("KIT-211 L1: end", format(Sys.time(), "%H:%M:%S"), "\n")
if (!all(checks$pass)) quit(status = 1L)
