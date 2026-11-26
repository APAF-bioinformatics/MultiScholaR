# KIT-281 parse-and-self-test send (script audit S09): the first contact of the three scripts with R, before L1. Evidence code under docs/qa: no package code
# is loaded by this driver (the one simulation task it starts loads the bound tree, as every simulation does).
#
# Usage, from the package root of the bound commit on a run-tests runner (full git history):
#   Rscript docs/qa/workflow-execution/kit/kit-281/smoke/smoke-k281.R <public out dir>
# It (1) parses sim-k281.R, score-k281.R and l1-k281.R; (2) runs the scorer self-test (selftest/run-selftest.sh) and records its SELFTEST line, the stop kind of
# the probe case and the completion of cal105; (3) runs ONE simulation task, smoke-flagged (K281_SMOKE=1, K211_LIFECYCLE=1, cell A, r = 9301, a seed range used by
# nothing else), into tempdir() and uploads only checks-sim.tsv and the start provenance; (4) runs the scorer on that directory and expects it to refuse (smoke flag).
# Nothing simulated leaves the runner except those two files. The smoke-flagged rows can never be scored as a result.
options(warnPartialMatchDollar = TRUE)
globalCallingHandlers(warning = function(w) {
    if (grepl("partial match of", tolower(conditionMessage(w)), fixed = TRUE)) stop(paste("partial-match guard:", conditionMessage(w)), call. = FALSE)
})
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
public <- args[[1L]]
dir.create(public, recursive = TRUE, showWarnings = FALSE)
here <- "docs/qa/workflow-execution/kit/kit-281"
work <- file.path(tempdir(), "k281-smoke")
dir.create(work, showWarnings = FALSE)
rscript <- file.path(R.home("bin"), "Rscript")
path_env <- paste0("PATH=", R.home("bin"), ":", Sys.getenv("PATH"))
results <- list()
addResult <- function(id, pass, detail) {
    results[[length(results) + 1L]] <<- data.frame(id = id, pass = isTRUE(pass), detail = detail, stringsAsFactors = FALSE)
    cat(sprintf("%s\t%s\t%s\n", id, isTRUE(pass), detail))
}
readText <- function(path) utils::read.delim(path, colClasses = "character", na.strings = character(0), check.names = FALSE, stringsAsFactors = FALSE)
cat("KIT-281 smoke: start", format(Sys.time(), "%H:%M:%S"), "\n")

# ---- 1. three parses ------------------------------------------------------------------------------------------------------------------------------------
for (script in c("sim-k281.R", "score-k281.R", "l1-k281.R")) {
    outcome <- tryCatch({
        invisible(parse(file = file.path(here, script)))
        "parses"
    }, error = function(e) paste("does not parse:", substr(conditionMessage(e), 1L, 200L)))
    addResult(id = paste("parse", script), pass = identical(outcome, "parses"), detail = outcome)
}

# ---- 2. the scorer self-test (with the S01 and S07 cases) -----------------------------------------------------------------------------------------------------
st_dir <- normalizePath(file.path(here, "selftest"))
st_out <- file.path(work, "selftest-out")
status_st <- system2("bash", c(shQuote(file.path(st_dir, "run-selftest.sh")), shQuote(st_dir), shQuote(st_out)), env = path_env,
    stdout = file.path(work, "selftest.log"), stderr = file.path(work, "selftest.log"))
st_log <- if (file.exists(file.path(work, "selftest.log"))) readLines(file.path(work, "selftest.log"), warn = FALSE) else character()
st_line <- grep("^SELFTEST: ", st_log, value = TRUE)
writeLines(c(st_line, grep("DIFFERS", st_log, value = TRUE)), file.path(public, "selftest-summary.txt"))
st_zero <- length(st_line) == 1L && grepl(", 0 differ$", st_line)
addResult(id = "selftest", pass = status_st == 0L && st_zero, detail = if (length(st_line)) st_line[[length(st_line)]] else "no SELFTEST line")
probe_row <- if (file.exists(file.path(st_out, "probe", "outcome.tsv"))) readText(file.path(st_out, "probe", "outcome.tsv"))
addResult(id = "probe stopped as instrument", pass = !is.null(probe_row) && identical(probe_row[["outcome"]], "Stop") && identical(probe_row[["kind"]], "instrument"),
    detail = if (is.null(probe_row)) "no outcome.tsv" else paste(probe_row[["outcome"]], probe_row[["kind"]]))
cal_row <- if (file.exists(file.path(st_out, "cal105", "outcome.tsv"))) readText(file.path(st_out, "cal105", "outcome.tsv"))
addResult(id = "cal105 completed", pass = !is.null(cal_row) && identical(cal_row[["outcome"]], "Calibrated"), detail = if (is.null(cal_row)) "no outcome.tsv" else cal_row[["outcome"]])

# ---- 3. one smoke-flagged simulation task (cell A, r = 9301, three masks) -------------------------------------------------------------------------------------
sim_root <- file.path(work, "rows")
sim_job <- file.path(sim_root, "r9301")
status_sim <- system2(rscript, c(shQuote(file.path(here, "sim-k281.R")), shQuote(sim_job), "9301", "9301"),
    env = c("K281_SMOKE=1", "K211_LIFECYCLE=1", "K211_CELLS=A"), stdout = file.path(work, "sim.log"), stderr = file.path(work, "sim.log"))
addResult(id = "simulation task exit", pass = status_sim == 0L, detail = sprintf("exit %d", status_sim))
checks_file <- file.path(sim_job, "checks-sim.tsv")
if (file.exists(checks_file)) {
    sim_checks <- readText(checks_file)
    file.copy(checks_file, file.path(public, "checks-sim.tsv"), overwrite = TRUE)
    for (id in c("E2", "E3", "E7", "E10")) {
        row <- sim_checks[sim_checks[["check"]] == id, , drop = FALSE]
        addResult(id = paste("sim", id), pass = nrow(row) == 1L && identical(row[["units"]], "3") && identical(row[["passed"]], "3"),
            detail = if (nrow(row) == 1L) paste(row[["passed"]], "of", row[["units"]]) else "absent")
    }
} else addResult(id = "sim checks", pass = FALSE, detail = "checks-sim.tsv absent")
start_file <- file.path(sim_job, "provenance-sim-start.json")
if (file.exists(start_file)) file.copy(start_file, file.path(public, "provenance-sim-start.json"), overwrite = TRUE)
final_file <- file.path(sim_job, "provenance-sim.json")
final_text <- if (file.exists(final_file)) paste(readLines(final_file, warn = FALSE), collapse = "\n") else ""
addResult(id = "provenance smoke true and complete", pass = grepl("\"smoke\"[[:space:]]*:[[:space:]]*true", final_text) && grepl("\"status\"[[:space:]]*:[[:space:]]*\"complete\"", final_text),
    detail = "provenance-sim.json (the marker assertion ran at load: exit 0)")

# ---- 4. the scorer refuses the smoke directory ------------------------------------------------------------------------------------------------------------------
status_score <- system2(rscript, c(shQuote(file.path(here, "score-k281.R")), shQuote(sim_root), shQuote(file.path(work, "score-smoke"))),
    env = c("K281_RANGE=9301,9301", "K211_B_BOOT=1000"), stdout = file.path(work, "score.log"), stderr = file.path(work, "score.log"))
score_log <- if (file.exists(file.path(work, "score.log"))) readLines(file.path(work, "score.log"), warn = FALSE) else character()
addResult(id = "scorer refuses the smoke directory", pass = status_score == 1L && any(grepl("smoke flag", score_log, fixed = TRUE)),
    detail = sprintf("exit %d", status_score))

# ---- public record -------------------------------------------------------------------------------------------------------------------------------------------------------
results <- do.call(rbind, results)
utils::write.table(results, file.path(public, "smoke-checks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("KIT-281 smoke: end", format(Sys.time(), "%H:%M:%S"), "\n")
if (!all(results[["pass"]])) quit(status = 1L)
