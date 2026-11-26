# KIT-275 E4: runs the scorer on the constructed cases and compares it with the predictions derived in Python. Run from the package root
# by the run-tests `scripts` mechanism (Rscript). Steps: parse every script; build the cases and predictions (python3 -I make-selftest.py);
# run score-k275.R on each case, four at a time; compare (python3 -I check-selftest.py). Output directory: K275_SELFTEST_OUT (default
# k275-selftest). Exit status 1 when a script does not parse or a prediction differs.
options(warnPartialMatchDollar = TRUE)
globalCallingHandlers(warning = function(w) {
    if (grepl("partial match of", tolower(conditionMessage(w)), fixed = TRUE)) stop(paste("partial-match guard:", conditionMessage(w)), call. = FALSE)
})
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])))
kit <- dirname(here)
work <- Sys.getenv("K275_SELFTEST_OUT", "k275-selftest")
dir.create(work, recursive = TRUE, showWarnings = FALSE)
work <- normalizePath(work)
scripts <- file.path(kit, c("sim-k275.R", "score-k275.R", "fixture-k275.R"))
parsed <- vapply(c(scripts, file.path(here, "run-selftest.R")), function(f) {
    ok <- tryCatch({
        parse(file = f)
        TRUE
    }, error = function(e) {
        cat(sprintf("parse: %s FAILED: %s\n", basename(f), conditionMessage(e)))
        FALSE
    })
    if (ok) cat(sprintf("parse: %s parses\n", basename(f)))
    ok
}, logical(1L))
if (!all(parsed)) quit(save = "no", status = 1L)
cat("parse: all scripts parse\n")
# E0 (local and git only): the generator copy against 0866982f lists exactly the marked changes
e0 <- system2("python3", c("-I", shQuote(file.path(kit, "e0-diff.py"))), stdout = TRUE, stderr = TRUE)
cat(e0, sep = "\n")
if (!is.null(attr(e0, "status")) && attr(e0, "status") != 0L) quit(save = "no", status = 1L)
made <- system2("python3", c("-I", shQuote(file.path(here, "make-selftest.py")), shQuote(work)), stdout = TRUE, stderr = TRUE)
cat(made, sep = "\n")
index <- utils::read.delim(file.path(work, "index.tsv"), stringsAsFactors = FALSE, colClasses = "character", na.strings = character())
runCase <- function(i) {
    name <- index[["case"]][[i]]
    mode <- index[["mode"]][[i]]
    folder <- file.path(work, "cases", name)
    out <- file.path(work, "out", name)
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    scorer <- file.path(kit, "score-k275.R")
    env <- strsplit(index[["env"]][[i]], ";", fixed = TRUE)[[1L]]
    command <- switch(mode,
        score = c(shQuote(scorer), shQuote(file.path(folder, "rows")), shQuote(out)),
        `score-nog1` = c(shQuote(scorer), shQuote(file.path(folder, "rows")), shQuote(out)),
        blind = c(shQuote(scorer), "--blind", shQuote(file.path(folder, "rows")), shQuote(out), name),
        labels = c(shQuote(scorer), "--labels", shQuote(file.path(folder, "table.tsv")), shQuote(out)),
        controls = c(shQuote(scorer), shQuote(file.path(folder, "rows")), shQuote(out)),
        g1 = c(shQuote(scorer), "--g1", shQuote(file.path(folder, "g1.tsv")), shQuote(out)))
    if (identical(mode, "score")) env <- c(env, sprintf("K275_G1=%s", file.path(folder, "g1.tsv")))
    code <- system2("Rscript", command, env = env, stdout = file.path(out, "stdout.txt"), stderr = file.path(out, "stderr.txt"))
    writeLines(as.character(code), file.path(out, "exit.txt"))
    code
}
codes <- parallel::mclapply(seq_len(nrow(index)), runCase, mc.cores = 4L, mc.preschedule = FALSE)
cat(sprintf("selftest: %d cases run, scorer exit codes %s\n", nrow(index), paste(unlist(codes), collapse = " ")))
checked <- system2("python3", c("-I", shQuote(file.path(here, "check-selftest.py")), shQuote(work)), stdout = TRUE, stderr = TRUE)
cat(checked, sep = "\n")
status <- attr(checked, "status")
quit(save = "no", status = if (is.null(status)) 0L else status)
