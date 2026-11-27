# KIT-275 r2.1 pilot self-test: runs on a runner from the package root (run-tests `scripts`). Steps: the bound scripts' SHA-256 (first 16 hex digits)
# against the declared values; parse the pilot's scripts; build the constructed cases and predictions (python3 -I pilot-selftest.py make);
# run score-pilot.R on each case; compare (python3 -I pilot-selftest.py check). Output directory: K275_PILOT_SELFTEST_OUT (default k275-pilot-selftest).
# Exit status 1 when a hash differs, a script does not parse or a prediction differs.
options(warnPartialMatchDollar = TRUE)
globalCallingHandlers(warning = function(w) {
    if (grepl("partial match of", tolower(conditionMessage(w)), fixed = TRUE)) stop(paste("partial-match guard:", conditionMessage(w)), call. = FALSE)
})
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])))
diag_dir <- dirname(here)
kit <- dirname(diag_dir)
work <- Sys.getenv("K275_PILOT_SELFTEST_OUT", "k275-pilot-selftest")
dir.create(work, recursive = TRUE, showWarnings = FALSE)
work <- normalizePath(work)
bound <- c("sim-k275.R" = "2cc87ea569eef7fd", "score-k275.R" = "a43a2c10c9a41ec2", "fixture-k275.R" = "407949a9a6c0f90f")
hashes_ok <- vapply(names(bound), function(f) {
    sha <- strsplit(system2("sha256sum", shQuote(file.path(kit, f)), stdout = TRUE), " ")[[1L]][[1L]]
    ok <- startsWith(sha, bound[[f]])
    cat(sprintf("bound: %s %s (%s)\n", f, substr(sha, 1L, 16L), if (ok) "as declared" else "DIFFERS"))
    ok
}, logical(1L))
if (!all(hashes_ok)) quit(save = "no", status = 1L)
scripts <- c(file.path(diag_dir, c("pilot-k275.R", "score-pilot.R")), file.path(here, "run-pilot-selftest.R"))
parsed <- vapply(scripts, function(f) {
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
made <- system2("python3", c("-I", shQuote(file.path(here, "pilot-selftest.py")), "make", shQuote(work)), stdout = TRUE, stderr = TRUE)
cat(made, sep = "\n")
index <- utils::read.delim(file.path(work, "index.tsv"), stringsAsFactors = FALSE, colClasses = "character", na.strings = character())
runCase <- function(i) {
    name <- index[["case"]][[i]]
    out <- file.path(work, "out", name)
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    code <- system2("Rscript", c(shQuote(file.path(diag_dir, "score-pilot.R")), shQuote(file.path(work, "cases", name, "rows")), shQuote(out)),
        stdout = file.path(out, "stdout.txt"), stderr = file.path(out, "stderr.txt"))
    writeLines(as.character(code), file.path(out, "exit.txt"))
    code
}
codes <- lapply(seq_len(nrow(index)), runCase)
cat(sprintf("selftest: %d cases run, scorer exit codes %s\n", nrow(index), paste(unlist(codes), collapse = " ")))
checked <- system2("python3", c("-I", shQuote(file.path(here, "pilot-selftest.py")), "check", shQuote(work)), stdout = TRUE, stderr = TRUE)
cat(checked, sep = "\n")
status <- attr(checked, "status")
quit(save = "no", status = if (is.null(status)) 0L else status)
