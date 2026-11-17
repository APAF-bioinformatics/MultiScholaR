# KIT-124 on hosted runners: a simulator's out dir, checked before it is uploaded.
#
# Usage: Rscript --vanilla check-output.R <sim: s1..s4> <out dir> <cell> <runs>
#
# Two checks; either failing gives exit status 1, so the job fails and uploads nothing.
#   Complete: s<i>-runs.tsv holds the cell only, with runs 1..<runs>, every run with the same number of rows, and not
#   every row crashed (a table of nothing but crashes is a broken environment, not a result).
#   Clean: no file of the out dir (text read by line, .rds by its character values) holds compiler or build output,
#   a source location of the internal packages, or anything shaped like a token. A hit is reported by file and
#   pattern, never by its text.
# Prints the rows per arm and status and the seconds per run (the arms' own seconds summed over a run).

FORBIDDEN <- c(
    rust_location = "\\.rs:[0-9]", rust_tree = "src/rust|\\.cargo|Cargo\\.(toml|lock)", rust_build = "rustc|cargo |panicked at",
    rust_messages = "Compiling [a-z_]+ v[0-9]|error\\[E[0-9]+\\]|warning: unused",
    c_compiler = "(^|[ /])(gcc|g\\+\\+|clang|cc1|ld)( |:)|-fPIC|-std=gnu", make = "make(\\[[0-9]+\\])?: \\*\\*\\*",
    c_location = "\\.(c|cc|cpp|h|hpp|o):[0-9]+", internal_repository = "fast_limpa_dpc|fast_duckdb|internal-build-logs",
    token = "gh[pousr]_[A-Za-z0-9]{20,}|github_pat_|Bearer "
)

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 4L) stop("usage: Rscript --vanilla check-output.R <sim> <out dir> <cell> <runs>")
sim <- arguments[[1L]]; out_dir <- arguments[[2L]]; cell <- arguments[[3L]]; runs_wanted <- as.integer(arguments[[4L]])
problems <- character()

runs_path <- file.path(out_dir, sprintf("%s-runs.tsv", sim))
if (!file.exists(runs_path)) {
    problems <- c(problems, sprintf("%s was not written", basename(runs_path)))
} else {
    rows <- utils::read.delim(runs_path, colClasses = "character", quote = "")
    if (!nrow(rows) || !all(rows$cell == cell)) problems <- c(problems, sprintf("the table does not hold cell %s alone", cell))
    per_run <- table(as.integer(rows$run))
    if (!identical(sort(as.integer(names(per_run))), seq_len(runs_wanted))) {
        problems <- c(problems, sprintf("the table's runs are not 1 to %d", runs_wanted))
    }
    if (length(unique(per_run)) > 1L) problems <- c(problems, "runs differ in their number of rows")
    if (nrow(rows) && all(rows$status == "crashed")) problems <- c(problems, "every row crashed")
    cat(sprintf("%s %s: %d rows, %d runs\n", sim, cell, nrow(rows), length(per_run)))
    print(table(arm = rows$arm, status = rows$status))
    first_subset <- !duplicated(rows[c("run", "arm")])
    seconds <- tapply(as.numeric(rows$seconds[first_subset]), rows$run[first_subset], sum, na.rm = TRUE)
    cat(sprintf("seconds per run, arms only: mean %.2f, max %.2f\n", mean(seconds), max(seconds)))
}

characterValues <- function(path) {
    if (grepl("\\.rds$", path)) {
        return(unlist(rapply(list(readRDS(path)), \(x) as.character(x), classes = c("character", "factor"), how = "list"),
            use.names = FALSE))
    }
    readLines(path, warn = FALSE)
}
files <- list.files(out_dir, recursive = TRUE, full.names = TRUE)
for (path in files) {
    values <- characterValues(path)
    for (name in names(FORBIDDEN)) {
        hits <- sum(grepl(FORBIDDEN[[name]], values, useBytes = TRUE))
        if (hits) problems <- c(problems, sprintf("%s: %d value(s) match the forbidden pattern `%s`",
            substring(path, nchar(out_dir) + 2L), hits, name))
    }
}
cat(sprintf("%d files scanned for build output, internal source locations and tokens\n", length(files)))

if (length(problems)) {
    cat(sprintf("::error::output check failed: %s\n", problems), sep = "")
    quit(save = "no", status = 1L)
}
cat("output check passed\n")
