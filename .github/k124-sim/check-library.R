# KIT-124 on hosted runners: the runner's library against the laptop's, before a simulator starts.
#
# Usage: Rscript --vanilla check-library.R <laptop-versions.tsv> <out.tsv>
#
# <laptop-versions.tsv> (package, version, library, check) lists the laptop's versions of the packages the simulators
# load and their dependencies, as installed in the library paths the laptop runs used. Every listed package is looked
# up in this R's library. Writes <out.tsv> (package, laptop, runner, check, verdict) and prints the rows that differ.
# The exit status is 1 when R is not 4.6.1, when a package marked `hard` is missing or at another version, or when an
# internal or GitHub package was not built from the pinned commit; other differences are reported only.

R_WANTED <- "4.6.1"
PINNED_SHA <- c(
    fastlimpadpc = "87a90b894b2057744fc9c315a4fb7d7db3ce3af7",
    fastduckdb = "e724f90bdbdf1073635b7112c80f8dd406f2c79b",
    RUVIIIC = "c534151f6b2265b36294ff630f4dd13c2e86a9fa"
)

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) stop("usage: Rscript --vanilla check-library.R <laptop-versions.tsv> <out.tsv>")
laptop <- utils::read.delim(arguments[[1L]], stringsAsFactors = FALSE)
installed <- utils::installed.packages(lib.loc = .Library, fields = "RemoteSha")

runner <- unname(installed[match(laptop$package, rownames(installed)), "Version"])
same <- !is.na(runner) & runner == laptop$version
result <- data.frame(package = laptop$package, laptop = laptop$version, runner = runner, check = laptop$check,
    verdict = ifelse(same, "same", ifelse(is.na(runner), "missing", "differs")), stringsAsFactors = FALSE)

failures <- character()
if (!identical(as.character(getRversion()), R_WANTED)) failures <- c(failures, sprintf("R is %s, not %s", getRversion(), R_WANTED))
hard_off <- result[result$check == "hard" & result$verdict != "same", ]
if (nrow(hard_off)) {
    failures <- c(failures, sprintf("%s: laptop %s, runner %s", hard_off$package, hard_off$laptop, hard_off$runner))
}
for (package in names(PINNED_SHA)) {
    sha <- if (package %in% rownames(installed)) unname(installed[package, "RemoteSha"]) else NA_character_
    if (!identical(sha, PINNED_SHA[[package]])) {
        failures <- c(failures, sprintf("%s: built from %s, not %s", package, sha, PINNED_SHA[[package]]))
    }
    result <- rbind(result, data.frame(package = paste0(package, " (commit)"), laptop = "", runner = sha,
        check = "hard", verdict = if (identical(sha, PINNED_SHA[[package]])) "pinned" else "differs"))
}
extra <- data.frame(package = c("R", "BLAS", "LAPACK"), laptop = c(R_WANTED, "", ""),
    runner = c(as.character(getRversion()), basename(extSoftVersion()[["BLAS"]]), basename(La_library())),
    check = c("hard", "report", "report"), verdict = "", stringsAsFactors = FALSE)
extra$verdict[[1L]] <- if (identical(extra$runner[[1L]], R_WANTED)) "same" else "differs"
result <- rbind(extra, result)

dir.create(dirname(arguments[[2L]]), recursive = TRUE, showWarnings = FALSE)
utils::write.table(result, arguments[[2L]], sep = "\t", quote = FALSE, row.names = FALSE, na = "")
options(width = 200)
cat(sprintf("%d packages listed: %d same, %d differ, %d missing\n", nrow(laptop), sum(same),
    sum(!same & !is.na(runner)), sum(is.na(runner))))
print(result[result$check == "hard" | !result$verdict %in% c("same", ""), ], row.names = FALSE)
if (length(failures)) {
    cat(sprintf("::error::library check failed: %s\n", failures), sep = "")
    quit(save = "no", status = 1L)
}
cat("library check passed: R, the hard-checked packages and the pinned commits are the laptop's\n")
