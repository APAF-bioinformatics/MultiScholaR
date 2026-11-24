# KIT-212 control C (ticket section 10; no expected value changes, nothing is re-frozen here): the W1 Summary captures
# (`omicSummaryW1Frozen`, tests/testthat/test-omic-summary-module.R) of the nine cases dumped (MULTISCHOLAR_KIT108_W1_DUMP)
# on the base commit's tree once and on the head tree twice, each in its own process. Evidence code under docs/qa: no
# package code. Runs on a GitHub-hosted runner (run-tests `scripts`) from the package root of the head commit.
#
# Usage: Rscript docs/qa/workflow-execution/kit/kit-212/control-w1.R <out dir> <base commit>
#
# A leaf is one field of one step of one case (the unit the file's test compares with expect_identical). Changed leaves:
# base dump and a head dump differ. Checks: (i) at every changed leaf the base dump equals the stored value; (ii) the two
# head dumps are identical; (iii) every changed leaf is in one of the nine cases at its Copy step, its Report step or a
# later step (from the first step whose name contains "copy" or "report"). Leaves where the base dump differs from the
# stored value and that do not change are recorded as pre-existing. Writes leaves.tsv, verdict.txt, the dumps and logs.

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
dir.create(args[[1L]], showWarnings = FALSE, recursive = TRUE)
out <- normalizePath(args[[1L]], mustWork = TRUE)
base <- args[[2L]]
CASES <- c("lipid_main", "lipid_variants", "lipid_github_fallback", "lipid_global_another", "metab_main", "metab_variants",
    "metab_github_fallback", "prot_owner", "prot_owner_variants")
work <- file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "k212-control")
src <- file.path(work, "base-src")
dir.create(src, recursive = TRUE, showWarnings = FALSE)
stopifnot(system(sprintf("git archive %s | tar -x -C %s", shQuote(base), shQuote(src))) == 0L)
rscript <- file.path(R.home("bin"), "Rscript")
FILE <- file.path("tests", "testthat", "test-omic-summary-module.R")

dump <- function(tree, label) {
    path <- file.path(out, paste0("dump-", label, ".R"))
    unlink(path)
    descs <- sprintf("(b) %s: the outputs, notifications, progress, reactive values, calls, global and files of each step are today's",
        CASES)
    code <- paste0("suppressMessages(pkgload::load_all(", deparse(tree), ", export_all = TRUE, helpers = TRUE, quiet = TRUE)); ",
        "for (d in ", paste(deparse(descs), collapse = ""), ") testthat::test_file(file.path(", deparse(tree), ", ",
        deparse(FILE), "), desc = d, reporter = 'summary', stop_on_failure = FALSE)")
    log <- file.path(out, paste0("dump-", label, ".log.txt"))
    exit <- system2(rscript, c("-e", shQuote(code)), stdout = log, stderr = log,
        env = c(paste0("MULTISCHOLAR_KIT108_W1_DUMP=", path), paste0("MULTISCHOLAR_KIT108_W1_CASES=", paste(CASES, collapse = ",")),
            "NOT_CRAN=true"))
    list(exit = exit, path = path)
}
runs <- list(base = dump(src, "base"), head_1 = dump(normalizePath("."), "head-1"), head_2 = dump(normalizePath("."), "head-2"))

evaluate <- function(code) {
    env <- new.env()
    assign("omicSummaryW1Frozen", list(), envir = env)
    ok <- tryCatch({ eval(parse(text = code), envir = env); TRUE }, error = function(e) {
        message("parse failed: ", conditionMessage(e))
        FALSE
    })
    if (ok) get("omicSummaryW1Frozen", envir = env) else NULL
}
dumped <- lapply(runs, function(run) if (file.exists(run$path)) evaluate(paste(readLines(run$path, warn = FALSE), collapse = "\n")) else NULL)
lines <- readLines(file.path(src, FILE), warn = FALSE)
at <- c(grep("^# OMIC_SUMMARY_FROZEN_START$", lines), grep("^# OMIC_SUMMARY_FROZEN_END$", lines))
stored <- evaluate(paste(lines[(at[[1L]] + 1L):(at[[2L]] - 1L)], collapse = "\n"))

short <- function(x) {
    text <- paste(deparse(x, width.cutoff = 500L), collapse = " ")
    text <- gsub("[\t\n]", " ", text)
    if (nchar(text) > 300L) paste0(substr(text, 1L, 300L), "...") else text
}
rows <- list()
for (case in CASES) {
    steps <- unique(c(names(stored[[case]]), names(dumped$base[[case]]), names(dumped$head_1[[case]])))
    first <- match(TRUE, grepl("copy|report", names(stored[[case]])))
    for (step in steps) {
        fields <- unique(c(names(stored[[case]][[step]]), names(dumped$base[[case]][[step]]), names(dumped$head_1[[case]][[step]])))
        index <- match(step, names(stored[[case]]))
        for (field in fields) {
            s <- stored[[case]][[step]][[field]]
            b <- dumped$base[[case]][[step]][[field]]
            h1 <- dumped$head_1[[case]][[step]][[field]]
            h2 <- dumped$head_2[[case]][[step]][[field]]
            changed <- !identical(b, h1) || !identical(b, h2)
            pre <- !identical(b, s)
            if (!changed && !pre) next
            rows[[length(rows) + 1L]] <- data.frame(case = case, step = step, field = field,
                stored = short(s), base_dump = short(b), head_dump = short(h1),
                kind = if (changed) "changed" else "pre-existing",
                i_base_equals_stored = if (changed) !pre else NA,
                iii_copy_report_or_later = if (changed) isTRUE(!is.na(first) && !is.na(index) && index >= first) else NA,
                stringsAsFactors = FALSE)
        }
    }
}
leaves <- if (length(rows)) do.call(rbind, rows) else data.frame(case = character(), step = character(), field = character(),
    stored = character(), base_dump = character(), head_dump = character(), kind = character(), i_base_equals_stored = logical(),
    iii_copy_report_or_later = logical())
utils::write.table(leaves, file.path(out, "leaves.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
changed <- leaves[leaves$kind == "changed", , drop = FALSE]
all_cases <- vapply(c("base", "head_1", "head_2"), function(run) all(CASES %in% names(dumped[[run]])), logical(1))
verdict <- c(
    sprintf("exits: base=%s head_1=%s head_2=%s", runs$base$exit, runs$head_1$exit, runs$head_2$exit),
    sprintf("every case dumped: base=%s head_1=%s head_2=%s", all_cases[["base"]], all_cases[["head_1"]], all_cases[["head_2"]]),
    sprintf("changed leaves: %d; pre-existing leaves: %d", nrow(changed), sum(leaves$kind == "pre-existing")),
    sprintf("(i) every changed leaf's base dump equals the stored value: %s", all(changed$i_base_equals_stored)),
    sprintf("(ii) the two head dumps are identical: %s", identical(dumped$head_1, dumped$head_2)),
    sprintf("(iii) every changed leaf is at a Copy, Report or later step of the nine cases: %s", all(changed$iii_copy_report_or_later)),
    sprintf("cases with changed leaves: %s", paste(unique(changed$case), collapse = ",")))
writeLines(verdict, file.path(out, "verdict.txt"))
writeLines(verdict)
