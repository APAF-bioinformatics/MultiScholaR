# KIT-213 control C (ticket 10; A3; stop rule 3), KIT-203 OD5's two-dump form. Evidence code under docs/qa: no package
# code. Runs on a GitHub-hosted runner (run-tests `scripts`) from the package root of the commit under test (head).
#
# Usage: Rscript docs/qa/workflow-execution/kit/kit-213/control-w1.R <out dir> <base commit> <case,case,...>
#
# Three separate processes, each pkgload::load_all() of a source tree and, per listed case, the case's W1 test of that
# tree's tests/testthat/test-omic-summary-module.R (testthat::test_file, `desc` the case's test, MULTISCHOLAR_KIT108_W1_CASES
# the case, MULTISCHOLAR_KIT108_W1_DUMP the side's dump): once on the base commit's tree (git archive), twice on head.
# In each process the package's createWorkflowArgsFromConfig is wrapped (the test's recorder calls the wrapper, which
# calls the real writer) to keep every study_parameters.txt the case writes, masked as the W1 capture masks it
# (omicSummaryW1MaskVersion, then the entry suite's guiEntryMask). The script then checks:
#   text: base against head, only L lines (`pairs`, `logcpm`, `grid`, `results` [input]) differ; head against head,
#         identical;
#   (i)   at every changed leaf (base dump against head dump), the base dump equals the stored omicSummaryW1Frozen value
#         (the base tree's frozen block);
#   (ii)  the two head dumps are byte-identical;
#   (iii) changed leaves fall only among the `study_parameters` leaves of the listed cases; nothing else moves.
# A leaf whose base dump already differs from the stored value is written as pre-existing (left untouched unless it is
# also a changed leaf, which fails (i)). Writes <out dir>/leaves.tsv (case, step, field, stored, base dump, head dump,
# status), text.tsv, verdict.txt, the dumps and the side logs. Exit status 0 (a difference is a result: the verdict
# names it).

args <- commandArgs(trailingOnly = TRUE)
L_PATTERN <- "^  ([^ ]+\\.)?(pairs|logcpm|grid|results) \\[input\\] = "
W1_DESC <- "(b) %s: the outputs, notifications, progress, reactive values, calls, global and files of each step are today's"

# --- one side (a separate process) -----------------------------------------------------------------------------------
if (length(args) && identical(args[[1L]], "--side")) {
    tree <- args[[2L]]
    side_out <- args[[3L]]
    cases <- strsplit(args[[4L]], ",", fixed = TRUE)[[1L]]
    dump <- args[[5L]]
    dir.create(side_out, recursive = TRUE, showWarnings = FALSE)
    suppressMessages(pkgload::load_all(tree, export_all = TRUE, helpers = TRUE, quiet = TRUE))
    real <- get("createWorkflowArgsFromConfig", envir = asNamespace("MultiScholaR"))
    mask <- get("guiEntryMask", envir = asNamespace("MultiScholaR"), inherits = TRUE)
    written <- new.env()
    written$case <- ""
    written$n <- 0L
    testthat::local_mocked_bindings(createWorkflowArgsFromConfig = function(...) {
        path <- real(...)
        written$n <- written$n + 1L
        lines <- readLines(path, warn = FALSE)
        start <- which(lines == "Version Information:")
        if (length(start)) {
            start <- start[[1L]]
            blank <- which(!nzchar(lines) & seq_along(lines) > start)
            end <- if (length(blank)) blank[[1L]] - 1L else length(lines)
            lines <- c(lines[seq_len(start - 1L)], "Version Information: <masked>", lines[-seq_len(end)])
        }
        masked <- strsplit(mask(lines), "\n", fixed = TRUE)[[1L]]
        writeLines(masked, file.path(side_out, sprintf("%s-%02d.txt", written$case, written$n)))
        path
    }, .package = "MultiScholaR", .env = globalenv())
    for (case in cases) {
        written$case <- case
        written$n <- 0L
        Sys.setenv(MULTISCHOLAR_KIT108_W1_CASES = case, MULTISCHOLAR_KIT108_W1_DUMP = dump, NOT_CRAN = "true")
        testthat::test_file(file.path(tree, "tests", "testthat", "test-omic-summary-module.R"), desc = sprintf(W1_DESC, case),
            reporter = "summary", stop_on_failure = FALSE)
    }
    quit(save = "no", status = 0L)
}

# --- the control -----------------------------------------------------------------------------------------------------
stopifnot(length(args) == 3L)
dir.create(args[[1L]], showWarnings = FALSE, recursive = TRUE)
out <- normalizePath(args[[1L]], mustWork = TRUE)
base <- args[[2L]]
cases <- strsplit(args[[3L]], ",", fixed = TRUE)[[1L]]
work <- file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "k213-control")
src <- file.path(work, "base-src")
dir.create(src, recursive = TRUE, showWarnings = FALSE)
stopifnot(system(sprintf("git archive %s | tar -x -C %s", shQuote(base), shQuote(src))) == 0L)
rscript <- file.path(R.home("bin"), "Rscript")
self <- normalizePath("docs/qa/workflow-execution/kit/kit-213/control-w1.R", mustWork = TRUE)

side <- function(tree, label) {
    dump <- file.path(out, paste0("dump-", label, ".R"))
    unlink(dump)
    log <- file.path(out, paste0("side-", label, ".log.txt"))
    exit <- system2(rscript, c(shQuote(self), "--side", shQuote(tree), shQuote(file.path(out, paste0("text-", label))),
        shQuote(paste(cases, collapse = ",")), shQuote(dump)), stdout = log, stderr = log)
    list(exit = exit, dump = dump, text = file.path(out, paste0("text-", label)))
}
runs <- list(base = side(src, "base"), head_1 = side(normalizePath("."), "head-1"), head_2 = side(normalizePath("."), "head-2"))

# The captures of a dump, or of the base tree's frozen block, by case.
frozen <- function(code) {
    env <- new.env()
    assign("omicSummaryW1Frozen", list(), envir = env)
    ok <- tryCatch({ eval(parse(text = code), envir = env); TRUE }, error = function(e) conditionMessage(e))
    if (!isTRUE(ok)) stop("unparsable captures: ", ok)
    get("omicSummaryW1Frozen", envir = env)
}
fileText <- function(f) if (file.exists(f)) paste(readLines(f, warn = FALSE), collapse = "\n") else ""
stored_lines <- readLines(file.path(src, "tests", "testthat", "test-omic-summary-module.R"), warn = FALSE)
at <- c(grep("^# OMIC_SUMMARY_FROZEN_START$", stored_lines), grep("^# OMIC_SUMMARY_FROZEN_END$", stored_lines))
stored <- frozen(paste(stored_lines[(at[[1L]] + 1L):(at[[2L]] - 1L)], collapse = "\n"))
dumps <- lapply(runs, function(r) frozen(fileText(r$dump)))

show <- function(x) if (is.null(x)) "<missing>" else gsub("[\t\n]", " ", paste(deparse(x, width.cutoff = 500L), collapse = " "))
rows <- list()
for (case in cases) {
    steps <- unique(c(names(stored[[case]]), names(dumps$base[[case]]), names(dumps$head_1[[case]])))
    for (step in steps) {
        fields <- unique(c(names(stored[[case]][[step]]), names(dumps$base[[case]][[step]]), names(dumps$head_1[[case]][[step]])))
        for (field in fields) {
            s <- stored[[case]][[step]][[field]]
            b <- dumps$base[[case]][[step]][[field]]
            h <- dumps$head_1[[case]][[step]][[field]]
            changed <- !identical(b, h)
            pre <- !identical(b, s)
            if (!changed && !pre) next
            rows[[length(rows) + 1L]] <- data.frame(case = case, step = step, field = field, stored = show(s), base_dump = show(b),
                head_dump = show(h), status = if (changed && pre) "changed+pre-existing" else if (changed) "changed" else "pre-existing")
        }
    }
}
leaves <- if (length(rows)) do.call(rbind, rows) else data.frame(case = character(), step = character(), field = character(),
    stored = character(), base_dump = character(), head_dump = character(), status = character())
utils::write.table(leaves, file.path(out, "leaves.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# The masked study parameters, side by side.
texts <- function(dir) { f <- list.files(dir, pattern = "\\.txt$"); stats::setNames(lapply(file.path(dir, f), readLines, warn = FALSE), f) }
tb <- texts(runs$base$text)
th1 <- texts(runs$head_1$text)
th2 <- texts(runs$head_2$text)
text_rows <- list()
for (name in union(names(tb), names(th1))) {
    a <- tb[[name]] %||% character()
    b <- th1[[name]] %||% character()
    n <- max(length(a), length(b))
    a <- c(a, rep("<no line>", n - length(a)))
    b <- c(b, rep("<no line>", n - length(b)))
    for (i in which(a != b)) text_rows[[length(text_rows) + 1L]] <- data.frame(file = name, line = i,
        l_line = grepl(L_PATTERN, a[[i]]) && grepl(L_PATTERN, b[[i]]), base = a[[i]], head = b[[i]])
}
text_table <- if (length(text_rows)) do.call(rbind, text_rows) else data.frame(file = character(), line = integer(),
    l_line = logical(), base = character(), head = character())
utils::write.table(text_table, file.path(out, "text.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

changed <- leaves[grepl("changed", leaves$status), , drop = FALSE]
verdict <- c(
    sprintf("sides: exit base=%s head-1=%s head-2=%s; files base=%d head-1=%d head-2=%d", runs$base$exit, runs$head_1$exit,
        runs$head_2$exit, length(tb), length(th1), length(th2)),
    sprintf("cases run: %s; in each dump: base=%s head-1=%s head-2=%s", paste(cases, collapse = ","),
        all(cases %in% names(dumps$base)), all(cases %in% names(dumps$head_1)), all(cases %in% names(dumps$head_2))),
    sprintf("text: base vs head differing lines %d, all L lines: %s; head-1 vs head-2 identical: %s", nrow(text_table),
        all(text_table$l_line), identical(th1, th2)),
    sprintf("(i) base dump equals stored at every changed leaf: %s (%d changed leaves)",
        !any(changed$status == "changed+pre-existing"), nrow(changed)),
    sprintf("(ii) the two head dumps byte-identical: %s", identical(fileText(runs$head_1$dump), fileText(runs$head_2$dump)) &&
        nzchar(fileText(runs$head_1$dump))),
    sprintf("(iii) changed leaves only study_parameters leaves of the listed cases: %s", all(changed$field == "study_parameters")),
    sprintf("pre-existing leaves (base dump differs from stored, not changed): %d", sum(leaves$status == "pre-existing")))
writeLines(verdict, file.path(out, "verdict.txt"))
cat(verdict, sep = "\n")
