# KIT-205 audit K205-I01 control (no expected value changes; nothing is re-frozen): dumps of the enrichment tab's UI
# capture `omicEnrichW1FrozenUi` (test (a), MULTISCHOLAR_KIT107_W1_DUMP) and of the W113 captures (test "(KIT-113 W1)",
# MULTISCHOLAR_KIT113_W1_DUMP) of tests/testthat/test-omic-enrich-module.R. Evidence code under docs/qa: no package code.
# Runs on a GitHub-hosted runner (run-tests `scripts`) from the package root of the commit under test (head).
#
# Usage: Rscript docs/qa/workflow-execution/kit/kit-205/dumps-k205.R <out dir> <base commit>
#
# Three separate processes, each pkgload::load_all() of a source tree and the two tests of that tree's file: once on the
# base commit's tree (git archive), twice on head. Per capture (the UI block; the W113 cases `gprofiler`,
# `kernel_gprofiler`, `restore`, `organism_filter/no_mapping`, `organism_filter/mapped`) it records whether the base dump
# equals the base tree's stored block (text and value), whether the two head dumps are byte-identical, and which leaves
# differ between the base and head dumps. Writes <out dir>/dump-{ui,w113}-{base,head-1,head-2}.R, dumps.tsv and the logs.

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
dir.create(args[[1L]], showWarnings = FALSE, recursive = TRUE)
out <- normalizePath(args[[1L]], mustWork = TRUE)
base <- args[[2L]]
work <- file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "k205-dumps")
src <- file.path(work, "base-src")
dir.create(src, recursive = TRUE, showWarnings = FALSE)
rscript <- file.path(R.home("bin"), "Rscript")
archive_exit <- system(sprintf("git archive %s | tar -x -C %s", shQuote(base), shQuote(src)))
DESCS <- c("(a) the enrichment UI keeps its HTML, ids, test ids, tab titles, labels, help texts, choices and defaults",
    "(KIT-113 W1) each W1 (b) session's committed record, results, product, provenance and analysis record step")
CASES <- c("gprofiler", "kernel_gprofiler", "restore", "organism_filter/no_mapping", "organism_filter/mapped")

dump <- function(tree, label) {
    ui <- file.path(out, paste0("dump-ui-", label, ".R"))
    w113 <- file.path(out, paste0("dump-w113-", label, ".R"))
    unlink(c(ui, w113))
    code <- paste0("suppressMessages(pkgload::load_all(", deparse(tree), ", export_all = TRUE, helpers = TRUE, ",
        "quiet = TRUE)); for (d in ", paste(deparse(DESCS), collapse = ""), ") testthat::test_file(file.path(",
        deparse(tree), ", 'tests', 'testthat', 'test-omic-enrich-module.R'), desc = d, reporter = 'summary', ",
        "stop_on_failure = FALSE)")
    log <- file.path(out, paste0("dump-", label, ".log.txt"))
    exit <- system2(rscript, c("-e", shQuote(code)), stdout = log, stderr = log,
        env = c(paste0("MULTISCHOLAR_KIT107_W1_DUMP=", ui), paste0("MULTISCHOLAR_KIT113_W1_DUMP=", w113), "NOT_CRAN=true"))
    list(exit = exit, ui = ui, w113 = w113)
}
runs <- list(base = dump(src, "base"), head_1 = dump(normalizePath("."), "head-1"),
    head_2 = dump(normalizePath("."), "head-2"))

text <- function(f) if (file.exists(f)) sub("\\s+$", "", paste(readLines(f, warn = FALSE), collapse = "\n")) else ""
bytes <- function(f) if (file.exists(f)) readBin(f, "raw", file.info(f)$size) else raw()
block <- function(tree, start, end) {
    lines <- readLines(file.path(tree, "tests", "testthat", "test-omic-enrich-module.R"), warn = FALSE)
    at <- c(grep(paste0("^# ", start, "$"), lines), grep(paste0("^# ", end, "$"), lines))
    paste(lines[(at[[1L]] + 1L):(at[[2L]] - 1L)], collapse = "\n")
}
# Evaluate a block or dump and return the object it defines (`omicEnrichW1FrozenUi` or the W113 cases by key).
values <- function(code, name) {
    env <- new.env()
    if (identical(name, "omicEnrichW113Frozen")) assign(name, list(), envir = env)
    ok <- tryCatch({ eval(parse(text = code), envir = env); TRUE }, error = function(e) FALSE)
    if (ok && exists(name, envir = env)) get(name, envir = env) else NULL
}
# W113 dump text per case (the dump appends `omicEnrichW113Frozen[["key"]] <- ...` blocks)
caseText <- function(code, key) {
    parts <- strsplit(code, "\n(?=omicEnrichW113Frozen\\[\\[)", perl = TRUE)[[1L]]
    hit <- parts[startsWith(parts, sprintf("omicEnrichW113Frozen[[\"%s\"]] <- ", key))]
    if (length(hit)) sub("\\s+$", "", hit[[1L]]) else ""
}
leaves <- function(x, path = character()) {
    if (is.list(x) && length(x)) {
        keys <- names(x) %||% rep("", length(x))
        keys[!nzchar(keys)] <- paste0("[", which(!nzchar(keys)), "]")
        return(unlist(lapply(seq_along(x), function(i) leaves(x[[i]], c(path, keys[[i]])))))
    }
    if (!is.list(x) && length(x) > 1L && !is.null(names(x))) {
        return(unlist(lapply(seq_along(x), function(i) leaves(unname(x[[i]]), c(path, names(x)[[i]])))))
    }
    stats::setNames(digest::digest(x, algo = "sha256"), paste(path, collapse = "/"))
}
differ <- function(a, b) {
    la <- leaves(a)
    lb <- leaves(b)
    keys <- sort(union(names(la), names(lb)))
    keys[vapply(keys, function(k) !identical(la[k][[1L]], lb[k][[1L]]), logical(1))]
}
rows <- list()
ui_base_stored <- block(src, "OMIC_ENRICH_FROZEN_UI_START", "OMIC_ENRICH_FROZEN_UI_END")
ui_b <- text(runs$base$ui)
rows[[1L]] <- data.frame(capture = "omicEnrichW1FrozenUi",
    base_dump_equals_stored_text = nzchar(ui_b) && identical(ui_b, ui_base_stored),
    base_dump_equals_stored_value = identical(values(ui_b, "omicEnrichW1FrozenUi"),
        values(ui_base_stored, "omicEnrichW1FrozenUi")),
    head_dumps_byte_identical = length(bytes(runs$head_1$ui)) > 0L && identical(bytes(runs$head_1$ui), bytes(runs$head_2$ui)),
    head_sha256 = if (file.exists(runs$head_1$ui)) digest::digest(file = runs$head_1$ui, algo = "sha256") else NA,
    leaves_differing_base_to_head = paste(differ(values(ui_b, "omicEnrichW1FrozenUi"),
        values(text(runs$head_1$ui), "omicEnrichW1FrozenUi")), collapse = "; "), stringsAsFactors = FALSE)
w113_stored <- block(src, "OMIC_ENRICH_W113_FROZEN_START", "OMIC_ENRICH_W113_FROZEN_END")
stored_values <- values(w113_stored, "omicEnrichW113Frozen")
base_values <- values(text(runs$base$w113), "omicEnrichW113Frozen")
head_values <- values(text(runs$head_1$w113), "omicEnrichW113Frozen")
for (key in CASES) {
    rows[[length(rows) + 1L]] <- data.frame(capture = paste0("omicEnrichW113Frozen[[\"", key, "\"]]"),
        base_dump_equals_stored_text = nzchar(caseText(text(runs$base$w113), key)) &&
            identical(caseText(text(runs$base$w113), key), caseText(w113_stored, key)),
        base_dump_equals_stored_value = !is.null(base_values[[key]]) && identical(base_values[[key]], stored_values[[key]]),
        head_dumps_byte_identical = nzchar(caseText(text(runs$head_1$w113), key)) &&
            identical(caseText(text(runs$head_1$w113), key), caseText(text(runs$head_2$w113), key)),
        head_sha256 = digest::digest(caseText(text(runs$head_1$w113), key), algo = "sha256"),
        leaves_differing_base_to_head = paste(differ(base_values[[key]], head_values[[key]]), collapse = "; "),
        stringsAsFactors = FALSE)
}
table <- do.call(rbind, rows)
table$note <- c(sprintf("archive %d; exits base %d head-1 %d head-2 %d; W113 dump files byte-identical: %s", archive_exit,
    runs$base$exit, runs$head_1$exit, runs$head_2$exit, identical(bytes(runs$head_1$w113), bytes(runs$head_2$w113))),
    rep("", nrow(table) - 1L))
utils::write.table(table, file.path(out, "dumps.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(table)
