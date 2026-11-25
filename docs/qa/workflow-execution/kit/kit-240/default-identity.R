# KIT-240 (row I of the ticket, section 5): the integration tabs are on by the package default. Proves byte for byte:
#   (i)   head with `MultiScholaR.integration_gui` unset renders what base renders with it TRUE, for `app_ui(NULL)` and
#         `mod_integration_ui("integration_workflow")`, in production and test mode; the only differences allowed are the
#         reworded notes and marks of KIT-240 item 8 in the module (VAL-221 (4)): base's strings of the wording keys
#         `settings_provisional` and `results_provisional`, and base's "(<mark_provisional>)", replaced by head's
#         strings of the same keys and head's "(<mark_decided>)" (the mark every Model and Results setting carries at
#         head). The pairs are read from the two wording tables, never typed here; each is compared in its HTML-escaped
#         form; each old string must occur at least once in base's render and never in head's (counts printed and
#         checked; a pair that finds nothing is an instrument failure). `app_ui(NULL)` takes no substitution;
#   (ii)  head with the option FALSE renders exactly what base renders with it unset (no substitution), and the head's
#         module is then the "not yet implemented" panel;
#   (iii) the head renders the same twice, in two R processes.
# Each side is loaded with pkgload::load_all in its own R process (base from a git worktree of <base>, never a stash),
# as kit-138's gate-off-identity.R does. shiny's tabset ids are random draws from shiny's private seed; it is set to the
# same value before every render on both sides, so that equal pages are equal bytes.
#
# Usage (from the package root; full git history): Rscript default-identity.R <base commit> <out dir>
# Writes <out dir>/default-identity.tsv (side, option, mode, object, sha256, bytes), substitutions.tsv (pair, old, new,
# count in base, count in head), the renders, and default-identity-verdict.txt; exits 1 when a check fails.

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 2L)
base <- arguments[[1L]]
out_dir <- arguments[[2L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir)
head <- normalizePath(".")
tree <- file.path(tempdir(), "default-identity-base")
stopifnot(system2("git", c("worktree", "add", "--detach", tree, base)) == 0L)

child <- '
arguments <- commandArgs(trailingOnly = TRUE)
package <- arguments[[1L]]
out_dir <- arguments[[2L]]
side <- arguments[[3L]]
options_wanted <- strsplit(arguments[[5L]], ",", fixed = TRUE)[[1L]]
suppressMessages(pkgload::load_all(package, quiet = TRUE, export_all = TRUE, helpers = FALSE))
globals <- get(".globals", envir = asNamespace("shiny"))
render <- function(tag) {
    set.seed(20261125L)
    assign("ownSeed", get(".Random.seed", envir = globalenv()), envir = globals)
    as.character(htmltools::renderTags(tag)$html)
}
unavailable <- shiny::div(class = "alert alert-warning", shiny::h4("Module not yet implemented"),
    shiny::p(paste("The", "integration", "workflow module is under development.")))
rows <- list()
for (option in options_wanted) {
    options(MultiScholaR.integration_gui = switch(option, unset = NULL, true = TRUE, false = FALSE))
    for (mode in c("production", "test")) {
        options(multischolar.test_mode = identical(mode, "test"))
        pages <- list(app_ui = render(app_ui(NULL)), module = render(mod_integration_ui("integration_workflow")))
        for (object in names(pages)) {
            page <- pages[[object]]
            writeLines(page, file.path(out_dir, sprintf("%s-%s-%s-%s.html", side, option, mode, object)), useBytes = TRUE)
            rows[[length(rows) + 1L]] <- data.frame(side = side, option = option, mode = mode, object = object,
                sha256 = digest::digest(page, algo = "sha256", serialize = FALSE), bytes = nchar(page, type = "bytes"),
                module_is_unavailable_panel = identical(mod_integration_ui("integration_workflow"), unavailable),
                stringsAsFactors = FALSE)
        }
    }
}
utils::write.table(do.call(rbind, rows), arguments[[4L]], sep = "\\t", quote = FALSE, row.names = FALSE)
texts <- integrationGuiTexts()
keys <- c("settings_provisional", "results_provisional", "mark_provisional", "mark_decided")
writeLines(jsonlite::toJSON(as.list(texts[intersect(keys, names(texts))]), auto_unbox = TRUE),
    file.path(out_dir, sprintf("%s-texts.json", side)), useBytes = TRUE)
'
script <- file.path(tempdir(), "default-identity-child.R")
writeLines(child, script)

runs <- list(base = list(package = tree, options = "true,unset"), head = list(package = head, options = "unset,false"),
    `head-again` = list(package = head, options = "unset"))
results <- list()
for (side in names(runs)) {
    target <- file.path(tempdir(), sprintf("%s.tsv", side))
    log <- file.path(out_dir, sprintf("%s.log", side))
    status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", script, runs[[side]]$package, out_dir, side, target,
        runs[[side]]$options), stdout = log, stderr = log)
    if (status != 0L || !file.exists(target)) stop(sprintf("the %s render failed (status %d); see %s.log", side, status, side))
    results[[side]] <- utils::read.delim(target, stringsAsFactors = FALSE)
}
rows <- do.call(rbind, results)
rownames(rows) <- NULL
utils::write.table(rows, file.path(out_dir, "default-identity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
page <- function(side, option, mode, object) {
    path <- file.path(out_dir, sprintf("%s-%s-%s-%s.html", side, option, mode, object))
    paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}
sha <- function(x) digest::digest(x, algo = "sha256", serialize = FALSE)
count <- function(text, pattern) {
    at <- gregexpr(pattern, text, fixed = TRUE)[[1L]]
    if (identical(at[[1L]], -1L)) 0L else length(at)
}

# The substitution pairs of (i), from the two wording tables.
base_texts <- jsonlite::read_json(file.path(out_dir, "base-texts.json"))
head_texts <- jsonlite::read_json(file.path(out_dir, "head-texts.json"))
pairs <- data.frame(pair = c("settings_provisional", "results_provisional", "mark"),
    old = c(base_texts$settings_provisional, base_texts$results_provisional, sprintf("(%s)", base_texts$mark_provisional)),
    new = c(head_texts$settings_provisional, head_texts$results_provisional, sprintf("(%s)", head_texts$mark_decided)),
    stringsAsFactors = FALSE)
pairs$old <- htmltools::htmlEscape(pairs$old)
pairs$new <- htmltools::htmlEscape(pairs$new)

checks <- character()
check <- function(name, ok) {
    checks[[name]] <<- if (isTRUE(ok)) "TRUE" else "FALSE"
    invisible(ok)
}
substitutions <- list()
for (mode in c("production", "test")) {
    # (i) app_ui: no substitution.
    check(sprintf("(i) app_ui, %s: head unset = base TRUE", mode),
        identical(sha(page("head", "unset", mode, "app_ui")), sha(page("base", "true", mode, "app_ui"))))
    # (i) the module: base's reworded strings replaced by head's, each counted first.
    base_module <- page("base", "true", mode, "module")
    head_module <- page("head", "unset", mode, "module")
    for (k in seq_len(nrow(pairs))) {
        in_base <- count(base_module, pairs$old[[k]])
        in_head <- count(head_module, pairs$old[[k]])
        substitutions[[length(substitutions) + 1L]] <- data.frame(mode = mode, pair = pairs$pair[[k]], old = pairs$old[[k]],
            new = pairs$new[[k]], count_in_base = in_base, count_in_head = in_head, stringsAsFactors = FALSE)
        check(sprintf("(i) module, %s: %s found in base (%d) and not in head (%d)", mode, pairs$pair[[k]], in_base, in_head),
            in_base >= 1L && in_head == 0L)
        base_module <- gsub(pairs$old[[k]], pairs$new[[k]], base_module, fixed = TRUE)
    }
    check(sprintf("(i) module, %s: head unset = base TRUE after the substitutions", mode),
        identical(sha(head_module), sha(base_module)))
    if (!identical(head_module, base_module)) {
        writeLines(base_module, file.path(out_dir, sprintf("base-true-%s-module-substituted.html", mode)), useBytes = TRUE)
        a <- strsplit(base_module, "")[[1L]]
        b <- strsplit(head_module, "")[[1L]]
        n <- min(length(a), length(b))
        first <- which(a[seq_len(n)] != b[seq_len(n)])[1L]
        if (is.na(first)) first <- n + 1L
        message(sprintf("(i) %s: first difference at character %d; base: %s | head: %s", mode, first,
            paste(a[seq(first, min(length(a), first + 120L))], collapse = ""),
            paste(b[seq(first, min(length(b), first + 120L))], collapse = "")))
    }
    # (ii) head FALSE = base unset, no substitution.
    for (object in c("app_ui", "module")) {
        check(sprintf("(ii) %s, %s: head FALSE = base unset", object, mode),
            identical(sha(page("head", "false", mode, object)), sha(page("base", "unset", mode, object))))
    }
    # (iii) the head twice.
    for (object in c("app_ui", "module")) {
        check(sprintf("(iii) %s, %s: head unset the same in two processes", object, mode),
            identical(sha(page("head", "unset", mode, object)), sha(page("head-again", "unset", mode, object))))
    }
}
check("(ii) head FALSE: the module is the unavailable panel",
    all(rows$module_is_unavailable_panel[rows$side == "head" & rows$option == "false"]))
check("head unset: the module is not the unavailable panel",
    !any(rows$module_is_unavailable_panel[rows$side %in% c("head", "head-again") & rows$option == "unset"]))
substitutions <- do.call(rbind, substitutions)
utils::write.table(substitutions, file.path(out_dir, "substitutions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

verdict <- c(paste("base:", base), paste("head:", system2("git", c("rev-parse", "HEAD"), stdout = TRUE)), R.version.string,
    "substitutions (old -> new, count in base, count in head):",
    sprintf("  %s %s: %s -> %s (%d, %d)", substitutions$mode, substitutions$pair, substitutions$old, substitutions$new,
        substitutions$count_in_base, substitutions$count_in_head),
    sprintf("%s: %s", names(checks), checks), paste("all checks pass:", all(checks == "TRUE")))
writeLines(verdict, file.path(out_dir, "default-identity-verdict.txt"))
writeLines(verdict)
system2("git", c("worktree", "remove", "--force", tree))
if (!all(checks == "TRUE")) quit(status = 1L)
