# KIT-138: with the integration gate off (`MultiScholaR.integration_gui` unset, the default) the app's page is
# byte-identical to the base commit's. The base tree and the head are each loaded in their own R process
# (pkgload::load_all, as the runner's per-file runs load the package) and `app_ui(NULL)` rendered in production mode and
# in test mode (`multischolar.test_mode`), twice on the head (a page that differs between two renders of one tree could
# not be compared). The integration module's UI with the gate off is compared with the panel app_server showed for an
# omic without a module.
#
# Usage (from the package root; full git history): Rscript gate-off-identity.R <base commit> <out dir>
# Writes <out dir>/gate-off-identity.tsv (side, mode, render, sha256, bytes) and gate-off-verdict.txt; exits 1 when a
# gate-off page differs.

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 2L)
base <- arguments[[1L]]
out_dir <- arguments[[2L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir)
head <- normalizePath(".")
tree <- file.path(tempdir(), "gate-off-base")
stopifnot(system2("git", c("worktree", "add", "--detach", tree, base)) == 0L)

child <- '
arguments <- commandArgs(trailingOnly = TRUE)
suppressMessages(pkgload::load_all(arguments[[1L]], quiet = TRUE, export_all = TRUE, helpers = FALSE))
options(MultiScholaR.integration_gui = NULL)
rows <- list()
for (mode in c("production", "test")) {
    options(multischolar.test_mode = identical(mode, "test"))
    page <- as.character(htmltools::renderTags(app_ui(NULL))$html)
    file <- file.path(arguments[[2L]], sprintf("%s-%s.html", arguments[[3L]], mode))
    writeLines(page, file, useBytes = TRUE)
    rows[[mode]] <- data.frame(side = arguments[[3L]], mode = mode, sha256 = digest::digest(page, algo = "sha256", serialize = FALSE),
        bytes = nchar(page, type = "bytes"), stringsAsFactors = FALSE)
}
module <- if (exists("mod_integration_ui")) {
    identical(mod_integration_ui("integration_workflow"), shiny::div(class = "alert alert-warning",
        shiny::h4("Module not yet implemented"), shiny::p(paste("The", "integration", "workflow module is under development."))))
} else NA
rows <- do.call(rbind, rows)
rows$module_is_unavailable_panel <- module
utils::write.table(rows, arguments[[4L]], sep = "\\t", quote = FALSE, row.names = FALSE)
'
script <- file.path(tempdir(), "gate-off-child.R")
writeLines(child, script)

results <- list()
for (side in c("base", "head", "head-again")) {
    package <- if (side == "base") tree else head
    target <- file.path(tempdir(), sprintf("%s.tsv", side))
    status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", script, package, out_dir, side, target),
        stdout = file.path(out_dir, sprintf("%s.log", side)), stderr = file.path(out_dir, sprintf("%s.log", side)))
    if (status != 0L || !file.exists(target)) stop(sprintf("the %s render failed (status %d); see %s.log", side, status, side))
    results[[side]] <- utils::read.delim(target, stringsAsFactors = FALSE)
}
rows <- do.call(rbind, results)
utils::write.table(rows, file.path(out_dir, "gate-off-identity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
sha <- function(side, mode) rows$sha256[rows$side == side & rows$mode == mode]
deterministic <- all(vapply(c("production", "test"), function(mode) identical(sha("head", mode), sha("head-again", mode)), logical(1)))
identical_pages <- all(vapply(c("production", "test"), function(mode) identical(sha("base", mode), sha("head", mode)), logical(1)))
module <- isTRUE(all(rows$module_is_unavailable_panel[rows$side != "base"]))
verdict <- c(paste("base:", base), paste("head:", system2("git", c("rev-parse", "HEAD"), stdout = TRUE)), R.version.string,
    paste("head renders deterministic:", deterministic), paste("gate-off pages identical to base:", identical_pages),
    paste("gate-off integration module is the unavailable panel:", module))
writeLines(verdict, file.path(out_dir, "gate-off-verdict.txt"))
writeLines(verdict)
system2("git", c("worktree", "remove", "--force", tree))
if (!(deterministic && identical_pages && module)) quit(status = 1L)
