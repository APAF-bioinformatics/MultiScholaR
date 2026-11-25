# KIT-167 section 8.8 (A11): with the reference mapping gate off (`MultiScholaR.reference_mapping` unset, the default)
# the pages and module configurations the ticket changes are identical to the base commit's. The kit-138 pattern: the
# base tree and the head each loaded in their own R process (pkgload::load_all), the head twice (a value that differs
# between two loads of one tree could not be compared).
#
# Hashed per side: app_ui(NULL); mod_prot_import_ui("x"); mod_omic_import_ui("x", transcriptomics configuration);
# mod_prot_enrich_ui("x"); mod_omic_enrich_ui("x", the proteomics configuration's enrichment group); omicGuiConfig() of proteomics and of
# transcriptomics, as values (every function replaced by a marker) and, separately, with every function deparsed. The
# verdict reads the pages and the values; the deparsed functions are reported, not judged: the ticket changes the bodies
# of functions the proteomics configuration holds (the annotation resolver, 8.3 item 3; the enrichment declaration,
# 8.3 item 5), whatever the gate says, so their deparsed text differs by design. A UI function that does not exist or fails is hashed as its error text
# (the same on both sides when neither changed). Table check: the proteomics FASTA step's file is unchanged between
# base and head (`aa_seq_tbl_final` is built there); the enrichment identity without the choice and the integration
# features without a mapping are tests of test-reference-mapping-enrichment.R.
#
# Usage (from the package root; full git history): Rscript gate-off-identity.R <base commit> <out dir>
# Writes <out dir>/gate-off-identity.tsv (side, item, sha256) and gate-off-verdict.txt; exits 1 when an item differs.

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 2L)
base <- arguments[[1L]]
out_dir <- arguments[[2L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir)
head <- normalizePath(".")
tree <- file.path(tempdir(), "k167-gate-off-base")
stopifnot(system2("git", c("worktree", "add", "--detach", tree, base)) == 0L)

child <- '
arguments <- commandArgs(trailingOnly = TRUE)
suppressMessages(pkgload::load_all(arguments[[1L]], quiet = TRUE, export_all = TRUE, helpers = FALSE))
options(MultiScholaR.reference_mapping = NULL)
# An item whose rendering fails is recorded as its error text, so the verdict counts it (run 45347302450 hashed the
# error text of the enrichment page, called with the whole configuration instead of its enrichment group).
sha <- function(text) if (length(text) && startsWith(text[[1L]], "error:")) gsub("[\\t\\n]", " ", text[[1L]]) else
    digest::digest(paste(text, collapse = "\\n"), algo = "sha256", serialize = FALSE)
# Shiny gives each tabset a random id (`data-tabsetid`, `tab-<n>-<m>`; run 45339028522 found mod_prot_enrich_ui
# different between two renders of the head): those numbers are masked before hashing.
page <- function(expr) {
    html <- tryCatch(as.character(htmltools::renderTags(expr)$html), error = function(e) paste("error:", conditionMessage(e)))
    html <- gsub("data-tabsetid=.[0-9]+.", "data-tabsetid=ID", html)
    gsub("tab-[0-9]+-", "tab-ID-", html)
}
plain <- function(value, functions) {
    if (is.function(value)) return(if (functions) paste(deparse(value), collapse = "\\n") else "<function>")
    if (is.list(value)) return(lapply(value, plain, functions = functions))
    value
}
config <- function(omic, functions = FALSE) tryCatch(digest::digest(plain(omicGuiConfig(omic), functions),
    algo = "sha256"), error = function(e) paste("error:", conditionMessage(e)))
items <- list(
    app_ui = sha(page(app_ui(NULL))),
    mod_prot_import_ui = sha(page(mod_prot_import_ui("x"))),
    mod_omic_import_ui_transcriptomics = sha(page(mod_omic_import_ui("x", omicGuiConfig("transcriptomics")))),
    mod_prot_enrich_ui = sha(page(mod_prot_enrich_ui("x"))),
    mod_omic_enrich_ui_proteomics = sha(page(mod_omic_enrich_ui("x", omicGuiConfig("proteomics")$enrichment))),
    omic_gui_config_proteomics = config("proteomics"),
    omic_gui_config_transcriptomics = config("transcriptomics"),
    functions_omic_gui_config_proteomics = config("proteomics", TRUE),
    functions_omic_gui_config_transcriptomics = config("transcriptomics", TRUE))
rows <- data.frame(side = arguments[[2L]], item = names(items), sha256 = unlist(items), stringsAsFactors = FALSE)
utils::write.table(rows, arguments[[3L]], sep = "\\t", quote = FALSE, row.names = FALSE)
'
script <- file.path(tempdir(), "k167-gate-off-child.R")
writeLines(child, script)

results <- list()
for (side in c("base", "head", "head-again")) {
    package <- if (side == "base") tree else head
    target <- file.path(tempdir(), sprintf("k167-%s.tsv", side))
    status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", script, package, side, target),
        stdout = file.path(out_dir, sprintf("%s.log", side)), stderr = file.path(out_dir, sprintf("%s.log", side)))
    if (status != 0L || !file.exists(target)) stop(sprintf("the %s render failed (status %d); see %s.log", side, status, side))
    results[[side]] <- utils::read.delim(target, stringsAsFactors = FALSE)
}
rows <- do.call(rbind, results)
utils::write.table(rows, file.path(out_dir, "gate-off-identity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
value <- function(side, judged = TRUE) {
    keep <- rows$side == side & (startsWith(rows$item, "functions_") != judged)
    stats::setNames(rows$sha256[keep], rows$item[keep])
}
deterministic <- identical(value("head"), value("head-again"))
same <- identical(value("base"), value("head"))
functions_changed <- names(which(value("base", FALSE) != value("head", FALSE)))
errors <- sum(grepl("^error:", rows$sha256))
fasta_step <- system2("git", c("diff", "--quiet", base, "HEAD", "--", "R/func_prot_annotation_fasta_processing_helpers.R"))
verdict <- c(paste("base:", base), paste("head:", system2("git", c("rev-parse", "HEAD"), stdout = TRUE)), R.version.string,
    paste("head values deterministic:", deterministic), paste("gate-off values identical to base:", same),
    paste("items hashed as an error text (both sides):", errors),
    paste("deparsed configuration functions changed (reported, not judged):",
        if (length(functions_changed)) paste(functions_changed, collapse = ", ") else "none"),
    paste("FASTA step file unchanged (aa_seq_tbl_final):", identical(fasta_step, 0L)),
    if (!same) paste("differing items:", paste(names(which(value("base") != value("head"))), collapse = ", ")))
writeLines(verdict, file.path(out_dir, "gate-off-verdict.txt"))
writeLines(verdict)
system2("git", c("worktree", "remove", "--force", tree))
if (!(deterministic && same && errors == 0L && identical(fasta_step, 0L))) quit(status = 1L)
