# KIT-063 item A: without an uploaded sample match, the model's parameters, its scientific identity and its aligned
# views are those of the base commit. The base tree and the head are each loaded in their own R process
# (pkgload::load_all, as the runner's per-file runs load the package); each imports the tier-1 corpus into a new
# integration store and writes the digests of: the resolved parameters; the identity's parameter and selector
# components (the input and design components name the store's entity ids, written for information); the aligned views
# and sample table; and the product of a model run with a stand-in trainer that returns the stored double-precision
# reference (tests/testdata/integration/mofa_reference), so that the whole model step is compared without MOFA2. On the
# head the same is written again with the downloaded match table uploaded unchanged.
#
# Usage (from the package root; full git history): Rscript match-identity.R <base commit> <out dir>
# Writes <out dir>/match-identity.tsv (side, item, digest) and match-identity-verdict.txt; exits 1 when a compared
# digest differs.

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 2L)
base <- arguments[[1L]]
out_dir <- arguments[[2L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir)
head <- normalizePath(".")
tree <- file.path(tempdir(), "match-identity-base")
stopifnot(system2("git", c("worktree", "add", "--detach", tree, base)) == 0L)

child <- '
arguments <- commandArgs(trailingOnly = TRUE)
suppressMessages(pkgload::load_all(arguments[[1L]], quiet = TRUE, export_all = TRUE, helpers = FALSE))
side <- arguments[[2L]]
data <- arguments[[3L]]
corpus <- file.path(data, "corpus")
reference <- file.path(data, "mofa_reference")
root <- tempfile("match-identity-")
integration <- file.path(root, "integration")
for (omic in c("lipidomics", "metabolomics", "proteomics", "transcriptomics")) {
    dir.create(file.path(integration, omic), recursive = TRUE)
    stopifnot(all(file.copy(list.files(file.path(corpus, omic), full.names = TRUE), file.path(integration, omic))))
}
stored <- function(name) utils::read.delim(file.path(reference, paste0(name, ".tsv")), check.names = FALSE, stringsAsFactors = FALSE)
trainer <- function(views, options) {
    factors <- stored("factors")
    f <- as.matrix(factors[-1L])
    rownames(f) <- factors$biological_sample
    f <- f[colnames(views[[1L]]), , drop = FALSE]
    w <- stored("weights")
    weights <- lapply(stats::setNames(names(views), names(views)), function(view) {
        rows <- w[w$view == view, , drop = FALSE]
        m <- as.matrix(rows[-(1:2)])
        rownames(m) <- paste(view, rows$feature_id, sep = ":")
        m[rownames(views[[view]]), , drop = FALSE]
    })
    v <- stored("variance_explained")
    r2 <- matrix(v$r2_percent, nrow = length(unique(v$factor)), dimnames = list(unique(v$factor), unique(v$view)))
    r2 <- r2[, names(views), drop = FALSE]
    list(factors = f, weights = weights, r2_per_factor = r2, r2_total = colSums(r2), elbo = c(-2, -1), iterations = 2L,
        options = options, versions = list(MOFA2 = "stored reference"), messages = character(),
        model_file = list(bytes = as.raw(1:3), sha256 = strrep("c", 64L), size = 3))
}
rows <- list()
record <- function(item, value) rows[[length(rows) + 1L]] <<- data.frame(side = side, item = item,
    digest = if (is.character(value) && length(value) == 1L && grepl("^[0-9a-f]{64}$", value)) value else
        digest::digest(value, algo = "sha256"), stringsAsFactors = FALSE)
measure <- function(label, parameters) {
    manager <- newIntegrationStateOwner(file.path(root, paste0("store-", label)))
    on.exit(closeWorkflowStateOwner(manager), add = TRUE)
    scan <- integrationScanProject(integration)
    integrationCommitImport(manager, integrationImportRequest(scan, scan$views$view))
    object <- integrationImportedObject(manager)
    if (is.function(parameters)) parameters <- parameters(object)
    resolved <- integrationResolveMofaParameters(object, parameters)
    record(paste0(label, ": resolved parameters"), resolved)
    declaration <- workflowOperationEntry(INTEGRATION_MODEL_OPERATION)$declaration
    identity <- withWorkflowStateScope(manager, function() {
        source <- getWorkflowStateReference(manager, INTEGRATION_INPUTS_STATE)
        workflowComponentScientificIdentity(manager, source, declaration$name, resolved[declaration$parameters])
    })
    record(paste0(label, ": identity, parameters"), identity$parameters)
    record(paste0(label, ": identity, selectors"), identity$selectors)
    record(paste0(label, ": identity, software (information)"), identity$software)
    record(paste0(label, ": identity, input entity (information)"), identity$input)
    aligned <- integrationAlignViews(object, resolved)
    record(paste0(label, ": aligned views"), aligned$views)
    record(paste0(label, ": aligned samples"), aligned$samples)
    run <- runIntegrationModel(manager, parameters, trainer = trainer)
    value <- run$value
    for (name in c("factors", "weights", "variance_explained", "variance_total", "samples", "view_summary", "notes"))
        record(paste0(label, ": product ", name), value[[name]])
    record(paste0(label, ": export parameters"), manager$.context$snapshot$exports$integration_mofa_model$parameters)
}
measure("no table", list())
if (identical(side, "head")) measure("table uploaded unchanged", function(object) list(sample_match_table = integrationSampleMatchTemplate(object)))
utils::write.table(do.call(rbind, rows), arguments[[4L]], sep = "\\t", quote = FALSE, row.names = FALSE)
'
script <- file.path(tempdir(), "match-identity-child.R")
writeLines(child, script)
data <- file.path(head, "tests", "testdata", "integration")

results <- list()
for (side in c("base", "head")) {
    package <- if (side == "base") tree else head
    target <- file.path(tempdir(), sprintf("%s.tsv", side))
    status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", script, package, side, data, target),
        stdout = file.path(out_dir, sprintf("%s.log", side)), stderr = file.path(out_dir, sprintf("%s.log", side)))
    if (status != 0L || !file.exists(target)) stop(sprintf("the %s side failed (status %d); see %s.log", side, status, side))
    results[[side]] <- utils::read.delim(target, stringsAsFactors = FALSE)
}
rows <- do.call(rbind, results)
utils::write.table(rows, file.path(out_dir, "match-identity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
digest_of <- function(side, item) rows$digest[rows$side == side & rows$item == item]
compared <- unique(rows$item[rows$side == "base" & !grepl("(information)", rows$item, fixed = TRUE)])
same_as_base <- vapply(compared, function(item) identical(digest_of("base", item), digest_of("head", item)), logical(1))
unchanged_items <- sub("^no table: ", "table uploaded unchanged: ", compared)
same_unchanged <- vapply(seq_along(compared), function(i) identical(digest_of("head", compared[[i]]),
    digest_of("head", unchanged_items[[i]])), logical(1))
verdict <- c(paste("base:", base), paste("head:", system2("git", c("rev-parse", "HEAD"), stdout = TRUE)), R.version.string,
    sprintf("%s: head %s base", compared, ifelse(same_as_base, "identical to", "DIFFERS from")),
    sprintf("%s: the unchanged uploaded table %s", unchanged_items, ifelse(same_unchanged, "identical", "DIFFERS")),
    paste("every compared digest identical:", all(same_as_base) && all(same_unchanged)))
writeLines(verdict, file.path(out_dir, "match-identity-verdict.txt"))
writeLines(verdict)
system2("git", c("worktree", "remove", "--force", tree))
if (!(all(same_as_base) && all(same_unchanged))) quit(status = 1L)
