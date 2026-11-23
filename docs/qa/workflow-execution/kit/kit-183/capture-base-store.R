# KIT-183: capture the base store (tests/testdata/integration/kit183-base-store/) with the package as it was before
# KIT-183 (run on a commit whose R/ equals the branch point). Re-captured with this script when a later ticket changes
# the integration store or the enrichment record (ticket 8.5).
#
# Usage, from the package root with MultiScholaR installed: Rscript capture-base-store.R <out dir>
# Writes <out dir>/store/ (an integration store holding the corpus, the stored reference model and one KEGG enrichment
# of Factor1 on metabolomics_GCMS made by the pre-ranked test, every setting at its default), base-terms.tsv and
# base-summary.txt (that enrichment's terms and summary), gsea-deparse.txt and gsea-fixture-value.rds
# (integrationEnrichmentGsea's formals and body, and its value on the determinism fixture), provenance.tsv.
# The recipe is the one test-integration-enrichment-tests.R builds (intEtProject, intEtModel, intEtCorpusSets).

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 1L)
# Prefixed names: the evaluation environment's parent is the package namespace, whose imports hold functions such as
# utils::data, which a plain `data` would find first.
k183_out <- arguments[[1L]]
if (dir.exists(k183_out)) stop("the output folder exists: ", k183_out)
dir.create(k183_out, recursive = TRUE)
k183_out <- normalizePath(k183_out)
k183_data <- normalizePath(file.path("tests", "testdata", "integration"), mustWork = TRUE)

local(envir = new.env(parent = asNamespace("MultiScholaR")), {
    out <- k183_out
    data <- k183_data
    stored <- function(name) utils::read.delim(file.path(data, "mofa_reference", paste0(name, ".tsv")),
        check.names = FALSE, stringsAsFactors = FALSE)
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
        list(factors = f, weights = weights, r2_per_factor = r2, r2_total = colSums(r2), elbo = c(-2, -1), iterations = 186L,
            options = options, versions = list(MOFA2 = "stored reference"), messages = character(),
            model_file = list(bytes = as.raw(1:3), sha256 = strrep("d", 64L), size = 3))
    }
    kegg_ids <- function(n) {
        ids <- sprintf("C%05d", 10000L + seq_len(n))
        pairs <- seq(2L, 20L, by = 2L)
        ids[pairs] <- ids[pairs - 1L]
        ids
    }
    # The project: the corpus, metabolomics_GCMS given a kegg_id column.
    root <- tempfile("k183-capture-")
    integration_dir <- file.path(root, "integration")
    for (omic in c("lipidomics", "metabolomics", "proteomics", "transcriptomics")) {
        dir <- file.path(integration_dir, omic)
        dir.create(dir, recursive = TRUE)
        stopifnot(all(file.copy(list.files(file.path(data, "corpus", omic), full.names = TRUE), dir)))
    }
    dir <- file.path(integration_dir, "metabolomics")
    # tests/testthat/helper-integration-input.R's readIntegrationInputWithMatrices (no corrected matrices here).
    input <- readIntegrationInput(dir)
    for (name in names(input$views)) input$views[[name]]$matrix <- openIntegrationInputView(dir, name, input)$read()
    read <- input$manifest
    manifest <- list(omic_type = read$omic_type, workflow = read$workflow, matrix_role = read$matrix_role, scale = read$scale,
        correction = NULL, source_state = list(name = "kit183-kegg"), sample_id_column = read$sample_id_column,
        group_column = read$group_column, technical_replicate_column = read$technical_replicate_column,
        biological_sample_column = read$biological_sample_column)
    features <- lapply(input$views, `[[`, "features")
    features$metabolomics_GCMS$kegg_id <- kegg_ids(nrow(features$metabolomics_GCMS))
    writeIntegrationInput(dir, manifest, input$samples, lapply(input$views, `[[`, "matrix"), features)

    manager <- newIntegrationStateOwner(file.path(out, "store"))
    scan <- integrationScanProject(integration_dir)
    integrationCommitImport(manager, integrationImportRequest(scan, scan$views$view))
    run <- runIntegrationModel(manager, list(), trainer = trainer)
    value <- run$value

    # The sets: complete features on single-feature identifiers (intEtCorpusSets).
    view <- "metabolomics_GCMS"
    object <- integrationImportedObject(manager)
    table <- object@views[[view]]
    m <- as.matrix(table[-1L])
    rownames(m) <- as.character(table$feature_id)
    design <- object@design_matrix
    design <- design[design$omic == integrationViewOmics(object)[[view]], , drop = FALSE]
    key <- design$biological_sample[match(colnames(m), design$run)]
    samples <- as.character(value$factors$biological_sample)
    samples <- samples[samples %in% key]
    x <- vapply(samples, function(s) {
        means <- rowMeans(m[, key == s, drop = FALSE], na.rm = TRUE)
        means[is.nan(means)] <- NA_real_
        means
    }, numeric(nrow(m)))
    rownames(x) <- rownames(m)
    annotation <- object@features[[view]]
    ids <- stats::setNames(annotation$kegg_id, annotation$feature_id)
    complete <- rownames(x)[stats::complete.cases(x)]
    pool <- intersect(names(ids)[!duplicated(ids) & !duplicated(ids, fromLast = TRUE)], complete)
    sets <- list(first = unname(ids[pool[1:10]]), second = unname(ids[pool[11:20]]), third = unname(ids[pool[21:30]]))
    pathway <- paste0("path:hsa", sprintf("%05d", seq_along(sets)))
    mapping <- list(database = "KEGG", release = "KEGG Release fixture", sha256 = strrep("e", 64L), url = "fixture",
        downloaded_at = "2026-11-23T00:00:00+0000",
        table = data.frame(pathway = rep(pathway, lengths(sets)), compound = paste0("cpd:", unlist(sets, use.names = FALSE)),
            name = rep(names(sets), lengths(sets)), stringsAsFactors = FALSE))
    enrichment <- runIntegrationEnrichment(manager, list(factor = "Factor1", view = view, database = "kegg",
        id_column = "kegg_id"), mapping = mapping)
    utils::write.table(as.data.frame(enrichment$run$value$terms), file.path(out, "base-terms.tsv"), sep = "\t",
        quote = FALSE, row.names = FALSE)
    writeLines(integrationEnrichmentSummary(enrichment$run$value), file.path(out, "base-summary.txt"))
    closeWorkflowStateOwner(manager)
    unlink(root, recursive = TRUE)

    # integrationEnrichmentGsea as at base: formals and body, and its value on the determinism fixture.
    gsea <- utils::removeSource(integrationEnrichmentGsea)
    writeLines(c(deparse(formals(gsea)), deparse(body(gsea))), file.path(out, "gsea-deparse.txt"))
    gene_ids <- sprintf("G%02d", 1:60)
    ranked <- data.frame(identifier = gene_ids, value = seq(3, -3, length.out = 60), stringsAsFactors = FALSE)
    term2gene <- data.frame(term = rep(c("top", "bottom", "spread", "middle"), c(10, 10, 10, 15)),
        gene = c(gene_ids[1:10], gene_ids[51:60], gene_ids[seq(1, 60, by = 6)], gene_ids[21:35]), stringsAsFactors = FALSE)
    term2name <- data.frame(term = c("top", "bottom", "spread", "middle"), name = c("Top set", "Bottom set",
        "Spread set", "Middle set"), stringsAsFactors = FALSE)
    saveRDS(integrationEnrichmentGsea(ranked, term2gene, term2name, seed = 11L), file.path(out, "gsea-fixture-value.rds"))

    commit <- tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE), error = function(e) NA_character_)
    files <- list.files(out, recursive = TRUE, all.files = TRUE)
    provenance <- data.frame(item = c("commit", "R", "MultiScholaR", "clusterProfiler", "limma", "captured_at", files),
        value = c(commit, R.version.string, as.character(utils::packageVersion("MultiScholaR")),
            as.character(utils::packageVersion("clusterProfiler")), as.character(utils::packageVersion("limma")),
            format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
            vapply(file.path(out, files), function(path) digest::digest(file = path, algo = "sha256"), character(1))),
        stringsAsFactors = FALSE)
    utils::write.table(provenance, file.path(out, "provenance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("captured %d files in %s\n", length(files), out))
})
