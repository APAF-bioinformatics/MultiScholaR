# KIT-184: the base stores, made with the package as it was before KIT-184 (run on the branch point's tree, loaded in
# its own process by base-store-job.R; VAL-64 (2): the stores are built and used inside one runner job and never kept).
#
# Usage, with the base MultiScholaR loaded (pkgload::load_all) and the base tree as working directory:
#   source("capture-base-store.R") with commandArgs(trailingOnly = TRUE) = <out dir>
# Writes under <out dir>:
#   plain/lipidomics/          the corpus's lipidomics input rewritten by the base's writer (no correction);
#   corrected/integration/     the corpus as a project whose metabolomics and proteomics inputs hold a corrected matrix
#                              and a correction of k = 2 (the recipe of tests/testthat/test-integration-removed-factors.R,
#                              intRfProject), written by the base's writer: no removed factors;
#   corrected/store/           a store on it, the stand-in model trained on the corrected values of its three corrected
#                              data sets, with a KEGG camera entry and a KEGG gsea entry of Factor1 on metabolomics_GCMS;
#   default/store/             a store on the same inputs, the default model, with a KEGG camera entry;
#   provenance.tsv             each file's sha256.

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 1L)
k184_out <- arguments[[1L]]
if (dir.exists(k184_out)) stop("the output folder exists: ", k184_out)
dir.create(k184_out, recursive = TRUE)
k184_out <- normalizePath(k184_out)
k184_data <- normalizePath(file.path("tests", "testdata", "integration"), mustWork = TRUE)

local(envir = new.env(parent = asNamespace("MultiScholaR")), {
    out <- k184_out
    data <- k184_data
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
    read_whole <- function(dir) {
        input <- readIntegrationInput(dir)
        for (view in names(input$views)) input$views[[view]]$matrix <- openIntegrationInputView(dir, view, input)$read()
        input
    }
    fields <- c("omic_type", "workflow", "matrix_role", "scale", "sample_id_column", "group_column",
        "technical_replicate_column", "biological_sample_column")

    # 1. The lipidomics input without a correction, rewritten by the base's writer.
    plain <- read_whole(file.path(data, "corpus", "lipidomics"))
    manifest <- plain$manifest[fields]
    manifest$source_state <- list(name = "kit184-plain")
    writeIntegrationInput(file.path(out, "plain", "lipidomics"), manifest, plain$samples,
        lapply(plain$views, `[[`, "matrix"), lapply(plain$views, `[[`, "features"))

    # 2. The project with corrected matrices (intRfProject's recipe).
    integration_dir <- file.path(out, "corrected", "integration")
    for (omic in c("lipidomics", "metabolomics", "proteomics", "transcriptomics")) {
        dir <- file.path(integration_dir, omic)
        dir.create(dir, recursive = TRUE)
        stopifnot(all(file.copy(list.files(file.path(data, "corpus", omic), full.names = TRUE), dir)))
    }
    for (omic in c("metabolomics", "proteomics")) {
        dir <- file.path(integration_dir, omic)
        input <- read_whole(dir)
        manifest <- input$manifest[fields]
        manifest$correction <- list(state = "ruv_corrected", method = "kit184-test", k = 2L, replicate_structure = "group",
            grouping_variable = "group")
        manifest$source_state <- list(name = "kit184")
        features <- lapply(input$views, `[[`, "features")
        matrices <- lapply(input$views, `[[`, "matrix")
        if (omic == "metabolomics") for (view in names(features)) features[[view]]$kegg_id <- kegg_ids(nrow(features[[view]]))
        replicate <- as.integer(sub("^.*_", "", as.character(input$samples$tech_rep_group)))
        w <- cbind(ifelse(replicate %% 2L == 0L, 1, -1), replicate - 5.5)
        corrected <- lapply(matrices, function(table) {
            m <- as.matrix(table[-1L])
            i <- seq_len(nrow(m))
            a <- cbind(0.30 + 0.01 * (i %% 13L), 0.02 + 0.05 * ((i %% 5L) - 2))
            m <- m - a %*% t(w) - ifelse(i %% 2L == 0L, 0.5, 0)
            table[-1L] <- as.data.frame(m)
            table
        })
        writeIntegrationInput(dir, manifest, input$samples, matrices, features, corrected = corrected)
    }

    # 3. The two stores, each with its entries.
    sets_of <- function(manager, value, view) {
        object <- integrationImportedObject(manager)
        table <- object@views[[view]]
        m <- as.matrix(table[-1L])
        rownames(m) <- as.character(table$feature_id)
        design <- object@design_matrix
        design <- design[design$omic == integrationViewOmics(object)[[view]], , drop = FALSE]
        key <- design$biological_sample[match(colnames(m), design$run)]
        samples <- as.character(value$factors$biological_sample)
        samples <- samples[samples %in% key]
        x <- vapply(samples, function(s) rowMeans(m[, key == s, drop = FALSE], na.rm = TRUE), numeric(nrow(m)))
        rownames(x) <- rownames(m)
        ids <- stats::setNames(object@features[[view]]$kegg_id, object@features[[view]]$feature_id)
        complete <- rownames(x)[stats::complete.cases(x)]
        pool <- intersect(names(ids)[!duplicated(ids) & !duplicated(ids, fromLast = TRUE)], complete)
        list(first = unname(ids[pool[1:10]]), second = unname(ids[pool[11:20]]), third = unname(ids[pool[21:30]]))
    }
    mapping_of <- function(sets) {
        pathway <- paste0("path:hsa", sprintf("%05d", seq_along(sets)))
        list(database = "KEGG", release = "KEGG Release fixture", sha256 = strrep("e", 64L), url = "fixture",
            downloaded_at = "2026-11-24T00:00:00+0000",
            table = data.frame(pathway = rep(pathway, lengths(sets)), compound = paste0("cpd:", unlist(sets, use.names = FALSE)),
                name = rep(names(sets), lengths(sets)), stringsAsFactors = FALSE))
    }
    view <- "metabolomics_GCMS"
    request <- list(factor = "Factor1", view = view, database = "kegg", id_column = "kegg_id")
    for (part in c("corrected", "default")) {
        manager <- newIntegrationStateOwner(file.path(out, part, "store"))
        scan <- integrationScanProject(integration_dir)
        integrationCommitImport(manager, integrationImportRequest(scan, scan$views$view))
        parameters <- if (part == "corrected") list(input_matrix = "ruv_corrected",
            views = c("proteomics", "metabolomics_LCMS_Pos", "metabolomics_GCMS")) else list()
        run <- runIntegrationModel(manager, parameters, trainer = trainer)
        mapping <- mapping_of(sets_of(manager, run$value, view))
        camera <- runIntegrationEnrichment(manager, c(request, list(test = "camera")), mapping = mapping)
        stopifnot(identical(camera$run$value$design_covariates, "none"))
        if (part == "corrected") runIntegrationEnrichment(manager, c(request, list(test = "gsea")), mapping = mapping)
        closeWorkflowStateOwner(manager)
    }
    files <- list.files(out, recursive = TRUE, all.files = TRUE, no.. = TRUE)
    files <- files[!grepl("(^|/)[.]workflow-entities/", files)]
    utils::write.table(data.frame(file = files, sha256 = vapply(file.path(out, files), function(path) digest::digest(file = path,
        algo = "sha256"), character(1)), stringsAsFactors = FALSE), file.path(out, "provenance.tsv"), sep = "\t", quote = FALSE,
        row.names = FALSE)
})
