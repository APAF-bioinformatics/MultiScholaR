makeSavedConsumerDiaOutputs <- function() {
    set.seed(10010L)
    samples <- sprintf("S%02d", seq_len(8L))
    features <- sprintf("P%03d", seq_len(200L))
    design <- data.frame(
        group = factor(rep(c("control", "treated"), each = 4L),
            levels = c("control", "treated")),
        row.names = samples
    )
    quantities <- matrix(
        stats::rnorm(length(features) * length(samples)),
        nrow = length(features),
        dimnames = list(features, samples)
    )
    changed <- seq.int(10L, length(features), by = 10L)
    quantities[changed, design$group == "treated"] <-
        quantities[changed, design$group == "treated"] + 2
    result <- suppressMessages(runTestsContrasts(
        quantities,
        "grouptreated-groupcontrol",
        design,
        "~0+group",
        weights = NA,
        treat_lfc_cutoff = NA
    ))
    list(samples = samples, features = features, design = design,
        quantities = quantities, result = result)
}

makeSavedConsumerProteinProduct <- function(outputs) {
    table <- data.frame(
        Protein.Ids = outputs$features,
        outputs$quantities,
        check.names = FALSE
    )
    rownames(table) <- NULL
    design <- data.frame(
        Run = outputs$samples,
        group = outputs$design$group,
        row.names = outputs$samples
    )
    object <- methods::new(
        "ProteinQuantitativeData",
        protein_quant_table = table,
        design_matrix = design,
        sample_id = "Run",
        group_id = "group",
        args = list(normaliseBetweenSamples = list(normalisation_method = "none"))
    )
    capture.output(normalized <- normaliseBetweenSamples(
        object,
        normalisation_method = "none"
    ))
    normalized
}

renderSavedConsumerVolcano <- function(table) {
    path <- withr::local_tempfile(fileext = ".png")
    grDevices::png(path, width = 400L, height = 300L)
    on.exit(grDevices::dev.off(), add = TRUE)
    graphics::plot(
        table$logFC,
        -log10(table$raw_pvalue),
        xlab = "log2 fold change",
        ylab = "-log10 p-value"
    )
    grDevices::dev.off()
    on.exit(NULL, add = FALSE)
    readBin(path, "raw", n = file.info(path)$size)
}

consumerScientificIdentity <- function(outputs, selector) {
    newSavedWorkflowScientificIdentity(
        digest::digest(outputs$quantities, algo = "sha256"),
        digest::digest(outputs$design, algo = "sha256"),
        digest::digest(list(formula = "~0+group", weights = NA), algo = "sha256"),
        selector,
        paste0("MultiScholaR-", as.character(utils::packageVersion("MultiScholaR")))
    )
}

test_that("real DIA products fits and tables use explicit schemas and exact keys", {
    producer_calls <- 0L
    produce <- function() {
        producer_calls <<- producer_calls + 1L
        makeSavedConsumerDiaOutputs()
    }
    outputs <- produce()
    product_value <- makeSavedConsumerProteinProduct(outputs)
    product <- newSavedWorkflowScientificProduct(
        product_value,
        "between_sample_normalization",
        provenance = list(route = "iq", method = "none")
    )
    fit <- newSavedWorkflowDaFit(
        outputs$result$fit.eb,
        rownames(outputs$result$fit.eb$coefficients),
        rownames(outputs$result$fit.eb$design),
        colnames(outputs$result$fit.eb$design)
    )
    result_table <- newSavedWorkflowResultTable(
        outputs$result$results[[1L]],
        ".rownames",
        "grouptreated-groupcontrol",
        "fdr_qvalue"
    )
    store <- newWorkflowEntityStore(
        withr::local_tempdir(),
        "consumer-project",
        "DIA",
        max_cache_bytes = 67108864
    )
    product_identity <- consumerScientificIdentity(outputs, "normalized-product")
    fit_identity <- consumerScientificIdentity(outputs, "base-fit")
    table_identity <- consumerScientificIdentity(outputs, "treated-control-table")
    product_ref <- saveWorkflowConsumerEntity(
        store,
        product,
        "normalized-product",
        product_identity,
        "scientific_product",
        factors = list(route = "iq", method = "none")
    )
    fit_ref <- saveWorkflowConsumerEntity(
        store,
        fit,
        "da-fit",
        fit_identity,
        "da_fit",
        factors = list(route = "iq", formula = "~0+group")
    )
    table_ref <- saveWorkflowConsumerEntity(
        store,
        result_table,
        "da-table",
        table_identity,
        "result_table",
        factors = list(route = "iq", contrast = "treated-control")
    )
    expect_identical(
        c(product_ref$kind, fit_ref$kind, table_ref$kind),
        c("scientific_product", "da_fit", "result_table")
    )

    product_hit <- resolveSavedWorkflowEntity(
        store,
        "scientific_product",
        product_identity,
        factors = list(route = "iq", method = "none")
    )
    fit_hit <- resolveSavedWorkflowEntity(
        store,
        "da_fit",
        fit_identity,
        factors = list(route = "iq", formula = "~0+group")
    )
    table_hit <- resolveSavedWorkflowEntity(
        store,
        "result_table",
        table_identity,
        factors = list(route = "iq", contrast = "treated-control")
    )
    expect_identical(product_hit$value, product)
    expect_identical(fit_hit$value, fit)
    expect_identical(table_hit$value, result_table)
    expect_identical(producer_calls, 1L)
    expect_false(product_hit$capabilities$column_projection)
    expect_false(table_hit$capabilities$predicate_pushdown)

    stale <- table_identity
    stale$software <- paste0(stale$software, "-changed")
    miss <- resolveSavedWorkflowEntity(
        store,
        "result_table",
        stale,
        factors = list(route = "iq", contrast = "treated-control")
    )
    expect_identical(miss$status, "miss")
    expect_null(miss$value)
    expect_identical(producer_calls, 1L)
    expect_error(
        hydrateSavedWorkflowProjection(store, table_ref, c("logFC", "fdr_qvalue")),
        class = "saved_entity_unsupported"
    )
})

test_that("render identities invalidate style only and asset hits do not hydrate science", {
    outputs <- makeSavedConsumerDiaOutputs()
    result_table <- newSavedWorkflowResultTable(
        outputs$result$results[[1L]],
        ".rownames",
        "grouptreated-groupcontrol",
        "fdr_qvalue"
    )
    store <- newWorkflowEntityStore(
        withr::local_tempdir(),
        "asset-project",
        "DIA",
        max_cache_bytes = 67108864
    )
    table_identity <- consumerScientificIdentity(outputs, "volcano-table")
    table_ref <- saveWorkflowConsumerEntity(
        store,
        result_table,
        "volcano-table",
        table_identity,
        "result_table"
    )
    renderer_calls <- 0L
    render <- function() {
        renderer_calls <<- renderer_calls + 1L
        renderSavedConsumerVolcano(result_table$table)
    }
    asset <- newSavedWorkflowRenderedAsset(
        render(),
        "image/png",
        c(400L, 300L)
    )
    render_identity <- newSavedWorkflowRenderIdentity(
        table_ref,
        aesthetics = list(color = "navy", alpha = 0.7),
        format = "png",
        dimensions = c(400L, 300L),
        renderer = "base-graphics-v1"
    )
    asset_ref <- saveWorkflowConsumerEntity(
        store,
        asset,
        "volcano-asset",
        render_identity,
        "rendered_asset"
    )
    pending_input <- file.path(store$directory, ".pending", "publication-input")
    writeBin(as.raw(c(9L, 8L, 7L)), pending_input)
    releaseWorkflowEntityCache(store)
    expect_true(file.exists(pending_input))
    before <- workflowEntityCacheInfo(store)
    hit <- resolveSavedWorkflowAsset(store, table_ref, render_identity)
    after <- workflowEntityCacheInfo(store)
    expect_identical(hit$status, "hit")
    expect_identical(hit$value, asset)
    expect_identical(renderer_calls, 1L)
    expect_identical(after$payload_reads, before$payload_reads + 1L)
    expect_identical(after$payload_decodes, before$payload_decodes + 1L)
    expect_identical(after$entries, 1L)
    second <- resolveSavedWorkflowAsset(store, table_ref, render_identity)
    expect_identical(second$value, asset)
    expect_identical(renderer_calls, 1L)
    expect_identical(workflowEntityCacheInfo(store)$cache_hits, after$cache_hits + 1L)

    changed_style <- newSavedWorkflowRenderIdentity(
        table_ref,
        aesthetics = list(color = "firebrick", alpha = 0.7),
        format = "png",
        dimensions = c(400L, 300L),
        renderer = "base-graphics-v1"
    )
    expect_false(identical(changed_style, render_identity))
    expect_identical(
        resolveSavedWorkflowAsset(store, table_ref, changed_style)$status,
        "miss"
    )
    expect_identical(renderer_calls, 1L)
    expect_identical(
        resolveSavedWorkflowEntity(store, "result_table", table_identity)$status,
        "hit"
    )

    path <- file.path(store$directory, asset_ref$payload$relative_path)
    original <- readBin(path, "raw", n = file.info(path)$size)
    changed <- original
    changed[[length(changed)]] <- as.raw(bitwXor(as.integer(changed[[length(changed)]]), 1L))
    writeBin(changed, path)
    expect_error(
        resolveSavedWorkflowAsset(store, table_ref, render_identity),
        class = "saved_entity_corrupt"
    )
    writeBin(original, path)
    releaseWorkflowEntityCache(store)
    expect_true(file.exists(pending_input))
    expect_identical(
        resolveSavedWorkflowAsset(store, table_ref, render_identity)$value,
        asset
    )

    foreign <- newWorkflowEntityStore(
        withr::local_tempdir(),
        "foreign-project",
        "DIA"
    )
    expect_error(
        resolveSavedWorkflowAsset(foreign, table_ref, render_identity),
        class = "saved_entity_scope"
    )
})

test_that("unavailable panels remain unavailable and malformed envelopes fail", {
    outputs <- makeSavedConsumerDiaOutputs()
    table <- newSavedWorkflowResultTable(
        outputs$result$results[[1L]],
        ".rownames",
        "grouptreated-groupcontrol",
        "fdr_qvalue"
    )
    store <- newWorkflowEntityStore(withr::local_tempdir(), "unavailable-project", "DIA")
    table_identity <- consumerScientificIdentity(outputs, "unavailable-table")
    table_ref <- saveWorkflowConsumerEntity(
        store,
        table,
        "unavailable-table",
        table_identity,
        "result_table"
    )
    identity <- newSavedWorkflowRenderIdentity(
        table_ref,
        aesthetics = list(panel = "limpa_uncertainty"),
        format = "png",
        dimensions = c(400L, 300L),
        renderer = "base-graphics-v1"
    )
    unavailable <- newSavedWorkflowRenderedAsset(
        raw(),
        NA_character_,
        c(NA_integer_, NA_integer_),
        status = "unavailable",
        reason = "uncertainty panel is unavailable for the IQ route"
    )
    saveWorkflowConsumerEntity(
        store,
        unavailable,
        "unavailable-panel",
        identity,
        "rendered_asset"
    )
    hit <- resolveSavedWorkflowAsset(store, table_ref, identity)
    expect_identical(hit$value$status, "unavailable")
    expect_length(hit$value$bytes, 0L)
    expect_match(hit$value$reason, "IQ route")
    expect_error(
        newSavedWorkflowRenderedAsset(raw(), "image/png", c(400L, 300L)),
        class = "saved_entity_invalid"
    )
    expect_error(
        newSavedWorkflowScientificProduct(NULL, "bad", status = "available"),
        class = "saved_entity_invalid"
    )
    expect_error(
        saveWorkflowConsumerEntity(store, list(), "bad", table_identity, "da_fit"),
        class = "saved_entity_invalid"
    )
    expect_error(
        newSavedWorkflowResultTable(data.frame(x = 1), "missing", "contrast", "q"),
        class = "saved_entity_invalid"
    )
})
