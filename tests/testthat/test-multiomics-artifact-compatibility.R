test_that("phosphosite and multiomics APIs are not false artifact workflows", {
    capabilities <- workflowCapabilityCatalogue()
    identities <- lapply(capabilities, `[[`, "identity")
    expect_false(any(vapply(
        identities,
        \(identity) identity$omic_type %in% c(
            "phosphoproteomics",
            "multiomics"
        ),
        logical(1)
    )))
    expect_true(is.function(processMultisiteEvidence))
    expect_true(is.function(runMetabolomicsEnrichmentAnalysis))
    expect_true(is.function(plotMofaWeights))
    decisions <- c(
        artifactProteomicsNonDiaCloseoutDecisions(),
        artifactMetabolomicsCloseoutDecisions(),
        artifactLipidomicsCloseoutDecisions()
    )
    expect_true(all(vapply(
        decisions,
        \(decision) identical(decision$promotion_status, "withheld"),
        logical(1)
    )))
    dia <- artifactDiaWorkflowDescriptor()
    expect_false(dia$certification$auto_eligible)
    expect_identical(dia$certification$status, "dual_write")
})

test_that("multiomics enrichment plotting remains deterministic", {
    skip_if_not_installed("tidytext")
    input <- data.frame(
        comparison = c("B_vs_A", "B_vs_A", "C_vs_A"),
        category = c("Pathway", "Pathway", "Process"),
        termDescription = c(
            "Alpha beta pathway",
            "Gamma pathway",
            "Delta process"
        ),
        enrichmentScore = c(3, 2, 4),
        falseDiscoveryRate = c(0.01, 0.02, 0.005),
        direction = c("top", "bottom", "both ends"),
        genesMapped = c(5, 3, 7),
        stringsAsFactors = FALSE
    )
    first <- printStringDbFunctionalEnrichmentBarGraph(input, word_limit = 2)
    second <- printStringDbFunctionalEnrichmentBarGraph(input, word_limit = 2)
    expect_s3_class(first, "ggplot")
    expect_identical(
        ggplot2::ggplot_build(first)$data,
        ggplot2::ggplot_build(second)$data
    )
    labels <- first$data$termDescriptionAbbrev
    expect_true("Alpha beta ..." %in% labels)
    expect_identical(
        levels(first$data$direction),
        c("top", "bottom", "both ends")
    )
})

test_that("MOFA plotting preserves class ordering and file-product contract", {
    weights <- matrix(
        c(0.5, -0.7, 0.2, -0.1, 0.9, -0.4),
        ncol = 1L,
        dimnames = list(
            c(
                "GENE1_proteome", "GENE2_proteome", "GENE3_proteome",
                "GENE4_proteome", "GENE5_proteome", "GENE6_proteome"
            ),
            "Factor1"
        )
    )
    output_dir <- withr::local_tempdir()
    products <- character()
    environment <- new.env(parent = asNamespace("MultiScholaR"))
    environment$fakeGetWeights <- \(model) model$weights
    environment$omic_type <- "multiomics"
    environment$experiment_label <- "compatibility"
    environment$getProjectPaths <- \(...) list(mofa_plots_dir = output_dir)
    environment$savePlot <- function(
        plot,
        base_path,
        plot_name,
        formats
    ) {
        products <<- file.path(base_path, paste0(plot_name, ".", formats))
        invisible(lapply(
            products,
            \(path) writeLines("deterministic plot", path)
        ))
        invisible(products)
    }
    plot_function <- plotMofaWeights
    body_text <- paste(deparse(body(plot_function)), collapse = "\n")
    body_text <- gsub(
        "MOFA2::get_weights",
        "fakeGetWeights",
        body_text,
        fixed = TRUE
    )
    body(plot_function) <- parse(text = body_text)[[1L]]
    environment(plot_function) <- environment
    plot <- plot_function(
        list(weights = list(proteomics = weights)),
        "proteomics",
        factor_level = Factor1
    )
    expect_s3_class(plot, "ggplot")
    expect_identical(nrow(plot$data), 6L)
    expect_false(any(grepl("_proteome$", as.character(plot$data$gene_symbol))))
    expect_identical(levels(plot$data$Direction), c("Positive", "Negative"))
    expect_setequal(tools::file_ext(products), c("png", "pdf"))
    expect_true(all(file.exists(products)))
})
