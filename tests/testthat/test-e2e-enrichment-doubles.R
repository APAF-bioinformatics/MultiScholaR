library(testthat)

# Tests for helper-e2e-enrichment-doubles.R
# Verifies that fake_gost() and fake_enricher() return stable, correct structures.

# --------------------------------------------------------------------------- #
# fake_gost()
# --------------------------------------------------------------------------- #

test_that("fake_gost returns list with $result and $meta", {
    out <- fake_gost()

    expect_type(out, "list")
    expect_true("result" %in% names(out))
    expect_true("meta" %in% names(out))
    expect_s3_class(out$result, "data.frame")
    expect_type(out$meta, "list")
})

test_that("fake_gost $result has required gprofiler2 columns", {
    required_cols <- c(
        "source", "term_name", "term_id", "term_size"
        , "query_size", "intersection_size", "p_value"
        , "precision", "recall", "effective_domain_size"
        , "source_order", "parents", "significant"
        , "query", "intersection"
    )
    result <- fake_gost()$result

    expect_true(
        all(required_cols %in% names(result))
        , info = paste("Missing columns:", paste(setdiff(required_cols, names(result)), collapse = ", "))
    )
})

test_that("fake_gost returns 5 terms with no sources filter", {
    result <- fake_gost()$result
    expect_equal(nrow(result), 5L)
})

test_that("fake_gost sources filter subsets correctly", {
    go_only <- fake_gost(sources = c("GO:BP"))$result
    expect_true(all(go_only$source == "GO:BP"))
    expect_true(nrow(go_only) > 0L)

    kegg_only <- fake_gost(sources = c("KEGG"))$result
    expect_equal(nrow(kegg_only), 1L)
    expect_equal(kegg_only$term_id, "KEGG:04110")
})

test_that("fake_gost contains recognizable term IDs", {
    result <- fake_gost()$result
    expect_true("GO:0006915" %in% result$term_id)
    expect_true("GO:0007049" %in% result$term_id)
    expect_true("KEGG:04110" %in% result$term_id)
    expect_true("REAC:R-HSA-109582" %in% result$term_id)
})

test_that("fake_gost $result is fully deterministic (query arg does not affect result rows)", {
    r1 <- fake_gost(query = c("P1", "P2"), sources = c("GO:BP", "KEGG"))
    r2 <- fake_gost(query = c("P9", "P10"), sources = c("GO:BP", "KEGG"))
    # $result must be identical; $meta intentionally reflects the query arg
    expect_identical(r1$result, r2$result)
})

test_that("fake_gost $meta contains query_metadata with user_threshold", {
    out <- fake_gost(user_threshold = 0.01)
    expect_equal(out$meta$query_metadata$user_threshold, 0.01)
})

test_that("fake_gost parents column is a list column", {
    result <- fake_gost()$result
    expect_true(is.list(result$parents))
})

# --------------------------------------------------------------------------- #
# fake_enricher()
# --------------------------------------------------------------------------- #

test_that("fake_enricher returns S4 object of class e2eEnrichResult", {
    out <- fake_enricher()
    expect_s4_class(out, "e2eEnrichResult")
})

test_that("fake_enricher @result has required clusterProfiler columns", {
    required_cols <- c(
        "ID", "Description", "GeneRatio", "BgRatio"
        , "pvalue", "p.adjust", "qvalue", "geneID", "Count"
    )
    result <- fake_enricher()@result

    expect_s3_class(result, "data.frame")
    expect_true(
        all(required_cols %in% names(result))
        , info = paste("Missing columns:", paste(setdiff(required_cols, names(result)), collapse = ", "))
    )
})

test_that("fake_enricher returns 3 terms", {
    result <- fake_enricher()@result
    expect_equal(nrow(result), 3L)
})

test_that("fake_enricher result contains GO term IDs", {
    result <- fake_enricher()@result
    expect_true("GO:0006915" %in% result$ID)
    expect_true("GO:0007049" %in% result$ID)
})

test_that("fake_enricher is fully deterministic", {
    r1 <- fake_enricher(gene = c("P1", "P2"))
    r2 <- fake_enricher(gene = c("P9", "P10"))
    expect_identical(r1@result, r2@result)
})

test_that("fake_enricher p-values are realistic (0 < pvalue < 0.05)", {
    result <- fake_enricher()@result
    expect_true(all(result$pvalue > 0))
    expect_true(all(result$pvalue < 0.05))
    expect_true(all(result$p.adjust >= result$pvalue))
})

# --------------------------------------------------------------------------- #
# install_enrichment_doubles()
# --------------------------------------------------------------------------- #

test_that("install_enrichment_doubles runs without error", {
    # Should run silently even when packages are not installed (skips gracefully)
    expect_no_error(install_enrichment_doubles())
})

test_that("install_enrichment_doubles installs gprofiler2 double when package present", {
    skip_if_not(
        requireNamespace("gprofiler2", quietly = TRUE)
        , "gprofiler2 not installed"
    )

    install_enrichment_doubles()
    result <- gprofiler2::gost(query = c("P1", "P2"), organism = "hsapiens")

    expect_type(result, "list")
    expect_true("result" %in% names(result))
    expect_equal(nrow(result$result), 5L)
})

test_that("install_enrichment_doubles installs fake_gost that returns expected structure", {
    skip_if_not(
        requireNamespace("gprofiler2", quietly = TRUE)
        , "gprofiler2 not installed"
    )

    # install_enrichment_doubles scopes to the withr local block below
    local({
        install_enrichment_doubles(.local_envir = parent.frame())
        result <- gprofiler2::gost(query = c("P1", "P2"), organism = "hsapiens")
        expect_equal(nrow(result$result), 5L)
        expect_true("GO:0006915" %in% result$result$term_id)
    })
})
