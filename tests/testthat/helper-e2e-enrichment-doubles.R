# E2E test doubles for gprofiler2::gost() and clusterProfiler::enricher().
# Provides deterministic doubles that avoid network access during browser tests.
# Use install_enrichment_doubles() inside an E2E test setup to activate.

# --------------------------------------------------------------------------- #
# Mock S4 class
# --------------------------------------------------------------------------- #

if (!methods::isClass("e2eEnrichResult")) {
    methods::setClass(
        "e2eEnrichResult"
        , slots = c(result = "data.frame")
    )
}

# --------------------------------------------------------------------------- #
# Static enrichment data
# --------------------------------------------------------------------------- #

.e2e_gost_result_df <- data.frame(
    source = c("GO:BP", "GO:BP", "GO:CC", "KEGG", "REAC")
    , term_name = c(
        "apoptotic process"
        , "cell cycle"
        , "mitochondrion"
        , "Cell cycle"
        , "Hemostasis"
    )
    , term_id = c(
        "GO:0006915"
        , "GO:0007049"
        , "GO:0005739"
        , "KEGG:04110"
        , "REAC:R-HSA-109582"
    )
    , term_size = c(1392L, 1780L, 892L, 131L, 418L)
    , query_size = c(10L, 10L, 10L, 10L, 10L)
    , intersection_size = c(4L, 5L, 3L, 3L, 2L)
    , p_value = c(0.0012, 0.0034, 0.0085, 0.0127, 0.0243)
    , precision = c(0.400, 0.500, 0.300, 0.300, 0.200)
    , recall = c(0.00287, 0.00281, 0.00336, 0.02290, 0.00478)
    , effective_domain_size = c(20717L, 20717L, 20717L, 20717L, 20717L)
    , source_order = c(1L, 2L, 3L, 4L, 5L)
    , parents = I(list(
        c("GO:0008219")
        , c("GO:0022402")
        , c("GO:0043231")
        , character(0)
        , character(0)
    ))
    , significant = c(TRUE, TRUE, TRUE, TRUE, TRUE)
    , query = c("query_1", "query_1", "query_1", "query_1", "query_1")
    , intersection = c(
        "P04637/P10600/Q07817/P55957"
        , "P04637/P10600/Q07817/P55957/O14757"
        , "P04637/P10600/Q07817"
        , "P04637/P10600/Q07817"
        , "P04637/P10600"
    )
    , stringsAsFactors = FALSE
)

.e2e_enricher_result_df <- data.frame(
    ID = c("GO:0006915", "GO:0007049", "GO:0005739")
    , Description = c("apoptotic process", "cell cycle", "mitochondrion")
    , GeneRatio = c("4/10", "5/10", "3/10")
    , BgRatio = c("1392/20717", "1780/20717", "892/20717")
    , pvalue = c(0.0012, 0.0034, 0.0085)
    , p.adjust = c(0.0060, 0.0085, 0.0142)
    , qvalue = c(0.0058, 0.0083, 0.0139)
    , geneID = c(
        "P04637/P10600/Q07817/P55957"
        , "P04637/P10600/Q07817/P55957/O14757"
        , "P04637/P10600/Q07817"
    )
    , Count = c(4L, 5L, 3L)
    , stringsAsFactors = FALSE
)

# --------------------------------------------------------------------------- #
# fake_gost()
# --------------------------------------------------------------------------- #

#' Deterministic double for gprofiler2::gost()
#'
#' @param query Gene/protein IDs (accepted, ignored).
#' @param organism Organism identifier (accepted, ignored).
#' @param ordered_query Logical (accepted, ignored).
#' @param multi_query Logical (accepted, ignored).
#' @param significant Filter to significant results (accepted, ignored).
#' @param exclude_iea Exclude IEA evidence (accepted, ignored).
#' @param measure_underrepresentation Logical (accepted, ignored).
#' @param evcodes Include evidence codes (accepted, ignored).
#' @param user_threshold P-value threshold (accepted, ignored).
#' @param correction_method Multiple testing method (accepted, ignored).
#' @param domain_scope Annotation scope (accepted, ignored).
#' @param custom_bg Custom background (accepted, ignored).
#' @param numeric_ns Numeric namespace (accepted, ignored).
#' @param sources Ontology sources to filter (used to subset returned terms).
#' @param as_short_link Return short link (accepted, ignored).
#' @param highlight Highlight driver terms (accepted, ignored).
#' @return List with $result data.frame and $meta matching gprofiler2 output.
fake_gost <- function(
    query = character()
    , organism = "hsapiens"
    , ordered_query = FALSE
    , multi_query = FALSE
    , significant = TRUE
    , exclude_iea = FALSE
    , measure_underrepresentation = FALSE
    , evcodes = FALSE
    , user_threshold = 0.05
    , correction_method = "g_SCS"
    , domain_scope = "annotated"
    , custom_bg = NULL
    , numeric_ns = ""
    , sources = NULL
    , as_short_link = FALSE
    , highlight = FALSE
) {
    result <- .e2e_gost_result_df
    if (!is.null(sources) && length(sources) > 0) {
        result <- result[result$source %in% sources, ]
    }
    list(
        result = result
        , meta = list(
            query_metadata = list(
                queries = list(query_1 = query)
                , user_threshold = user_threshold
                , organism = organism
                , sources = if (is.null(sources)) {
                    c("GO:BP", "GO:CC", "GO:MF", "KEGG", "REAC")
                } else {
                    sources
                }
            )
        )
    )
}

# --------------------------------------------------------------------------- #
# fake_enricher()
# --------------------------------------------------------------------------- #

#' Deterministic double for clusterProfiler::enricher()
#'
#' @param gene Gene/protein IDs (accepted, ignored).
#' @param universe Background gene set (accepted, ignored).
#' @param TERM2GENE Term-to-gene mapping (accepted, ignored).
#' @param TERM2NAME Term-to-name mapping (accepted, ignored).
#' @param pvalueCutoff P-value cutoff (accepted, ignored).
#' @param pAdjustMethod Adjustment method (accepted, ignored).
#' @param ... Additional arguments (accepted, ignored).
#' @return S4 object of class e2eEnrichResult with @result slot.
fake_enricher <- function(
    gene = character()
    , universe = NULL
    , TERM2GENE = NULL
    , TERM2NAME = NULL
    , pvalueCutoff = 0.05
    , pAdjustMethod = "BH"
    , ...
) {
    methods::new("e2eEnrichResult", result = .e2e_enricher_result_df)
}

# --------------------------------------------------------------------------- #
# install_enrichment_doubles()
# --------------------------------------------------------------------------- #

#' Install deterministic doubles for enrichment functions
#'
#' Replaces gprofiler2::gost() and clusterProfiler::enricher() with
#' deterministic doubles that avoid network access in E2E browser tests.
#' Automatically restored when `.local_envir` exits (withr::defer).
#'
#' If a package is not installed, that double is silently skipped.
#'
#' @param .local_envir Environment to scope teardown into (default: caller frame).
install_enrichment_doubles <- function(.local_envir = parent.frame()) {
    if (requireNamespace("gprofiler2", quietly = TRUE)) {
        local_mocked_bindings(
            gost = fake_gost
            , .package = "gprofiler2"
            , .env = .local_envir
        )
    }

    if (requireNamespace("clusterProfiler", quietly = TRUE)) {
        .install_ns_binding(
            env = asNamespace("clusterProfiler")
            , name = "enricher"
            , value = fake_enricher
            , .local_envir = .local_envir
        )
    }

    invisible(NULL)
}

# Inline localNamespaceBinding pattern (self-contained, no dependency on
# test-general-enrichment-shared.R which is a test file, not a helper).
.install_ns_binding <- function(env, name, value, .local_envir = parent.frame()) {
    had_binding <- exists(name, envir = env, inherits = FALSE)
    old_value <- if (had_binding) get(name, envir = env, inherits = FALSE) else NULL
    was_locked <- had_binding && bindingIsLocked(name, env)

    if (was_locked) unlockBinding(name, env)
    assign(name, value, envir = env)
    if (was_locked) lockBinding(name, env)

    withr::defer({
        if (exists(name, envir = env, inherits = FALSE) &&
                bindingIsLocked(name, env)) {
            unlockBinding(name, env)
        }
        if (had_binding) {
            assign(name, old_value, envir = env)
        } else if (exists(name, envir = env, inherits = FALSE)) {
            rm(list = name, envir = env)
        }
        if (was_locked && exists(name, envir = env, inherits = FALSE)) {
            lockBinding(name, env)
        }
    }, envir = .local_envir)
}
