.aux077RepoPath <- function(...) {
    file.path(
        normalizePath(testthat::test_path("..", ".."), mustWork = TRUE),
        ...
    )
}

for (.aux077_source in c(
    "omics_publication_protocol.R",
    "omics_publication_auxiliary_contracts.R",
    "omics_publication_auxiliary_model.R",
    "omics_publication_auxiliary_responses.R",
    "omics_publication_auxiliary_sources.R",
    "omics_publication_auxiliary_governance.R",
    "omics_publication_auxiliary_truth.R",
    "omics_publication_workload_auxiliary.R"
)) {
    sys.source(
        .aux077RepoPath("tools", "profiling", .aux077_source),
        envir = environment()
    )
}
rm(.aux077_source)

test_that("phosphosite optimized primitives preserve edge contracts", {
    expect_null(getMaxProb("PEPTIDE", 1L))
    expect_identical(getBestPosition("PEPTIDE", 1L), numeric())
    expect_identical(
        getMaxProb("AA(0.9)BB(0.9)CC(0.2)", 1L),
        c(0.9)
    )
    expect_identical(
        getBestPosition("AA(0.9)BB(0.9)CC(0.2)", 1L),
        c(2L, 4L)
    )
    expect_error(
        getBestPosition("AApBB(0.9)", 1L),
        "should not have little 'p'"
    )
    expect_identical(getPosString(c(10L, 20L), c(2L, 4L)), "(11;13)|(21;23)")
})

test_that("phosphosite truth separates accession and paralog cardinality", {
    workload_id <- "auxiliary.phosphosite_stages.fixture_correctness.v1"
    contract <- auxPublicationReadContract(.aux077RepoPath(
        "tests",
        "testdata",
        "omics-performance",
        "auxiliary",
        "workloads",
        paste0(workload_id, ".json")
    ))
    truth <- publicationReadJson(.aux077RepoPath(
        "tests",
        "testdata",
        "omics-performance",
        "auxiliary",
        "truth",
        paste0(workload_id, ".json")
    ))
    truth_path <- .aux077RepoPath(
        "tests",
        "testdata",
        "omics-performance",
        "auxiliary",
        "truth",
        paste0(workload_id, ".json")
    )
    expect_identical(
        publicationFileDigest(truth_path),
        contract$expected_digests$truth_sha256
    )
    expect_identical(
        truth$contract_basis_sha256,
        auxPublicationContractBasis(contract)
    )
    recalculated <- auxPublicationBuildTruth(
        contract,
        auxPublicationBuildPayload(contract),
        contract$expected_digests$payload_sha256
    )
    expect_equal(
        recalculated$facts$expected_paralog_rows,
        truth$facts$expected_paralog_rows
    )
    expect_lte(truth$facts$expected_paralog_rows, truth$facts$expected_long_rows)
})

test_that("pathway lookups retain first-match and assay fallback semantics", {
    names_mapping <- data.frame(
        metabolite = c("LC first", "LC duplicate", "GC first", "Fallback"),
        chebi_id = c("CHEBI:1", "CHEBI:1", "CHEBI:1", "CHEBI:2"),
        assay = c("LC-MS", "LC-MS", "GC-MS", "LC-MS"),
        stringsAsFactors = FALSE
    )
    kegg_mapping <- data.frame(
        kegg_id = c("cpd:C1", "cpd:C1", "cpd:C2"),
        chebi_id = c("CHEBI:1", "CHEBI:2", "CHEBI:9"),
        stringsAsFactors = FALSE
    )
    lookups <- multiomicsPathwayNameLookups(names_mapping, kegg_mapping)

    expect_identical(
        multiomicsMapPathwayIds("cpd:C1/CHEBI:2/raw", "LC-MS", lookups),
        "LC first, Fallback, raw"
    )
    expect_identical(
        multiomicsMapPathwayIds("CHEBI:1", "GC-MS", lookups),
        "GC first"
    )
    expect_identical(
        multiomicsMapPathwayIds("CHEBI:2", "GC-MS", lookups),
        "Fallback"
    )
    expect_identical(
        multiomicsMapPathwayIds("cpd:C2", "LC-MS", lookups),
        "cpd:C2"
    )
    expect_identical(
        multiomicsMapPathwayIds("", "LC-MS", lookups),
        NA_character_
    )
})

test_that("auxiliary APIs remain outside workflow and backend ownership", {
    functions <- c(
        "processMultisiteEvidence",
        "plotMofaWeights",
        "runMetabolomicsEnrichmentAnalysis",
        "runMetabolomicsPathwayEnrichment",
        "runOneStringDbRankEnrichmentMofa"
    )
    forbidden <- c(
        "workflowData", "workflow_data", "WorkflowState",
        "ArtifactWorkflowState", "storage_policy", "requested_backend"
    )
    for (name in functions) {
        fun <- getExportedValue("MultiScholaR", name)
        expect_false(any(names(formals(fun)) %in% forbidden), info = name)
    }
})
