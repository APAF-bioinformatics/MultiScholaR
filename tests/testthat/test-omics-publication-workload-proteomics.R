.publicationProteomicsRepoPath <- function(...) {
    file.path(
        normalizePath(testthat::test_path("..", ".."), mustWork = TRUE),
        ...
    )
}

for (.publication_proteomics_source in c(
    "omics_workload_contract.R",
    "omics_publication_protocol.R",
    "omics_publication_proteomics_contracts.R",
    "omics_publication_proteomics_model.R",
    "omics_publication_proteomics_serializers.R",
    "omics_publication_proteomics_truth.R",
    "omics_publication_proteomics_sources.R",
    "omics_publication_proteomics_governance.R",
    "omics_publication_proteomics_negative.R",
    "omics_publication_proteomics_cleanup.R",
    "omics_publication_proteomics_pilot.R",
    "omics_publication_workload_proteomics.R"
)) {
    sys.source(
        .publicationProteomicsRepoPath(
            "tools",
            "profiling",
            .publication_proteomics_source
        ),
        envir = environment()
    )
}
rm(.publication_proteomics_source)

.publicationProteomicsParameters <- function() {
    publicationReadJson(.publicationProteomicsRepoPath(
        "tests",
        "testdata",
        "omics-performance",
        "proteomics",
        "parameters-v1.json"
    ))
}

.publicationProteomicsBinding <- function(path) {
    list(
        path = path,
        sha256 = publicationFileDigest(path)
    )
}

.publicationProteomicsDimensions <- function(format) {
    switch(
        format,
        diann = list(
            protein_count = 4L,
            peptide_count = 12L,
            precursor_count = 16L,
            sample_count = 6L,
            quantity_count = 96L,
            input_row_count = 96L,
            input_column_count = 15L,
            assay_count = 1L
        ),
        maxquant = list(
            protein_count = 12L,
            peptide_count = 0L,
            precursor_count = 0L,
            sample_count = 6L,
            quantity_count = 72L,
            input_row_count = 12L,
            input_column_count = 10L,
            assay_count = 1L
        ),
        fragpipe = list(
            protein_count = 12L,
            peptide_count = 0L,
            precursor_count = 0L,
            sample_count = 6L,
            quantity_count = 72L,
            input_row_count = 12L,
            input_column_count = 8L,
            assay_count = 1L
        ),
        pd_tmt = list(
            protein_count = 4L,
            peptide_count = 16L,
            precursor_count = 0L,
            sample_count = 6L,
            quantity_count = 96L,
            input_row_count = 16L,
            input_column_count = 8L,
            assay_count = 1L
        )
    )
}

.publicationProteomicsFormatSchema <- function(format) {
    orientation <- c(
        diann = "long",
        maxquant = "wide",
        fragpipe = "wide",
        pd_tmt = "wide_peptide_input"
    )[[format]]
    list(
        schema_id = paste("proteomics", format, "fixture", sep = "."),
        schema_version = "1.0.0",
        orientation = unname(orientation),
        delimiter = "tab",
        required_columns = list("synthetic_required_column"),
        quantity_column_rule = "format_specific_fixture",
        line_ending = "LF",
        encoding = "UTF-8"
    )
}

.publicationProteomicsContract <- function(format, seed = 210001L) {
    capability <- c(
        list(omic_type = "proteomics"),
        proteomicsPublicationCapabilities()[[format]]
    )
    parameter_path <- paste0(
        "tests/testdata/omics-performance/proteomics/parameters-v1.json"
    )
    source_paths <- c(
        "tools/profiling/omics_publication_proteomics_contracts.R",
        "tools/profiling/omics_publication_proteomics_model.R"
    )
    streams <- as.list(stats::setNames(
        seed + seq_len(8L) * 100L,
        c(
            "hierarchy", "multiplicity", "peptide_offsets", "batch",
            "residual", "mar", "mnar", "annotations"
        )
    ))
    list(
        schema = "multischolar.omics_publication_proteomics_workload",
        schema_version = "1.0.0",
        contract_id = paste("test", format, seed, sep = "."),
        owner_ticket_id = "OMICS-ART-065",
        status = "frozen",
        workload_id = paste("proteomics", format, "representative", sep = "."),
        workload_class = "representative",
        capability = capability,
        format_schema = .publicationProteomicsFormatSchema(format),
        dimensions = .publicationProteomicsDimensions(format),
        model_profile_id = "proteomics.publication.declared.v1",
        parameter_authority = .publicationProteomicsBinding(parameter_path),
        source_authority = .publicationProteomicsBinding(parameter_path),
        split_authority = .publicationProteomicsBinding(parameter_path),
        generator = list(
            mode = "generated",
            source_bindings = lapply(
                as.list(source_paths),
                .publicationProteomicsBinding
            ),
            chunk_rows = 100L,
            output_filename = paste0(format, ".tsv"),
            truth_filename = paste0(format, "-truth.json"),
            fixture_payload = NULL,
            fixture_truth = NULL
        ),
        rng = list(
            kind = "L'Ecuyer-CMRG",
            normal_kind = "Inversion",
            sample_kind = "Rejection",
            seed_family_id = "generated.pilot.200000-299999.v1",
            seed = seed,
            streams = streams
        ),
        truth_contract = list(
            schema_id = "multischolar.omics_publication_proteomics_truth",
            oracle_type = "generator_derived_latent_truth",
            independent_from_importer = TRUE,
            source_bindings = list(.publicationProteomicsBinding(
                "tools/profiling/omics_publication_proteomics_truth.R"
            ))
        ),
        execution = list(
            preparation_processes = 2L,
            generation_inside_measured_worker = FALSE,
            generated_payload_committed = FALSE,
            timeout_seconds = 1800L,
            maximum_temporary_bytes = 53687091200
        ),
        privacy = list(
            classification = "public_generated",
            private_source = FALSE,
            private_values_retained = FALSE,
            source_project_identity_reused = FALSE
        ),
        claim_scope = list(
            evidence_class = "public_generated",
            scientific_authority = FALSE,
            performance_authority = TRUE,
            cross_project_authority = FALSE,
            promotion_authority = FALSE,
            limitations = list("Declared synthetic workload only")
        ),
        expected_digests = list(
            payload_sha256 = strrep("a", 64L),
            truth_sha256 = strrep("b", 64L)
        ),
        publication_authority = FALSE
    )
}

test_that("proteomics parameter and exclusion authorities are exact", {
    parameters <- .publicationProteomicsParameters()
    exclusions <- publicationReadJson(.publicationProteomicsRepoPath(
        "tests",
        "testdata",
        "omics-performance",
        "proteomics",
        "exclusions-v1.json"
    ))
    expect_silent(proteomicsPublicationValidateParameters(parameters))
    expect_silent(proteomicsPublicationValidateExclusions(exclusions))
    expect_length(parameters$parameters, 22L)

    changed <- publicationGovernanceCopy(parameters)
    changed$parameters[[11L]]$value[[1L]] <- 0.09
    expect_error(
        proteomicsPublicationValidateParameters(changed),
        class = "multischolar_publication_proteomics_contract_error"
    )

    changed <- publicationGovernanceCopy(parameters)
    changed$parameters[[1L]]$allowed_claim_vocabulary <- list("realistic")
    expect_error(
        proteomicsPublicationValidateParameters(changed),
        class = "multischolar_publication_proteomics_contract_error"
    )

    changed_exclusions <- publicationGovernanceCopy(exclusions)
    changed_exclusions$excluded_capabilities[[1L]]$serializer_dispatch <- TRUE
    expect_error(
        proteomicsPublicationValidateExclusions(changed_exclusions),
        class = "multischolar_publication_proteomics_contract_error"
    )
})

test_that("all four generated format contracts validate exact dimensions", {
    for (format in names(proteomicsPublicationCapabilities())) {
        contract <- .publicationProteomicsContract(format)
        expect_silent(proteomicsPublicationValidateWorkload(contract))
    }

    changed <- .publicationProteomicsContract("diann")
    changed$dimensions$quantity_count <- 95L
    expect_error(
        proteomicsPublicationValidateWorkload(changed),
        class = "multischolar_publication_proteomics_contract_error"
    )

    changed <- .publicationProteomicsContract("maxquant")
    changed$capability$capability_id <-
        "proteomics.fragpipe.protein.lfq.v1"
    expect_error(
        proteomicsPublicationValidateWorkload(changed),
        class = "multischolar_publication_proteomics_contract_error"
    )
})

test_that("hierarchical model is deterministic and mechanism labelled", {
    parameters <- .publicationProteomicsParameters()
    for (format in names(proteomicsPublicationCapabilities())) {
        contract <- .publicationProteomicsContract(format)
        first <- proteomicsPublicationModelPlan(contract, parameters)
        second <- proteomicsPublicationModelPlan(contract, parameters)
        entity_count <- if (identical(format, "diann")) {
            contract$dimensions$precursor_count
        } else if (identical(format, "pd_tmt")) {
            contract$dimensions$peptide_count
        } else {
            contract$dimensions$protein_count
        }
        first_block <- proteomicsPublicationGenerateBlock(
            first,
            seq_len(entity_count)
        )
        second_block <- proteomicsPublicationGenerateBlock(
            second,
            seq_len(entity_count)
        )

        expect_identical(first$proteins, second$proteins, info = format)
        expect_identical(first$peptides, second$peptides, info = format)
        expect_identical(first$precursors, second$precursors, info = format)
        expect_identical(first_block, second_block, info = format)
        expect_equal(
            dim(first_block$values),
            c(entity_count, contract$dimensions$sample_count),
            info = format
        )
        expect_false(any(
            first_block$mar_missing & first_block$mnar_missing
        ), info = format)
        expect_true(all(
            is.na(first_block$values) ==
                (first_block$mar_missing | first_block$mnar_missing)
        ), info = format)
        expect_gt(first$correlation$minimum_eigenvalue, 0)
        expect_true(all(first_block$detection_probability >= 0), info = format)
        expect_true(all(first_block$detection_probability <= 1), info = format)
        expect_true(all(first_block$residual_sigma > 0), info = format)
        expect_setequal(
            unique(first$proteins$effect_class),
            c("up", "down", "unaffected")
        )
    }
})

test_that("seed and hierarchy changes invalidate model identity", {
    parameters <- .publicationProteomicsParameters()
    first <- .publicationProteomicsContract("diann", seed = 210001L)
    second <- .publicationProteomicsContract("diann", seed = 210002L)
    first_plan <- proteomicsPublicationModelPlan(first, parameters)
    second_plan <- proteomicsPublicationModelPlan(second, parameters)
    expect_false(identical(first_plan$proteins, second_plan$proteins))
    expect_false(identical(first_plan$peptides, second_plan$peptides))

    changed <- publicationGovernanceCopy(first)
    changed$dimensions$peptide_count <- 3L
    changed$dimensions$precursor_count <- 16L
    expect_error(
        proteomicsPublicationModelPlan(changed, parameters),
        class = "multischolar_publication_proteomics_contract_error"
    )
})

test_that("format serializers reproduce bytes truth and importer oracles", {
    for (format in names(proteomicsPublicationCapabilities())) {
        contract <- .publicationProteomicsContract(format)
        first_root <- tempfile(paste0("proteomics-publication-", format, "-a-"))
        second_root <- tempfile(paste0("proteomics-publication-", format, "-b-"))
        first <- proteomicsPublicationPrepareGenerated(
            contract,
            first_root,
            verify_expected = FALSE
        )
        second <- proteomicsPublicationPrepareGenerated(
            contract,
            second_root,
            verify_expected = FALSE
        )

        expect_identical(first$payload$sha256, second$payload$sha256)
        expect_identical(first$truth$sha256, second$truth$sha256)
        expect_identical(first$receipt$sha256, second$receipt$sha256)
        result <- proteomicsPublicationRunImported(
            contract,
            first$payload$path,
            first$truth$path
        )
        expect_true(result$workflow_evidence$truth_valid, info = format)
        expect_true(
            result$workflow_evidence$differential_direction$valid,
            info = format
        )
        expect_false(result$workflow_evidence$promotion_authority)

        frozen <- contract
        frozen$expected_digests <- list(
            payload_sha256 = first$payload$sha256,
            truth_sha256 = first$truth$sha256
        )
        third <- proteomicsPublicationPrepareGenerated(
            frozen,
            tempfile(paste0("proteomics-publication-", format, "-c-")),
            verify_expected = TRUE
        )
        expect_identical(third$payload$sha256, first$payload$sha256)
        expect_identical(third$truth$sha256, first$truth$sha256)
    }
})

test_that("private DIA and FASTA inspection retains only approved aggregates", {
    root <- withr::local_tempdir(pattern = "proteomics-private-")
    salt_path <- file.path(root, "private", "salt")
    tsv_path <- file.path(root, "private-report.tsv")
    fasta_path <- file.path(root, "reference.crdownload")
    writeLines(c(
        paste(c(
            "Protein.Group", "Stripped.Sequence", "Run",
            "Precursor.Quantity", "Precursor.Normalised"
        ), collapse = "\t"),
        "P1\tPEPTIDEA\tS1\t100\t90",
        "P1\tPEPTIDEB\tS2\t110\t100",
        "P2\tPEPTIDEC\tS1\t120\t110"
    ), tsv_path, useBytes = TRUE)
    writeLines(c(
        ">sp|P1|Synthetic one",
        "MPEPTIDEK",
        ">sp|P2|Synthetic two",
        "ACDEFGHIKLMNPQRSTVWY"
    ), fasta_path, useBytes = TRUE)
    proteomicsPublicationCreateSalt(salt_path)
    expect_identical(
        as.integer(file.info(salt_path)$mode) %% 512L,
        strtoi("600", base = 8L)
    )
    envelope <- proteomicsPublicationPrivateEnvelope(
        tsv_path,
        fasta_path,
        salt_path
    )
    expect_silent(proteomicsPublicationValidatePrivateEnvelope(envelope))
    expect_identical(envelope$report$row_count, 3L)
    expect_identical(envelope$report$column_count, 5L)
    expect_identical(envelope$report$unique_run_count, 2L)
    expect_identical(envelope$fasta$record_count, 2L)
    expect_true(envelope$fasta$syntax_valid)
    expect_identical(
        envelope$fasta$completeness_status,
        "unconfirmed_download"
    )
    encoded <- jsonlite::toJSON(envelope, auto_unbox = TRUE, null = "null")
    expect_false(grepl(root, encoded, fixed = TRUE))
    expect_false(grepl("PEPTIDEA", encoded, fixed = TRUE))
    expect_false(grepl("sp|P1", encoded, fixed = TRUE))

    leaked <- publicationGovernanceCopy(envelope)
    leaked$report$path <- tsv_path
    expect_error(
        proteomicsPublicationValidatePrivateEnvelope(leaked),
        class = "multischolar_publication_proteomics_privacy_error"
    )

    expect_error(
        proteomicsPublicationReadSalt(file.path(root, "missing-salt")),
        class = "multischolar_publication_proteomics_privacy_error"
    )
})

test_that("private FASTA syntax cannot establish download completeness", {
    root <- withr::local_tempdir(pattern = "proteomics-fasta-")
    salt_path <- file.path(root, "salt")
    fasta_path <- file.path(root, "bad.crdownload")
    proteomicsPublicationCreateSalt(salt_path)
    writeLines(c(">record", "MPEPTIDE1"), fasta_path)
    record <- proteomicsPublicationInspectPrivateFasta(
        fasta_path,
        proteomicsPublicationReadSalt(salt_path)
    )
    expect_false(record$syntax_valid)
    expect_identical(record$invalid_sequence_line_count, 1L)
    expect_identical(record$completeness_status, "unconfirmed_download")
})

test_that("source split and governance successors fail closed", {
    root <- .publicationProteomicsRepoPath(
        "tests",
        "testdata",
        "omics-performance",
        "proteomics"
    )
    acquisition <- publicationReadJson(file.path(
        root,
        "public-acquisition-v1.json"
    ))
    sources <- publicationReadJson(file.path(root, "sources-v1.json"))
    splits <- publicationReadJson(file.path(root, "splits-v1.json"))
    governance <- publicationReadJson(file.path(
        root,
        "governance-successor-v1.json"
    ))
    expect_silent(proteomicsPublicationValidateAcquisition(acquisition))
    expect_silent(proteomicsPublicationValidateSources(sources))
    expect_silent(proteomicsPublicationValidateSplits(splits))
    expect_silent(proteomicsPublicationValidateGovernanceSuccessor(governance))
    expect_identical(length(acquisition$candidates), 17L)
    expect_true(all(vapply(
        sources$capability_decisions,
        \(decision) !isTRUE(decision$promotion_eligible),
        logical(1)
    )))

    synthetic_substitution <- publicationGovernanceCopy(acquisition)
    synthetic_substitution$candidates[[1L]]$counts_toward_cross_project <- TRUE
    expect_error(
        proteomicsPublicationValidateAcquisition(synthetic_substitution),
        class = "multischolar_publication_proteomics_contract_error"
    )

    promoted <- publicationGovernanceCopy(sources)
    promoted$capability_decisions[[1L]]$promotion_eligible <- TRUE
    expect_error(
        proteomicsPublicationValidateSources(promoted),
        class = "multischolar_publication_proteomics_contract_error"
    )

    leaked <- publicationGovernanceCopy(splits)
    leaked$assignments[[1L]]$source_id <-
        leaked$pilot_calibration_assignments[[1L]]$source_id
    expect_error(
        proteomicsPublicationValidateSplits(leaked),
        class = "multischolar_publication_proteomics_contract_error"
    )

    predecessor_drift <- publicationGovernanceCopy(governance)
    predecessor_drift$projects_predecessor$sha256 <- strrep("0", 64L)
    expect_error(
        proteomicsPublicationValidateGovernanceSuccessor(predecessor_drift),
        class = "multischolar_publication_proteomics_contract_error"
    )
})

test_that("hand-reviewed fixtures replay independently of generation", {
    paths <- sort(list.files(
        .publicationProteomicsRepoPath(
            "tests",
            "testdata",
            "omics-performance",
            "proteomics",
            "workloads"
        ),
        pattern = "^fixture-.*[.]json$",
        full.names = TRUE
    ))
    expect_length(paths, 4L)
    for (path in paths) {
        contract <- publicationReadJson(path)
        expect_silent(proteomicsPublicationValidateWorkload(contract))
        first <- proteomicsPublicationPrepareFixture(
            contract,
            tempfile("proteomics-fixture-a-")
        )
        second <- proteomicsPublicationPrepareFixture(
            contract,
            tempfile("proteomics-fixture-b-")
        )
        expect_identical(first$payload$sha256, second$payload$sha256)
        expect_identical(first$truth$sha256, second$truth$sha256)
        expect_identical(first$receipt$sha256, second$receipt$sha256)
        result <- proteomicsPublicationRunImported(
            contract,
            first$payload$path,
            first$truth$path
        )
        expect_true(result$workflow_evidence$truth_valid)
        expect_identical(
            result$workflow_evidence$differential_direction$oracle,
            "hand_reviewed_fixture_e2e"
        )
        expect_false(result$workflow_evidence$promotion_authority)
    }
})

test_that("format-specific malformed and edge contracts match current behavior", {
    authority <- publicationReadJson(.publicationProteomicsRepoPath(
        "tests",
        "testdata",
        "omics-performance",
        "proteomics",
        "negative",
        "contracts-v1.json"
    ))
    expect_silent(proteomicsPublicationValidateNegativeAuthority(authority))
    expect_length(authority$cases, 16L)
    for (case in authority$cases) {
        prepared <- proteomicsPublicationPrepareNegative(
            case,
            tempfile(paste0("proteomics-negative-", case$format, "-"))
        )
        evidence <- proteomicsPublicationEvaluateNegative(case, prepared$path)
        expect_identical(evidence$observed_outcome, case$expected_outcome)
        expect_false(evidence$performance_authority)
        expect_false(evidence$promotion_authority)
    }

    changed <- publicationGovernanceCopy(authority)
    changed$cases[[1L]]$performance_authority <- TRUE
    expect_error(
        proteomicsPublicationValidateNegativeAuthority(changed),
        class = "multischolar_publication_proteomics_contract_error"
    )
})

test_that("generated payload cleanup is owned archived and symlink safe", {
    sandbox <- withr::local_tempdir(pattern = "proteomics-cleanup-")
    withr::local_dir(sandbox)
    root <- file.path(".omics-publication-workloads", "proteomics")
    generated <- file.path(root, "generated", "attempt-1")
    protected <- file.path(root, "private")
    dir.create(file.path(generated, "logs"), recursive = TRUE)
    dir.create(file.path(generated, "payload"), recursive = TRUE)
    dir.create(protected, recursive = TRUE)
    writeLines("log", file.path(generated, "logs", "run.log"))
    writeLines("payload", file.path(generated, "payload", "large.tsv"))
    writeLines("secret", file.path(protected, "salt"))
    expect_true(file.symlink(
        normalizePath(file.path(protected, "salt")),
        file.path(generated, "payload", "protected-link")
    ))
    plan <- list(
        schema = "multischolar.omics_publication_proteomics_cleanup_plan",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-065",
        status = "approved_for_execution",
        root = root,
        archive_root = "cleanup-evidence/test",
        protected_paths = list("private"),
        removals = list(list(
            path = "generated/attempt-1",
            reason = "test generated payload"
        )),
        publication_authority = FALSE
    )
    dry_run <- proteomicsPublicationRunCleanup(plan, execute = FALSE)
    expect_identical(dry_run$status, "dry_run")
    expect_true(dir.exists(generated))
    expect_identical(dry_run$removals[[1L]]$before$symlink_count, 1L)

    result <- proteomicsPublicationRunCleanup(plan, execute = TRUE)
    expect_identical(result$status, "passed")
    expect_false(dir.exists(generated))
    expect_true(file.exists(file.path(protected, "salt")))
    expect_true(file.exists(file.path(
        root,
        "cleanup-evidence/test/generated__attempt-1/retained/logs/run.log"
    )))
})

test_that("heavy qualification uses measured historical thresholds only", {
    measurement <- list(
        status = "passed",
        publication_certifiable = TRUE,
        timed_out = FALSE,
        safety_aborted = FALSE,
        phase_evidence = list(valid = TRUE),
        safety_evidence = publicationTestSafetyEvidence(),
        metrics = list(
            peak_charged_memory_bytes = 4 * 1024^3,
            elapsed_seconds = 30
        )
    )
    qualified <- proteomicsPublicationPilotQualification(measurement)
    expect_true(qualified$qualified)
    expect_identical(
        proteomicsPublicationPilotStatus(
            list(status = "failed", safety_aborted = TRUE, timed_out = FALSE),
            list(qualified = FALSE)
        ),
        "safety_aborted_no_dimension_decision"
    )

    measurement$metrics$peak_charged_memory_bytes <- 4 * 1024^3 - 1
    expect_false(proteomicsPublicationPilotQualification(measurement)$qualified)
    measurement$metrics$elapsed_seconds <- 60
    expect_true(proteomicsPublicationPilotQualification(measurement)$qualified)
    measurement$safety_aborted <- TRUE
    expect_false(proteomicsPublicationPilotQualification(measurement)$qualified)

    environment <- proteomicsPublicationPilotThreadEnvironment(
        "/tmp/home",
        "/tmp/temp",
        "/tmp/package",
        "/tmp/dependency",
        "/tmp/site"
    )
    expect_identical(environment[["OMP_NUM_THREADS"]], "1")
    expect_identical(environment[["R_LIBS_SITE"]], "/tmp/site")
    expect_match(environment[["R_LIBS"]], "/tmp/package", fixed = TRUE)
})

test_that("metabolomics and lipidomics adapters reject proteomics identity", {
    contract <- omicsWorkloadReadContract(.publicationProteomicsRepoPath(
        "tests",
        "testdata",
        "omics-parity",
        "proteomics",
        "workloads",
        "maxquant-protein-lfq-public-v1.json"
    ))
    for (name in c("metabolomics", "lipidomics")) {
        path <- .publicationProteomicsRepoPath(
            "tools",
            "profiling",
            paste0("omics_workload_", name, ".R")
        )
        environment <- new.env(parent = globalenv())
        sys.source(path, envir = environment)
        adapter <- environment$omicsWorkloadAdapter()
        changed <- contract
        changed$adapter <- list(
            adapter_id = adapter$adapter_id,
            adapter_version = adapter$adapter_version,
            source_sha256 = omicsWorkloadFileDigest(path)
        )
        expect_error(
            omicsWorkloadLoadAdapter(path, changed),
            class = "omics_workload_binding_error",
            info = name
        )
    }
})
