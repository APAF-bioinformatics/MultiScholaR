.settlementFixture <- function(
    fail_at = NULL,
    hydration_active = FALSE,
    descriptor = artifactDiaWorkflowDescriptor()
) {
    state <- new.env(parent = emptyenv())
    state$current <- list(generation_id = "generation-old", marker = "old")
    state$encode_calls <- 0L
    state$rollback_calls <- 0L
    state$cleanup_calls <- 0L
    state$parity_calls <- 0L
    state$release_calls <- 0L
    state$close_calls <- 0L
    maybe_fail <- function(stage) {
        if (identical(fail_at, stage)) {
            rlang::abort(
                paste("injected", stage, "failure"),
                class = paste0("settlement_test_", stage)
            )
        }
    }
    roles <- c("data_tbl", "data_cln")
    operations <- list(
        current = function() state$current,
        encode = function(sources) {
            state$encode_calls <- state$encode_calls + 1L
            maybe_fail("encode")
            list(
                intent_id = "intent-new",
                generation_id = "generation-new",
                refs = stats::setNames(lapply(roles, function(role) {
                    list(
                        artifact_id = paste0("artifact-", role),
                        semantic_digest = strrep(substr(role, 1L, 1L), 64L)
                    )
                }), roles),
                bounded_metadata = list(row_count = 4L, column_count = 3L),
                temporary_paths = character(),
                complete_payload_hydrated = FALSE
            )
        },
        bounded_validate = function(prepared) {
            maybe_fail("bounded_validate")
            list(
                valid = TRUE,
                refs_validated = TRUE,
                bounded_metadata_validated = TRUE,
                complete_payload_hydrated = FALSE
            )
        },
        publish = function(prepared, previous) {
            state$current <- list(
                generation_id = prepared$generation_id,
                marker = "new"
            )
            maybe_fail("publish")
            list(
                current = state$current,
                generation_id = prepared$generation_id,
                refs = prepared$refs,
                published_atomically = TRUE
            )
        },
        lean_snapshot = function(publication, bounded_metadata) {
            maybe_fail("lean_snapshot")
            list(
                lineage = list(parent = "generation-old"),
                current = publication$current,
                bounded_metadata = bounded_metadata,
                audit_state = list(status = "settled"),
                config_state = list(mode = "artifact"),
                descriptor_pin = list(
                    descriptor_id = descriptor$descriptor_id,
                    descriptor_digest = descriptor$descriptor_digest
                ),
                rollback_checkpoint = list(generation_id = "generation-old"),
                locators = lapply(publication$refs, function(ref) {
                    list(
                        artifact_id = ref$artifact_id,
                        semantic_digest = ref$semantic_digest
                    )
                })
            )
        },
        release_sources = function(workflow_data, source_fields) {
            state$release_calls <- state$release_calls + 1L
            workflow_data[[source_fields[[1L]]]] <- NULL
            maybe_fail("release_sources")
            workflow_data[[source_fields[[2L]]]] <- NULL
            invisible(TRUE)
        },
        close_resources = function(workflow_data) {
            state$close_calls <- state$close_calls + 1L
            maybe_fail("close_resources")
            list(
                registry_connections = 0L,
                query_handles = 0L,
                results = 0L,
                writers = 0L,
                locks = 0L
            )
        },
        rollback = function(previous, publication, prepared) {
            state$rollback_calls <- state$rollback_calls + 1L
            state$current <- previous
            invisible(TRUE)
        },
        cleanup_prepared = function(prepared) {
            state$cleanup_calls <- state$cleanup_calls + 1L
            invisible(TRUE)
        },
        hydration_active = function(workflow_data) hydration_active
    )
    parity <- list(
        worker_id = "settlement.parity.worker.v1",
        process_isolation = "fresh_R_process",
        required = TRUE,
        input_binding_ids = c(
            "source_input", "artifact_refs", "descriptor", "package_revision"
        ),
        output_digest_required = TRUE,
        candidate_process_allowed = FALSE
    )
    contract <- newArtifactSettlementContract(
        contract_id = "test.artifact.settlement.v1",
        descriptor = descriptor,
        capability_id = descriptor$descriptor_id,
        source_fields = roles,
        stage_roles = roles,
        compatibility_strategy = "read_only_artifact_reconstruction",
        consumer_inventory_digest = function() strrep("a", 64L),
        resolver = function(ref) ref,
        parity_worker_contract = parity,
        operations = operations,
        max_lean_bytes = 1024^2
    )
    workflow_data <- new.env(parent = emptyenv())
    workflow_data$data_tbl <- data.frame(
        feature = c("A", "B"),
        sample = c("S1", "S2"),
        value = c(1, 2),
        stringsAsFactors = FALSE
    )
    workflow_data$data_cln <- workflow_data$data_tbl
    list(
        contract = contract,
        state = state,
        workflow_data = workflow_data
    )
}

test_that("settlement encodes once and reaches payload-free rest", {
    fixture <- .settlementFixture()

    result <- settleArtifactWorkflow(
        fixture$workflow_data,
        fixture$contract
    )
    expect_identical(result$status, "settled")
    expect_identical(
        unlist(result$transitions, use.names = FALSE),
        c(
            "idle", "prepared", "bounded_validated", "published",
            "source_released", "settled"
        )
    )
    expect_identical(fixture$state$encode_calls, 1L)
    expect_identical(fixture$state$rollback_calls, 0L)
    expect_identical(fixture$state$cleanup_calls, 0L)
    expect_null(fixture$workflow_data$data_tbl)
    expect_null(fixture$workflow_data$data_cln)
    expect_true(result$source_fields_released)
    expect_true(result$resources_closed)
    expect_false(result$complete_parity_executed_in_commit_process)
    expect_true(result$parity_worker_contract$required)
    expect_identical(fixture$state$parity_calls, 0L)
    expect_silent(artifactSettlementValidateLeanSnapshot(
        result$snapshot,
        1024^2
    ))
})

test_that("shared settlement accepts every descriptor without science dispatch", {
    descriptors <- artifactDescriptorCatalogueValues(
        artifactWorkflowDescriptorCatalogue()
    )
    for (descriptor in descriptors) {
        fixture <- .settlementFixture(descriptor = descriptor)
        result <- settleArtifactWorkflow(
            fixture$workflow_data,
            fixture$contract
        )
        expect_identical(
            result$descriptor_id,
            descriptor$descriptor_id,
            info = descriptor$descriptor_id
        )
        expect_identical(
            result$capability_id,
            descriptor$descriptor_id
        )
        expect_null(fixture$workflow_data$data_tbl)
        expect_null(fixture$workflow_data$data_cln)
        expect_false(result$complete_parity_executed_in_commit_process)
    }
})

test_that("every settlement failure restores previous current and sources", {
    stages <- c(
        "encode", "bounded_validate", "publish", "lean_snapshot",
        "release_sources", "close_resources"
    )
    for (stage in stages) {
        fixture <- .settlementFixture(fail_at = stage)
        original_tbl <- fixture$workflow_data$data_tbl
        original_cln <- fixture$workflow_data$data_cln

        error <- rlang::catch_cnd(settleArtifactWorkflow(
            fixture$workflow_data,
            fixture$contract
        ))
        expect_s3_class(
            error,
            "multischolar_artifact_settlement_rolled_back"
        )
        expect_identical(
            fixture$state$current,
            list(generation_id = "generation-old", marker = "old"),
            info = stage
        )
        expect_identical(fixture$workflow_data$data_tbl, original_tbl)
        expect_identical(fixture$workflow_data$data_cln, original_cln)
        expect_identical(fixture$state$rollback_calls, 1L)
        expect_identical(fixture$state$cleanup_calls, 1L)
        expect_identical(
            unlist(error$transitions, use.names = FALSE)[[
                length(error$transitions)
            ]],
            "rolled_back"
        )
    }
})

test_that("settlement validators reject tampered transition outputs", {
    mutations <- list(
        complete_payload = function(fixture) {
            operation <- fixture$contract$operations$encode
            fixture$contract$operations$encode <- function(sources) {
                prepared <- operation(sources)
                prepared$complete_payload_hydrated <- TRUE
                prepared
            }
            fixture
        },
        bounded_validation = function(fixture) {
            fixture$contract$operations$bounded_validate <- function(prepared) {
                list(
                    valid = FALSE,
                    refs_validated = TRUE,
                    bounded_metadata_validated = TRUE,
                    complete_payload_hydrated = FALSE
                )
            }
            fixture
        },
        publication_digest = function(fixture) {
            operation <- fixture$contract$operations$publish
            fixture$contract$operations$publish <- function(prepared, previous) {
                publication <- operation(prepared, previous)
                publication$refs$data_tbl$semantic_digest <- strrep("0", 64L)
                publication
            }
            fixture
        },
        incomplete_release = function(fixture) {
            fixture$contract$operations$release_sources <- function(
                workflow_data,
                source_fields
            ) {
                workflow_data[[source_fields[[1L]]]] <- NULL
                invisible(TRUE)
            }
            fixture
        }
    )
    for (name in names(mutations)) {
        fixture <- mutations[[name]](.settlementFixture())
        error <- rlang::catch_cnd(settleArtifactWorkflow(
            fixture$workflow_data,
            fixture$contract
        ))
        expect_s3_class(
            error,
            "multischolar_artifact_settlement_rolled_back"
        )
        expect_identical(
            fixture$state$current,
            list(generation_id = "generation-old", marker = "old"),
            info = name
        )
        expect_false(is.null(fixture$workflow_data$data_tbl))
        expect_false(is.null(fixture$workflow_data$data_cln))
    }
})

test_that("every declared resource must be closed before settlement", {
    resources <- c(
        "registry_connections", "query_handles", "results", "writers", "locks"
    )
    for (resource in resources) {
        fixture <- .settlementFixture()
        operation <- fixture$contract$operations$close_resources
        fixture$contract$operations$close_resources <- function(workflow_data) {
            info <- operation(workflow_data)
            info[[resource]] <- 1L
            info
        }
        expect_error(
            settleArtifactWorkflow(fixture$workflow_data, fixture$contract),
            class = "multischolar_artifact_settlement_rolled_back",
            info = resource
        )
        expect_identical(
            fixture$state$current,
            list(generation_id = "generation-old", marker = "old")
        )
        expect_false(is.null(fixture$workflow_data$data_tbl))
        expect_false(is.null(fixture$workflow_data$data_cln))
    }
})

test_that("mutated settlement authorities fail before encoding", {
    mutations <- list(
        parity = function(contract) {
            contract$parity_worker_contract$required <- FALSE
            contract
        },
        consumer_digest = function(contract) {
            contract$consumer_inventory_digest <- function() "invalid"
            contract
        },
        operation_order = function(contract) {
            contract$operations <- rev(contract$operations)
            contract
        }
    )
    for (name in names(mutations)) {
        fixture <- .settlementFixture()
        fixture$contract <- mutations[[name]](fixture$contract)
        expect_error(
            settleArtifactWorkflow(fixture$workflow_data, fixture$contract),
            class = "multischolar_invalid_artifact_settlement_contract",
            info = name
        )
        expect_identical(fixture$state$encode_calls, 0L)
        expect_identical(
            fixture$state$current,
            list(generation_id = "generation-old", marker = "old")
        )
    }
})

test_that("rollback and cleanup failures remain explicit and data recoverable", {
    fixture <- .settlementFixture(fail_at = "bounded_validate")
    fixture$contract$operations$rollback <- function(...) {
        stop("injected rollback failure")
    }
    fixture$contract$operations$cleanup_prepared <- function(...) {
        stop("injected cleanup failure")
    }
    error <- rlang::catch_cnd(settleArtifactWorkflow(
        fixture$workflow_data,
        fixture$contract
    ))
    expect_s3_class(error, "multischolar_artifact_settlement_rolled_back")
    expect_s3_class(error$rollback_error, "error")
    expect_s3_class(error$cleanup_error, "error")
    expect_identical(
        fixture$state$current,
        list(generation_id = "generation-old", marker = "old")
    )
    expect_false(is.null(fixture$workflow_data$data_tbl))
    expect_false(is.null(fixture$workflow_data$data_cln))
})

test_that("settlement rejects active hydration before encoding", {
    fixture <- .settlementFixture(hydration_active = TRUE)

    expect_error(
        settleArtifactWorkflow(fixture$workflow_data, fixture$contract),
        class = "multischolar_overlapping_artifact_hydration"
    )
    expect_identical(fixture$state$encode_calls, 0L)
    expect_identical(fixture$state$rollback_calls, 0L)
    expect_false(is.null(fixture$workflow_data$data_tbl))
})

test_that("lean snapshots reject complete and process-bound ownership", {
    fixture <- .settlementFixture()
    result <- settleArtifactWorkflow(fixture$workflow_data, fixture$contract)

    data_snapshot <- result$snapshot
    data_snapshot$bounded_metadata$payload <- data.frame(value = 1)
    expect_error(
        artifactSettlementValidateLeanSnapshot(data_snapshot, 1024^2),
        class = "multischolar_nonlean_artifact_settlement"
    )

    environment_snapshot <- result$snapshot
    environment_snapshot$audit_state$owner <- new.env(parent = emptyenv())
    expect_error(
        artifactSettlementValidateLeanSnapshot(environment_snapshot, 1024^2),
        class = "multischolar_nonlean_artifact_settlement"
    )

    oversized <- result$snapshot
    oversized$config_state$annotation <- paste(rep("x", 2000L), collapse = "")
    expect_error(
        artifactSettlementValidateLeanSnapshot(oversized, 100L),
        class = "multischolar_nonlean_artifact_settlement"
    )
})

test_that("settlement code does not depend on garbage collection", {
    source <- readLines(
        testthat::test_path(
            "..",
            "..",
            "R",
            "utils_artifact_settlement_contract.R"
        ),
        warn = FALSE
    )
    expect_false(any(grepl("gc(", source, fixed = TRUE)))
    expect_false(any(grepl("artifactReleaseTransientMemory", source, fixed = TRUE)))
})
