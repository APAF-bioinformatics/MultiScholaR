nondiaCloseout028SkipDependencies <- function() {
    for (package in c("DBI", "duckdb", "filelock", "ps")) {
        testthat::skip_if_not_installed(package)
    }
}

nondiaCloseout028RepoPath <- function(...) {
    root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
    file.path(root, ...)
}

nondiaCloseout028ReadEvidence <- function(decision) {
    jsonlite::read_json(
        nondiaCloseout028RepoPath(
            decision$performance_evidence$memory_baseline$evidence_id
        ),
        simplifyVector = FALSE
    )
}

nondiaCloseout028PreserveGlobals <- function() {
    names <- c(
        "data_tbl", "data_cln", "design_matrix", "contrasts_tbl",
        "uniprot_dat_cln", "filtering_progress", "project_dirs",
        "workflow_context"
    )
    existing <- names[vapply(
        names,
        exists,
        logical(1),
        envir = .GlobalEnv,
        inherits = FALSE
    )]
    values <- mget(existing, envir = .GlobalEnv, inherits = FALSE)
    if (length(existing) > 0L) rm(list = existing, envir = .GlobalEnv)
    withr::defer({
        present <- intersect(names, ls(envir = .GlobalEnv))
        if (length(present) > 0L) rm(list = present, envir = .GlobalEnv)
        list2env(values, envir = .GlobalEnv)
    })
    names
}

nondiaCloseout028ResourceMetrics <- function() {
    process <- ps::ps_handle()
    list(
        rss_bytes = as.numeric(ps::ps_memory_info(process)[["rss"]]),
        open_handles = as.integer(ps::ps_num_fds(process)),
        child_count = length(ps::ps_children(process, recursive = TRUE))
    )
}

nondiaCloseout028TempCount <- function(registry) {
    path <- projectRegistryPath(registry, "temporary")
    if (!dir.exists(path)) return(0L)
    length(list.files(path, all.files = TRUE, no.. = TRUE))
}

nondiaCloseout028Cycle <- function(registry, cycle, interrupt = FALSE) {
    before <- nondiaCloseout028ResourceMetrics()
    session <- initializeProjectRegistry(registry)
    open <- projectRegistrySessionResourceInfo(session)
    connection <- projectRegistrySessionConnection(session)
    result_events <- character()
    if (isTRUE(interrupt)) {
        error <- rlang::catch_cnd(projectRegistryFetchBound(
            connection,
            "SELECT * FROM closeout_missing_table WHERE cycle = ?",
            list(cycle),
            result_observer = \(event) {
                result_events <<- c(result_events, event)
            }
        ))
        stopifnot(inherits(error, "error"))
    } else {
        projectRegistryFetchBound(
            connection,
            "SELECT ? AS cycle",
            list(cycle),
            result_observer = \(event) {
                result_events <<- c(result_events, event)
            }
        )
    }
    live_results <- sum(result_events == "opened") -
        sum(result_events == "cleared")
    closeProjectRegistry(session)
    closed <- projectRegistrySessionResourceInfo(session)
    after <- nondiaCloseout028ResourceMetrics()
    list(
        before = before,
        after = after,
        open = open,
        closed = closed,
        live_results = live_results,
        writer_guard = file.exists(projectRegistryPath(registry, "owner")),
        temporary_entries = nondiaCloseout028TempCount(registry),
        interrupted = interrupt
    )
}

test_that("non-DIA closeout withholds each tuple on exact missing evidence", {
    decisions <- artifactProteomicsNonDiaCloseoutDecisions()
    supported_ids <- names(artifactProteomicsNonDiaWorkflowDescriptors())
    spectronaut_ids <- c(
        "proteomics.spectronaut.protein.lfq.v1",
        "proteomics.spectronaut.peptide.lfq.v1"
    )
    expect_setequal(names(decisions), c(supported_ids, spectronaut_ids))

    descriptors <- artifactProteomicsNonDiaWorkflowDescriptors()
    for (capability_id in supported_ids) {
        decision <- decisions[[capability_id]]
        descriptor <- descriptors[[capability_id]]
        evidence <- nondiaCloseout028ReadEvidence(decision)
        expect_identical(decision$capability_id, descriptor$descriptor_id)
        expect_identical(
            decision$capability_version,
            descriptor$descriptor_version
        )
        expect_identical(
            decision$descriptor_digest,
            descriptor$descriptor_digest
        )
        expect_identical(decision$promotion_status, "withheld")
        expect_identical(
            decision$reason_code,
            "paired_artifact_performance_evidence_absent"
        )
        expect_identical(decision$effective_default_backend, "memory")
        expect_true(decision$forced_artifact_canary)
        expect_false(decision$certification$auto_eligible)
        expect_identical(decision$certification$status, "dual_write")
        expect_identical(decision$maximum_rollout, "dual_write")
        expect_identical(decision$rollback$target_backend, "memory")
        expect_identical(evidence$status, "passed")
        expect_identical(evidence$summary$completed, 5L)
        expect_length(evidence$runs, 5L)
        expect_true(all(vapply(evidence$runs, \(run) {
            identical(run$status, "passed") &&
                identical(run$worker$session_evidence$status, "memory_only") &&
                identical(run$metrics$peak_artifact_disk_bytes, 0L)
        }, logical(1))))
        expect_identical(
            decision$performance_evidence$paired_artifact$valid_pairs,
            0L
        )
        expect_identical(
            decision$performance_evidence$paired_artifact$required_pairs,
            5L
        )
    }

    for (capability_id in spectronaut_ids) {
        decision <- decisions[[capability_id]]
        expect_identical(decision$support_status, "advertised_unverified")
        expect_identical(
            decision$reason_code,
            "scientific_support_evidence_absent"
        )
        expect_false(decision$forced_artifact_canary)
        expect_identical(decision$maximum_rollout, "none")
        expect_identical(decision$effective_default_backend, "memory")
    }

    first <- decisions[[supported_ids[[1L]]]]
    forged <- first
    forged$performance_evidence$paired_artifact$status <- "passed"
    forged$performance_evidence$paired_artifact$valid_pairs <- 5L
    expect_error(
        validateArtifactProteomicsNonDiaCloseout(forged),
        class = "multischolar_incomplete_prot_nondia_closeout_evidence"
    )
    expect_identical(
        decisions[[supported_ids[[2L]]]]$reason_code,
        "paired_artifact_performance_evidence_absent"
    )
})

test_that("closeout policy preserves memory and bounded forced canaries", {
    decisions <- artifactProteomicsNonDiaCloseoutDecisions()
    descriptors <- artifactProteomicsNonDiaWorkflowDescriptors()
    capabilities <- mergeWorkflowDescriptorCapabilities()
    for (descriptor in descriptors) {
        identity <- descriptor$identity
        automatic <- resolveWorkflowBackend(
            identity,
            requested_backend = "auto",
            project_state = "new",
            capabilities = capabilities
        )
        expect_identical(automatic$effective_backend, "memory")
        expect_identical(automatic$reason_code, "auto_capability_not_promoted")
        legacy <- resolveWorkflowBackend(
            identity,
            requested_backend = "auto",
            project_state = "legacy_memory",
            capabilities = capabilities
        )
        expect_identical(legacy$effective_backend, "memory")
        expect_identical(legacy$reason_code, "auto_preserve_legacy_memory")
        memory <- resolveWorkflowBackend(
            identity,
            requested_backend = "memory",
            project_state = "new",
            capabilities = capabilities
        )
        expect_identical(memory$effective_backend, "memory")
        forced <- resolveWorkflowBackend(
            identity,
            requested_backend = "artifact",
            requested_rollout = "dual_write",
            project_state = "new",
            migration_requested = TRUE,
            capabilities = capabilities
        )
        expect_identical(forced$effective_backend, "artifact")
        expect_identical(forced$effective_rollout, "dual_write")
        expect_error(
            resolveWorkflowBackend(
                identity,
                requested_backend = "artifact",
                requested_rollout = "read_through",
                project_state = "new",
                migration_requested = TRUE,
                capabilities = capabilities
            ),
            class = "multischolar_artifact_rollout_not_certified"
        )
        reopened <- resolveWorkflowBackend(
            identity,
            requested_backend = "auto",
            requested_rollout = "evict",
            project_state = "artifact_valid",
            capabilities = capabilities
        )
        expect_identical(reopened$effective_backend, "artifact")
        expect_identical(reopened$effective_rollout, "dual_write")
        decision <- decisions[[descriptor$descriptor_id]]
        expect_identical(
            reopened$capability_version,
            decision$capability_version
        )
    }

    spectronaut <- Filter(
        \(capability) identical(capability$identity$input_format, "spectronaut"),
        workflowCapabilityCatalogue()
    )
    for (capability in spectronaut) {
        expect_error(
            resolveWorkflowBackend(
                capability$identity,
                requested_backend = "artifact",
                project_state = "new",
                migration_requested = TRUE,
                capabilities = capabilities
            ),
            class = "multischolar_artifact_not_certified"
        )
    }
})

test_that("exact tuple registry cycles close without leaks or contamination", {
    nondiaCloseout028SkipDependencies()
    global_names <- nondiaCloseout028PreserveGlobals()
    descriptors <- artifactProteomicsNonDiaWorkflowDescriptors()
    specs <- artifactProteomicsNonDiaDescriptorSpecs()
    roots <- character()
    database_paths <- character()
    for (descriptor in descriptors) {
        root <- withr::local_tempdir(
            pattern = paste0(
                "nondia-028-lifecycle-",
                descriptor$identity$input_format,
                "-"
            )
        )
        roots <- c(roots, normalizePath(root, winslash = "/"))
        paths <- buildArtifactPaths(root, list(
            omic_type = descriptor$identity$omic_type,
            omic_label = paste0(
                "closeout_",
                descriptor$identity$input_format
            ),
            workflow_slug = descriptor$identity$workflow_slug
        ))
        registry <- newProjectRegistry(
            paths,
            paste0("closeout-", descriptor$identity$input_format)
        )
        database_paths <- c(
            database_paths,
            normalizePath(
                projectRegistryPath(registry, "database"),
                winslash = "/",
                mustWork = FALSE
            )
        )
        cycles <- lapply(seq_len(5L), \(cycle) {
            nondiaCloseout028Cycle(
                registry,
                cycle,
                interrupt = identical(cycle, 3L)
            )
        })
        for (cycle in cycles) {
            expect_identical(cycle$open$connection_count, 1L)
            expect_identical(cycle$open$writer_guard_count, 1L)
            expect_identical(cycle$closed$connection_count, 0L)
            expect_identical(cycle$closed$writer_guard_count, 0L)
            expect_identical(cycle$live_results, 0L)
            expect_false(cycle$writer_guard)
            expect_identical(cycle$temporary_entries, 0L)
            expect_identical(cycle$after$child_count, 0L)
        }
        rss <- vapply(
            cycles[-1L],
            \(cycle) cycle$after$rss_bytes,
            numeric(1)
        )
        plateau <- diff(range(rss))
        limit <- specs[[descriptor$descriptor_id]]$thresholds[[
            "maximum_retained_tree_rss_bytes"
        ]]
        expect_lte(plateau, limit)
        handles <- vapply(
            cycles,
            \(cycle) cycle$after$open_handles,
            integer(1)
        )
        expect_lte(diff(range(handles)), 1L)
    }
    expect_identical(anyDuplicated(roots), 0L)
    expect_identical(anyDuplicated(database_paths), 0L)
    expect_true(all(vapply(seq_along(roots), \(index) {
        artifactPathIsContained(database_paths[[index]], roots[[index]])
    }, logical(1))))
    expect_false(any(vapply(
        global_names,
        exists,
        logical(1),
        envir = .GlobalEnv,
        inherits = FALSE
    )))
})

test_that("closeout rollback contracts remain exact and data only", {
    decisions <- artifactProteomicsNonDiaCloseoutDecisions()
    descriptors <- artifactProteomicsNonDiaWorkflowDescriptors()
    for (descriptor in descriptors) {
        decision <- decisions[[descriptor$descriptor_id]]
        contract <- protNonDiaPayloadEvictionContract(descriptor)
        expect_identical(
            contract$descriptor$descriptor_id,
            descriptor$descriptor_id
        )
        expect_identical(
            decision$rollback$strategy_id,
            "force_memory_ignore_tuple_generations"
        )
        expect_identical(decision$rollback$target_backend, "memory")
        expect_identical(
            descriptor$migration$from_backend,
            "memory"
        )
        expect_identical(
            descriptor$migration$to_backend,
            "artifact"
        )
        expect_silent(artifactResourceDataOnly(decision))
        expect_silent(artifactResourceDataOnly(contract))
    }
    expect_identical(protNonDiaSessionMode("restore"), "enabled")
    expect_true(is.function(reconstructProtNonDiaCompatibilitySession))
    expect_true(is.function(writeProtNonDiaCompatibilitySession))
})
