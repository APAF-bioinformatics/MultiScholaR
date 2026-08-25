test_that("retention and trash remain exact-workflow scoped", {
    operationalArtifactSkipDependencies()
    root <- withr::local_tempdir(pattern = "omics-art-048-retention-")
    specifications <- list(
        dia = c("proteomics", "proteomics-study", "prot_dia"),
        nondia = c("proteomics", "nondia-study", "prot_lfq"),
        metabolomics = c("metabolomics", "metabolomics-study", "metab_standard"),
        lipidomics = c("lipidomics", "lipidomics-study", "lipid_standard")
    )
    stores <- list()
    current_refs <- list()
    historical_refs <- list()

    for (name in names(specifications)) {
        specification <- specifications[[name]]
        current <- operationalArtifactFixture(
            root,
            specification[[1L]],
            specification[[2L]],
            specification[[3L]],
            paste0("current-", name)
        )
        historical <- current
        historical$logical_key$generation_id <- paste0("historical-", name)
        stores[[name]] <- current$store
        current_refs[[name]] <- artifactStoreWriteParquet(
            current$store,
            current$encoded,
            current$logical_key
        )
        historical_refs[[name]] <- artifactStoreWriteParquet(
            historical$store,
            historical$encoded,
            historical$logical_key
        )
    }

    for (name in names(stores)) {
        store <- stores[[name]]
        current <- current_refs[[name]]
        historical <- historical_refs[[name]]
        expect_error(
            artifactStoreTrash(
                store,
                current,
                current_generation_ids = current$logical_key$generation_id
            ),
            class = "multischolar_protected_artifact_generation"
        )
        expect_error(
            artifactStoreTrash(
                store,
                historical,
                referenced_generation_ids = historical$logical_key$generation_id
            ),
            class = "multischolar_protected_artifact_generation"
        )
        other_name <- setdiff(names(stores), name)[[1L]]
        expect_error(
            artifactStoreTrash(stores[[other_name]], historical),
            class = "multischolar_cross_project_artifact"
        )

        trashed <- artifactStoreTrash(store, historical)
        trash_path <- artifactStoreResolveFile(
            store,
            trashed$trash_relative_path,
            must_exist = TRUE
        )
        expect_true(artifactPathIsContained(trash_path, store$project_root))
        expect_true(file.exists(file.path(trash_path, "artifact.json")))
        expect_true(dir.exists(file.path(trash_path, "payload")))
        current_paths <- artifactStoreManagedPaths(
            store,
            current$logical_key,
            current$artifact_id
        )
        expect_true(file.exists(artifactStoreResolveFile(
            store,
            current_paths$sidecar,
            must_exist = TRUE
        )))
    }
})

test_that("unsupported destructive APIs cannot be inferred from closeout", {
    closeout <- operationalArtifactReadCloseout()
    contract <- closeout$operational_contract
    expect_setequal(
        unlist(contract$protected_reachability, use.names = FALSE),
        c(
            "current_state", "active_lineage", "retained_history",
            "artifact_dependency_edges", "reports", "exports",
            "portable_sessions"
        )
    )
    expect_identical(contract$tombstone_api, "unsupported")
    expect_identical(contract$deliberate_purge_api, "unsupported")
    namespace <- asNamespace("MultiScholaR")
    expect_false(exists("artifactStoreTombstone", envir = namespace))
    expect_false(exists("artifactStorePurge", envir = namespace))
    expect_identical(closeout$rollback$delete_immutable_generations, FALSE)
})
