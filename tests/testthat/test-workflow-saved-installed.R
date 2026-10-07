installedCompatibilitySnapshot <- function(library_path) {
    callr::r(
        function(library_path) {
            .libPaths(c(library_path, .libPaths()))
            suppressPackageStartupMessages(library(MultiScholaR))
            namespace <- asNamespace("MultiScholaR")
            class_snapshot <- function(class_name) {
                definition <- methods::getClassDef(class_name)
                list(
                    slots = methods::slotNames(definition),
                    slot_classes = unname(definition@slots)
                )
            }
            generic_snapshot <- function(generic_name) {
                definition <- methods::getGeneric(generic_name)
                list(
                    formals = names(formals(definition)),
                    signatures = sort(names(methods::findMethods(generic_name)))
                )
            }
            list(
                package_path = normalizePath(find.package("MultiScholaR")),
                version = as.character(utils::packageVersion("MultiScholaR")),
                exports = sort(getNamespaceExports("MultiScholaR")),
                loaded_namespaces = sort(loadedNamespaces()),
                workflow_state_methods = sort(names(MultiScholaR::WorkflowState$public_methods)),
                workflow_state_fields = sort(names(MultiScholaR::WorkflowState$public_fields)),
                public_formals = lapply(
                    c("runTestsContrasts", "run_app"),
                    \(name) names(formals(getExportedValue("MultiScholaR", name)))
                ),
                classes = lapply(
                    c("PeptideQuantitativeData", "ProteinQuantitativeData"),
                    class_snapshot
                ),
                generics = lapply(
                    c("normaliseBetweenSamples", "differentialAbundanceAnalysis"),
                    generic_snapshot
                ),
                package_dynlibs = unname(getNamespaceInfo(namespace, "dynlibs"))
            )
        },
        args = list(library_path)
    )
}

installedSavedEntityProbe <- function(library_path, project_root, allow_native_runtime = FALSE) {
    callr::r(
        function(library_path, project_root, allow_native_runtime) {
            .libPaths(c(library_path, .libPaths()))
            suppressPackageStartupMessages(library(MultiScholaR))
            namespace <- asNamespace("MultiScholaR")
            api <- \(name) get(name, envir = namespace, inherits = FALSE)
            expected_path <- normalizePath(file.path(library_path, "MultiScholaR"))
            stopifnot(identical(normalizePath(find.package("MultiScholaR")), expected_path))

            object <- methods::new(
                "ProteinQuantitativeData",
                protein_quant_table = data.frame(
                    Protein.Ids = c("P1", "P2"),
                    S1 = c(NA_real_, Inf),
                    S2 = c(NaN, -Inf)
                ),
                design_matrix = data.frame(
                    Run = c("S1", "S2"),
                    group = ordered(c("control", "treated")),
                    replicates = c("1", "1")
                ),
                sample_id = "Run",
                group_id = "group",
                technical_replicate_id = "replicates",
                args = list(
                    normaliseBetweenSamples = list(normalisation_method = "none")
                )
            )
            methods::validObject(object)
            normalized <- normaliseBetweenSamples(object, normalisation_method = "none")
            stopifnot(
                methods::is(normalized, "ProteinQuantitativeData"),
                identical(normalized@protein_quant_table$Protein.Ids, c("P1", "P2")),
                identical(is.na(normalized@protein_quant_table$S1), c(TRUE, TRUE)),
                identical(is.na(normalized@protein_quant_table$S2), c(TRUE, TRUE))
            )

            manager <- MultiScholaR::WorkflowState$new()
            manager$saveState("raw", object, list(route = "iq"), "installed raw")
            manager$saveState("normalized", normalized, list(route = "iq"), "installed normalized")
            manager$states$raw$description <- "externally assigned description"
            stopifnot(
                identical(manager$getState("raw"), object),
                identical(manager$revertToState("raw"), object),
                identical(manager$states$raw$description, "externally assigned description")
            )
            legacy_path <- file.path(project_root, "legacy-session.rds")
            saveRDS(manager, legacy_path)
            legacy <- readRDS(legacy_path)
            stopifnot(
                identical(legacy$getState("raw"), object),
                identical(legacy$state_history, manager$state_history)
            )

            history <- api("prepareDIAWorkflowHistory")(manager$states)
            identity <- list(
                input = "installed-small-v1",
                design = "two-sample-v1",
                parameters = "normalisation-none-v1",
                selectors = "ordered-samples-v1",
                software = "installed-candidate-v1"
            )
            metadata <- list(
                class = class(history),
                shape = as.integer(dim(history)),
                factors = list(group = ordered(c("control", "treated")))
            )
            store <- api("newWorkflowEntityStore")(
                project_root,
                "installed-project",
                "DIA",
                max_cache_bytes = 67108864
            )
            reference <- api("saveWorkflowCheckpoint")(
                store,
                history,
                "installed-history",
                identity,
                metadata
            )
            restored <- api("hydrateSavedWorkflowCheckpoint")(store, reference)
            restored_again <- api("hydrateSavedWorkflowCheckpoint")(store, reference)
            stopifnot(
                identical(serialize(restored, NULL), serialize(history, NULL)),
                identical(restored_again, restored)
            )
            resident <- api("workflowEntityCacheInfo")(store)
            stopifnot(
                resident$entries == 1L,
                resident$cache_hits == 1L,
                resident$payload_reads == 1L,
                resident$payload_decodes == 1L,
                resident$bytes <= 67108864,
                resident$open_resources == 0L
            )
            api("releaseWorkflowEntityCache")(store)
            stopifnot(api("workflowEntityCacheInfo")(store)$entries == 0L)
            api("closeWorkflowEntityStore")(store)
            api("closeWorkflowEntityStore")(store)

            reopened <- api("newWorkflowEntityStore")(project_root, "installed-project", "DIA")
            recovered <- api("recoverSavedWorkflowRequest")(reopened, "installed-history")
            stopifnot(
                identical(recovered, reference),
                identical(
                    serialize(api("hydrateSavedWorkflowCheckpoint")(reopened, recovered), NULL),
                    serialize(history, NULL)
                ),
                api("workflowEntityCacheInfo")(reopened)$max_bytes == 0
            )

            foreign_root <- tempfile("foreign-project-")
            dir.create(foreign_root)
            foreign <- api("newWorkflowEntityStore")(
                foreign_root,
                "foreign-project",
                "DIA"
            )
            scope_error <- tryCatch(
                {
                    api("hydrateSavedWorkflowCheckpoint")(foreign, reference)
                    NULL
                },
                saved_entity_scope = \(error) error
            )
            connection <- file(file.path(project_root, "unsupported.txt"), "w")
            on.exit(close(connection), add = TRUE)
            unsupported_error <- tryCatch(
                {
                    api("saveWorkflowCheckpoint")(
                        reopened,
                        list(connection = connection),
                        "unsupported",
                        identity,
                        list(class = "list", shape = integer(), factors = list())
                    )
                    NULL
                },
                saved_entity_unsupported = \(error) error
            )
            stopifnot(
                inherits(scope_error, "saved_entity_scope"),
                inherits(unsupported_error, "saved_entity_unsupported"),
                inherits(get("digest_impl", asNamespace("digest")), "NativeSymbolInfo")
            )
            native_routines <- if (allow_native_runtime) {
                dll <- getLoadedDLLs()[["MultiScholaR"]]
                stopifnot(!is.null(dll), !dll[["dynamicLookup"]],
                    identical(normalizePath(dirname(dirname(dll[["path"]]))), expected_path))
                names(getDLLRegisteredRoutines(dll)$.Call)
            } else {
                stopifnot(length(getNamespaceInfo(namespace, "dynlibs")) == 0L)
                character()
            }
            api("closeWorkflowEntityStore")(reopened)
            api("closeWorkflowEntityStore")(foreign)

            files <- list.files(project_root, recursive = TRUE, all.files = TRUE,
                full.names = TRUE, include.dirs = FALSE)
            default_root <- tempfile("default-cache-")
            dir.create(default_root)
            default_store <- api("newWorkflowEntityStore")(
                default_root,
                "default-cache",
                "DIA"
            )
            default_cache_bytes <- api("workflowEntityCacheInfo")(default_store)$max_bytes
            api("closeWorkflowEntityStore")(default_store)
            list(
                package_path = expected_path,
                version = as.character(utils::packageVersion("MultiScholaR")),
                s4_dispatch = TRUE,
                public_workflow_state = TRUE,
                legacy_round_trip = TRUE,
                digest_native_symbol = TRUE,
                package_native_routines = length(native_routines),
                native_routines = native_routines,
                cache = resident,
                default_cache_bytes = default_cache_bytes,
                artifact_files = length(files),
                artifact_bytes = sum(file.info(files)$size),
                exact_round_trip = TRUE,
                explicit_rejections = c("saved_entity_scope", "saved_entity_unsupported")
            )
        },
        args = list(library_path, project_root, allow_native_runtime)
    )
}

test_that("fresh installed saved-entity utility preserves package compatibility", {
    candidate_library <- Sys.getenv("MULTISCHOLAR_DIA008_LIBRARY", unset = "")
    baseline_library <- Sys.getenv("MULTISCHOLAR_DIA008_B0_LIBRARY", unset = "")
    receipt_path <- Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT", unset = "")
    native_runtime <- identical(Sys.getenv("MULTISCHOLAR_PERF108_COMPATIBILITY"), "true")
    skip_if(!nzchar(candidate_library) || !nzchar(baseline_library),
        "DIA-008 installed libraries are not bound")
    expect_true(dir.exists(candidate_library))
    expect_true(dir.exists(baseline_library))

    baseline <- installedCompatibilitySnapshot(baseline_library)
    candidate <- installedCompatibilitySnapshot(candidate_library)
    expect_identical(candidate$exports, baseline$exports)
    expect_identical(candidate$loaded_namespaces, baseline$loaded_namespaces)
    expect_identical(candidate$workflow_state_methods, baseline$workflow_state_methods)
    expect_identical(candidate$workflow_state_fields, baseline$workflow_state_fields)
    expect_identical(candidate$public_formals, baseline$public_formals)
    expect_identical(candidate$classes, baseline$classes)
    expect_identical(candidate$generics, baseline$generics)
    if (!native_runtime) {
        expect_identical(candidate$package_dynlibs, baseline$package_dynlibs)
        expect_length(candidate$package_dynlibs, 0L)
    }

    project_root <- withr::local_tempdir(pattern = "dia008-installed-")
    probe <- installedSavedEntityProbe(candidate_library, project_root, native_runtime)
    expect_true(probe$s4_dispatch)
    expect_true(probe$public_workflow_state)
    expect_true(probe$legacy_round_trip)
    expect_true(probe$digest_native_symbol)
    if (native_runtime) {
        expect_setequal(probe$native_routines, c("workflow_reclaim_pages",
            "workflow_release_file_pages", "workflow_allocator_snapshot",
            "workflow_dia_ordered_new", "workflow_dia_ordered_push", "workflow_dia_ordered_close"))
    } else expect_identical(probe$package_native_routines, 0L)
    expect_identical(probe$default_cache_bytes, 0)
    expect_lte(probe$cache$bytes, 67108864)
    expect_identical(probe$cache$open_resources, 0L)
    expect_true(probe$exact_round_trip)
    expect_setequal(
        probe$explicit_rejections,
        c("saved_entity_scope", "saved_entity_unsupported")
    )
    expect_lt(probe$artifact_bytes, 1048576)

    if (nzchar(receipt_path)) {
        jsonlite::write_json(
            list(
                schema = "multischolar.dia008.installed_qualification.v1",
                status = "passed",
                candidate = candidate,
                baseline = baseline,
                probe = probe,
                source_test = normalizePath(testthat::test_path(
                    "test-workflow-saved-installed.R"
                )),
                performance_claim = FALSE
            ),
            receipt_path,
            auto_unbox = TRUE,
            pretty = TRUE
        )
    }
})
