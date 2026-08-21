omics019ArtifactWorkflow <- function(root) {
    paths <- list(
        base_dir = root,
        project_id = "omics-019-project",
        omic_type = "proteomics",
        omic_label = "dia-explicit-context",
        source_dir = file.path(root, "source"),
        results_dir = file.path(root, "results")
    )
    dir.create(paths$source_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(paths$results_dir, recursive = TRUE, showWarnings = FALSE)
    workflow <- new.env(parent = emptyenv())
    workflow$workflow_context <- createWorkflowContext(
        paths,
        "proteomics",
        paths$omic_label,
        storage_policy = list(
            requested_backend = "artifact",
            requested_rollout = "dual_write",
            migration_requested = TRUE,
            project_id = paths$project_id
        )
    )
    workflow$state_manager <- WorkflowState$new()
    workflow$state_manager$setWorkflowType("DIA")
    workflow$data_format <- "diann"
    workflow$data_type <- "peptide"
    prepareProtDiaArtifactContext(workflow)
    workflow
}

omics019PreserveGlobals <- function(global_names) {
    previous <- lapply(global_names, function(name) {
        present <- exists(name, envir = .GlobalEnv, inherits = FALSE)
        list(
            present = present,
            value = if (present) get(name, envir = .GlobalEnv) else NULL
        )
    })
    names(previous) <- global_names
    withr::defer({
        for (name in names(previous)) {
            if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
                rm(list = name, envir = .GlobalEnv)
            }
            if (previous[[name]]$present) {
                assign(name, previous[[name]]$value, envir = .GlobalEnv)
            }
        }
    }, envir = parent.frame())
    invisible(previous)
}

test_that("proteomics global usage is exhaustively classified", {
    dependency_inventory <- protContextDependencyInventory()
    usage_inventory <- protContextGlobalUsageInventory()

    expect_setequal(
        dependency_inventory$dependency,
        c(
            "project_paths", "design_matrix", "contrasts", "config",
            "annotations", "sequences", "experiment_paths",
            "filtering_progress"
        )
    )
    required_fields <- c(
        "classification", "explicit_owner", "legacy_consumers",
        "warning_policy", "removal_prerequisite"
    )
    for (field in required_fields) {
        expect_true(all(nzchar(dependency_inventory[[field]])), info = field)
    }
    for (field in c(
        "functions", "dependencies", "classification", "legacy_consumers",
        "artifact_policy", "warning_policy", "removal_prerequisite"
    )) {
        expect_true(all(nzchar(usage_inventory[[field]])), info = field)
    }

    source_dir <- testthat::test_path("..", "..", "R")
    skip_if_not(dir.exists(source_dir), "package R sources are unavailable")
    files <- list.files(
        source_dir,
        pattern = "^(func_prot|mod_prot).*\\.R$",
        full.names = TRUE
    )
    global_files <- files[vapply(files, function(path) {
        any(grepl(
            "\\.GlobalEnv|globalenv\\s*\\(",
            readLines(path, warn = FALSE)
        ))
    }, logical(1L))]
    expect_setequal(basename(global_files), usage_inventory$file)
})

test_that("artifact context always wins without inspecting stale globals", {
    workflow <- omics019ArtifactWorkflow(withr::local_tempdir())
    expect_true(protContextArtifactOwned(workflow))

    dependencies <- protContextDependencyInventory()$dependency
    stale_globals <- new.env(parent = emptyenv())
    for (dependency in dependencies) {
        specification <- protContextSpecification(dependency)
        assign(
            specification$global_name,
            structure(dependency, owner = "stale_global"),
            envir = stale_globals
        )
        explicit <- structure(dependency, owner = "explicit_context")
        resolved <- resolveProtContextDependency(
            dependency,
            workflow_data = workflow,
            explicit_value = explicit,
            required = TRUE,
            global_env = stale_globals,
            exists_fn = function(...) stop("artifact inspected a global"),
            get_fn = function(...) stop("artifact read a global")
        )
        expect_identical(resolved$value, explicit, info = dependency)
        expect_true(resolved$artifact_owned, info = dependency)
        expect_false(resolved$legacy_adapter_used, info = dependency)
    }

    writes <- 0L
    for (dependency in dependencies) {
        published <- publishProtContextLegacyGlobal(
            dependency,
            "forbidden",
            workflow_data = workflow,
            global_env = stale_globals,
            assign_fn = function(...) writes <<- writes + 1L
        )
        expect_false(published, info = dependency)
    }
    expect_identical(writes, 0L)
})

test_that("missing artifact dependencies fail with an ownership diagnosis", {
    workflow <- omics019ArtifactWorkflow(withr::local_tempdir())
    expect_error(
        resolveProtContextDependency(
            "contrasts",
            workflow_data = workflow,
            required = TRUE
        ),
        class = "multischolar_missing_prot_context_dependency"
    )
})

test_that("legacy adapters preserve precedence and can be disabled", {
    legacy_globals <- new.env(parent = emptyenv())
    assign("contrasts_tbl", "legacy", envir = legacy_globals)
    resolved <- resolveProtContextDependency(
        "contrasts",
        explicit_value = "explicit",
        global_env = legacy_globals
    )
    expect_identical(resolved$value, "legacy")
    expect_true(resolved$legacy_adapter_used)

    withr::local_options(list(
        multischolar.prot.legacy_global_adapters = "disabled"
    ))
    resolved <- resolveProtContextDependency(
        "contrasts",
        explicit_value = "explicit",
        global_env = legacy_globals,
        exists_fn = function(...) stop("disabled adapter inspected a global"),
        get_fn = function(...) stop("disabled adapter read a global")
    )
    expect_identical(resolved$value, "explicit")
    expect_false(resolved$legacy_adapter_used)
    expect_false(publishProtContextLegacyGlobal(
        "contrasts",
        "replacement",
        global_env = legacy_globals,
        assign_fn = function(...) stop("disabled adapter wrote a global")
    ))
    expect_length(captureProtContextLegacyGlobals("contrasts"), 0L)
    expect_null(restoreProtContextLegacyGlobals(list(
        contrasts = list(exists = TRUE, value = "replacement")
    )))
    expect_identical(get("contrasts_tbl", envir = legacy_globals), "legacy")
})

test_that("artifact session application does not mutate scientific globals", {
    workflow <- omics019ArtifactWorkflow(withr::local_tempdir())
    session_data <- module_ci_prot_da_session_payload()
    restored_manager <- WorkflowState$new()
    restored_manager$setWorkflowType("DIA")
    restored_manager$saveState(
        "correlation_filtered",
        session_data$current_s4_object,
        session_data$config_list,
        "OMICS-ART-019 explicit context fixture"
    )
    workflow$tab_status <- list(
        normalization = "pending",
        differential_expression = "locked"
    )
    workflow$state_update_trigger <- NULL
    da_data <- new.env(parent = emptyenv())
    da_data$current_s4_object <- NULL
    da_data$contrasts_available <- NULL
    da_data$formula_from_s4 <- NULL
    bundle <- list(
        format = "artifact",
        session_data = session_data,
        state_manager = restored_manager,
        workflow_context = workflow$workflow_context,
        manifest = list(schema_version = "test"),
        uniprot_dat_cln = data.frame(marker = "explicit"),
        aa_seq_tbl_final = data.frame(sequence = "PEPTIDE")
    )

    global_names <- c(
        "contrasts_tbl", "config_list", "uniprot_dat_cln",
        "aa_seq_tbl_final", "project_dirs"
    )
    omics019PreserveGlobals(global_names)
    for (name in global_names) {
        assign(name, structure(name, owner = "stale_global"), envir = .GlobalEnv)
    }
    before <- lapply(global_names, get, envir = .GlobalEnv)

    applied <- applyProtDaSessionBundle(workflow, da_data, bundle)

    after <- lapply(global_names, get, envir = .GlobalEnv)
    expect_identical(after, before)
    expect_identical(workflow$contrasts_tbl, session_data$contrasts_tbl)
    expect_identical(workflow$config_list, session_data$config_list)
    expect_identical(workflow$uniprot_dat_cln, bundle$uniprot_dat_cln)
    expect_identical(workflow$aa_seq_tbl_final, bundle$aa_seq_tbl_final)
    expect_identical(applied$annotations, bundle$uniprot_dat_cln)
})

test_that("proteomics public module signature remains stable", {
    expect_identical(
        names(formals(mod_proteomics_server)),
        c(
            "id", "project_dirs", "omic_type", "experiment_label",
            "volumes", "storage_policy"
        )
    )
})
