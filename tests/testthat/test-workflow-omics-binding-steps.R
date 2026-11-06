# KIT-094: the shared omics core reads each binding's own steps. A stub third binding, injected by mocking
# `omicsDeclarationBindings` (no registration API), has its steps, worker steps, default variants and publications
# registered, its filtering, native and dependency steps read by the worker's compute handler, its steps replayed, its
# record outcomes hook called, its input role checked, and its workflow type accepted with its route. The lists the
# bindings now own are characterised for lipidomics and metabolomics: the characterisation equals the snapshot taken on
# the base commit before any edit (tests/testthat/fixtures/kit094-characterisation.rds, written by
# docs/qa/workflow-execution/kit/kit-094/characterise.R from the functions below).
library(testthat)

#' The workflow types the allowlists and routes are characterised on.
kit094Types <- c("DIA", "DIA_LIMPA", "LFQ", "TMT", "LIPIDOMICS", "METABOLOMICS", "TRANSCRIPTOMICS", "PROTEOMICS",
    "lipidomics", "STUBOMICS")

#' The record input roles put through the replay's input check.
kit094Roles <- c("diann_report", "maxquant_report", "fragpipe_report", "pd_tmt_report", "spectronaut_report", "fasta",
    "design_file", "lipid_assay", "metab_assay", "lipid_assay_1", "lipid_assay_12", "metab_assay_3", "metab_assay_",
    "lipid_assay1", "lipid_assayx", "xlipid_assay_1", "LIPID_ASSAY_1", "lipid", "assay_1", "trx_assay_1",
    "stub_assay_1", "other")

#' Whether the replay checks a record input of this role: a checked role with no file is refused with
#' `workflow_replay_inputs` before any step; an unchecked one is skipped and the replay fails later, on the store.
kit094InputRoleChecked <- function(role) {
    record <- list(inputs = list(list(role = role, sha256 = "0")), steps = list())
    testthat::with_mocked_bindings(
        tryCatch({
            replayWorkflowAnalysisRecord(record, tempfile("kit094-"), character(), list(), tempfile("kit094-"))
            "unchecked"
        }, workflow_replay_inputs = function(error) "checked", error = function(error) "unchecked"),
        validateWorkflowAnalysisRecord = function(record) list(valid = TRUE), .package = "MultiScholaR")
}

#' What the reopen, the session's `setWorkflowType` and the route lookup do with each type.
kit094TypeSets <- function(root) {
    # The message's whitespace is normalised: cli wraps it at the console's width.
    outcome <- function(expr) tryCatch(expr, error = function(error) paste("error", class(error)[[1L]],
        gsub("\\s+", " ", conditionMessage(error))))
    reopen <- vapply(kit094Types, function(type) {
        directory <- file.path(root, "reopen", type)
        dir.create(directory, recursive = TRUE)
        closeWorkflowStateOwner(newWorkflowStateOwner(newWorkflowEntityStore(directory, "kit094", type)))
        owner <- newWorkflowStateOwner(newWorkflowEntityStore(directory, "kit094", type), reopen = TRUE)
        on.exit(closeWorkflowStateOwner(owner))
        owner$.context$omic_type %||% NA_character_
    }, character(1))
    directory <- file.path(root, "session")
    dir.create(directory, recursive = TRUE)
    owner <- newWorkflowStateOwner(newWorkflowEntityStore(directory, "kit094", "DIA"))
    on.exit(closeWorkflowStateOwner(owner))
    session <- vapply(kit094Types, function(type) outcome({
        owner$setWorkflowType(type)
        owner$workflow_type
    }), character(1))
    list(reopen = reopen, session = session, route = vapply(kit094Types, function(type) outcome(workflowRouteForType(type)), character(1)))
}

#' The registry fields this ticket touches, for every operation of one omic.
kit094Registry <- function(omic_type) {
    names <- Filter(function(name) identical(workflowOperationEntry(name)$omic_type, omic_type), workflowOperationNames())
    stats::setNames(lapply(names, function(name) {
        entry <- workflowOperationEntry(name)
        c(entry[c("compute", "handler", "validate", "completion", "publication", "variants", "default_variant", "step_label")],
            list(sectionless = "section" %in% names(entry)))
    }), names)
}

#' One fixture journey's analysis record: its steps, parameters and outcomes, its inputs and its DA, with the software
#' identity removed (the entity ids, which carry it, read `<entity>`).
kit094Record <- function(fixture, journey) {
    on.exit(closeDIAStateAdapter(fixture$owner), add = TRUE)
    journey(fixture)
    record <- buildWorkflowAnalysisRecord(fixture$owner)
    entity <- function(x) {
        if (is.list(x)) {
            x[] <- lapply(x, entity)
            return(x)
        }
        if (is.character(x)) x[grepl("^[0-9a-f]{64}$", x)] <- "<entity>"
        x
    }
    keep <- c("order", "state", "operation", "variant", "preceded_by", "step_label", "dimensions", "outcomes", "description",
        "parameters", "stored_as_object")
    list(omic_type = record$omic_type, inputs = lapply(record$inputs, `[`, c("role", "sha256")),
        steps = entity(lapply(record$steps, `[`, keep)),
        da = entity(record$analyses$differential_abundance[c("operation", "parameters", "contrasts")]))
}

#' Everything the frozen dump does not cover of the lists that move into the bindings.
kit094Characterisation <- function(root) {
    bindings <- omicsDeclarationBindings()
    omic_types <- vapply(bindings, `[[`, "", "omic_type")
    list(
        report_sources = workflowRecordReportSources,
        types = kit094TypeSets(root),
        input_roles = vapply(kit094Roles, kit094InputRoleChecked, character(1)),
        steps = lapply(bindings, function(binding) names(omicsWorkflowOperations(binding))),
        registry = lapply(stats::setNames(omic_types, omic_types), kit094Registry),
        records = list(
            lipidomics = kit094Record(lipidDeclarationsOwner(file.path(root, "lipid")), lipidDeclarationsJourney),
            metabolomics = kit094Record(metabDeclarationsOwner(file.path(root, "metab")), metabDeclarationsJourney)))
}

#' A third binding with steps of its own; `calls` records what replay and the record hook were handed.
kit094StubBinding <- function(calls) {
    declare <- function(binding, name, output_class = binding$class, parameters = character(), partitioned = FALSE) {
        declaration <- newOmicsOperation(binding, name, "whole_matrix", parameters = parameters,
            product_type = paste0("stub.", name), output_class = output_class,
            kernel = function(object, parameters, store) list())
        if (partitioned) declaration$partitioned <- function(parameters, store) list()
        declaration
    }
    list(omic_type = "stubomics", prefix = "stub", class = "LipidomicsAssayData", fields = c("lipid_data", "design_matrix"),
        design_input_prefix = "stub_assay_", step_labels = c(stub_filter = "1_Stub_Filtered"),
        steps = function(binding) list(declare(binding, "stub_intensity_filter", partitioned = TRUE),
            declare(binding, "stub_filter", partitioned = TRUE), declare(binding, "stub_native", output_class = NULL),
            declare(binding, "stub_dependent", parameters = "anchor"), omicsStoreOperation(binding, "stub_store")),
        worker_steps = c("stub_intensity_filter", "stub_filter", "stub_native", "stub_dependent"),
        default_variants = c(stub_filter = "partitioned"),
        publications = list(c("figures", "computeWorkflowOmicsQcFiguresAction")),
        filtering_steps = "stub_filter", native_steps = "stub_native", dependency_parameters = "anchor",
        replay = list(store = function(binding, manager, step, parameters, config, description) {
            calls$replay <- list(state = step$state, config = config, description = description, parameters = parameters)
            invisible(NULL)
        }),
        record_outcomes = function(binding, context, operation, record, outcomes) {
            outcomes$stub_hook <- operation
            outcomes
        })
}

#' Inject the stub as a third binding for the calling test, the registry rebuilt with it and rebuilt again without it
#' when the test ends.
kit094LocalStub <- function(env = parent.frame()) {
    calls <- new.env(parent = emptyenv())
    stub <- kit094StubBinding(calls)
    # Deferred first, so it runs after the mock is undone.
    withr::defer(registerBuiltinWorkflowOperations(), envir = env)
    testthat::local_mocked_bindings(omicsDeclarationBindings = function() list(lipidomics = lipidomicsDeclarationBinding,
        metabolomics = metabolomicsDeclarationBinding, stubomics = stub), .package = "MultiScholaR", .env = env)
    registerBuiltinWorkflowOperations()
    list(binding = stub, calls = calls)
}

#' KIT-095: the characterisation with exactly the transcriptomics entries taken out, for comparison with the snapshot
#' taken before transcriptomics existed. The snapshot proves lipid and metab unchanged; a third omic adds its own
#' entries and nothing else. The supported-types message is the one allowed textual change: it now lists TRANSCRIPTOMICS.
kit095WithoutTranscriptomics <- function(current, expected) {
    stripped <- current
    stripped$report_sources[c("ruv_mode", "ruviii_m_source")] <- NULL
    for (name in c("ruv_applied", "best_k", "best_percentage", "num_neg_ctrl")) {
        stripped$report_sources[[name]]$record <- Filter(function(r) !identical(r[[2L]], "trx_ruvs"), stripped$report_sources[[name]]$record)
    }
    stripped$steps$transcriptomics <- NULL
    stripped$registry$transcriptomics <- NULL
    for (set in names(stripped$types)) {
        stripped$types[[set]][["TRANSCRIPTOMICS"]] <- expected$types[[set]][["TRANSCRIPTOMICS"]]
        stripped$types[[set]] <- vapply(stripped$types[[set]], function(text) sub("LIPIDOMICS, METABOLOMICS and TRANSCRIPTOMICS",
            "LIPIDOMICS and METABOLOMICS", text, fixed = TRUE), character(1))
    }
    stripped
}

test_that("KIT-094: the characterisation of the lists the bindings own equals the snapshot taken on the base commit", {
    skip_if_not(workflowOwnerCapable())
    snapshot_file <- test_path("fixtures", "kit094-characterisation.rds")
    expect_true(file.exists(snapshot_file))
    expected <- readRDS(snapshot_file)
    root <- withr::local_tempdir()
    current <- suppressWarnings(suppressMessages(kit094Characterisation(root)))
    expect_identical(names(current), names(expected))
    # KIT-095: transcriptomics is the third binding; lipid and metab must still equal the snapshot exactly.
    stripped <- kit095WithoutTranscriptomics(current, expected)
    for (part in names(expected)) expect_identical(stripped[[part]], expected[[part]], info = part)
    # And transcriptomics is accepted: the reopen, the session and the route serve it.
    expect_identical(unname(current$types$reopen[["TRANSCRIPTOMICS"]]), "transcriptomics")
    expect_identical(unname(current$types$route[["TRANSCRIPTOMICS"]]), "transcriptomics")
    expect_true("transcriptomics" %in% names(current$steps))
})

test_that("KIT-094: a stub third binding's own steps, worker steps, default variants and publications register", {
    kit094LocalStub()
    names <- Filter(function(name) identical(workflowOperationEntry(name)$omic_type, "stubomics"), workflowOperationNames())
    expect_setequal(names, c("stub_intensity_filter", "stub_filter", "stub_native", "stub_dependent", "stub_store", "stub_figures"))
    for (name in c("stub_intensity_filter", "stub_filter", "stub_native", "stub_dependent")) {
        entry <- workflowOperationEntry(name)
        expect_identical(entry$compute, "omics", info = name)
        expect_identical(entry$handler, "computeWorkflowOmicsOperationAction", info = name)
    }
    expect_identical(workflowOperationEntry("stub_store")$compute, "in_process")
    expect_null(workflowOperationEntry("stub_store")$handler)
    # The default variant is the binding's, not a name rule: the stub's `_intensity_filter` stays `whole`.
    expect_identical(workflowOperationEntry("stub_filter")$default_variant, "partitioned")
    expect_identical(workflowOperationEntry("stub_intensity_filter")$default_variant, "whole")
    expect_identical(workflowOperationEntry("stub_filter")$step_label, "1_Stub_Filtered")
    # The two real bindings are untouched by a third.
    expect_identical(workflowOperationEntry("lipid_intensity_filter")$default_variant, "partitioned")
    expect_identical(workflowOperationEntry("metab_qc_figures")$handler, "computeWorkflowOmicsQcFiguresAction")
    # The stub's one publication, under its own handler.
    expect_true("stub_figures" %in% workflowOperationNames())
    expect_false("stub_qc_figures" %in% workflowOperationNames())
    expect_identical(workflowOperationEntry("stub_figures")$handler, "computeWorkflowOmicsQcFiguresAction")
    expect_identical(workflowOperationEntry("stub_figures")$input_class, "LipidomicsAssayData")
})

test_that("KIT-094: the worker's compute handler reads the stub's filtering, native and dependency steps", {
    kit094LocalStub()
    seen <- new.env(parent = emptyenv())
    local_mocked_bindings(
        runSavedOperation = function(manager, declaration, source, parameters, dependencies = list()) {
            seen[[declaration$name]] <- dependencies
            list(status = "computed", value = list(final = "final-reference"), reference = "science-reference")
        },
        readWorkflowStateMetadata = function(manager, reference = NULL) list(data = list(products = list())),
        omicsFilteringMetricsReference = function(binding, manager, state) paste("metrics of", binding$prefix),
        collectWorkflowNativeReferences = function(value) "native-references")
    manager <- list(.context = list(snapshot = list(products = list(p1 = "anchor-reference"))))
    request <- function(operation, parameters = list()) list(operation = operation, input_references = list(state = "source"),
        parameters = parameters, worker_cap = 1)
    expect_identical(computeWorkflowOmicsOperationAction(manager, request("stub_filter"), "whole")$display$filtering_metrics,
        "metrics of stub")
    expect_null(computeWorkflowOmicsOperationAction(manager, request("stub_intensity_filter"), "whole")$display)
    expect_identical(computeWorkflowOmicsOperationAction(manager, request("stub_native"), "whole")$native, "native-references")
    expect_identical(seen$stub_filter, list())
    computeWorkflowOmicsOperationAction(manager, request("stub_dependent", list(anchor = "p1")), "whole")
    expect_identical(seen$stub_dependent, list("anchor-reference"))
})

test_that("KIT-094: replay dispatches the stub's steps through its own entries", {
    stub <- kit094LocalStub()
    local_mocked_bindings(readWorkflowStateMetadata = function(manager, reference = NULL) list(config = list(k = 1)))
    expect_no_error(replayWorkflowOmicsStep(stub$binding, NULL, list(operation = "stub_store", state = "stub_stored",
        description = "stored"), list(a = 1)))
    expect_identical(stub$calls$replay, list(state = "stub_stored", config = list(k = 1), description = "stored",
        parameters = list(a = 1)))
    expect_error(replayWorkflowOmicsStep(stub$binding, NULL, list(operation = "stub_filter", state = "x"), list()),
        class = "workflow_replay_unsupported")
    expect_error(replayWorkflowOmicsStep(stub$binding, NULL, list(operation = "lipid_log2", state = "x"), list()),
        class = "workflow_replay_unsupported")
})

test_that("KIT-094: the record calls the stub's outcomes hook for its own operations only", {
    kit094LocalStub()
    expect_identical(workflowRecordStepOutcomes(NULL, "stub_stored", "stub_store", list(config = list()), list(rows = 3L, columns = 4L)),
        list(rows = 3L, columns = 4L, stub_hook = "stub_store"))
    expect_identical(workflowRecordStepOutcomes(NULL, "x", "lipid_log2", list(config = list()), list(rows = 1L, columns = 2L)),
        list(rows = 1L, columns = 2L))
})

test_that("KIT-094: the replay checks the stub's input roles, from its input prefix", {
    kit094LocalStub()
    expect_identical(kit094InputRoleChecked("stub_assay_1"), "checked")
    expect_identical(kit094InputRoleChecked("stub_assay"), "checked")
    expect_identical(kit094InputRoleChecked("stub_assayx"), "unchecked")
})

test_that("KIT-094: the stub's workflow type is accepted by the reopen and the session, and its route is returned", {
    skip_if_not(workflowOwnerCapable())
    kit094LocalStub()
    expect_identical(tryCatch(workflowRouteForType("STUBOMICS"), error = conditionMessage), "stubomics")
    expect_identical(workflowRouteForType("LIPIDOMICS"), "lipidomics")
    expect_identical(workflowRouteForType("DIA"), "iq")
    directory <- withr::local_tempdir()
    closeWorkflowStateOwner(newWorkflowStateOwner(newWorkflowEntityStore(directory, "kit094", "STUBOMICS")))
    owner <- newWorkflowStateOwner(newWorkflowEntityStore(directory, "kit094", "STUBOMICS"), reopen = TRUE)
    withr::defer(closeWorkflowStateOwner(owner))
    expect_identical(owner$.context$omic_type, "stubomics")
    expect_no_error(owner$setWorkflowType("STUBOMICS"))
    expect_identical(owner$workflow_type, "STUBOMICS")
})
