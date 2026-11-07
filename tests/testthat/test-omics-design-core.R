# KIT-103 W2: the design core of lipid, metab and trx (R/func_omics_design_shared.R,
# R/utils_workflow_omics_design.R, the objects in R/utils_workflow_{lipid,metab}_declarations.R and
# the trx copies in R/func_trx_design.R). Each shared function is checked against today's lipid,
# metab and trx function on the module CI fixtures (MCI-013, MCI-019, MCI-027): the validator's
# results and messages, the initial state and its run order, the edits, the save payload, the
# commit's bus, artifact bytes, state and provenance, and the design import. The comparison blocks
# call today's `R/mod_*` design functions; when KIT-103 W4/W5 delete them, those blocks are
# quarantined or kept against frozen values (as KIT-102 did for the import core).

# The configuration, with the function fields KIT-103 W7 wired.
omicsDesignTestConfig <- function(omic_type) {
    validateOmicGuiConfig(omicGuiConfig(omic_type))
}

quietDesign <- function(...) invisible(NULL)

# Paths under a run's root, made relative, so two runs in two roots compare.
omicsDesignRelative <- function(x, root) {
    rapply(x, function(value) gsub(root, "<root>", value, fixed = TRUE), classes = "character",
        how = "replace")
}

# The files under a directory, as relative path -> md5.
omicsDesignArtifactBytes <- function(dir) {
    files <- sort(list.files(dir, recursive = TRUE))
    stats::setNames(unname(tools::md5sum(file.path(dir, files))), files)
}

# A record without its commit times.
omicsDesignWithoutTimes <- function(x) {
    if (!is.list(x)) return(x)
    x <- x[!names(x) %in% c("committed_at", "created_at", "recorded_at", "timestamp", "time")]
    if (length(x)) x[] <- lapply(x, omicsDesignWithoutTimes)
    x
}

# What the owner recorded: the history, the step provenance (without its times), the workflow's
# inputs and each state's object.
omicsDesignOwnerRecord <- function(manager) {
    snapshot <- manager$.context$snapshot
    states <- setdiff(manager$getHistory(), "initial")
    list(history = manager$getHistory(), step_provenance = omicsDesignWithoutTimes(snapshot$step_provenance),
        file_dependencies = snapshot$file_dependencies,
        states = lapply(stats::setNames(states, states), function(name) manager$getState(name)))
}

# The bus fields the design step writes.
omicsDesignBus <- function(workflow) {
    list(design_matrix = workflow$design_matrix, contrasts_tbl = workflow$contrasts_tbl,
        config_list = workflow$config_list, column_mapping = workflow$column_mapping,
        data_cln = getWorkflowOmicsBusTables(workflow, "data_cln"),
        data_tbl = getWorkflowOmicsBusTables(workflow, "data_tbl"), tab_status = workflow$tab_status)
}

# Notifications captured from today's shells (shiny's functions mocked in the calling test).
omicsDesignCaptureShiny <- function(env = parent.frame()) {
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    testthat::local_mocked_bindings(showNotification = function(ui, ..., type = "default") {
        seen$notes <- c(seen$notes, paste(as.character(ui), collapse = " "))
        invisible(NULL)
    }, showModal = function(...) invisible(NULL), removeModal = function(...) invisible(NULL),
    removeNotification = function(...) invisible(NULL), .package = "shiny", .env = env)
    seen
}

# Today's metab builder commit (the observer's flow in mod_metab_design.R), run once.
runTodayMetabDesignCommit <- function(results, workflow, paths) {
    registerMetabDesignBuilderResultsObserver(builderResultsRv = function() results,
        workflowData = workflow, experimentPaths = paths, qcTrigger = function(value) NULL,
        observeEventFn = function(eventExpr, handlerExpr, ...) handlerExpr)
}

# The fixtures of one omic: payload, workflow and paths (MCI-019 lipid, MCI-013 metab).
omicsDesignFixture <- function(omic_type, layout = "lc_pair", ...) {
    if (identical(omic_type, "lipidomics")) {
        payload <- module_ci_lipid_design_payload(layout = layout, ...)
        list(payload = payload, workflow = module_ci_lipid_design_workflow(payload),
            paths = module_ci_lipid_design_paths())
    } else {
        payload <- module_ci_metab_design_payload(layout = layout, ...)
        list(payload = payload, workflow = module_ci_metab_design_workflow(payload),
            paths = module_ci_metab_design_paths())
    }
}

# ------------------------------------------------------------------------------------------------
# C1: the validator.
# ------------------------------------------------------------------------------------------------

test_that("C1 the shared validator is today's lipid and metab validator, labels and messages included", {
    skip("quarantined: Q-35 (KIT-103 W4) it compared the core with today's lipid design functions, deleted by W4 (and metab's, which W5 deletes); the comparison passed at fe1a9742 (kit-103 W2); the core is covered by the live blocks of this file, test-omic-design-module.R, MCI-013, MCI-019, MCI-027 and the entry lanes")
    families <- list(
        lipidomics = list(rules = omicsDesignRules("lipidomics"), alignment = validateLipidDesignAlignment,
            contrasts = validateLipidDesignContrasts, preflight = validateLipidDesignDaPreflight,
            patterns = lipidDesignNonSamplePatterns(), payload = module_ci_lipid_design_payload,
            assays = module_ci_lipid_design_assays, matrix = module_ci_lipid_design_matrix,
            mapping = module_ci_lipid_design_column_mapping,
            contrast = module_ci_lipid_design_contrasts),
        metabolomics = list(rules = omicsDesignRules("metabolomics"), alignment = validateMetabDesignAlignment,
            contrasts = validateMetabDesignContrasts, preflight = validateMetabDesignDaPreflight,
            patterns = metabDesignNonSamplePatterns(), payload = module_ci_metab_design_payload,
            assays = module_ci_metab_design_assays, matrix = module_ci_metab_design_matrix,
            mapping = module_ci_metab_design_column_mapping,
            contrast = module_ci_metab_design_contrasts))
    for (omic in names(families)) {
        f <- families[[omic]]
        # Both pattern lists are kept exactly; the id key is the configuration's.
        expect_identical(f$rules$non_sample_patterns, f$patterns, info = omic)
        expect_identical(f$rules$id_key, omicGuiConfig(omic)$design$id_mapping_key, info = omic)
        design <- f$matrix("two_group")
        mapping <- f$mapping(design$Run)
        no_samples_mapping <- mapping
        no_samples_mapping$sample_columns <- NULL
        odd_assay <- f$assays("lc_pair", samples = design$Run)
        # Columns that one list drops and the other keeps, and a non-numeric one.
        odd_assay$LCMS_Pos$Name <- 1
        odd_assay$LCMS_Pos$Class <- 2
        odd_assay$LCMS_Pos$Annotation_extra <- 3
        odd_assay$LCMS_Pos$`Ion mode` <- 4
        odd_assay$LCMS_Pos$`Peak ID` <- 5
        alignment_cases <- list(
            lc_pair = list(design, f$assays("lc_pair", samples = design$Run), mapping),
            gc = list(design, f$assays("gc", samples = design$Run), mapping),
            reordered = list(design, f$assays("combined", samples = design$Run, reorder_secondary = TRUE), mapping),
            missing = list(design, f$assays("lc_pair", samples = design$Run, missing_sample = "WT_2"), mapping),
            extra = list(design, f$assays("lc_pair", samples = design$Run, extra_sample = "EXTRA_1"), mapping),
            case = list(design, f$assays("lc_pair", samples = design$Run, case_varied = TRUE), mapping),
            inferred = list(design, odd_assay, no_samples_mapping),
            inferred_null_mapping = list(design, odd_assay, NULL),
            unnamed = list(design, unname(f$assays("lc_pair", samples = design$Run)), mapping),
            duplicate_names = list(design, stats::setNames(f$assays("lc_pair", samples = design$Run), c("A", "A")), mapping),
            not_frame = list(design, list(A = "x"), mapping),
            empty_list = list(design, list(), mapping),
            empty_design = list(design[0, ], f$assays("lc_pair", samples = design$Run), mapping),
            no_run = list(design[setdiff(names(design), "Run")], f$assays("lc_pair", samples = design$Run), mapping),
            blank_group = list(f$matrix("missing_group"), f$assays("lc_pair", samples = design$Run), mapping),
            no_replicates = list(design[setdiff(names(design), "replicates")], f$assays("lc_pair", samples = design$Run), mapping))
        for (case in names(alignment_cases)) {
            args <- alignment_cases[[case]]
            expect_identical(validateOmicsDesignAlignment(args[[1]], args[[2]], f$rules, columnMapping = args[[3]]),
                f$alignment(args[[1]], args[[2]], columnMapping = args[[3]]), info = paste(omic, case))
        }
        expect_identical(
            validateOmicsDesignAlignment(design, f$assays("lc_pair", samples = design$Run), f$rules, mapping,
                requireGroups = FALSE, requireReplicates = FALSE),
            f$alignment(design, f$assays("lc_pair", samples = design$Run), mapping, requireGroups = FALSE,
                requireReplicates = FALSE), info = omic)
        contrast_cases <- c("raw_terms", "friendly", "reversed", "duplicate", "invalid_term", "empty", "no_contrasts",
            "valid")
        for (kind in contrast_cases) for (allow in c(FALSE, TRUE)) for (formula in c("~ 0 + group", "~ group", "", "~ 0 + missing")) {
            expect_identical(validateOmicsDesignContrasts(design, f$contrast(kind), f$rules, formula, allowEmpty = allow),
                f$contrasts(design, f$contrast(kind), formula, allowEmpty = allow),
                info = paste(omic, kind, allow, formula))
        }
        friendly_dup <- data.frame(contrasts = c("groupKO-groupWT", "groupWT-groupKO"),
            friendly_names = c("X", "X"), stringsAsFactors = FALSE)
        eq_input <- data.frame(contrasts = c("KO_vs_WT=groupKO-groupWT", " "), stringsAsFactors = FALSE)
        for (tbl in list(friendly_dup, eq_input, data.frame(contrasts = "groupKO - groupWT"))) {
            expect_identical(validateOmicsDesignContrasts(design, tbl, f$rules), f$contrasts(design, tbl), info = omic)
        }
        for (kind in c("two_group", "three_group", "unbalanced", "batch_aware", "extra_metadata", "missing_group")) {
            payload <- f$payload(kind = kind, layout = "combined")
            for (require in c(TRUE, FALSE)) for (contrasts in list(payload$contrasts_tbl, NULL)) {
                expect_identical(
                    validateOmicsDesignDaPreflight(payload$design_matrix, payload$data_cln, f$rules,
                        contrastsTbl = contrasts, formulaString = "~ 0 + group", columnMapping = payload$column_mapping,
                        requireContrasts = require),
                    f$preflight(payload$design_matrix, payload$data_cln, contrastsTbl = contrasts,
                        formulaString = "~ 0 + group", columnMapping = payload$column_mapping, requireContrasts = require),
                    info = paste(omic, kind, require))
            }
        }
    }
    # The labels are the omic's.
    lipid <- validateOmicsDesignAlignment(design, list(), omicsDesignRules("lipidomics"))
    metab <- validateOmicsDesignAlignment(design, list(), omicsDesignRules("metabolomics"))
    expect_true("Lipidomics design requires a non-empty assay list." %in% lipid$errors)
    expect_true("Metabolomics design requires a non-empty assay list." %in% metab$errors)
})

test_that("C1 the save preflight is today's commit preflight (lipid, metab)", {
    for (omic in c("lipidomics", "metabolomics")) {
        fixture <- omicsDesignFixture(omic, layout = "combined")
        old <- if (identical(omic, "lipidomics")) validateLipidDesignDaPreflight else validateMetabDesignDaPreflight
        drifted <- fixture$payload
        drifted$data_cln[[1]]$WT_1 <- NULL
        no_formula <- fixture$payload
        no_formula$config_list$deAnalysisParameters$formula_string <- ""
        for (results in list(fixture$payload, drifted, no_formula)) {
            expect_identical(validateOmicsDesignSave(results, results$column_mapping, omic),
                old(designMatrix = results$design_matrix, assayList = results$data_cln,
                    contrastsTbl = results$contrasts_tbl, formulaString = "~ 0 + group",
                    columnMapping = results$column_mapping, requireContrasts = FALSE), info = omic)
        }
    }
})

# ------------------------------------------------------------------------------------------------
# C2: the formula.
# ------------------------------------------------------------------------------------------------

test_that("C2 the formula resolution is today's (lipid's function, metab's inline rule)", {
    skip("quarantined: Q-35 (KIT-103 W4) it compared the core with today's lipid design functions, deleted by W4 (and metab's, which W5 deletes); the comparison passed at fe1a9742 (kit-103 W2); the core is covered by the live blocks of this file, test-omic-design-module.R, MCI-013, MCI-019, MCI-027 and the entry lanes")
    metab_inline <- function(config_list) {
        formula_string <- config_list[["deAnalysisParameters"]][["formula_string"]]
        if (!metabDesignScalarString(formula_string)) formula_string <- "~ 0 + group"
        formula_string
    }
    configs <- list(NULL, list(), list(deAnalysisParameters = list()),
        list(deAnalysisParameters = list(formula_string = "")),
        list(deAnalysisParameters = list(formula_string = "  ")),
        list(deAnalysisParameters = list(formula_string = NA_character_)),
        list(deAnalysisParameters = list(formula_string = c("~ group", "~ batch"))),
        list(deAnalysisParameters = list(formula_string = 1)),
        list(deAnalysisParameters = list(formula_string = "~ group + batch")))
    for (i in seq_along(configs)) {
        expect_identical(resolveOmicsDesignFormula(configs[[i]]), lipidDesignResolveFormulaString(configs[[i]]),
            info = i)
        expect_identical(resolveOmicsDesignFormula(configs[[i]]), metab_inline(configs[[i]]), info = i)
    }
    expect_identical(resolveOmicsDesignFormula(list(), default = "~ batch"), "~ batch")
    expect_identical(resolveOmicsDesignFormula("not a list"), "~ 0 + group")
})

# ------------------------------------------------------------------------------------------------
# C3: the initial state and the run order.
# ------------------------------------------------------------------------------------------------

test_that("C3 the sample discovery and initial state are today's (lipid, trx, metab)", {
    skip("quarantined: Q-35 (KIT-103 W4) it compared the core with today's lipid design functions, deleted by W4 (and metab's, which W5 deletes); the comparison passed at fe1a9742 (kit-103 W2); the core is covered by the live blocks of this file, test-omic-design-module.R, MCI-013, MCI-019, MCI-027 and the entry lanes")
    samples <- c("WT_2", "KO_10", "WT_1", "KO_2")
    lipid_assays <- module_ci_lipid_design_assays("lc_pair", samples = samples)
    metab_assays <- module_ci_metab_design_assays("lc_pair", samples = samples)
    lipid_mapping <- module_ci_lipid_design_column_mapping(samples)
    metab_mapping <- module_ci_metab_design_column_mapping(samples)
    strip <- function(mapping) {
        mapping$sample_columns <- NULL
        mapping
    }
    existing <- module_ci_lipid_design_matrix("two_group", samples = c("WT_2", "WT_1", "KO_2", "KO_10"))
    contrasts <- data.frame(contrasts = c("groupKO-groupWT", "KO-WT"), stringsAsFactors = FALSE)
    config_list <- module_ci_lipid_design_config()
    unassigned <- existing
    unassigned$group <- NA_character_
    lipid_config <- omicsDesignTestConfig("lipidomics")
    trx_config <- omicsDesignTestConfig("transcriptomics")
    metab_config <- omicsDesignTestConfig("metabolomics")
    lipid_name_assays <- lipid_assays
    lipid_name_assays$LCMS_Pos$Lipid.name <- 1
    lipid_name_assays$LCMS_Pos$Metabolite.name <- 2
    for (mapping in list(lipid_mapping, strip(lipid_mapping), NULL)) for (assays in list(lipid_assays, lipid_name_assays)) {
        expect_identical(omicsDesignSampleColumns(assays, mapping, lipid_config),
            getLipidDesignSampleColumns(assays, mapping, messageFn = quietDesign))
        expect_identical(omicsDesignSampleColumns(assays, mapping, trx_config),
            getLipidDesignSampleColumns(assays, mapping, messageFn = quietDesign))
        for (existing_design in list(NULL, existing, function() existing, unassigned, function() stop("gone"))) {
            for (existing_contrasts in list(NULL, contrasts, function() contrasts)) {
                expect_identical(
                    buildOmicsDesignInitialState(assays, config_list, mapping, lipid_config, existing_design, existing_contrasts),
                    buildLipidDesignInitialState(assays, config_list, mapping, existing_design, existing_contrasts,
                        logErrorFn = quietDesign, messageFn = quietDesign))
                expect_identical(
                    buildOmicsDesignInitialState(assays, config_list, mapping, trx_config, existing_design, existing_contrasts),
                    buildLipidDesignInitialState(assays, config_list, mapping, existing_design, existing_contrasts,
                        mixedsortFn = trxOmicProfile()$run_order, logErrorFn = quietDesign, messageFn = quietDesign))
            }
        }
    }
    expect_null(buildOmicsDesignInitialState(list(), config_list, lipid_mapping, lipid_config))
    expect_null(buildOmicsDesignInitialState(list(A = "x"), config_list, lipid_mapping, lipid_config))
    # Metab: the builder server's closure (get_sample_columns, initial_state), without an owner.
    metab_name_assays <- metab_assays
    metab_name_assays$LCMS_Pos$Lipid.name <- 1
    metab_name_assays$LCMS_Pos$Metabolite.name <- 2
    for (mapping in list(metab_mapping, strip(metab_mapping), NULL)) for (assays in list(metab_assays, metab_name_assays)) {
        for (existing_design in list(NULL, existing, unassigned)) for (existing_contrasts in list(NULL, contrasts)) {
            today <- new.env()
            shiny::testServer(mod_metab_design_builder_server, args = list(data_tbl = function() assays,
                config_list = shiny::reactive(config_list), column_mapping = shiny::reactive(mapping),
                existing_design_matrix = if (!is.null(existing_design)) shiny::reactive(existing_design),
                existing_contrasts = if (!is.null(existing_contrasts)) shiny::reactive(existing_contrasts)), {
                today$columns <- suppressMessages(get_sample_columns(assays, mapping))
                today$state <- suppressMessages(shiny::isolate(initial_state()))
            })
            expect_identical(omicsDesignSampleColumns(assays, mapping, metab_config), today$columns)
            expect_identical(buildOmicsDesignInitialState(assays, config_list, mapping, metab_config, existing_design,
                existing_contrasts), today$state)
        }
    }
})

test_that("C3 the run order per omic: mixedsort for lipid and metab, input order for trx, file order for an imported design", {
    samples <- c("S10", "S2", "S1", "S3")
    assays <- module_ci_lipid_design_assays("lc_pair", samples = samples)
    mapping <- module_ci_lipid_design_column_mapping(samples)
    metab_assays <- module_ci_metab_design_assays("lc_pair", samples = samples)
    metab_mapping <- module_ci_metab_design_column_mapping(samples)
    config_list <- module_ci_lipid_design_config()
    fresh <- function(config, assays, mapping, existing = NULL) {
        buildOmicsDesignInitialState(assays, config_list, mapping, config, existing)$design_matrix$Run
    }
    expect_identical(fresh(omicsDesignTestConfig("lipidomics"), assays, mapping), c("S1", "S2", "S3", "S10"))
    expect_identical(fresh(omicsDesignTestConfig("metabolomics"), metab_assays, metab_mapping), c("S1", "S2", "S3", "S10"))
    expect_identical(fresh(omicsDesignTestConfig("transcriptomics"), assays, mapping), samples)
    # An imported design keeps its file order, for every omic.
    imported <- data.frame(Run = c("S3", "S10", "S1", "S2"), group = c("B", "B", "A", "A"), stringsAsFactors = FALSE)
    for (omic in c("lipidomics", "metabolomics", "transcriptomics")) {
        a <- if (identical(omic, "metabolomics")) metab_assays else assays
        m <- if (identical(omic, "metabolomics")) metab_mapping else mapping
        expect_identical(fresh(omicsDesignTestConfig(omic), a, m, imported), imported$Run, info = omic)
    }
    # The rule is the configuration's: none, none applied.
    config <- omicsDesignTestConfig("lipidomics")
    config$design$run_order <- NULL
    expect_error(buildOmicsDesignInitialState(assays, config_list, mapping, config), class = "omics_design_config_invalid")
    config$design$run_order <- rev
    expect_identical(fresh(config, assays, mapping), rev(samples))
})

test_that("C5 the saved assay columns follow the design's run order", {
    samples <- c("S1", "S2", "S3", "S10")
    assays <- module_ci_lipid_design_assays("lc_pair", samples = samples)
    mapping <- module_ci_lipid_design_column_mapping(samples)
    design <- data.frame(Run = c("S10", "S2", "S1", "S3"), group = c("B", "A", "A", "B"), factor1 = NA_character_,
        factor2 = NA_character_, factor3 = NA_character_, batch = NA_character_, replicates = c(1L, 1L, 2L, 2L),
        tech_reps = NA_integer_, stringsAsFactors = FALSE)
    payload <- buildOmicsDesignSavePayload(design, character(), assays, mapping, NULL, list(), "~ 0 + group",
        "lipid_id_col")
    for (assay in payload$data_cln) expect_identical(names(assay), c("LipidName", "LipidClass", design$Run))
    removed <- buildOmicsDesignSavePayload(design, "S2", assays, mapping, NULL, list(), "~ 0 + group", "lipid_id_col")
    expect_identical(names(removed$data_cln$LCMS_Pos), c("LipidName", "LipidClass", "S10", "S1", "S3"))
})

# ------------------------------------------------------------------------------------------------
# C4: the edits.
# ------------------------------------------------------------------------------------------------

test_that("C4 the design edits are today's lipid pure functions", {
    skip("quarantined: Q-35 (KIT-103 W4) it compared the core with today's lipid design functions, deleted by W4 (and metab's, which W5 deletes); the comparison passed at fe1a9742 (kit-103 W2); the core is covered by the live blocks of this file, test-omic-design-module.R, MCI-013, MCI-019, MCI-027 and the entry lanes")
    samples <- c("WT_1", "WT_2", "WT_3", "KO_1", "KO_2", "KO_3")
    assays <- module_ci_lipid_design_assays("lc_pair", samples = samples)
    design <- buildOmicsDesignInitialState(assays, list(), module_ci_lipid_design_column_mapping(samples),
        omicsDesignTestConfig("lipidomics"))$design_matrix
    rename <- c(KO_3 = "KO_3x", WT_1 = "WT_1x")
    expect_identical(renameOmicsDesignRuns(design, rename), renameLipidDesignDesignMatrixRuns(design, rename))
    expect_identical(renameOmicsDesignAssayColumns(assays, rename), renameLipidDesignAssayColumns(assays, rename))
    expect_identical(renameOmicsDesignTrackedSamples(samples, rename), renameLipidDesignTrackedSamples(samples, rename))
    expect_identical(applyOmicsDesignSampleRenameMap(design, assays, samples, rename),
        applyLipidDesignSampleRenameMap(design, assays, samples, rename))
    # The last range is out of bounds: NA names with extract_experiment's warning, in both.
    for (mode in list(list("before_underscore"), list("after_underscore"), list("range", 1, 2), list("range", 2, 2),
        list("range", 2, 4))) {
        expect_identical(suppressWarnings(buildOmicsDesignBulkRenameMap(samples, mode[[1]], mode[[2]] %||% NULL,
            mode[[3]] %||% NULL)), suppressWarnings(buildLipidDesignBulkRenameMap(samples, mode[[1]],
            mode[[2]] %||% NULL, mode[[3]] %||% NULL)))
    }
    expect_error(buildOmicsDesignBulkRenameMap(samples, "other"), "Unsupported transform mode.", fixed = TRUE)
    expect_error(buildOmicsDesignBulkRenameMap(samples, "range", NULL, 2), class = "shiny.silent.error")
    for (case in list(list(character(), "WT"), list("WT", " WT "), list("WT", ""), list("WT", "KO"))) {
        expect_identical(appendOmicsDesignFactorName(case[[1]], case[[2]]), appendLipidDesignFactorName(case[[1]], case[[2]]))
    }
    for (start in list(NULL, NA, 3L)) {
        expect_identical(buildOmicsDesignReplicateNumbers(samples[1:3], start),
            buildLipidDesignMetadataReplicateNumbers(samples[1:3], start))
    }
    for (parts in list(list("WT", "B1", "X"), list("WT", "", NULL), list("", NULL, NULL), list("WT", NULL, "X"))) {
        expect_identical(do.call(buildOmicsDesignGroupName, parts), do.call(buildLipidDesignMetadataGroupName, parts))
    }
    assigned <- design
    for (step in list(list(samples[1:3], "WT", "", "", 1L), list(samples[4:6], "KO", "B1", "X", NA),
        list(samples[4:5], "KO", "", NULL, 2L))) {
        next_design <- applyOmicsDesignMetadataAssignment(assigned, step[[1]], step[[2]], step[[3]], step[[4]], step[[5]])
        expect_identical(next_design, applyLipidDesignMetadataAssignment(assigned, step[[1]], step[[2]], step[[3]], step[[4]],
            step[[5]]))
        assigned <- next_design
    }
    expect_identical(listOmicsDesignAssignedGroups(assigned), listLipidDesignAssignedGroups(assigned))
    # A NULL second factor fails in both; the builder never sends one (its input is "").
    expect_error(applyLipidDesignMetadataAssignment(assigned, "WT_1", "WT", NULL), "replacement has length zero")
    expect_error(applyOmicsDesignMetadataAssignment(assigned, "WT_1", "WT", NULL), "replacement has length zero")
    assigned <- applyOmicsDesignMetadataAssignment(assigned, samples[4:6], "KO", "", "", 1L)
    for (case in list(list("WT_1", "lowest"), list(c("WT_1", "KO_1"), "lowest"), list(c("WT_2", "WT_3"), "lowest"),
        list(c("WT_3", "WT_1"), "first"), list(c("KO_1", "KO_2"), "manual", 5L))) {
        expect_identical(applyOmicsDesignTechRepAssignment(assigned, case[[1]], case[[2]], case[[3]] %||% NULL),
            applyLipidDesignTechRepAssignment(assigned, case[[1]], case[[2]], case[[3]] %||% NULL))
    }
    expect_error(applyOmicsDesignTechRepAssignment(assigned, c("WT_2", "WT_3"), "other"),
        "Unsupported technical replicate assignment mode.", fixed = TRUE)
    contrasts <- buildOmicsDesignContrastState(NULL)
    for (pair in list(c("KO", "WT"), c("KO", "WT"), c("WT", "WT"), c("", "WT"), c("WT", "KO"))) {
        expect_identical(appendOmicsDesignContrast(contrasts, pair[[1]], pair[[2]]),
            appendLipidDesignContrast(contrasts, pair[[1]], pair[[2]]))
        contrasts <- appendOmicsDesignContrast(contrasts, pair[[1]], pair[[2]])
    }
    expect_identical(appendOmicsDesignContrast(contrasts, NULL, "WT"), appendLipidDesignContrast(contrasts, NULL, "WT"))
    expect_identical(appendOmicsDesignRemovedSamples("A", c("B", "A")), appendLipidDesignRemovedSamples("A", c("B", "A")))
    existing <- data.frame(contrasts = c("groupKO-groupWT", "KO-WT", "groupA.1-groupB"), stringsAsFactors = FALSE)
    for (value in list(NULL, existing, existing[0, , drop = FALSE], "x")) {
        expect_identical(buildOmicsDesignContrastState(value), buildLipidDesignContrastState(value))
    }
})

# ------------------------------------------------------------------------------------------------
# C5: the save payload, and a whole builder session against today's builder of each omic.
# ------------------------------------------------------------------------------------------------

test_that("C5 the save payload is today's lipid payload", {
    skip("quarantined: Q-35 (KIT-103 W4) it compared the core with today's lipid design functions, deleted by W4 (and metab's, which W5 deletes); the comparison passed at fe1a9742 (kit-103 W2); the core is covered by the live blocks of this file, test-omic-design-module.R, MCI-013, MCI-019, MCI-027 and the entry lanes")
    samples <- c("WT_1", "WT_2", "KO_1", "KO_2")
    assays <- module_ci_lipid_design_assays("lc_pair", samples = samples)
    design <- module_ci_lipid_design_matrix("two_group")
    contrasts <- appendOmicsDesignContrast(buildOmicsDesignContrastState(NULL), "KO", "WT")
    contrasts <- appendOmicsDesignContrast(contrasts, "WT.1", "KO")
    mappings <- list(module_ci_lipid_design_column_mapping(samples),
        list(lipid_id_col = "LipidName", annotation_col = NA_character_),
        list(lipid_id_col = "LipidName", annotation_col = ""), list(lipid_id_col = "Absent"), NULL)
    unassigned <- design
    unassigned$group <- NA_character_
    for (mapping in mappings) for (removed in list(character(), "KO_1")) for (formula in c("~ 0 + group", "~ group")) {
        for (contrast_data in list(contrasts, contrasts[0, ], NULL)) for (d in list(design, unassigned)) {
            expect_identical(
                buildOmicsDesignSavePayload(d, removed, assays, mapping, contrast_data, module_ci_lipid_design_config(),
                    formula, "lipid_id_col"),
                buildLipidDesignSaveResultsPayload(d, removed, assays, mapping, contrast_data,
                    module_ci_lipid_design_config(), formula))
        }
    }
    # The contract's contrast format.
    table <- buildOmicsDesignSaveContrastsTable(contrasts, "~ 0 + group")
    expect_identical(table$contrasts, c("groupKO-groupWT", "groupWT.1-groupKO"))
    expect_identical(table$friendly_names, c("KO_vs_WT", "WT_1_vs_KO"))
    expect_identical(buildOmicsDesignSaveContrastsTable(contrasts, "~ group")$contrasts, c("KO-WT", "WT.1-KO"))
})

# One builder session on today's builder server: factors, assignment, a rename, technical
# replicates, a removal, a contrast, then save.
omicsDesignBuilderSession <- function(server, assays, mapping, extra_args = list()) {
    out <- new.env()
    shiny::testServer(server, args = c(list(data_tbl = function() assays,
        config_list = shiny::reactive(list(deAnalysisParameters = list(formula_string = "~ 0 + group"))),
        column_mapping = shiny::reactive(mapping)), extra_args), {
        session$setInputs(new_factor = "WT")
        session$setInputs(add_factor = 1)
        session$setInputs(new_factor = "KO")
        session$setInputs(add_factor = 2)
        session$setInputs(sample_to_rename = "KO_3", new_sample_name = "KO_3x")
        session$setInputs(rename_sample = 1)
        session$setInputs(selected_runs = c("WT_1", "WT_2", "WT_3"), factor1_select = "WT", factor2_select = "",
            factor3_select = "", replicate_start = 1L)
        session$setInputs(assign_metadata = 1)
        session$setInputs(selected_runs = c("KO_1", "KO_2", "KO_3x"), factor1_select = "KO", replicate_start = 1L)
        session$setInputs(assign_metadata = 2)
        session$setInputs(tech_rep_samples = c("WT_2", "WT_3"), tech_rep_assignment_mode = "lowest")
        session$setInputs(assign_tech_reps = 1)
        session$setInputs(samples_to_remove = "KO_1")
        session$setInputs(remove_samples = 1)
        session$setInputs(contrast_group1 = "KO", contrast_group2 = "WT", formula_string = "~ 0 + group")
        session$setInputs(add_contrast = 1)
        session$setInputs(save_results = 1)
        out$result <- shiny::isolate(session$returned())
        out$factors <- shiny::isolate(factors())
    })
    out
}

# The same session through the core.
omicsDesignCoreSession <- function(assays, mapping, config) {
    state <- buildOmicsDesignInitialState(assays, list(deAnalysisParameters = list(formula_string = "~ 0 + group")),
        mapping, config)
    factors <- appendOmicsDesignFactorName(appendOmicsDesignFactorName(state$factors, "WT"), "KO")
    renamed <- applyOmicsDesignSampleRenameMap(state$design_matrix, state$data_cln, state$sample_names,
        c(KO_3 = "KO_3x"))
    design <- applyOmicsDesignMetadataAssignment(renamed$designMatrix, c("WT_1", "WT_2", "WT_3"), "WT", "", "", 1L)
    design <- applyOmicsDesignMetadataAssignment(design, c("KO_1", "KO_2", "KO_3x"), "KO", "", "", 1L)
    design <- applyOmicsDesignTechRepAssignment(design, c("WT_2", "WT_3"), "lowest")$designMatrix
    contrasts <- appendOmicsDesignContrast(state$contrasts, "KO", "WT")
    removed <- appendOmicsDesignRemovedSamples(character(), "KO_1")
    payload <- buildOmicsDesignSavePayload(design, removed, renamed$assayList, mapping, contrasts,
        list(deAnalysisParameters = list(formula_string = "~ 0 + group")), "~ 0 + group", config$design$id_mapping_key)
    list(payload = payload, factors = factors, removed = removed)
}

test_that("C3 to C5 a builder session through the core gives today's result (lipid, metab, trx)", {
    samples <- c("WT_2", "KO_1", "WT_1", "KO_3", "WT_3", "KO_2")
    lipid_assays <- module_ci_lipid_design_assays("lc_pair", samples = samples)
    lipid_mapping <- module_ci_lipid_design_column_mapping(samples)
    today <- omicsDesignBuilderSession(mod_lipid_design_builder_server, lipid_assays, lipid_mapping)
    core <- omicsDesignCoreSession(lipid_assays, lipid_mapping, omicsDesignTestConfig("lipidomics"))
    expect_identical(core$payload, today$result)
    expect_identical(core$payload$design_matrix$Run, c("KO_2", "KO_3x", "WT_1", "WT_2", "WT_3"))

    metab_assays <- module_ci_metab_design_assays("lc_pair", samples = samples)
    metab_mapping <- module_ci_metab_design_column_mapping(samples)
    today <- suppressMessages(omicsDesignBuilderSession(mod_metab_design_builder_server, metab_assays, metab_mapping))
    core <- omicsDesignCoreSession(metab_assays, metab_mapping, omicsDesignTestConfig("metabolomics"))
    expect_identical(core$payload, today$result)

    trx <- omicsDesignTestConfig("transcriptomics")
    today <- omicsDesignBuilderSession(mod_lipid_design_builder_server, lipid_assays, lipid_mapping,
        list(profile = trxOmicProfile()))
    core <- omicsDesignCoreSession(lipid_assays, lipid_mapping, trx)
    # The builder adds the level order and the removed runs (`record_excluded_runs`).
    expect_true(trx$design$record_excluded_runs)
    core$payload$group_levels <- trx$design$group_levels(core$payload$design_matrix, core$factors)
    core$payload$excluded_runs <- core$removed
    expect_identical(core$payload, today$result)
    expect_identical(core$payload$design_matrix$Run, c("WT_2", "WT_1", "KO_3x", "WT_3", "KO_2"))
})

# ------------------------------------------------------------------------------------------------
# C7: the objects.
# ------------------------------------------------------------------------------------------------

test_that("C7 the lipid and metab objects hold today's argument mapping", {
    lipid <- module_ci_lipid_design_payload(layout = "lc_pair")
    lipid$design_matrix$tech_rep_group <- paste(lipid$design_matrix$group, lipid$design_matrix$replicates, sep = "_")
    metab <- module_ci_metab_design_payload(layout = "lc_pair")
    metab$design_matrix$tech_rep_group <- paste(metab$design_matrix$group, metab$design_matrix$replicates, sep = "_")
    pick <- function(value) if (!is.null(value) && !is.na(value) && nzchar(value)) value else NA_character_
    pattern <- function(value) if (!is.null(value) && !is.na(value)) value else NA_character_
    variants <- list(identity, function(m) {
        m$annotation_col <- NA_character_
        m
    }, function(m) {
        m$annotation_col <- NULL
        m$is_pattern <- "^LC"
        m
    }, function(m) {
        m$annotation_col <- ""
        m$is_pattern <- NULL
        m
    })
    for (variant in variants) {
        mapping <- variant(lipid$column_mapping)
        expect_identical(lipidDesignObject(lipid, mapping), createLipidomicsAssayData(lipid_data = lipid$data_cln,
            design_matrix = lipid$design_matrix, lipid_id_column = mapping$lipid_id_col,
            annotation_id_column = pick(mapping$annotation_col), sample_id = "Run", group_id = "group",
            technical_replicate_id = "tech_rep_group", database_identifier_type = "Unknown",
            internal_standard_regex = pattern(mapping$is_pattern), args = lipid$config_list))
        mapping <- variant(metab$column_mapping)
        expect_identical(metabDesignObject(metab, mapping), createMetaboliteAssayData(metabolite_data = metab$data_cln,
            design_matrix = metab$design_matrix, metabolite_id_column = mapping$metabolite_id_col,
            annotation_id_column = pick(mapping$annotation_col), sample_id = "Run", group_id = "group",
            technical_replicate_id = "tech_rep_group", database_identifier_type = "Unknown",
            internal_standard_regex = pattern(mapping$is_pattern), args = metab$config_list))
    }
})

# ------------------------------------------------------------------------------------------------
# C8: the commit.
# ------------------------------------------------------------------------------------------------

test_that("C8 the commit is today's for lipid and metab: bus, artifact bytes, state and provenance", {
    skip("quarantined: Q-37 (KIT-103 W5) it compared the core commit with today's metab builder-results observer (registerMetabDesignBuilderResultsObserver), deleted by W5, and today's lipid commit, now the shared commit itself (W4); the comparison passed at fe1a9742 (kit-103 W2); the commit is covered by MCI-013.5, MCI-019.5, MCI-027, test-omic-design-module.R and the entry lanes' design step")
    for (omic in c("lipidomics", "metabolomics")) {
        old <- omicsDesignFixture(omic)
        new <- omicsDesignFixture(omic)
        seen <- omicsDesignCaptureShiny()
        shiny::withReactiveDomain(shiny::MockShinySession$new(), {
            if (identical(omic, "lipidomics")) {
                runLipidDesignBuilderObserverShell(results = old$payload, workflowData = old$workflow,
                    experimentPaths = old$paths, qcTrigger = function(value) NULL)
            } else {
                runTodayMetabDesignCommit(old$payload, old$workflow, old$paths)
            }
            committed <- commitOmicsDesign(new$workflow, new$payload, new$paths$source_dir, omicsDesignTestConfig(omic))
        })
        expect_identical(seen$notes, "Design matrix and contrasts saved successfully!", info = omic)
        expect_s4_class(committed$object, if (identical(omic, "lipidomics")) "LipidomicsAssayData" else "MetaboliteAssayData")
        old_bus <- omicsDesignBus(old$workflow)
        new_bus <- omicsDesignBus(new$workflow)
        # The GUI sets the tab status after the commit (KIT-102's split).
        expect_identical(old_bus$tab_status$design_matrix, "complete")
        expect_identical(new_bus$tab_status$design_matrix, "pending")
        old_bus$tab_status <- new_bus$tab_status <- NULL
        expect_identical(new_bus, old_bus, info = omic)
        expect_identical(omicsDesignArtifactBytes(new$paths$source_dir), omicsDesignArtifactBytes(old$paths$source_dir),
            info = omic)
        expect_length(omicsDesignArtifactBytes(new$paths$source_dir), 8L)
        expect_identical(omicsDesignOwnerRecord(new$workflow$state_manager),
            omicsDesignOwnerRecord(old$workflow$state_manager), info = omic)
    }
})

test_that("C8 the commit is today's for trx (MCI-027's corpus through the builder)", {
    skip_if_not(workflowOwnerCapable())
    design <- module_ci_trx_corpus_design()
    old <- module_ci_trx_fixture_owner()
    new <- module_ci_trx_fixture_owner()
    for (fixture in list(old, new)) module_ci_trx_import(fixture, file.path(module_ci_trx_corpus(), "counts.tsv"))
    old_result <- module_ci_trx_build_design(old, design, c("RPMI", "Sera"))
    new_result <- module_ci_trx_build_design(new, design, c("RPMI", "Sera"))
    module_ci_trx_commit(old, old_result)
    shiny::withReactiveDomain(new$session, commitOmicsDesign(new$workflow, new_result, new$paths$source_dir,
        omicsDesignTestConfig("transcriptomics")))
    old_bus <- omicsDesignBus(old$workflow)
    new_bus <- omicsDesignBus(new$workflow)
    old_bus$tab_status <- new_bus$tab_status <- NULL
    expect_identical(omicsDesignRelative(new_bus, new$paths$base_dir), omicsDesignRelative(old_bus, old$paths$base_dir))
    expect_identical(omicsDesignArtifactBytes(new$paths$source_dir), omicsDesignArtifactBytes(old$paths$source_dir))
    expect_identical(omicsDesignOwnerRecord(new$workflow$state_manager), omicsDesignOwnerRecord(old$workflow$state_manager))
    expect_identical(new$workflow$state_manager$getHistory()[[2L]], "trx_raw_data_s4")
})

test_that("C8 refusals come before the bus and the files; a later failure carries its message", {
    skip("quarantined: Q-37 (KIT-103 W5) it took today's refusal text from metab's builder-results observer (registerMetabDesignBuilderResultsObserver), deleted by W5; the comparison passed at fe1a9742 (kit-103 W2); the refusals are covered by KIT-101 F-01's two tests and the W4 commit blocks in test-omic-design-module.R, and MCI-019.6")
    for (omic in c("lipidomics", "metabolomics")) {
        config <- omicsDesignTestConfig(omic)
        # A missing source directory, as today's text.
        fixture <- omicsDesignFixture(omic)
        before <- omicsDesignBus(fixture$workflow)
        error <- tryCatch(commitOmicsDesign(fixture$workflow, fixture$payload, file.path(fixture$paths$base_dir, "absent"),
            config), error = identity)
        expect_s3_class(error, "omics_design_refused")
        expect_identical(conditionMessage(error), "Could not find source directory to save files.")
        expect_identical(omicsDesignBus(fixture$workflow), before)
        error <- tryCatch(commitOmicsDesign(fixture$workflow, fixture$payload, NULL, config), error = identity)
        expect_s3_class(error, "omics_design_refused")
        # A failed preflight: today's notification text, nothing written, no state.
        drifted <- omicsDesignFixture(omic)
        drifted$payload$data_cln[[1]]$WT_2 <- NULL
        today <- omicsDesignFixture(omic)
        today$payload$data_cln[[1]]$WT_2 <- NULL
        seen <- omicsDesignCaptureShiny()
        shiny::withReactiveDomain(shiny::MockShinySession$new(), {
            if (identical(omic, "lipidomics")) {
                runLipidDesignBuilderObserverShell(results = today$payload, workflowData = today$workflow,
                    experimentPaths = today$paths)
            } else {
                runTodayMetabDesignCommit(today$payload, today$workflow, today$paths)
            }
        })
        before <- omicsDesignBus(drifted$workflow)
        error <- tryCatch(commitOmicsDesign(drifted$workflow, drifted$payload, drifted$paths$source_dir, config),
            error = identity)
        expect_s3_class(error, "omics_design_refused")
        expect_identical(conditionMessage(error), seen$notes, info = omic)
        expect_match(conditionMessage(error), paste0("^Invalid ", config$design$label_lower, " design: "))
        expect_identical(omicsDesignBus(drifted$workflow), before)
        expect_length(list.files(drifted$paths$source_dir), 0L)
        expect_identical(drifted$workflow$state_manager$getHistory(), "initial")
        # A failure after the bus write: the original message, classed.
        failing <- omicsDesignFixture(omic)
        config$design$build_object <- function(results, columnMapping) stop("object refused")
        error <- tryCatch(shiny::withReactiveDomain(shiny::MockShinySession$new(),
            commitOmicsDesign(failing$workflow, failing$payload, failing$paths$source_dir, config)), error = identity)
        expect_s3_class(error, "omics_design_save_failed")
        expect_identical(conditionMessage(error), "object refused")
        expect_identical(failing$workflow$design_matrix, failing$payload$design_matrix)
        expect_true(file.exists(file.path(failing$paths$source_dir, "design_matrix.tab")))
        expect_identical(failing$workflow$state_manager$getHistory(), "initial")
    }
    # Trx's refusal text uses its own word.
    trx <- omicsDesignTestConfig("transcriptomics")
    trx$design$validate <- function(results, columnMapping) list(valid = FALSE, errors = c("a", "b"))
    fixture <- omicsDesignFixture("lipidomics")
    error <- tryCatch(commitOmicsDesign(fixture$workflow, fixture$payload, fixture$paths$source_dir, trx), error = identity)
    expect_identical(conditionMessage(error), "Invalid transcriptomics design: a; b")
})

test_that("C8 the artifacts are today's seven files, without contrasts or a column mapping", {
    results <- module_ci_lipid_design_payload(layout = "combined")
    dir <- withr::local_tempdir()
    paths <- writeOmicsDesignArtifacts(dir, results, results$column_mapping, results$config_list)
    expect_identical(sort(basename(unlist(paths, use.names = FALSE))), sort(c("design_matrix.tab", "contrast_strings.tab",
        "data_cln_LCMS_Pos.tab", "data_cln_GCMS.tab", "assay_manifest.txt", "column_mapping.json", "manifest.json",
        "config.ini")))
    expect_identical(readLines(file.path(dir, "contrast_strings.tab")), results$contrasts_tbl$contrasts)
    expect_identical(jsonlite::read_json(file.path(dir, "manifest.json")), list(data_path = "assay_manifest.txt",
        design_matrix_path = "design_matrix.tab", contrast_strings_path = "contrast_strings.tab",
        column_mapping_path = "column_mapping.json"))
    bare <- withr::local_tempdir()
    results$contrasts_tbl <- results$contrasts_tbl[0, ]
    writeOmicsDesignArtifacts(bare, results, NULL, NULL)
    expect_setequal(list.files(bare), c("design_matrix.tab", "data_cln_LCMS_Pos.tab", "data_cln_GCMS.tab",
        "assay_manifest.txt", "manifest.json"))
})

# ------------------------------------------------------------------------------------------------
# C10: the design import.
# ------------------------------------------------------------------------------------------------

# Today's import of one omic onto `workflow` (lipid's confirmation shell; metab's observer chain).
runTodayDesignImport <- function(omic, workflow, paths, import_root, baseline) {
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    notify <- function(ui, ..., type = NULL) seen$notes <- c(seen$notes, paste(as.character(ui), collapse = " "))
    seen$qc <- FALSE
    shiny::withReactiveDomain(omicsOwnerSession(), if (identical(omic, "lipidomics")) {
        runLipidDesignImportConfirmationShell(workflowData = workflow, experimentPaths = paths, importPath = import_root,
            qcTrigger = function(value) seen$qc <- value, removeModalFn = quietDesign, removeNotificationFn = quietDesign,
            showNotificationFn = notify, updateLipidFilteringFn = baseline)
    } else {
        registerMetabDesignImportObserverShell(input = list(import_dir = import_root), resolvedVolumes = NULL,
            workflowData = workflow, experimentPaths = paths, qcTrigger = function(value) seen$qc <- value,
            observeEventFn = function(eventExpr, handlerExpr, ...) handlerExpr,
            resolveImportPreflight = function(input, resolvedVolumes) {
                resolveMetabDesignImportPreflight(input, resolvedVolumes, parseDirPathFn = function(v, d) d)
            },
            initializeImportedQcBaseline = function(s4Object) {
                initializeMetabDesignImportedQcBaseline(s4Object, updateFilteringFn = baseline)
            },
            completeImportedPostCheckpoint = function(...) {
                completeMetabDesignImportedPostCheckpoint(..., removeNotification = quietDesign,
                    showNotification = notify)
            },
            removeModal = quietDesign, showNotification = notify, removeNotification = quietDesign)
    })
    seen
}

# The core's import, as the shared frame will run it: the file check, the reader, the commit.
runCoreDesignImport <- function(omic, workflow, paths, import_root, baseline) {
    config <- omicsDesignTestConfig(omic)
    config$design$design_import$qc_baseline <- baseline
    checked <- checkOmicsDesignImportFiles(import_root)
    if (!isTRUE(checked$ok)) return(list(error = checked$message))
    tryCatch(shiny::withReactiveDomain(omicsOwnerSession(), {
        imported <- readOmicsDesignImport(import_root, config, paths$source_dir)
        commitOmicsDesignImport(workflow, imported, paths$source_dir, config)
    }), error = function(error) list(error = paste("Error during import:", conditionMessage(error)), condition = error))
}

# Today's imports hold vroom's tables, which the owner refuses (W2 finding; the reader test below):
# for the comparisons, today's and the core's reads are plain `read.delim` tables, as MCI-019.4.
omicsDesignPlainReads <- function(env = parent.frame()) {
    testthat::local_mocked_bindings(vroom = function(file, show_col_types = FALSE, ...) {
        utils::read.delim(file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
    }, .package = "vroom", .env = env)
}

omicsDesignImportCase <- function(omic, edit = identity) {
    fixture <- omicsDesignFixture(omic, layout = "combined")
    root <- tempfile("omics-design-import-")
    if (identical(omic, "lipidomics")) {
        module_ci_lipid_design_write_import_pack(root, fixture$payload)
    } else {
        module_ci_metab_design_write_import_pack(root, fixture$payload)
    }
    edit(root)
    root
}

test_that("C10 the design import is today's for lipid and metab: bus, state, provenance, baseline and texts", {
    skip("quarantined: Q-35 (KIT-103 W4) it compared the core with today's lipid design functions, deleted by W4 (and metab's, which W5 deletes); the comparison passed at fe1a9742 (kit-103 W2); the core is covered by the live blocks of this file, test-omic-design-module.R, MCI-013, MCI-019, MCI-027 and the entry lanes")
    omicsDesignPlainReads()
    edits <- list(
        full = identity,
        inferred_mapping = function(root) unlink(file.path(root, "column_mapping.json")),
        no_contrasts = function(root) unlink(file.path(root, "contrast_strings.tab")),
        empty_contrasts = function(root) writeLines(character(), file.path(root, "contrast_strings.tab")),
        no_config = function(root) unlink(file.path(root, "config.ini")),
        reordered_design = function(root) {
            design <- utils::read.delim(file.path(root, "design_matrix.tab"), check.names = FALSE)
            utils::write.table(design[rev(seq_len(nrow(design))), ], file.path(root, "design_matrix.tab"), sep = "\t",
                quote = FALSE, row.names = FALSE)
        })
    for (omic in c("lipidomics", "metabolomics")) for (case in names(edits)) {
        root <- omicsDesignImportCase(omic, edits[[case]])
        old <- omicsDesignFixture(omic)
        new <- omicsDesignFixture(omic)
        calls <- new.env()
        calls$old <- calls$new <- list()
        today <- runTodayDesignImport(omic, old$workflow, old$paths, root, function(...) calls$old <- c(calls$old, list(list(...))))
        core <- runCoreDesignImport(omic, new$workflow, new$paths, root, function(...) calls$new <- c(calls$new, list(list(...))))
        info <- paste(omic, case)
        expect_null(core$error, info = info)
        expect_identical(sprintf("Design imported successfully! Loaded %d assays with %d samples.", core$n_assays,
            core$n_samples), today$notes[[length(today$notes)]], info = info)
        expect_true(today$qc, info = info)
        old_bus <- omicsDesignBus(old$workflow)
        new_bus <- omicsDesignBus(new$workflow)
        expect_identical(old_bus$tab_status$design_matrix, "complete", info = info)
        old_bus$tab_status <- new_bus$tab_status <- NULL
        expect_identical(new_bus, old_bus, info = info)
        expect_identical(omicsDesignOwnerRecord(new$workflow$state_manager),
            omicsDesignOwnerRecord(old$workflow$state_manager), info = info)
        expect_identical(calls$new, calls$old, info = info)
        expect_length(calls$new, 1L)
        # The imported design keeps its file order.
        expect_identical(new$workflow$design_matrix$Run,
            utils::read.delim(file.path(root, "design_matrix.tab"), check.names = FALSE)$Run, info = info)
    }
})

test_that("C10 the import's refusals are today's texts and write nothing", {
    skip("quarantined: Q-35 (KIT-103 W4) it compared the core with today's lipid design functions, deleted by W4 (and metab's, which W5 deletes); the comparison passed at fe1a9742 (kit-103 W2); the core is covered by the live blocks of this file, test-omic-design-module.R, MCI-013, MCI-019, MCI-027 and the entry lanes")
    omicsDesignPlainReads()
    edits <- list(
        no_design = function(root) unlink(file.path(root, "design_matrix.tab")),
        no_manifest = function(root) unlink(file.path(root, "assay_manifest.txt")),
        empty_manifest = function(root) writeLines(c("", ""), file.path(root, "assay_manifest.txt")),
        missing_assay = function(root) unlink(file.path(root, "data_cln_GCMS.tab")),
        drifted = function(root) {
            design <- utils::read.delim(file.path(root, "design_matrix.tab"), check.names = FALSE)
            design$Run[[1]] <- "ZZ_9"
            utils::write.table(design, file.path(root, "design_matrix.tab"), sep = "\t", quote = FALSE, row.names = FALSE)
        },
        bad_contrast = function(root) writeLines("groupMISSING-groupWT", file.path(root, "contrast_strings.tab")))
    for (omic in c("lipidomics", "metabolomics")) for (case in names(edits)) {
        root <- omicsDesignImportCase(omic, edits[[case]])
        old <- omicsDesignFixture(omic)
        new <- omicsDesignFixture(omic)
        today <- runTodayDesignImport(omic, old$workflow, old$paths, root, quietDesign)
        before <- omicsDesignBus(new$workflow)
        core <- runCoreDesignImport(omic, new$workflow, new$paths, root, quietDesign)
        info <- paste(omic, case)
        if (identical(omic, "metabolomics") && case %in% c("drifted", "bad_contrast")) {
            # The metab import runs no DA preflight (design_import$preflight FALSE): an unknown contrast term is
            # imported, and a drifted design fails later, at the object, after the bus is written, in both.
            expect_false(omicsDesignTestConfig(omic)$design$design_import$preflight)
            if (identical(case, "bad_contrast")) {
                expect_null(core$error, info = info)
            } else {
                expect_identical(core$error, today$notes[[length(today$notes)]], info = info)
                expect_match(core$error, "invalid class", fixed = TRUE)
                expect_identical(omicsDesignBus(new$workflow), omicsDesignBus(old$workflow), info = info)
            }
            next
        }
        expect_identical(core$error, today$notes[[length(today$notes)]], info = info)
        if (case %in% c("drifted", "bad_contrast")) expect_s3_class(core$condition, "omics_design_refused")
        if (identical(case, "drifted")) expect_match(core$error, "Invalid lipidomics design import: ", fixed = TRUE)
        expect_identical(omicsDesignBus(new$workflow), before, info = info)
        expect_identical(new$workflow$state_manager$getHistory(), "initial", info = info)
    }
})

test_that("C10 the reader infers the column mapping as today, by the omic's id key and candidates", {
    for (omic in c("lipidomics", "metabolomics")) {
        root <- omicsDesignImportCase(omic, function(root) unlink(file.path(root, "column_mapping.json")))
        config <- omicsDesignTestConfig(omic)
        imported <- readOmicsDesignImport(root, config, NULL)
        key <- config$design$id_mapping_key
        expect_identical(names(imported$column_mapping), c(key, "annotation_col", "sample_columns", "is_pattern"))
        expect_identical(imported$column_mapping$sample_columns, imported$design_matrix$Run)
        expect_identical(imported$design_file, file.path(root, "design_matrix.tab"))
        expect_identical(imported$assay_names, c("LCMS_Pos", "GCMS"))
    }
    # A Peak ID column is the id; no annotation candidate leaves the annotation NULL.
    root <- omicsDesignImportCase("lipidomics", function(root) {
        unlink(file.path(root, "column_mapping.json"))
        for (assay in c("LCMS_Pos", "GCMS")) {
            path <- file.path(root, paste0("data_cln_", assay, ".tab"))
            table <- utils::read.delim(path, check.names = FALSE)
            names(table)[names(table) == "LipidClass"] <- "Peak ID"
            utils::write.table(table, path, sep = "\t", quote = FALSE, row.names = FALSE)
        }
    })
    mapping <- readOmicsDesignImport(root, omicsDesignTestConfig("lipidomics"), NULL)$column_mapping
    expect_identical(mapping$lipid_id_col, "Peak ID")
    expect_true("annotation_col" %in% names(mapping))
    expect_null(mapping$annotation_col)
    # The contrasts' friendly names as today's import (P-09 keeps the two rules).
    root <- omicsDesignImportCase("lipidomics", function(root) {
        writeLines(c("groupKO-groupWT", "groupA.1-groupB"), file.path(root, "contrast_strings.tab"))
    })
    contrasts <- readOmicsDesignImport(root, omicsDesignTestConfig("lipidomics"), NULL)$contrasts_tbl
    expect_identical(contrasts$friendly_names, c("KO_vs_WT", "A.1_vs_B"))
    expect_identical(contrasts$full_format, c("KO_vs_WT=groupKO-groupWT", "A.1_vs_B=groupA.1-groupB"))
})

test_that("W2 finding: the design import holds vroom's tables (today's fails: 'payload contains a nonportable resource')", {
    skip("quarantined: Q-35 (KIT-103 W4) it compared the core with today's lipid design functions, deleted by W4 (and metab's, which W5 deletes); the comparison passed at fe1a9742 (kit-103 W2); the core is covered by the live blocks of this file, test-omic-design-module.R, MCI-013, MCI-019, MCI-027 and the entry lanes")
    for (omic in c("lipidomics", "metabolomics")) {
        root <- omicsDesignImportCase(omic)
        # Today: the owner refuses the tables' live parser pointer when they are held (KIT-075).
        old <- omicsDesignFixture(omic)
        today <- runTodayDesignImport(omic, old$workflow, old$paths, root, quietDesign)
        expect_identical(today$notes[[length(today$notes)]], "Error during import: payload contains a nonportable resource",
            info = omic)
        expect_identical(old$workflow$state_manager$getHistory(), "initial")
        # The core: the same read, held, committed and the same values as a plain read.
        fixture <- omicsDesignFixture(omic)
        core <- runCoreDesignImport(omic, fixture$workflow, fixture$paths, root, quietDesign)
        expect_null(core$error, info = omic)
        expect_identical(fixture$workflow$state_manager$getHistory()[[2L]], omicGuiConfig(omic)$design$state_name)
        imported <- readOmicsDesignImport(root, omicsDesignTestConfig(omic), NULL)
        for (table in c(list(imported$design_matrix), imported$assay_list)) {
            expect_s3_class(table, "spec_tbl_df")
            expect_null(attr(table, "problems", exact = TRUE))
        }
        held <- getWorkflowOmicsBusTables(fixture$workflow, "data_cln")
        for (assay in imported$assay_names) {
            plain <- utils::read.delim(file.path(root, paste0("data_cln_", assay, ".tab")), check.names = FALSE,
                stringsAsFactors = FALSE)
            expect_equal(as.data.frame(held[[assay]]), plain, ignore_attr = TRUE, info = paste(omic, assay))
        }
    }
})
