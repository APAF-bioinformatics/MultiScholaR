# KIT-224 (decision M8, confirmed by the owner 2026-11-25): the legacy lipidomics and metabolomics RUV-III-C
# correction (`ruvIII_C_Varying`, the shipped Automatic default and Manual) keeps each feature's level, as proteomics
# does since KIT-123 R2; the transcriptomics RUVIII export file and integration input keep each gene's TMM logCPM
# level (OD-2, T-a'). Analyses recorded before replay with the shift (replay rule marker `ruv_feature_levels`).
# E1 to E6 and G1, G2 of the ticket's section 10. The option and its tokens are written as text so that the file runs
# on the commit before the change (where the marker does not exist).

kit224Quiet <- function(expr) {
    value <- NULL
    utils::capture.output(value <- suppressWarnings(suppressMessages(expr)))
    value
}

KIT224_OPTION <- "MultiScholaR.ruv_feature_levels"

# KIT-185's lipid shape (ticket section 4): two assays, two groups of six runs, levels about 20 log2, one sample-wise
# unwanted factor uncorrelated with the groups. `missing`: the share of cells removed, value-dependent (drawn among the
# lowest 20 % of the values, at most two per feature and group).
kit224Mats <- function(seed = 224000001L, missing = 0, p = 60L) {
    set.seed(seed)
    design <- data.frame(Sample_ID = sprintf("S%02d", 1:12), group = rep(c("A", "B"), each = 6L), stringsAsFactors = FALSE)
    f <- stats::residuals(stats::lm(stats::rnorm(12L) ~ factor(design$group)))
    f <- f / sqrt(mean(f^2))
    mats <- lapply(c(pos = "pos", neg = "neg"), function(assay) {
        effect <- c(rep(1, 6L), rep(-1, 6L), rep(0, p - 12L))
        y <- stats::rnorm(p, 20, 2) + outer(effect, as.numeric(design$group == "B")) + outer(stats::rnorm(p, sd = 0.6), f) +
            matrix(stats::rnorm(p * 12L, sd = 0.5), p)
        dimnames(y) <- list(sprintf("%s_f%03d", assay, seq_len(p)), design$Sample_ID)
        if (missing > 0) {
            low <- which(y <= stats::quantile(y, 0.2), arr.ind = TRUE)
            low <- low[sample.int(nrow(low)), , drop = FALSE]
            key <- paste(low[, 1L], design$group[low[, 2L]])
            low <- low[stats::ave(seq_along(key), key, FUN = seq_along) <= 2L, , drop = FALSE]
            y[low[seq_len(min(nrow(low), round(missing * length(y)))), , drop = FALSE]] <- NA_real_
        }
        y
    })
    list(design = design, mats = mats)
}

kit224AssayTable <- function(m) {
    table <- tibble::as_tibble(as.data.frame(m, check.names = FALSE))
    table <- tibble::add_column(table, database_identifier = rownames(m), .before = 1L)
    tibble::add_column(table, name = rownames(m), .after = 1L)
}

kit224Object <- function(fixture, omic = "lipid") {
    tables <- lapply(fixture$mats, kit224AssayTable)
    if (identical(omic, "lipid")) {
        methods::new("LipidomicsAssayData", lipid_data = tables, lipid_id_column = "database_identifier",
            annotation_id_column = "name", design_matrix = fixture$design, sample_id = "Sample_ID", group_id = "group", args = list())
    } else {
        methods::new("MetaboliteAssayData", metabolite_data = tables, metabolite_id_column = "database_identifier",
            annotation_id_column = "name", design_matrix = fixture$design, sample_id = "Sample_ID", group_id = "group", args = list())
    }
}

kit224Slot <- function(object) if (methods::is(object, "LipidomicsAssayData")) "lipid_data" else "metabolite_data"

kit224Values <- function(table, samples) {
    samples <- intersect(samples, names(table))
    m <- as.matrix(as.data.frame(table[, samples], check.names = FALSE))
    rownames(m) <- table[["database_identifier"]]
    m
}

# The legacy correction as the application calls it: every feature a potential control, k = 1 per assay.
kit224Legacy <- function(object, fixture) {
    kit224Quiet(ruvIII_C_Varying(object, ruv_grouping_variable = "group", ruv_number_k = list(pos = 1L, neg = 1L),
        ctrl = lapply(fixture$mats, rownames)))
}

# Per assay, the largest change of a feature's mean (over the samples finite before and after) made by the correction.
kit224LevelChange <- function(corrected, fixture) {
    samples <- fixture$design$Sample_ID
    vapply(names(fixture$mats), function(assay) {
        after <- kit224Values(methods::slot(corrected, kit224Slot(corrected))[[assay]], samples)
        before <- fixture$mats[[assay]][rownames(after), colnames(after), drop = FALSE]
        both <- is.finite(after) & is.finite(before)
        mean_in <- function(m) rowSums(replace(m, !both, 0)) / rowSums(both)
        max(abs(mean_in(after) - mean_in(before)), na.rm = TRUE)
    }, numeric(1))
}

# --- E1, E1m, E2: the S4 methods keep each feature's level ---------------------------------------------------------------

test_that("KIT-224 E1: the lipid ruvIII_C_Varying keeps each feature's mean (complete data, k = 1)", {
    fixture <- kit224Mats(224000001L)
    change <- kit224LevelChange(kit224Legacy(kit224Object(fixture, "lipid"), fixture), fixture)
    for (assay in names(change)) {
        expect_lt(change[[assay]], 1e-9, label = sprintf("E1 lipid %s: largest |mean after - mean before| = %.6g", assay,
            change[[assay]]))
    }
})

test_that("KIT-224 E1m: the lipid ruvIII_C_Varying keeps each feature's mean over its corrected samples (5 % missing)", {
    fixture <- kit224Mats(224000002L, missing = 0.05)
    expect_true(anyNA(fixture$mats$pos) && anyNA(fixture$mats$neg))
    change <- kit224LevelChange(kit224Legacy(kit224Object(fixture, "lipid"), fixture), fixture)
    for (assay in names(change)) {
        expect_lt(change[[assay]], 1e-9, label = sprintf("E1m lipid %s: largest |mean after - mean before| = %.6g", assay,
            change[[assay]]))
    }
})

test_that("KIT-224 E2: the metabolomics ruvIII_C_Varying keeps each feature's mean, as the lipid method", {
    fixture <- kit224Mats(224000001L)
    metab <- kit224Legacy(kit224Object(fixture, "metab"), fixture)
    change <- kit224LevelChange(metab, fixture)
    for (assay in names(change)) {
        expect_lt(change[[assay]], 1e-9, label = sprintf("E2 metab %s: largest |mean after - mean before| = %.6g", assay,
            change[[assay]]))
    }
    lipid <- kit224Legacy(kit224Object(fixture, "lipid"), fixture)
    expect_identical(metab@metabolite_data, lipid@lipid_data)
})

# --- E4 (i), (ii): the marker and the earlier rule ---------------------------------------------------------------------

test_that("KIT-224 E4 (i): the replay marker ruv_feature_levels exists; the legacy declarations fix it", {
    expect_identical(workflowReplayRuleMarkers$ruv_feature_levels,
        list(option = KIT224_OPTION, current = "restored", earlier = "shifted"))
    marked <- list(parameters = list(list(name = "ruv_feature_levels", source = "fixed",
        value = encodeWorkflowRecordValue("restored"))))
    registry <- list(parameters = list(list(name = "ruv_feature_levels", source = "registry",
        value = encodeWorkflowRecordValue("restored"))))
    expect_identical(workflowReplayRule(marked, "ruv_feature_levels"), "restored")
    expect_identical(workflowReplayRule(registry, "ruv_feature_levels"), "shifted")
    expect_identical(workflowReplayRule(list(parameters = list()), "ruv_feature_levels"), "shifted")
    for (binding in list(lipidomicsDeclarationBinding, metabolomicsDeclarationBinding)) {
        expect_identical(omicsRuvOperation(binding)$fixed$ruv_feature_levels, "restored", info = binding$prefix)
        # the declared design restores in its own step; its declaration fixes no marker (the replay leaves it as computed)
        expect_null(omicsRuvOperation(binding, "declared")$fixed$ruv_feature_levels, info = binding$prefix)
    }
})

test_that("KIT-224 E4 (ii): under the earlier rule the S4 methods return RUVIIIC::RUVIII_C_Varying's values exactly", {
    fixture <- kit224Mats(224000001L)
    design <- fixture$design
    for (omic in c("lipid", "metab")) {
        shifted <- withr::with_options(stats::setNames(list("shifted"), KIT224_OPTION),
            kit224Legacy(kit224Object(fixture, omic), fixture))
        for (assay in names(fixture$mats)) {
            after <- kit224Values(methods::slot(shifted, kit224Slot(shifted))[[assay]], design$Sample_ID)
            Y <- t(fixture$mats[[assay]][, design$Sample_ID])
            M <- ruv::replicate.matrix(design$group)
            rownames(M) <- rownames(Y)
            reference <- t(RUVIIIC::RUVIII_C_Varying(k = 1L, Y = Y, M = M, toCorrect = colnames(Y),
                potentialControls = colnames(Y), progress = FALSE))
            expect_identical(after, reference, info = paste(omic, assay))
        }
        expect_gt(max(kit224LevelChange(shifted, fixture)), 0.01, label = paste("E4 (ii) the earlier rule's level shift,", omic))
    }
})

# --- G1, G2: fold changes; the declared path ---------------------------------------------------------------------------

test_that("KIT-224 G1: the legacy DA's fold changes are the same with the levels kept and with the earlier shift", {
    fixture <- kit224Mats(224000001L)
    contrasts <- data.frame(contrasts = "groupB-groupA", friendly_names = "B vs A", stringsAsFactors = FALSE)
    da <- function(corrected) {
        kit224Quiet(runLipidsDA(corrected, contrasts, "~ 0 + group", 0.05, 0, TRUE, TRUE))$da_lipids_long
    }
    object <- kit224Object(fixture, "lipid")
    kept <- da(kit224Legacy(object, fixture))
    shifted <- da(withr::with_options(stats::setNames(list("shifted"), KIT224_OPTION), kit224Legacy(object, fixture)))
    key <- function(d) paste(d$assay, d$comparison, d$lipid_id)
    expect_setequal(key(kept), key(shifted))
    at <- match(key(shifted), key(kept))
    expect_lt(max(abs(kept$logFC[at] - shifted$logFC)), 1e-10, label = "G1 largest |logFC difference|")
})

test_that("KIT-224 G2: the declared designs restore once: their corrected matrices do not depend on the option", {
    fixture <- kit224Mats(224000001L)
    supplied <- lapply(fixture$mats, function(m) rownames(m)[13:40])
    for (omic in c("lipid", "metab")) {
        object <- kit224Object(fixture, omic)
        run <- function() list(
            groups = kit224Quiet(runOmicsRuvFactors(object, "groups_all_features", 1L, "user")),
            supplied = kit224Quiet(runOmicsRuvFactors(object, "supplied_controls", 1L, "user", ruv_control_ids = supplied)))
        default <- run()
        shifted <- withr::with_options(stats::setNames(list("shifted"), KIT224_OPTION), run())
        for (design in names(default)) {
            expect_identical(default[[design]]$corrected, shifted[[design]]$corrected, info = paste(omic, design))
            expect_identical(default[[design]]$factors, shifted[[design]]$factors, info = paste(omic, design))
        }
        expect_lt(max(kit224LevelChange(default$groups$corrected, fixture)), 1e-9, label = paste("G2 declared levels,", omic))
    }
})

# --- E3, E4 (iii): the declared legacy step on the application's path ---------------------------------------------------

test_that("KIT-224 E3 and E4 (iii): the legacy step lipid_ruv (Automatic) keeps each feature's level; its record lists the marker", {
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    values <- kit172Values(172101L, planted = TRUE, missing = FALSE)$values
    journey <- kit172Journey(withr::local_tempdir(), kit172Object(values, "lipid"), "lipid", itsd = TRUE,
        method = "cyclicloess", offset = 0, ruv = kit172LegacyRuv("automatic"), da = FALSE)
    withr::defer(closeDIAStateAdapter(journey$owner))
    after <- kit172StateMatrix(journey, "lipid_ruv_corrected")
    before <- kit172StateMatrix(journey, "lipid_normalized")[rownames(after), colnames(after), drop = FALSE]
    both <- is.finite(after) & is.finite(before)
    mean_in <- function(m) rowSums(replace(m, !both, 0)) / rowSums(both)
    change <- max(abs(mean_in(after) - mean_in(before)), na.rm = TRUE)
    expect_lt(change, 1e-9, label = sprintf("E3 lipid_ruv Automatic: largest |mean after - mean before| = %.6g", change))
    record <- buildWorkflowAnalysisRecord(journey$owner)
    step <- Filter(function(step) identical(step$operation, "lipid_ruv"), record$steps)
    expect_length(step, 1L)
    marker <- Filter(function(p) identical(p$name, "ruv_feature_levels"), step[[1L]]$parameters)
    expect_length(marker, 1L)
    if (length(marker)) {
        expect_identical(marker[[1L]]$source, "fixed")
        expect_identical(decodeWorkflowRecordValue(marker[[1L]]$value), "restored")
    }
})

# --- E4 (iv): save, reopen and replay under each rule ------------------------------------------------------------------

# The replay project and journey of test-kit129-run-path.R (kit129ReplayProject, kit129ReplayInputs,
# kit129SaveReopenReplay), copied.
kit224ReplayProject <- function(source_dir, preconditions, env = parent.frame()) {
    harness <- new.env(parent = globalenv())
    sys.source(file.path(testthat::test_path("..", ".."), "tools", "profiling", "workflow_kit.R"), harness)
    api <- function(name) get(name, asNamespace("MultiScholaR"))
    copy <- withr::local_tempdir(.local_envir = env)
    stopifnot(system2("cp", c("-a", file.path(source_dir, "."), copy)) == 0L)
    store_root <- dirname(dirname(list.files(source_dir, "^store[.]json$", recursive = TRUE, all.files = TRUE, full.names = TRUE)[[1L]]))
    copy_root <- file.path(copy, substring(store_root, nchar(source_dir) + 2L))
    list(prepared = harness$kitPrepareReplayProject(copy_root, preconditions, api), origin = store_root, harness = harness, api = api)
}

kit224ReplayInputs <- function(record, source_dir) {
    files <- list.files(source_dir, recursive = TRUE, full.names = TRUE)
    inputs <- Filter(function(input) is.character(input$basename) && any(basename(files) == input$basename), record$inputs)
    stats::setNames(vapply(inputs, function(input) files[basename(files) == input$basename][[1L]], ""),
        vapply(inputs, function(input) input$role, ""))
}

kit224SaveReopenReplay <- function(owner, source_dir, preconditions, paths = list(), env = parent.frame()) {
    json <- file.path(withr::local_tempdir(.local_envir = env), "analysis_record.json")
    writeWorkflowAnalysisRecord(owner, json)
    record <- jsonlite::read_json(json, simplifyVector = FALSE)
    project <- kit224ReplayProject(source_dir, preconditions, env = env)
    replay <- kit224Quiet(project$harness$kitReplayJourney(record, project$prepared, project$origin,
        kit224ReplayInputs(record, source_dir), paths, withr::local_tempdir(.local_envir = env), project$api))
    list(record = record, replay = replay)
}

kit224RuvMarkers <- function(record, prefix) {
    step <- Filter(function(step) identical(step$operation, paste0(prefix, "_ruv")), record$steps)[[1L]]
    Filter(function(p) identical(p$name, "ruv_feature_levels"), step$parameters)
}

test_that("KIT-224 E4 (iv): a record made under each rule replays without a mismatch; only the current rule lists the marker", {
    for (omic in c("lipidomics", "metabolomics")) {
        binding <- omicsDeclarationBindings()[[omic]]
        slot_name <- if (identical(omic, "lipidomics")) "lipid_data" else "metabolite_data"
        object <- kit224Object(kit224Mats(224000003L), if (identical(omic, "lipidomics")) "lipid" else "metab")
        owners <- list()
        for (rule in c("restored", "shifted")) {
            root <- withr::local_tempdir()
            fixture <- if (identical(omic, "lipidomics")) lipidDeclarationsOwner(root, object = object) else
                metabDeclarationsOwner(root, object = object)
            withr::defer(closeDIAStateAdapter(fixture$owner))
            options <- if (identical(rule, "shifted")) stats::setNames(list("shifted"), KIT224_OPTION) else list()
            withr::with_options(options, kit224Quiet(runWorkflowOmicsReferenceRuv(binding, fixture$owner,
                lipidDeclarationsRuvParameters(), fixture$workflow$config_list)))
            owners[[rule]] <- fixture
        }
        corrected <- lapply(owners, function(fixture) {
            methods::slot(fixture$owner$getState(paste0(binding$prefix, "_ruv_corrected")), slot_name)
        })
        # non-vacuity: the two rules' corrected tables differ
        expect_false(identical(corrected$restored, corrected$shifted), info = paste(omic, "the two rules' corrected tables differ"))
        for (rule in names(owners)) {
            fixture <- owners[[rule]]
            preconditions <- unlist(fixture$owner$.context$snapshot$state_history)[1:2]
            saved <- kit224SaveReopenReplay(fixture$owner, fixture$source_dir, preconditions, fixture$paths)
            expect_identical(saved$replay$states, unlist(fixture$owner$.context$snapshot$state_history), info = paste(omic, rule))
            expect_identical(as.character(saved$replay$mismatches), character(), info = paste(omic, rule))
            expect_length(kit224RuvMarkers(saved$record, binding$prefix), if (identical(rule, "restored")) 1L else 0L)
        }
        expect_null(getOption(KIT224_OPTION))
    }
})

# --- E5: the transcriptomics RUVIII export file and integration input keep each gene's level (OD-2, T-a') -------------

# test-integration-input-transcriptomics.R's owner and save (integrationTrxOwner, integrationTrxSave), copied: the
# KIT-040 corpus trx_edge_12 on a session owner, the expression filter and QC finalise done.
kit224TrxOwner <- function(env = parent.frame()) {
    skipWithoutTrxDeclReference()
    reference <- trxDeclReference("trx_edge_12", env = env)
    fixture <- trxDeclarationsOwner(withr::local_tempdir(.local_envir = env), reference$dir)
    withr::defer(closeDIAStateAdapter(fixture$owner), envir = env)
    kit224Quiet({
        runTrxReferenceExpressionFilter(fixture$owner, list())
        commitWorkflowOmicsStoreStep(fixture$owner, "trx_qc_finalize", "trx_qc_complete", list(), "QC processing complete")
    })
    fixture
}

kit224TrxSave <- function(fixture, root, label = "kit224") {
    paths <- list(source_dir = file.path(root, "scripts", "transcriptomics"), integration_dir = file.path(root, "integration"),
        base_dir = root)
    dir.create(paths$source_dir, recursive = TRUE, showWarnings = FALSE)
    failure <- NULL
    kit224Quiet(saveWorkflowSummaryParameters(omicGuiConfig("transcriptomics")$summary, fixture$workflow,
        projectDirs = list(transcriptomics = paths), omicType = "transcriptomics", label = label, description = "d",
        notify = function(message, ...) NULL, done = function(value) NULL,
        failed = function(error) failure <<- error, createWorkflowArgsFn = function(...) "study_parameters.txt",
        writeAnalysisRecordFn = function(...) NULL, catFn = function(...) NULL))
    if (!is.null(failure)) stop(failure)
    list(dir = file.path(paths$integration_dir, "transcriptomics"),
        tsv = file.path(paths$integration_dir, "transcriptomics_normalized_data.tsv"))
}

test_that("KIT-224 E5: the RUVIII export file and the corrected integration input keep each gene's TMM logCPM mean", {
    fixture <- kit224TrxOwner()
    owner <- fixture$owner
    ruvs <- kit224Quiet(runTrxReferenceRuvs(owner, trxRuvGridParameters(), list()))
    export <- kit224Quiet(runTrxReferenceRuviiiExport(owner, "group"))
    log_cpm <- withWorkflowStateScope(owner, function() trxSavedMatrix(owner, ruvs$logcpm$value$log_cpm))
    new_y <- t(export$value$newY)
    saved <- kit224TrxSave(fixture, withr::local_tempdir())
    # (a) the export file (samples by genes, a Run column)
    exported <- as.data.frame(vroom::vroom(saved$tsv, show_col_types = FALSE))
    file_values <- t(as.matrix(exported[-1L]))
    colnames(file_values) <- exported$Run
    # (b) the integration input's corrected matrix (genes by runs)
    input <- readIntegrationInputWithMatrices(saved$dir)
    corrected <- input$views$transcriptomics$corrected_matrix
    input_values <- as.matrix(corrected[-1L])
    rownames(input_values) <- corrected$feature_id
    for (case in list(list(name = "file", values = file_values), list(name = "integration input", values = input_values))) {
        values <- case$values
        genes <- rownames(values)
        runs <- colnames(values)
        level <- max(abs(rowMeans(values) - rowMeans(log_cpm[genes, runs, drop = FALSE])))
        expect_lt(level, 1e-9, label = sprintf("E5 %s: largest |gene mean - TMM logCPM mean| = %.6g", case$name, level))
        # only a constant per gene is added to RUVIII's values
        within <- values - new_y[genes, runs, drop = FALSE]
        expect_lt(max(apply(within, 1L, function(x) diff(range(x)))), 1e-9, label = paste("E5", case$name, "within-gene change"))
    }
    expect_identical(input$manifest$correction$levels, "restored_to_tmm_logcpm")
})

# --- E6: the repaired legacy route does not depend on the unit (KIT-172 T1b) -------------------------------------------

test_that("KIT-224 E6: the legacy Manual route's calls do not change when every value is multiplied by 2^10 (T1b, seed 172102)", {
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    values <- kit172Values(172102L, planted = TRUE, missing = FALSE)$values
    run <- function(input) {
        root <- tempfile("kit224-e6-")
        journey <- kit172Quiet(kit172Journey(root, kit172Object(input, "lipid"), "lipid", itsd = TRUE,
            method = "cyclicloess", filter = TRUE, offset = 0, ruv = kit172LegacyRuv("manual")))
        on.exit({
            closeDIAStateAdapter(journey$owner)
            unlink(root, recursive = TRUE)
        }, add = TRUE)
        kit172DaTable(journey)
    }
    reference <- run(values)
    scaled <- run(values * 2^10)
    calls <- function(da) ifelse(!is.na(da$fdr_qvalue) & da$fdr_qvalue < 0.05, sign(da$logFC), 0)
    key <- function(d) paste(d$assay, d$comparison, d$feature)
    expect_setequal(key(scaled), key(reference))
    at <- match(key(reference), key(scaled))
    changes <- sum(calls(scaled)[at] != calls(reference), na.rm = TRUE)
    max_dt <- max(abs(scaled$t[at] - reference$t), na.rm = TRUE)
    max_dlogfc <- max(abs(scaled$logFC[at] - reference$logFC), na.rm = TRUE)
    expect_identical(changes, 0L, label = sprintf("E6 calls changed (of %d; max |delta t| %.4g)", nrow(reference), max_dt))
    expect_lt(max_dt, 1e-6, label = sprintf("E6 max |delta t| = %.6g", max_dt))
    expect_lt(max_dlogfc, 1e-9, label = sprintf("E6 max |delta logFC| = %.6g", max_dlogfc))
})
