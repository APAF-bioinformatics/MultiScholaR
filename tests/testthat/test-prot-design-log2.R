# KIT-177: plain log2 at the label-free and TMT protein design step (ticket sections 3 and 5: T1a, T2, T3, T3b, T4, T6,
# T7, T8, T9). The design state is made by the step under test, `buildProtDesignStateCheckpoint`; where a test needs the
# earlier state it makes it here, the old way (the wide table, log2(x + 1), `ProteinQuantitativeData`). The route's
# later steps run as declared operations in process, as the tabs run them.

# --- Fixture F177 and the route (fixed in the ticket, section 5) ------------------------------------------------------

#' F177: 400 proteins, two groups of 5 runs, levels log-uniform over `range`, log2-normal noise SD 0.25 in every cell;
#' per protein and group 0 to 2 cells missing at random (10 % on average) and, among the observed, cells set to exactly
#' 1 with probability 2 % (at most 2 per protein and group); no duplicate accession.
k177F177 <- function(seed, n = 400L, range = c(0.5, 1e7)) {
    set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
    level <- exp(stats::runif(n, log(range[[1L]]), log(range[[2L]])))
    values <- level * 2^matrix(stats::rnorm(n * 10L, 0, 0.25), n, 10L)
    for (i in seq_len(n)) for (g in list(1:5, 6:10)) {
        missing <- sample(g, sample(0:2, 1L, prob = c(0.6, 0.3, 0.1)))
        values[i, missing] <- NA_real_
        observed <- setdiff(g, missing)
        ones <- observed[stats::runif(length(observed)) < 0.02]
        values[i, utils::head(ones, 2L)] <- 1
    }
    runs <- sprintf("S%02d", 1:10)
    list(ids = sprintf("P%05d", seq_len(n)), runs = runs, values = values,
        design = data.frame(Run = runs, group = rep(c("control", "treated"), each = 5L), replicates = rep(1:5, 2L),
            stringsAsFactors = FALSE))
}

#' The long table a reader returns for one shape (column names as the MaxQuant and Proteome Discoverer readers write).
k177Long <- function(sim, shape, scale = 1) {
    quantity <- if (identical(shape, "TMT")) "Abundance" else "Intensity"
    long <- data.frame(Protein.Ids = rep(sim$ids, times = length(sim$runs)), Run = rep(sim$runs, each = length(sim$ids)),
        stringsAsFactors = FALSE)
    long[[quantity]] <- as.vector(sim$values) * scale
    list(data = long, mapping = list(protein_col = "Protein.Ids", run_col = "Run", quantity_col = quantity))
}

#' Close an owner when the test ends (the value captured now, so a loop variable reassigned later is not reread).
k177Defer <- function(manager, env = parent.frame()) {
    force(manager)
    withr::defer(closeDIAStateAdapter(manager), envir = env)
}

k177Owner <- function(root, type) {
    dir.create(root, recursive = TRUE, showWarnings = FALSE)
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, paste0("k177-", basename(root)), type),
        WorkflowState$new(audit_enabled = FALSE))
    manager$.context$reference_states <- TRUE
    manager$setWorkflowType(type)
    paths <- list(source_dir = file.path(root, "src"), protein_qc_dir = file.path(root, "pqc"),
        da_output_dir = file.path(root, "da"), publication_graphs_dir = file.path(root, "pub"), time_dir = file.path(root, "time"),
        results_dir = file.path(root, "res"))
    for (p in paths) dir.create(p, recursive = TRUE, showWarnings = FALSE)
    manager$.context$experiment_paths <- paths
    manager
}

#' The design state by the step under test on a fresh owner.
k177Design <- function(root, type, long, design) {
    manager <- k177Owner(root, type)
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    workflow$data_cln <- long$data
    workflow$column_mapping <- long$mapping
    workflow$design_matrix <- design
    workflow$config_list <- list(globalParameters = list(workflow_type = type))
    withDIAStateScope(manager, function() buildProtDesignStateCheckpoint(workflow, type, "Save"))
    list(manager = manager, workflow = workflow)
}

#' The design state the old way (log2(x + 1)), as KIT-177 T5 makes it.
k177DesignOld <- function(root, type, long, design) {
    manager <- k177Owner(root, type)
    id <- long$mapping$protein_col
    wide <- tidyr::pivot_wider(long$data, id_cols = dplyr::all_of(id), names_from = dplyr::all_of(long$mapping$run_col),
        values_from = dplyr::all_of(long$mapping$quantity_col))
    for (run in setdiff(names(wide), id)) wide[[run]] <- log2(wide[[run]] + 1)
    design$tech_rep_group <- paste(design$group, design$replicates, sep = "_")
    object <- ProteinQuantitativeData(protein_quant_table = wide, protein_id_column = id,
        protein_id_table = dplyr::distinct(wide, !!rlang::sym(id)), design_matrix = design, sample_id = "Run",
        group_id = "group", technical_replicate_id = "tech_rep_group", args = list(globalParameters = list(workflow_type = type)))
    setDIAStepProvenance(manager, "design_builder", fixed = list(pivot = "protein x run", transform = "log2(x + 1)"))
    manager$saveState("protein_s4_initial", object, list(), "design")
    list(manager = manager)
}

k177Matrix <- function(object, ids, runs) {
    table <- object@protein_quant_table
    as.matrix(as.data.frame(table)[match(ids, table[[object@protein_id_column]]), runs, drop = FALSE])
}

#' The route after the design state: QC at the given settings, normalisation (RUV skipped), correlation filter, DA.
#' Returns the design and normalised objects, the intensity-filtered protein ids and the DA columns (id, logFC, t, p).
k177Route <- function(manager, norm_method = "cyclicloess", intensity = TRUE,
    qc = list(min_reps_per_group = 2, min_groups = 2, percentile = 1), contrast = c("treated", "control")) {
    workflow <- list2env(list(state_manager = manager, qc_params = list(protein_qc = list()), config_list = list()))
    paths <- manager$.context$experiment_paths
    out <- withDIAStateScope(manager, function() {
        design <- manager$getState("protein_s4_initial")
        finish <- function(value) finishDIAReferenceProteinStep(workflow, value, publish = FALSE)
        finish(computeDIAReferenceAccessionCleanup(manager, getDIAStateReference(manager),
            list(delimiter = ";", aggregation_method = "mean")))
        if (isTRUE(intensity)) finish(computeDIAReferenceProteinIntensityFilter(manager, getDIAStateReference(manager),
            list(strict_mode = FALSE, min_reps_per_group = qc$min_reps_per_group, min_groups = qc$min_groups,
                proteins_intensity_cutoff_percentile = qc$percentile)))
        filtered <- manager$getState(manager$current_state)
        finish(computeDIAReferenceProteinDuplicateRemoval(manager, getDIAStateReference(manager), list(aggregation_method = "mean")))
        completeDIAProteinAction(workflow, computeDIAReferenceProteinReplicateFilter(manager, getDIAStateReference(manager),
            list(grouping_variable = "group", parallel_cores = 1L), variant = "rows"), publish = FALSE)
        parameters <- list(norm_method = norm_method, ruv_mode = "skip", auto_percentage_min = NULL, auto_percentage_max = NULL,
            separation_metric = NULL, k_penalty_weight = NULL, adaptive_k_penalty = NULL, ruv_percentage = NULL, ruv_k = NULL,
            grouping = "group")
        run <- computeDIAReferenceNormalization(manager, getDIAStateReference(manager), parameters)
        commit <- diaNormalizationCommit(manager, parameters, run)
        if (!is.null(commit)) commitDIAReferenceState(manager, commit$state_name, commit$checkpoint, commit$config, commit$description)
        normalised <- manager$getState("normalized")
        correlated <- computeDIAReferenceCorrelation(manager, run$value$final, "group", 0.5)
        commit <- diaCorrelationCommit(manager, correlated$final, 0.5, FALSE)
        if (!is.null(commit)) commitDIAReferenceState(manager, commit$state_name, commit$checkpoint, commit$config, commit$description)
        name <- paste0(contrast[[1L]], ".vs.", contrast[[2L]])
        contrasts <- data.frame(contrasts = name, full_format = sprintf("%s=group%s-group%s", name, contrast[[1L]], contrast[[2L]]),
            friendly_names = name, stringsAsFactors = FALSE)
        da_data <- new.env(parent = emptyenv())
        da_data$current_s4_object <- getDIAStateReference(manager)
        da_data$contrasts_available <- contrasts$contrasts
        workflow$contrasts_tbl <- contrasts
        manifest <- runDIAReferenceDaAnalysis(list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0),
            da_data, workflow, paths, contrasts = contrasts, publish = FALSE)
        tables <- readDIAContrastOutput(manifest$results[[1L]], c("da_proteins_long", "contrasts_results_table"))
        list(design = design, normalised = normalised, filtered_ids = as.character(filtered@protein_quant_table[[filtered@protein_id_column]]),
            correlated = manager$getState("correlation_filtered"), long = as.data.frame(tables$da_proteins_long),
            results = as.data.frame(tables$contrasts_results_table))
    })
    out$da <- k177DaTable(out, out$design@protein_id_column)
    out
}

k177DaTable <- function(route, id) {
    pick <- function(table, names) {
        hit <- names[names %in% names(table)]
        if (length(hit)) table[[hit[[1L]]]] else NULL
    }
    long <- route$long
    ids <- as.character(pick(long, c(id, "uniprot_acc")))
    out <- data.frame(id = ids, logFC = as.numeric(pick(long, c("log2FC", "logFC"))), p = as.numeric(pick(long, c("raw_pvalue",
        "P.Value"))), stringsAsFactors = FALSE)
    t <- pick(long, "t")
    if (is.null(t)) {
        rt <- route$results
        rid <- if (id %in% names(rt)) as.character(rt[[id]]) else rownames(rt)
        t <- as.numeric(pick(rt, "t"))[match(ids, rid)]
    }
    out$t <- as.numeric(t)
    out
}

k177Calls <- function(da) {
    tested <- !is.na(da$p)
    q <- rep(NA_real_, nrow(da))
    q[tested] <- stats::p.adjust(da$p[tested], "BH")
    da$id[tested & q < 0.05]
}

# --- T1a unit equivariance (asserted) --------------------------------------------------------------------------------

test_that("KIT-177 T1a: multiplying the table by s adds log2(s) at the design step and changes no fold change", {
    skip_if_not(workflowOwnerCapable())
    sim <- k177F177(177001L)
    for (shape in c("TMT", "LFQ")) for (norm_method in c("none", "quantile", "cyclicloess")) {
        runAt <- function(s, intensity = TRUE, keep = NULL) {
            long <- k177Long(sim, shape, s)
            if (!is.null(keep)) long$data <- long$data[long$data$Protein.Ids %in% keep, , drop = FALSE]
            made <- k177Design(tempfile("k177-t1a-"), shape, long, sim$design)
            on.exit(closeDIAStateAdapter(made$manager), add = TRUE)
            utils::capture.output(route <- suppressMessages(suppressWarnings(k177Route(made$manager, norm_method, intensity))))
            route
        }
        one <- runAt(1)
        one_unfiltered <- runAt(1, intensity = FALSE, keep = one$filtered_ids)
        for (s in c(2^-10, 2^10, 0.001, 1000)) {
            exact <- s %in% c(2^-10, 2^10)
            reference <- if (exact) one else one_unfiltered
            other <- if (exact) runAt(s) else runAt(s, intensity = FALSE, keep = one$filtered_ids)
            label <- sprintf("%s %s s=%g", shape, norm_method, s)
            ids <- reference$design@protein_quant_table$Protein.Ids
            d1 <- k177Matrix(reference$design, ids, sim$runs)
            ds <- k177Matrix(other$design, ids, sim$runs)
            expect_identical(is.na(ds), is.na(d1), label = paste(label, "(1) missing cells"))
            expect_lte(max(abs(ds - d1 - log2(s)), na.rm = TRUE), 1e-12, label = paste(label, "(1) design-step shift"))
            nids <- reference$normalised@protein_quant_table$Protein.Ids
            runs <- intersect(sim$runs, names(reference$normalised@protein_quant_table))
            n1 <- k177Matrix(reference$normalised, nids, runs)
            ns <- k177Matrix(other$normalised, nids, runs)
            expect_lte(max(abs(ns - n1 - log2(s)), na.rm = TRUE), 1e-9, label = paste(label, "(2) normalised shift"))
            both <- merge(reference$da, other$da, by = "id", suffixes = c(".1", ".s"))
            expect_identical(nrow(both), nrow(reference$da), label = paste(label, "same proteins tested"))
            expect_lte(max(abs(both$logFC.s - both$logFC.1), na.rm = TRUE), 1e-9, label = paste(label, "(3) logFC"))
            rel_t <- max(abs(both$t.s - both$t.1) / pmax(abs(both$t.1), .Machine$double.xmin), na.rm = TRUE)
            zero <- both$p.1 == 0 & both$p.s == 0
            log_p <- abs(log(both$p.s) - log(both$p.1)) - (1e-6 * abs(log(both$p.1)) + 1e-12)
            excess <- max(c(-Inf, log_p[!zero & !is.na(log_p)]))
            # Exceedance rule: with (1) to (3) holding, an excess in t or p alone is reported as a numerical finding.
            if (rel_t > 1e-6 || excess > 0)
                message(sprintf("KIT177-T1A-TP %s max_rel_t=%.3g max_logp_excess=%.3g (limma treat/eBayes)", label, rel_t, excess))
            expect_true(all(!is.na(both$t.1)), label = paste(label, "t reported"))
        }
    }
})

# --- T2 fold change --------------------------------------------------------------------------------------------------

test_that("KIT-177 T2: a protein at 10 against 80 reports a log2 fold change of 3", {
    skip_if_not(workflowOwnerCapable())
    sim <- k177F177(177001L)
    sim$values[1L, ] <- c(rep(10, 5L), rep(80, 5L))
    made <- k177Design(withr::local_tempdir(), "TMT", k177Long(sim, "TMT"), sim$design)
    k177Defer(made$manager)
    utils::capture.output(route <- suppressMessages(suppressWarnings(k177Route(made$manager, "none", intensity = FALSE))))
    expect_equal(route$da$logFC[route$da$id == sim$ids[[1L]]], 3, tolerance = 1e-12)
})

# --- T3 the rule per cell (direct call) ------------------------------------------------------------------------------

k177CellCase <- function(value) {
    long <- data.frame(Protein.Ids = c("P1", "P1", "P2", "P2", "P3", "P4"), Run = c("A1", "B1", "A1", "B1", "A1", "B1"),
        Intensity = c(value, 100, 200, 300, 400, 500), stringsAsFactors = FALSE)
    list(data = long, mapping = list(protein_col = "Protein.Ids", run_col = "Run", quantity_col = "Intensity"))
}

test_that("KIT-177 T3: 0, NaN and NA become missing, counted apart; a negative or infinite quantity is refused", {
    skip_if_not(workflowOwnerCapable())
    design <- data.frame(Run = c("A1", "B1"), group = c("A", "B"), replicates = c(1L, 1L), stringsAsFactors = FALSE)
    long <- k177CellCase(0)
    long$data$Intensity[long$data$Protein.Ids == "P3"] <- NaN
    made <- expect_no_warning(k177Design(withr::local_tempdir(), "LFQ", long, design))
    k177Defer(made$manager)
    object <- made$manager$getState("protein_s4_initial")
    m <- k177Matrix(object, c("P1", "P2", "P3", "P4"), c("A1", "B1"))
    expect_true(is.na(m[1L, 1L]))                 # 0
    expect_true(is.na(m[3L, 1L]))                 # NaN
    expect_true(is.na(m[4L, 1L]))                 # NA filled by the pivot
    expect_false(any(is.nan(m)))
    expect_false(any(is.infinite(m)))
    expect_equal(m[2L, ], c(A1 = log2(200), B1 = log2(300)))
    fixed <- made$manager$.context$snapshot$step_provenance$protein_s4_initial$fixed
    expect_identical(as.integer(fixed$nonpositive_set_to_na), 1L)
    expect_identical(as.integer(fixed$nonfinite_set_to_na), 1L)
    text <- protDesignLog2NotificationText(made$workflow$prot_design_log2_counts)
    expect_identical(text, paste("1 values of 0 and 1 values that were not a number were read as not measured (missing).",
        "They are not used in the analysis."))
    for (value in c(-5, Inf, -Inf)) {
        manager <- k177Owner(withr::local_tempdir(), "LFQ")
        k177Defer(manager)
        workflow <- new.env(parent = emptyenv())
        workflow$state_manager <- manager
        workflow$data_cln <- k177CellCase(value)$data
        workflow$column_mapping <- k177CellCase(value)$mapping
        workflow$design_matrix <- design
        workflow$config_list <- list()
        before <- list(head = manager$current_state, history = manager$.context$snapshot$state_history)
        error <- tryCatch(withDIAStateScope(manager, function() buildProtDesignStateCheckpoint(workflow, "LFQ", "Save")),
            error = function(e) e)
        expect_s3_class(error, "prot_design_quantity_invalid")
        expect_match(conditionMessage(error), "Run 'A1' contains negative or infinite values (for example", fixed = TRUE)
        expect_match(conditionMessage(error), format(value), fixed = TRUE)
        expect_identical(manager$current_state, before$head)
        expect_identical(manager$.context$snapshot$state_history, before$history)
    }
})

# --- T3b through the route -------------------------------------------------------------------------------------------

makeK177FixtureSession <- function() {
    session <- new.env(parent = emptyenv())
    session$userData <- new.env(parent = emptyenv())
    session$onFlush <- session$onFlushed <- session$onSessionEnded <- function(...) function() invisible(NULL)
    session
}

test_that("KIT-177 T3b: zeros from a Spectronaut table and from a stored data_cln.tab are missing, recorded and notified", {
    skip_if_not(workflowOwnerCapable())
    notes <- list()
    testthat::local_mocked_bindings(showNotification = function(ui, ..., type = "default") {
        notes[[length(notes) + 1L]] <<- list(text = as.character(ui), type = type)
        invisible(NULL)
    }, .package = "shiny")
    root <- withr::local_tempdir()
    runs <- c("WT_1", "WT_2", "WT_3", "KO_1", "KO_2", "KO_3")
    ids <- sprintf("P%05d", 1:8)
    values <- matrix(round(2^(10 + seq_len(48) / 10), 1), 8, 6, dimnames = list(ids, runs))
    values[cbind(c(2, 5, 7), c(1, 4, 6))] <- 0
    design <- data.frame(Run = runs, group = rep(c("WT", "KO"), each = 3L), replicates = rep(1:3, 2L), stringsAsFactors = FALSE)
    # (i) Spectronaut, read by its reader, the design saved through the builder's save helper.
    report <- file.path(root, "spectronaut.tsv")
    utils::write.table(data.frame(PG.ProteinGroups = rep(ids, 6L), R.FileName = rep(runs, each = 8L),
        PG.Quantity = as.vector(values), PG.Qvalue = 0.001), report, sep = "\t", quote = FALSE, row.names = FALSE)
    fasta <- file.path(root, "proteins.fasta")
    writeLines(paste0(">sp|", ids, "|PROT", seq_along(ids), "\nMPEPTIDEK"), fasta)
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- WorkflowState$new(audit_enabled = FALSE)
    workflow$config_list <- list(globalParameters = list(workflow_type = "LFQ"))
    paths <- list(source_dir = file.path(root, "src"), protein_qc_dir = file.path(root, "pqc"), peptide_qc_dir = file.path(root, "peq"))
    for (p in paths) dir.create(p, recursive = TRUE, showWarnings = FALSE)
    initializeDIANativeImport(workflow, paths, makeK177FixtureSession(), "spectronaut", "iq", data_type = "protein")
    manager <- workflow$state_manager
    k177Defer(manager)
    testthat::local_mocked_bindings(persistProtDesignBuilderArtifacts = function(...) invisible(NULL),
        completeProtDesignPostCheckpoint = function(...) invisible(NULL), .package = "MultiScholaR")
    withDIAStateScope(manager, function() {
        registerDIAInputFiles(manager, c(spectronaut_report = report, fasta = fasta))
        imported <- suppressMessages(importSpectronautData(report))
        applyProtImportResultToWorkflow(workflow, imported, "spectronaut", fasta, logInfo = function(...) NULL,
            showNotification = function(...) NULL)
        completeDIAInputFiles(manager)
        workflow$design_matrix <- design
    })
    saved <- withDIAStateScope(manager, function() suppressMessages(runProtDesignBuilderSaveFlow(list(design_matrix = design,
        contrasts_tbl = NULL, config_list = list()), workflow, paths, session = NULL)))
    expect_true(saved)
    object <- manager$getState("protein_s4_initial")
    m <- k177Matrix(object, ids, runs)
    expect_identical(which(is.na(m)), which(values == 0))
    record <- buildWorkflowAnalysisRecord(manager)
    step <- Filter(function(s) identical(s$state, "protein_s4_initial"), record$steps)[[1L]]
    counted <- Filter(function(p) identical(p$name, "nonpositive_set_to_na"), step$parameters)
    expect_length(counted, 1L)
    expect_identical(as.integer(counted[[1L]]$value$value), 3L)
    texts <- vapply(notes, `[[`, "", "text")
    expect_true("3 values of 0 were read as not measured (missing). They are not used in the analysis." %in% texts)
    # (ii) A stored data_cln.tab with zeros, read by the "Import existing design" helpers.
    notes <- list()
    import_dir <- file.path(root, "stored")
    dir.create(import_dir)
    utils::write.table(data.frame(Protein.Ids = rep(ids, 6L), Run = rep(runs, each = 8L), Intensity = as.vector(values)),
        file.path(import_dir, "data_cln.tab"), sep = "\t", quote = FALSE, row.names = FALSE)
    utils::write.table(design, file.path(import_dir, "design_matrix.tab"), sep = "\t", quote = FALSE, row.names = FALSE)
    stored_paths <- list(source_dir = file.path(root, "stored-src"))
    dir.create(stored_paths$source_dir)
    writeLines(c("[globalParameters]", "workflow_type = LFQ"), file.path(stored_paths$source_dir, "config.ini"))
    owner <- k177Owner(file.path(root, "stored-owner"), "LFQ")
    k177Defer(owner)
    stored <- new.env(parent = emptyenv())
    stored$state_manager <- owner
    withDIAStateScope(owner, function() {
        artifacts <- suppressMessages(loadProtDesignImportedConfigAndTables(stored, stored_paths,
            file.path(import_dir, "design_matrix.tab"), file.path(import_dir, "data_cln.tab"),
            file.path(import_dir, "contrast_strings.tab")))
        suppressMessages(runProtDesignImportConfirmationFlow(stored, stored_paths, list(importPath = import_dir,
            designFile = file.path(import_dir, "design_matrix.tab"), fastaPath = NULL), artifacts, 9606, "human", NULL,
            hydrateFastaSidecar = function(...) NULL, hydrateUniprotSidecar = function(...) NULL,
            completePostCheckpointFn = function(...) NULL))
    })
    object <- owner$getState("protein_s4_initial")
    expect_identical(which(is.na(k177Matrix(object, ids, runs))), which(values == 0))
    expect_identical(as.integer(owner$.context$snapshot$step_provenance$protein_s4_initial$fixed$nonpositive_set_to_na), 3L)
    expect_true("3 values of 0 were read as not measured (missing). They are not used in the analysis." %in%
        vapply(notes, `[[`, "", "text"))
})

# --- T4 provenance ---------------------------------------------------------------------------------------------------

test_that("KIT-177 T4: the design step records log2(x); 0 missing and both counts (0 on a fixture without zeros)", {
    skip_if_not(workflowOwnerCapable())
    fixture <- testthat::test_path("..", "testdata", "e2e", "prot_lfq_maxquant")
    imported <- suppressMessages(importMaxQuantData(file.path(fixture, "proteinGroups.txt")))
    runs <- unique(imported$data$Run)
    design <- data.frame(Run = runs, group = sub("_[0-9]+$", "", runs), replicates = as.integer(sub("^.*_", "", runs)),
        stringsAsFactors = FALSE)
    made <- k177Design(withr::local_tempdir(), "LFQ", list(data = imported$data, mapping = imported$column_mapping), design)
    k177Defer(made$manager)
    fixed <- made$manager$.context$snapshot$step_provenance$protein_s4_initial$fixed
    expect_identical(fixed$transform, "log2(x); 0 missing")
    expect_identical(as.integer(fixed$nonpositive_set_to_na), 0L)
    expect_identical(as.integer(fixed$nonfinite_set_to_na), 0L)
    expect_null(protDesignLog2NotificationText(made$workflow$prot_design_log2_counts))
    expect_true(validateWorkflowAnalysisRecord(buildWorkflowAnalysisRecord(made$manager))$valid)
})

# --- T6 intensity-scale data -----------------------------------------------------------------------------------------

test_that("KIT-177 T6: on the five e2e fixtures the design step moves each value by log2(1 + 1/v) at most", {
    skip_if_not(workflowOwnerCapable())
    fixtures <- list(prot_lfq = c("proteinGroups.txt", "importMaxQuantData", "LFQ", 2.761e-3),
        prot_tmt = c("pd_tmt_peptides.tsv", "importProteomeDiscovererTMTData", "TMT", 2.761e-3),
        prot_lfq_fragpipe = c("combined_protein.tsv", "importFragPipeData", "LFQ", 4.303e-4),
        prot_tmt_pd = c("pd_tmt_peptides.tsv", "importProteomeDiscovererTMTData", "TMT", 2.799e-4),
        prot_lfq_maxquant = c("proteinGroups.txt", "importMaxQuantData", "LFQ", 2.492e-4))
    for (lane in names(fixtures)) {
        spec <- fixtures[[lane]]
        dir <- testthat::test_path("..", "testdata", "e2e", lane)
        imported <- suppressMessages(get(spec[[2L]])(file.path(dir, spec[[1L]])))
        long <- list(data = as.data.frame(imported$data), mapping = imported$column_mapping)
        runs <- unique(as.character(long$data[[long$mapping$run_col]]))
        sample <- sub("^.*((WT|KO)_[0-9]+)$", "\\1", runs)
        design <- data.frame(Run = runs, group = sub("_[0-9]+$", "", sample), replicates = as.integer(sub("^.*_", "", sample)),
            stringsAsFactors = FALSE)
        new <- k177Design(withr::local_tempdir(), spec[[3L]], long, design)
        k177Defer(new$manager)
        old <- k177DesignOld(withr::local_tempdir(), spec[[3L]], long, design)
        k177Defer(old$manager)
        a <- new$manager$getState("protein_s4_initial")
        b <- old$manager$getState("protein_s4_initial")
        ids <- a@protein_quant_table[[a@protein_id_column]]
        mn <- k177Matrix(a, ids, runs)
        mo <- k177Matrix(b, ids, runs)
        quantities <- long$data[[long$mapping$quantity_col]]
        bound <- log2(1 + 1 / min(quantities[!is.na(quantities) & quantities > 0]))
        expect_equal(bound, as.numeric(spec[[4L]]), tolerance = 1e-3, label = paste(lane, "bound"))
        difference <- mo - mn
        expect_lte(max(difference, na.rm = TRUE), bound + 1e-12, label = paste(lane, "largest change"))
        expect_lte(max(abs(difference - log2(1 + 2^-mn)), na.rm = TRUE), 1e-12, label = paste(lane, "each cell log2(1 + 1/v)"))
        # Changed calls at the GUI-entry lane settings (helper-gui-entry.R guiEntryProtQcSettings; cyclic loess, RUV skipped).
        qc <- list(min_reps_per_group = 1, min_groups = 1, percentile = 0.1)
        routes <- tryCatch({
            utils::capture.output(route_new <- suppressMessages(suppressWarnings(k177Route(new$manager, "cyclicloess", qc = qc,
                contrast = c("KO", "WT")))))
            utils::capture.output(route_old <- suppressMessages(suppressWarnings(k177Route(old$manager, "cyclicloess", qc = qc,
                contrast = c("KO", "WT")))))
            list(new = route_new, old = route_old)
        }, error = function(e) e)
        if (inherits(routes, "error")) {
            message(sprintf("KIT177-T6\t%s\t%d\t%.6g\t%.6g\troute stopped: %s", lane, length(ids), bound,
                max(difference, na.rm = TRUE), gsub("[\t\r\n]+", " ", conditionMessage(routes))))
            next
        }
        route_new <- routes$new
        route_old <- routes$old
        a_calls <- k177Calls(route_old$da)
        b_calls <- k177Calls(route_new$da)
        both <- merge(route_old$da, route_new$da, by = "id", suffixes = c(".old", ".new"))
        message(sprintf("KIT177-T6\t%s\t%d\t%.6g\t%.6g\t%d\t%d\t%d\t%.6g", lane, length(ids), bound,
            max(difference, na.rm = TRUE), length(a_calls), length(b_calls),
            length(union(setdiff(a_calls, b_calls), setdiff(b_calls, a_calls))),
            if (nrow(both)) max(abs(both$logFC.old - both$logFC.new), na.rm = TRUE) else NA_real_))
    }
})

# --- T7 downstream on small values -----------------------------------------------------------------------------------

test_that("KIT-177 T7: on values of 0.5 to 50 every later step runs and the counts are kept", {
    skip_if_not(workflowOwnerCapable())
    sim <- k177F177(177001L, range = c(0.5, 50))
    for (shape in c("TMT", "LFQ")) {
        long <- k177Long(sim, shape)
        made <- k177Design(withr::local_tempdir(), shape, long, sim$design)
        k177Defer(made$manager)
        design_state <- made$manager$getState("protein_s4_initial")
        q <- long$data[[long$mapping$quantity_col]]
        expect_identical(sum(!is.na(k177Matrix(design_state, sim$ids, sim$runs))), sum(!is.na(q) & q != 0), label = shape)
        missing_per_run <- colSums(is.na(k177Matrix(design_state, sim$ids, sim$runs)))
        expected <- vapply(sim$runs, function(run) sum(is.na(q[long$data$Run == run]) | q[long$data$Run == run] %in% 0),
            numeric(1))
        expect_equal(unname(missing_per_run), unname(expected), label = paste(shape, "missing per run"))
        utils::capture.output(route <- suppressMessages(suppressWarnings(k177Route(made$manager, "cyclicloess"))))
        expect_gt(nrow(route$da), 0L)
        kept <- route$correlated@protein_quant_table
        runs <- intersect(sim$runs, names(kept))
        ids <- kept[[route$correlated@protein_id_column]]
        entering <- !is.na(k177Matrix(route$correlated, ids, runs))
        leaving <- !is.na(k177Matrix(design_state, ids, runs))
        expect_identical(unname(entering), unname(leaving), label = paste(shape, "observed cells of the kept proteins"))
    }
})

# --- T8 reopened project ---------------------------------------------------------------------------------------------

k177ReopenCase <- function(type, state, fixed) {
    manager <- k177Owner(tempfile("k177-t8-"), type)
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    if (is.null(fixed)) setDIAStepProvenance(manager, "design_builder") else setDIAStepProvenance(manager, "design_builder", fixed = fixed)
    manager$saveState(state, object, list(), "design")
    manager
}

test_that("KIT-177 T8: an earlier label-free or TMT design state is told once; never a head, DIA or other state", {
    skip_if_not(workflowOwnerCapable())
    reopen <- omicGuiProtLog2Texts()[["reopen"]]
    old <- list(pivot = "protein x run", transform = "log2(x + 1)")
    cases <- list(
        list(type = "LFQ", state = "protein_s4_initial", fixed = old, expected = 1L),
        list(type = "TMT", state = "protein_s4_initial", fixed = NULL, expected = 1L),
        list(type = "LFQ", state = "protein_s4_initial", fixed = list(pivot = "protein x run", transform = "log2(x); 0 missing",
            nonpositive_set_to_na = 0L, nonfinite_set_to_na = 0L), expected = 0L),
        list(type = "DIA", state = "raw_data_s4", fixed = NULL, expected = 0L),
        list(type = "TMT", state = "protein_s4_initial", fixed = list(pivot = "protein x run", transform = "log2(x + 0.5)"),
            expected = 0L))
    for (case in cases) {
        manager <- k177ReopenCase(case$type, case$state, case$fixed)
        k177Defer(manager)
        label <- paste(case$type, case$state, if (is.null(case$fixed)) "<none>" else case$fixed$transform)
        # The design tab's start (the module server), then the load-session site's call on the same owner.
        shown <- list()
        notify <- function(ui, ..., type = "default") { shown[[length(shown) + 1L]] <<- list(text = as.character(ui), type = type) }
        testthat::local_mocked_bindings(showNotification = notify, .package = "shiny")
        testthat::local_mocked_bindings(registerProtDesignServerShells = function(...) invisible(NULL), .package = "MultiScholaR")
        shiny::testServer(mod_prot_design_server, args = list(workflow_data = shiny::reactiveValues(state_manager = manager),
            experiment_paths = list()), { session$flushReact() })
        notifyProtDesignEarlierTransform(manager)
        texts <- vapply(shown, `[[`, "", "text")
        expect_identical(sum(texts == reopen), case$expected, label = label)
        if (case$expected) expect_identical(shown[[which(texts == reopen)]]$type, "message", label = label)
    }
    # The load-session handler calls the same notifier on the restored owner.
    expect_true(any(grepl("notifyProtDesignEarlierTransform(manager)", deparse(body(da_server_load_session_handler)), fixed = TRUE)))
})

# --- T9 the scale warning (D6) ---------------------------------------------------------------------------------------

k177ScaleCase <- function(table, ids, runs, design, settings, mock_warning = FALSE) {
    root <- tempfile("k177-t9-")
    dir.create(root)
    paths <- list(source_dir = root, protein_qc_dir = file.path(root, "qc"), time_dir = file.path(root, "time"),
        publication_graphs_dir = file.path(root, "graphs"))
    for (path in paths[-1]) dir.create(path)
    notes <- list()
    logs <- character()
    normalised <- NULL
    if (mock_warning) testthat::local_mocked_bindings(warnProtNormScaleNegativeLogs = function(...) invisible(NA_integer_),
        .package = "MultiScholaR")
    shiny::testServer(function(input, output, session) {
        manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "k177-scale", "LFQ"), WorkflowState$new(audit_enabled = FALSE))
        manager$.context$reference_states <- TRUE
        manager$.context$experiment_paths <- paths
        workflow <- shiny::reactiveValues(state_manager = manager, config_list = list())
        attachDIAStateSession(manager, session, workflow)
        design$tech_rep_group <- paste(design$group, design$replicates, sep = "_")
        object <- ProteinQuantitativeData(protein_quant_table = table, protein_id_column = "Protein.Ids",
            protein_id_table = table[, "Protein.Ids", drop = FALSE], design_matrix = design, sample_id = "Run", group_id = "group",
            technical_replicate_id = "tech_rep_group", args = list())
        manager$saveState("protein_replicate_filtered", object, list(), "source")
        norm <- shiny::reactiveValues()
    }, {
        completed <- shiny::isolate(runProtNormNormalizationWorkflow(settings, workflow, norm, paths, "proteomics", "k177",
            generatePostNormalizationQcFn = function(...) NULL, generateRuvCorrectedQcFn = function(...) NULL,
            getRuvGroupingVariableFn = function() "group",
            showNotificationFn = function(message, ..., type = "default") notes[[length(notes) + 1L]] <<- list(text = message, type = type),
            handleErrorFn = function(error, ...) stop(error), messageFn = function(...) logs <<- c(logs, paste0(...))))
        expect_true(completed)
        # The managed path finishes through the same function: called again on the delivered result.
        managed <- list()
        finishProtNormNormalization(TRUE, settings, workflow, norm, function(mode) "done",
            function(message, ..., type = "default") managed[[length(managed) + 1L]] <<- list(text = message, type = type),
            function(...) NULL)
        normalised <<- list(table = withDIAStateScope(manager, function() manager$getState(manager$current_state))@protein_quant_table,
            managed = managed)
        closeDIAStateAdapter(manager)
    })
    list(notes = notes, logs = logs, table = normalised$table, managed = normalised$managed)
}

test_that("KIT-177 T9: scale on log2 values below 0 warns once with their count, and its result is unchanged", {
    skip_if_not(workflowOwnerCapable())
    sim <- k177F177(177001L)
    table <- tibble::as_tibble(cbind(data.frame(Protein.Ids = sim$ids, stringsAsFactors = FALSE),
        as.data.frame(stats::setNames(as.data.frame(log2(sim$values)), sim$runs))))
    n <- sum(log2(sim$values) < 0, na.rm = TRUE)
    expect_gt(n, 0L)
    texts <- omicGuiProtLog2Texts()
    scale <- list(norm_method = "scale", ruv_mode = "skip")
    case <- k177ScaleCase(table, sim$ids, sim$runs, sim$design, scale)
    warnings <- Filter(function(note) identical(note$type, "warning"), case$notes)
    expect_length(warnings, 1L)
    expect_identical(warnings[[1L]]$text, sprintf(texts[["scale_warning"]], n))
    expect_identical(sum(case$logs == sprintf(texts[["scale_log"]], n)), 1L)
    expect_identical(sum(vapply(case$managed, function(note) identical(note$type, "warning"), logical(1))), 1L)
    without <- k177ScaleCase(table, sim$ids, sim$runs, sim$design, scale, mock_warning = TRUE)
    expect_identical(case$table, without$table)
    loess <- k177ScaleCase(table, sim$ids, sim$runs, sim$design, list(norm_method = "cyclicloess", ruv_mode = "skip"))
    expect_false(any(vapply(loess$notes, function(note) identical(note$type, "warning"), logical(1))))
    # prot_lfq_maxquant (minimum 5,788.1): no log2 value below 0, no warning.
    imported <- suppressMessages(importMaxQuantData(testthat::test_path("..", "testdata", "e2e", "prot_lfq_maxquant",
        "proteinGroups.txt")))
    wide <- tidyr::pivot_wider(imported$data, id_cols = "Protein.Ids", names_from = "Run", values_from = "Intensity")
    runs <- setdiff(names(wide), "Protein.Ids")
    for (run in runs) wide[[run]] <- log2(wide[[run]])
    design <- data.frame(Run = runs, group = sub("_[0-9]+$", "", runs), replicates = as.integer(sub("^.*_", "", runs)),
        stringsAsFactors = FALSE)
    fixture <- k177ScaleCase(wide, wide$Protein.Ids, runs, design, scale)
    expect_false(any(vapply(fixture$notes, function(note) identical(note$type, "warning"), logical(1))))
})
