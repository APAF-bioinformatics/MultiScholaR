# KIT-255 (VAL-283): the lipidomics and metabolomics sample-correlation filter judges each sample within each assay. A
# sample below the threshold in one assay only has its values set to missing in that assay (the state's table, the
# uncorrected copy the new RUV design tests, the matrix the integration input takes for an RUV-corrected view), stays in
# the design, and leaves only when no assay keeps a value for it. Written before the fix: on the tests-only commit each E
# row fails at its stated base value; the G rows pass on both commits (ticket r2.2 section 10.2). Fixtures F172
# (helper-kit172.R) and K2 (G0's draw: B_6 passes in LCMS_Pos, fails in LCMS_Neg).

library(testthat)

# ---- fixtures ----

kit255S <- "B_6"
kit255Runs <- kit172Design()$Run

#' K2's values (`missing = FALSE` for the RUV objects).
kit255Values <- function(missing = TRUE) {
    pos <- kit172Values(255001L, missing = missing)$values
    neg <- kit172Values(255003L, missing = missing)$values
    logs <- log2(neg)
    logs <- logs[is.finite(logs)]
    m <- mean(logs)
    s <- stats::sd(logs)
    set.seed(255002L, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
    observed <- which(!is.na(neg[, kit255S]))
    neg[observed, kit255S] <- 2^stats::rnorm(length(observed), m, s)
    list(LCMS_Pos = pos, LCMS_Neg = neg)
}

kit255Slot <- function(class) if (identical(class, "lipid")) "lipid_data" else "metabolite_data"
kit255Assays <- function(object, class) methods::slot(object, kit255Slot(class))
kit255Design <- function(object) as.character(as.data.frame(object@design_matrix)[["Run"]])

#' Within-group pairs of `runs`; every pair of a sample in `failing[[assay]]` at `best`, the others at `high`.
kit255Pairs <- function(failing = list(LCMS_Neg = kit255S), best = 0.1, high = 0.9, assays = c("LCMS_Pos", "LCMS_Neg"),
    runs = kit255Runs) {
    groups <- split(runs, sub("_.*$", "", runs))
    base <- do.call(rbind, lapply(groups, function(g) {
        p <- utils::combn(g, 2L)
        data.frame(Run.x = p[1L, ], Run.y = p[2L, ], stringsAsFactors = FALSE)
    }))
    rownames(base) <- NULL
    stats::setNames(lapply(assays, function(assay) {
        low <- base$Run.x %in% failing[[assay]] | base$Run.y %in% failing[[assay]]
        out <- base
        out$pearson_correlation <- ifelse(low, best, high)
        out
    }), assays)
}

kit255Filter <- function(object, class, pairs, threshold = 0.5) {
    run <- if (identical(class, "lipid")) filterSamplesByLipidCorrelationThreshold else filterSamplesByMetaboliteCorrelationThreshold
    suppressWarnings(suppressMessages(run(object, pairs, threshold)))
}

kit255Blank <- function(values, assay, runs) {
    values[[assay]][, runs] <- NA_real_
    values
}
kit255Without <- function(values, runs) lapply(values, function(m) m[, setdiff(colnames(m), runs), drop = FALSE])
kit255Design1 <- function(runs) {
    design <- kit172Design()
    design[!(design$Run %in% runs), , drop = FALSE]
}
kit255Obj <- function(values, class, design = NULL) {
    object <- kit172Object(values, class)
    if (!is.null(design)) object@design_matrix <- design
    object
}

kit255Kernel <- function(class, object, factors = NULL) {
    run <- if (identical(class, "lipid")) runLipidsDA else runMetabolitesDA
    suppressWarnings(suppressMessages(utils::capture.output(
        value <- tryCatch(if (is.null(factors)) run(object, kit172Contrasts(), "~ 0 + group", 0.05, 0) else
            run(object, kit172Contrasts(), "~ 0 + group", 0.05, 0, ruv_factors = factors), error = function(e) e))))
    value
}

#' Whether two long tables agree on every common numeric column that is not a run column, to 1e-10.
kit255TablesEqual <- function(x, y) {
    x <- as.data.frame(x)
    y <- as.data.frame(y)
    id <- intersect(c("lipid_id", "metabolite_id"), names(x))[[1L]]
    x <- x[order(as.character(x[[id]])), , drop = FALSE]
    y <- y[order(as.character(y[[id]])), , drop = FALSE]
    if (!identical(as.character(x[[id]]), as.character(y[[id]]))) return(FALSE)
    common <- setdiff(intersect(names(x), names(y)), kit255Runs)
    numeric <- common[vapply(common, function(n) is.numeric(x[[n]]) && is.numeric(y[[n]]), logical(1))]
    all(vapply(numeric, function(n) isTRUE(all.equal(x[[n]], y[[n]], tolerance = 1e-10)), logical(1)))
}

kit255AllNa <- function(x) is.numeric(x) && all(is.na(x))

# The fixture of test-metab-s4-qc-shared.R (newMetabS4QcObject), copied: Sample_3 and Sample_4 correlate at -1 in LCMS_Pos
# and +1 in LCMS_Neg.
kit255MetabQc <- function() {
    methods::new("MetaboliteAssayData",
        metabolite_data = list(
            LCMS_Pos = tibble::tibble(feature = c("M1", "M2", "M3"), annotation_id = c("alpha", "beta", "gamma"),
                Sample_1 = c(1, 2, 3), Sample_2 = c(1, 2, 3), Sample_3 = c(3, 2, 1), Sample_4 = c(1, 2, 3)),
            LCMS_Neg = tibble::tibble(feature = c("N1", "N2", "N3"), annotation_id = c("delta", "epsilon", "zeta"),
                Sample_1 = c(5, 6, 7), Sample_2 = c(5, 6, 7), Sample_3 = c(10, 9, 8), Sample_4 = c(10, 9, 8))),
        metabolite_id_column = "feature", annotation_id_column = "annotation_id", internal_standard_regex = "^ISTD",
        design_matrix = data.frame(Run = c("Sample_1", "Sample_2", "Sample_3", "Sample_4"),
            TechRep = c("pair_a", "pair_a", "pair_b", "pair_b"), AltGroup = c("grp1", "grp1", "grp2", "grp2"), stringsAsFactors = FALSE),
        sample_id = "Run", technical_replicate_id = "TechRep", args = list(pearsonCorForSamplePairs = list()))
}

kit255Notes <- function(env = parent.frame()) {
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    testthat::local_mocked_bindings(
        showNotification = function(ui, ..., type = "default") {
            seen$notes <- c(seen$notes, paste0("[", type, "] ", paste(as.character(ui), collapse = " ")))
            invisible(NULL)
        },
        removeNotification = function(...) invisible(NULL),
        withProgress = function(expr, ...) expr,
        incProgress = function(...) invisible(NULL),
        updateSelectInput = function(...) invisible(NULL),
        .package = "shiny", .env = env)
    seen
}

# A session owner holding `object` as its design state (closed when the calling test ends).
kit255Owner <- function(object, class, env = parent.frame()) {
    fx <- kit172Quiet(kit172Owner(file.path(withr::local_tempdir(.local_envir = env), "project"), object, class))
    withr::defer(try(closeDIAStateAdapter(fx$owner), silent = TRUE), envir = env)
    fx
}

# ---- E1, E2: the kernels ----

for (class in c("lipid", "metab")) {
    test_that(paste("KIT-255 E1:", class, "a sample failing in LCMS_Neg only is set missing there and stays in the design"), {
        input <- kit172Object(kit255Values(), class)
        out <- kit255Filter(input, class, kit255Pairs())
        before <- kit255Assays(input, class)
        after <- kit255Assays(out, class)
        expect_identical(kit255Design(out), kit255Runs)
        expect_identical(after$LCMS_Pos, before$LCMS_Pos)
        expect_true(kit255AllNa(after$LCMS_Neg[[kit255S]]))
        expect_identical(names(after$LCMS_Neg), names(before$LCMS_Neg))
        other <- setdiff(names(before$LCMS_Neg), kit255S)
        expect_identical(after$LCMS_Neg[other], before$LCMS_Neg[other])
        record <- out@args$correlation_exclusions
        expect_identical(record$threshold, 0.5)
        expect_identical(names(record$excluded), "LCMS_Neg")
        expect_equal(record$excluded$LCMS_Neg, c(B_6 = 0.1))
        expect_identical(out@args$deAnalysisParameters, input@args$deAnalysisParameters)
    })
}

test_that("KIT-255 E2: metabolomics, the repository fixture keeps four samples and blanks Sample_3 and Sample_4 in LCMS_Pos", {
    input <- kit255MetabQc()
    pairs <- suppressWarnings(suppressMessages(pearsonCorForSamplePairs(input)))
    out <- kit255Filter(input, "metab", pairs)
    expect_identical(kit255Design(out), c("Sample_1", "Sample_2", "Sample_3", "Sample_4"))
    pos <- out@metabolite_data$LCMS_Pos
    expect_identical(names(pos), names(input@metabolite_data$LCMS_Pos))
    expect_true(kit255AllNa(pos$Sample_3) && kit255AllNa(pos$Sample_4))
    expect_identical(pos$Sample_1, input@metabolite_data$LCMS_Pos$Sample_1)
    expect_identical(out@metabolite_data$LCMS_Neg, input@metabolite_data$LCMS_Neg)
    expect_identical(names(out@args$correlation_exclusions$excluded), "LCMS_Pos")
})

# ---- E3: the validators ----

test_that("KIT-255 E3: the object the filter makes passes the design validators and the class validity", {
    out <- kit255Filter(kit172Object(kit255Values(), "lipid"), "lipid", kit255Pairs())
    expect_true(validateLipidDesignAlignment(out@design_matrix, out@lipid_data, sampleIdCol = "Run", groupCol = "group",
        requireReplicates = FALSE)$valid)
    expect_true(isTRUE(methods::validObject(out, test = TRUE)))
    pairs <- suppressWarnings(suppressMessages(pearsonCorForSamplePairs(kit255MetabQc())))
    metab <- kit255Filter(kit255MetabQc(), "metab", pairs)
    expect_true(validateMetabDesignAlignment(metab@design_matrix, metab@metabolite_data, sampleIdCol = "Run", groupCol = "AltGroup",
        requireReplicates = FALSE)$valid)
    expect_true(isTRUE(methods::validObject(metab, test = TRUE)))
})

# ---- E4, E5: the DA on the object, in process and the tab's preflight ----

for (class in c("lipid", "metab")) {
    test_that(paste("KIT-255 E4, E5:", class, "the DA runs on the filtered object, LCMS_Neg as without the sample; the tab's preflight passes"), {
        skip_if_not(workflowOwnerCapable())
        values <- kit255Values()
        out <- kit255Filter(kit172Object(values, class), class, kit255Pairs())
        binding <- kit172Binding(class)
        fx <- kit255Owner(out, class)
        da <- tryCatch(kit172Quiet(runWorkflowOmicsReferenceDaAnalysis(binding, fx$owner, getWorkflowStateReference(fx$owner),
            list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0), kit172Contrasts())), error = function(e) e)
        expect_false(inherits(da, "error"), info = if (inherits(da, "error")) conditionMessage(da))
        if (!inherits(da, "error")) {
            entry <- da$manifest$results[[1L]]
            expect_identical(sort(names(entry$assays)), c("LCMS_Neg", "LCMS_Pos"))
            expect_length(entry$errors, 0L)
        }
        b <- kit255Kernel(class, out)
        a <- kit255Kernel(class, kit255Obj(kit255Without(values, kit255S), class, kit255Design1(kit255S)))
        expect_false(inherits(b, "error"))
        expect_false(inherits(a, "error"))
        if (!inherits(b, "error") && !inherits(a, "error")) {
            expect_true(kit255TablesEqual(b$per_assay_results$LCMS_Neg, a$per_assay_results$LCMS_Neg))
            expect_false(kit255TablesEqual(b$per_assay_results$LCMS_Pos, a$per_assay_results$LCMS_Pos))
        }
        seen <- kit255Notes()
        testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) NULL, .package = "MultiScholaR")
        daConfig <- if (identical(class, "lipid")) omicGuiLipidDa() else omicGuiMetabDa()
        session <- shiny::MockShinySession$new()
        value <- shiny::withReactiveDomain(session, omicDaPreflight(fx$workflow,
            list(reference = getWorkflowStateReference(fx$owner)), kit172Contrasts(), list(formula_string = "~ 0 + group"), daConfig,
            binding, running = FALSE))
        expect_null(value, info = if (!is.null(value)) value$errorMessage)
    })
}

# ---- E12: the marker's earlier token reproduces the base body ----

for (class in c("lipid", "metab")) {
    test_that(paste("KIT-255 E12, E22:", class, "no option: the current rule; the earlier token: the base body (stray column, pooled design, recycled table)"), {
        input <- kit172Object(kit255Values(), class)
        out <- kit255Filter(input, class, kit255Pairs())
        expect_identical(kit255Design(out), kit255Runs)
        old <- withr::with_options(list(MultiScholaR.correlation_exclusion = "design_pooled"),
            kit255Filter(input, class, kit255Pairs()))
        expect_identical(kit255Design(old), setdiff(kit255Runs, kit255S))
        expect_true(kit255S %in% names(kit255Assays(old, class)$LCMS_Pos))
        expect_false(kit255S %in% names(kit255Assays(old, class)$LCMS_Neg))
        expect_null(old@args$correlation_exclusions)
        # E22: the pairs list holds only LCMS_Neg (named)
        only <- kit255Pairs()["LCMS_Neg"]
        new <- kit255Filter(input, class, only)
        expect_identical(kit255Design(new), kit255Runs)
        expect_identical(kit255Assays(new, class)$LCMS_Pos, kit255Assays(input, class)$LCMS_Pos)
        expect_true(kit255AllNa(kit255Assays(new, class)$LCMS_Neg[[kit255S]]))
        recycled <- withr::with_options(list(MultiScholaR.correlation_exclusion = "design_pooled"), kit255Filter(input, class, only))
        expect_identical(kit255Design(recycled), setdiff(kit255Runs, kit255S))
        expect_false(kit255S %in% names(kit255Assays(recycled, class)$LCMS_Pos))
        expect_false(kit255S %in% names(kit255Assays(recycled, class)$LCMS_Neg))
    })
}

# ---- E18, E19, E25: empty samples, untested assays ----

for (class in c("lipid", "metab")) {
    test_that(paste("KIT-255 E18:", class, "(vii) an empty sample's pairs are not counted: it stays, is not recorded, and its partner keeps its values"), {
        values <- kit255Blank(kit255Values(), "LCMS_Neg", kit255S)
        input <- kit172Object(values, class)
        pairs <- kit255Pairs(failing = list())
        neg <- pairs$LCMS_Neg
        keep <- !(neg$Run.x %in% c("B_5", "B_6") | neg$Run.y %in% c("B_5", "B_6"))
        pairs$LCMS_Neg <- rbind(neg[keep, ], data.frame(Run.x = "B_5", Run.y = "B_6", pearson_correlation = NA_real_))
        out <- kit255Filter(input, class, pairs)
        expect_identical(kit255Design(out), kit255Runs)
        expect_identical(kit255Assays(out, class)$LCMS_Neg[["B_5"]], kit255Assays(input, class)$LCMS_Neg[["B_5"]])
        expect_false(kit255S %in% names(out@args$correlation_exclusions$excluded$LCMS_Neg))
    })
    test_that(paste("KIT-255 E19:", class, "a sample failing in LCMS_Neg and untested in LCMS_Pos (an empty table) stays with its LCMS_Pos values"), {
        input <- kit172Object(kit255Values(), class)
        pairs <- kit255Pairs()
        pairs$LCMS_Pos <- pairs$LCMS_Pos[0L, ]
        out <- kit255Filter(input, class, pairs)
        expect_identical(kit255Design(out), kit255Runs)
        expect_identical(kit255Assays(out, class)$LCMS_Pos, kit255Assays(input, class)$LCMS_Pos)
        expect_true(kit255AllNa(kit255Assays(out, class)$LCMS_Neg[[kit255S]]))
    })
    test_that(paste("KIT-255 E25, G7:", class, "a sample empty in every assay leaves the design, with or without a replicate pair"), {
        values <- kit255Blank(kit255Blank(kit255Values(), "LCMS_Neg", kit255S), "LCMS_Pos", kit255S)
        input <- kit172Object(values, class)
        # E25: alone in its group (no pair in any assay)
        alone <- kit255Pairs(failing = list())
        alone <- lapply(alone, function(t) t[!(t$Run.x == kit255S | t$Run.y == kit255S), ])
        out <- kit255Filter(input, class, alone)
        expect_identical(kit255Design(out), setdiff(kit255Runs, kit255S))
        expect_false(any(vapply(kit255Assays(out, class), function(t) kit255S %in% names(t), logical(1))))
        # G7: with a replicate partner (its pairs hold NA)
        paired <- lapply(kit255Pairs(failing = list()), function(t) {
            t$pearson_correlation[t$Run.x == kit255S | t$Run.y == kit255S] <- NA_real_
            t
        })
        out <- kit255Filter(input, class, paired)
        expect_identical(kit255Design(out), setdiff(kit255Runs, kit255S))
        expect_null(out@args$correlation_exclusions)
    })
}

# ---- G rows ----

for (class in c("lipid", "metab")) {
    test_that(paste("KIT-255 G1, G2, G3, G5:", class, "failing everywhere, nowhere, one assay, twice, no cascade"), {
        input <- kit172Object(kit255Values(), class)
        # G1: a sample failing in every assay leaves the design and every assay; no args entry
        out <- kit255Filter(input, class, kit255Pairs(failing = list(LCMS_Pos = kit255S, LCMS_Neg = kit255S)))
        expect_identical(kit255Design(out), setdiff(kit255Runs, kit255S))
        expect_false(any(vapply(kit255Assays(out, class), function(t) kit255S %in% names(t), logical(1))))
        expect_identical(out@args, input@args)
        # G2: no sample failing: the object is the input
        expect_identical(kit255Filter(input, class, kit255Pairs(failing = list())), input)
        # G3: one assay
        one <- kit172Object(kit255Values()["LCMS_Neg"], class)
        out <- kit255Filter(one, class, kit255Pairs(assays = "LCMS_Neg"))
        expect_identical(kit255Design(out), setdiff(kit255Runs, kit255S))
        expect_false(kit255S %in% names(kit255Assays(out, class)$LCMS_Neg))
        expect_identical(out@args, one@args)
        # G5: applied twice, the second state is the first
        first <- kit255Filter(input, class, kit255Pairs())
        expect_identical(kit255Filter(first, class, kit255Pairs()), first)
    })

    test_that(paste("KIT-255 G6:", class, "A_1 failing in LCMS_Pos only leaves the design as it was, one press or two (no cascade)"), {
        input <- kit172Object(kit255Values(), class)
        # A_1 fails in LCMS_Pos only; A_2 passes in LCMS_Neg only through A_1: no cascade, one press or two
        cascade <- kit255Pairs(failing = list(LCMS_Pos = "A_1"))
        neg <- cascade$LCMS_Neg
        in_a <- startsWith(neg$Run.x, "A_")
        # A_1-A_2 and the pairs among A_3 to A_6 pass; every other A pair is low: A_2 passes through A_1 alone.
        keep <- paste(neg$Run.x, neg$Run.y) %in% c("A_1 A_2", "A_3 A_4", "A_3 A_5", "A_3 A_6", "A_4 A_5", "A_4 A_6", "A_5 A_6")
        neg$pearson_correlation <- ifelse(!in_a | keep, 0.99, 0.1)
        cascade$LCMS_Neg <- neg
        once <- kit255Filter(input, class, cascade)
        expect_identical(kit255Design(once), kit255Runs)
        expect_identical(kit255Design(kit255Filter(once, class, cascade)), kit255Runs)
    })

    test_that(paste("KIT-255 G4:", class, "the per-assay decision is the base body's, on random pair tables"), {
        set.seed(255004L)
        for (round in 1:12) {
            assays <- c("LCMS_Pos", "LCMS_Neg", "GCMS")[seq_len(sample(2:3, 1L))]
            values <- stats::setNames(lapply(assays, function(a) kit172Values(255100L + round, missing = FALSE, n = 45L, n_itsd = 5L)$values), assays)
            input <- kit172Object(values, class)
            t <- sample(c(0.3, 0.5, 0.7), 1L)
            pairs <- stats::setNames(lapply(assays, function(a) {
                table <- kit255Pairs(failing = list(), assays = a)[[1L]]
                table$pearson_correlation <- stats::runif(nrow(table), 0, 1)
                table
            }), assays)
            out <- kit255Filter(input, class, pairs, t)
            removed <- setdiff(kit255Runs, kit255Design(out))
            for (a in assays) {
                table <- pairs[[a]]
                everyone <- unique(c(table$Run.x, table$Run.y))
                passing <- table[table$pearson_correlation >= t, ]
                base_failing <- setdiff(everyone, unique(c(passing$Run.x, passing$Run.y)))
                recorded <- names(out@args$correlation_exclusions$excluded[[a]])
                expect_setequal(base_failing, c(recorded, removed))
            }
        }
    })
}

# ---- E23: the heatmaps drop a blank column ----

for (class in c("lipid", "metab")) {
    test_that(paste("KIT-255 E23:", class, "the DA heatmap does not draw a sample blank in the assay as a column of row means"), {
        skip_if_not(requireNamespace("ComplexHeatmap", quietly = TRUE))
        out <- kit255Filter(kit172Object(kit255Values(), class), class, kit255Pairs())
        da <- kit255Kernel(class, out)
        skip_if(inherits(da, "error"))
        heat <- if (identical(class, "lipid")) generateLipidDAHeatmap else generateMetabDAHeatmap
        run <- function(view) suppressWarnings(suppressMessages(utils::capture.output(
            value <- heat(da, selected_contrast = "groupB-groupA", selected_assay = view, top_n = 50, da_q_val_thresh = 2))))
        neg <- NULL
        utils::capture.output(neg <- heat(da, selected_contrast = "groupB-groupA", selected_assay = "LCMS_Neg", top_n = 50, da_q_val_thresh = 2))
        expect_false(is.null(neg))
        expect_false(kit255S %in% colnames(neg$plot@matrix))
        combined <- NULL
        utils::capture.output(combined <- heat(da, selected_contrast = "groupB-groupA", selected_assay = "Combined", top_n = 50, da_q_val_thresh = 2))
        expect_false(is.null(combined))
        expect_true(kit255S %in% colnames(combined$plot@matrix))
    })
}

# ---- E24: the removed factors of a view holding a blank run ----

test_that("KIT-255 E24: a corrected view with a blank run has no removed factors, and the reason names the run", {
    runs <- c("R1", "R2", "R3", "R4", "R5")
    u <- data.frame(feature_id = sprintf("F%d", 1:6), matrix(stats::rnorm(30), 6, dimnames = list(NULL, runs)), check.names = FALSE)
    v <- u
    v[2:6] <- v[2:6] + matrix(stats::rnorm(5), 6, 5, byrow = TRUE)
    u$R3 <- NA_real_
    v$R3 <- NA_real_
    out <- integrationInputRemovedFactors(u, v, 1L, runs)
    expect_identical(out$record$status, "not_determined")
    expect_identical(out$record$features, 0L)
    expect_identical(out$record$reason, "run R3 has no value in this data set, and no feature is then measured in every run")
    u$R4 <- NA_real_
    v$R4 <- NA_real_
    expect_identical(integrationInputRemovedFactors(u, v, 1L, runs)$record$reason,
        "runs R3, R4 have no value in this data set, and no feature is then measured in every run")
})

# ---- the tab's fixture, the in-process run and the managed worker journey (G0 phase 0's, G1's, copied) ----

kit255Omic <- function(class) c(lipid = "lipidomics", metab = "metabolomics")[[class]]

kit255Fixture <- function(root, values, class) {
    fixture <- kit172Owner(root, kit172Object(values, class), class)
    binding <- kit172Binding(class)
    kit172Quiet({
        runWorkflowOmicsReferenceIntensityFilter(binding, fixture$owner, 1, 0.5, fixture$workflow$config_list)
        commitWorkflowOmicsStoreStep(fixture$owner, paste0(binding$prefix, "_qc_finalize"), paste0(binding$prefix, "_qc_complete"),
            fixture$workflow$config_list, "QC complete")
    })
    fixture$class <- class
    fixture$omic <- kit255Omic(class)
    fixture$binding <- binding
    fixture$prefix <- binding$prefix
    config <- omicGuiConfig(fixture$omic)
    config$binding$binding$qc_plots <- function(...) invisible(NULL)
    fixture$config <- config
    log <- new.env(parent = emptyenv())
    log$lines <- character()
    fixture$log <- log
    fixture$addLog <- function(entry) log$lines <- c(log$lines, entry)
    fixture
}

kit255Input <- function(threshold = 0.5) {
    list(apply_itsd = FALSE, itsd_aggregation = "median", log_offset = 0, norm_method = "quantile", ruv_mode = "skip",
        auto_percentage_min = 1, auto_percentage_max = 20, ruv_grouping_variable = "group", separation_metric = "max_difference",
        k_penalty_weight = 0.5, adaptive_k_penalty = TRUE, ruv_k = 1L, ruv_percentage = 20, min_pearson_correlation_threshold = threshold)
}

kit255RunTab <- function(fixture, nd, input = kit255Input()) {
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) NULL, .package = "MultiScholaR")
    kit172Quiet(shiny::isolate(runOmicNormPipeline(input, fixture$workflow, list(), nd, fixture$addLog,
        function(...) list(color_var = "group", shape_var = NULL), fixture$config)))
}
kit255ApplyTab <- function(fixture, nd, input = kit255Input()) {
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) NULL, .package = "MultiScholaR")
    kit172Quiet(shiny::isolate(applyOmicNormCorrelation(input, fixture$workflow, nd, fixture$addLog, fixture$config)))
}

kit255Server <- function(fx, probe) {
    force(fx)
    force(probe)
    function(input, output, session) {
        workflow <- fx$workflow
        runner <- installWorkflowOmicsActions(session, workflow, fx$omic, "kit255", fx$paths)
        runner$progress <- function(message) list(set = function(...) invisible(NULL), inc = function(...) invisible(NULL),
            close = function(...) invisible(NULL))
        probe$runner <- runner
        nd <- omicNormReactiveState()
        nd$assay_names <- c("LCMS_Pos", "LCMS_Neg")
        probe$nd <- nd
        aesthetics <- function(...) list(color_var = "group", shape_var = NULL)
        config <- omicGuiConfig(fx$omic)
        shiny::observeEvent(input$run_normalization, runOmicNormPipeline(input, workflow, fx$paths, nd, fx$addLog, aesthetics, config))
        shiny::observeEvent(input$apply_correlation_filter, applyOmicNormCorrelation(input, workflow, nd, fx$addLog, config))
    }
}

kit255Settle <- function(session, probe, rounds = 40L) {
    for (i in seq_len(rounds)) {
        if (is.null(shiny::isolate(probe$runner$step$active()))) break
        session$elapse(100)
    }
    is.null(shiny::isolate(probe$runner$step$active()))
}

#' The managed journey on K2 at QC complete: the tab's run, then Apply(0.5), under the worker-journey mocks. Returns the
#' fixture (owner open) and the probe.
kit255Managed <- function(class, values, env = parent.frame()) {
    fx <- kit255Fixture(file.path(withr::local_tempdir(.local_envir = env), "project"), values, class)
    withr::defer(try(closeDIAStateAdapter(fx$owner), silent = TRUE), envir = env)
    probe <- new.env()
    omicsManagedChainMocks(fx, probe, env = env)
    suppressWarnings(suppressMessages(utils::capture.output(shiny::testServer(kit255Server(fx, probe), {
        do.call(session$setInputs, kit255Input())
        session$setInputs(run_normalization = 1L)
        session$flushReact()
        probe$settled_run <- kit255Settle(session, probe)
        session$setInputs(min_pearson_correlation_threshold = 0.5)
        session$setInputs(apply_correlation_filter = 1L)
        session$flushReact()
        probe$settled_apply <- kit255Settle(session, probe)
    }))))
    list(fx = fx, probe = probe)
}

kit255Head <- function(fx) readWorkflowOmicsState(fx$binding, fx$owner, NULL)

# ---- E6, G12: the managed worker and the tab in process make the same object, and the DA runs on it ----

for (class in c("lipid", "metab")) {
    test_that(paste("KIT-255 E6, G12:", class, "managed and in-process Apply give the same state; the DA runs on it"), {
        skip_if_not(workflowOwnerCapable())
        seen <- kit255Notes()
        values <- kit255Values()
        fx <- kit255Fixture(file.path(withr::local_tempdir(), "project"), values, class)
        withr::defer(try(closeDIAStateAdapter(fx$owner), silent = TRUE))
        nd <- omicNormReactiveState()
        kit255RunTab(fx, nd)
        kit255ApplyTab(fx, nd)
        expect_identical(fx$owner$current_state, paste0(fx$prefix, "_correlation_filtered"))
        in_process <- kit255Head(fx)
        managed <- kit255Managed(class, values)
        expect_true(isTRUE(managed$probe$settled_run) && isTRUE(managed$probe$settled_apply))
        worker <- kit255Head(managed$fx)
        slot <- kit255Slot(class)
        expect_identical(worker@design_matrix, in_process@design_matrix)
        expect_identical(methods::slot(worker, slot), methods::slot(in_process, slot))
        expect_identical(worker@args, in_process@args)
        # E6: under the new rule the sample stays, so the DA of the managed state completes for both assays
        expect_identical(kit255Design(worker), kit255Runs)
        da <- tryCatch(kit172Quiet(runWorkflowOmicsReferenceDaAnalysis(managed$fx$binding, managed$fx$owner,
            getWorkflowStateReference(managed$fx$owner), list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0),
            kit172Contrasts())), error = function(e) e)
        expect_false(inherits(da, "error"), info = if (inherits(da, "error")) conditionMessage(da))
    })
}

# ---- E15: the tab's log, notification and summary ----

for (class in c("lipid", "metab")) {
    test_that(paste("KIT-255 E15:", class, "the tab's log, notification and summary name the sample set missing, the assay, the threshold and the correlation"), {
        skip_if_not(workflowOwnerCapable())
        seen <- kit255Notes()
        fx <- kit255Fixture(file.path(withr::local_tempdir(), "project"), kit255Values(), class)
        withr::defer(try(closeDIAStateAdapter(fx$owner), silent = TRUE))
        nd <- omicNormReactiveState()
        kit255RunTab(fx, nd)
        kit255ApplyTab(fx, nd)
        text <- paste(c(fx$log$lines, seen$notes), collapse = " | ")
        expect_match(text, "B_6", fixed = TRUE)
        expect_match(text, "LCMS_Neg", fixed = TRUE)
        expect_match(text, "0.5", fixed = TRUE)
        expect_match(text, "set to missing", fixed = TRUE)
        summary <- NULL
        shiny::testServer(function(input, output, session) {
            output$summary <- omicNormCorrelationSummary(nd, fx$workflow, fx$config$binding$read_metadata)
        }, {
            summary <<- output$summary
        })
        expect_match(summary, "Removed: 0 samples", fixed = TRUE)
        expect_match(summary, "Excluded in one assay only: B_6 in LCMS_Neg", fixed = TRUE)
        # with no sample set missing the summary has no such line
        clean <- kit255Fixture(file.path(withr::local_tempdir(), "project"), kit172Values(255001L)$values, class)
        withr::defer(try(closeDIAStateAdapter(clean$owner), silent = TRUE))
        nd2 <- omicNormReactiveState()
        kit255RunTab(clean, nd2)
        kit255ApplyTab(clean, nd2)
        summary2 <- NULL
        shiny::testServer(function(input, output, session) {
            output$summary <- omicNormCorrelationSummary(nd2, clean$workflow, clean$config$binding$read_metadata)
        }, {
            summary2 <<- output$summary
        })
        expect_false(grepl("Excluded in one assay only", summary2, fixed = TRUE))
    })
}

# ---- E17: the final plots ----

for (class in c("lipid", "metab")) {
    test_that(paste("KIT-255 E17:", class, "PCA and density omit the sample in the assay where it is blank, and show it where it is not"), {
        out <- kit255Filter(kit172Object(kit255Values(), class), class, kit255Pairs())
        samples <- function(p) vapply(p, function(el) {
            d <- as.data.frame(el$data)
            col <- names(d)[vapply(d, function(x) any(as.character(x) %in% kit255Runs), logical(1))][1L]
            length(unique(as.character(d[[col]])))
        }, integer(1))
        pca <- suppressWarnings(suppressMessages(plotPca(out, grouping_variable = "group")))
        expect_equal(unname(samples(pca)), c(12L, 11L))
        density <- suppressWarnings(suppressMessages(plotDensity(out, grouping_variable = "group")))
        expect_equal(unname(samples(density)), c(12L, 11L))
    })
}

# ---- E20: a replicate set failing in one assay ----

for (class in c("lipid", "metab")) {
    test_that(paste("KIT-255 E20:", class, "a group blank in one assay is not testable there and tested in the other"), {
        values <- kit255Values()
        out <- kit255Filter(kit172Object(values, class), class, kit255Pairs(failing = list(LCMS_Neg = c("B_1", "B_2", "B_3", "B_4", "B_5", "B_6"))))
        expect_identical(kit255Design(out), kit255Runs)
        da <- kit255Kernel(class, out)
        expect_false(inherits(da, "error"), info = if (inherits(da, "error")) conditionMessage(da))
        if (!inherits(da, "error")) {
            expect_equal(sum(is.finite(da$per_assay_results$LCMS_Neg$logFC)), 0L)
            expect_gt(sum(is.finite(da$per_assay_results$LCMS_Pos$logFC)), 0L)
        }
    })
}

# ---- the RUV objects: the declared journey to the ruv_corrected head, the filter with the pairs bound (RUV removes a
# one-sample anomaly, so real pairs after RUV would not fail B_6: G1 run 45691779984) ----

kit255MockPairs <- function(env = parent.frame()) {
    testthat::local_mocked_bindings(pearsonCorForSamplePairs = function(theObject, ...) {
        runs <- as.character(methods::slot(theObject, "design_matrix")$Run)
        slot <- if (inherits(theObject, "LipidomicsAssayData")) "lipid_data" else "metabolite_data"
        assays <- names(methods::slot(theObject, slot))
        groups <- split(sort(runs), sub("_.*$", "", sort(runs)))
        rows <- do.call(rbind, lapply(groups, function(g) {
            p <- utils::combn(g, 2L)
            data.frame(Run.x = p[1L, ], Run.y = p[2L, ], stringsAsFactors = FALSE)
        }))
        rownames(rows) <- NULL
        stats::setNames(lapply(assays, function(assay) {
            out <- rows
            low <- identical(assay, "LCMS_Neg") & (out$Run.x == kit255S | out$Run.y == kit255S)
            out$pearson_correlation <- ifelse(low, -0.015, 0.95)
            out
        }), assays)
    }, .package = "MultiScholaR", .env = env)
    invisible(TRUE)
}

kit255RuvNew <- list(ruv_design = "groups_all_features", ruv_grouping_variable = "group", ruv_k = 1L, ruv_k_rule = "user", ruv_k_seed = NULL,
    ruv_control_ids = NULL)

#' An owner taken to the ruv_corrected head (QC complete, log2, quantile, RUV) on K2 without missing cells.
kit255BuiltRuv <- function(kind, class, env = parent.frame()) {
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE, .local_envir = env)
    fx <- kit172Quiet(kit172Owner(file.path(withr::local_tempdir(.local_envir = env), "project"), kit172Object(kit255Values(missing = FALSE), class), class))
    withr::defer(try(closeDIAStateAdapter(fx$owner), silent = TRUE), envir = env)
    binding <- kit172Binding(class)
    prefix <- binding$prefix
    config <- fx$workflow$config_list
    kit172Quiet({
        commitWorkflowOmicsStoreStep(fx$owner, paste0(prefix, "_qc_finalize"), paste0(prefix, "_qc_complete"), config, "QC complete")
        runWorkflowOmicsReferenceLog2(binding, fx$owner, 0, config)
        runWorkflowOmicsReferenceBetweenSampleNormalisation(binding, fx$owner, "quantile", config)
        runWorkflowOmicsReferenceRuv(binding, fx$owner, if (identical(kind, "new")) kit255RuvNew else
            utils::modifyList(kit172LegacyRuv("manual"), list(manual_k = 1L, manual_percentage = 10)), config)
    })
    list(fx = fx, binding = binding, prefix = prefix, config = config)
}

kit255Correlate <- function(built, env = parent.frame()) {
    kit255MockPairs(env)
    kit172Quiet(runWorkflowOmicsReferenceCorrelation(built$binding, built$fx$owner, "group", 0.5, built$config))
}

for (class in c("lipid", "metab")) {
    for (kind in c("new", "leg")) {
        test_that(paste("KIT-255 E7 to E10:", class, "RUV", kind, "the failing column is blank in the state, the uncorrected copy and the integration matrices"), {
            skip_if_not(workflowOwnerCapable())
            built <- kit255BuiltRuv(kind, class)
            owner <- built$fx$owner
            slot <- built$binding$data_slot
            before <- readWorkflowOmicsState(built$binding, owner, NULL)
            kit255Correlate(built)
            expect_identical(owner$current_state, paste0(built$prefix, "_correlation_filtered"))
            after <- readWorkflowOmicsState(built$binding, owner, NULL)
            expect_identical(kit255Design(after), kit255Runs)
            expect_true(kit255AllNa(methods::slot(after, slot)$LCMS_Neg[[kit255S]]))
            if (identical(kind, "new")) {
                # E7: place 2
                copy <- after@args$ruv_uncorrected_assays
                expect_true(kit255AllNa(as.numeric(as.data.frame(copy$LCMS_Neg)[[kit255S]])))
                expect_identical(copy$LCMS_Pos, before@args$ruv_uncorrected_assays$LCMS_Pos)
                other <- setdiff(names(before@args$ruv_uncorrected_assays$LCMS_Neg), kit255S)
                expect_identical(as.data.frame(copy$LCMS_Neg)[other], as.data.frame(before@args$ruv_uncorrected_assays$LCMS_Neg)[other])
                # E8: the RUV design's DA completes for both assays
                da <- tryCatch(withr::with_options(list(MultiScholaR.ruv_design_panel = TRUE), kit172Quiet(runWorkflowOmicsReferenceDaAnalysis(
                    built$binding, owner, getWorkflowStateReference(owner), list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05,
                        treat_lfc_cutoff = 0), kit172Contrasts()))), error = function(e) e)
                expect_false(inherits(da, "error"), info = if (inherits(da, "error")) conditionMessage(da))
                if (!inherits(da, "error")) {
                    entry <- da$manifest$results[[1L]]
                    expect_identical(sort(names(entry$assays)), c("LCMS_Neg", "LCMS_Pos"))
                    expect_length(entry$errors, 0L)
                }
            }
            # E9, E10: place 3 (the writer's matrices); the normalised state itself is untouched
            source <- withWorkflowStateScope(owner, function() omicsIntegrationInputSource(built$binding, owner))
            expect_identical(nrow(source$samples), 12L)
            view <- integrationInputViewName(built$binding$omic_type, "LCMS_Neg")
            expect_identical(ncol(source$matrices[[view]]), 13L)
            expect_true(kit255AllNa(source$matrices[[view]][[kit255S]]))
            expect_true(kit255AllNa(source$corrected[[view]][[kit255S]]))
            expect_false(kit255S %in% names(source$features[[integrationInputViewName(built$binding$omic_type, "LCMS_Pos")]]))
            normalised <- readWorkflowOmicsState(built$binding, owner, paste0(built$prefix, "_normalized"))
            expect_true(all(is.finite(methods::slot(normalised, slot)$LCMS_Neg[[kit255S]])))
        })
    }
}

# ---- E13, E14, T-OLD: the marker and the record; older records replay as computed ----

#' The record of `fx` replayed on a copy of its store cut to the design state; the assay inputs matched to the record's roles by
#' their SHA-256 (as KIT-254's replay). Returns the replay (`states`, `mismatches`) or its error.
kit255Replay <- function(source_dir, store_root, record, paths, binding) {
    copy <- withr::local_tempdir(.local_envir = parent.frame())
    stopifnot(system2("cp", c("-a", file.path(source_dir, "."), copy)) == 0L)
    harness <- new.env(parent = globalenv())
    sys.source(file.path(testthat::test_path("..", ".."), "tools", "profiling", "workflow_kit.R"), harness)
    api <- function(name) get(name, asNamespace("MultiScholaR"))
    copy_root <- file.path(copy, substring(store_root, nchar(source_dir) + 2L))
    prepared <- harness$kitPrepareReplayProject(copy_root, c("initial", paste0(binding$prefix, "_raw_data_s4")), api)
    sources <- list.files(file.path(source_dir, "import_sources"), full.names = TRUE)
    digests <- vapply(sources, function(path) digest::digest(file = path, algo = "sha256"), character(1))
    roles <- Filter(function(input) startsWith(input$role, binding$design_input_prefix), record$inputs)
    inputs <- c(stats::setNames(unname(sources[match(vapply(roles, function(input) input$sha256, ""), digests)]),
        vapply(roles, function(input) input$role, "")), design_file = file.path(source_dir, "design_matrix.tab"))
    tryCatch(kit172Quiet(harness$kitReplayJourney(record, prepared, store_root, inputs, paths, withr::local_tempdir(.local_envir = parent.frame()), api)),
        error = function(e) e)
}

kit255StoreRoot <- function(source_dir) dirname(dirname(list.files(source_dir, "^store[.]json$", recursive = TRUE, all.files = TRUE, full.names = TRUE)[[1L]]))

for (class in c("lipid", "metab")) {
    test_that(paste("KIT-255 E13, E14:", class, "a record made now lists the marker and the exclusions, and replays exactly"), {
        skip_if_not(workflowOwnerCapable())
        binding <- kit172Binding(class)
        prefix <- binding$prefix
        fx <- kit255Owner(kit172Object(kit255Values(), class), class)
        config <- fx$workflow$config_list
        kit172Quiet({
            commitWorkflowOmicsStoreStep(fx$owner, paste0(prefix, "_qc_finalize"), paste0(prefix, "_qc_complete"), config, "QC complete")
            runWorkflowOmicsReferenceCorrelation(binding, fx$owner, "group", 0.5, config)
        })
        json <- tempfile(fileext = ".json")
        writeWorkflowAnalysisRecord(fx$owner, json)
        record <- jsonlite::read_json(json, simplifyVector = FALSE)
        step <- Filter(function(s) identical(s$operation, paste0(prefix, "_pearson_selection")), record$steps)[[1L]]
        marker <- Filter(function(p) identical(p$name, "correlation_exclusion"), step$parameters)
        expect_length(marker, 1L)
        expect_identical(unlist(lapply(marker, function(p) p$source)), "fixed")
        expect_identical(workflowReplayRule(step, "correlation_exclusion"), "within_assay")
        # E14
        excluded <- step$outcomes$correlation_excluded_per_assay
        expect_identical(names(excluded), "LCMS_Neg")
        expect_identical(names(excluded$LCMS_Neg), kit255S)
        best <- unlist(excluded$LCMS_Neg[[kit255S]]$best_correlation)
        expect_true(is.numeric(best) && identical(round(best, 3), best))
        expect_identical(as.character(unlist(step$outcomes$correlation_removed_from_design)), character())
        # E13: it replays on head with the states identical
        replay <- kit255Replay(fx$source_dir, kit255StoreRoot(fx$source_dir), record, fx$paths, binding)
        expect_false(inherits(replay, "error"), info = if (inherits(replay, "error")) conditionMessage(replay))
        if (!inherits(replay, "error")) expect_identical(as.character(replay$mismatches), character())
    })

    test_that(paste("KIT-255 T-OLD:", class, "a project saved by the base commit replays as computed: the sample out of the design, its column in the other assay"), {
        skip_if_not(workflowOwnerCapable())
        archive <- testthat::test_path("..", "testdata", "kit-255", "base-record.tar.gz")
        expect_true(file.exists(archive))
        unpacked <- withr::local_tempdir()
        utils::untar(archive, exdir = unpacked)
        root <- file.path(unpacked, class)
        binding <- kit172Binding(class)
        source_dir <- file.path(root, "source")
        store_root <- kit255StoreRoot(source_dir)
        record <- jsonlite::read_json(file.path(root, "analysis_record.json"), simplifyVector = FALSE)
        step <- Filter(function(s) identical(s$operation, paste0(binding$prefix, "_pearson_selection")), record$steps)[[1L]]
        expect_false("correlation_exclusion" %in% vapply(step$parameters, function(p) as.character(p$name), ""))
        expect_null(step$outcomes$correlation_excluded_per_assay)
        paths <- stats::setNames(list(source_dir, file.path(root, "qc"), file.path(root, "da"), file.path(root, "pub")),
            c("source_dir", if (identical(class, "lipid")) "lipid_qc_dir" else "metabolite_qc_dir", "da_output_dir", "publication_graphs_dir"))
        replay <- kit255Replay(source_dir, store_root, record, paths, binding)
        expect_false(inherits(replay, "error"), info = if (inherits(replay, "error")) conditionMessage(replay))
        if (!inherits(replay, "error")) expect_identical(as.character(replay$mismatches), character())
        identity <- readSavedStoreIdentity(file.path(store_root, ".workflow-entities", "store.json"))
        owner <- newDIAStateAdapter(newWorkflowEntityStore(store_root, identity$project_id, identity$workflow_id), reopen = TRUE)
        withr::defer(closeDIAStateAdapter(owner))
        object <- readWorkflowOmicsState(binding, owner, NULL)
        expect_false(kit255S %in% kit255Design(object))
        expect_true(kit255S %in% names(methods::slot(object, binding$data_slot)$LCMS_Pos))
    })
}
