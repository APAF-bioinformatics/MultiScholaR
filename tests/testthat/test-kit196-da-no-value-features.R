# KIT-196 (VAL-106 (1), VAL-116): a lipidomics or metabolomics assay holding a feature with no value in any fitted sample
# stopped limma's empirical-Bayes step with the trend on ("prior.weights contain NA values" when the positive residual df
# are unequal, "covariate contains NA or infinite values" when they are equal). Such a feature is removed before the
# moderated test and kept in every contrast table as a row without statistics (AveExpr missing); the other features equal
# the same DA on the data without it. Ticket section 9, T1 to T9 (T10 is in test-kit128-da-assay-failure.R).

kit196Contrasts <- function() {
    contrasts <- omicsDaContrastsTable()
    contrasts$full_format <- paste0(contrasts$friendly_names, "=", contrasts$contrasts)
    contrasts
}
kit196Settings <- function() list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0)
kit196Slot <- function(omic) if (identical(omic, "lipid")) "lipid_data" else "metabolite_data"
kit196Id <- function(omic) if (identical(omic, "lipid")) "lipid_id" else "metabolite_id"
kit196LongTable <- function(omic) if (identical(omic, "lipid")) "da_lipids_long" else "da_metabolites_long"
kit196Binding <- function(omic) if (identical(omic, "lipid")) lipidomicsDeclarationBinding else metabolomicsDeclarationBinding
kit196Statistics <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "fdr_qvalue", "fdr_value_bh", "raw_pvalue")
kit196Equal <- c("logFC", "AveExpr", "t", "P.Value", "fdr_qvalue", "fdr_value_bh", "raw_pvalue", "significant")

#' The fixture (omicsDaContrastsFixture) with, in assay `pos`: the scatter (NA at F1_2 S1, F1_3 S5 and S9, F1_4 S12),
#' the empty feature (F1_empty, NA in every sample), or every value missing.
kit196Object <- function(omic, empty = TRUE, scatter = FALSE, all_missing = FALSE, edit = NULL) {
    object <- omicsDaContrastsFixture(omic)
    slot_name <- kit196Slot(omic)
    id <- kit196Id(omic)
    tables <- methods::slot(object, slot_name)
    pos <- tables$pos
    samples <- object@design_matrix$Run
    if (scatter) {
        pos[pos[[id]] == "F1_2", "S1"] <- NA_real_
        pos[pos[[id]] == "F1_3", c("S5", "S9")] <- NA_real_
        pos[pos[[id]] == "F1_4", "S12"] <- NA_real_
    }
    if (is.function(edit)) pos <- edit(pos, id, samples)
    if (empty) {
        row <- pos[1L, , drop = FALSE]
        row[[id]] <- "F1_empty"
        row[[setdiff(names(pos), c(id, samples))[[1L]]]] <- "Name empty"
        row[, samples] <- NA_real_
        pos <- rbind(pos, row)
        rownames(pos) <- NULL
    }
    if (all_missing) pos[, samples] <- NA_real_
    tables$pos <- pos
    methods::slot(object, slot_name) <- tables
    object
}

#' The orchestrator (runLipidsDA / runMetabolitesDA) with the three contrasts: its result, its warnings and its error.
kit196Run <- function(omic, object, eBayes_trend = TRUE, eBayes_robust = TRUE, contrasts = omicsDaContrastsTable()) {
    warnings <- list()
    error <- NULL
    result <- NULL
    utils::capture.output(result <- tryCatch(withCallingHandlers(suppressMessages(omicsDaKernel(omic)(theObject = object,
        contrasts_tbl = contrasts, formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0,
        eBayes_trend = eBayes_trend, eBayes_robust = eBayes_robust)), warning = function(w) {
        warnings[[length(warnings) + 1L]] <<- w
        invokeRestart("muffleWarning")
    }), error = function(e) {
        error <<- e
        NULL
    }))
    list(result = result, warnings = warnings, error = error)
}

#' The rows of a result table, keyed and ordered (a missing table or column gives NA; the base side fails by assertion).
kit196Key <- function(table, id, columns = kit196Equal) {
    table <- as.data.frame(table %||% data.frame())
    need <- c(id, "assay", "comparison", columns)
    for (column in setdiff(need, names(table))) table[[column]] <- rep(NA, nrow(table))
    table <- table[, need, drop = FALSE]
    table <- table[order(table$assay, table$comparison, table[[id]]), , drop = FALSE]
    rownames(table) <- NULL
    table
}

#' "Completes": no failed-assay warning, no failed assay, no error.
kit196ExpectCompletes <- function(run, info = NULL) {
    expect_identical(run$result$failed_assays$pos %||% "<none>", "<none>", info = info)
    expect_null(run$result$failed_assays, info = info)
    expect_false(any(vapply(run$warnings, inherits, logical(1L), "omics_da_assay_failed")), info = info)
    expect_null(run$error, info = info)
}

#' The empty feature's rows: one per contrast, every statistic NA, AveExpr NA, significant NA; per contrast the rows with
#' AveExpr missing are the features without values (`pos` 1, `neg` 0).
kit196ExpectUntested <- function(long, id, contrasts = omicsDaContrastsTable()$contrasts, info = NULL) {
    long <- kit196Key(long, id, unique(c(kit196Statistics, "significant")))
    empty <- long[long[[id]] %in% "F1_empty", , drop = FALSE]
    expect_identical(nrow(empty), length(contrasts), info = info)
    expect_identical(sort(as.character(empty$comparison)), sort(contrasts), info = info)
    for (column in c(kit196Statistics, "significant")) expect_true(nrow(empty) > 0L && all(is.na(empty[[column]])), info = paste(info, column))
    for (contrast in contrasts) {
        for (assay in c("pos", "neg")) {
            rows <- long[long$comparison %in% contrast & long$assay %in% assay, , drop = FALSE]
            expect_identical(sum(is.na(rows$AveExpr)), if (identical(assay, "pos")) 1L else 0L, info = paste(info, contrast, assay))
            expect_identical(as.character(rows[[id]][is.na(rows$AveExpr)]), if (identical(assay, "pos")) "F1_empty" else character(),
                info = paste(info, contrast, assay))
        }
    }
}

#' The tested rows (every row but F1_empty) against the reference's: identical, or (T1, a complete reference) logFC,
#' AveExpr and significant identical and the statistics equal to a relative 1e-12, the largest relative difference printed.
kit196ExpectTestedEqual <- function(long, ref_long, id, rounding = FALSE, info = NULL) {
    got <- kit196Key(long, id)
    got <- got[!got[[id]] %in% "F1_empty", , drop = FALSE]
    want <- kit196Key(ref_long, id)
    rownames(got) <- rownames(want) <- NULL
    expect_identical(nrow(got), nrow(want), info = info)
    if (!rounding) {
        expect_identical(got, want, info = info)
        return(invisible())
    }
    expect_identical(got[c(id, "assay", "comparison", "logFC", "AveExpr", "significant")],
        want[c(id, "assay", "comparison", "logFC", "AveExpr", "significant")], info = info)
    for (column in c("t", "P.Value", "fdr_qvalue", "fdr_value_bh", "raw_pvalue")) {
        a <- as.numeric(got[[column]])
        b <- as.numeric(want[[column]])
        largest <- if (length(a) && length(a) == length(b)) max(abs(a - b) / pmax(abs(b), .Machine$double.xmin), 0, na.rm = TRUE) else NA
        cat(sprintf("KIT-196 T1 %s %s: largest relative difference %s\n", info %||% "", column, format(largest, digits = 3)))
        expect_equal(a, b, tolerance = 1e-12, info = paste(info, column))
    }
}

#' The pos matrix and design of an object (for the direct limma and contrast-half calls).
kit196Matrix <- function(omic, object) {
    pos <- methods::slot(object, kit196Slot(omic))$pos
    samples <- object@design_matrix$Run
    m <- as.matrix(pos[, samples])
    rownames(m) <- pos[[kit196Id(omic)]]
    m
}

# ---- T6 (b)'s builders, copied from test-kit129-run-path.R:6-61 (a test file; KIT-196 R06), taking the object -----------
kit196Quiet <- function(expr) {
    value <- NULL
    utils::capture.output(value <- suppressWarnings(suppressMessages(expr)))
    value
}
kit196OmicsOwner <- function(omic, object, env = parent.frame()) {
    root <- withr::local_tempdir(.local_envir = env)
    fixture <- if (identical(omic, "lipid")) lipidDeclarationsOwner(root, object = object) else metabDeclarationsOwner(root, object = object)
    withr::defer(closeDIAStateAdapter(fixture$owner), envir = env)
    fixture$object <- object
    fixture
}
kit196OmicsRuvParameters <- function(ruv_design = "groups_all_features", ruv_k = 1L, ruv_k_rule = "user", grouping = "group") {
    list(ruv_design = ruv_design, ruv_grouping_variable = grouping, ruv_k = ruv_k, ruv_k_rule = ruv_k_rule, ruv_k_seed = NULL,
        ruv_control_ids = NULL)
}

# ---- T7's managed DA tab, built as test-omic-da-module.R:756-797 builds lipid_managed (copied, not called) ----------------
kit196ManagedDa <- function(omic, object, env) {
    omic_type <- if (identical(omic, "lipid")) "lipidomics" else "metabolomics"
    prefix <- kit196Binding(omic)$prefix
    root <- normalizePath(withr::local_tempdir(pattern = paste0("kit196-", omic, "-"), .local_envir = env))
    paths <- list(source_dir = file.path(root, "source"), da_output_dir = file.path(root, "da_output"),
        publication_graphs_dir = file.path(root, "publication_graphs"))
    for (path in paths) dir.create(path, recursive = TRUE)
    bus <- list2env(list(contrasts_tbl = kit196Contrasts(), design_matrix = object@design_matrix,
        config_list = list(da = list(mode = "module-ci")), tab_status = list(setup_import = "complete",
            design_matrix = "complete", quality_control = "complete", normalization = "complete",
            differential_analysis = "pending")), parent = emptyenv())
    bus <- omicsTestOwner(object, paste0(prefix, "_correlation_filtered"), workflow = bus, omic_type = omic_type,
        experiment_paths = paths, env = env)
    probe <- new.env()
    omicsManagedChainMocks(list(owner = bus$state_manager, workflow = bus, paths = paths, source_dir = paths$source_dir), probe, env = env)
    server_fn <- if (identical(omic, "lipid")) mod_lipid_da_server else mod_metab_da_server
    ui_fn <- if (identical(omic, "lipid")) mod_lipid_da_ui else mod_metab_da_ui
    errors <- character()
    server <- function(input, output, session) {
        runner <- installWorkflowOmicsActions(session, bus, omic_type, "kit196", paths)
        runner$progress <- function(message) list(inc = function(...) invisible(NULL), close = function() invisible(NULL))
        session$userData$da <- server_fn("de", bus, paths, omic_type, "kit196")
    }
    suppressWarnings(suppressMessages(utils::capture.output(shiny::testServer(server, {
        defaults <- omicsUiInputDefaults(ui_fn("x"))
        names(defaults) <- sub("^x-", "de-", names(defaults))
        tryCatch({
            do.call(session$setInputs, defaults)
            session$setInputs(`de-formula_string` = "~ 0 + group", `de-da_q_val_thresh` = 0.05, `de-treat_lfc_cutoff` = 0)
            session$setInputs(`de-run_da_analysis` = 1L)
            session$flushReact()
        }, error = function(e) errors <<- c(errors, conditionMessage(e)))
    }))))
    list(manager = bus$state_manager, probe = probe, errors = errors)
}

for (omic in c("lipid", "metab")) {
    id <- kit196Id(omic)
    long_table <- kit196LongTable(omic)

    test_that(paste("KIT-196 T1:", omic, "a feature without values (equal df) is not tested and the others equal the DA without it"), {
        run <- kit196Run(omic, kit196Object(omic))
        ref <- kit196Run(omic, kit196Object(omic, empty = FALSE))
        kit196ExpectCompletes(run, "T1")
        expect_identical(run$result$features_without_values, list(pos = "F1_empty"))
        expect_null(ref$result$features_without_values)
        kit196ExpectUntested(run$result$per_assay_results$pos %||% data.frame(), id, info = "T1 pos")
        long <- run$result[[long_table]] %||% data.frame()
        kit196ExpectUntested(long, id, info = "T1 long")
        kit196ExpectTestedEqual(run$result$per_assay_results$pos, ref$result$per_assay_results$pos, id, rounding = TRUE, info = "T1 pos")
        expect_identical(run$result$per_assay_results$neg, ref$result$per_assay_results$neg)
    })

    test_that(paste("KIT-196 T2:", omic, "a feature without values (unequal df, KIT-188's message) is not tested and the others equal the DA without it"), {
        run <- kit196Run(omic, kit196Object(omic, scatter = TRUE))
        ref <- kit196Run(omic, kit196Object(omic, empty = FALSE, scatter = TRUE))
        kit196ExpectCompletes(run, "T2")
        expect_identical(run$result$features_without_values, list(pos = "F1_empty"))
        expect_null(ref$result$features_without_values)
        kit196ExpectUntested(run$result$per_assay_results$pos %||% data.frame(), id, info = "T2 pos")
        kit196ExpectUntested(run$result[[long_table]] %||% data.frame(), id, info = "T2 long")
        kit196ExpectTestedEqual(run$result$per_assay_results$pos, ref$result$per_assay_results$pos, id, info = "T2 pos")
        expect_identical(run$result$per_assay_results$neg, ref$result$per_assay_results$neg)
    })

    test_that(paste("KIT-196 T3:", omic, "limma's mechanism: an NA trend covariate stops eBayes with the df case's message"), {
        cat(sprintf("KIT-196 T3: limma %s\n", as.character(utils::packageVersion("limma"))))
        design_table <- kit196Object(omic)@design_matrix
        design <- stats::model.matrix(~ 0 + group, data = design_table)
        rownames(design) <- design_table$Run
        contrasts <- limma::makeContrasts(contrasts = omicsDaContrastsTable()$contrasts, levels = colnames(design))
        moderate <- function(m, trend = TRUE) tryCatch(limma::eBayes(limma::contrasts.fit(limma::lmFit(m, design), contrasts),
            trend = trend, robust = TRUE), error = function(e) e)
        unequal <- kit196Matrix(omic, kit196Object(omic, scatter = TRUE))
        fit <- limma::lmFit(unequal, design)
        expect_true(is.nan(fit$Amean[["F1_empty"]]))
        positive <- fit$df.residual[fit$df.residual > 0]
        expect_false(identical(min(positive), max(positive)))
        e <- moderate(unequal)
        expect_s3_class(e, "error")
        expect_identical(conditionMessage(e), "prior.weights contain NA values")
        expect_identical(deparse(conditionCall(e)[[1L]]), "fitFDistUnequalDF1")
        equal <- kit196Matrix(omic, kit196Object(omic))
        e <- moderate(equal)
        expect_s3_class(e, "error")
        expect_identical(conditionMessage(e), "covariate contains NA or infinite values")
        for (m in list(unequal, equal)) {
            expect_s4_class(moderate(m, trend = FALSE), "MArrayLM")
            expect_s4_class(moderate(m[rownames(m) != "F1_empty", , drop = FALSE]), "MArrayLM")
        }
    })

    test_that(paste("KIT-196 T4:", omic, "guards: a feature without a value in one group, and one with zero residual df, are tested"), {
        no_group <- kit196Object(omic, empty = FALSE, scatter = TRUE, edit = function(pos, id, samples) {
            pos[pos[[id]] == "F1_5", samples[9:12]] <- NA_real_
            pos
        })
        run <- kit196Run(omic, no_group)
        kit196ExpectCompletes(run, "T4 (a)")
        f5 <- kit196Key(run$result$per_assay_results$pos, id)
        f5 <- f5[f5[[id]] %in% "F1_5", , drop = FALSE]
        expect_identical(nrow(f5), 3L)
        expect_true(all(is.na(f5$logFC[f5$comparison %in% c("groupC-groupA", "groupC-groupB")])))
        expect_true(all(is.finite(f5$logFC[f5$comparison %in% "groupB-groupA"])))
        zero_df <- kit196Object(omic, empty = FALSE, scatter = TRUE, edit = function(pos, id, samples) {
            pos[pos[[id]] == "F1_6", samples[-c(1L, 5L, 9L)]] <- NA_real_
            pos
        })
        run <- kit196Run(omic, zero_df)
        kit196ExpectCompletes(run, "T4 (b)")
        f6 <- kit196Key(run$result$per_assay_results$pos, id)
        f6 <- f6[f6[[id]] %in% "F1_6", , drop = FALSE]
        expect_identical(nrow(f6), 3L)
        expect_true(all(is.finite(f6$logFC)) && all(is.finite(f6$P.Value)))
    })

    test_that(paste("KIT-196 T5:", omic, "an assay in which no feature has a value is refused with a classed message; the other assay is analysed"), {
        object <- kit196Object(omic, empty = FALSE, all_missing = TRUE)
        run <- kit196Run(omic, object)
        expect_identical(run$result$failed_assays$pos %||% "<none>", "no feature has a value in any sample, so nothing can be tested")
        expect_identical(names(run$result$per_assay_results), "neg")
        expect_true(any(vapply(run$warnings, inherits, logical(1L), "omics_da_assay_failed")))
        fit_half <- if (identical(omic, "lipid")) fitLipidDaAssayModel else fitMetabDaAssayModel
        contrast_half <- if (identical(omic, "lipid")) applyLipidDaContrasts else applyMetabDaContrasts
        e <- tryCatch(suppressMessages(contrast_half(fit_half(kit196Matrix(omic, object), object@design_matrix, "~ 0 + group", "Run"),
            omicsDaContrastsTable()$contrasts, NA, TRUE, TRUE)), error = function(e) e)
        expect_s3_class(e, "omics_da_no_observed_features")
    })

    # The DA tab's paths: their results, the manifest's errors per contrast and the saved contrast fits.
    kit196ExpectTab <- function(errors, results, ref_results, label) {
        contrasts <- omicsDaContrastsTable()$contrasts
        expect_identical(errors, stats::setNames(rep(list(list()), length(contrasts)), contrasts), info = label)
        expect_length(results$failed_contrasts, 0L)
        long <- results[[long_table]] %||% data.frame()
        kit196ExpectUntested(long, id, info = label)
        kit196ExpectTestedEqual(long, ref_results[[long_table]], id, info = label)
        for (contrast in contrasts) {
            rows <- long[long$comparison %in% contrast & long$assay %in% "pos", , drop = FALSE]
            n <- sum(is.na(rows$logFC))
            expect_true("F1_empty" %in% rows[[id]][is.na(rows$logFC)], info = paste(label, contrast))
            text <- if (nrow(long)) omicDaSummaryText(long, contrast, "pos", 0.05, list(noun = "x", significance_label = "Q",
                not_tested = TRUE)) else ""
            expect_true(grepl(sprintf("Not tested (no value in a compared group): %d", n), text, fixed = TRUE), info = paste(label, contrast))
        }
        fitted <- rownames(results$contrasts_results$pos$fit.eb$coefficients)
        expected <- setdiff(as.character(methods::slot(kit196Object(omic, scatter = TRUE), kit196Slot(omic))$pos[[id]]), "F1_empty")
        expect_identical(sort(fitted %||% character()), sort(expected), info = label)
    }
    kit196Errors <- function(manifest) lapply(manifest$results, function(entry) entry$errors)

    test_that(paste("KIT-196 T6:", omic, "the DA tab path in process completes and reports the feature without values as not tested"), {
        test_env <- environment()
        tab <- function(object) {
            workflow <- omicsTestOwner(object, "da_input", env = test_env)
            manager <- workflow$state_manager
            wrapper <- if (identical(omic, "lipid")) runLipidReferenceDaAnalysis else runMetabReferenceDaAnalysis
            out <- NULL
            suppressWarnings(suppressMessages(utils::capture.output(out <- tryCatch(wrapper(manager, getWorkflowStateReference(manager),
                kit196Settings(), kit196Contrasts()), error = function(e) list(error = conditionMessage(e))))))
            out
        }
        out <- tab(kit196Object(omic, scatter = TRUE))
        ref <- tab(kit196Object(omic, empty = FALSE, scatter = TRUE))
        expect_null(out$error)
        kit196ExpectTab(kit196Errors(out$manifest), out$results, ref$results, "T6")
    })

    test_that(paste("KIT-196 T6 (b):", omic, "the declared DA without RUV factors (design none) completes and reports the feature as not tested"), {
        withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
        test_env <- environment()
        binding <- kit196Binding(omic)
        declared <- function(object) {
            fixture <- kit196OmicsOwner(omic, object, env = test_env)
            none <- kit196OmicsRuvParameters("none", NULL, NULL)
            ruv <- kit196Quiet(runWorkflowOmicsReferenceRuv(binding, fixture$owner, none, list()))
            # The precondition (an instrument check, not T6 (b)'s value): the RUV step ran as the option None.
            expect_identical(ruv$value$ruv_design, "none")
            expect_identical(unname(unlist(ruv$value$best_k)), c(0L, 0L))
            da <- kit196Quiet(tryCatch(runWorkflowOmicsReferenceDaAnalysis(binding, fixture$owner, ruv$reference, kit196Settings(),
                kit196Contrasts()), error = function(e) list(error = conditionMessage(e))))
            expect_false(is.null(da$manifest$results[[1L]]$model$ruv_da_method))
            da
        }
        out <- declared(kit196Object(omic, scatter = TRUE))
        ref <- declared(kit196Object(omic, empty = FALSE, scatter = TRUE))
        expect_null(out$error)
        kit196ExpectTab(kit196Errors(out$manifest), out$results, ref$results, "T6 (b)")
    })

    test_that(paste("KIT-196 T7:", omic, "the DA tab path managed (worker) completes and reports the feature as not tested"), {
        managed <- function(object) {
            # The owner and the worker's mocks live for this call only.
            run <- kit196ManagedDa(omic, object, environment())
            manager <- run$manager
            reference <- manager$.context$snapshot$da_results
            manifest <- if (!is.null(reference)) hydrateSavedWorkflowEntity(manager$.context$store, reference)$value
            list(manifest = manifest, results = if (!is.null(manifest)) readWorkflowOmicsDaResults(kit196Binding(omic), manager, manifest),
                operations = run$probe$operations, errors = run$errors)
        }
        out <- managed(kit196Object(omic, scatter = TRUE))
        ref <- managed(kit196Object(omic, empty = FALSE, scatter = TRUE))
        # The precondition (an instrument check): the DA ran on the worker.
        expect_true(paste0(kit196Binding(omic)$prefix, "_da_contrast") %in% ref$operations)
        expect_false(is.null(out$manifest))
        kit196ExpectTab(kit196Errors(out$manifest), out$results, ref$results, "T7")
    })
}

# ---- T8: seeded sparse data modelled on KIT-188's M12 cell (sim-k188.R:62, :121-137, :251-265) -------------------------
kit196Sparse <- function(seed) {
    withr::with_seed(seed, .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection", {
        groups <- LETTERS[1:4]
        group <- rep(groups, each = 3L)
        n <- 150L
        mu <- stats::rnorm(n, 2.2, 3)
        x <- 10 / stats::rchisq(n, 10)
        s2 <- 0.25 * x / stats::median(x)
        n_eff <- round(0.2 * n)
        where <- sample.int(n, n_eff)
        signs <- sample(rep(c(1, -1), c(round(n_eff / 2), n_eff - round(n_eff / 2))))
        target <- sample(2:4, n_eff, replace = TRUE)
        beta <- matrix(0, n, 4L)
        beta[cbind(where, target)] <- signs
        y <- mu + beta[, match(group, groups), drop = FALSE] + matrix(stats::rnorm(n * 12L), n, 12L) * sqrt(s2)
        lod <- stats::uniroot(function(l) mean(stats::plogis(l - y)) - 0.20, c(min(y) - 60, max(y) + 60), tol = 1e-10)$root
        missing <- matrix(stats::runif(n * 12L), n, 12L) < stats::plogis(lod - y)
    })
    y[missing] <- NA_real_
    dimnames(y) <- list(sprintf("F%04d", seq_len(n)), sprintf("R%03d", seq_len(12L)))
    list(y = y, design = data.frame(Run = colnames(y), group = group, stringsAsFactors = FALSE))
}
#' The three preconditions, from the matrix and limma::lmFit only (R04).
kit196SparseQualifies <- function(sim) {
    without <- rowSums(!is.na(sim$y)) == 0L
    design <- stats::model.matrix(~ 0 + group, data = sim$design)
    df <- limma::lmFit(sim$y, design)$df.residual
    positive <- df[df > 0]
    any(without) && length(positive) > 0L && !identical(min(positive), max(positive)) && anyNA(sim$y[!without, , drop = FALSE])
}
kit196SparseObject <- function(sim, keep = rownames(sim$y)) {
    y <- sim$y[keep, , drop = FALSE]
    table <- data.frame(metabolite_id = rownames(y), metabolite_name = rownames(y), stringsAsFactors = FALSE)
    table <- cbind(table, as.data.frame(y, check.names = FALSE))
    methods::new("MetaboliteAssayData", metabolite_data = list(LCMS_Pos = table), design_matrix = sim$design, sample_id = "Run",
        group_id = "group", metabolite_id_column = "metabolite_id", annotation_id_column = "metabolite_name")
}

test_that("KIT-196 T8: seeded sparse data (KIT-188's M12 cell, VD 20 %) complete; the tested rows equal the DA without the features without values", {
    seed <- NULL
    for (candidate in 196001:196010) {
        if (kit196SparseQualifies(kit196Sparse(candidate))) {
            seed <- candidate
            break
        }
    }
    expect_true(!is.null(seed))
    skip_if(is.null(seed), "no qualifying seed (a result, recorded above)")
    sim <- kit196Sparse(seed)
    without <- rownames(sim$y)[rowSums(!is.na(sim$y)) == 0L]
    cat(sprintf("KIT-196 T8: seed %d, features without values %d (%s)\n", seed, length(without), paste(without, collapse = ",")))
    contrasts <- omicsDaContrastsTable(c("groupB-groupA", "groupC-groupA", "groupD-groupA"))
    run <- kit196Run("metab", kit196SparseObject(sim), contrasts = contrasts)
    message <- if (!is.null(run$error)) conditionMessage(run$error)
    expect_identical(message %||% "<none>", "<none>")
    expect_identical(run$result$features_without_values, list(LCMS_Pos = without))
    ref <- kit196Run("metab", kit196SparseObject(sim, setdiff(rownames(sim$y), without)), contrasts = contrasts)
    got <- kit196Key(run$result$da_metabolites_long, "metabolite_id")
    got <- got[!got$metabolite_id %in% without, , drop = FALSE]
    rownames(got) <- NULL
    expect_identical(got, kit196Key(ref$result$da_metabolites_long, "metabolite_id"))
})

test_that("KIT-196 T9: with the trend off, the tested rows equal the DA without the feature without values (D2)", {
    for (omic in c("lipid", "metab")) {
        id <- kit196Id(omic)
        for (robust in c(FALSE, TRUE)) {
            run <- kit196Run(omic, kit196Object(omic, scatter = TRUE), eBayes_trend = FALSE, eBayes_robust = robust)
            ref <- kit196Run(omic, kit196Object(omic, empty = FALSE, scatter = TRUE), eBayes_trend = FALSE, eBayes_robust = robust)
            label <- sprintf("T9 (%s) %s", if (robust) "b" else "a", omic)
            got <- kit196Key(run$result$per_assay_results$pos, id)
            got <- got[!got[[id]] %in% "F1_empty", , drop = FALSE]
            rownames(got) <- NULL
            want <- kit196Key(ref$result$per_assay_results$pos, id)
            largest <- if (nrow(got) == nrow(want)) max(abs(as.numeric(got$P.Value) - as.numeric(want$P.Value)), 0, na.rm = TRUE) else NA
            cat(sprintf("KIT-196 %s: tested rows equal the reference %s; largest absolute difference in P.Value %s\n", label,
                identical(got, want), format(largest, digits = 3)))
            expect_identical(got, want, info = label)
        }
    }
})
