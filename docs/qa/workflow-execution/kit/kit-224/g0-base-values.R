# KIT-224 G0: on the base commit, the values of E1 to E6 (acceptance criterion A1; ticket section 10).
#
# Usage, from the package root on a run-tests runner (the tests-only commit: nothing under R/ differs from the base):
#   Rscript docs/qa/workflow-execution/kit/kit-224/g0-base-values.R [<out dir>]
#   <out dir>  default docs/qa/workflow-execution/kit/kit-224/g0
#
# Under pkgload::load_all(".", export_all = TRUE, helpers = TRUE). The fixtures are the test file's
# (tests/testthat/test-kit224-ruv-feature-levels.R; `kit224Mats` and `kit224Object` copied below) and KIT-172's
# (tests/testthat/helper-kit172.R). Recorded:
#   E1, E1m, E2  per assay and feature, the mean before and after the legacy ruvIII_C_Varying (k = 1, every feature a
#                potential control) over the samples finite in both, and their difference (levels.tsv);
#   DA           the legacy DA (runLipidsDA, trend and robust) on E1's corrected lipid object at fold-change cutoff 0 and
#                0.5: logFC, t, P.Value, fdr_qvalue, fdr_value_bh, AveExpr per feature (da.tsv);
#   CF           the samples the correlation filter (Pearson pairs within `group`, threshold 0.5) keeps on E1's
#                corrected lipid object;
#   E3           the legacy step lipid_ruv (Automatic) on KIT-172 seed 172101 (complete): the largest change of a feature's
#                mean between lipid_normalized and lipid_ruv_corrected; whether its record lists ruv_feature_levels;
#   E5           transcriptomics trx_edge_12: the largest |gene mean - TMM logCPM mean| of RUVIII's newY, of the export
#                file and of the integration input's corrected matrix; the manifest's correction$levels;
#   E6           KIT-172 T1b legacy_manual, lipid seed 172102, complete data, every value x 2^10 against x 1: calls changed
#                (q < 0.05, sign of logFC), max |delta t|, max |delta logFC|.
# Written: <out>/g0-base-values.tsv (id, key, value), <out>/levels.tsv, <out>/da.tsv. An error at a site is recorded as
# its value ("error: <message>"); the exit status is non-zero if any site failed.

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1L) args[[1L]] else "docs/qa/workflow-execution/kit/kit-224/g0"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
cat("KIT-224 G0 at", commit, "out", out, "\n")

suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = TRUE))
g0 <- new.env(parent = asNamespace("MultiScholaR"))

local(envir = g0, {
    quiet <- function(expr) {
        value <- NULL
        utils::capture.output(value <- suppressWarnings(suppressMessages(expr)))
        value
    }
    rows <- list()
    failed <- FALSE
    put <- function(id, key, value) {
        rows[[length(rows) + 1L]] <<- data.frame(id = id, key = key, value = as.character(value), stringsAsFactors = FALSE)
    }
    site <- function(id, expr) {
        tryCatch(expr, error = function(e) {
            failed <<- TRUE
            put(id, "error", paste("error:", conditionMessage(e)))
        })
    }

    # --- the test file's fixture (copied) ---
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
    slotOf <- function(object) if (methods::is(object, "LipidomicsAssayData")) "lipid_data" else "metabolite_data"
    valuesOf <- function(table, samples) {
        samples <- intersect(samples, names(table))
        m <- as.matrix(as.data.frame(table[, samples], check.names = FALSE))
        rownames(m) <- table[["database_identifier"]]
        m
    }
    legacy <- function(object, fixture) {
        quiet(ruvIII_C_Varying(object, ruv_grouping_variable = "group", ruv_number_k = list(pos = 1L, neg = 1L),
            ctrl = lapply(fixture$mats, rownames)))
    }

    # --- E1, E1m, E2 ---
    level_rows <- list()
    corrected_e1 <- NULL
    for (case in list(list(id = "E1", omic = "lipid", seed = 224000001L, missing = 0),
                      list(id = "E1m", omic = "lipid", seed = 224000002L, missing = 0.05),
                      list(id = "E2", omic = "metab", seed = 224000001L, missing = 0))) {
        site(case$id, {
            fixture <- kit224Mats(case$seed, case$missing)
            corrected <- legacy(kit224Object(fixture, case$omic), fixture)
            if (identical(case$id, "E1")) corrected_e1 <- corrected
            for (assay in names(fixture$mats)) {
                after <- valuesOf(methods::slot(corrected, slotOf(corrected))[[assay]], fixture$design$Sample_ID)
                before <- fixture$mats[[assay]][rownames(after), colnames(after), drop = FALSE]
                both <- is.finite(after) & is.finite(before)
                mean_in <- function(m) rowSums(replace(m, !both, 0)) / rowSums(both)
                d <- mean_in(after) - mean_in(before)
                level_rows[[length(level_rows) + 1L]] <- data.frame(id = case$id, assay = assay, feature = rownames(after),
                    mean_before = mean_in(before), mean_after = mean_in(after), difference = d, stringsAsFactors = FALSE)
                put(case$id, paste0(assay, ":max_abs_level_change"), sprintf("%.10g", max(abs(d), na.rm = TRUE)))
                put(case$id, paste0(assay, ":median_abs_level_change"), sprintf("%.10g", stats::median(abs(d), na.rm = TRUE)))
                put(case$id, paste0(assay, ":features"), nrow(after))
            }
        })
    }
    utils::write.table(do.call(rbind, level_rows), file.path(out, "levels.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

    # --- DA and the correlation filter on E1's corrected lipid object ---
    site("DA", {
        contrasts <- data.frame(contrasts = "groupB-groupA", friendly_names = "B vs A", stringsAsFactors = FALSE)
        da_rows <- list()
        for (cutoff in c(0, 0.5)) {
            long <- quiet(runLipidsDA(corrected_e1, contrasts, "~ 0 + group", 0.05, cutoff, TRUE, TRUE))$da_lipids_long
            cols <- intersect(c("lipid_id", "assay", "comparison", "logFC", "t", "P.Value", "fdr_qvalue", "fdr_value_bh",
                "AveExpr", "significant"), names(long))
            frame <- as.data.frame(long)[cols]
            frame$cutoff <- cutoff
            da_rows[[length(da_rows) + 1L]] <- frame
            put("DA", sprintf("cutoff_%s:significant", cutoff), sum(long$significant != "NS"))
        }
        utils::write.table(do.call(rbind, da_rows), file.path(out, "da.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    })
    site("CF", {
        pairs <- quiet(pearsonCorForSamplePairs(theObject = corrected_e1, correlation_group = "group"))
        filtered <- quiet(lipidomicsDeclarationBinding$filter_correlation(corrected_e1, pairs, 0.5))
        put("CF", "kept_samples", paste(as.character(filtered@design_matrix$Sample_ID), collapse = ","))
    })

    # --- E3 and E4 (iii) ---
    site("E3", {
        options(MultiScholaR.ruv_design_panel = TRUE)
        values <- kit172Values(172101L, planted = TRUE, missing = FALSE)$values
        root <- tempfile("kit224-g0-e3-")
        journey <- kit172Journey(root, kit172Object(values, "lipid"), "lipid", itsd = TRUE, method = "cyclicloess", offset = 0,
            ruv = kit172LegacyRuv("automatic"), da = FALSE)
        after <- kit172StateMatrix(journey, "lipid_ruv_corrected")
        before <- kit172StateMatrix(journey, "lipid_normalized")[rownames(after), colnames(after), drop = FALSE]
        both <- is.finite(after) & is.finite(before)
        mean_in <- function(m) rowSums(replace(m, !both, 0)) / rowSums(both)
        put("E3", "max_abs_level_change", sprintf("%.10g", max(abs(mean_in(after) - mean_in(before)), na.rm = TRUE)))
        record <- buildWorkflowAnalysisRecord(journey$owner)
        step <- Filter(function(step) identical(step$operation, "lipid_ruv"), record$steps)[[1L]]
        put("E4_iii", "record_lists_ruv_feature_levels",
            any(vapply(step$parameters, function(p) identical(p$name, "ruv_feature_levels"), logical(1))))
        closeDIAStateAdapter(journey$owner)
        unlink(root, recursive = TRUE)
    })

    # --- E5 ---
    e5 <- function() {
        reference <- trxDeclReference("trx_edge_12", env = environment())
        fixture <- trxDeclarationsOwner(tempfile("kit224-g0-e5-"), reference$dir)
        owner <- fixture$owner
        quiet({
            runTrxReferenceExpressionFilter(owner, list())
            commitWorkflowOmicsStoreStep(owner, "trx_qc_finalize", "trx_qc_complete", list(), "QC processing complete")
        })
        ruvs <- quiet(runTrxReferenceRuvs(owner, trxRuvGridParameters(), list()))
        export <- quiet(runTrxReferenceRuviiiExport(owner, "group"))
        log_cpm <- withWorkflowStateScope(owner, function() trxSavedMatrix(owner, ruvs$logcpm$value$log_cpm))
        new_y <- t(export$value$newY)
        root <- tempfile("kit224-g0-e5-save-")
        paths <- list(source_dir = file.path(root, "scripts", "transcriptomics"), integration_dir = file.path(root, "integration"),
            base_dir = root)
        dir.create(paths$source_dir, recursive = TRUE, showWarnings = FALSE)
        quiet(saveWorkflowSummaryParameters(omicGuiConfig("transcriptomics")$summary, fixture$workflow,
            projectDirs = list(transcriptomics = paths), omicType = "transcriptomics", label = "kit224", description = "d",
            notify = function(message, ...) NULL, done = function(value) NULL, failed = function(error) stop(error),
            createWorkflowArgsFn = function(...) "study_parameters.txt", writeAnalysisRecordFn = function(...) NULL,
            catFn = function(...) NULL))
        exported <- as.data.frame(vroom::vroom(file.path(paths$integration_dir, "transcriptomics_normalized_data.tsv"),
            show_col_types = FALSE))
        file_values <- t(as.matrix(exported[-1L]))
        colnames(file_values) <- exported$Run
        input <- readIntegrationInputWithMatrices(file.path(paths$integration_dir, "transcriptomics"))
        corrected <- input$views$transcriptomics$corrected_matrix
        input_values <- as.matrix(corrected[-1L])
        rownames(input_values) <- corrected$feature_id
        level <- function(values) max(abs(rowMeans(values) - rowMeans(log_cpm[rownames(values), colnames(values), drop = FALSE])))
        put("E5", "newY:max_abs_level_change", sprintf("%.10g", level(new_y)))
        put("E5", "file:max_abs_level_change", sprintf("%.10g", level(file_values)))
        put("E5", "integration_input:max_abs_level_change", sprintf("%.10g", level(input_values)))
        put("E5", "manifest_correction_levels", if (is.null(input$manifest$correction$levels)) "<absent>" else
            input$manifest$correction$levels)
        put("E5", "file_sha256", digest::digest(file = file.path(paths$integration_dir, "transcriptomics_normalized_data.tsv"),
            algo = "sha256"))
        closeDIAStateAdapter(owner)
    }
    site("E5", e5())

    # --- E6 ---
    site("E6", {
        options(MultiScholaR.ruv_design_panel = TRUE)
        values <- kit172Values(172102L, planted = TRUE, missing = FALSE)$values
        run <- function(input) {
            root <- tempfile("kit224-g0-e6-")
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
        at <- match(key(reference), key(scaled))
        put("E6", "rows", nrow(reference))
        put("E6", "call_changes", sum(calls(scaled)[at] != calls(reference), na.rm = TRUE))
        put("E6", "max_abs_dt", sprintf("%.10g", max(abs(scaled$t[at] - reference$t), na.rm = TRUE)))
        put("E6", "max_abs_dlogfc", sprintf("%.10g", max(abs(scaled$logFC[at] - reference$logFC), na.rm = TRUE)))
        put("E6", "calls_x1", sum(calls(reference) != 0))
    })

    table <- do.call(rbind, rows)
    utils::write.table(table, file.path(out, "g0-base-values.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    print(table, right = FALSE)
    if (failed) quit(status = 1L)
})
