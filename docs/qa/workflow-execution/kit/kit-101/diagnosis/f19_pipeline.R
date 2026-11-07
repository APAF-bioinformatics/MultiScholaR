# KIT-101 revision 2.3, F-19 diagnosis: which step of the DIA_LIMPA route makes the false positives.
# Diagnosis only: no package code is changed. Runs in one R process (load_all, the kit_100 library on R_LIBS):
#   Phase 1: the entry-suite lane `prot_dia_limpa` (tests/testthat/helper-gui-entry.R, mod_proteomics_server under
#            shiny::testServer, the production step functions) from import to its DA, then a dump step that reads every
#            committed state through the owner (manager$getState) and copies the DA tables the GUI wrote.
#            The same lane driver on the DIA route (`prot_dia` steps, use_limpa FALSE) on the prot_dia_limpa fixture
#            (variant e).
#   Phase 2: DA variants on the dumped states, each through the DA kernel the GUI's DA declaration calls
#            (`differentialAbundanceAnalysis`, diaContrastOperation's arguments): RUV none / k = 1 / k = 2 (the RUV
#            declaration's sweep and computeRuvCovariates), RUV-III correction mode, limpa DPC matrix with plain limma,
#            a plain peptide mean summarisation with limma, and the same with a per-sample median centring.
# Usage (repository root): R_LIBS=<kit_100 library> Rscript docs/qa/workflow-execution/kit/kit-101/diagnosis/f19_pipeline.R
# Optional: F19_OUT=<dir> (default tempdir), F19_PHASE=1|2|both (default both; 2 reads F19_OUT's dumps).

suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = FALSE))
OUT <- Sys.getenv("F19_OUT", file.path(tempdir(), "f19"))
PHASE <- Sys.getenv("F19_PHASE", "both")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
cat("output:", OUT, "\n")
planted <- sprintf("P%05d", 1:40)
pid <- function(x) { m <- regmatches(x, regexpr("P[0-9]{5}", x)); out <- rep(NA_character_, length(x)); out[grepl("P[0-9]{5}", x)] <- m; out }

# ------------------------------------------------------------------------------------------------------------------
# Phase 1: the lanes.
# ------------------------------------------------------------------------------------------------------------------
if (PHASE %in% c("1", "both")) {
    suppressMessages(testthat::source_test_helpers("tests/testthat", env = globalenv()))
    original_steps <- guiEntrySteps
    original_lane <- guiEntryLane
    dumpStep <- function(tag) guiEntryStep("diag_dump", "da", function(api) {
        wd <- api$workflow_data
        manager <- shiny::isolate(wd$state_manager)
        history <- unique(guiEntryHistory(manager))
        states <- list()
        for (name in history) states[[name]] <- tryCatch(withWorkflowStateScope(manager, function() manager$getState(name)),
            error = function(e) paste("<error>", conditionMessage(e)))
        saveRDS(list(history = history, head = manager$current_state, workflow_type = manager$workflow_type,
            states = states, contrasts_tbl = shiny::isolate(wd$contrasts_tbl), config_list = shiny::isolate(wd$config_list)),
            file.path(OUT, paste0(tag, "_states.rds")))
        files <- list.files(api$inputs$paths$base_dir, recursive = TRUE, full.names = TRUE)
        keep <- files[grepl("\\.(tsv|tab|csv)$", files) & grepl("de_|da_|differential|ruv", basename(files), ignore.case = TRUE)]
        target <- file.path(OUT, paste0(tag, "_gui_files")); dir.create(target, showWarnings = FALSE)
        file.copy(keep, target, overwrite = TRUE)
        writeLines(sub(api$inputs$paths$base_dir, "<base>", keep, fixed = TRUE), file.path(target, "index.txt"))
    })
    guiEntrySteps <- function(lane) {
        steps <- original_steps(lane)
        ids <- vapply(steps, `[[`, "", "id")
        c(steps[seq_len(which(ids == "da"))], list(dumpStep(lane$diag_tag)))
    }
    guiEntryLane <- function(lane_id) {
        if (identical(lane_id, "diag_dia_on_limpa_fixture")) {
            lane <- original_lane("prot_dia")
            lane$lane_id <- lane_id
            lane$fixture_dir <- "prot_dia_limpa"
            lane$diag_tag <- "dia_route"
            return(lane)
        }
        lane <- original_lane(lane_id)
        lane$diag_tag <- "limpa_route"
        lane
    }
    for (id in c("prot_dia_limpa", "diag_dia_on_limpa_fixture")) {
        t0 <- Sys.time()
        record <- guiEntryRunLane(id)
        cat(sprintf("lane %s: %d steps, %.0f s; heads: %s\n", id, length(record$steps), record$seconds,
            paste(unique(vapply(record$steps, `[[`, "", "head")), collapse = " > ")))
        notes <- unlist(lapply(record$steps, function(s) vapply(s$notifications, function(n) paste0("[", s$step, "/", n$type, "] ", substr(n$message, 1, 160)), "")))
        cat(paste(" ", grep("error|warn", notes, value = TRUE, ignore.case = TRUE), collapse = "\n"), "\n")
    }
}

# ------------------------------------------------------------------------------------------------------------------
# Phase 2: the variants.
# ------------------------------------------------------------------------------------------------------------------
if (PHASE %in% c("2", "both")) {
    limpa <- readRDS(file.path(OUT, "limpa_route_states.rds"))
    dia <- readRDS(file.path(OUT, "dia_route_states.rds"))
    cat("\nDIA_LIMPA history:", paste(limpa$history, collapse = " > "), "| head", limpa$head, "\n")
    cat("DIA history:      ", paste(dia$history, collapse = " > "), "| head", dia$head, "\n")
    contrast <- as.character(limpa$contrasts_tbl[[1L]][[1L]])
    cat("contrast:", contrast, "\n")

    # The DA declaration's kernel (diaContrastOperation), outside the store.
    daCalls <- function(object, formula = "~ 0 + group") {
        single <- data.frame(contrasts = contrast, stringsAsFactors = FALSE)
        arguments <- list(contrasts_tbl = single, formula_string = formula, da_q_val_thresh = 0.05, treat_lfc_cutoff = 0,
            eBayes_trend = TRUE, eBayes_robust = TRUE, args_group_pattern = "(\\d+)", args_row_id = object@protein_id_column,
            qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue")
        for (section in c("differentialAbundanceAnalysis", "differentialAbundanceAnalysisHelper"))
            if (is.null(object@args[[section]])) object@args[[section]] <- arguments
        result <- NULL
        invisible(utils::capture.output(result <- suppressMessages(suppressWarnings(differentialAbundanceAnalysis(object,
            contrasts_tbl = single, formula_string = formula, da_q_val_thresh = 0.05, treat_lfc_cutoff = 0,
            qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue")))))
        tab <- result$contrasts_results_table
        list(table = data.frame(protein = pid(tab[[object@protein_id_column]]), logFC = tab$logFC, p = tab$raw_pvalue,
            q = tab$fdr_qvalue, bh = tab$fdr_value_bh_adjustment, stringsAsFactors = FALSE), model = result$da_model)
    }
    rows <- list()
    summariseCalls <- function(variant, calls, note = "") {
        tab <- calls$table
        sig <- tab$protein[!is.na(tab$q) & tab$q < 0.05]
        sig_bh <- tab$protein[!is.na(tab$bh) & tab$bh < 0.05]
        row <- data.frame(variant = variant, model = calls$model$model %||% "?", formula = calls$model$formula_used %||% "?",
            proteins = nrow(tab), sig = length(sig), planted_recall = sum(sig %in% planted), false_pos = sum(!sig %in% planted),
            sig_bh = length(sig_bh), false_pos_bh = sum(!sig_bh %in% planted),
            median_logFC_nonplanted = round(stats::median(tab$logFC[!tab$protein %in% planted]), 3),
            p_lt_0.05_nonplanted = sum(tab$p[!tab$protein %in% planted] < 0.05, na.rm = TRUE),
            note = note, stringsAsFactors = FALSE)
        rows[[length(rows) + 1L]] <<- row
        cat(sprintf("%-58s %-6s sig %3d planted %2d FP %2d | BH sig %3d FP %2d | median logFC null %6.3f | null p<0.05 %3d  %s\n",
            variant, row$model, row$sig, row$planted_recall, row$false_pos, row$sig_bh, row$false_pos_bh,
            row$median_logFC_nonplanted, row$p_lt_0.05_nonplanted, note))
        invisible(row)
    }

    states <- limpa$states
    ruv_state <- states[["ruv_covariates"]]
    pre_name <- rev(setdiff(limpa$history, c("ruv_covariates", "ruv_corrected")))[[1L]]
    pre <- states[[pre_name]]
    cat("pre-RUV state:", pre_name, "| class", class(pre), "| proteins", nrow(pre@protein_quant_table), "\n")
    cat("ruv_covariates design columns:", paste(names(ruv_state@design_matrix), collapse = ", "), "\n")
    print(ruv_state@args$ruvCovariates[c("k", "percentage_as_neg_ctrl", "control_count", "replicate_structure", "block_reason")])
    elist <- pre@args$limpa_dpc_quant_results$quantified_elist
    cat("EList present:", !is.null(elist), "| dim", paste(dim(elist$E), collapse = " x "), "\n")
    group <- setNames(as.character(pre@design_matrix[[pre@group_id]]), pre@design_matrix[[pre@sample_id]])
    med <- apply(elist$E, 2, stats::median)
    cat("per-sample median of the DPC-quant EList E (log2):\n"); print(round(med, 3))
    cat(sprintf("mean KO - mean WT of the per-sample medians: %.3f\n", mean(med[group[names(med)] == "KO"]) - mean(med[group[names(med)] == "WT"])))

    cat("\nVariant table (q = fdr_qvalue as the DA writes it; BH = fdr_value_bh_adjustment):\n")
    # GUI-equivalent: the committed ruv_covariates state as the DA loads it.
    summariseCalls("GUI state ruv_covariates (k=2 covariates) through DA kernel", daCalls(ruv_state))
    gui_files <- list.files(file.path(OUT, "limpa_route_gui_files"), pattern = "\\.(tsv|tab|csv)$", full.names = TRUE)
    for (f in gui_files) {
        t <- tryCatch(utils::read.delim(f, stringsAsFactors = FALSE), error = function(e) NULL)
        if (is.null(t) || !"fdr_qvalue" %in% names(t)) next
        idcol <- intersect(c("Protein.Ids", "uniprot_acc", "Protein.Group", "protein_id"), names(t))[1L]
        if (is.na(idcol)) idcol <- names(t)[[1L]]
        sig <- pid(t[[idcol]])[t$fdr_qvalue < 0.05]
        cat(sprintf("  GUI file %s: rows %d, q<0.05 %d, planted %d, FP %d\n", basename(f), nrow(t), length(sig), sum(sig %in% planted), sum(!sig %in% planted)))
    }

    # (a) no RUV.
    summariseCalls("(a) no RUV: dpcDE ~0+group", daCalls(pre))

    # The RUV declaration's sweep on the pre-RUV state, then k = 1 and k = 2 with its controls.
    parameters <- resolveDIARuvCovariatesParameters(limpa$config_list, "group", "covariates")
    sweep <- suppressMessages(findBestNegCtrlPercentageCapped(pre, percentage_range = seq(parameters$percentage_min, parameters$percentage_max, by = 1),
        num_components_to_impute = parameters$num_components_to_impute, ruv_grouping_variable = parameters$ruv_grouping_variable,
        ruv_qval_cutoff = parameters$ruv_qval_cutoff, ruv_fdr_method = parameters$ruv_fdr_method, separation_metric = parameters$separation_metric,
        k_penalty_weight = parameters$k_penalty_weight, adaptive_k_penalty = isTRUE(parameters$adaptive_k_penalty), verbose = FALSE))
    controls <- sweep$best_control_genes_index
    ctl_ids <- pid(pre@protein_quant_table[[pre@protein_id_column]])[controls]
    cat(sprintf("\nRUV sweep: best_k %d, percentage %s, controls %d (planted among controls: %d), max_acceptable_k %s\n",
        sweep$best_k, sweep$best_percentage, sum(controls), sum(ctl_ids %in% planted), sweep$max_acceptable_k))
    withW <- function(object, k, ctl = controls, mode = "covariates") {
        f <- computeRuvCovariates(object, k, ctl, "group")
        samples <- as.character(object@design_matrix[[object@sample_id]])
        if (identical(mode, "correction")) {
            r <- object@args$limpa_dpc_quant_results
            E <- r$quantified_elist$E
            r$quantified_elist$E <- f$adjusted[rownames(E), colnames(E), drop = FALSE]
            object@args$limpa_dpc_quant_results <- r
            object@args$ruvCovariates <- list(ruv_mode = "correction", k = k)
        } else {
            for (i in seq_len(k)) object@design_matrix[[paste0("W", i)]] <- unname(f$W[samples, i])
            object@args$ruvCovariates <- list(ruv_mode = "covariates", k = k)
        }
        attr(object, "W") <- f$W
        object
    }
    k1 <- withW(pre, 1L); k2 <- withW(pre, 2L)
    W2 <- attr(k2, "W")
    gui_W <- as.matrix(ruv_state@design_matrix[, grep("^W[0-9]+$", names(ruv_state@design_matrix)), drop = FALSE])
    cat(sprintf("recomputed k=2 W equals the GUI state's W: max abs diff %.3g\n", max(abs(unname(W2) - unname(gui_W)))))
    g01 <- as.numeric(group[rownames(W2)] == "KO")
    cat("correlation of W with the group indicator (KO = 1):", paste(sprintf("W%d %.3f", seq_len(ncol(W2)), stats::cor(W2, g01)), collapse = ", "), "\n")
    cat("W (k=2):\n"); print(round(cbind(W2, KO = g01), 3))
    summariseCalls("(b) RUV k=1 covariates: dpcDE ~0+group+W1", daCalls(k1))
    summariseCalls("(c) RUV k=2 covariates: dpcDE ~0+group+W1+W2 (recomputed)", daCalls(k2))
    summariseCalls("(c') RUV k=2 correction mode (E <- RUV-III newY)", daCalls(withW(pre, 2L, mode = "correction")))

    # (d) limpa DPC matrix with plain limma, and a plain peptide summarisation.
    plainLimma <- function(object) { object@args$limpa_dpc_quant_results <- NULL; object }
    summariseCalls("(d1) limpa DPC matrix, plain limma (no dpcDE weights)", daCalls(plainLimma(pre)))
    y <- pre@args$limpa_dpc_quant_results$y_peptide_for_dpc
    fmap <- pre@args$peptide_feature_key_map
    prot <- fmap$active_protein_id[match(rownames(y), fmap$feature_key)]
    cat(sprintf("peptide matrix: %d x %d, missing %.1f%%, protein map NA %d\n", nrow(y), ncol(y), 100 * mean(is.na(y)), sum(is.na(prot))))
    pm <- rowsum(ifelse(is.na(y), 0, y), prot) / rowsum(1 * !is.na(y), prot)
    pm[!is.finite(pm)] <- NA
    meanState <- function(M) {
        object <- plainLimma(pre)
        tab <- data.frame(id = rownames(M), M, check.names = FALSE, stringsAsFactors = FALSE)
        names(tab)[1L] <- object@protein_id_column
        keep <- tab[[object@protein_id_column]] %in% object@protein_quant_table[[object@protein_id_column]]
        object@protein_quant_table <- tab[keep, c(object@protein_id_column, colnames(object@protein_quant_table)[-1L]), drop = FALSE]
        rownames(object@protein_quant_table) <- NULL
        object
    }
    stopifnot(identical(colnames(pre@protein_quant_table)[1L], pre@protein_id_column))
    summariseCalls("(d2) plain mean of log2 peptides per protein, limma", daCalls(meanState(pm)))

    # (f) the normalisation check: per-sample median centring of the quantified matrix (E and the protein table).
    centre <- function(object) {
        r <- object@args$limpa_dpc_quant_results
        m <- apply(r$quantified_elist$E, 2, stats::median)
        r$quantified_elist$E <- sweep(r$quantified_elist$E, 2, m)
        object@args$limpa_dpc_quant_results <- r
        samples <- names(m)
        object@protein_quant_table[samples] <- sweep(as.matrix(object@protein_quant_table[samples]), 2, m)
        object
    }
    pre_c <- centre(pre)
    summariseCalls("(f1) median-centred E, no RUV: dpcDE ~0+group", daCalls(pre_c))
    sweep_c <- suppressMessages(findBestNegCtrlPercentageCapped(pre_c, percentage_range = seq(parameters$percentage_min, parameters$percentage_max, by = 1),
        num_components_to_impute = parameters$num_components_to_impute, ruv_grouping_variable = parameters$ruv_grouping_variable,
        ruv_qval_cutoff = parameters$ruv_qval_cutoff, ruv_fdr_method = parameters$ruv_fdr_method, separation_metric = parameters$separation_metric,
        k_penalty_weight = parameters$k_penalty_weight, adaptive_k_penalty = isTRUE(parameters$adaptive_k_penalty), verbose = FALSE))
    cat(sprintf("RUV sweep on the centred matrix: best_k %d, percentage %s, controls %d\n", sweep_c$best_k, sweep_c$best_percentage, sum(sweep_c$best_control_genes_index)))
    summariseCalls(sprintf("(f2) median-centred E, RUV k=%d covariates (its own sweep)", sweep_c$best_k),
        daCalls(withW(pre_c, sweep_c$best_k, sweep_c$best_control_genes_index)))
    summariseCalls("(f3) median-centred E, RUV k=2 covariates (centred sweep's controls)",
        daCalls(withW(pre_c, 2L, sweep_c$best_control_genes_index)))
    pm_c <- sweep(pm, 2, apply(pm, 2, stats::median, na.rm = TRUE))
    summariseCalls("(f4) plain mean of log2 peptides, median-centred, limma", daCalls(meanState(pm_c)))

    # (e) the DIA route on the same fixture: the head state the DA loaded, through the same kernel.
    dia_head <- dia$states[[dia$head]]
    summariseCalls(sprintf("(e) DIA route (IQ MaxLFQ, %s) on the same fixture", dia$head), daCalls(dia_head))
    dia_files <- list.files(file.path(OUT, "dia_route_gui_files"), pattern = "\\.(tsv|tab|csv)$", full.names = TRUE)
    for (f in dia_files) {
        t <- tryCatch(utils::read.delim(f, stringsAsFactors = FALSE), error = function(e) NULL)
        if (is.null(t) || !"fdr_qvalue" %in% names(t)) next
        idcol <- intersect(c("Protein.Ids", "uniprot_acc", "Protein.Group", "protein_id"), names(t))[1L]
        if (is.na(idcol)) idcol <- names(t)[[1L]]
        sig <- pid(t[[idcol]])[t$fdr_qvalue < 0.05]
        cat(sprintf("  DIA GUI file %s: rows %d, q<0.05 %d, planted %d, FP %d\n", basename(f), nrow(t), length(sig), sum(sig %in% planted), sum(!sig %in% planted)))
    }
    dmed <- apply(as.matrix(dia_head@protein_quant_table[, -1L]), 2, stats::median, na.rm = TRUE)
    cat("DIA route: per-sample median of its normalised protein matrix:\n"); print(round(dmed, 3))
    dia_norm <- dia$config_list$normaliseBetweenSamples
    if (!is.null(dia_norm)) { cat("DIA route normaliseBetweenSamples config:\n"); str(dia_norm) }
    cat("DIA_LIMPA limpaProteinNormalisation config:\n"); str(limpa$config_list$limpaProteinNormalisation)

    table <- do.call(rbind, rows)
    utils::write.table(table, file.path(OUT, "f19_variants.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    cat("\nwrote", file.path(OUT, "f19_variants.tsv"), "\n")
}
