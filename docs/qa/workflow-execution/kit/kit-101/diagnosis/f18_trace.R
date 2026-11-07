# KIT-101 revision 2.3, F-18 diagnosis: what the MaxQuant and FragPipe LFQ imports do with a zero intensity, traced
# through the proteomics LFQ route with the production code. Diagnosis only: no package code is changed.
#   Part A: the two readers on a tiny constructed table (3 proteins x 4 samples, one 0 and one empty cell each).
#   Part B: the entry-suite lane driver (tests/testthat/helper-gui-entry.R; mod_proteomics_server under testServer)
#           on the prot_lfq / prot_lfq_fragpipe lanes, with fixtures written by the committed generator at the KIT-021
#           recipe's own 8% missing (the lanes' committed fixtures use 0%), in two arms that differ only in how a
#           missing cell is written: "zero" (0, as MaxQuant/FragPipe write it) and "empty" (an empty field, read as NA).
#           Each run stops after the DA and dumps every committed state; one tracer cell (a non-planted protein with
#           exactly one missing cell) is followed through the import table, every state and the DA.
#   Part C: aggregate zero counts of the synthetic LFQ corpora the KIT-021 LFQ oracles were recorded on (counts only).
# Usage (repository root): R_LIBS=<kit_100 library> Rscript docs/qa/workflow-execution/kit/kit-101/diagnosis/f18_trace.R
# Optional: F18_OUT=<dir> (default tempdir); F18_ARMS=<key,key> to run a subset of the four lane arms.

suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = FALSE))
OUT <- Sys.getenv("F18_OUT", file.path(tempdir(), "f18"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
cat("output:", OUT, "\n")
planted <- sprintf("P%05d", 1:40)
pid <- function(x) { out <- rep(NA_character_, length(x)); hit <- grepl("P[0-9]{5}", x); out[hit] <- regmatches(x, regexpr("P[0-9]{5}", x)); out }

# ------------------------------------------------------------------------------------------------------------------
# Part A: the readers on a tiny constructed input.
# ------------------------------------------------------------------------------------------------------------------
cat("\n== Part A: readers ==\n")
tiny <- file.path(OUT, "tiny"); dir.create(tiny, showWarnings = FALSE)
mq <- file.path(tiny, "proteinGroups.txt")
writeLines(c(
    paste("Protein IDs", "Majority protein IDs", "Gene names", "LFQ intensity WT_1", "LFQ intensity WT_2", "LFQ intensity KO_1",
        "LFQ intensity KO_2", "Potential contaminant", "Reverse", sep = "\t"),
    paste("P00001", "P00001", "G1", "1000", "0", "1200", "1100", "", "", sep = "\t"),
    paste("P00002", "P00002", "G2", "2000", "", "2100", "1900", "", "", sep = "\t"),
    paste("P00003", "P00003", "G3", "3000", "3100", "2900", "3050", "", "", sep = "\t")), mq)
fp <- file.path(tiny, "combined_protein.tsv")
writeLines(c(
    paste("Protein", "Protein ID", "Gene", "WT_1 MaxLFQ Intensity", "WT_2 MaxLFQ Intensity", "KO_1 MaxLFQ Intensity",
        "KO_2 MaxLFQ Intensity", sep = "\t"),
    paste("sp|P00001|A", "P00001", "G1", "1000", "0", "1200", "1100", sep = "\t"),
    paste("sp|P00002|B", "P00002", "G2", "2000", "", "2100", "1900", sep = "\t"),
    paste("sp|P00003|C", "P00003", "G3", "3000", "3100", "2900", "3050", sep = "\t")), fp)
show <- function(label, res) {
    d <- as.data.frame(res$data)
    cell <- function(p, r) d$Intensity[d$Protein.Ids == p & d$Run == r]
    cat(sprintf("%-9s P00001/WT_2 (file '0') -> %s | P00002/WT_2 (file empty) -> %s | quantity_col %s\n", label,
        format(cell("P00001", "WT_2")), format(cell("P00002", "WT_2")), res$column_mapping$quantity_col))
    cat(sprintf("          the design step's transform (mod_prot_design_state_helpers.R:74, log2(.x + 1)): 0 -> %s, NA -> %s\n",
        format(log2(cell("P00001", "WT_2") + 1)), format(log2(cell("P00002", "WT_2") + 1))))
}
show("MaxQuant", suppressMessages(importMaxQuantData(mq)))
show("FragPipe", suppressMessages(importFragPipeData(fp)))

# ------------------------------------------------------------------------------------------------------------------
# Part B: the LFQ lanes, zero against empty.
# ------------------------------------------------------------------------------------------------------------------
cat("\n== Part B: LFQ route, zero against empty ==\n")
gen_env <- new.env(parent = globalenv())
for (e in parse("tests/testdata/e2e/generate_fixtures.R")) {
    if (is.call(e) && identical(e[[1L]], as.name("<-")) && is.call(e[[3L]]) && identical(e[[3L]][[1L]], as.name("function"))) eval(e, gen_env)
}
assign("FIXTURE_ROOT", normalizePath("tests/testdata/e2e"), gen_env)
suppressMessages(testthat::source_test_helpers("tests/testthat", env = globalenv()))
original_steps <- guiEntrySteps
original_lane <- guiEntryLane
original_fixture_dir <- guiEntryFixtureDir
arms <- list()
for (lane_id in c("prot_lfq", "prot_lfq_fragpipe")) {
    base <- original_lane(lane_id)
    base$fixture_recipe$missing_fraction <- 0.08
    truth <- gen_env$protFixtureMatrix(base, gen_env$kitProteinTools())
    for (arm in c("zero", "empty")) {
        dir <- file.path(OUT, paste0(lane_id, "_", arm)); dir.create(dir, showWarnings = FALSE)
        generator <- if (identical(base$import_tool, "maxquant")) gen_env$generateProtMaxQuantFixture else gen_env$generateProtFragPipeFixture
        invisible(capture.output(generator(base, dir)))
        path <- file.path(dir, base$seed_file)
        if (identical(arm, "empty")) {
            t <- utils::read.delim(path, colClasses = "character", check.names = FALSE, quote = "", na.strings = character())
            cols <- grep("Intensity$|^LFQ\\.intensity\\.", names(t), value = TRUE)
            for (col in cols) t[[col]][t[[col]] == "0"] <- ""
            utils::write.table(t, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
        }
        t <- utils::read.delim(path, colClasses = "character", check.names = FALSE, quote = "", na.strings = character())
        cols <- grep("Intensity$|^LFQ\\.intensity\\.", names(t), value = TRUE)
        cat(sprintf("%s/%s fixture: %d rows, %d intensity cells, %d written as 0, %d empty\n", lane_id, arm, nrow(t),
            nrow(t) * length(cols), sum(as.matrix(t[cols]) == "0"), sum(as.matrix(t[cols]) == "")))
        arms[[paste0(lane_id, "_", arm)]] <- list(lane_id = lane_id, arm = arm, dir = dir, truth = truth)
    }
}
# The tracer: per lane, the first non-planted protein with exactly one missing cell.
tracerOf <- function(truth) {
    ids <- sprintf("P%05d", seq_len(nrow(truth$missing)))
    one <- which(rowSums(truth$missing) == 1L & !truth$differential)[[1L]]
    list(protein = ids[[one]], sample = truth$design$sample[which(truth$missing[one, ])])
}
dumpStep <- function(tag) guiEntryStep("diag_dump", "da", function(api) {
    wd <- api$workflow_data
    manager <- shiny::isolate(wd$state_manager)
    history <- unique(guiEntryHistory(manager))
    states <- list()
    for (name in history) states[[name]] <- tryCatch(withWorkflowStateScope(manager, function() manager$getState(name)),
        error = function(e) paste("<error>", conditionMessage(e)))
    data_cln <- tryCatch(shiny::isolate(getDIAOwnedPayload(wd, "data_cln")), error = function(e) paste("<error>", conditionMessage(e)))
    saveRDS(list(history = history, head = manager$current_state, states = states, data_cln = data_cln,
        contrasts_tbl = shiny::isolate(wd$contrasts_tbl)), file.path(OUT, paste0(tag, "_states.rds")))
})
guiEntrySteps <- function(lane) {
    steps <- original_steps(lane)
    ids <- vapply(steps, `[[`, "", "id")
    c(steps[seq_len(which(ids == "da"))], list(dumpStep(lane$diag_tag)))
}
guiEntryLane <- function(lane_id) {
    a <- arms[[lane_id]]
    lane <- original_lane(a$lane_id)
    lane$diag_tag <- lane_id
    lane$diag_dir <- a$dir
    lane
}
guiEntryFixtureDir <- function(lane) if (!is.null(lane$diag_dir)) normalizePath(lane$diag_dir) else original_fixture_dir(lane)

daCalls <- function(object, contrast) {
    single <- data.frame(contrasts = contrast, stringsAsFactors = FALSE)
    arguments <- list(contrasts_tbl = single, formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0,
        eBayes_trend = TRUE, eBayes_robust = TRUE, args_group_pattern = "(\\d+)", args_row_id = object@protein_id_column,
        qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue")
    for (section in c("differentialAbundanceAnalysis", "differentialAbundanceAnalysisHelper"))
        if (is.null(object@args[[section]])) object@args[[section]] <- arguments
    result <- NULL
    invisible(utils::capture.output(result <- suppressMessages(suppressWarnings(differentialAbundanceAnalysis(object,
        contrasts_tbl = single, formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0,
        qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue")))))
    tab <- result$contrasts_results_table
    data.frame(protein = pid(tab[[object@protein_id_column]]), logFC = tab$logFC, q = tab$fdr_qvalue,
        bh = tab$fdr_value_bh_adjustment, stringsAsFactors = FALSE)
}

summary_rows <- list()
# F18_ARMS (comma-separated keys) reruns a subset of the arms (the correction pass after a DA failure on one arm).
run_arms <- if (nzchar(Sys.getenv("F18_ARMS"))) strsplit(Sys.getenv("F18_ARMS"), ",")[[1L]] else names(arms)
for (key in run_arms) {
    a <- arms[[key]]
    tr <- tracerOf(a$truth)
    record <- guiEntryRunLane(key)
    notes <- unlist(lapply(record$steps, function(s) vapply(s$notifications, function(n) paste0("[", s$step, "/", n$type, "] ", substr(n$message, 1, 140)), "")))
    cat(sprintf("\n-- %s: lane ran %d steps in %.0f s; head %s\n", key, length(record$steps), record$seconds, record$steps[[length(record$steps)]]$head))
    bad <- grep("error|warn", notes, value = TRUE, ignore.case = TRUE)
    if (length(bad)) cat(paste(" ", bad, collapse = "\n"), "\n")
    d <- readRDS(file.path(OUT, paste0(key, "_states.rds")))
    cat(sprintf("tracer: %s / %s (true value %.1f, log2 %.2f)\n", tr$protein, tr$sample,
        a$truth$values[as.integer(sub("P", "", tr$protein)), match(tr$sample, a$truth$design$sample)],
        log2(a$truth$values[as.integer(sub("P", "", tr$protein)), match(tr$sample, a$truth$design$sample)])))
    if (is.data.frame(d$data_cln)) {
        dc <- as.data.frame(d$data_cln)
        v <- dc$Intensity[pid(dc$Protein.Ids) %in% tr$protein & tolower(dc$Run) == tolower(tr$sample)]
        cat(sprintf("  %-28s Intensity = %s   (import table rows %d, Intensity == 0: %d, NA: %d)\n", "import (data_cln)",
            paste(format(v), collapse = ","), nrow(dc), sum(dc$Intensity == 0, na.rm = TRUE), sum(is.na(dc$Intensity))))
    }
    first_samples <- NULL
    for (name in d$history) {
        s <- d$states[[name]]
        if (!methods::is(s, "ProteinQuantitativeData")) next
        tab <- s@protein_quant_table
        samples <- setdiff(names(tab), s@protein_id_column)
        if (is.null(first_samples)) first_samples <- samples
        row <- which(pid(tab[[s@protein_id_column]]) == tr$protein)
        col <- samples[tolower(samples) == tolower(tr$sample)]
        m <- as.matrix(tab[samples])
        value <- if (length(row) && length(col)) format(round(tab[[col]][row[[1L]]], 4)) else "<row or sample removed>"
        cat(sprintf("  %-28s value = %-22s proteins %4d samples %2d | cells == 0: %4d, NA: %4d, min %.2f\n", name, value,
            nrow(tab), length(samples), sum(m == 0, na.rm = TRUE), sum(is.na(m)), min(m, na.rm = TRUE)))
    }
    head <- d$states[[d$head]]
    final_samples <- setdiff(names(head@protein_quant_table), head@protein_id_column)
    removed <- setdiff(first_samples, final_samples)
    calls <- tryCatch(daCalls(head, as.character(d$contrasts_tbl[[1L]][[1L]])), error = function(e) {
        cat("  DA failed:", conditionMessage(e), "\n")
        data.frame(protein = character(), logFC = numeric(), q = numeric(), bh = numeric())
    })
    sig <- calls$protein[!is.na(calls$q) & calls$q < 0.05]
    trow <- calls[calls$protein == tr$protein, ]
    cat(sprintf("  samples removed before DA: %s\n", if (length(removed)) paste(removed, collapse = ", ") else "none"))
    cat(sprintf("  DA (head %s): proteins %d, q<0.05 %d, planted %d/40, false positives %d; tracer logFC %s q %s\n", d$head,
        nrow(calls), length(sig), sum(sig %in% planted), sum(!sig %in% planted),
        if (nrow(trow)) format(round(trow$logFC, 3)) else "NA", if (nrow(trow)) format(signif(trow$q, 3)) else "NA"))
    summary_rows[[key]] <- data.frame(arm = key, samples_in = length(first_samples), samples_removed = length(removed),
        removed = paste(removed, collapse = ";"), da_proteins = nrow(calls), sig = length(sig), planted_recall = sum(sig %in% planted),
        false_pos = sum(!sig %in% planted), tracer = paste(tr$protein, tr$sample), tracer_logFC = if (nrow(trow)) trow$logFC else NA,
        tracer_q = if (nrow(trow)) trow$q else NA, stringsAsFactors = FALSE)
}
summary <- do.call(rbind, summary_rows)
utils::write.table(summary, file.path(OUT, "f18_arms.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n"); print(summary, row.names = FALSE)

# ------------------------------------------------------------------------------------------------------------------
# Part C: zero counts of the synthetic LFQ corpora (aggregate counts only; no real dataset is opened).
# ------------------------------------------------------------------------------------------------------------------
cat("\n== Part C: synthetic LFQ corpora, zero counts ==\n")
corpora <- c(
    lfq_300 = "/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/kit/corpus/lfq-300/proteinGroups.txt",
    lfq_fixture = "/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/kit/corpus/lfq-fixture/proteinGroups.txt",
    e2e_prot_lfq = "tests/testdata/e2e/prot_lfq/proteinGroups.txt",
    e2e_prot_lfq_seed_fragpipe = "tests/testdata/e2e/prot_lfq/seed_combined_protein.tsv",
    e2e_prot_lfq_maxquant = "tests/testdata/e2e/prot_lfq_maxquant/proteinGroups.txt",
    e2e_prot_lfq_fragpipe = "tests/testdata/e2e/prot_lfq_fragpipe/combined_protein.tsv")
for (name in names(corpora)) {
    path <- corpora[[name]]
    if (!file.exists(path)) { cat(name, ": absent\n"); next }
    t <- utils::read.delim(path, colClasses = "character", check.names = FALSE, quote = "", na.strings = character())
    cols <- grep("Intensity$|^LFQ\\.intensity\\.|^LFQ intensity ", names(t), value = TRUE)
    m <- as.matrix(t[cols])
    cat(sprintf("%-28s sha256 %s… rows %5d, intensity columns %3d, cells %8d, '0' %7d (%.2f%%), empty %d\n", name,
        substr(digest::digest(file = path, algo = "sha256"), 1, 16), nrow(t), length(cols), length(m), sum(m == "0"),
        100 * mean(m == "0"), sum(m == "")))
}
