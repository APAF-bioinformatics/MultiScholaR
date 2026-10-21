# KIT-053 script-equivalence check (executable): the limpa script's chunks 11-22 run with the package's own functions in
# process on a bound DIA-NN corpus with `engine = limpa` (the script's engine), starting from the rollup state of a
# route store recorded by kit-053/oracle-record.R, and the contrast tables compared with that store's owner-path DA
# (the route's output). The replicate structure the route resolves (KIT-035's rule) is recorded beside the script's
# fixed `group`; on a design without technical replicates the two coincide.
# Usage: [KIT053_INSTALLED=1] R_LIBS=<library> Rscript docs/qa/workflow-execution/kit/limpa-workflow/script-equivalence.R <corpus key> [label]
if (identical(Sys.getenv("KIT053_INSTALLED"), "1")) suppressMessages(library(MultiScholaR)) else suppressMessages(devtools::load_all(quiet = TRUE))
started <- Sys.time()
args <- commandArgs(TRUE)
corpus_key <- args[[1L]]
label <- if (length(args) > 1L) args[[2L]] else "tree"
repo <- normalizePath(".")
`%||%` <- function(a, b) if (is.null(a)) b else a
ns <- asNamespace("MultiScholaR"); api <- function(name) get(name, ns)
for (f in ls(ns, all.names = FALSE)) if (!exists(f, envir = globalenv(), inherits = TRUE) && is.function(get(f, ns))) assign(f, get(f, ns), envir = globalenv())
harness <- new.env(parent = globalenv()); sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), harness)
recorded <- readRDS(file.path(repo, "docs/qa/workflow-execution/kit/kit-053", paste0("oracles-limpa_", corpus_key, ".rds")))
stopifnot(identical(recorded$route, "DIA_LIMPA"), dir.exists(recorded$store))
store_root <- dirname(recorded$store)
identity <- readSavedStoreIdentity(file.path(store_root, ".workflow-entities", "store.json"))
manager <- newDIAStateAdapter(newWorkflowEntityStore(store_root, identity$project_id, identity$workflow_id), reopen = TRUE)
sha <- function(x) digest::digest(x, algo = "sha256")
table_sha <- function(df) { df <- as.data.frame(df); rownames(df) <- NULL; sha(list(names(df), lapply(df, function(col) { attributes(col) <- NULL; col }))) }
out <- list(schema = "kit.script_equivalence.v1", ticket = "KIT-053", corpus = corpus_key, label = label, date = format(Sys.Date()),
    head_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE), store = recorded$store, store_build = recorded$build,
    script = list(repository = "APAF-bioinformatics/e_35036_MorpheusTex_GarryNunn_20260329_b", commit = "4aa97fc3b85e897ee457052682d8bba2a7afbbc4",
        file = "scripts/proteomics_isoseq/DIA_workflow_limpa_isoseq_tissue.Rmd", chunks = "11-22 with package functions; findBestNegCtrlPercentageCustom vendored"),
    engine = "limpa", route_dpc_engine = recorded$settings$dpc_engine)
vendored <- new.env(parent = globalenv())
sys.source(file.path(repo, "docs/qa/workflow-execution/kit/limpa-workflow/findBestNegCtrlPercentageCustom.R"), vendored)
withDIAStateScope(manager, function() {
    # The script's chunks from the rollup state (chunk 12) of the route's store, in memory.
    rolled <- manager$getState("precursor_rollup")
    fasta_path <- list.files(dirname(dirname(store_root)), pattern = "^aa_seq_tbl_final[.]RDS$", recursive = TRUE, full.names = TRUE)
    stopifnot(length(fasta_path) >= 1L)
    fasta <- readRDS(fasta_path[[1L]])
    t0 <- Sys.time()
    log2_obj <- log2TransformPeptideMatrix(calcPeptideMatrix(rolled))                                     # chunk 13
    seqinr <- dplyr::rename(fasta, uniprot_acc = "database_id")
    acc_obj <- suppressMessages(chooseBestProteinAccession(theObject = log2_obj, delim = ";", seqinr_obj = seqinr,   # chunk 16
        seqinr_accession_column = "uniprot_acc", aggregation_method = "mean"))
    dpc_obj <- suppressMessages(proteinMissingValueImputationLimpa(acc_obj, dpc_slope = 0.8, verbose = FALSE, chunk = 1000, engine = "limpa")) # chunk 17
    normalised <- dpc_obj                                                                                  # chunk 17b (skipped in the script)
    filtered <- suppressMessages(filterMinNumPeptidesPerProtein(normalised, num_peptides_per_protein_thresh = 1, num_peptidoforms_per_protein_thresh = 2, verbose = FALSE)) # chunk 18
    pairs <- suppressMessages(pearsonCorForSamplePairs(filtered, tech_rep_remove_regex = "pool", correlation_group = "group"))        # chunk 19a
    correlated <- suppressMessages(filterSamplesByProteinCorrelationThreshold(filtered, pearson_correlation_per_pair = pairs, min_pearson_correlation_threshold = 0.5))
    correlated <- subsetLimpaResultsToSamples(correlated)
    sweep <- suppressMessages(suppressWarnings(vendored$findBestNegCtrlPercentageCustom(correlated, percentage_range = seq(1, 30, by = 1), # chunk 19b
        num_components_to_impute = 5, ruv_grouping_variable = "group", ruv_fdr_method = "BH", ruv_qval_cutoff = 0.05,
        separation_metric = "max_difference", k_penalty_weight = 0.5, adaptive_k_penalty = TRUE, verbose = FALSE)))
    best_k <- as.integer(sweep$best_k); controls <- sweep$best_control_genes_index
    # The script's replicate matrix is the design's `group`; the route's rule (KIT-035) is the technical block when the
    # design declares one, else the grouping variable. Both are computed; the comparison follows the route's rule.
    design <- correlated@design_matrix
    samples <- as.character(design[[correlated@sample_id]])
    Y <- t(as.matrix(tibble::column_to_rownames(as.data.frame(correlated@protein_quant_table), correlated@protein_id_column))[, samples])
    M_script <- getRuvIIIReplicateMatrixHelper(design, !!rlang::sym(correlated@sample_id), !!rlang::sym("group"))[samples, , drop = FALSE]
    blocks <- resolveProtRuvReplicateBlocks(design, correlated@sample_id, correlated@group_id, correlated@technical_replicate_id, "group")
    M_route <- getRuvIIIReplicateMatrixHelper(blocks$design, !!rlang::sym(correlated@sample_id), .ruv_replicate_block)[samples, , drop = FALSE]
    out$replicate_structure <<- list(route = blocks$structure, script = "group", coincide = identical(unname(M_script[, order(colnames(M_script))]), unname(M_route[, order(colnames(M_route))])))
    ruv_res <- ruv::RUVIII(Y = Y, M = M_route, ctl = controls, k = best_k, return.info = TRUE)
    alpha <- ruv_res$fullalpha[seq_len(best_k), , drop = FALSE]; ac <- alpha[, controls, drop = FALSE]
    W <- Y[, controls, drop = FALSE] %*% t(ac) %*% solve(ac %*% t(ac)); colnames(W) <- paste0("W", seq_len(best_k))
    with_w <- correlated
    for (i in seq_len(best_k)) with_w@design_matrix[[paste0("W", i)]] <- unname(W[samples, i])
    # chunk 22: dpcDE with the formula plus the W terms on the quantified EList (the DA helper's branch does the same).
    da <- suppressMessages(differentialAbundanceAnalysis(with_w, contrasts_tbl = data.frame(contrasts = recorded_contrasts <- {
        m <- resolveDIAStateValue(manager$.context, manager$.context$snapshot$da_results); m$contrasts$full_format }, stringsAsFactors = FALSE),
        formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0, args_row_id = with_w@protein_id_column,
        qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue"))
    out$seconds_script_chain <<- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    out$script_chain <<- list(dpc = unname(as.numeric(dpc_obj@args$limpa_dpc_quant_results$dpc_parameters_used)), proteins = nrow(filtered@protein_quant_table),
        samples = length(samples), best_percentage = sweep$best_percentage, best_k = best_k, controls = sum(controls), da_model = da$da_model$model,
        formula_used = da$da_model$formula_used)
    # The route's owner-path DA tables from the store's manifest.
    manifest <- resolveDIAStateValue(manager$.context, manager$.context$snapshot$da_results)
    owner_long <- materializeDIATable(manager$.context, manifest$results[[1L]]$tables$da_proteins_long)
    owner_wide <- materializeDIATable(manager$.context, manifest$results[[1L]]$tables$da_proteins_wide)
    id <- with_w@protein_id_column
    key <- c(id, "comparison", "log2FC", "raw_pvalue", "fdr_qvalue")
    memory_long <- as.data.frame(da$da_proteins_long)
    out$compare <<- list(owner_rows = nrow(owner_long), memory_rows = nrow(memory_long), id_column = id,
        key_columns = key[key %in% intersect(names(owner_long), names(memory_long))],
        long_key_sha_owner = table_sha(owner_long[order(owner_long[[id]]), intersect(key, names(owner_long))]),
        long_key_sha_memory = table_sha(memory_long[order(memory_long[[id]]), intersect(key, names(memory_long))]),
        wide_sha_owner = table_sha(owner_wide), wide_sha_memory = table_sha(as.data.frame(da$da_proteins_wide)),
        route_model = manifest$model[c("model", "formula_used", "k", "block_reason", "route")])
    out$compare$long_identical <<- identical(out$compare$long_key_sha_owner, out$compare$long_key_sha_memory)
    out$compare$wide_identical <<- identical(out$compare$wide_sha_owner, out$compare$wide_sha_memory)
    if (!out$compare$long_identical) {
        a <- owner_long[order(owner_long[[id]]), ]; b <- memory_long[order(memory_long[[id]]), ]
        common <- intersect(a[[id]], b[[id]])
        out$compare$max_abs_dlogFC <<- max(abs(a$log2FC[match(common, a[[id]])] - b$log2FC[match(common, b[[id]])]))
        out$compare$max_abs_dp <<- max(abs(a$raw_pvalue[match(common, a[[id]])] - b$raw_pvalue[match(common, b[[id]])]))
    }
})
closeDIAStateAdapter(manager)
out$seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
cat("SCRIPT CHAIN:", paste(names(out$script_chain), unlist(lapply(out$script_chain, function(x) paste(if (is.numeric(x)) signif(x, 6) else x, collapse = "/"))), collapse = " "), "\n")
cat("REPLICATE STRUCTURE: route", out$replicate_structure$route, "| script group | coincide", out$replicate_structure$coincide, "\n")
cat("DA TABLES: long identical", out$compare$long_identical, "| wide identical", out$compare$wide_identical, "| rows", out$compare$owner_rows, out$compare$memory_rows,
    "| route model", paste(unlist(out$compare$route_model), collapse = "/"), if (!is.null(out$compare$max_abs_dlogFC)) paste("| max|dlogFC|", signif(out$compare$max_abs_dlogFC, 3), "max|dp|", signif(out$compare$max_abs_dp, 3)) else "", "\n")
jsonlite::write_json(out, file.path(repo, "docs/qa/workflow-execution/kit/kit-053", paste0("script-equivalence-", corpus_key, "-", label, ".json")), auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = NA)
cat("TOTAL", round(out$seconds), "s\n")
