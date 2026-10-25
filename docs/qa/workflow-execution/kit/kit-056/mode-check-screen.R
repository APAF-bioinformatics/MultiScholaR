# KIT-056 mode check (in process, the screen corpus): two recorder stores of the DIA_LIMPA route on the same build,
# one in the covariates mode (kit-057's record) and one in the correction mode (kit-056's record), compared:
#   (a) the same sweep: percentage, k, control set; the covariates design has W1..Wk, the correction design none;
#   (b) the correction state's matrix == computeRuvCovariates()$adjusted == t(ruv::RUVIII(Y, M, ctl, k)$newY) (exact),
#       max |delta| against RUVIIIC::RUVIII_C_Varying on the same inputs (characterisation); the owned EList's E equals
#       the matrix; standard errors and observation counts identical to the source's; the owned entity differs;
#   (c) the correction-mode DA (the route's runner, in the store's manifest) == a direct limpa::dpcDE on the corrected
#       EList with the plain design -> convertDpcDAToStandardFormat (long table identical);
#   (d) covariates-mode vs correction-mode DA tables: max |dlogFC|, max |dp|, overlap of the q < 0.05 sets (descriptive).
# Usage: [KIT053_INSTALLED=1] R_LIBS=<library> Rscript docs/qa/workflow-execution/kit/kit-056/mode-check-screen.R <covariates project> <correction project> <label>
if (identical(Sys.getenv("KIT053_INSTALLED"), "1")) suppressMessages(library(MultiScholaR)) else suppressMessages(devtools::load_all(quiet = TRUE))
`%||%` <- function(a, b) if (is.null(a)) b else a
ns <- asNamespace("MultiScholaR"); api <- function(name) get(name, ns)
h <- new.env(parent = globalenv()); sys.source("tools/profiling/workflow_kit.R", h)
args <- commandArgs(TRUE); cov_project <- args[[1L]]; corr_project <- args[[2L]]; label <- if (length(args) > 2L) args[[3L]] else "tree"
open <- function(project) { store <- h$kitFindStore(project); root <- dirname(store); id <- api("readSavedStoreIdentity")(file.path(root, ".workflow-entities", "store.json"))
    st <- api("newWorkflowEntityStore")(root, id$project_id, id$workflow_id); list(m = api("newDIAStateAdapter")(st, reopen = TRUE), store = st) }
A <- open(cov_project); B <- open(corr_project)
out <- list(schema = "kit.mode_check.v1", ticket = "KIT-056", label = label, date = format(Sys.Date()), head_commit = system("git rev-parse HEAD", intern = TRUE))
read_state <- function(x, state) api("withDIAStateScope")(x$m, function() api("readDIAExecutionComponents")(x$m, api("getDIAStateReference")(x$m, state)))
cov <- read_state(A, "ruv_covariates"); corr <- read_state(B, "ruv_corrected"); src <- read_state(B, "correlation_filtered")
# (a)
a <- cov@args$ruvCovariates; b <- corr@args$ruvCovariates
out$a_same_sweep <- list(k = list(covariates = a$k, correction = b$k, identical = identical(a$k, b$k)),
    percentage = list(covariates = a$percentage_as_neg_ctrl, correction = b$percentage_as_neg_ctrl, identical = identical(a$percentage_as_neg_ctrl, b$percentage_as_neg_ctrl)),
    control_count = list(covariates = a$control_count, correction = b$control_count, identical = identical(a$control_count, b$control_count)),
    modes = list(covariates = a$ruv_mode, correction = b$ruv_mode),
    w_in_design = list(covariates = grep("^W[0-9]+$", names(cov@design_matrix), value = TRUE), correction = grep("^W[0-9]+$", names(corr@design_matrix), value = TRUE)))
# (b) recompute on the correction store's correlation state
sweep <- suppressMessages(suppressWarnings(findBestNegCtrlPercentageCapped(src, verbose = FALSE)))
factors <- computeRuvCovariates(src, sweep$best_k, sweep$best_control_genes_index, a$ruv_grouping_variable %||% "group")
key <- corr@protein_id_column; samples <- as.character(corr@design_matrix[[corr@sample_id]])
state_matrix <- as.matrix(corr@protein_quant_table[, samples]); rownames(state_matrix) <- corr@protein_quant_table[[key]]
adj <- factors$adjusted[rownames(state_matrix), samples]
Y <- t(as.matrix(src@protein_quant_table[, samples])); colnames(Y) <- src@protein_quant_table[[key]]
blocks <- api("resolveProtRuvReplicateBlocks")(src@design_matrix, src@sample_id, src@group_id, src@technical_replicate_id, a$ruv_grouping_variable %||% "group")
M <- api("getRuvIIIReplicateMatrixHelper")(blocks$design, !!rlang::sym(src@sample_id), .ruv_replicate_block); M <- M[samples, , drop = FALSE]
direct <- t(ruv::RUVIII(Y = Y, M = M, ctl = sweep$best_control_genes_index, k = sweep$best_k))[rownames(state_matrix), samples]
ruviiic <- tryCatch({ r <- RUVIIIC::RUVIII_C_Varying(k = sweep$best_k, Y = Y, M = M, toCorrect = colnames(Y), potentialControls = colnames(Y)[sweep$best_control_genes_index]); t(r)[rownames(state_matrix), samples] }, error = function(e) NULL)
E <- corr@args$limpa_dpc_quant_results$quantified_elist$E[rownames(state_matrix), samples]
out$b_correction_matrix <- list(sweep_k = sweep$best_k, state_equals_computeRuvCovariates = isTRUE(all.equal(unname(state_matrix), unname(adj), tolerance = 1e-12)),
    state_equals_ruv_RUVIII_newY = isTRUE(all.equal(unname(state_matrix), unname(direct), tolerance = 1e-12)),
    max_abs_delta_vs_RUVIIIC = if (is.null(ruviiic)) NA_real_ else max(abs(state_matrix - ruviiic)),
    owned_E_equals_state_matrix = isTRUE(all.equal(unname(E), unname(state_matrix), tolerance = 1e-12)),
    standard_errors_identical_to_source = identical(corr@args$limpa_dpc_quant_results$standard_errors, src@args$limpa_dpc_quant_results$standard_errors),
    n_observations_identical_to_source = identical(corr@args$limpa_dpc_quant_results$n_observations, src@args$limpa_dpc_quant_results$n_observations))
owned_id <- function(x, state) api("withDIAStateScope")(x$m, function() api("readDIAStateMetadata")(x$m, api("getDIAStateReference")(x$m, state))$data$slots$args$limpa_dpc_quant_results$entity_id)
out$b_correction_matrix$owned_entity_changed <- !identical(owned_id(B, "ruv_corrected"), owned_id(B, "correlation_filtered"))
# (c) the runner's DA (the store's manifest) against a direct dpcDE on the corrected EList
manifest_of <- function(x) api("withDIAStateScope")(x$m, function() api("resolveDIAStateValue")(x$m$.context, x$m$.context$snapshot$da_results))
table_of <- function(x, ref) api("withDIAStateScope")(x$m, function() as.data.frame(api("materializeDIATable")(x$m$.context, ref)))
mb <- manifest_of(B); ma <- manifest_of(A)
long_of <- function(x, manifest) do.call(rbind, lapply(manifest$results, function(r) table_of(x, r$tables$da_proteins_long)))
long_b <- long_of(B, mb); long_a <- long_of(A, ma)
design <- model.matrix(as.formula(mb$model$formula_used), corr@design_matrix); rownames(design) <- samples
replicate <- api("resolveProtTechnicalReplicateBlock")(corr@design_matrix, corr@group_id, corr@technical_replicate_id)
elist <- corr@args$limpa_dpc_quant_results$quantified_elist[, samples]
fit <- limpa::dpcDE(elist, design, plot = FALSE, block = if (isTRUE(replicate$technical)) replicate$block else NULL)
contrast <- as.list(mb$contrasts$full_format %||% mb$contrasts$contrasts)
out$c_da_direct <- list(formula = mb$model$formula_used, model = mb$model$model, ruv_mode = mb$model$ruv_mode, k = mb$model$k, contrasts = contrast,
    long_rows = nrow(long_b), long_columns = names(long_b), note = "the runner's long table is compared with a direct dpcDE -> convertDpcDAToStandardFormat on the same contrasts")
protein_col <- if (key %in% names(long_b)) key else intersect(c("Protein.Ids", "uniprot_acc", "protein_id", mb$protein_id_column), names(long_b))[1]
p_col <- intersect(c("raw_pvalue", "P.Value", "p_value", "pvalue"), names(long_b))[1]; q_col <- intersect(c("fdr_qvalue", "adj.P.Val", "q_value", "qvalue"), names(long_b))[1]
lfc_col <- intersect(c("log2FC", "logFC"), names(long_b))[1]
out$c_da_direct$columns_used <- list(protein = protein_col, p = p_col, q = q_col)
# The runner's eBayes settings (the DIA DA runner fixes trend and robust to TRUE; the manifest's parameters when recorded).
trend <- isTRUE(mb$parameters$eBayes_trend %||% TRUE); robust <- isTRUE(mb$parameters$eBayes_robust %||% TRUE)
out$c_da_direct$eBayes <- list(trend = trend, robust = robust)
direct_da <- tryCatch({ r <- convertDpcDAToStandardFormat(dpc_fit = fit, contrast_strings = unlist(contrast), design_matrix = design, eBayes_trend = trend, eBayes_robust = robust); r }, error = function(e) list(error = conditionMessage(e)))
if (is.null(direct_da$error)) { d <- do.call(rbind, lapply(names(direct_da$results), function(nm) { t <- direct_da$results[[nm]]; t$comparison <- nm; t }))
    out$c_da_direct$direct_columns <- names(d)
    d_protein <- as.character(d$uniprot_acc %||% d[[protein_col]] %||% rownames(d)); d_comp <- as.character(d$comparison)
    key_b <- paste(long_b[[protein_col]], long_b$comparison); key_d <- paste(d_protein, d_comp); m <- match(key_b, key_d)
    if (all(is.na(m)) && length(unique(long_b$comparison)) == 1L) m <- match(as.character(long_b[[protein_col]]), d_protein)
    out$c_da_direct$sample_keys <- list(runner = head(key_b, 2), direct = head(key_d, 2))
    d_p <- d$raw_pvalue %||% d$P.Value %||% d$p_value; d_lfc <- d$logFC %||% d$log2FC
    out$c_da_direct$matched <- sum(!is.na(m)); out$c_da_direct$max_abs_dlogFC <- max(abs(long_b[[lfc_col]] - d_lfc[m]), na.rm = TRUE); out$c_da_direct$max_abs_dp <- if (is.null(d_p)) NA_real_ else max(abs(long_b[[p_col]] - d_p[m]), na.rm = TRUE) } else out$c_da_direct$error <- direct_da$error
# (d) covariates vs correction DA (descriptive)
ka <- paste(long_a[[protein_col]], long_a$comparison); kb <- paste(long_b[[protein_col]], long_b$comparison); m <- match(kb, ka)
sig <- function(t) paste(t[[protein_col]], t$comparison)[!is.na(t[[q_col]]) & t[[q_col]] < 0.05]
out$d_modes_compared <- list(rows = c(covariates = nrow(long_a), correction = nrow(long_b)), matched = sum(!is.na(m)), max_abs_dlogFC = max(abs(long_b[[lfc_col]] - long_a[[lfc_col]][m]), na.rm = TRUE),
    max_abs_dp = max(abs(long_b[[p_col]] - long_a[[p_col]][m]), na.rm = TRUE), significant = c(covariates = length(sig(long_a)), correction = length(sig(long_b)), common = length(intersect(sig(long_a), sig(long_b)))),
    formula = c(covariates = ma$model$formula_used, correction = mb$model$formula_used))
jsonlite::write_json(out, file.path("docs/qa/workflow-execution/kit/kit-056", paste0("mode-check-screen-", label, ".json")), auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null", na = "null")
cat("MODE CHECK", label, "| same sweep", out$a_same_sweep$k$identical && out$a_same_sweep$percentage$identical, "| matrix == RUVIII newY", out$b_correction_matrix$state_equals_ruv_RUVIII_newY,
    "| E == matrix", out$b_correction_matrix$owned_E_equals_state_matrix, "| SE identical", out$b_correction_matrix$standard_errors_identical_to_source, "| RUVIIIC max |d|", signif(out$b_correction_matrix$max_abs_delta_vs_RUVIIIC, 3),
    "| DA direct max |dlogFC|", signif(out$c_da_direct$max_abs_dlogFC %||% NA, 3), "| modes max |dlogFC|", signif(out$d_modes_compared$max_abs_dlogFC, 3), "| sig", paste(out$d_modes_compared$significant, collapse = "/"), "\n")
