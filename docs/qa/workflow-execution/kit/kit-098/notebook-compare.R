# KIT-098 notebook comparator (kit-098/plan.md section 6; protocol v1.8 section 10 item 2): a journey member's store against
# the second notebook reference (the unchanged notebook on design_matrix.tab in the counts file's run order), `identical()`
# throughout, no alignment and no tolerance. For the technical-replicate journey (098b) the RUVIII product is compared with
# the direct ruv::RUVIII call of KIT-093 on the reference's inputs. Reads stores only: the member's project is reflinked to a
# scratch directory and opened there, so the member's runtime directory is never written. Every mapped product is reported
# as compared (identical or not) or excluded with its reason.
# Usage: R_LIBS=<the member's build library> Rscript docs/qa/workflow-execution/kit/kit-098/notebook-compare.R <ticket 098a|098b>
#   <member dir (holds journey-result.rds) or a recorder's oracle rds> <scratch dir> <output prefix>
suppressMessages(library(MultiScholaR))
args <- commandArgs(TRUE)
ticket <- args[[1L]]; source_path <- normalizePath(args[[2L]]); scratch <- args[[3L]]; out <- args[[4L]]
repo <- normalizePath(".")
`%||%` <- function(a, b) if (is.null(a)) b else a
ns <- asNamespace("MultiScholaR"); api <- function(name) get(name, ns)
bindings <- jsonlite::read_json(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
spec <- bindings$tickets[[ticket]]
reference_spec <- spec$notebook_reference
stopifnot(!is.null(reference_spec$path), identical(digest::digest(file = reference_spec$path, algo = "sha256"), reference_spec$sha256))
reference <- readRDS(reference_spec$path)$products
technical <- identical(spec$driver_settings$ruviii_m_source, "technical_replicate")
# The store: a member's (journey-result.rds names it) or, for the comparator's dry run, a recorder's (its record names it).
recorded <- if (dir.exists(source_path)) readRDS(file.path(source_path, "journey-result.rds")) else readRDS(source_path)
store_dir <- recorded$store
stopifnot(is.character(store_dir), dir.exists(store_dir))
project <- sub("^(.*/project)/.*$", "\\1", store_dir)
stopifnot(startsWith(store_dir, project))
stopifnot(!dir.exists(scratch))
dir.create(scratch, recursive = TRUE)
stopifnot(identical(system2("cp", c("-a", "--reflink=auto", shQuote(project), shQuote(file.path(scratch, "project")))), 0L))
copy <- file.path(scratch, "project", substring(store_dir, nchar(project) + 2L))
identity <- api("readSavedStoreIdentity")(file.path(copy, "store.json"))
owner <- api("newDIAStateAdapter")(api("newWorkflowEntityStore")(dirname(copy), identity$project_id, identity$workflow_id), reopen = TRUE)
on.exit(api("closeDIAStateAdapter")(owner), add = TRUE)
store <- owner$.context$store
scoped <- function(f) api("withWorkflowStateScope")(owner, f)
state_inputs <- function(name) scoped(function() api("trxAssayInputs")(api("readTrxWorkflowState")(owner, name)))
hydrate <- function(reference) api("hydrateSavedWorkflowEntity")(store, reference)
products <- lapply(owner$.context$snapshot$products, hydrate)
by_type <- function(type) {
    hit <- Filter(function(p) identical(p$product_type, type), products)
    if (length(hit) != 1L) stop(sprintf("%d products of type %s in the store (one expected)", length(hit), type))
    hit[[1L]]$value
}
matrix_of <- function(table) scoped(function() api("trxSavedMatrix")(owner, table))
filter <- by_type("trx.expression_filter_result")
logcpm <- by_type("trx.tmm_logcpm")
grid <- by_type("trx.ruv_grid")
ruvs <- by_type("trx.ruvs_result")
ruviii <- by_type("trx.ruviii")
manifest <- scoped(function() api("resolveWorkflowStateValue")(owner$.context, owner$.context$snapshot$da_results))
contrast <- manifest$results[[1L]]
applied <- hydrate(contrast$products$da_contrast)$value$assays[[1L]]$results
fit <- hydrate(contrast$assays[[1L]]$fit)$fit
members <- grid$grid_members
pct <- vapply(members, function(r) paste0("pct_", r$percentage), "")
grid_members <- do.call(rbind, lapply(seq(5, 30, by = 5), function(p) {
    ok <- stats::setNames(members, vapply(members, function(r) as.character(r$percentage), ""))[[as.character(p)]]
    if (is.null(ok)) return(data.frame(percentage = p, status = "failed", n_controls = NA_real_, best_k = NA_integer_,
        separation_score = NA_real_, k_penalty = NA_real_, composite_score = NA_real_))
    data.frame(percentage = ok$percentage, status = "complete", n_controls = ok$n_controls, best_k = ok$best_k,
        separation_score = ok$separation_score, k_penalty = ok$k_penalty, composite_score = ok$composite_score)
}))
journey <- list(
    count_matrix_filtered = list(state_inputs("trx_expression_filtered")$counts, "state trx_expression_filtered, trxAssayInputs()$counts"),
    filter_effective_lib_size = list(filter$effective_lib_size, "trx.expression_filter_result $effective_lib_size"),
    filter_min_sample_size = list(filter$min_sample_size, "trx.expression_filter_result $min_sample_size"),
    filter_median_lib_size = list(filter$median_lib_size, "trx.expression_filter_result $median_lib_size"),
    dge_filtered_samples = list(logcpm$samples, "trx.tmm_logcpm $samples"),
    log_cpm_tmm = list(matrix_of(logcpm$log_cpm), "trx.tmm_logcpm $log_cpm"),
    ordered_features_idx = list(grid$prelim_order, "trx.ruv_grid $prelim_order"),
    max_acceptable_k = list(grid$max_acceptable_k, "trx.ruv_grid $max_acceptable_k"),
    grid_members = list(grid_members, "trx.ruv_grid $grid_members (KIT-095's shape)"),
    grid_cancor_data = list(stats::setNames(lapply(members, `[[`, "cancor_data"), pct), "trx.ruv_grid members' cancor_data"),
    grid_control_genes = list(stats::setNames(lapply(members, `[[`, "control_genes"), pct), "trx.ruv_grid members' control_genes"),
    df_summary = list(grid$grid_summary, "trx.ruv_grid $grid_summary"),
    best_res = list(grid$best[c("percentage", "n_controls", "best_k", "separation_score", "k_penalty", "composite_score", "ctl_idx", "control_genes")],
        "trx.ruv_grid $best (the notebook's fields)"),
    best_control_genes = list(grid$control_genes, "trx.ruv_grid $control_genes"),
    best_k = list(ruviii$k, "trx.ruviii $k"),
    w_matrix = list(ruvs$W, "trx.ruvs_result $W"),
    ruvs_counts = list(state_inputs("trx_ruvs")$counts, "state trx_ruvs, trxAssayInputs()$counts"),
    log_cpm_corrected = list(matrix_of(ruvs$corrected), "trx.ruvs_result $corrected"),
    design_final = list(fit$design, "DA fit $design"),
    fit_final = list(fit, "DA fit"),
    design_formula_str = list(manifest$parameters$formula, "DA manifest $parameters$formula"),
    contrast_entry = list(names(manifest$results)[[1L]], "DA manifest, first result's name"),
    contrasts_matrix = list(applied$contrasts_matrix, "first contrast's applied product $contrasts_matrix"),
    qlf = list(applied$qlf, "first contrast's applied product $qlf"),
    qlf_table = list(applied$qlf$table, "first contrast's applied product $qlf$table"),
    top_genes = list(applied$top_table, "first contrast's applied product $top_table"),
    ruviii_ctl_mask = list(ruviii$ctl_mask, "trx.ruviii $ctl_mask"))
rows <- list()
row <- function(product, source, compared, identical_, note = "") rows[[length(rows) + 1L]] <<- data.frame(product = product,
    source = source, status = if (compared) "compared" else "excluded", identical = if (compared) identical_ else NA, note = note, stringsAsFactors = FALSE)
for (name in names(journey)) row(name, journey[[name]][[2L]], TRUE, identical(journey[[name]][[1L]], reference[[name]]))
stage <- function(name, value, expected) row(name, "trx.expression_filter_result $counts", TRUE, identical(as.integer(value), as.integer(expected)))
stage("filter_count_before", filter$counts$before, length(reference$keep_expr))
stage("filter_count_after_filter_by_expr", filter$counts$after_filter_by_expr, sum(reference$keep_expr))
stage("filter_count_after_zero_variance", filter$counts$after_zero_variance, sum(reference$gene_vars_positive))
if (!technical) {
    row("ruviii_newY", "trx.ruviii $newY", TRUE, identical(ruviii$newY, reference$ruviii_newY))
    row("ruviii_M", "trx.ruviii $M", TRUE, identical(ruviii$M, reference$ruviii_M))
} else {
    # KIT-093's direct oracle on the reference's inputs; M from the journey design's technical-replicate sets, in run order,
    # whose partition must be design.tsv's.
    x <- state_inputs("trx_ruvs")
    object <- scoped(function() api("readTrxWorkflowState")(owner, "trx_ruvs"))
    sets <- as.character(x$design[[object@technical_replicate_id]])
    corpus <- bindings$corpora[[spec$corpus]]
    tsv <- utils::read.delim(corpus$files[[corpus$design_file]]$path, stringsAsFactors = FALSE)
    truth <- tsv$tech_rep_group[match(x$design$Run, tsv$sample)]
    partition <- function(v) unname(sort(vapply(split(seq_along(v), v), function(i) paste(i, collapse = ","), "")))
    row("technical_replicate_sets", "design state's technical replicate column", TRUE, identical(partition(sets), partition(truth)),
        "the partition of runs equals design.tsv's tech_rep_group")
    M <- ruv::replicate.matrix(sets)
    direct <- ruv::RUVIII(Y = t(reference$log_cpm_tmm), M = M, ctl = reference$ruviii_ctl_mask, k = reference$best_k, return.info = TRUE)$newY
    row("ruviii_newY", "trx.ruviii $newY", TRUE, identical(ruviii$newY, direct), "against the direct ruv::RUVIII call (KIT-093), not the notebook")
    row("ruviii_M", "trx.ruviii $M", TRUE, identical(ruviii$M, M), "against ruv::replicate.matrix(technical replicate sets)")
    row("ruviii_m_source", "trx.ruviii $m_source", TRUE, identical(ruviii$m_source, "technical_replicate"))
}
excluded <- c(keep_expr = "not stored by the journey; the kept gene set is compared exactly through count_matrix_filtered's row names and the stage counts",
    gene_vars = "not stored by the journey (as keep_expr)", gene_vars_positive = "not stored by the journey (as keep_expr)",
    log_cpm_raw = "not stored (KIT-095: trx_qc_figures computes it for stage raw)",
    tt_prelim = "a notebook intermediate the journey does not produce", dge_raw_samples = "a notebook intermediate the journey does not produce",
    ruvs_counts_identical = "a notebook intermediate (KIT-093: RUVs counts equal the filtered counts; ruvs_counts compared)",
    glimma_html_exists = "a file check of the notebook's render", feature_meta = "not a head product",
    top_genes_annotated = "a publication table, not a head product", top_up = "a publication table, not a head product",
    top_down = "a publication table, not a head product", top_de_combined = "a publication table, not a head product")
for (name in names(excluded)) row(name, "-", FALSE, NA, excluded[[name]])
stats <- grep("^(pca|rle)_|^(cca|norm_grid|ruviii_grid|volcano|bar|glimma_anno|summary_text)", names(reference), value = TRUE)
for (name in stats) row(name, "-", FALSE, NA, "renderer and PCA/RLE statistics outputs: a display value outside the head's products; their input matrices are compared")
table <- do.call(rbind, rows)
missing <- setdiff(names(reference), table$product)
for (name in missing) row(name, "-", FALSE, NA, "not mapped")
table <- do.call(rbind, rows)
utils::write.table(table, paste0(out, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
compared <- table[table$status == "compared", ]
lines <- c(sprintf("KIT-098 notebook comparison, ticket %s, store %s", ticket, store_dir),
    sprintf("reference %s (sha256 %s)", reference_spec$path, reference_spec$sha256),
    sprintf("compared %d: identical %d, not identical %d (%s); excluded with a reason %d; unmapped reference products %d",
        nrow(compared), sum(compared$identical), sum(!compared$identical), paste(compared$product[!compared$identical], collapse = ", "),
        sum(table$status == "excluded" & table$note != "not mapped"), sum(table$note == "not mapped")),
    sprintf("NOTEBOOK %s", if (all(compared$identical) && !any(table$note == "not mapped")) "identical" else "mismatch"))
writeLines(lines, paste0(out, ".txt"))
cat(lines, sep = "\n")
