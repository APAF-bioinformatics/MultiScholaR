# KIT-217 control S: the lipidomics and metabolomics intensity filter's per-group count, base against head, on the S4
# methods and the declared step (ticket KIT-217 r1.6 section 5.8, acceptance criterion A4). Runs in R2 on the fix commit.
#
# Usage, from the package root of the head commit on a run-tests runner (full git history, `fetch_refs:
# ["kit-217-measure"]`):
#   Rscript docs/qa/workflow-execution/kit/kit-217/filter-missing-control.R <base commit> [<out dir>]
#   <out dir>  default docs/qa/workflow-execution/kit/kit-217/s (made absolute)
#
# KIT-180's mechanism (as KIT-202 S1): `git worktree add <tmp> <base>`, then per side a fresh Rscript process (this
# script, `--side`) whose working directory is that tree's root, with
# pkgload::load_all(<tree>, export_all = TRUE, helpers = TRUE). Each side writes <out>/<side>-rows.tsv; the comparison
# is made here.
#
# Inputs (written once by the parent to <out>/inputs.rds, read by both sides):
#   complete  the four lipid/metab GUI-entry lanes' design objects (KIT-202's import path, copied from KIT-202 S1,
#             docs/qa/workflow-execution/kit/kit-202/filter-samples-control.R at 185ebf75, function `run_input`, built
#             in each side) and M's tables of the mask `none` for cells L2 and M2, replicates 1 to 5;
#   sparse    M's tables of the masks vd10 and vd25 for L2 and M2, replicates 1 to 5;
#   module-CI the module-CI lipid and metabolite QC objects (each assay holds one missing value, so they are complete
#             only at g = 1; ticket r1.6): in (i) at g = 1, in (ii) at g = 1 and 2.
#   M's tables come from KIT-217's measurement generator, docs/qa/workflow-execution/kit/kit-217/m/sim-k217.R on branch
#   kit-217-measure at 30f57029 (read with `git show`; its functions GROUPS, cell, MED_VAR_PROT, MED_VAR_OMIC, maskMech,
#   maskRate, scale1, generateReplicate, solveLod, maskMatrix and CELLS evaluated, nothing else), with strength 0 as M
#   ran them; raw scale 2^y, missing values NA, one table per assay as M's omicObject builds them.
# Criteria (written to <out>/compare.tsv):
#   (i)   complete inputs: each S4 method at percentiles {1, 10, 20, 50} x the slider's 19 proportions x g in {1, 2}
#         (module-CI: g = 1): head equal to base (kept-id digest, filtered-table digest, rows kept);
#   (ii)  sparse inputs and module-CI at percentiles {1, 10} x proportions {0.2, 0.5}: head equal to the oracle at m = 2
#         and g = 1 and 2 (L2, module-CI) or g = 1 to 4 (M2); base equal to the oracle of today's rule (C). The oracle is
#         written here from the rule (ticket 5.1 and L-1 to L-3), not from package code;
#   (iii) the declared step at the defaults (10, 0.5, m 2, g 1), `whole` and `partitioned`, one owner per variant, on
#         every input of (i) and (ii): head `whole` equal to head `partitioned` (kept-id digest, table digest, rows kept,
#         rows removed for missing values and for low intensity);
#   (iv)  the head under `MultiScholaR.intensity_filter_missing = "not_counted"` equal to the base on every row of (ii).
# Written: <out>/filter-identity.tsv (every side's rows) and <out>/compare.tsv (criterion, input, assay, setting, equal).
# Exit status non-zero when a side fails or any criterion row is not equal.

args <- commandArgs(trailingOnly = TRUE)

# ---- Child: one side --------------------------------------------------------------------------------------------------
childMain <- function(side, out) {
    inputs <- readRDS(file.path(out, "inputs.rds"))
    e2e <- normalizePath(file.path("tests", "testdata", "e2e"), mustWork = TRUE)
    manifest <- jsonlite::fromJSON(file.path(e2e, "manifest.json"), simplifyVector = FALSE)
    lane_of <- function(id) {
        for (lane in manifest$lanes) if (identical(lane$lane_id, id)) return(lane)
        stop("lane not in the manifest: ", id)
    }
    run_input <- function(omic, vendor_format, assay1_file, assay1_name, assay2_file = NULL, assay2_name = NULL,
                          custom = NULL, sanitize = TRUE, design = NULL) {
        config <- omicGuiConfig(omic)
        import <- config$import
        preview <- config$import$preview(assay1_file, import = import)
        headers <- names(preview$assay_data)
        if (identical(vendor_format, "custom")) {
            id_col <- resolveOmicsImportColumn(headers, vendor_format, NULL, custom$id)
            annotation_col <- resolveOmicsImportColumn(headers, vendor_format, NULL, custom$annotation)
            sample_pattern <- custom$samples
            is_pattern <- custom$is_pattern
        } else {
            id_col <- preview$updates$id$selected
            annotation_col <- preview$updates$annotation$selected
            sample_pattern <- ""
            is_pattern <- preview$updates$is_pattern %||% ""
        }
        exclude <- if (identical(omic, "lipidomics")) OMICS_IMPORT_NORMALISED_COLUMNS
        sample_columns <- resolveOmicsImportSampleColumns(preview$assay_data, vendor_format, sample_pattern,
            preview$import_result, exclude)
        source_dir <- file.path(tempfile("kit217-s-"), "source")
        dir.create(source_dir, recursive = TRUE)
        request <- list(assay1_name = assay1_name, assay1_data = preview$assay_data, assay1_file = assay1_file,
            assay2_file = assay2_file, assay2_name = assay2_name, vendor_format = vendor_format,
            detected_format = preview$detected_format, id_col = id_col, annotation_col = annotation_col,
            sample_columns = sample_columns, sanitize_names = sanitize, is_pattern = is_pattern,
            experiment_paths = list(source_dir = source_dir))
        imported <- importOmicsTables(request, config)
        data_cln <- buildOmicsDesignSaveDataList(imported$assays, design$sample, imported$column_mapping,
            import$id_mapping_key)
        design_matrix <- data.frame(Run = design$sample, group = design$group,
            tech_rep_group = paste(design$group, ave(design$sample, design$group, FUN = seq_along), sep = "_"),
            stringsAsFactors = FALSE)
        config$design$build_object(list(data_cln = data_cln, design_matrix = design_matrix, config_list = list()),
            imported$column_mapping)
    }
    custom_inputs <- list(id = "Feature.Name", annotation = "", samples = "^(WT|KO)_", is_pattern = "^ITSD|_d[0-9]+$")
    lane_input <- function(lane_id) {
        lane <- lane_of(lane_id)
        dir <- file.path(e2e, lane$fixture_dir)
        assays <- unlist(lane$assays)
        design <- utils::read.delim(file.path(dir, "design.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
        second <- !is.null(lane$assay2_file) && length(assays) >= 2L
        run_input(lane$omic_type, lane$import_tool, file.path(dir, lane$seed_file), assays[[1L]],
            assay2_file = if (second) file.path(dir, lane$assay2_file), assay2_name = if (second) assays[[2L]],
            custom = if (identical(lane$import_tool, "custom")) custom_inputs, sanitize = FALSE, design = design)
    }
    # M's object (sim-k217.R `assayTable` and `omicObject`, copied): one table per assay, id column `feature_id`.
    m_object <- function(item) {
        raw <- item$raw
        tables <- lapply(split(rownames(raw), factor(item$assay[rownames(raw)], levels = unique(item$assay))), function(ids) {
            table <- tibble::as_tibble(as.data.frame(raw[ids, , drop = FALSE], check.names = FALSE))
            table <- tibble::add_column(table, feature_id = ids, .before = 1L)
            tibble::add_column(table, name = ids, .after = 1L)
        })
        if (identical(item$omic, "lipid")) {
            methods::new("LipidomicsAssayData", lipid_data = tables, lipid_id_column = "feature_id", annotation_id_column = "name",
                design_matrix = item$design, sample_id = "Sample_ID", group_id = "group", technical_replicate_id = "extract",
                args = list())
        } else {
            methods::new("MetaboliteAssayData", metabolite_data = tables, metabolite_id_column = "feature_id",
                annotation_id_column = "name", design_matrix = item$design, sample_id = "Sample_ID", group_id = "group",
                technical_replicate_id = "extract", args = list())
        }
    }
    quiet <- function(expr) suppressWarnings(suppressMessages(expr))
    sha <- function(x) digest::digest(x, algo = "sha256")
    is_lipid <- function(object) methods::is(object, "LipidomicsAssayData")
    binding_of <- function(object) if (is_lipid(object)) lipidomicsDeclarationBinding else metabolomicsDeclarationBinding
    names_of <- function(object) {
        prefix <- if (is_lipid(object)) "lipids" else "metabolites"
        list(percentile = paste0(prefix, "_intensity_cutoff_percentile"),
            proportion = paste0(prefix, "_proportion_of_samples_below_cutoff"), m = paste0(prefix, "_min_values_per_group"),
            g = paste0(prefix, "_min_groups"))
    }
    has_counts <- function(object) names_of(object)$g %in% names(formals(if (is_lipid(object)) lipidIntensityFiltering else
        metaboliteIntensityFiltering))
    s4 <- function(object, percentile, proportion, g = NULL) {
        n <- names_of(object)
        call_args <- stats::setNames(list(object, percentile, proportion), c("theObject", n$percentile, n$proportion))
        if (!is.null(g) && has_counts(object)) call_args[[n$g]] <- g
        quiet(do.call(if (is_lipid(object)) "lipidIntensityFiltering" else "metaboliteIntensityFiltering", call_args))
    }
    rows <- list()
    add <- function(...) rows[[length(rows) + 1L]] <<- data.frame(..., stringsAsFactors = FALSE)
    record_object <- function(input, object, filtered, route, variant, percentile, proportion, g, rule, value = NULL) {
        binding <- binding_of(object)
        id <- methods::slot(object, binding$id_slot)
        for (assay in names(methods::slot(filtered, binding$data_slot))) {
            table <- methods::slot(filtered, binding$data_slot)[[assay]]
            add(input = input, omic = binding$prefix, route = route, variant = variant, percentile = percentile,
                proportion = proportion, g = g, side = side, rule = rule, assay = assay,
                kept_digest = sha(as.character(table[[id]])), table_digest = sha(as.data.frame(table)), n_kept = nrow(table),
                removed_missing = if (is.null(value$removed_missing)) NA_integer_ else value$removed_missing[[assay]],
                removed_low = if (is.null(value$removed_low)) NA_integer_ else value$removed_low[[assay]])
        }
    }
    # The oracle (ticket 5.1, L-1 to L-3), written from the definition: per assay the type-7 percentile of the design
    # samples' values (missing removed); below = a value strictly below it; share = below / samples; present = finite;
    # per group (labels not missing, not blank, not named like "pool") need = min(group size, m); a feature is kept when
    # its share is under the proportion and, for rule A, at least g groups have `need` present values.
    oracle <- function(input, object, percentiles, proportions, gs, rule) {
        binding <- binding_of(object)
        id <- methods::slot(object, binding$id_slot)
        design <- object@design_matrix
        samples <- as.character(design[[object@sample_id]])
        labels <- as.character(design[[object@group_id]])
        counted <- !is.na(labels) & nzchar(trimws(labels)) & !grepl("pool", labels, ignore.case = TRUE)
        for (assay in names(methods::slot(object, binding$data_slot))) {
            table <- as.data.frame(methods::slot(object, binding$data_slot)[[assay]])
            columns <- intersect(names(table), samples)
            x <- as.matrix(table[, columns, drop = FALSE])
            group_of <- labels[match(columns, samples)]
            in_count <- counted[match(columns, samples)]
            for (percentile in percentiles) for (proportion in proportions) for (g in gs) {
                tau <- unname(stats::quantile(as.vector(x), probs = percentile / 100, na.rm = TRUE, type = 7))
                below <- !is.na(x) & x < tau
                share_ok <- rowSums(below) / ncol(x) < proportion
                keep <- share_ok
                if (identical(rule, "oracle_A")) {
                    met <- rep(0L, nrow(x))
                    for (grp in unique(group_of[in_count])) {
                        cols <- which(in_count & group_of == grp)
                        met <- met + as.integer(rowSums(is.finite(x[, cols, drop = FALSE])) >= min(length(cols), 2L))
                    }
                    keep <- share_ok & met >= g
                }
                add(input = input, omic = binding$prefix, route = "oracle", variant = "", percentile = percentile,
                    proportion = proportion, g = g, side = side, rule = rule, assay = assay,
                    kept_digest = sha(as.character(table[[id]][keep])), table_digest = "", n_kept = sum(keep),
                    removed_missing = NA_integer_, removed_low = NA_integer_)
            }
        }
    }
    step <- function(input, object, variant) {
        binding <- binding_of(object)
        root <- file.path(tempfile("kit217-s-owner-"))
        dir.create(root)
        fixture <- if (is_lipid(object)) lipidDeclarationsOwner(root, object) else metabDeclarationsOwner(root, object)
        on.exit(closeDIAStateAdapter(fixture$owner), add = TRUE)
        filter_step <- omicsIntensityFilterStep(binding, 10, 0.5)
        run <- quiet(runWorkflowOmicsReferenceOperation(fixture$owner, workflowOperationEntry(filter_step$operation)$declaration,
            c(filter_step$parameters, list(variant = variant)), filter_step$state_name, fixture$workflow$config_list,
            filter_step$description))
        filtered <- fixture$owner$getState(filter_step$state_name)
        record_object(input, object, filtered, "declared_step", variant, 10, 0.5, 1L, "shipped", run$value)
    }
    percentiles_i <- c(1, 10, 20, 50)
    proportions_i <- round(seq(0.1, 1, by = 0.05), 2)
    failed <- FALSE
    builders <- c(
        list(lipid_canonical = function() lane_input("lipid_canonical"), metab_lc = function() lane_input("metab_lc"),
            metab_gc = function() lane_input("metab_gc"), metab_combined = function() lane_input("metab_combined"),
            module_ci_lipid = function() module_ci_lipid_qc_object(), module_ci_metab = function() module_ci_metab_qc_object()),
        lapply(inputs, function(item) function() m_object(item)))
    kind_of <- function(input) if (input %in% c("module_ci_lipid", "module_ci_metab")) "module_ci" else
        if (input %in% names(inputs)) inputs[[input]]$kind else "complete"
    for (input in names(builders)) {
        cat("KIT-217 S", side, "input", input, "\n")
        ok <- tryCatch({
            object <- quiet(builders[[input]]())
            kind <- kind_of(input)
            groups_in_design <- length(unique(object@design_matrix[[object@group_id]]))
            # (i): complete inputs (module-CI at g = 1), every setting
            if (kind %in% c("complete", "module_ci")) {
                gs <- if (identical(kind, "module_ci")) 1L else 1:2
                for (percentile in percentiles_i) for (proportion in proportions_i) {
                    if (identical(side, "base")) {
                        record_object(input, object, s4(object, percentile, proportion), "s4_i", "", percentile, proportion,
                            NA_integer_, "shipped")
                    } else for (g in gs) {
                        record_object(input, object, s4(object, percentile, proportion, g), "s4_i", "", percentile, proportion,
                            g, "shipped")
                    }
                }
            }
            # (ii) and (iv): sparse inputs and module-CI
            if (kind %in% c("sparse", "module_ci")) {
                gs <- if (groups_in_design > 2L) seq_len(groups_in_design) else 1:2
                for (percentile in c(1, 10)) for (proportion in c(0.2, 0.5)) {
                    if (identical(side, "base")) {
                        record_object(input, object, s4(object, percentile, proportion), "s4_ii", "", percentile, proportion,
                            NA_integer_, "shipped")
                    } else {
                        for (g in gs) record_object(input, object, s4(object, percentile, proportion, g), "s4_ii", "",
                            percentile, proportion, g, "shipped")
                        withr::with_options(list(MultiScholaR.intensity_filter_missing = "not_counted"),
                            record_object(input, object, s4(object, percentile, proportion), "s4_ii", "", percentile,
                                proportion, NA_integer_, "earlier_option"))
                    }
                }
                if (identical(side, "head")) {
                    oracle(input, object, c(1, 10), c(0.2, 0.5), gs, "oracle_A")
                    oracle(input, object, c(1, 10), c(0.2, 0.5), NA_integer_, "oracle_C")
                }
            }
            # (iii): the declared step at the defaults, both variants (head)
            if (identical(side, "head")) for (variant in c("whole", "partitioned")) step(input, object, variant)
            TRUE
        }, error = function(e) {
            cat("KIT-217 S", side, "input", input, "ERROR", class(e)[[1L]], conditionMessage(e), "\n")
            FALSE
        })
        if (!ok) failed <- TRUE
    }
    utils::write.table(do.call(rbind, rows), file.path(out, paste0(side, "-rows.tsv")), sep = "\t", row.names = FALSE,
        quote = FALSE)
    quit(status = if (failed) 1L else 0L)
}

if (length(args) >= 1L && identical(args[[1L]], "--side")) {
    setwd(args[[2L]])
    suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = TRUE))
    environment(childMain) <- new.env(parent = asNamespace("MultiScholaR"))
    childMain(args[[3L]], args[[4L]])
}

# ---- Parent: the inputs, both sides, then the comparison ---------------------------------------------------------------
stopifnot(length(args) >= 1L, grepl("^[0-9a-f]{40}$", args[[1L]]))
base_commit <- args[[1L]]
out <- if (length(args) >= 2L) args[[2L]] else "docs/qa/workflow-execution/kit/kit-217/s"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
head_root <- normalizePath(".")
head_commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
script <- normalizePath(file.path("docs", "qa", "workflow-execution", "kit", "kit-217", "filter-missing-control.R"))
cat("KIT-217 S: head", head_commit, "base", base_commit, "out", out, "\n")

# M's generator, by name, from the measurement branch (nothing else of the file is evaluated).
generator <- new.env(parent = globalenv())
sim_text <- system2("git", c("show", "origin/kit-217-measure:docs/qa/workflow-execution/kit/kit-217/m/sim-k217.R"),
    stdout = TRUE)
stopifnot(length(sim_text) > 100L)
wanted <- c("GROUPS", "cell", "MED_VAR_PROT", "MED_VAR_OMIC", "maskMech", "maskRate", "scale1", "generateReplicate",
    "solveLod", "maskMatrix", "CELLS")
for (expr in parse(text = sim_text, keep.source = FALSE)) {
    if (is.call(expr) && as.character(expr[[1L]]) %in% c("<-", "=") && is.symbol(expr[[2L]]) &&
        as.character(expr[[2L]]) %in% wanted) eval(expr, envir = generator)
}
stopifnot(all(wanted %in% ls(generator)))
m_inputs <- list()
for (cellname in c("L2", "M2")) for (r in 1:5) {
    sim <- generator$generateReplicate(cellname, r, strength = 0)
    for (mask in c("none", "vd10", "vd25")) {
        masked <- sim$y
        masked[generator$maskMatrix(sim, mask)] <- NA_real_
        m_inputs[[sprintf("%s_%s_r%d", cellname, mask, r)]] <- list(omic = sim$p$omic, raw = 2^masked, assay = sim$assay,
            design = sim$design, kind = if (identical(mask, "none")) "complete" else "sparse")
    }
}
saveRDS(m_inputs, file.path(out, "inputs.rds"))
cat("KIT-217 S: M inputs", length(m_inputs), "\n")

rscript <- file.path(R.home("bin"), "Rscript")
runSide <- function(tree, side) {
    console <- suppressWarnings(system2(rscript, c(shQuote(script), "--side", shQuote(tree), side, shQuote(out)),
        stdout = TRUE, stderr = TRUE))
    writeLines(console, file.path(out, paste0(side, "-console.txt")))
    cat(grep("^KIT-217 S", console, value = TRUE), sep = "\n")
    cat(utils::tail(console, 10L), sep = "\n")
    attr(console, "status") %||% 0L
}
base_tree <- file.path(tempdir(), "kit217-s-base")
invisible(system2("git", c("worktree", "add", "--detach", shQuote(base_tree), base_commit), stdout = TRUE, stderr = TRUE))
stopifnot(dir.exists(file.path(base_tree, "tests", "testthat")))
status <- c(base = runSide(base_tree, "base"), head = runSide(head_root, "head"))

read_rows <- function(side) {
    path <- file.path(out, paste0(side, "-rows.tsv"))
    if (!file.exists(path)) return(NULL)
    utils::read.delim(path, colClasses = "character", na.strings = "NA")
}
rows <- rbind(read_rows("base"), read_rows("head"))
utils::write.table(rows, file.path(out, "filter-identity.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

key <- function(x) paste(x$input, x$assay, x$percentile, x$proportion, x$g, x$variant, sep = "|")
pick <- function(route, side, rule) rows[rows$route == route & rows$side == side & rows$rule == rule, , drop = FALSE]
compare <- function(criterion, a, b, fields) {
    keys <- union(key(a), key(b))
    if (!length(keys)) return(data.frame(criterion = criterion, input = "<no rows>", assay = "", setting = "", equal = FALSE,
        stringsAsFactors = FALSE))
    a <- a[match(keys, key(a)), , drop = FALSE]
    b <- b[match(keys, key(b)), , drop = FALSE]
    equal <- !is.na(a$input) & !is.na(b$input)
    for (field in fields) equal <- equal & ((a[[field]] == b[[field]]) %in% TRUE | (is.na(a[[field]]) & is.na(b[[field]])))
    parts <- strsplit(keys, "|", fixed = TRUE)
    data.frame(criterion = criterion, input = vapply(parts, `[`, "", 1L), assay = vapply(parts, `[`, "", 2L),
        setting = vapply(parts, function(p) paste(p[3:6], collapse = " "), ""), equal = equal, stringsAsFactors = FALSE)
}
# The base ran (i) and (ii) once per setting (no g); its row stands for every g of the head.
with_g <- function(base, head) {
    expanded <- lapply(seq_len(nrow(head)), function(i) {
        match_row <- base[base$input == head$input[[i]] & base$assay == head$assay[[i]] &
            base$percentile == head$percentile[[i]] & base$proportion == head$proportion[[i]], , drop = FALSE]
        if (!nrow(match_row)) return(NULL)
        match_row$g <- head$g[[i]]
        match_row[1L, , drop = FALSE]
    })
    do.call(rbind, expanded)
}
drop_variant <- function(x) { x$variant <- ""; x }
step_rows <- function(variant) drop_variant(rows[rows$route == "declared_step" & rows$side == "head" & rows$variant == variant, ,
    drop = FALSE])
drop_g <- function(x) { x$g <- NA_character_; x }
head_i <- pick("s4_i", "head", "shipped")
head_ii <- pick("s4_ii", "head", "shipped")
result <- rbind(
    compare("i_base_equals_head", with_g(pick("s4_i", "base", "shipped"), head_i), head_i,
        c("kept_digest", "table_digest", "n_kept")),
    compare("ii_head_equals_oracle_A", head_ii, pick("oracle", "head", "oracle_A"), c("kept_digest", "n_kept")),
    compare("ii_base_equals_oracle_C", drop_g(pick("s4_ii", "base", "shipped")), drop_g(pick("oracle", "head", "oracle_C")),
        c("kept_digest", "n_kept")),
    compare("iii_head_whole_equals_partitioned", step_rows("whole"), step_rows("partitioned"),
        c("kept_digest", "table_digest", "n_kept", "removed_missing", "removed_low")),
    compare("iv_head_earlier_equals_base", drop_g(pick("s4_ii", "head", "earlier_option")), drop_g(pick("s4_ii", "base", "shipped")),
        c("kept_digest", "table_digest", "n_kept")))
utils::write.table(result, file.path(out, "compare.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
summary <- stats::aggregate(equal ~ criterion, result, function(x) paste0(sum(x), "/", length(x)))
print(summary)
cat("side status:", paste(names(status), status, collapse = ", "), "\n")
quit(status = if (any(status != 0L) || !all(result$equal)) 1L else 0L)
