# KIT-202 script S1: the lipidomics and metabolomics intensity filter, base against head, on the S4 methods and the
# declared step (ticket section 9 "S1", acceptance criterion A2). Runs in R2 on the fix commit.
#
# Usage, from the package root of the head commit on a run-tests runner (full git history):
#   Rscript docs/qa/workflow-execution/kit/kit-202/filter-samples-control.R <base commit> [<out dir>]
#   <out dir>  default docs/qa/workflow-execution/kit/kit-202/s1 (made absolute)
#
# KIT-180's mechanism (as KIT-198 S1): `git worktree add <tmp> <base>`, then per side a fresh Rscript process (this
# script, `--side`) whose working directory is that tree's root, with
# pkgload::load_all(<tree>, export_all = TRUE, helpers = TRUE). Each side writes <out>/<side>-rows.tsv; the comparison
# is made here.
# Inputs. (A) Unaffected: the four GUI-entry lanes' design objects (KIT-198 G0's import path, copied below from
#   docs/qa/workflow-execution/kit/kit-198/g0-numeric-columns.R at b0e913b8, function `run_input`), the module-CI lipid
#   and metabolite QC objects. (B) Affected: the two MS-DIAL objects of G0 (metabolite seed; module-CI lipid MS-DIAL
#   table with `Peak ID` = 1:3) and the metabolomics workbook's shape (createMetaboliteAssayData on the seed, every
#   column kept).
# Criteria (written to <out>/compare.tsv):
#   (i)   (A): each S4 method at percentiles {1, 10, 20, 50} x the slider's 19 proportions: base equal to head (kept-id
#         digest, filtered-table digest, threshold);
#   (ii)  (B): the head equal to the oracle on the design's samples, and the base equal to the same oracle on every
#         numeric column (kept-id digest and threshold); the oracle is written here from the rule (ticket 5.3), not from
#         package code: type-7 quantile of the values, missing removed; share = values strictly below it (missing not
#         counted) over the number of columns counted; kept when the share is below the proportion;
#   (iii) the declared step, `whole` and `partitioned`, on every (A) and (B) object at (1, 0.5), (10, 0.5), (10, 0.2),
#         one owner (store root) per variant and per side: head `whole` equal to head `partitioned`; (A) base equal to
#         head;
#   (iv)  the head under `MultiScholaR.intensity_filter_samples = "numeric_columns"` equal to the base on every (B) row
#         of (i).
# Written: <out>/filter-identity.tsv (input, omic, route, variant, percentile, proportion, side, rule, assay, kept
# digest, table digest, threshold, n kept) and <out>/compare.tsv (criterion, input, assay, setting, equal).
# Exit status non-zero when a side fails or any criterion row is not equal.

args <- commandArgs(trailingOnly = TRUE)

# ---- Child: one side --------------------------------------------------------------------------------------------------
childMain <- function(side, out) {
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
        source_dir <- file.path(tempfile("kit202-s1-"), "source")
        dir.create(source_dir, recursive = TRUE)
        request <- list(assay1_name = assay1_name, assay1_data = preview$assay_data, assay1_file = assay1_file,
            assay2_file = assay2_file, assay2_name = assay2_name, vendor_format = vendor_format,
            detected_format = preview$detected_format, id_col = id_col, annotation_col = annotation_col,
            sample_columns = sample_columns, sanitize_names = sanitize, is_pattern = is_pattern,
            experiment_paths = list(source_dir = source_dir))
        imported <- importOmicsTables(request, config)
        if (is.null(design)) {
            id_and_annotation <- c(imported$column_mapping[[import$id_mapping_key]], imported$column_mapping$annotation_col)
            samples <- setdiff(imported$sample_columns, id_and_annotation)
            design <- data.frame(sample = samples, group = sub("_.*$", "", samples), stringsAsFactors = FALSE)
        }
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
    lipid_msdial <- function() {
        assay <- module_ci_lipid_msdial_assay()
        assay[["Peak ID"]] <- seq_len(nrow(assay))
        path <- file.path(tempfile("kit202-s1-lipid-"), "msdial_lipid_pos.csv")
        dir.create(dirname(path), recursive = TRUE)
        utils::write.csv(assay, path, row.names = FALSE)
        run_input("lipidomics", "msdial", path, "LCMS_Pos")
    }
    workbook <- function() {
        seed <- utils::read.csv(file.path(e2e, "metab_lc", "seed_lcms_pos.csv"), check.names = FALSE,
            stringsAsFactors = FALSE)
        design <- data.frame(Run = c("ctrl_1", "ctrl_2", "treat_1", "treat_2"),
            group = c("ctrl", "ctrl", "treat", "treat"), stringsAsFactors = FALSE)
        createMetaboliteAssayData(metabolite_data = list(LCMS_Pos = seed), design_matrix = design,
            metabolite_id_column = "Alignment ID", annotation_id_column = "Name", sample_id = "Run", group_id = "group")
    }
    builders <- list(
        lipid_canonical = function() lane_input("lipid_canonical"), metab_lc = function() lane_input("metab_lc"),
        metab_gc = function() lane_input("metab_gc"), metab_combined = function() lane_input("metab_combined"),
        module_ci_lipid = function() module_ci_lipid_qc_object(), module_ci_metab = function() module_ci_metab_qc_object(),
        lipid_msdial = lipid_msdial,
        metab_msdial = function() run_input("metabolomics", "msdial", file.path(e2e, "metab_lc", "seed_lcms_pos.csv"),
            "LCMS_Pos"),
        workbook = workbook)
    affected <- c("lipid_msdial", "metab_msdial", "workbook")
    quiet <- function(expr) suppressWarnings(suppressMessages(expr))
    sha <- function(x) digest::digest(x, algo = "sha256")
    is_lipid <- function(object) methods::is(object, "LipidomicsAssayData")
    binding_of <- function(object) if (is_lipid(object)) lipidomicsDeclarationBinding else metabolomicsDeclarationBinding
    percentiles <- c(1, 10, 20, 50)
    proportions <- round(seq(0.1, 1, by = 0.05), 2)
    rows <- list()
    add <- function(...) rows[[length(rows) + 1L]] <<- data.frame(..., stringsAsFactors = FALSE)
    failed <- FALSE
    # One filtered object's rows: per assay, the kept ids' digest, the filtered table's digest (as a data frame, so the
    # digest reads values, names and order), the threshold the method's message prints (15 significant digits).
    record_object <- function(input, object, filtered, thresholds, route, variant, percentile, proportion, rule) {
        binding <- binding_of(object)
        id <- methods::slot(object, binding$id_slot)
        for (assay in names(methods::slot(filtered, binding$data_slot))) {
            table <- methods::slot(filtered, binding$data_slot)[[assay]]
            add(input = input, omic = binding$prefix, route = route, variant = variant, percentile = percentile,
                proportion = proportion, side = side, rule = rule, assay = assay,
                kept_digest = sha(as.character(table[[id]])), table_digest = sha(as.data.frame(table)),
                threshold = if (is.na(thresholds[assay])) NA_character_ else thresholds[[assay]], n_kept = nrow(table))
        }
    }
    s4 <- function(object, percentile, proportion) {
        messages <- character()
        filtered <- withCallingHandlers(suppressWarnings(if (is_lipid(object)) {
            lipidIntensityFiltering(object, lipids_intensity_cutoff_percentile = percentile,
                lipids_proportion_of_samples_below_cutoff = proportion)
        } else {
            metaboliteIntensityFiltering(object, metabolites_intensity_cutoff_percentile = percentile,
                metabolites_proportion_of_samples_below_cutoff = proportion)
        }), message = function(m) {
            messages <<- c(messages, conditionMessage(m))
            invokeRestart("muffleMessage")
        })
        lines <- trimws(grep("^Calculated minimum intensity threshold for assay '", messages, value = TRUE))
        list(filtered = filtered, thresholds = stats::setNames(sub(".*': ", "", lines),
            sub("^Calculated minimum intensity threshold for assay '([^']*)'.*$", "\\1", lines)))
    }
    oracle <- function(input, object, columns_of, rule) {
        binding <- binding_of(object)
        id <- methods::slot(object, binding$id_slot)
        for (assay in names(methods::slot(object, binding$data_slot))) {
            table <- as.data.frame(methods::slot(object, binding$data_slot)[[assay]])
            columns <- columns_of(table, object, id)
            values <- as.matrix(table[, columns, drop = FALSE])
            for (percentile in percentiles) for (proportion in proportions) {
                threshold <- unname(stats::quantile(as.vector(values), probs = percentile / 100, na.rm = TRUE, type = 7))
                keep <- rowSums(values < threshold, na.rm = TRUE) / length(columns) < proportion
                add(input = input, omic = binding$prefix, route = "oracle", variant = "", percentile = percentile,
                    proportion = proportion, side = side, rule = rule, assay = assay,
                    kept_digest = sha(as.character(table[[id]][keep])), table_digest = "", threshold = as.character(threshold),
                    n_kept = sum(keep))
            }
        }
    }
    design_columns <- function(table, object, id) {
        setdiff(intersect(names(table), as.character(object@design_matrix[[object@sample_id]])), id)
    }
    numeric_columns <- function(table, object, id) names(table)[vapply(table, is.numeric, logical(1))]
    step <- function(input, object, percentile, proportion, variant) {
        binding <- binding_of(object)
        root <- file.path(tempfile("kit202-s1-owner-"))
        dir.create(root)
        fixture <- if (is_lipid(object)) lipidDeclarationsOwner(root, object) else metabDeclarationsOwner(root, object)
        on.exit(closeDIAStateAdapter(fixture$owner), add = TRUE)
        filter_step <- omicsIntensityFilterStep(binding, percentile, proportion)
        quiet(runWorkflowOmicsReferenceOperation(fixture$owner, workflowOperationEntry(filter_step$operation)$declaration,
            c(filter_step$parameters, list(variant = variant)), filter_step$state_name, fixture$workflow$config_list,
            filter_step$description))
        filtered <- fixture$owner$getState(filter_step$state_name)
        record_object(input, object, filtered, character(), "declared_step", variant, percentile, proportion, "shipped")
    }
    for (input in names(builders)) {
        cat("KIT-202 S1", side, "input", input, "\n")
        ok <- tryCatch({
            object <- quiet(builders[[input]]())
            rules <- if (identical(side, "head") && input %in% affected) c("shipped", "earlier_option") else "shipped"
            for (rule in rules) {
                options <- if (identical(rule, "earlier_option")) list(MultiScholaR.intensity_filter_samples = "numeric_columns") else list()
                withr::with_options(options, for (percentile in percentiles) for (proportion in proportions) {
                    # (Named `res`: `out` is this side's output directory; R2's run 45469464049 overwrote it here.)
                    res <- s4(object, percentile, proportion)
                    record_object(input, object, res$filtered, res$thresholds, "s4", "", percentile, proportion, rule)
                })
            }
            if (identical(side, "head") && input %in% affected) {
                oracle(input, object, design_columns, "oracle_design")
                oracle(input, object, numeric_columns, "oracle_numeric_columns")
            }
            for (setting in list(c(1, 0.5), c(10, 0.5), c(10, 0.2))) for (variant in c("whole", "partitioned")) {
                step(input, object, setting[[1L]], setting[[2L]], variant)
            }
            TRUE
        }, error = function(e) {
            cat("KIT-202 S1", side, "input", input, "ERROR", class(e)[[1L]], conditionMessage(e), "\n")
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

# ---- Parent: both sides, then the comparison ---------------------------------------------------------------------------
stopifnot(length(args) >= 1L, grepl("^[0-9a-f]{40}$", args[[1L]]))
base_commit <- args[[1L]]
out <- if (length(args) >= 2L) args[[2L]] else "docs/qa/workflow-execution/kit/kit-202/s1"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
head_root <- normalizePath(".")
head_commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
script <- normalizePath(file.path("docs", "qa", "workflow-execution", "kit", "kit-202", "filter-samples-control.R"))
cat("KIT-202 S1: head", head_commit, "base", base_commit, "out", out, "\n")
rscript <- file.path(R.home("bin"), "Rscript")
runSide <- function(tree, side) {
    console <- suppressWarnings(system2(rscript, c(shQuote(script), "--side", shQuote(tree), side, shQuote(out)),
        stdout = TRUE, stderr = TRUE))
    writeLines(console, file.path(out, paste0(side, "-console.txt")))
    cat(grep("^KIT-202 S1", console, value = TRUE), sep = "\n")
    cat(utils::tail(console, 10L), sep = "\n")
    attr(console, "status") %||% 0L
}
base_tree <- file.path(tempdir(), "kit202-s1-base")
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

affected <- c("lipid_msdial", "metab_msdial", "workbook")
key <- function(x) paste(x$input, x$assay, x$percentile, x$proportion, x$variant, sep = "|")
pick <- function(route, side, rule, inputs = NULL) {
    x <- rows[rows$route == route & rows$side == side & rows$rule == rule, , drop = FALSE]
    if (!is.null(inputs)) x <- x[x$input %in% inputs, , drop = FALSE]
    x
}
# Thresholds: both absent (the declared step's rows) is equal; otherwise the same text, or with `numeric_threshold`
# the same number to 12 digits (the oracle against the method's message).
same_threshold <- function(a, b, numeric_threshold) {
    both_absent <- is.na(a) & is.na(b)
    same <- if (numeric_threshold) {
        mapply(function(x, y) isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = 1e-12)), a, b)
    } else (a == b) %in% TRUE
    both_absent | (!is.na(a) & !is.na(b) & same)
}
compare <- function(criterion, a, b, fields, numeric_threshold = FALSE) {
    keys <- union(key(a), key(b))
    if (!length(keys)) return(data.frame(criterion = criterion, input = "<no rows>", assay = "", setting = "", equal = FALSE,
        stringsAsFactors = FALSE))
    a <- a[match(keys, key(a)), , drop = FALSE]
    b <- b[match(keys, key(b)), , drop = FALSE]
    equal <- !is.na(a$input) & !is.na(b$input)
    for (field in fields) equal <- equal & (a[[field]] == b[[field]]) %in% TRUE
    equal <- equal & same_threshold(a$threshold, b$threshold, numeric_threshold)
    parts <- strsplit(keys, "|", fixed = TRUE)
    data.frame(criterion = criterion, input = vapply(parts, `[`, "", 1L), assay = vapply(parts, `[`, "", 2L),
        setting = vapply(parts, function(p) paste(p[3:5], collapse = " "), ""), equal = equal, stringsAsFactors = FALSE)
}
unaffected <- setdiff(unique(rows$input), affected)
drop_variant <- function(x) { x$variant <- ""; x }
step_rows <- function(side, variant, inputs = NULL) {
    x <- pick("declared_step", side, "shipped", inputs)
    drop_variant(x[x$variant == variant, , drop = FALSE])
}
result <- rbind(
    compare("i_base_equals_head", pick("s4", "base", "shipped", unaffected), pick("s4", "head", "shipped", unaffected),
        c("kept_digest", "table_digest", "n_kept")),
    compare("ii_head_equals_oracle_design", pick("s4", "head", "shipped", affected),
        drop_variant(pick("oracle", "head", "oracle_design", affected)), c("kept_digest", "n_kept"), numeric_threshold = TRUE),
    compare("ii_base_equals_oracle_numeric_columns", pick("s4", "base", "shipped", affected),
        drop_variant(pick("oracle", "head", "oracle_numeric_columns", affected)), c("kept_digest", "n_kept"),
        numeric_threshold = TRUE),
    compare("iii_head_whole_equals_partitioned", step_rows("head", "whole"), step_rows("head", "partitioned"),
        c("kept_digest", "table_digest", "n_kept")),
    compare("iii_base_whole_equals_partitioned", step_rows("base", "whole"), step_rows("base", "partitioned"),
        c("kept_digest", "table_digest", "n_kept")),
    compare("iii_unaffected_base_equals_head_whole", step_rows("base", "whole", unaffected),
        step_rows("head", "whole", unaffected), c("kept_digest", "table_digest", "n_kept")),
    compare("iii_unaffected_base_equals_head_partitioned", step_rows("base", "partitioned", unaffected),
        step_rows("head", "partitioned", unaffected), c("kept_digest", "table_digest", "n_kept")),
    compare("iv_head_earlier_equals_base", pick("s4", "head", "earlier_option", affected),
        pick("s4", "base", "shipped", affected), c("kept_digest", "table_digest", "n_kept")))
utils::write.table(result, file.path(out, "compare.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
summary <- stats::aggregate(equal ~ criterion, result, function(x) paste0(sum(x), "/", length(x)))
print(summary)
cat("side status:", paste(names(status), status, collapse = ", "), "\n")
quit(status = if (any(status != 0L) || !all(result$equal)) 1L else 0L)
