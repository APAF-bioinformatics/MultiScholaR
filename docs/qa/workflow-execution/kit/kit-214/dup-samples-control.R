# KIT-214 script S1: lipidomics and metabolomics duplicate resolution (C1), base against head, on the declared step
# (ticket section 10 "S1", acceptance criterion A3). Committed with the tests; runs in R2 on the fix commit.
#
# Usage, from the package root of the head commit on a run-tests runner (full git history):
#   Rscript docs/qa/workflow-execution/kit/kit-214/dup-samples-control.R <base commit> [<out dir>]
#   <out dir>  default docs/qa/workflow-execution/kit/kit-214/s1 (made absolute)
#
# KIT-202 S1's mechanism: `git worktree add <tmp> <base>`, then per side a fresh Rscript process (this script, `--side`)
# whose working directory is that tree's root, with pkgload::load_all(<tree>, export_all = TRUE, helpers = TRUE). Each
# side writes <out>/<side>-rows.tsv and <out>/<side>-kept.tsv; the comparison is made here.
# Inputs (the import path is G0's, copied below). Every object's assay tables get a text column `.row` (the row's input
# position as text), added after the import and design steps (the design save drops other columns) and before the
# owner is built; it is neither numeric nor a design sample, so neither rule averages it.
#   (A)  unaffected: the four GUI-entry lanes' design objects and the module-CI lipid and metabolite QC objects (G1);
#   (A2) unaffected and real-shaped: lipid_canonical's LCMS_Pos table and metab_lc's lcms_pos table, each with its first
#        row copied twice (intensities x 0.5 and x 2), through the import path;
#   (B)  affected: 2.2's two tables (numeric id duplicated; numeric annotation with a duplicated text id), lipid and
#        metabolite, built in code; module_ci_lipid_msdial_assay(duplicate_feature = TRUE) built in code (its numeric
#        `RT (min)`, `Precursor m/z`, `Total Score` kept); the metabolite seed with `Alignment ID` 1 duplicated (a copy of
#        its row, intensities x 0.5) through the MS-DIAL import;
#   one input classed on the head at run time: module_ci_lipid_msdial_assay(duplicate_feature = TRUE) through the lipid
#        MS-DIAL import, affected when an assay's numeric columns differ from its design samples less the id (G0 records
#        whether it keeps its vendor columns), unaffected otherwise.
# Routes: the declared duplicate step, `whole` and `partitioned`, one owner (store root) per input, variant, rule and
# side. On the head the (B) and run-time-classed inputs also run under
# `MultiScholaR.duplicate_resolution_samples = "numeric_columns"` (rule `earlier_option`).
# Oracle (head side; base R only, no dplyr and no package function): for each id with more than one row, m_r =
# mean(x_r[!is.na(x_r)]) over the rule's columns in the table's order (-Inf when none); the first row in input order
# among those with the largest m_r is kept. Rule `oracle_design`: D = setdiff(intersect(columns, design samples), id);
# rule `oracle_numeric_columns`: every numeric column. G0 records the base's choice on an exact tie.
# Criteria (written to <out>/compare.tsv; criterion, input, assay, setting, equal):
#   (i)   (A), (A2): base equal to head (kept-pairs digest, table digest, n kept), each variant; failed as vacuous unless
#         at least one lipid and one metabolite (A)/(A2) input has a duplicated id (row `i_non_vacuous`);
#   (ii)  (B): head equal to the oracle over D, and base equal to the oracle over every numeric column (kept pairs);
#   (iii) head `whole` equal to head `partitioned` on every input (kept pairs, table digest, n kept);
#   (iv)  head under the earlier option equal to base on every (B) input, each variant.
# Written: <out>/dup-identity.tsv (input, omic, class, variant, side, rule, assay, n_duplicated_ids, kept-pairs digest,
# table digest, n kept), <out>/compare.tsv, <out>/changed-rows.tsv (input, assay, id, base .row, head .row: the (B) ids
# whose kept row differs between base and head, `whole`). Exit status non-zero when a side fails or any criterion row
# is not equal.

args <- commandArgs(trailingOnly = TRUE)

# ---- Child: one side --------------------------------------------------------------------------------------------------
childMain <- function(side, out) {
    e2e <- normalizePath(file.path("tests", "testdata", "e2e"), mustWork = TRUE)
    manifest <- jsonlite::fromJSON(file.path(e2e, "manifest.json"), simplifyVector = FALSE)
    lane_of <- function(id) {
        for (lane in manifest$lanes) if (identical(lane$lane_id, id)) return(lane)
        stop("lane not in the manifest: ", id)
    }
    quiet <- function(expr) suppressWarnings(suppressMessages(expr))
    sha <- function(x) digest::digest(x, algo = "sha256")
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
        source_dir <- file.path(tempfile("kit214-s1-"), "source")
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
    write_table <- function(table, file_name, csv = TRUE) {
        path <- file.path(tempfile("kit214-s1-input-"), file_name)
        dir.create(dirname(path), recursive = TRUE)
        if (csv) utils::write.csv(table, path, row.names = FALSE) else {
            utils::write.table(table, path, sep = "\t", quote = FALSE, row.names = FALSE)
        }
        path
    }
    custom_inputs <- list(id = "Feature.Name", annotation = "", samples = "^(WT|KO)_", is_pattern = "^ITSD|_d[0-9]+$")
    lane_input <- function(lane_id, seed_file = NULL) {
        lane <- lane_of(lane_id)
        dir <- file.path(e2e, lane$fixture_dir)
        assays <- unlist(lane$assays)
        design <- utils::read.delim(file.path(dir, "design.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
        second <- !is.null(lane$assay2_file) && length(assays) >= 2L
        run_input(lane$omic_type, lane$import_tool, seed_file %||% file.path(dir, lane$seed_file), assays[[1L]],
            assay2_file = if (second) file.path(dir, lane$assay2_file), assay2_name = if (second) assays[[2L]],
            custom = if (identical(lane$import_tool, "custom")) custom_inputs, sanitize = FALSE, design = design)
    }
    lane_duplicated <- function(lane_id) {
        lane <- lane_of(lane_id)
        source <- file.path(e2e, lane$fixture_dir, lane$seed_file)
        table <- utils::read.delim(source, check.names = FALSE, stringsAsFactors = FALSE)
        samples <- grep("^(WT|KO)_", names(table), value = TRUE)
        half <- table[1L, , drop = FALSE]
        half[samples] <- half[samples] * 0.5
        double <- table[1L, , drop = FALSE]
        double[samples] <- double[samples] * 2
        lane_input(lane_id, write_table(rbind(table, half, double), basename(source), csv = FALSE))
    }
    seed <- function() utils::read.csv(file.path(e2e, "metab_lc", "seed_lcms_pos.csv"), check.names = FALSE,
        stringsAsFactors = FALSE)
    seed_samples <- c("ctrl_1", "ctrl_2", "treat_1", "treat_2")
    msdial <- function(omic, table, file_name) run_input(omic, "msdial", write_table(table, file_name), "LCMS_Pos")
    design4 <- function(samples = c("S1", "S2", "S3", "S4")) {
        groups <- rep(c("A", "B"), each = 2L)
        data.frame(Run = samples, group = groups, tech_rep_group = paste(groups, c(1L, 2L, 1L, 2L), sep = "_"),
            stringsAsFactors = FALSE)
    }
    build <- function(omic, table, id_column, annotation_column, design = design4()) {
        if (identical(omic, "lipidomics")) {
            createLipidomicsAssayData(lipid_data = list(LCMS_Pos = table), design_matrix = design,
                lipid_id_column = id_column, annotation_id_column = annotation_column, sample_id = "Run",
                group_id = "group", technical_replicate_id = "tech_rep_group")
        } else {
            createMetaboliteAssayData(metabolite_data = list(LCMS_Pos = table), design_matrix = design,
                metabolite_id_column = id_column, annotation_id_column = annotation_column, sample_id = "Run",
                group_id = "group", technical_replicate_id = "tech_rep_group")
        }
    }
    e1_table <- function() data.frame(feature_id = c(1000L, 1000L, 2000L), name = c("A", "B", "C"), S1 = c(50, 60, 100),
        S2 = c(NA, 60, 110), S3 = c(NA, 60, 120), S4 = c(NA, 60, 130), stringsAsFactors = FALSE)
    e2_table <- function() data.frame(feature_id = c("M1", "M1", "M2"), mz = c(100, 10, 50), S1 = c(10, 12, 30),
        S2 = c(10, 12, 31), S3 = c(10, 12, 32), S4 = c(10, 12, 33), stringsAsFactors = FALSE)
    builders <- list(
        lipid_canonical = function() lane_input("lipid_canonical"),
        metab_lc = function() lane_input("metab_lc"),
        metab_gc = function() lane_input("metab_gc"),
        metab_combined = function() lane_input("metab_combined"),
        module_ci_lipid = function() module_ci_lipid_qc_object(),
        module_ci_metab = function() module_ci_metab_qc_object(),
        a2_lipid_canonical = function() lane_duplicated("lipid_canonical"),
        a2_metab_lc = function() lane_duplicated("metab_lc"),
        e1_lipid = function() build("lipidomics", e1_table(), "feature_id", "name"),
        e1_metab = function() build("metabolomics", e1_table(), "feature_id", "name"),
        e2_lipid = function() build("lipidomics", e2_table(), "feature_id", "mz"),
        e2_metab = function() build("metabolomics", e2_table(), "feature_id", "mz"),
        lipid_msdial_dup_code = function() build("lipidomics", module_ci_lipid_msdial_assay(duplicate_feature = TRUE),
            "Peak ID", "Name", design4(c("WT_1", "WT_2", "KO_1", "KO_2"))),
        metab_msdial_dup = function() {
            table <- seed()
            copy <- table[1L, , drop = FALSE]
            copy[seed_samples] <- copy[seed_samples] * 0.5
            msdial("metabolomics", rbind(table, copy), "seed_lcms_pos_dup.csv")
        },
        lipid_msdial_dup_import = function() msdial("lipidomics", module_ci_lipid_msdial_assay(duplicate_feature = TRUE),
            "msdial_lipid_pos_dup.csv"))
    classes <- c(lipid_canonical = "A", metab_lc = "A", metab_gc = "A", metab_combined = "A", module_ci_lipid = "A",
        module_ci_metab = "A", a2_lipid_canonical = "A2", a2_metab_lc = "A2", e1_lipid = "B", e1_metab = "B",
        e2_lipid = "B", e2_metab = "B", lipid_msdial_dup_code = "B", metab_msdial_dup = "B",
        lipid_msdial_dup_import = "run_time")
    is_lipid <- function(object) methods::is(object, "LipidomicsAssayData")
    binding_of <- function(object) if (is_lipid(object)) lipidomicsDeclarationBinding else metabolomicsDeclarationBinding
    assays_of <- function(object) methods::slot(object, binding_of(object)$data_slot)
    design_columns <- function(table, object, id) {
        setdiff(intersect(names(table), as.character(object@design_matrix[[object@sample_id]])), id)
    }
    numeric_columns <- function(table) names(table)[vapply(table, is.numeric, logical(1))]
    with_row <- function(object) {
        assays <- lapply(assays_of(object), function(table) {
            table$.row <- as.character(seq_len(nrow(table)))
            table
        })
        methods::slot(object, binding_of(object)$data_slot, check = FALSE) <- assays
        object
    }
    oracle_kept <- function(table, id, columns) {
        values <- do.call(cbind, lapply(columns, function(column) as.numeric(table[[column]])))
        means <- apply(values, 1L, function(x) {
            x <- x[!is.na(x)]
            if (length(x)) mean(x) else -Inf
        })
        keys <- as.character(table[[id]])
        kept <- vapply(unique(keys), function(key) {
            rows <- which(keys == key)
            rows[[which.max(means[rows])]]
        }, integer(1))
        data.frame(id = keys[kept], row = table$.row[kept], stringsAsFactors = FALSE)
    }
    pairs_digest <- function(id, row) sha(sort(paste(id, row, sep = "|")))
    rows <- list()
    kept_rows <- list()
    add <- function(...) rows[[length(rows) + 1L]] <<- data.frame(..., stringsAsFactors = FALSE)
    record <- function(input, object, class, variant, rule, assay, input_table, table) {
        id <- methods::slot(object, binding_of(object)$id_slot)
        keys <- as.character(input_table[[id]])
        duplicated_ids <- unique(keys[duplicated(keys)])
        out_keys <- as.character(table[[id]])
        frame <- as.data.frame(table)
        rownames(frame) <- NULL
        add(input = input, omic = binding_of(object)$prefix, class = class, variant = variant, side = side, rule = rule,
            assay = assay, n_duplicated_ids = length(duplicated_ids), kept_digest = pairs_digest(out_keys, table$.row),
            table_digest = sha(frame), n_kept = nrow(table))
        dup <- out_keys %in% duplicated_ids
        if (any(dup)) {
            kept_rows[[length(kept_rows) + 1L]] <<- data.frame(input = input, assay = assay, variant = variant, side = side,
                rule = rule, id = out_keys[dup], row = table$.row[dup], stringsAsFactors = FALSE)
        }
    }
    step <- function(input, object, class, variant, rule) {
        binding <- binding_of(object)
        root <- tempfile("kit214-s1-owner-")
        dir.create(root)
        fixture <- if (is_lipid(object)) lipidDeclarationsOwner(root, object) else metabDeclarationsOwner(root, object)
        on.exit(closeDIAStateAdapter(fixture$owner), add = TRUE)
        dup_step <- omicsDuplicateResolutionStep(binding)
        run <- function() quiet(runWorkflowOmicsReferenceOperation(fixture$owner,
            workflowOperationEntry(dup_step$operation)$declaration, c(dup_step$parameters, list(variant = variant)),
            dup_step$state_name, fixture$workflow$config_list, dup_step$description))
        if (identical(rule, "earlier_option")) {
            withr::with_options(list(MultiScholaR.duplicate_resolution_samples = "numeric_columns"), run())
        } else run()
        resolved <- assays_of(quiet(fixture$owner$getState(dup_step$state_name)))
        for (assay in names(resolved)) {
            record(input, object, class, variant, rule, assay, assays_of(object)[[assay]], resolved[[assay]])
        }
    }
    failed <- FALSE
    for (input in names(builders)) {
        cat("KIT-214 S1", side, "input", input, "\n")
        ok <- tryCatch({
            object <- with_row(quiet(builders[[input]]()))
            id <- methods::slot(object, binding_of(object)$id_slot)
            class <- classes[[input]]
            if (identical(class, "run_time")) {
                differs <- vapply(assays_of(object), function(table) {
                    !setequal(numeric_columns(table), design_columns(table, object, id))
                }, logical(1))
                class <- if (any(differs)) "B" else "A2"
            }
            rules <- if (identical(side, "head") && identical(class, "B")) c("shipped", "earlier_option") else "shipped"
            for (rule in rules) for (variant in c("whole", "partitioned")) step(input, object, class, variant, rule)
            if (identical(side, "head")) {
                for (assay in names(assays_of(object))) {
                    table <- assays_of(object)[[assay]]
                    for (oracle in c("oracle_design", "oracle_numeric_columns")) {
                        columns <- if (identical(oracle, "oracle_design")) design_columns(table, object, id) else {
                            numeric_columns(table)
                        }
                        kept <- oracle_kept(table, id, columns)
                        keys <- as.character(table[[id]])
                        add(input = input, omic = binding_of(object)$prefix, class = class, variant = "", side = side,
                            rule = oracle, assay = assay, n_duplicated_ids = length(unique(keys[duplicated(keys)])),
                            kept_digest = pairs_digest(kept$id, kept$row), table_digest = "", n_kept = nrow(kept))
                    }
                }
            }
            TRUE
        }, error = function(e) {
            cat("KIT-214 S1", side, "input", input, "ERROR", class(e)[[1L]], conditionMessage(e), "\n")
            FALSE
        })
        if (!ok) failed <- TRUE
    }
    utils::write.table(do.call(rbind, rows), file.path(out, paste0(side, "-rows.tsv")), sep = "\t", row.names = FALSE,
        quote = FALSE)
    utils::write.table(do.call(rbind, kept_rows), file.path(out, paste0(side, "-kept.tsv")), sep = "\t",
        row.names = FALSE, quote = FALSE)
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
out <- if (length(args) >= 2L) args[[2L]] else "docs/qa/workflow-execution/kit/kit-214/s1"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
head_root <- normalizePath(".")
head_commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
script <- normalizePath(file.path("docs", "qa", "workflow-execution", "kit", "kit-214", "dup-samples-control.R"))
cat("KIT-214 S1: head", head_commit, "base", base_commit, "out", out, "\n")
rscript <- file.path(R.home("bin"), "Rscript")
runSide <- function(tree, side) {
    console <- suppressWarnings(system2(rscript, c(shQuote(script), "--side", shQuote(tree), side, shQuote(out)),
        stdout = TRUE, stderr = TRUE))
    writeLines(console, file.path(out, paste0(side, "-console.txt")))
    cat(grep("^KIT-214 S1", console, value = TRUE), sep = "\n")
    cat(utils::tail(console, 10L), sep = "\n")
    attr(console, "status") %||% 0L
}
base_tree <- file.path(tempdir(), "kit214-s1-base")
invisible(system2("git", c("worktree", "add", "--detach", shQuote(base_tree), base_commit), stdout = TRUE, stderr = TRUE))
stopifnot(dir.exists(file.path(base_tree, "tests", "testthat")))
status <- c(base = runSide(base_tree, "base"), head = runSide(head_root, "head"))

read_side <- function(side, kind) {
    path <- file.path(out, paste0(side, "-", kind, ".tsv"))
    if (!file.exists(path) || !file.size(path)) return(NULL)
    utils::read.delim(path, colClasses = "character", na.strings = "NA")
}
rows <- rbind(read_side("base", "rows"), read_side("head", "rows"))
utils::write.table(rows, file.path(out, "dup-identity.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# The classes are the head's (the run-time-classed input is classed there).
head_classes <- unique(rows[rows$side == "head", c("input", "class")])
class_of <- stats::setNames(head_classes$class, head_classes$input)
inputs_in <- function(classes) names(class_of)[class_of %in% classes]
unaffected <- inputs_in(c("A", "A2"))
affected <- inputs_in("B")

key <- function(x) paste(x$input, x$assay, x$variant, sep = "|")
pick <- function(side, rule, inputs, variant = NULL) {
    x <- rows[rows$side == side & rows$rule == rule & rows$input %in% inputs, , drop = FALSE]
    if (!is.null(variant)) {
        x <- x[x$variant == variant, , drop = FALSE]
        x$variant <- ""
    }
    x
}
compare <- function(criterion, a, b, fields) {
    keys <- union(key(a), key(b))
    a <- a[match(keys, key(a)), , drop = FALSE]
    b <- b[match(keys, key(b)), , drop = FALSE]
    equal <- !is.na(a$input) & !is.na(b$input)
    for (field in fields) equal <- equal & (a[[field]] == b[[field]]) %in% TRUE
    parts <- strsplit(keys, "|", fixed = TRUE)
    data.frame(criterion = criterion, input = vapply(parts, `[`, "", 1L), assay = vapply(parts, `[`, "", 2L),
        setting = vapply(parts, function(p) if (length(p) >= 3L) p[[3L]] else "", ""), equal = equal,
        stringsAsFactors = FALSE)
}
step_fields <- c("kept_digest", "table_digest", "n_kept")
oracle_fields <- c("kept_digest", "n_kept")
oracle_of <- function(rule, inputs, variant) {
    x <- pick("head", rule, inputs)
    x$variant <- variant
    x
}
result <- list(
    compare("i_base_equals_head", pick("base", "shipped", unaffected), pick("head", "shipped", unaffected), step_fields))
# Non-vacuity of (i): at least one lipid and one metabolite (A)/(A2) input with a duplicated id.
duplicated_unaffected <- rows[rows$side == "head" & rows$rule == "shipped" & rows$input %in% unaffected &
    as.integer(rows$n_duplicated_ids) >= 1L, , drop = FALSE]
for (omic in c("lipid", "metab")) {
    result[[length(result) + 1L]] <- data.frame(criterion = "i_non_vacuous", input = omic, assay = "", setting = "",
        equal = any(duplicated_unaffected$omic == omic), stringsAsFactors = FALSE)
}
for (variant in c("whole", "partitioned")) {
    result[[length(result) + 1L]] <- compare("ii_head_equals_oracle_design", pick("head", "shipped", affected, variant),
        pick("head", "oracle_design", affected), oracle_fields)
    result[[length(result) + 1L]] <- compare("ii_base_equals_oracle_numeric_columns",
        pick("base", "shipped", affected, variant), pick("head", "oracle_numeric_columns", affected), oracle_fields)
}
all_inputs <- names(class_of)
result[[length(result) + 1L]] <- compare("iii_head_whole_equals_partitioned",
    pick("head", "shipped", all_inputs, "whole"), pick("head", "shipped", all_inputs, "partitioned"), step_fields)
result[[length(result) + 1L]] <- compare("iv_head_earlier_equals_base", pick("head", "earlier_option", affected),
    pick("base", "shipped", affected), step_fields)
result <- do.call(rbind, result)
utils::write.table(result, file.path(out, "compare.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# The (B) ids whose kept row differs between base and head (`whole`, the shipped rule on each side).
kept <- rbind(read_side("base", "kept"), read_side("head", "kept"))
changed <- data.frame(input = character(), assay = character(), id = character(), base_row = character(),
    head_row = character(), stringsAsFactors = FALSE)
if (!is.null(kept)) {
    kept <- kept[kept$variant == "whole" & kept$rule == "shipped" & kept$input %in% affected, , drop = FALSE]
    base_kept <- kept[kept$side == "base", , drop = FALSE]
    head_kept <- kept[kept$side == "head", , drop = FALSE]
    both <- merge(base_kept[c("input", "assay", "id", "row")], head_kept[c("input", "assay", "id", "row")],
        by = c("input", "assay", "id"), all = TRUE, suffixes = c("_base", "_head"))
    differs <- !(both$row_base == both$row_head) %in% TRUE
    changed <- data.frame(input = both$input[differs], assay = both$assay[differs], id = both$id[differs],
        base_row = both$row_base[differs], head_row = both$row_head[differs], stringsAsFactors = FALSE)
}
utils::write.table(changed, file.path(out, "changed-rows.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
summary <- stats::aggregate(equal ~ criterion, result, function(x) paste0(sum(x), "/", length(x)))
print(summary)
cat("changed rows:", nrow(changed), "\n")
cat("side status:", paste(names(status), status, collapse = ", "), "\n")
quit(status = if (any(status != 0L) || !all(result$equal)) 1L else 0L)
