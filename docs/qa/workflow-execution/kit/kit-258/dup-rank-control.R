# KIT-258 script S1: lipidomics and metabolomics duplicate resolution, base against head, on the declared step
# (ticket section 10 "S1" as r1.2 re-specified it (Re-check 1, R1-02); acceptance criterion A4 (a) to (e)).
#
# Usage, from the package root of the head commit on a run-tests runner (full git history):
#   Rscript docs/qa/workflow-execution/kit/kit-258/dup-rank-control.R <base commit> [<out dir>]
#   <out dir>  default docs/qa/workflow-execution/kit/kit-258/s1f (made absolute; s1 held the first run, before the
#              owner's VAL-284, at 81fdcf7c with "measured" = not NA)
#
# KIT-214's `dup-samples-control.R` (docs/qa/workflow-execution/kit/kit-214/), its inputs and variants kept; its classes
# and criteria (i) to (iv) replaced by A4's (a) to (e). Mechanism as there: `git worktree add <tmp> <base>`, then per
# side a fresh Rscript process (this script, `--side`) whose working directory is that tree's root, with
# pkgload::load_all(<tree>, export_all = TRUE, helpers = TRUE); each side writes <out>/<side>-rows.tsv and
# <out>/<side>-kept.tsv (the base also <out>/base-groups.tsv); the comparison is made here.
#
# Inputs (KIT-214's 15, built as there; every assay gets a text column `.row`, the row's input position): the four
# GUI-entry lanes, the module-CI lipid and metabolite QC objects, the two lane tables with a duplicated first row
# (x 0.5, x 2), section 2.2's two synthetic tables per omic, module_ci_lipid_msdial_assay(duplicate_feature = TRUE)
# built in code and through the lipid MS-DIAL import, the metabolite seed with a duplicated row through the import.
# Routes: the declared duplicate step, `whole` and `partitioned`, one owner (store root) per input, variant, rule, side.
# Rules: base `shipped` (KIT-214: design samples, highest mean) and `earlier_samples`
# (MultiScholaR.duplicate_resolution_samples = "numeric_columns"); head `shipped` (KIT-258), `earlier_rank`
# (MultiScholaR.duplicate_resolution_rank = "mean") and `earlier_both` (both earlier options).
# Classes, on the base side, per input and assay, over D = setdiff(intersect(columns, design samples), id), a value
# "measured" when neither NA nor 0 (VAL-284 (1)): `a` no duplicated id; `b` every duplicated id's rows equally
# measured; `e` the rest. Groups: every row of a duplicated id with its measured count, the mean of its measured values,
# the mean of its non-missing values (KIT-214's) and its annotation (<out>/base-groups.tsv).
# Oracles (head side; base R only): `oracle_kit214` the highest mean of the non-missing values over D (zeros averaged),
# first in input order on a tie (base's rule); `oracle_kit258` the most measured values over D (neither NA nor 0), then
# the highest mean of the measured values, then the first (head's rule, VAL-254 and VAL-284).
# Criteria (<out>/compare.tsv; criterion, input, assay, setting, equal):
#   (a)  class `a`: base shipped equal to head shipped (kept-pairs digest, table digest, n kept), each variant;
#   (b)  class `b`: the same; and `b_group`: every equally measured duplicated id of any assay keeps the same row on
#        base and head, each variant; non-vacuous rows `a_non_vacuous`, `b_non_vacuous` (lipid and metabolite);
#   (c)  head shipped equal to `oracle_kit258` (kept pairs), each variant; head `whole` equal to head `partitioned`
#        (kept pairs, table digest, n kept); and base shipped equal to `oracle_kit214` (instrument check), each variant;
#   (d1) head `earlier_rank` equal to base `shipped`, every input, assay and variant (records made under KIT-214);
#   (d2) head `earlier_both` equal to base `earlier_samples`, every input, assay and variant (records made before);
#   (e)  <out>/changed-rows.tsv: every id whose kept row differs between base shipped and head shipped (each variant),
#        with both rows' measured counts, means and annotations; criterion `e_as_predicted`: exactly one id,
#        lipid_msdial_dup_import / LCMS_Pos / P001, base row 2, head row 1, in each variant (ticket section 6;
#        stop rule 1 on any other).
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
    quiet <- function(expr) suppressWarnings(suppressMessages(expr))
    sha <- function(x) digest::digest(x, algo = "sha256")
    # ---- Inputs: KIT-214's, copied unchanged from dup-samples-control.R:54-171 ----
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
        source_dir <- file.path(tempfile("kit258-s1-"), "source")
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
        path <- file.path(tempfile("kit258-s1-input-"), file_name)
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
    # ---- Reading objects ----
    is_lipid <- function(object) methods::is(object, "LipidomicsAssayData")
    binding_of <- function(object) if (is_lipid(object)) lipidomicsDeclarationBinding else metabolomicsDeclarationBinding
    assays_of <- function(object) methods::slot(object, binding_of(object)$data_slot)
    design_columns <- function(table, object, id) {
        setdiff(intersect(names(table), as.character(object@design_matrix[[object@sample_id]])), id)
    }
    with_row <- function(object) {
        assays <- lapply(assays_of(object), function(table) {
            table$.row <- as.character(seq_len(nrow(table)))
            table
        })
        methods::slot(object, binding_of(object)$data_slot, check = FALSE) <- assays
        object
    }
    # Per row over `columns`: the measured count and the mean of the measured values (-Inf when none), a value measured
    # when neither NA nor 0 (VAL-284 (1)); and the mean of the non-missing values, zeros averaged (KIT-214's rule).
    row_stats <- function(table, columns) {
        values <- do.call(cbind, lapply(columns, function(column) as.numeric(table[[column]])))
        mean_of <- function(x) if (length(x)) mean(x) else -Inf
        list(counts = apply(values, 1L, function(x) sum(!is.na(x) & x != 0)),
            means = apply(values, 1L, function(x) mean_of(x[!is.na(x) & x != 0])),
            nonmissing_means = apply(values, 1L, function(x) mean_of(x[!is.na(x)])))
    }
    # The kept row per id: `by_count` FALSE is KIT-214's rule, TRUE KIT-258's.
    oracle_kept <- function(table, id, columns, by_count) {
        stats <- row_stats(table, columns)
        keys <- as.character(table[[id]])
        kept <- vapply(unique(keys), function(key) {
            rows <- which(keys == key)
            if (!by_count) return(rows[[which.max(stats$nonmissing_means[rows])]])
            rows <- rows[stats$counts[rows] == max(stats$counts[rows])]
            rows[[which.max(stats$means[rows])]]
        }, integer(1))
        data.frame(id = keys[kept], row = table$.row[kept], stringsAsFactors = FALSE)
    }
    pairs_digest <- function(id, row) sha(sort(paste(id, row, sep = "|")))
    rows <- list()
    kept_rows <- list()
    groups <- list()
    add <- function(...) rows[[length(rows) + 1L]] <<- data.frame(..., stringsAsFactors = FALSE)
    record <- function(input, object, variant, rule, assay, input_table, table) {
        id <- methods::slot(object, binding_of(object)$id_slot)
        keys <- as.character(input_table[[id]])
        duplicated_ids <- unique(keys[duplicated(keys)])
        out_keys <- as.character(table[[id]])
        frame <- as.data.frame(table)
        rownames(frame) <- NULL
        add(input = input, omic = binding_of(object)$prefix, class = "", variant = variant, side = side, rule = rule,
            assay = assay, n_duplicated_ids = length(duplicated_ids), kept_digest = pairs_digest(out_keys, table$.row),
            table_digest = sha(frame), n_kept = nrow(table))
        dup <- out_keys %in% duplicated_ids
        if (any(dup)) {
            kept_rows[[length(kept_rows) + 1L]] <<- data.frame(input = input, assay = assay, variant = variant,
                side = side, rule = rule, id = out_keys[dup], row = table$.row[dup], stringsAsFactors = FALSE)
        }
    }
    rule_options <- list(
        shipped = list(),
        earlier_samples = list(MultiScholaR.duplicate_resolution_samples = "numeric_columns"),
        earlier_rank = list(MultiScholaR.duplicate_resolution_rank = "mean"),
        earlier_both = list(MultiScholaR.duplicate_resolution_samples = "numeric_columns",
            MultiScholaR.duplicate_resolution_rank = "mean"))
    side_rules <- if (identical(side, "base")) c("shipped", "earlier_samples") else c("shipped", "earlier_rank", "earlier_both")
    step <- function(input, object, variant, rule) {
        binding <- binding_of(object)
        root <- tempfile("kit258-s1-owner-")
        dir.create(root)
        fixture <- if (is_lipid(object)) lipidDeclarationsOwner(root, object) else metabDeclarationsOwner(root, object)
        on.exit(closeDIAStateAdapter(fixture$owner), add = TRUE)
        dup_step <- omicsDuplicateResolutionStep(binding)
        withr::with_options(rule_options[[rule]], quiet(runWorkflowOmicsReferenceOperation(fixture$owner,
            workflowOperationEntry(dup_step$operation)$declaration, c(dup_step$parameters, list(variant = variant)),
            dup_step$state_name, fixture$workflow$config_list, dup_step$description)))
        resolved <- assays_of(quiet(fixture$owner$getState(dup_step$state_name)))
        for (assay in names(resolved)) {
            record(input, object, variant, rule, assay, assays_of(object)[[assay]], resolved[[assay]])
        }
    }
    failed <- FALSE
    for (input in names(builders)) {
        cat("KIT-258 S1", side, "input", input, "\n")
        ok <- tryCatch({
            object <- with_row(quiet(builders[[input]]()))
            id <- methods::slot(object, binding_of(object)$id_slot)
            annotation <- object@annotation_id_column
            for (rule in side_rules) for (variant in c("whole", "partitioned")) step(input, object, variant, rule)
            for (assay in names(assays_of(object))) {
                table <- assays_of(object)[[assay]]
                columns <- design_columns(table, object, id)
                keys <- as.character(table[[id]])
                duplicated_ids <- unique(keys[duplicated(keys)])
                if (identical(side, "base")) {
                    # Classes and groups, on the base side.
                    stats <- row_stats(table, columns)
                    equal <- vapply(duplicated_ids, function(key) {
                        length(unique(stats$counts[keys == key])) == 1L
                    }, logical(1))
                    class <- if (!length(duplicated_ids)) "a" else if (all(equal)) "b" else "e"
                    add(input = input, omic = binding_of(object)$prefix, class = class, variant = "", side = side,
                        rule = "class", assay = assay, n_duplicated_ids = length(duplicated_ids), kept_digest = "",
                        table_digest = "", n_kept = NA_integer_)
                    in_group <- which(keys %in% duplicated_ids)
                    if (length(in_group)) {
                        names_of <- if (length(annotation) == 1L && !is.na(annotation) && annotation %in% names(table)) {
                            as.character(table[[annotation]][in_group])
                        } else rep("", length(in_group))
                        groups[[length(groups) + 1L]] <- data.frame(input = input, assay = assay, id = keys[in_group],
                            row = table$.row[in_group], n_measured = stats$counts[in_group],
                            mean_measured = as.character(round(stats$means[in_group], 6)),
                            mean_nonmissing = as.character(round(stats$nonmissing_means[in_group], 6)),
                            annotation = names_of,
                            equally_measured = equal[match(keys[in_group], duplicated_ids)], stringsAsFactors = FALSE)
                    }
                } else {
                    for (oracle in c("oracle_kit214", "oracle_kit258")) {
                        kept <- oracle_kept(table, id, columns, by_count = identical(oracle, "oracle_kit258"))
                        add(input = input, omic = binding_of(object)$prefix, class = "", variant = "", side = side,
                            rule = oracle, assay = assay, n_duplicated_ids = length(duplicated_ids),
                            kept_digest = pairs_digest(kept$id, kept$row), table_digest = "", n_kept = nrow(kept))
                    }
                }
            }
            TRUE
        }, error = function(e) {
            cat("KIT-258 S1", side, "input", input, "ERROR", class(e)[[1L]], conditionMessage(e), "\n")
            FALSE
        })
        if (!ok) failed <- TRUE
    }
    utils::write.table(do.call(rbind, rows), file.path(out, paste0(side, "-rows.tsv")), sep = "\t", row.names = FALSE,
        quote = FALSE)
    utils::write.table(do.call(rbind, kept_rows), file.path(out, paste0(side, "-kept.tsv")), sep = "\t",
        row.names = FALSE, quote = FALSE)
    if (identical(side, "base")) {
        utils::write.table(do.call(rbind, groups), file.path(out, "base-groups.tsv"), sep = "\t", row.names = FALSE,
            quote = FALSE)
    }
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
out <- if (length(args) >= 2L) args[[2L]] else "docs/qa/workflow-execution/kit/kit-258/s1f"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
head_root <- normalizePath(".")
head_commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
script <- normalizePath(file.path("docs", "qa", "workflow-execution", "kit", "kit-258", "dup-rank-control.R"))
cat("KIT-258 S1: head", head_commit, "base", base_commit, "out", out, "\n")
rscript <- file.path(R.home("bin"), "Rscript")
runSide <- function(tree, side) {
    console <- suppressWarnings(system2(rscript, c(shQuote(script), "--side", shQuote(tree), side, shQuote(out)),
        stdout = TRUE, stderr = TRUE))
    writeLines(console, file.path(out, paste0(side, "-console.txt")))
    cat(grep("^KIT-258 S1", console, value = TRUE), sep = "\n")
    cat(utils::tail(console, 10L), sep = "\n")
    attr(console, "status") %||% 0L
}
base_tree <- file.path(tempdir(), "kit258-s1-base")
invisible(system2("git", c("worktree", "add", "--detach", shQuote(base_tree), base_commit), stdout = TRUE, stderr = TRUE))
stopifnot(dir.exists(file.path(base_tree, "tests", "testthat")))
status <- c(base = runSide(base_tree, "base"), head = runSide(head_root, "head"))

read_side <- function(name) {
    path <- file.path(out, paste0(name, ".tsv"))
    if (!file.exists(path) || !file.size(path)) return(NULL)
    utils::read.delim(path, colClasses = "character", na.strings = "NA")
}
rows <- rbind(read_side("base-rows"), read_side("head-rows"))
utils::write.table(rows, file.path(out, "dup-identity.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
groups <- read_side("base-groups")

# Classes per input and assay, from the base side.
classes <- rows[rows$side == "base" & rows$rule == "class", c("input", "omic", "assay", "class")]
class_key <- function(x) paste(x$input, x$assay, sep = "|")
in_class <- function(x, class) x[class_key(x) %in% class_key(classes[classes$class == class, , drop = FALSE]), , drop = FALSE]

key <- function(x) paste(x$input, x$assay, x$variant, sep = "|")
pick <- function(side, rule, variant = NULL) {
    x <- rows[rows$side == side & rows$rule == rule, , drop = FALSE]
    if (!is.null(variant)) {
        x <- x[x$variant == variant, , drop = FALSE]
        x$variant <- variant
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
oracle_as <- function(rule, variant) {
    x <- rows[rows$side == "head" & rows$rule == rule, , drop = FALSE]
    x$variant <- variant
    x
}
flag <- function(criterion, input, equal) data.frame(criterion = criterion, input = input, assay = "", setting = "",
    equal = equal, stringsAsFactors = FALSE)
result <- list()
# (a), (b): assays without a duplicated id, and assays whose duplicated ids are all equally measured.
for (class in c("a", "b")) {
    result[[length(result) + 1L]] <- compare(paste0(class, "_base_equals_head"), in_class(pick("base", "shipped"), class),
        in_class(pick("head", "shipped"), class), step_fields)
    for (omic in c("lipid", "metab")) {
        result[[length(result) + 1L]] <- flag(paste0(class, "_non_vacuous"), omic,
            any(classes$class == class & classes$omic == omic))
    }
}
# (b) per group: every equally measured duplicated id keeps the same row on base and head, each variant.
kept <- rbind(read_side("base-kept"), read_side("head-kept"))
kept <- kept[kept$rule == "shipped", , drop = FALSE]
both <- merge(kept[kept$side == "base", c("input", "assay", "variant", "id", "row")],
    kept[kept$side == "head", c("input", "assay", "variant", "id", "row")],
    by = c("input", "assay", "variant", "id"), all = TRUE, suffixes = c("_base", "_head"))
equal_groups <- unique(groups[groups$equally_measured == "TRUE", c("input", "assay", "id")])
group_key <- function(x) paste(x$input, x$assay, x$id, sep = "|")
b_groups <- both[group_key(both) %in% group_key(equal_groups), , drop = FALSE]
result[[length(result) + 1L]] <- data.frame(criterion = "b_group_base_equals_head", input = b_groups$input,
    assay = b_groups$assay, setting = paste(b_groups$variant, b_groups$id, sep = "|"),
    equal = (b_groups$row_base == b_groups$row_head) %in% TRUE, stringsAsFactors = FALSE)
result[[length(result) + 1L]] <- flag("b_group_non_vacuous", "", nrow(b_groups) >= 1L &&
    nrow(b_groups) == 2L * nrow(equal_groups))
# (c): head = KIT-258's oracle, whole = partitioned; base = KIT-214's oracle (instrument check).
for (variant in c("whole", "partitioned")) {
    result[[length(result) + 1L]] <- compare("c_head_equals_oracle_kit258", pick("head", "shipped", variant),
        oracle_as("oracle_kit258", variant), oracle_fields)
    result[[length(result) + 1L]] <- compare("c_base_equals_oracle_kit214", pick("base", "shipped", variant),
        oracle_as("oracle_kit214", variant), oracle_fields)
}
whole <- pick("head", "shipped", "whole")
whole$variant <- ""
partitioned <- pick("head", "shipped", "partitioned")
partitioned$variant <- ""
result[[length(result) + 1L]] <- compare("c_head_whole_equals_partitioned", whole, partitioned, step_fields)
# (d1), (d2): the earlier rules on the head reproduce the base.
result[[length(result) + 1L]] <- compare("d1_head_earlier_rank_equals_base", pick("head", "earlier_rank"),
    pick("base", "shipped"), step_fields)
result[[length(result) + 1L]] <- compare("d2_head_earlier_both_equals_base_earlier_samples", pick("head", "earlier_both"),
    pick("base", "earlier_samples"), step_fields)
# (e): the ids whose kept row changed, with both rows' measured counts, means and annotations.
changed <- both[!(both$row_base == both$row_head) %in% TRUE, , drop = FALSE]
describe <- function(input, assay, id, row) {
    g <- groups[groups$input == input & groups$assay == assay & groups$id == id & groups$row == row, , drop = FALSE]
    if (nrow(g) != 1L) return(c(n_measured = NA, mean_measured = NA, mean_nonmissing = NA, annotation = NA))
    c(n_measured = g$n_measured, mean_measured = g$mean_measured, mean_nonmissing = g$mean_nonmissing,
        annotation = g$annotation)
}
base_info <- if (nrow(changed)) t(mapply(describe, changed$input, changed$assay, changed$id, changed$row_base)) else NULL
head_info <- if (nrow(changed)) t(mapply(describe, changed$input, changed$assay, changed$id, changed$row_head)) else NULL
changed_out <- data.frame(input = changed$input, assay = changed$assay, variant = changed$variant, id = changed$id,
    base_row = changed$row_base, base_n_measured = if (nrow(changed)) base_info[, "n_measured"] else character(),
    base_mean = if (nrow(changed)) base_info[, "mean_measured"] else character(),
    base_mean_nonmissing = if (nrow(changed)) base_info[, "mean_nonmissing"] else character(),
    base_annotation = if (nrow(changed)) base_info[, "annotation"] else character(),
    head_row = changed$row_head, head_n_measured = if (nrow(changed)) head_info[, "n_measured"] else character(),
    head_mean = if (nrow(changed)) head_info[, "mean_measured"] else character(),
    head_mean_nonmissing = if (nrow(changed)) head_info[, "mean_nonmissing"] else character(),
    head_annotation = if (nrow(changed)) head_info[, "annotation"] else character(), stringsAsFactors = FALSE)
rownames(changed_out) <- NULL
utils::write.table(changed_out, file.path(out, "changed-rows.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
for (variant in c("whole", "partitioned")) {
    x <- changed_out[changed_out$variant == variant, , drop = FALSE]
    as_predicted <- nrow(x) == 1L && identical(x$input, "lipid_msdial_dup_import") && identical(x$assay, "LCMS_Pos") &&
        identical(x$id, "P001") && identical(x$base_row, "2") && identical(x$head_row, "1")
    result[[length(result) + 1L]] <- data.frame(criterion = "e_as_predicted", input = "", assay = "", setting = variant,
        equal = as_predicted, stringsAsFactors = FALSE)
}
result <- do.call(rbind, result)
utils::write.table(result, file.path(out, "compare.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(classes, file.path(out, "classes.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
summary <- stats::aggregate(equal ~ criterion, result, function(x) paste0(sum(x), "/", length(x)))
print(summary)
cat("changed ids (both variants):", nrow(changed_out), "\n")
cat("side status:", paste(names(status), status, collapse = ", "), "\n")
quit(status = if (any(status != 0L) || !all(result$equal)) 1L else 0L)
