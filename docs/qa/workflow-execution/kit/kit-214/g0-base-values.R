# KIT-214 G0: on the base commit, each site's value on every fixture of the ticket's section 4 (acceptance criterion A1).
#
# Usage, from the package root on a run-tests runner (the tests-only commit: nothing under R/ differs from the base):
#   Rscript docs/qa/workflow-execution/kit/kit-214/g0-base-values.R [<out dir>]
#   <out dir>  default docs/qa/workflow-execution/kit/kit-214/g0
#
# Under pkgload::load_all(".", export_all = TRUE, helpers = TRUE). Sites (ticket 2.1):
#   columns  per assay, every column with its class (whether a fixture keeps its vendor columns);
#   C1       the declared duplicate step, `whole` and `partitioned`, each on its own owner in its own store root, on the
#            object with a text column `.row` (the row's input position) added to each assay after the import and design
#            steps: per duplicated id the kept row's `.row`; the number of duplicated ids; the kept pairs' digest;
#   C2       the QC summary's "Number of Samples:" (omicQcDataSummaryUi);
#   C3       per assay `Samples` and `Missingness`, from the whole object (omicQcAssayStatsFrame) and from the store
#            (readWorkflowOmicsShell with its streamed missingness);
#   C4       analyzeWorkflowOmicsItsd on an owner holding the object: per standard mean and CV, the CV > 30 % count, the
#            long table's `Sample` values, or the error;
#   C5       the normalisation's ITSD selection table through the configuration's `norm$itsd_table` (three arguments, as
#            the panel calls it): per row the feature id, mean and CV.
# Fixtures (ticket 4, 10): the metabolite MS-DIAL seed (`Alignment ID` 1 to 8) and the seed with `Alignment ID` 0 to 7,
#   module_ci_lipid_msdial_assay() with `Peak ID` = 1:3, module_ci_lipid_msdial_assay(duplicate_feature = TRUE) through
#   the lipid MS-DIAL import and built in code, the module-CI QC objects, 2.2's two synthetic tables (lipid and
#   metabolite), the four GUI-entry lanes, S1's (A2) inputs (lipid_canonical LCMS_Pos and metab_lc lcms_pos with one id
#   duplicated by two copied rows, intensities x 0.5 and x 2), the seed with a duplicated `Alignment ID`, E9's table (a
#   numeric annotation column), G7's objects (a text design sample) and one deliberate exact tie per omic.
# Written: <out>/g0-base-values.tsv (fixture, site, variant, assay, key, value). Exit status non-zero when a fixture
# cannot be built; an error at a site is recorded as its value ("error: <message>").
# Synthetic and repository fixtures only.

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1L) args[[1L]] else "docs/qa/workflow-execution/kit/kit-214/g0"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
cat("KIT-214 G0 at", commit, "out", out, "\n")

suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = TRUE))
g0 <- new.env(parent = asNamespace("MultiScholaR"))

local(envir = g0, {
    e2e <- normalizePath(file.path("tests", "testdata", "e2e"), mustWork = TRUE)
    manifest <- jsonlite::fromJSON(file.path(e2e, "manifest.json"), simplifyVector = FALSE)
    lane_of <- function(id) {
        for (lane in manifest$lanes) if (identical(lane$lane_id, id)) return(lane)
        stop("lane not in the manifest: ", id)
    }
    quiet <- function(expr) suppressWarnings(suppressMessages(expr))
    sha <- function(x) digest::digest(x, algo = "sha256")
    num <- function(x) format(x, digits = 15)

    # The application's import path, copied from KIT-202's tests (branch kit-202-intensity-filter-samples at f9d7d910,
    # `k202RunInput`), itself from KIT-198's G0 (`run_input`, b0e913b8).
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
        source_dir <- file.path(tempfile("kit214-g0-"), "source")
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
        path <- file.path(tempfile("kit214-g0-input-"), file_name)
        dir.create(dirname(path), recursive = TRUE)
        if (csv) utils::write.csv(table, path, row.names = FALSE) else {
            utils::write.table(table, path, sep = "\t", quote = FALSE, row.names = FALSE)
        }
        path
    }
    custom_inputs <- list(id = "Feature.Name", annotation = "", samples = "^(WT|KO)_", is_pattern = "^ITSD|_d[0-9]+$")
    # A lane; `seed_file` replaces the lane's first assay file (S1's (A2) inputs).
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
    # S1's (A2): the lane's first assay with its first row copied twice, intensities x 0.5 and x 2.
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
    tie_table <- function() data.frame(feature_id = c("T", "T", "U"), name = c("first", "second", "single"),
        S1 = c(10, 40, 5), S2 = c(20, 30, 5), S3 = c(30, 20, 5), S4 = c(40, 10, 5), stringsAsFactors = FALSE)
    is_lipid <- function(object) methods::is(object, "LipidomicsAssayData")
    binding_of <- function(object) if (is_lipid(object)) lipidomicsDeclarationBinding else metabolomicsDeclarationBinding
    omic_of <- function(object) if (is_lipid(object)) "lipidomics" else "metabolomics"
    assays_of <- function(object) methods::slot(object, binding_of(object)$data_slot)
    set_assays <- function(object, assays) {
        methods::slot(object, binding_of(object)$data_slot, check = FALSE) <- assays
        object
    }
    text_sample <- function(object) {
        assays <- assays_of(object)
        sample <- as.character(object@design_matrix[[object@sample_id]])[[1L]]
        assays[[1L]][[sample]] <- as.character(assays[[1L]][[sample]])
        set_assays(object, assays)
    }
    with_row <- function(object) set_assays(object, lapply(assays_of(object), function(table) {
        table$.row <- as.character(seq_len(nrow(table)))
        table
    }))

    builders <- list(
        metab_msdial = function() run_input("metabolomics", "msdial", file.path(e2e, "metab_lc", "seed_lcms_pos.csv"),
            "LCMS_Pos"),
        metab_msdial_id0 = function() {
            table <- seed()
            table[["Alignment ID"]] <- table[["Alignment ID"]] - 1L
            msdial("metabolomics", table, "seed_lcms_pos_id0.csv")
        },
        metab_msdial_dup = function() {
            table <- seed()
            copy <- table[1L, , drop = FALSE]
            copy[seed_samples] <- copy[seed_samples] * 0.5
            msdial("metabolomics", rbind(table, copy), "seed_lcms_pos_dup.csv")
        },
        lipid_msdial = function() {
            assay <- module_ci_lipid_msdial_assay()
            assay[["Peak ID"]] <- seq_len(nrow(assay))
            msdial("lipidomics", assay, "msdial_lipid_pos.csv")
        },
        lipid_msdial_dup_import = function() msdial("lipidomics", module_ci_lipid_msdial_assay(duplicate_feature = TRUE),
            "msdial_lipid_pos_dup.csv"),
        lipid_msdial_dup_code = function() build("lipidomics", module_ci_lipid_msdial_assay(duplicate_feature = TRUE),
            "Peak ID", "Name", design4(c("WT_1", "WT_2", "KO_1", "KO_2"))),
        module_ci_lipid = function() module_ci_lipid_qc_object(),
        module_ci_metab = function() module_ci_metab_qc_object(),
        e1_lipid = function() build("lipidomics", e1_table(), "feature_id", "name"),
        e1_metab = function() build("metabolomics", e1_table(), "feature_id", "name"),
        e2_lipid = function() build("lipidomics", e2_table(), "feature_id", "mz"),
        e2_metab = function() build("metabolomics", e2_table(), "feature_id", "mz"),
        tie_lipid = function() build("lipidomics", tie_table(), "feature_id", "name"),
        tie_metab = function() build("metabolomics", tie_table(), "feature_id", "name"),
        e9_metab = function() build("metabolomics", seed()[, c("Name", "Average Mz", seed_samples)], "Name", "Average Mz",
            design4(seed_samples)),
        text_sample_lipid = function() text_sample(module_ci_lipid_qc_object()),
        text_sample_metab = function() text_sample(module_ci_metab_qc_object()),
        lipid_canonical = function() lane_input("lipid_canonical"),
        metab_lc = function() lane_input("metab_lc"),
        metab_gc = function() lane_input("metab_gc"),
        metab_combined = function() lane_input("metab_combined"),
        a2_lipid_canonical = function() lane_duplicated("lipid_canonical"),
        a2_metab_lc = function() lane_duplicated("metab_lc"))
    patterns <- c(metab_msdial = "ITSD", metab_msdial_id0 = "ITSD", metab_msdial_dup = "ITSD", e9_metab = "ITSD",
        module_ci_lipid = "^IS_", module_ci_metab = "^IS_", text_sample_lipid = "^IS_", text_sample_metab = "^IS_")

    rows <- list()
    add <- function(fixture, site, variant, assay, key, value) {
        rows[[length(rows) + 1L]] <<- data.frame(fixture = fixture, site = site, variant = variant, assay = assay,
            key = key, value = as.character(value), stringsAsFactors = FALSE)
    }
    site <- function(fixture, name, variant, expr) {
        tryCatch(expr, error = function(e) {
            add(fixture, name, variant, "", "error", paste("error:", gsub("[\t\n]+", " ", conditionMessage(e))))
            cat("KIT-214 G0", fixture, name, variant, "error:", conditionMessage(e), "\n")
        })
    }
    owner_of <- function(object) {
        root <- tempfile("kit214-g0-owner-")
        dir.create(root)
        if (is_lipid(object)) lipidDeclarationsOwner(root, object) else metabDeclarationsOwner(root, object)
    }

    record_c1 <- function(fixture, object, variant) {
        marked <- with_row(object)
        binding <- binding_of(object)
        id <- methods::slot(object, binding$id_slot)
        owner <- owner_of(marked)
        on.exit(closeDIAStateAdapter(owner$owner), add = TRUE)
        step <- omicsDuplicateResolutionStep(binding)
        quiet(runWorkflowOmicsReferenceOperation(owner$owner, workflowOperationEntry(step$operation)$declaration,
            c(step$parameters, list(variant = variant)), step$state_name, owner$workflow$config_list, step$description))
        resolved <- assays_of(quiet(owner$owner$getState(step$state_name)))
        for (assay in names(resolved)) {
            input <- assays_of(marked)[[assay]]
            keys <- as.character(input[[id]])
            duplicated_ids <- unique(keys[duplicated(keys)])
            table <- resolved[[assay]]
            add(fixture, "C1", variant, assay, "n_duplicated_ids", length(duplicated_ids))
            add(fixture, "C1", variant, assay, "n_rows", nrow(table))
            add(fixture, "C1", variant, assay, "kept_pairs_digest",
                sha(sort(paste(as.character(table[[id]]), table$.row, sep = "|"))))
            for (key in duplicated_ids) {
                kept <- table$.row[as.character(table[[id]]) == key]
                add(fixture, "C1", variant, assay, paste0("kept_row:", key), paste(kept, collapse = ";"))
            }
        }
    }

    record_displays <- function(fixture, object) {
        binding <- binding_of(object)
        config <- omicGuiConfig(omic_of(object))
        id <- methods::slot(object, binding$id_slot)
        for (assay in names(assays_of(object))) {
            table <- assays_of(object)[[assay]]
            add(fixture, "columns", "", assay, "id_column", id)
            add(fixture, "columns", "", assay, "annotation_column", object@annotation_id_column)
            add(fixture, "columns", "", assay, "columns",
                paste(sprintf("%s:%s", names(table), vapply(table, function(x) class(x)[[1L]], "")), collapse = ";"))
            add(fixture, "columns", "", assay, "design_samples",
                paste(as.character(object@design_matrix[[object@sample_id]]), collapse = ";"))
        }
        site(fixture, "C2", "", {
            html <- paste(as.character(omicQcDataSummaryUi(object, config)), collapse = "")
            found <- regmatches(html, regexec("Number of Samples:</strong>\\s*</td>\\s*<td>([^<]*)</td>", html))[[1L]]
            add(fixture, "C2", "", "", "number_of_samples", if (length(found) >= 2L) found[[2L]] else "not found")
        })
        site(fixture, "C3", "frame", {
            stats <- omicQcAssayStatsFrame(object, config)
            for (i in seq_len(nrow(stats))) {
                add(fixture, "C3", "frame", stats$Assay[[i]], "Samples", stats$Samples[[i]])
                add(fixture, "C3", "frame", stats$Assay[[i]], "Missingness", num(stats$Missingness[[i]]))
            }
        })
        for (assay in names(assays_of(object))) {
            table <- assays_of(object)[[assay]]
            annotation <- object@annotation_id_column
            if (!(length(annotation) == 1L && !is.na(annotation) && annotation %in% names(table))) annotation <- NULL
            site(fixture, "C5", assay, {
                built <- quiet(config$norm$itsd_table(table, id, annotation))
                for (i in seq_len(nrow(built))) {
                    add(fixture, "C5", "", assay, paste0("row:", built$feature_id[[i]]),
                        paste(num(built$mean_intensity[[i]]), num(built$cv_percent[[i]]), sep = "|"))
                }
            })
        }
        owner <- tryCatch(owner_of(object), error = function(e) {
            add(fixture, "owner", "", "", "error", paste("error:", gsub("[\t\n]+", " ", conditionMessage(e))))
            NULL
        })
        if (is.null(owner)) return(invisible(NULL))
        on.exit(closeDIAStateAdapter(owner$owner), add = TRUE)
        site(fixture, "C3", "shell", {
            shell <- readWorkflowOmicsShell(binding, owner$owner, measures = c("features", "missing"))
            stats <- omicQcAssayStatsFrame(shell, config)
            for (i in seq_len(nrow(stats))) {
                add(fixture, "C3", "shell", stats$Assay[[i]], "Samples", stats$Samples[[i]])
                add(fixture, "C3", "shell", stats$Assay[[i]], "Missingness", num(stats$Missingness[[i]]))
            }
        })
        pattern <- if (fixture %in% names(patterns)) patterns[[fixture]] else if (startsWith(fixture, "metab_") ||
            startsWith(fixture, "lipid_canonical") || startsWith(fixture, "a2_")) "^ITSD|_d[0-9]+$" else "ITSD"
        site(fixture, "C4", pattern, {
            analysis <- quiet(analyzeWorkflowOmicsItsd(binding, owner$owner, pattern, config$qc$itsd))
            metrics <- as.data.frame(analysis$metrics)
            for (i in seq_len(nrow(metrics))) {
                add(fixture, "C4", pattern, metrics$assay[[i]], paste0("mean:", metrics$is_id[[i]]),
                    num(metrics$mean_intensity[[i]]))
                add(fixture, "C4", pattern, metrics$assay[[i]], paste0("cv:", metrics$is_id[[i]]), num(metrics$cv[[i]]))
            }
            add(fixture, "C4", pattern, "", "cv_over_30", sum(metrics$cv > 30, na.rm = TRUE))
            add(fixture, "C4", pattern, "", "result_text_cv_line",
                grep("IS with CV > 30%", strsplit(analysis$result_text, "\n", fixed = TRUE)[[1L]], value = TRUE))
            add(fixture, "C4", pattern, "", "long_samples",
                if (is.null(analysis$long_data)) "none" else paste(sort(unique(analysis$long_data$Sample)), collapse = ";"))
        })
    }

    failed <- FALSE
    for (fixture in names(builders)) {
        cat("KIT-214 G0 fixture", fixture, "\n")
        object <- tryCatch(quiet(builders[[fixture]]()), error = function(e) {
            cat("KIT-214 G0 fixture", fixture, "BUILD ERROR", conditionMessage(e), "\n")
            NULL
        })
        if (is.null(object)) {
            failed <- TRUE
            add(fixture, "build", "", "", "error", "build failed")
            next
        }
        record_displays(fixture, object)
        for (variant in c("whole", "partitioned")) site(fixture, "C1", variant, record_c1(fixture, object, variant))
    }

    table <- do.call(rbind, rows)
    utils::write.table(table, file.path(out, "g0-base-values.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    cat("KIT-214 G0: rows", nrow(table), "fixtures", length(builders), "failed", failed, "\n")
    status <- if (failed) 1L else 0L
})
quit(status = g0$status)
