# KIT-073: the shared declaration core of the lipidomics and metabolomics routes (KIT-031/032 wrote them as two
# near-identical files). One `newWorkflowOperation` per step, every one in process (`compute = "in_process"`, run
# through `runSavedOperation` on the `whole` variant), built from a per-omic binding: the class, the stored slots,
# the names and strings each omic records, and the small kernel adaptors where the two omics call their kernels
# differently (utils_workflow_lipid_declarations.R, utils_workflow_metab_declarations.R). The kernels are unchanged and
# called under the names the modules call them by, so the recorded `@args` keep their keys. The wrappers are what the
# modules and the record's replay call: the operation on the head (or a named source), its products on the checkpoint,
# the commit under the module's state name with the module's configuration and description. Nothing here reads GUI
# inputs, the global environment or a Shiny session.

#' One omics declaration: the binding's class in, both stored components hydrated, the merge rules of a partitioned
#' unit supplied.
#' @noRd
newOmicsOperation <- function(binding, name, unit, parameters, product_type, kernel, output_class = binding$class,
                              global_statistics = NULL, fixed = NULL) {
    partitioned <- unit %in% c("row", "complete_sample", "feature_group")
    newWorkflowOperation(name, binding$class, fields = binding$fields, unit = unit, kernel = kernel,
        output_class = output_class, parameters = parameters, product_type = product_type,
        global_statistics = global_statistics, merge = if (partitioned) binding$merge else NULL, fixed = fixed)
}

#' A declared numeric parameter as a kernel reads it: the text before any comment, as a number.
#' @noRd
omicsCleanNumericParameter <- function(value) as.numeric(trimws(sub("#.*$", "", value)))

#' Unique feature ids per assay from saved tables, as `omicsFeatureCounts` counts them, reading the id column alone (KIT-079).
#' @noRd
omicsPartitionedFeatureCounts <- function(binding, store, tables, id_column) {
    context <- store$manager$.context
    vapply(tables, function(table) {
        scope <- newWorkflowTableReadScope(context$store, table, metadata_only = TRUE)
        columns <- scope$metadata$columns
        releaseWorkflowTableReadScope(scope)
        if (id_column %in% columns) length(unique(materializeWorkflowTable(context, table, id_column)[[id_column]])) else as.double(table$metadata$shape[[1L]])
    }, numeric(1))
}

#' Declared row and feature-group working bytes per column of a partitioned omics step (the proteomics row filter's figure).
#' @noRd
omicsPartitionedWorkingBytesPerColumn <- 48

#' Run a partitioned variant (KIT-079): the kernel on the object with empty assays gives what the step records and the
#' design it carries; each assay's output table is written from its pages or groups by `write_assay`; the checkpoint is
#' the one the kernel would have produced, with those tables as its list slot's elements.
#' @param on_template `function(template)`: the unchanged kernel on the object with zero-row assays.
#' @param label The checkpoint's owned-value label.
#' @param write_assay `function(assay, component, prototype, template, partition)` returning the saved table's reference.
#'   `prototype` is the kernel's zero-row output for the assay; `partition` the assay's readers and writers.
#' @return `list(final, tables, source_tables, id_column)`.
#' @noRd
omicsRunPartitioned <- function(binding, store, on_template, label, write_assay) {
    manager <- store$manager
    context <- manager$.context
    source <- store$source
    spec <- readWorkflowStateMetadata(manager, source)$data
    worker_cap <- store$worker_cap %||% (12 * 1024^3)
    source_tables <- spec$components[[binding$data_slot]]
    # KIT-095: a binding whose class refuses an object with empty assays reads its template its own way.
    template <- if (is.function(binding$read_template)) binding$read_template(manager, source) else
        readWorkflowExecutionComponents(manager, source, fields = binding$design_field)
    result <- suppressWarnings(suppressMessages(on_template(template)))
    tables <- source_tables
    for (assay in names(source_tables)) {
        component <- source_tables[[assay]]
        if (!isWorkflowTableReference(component)) next
        scope <- newWorkflowTableReadScope(context$store, component, reader = if (component$metadata$shape[[2L]] > 64L) "duckdb" else "arrow")
        partition <- omicsAssayPartition(scope, context, component, label, assay, worker_cap)
        tables[[assay]] <- tryCatch(suppressWarnings(suppressMessages(write_assay(assay, component,
            methods::slot(result, binding$data_slot)[[assay]], template, partition))), finally = releaseWorkflowTableReadScope(scope))
    }
    result <- inheritWorkflowExecutionSource(result, manager, source)
    # The assays are written above; only the design comes from this object.
    attr(result, "dia_execution_fields") <- binding$design_field
    checkpoint <- checkpointWorkflowExecutionResult(manager, result)
    checkpoint$components[[binding$data_slot]] <- tables
    checkpoint$dimensions$rows <- sum(vapply(tables, function(table) as.integer(table$metadata$shape[[1L]]), integer(1)))
    list(final = store$checkpoint_record(checkpoint, label), tables = tables, source_tables = source_tables,
        id_column = methods::slot(template, binding$id_slot))
}

#' One assay's readers and writers for a partitioned variant.
#' @noRd
omicsAssayPartition <- function(scope, context, component, label, assay, worker_cap) {
    columns <- scope$metadata$columns
    prototype <- readSavedWorkflowTablePrototype(scope)
    identity_of <- function(kind) newSavedWorkflowScientificIdentity(component$entity_id, label, digest::digest(list(assay, kind), algo = "sha256"),
        "omics-partitioned-table-v1", context$software)
    working <- omicsPartitionedWorkingBytesPerColumn * length(columns)
    overhead <- 2 * 67108864 + 2 * 8388608
    list(columns = columns, prototype = prototype,
        numeric_columns = function() names(prototype)[vapply(prototype, is.numeric, logical(1))],
        # Every value of the columns, for a pre-pass that needs them all.
        all_values = function(selected) {
            values <- list()
            streamWorkflowTableRows(scope, selected, function(page, offset, row_names) {
                values[[length(values) + 1L]] <<- unlist(page, use.names = FALSE)
                invisible(NULL)
            })
            unlist(values, use.names = FALSE)
        },
        has_duplicates = function(id_column) {
            if (!scope$metadata$rows) return(FALSE)
            cursor <- openWorkflowTableGroups(scope, id_column, id_column)
            on.exit(cursor$close(), add = TRUE)
            any(cursor$group_sizes > 1L)
        },
        write_rows = function(output_prototype, page_kernel) saveWorkflowPartitionedTable(context, output_prototype, identity_of("rows"), function(append) {
            if (!scope$metadata$rows) return(invisible(NULL))
            partitionByRows(scope, columns, kernel = function(page, index, offset) page_kernel(page), consume = function(value, index) append(value),
                worker_cap = worker_cap, working_bytes_per_row = working, overhead_bytes = overhead, policy = "omics-row-pages-v1")
            invisible(NULL)
        }),
        write_groups = function(output_prototype, id_column, batch_kernel) saveWorkflowPartitionedTable(context, output_prototype, identity_of("groups"), function(append) {
            partitionByFeatureGroup(scope, columns, id_column, kernel = function(batch, index) batch_kernel(batch), consume = function(value, index) append(value),
                worker_cap = worker_cap, working_bytes_per_row = working, overhead_bytes = overhead, policy = "omics-feature-groups-v1", group_order = "key")
            invisible(NULL)
        }))
}

#' Unique feature ids per assay (the QC tabs' count), by the id column when the table has it.
#' @noRd
omicsFeatureCounts <- function(binding, object) {
    id_column <- methods::slot(object, binding$id_slot)
    vapply(methods::slot(object, binding$data_slot), function(assay) {
        if (id_column %in% names(assay)) length(unique(assay[[id_column]])) else nrow(assay)
    }, numeric(1))
}

#' Refuse sample columns the lipid and metabolite intensity filter cannot count (KIT-202, D3, VAL-130).
#' @param assay The assay's name, or NULL for a table given to an exported helper.
#' @param invalid The design samples (or given columns) that are absent from the table or not numeric.
#' @noRd
omicsIntensityFilterRefuse <- function(assay, invalid) {
    shown <- paste(utils::head(invalid, 3L), collapse = ", ")
    message <- if (!length(invalid)) {
        sprintf("Assay '%s': no sample of the design is left to count. The intensity filter counts only the design's samples.",
            assay %||% "assay table")
    } else if (is.null(assay)) {
        sprintf(paste("%d sample column(s) given are missing from the table or not numeric: %s. A run with no value at all",
            "is read as not numeric: remove it in the design."), length(invalid), shown)
    } else {
        sprintf(paste("Assay '%s': %d sample(s) of the design are missing from the table or not numeric: %s. The intensity",
            "filter counts only the design's samples. A run with no value at all is read as not numeric: remove it in the",
            "design."), assay, length(invalid), shown)
    }
    rlang::abort(message, class = "omics_intensity_filter_samples_invalid", assay = assay, columns = invalid)
}

#' The columns the lipid and metabolite intensity filter counts as samples in one assay (KIT-202, VAL-118 (1)): the
#' design's samples less the id column, in the table's order; under the earlier rule of a record made before
#' (`workflowReplayRuleMarkers$intensity_filter_samples`), every numeric column, as then. A design sample absent from
#' the table or not numeric, or no sample left, is refused (D3).
#' @param table_columns The assay table's column names, in order.
#' @param numeric_columns The names of its numeric columns.
#' @param design_samples The design's sample ids.
#' @param id_column The assay's identifier column.
#' @param assay The assay's name (for the refusal).
#' @return The columns to count.
#' @noRd
omicsIntensityFilterSampleColumns <- function(table_columns, numeric_columns, design_samples, id_column, assay) {
    marker <- workflowReplayRuleMarkers$intensity_filter_samples
    if (identical(getOption(marker$option), marker$earlier)) return(numeric_columns)
    design_samples <- setdiff(unique(as.character(design_samples)), id_column)
    present <- intersect(design_samples, table_columns)
    invalid <- c(setdiff(design_samples, table_columns), setdiff(present, numeric_columns))
    samples <- intersect(table_columns, design_samples)
    if (length(invalid) || !length(samples)) omicsIntensityFilterRefuse(assay, invalid)
    samples
}

#' The columns an exported intensity filter helper counts (KIT-202, D1 (d), VAL-130): `sample_columns` as given, refused
#' when one is absent or not numeric; without it, every numeric column except the identifier, with a classed warning
#' naming them (none when no column is counted, where the helper keeps its earlier warning and return).
#' @noRd
omicsIntensityFilterHelperColumns <- function(assay_table, id_column, sample_columns = NULL) {
    numeric_columns <- names(assay_table)[vapply(assay_table, is.numeric, logical(1))]
    if (!is.null(sample_columns)) {
        sample_columns <- as.character(sample_columns)
        invalid <- c(setdiff(sample_columns, names(assay_table)),
            setdiff(intersect(sample_columns, names(assay_table)), numeric_columns))
        if (length(invalid)) omicsIntensityFilterRefuse(NULL, invalid)
        return(sample_columns)
    }
    counted <- setdiff(numeric_columns, id_column)
    if (length(counted)) {
        shown <- paste(utils::head(counted, 5L), collapse = ", ")
        if (length(counted) > 5L) shown <- paste0(shown, ", ...")
        rlang::warn(sprintf(paste("sample_columns not given: every numeric column except '%s' is counted as a sample: %s.",
            "Pass the design's samples when the table holds other numeric columns."), id_column, shown),
            class = "omics_intensity_filter_samples_guessed")
    }
    counted
}

#' Whether the lipid and metabolite intensity filter counts present values per group (KIT-217, VAL-173 (2)): FALSE
#' only while a record made before is replayed under its earlier rule (`workflowReplayRuleMarkers$intensity_filter_missing`).
#' @noRd
omicsIntensityFilterCountsInForce <- function() {
    marker <- workflowReplayRuleMarkers$intensity_filter_missing
    !identical(getOption(marker$option), marker$earlier)
}

#' Refuse a setting or a design the per-group count cannot use (KIT-217, 5.1; texts W-12 to W-15 of the ticket's 5.4).
#' @noRd
omicsIntensityFilterCountsRefuse <- function(message, setting = NULL, value = NULL) {
    rlang::abort(message, class = "omics_intensity_filter_counts_invalid", setting = setting, value = value)
}

#' A setting as a whole number, or NA: a config string has its comment and spaces removed first, as today's two.
#' @noRd
omicsIntensityFilterWholeNumber <- function(value) {
    if (is.character(value)) value <- trimws(sub("#.*$", "", value))
    number <- if (length(value) == 1L) suppressWarnings(as.numeric(value)) else NA_real_
    if (is.na(number) || !is.finite(number) || number != round(number)) NA_integer_ else as.integer(number)
}

#' What a value was entered as, for a refusal.
#' @noRd
omicsIntensityFilterShown <- function(value) {
    if (is.null(value) || !length(value)) return("nothing")
    paste(format(value), collapse = ", ")
}

#' The counted group of each sample (KIT-217, SD-11): its label, or NA for a missing or blank label (in no group, as
#' the QC finalise step reads blanks) and for a group named like "pool" (case-insensitive, the RUV ANOVA's rule).
#' @param labels Group labels, one per sample.
#' @noRd
omicsIntensityFilterCountedGroups <- function(labels) {
    labels <- as.character(labels)
    labels[is.na(labels) | !nzchar(trimws(labels)) | grepl("pool", labels, ignore.case = TRUE)] <- NA_character_
    labels
}

#' The per-group count's settings checked against the groups (KIT-217, 5.1, in the order the S4 methods check them):
#' values per group a whole number of at least 1; the design has the group column; at least one group is not a pool
#' group; groups required a whole number from 1 to the number of counted groups.
#' @param labels Group labels of the samples counted (the design's, or `group_labels`).
#' @param min_values_per_group,min_groups The settings as given.
#' @param group_id The design's group column, named in the refusal of a design without it.
#' @param has_column FALSE for a design without the group column.
#' @return `list(groups, min_values_per_group, min_groups)`, groups as `omicsIntensityFilterCountedGroups` gives them.
#' @noRd
omicsIntensityFilterCountSettings <- function(labels, min_values_per_group, min_groups, group_id = "group", has_column = TRUE) {
    m <- omicsIntensityFilterWholeNumber(min_values_per_group)
    if (is.na(m) || m < 1L) omicsIntensityFilterCountsRefuse(sprintf(
        "Values required per group must be a whole number of at least 1. You entered: %s.",
        omicsIntensityFilterShown(min_values_per_group)), "min_values_per_group", min_values_per_group)
    if (!isTRUE(has_column)) omicsIntensityFilterCountsRefuse(sprintf(paste("The design has no \"%s\" column, so missing",
        "values cannot be counted per group. Assign each sample to a group in the design and run the filter again."), group_id),
        "group_id", group_id)
    groups <- omicsIntensityFilterCountedGroups(labels)
    pools <- unique(as.character(labels)[!is.na(labels) & grepl("pool", as.character(labels), ignore.case = TRUE)])
    if (all(is.na(groups)) && length(pools)) omicsIntensityFilterCountsRefuse(sprintf(paste("No group is left to count",
        "missing values in: every group in the design is named like \"pool\" (%s), and those are not counted. Give at least",
        "one group of samples a name without \"pool\" in the design, then run the filter again."),
        paste(pools, collapse = ", ")), "groups", pools)
    counted <- length(unique(groups[!is.na(groups)]))
    g <- omicsIntensityFilterWholeNumber(min_groups)
    if (is.na(g) || g < 1L || g > counted) omicsIntensityFilterCountsRefuse(sprintf(paste("Groups required must be a whole",
        "number from 1 to %d, the number of groups in the design (groups named like \"pool\" are not counted). You entered: %s."),
        counted, omicsIntensityFilterShown(min_groups)), "min_groups", min_groups)
    list(groups = groups, min_values_per_group = m, min_groups = g)
}

#' The features that pass the per-group count (KIT-217, 5.1): at least `min_groups` counted groups each hold at least
#' min(n_i, `min_values_per_group`) present (finite) values; with `groups` NULL (an exported helper called without
#' group labels, SD-14), at least one present value in any sample.
#' @param values A matrix (or table) of the samples' values, one column per sample.
#' @param groups The counted group of each column (NA: in no counted group), or NULL.
#' @return Logical, one per row.
#' @noRd
omicsIntensityFilterGroupsPass <- function(values, groups, min_values_per_group, min_groups) {
    present <- is.finite(as.matrix(values))
    if (is.null(groups)) return(rowSums(present) >= 1L)
    met <- integer(nrow(present))
    for (group in unique(groups[!is.na(groups)])) {
        columns <- which(!is.na(groups) & groups == group)
        met <- met + (rowSums(present[, columns, drop = FALSE]) >= min(length(columns), min_values_per_group))
    }
    met >= min_groups
}

#' The lipid and metabolite intensity filter's decision per feature (KIT-217, 5.1), shared by the two exported helpers,
#' so the whole and partitioned variants decide alike: today's share (values strictly below the threshold over the
#' samples counted; a missing value is not below, `<` unchanged, VAL-112) under the proportion, and, unless a record
#' made before is replayed under its earlier rule, the per-group count of present values.
#' @param values The samples' values of one table or page, one column per sample.
#' @param groups The counted group of each column (NA for pool and unassigned samples), or NULL (no groups given).
#' @param threshold The assay's percentile threshold.
#' @param proportion The maximum share of samples below the threshold.
#' @param min_values_per_group,min_groups The checked settings (ignored with `groups` NULL).
#' @return Logical, one per row.
#' @noRd
omicsIntensityFilterKeep <- function(values, groups, threshold, proportion, min_values_per_group, min_groups) {
    values <- as.matrix(values)
    share <- rowSums(values < threshold, na.rm = TRUE) / ncol(values)
    keep <- share < proportion
    if (!omicsIntensityFilterCountsInForce()) return(keep)
    keep & omicsIntensityFilterGroupsPass(values, groups, min_values_per_group, min_groups)
}

#' An exported helper's groups and settings (KIT-217, 5.1, SD-14): `group_labels` aligned with the counted columns and
#' checked as the S4 methods check the design; without them, a classed warning and the no-value rule (settings ignored).
#' Under a record's earlier rule nothing is read or checked.
#' @param noun The feature noun's plural, for the warning ("lipids", "metabolites").
#' @return `list(groups, min_values_per_group, min_groups)`; `groups` NULL for the no-value rule.
#' @noRd
omicsIntensityFilterHelperSettings <- function(sample_columns, group_labels, min_values_per_group, min_groups, noun) {
    if (!omicsIntensityFilterCountsInForce()) return(list(groups = NULL, min_values_per_group = NULL, min_groups = NULL))
    if (is.null(group_labels)) {
        rlang::warn(sprintf(paste("No groups were given (group_labels), so missing values are not counted per group: only",
            "%s with no value in any sample are removed, and then the low-intensity rule is applied as before. Give the",
            "design's groups to apply the per-group rule."), noun), class = "omics_intensity_filter_groups_guessed")
        return(list(groups = NULL, min_values_per_group = NULL, min_groups = NULL))
    }
    if (length(group_labels) != length(sample_columns)) omicsIntensityFilterCountsRefuse(sprintf(paste("group_labels gives",
        "%d group(s) for %d sample column(s): give one group label per sample column, in the same order."),
        length(group_labels), length(sample_columns)), "group_labels", length(group_labels))
    omicsIntensityFilterCountSettings(group_labels, min_values_per_group %||% 2L, min_groups %||% 1L)
}

#' Features of one table that fail the per-group count (KIT-217: the step's `removed_missing`, counted first), or 0
#' under a record's earlier rule.
#' @noRd
omicsIntensityFilterRemovedMissing <- function(values, labels, min_values_per_group, min_groups) {
    if (!omicsIntensityFilterCountsInForce() || !nrow(values)) return(0L)
    settings <- omicsIntensityFilterCountSettings(labels, min_values_per_group, min_groups)
    as.integer(sum(!omicsIntensityFilterGroupsPass(values, settings$groups, settings$min_values_per_group,
        settings$min_groups)))
}

#' The intensity filter (row, after the per-assay quantile pre-pass; the user's 2026-10-28 decision removed the
#' `ceiling()` from that quantile in both kernels before any oracle was recorded). KIT-202: the threshold and the share
#' count the design's samples only (`omicsIntensityFilterSampleColumns`); the step records the rule marker. KIT-217:
#' the per-group count of present values (`omicsIntensityFilterKeep`); the step records its marker and, per assay, the
#' rows removed for missing values (`removed_missing`, counted first) and for low intensity (`removed_low`).
#' @noRd
omicsIntensityFilterOperation <- function(binding) {
    filter <- function(object, parameters) binding$filter_intensity(object, parameters)
    # KIT-217: the two settings as the step received them (the declared step always passes both; NULL only in the
    # replay of a record made before, where the earlier rule reads neither).
    count_setting <- function(parameters, i, default) parameters[[binding$intensity_parameters[[i]]]] %||% default
    declaration <- newOmicsOperation(binding, paste0(binding$prefix, "_intensity_filter"), "row",
        parameters = binding$intensity_parameters,
        product_type = paste0(binding$prefix, ".intensity_filter_result"),
        fixed = list(intensity_filter_samples = "design", intensity_filter_missing = "group_counts"),
        global_statistics = stats::setNames(
            "the percentile of the design's sample intensities of one assay table, computed once per assay before any row",
            binding$intensity_threshold),
        kernel = function(object, parameters, store) {
            before <- omicsFeatureCounts(binding, object)
            filtered <- filter(object, parameters)
            removed <- if (omicsIntensityFilterCountsInForce()) omicsIntensityFilterRemovedCounts(binding, object, filtered,
                count_setting(parameters, 3L, 2L), count_setting(parameters, 4L, 1L))
            object <- NULL
            value <- list(final = store$checkpoint(filtered, paste0(binding$prefix, "-intensity-filtered")),
                counts = list(before = before, after = omicsFeatureCounts(binding, filtered)))
            # KIT-217: the reasons, under the current rule only (a record made before replays to the value it saved).
            if (!is.null(removed)) value[c("removed_missing", "removed_low")] <- list(removed$missing, removed$low)
            value
        })
    # KIT-079: the percentile once per assay over its every intensity, then the kernel's own filter per page.
    declaration$partitioned <- function(parameters, store) {
        percentile <- omicsCleanNumericParameter(parameters[[binding$intensity_parameters[[1L]]]])
        proportion <- omicsCleanNumericParameter(parameters[[binding$intensity_parameters[[2L]]]])
        min_values_per_group <- count_setting(parameters, 3L, 2L)
        min_groups <- count_setting(parameters, 4L, 1L)
        removed <- new.env(parent = emptyenv())
        run <- omicsRunPartitioned(binding, store, function(template) filter(template, parameters), paste0(binding$prefix, "-intensity-filtered"),
            write_assay = function(assay, component, prototype, template, partition) {
                removed[[assay]] <- c(missing = 0L, low = 0L)
                id_column <- methods::slot(template, binding$id_slot)
                design <- methods::slot(template, "design_matrix")
                design_samples <- as.character(design[[methods::slot(template, "sample_id")]])
                # KIT-202: the design's samples of this assay (every numeric column under a record's earlier rule).
                sample_columns <- if (id_column %in% partition$columns) omicsIntensityFilterSampleColumns(
                    names(partition$prototype), partition$numeric_columns(), design_samples, id_column, assay)
                intensities <- if (length(sample_columns)) partition$all_values(sample_columns)
                if (!length(intensities) || all(is.na(intensities))) return(partition$write_rows(partition$prototype, identity))
                threshold <- stats::quantile(intensities, na.rm = TRUE, probs = c(percentile / 100))[1]
                intensities <- NULL
                # KIT-217: the design's group of each sample column (the S4 method on the template checked the design).
                group_id <- methods::slot(template, "group_id")
                group_labels <- if (group_id %in% names(design)) as.character(design[[group_id]])[match(sample_columns, design_samples)]
                page_filter <- function(page) {
                    kept <- binding$filter_intensity_page(page, threshold, proportion, id_column, sample_columns, group_labels,
                        min_values_per_group, min_groups)
                    missing <- omicsIntensityFilterRemovedMissing(as.data.frame(page)[, sample_columns, drop = FALSE],
                        group_labels, min_values_per_group, min_groups)
                    removed[[assay]] <- removed[[assay]] + c(missing, nrow(page) - nrow(kept) - missing)
                    kept
                }
                partition$write_rows(suppressWarnings(page_filter(partition$prototype)), page_filter)
            })
        per_assay <- function(kind) stats::setNames(vapply(names(run$source_tables), function(assay)
            if (is.null(removed[[assay]])) 0L else as.integer(removed[[assay]][[kind]]), integer(1)), names(run$source_tables))
        value <- list(final = run$final, counts = list(before = omicsPartitionedFeatureCounts(binding, store, run$source_tables, run$id_column),
            after = omicsPartitionedFeatureCounts(binding, store, run$tables, run$id_column)))
        # KIT-217: the reasons, under the current rule only (a record made before replays to the value it saved).
        if (omicsIntensityFilterCountsInForce()) value[c("removed_missing", "removed_low")] <- list(per_assay("missing"),
            per_assay("low"))
        value
    }
    declaration
}

#' The whole variant's rows removed per assay (KIT-217): for missing values, the rows failing the per-group count (as
#' the filter counts them: an assay the filter returns whole removes none); for low intensity, the other rows removed.
#' @return `list(missing, low)`, named integer vectors by assay.
#' @noRd
omicsIntensityFilterRemovedCounts <- function(binding, object, filtered, min_values_per_group, min_groups) {
    id_column <- methods::slot(object, binding$id_slot)
    design <- methods::slot(object, "design_matrix")
    design_samples <- as.character(design[[methods::slot(object, "sample_id")]])
    group_id <- methods::slot(object, "group_id")
    before <- methods::slot(object, binding$data_slot)
    after <- methods::slot(filtered, binding$data_slot)
    assays <- names(before) %||% as.character(seq_along(before))
    missing <- vapply(seq_along(before), function(i) {
        table <- before[[i]]
        if (!id_column %in% names(table)) return(0L)
        # The filter has checked these columns already (`omicsIntensityFilterSampleColumns`); the count reads the
        # design's samples of the table, the only columns it counts under the current rule.
        sample_columns <- intersect(names(table), setdiff(design_samples, id_column))
        if (!length(sample_columns)) return(0L)
        values <- as.data.frame(table)[, sample_columns, drop = FALSE]
        if (all(is.na(unlist(values, use.names = FALSE)))) return(0L)
        labels <- if (group_id %in% names(design)) as.character(design[[group_id]])[match(sample_columns, design_samples)]
        omicsIntensityFilterRemovedMissing(values, labels, min_values_per_group, min_groups)
    }, integer(1))
    removed <- vapply(seq_along(before), function(i) as.integer(nrow(before[[i]]) - nrow(after[[i]])), integer(1))
    list(missing = stats::setNames(missing, assays), low = stats::setNames(removed - missing, assays))
}

#' Duplicate resolution (feature group by the id column), through the omic's adaptor (`resolve_duplicates` returns
#' the resolved object and the step's statistics under the name the module reads).
#' @noRd
omicsDuplicateResolutionOperation <- function(binding) {
    declaration <- newOmicsOperation(binding, paste0(binding$prefix, "_duplicate_resolution"), "feature_group", parameters = character(),
        product_type = paste0(binding$prefix, ".duplicate_resolution_result"),
        kernel = function(object, parameters, store) {
            resolved <- binding$resolve_duplicates(object)
            c(list(final = store$checkpoint(resolved$object, paste0(binding$prefix, "-duplicates-resolved"))), resolved$value)
        })
    # KIT-079: an assay with no duplicate comes back as it is, in source order. One with any is resolved group by group in
    # the kernel's output order, its id's key order (`group_by`, C locale, missing last: the cursor's `key` order), with every
    # numeric column double: the kernel converts them all before it groups.
    declaration$partitioned <- function(parameters, store) {
        without_samples <- list()
        run <- omicsRunPartitioned(binding, store, function(template) binding$resolve_duplicates(template)$object,
            paste0(binding$prefix, "-duplicates-resolved"),
            write_assay = function(assay, component, prototype, template, partition) {
                id_column <- methods::slot(template, binding$id_slot)
                numeric_columns <- partition$numeric_columns()
                without_samples[[assay]] <<- !length(numeric_columns)
                if (!id_column %in% partition$columns || !length(numeric_columns) || !partition$has_duplicates(id_column)) {
                    return(partition$write_rows(partition$prototype, identity))
                }
                as_double <- function(table) dplyr::mutate(table, dplyr::across(dplyr::any_of(numeric_columns), as.numeric))
                partition$write_groups(tibble::as_tibble(as_double(partition$prototype)), id_column, function(batch) {
                    binding$resolve_duplicates_batch(as_double(batch), id_column, numeric_columns)
                })
            })
        rows <- function(tables) vapply(tables, function(table) as.integer(table$metadata$shape[[1L]]), integer(1))
        c(list(final = run$final), binding$duplicate_value(rows(run$source_tables), rows(run$tables), without_samples))
    }
    declaration
}

#' A store step: the current object under a new name (the QC finalise, the correlation skip, the reset). No kernel,
#' no product; the declaration names the step for the record and the replay.
#' @noRd
omicsStoreOperation <- function(binding, name) {
    newOmicsOperation(binding, name, "store", parameters = character(), product_type = paste0(binding$prefix, ".store_step"),
        kernel = function(object, parameters, store) list())
}

#' Refuse a manual internal-standard selection that names features the state no longer holds (P-19): the selection is
#' the ids of the rows ticked in the table the tab showed, and the state may have moved since (a second Run without
#' Reset removes the standards; a QC refilter removes features). The kernel would normalise by the standards that are
#' left, or skip the assay with a warning. An assay the state does not hold is left to the kernel, as before.
#' @param binding The omic's binding.
#' @param object The state the step runs on.
#' @param feature_ids The selection: feature ids per assay, or NULL for the pattern.
#' @return TRUE, invisibly; else an error of class `omics_itsd_selection_stale` naming the assay and up to three ids.
#' @noRd
omicsAssertManualItsdPresent <- function(binding, object, feature_ids) {
    if (!is.list(feature_ids) || !length(feature_ids)) return(invisible(TRUE))
    assays <- methods::slot(object, binding$data_slot)
    id_column <- methods::slot(object, binding$id_slot)
    for (assay in intersect(names(feature_ids), names(assays))) {
        absent <- setdiff(as.character(feature_ids[[assay]]), as.character(assays[[assay]][[id_column]]))
        if (!length(absent)) next
        rlang::abort(sprintf(paste("%d selected internal standard(s) are no longer in assay '%s': %s%s.",
            "Reset the normalisation, or untick them in the ITSD selection table, then run again."),
            length(absent), assay, paste(utils::head(absent, 3L), collapse = ", "), if (length(absent) > 3L) ", ..." else ""),
            class = "omics_itsd_selection_stale", call = NULL)
    }
    invisible(TRUE)
}

#' Internal-standard normalisation (complete sample after the mean-factor pre-pass: one aggregate standard signal per
#' sample, the mean over all samples, each sample scaled by mean / own aggregate; the standards removed afterwards).
#' @noRd
omicsItsdNormalisationOperation <- function(binding) {
    newOmicsOperation(binding, paste0(binding$prefix, "_itsd_normalisation"), "complete_sample",
        parameters = c("itsd_aggregation", "itsd_feature_ids"), product_type = paste0(binding$prefix, ".itsd_normalisation_result"),
        global_statistics = c(mean_standard_factor =
            "the mean over all samples of each assay's per-sample aggregate standard signal, computed once before any sample is scaled"),
        fixed = list(method = "ITSD"),
        kernel = function(object, parameters, store) {
            omicsAssertManualItsdPresent(binding, object, parameters$itsd_feature_ids)
            normaliseUntransformedDataFn <- normaliseUntransformedData
            normalised <- normaliseUntransformedDataFn(theObject = object, method = "ITSD",
                itsd_aggregation = parameters$itsd_aggregation, itsd_feature_ids = parameters$itsd_feature_ids)
            object <- NULL
            list(final = store$checkpoint(normalised, paste0(binding$prefix, "-itsd-normalised")))
        })
}

#' Per assay, the number of log2 values below 0 in the sample columns (KIT-172; read by the `scale` warning).
#' @noRd
omicsNegativeLogCounts <- function(binding, object) {
    samples <- as.character(methods::slot(object, "design_matrix")[[methods::slot(object, "sample_id")]])
    id_column <- methods::slot(object, binding$id_slot)
    vapply(methods::slot(object, binding$data_slot), function(assay) {
        columns <- setdiff(intersect(names(assay), samples), id_column)
        sum(vapply(columns, function(column) sum(assay[[column]] < 0, na.rm = TRUE), integer(1)))
    }, integer(1))
}

#' Log2 transformation (row: `log2(x)` per value at the default offset 0, `log2(x + offset)` above it).
#' KIT-172: at offset 0 the result value also holds `nonpositive` (per assay, the cells of 0 or below and the cells
#' that are not finite set to missing; summed over pages) and `negative_log` (per assay, the log2 values below 0). The
#' declared operation does not store the method's count field (`args$log_transform_nonpositive`): the partitioned
#' runner checkpoints the template's object, whose assays are empty, so the counts live in the result value only.
#' Above offset 0 the operation and its result are those of earlier releases.
#' @noRd
omicsLog2Operation <- function(binding) {
    transform <- function(object, parameters) {
        logTransformAssaysFn <- logTransformAssays
        logTransformAssaysFn(theObject = object, offset = parameters$offset)
    }
    plain <- function(parameters) isTRUE(parameters$offset == 0)
    # The count warning reaches no tab from a worker or from inside the runner's suppression; the counts are returned.
    quiet <- function(expr) withCallingHandlers(expr, omics_log2_nonpositive = function(w) invokeRestart("muffleWarning"))
    declaration <- newOmicsOperation(binding, paste0(binding$prefix, "_log2"), "row", parameters = "offset",
        product_type = paste0(binding$prefix, ".log2_result"),
        kernel = function(object, parameters, store) {
            transformed <- quiet(transform(object, parameters))
            object <- NULL
            if (!plain(parameters)) return(list(final = store$checkpoint(transformed, paste0(binding$prefix, "-log2"))))
            nonpositive <- transformed@args$log_transform_nonpositive
            transformed@args$log_transform_nonpositive <- NULL
            negative_log <- as.list(omicsNegativeLogCounts(binding, transformed))
            list(final = store$checkpoint(transformed, paste0(binding$prefix, "-log2")), nonpositive = nonpositive,
                negative_log = negative_log)
        })
    # KIT-079: the kernel per page, on the object with that page as its assay.
    declaration$partitioned <- function(parameters, store) {
        nonpositive <- list()
        negative_log <- list()
        run <- omicsRunPartitioned(binding, store, function(template) {
            x <- transform(template, parameters)
            x@args$log_transform_nonpositive <- NULL
            x
        }, paste0(binding$prefix, "-log2"),
            write_assay = function(assay, component, prototype, template, partition) {
                nonpositive[[assay]] <<- c(nonpositive = 0L, nonfinite = 0L)
                negative_log[[assay]] <<- 0L
                partition$write_rows(prototype, function(page) {
                    object <- template
                    methods::slot(object, binding$data_slot) <- stats::setNames(list(page), assay)
                    transformed <- withCallingHandlers(transform(object, parameters), omics_log2_nonpositive = function(w) {
                        nonpositive[[assay]] <<- nonpositive[[assay]] + w$counts[[assay]]
                        invokeRestart("muffleWarning")
                    })
                    if (plain(parameters)) negative_log[[assay]] <<- negative_log[[assay]] + omicsNegativeLogCounts(binding, transformed)[[assay]]
                    methods::slot(transformed, binding$data_slot)[[assay]]
                })
            })
        if (!plain(parameters)) return(list(final = run$final))
        list(final = run$final, nonpositive = nonpositive, negative_log = negative_log)
    }
    declaration
}

#' Between-sample normalisation (whole matrix; `none` is the module's identity, the state still committed).
#' @noRd
omicsBetweenSampleNormalisationOperation <- function(binding) {
    newOmicsOperation(binding, paste0(binding$prefix, "_between_sample_normalisation"), "whole_matrix",
        parameters = "normalisation_method", product_type = paste0(binding$prefix, ".between_sample_normalisation_result"),
        kernel = function(object, parameters, store) {
            normaliseBetweenSamplesFn <- normaliseBetweenSamples
            normalised <- if (identical(parameters$normalisation_method, "none")) object else
                normaliseBetweenSamplesFn(theObject = object, normalisation_method = parameters$normalisation_method)
            object <- NULL
            list(final = store$checkpoint(normalised, paste0(binding$prefix, "-normalised")),
                applied = !identical(parameters$normalisation_method, "none"))
        })
}

omicsRuvParameterNames <- c("ruv_mode", "ruv_grouping_variable", "percentage_min", "percentage_max", "separation_metric",
    "k_penalty_weight", "adaptive_k_penalty", "manual_k", "manual_percentage")

#' RUV-III (whole matrix): the per-assay optimisation (automatic sweep or manual), then `ruvIII_C_Varying` with the
#' best k and the control features per assay (the module's optimisation and correction in one declaration). The
#' optimisation results are a product of the step.
#' @param ruv_design KIT-129 integration item 1. `"legacy"` (the default): today's declaration, unchanged. `"declared"`:
#'   the new design's declaration (`omicsRuvDesignOperation`).
#' @param supplied KIT-129 fifth pass item 3: under `"declared"`, the declaration of the design `supplied_controls`.
#' @noRd
omicsRuvOperation <- function(binding, ruv_design = c("legacy", "declared"), supplied = FALSE) {
    if (identical(match.arg(ruv_design), "declared")) return(omicsRuvDesignOperation(binding, supplied = supplied))
    newOmicsOperation(binding, paste0(binding$prefix, "_ruv"), "whole_matrix", parameters = omicsRuvParameterNames,
        product_type = paste0(binding$prefix, ".ruv_result"),
        kernel = function(object, parameters, store) {
            ruvIII_C_VaryingFn <- ruvIII_C_Varying
            ruv_results <- binding$optimise_ruv(object, parameters$ruv_mode, parameters[setdiff(omicsRuvParameterNames, "ruv_mode")])
            best_k <- binding$extract_best_k(ruv_results)
            controls <- binding$extract_controls(ruv_results)
            corrected <- ruvIII_C_VaryingFn(theObject = object, ruv_grouping_variable = parameters$ruv_grouping_variable,
                ruv_number_k = best_k, ctrl = controls)
            object <- NULL
            final <- store$checkpoint(corrected, paste0(binding$prefix, "-ruv-corrected"))
            corrected <- NULL
            saved_results <- lapply(ruv_results, function(result) {
                for (field in intersect(c("cancor_plot", "best_cancor_plot"), names(result))) {
                    if (inherits(result[[field]], c("gg", "ggplot"))) result[[field]] <- stripProtNormQcPlotOwners(result[[field]])
                }
                result
            })
            list(final = final, optimisation = store$product(saved_results, "optimisation", paste0(binding$prefix, "-ruv-optimisation")),
                best_k = best_k,
                control_count = lapply(controls, function(x) if (is.null(x)) NA_integer_ else if (is.logical(x)) sum(x, na.rm = TRUE) else length(x)),
                percentage = lapply(ruv_results, function(x) x$best_percentage %||% NA_real_))
        })
}

#' The new design's RUV parameters (KIT-129 integration item 1), by the names of KIT-129 section 1e: the design, the
#' replicate-set column, k (a whole number, or one per assay by name; not given under the estimated count), its rule,
#' the count's seed and the supplied controls per assay.
#' @noRd
omicsRuvDesignParameterNames <- c("ruv_design", "ruv_grouping_variable", "ruv_k", "ruv_k_rule", "ruv_k_seed", "ruv_control_ids")

#' The new design's parameters today's declaration does not declare (its refusal reads these).
#' @noRd
omicsRuvDesignOnlyParameterNames <- setdiff(omicsRuvDesignParameterNames, omicsRuvParameterNames)

#' The designs the new declaration takes: the factors kernel's three (all features, supplied controls, first-pass
#' nulls), the technical-replicate design (KIT-129 third pass item 3: all features as controls, the declared units as
#' the replicate sets) and `none` (no RUV, the option set's None).
#' @noRd
OMICS_RUV_DESIGNS <- c("groups_all_features", "technical_replicates", "supplied_controls", "first_pass_nulls", "none")

#' The replicate sets of a declared lipidomics or metabolomics design (KIT-129 third pass item 3), as simulator S4's arms
#' define them: the technical-replicate design (arm 16) the declared technical-replicate column, the biological units,
#' refused when the design declares no technical replicates, when that column does not name each unit alone, or when
#' the request names another replicate-set column; every other design the request's column, or the groups (arms 7, 8,
#' 5 and 23).
#' @param object The normalised object.
#' @param design The declared design.
#' @param grouping The request's `ruv_grouping_variable`, or NULL.
#' @return `list(grouping, structure)`: the replicate-set column and the record's `ruv_replicate_structure` (`group` or
#'   `technical_replicates`, KIT-129 section 1e).
#' @noRd
omicsRuvDesignReplicateSets <- function(object, design, grouping) {
    invalid <- function(message) cli::cli_abort(message, class = "workflow_execution_invalid", call = NULL)
    if (!identical(design, "technical_replicates")) return(list(grouping = grouping %||% object@group_id, structure = "group"))
    column <- object@technical_replicate_id
    frame <- as.data.frame(object@design_matrix)
    technical <- resolveProtTechnicalReplicateBlock(frame, object@group_id, column)
    if (!isTRUE(technical$technical)) {
        invalid(sprintf("the design technical_replicates needs declared technical replicates; this design has none (%s)", technical$reason))
    }
    if (!is.null(grouping) && !identical(grouping, column)) {
        invalid(sprintf("the design technical_replicates takes the technical-replicate column `%s` as its replicate sets, not `%s`",
            column, grouping))
    }
    if (!ruvReplicateSetsAreUnits(frame[[column]], technical$block)) {
        invalid(sprintf(paste("the technical-replicate column `%s` does not name each biological unit alone (the same value",
            "in more than one group); give each unit its own value"), column))
    }
    list(grouping = column, structure = "technical_replicates")
}

#' Today's RUV declaration takes none of the new design's parameters (KIT-129 integration item 1): its identity reads its
#' nine parameters only, so `runSavedOperation` calls this before the identity, and a request that names the new design
#' is refused rather than served today's product (the legacy refusal of every omic, `workflow_execution_invalid`).
#' @param parameters A request's parameters.
#' @return `parameters`, invisibly.
#' @noRd
assertOmicsLegacyRuvParameters <- function(parameters) {
    given <- Filter(\(name) !is.null(parameters[[name, exact = TRUE]]), omicsRuvDesignOnlyParameterNames)
    if (length(given)) {
        cli::cli_abort(c("today's RUV declaration takes none of the new design's parameters; the request carries {.val {given}}",
            i = "the new design runs under omicsRuvOperation(binding, \"declared\")"), class = "workflow_execution_invalid", call = NULL)
    }
    invisible(parameters)
}

#' A value that is not given: NULL, or one NA.
#' @noRd
omicsRuvAbsent <- function(value) is.null(value) || (length(value) == 1L && !is.list(value) && is.na(value))

#' The new design's RUV step for lipidomics and metabolomics (KIT-129 integration item 1): `runOmicsRuvFactors` (KIT-127;
#' k by the shared resolver, KIT-129 item 3) per assay, or no RUV for the design `none`. The state is the corrected object
#' (figures, export and integration read it, decision M10) and carries what the DA tests: the uncorrected assays
#' (`args$ruv_uncorrected_assays`), the per-assay factor records (`args$ruv_factors`; both stored out of line, see
#' `workflowClassDeclarations`) and the declaration `args$ruvFactorsDeclaration` (the design, the k rule and seed, k per
#' assay), which names the DA's method (`workflowSourceRuvDaMethod`). A k is not given under the estimated count: the
#' count gives it, once, here.
#' @param supplied KIT-129 fifth pass item 3: the declaration of the design `supplied_controls`. Its identity names the
#'   list by its hash and count (`RUV_SUPPLIED_LIST_PARAMETERS`) in place of `ruv_control_ids`; the request's one list
#'   (`ruvDeclaredSuppliedList`) is applied to each assay as the listed features that assay holds, and an assay that
#'   holds none of them is refused (the factors kernel needs controls in every assay); listed ids in no assay are named
#'   in a warning (`ruv_control_ids_unmatched`) and in the product. Too few controls in an assay are the kernel's
#'   refusal (k below the controls observed in every sample, decision M32).
#' @noRd
omicsRuvDesignOperation <- function(binding, supplied = FALSE) {
    supplied <- isTRUE(supplied)
    declared_names <- if (supplied) c(setdiff(omicsRuvDesignParameterNames, "ruv_control_ids"), RUV_SUPPLIED_LIST_PARAMETERS) else
        omicsRuvDesignParameterNames
    newOmicsOperation(binding, paste0(binding$prefix, "_ruv"), "whole_matrix", parameters = declared_names,
        product_type = paste0(binding$prefix, ".ruv_result"),
        # KIT-129 fourth pass item 1: the replicate sets by design (`omicsRuvDesignReplicateSets`), not the groups alone.
        fixed = list(replicate_sets = paste("the groups (RUV-III); under the design technical_replicates the declared",
            "technical replicates, the biological units (RUV-III)"), test = "the uncorrected values with W1..Wk in the limma model"),
        kernel = function(object, parameters, store) {
            invalid <- function(message) cli::cli_abort(message, class = "workflow_execution_invalid", call = NULL)
            design <- parameters$ruv_design
            if (!checkmate::test_choice(design, OMICS_RUV_DESIGNS)) {
                invalid(sprintf("ruv_design must be one of %s", paste(OMICS_RUV_DESIGNS, collapse = ", ")))
            }
            slots <- omicsRuvAssaySlots(object)
            assays <- methods::slot(object, slots$data)
            sets <- omicsRuvDesignReplicateSets(object, design, parameters$ruv_grouping_variable)
            grouping <- sets$grouping
            rule <- parameters$ruv_k_rule
            seed <- if (omicsRuvAbsent(parameters$ruv_k_seed)) RUV_FACTOR_COUNT_SEED else parameters$ruv_k_seed
            # KIT-129 fifth pass item 3: the one supplied list, checked against the request's identity of it, then per
            # assay the listed features that assay holds.
            listed <- if (!identical(design, "none")) ruvDeclaredSuppliedList(parameters)
            if (!is.null(listed) != supplied) invalid("the design supplied_controls runs under its own declaration, with its list")
            unmatched <- character()
            separation_plots <- NULL
            if (supplied) {
                id_column <- methods::slot(object, slots$id)
                held <- lapply(assays, \(assay) as.character(assay[[id_column]]))
                for (assay in names(assays)) {
                    if (!any(listed %in% held[[assay]])) {
                        invalid(sprintf("assay '%s': none of the %d listed control identifiers is a feature of this assay", assay, length(listed)))
                    }
                }
                unmatched <- warnRuvSuppliedListUnmatched(listed, unlist(held, use.names = FALSE),
                    if (identical(binding$omic_type, "lipidomics")) "lipids" else "metabolites")
                parameters$ruv_control_ids <- lapply(held, \(ids) listed[listed %in% ids])
            }
            if (identical(design, "none")) {
                if (!omicsRuvAbsent(parameters$ruv_k) || !omicsRuvAbsent(rule) || !is.null(parameters$ruv_control_ids)) {
                    invalid("the design none takes no k, no k rule and no controls")
                }
                rule <- NA_character_
                factors <- lapply(stats::setNames(names(assays), names(assays)), function(assay) list(assay = assay, status = "none",
                    reason = "the design none: no correction for unwanted variation", reason_class = NA_character_, W = NULL, k = 0L,
                    k_requested = NA_integer_, ruv_k_rule = NA_character_, ruv_k_estimate = NA_integer_, k_estimate = NULL,
                    control_ids = character(), control_count = 0L, controls_left_out = character(), complete_features = NA_integer_,
                    not_tested = character(), ruv_design = "none", ruv_grouping_variable = grouping))
                corrected <- object
            } else {
                counted <- identical(rule, "num_sv_be")
                if (counted && !omicsRuvAbsent(parameters$ruv_k)) invalid("a k is not given with ruv_k_rule num_sv_be: the count gives k")
                # the technical-replicate design is the factors kernel's all-features design on the units (arm 16)
                kernel_design <- if (identical(design, "technical_replicates")) "groups_all_features" else design
                arguments <- list(theObject = object, ruv_design = kernel_design, ruv_k_rule = rule, ruv_control_ids = parameters$ruv_control_ids,
                    ruv_grouping_variable = grouping, ruv_k_seed = seed)
                if (!counted) arguments$ruv_k <- parameters$ruv_k
                # KIT-129 fifth pass item 4: the separation rule, computed here per assay on the design's own controls
                # (`omicsRuvSeparationRuleK`; arms 22 and 5 of simulator S4), when the request gives no k under it.
                if (identical(rule, "sweep_bounded") && omicsRuvAbsent(parameters$ruv_k) &&
                    design %in% c("supplied_controls", "first_pass_nulls")) {
                    separated <- omicsRuvSeparationRuleK(object, design, parameters$ruv_control_ids, grouping)
                    arguments$ruv_k <- separated$k
                    separation_plots <- lapply(separated$plots, function(plot) {
                        if (inherits(plot, c("gg", "ggplot"))) ggplot2::ggplotGrob(stripProtNormQcPlotOwners(plot))
                    })
                }
                result <- do.call(runOmicsRuvFactors, arguments)
                corrected <- result$corrected
                factors <- result$factors
                # KIT-129 integration item 3: under the user's k the count is recorded beside it, per assay, on the
                # values and with the units the count would read (`runOmicsRuvFactors`'s own).
                if (identical(rule, "user")) {
                    design_frame <- as.data.frame(object@design_matrix)
                    design_samples <- as.character(design_frame[[object@sample_id]])
                    for (assay in names(assays)) {
                        values <- omicsRuvAssayValues(assays[[assay]], methods::slot(object, slots$id), design_samples, assay)
                        rows <- design_frame[match(colnames(values), design_samples), , drop = FALSE]
                        technical <- resolveProtTechnicalReplicateBlock(rows, object@group_id, object@technical_replicate_id)
                        beside <- recordRuvFactorCountBesideUserK(values, rows[[grouping]], seed,
                            units = if (isTRUE(technical$technical)) technical$block)
                        factors[[assay]]$ruv_k_estimate <- beside$k_estimate
                        factors[[assay]]$k_estimate_reason <- beside$k_estimate_reason
                        factors[[assay]]$k_estimate <- beside$estimate
                    }
                }
            }
            # KIT-129 fourth pass item 4: the groups of the samples, for the factor-tied-to-groups check of each assay's W.
            group_design <- list(design = object@design_matrix, sample_id = object@sample_id, group_id = object@group_id)
            omic <- binding$omic_type
            object <- NULL
            # KIT-129 third pass item 3: each record names the declared design and its replicate structure (section 1e)
            records <- lapply(factors, function(record) {
                record$k_estimate <- NULL
                record$ruv_design <- design
                record$ruv_replicate_structure <- sets$structure
                record$ruv_group_association <- ruvGroupAssociationCheck(record$W, ruvFactorGroups(record$W, group_design$design,
                    group_design$sample_id, group_design$group_id), omic, design, rule, record$k %||% NA_integer_,
                    record$ruv_k_estimate %||% NA_integer_)
                record
            })
            k <- lapply(factors, `[[`, "k")
            corrected@args$ruv_uncorrected_assays <- assays
            corrected@args$ruv_factors <- records
            seed_used <- if (rule %in% c("num_sv_be", "user")) as.integer(seed) else NA_integer_
            corrected@args$ruvFactorsDeclaration <- list(ruv_design = design, ruv_k_rule = rule, ruv_k_seed = seed_used, k = k,
                tested_matrix = "uncorrected", ruv_replicate_structure = sets$structure, ruv_grouping_variable = grouping)
            final <- store$checkpoint(corrected, paste0(binding$prefix, "-ruv-corrected"))
            corrected <- assays <- NULL
            list(final = final, ruv_design = design, ruv_replicate_structure = sets$structure, ruv_grouping_variable = grouping,
                best_k = k, control_count = lapply(factors, `[[`, "control_count"),
                percentage = lapply(factors, function(record) NA_real_), ruv_k_rule = rule, ruv_k_seed = seed_used,
                ruv_k_estimate = lapply(factors, `[[`, "ruv_k_estimate"), k_estimate = lapply(factors, `[[`, "k_estimate"),
                ruv_applied = lapply(factors, function(record) identical(record$status, "ok")),
                ruv_not_applied_reason = lapply(factors, `[[`, "reason"), factors = records,
                parameters = parameters[declared_names],
                # KIT-129 fifth pass item 3: the supplied list by section 1e's names; per assay its controls are
                # `control_count` (those observed in every sample) and the records' `control_ids`.
                ruv_control_rule = if (supplied) "supplied_list", ruv_control_list_sha256 = if (supplied) parameters$ruv_control_list_sha256,
                ruv_control_list_count = if (supplied) length(listed), ruv_control_list_unmatched = if (supplied) unmatched,
                # KIT-129 fifth pass item 4: the separation rule's plot per assay, where the run computed the rule
                separation_plots = separation_plots)
        })
}

#' Pearson correlations per sample pair within the grouping variable (whole matrix; a product, no state).
#' @noRd
omicsPearsonPairsOperation <- function(binding) {
    newOmicsOperation(binding, paste0(binding$prefix, "_pearson_pairs"), "whole_matrix", parameters = "correlation_group",
        product_type = paste0(binding$prefix, "-pearson-pairs"), output_class = NULL,
        kernel = function(object, parameters, store) {
            pearsonCorForSamplePairsFn <- pearsonCorForSamplePairs
            pearsonCorForSamplePairsFn(theObject = object, correlation_group = parameters$correlation_group)
        })
}

#' The sample selection by the correlation threshold on the pairs product (whole matrix).
#' @noRd
omicsPearsonSelectionOperation <- function(binding) {
    newOmicsOperation(binding, paste0(binding$prefix, "_pearson_selection"), "whole_matrix",
        parameters = c("pairs", "min_pearson_correlation_threshold"),
        product_type = paste0(binding$prefix, ".pearson_selection_result"),
        kernel = function(object, parameters, store) {
            # KIT-078: the declared `pairs` entity, read through the store in either process (it was the undeclared
            # live `correlations`).
            manager <- store$manager
            pairs <- manager$.context$snapshot$products[[parameters$pairs]]
            if (is.null(pairs)) cli::cli_abort("the Pearson pairs product is not in the session", class = "workflow_operation_parameters")
            correlations <- hydrateSavedWorkflowEntity(manager$.context$store, pairs)$value
            filtered <- binding$filter_correlation(object, correlations, parameters$min_pearson_correlation_threshold)
            object <- NULL
            list(final = store$checkpoint(filtered, paste0(binding$prefix, "-correlation-filtered")))
        })
}

#' Differential abundance (whole matrix, product only): the omic's DA kernel once over every contrast and assay, its
#' long tables saved per assay and combined, the kernel's list (less the object) saved as the display product. The
#' GUI never builds the S4 result class; the manifest (`<prefix>.da_results.v1`) is the step's output.
#' @noRd
omicsDaOperation <- function(binding) {
    # The declared names are the DA names every route's record and report overlay read (`formula`, `qvalue`, `lfc`,
    # the proteomics `DA-contrast-output`'s); the kernel receives them under its own argument names.
    long_table <- binding$long_table
    newOmicsOperation(binding, paste0(binding$prefix, "_da_output"), "whole_matrix",
        parameters = c("formula", "qvalue", "lfc", "contrasts"),
        product_type = paste0(binding$prefix, "-da-references"), output_class = NULL,
        fixed = list(eBayes_trend = TRUE, eBayes_robust = TRUE),
        kernel = function(object, parameters, store) {
            result <- binding$run_da(object, parameters)
            object <- NULL
            result$theObject <- NULL
            accepted <- list(schema = paste0(binding$prefix, ".da_references.v1"), source = store$source, tables = list(),
                results = list(), products = list(), qvalue_warnings = result$qvalue_warnings)
            if (is.data.frame(result[[long_table]])) accepted$tables[[long_table]] <- store$table(as.data.frame(result[[long_table]]))
            for (assay in names(result$per_assay_results)) {
                table <- result$per_assay_results[[assay]]
                accepted$results[[assay]] <- list(assay = assay,
                    tables = if (is.data.frame(table)) stats::setNames(list(store$table(as.data.frame(table))), long_table) else list(),
                    rows = if (is.data.frame(table)) nrow(table) else 0L,
                    significant_counts = result$significant_counts[[assay]])
            }
            accepted$products$da_results <- store$product(result, "da_results", paste0(binding$prefix, "-da-results"))
            accepted$rows <- if (is.data.frame(result[[long_table]])) nrow(result[[long_table]]) else 0L
            accepted$columns <- names(result[[long_table]])
            # KIT-128: the assays the kernel could not compute, as it names them.
            if (length(result$failed_assays)) accepted$failed_assays <- result$failed_assays
            accepted$model <- c(list(model = "limma", formula_used = parameters$formula, eBayes_trend = TRUE,
                eBayes_robust = TRUE, treat_applied = isTRUE(parameters$lfc > 0)),
                omicsDaBlockRecord(list(list(
                    block_reason = lapply(result$contrasts_results, function(entry) entry$fit.eb$block_reason),
                    uses_duplicate_correlation = lapply(result$contrasts_results, function(entry) entry$fit.eb$uses_duplicate_correlation)))))
            accepted
        })
}

#' The DA model record's replicate fields (KIT-128), as the proteomics record names them: `block_reason` (the rule's
#' reason; the distinct ones joined when assays differ) and `uses_duplicate_correlation` (whether any assay was fitted
#' with the biological unit as block).
#' @param models Lists with `block_reason` and `uses_duplicate_correlation` per assay (a contrast's, or the kernel's).
#' @return The two fields, or an empty list when no fit recorded a reason (a transcriptomics fit, a fit saved earlier).
#' @noRd
omicsDaBlockRecord <- function(models) {
    reasons <- unique(as.character(unlist(lapply(models, `[[`, "block_reason"))))
    if (!length(reasons)) return(list())
    list(block_reason = paste(reasons, collapse = "; "),
        uses_duplicate_correlation = any(as.logical(unlist(lapply(models, `[[`, "uses_duplicate_correlation")))))
}

#' Identify one assay's uncontrasted limma fit (KIT-076): the values in order, the design and formula, the software.
#' New code: the proteomics `newWorkflowDiaBaseFitIdentity` hashes protein-only inputs.
#' @param data Feature-by-sample matrix.
#' @param design Design table of the run.
#' @param formula_string The model formula.
#' @param assay The assay's name.
#' @param software Installed software identity.
#' @param block The fit's block (KIT-128): the biological unit per sample when the design has technical replicates, else
#'   NULL, which leaves the identity as it was. A fit saved before KIT-128 on such a design is therefore not a hit.
#' @return The store's scientific identity.
#' @noRd
newWorkflowOmicsDaFitIdentity <- function(data, design, formula_string, assay, software, block = NULL, ruv = NULL) {
    newSavedWorkflowScientificIdentity(
        digest::digest(list(values = data, feature_order = rownames(data), sample_order = colnames(data)), algo = "sha256"),
        digest::digest(c(list(formula = formula_string, design = design), if (!is.null(block)) list(block = block)), algo = "sha256"),
        # KIT-129 integration item 1: a fit of a new RUV design names its method (NULL leaves today's identity as it was).
        digest::digest(c(list(layer = "omics_da_base_fit"), if (!is.null(ruv)) list(ruv = ruv)), algo = "sha256"),
        digest::digest(list(assay = assay), algo = "sha256"), software)
}

#' One assay's fit for the contrast operations (KIT-076): the saved `da_fit` when the identity has one, else the fit half
#' of the omic's kernel, saved once. A fit failure is returned as its message (KIT-128; before, the assay was dropped
#' silently) and the caller records it under the contrast's `errors`; ids the store's envelope refuses keep the fit
#' unsaved.
#' @param ruv,fitFn KIT-129 integration item 1: the DA method of a new RUV design, named in the fit's identity (NULL:
#'   today's identity), and the fit half (the binding's by default).
#' @param group_id,technical_replicate_id The object's declared columns (KIT-128), handed to the fit half; the default
#'   NULL keeps the fit and its identity as they were.
#' @return `list(fit, reference)`, or `list(error)` when the fit failed.
#' @noRd
omicsDaAssayFit <- function(binding, store, data, design, formula_string, sample_id, assay, group_id = NULL,
                            technical_replicate_id = NULL, ruv = NULL, fitFn = binding$da_fit) {
    manager <- store$manager
    # KIT-128: the block of the samples the fit reads, by the rule the fit half applies; only a technical one is in the
    # identity.
    block <- if (!is.null(technical_replicate_id) && sample_id %in% names(design)) {
        replicate <- resolveProtTechnicalReplicateBlock(as.data.frame(design)[design[[sample_id]] %in% colnames(data), , drop = FALSE],
            group_id, technical_replicate_id)
        if (isTRUE(replicate$technical)) replicate$block
    }
    identity <- newWorkflowOmicsDaFitIdentity(data, design, formula_string, assay, manager$.context$software, block = block, ruv = ruv)
    factors <- list(route = binding$prefix, layer = "base_fit")
    hit <- resolveSavedWorkflowEntity(manager$.context$store, "da_fit", identity, factors = factors)
    if (identical(hit$status, "hit")) return(list(fit = hit$value$fit, reference = hit$reference))
    failure <- NULL
    fit <- tryCatch(fitFn(data, design, formula_string, sample_id, group_id = group_id,
        technical_replicate_id = technical_replicate_id), error = function(e) {
        cat(sprintf("[ERROR] Processing assay %s: %s\n", assay, e$message))
        logger::log_error(sprintf("   ERROR processing assay %s: %s", assay, e$message))
        failure <<- conditionMessage(e)
        NULL
    })
    if (!is.null(failure)) return(list(error = failure))
    features <- rownames(fit$coefficients)
    if (anyDuplicated(features) || anyNA(features) || !all(nzchar(features))) return(list(fit = fit, reference = NULL))
    reference <- saveWorkflowConsumerEntity(manager$.context$store,
        newSavedWorkflowDaFit(fit, features, rownames(fit$design), colnames(fit$design)),
        basename(tempfile(paste0(binding$prefix, "-da-fit-"))), identity, "da_fit", factors = factors)
    rootWorkflowProduct(manager, reference, list(store$source))
    list(fit = fit, reference = reference)
}

#' Differential abundance for one contrast (KIT-076; whole matrix, product only): per assay the saved fit, the contrast
#' half of the omic's kernel, the orchestrator's table. A failure of this contrast on an assay is recorded under
#' `errors` and costs nothing else; store errors are not caught. The product keeps what the module's list is rebuilt
#' from (`readWorkflowOmicsDaResults`). KIT-095: the fixed settings and the design the kernel halves read are the
#' binding's (`da_fixed`, `da_design`); lipid and metab keep the limma settings and the object's design.
#' @param ruv_design KIT-129 integration item 1. `"legacy"` (the default): today's declaration, unchanged. `"declared"`:
#'   `omicsDaContrastDeclaredOperation` (lipidomics and metabolomics).
#' @noRd
omicsDaContrastOperation <- function(binding, ruv_design = c("legacy", "declared")) {
    if (identical(match.arg(ruv_design), "declared")) return(omicsDaContrastDeclaredOperation(binding))
    long_table <- binding$long_table
    da_design <- binding$da_design %||% function(object) object@design_matrix
    newOmicsOperation(binding, paste0(binding$prefix, "_da_contrast"), "whole_matrix",
        parameters = c("formula", "qvalue", "lfc", "contrast", "friendly_name"),
        product_type = paste0(binding$prefix, "-da-contrast-references"), output_class = NULL,
        fixed = binding$da_fixed %||% list(eBayes_trend = TRUE, eBayes_robust = TRUE),
        kernel = function(object, parameters, store) {
            contrast <- parameters$contrast
            accepted <- list(schema = paste0(binding$prefix, ".da_contrast_reference.v1"), source = store$source,
                contrast = contrast, comparison = parameters$friendly_name, assays = list(), errors = list(),
                qvalue_warnings = list(), products = list())
            assays <- methods::slot(object, binding$data_slot)
            design <- da_design(object)
            preflight <- tryCatch({
                binding$da_preflight(object, data.frame(contrasts = contrast, friendly_names = parameters$friendly_name,
                    stringsAsFactors = FALSE), parameters$formula)
                NULL
            }, error = function(e) conditionMessage(e))
            if (!is.null(preflight)) {
                accepted$errors <- stats::setNames(rep(list(preflight), length(assays)), names(assays))
                return(accepted)
            }
            products <- list()
            for (assay in names(assays)) {
                input <- binding$da_assay_input(assays[[assay]], assay, design, object@sample_id,
                    methods::slot(object, binding$id_slot), object@annotation_id_column)
                if (is.null(input)) next
                fit <- omicsDaAssayFit(binding, store, input$expr_matrix, design, parameters$formula, object@sample_id, assay,
                    group_id = object@group_id, technical_replicate_id = object@technical_replicate_id)
                # KIT-128: a failed fit is named in the errors channel the contrast half uses; the other assays go on.
                if (!is.null(fit$error)) {
                    accepted$errors[[assay]] <- fit$error
                    next
                }
                # KIT-128: which model the assay's fit is (the fit half records it; a trx fit carries none).
                if (!is.null(fit$fit$block_reason)) {
                    accepted$model$block_reason[[assay]] <- fit$fit$block_reason
                    accepted$model$uses_duplicate_correlation[[assay]] <- isTRUE(fit$fit$uses_duplicate_correlation)
                }
                outcome <- tryCatch({
                    applied <- binding$da_apply(fit$fit, contrast, if (parameters$lfc > 0) parameters$lfc else NA, TRUE, TRUE)
                    list(applied = applied, table = binding$da_format(applied$results[[contrast]], contrast, assay,
                        parameters$friendly_name, input$annotation_info, input$expr_matrix, design, object@sample_id,
                        object@group_id, parameters$qvalue, parameters$lfc))
                }, error = function(e) {
                    logger::log_error(sprintf("   ERROR in contrast %s, assay %s: %s", contrast, assay, conditionMessage(e)))
                    list(error = conditionMessage(e))
                })
                if (!is.null(outcome$error)) {
                    accepted$errors[[assay]] <- outcome$error
                    next
                }
                table <- as.data.frame(outcome$table)
                accepted$assays[[assay]] <- list(assay = assay, tables = stats::setNames(list(store$table(table)), long_table),
                    rows = nrow(table), fit = fit$reference, significant_counts = binding$da_count(outcome$table))
                if (length(outcome$applied$qvalue_warnings)) accepted$qvalue_warnings[[assay]] <- outcome$applied$qvalue_warnings
                products[[assay]] <- list(results = outcome$applied$results[[contrast]], table = outcome$table,
                    fit.eb = outcome$applied$fit.eb)
            }
            object <- NULL
            accepted$products$da_contrast <- store$product(list(assays = products), "da_contrast", paste0(binding$prefix, "-da-contrast"))
            accepted$rows <- sum(vapply(accepted$assays, function(entry) entry$rows, integer(1)))
            accepted$columns <- if (length(products)) names(products[[1L]]$table)
            accepted
        })
}

#' The DA of a new-design RUV state for one contrast (KIT-129 integration item 1; lipidomics and metabolomics): today's
#' contrast declaration (`omicsDaContrastOperation`), with the method declared as the parameter `ruv_da_method`:
#' per assay the uncorrected values (`args$ruv_uncorrected_assays`) with the assay's W in the model
#' (`applyOmicsRuvFactorsToDaInput`), each comparison's standard error exact (decision M29,
#' `omicsRuvExactContrastsFitFn`), the fit's identity naming the method; an assay without factors is tested as
#' `runLipidsDA` / `runMetabolitesDA` test it with its records (decision M10: the corrected matrix is for display).
#' A separate declaration so that today's kernel is unchanged.
#' @noRd
omicsDaContrastDeclaredOperation <- function(binding) {
    long_table <- binding$long_table
    da_design <- binding$da_design %||% function(object) object@design_matrix
    newOmicsOperation(binding, paste0(binding$prefix, "_da_contrast"), "whole_matrix",
        parameters = c("formula", "qvalue", "lfc", "contrast", "friendly_name", "ruv_da_method"),
        product_type = paste0(binding$prefix, "-da-contrast-references"), output_class = NULL,
        fixed = binding$da_fixed %||% list(eBayes_trend = TRUE, eBayes_robust = TRUE),
        kernel = function(object, parameters, store) {
            ruv_factors <- object@args$ruv_factors
            if (!is.list(ruv_factors) || !identical(workflowRuvDaMethodOfArgs(object@args), parameters$ruv_da_method)) {
                cli::cli_abort("the DA of a new RUV design needs a source that declares it, and its method",
                    class = "workflow_execution_invalid", call = NULL)
            }
            contrast <- parameters$contrast
            accepted <- list(schema = paste0(binding$prefix, ".da_contrast_reference.v1"), source = store$source,
                contrast = contrast, comparison = parameters$friendly_name, assays = list(), errors = list(),
                qvalue_warnings = list(), products = list())
            assays <- methods::slot(object, binding$data_slot)
            design <- da_design(object)
            preflight <- tryCatch({
                binding$da_preflight(object, data.frame(contrasts = contrast, friendly_names = parameters$friendly_name,
                    stringsAsFactors = FALSE), parameters$formula)
                NULL
            }, error = function(e) conditionMessage(e))
            if (!is.null(preflight)) {
                accepted$errors <- stats::setNames(rep(list(preflight), length(assays)), names(assays))
                return(accepted)
            }
            products <- list()
            for (assay in names(assays)) {
                tested <- object@args$ruv_uncorrected_assays[[assay]] %||% assays[[assay]]
                input <- binding$da_assay_input(tested, assay, design, object@sample_id,
                    methods::slot(object, binding$id_slot), object@annotation_id_column)
                if (is.null(input)) next
                # KIT-127's DA input: the supplied controls left out, W in the design and the formula (M28's checks).
                ruv_input <- tryCatch(applyOmicsRuvFactorsToDaInput(input$expr_matrix, design, parameters$formula,
                    object@sample_id, ruv_factors[[assay]], assay, object@group_id, object@technical_replicate_id),
                    error = function(e) list(error = conditionMessage(e)))
                if (!is.null(ruv_input$error)) {
                    accepted$errors[[assay]] <- ruv_input$error
                    next
                }
                fit <- omicsDaAssayFit(binding, store, ruv_input$expr_matrix, ruv_input$design_matrix, ruv_input$formula_string,
                    object@sample_id, assay, group_id = object@group_id, technical_replicate_id = object@technical_replicate_id,
                    ruv = parameters$ruv_da_method, fitFn = function(...) withOmicsRuvKOnBlockFailure(binding$da_fit(...), ruv_factors[[assay]]))
                # KIT-128: a failed fit is named in the errors channel the contrast half uses; the other assays go on.
                if (!is.null(fit$error)) {
                    accepted$errors[[assay]] <- fit$error
                    next
                }
                # KIT-128: which model the assay's fit is (the fit half records it; a trx fit carries none).
                if (!is.null(fit$fit$block_reason)) {
                    accepted$model$block_reason[[assay]] <- fit$fit$block_reason
                    accepted$model$uses_duplicate_correlation[[assay]] <- isTRUE(fit$fit$uses_duplicate_correlation)
                }
                # Decision M29: an assay with W in the model forms each comparison's standard error exactly.
                contrasts_fit <- if (is.null(ruv_factors[[assay]]$W)) limma::contrasts.fit else omicsRuvExactContrastsFitFn(fit$fit,
                    ruv_input$expr_matrix, ruv_input$design_matrix, object@sample_id, object@group_id, object@technical_replicate_id)
                outcome <- tryCatch({
                    applied <- binding$da_apply(fit$fit, contrast, if (parameters$lfc > 0) parameters$lfc else NA, TRUE, TRUE,
                        contrastsFitFn = contrasts_fit)
                    list(applied = applied, table = binding$da_format(applied$results[[contrast]], contrast, assay,
                        parameters$friendly_name, input$annotation_info, ruv_input$expr_matrix, design, object@sample_id,
                        object@group_id, parameters$qvalue, parameters$lfc))
                }, error = function(e) {
                    logger::log_error(sprintf("   ERROR in contrast %s, assay %s: %s", contrast, assay, conditionMessage(e)))
                    list(error = conditionMessage(e))
                })
                if (!is.null(outcome$error)) {
                    accepted$errors[[assay]] <- outcome$error
                    next
                }
                table <- as.data.frame(outcome$table)
                accepted$assays[[assay]] <- list(assay = assay, tables = stats::setNames(list(store$table(table)), long_table),
                    rows = nrow(table), fit = fit$reference, significant_counts = binding$da_count(outcome$table))
                if (length(outcome$applied$qvalue_warnings)) accepted$qvalue_warnings[[assay]] <- outcome$applied$qvalue_warnings
                products[[assay]] <- list(results = outcome$applied$results[[contrast]], table = outcome$table,
                    fit.eb = outcome$applied$fit.eb)
                # KIT-127: what the new design did not test, by assay.
                accepted$ruv_features_not_tested[[assay]] <- ruv_input$not_tested
            }
            accepted$model <- c(accepted$model, list(ruv_da_method = parameters$ruv_da_method, tested_matrix = "uncorrected",
                w_terms = lapply(ruv_factors, function(record) colnames(record$W) %||% character())))
            object <- NULL
            accepted$products$da_contrast <- store$product(list(assays = products), "da_contrast", paste0(binding$prefix, "-da-contrast"))
            accepted$rows <- sum(vapply(accepted$assays, function(entry) entry$rows, integer(1)))
            accepted$columns <- if (length(products)) names(products[[1L]]$table)
            accepted
        })
}

#' The contrast strings and friendly names of a contrasts table, as the orchestrators read them.
#' @noRd
omicsDaContrastNames <- function(contrasts) {
    strings <- if ("contrasts" %in% colnames(contrasts)) {
        as.character(contrasts$contrasts)
    } else if ("contrast_string" %in% colnames(contrasts)) {
        as.character(contrasts$contrast_string)
    } else {
        stop("contrasts_tbl must have a 'contrasts' or 'contrast_string' column")
    }
    list(strings = strings, friendly = if ("friendly_names" %in% colnames(contrasts)) as.character(contrasts$friendly_names) else strings)
}

#' The kernel's list rebuilt from a v2 manifest (KIT-076), in the kernel's order: assay-major, then contrast. The
#' exports, the Glimma pages, the volcanoes, the heatmaps, the selectors and the tables read this list. `fit.eb` here
#' is a coefficient view (one column per contrast): its consumers read the names and the row order; the full
#' per-contrast fit stays in the saved product. `failed_contrasts` names what the run could not compute.
#' @param manifest A `<prefix>.da_results.v2` manifest.
#' @return The list `runLipidsDA` / `runMetabolitesDA` return, without the source object.
#' @noRd
readWorkflowOmicsDaResults <- function(binding, manager, manifest) {
    # A contrast the preflight refused saved no product.
    products <- lapply(manifest$results, function(entry) {
        if (is.null(entry$products$da_contrast)) return(list())
        hydrateSavedWorkflowEntity(manager$.context$store, entry$products$da_contrast)$value$assays
    })
    per_assay <- list()
    contrasts_results <- list()
    warnings <- list()
    for (assay in manifest$assays) {
        present <- Filter(Negate(is.null), lapply(products, `[[`, assay))
        if (!length(present)) next
        per_assay[[assay]] <- dplyr::bind_rows(lapply(present, `[[`, "table"))
        contrasts_results[[assay]] <- list(results = lapply(present, `[[`, "results"),
            fit.eb = list(coefficients = do.call(cbind, lapply(unname(present), function(entry) entry$fit.eb$coefficients))),
            qvalue_warnings = do.call(c, c(list(list()), lapply(unname(manifest$results), function(entry) entry$qvalue_warnings[[assay]]))))
        if (length(contrasts_results[[assay]]$qvalue_warnings)) warnings[[assay]] <- contrasts_results[[assay]]$qvalue_warnings
    }
    failed <- Filter(length, lapply(manifest$results, `[[`, "errors"))
    c(list(contrasts_results = contrasts_results), stats::setNames(list(dplyr::bind_rows(per_assay)), binding$long_table),
        list(per_assay_results = per_assay, significant_counts = lapply(per_assay, binding$da_count), qvalue_warnings = warnings,
            da_q_val_thresh = manifest$parameters$qvalue, treat_lfc_cutoff = manifest$parameters$lfc, failed_contrasts = failed))
}

#' What a managed DA tab holds (KIT-082): the kernel's list without its per-sample columns and without the per-assay tables
#' and fits, which only the worker's volcano and export read. Everything the tab shows comes from these rows, unchanged.
#' @param binding The omic's binding.
#' @param results `readWorkflowOmicsDaResults`'s list.
#' @return The list, with `download_columns`: the whole table's column types by name, for the download.
#' @noRd
omicsDaTabResults <- function(binding, results) {
    long <- results[[binding$long_table]]
    keep <- !grepl("^intensity[.]", names(long))
    results$contrasts_results <- results$per_assay_results <- NULL
    results[[binding$long_table]] <- long[, keep, drop = FALSE]
    results$download_columns <- vapply(long, typeof, character(1))
    results
}

#' Write the DA download a product at a time (KIT-082), the same bytes as `write.csv` on the whole table, with one contrast
#' and assay held at once.
#' @param binding The omic's binding.
#' @param manager The session owner.
#' @param manifest The accepted DA manifest.
#' @param columns The whole table's column types by name: `bind_rows` gives every piece of a column one type, and the file's
#'   numbers are written by type.
#' @param file The download's path.
#' @noRd
writeWorkflowOmicsDaDownload <- function(binding, manager, manifest, columns, file) {
    first <- TRUE
    for (assay in manifest$assays) for (entry in manifest$results) {
        if (is.null(entry$products$da_contrast)) next
        table <- hydrateSavedWorkflowEntity(manager$.context$store, entry$products$da_contrast)$value$assays[[assay]][["table"]]
        if (is.null(table) || !nrow(table)) next
        for (absent in setdiff(names(columns), names(table))) table[[absent]] <- NA
        table <- as.data.frame(table)[, names(columns), drop = FALSE]
        for (column in names(columns)) storage.mode(table[[column]]) <- columns[[column]]
        utils::write.table(table, file, sep = ",", dec = ".", qmethod = "double", row.names = FALSE, col.names = first, append = !first)
        first <- FALSE
    }
    # No rows anywhere: the header, as `write.csv` writes for an empty table.
    if (first) writeLines(paste0('"', gsub('"', '""', names(columns)), '"', collapse = ","), file)
    invisible(file)
}

#' The steps of the lipidomics and metabolomics routes (a binding's `steps`, KIT-094), in registry order.
#' `<prefix>_da_output` is no longer run (KIT-076): it stays registered so v1 records and manifests keep resolving their
#' operation.
#' @param binding The omic's binding.
#' @return The declarations.
#' @noRd
omicsRouteSteps <- function(binding) {
    prefix <- binding$prefix
    list(omicsIntensityFilterOperation(binding), omicsDuplicateResolutionOperation(binding),
        omicsStoreOperation(binding, paste0(prefix, "_qc_finalize")), omicsItsdNormalisationOperation(binding),
        omicsLog2Operation(binding), omicsBetweenSampleNormalisationOperation(binding), omicsRuvOperation(binding),
        omicsPearsonPairsOperation(binding), omicsPearsonSelectionOperation(binding),
        omicsStoreOperation(binding, paste0(prefix, "_norm_complete")), omicsStoreOperation(binding, paste0(prefix, "_reset")),
        omicsDaOperation(binding), omicsDaContrastOperation(binding))
}

#' Every declaration of one omic, by name: the binding's own steps (KIT-094).
#' @noRd
omicsWorkflowOperations <- function(binding) {
    declarations <- binding$steps(binding)
    stats::setNames(declarations, vapply(declarations, `[[`, "", "name"))
}

#' The omics bindings the registry and the record read, by omic type.
#' @noRd
omicsDeclarationBindings <- function() {
    list(lipidomics = lipidomicsDeclarationBinding, metabolomics = metabolomicsDeclarationBinding,
        transcriptomics = transcriptomicsDeclarationBinding)
}

#' The bindings' omic types (KIT-094): the session's type allowlists and routes take an omic's from its binding
#' (`toupper` of it is the workflow type).
#' @noRd
omicsWorkflowTypes <- function() vapply(omicsDeclarationBindings(), function(binding) binding$omic_type, character(1))

#' The worker's publication operations of the lipidomics and metabolomics routes (a binding's `publications`, KIT-094):
#' `c(<name after the prefix>, <handler>)` each. KIT-078: what the steps' figures and exports read (the matrix, the
#' composite's rasters), in the worker; they commit nothing. KIT-082: what the tabs read from the whole matrix outside
#' their steps, in the worker too.
#' @noRd
omicsRoutePublications <- list(c("qc_figures", "computeWorkflowOmicsQcFiguresAction"),
    c("composite_figure", "computeWorkflowOmicsCompositeFigureAction"), c("da_publication", "computeWorkflowOmicsDaPublicationAction"),
    c("da_tab_results", "computeWorkflowOmicsDaTabResultsAction"), c("filtering_metrics", "computeWorkflowOmicsFilteringMetricsAction"),
    c("final_pca", "computeWorkflowOmicsFinalPcaAction"), c("volcano_page", "computeWorkflowOmicsVolcanoPageAction"),
    c("session_export", "computeWorkflowOmicsSessionExportAction"), c("integration_export", "computeWorkflowOmicsIntegrationExportAction"))

#' The realised outcomes the record keeps for a lipidomics or metabolomics step (a binding's `record_outcomes`,
#' KIT-094). KIT-031: the RUV step's realised result (k and the control count per assay, their totals) and the
#' internal-standard step's application, in the words the lipid report's RUV and ITSD sentences use. KIT-032: the
#' metabolomics RUV and ITSD steps record the same outcomes under their own operation names.
#' @param binding The omic's binding.
#' @param context The owner's context.
#' @param operation The step's operation, one of the binding's.
#' @param record The state's record.
#' @param outcomes The outcomes every step has.
#' @return `outcomes`, with the step's own.
#' @noRd
omicsRouteRecordOutcomes <- function(binding, context, operation, record, outcomes) {
    if (identical(operation, paste0(binding$prefix, "_ruv"))) {
        product <- record$data$products[[paste0(operation, "_science")]]
        saved <- if (isSavedWorkflowReference(product)) tryCatch({
            hydrated <- hydrateSavedWorkflowEntity(context$store, product)
            if (is.list(hydrated$value)) hydrated$value else hydrated
        }, error = function(error) NULL) else NULL
        k <- unlist(saved$best_k); controls <- unlist(saved$control_count); percentage <- unlist(saved$percentage)
        outcomes$ruv_applied <- TRUE
        outcomes$k <- if (length(k)) sum(k, na.rm = TRUE) else NA_real_
        outcomes$k_per_assay <- as.list(saved$best_k)
        outcomes$negative_control_count <- if (length(controls)) sum(controls, na.rm = TRUE) else NA_integer_
        outcomes$negative_control_count_per_assay <- as.list(saved$control_count)
        outcomes$negative_control_percentage <- if (length(percentage) && !all(is.na(percentage))) mean(percentage, na.rm = TRUE) else NA_real_
        outcomes$negative_control_percentage_per_assay <- as.list(saved$percentage)
        # KIT-129 integration item 3: a new design's product records its design and its k per assay: whether RUV was
        # applied, the rule, the estimate and the seed (today's product names none, and its outcomes are as before).
        if (!is.null(saved$ruv_design)) {
            outcomes$ruv_applied <- any(unlist(saved$ruv_applied))
            outcomes[c("ruv_design", "ruv_k_rule", "ruv_k_seed")] <- list(saved$ruv_design, saved$ruv_k_rule, saved$ruv_k_seed)
            # KIT-129 third pass item 3: the replicate structure and its column (products before it name neither)
            if (!is.null(saved$ruv_replicate_structure)) {
                outcomes[c("ruv_replicate_structure", "ruv_grouping_variable")] <- list(saved$ruv_replicate_structure, saved$ruv_grouping_variable)
            }
            outcomes$ruv_applied_per_assay <- as.list(saved$ruv_applied)
            outcomes$ruv_k_estimate_per_assay <- as.list(saved$ruv_k_estimate)
            outcomes$ruv_not_applied_reason_per_assay <- as.list(saved$ruv_not_applied_reason)
            # KIT-129 fourth pass item 4: the factor-tied-to-groups check per assay (statistic, threshold, outcome).
            for (field in ruvGroupAssociationRecordFields) {
                outcomes[[paste0(field, "_per_assay")]] <- lapply(saved$factors, function(record) record$ruv_group_association[[field]])
            }
            # KIT-129 sixth pass item 4: the new design's fields, one set in every omic; what differs by assay (k, the
            # estimate's permutations, the controls used, the list's found count) per assay.
            listed <- saved$ruv_control_list_count
            outcomes[RUV_DESIGN_RECORD_FIELDS] <- ruvDesignRecordFields(ruv_design = saved$ruv_design,
                ruv_replicate_structure = saved$ruv_replicate_structure, ruv_grouping_variable = saved$ruv_grouping_variable,
                ruv_k_rule = saved$ruv_k_rule, ruv_k_typed = if (identical(saved$ruv_k_rule, "user")) saved$parameters$ruv_k,
                ruv_k_seed = saved$ruv_k_seed, ruv_control_rule = saved$ruv_control_rule %||% switch(saved$ruv_design %||% "",
                    groups_all_features = , technical_replicates = "all_features", first_pass_nulls = "first_pass_drop_top_25", NULL),
                ruv_control_list_sha256 = saved$ruv_control_list_sha256, ruv_control_list_count = listed,
                ruv_control_list_found = if (length(listed)) listed - length(saved$ruv_control_list_unmatched),
                ruv_first_pass_rule = if (identical(saved$ruv_design, "first_pass_nulls")) "limma_smallest_pairwise_p",
                ruv_da_method = workflowRuvDaMethodOfArgs(record$data$slots$args))
            outcomes$ruv_k_per_assay <- as.list(saved$best_k)
            outcomes$ruv_control_count_per_assay <- as.list(saved$control_count)
            outcomes$ruv_k_permutations_per_assay <- lapply(saved$k_estimate, function(estimate) estimate$permutations %||% NA_integer_)
        }
    }
    if (identical(operation, paste0(binding$prefix, "_itsd_normalisation"))) outcomes$itsd_applied <- TRUE
    outcomes
}

# ---------------------------------------------------------------------------------------------------------------
# The wrappers: the operation on a source, its product on the checkpoint, the commit under the module's state name.
# ---------------------------------------------------------------------------------------------------------------

#' Run one state-producing omics declaration and commit its checkpoint (KIT-031; KIT-032: both omics).
#' @param manager The session owner (`DIAStateAdapter`).
#' @param declaration The declaration.
#' @param parameters Every declared parameter and any extra the kernel reads.
#' @param state_name,config,description The module's commit.
#' @param source Saved reference of the input state; the head when NULL.
#' @return `list(run, reference, value, object)`: the operation run, the committed state's reference, the kernel's
#'   value, and the committed object only when `hydrate` asks for it (KIT-075: the modules hold the reference and read
#'   per action; the argument stays for the harness's recorders).
#' @noRd
runWorkflowOmicsReferenceOperation <- function(manager, declaration, parameters, state_name, config, description, source = NULL,
    dependencies = list(), hydrate = FALSE) {
    if (!inherits(manager, "DIAStateAdapter")) cli::cli_abort("the step needs the session's state owner", class = "saved_entity_scope")
    withWorkflowStateScope(manager, function() {
        source <- source %||% getWorkflowStateReference(manager)
        run <- runSavedOperation(manager, declaration, source, parameters, dependencies)
        final <- readWorkflowStateMetadata(manager, run$value$final)$data
        final$products[[paste0(declaration$name, "_science")]] <- run$reference
        commitWorkflowReferenceState(manager, state_name, final, config, description)
        list(run = run, reference = getWorkflowStateReference(manager), value = run$value,
            object = if (hydrate) manager$getState(state_name) else NULL)
    })
}

#' Read an omics state for one action (KIT-075): the binding's components of the state, read in the owner's scope.
#' The caller uses the object within the action and holds only references between actions.
#' @param reference The state's saved reference or name; the head when NULL.
#' @return The transient object, or NULL when the owner has no such state.
#' @noRd
readWorkflowOmicsState <- function(binding, manager, reference = NULL) {
    if (!inherits(manager, "DIAStateAdapter")) cli::cli_abort(paste(binding$step_scope, "needs the session's state owner"), class = "saved_entity_scope")
    if (is.character(reference)) reference <- getWorkflowStateReference(manager, reference)
    reference <- reference %||% getWorkflowStateReference(manager)
    if (is.null(reference)) return(NULL)
    withWorkflowStateScope(manager, function() readWorkflowExecutionComponents(manager, reference, fields = binding$fields))
}

#' An omics state's metadata (KIT-075): the stored slots, the assay names and row counts (from the stored tables'
#' shapes) and the design (a sample table), read without any assay table.
#' @param reference The state's saved reference or name; the head when NULL.
#' @return `list(reference, class, slots, assay_names, assay_rows, design_matrix)`, or NULL when the owner has no such
#'   state.
#' @noRd
readWorkflowOmicsMetadata <- function(binding, manager, reference = NULL) {
    if (!inherits(manager, "DIAStateAdapter")) return(NULL)
    if (is.character(reference)) reference <- getWorkflowStateReference(manager, reference)
    reference <- reference %||% getWorkflowStateReference(manager)
    if (is.null(reference)) return(NULL)
    withWorkflowStateScope(manager, function() {
        spec <- readWorkflowStateMetadata(manager, reference)$data
        assays <- spec$components[[binding$data_slot]]
        list(reference = reference, class = spec$class, slots = spec$slots, assay_names = names(assays),
            assay_rows = vapply(assays, function(table) as.integer(table$metadata$shape[[1L]]), integer(1)),
            design_matrix = materializeWorkflowTable(manager$.context, spec$components$design_matrix))
    })
}

#' An omics state without its matrix (KIT-082): the slots and the design, each assay as its zero-row prototype, with the
#' stored row counts and, when asked, what the QC tabs count from the tables, read a page at a time with nothing whole kept.
#' @param reference The state's saved reference or name; the head when NULL.
#' @param measures Any of `features` (unique ids per assay) and `missing` (percent of cells NA or zero per assay).
#' @return The object, with the attributes `assay_rows`, `feature_counts` and `missing_percent`, or NULL with no such state.
#' @noRd
readWorkflowOmicsShell <- function(binding, manager, reference = NULL, measures = character()) {
    if (!inherits(manager, "DIAStateAdapter")) cli::cli_abort(paste(binding$step_scope, "needs the session's state owner"), class = "saved_entity_scope")
    if (is.character(reference)) reference <- getWorkflowStateReference(manager, reference)
    reference <- reference %||% getWorkflowStateReference(manager)
    if (is.null(reference)) return(NULL)
    withWorkflowStateScope(manager, function() {
        shell <- readWorkflowExecutionComponents(manager, reference, fields = binding$design_field)
        tables <- readWorkflowStateMetadata(manager, reference)$data$components[[binding$data_slot]]
        attr(shell, "assay_rows") <- vapply(tables, function(table) as.integer(table$metadata$shape[[1L]]), integer(1))
        id_column <- methods::slot(shell, binding$id_slot)
        if ("features" %in% measures) attr(shell, "feature_counts") <- vapply(omicsVisitAssays(binding, manager, reference, function(assay, scope) {
            if (id_column %in% scope$metadata$columns) length(unique(omicsReadRows(scope, id_column)[[id_column]])) else as.integer(scope$metadata$rows)
        }), identity, integer(1))
        if ("missing" %in% measures) attr(shell, "missing_percent") <- vapply(omicsVisitAssays(binding, manager, reference,
            function(assay, scope) omicsAssayMissingPercent(scope)), identity, numeric(1))
        shell
    })
}

#' Visit each assay table of a state through a bounded read scope.
#' @param visit `function(assay, scope)`.
#' @return The values, named by assay.
#' @noRd
omicsVisitAssays <- function(binding, manager, reference, visit, assays = NULL) {
    tables <- readWorkflowStateMetadata(manager, reference)$data$components[[binding$data_slot]]
    names_to_read <- if (is.null(assays)) names(tables) else intersect(assays, names(tables))
    stats::setNames(lapply(names_to_read, function(assay) {
        component <- tables[[assay]]
        scope <- newWorkflowTableReadScope(manager$.context$store, component, reader = if (component$metadata$shape[[2L]] > 64L) "duckdb" else "arrow")
        tryCatch(visit(assay, scope), finally = releaseWorkflowTableReadScope(scope))
    }), names_to_read)
}

#' The rows of a table that a filter keeps, read a page at a time: the filter's rows are all that is held.
#' @param keep `function(page)` returning a logical per row; every row when NULL.
#' @noRd
omicsReadRows <- function(scope, columns, keep = NULL) {
    pieces <- list()
    streamWorkflowTableRows(scope, columns, function(page, offset, row_names) {
        rownames(page) <- NULL
        pieces[[length(pieces) + 1L]] <<- if (is.function(keep)) page[keep(page), , drop = FALSE] else page
        invisible(NULL)
    })
    if (!length(pieces)) return(readSavedWorkflowTablePrototype(scope)[columns])
    rows <- do.call(rbind, pieces)
    rownames(rows) <- NULL
    rows
}

#' Percent of the numeric cells that are NA or zero, as the QC tabs' assay statistics count them.
#' @noRd
omicsAssayMissingPercent <- function(scope) {
    prototype <- readSavedWorkflowTablePrototype(scope)
    columns <- names(prototype)[vapply(prototype, is.numeric, logical(1))]
    if (!length(columns)) return(NA_real_)
    missing <- 0
    if (scope$metadata$rows) streamWorkflowTableRows(scope, columns, function(page, offset, row_names) {
        for (column in columns) missing <<- missing + sum(is.na(page[[column]]) | page[[column]] == 0)
        invisible(NULL)
    })
    round((missing / (scope$metadata$rows * length(columns))) * 100, 1)
}

#' An omics state whose assays hold only the rows a filter keeps (KIT-082), for a kernel that reads those rows only.
#' @param keep `function(page, id_column)` returning a logical per row.
#' @return The object, or NULL with no such state.
#' @noRd
readWorkflowOmicsRows <- function(binding, manager, reference = NULL, keep) {
    shell <- readWorkflowOmicsShell(binding, manager, reference)
    if (is.null(shell)) return(NULL)
    if (is.character(reference)) reference <- getWorkflowStateReference(manager, reference)
    reference <- reference %||% getWorkflowStateReference(manager)
    id_column <- methods::slot(shell, binding$id_slot)
    withWorkflowStateScope(manager, function() {
        methods::slot(shell, binding$data_slot, check = FALSE) <- omicsVisitAssays(binding, manager, reference, function(assay, scope) {
            omicsReadRows(scope, scope$metadata$columns, function(page) keep(page, id_column))
        })
        shell
    })
}

#' An omics state whose assays hold only the id column (KIT-082): duplicate detection and the study parameters read no other.
#' @return The object, or NULL with no such state.
#' @noRd
readWorkflowOmicsIds <- function(binding, manager, reference = NULL) {
    shell <- readWorkflowOmicsShell(binding, manager, reference)
    if (is.null(shell)) return(NULL)
    if (is.character(reference)) reference <- getWorkflowStateReference(manager, reference)
    reference <- reference %||% getWorkflowStateReference(manager)
    id_column <- methods::slot(shell, binding$id_slot)
    withWorkflowStateScope(manager, function() {
        methods::slot(shell, binding$data_slot, check = FALSE) <- omicsVisitAssays(binding, manager, reference, function(assay, scope) {
            if (id_column %in% scope$metadata$columns) return(omicsReadRows(scope, id_column))
            as.data.frame(matrix(nrow = scope$metadata$rows, ncol = 0L))
        })
        shell
    })
}

#' The ITSD candidate tables, built a page at a time by the unchanged builder and sorted as it sorts: its statistics are per
#' row, and a stable sort of stably sorted pages keeps the whole table's order.
#' @param build `function(page)`, the omic's builder.
#' @param reference The state; the head when NULL.
#' @return The tables, named by assay, for the assays the state holds.
#' @noRd
readWorkflowOmicsItsdTables <- function(binding, manager, assays, build, reference = NULL) {
    if (is.character(reference)) reference <- getWorkflowStateReference(manager, reference)
    reference <- reference %||% getWorkflowStateReference(manager)
    withWorkflowStateScope(manager, function() {
        omicsVisitAssays(binding, manager, reference, function(assay, scope) {
            pieces <- list()
            streamWorkflowTableRows(scope, scope$metadata$columns, function(page, offset, row_names) {
                rownames(page) <- NULL
                pieces[[length(pieces) + 1L]] <<- build(page)
                invisible(NULL)
            })
            if (!length(pieces)) return(build(readSavedWorkflowTablePrototype(scope)))
            if (length(pieces) == 1L) return(pieces[[1L]])
            dplyr::arrange(do.call(rbind, pieces), dplyr::desc(.data$is_candidate), .data$cv_percent)
        }, assays)
    })
}

#' The patterns the internal-standard kernels fall back to when the user's pattern finds nothing.
#' @noRd
omicsInternalStandardFallbackPatterns <- function() {
    c(
        "^IS[_-]", # IS_ or IS- prefix (e.g., IS_Caffeine, IS-Leucine)
        "^ITSD[_-]", # ITSD_ prefix (e.g., ITSD_A, ITSD_001)
        "ISTD", # ISTD anywhere (e.g., Caffeine-ISTD, ISTD_001)
        "ITSD", # ITSD anywhere (alternate spelling)
        "[_-]d\\d+$", # Deuterated suffix (e.g., Caffeine_d3, Leucine-d10)
        "[_-]d\\d+[_-]", # Deuterated mid-name (e.g., Caffeine_d3_IS)
        "\\(d\\d+\\)", # Deuterated in parens (e.g., Caffeine(d3))
        "(?i)internal.?standard", # "Internal Standard" / "Internal_Standard"
        "(?i)^is\\d+$", # IS followed by numbers (e.g., IS1, IS23)
        "13C[_-]?labeled", # 13C labeled (e.g., 13C-labeled_reference)
        "15N[_-]?labeled", # 15N labeled
        "-13C\\d*-", # 13C in middle (e.g., Glucose-13C6-IS)
        "-15N\\d*-" # 15N in middle
    )
}

#' The rows an internal-standard kernel can use: those matching the user's pattern or the fallback patterns in any column it
#' tries. A row whose value is NA is kept, because the kernels count a match of NA.
#' @param page One page of an assay.
#' @param pattern The user's pattern.
#' @param columns The columns the kernel tries.
#' @return A logical per row.
#' @noRd
omicsInternalStandardRows <- function(page, pattern, columns) {
    keep <- rep(FALSE, nrow(page))
    for (column in intersect(columns, names(page))) {
        values <- as.character(page[[column]])
        for (candidate in c(pattern, paste(omicsInternalStandardFallbackPatterns(), collapse = "|"))) {
            hit <- tryCatch(stringr::str_detect(values, candidate), error = function(error) rep(FALSE, length(values)))
            keep <- keep | is.na(hit) | hit
        }
    }
    keep
}

#' DA results whose heatmap reads the rows it has chosen (KIT-082): `read_source(ids)` gives the source state with only those rows.
#' @noRd
omicsDaResultsWithRowReader <- function(binding, results, daData, workflowData) {
    manager <- if (is.list(workflowData) || is.environment(workflowData)) workflowData$state_manager
    if (is.null(results) || is.null(manager)) return(results)
    state <- daData$current_state_name
    results$read_source <- function(ids) readWorkflowOmicsRows(binding, manager, state, function(page, id_column) page[[id_column]] %in% ids)
    results
}

#' Bind the omics session owner to its route at the import (KIT-075, the proteomics pattern): binding replaces the
#' pending store, so the import's tables can only be held by reference once the owner is bound. An owner that is
#' already bound, or a bus without an owner, is left as it is.
#' @noRd
bindWorkflowOmicsOwnerAtImport <- function(workflowData, workflow_type, format) {
    if (isDIASessionOwnerPending(workflowData$state_manager)) bindWorkflowOwnerTypeAtDesign(workflowData, workflow_type, format)
    invisible(workflowData$state_manager)
}

#' Hold an import or design table list on the bus by reference, in the omic's owner (KIT-075).
#' @noRd
setWorkflowOmicsBusTables <- function(workflowData, field, value) {
    setWorkflowOwnedPayload(workflowData, field, value, manager = workflowData$state_manager)
}

#' Read an import or design table list from the bus for one action (KIT-075).
#' @noRd
getWorkflowOmicsBusTables <- function(workflowData, field) {
    getWorkflowOwnedPayload(workflowData, field, manager = workflowData$state_manager)
}

#' A design builder's table-list holder (KIT-075): the get/set contract of a `reactiveVal`, holding a reference in the
#' owner `owner()` returns at each use (the builder starts before the import binds the owner); plain without an owner.
#' @noRd
omicsBuilderTableVal <- function(label, owner) {
    held <- shiny::reactiveVal()
    function(value) {
        if (missing(value)) return(resolveWorkflowOwnedValue(held(), owner()))
        held(saveWorkflowOwnedValue(value, label, owner()))
        invisible(NULL)
    }
}

#' One tab step, as a wrapper runs it in process and a managed action dispatches it (KIT-078): the declaration's name, its
#' parameters, and the module's state name and description for the commit.
#' KIT-217: the per-group count's two settings (values required per group, groups required; defaults 2 and 1) are
#' declared parameters too, so a step always records them; the replay of a record made before passes NULL for both.
#' @return `list(operation, parameters, state_name, description)`.
#' @noRd
omicsIntensityFilterStep <- function(binding, percentile, proportion, description = NULL, state_name = NULL,
    min_values_per_group = 2L, min_groups = 1L) {
    list(operation = paste0(binding$prefix, "_intensity_filter"),
        parameters = stats::setNames(list(percentile, proportion, min_values_per_group, min_groups), binding$intensity_parameters),
        state_name = state_name %||% paste0(binding$prefix, "_intensity_filtered"),
        description = description %||% sprintf(paste("Applied %s intensity filter (percentile: %d%%, proportion: %.2f,",
            "values required per group: %s, groups required: %s)"), binding$feature_noun, as.integer(percentile), proportion,
            omicsIntensityFilterShown(min_values_per_group), omicsIntensityFilterShown(min_groups)))
}

#' @noRd
omicsDuplicateResolutionStep <- function(binding, description = NULL, state_name = NULL) {
    list(operation = paste0(binding$prefix, "_duplicate_resolution"), parameters = list(),
        state_name = state_name %||% paste0(binding$prefix, "_duplicates_resolved"),
        description = description %||% sprintf("Resolved duplicate %s features by keeping highest intensity", binding$feature_noun))
}

#' @noRd
omicsItsdNormalisationStep <- function(binding, aggregation, feature_ids, description = NULL, state_name = NULL) {
    parameters <- list(itsd_aggregation = aggregation)
    parameters["itsd_feature_ids"] <- list(feature_ids)
    list(operation = paste0(binding$prefix, "_itsd_normalisation"), parameters = parameters,
        state_name = state_name %||% paste0(binding$prefix, "_itsd_norm"),
        description = description %||% paste("ITSD normalization (aggregation:", aggregation, ")"))
}

#' @noRd
omicsLog2Step <- function(binding, offset, description = NULL, state_name = NULL) {
    list(operation = paste0(binding$prefix, "_log2"), parameters = list(offset = offset),
        state_name = state_name %||% paste0(binding$prefix, "_log2"),
        description = description %||% paste("Log2 transformation (offset:", offset, ")"))
}

#' @noRd
omicsBetweenSampleNormalisationStep <- function(binding, method, description = NULL, state_name = NULL) {
    list(operation = paste0(binding$prefix, "_between_sample_normalisation"), parameters = list(normalisation_method = method),
        state_name = state_name %||% paste0(binding$prefix, "_normalized"),
        description = description %||% paste("Between-sample normalization (method:", method, ")"))
}

#' @param parameters The nine RUV parameters by their declared names (the tab's `ruv_params` plus the mode and grouping).
#' @noRd
omicsRuvStep <- function(binding, parameters, description = NULL, state_name = NULL) {
    # KIT-129 integration item 1: with the gate on, a request that names a new design keeps the new design's parameters.
    declared <- omicGuiRuvDesignEnabled() && !omicsRuvAbsent(parameters[["ruv_design", exact = TRUE]]) &&
        !identical(parameters$ruv_design, "legacy")
    list(operation = paste0(binding$prefix, "_ruv"), parameters = if (declared) stats::setNames(lapply(omicsRuvDesignParameterNames,
        \(name) parameters[[name, exact = TRUE]]), omicsRuvDesignParameterNames) else parameters[omicsRuvParameterNames],
        state_name = state_name %||% paste0(binding$prefix, "_ruv_corrected"),
        description = description %||% "RUV-III batch correction complete")
}

#' @noRd
omicsPearsonSelectionStep <- function(binding, pairs, threshold, description = NULL, state_name = NULL) {
    list(operation = paste0(binding$prefix, "_pearson_selection"),
        parameters = list(pairs = pairs, min_pearson_correlation_threshold = threshold),
        state_name = state_name %||% paste0(binding$prefix, "_correlation_filtered"),
        description = description %||% paste("Correlation filtering (threshold:", threshold, ")"))
}

#' Run one tab step through its wrapper's commit.
#' @noRd
runWorkflowOmicsReferenceStep <- function(manager, step, config, source = NULL, dependencies = list(), hydrate = FALSE) {
    entry <- workflowOperationEntry(step$operation)
    parameters <- step$parameters
    # KIT-089: a step whose default is not `whole` runs that variant in process too, as in the worker, so the record names
    # what ran. Other steps get no variant parameter.
    if (is.null(parameters$variant) && !identical(entry$default_variant %||% "whole", "whole")) parameters$variant <- entry$default_variant
    # KIT-129 integration item 1: the new design's declaration when the request names it, as on the worker.
    chosen <- workflowRequestOperation(manager, entry, source, parameters)
    runWorkflowOmicsReferenceOperation(manager, chosen$declaration, chosen$parameters, step$state_name,
        config, step$description, source = source, dependencies = dependencies, hydrate = hydrate)
}

#' @noRd
runWorkflowOmicsReferenceIntensityFilter <- function(binding, manager, percentile, proportion, config, description = NULL,
    state_name = NULL, hydrate = FALSE, min_values_per_group = 2L, min_groups = 1L) {
    runWorkflowOmicsReferenceStep(manager, omicsIntensityFilterStep(binding, percentile, proportion, description, state_name,
        min_values_per_group, min_groups), config, hydrate = hydrate)
}

#' @noRd
runWorkflowOmicsReferenceDuplicateResolution <- function(binding, manager, config, description = NULL, state_name = NULL,
    hydrate = FALSE) {
    runWorkflowOmicsReferenceStep(manager, omicsDuplicateResolutionStep(binding, description, state_name), config, hydrate = hydrate)
}

#' @noRd
runWorkflowOmicsReferenceItsdNormalisation <- function(binding, manager, aggregation, feature_ids, config, description = NULL,
    state_name = NULL, hydrate = FALSE) {
    runWorkflowOmicsReferenceStep(manager, omicsItsdNormalisationStep(binding, aggregation, feature_ids, description, state_name), config, hydrate = hydrate)
}

#' @noRd
runWorkflowOmicsReferenceLog2 <- function(binding, manager, offset, config, description = NULL, state_name = NULL, hydrate = FALSE) {
    runWorkflowOmicsReferenceStep(manager, omicsLog2Step(binding, offset, description, state_name), config, hydrate = hydrate)
}

#' @noRd
runWorkflowOmicsReferenceBetweenSampleNormalisation <- function(binding, manager, method, config, description = NULL,
    state_name = NULL, hydrate = FALSE) {
    runWorkflowOmicsReferenceStep(manager, omicsBetweenSampleNormalisationStep(binding, method, description, state_name), config, hydrate = hydrate)
}

#' @param parameters The nine RUV parameters by their declared names (the tab's `ruv_params` plus the mode and grouping).
#' @return The operation result with `optimisation`, the hydrated per-assay optimisation results.
#' @noRd
runWorkflowOmicsReferenceRuv <- function(binding, manager, parameters, config, description = NULL, state_name = NULL,
    hydrate = FALSE) {
    result <- runWorkflowOmicsReferenceStep(manager, omicsRuvStep(binding, parameters, description, state_name), config, hydrate = hydrate)
    # A new design's product has no sweep (KIT-129 integration item 1).
    if (!is.null(result$value$optimisation)) result$optimisation <- hydrateSavedWorkflowEntity(manager$.context$store, result$value$optimisation)$value
    result
}

#' The correlation filter: the pairs product, then the selection on it (the DIA two-operation pattern).
#' @return The selection's result with `correlations`, the per-assay correlation tables.
#' @noRd
runWorkflowOmicsReferenceCorrelation <- function(binding, manager, group, threshold, config, description = NULL,
    state_name = NULL, hydrate = FALSE) {
    if (!inherits(manager, "DIAStateAdapter")) {
        cli::cli_abort(paste(binding$step_scope, "needs the session's state owner"), class = "saved_entity_scope")
    }
    withWorkflowStateScope(manager, function() {
        source <- getWorkflowStateReference(manager)
        pairs <- runSavedOperation(manager, omicsPearsonPairsOperation(binding), source, list(correlation_group = group))
        # KIT-078: the selection reads the pairs through the head (a hit from a store whose head no longer roots it).
        if (is.null(manager$.context$snapshot$products[[pairs$reference$entity_id]])) rootWorkflowProduct(manager, pairs$reference, list(source))
        result <- runWorkflowOmicsReferenceStep(manager, omicsPearsonSelectionStep(binding, pairs$reference$entity_id, threshold, description, state_name),
            config, source = source, dependencies = list(pairs$reference), hydrate = hydrate)
        result$correlations <- pairs$value
        result
    })
}

#' Name a store step's provenance before the module's own `saveState` (the reset): the commit is attributed to the
#' declaration instead of `legacy`. A non-owner manager is refused (KIT-072).
#' @noRd
markWorkflowOmicsStoreStep <- function(manager, operation) {
    if (!inherits(manager, "DIAStateAdapter")) cli::cli_abort("the step needs the session's state owner", class = "saved_entity_scope")
    entry <- workflowOperationEntry(operation)
    setWorkflowStepProvenance(manager, operation, entry$default_variant %||% "whole", entry$parameters, entry$declaration$fixed)
    invisible(NULL)
}

#' Commit a store step on the owner (the QC finalise, the correlation skip, the reset): the head's record, or the
#' `source` state's, under the new name with the step's provenance, the prior products kept and no component re-serialised.
#' @noRd
commitWorkflowOmicsStoreStep <- function(manager, operation, state_name, config, description, source = NULL) {
    if (!inherits(manager, "DIAStateAdapter")) cli::cli_abort("the step needs the session's state owner", class = "saved_entity_scope")
    withWorkflowStateScope(manager, function() {
        markWorkflowOmicsStoreStep(manager, operation)
        # KIT-075: `source` (a saved state reference) recommits that state's record, as the reset does; the head's otherwise.
        commitWorkflowReferenceState(manager, state_name, readWorkflowStateMetadata(manager, source)$data, config, description)
        getWorkflowStateReference(manager)
    })
}

#' The design step on the owner (KIT-031/032): the `design_builder` provenance with the three design roles, and the
#' workflow's inputs (the import's source copies under `import_sources/` as `<input prefix><n>`, the design file).
#' @noRd
recordWorkflowOmicsDesignProvenance <- function(binding, manager, object, source_dir = NULL, design_file = NULL) {
    if (!inherits(manager, "DIAStateAdapter")) cli::cli_abort("the design step needs the session's state owner", class = "saved_entity_scope")
    setWorkflowStepProvenance(manager, "design_builder",
        parameter_names = c("sample_id", "group_id", "technical_replicate_id"),
        values = list(sample_id = object@sample_id, group_id = object@group_id,
            technical_replicate_id = object@technical_replicate_id))
    paths <- character()
    if (!is.null(source_dir) && dir.exists(file.path(source_dir, "import_sources"))) {
        sources <- sort(list.files(file.path(source_dir, "import_sources"), full.names = TRUE))
        if (length(sources)) paths[paste0(binding$design_input_prefix, seq_along(sources))] <- sources
    }
    if (is.null(design_file) && !is.null(source_dir)) design_file <- file.path(source_dir, "design_matrix.tab")
    if (!is.null(design_file) && file.exists(design_file)) paths[["design_file"]] <- design_file
    if (length(paths)) recordWorkflowInputs(manager, paths)
    invisible(NULL)
}

#' Run the omic's DA on a state, one declaration per contrast on the saved per-assay fits, and accept the manifest on
#' the head (the `runDIAReferenceDaAnalysis` twin; KIT-076, `<prefix>.da_results.v2`, keyed by contrast). A failed contrast
#' is recorded and the others complete.
#' @param manager The session owner.
#' @param source Saved reference of the state the DA runs on.
#' @param settings `formula_string`, `da_q_val_thresh`, `treat_lfc_cutoff`.
#' @param contrasts The contrasts table (`contrasts`, `friendly_names`, `full_format`).
#' @return `list(manifest, results, reference)`: `results` is the kernel's list for the module's display and files
#'   (`readWorkflowOmicsDaResults`; without the source object, whose readers take its components per action).
#' @noRd
runWorkflowOmicsReferenceDaAnalysis <- function(binding, manager, source, settings, contrasts) {
    if (!inherits(manager, "DIAStateAdapter")) {
        cli::cli_abort(paste(binding$da_scope, "needs the session's state owner"), class = "saved_entity_scope")
    }
    if (!isSavedWorkflowReference(source)) {
        cli::cli_abort(paste(binding$da_scope, "needs a saved state reference"), class = "saved_entity_scope")
    }
    withWorkflowStateScope(manager, function() {
        context <- manager$.context
        steps <- omicsDaContrastSteps(binding, settings, contrasts)
        results <- list()
        for (step in steps) {
            # KIT-129 integration item 1: a source of a new RUV design is tested by the DA declaration naming its method.
            chosen <- workflowRequestOperation(manager, workflowOperationEntry(step$operation), source, step$parameters)
            run <- runSavedOperation(manager, chosen$declaration, source, chosen$parameters)
            if (identical(run$status, "computed")) {
                snapshot <- context$snapshot
                snapshot$da_dependencies[[run$reference$entity_id]] <- collectWorkflowNativeReferences(run$value)
                publishWorkflowStateHead(context, snapshot)
            }
            results[[step$contrast]] <- run$value
        }
        acceptWorkflowOmicsDaManifest(binding, manager, source, settings, contrasts, results)
    })
}

#' One contrast's declaration parameters per row of the contrasts table (KIT-078: in process or a managed action each).
#' @return A list of `list(operation, contrast, parameters)`.
#' @noRd
omicsDaContrastSteps <- function(binding, settings, contrasts) {
    names_of <- omicsDaContrastNames(as.data.frame(contrasts, stringsAsFactors = FALSE))
    lapply(seq_along(names_of$strings), function(i) list(operation = paste0(binding$prefix, "_da_contrast"), contrast = names_of$strings[[i]],
        parameters = list(formula = settings$formula_string, qvalue = settings$da_q_val_thresh, lfc = settings$treat_lfc_cutoff,
            contrast = names_of$strings[[i]], friendly_name = names_of$friendly[[i]])))
}

#' Accept the per-contrast results' manifest on the head (KIT-078: the last step of the DA in process and under a managed
#' session), the DA step's provenance cleared (it commits no state). When no contrast completed the run stops with their
#' messages and the head keeps its previous result.
#' @param results The contrast operations' values, by contrast.
#' @return `list(manifest, results, reference)`, as `runWorkflowOmicsReferenceDaAnalysis`.
#' @noRd
acceptWorkflowOmicsDaManifest <- function(binding, manager, source, settings, contrasts, results, read_results = TRUE) {
    context <- manager$.context
    metadata <- readWorkflowOmicsMetadata(binding, manager, source)
    contrasts <- as.data.frame(contrasts, stringsAsFactors = FALSE)
    parameters <- list(formula = settings$formula_string, qvalue = settings$da_q_val_thresh, lfc = settings$treat_lfc_cutoff)
    clearWorkflowStepProvenance(manager)
    if (length(results) && all(vapply(results, function(entry) !length(entry$assays) && length(entry$errors) > 0L, logical(1)))) {
        stop(paste(unique(unlist(lapply(results, `[[`, "errors"))), collapse = "; "), call. = FALSE)
    }
    # KIT-095: the model record is the binding's (`da_model`); lipid and metab keep the limma one.
    model <- if (is.function(binding$da_model)) binding$da_model(parameters) else list(model = "limma",
        formula_used = parameters$formula, eBayes_trend = TRUE, eBayes_robust = TRUE, treat_applied = isTRUE(parameters$lfc > 0))
    model <- c(model, omicsDaBlockRecord(lapply(results, `[[`, "model")))
    manifest <- c(list(schema = paste0(binding$prefix, ".da_results.v2"), source = source),
        stats::setNames(list(metadata$slots[[binding$id_slot]]), binding$id_slot),
        list(sample_id = metadata$slots$sample_id, assays = metadata$assay_names, contrasts = contrasts, results = results,
            parameters = parameters, model = model))
    reference <- saveWorkflowStateValue(context, manifest, "da-result-references", codec = "rds")
    snapshot <- context$snapshot
    snapshot$da_results <- reference
    snapshot$da_dependencies[[reference$entity_id]] <- collectWorkflowNativeReferences(list(source, results))
    publishWorkflowStateHead(context, snapshot)
    list(manifest = manifest, results = if (read_results) readWorkflowOmicsDaResults(binding, manager, manifest), reference = reference)
}

# ---------------------------------------------------------------------------------------------------------------
# The worker (KIT-078): each action runs under `computeWorkflowOmicsAction`, in the owner's scope, and leaves the head
# to the controller's completion (`completeWorkflowOmicsAction`). Nothing here reads a Shiny session.
# ---------------------------------------------------------------------------------------------------------------

#' A directory the worker stages one role's files in; the controller publishes them (`publishWorkflowOmicsResult`). On the
#' worker it lies in the short folder `computeWorkflowOmicsAction` copies into `<request>/p` (KIT-052).
#' @noRd
omicsStagedPublicationDir <- function(request, role) {
    path <- file.path(request$publication_dir %||% file.path(request$directory, "omics-publication"), role)
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    path
}

#' The filtering step's metrics, saved for the controller's figures (the metabolite tables can pass the result's limit).
#' A failure costs the figures only, as it does in process.
#' @return The saved value's reference, or NULL.
#' @noRd
omicsFilteringMetricsReference <- function(binding, manager, state) {
    tryCatch({
        metrics <- binding$filtering_metrics(readWorkflowExecutionComponents(manager, state, fields = binding$fields))
        saveWorkflowStateValue(manager$.context, metrics, "omics-filtering-metrics", codec = "rds")
    }, error = function(error) {
        logger::log_warn(paste("Could not compute the filtering metrics:", conditionMessage(error)))
        NULL
    })
}

#' The worker's one compute handler: the registered declaration on the request's state, by the request's variant. A
#' state-producing declaration returns the checkpoint the completion commits; a product-only one returns the source's.
#' KIT-094: the binding names the parameters that are product dependencies, the steps whose product holds native
#' references and the steps whose figures need the filtering metrics.
#' @noRd
computeWorkflowOmicsOperationAction <- function(manager, request, variant) {
    entry <- workflowActionSpec(request$operation)
    binding <- omicsDeclarationBindings()[[entry$omic_type]]
    source <- request$input_references$state
    # KIT-129 integration item 1: the declaration the request runs, as in process (the new design's when it names it).
    chosen <- workflowRequestOperation(manager, entry, source, request$parameters)
    declaration <- chosen$declaration
    parameters <- chosen$parameters
    parameters$variant <- variant
    parameters$worker_cap <- request$worker_cap
    dependencies <- lapply(intersect(binding$dependency_parameters, declaration$parameters),
        function(parameter) manager$.context$snapshot$products[[parameters[[parameter]]]])
    run <- runSavedOperation(manager, declaration, source, parameters, dependencies)
    # An earlier worker's product whose action never completed has no root yet. KIT-111 F111-01: a hit roots what the
    # computed run rooted, the source, the dependencies and the saved references its value names; a value that is the
    # live science (a product-only declaration's) is never a root.
    if (identical(run$status, "hit")) {
        named <- tryCatch(collectWorkflowNativeReferences(run$value), saved_entity_invalid = function(error) list())
        rootWorkflowProduct(manager, run$reference, c(list(source), dependencies, named))
    }
    if (is.null(declaration$output_class)) {
        return(list(checkpoint = readWorkflowStateMetadata(manager, source)$data, reference = run$reference,
            native = if (declaration$name %in% binding$native_steps) collectWorkflowNativeReferences(run$value)))
    }
    final <- readWorkflowStateMetadata(manager, run$value$final)$data
    final$products[[paste0(declaration$name, "_science")]] <- run$reference
    filtering <- declaration$name %in% binding$filtering_steps
    list(checkpoint = final, reference = run$reference,
        display = if (filtering) list(filtering_metrics = omicsFilteringMetricsReference(binding, manager, run$value$final)))
}

#' One stage's per-assay QC figures from the request's state: the files staged, the assets saved.
#' @noRd
computeWorkflowOmicsQcFiguresAction <- function(manager, request, variant) {
    binding <- omicsDeclarationBindings()[[workflowActionSpec(request$operation)$omic_type]]
    source <- request$input_references$state
    binding$qc_plots(theObject = readWorkflowExecutionComponents(manager, source, fields = binding$fields),
        experiment_paths = stats::setNames(list(omicsStagedPublicationDir(request, binding$qc_dir_role)), binding$qc_dir_role),
        stage = request$parameters$stage, grouping_variable = request$parameters$grouping_variable,
        shape_variable = request$parameters$shape_variable)
    list(checkpoint = readWorkflowStateMetadata(manager, source)$data, reference = source)
}

#' The composite QC figure from the published per-assay files, staged and saved (its rasters are the application's largest
#' memory event when it is assembled there). A failure is reported to the tab, which warns as it does in process.
#' @noRd
computeWorkflowOmicsCompositeFigureAction <- function(manager, request, variant) {
    binding <- omicsDeclarationBindings()[[workflowActionSpec(request$operation)$omic_type]]
    source <- request$input_references$state
    composite <- tryCatch(binding$composite_figure(workflowRequestExperimentPaths(request)[[binding$qc_dir_role]],
        omicsStagedPublicationDir(request, binding$qc_dir_role), request$parameters$assay_names, request$parameters$ruv_mode),
        error = function(error) list(compositeSaved = FALSE, errorMessage = substr(conditionMessage(error), 1L, 1024L)))
    list(checkpoint = readWorkflowStateMetadata(manager, source)$data, reference = source,
        display = list(saved = isTRUE(composite$compositeSaved), message = composite$errorMessage))
}

#' The DA tables and figures from the accepted manifest and its source state: the files staged, the assets saved. A
#' failure is reported to the tab, which warns as it does in process; the science is already accepted.
#' @noRd
computeWorkflowOmicsDaPublicationAction <- function(manager, request, variant) {
    binding <- omicsDeclarationBindings()[[workflowActionSpec(request$operation)$omic_type]]
    source <- request$input_references$state
    manifest_reference <- manager$.context$snapshot$da_results
    if (is.null(manifest_reference)) cli::cli_abort("no accepted DA manifest", class = "workflow_execution_invalid")
    manifest <- resolveWorkflowStateValue(manager$.context, manifest_reference)
    results <- readWorkflowOmicsDaResults(binding, manager, manifest)
    results$theObject <- readWorkflowOmicsState(binding, manager, manifest$source)
    export <- tryCatch(list(success = isTRUE(binding$da_export(da_results_list = results,
        da_output_dir = omicsStagedPublicationDir(request, "da_output_dir"),
        publication_graphs_dir = omicsStagedPublicationDir(request, "publication_graphs_dir"),
        da_q_val_thresh = manifest$parameters$qvalue, lfc_threshold = manifest$parameters$lfc, heatmap_top_n = 50,
        heatmap_clustering = "both", heatmap_color_scheme = "RdBu"))), error = function(error) {
            list(success = FALSE, message = substr(conditionMessage(error), 1L, 1024L))
        })
    list(checkpoint = readWorkflowStateMetadata(manager, source)$data, reference = source, display = export)
}

#' The head's filtering metrics, for QC finalise's tracking figure (KIT-082).
#' @noRd
computeWorkflowOmicsFilteringMetricsAction <- function(manager, request, variant) {
    binding <- omicsDeclarationBindings()[[workflowActionSpec(request$operation)$omic_type]]
    source <- request$input_references$state
    list(checkpoint = readWorkflowStateMetadata(manager, source)$data, reference = source,
        display = list(filtering_metrics = omicsFilteringMetricsReference(binding, manager, source)))
}

#' The final-QC PCA of the request's state (KIT-082), saved as the drawn grobs: a plot's aesthetics keep the frames that built
#' it, which hold the state. A failure is the tab's to show, as it is in process.
#' @noRd
computeWorkflowOmicsDaTabResultsAction <- function(manager, request, variant) {
    binding <- omicsDeclarationBindings()[[workflowActionSpec(request$operation)$omic_type]]
    source <- request$input_references$state
    context <- manager$.context
    if (is.null(context$snapshot$da_results)) cli::cli_abort("no accepted DA manifest", class = "workflow_execution_invalid")
    manifest <- resolveWorkflowStateValue(context, context$snapshot$da_results)
    tab <- omicsDaTabResults(binding, readWorkflowOmicsDaResults(binding, manager, manifest))
    list(checkpoint = readWorkflowStateMetadata(manager, source)$data, reference = source,
        display = list(tab = saveWorkflowStateValue(context, tab, "omics-da-tab", codec = "rds")))
}

computeWorkflowOmicsFinalPcaAction <- function(manager, request, variant) {
    binding <- omicsDeclarationBindings()[[workflowActionSpec(request$operation)$omic_type]]
    source <- request$input_references$state
    parameters <- request$parameters
    display <- tryCatch(list(plot = saveWorkflowStateValue(manager$.context, omicsPlotGrob(binding$final_pca(
        readWorkflowOmicsState(binding, manager, parameters$state), parameters$grouping_variable, parameters$shape_variable)),
        "omics-final-pca", codec = "rds")), error = function(error) list(message = substr(conditionMessage(error), 1L, 1024L)))
    list(checkpoint = readWorkflowStateMetadata(manager, source)$data, reference = source, display = display)
}

#' A plot as its grobs, drawn at any size with `grid::grid.draw`.
#' @noRd
omicsPlotGrob <- function(plot) {
    # A bitmap device, as the tab's: the layout's text metrics are taken from the device the grobs are built on.
    path <- tempfile(fileext = ".png")
    grDevices::png(path)
    on.exit({
        grDevices::dev.off()
        unlink(path)
    }, add = TRUE)
    if (inherits(plot, "patchwork")) patchwork::patchworkGrob(plot) else ggplot2::ggplotGrob(plot)
}

#' One contrast and assay's interactive volcano as a self-contained page, saved as an asset and shown in an iframe (KIT-082,
#' the KIT-016 pattern). No page is `list(page = NULL)`: the widget reported nothing to draw.
#' @noRd
computeWorkflowOmicsVolcanoPageAction <- function(manager, request, variant) {
    binding <- omicsDeclarationBindings()[[workflowActionSpec(request$operation)$omic_type]]
    source <- request$input_references$state
    parameters <- request$parameters
    done <- function(display) list(checkpoint = readWorkflowStateMetadata(manager, source)$data, reference = source, display = display)
    context <- manager$.context
    manifest_reference <- context$snapshot$da_results
    if (is.null(manifest_reference)) cli::cli_abort("no accepted DA manifest", class = "workflow_execution_invalid")
    identity <- newSavedWorkflowScientificIdentity(manifest_reference$entity_id, "omics-volcano-page",
        digest::digest(list(parameters$contrast, parameters$assay, as.numeric(parameters$qvalue)), algo = "sha256"), "omics-volcano-page-v1", context$software)
    hit <- resolveSavedWorkflowEntity(context$store, "rendered_asset", identity, hydrate = FALSE)
    if (identical(hit$status, "hit")) {
        rootWorkflowProduct(manager, hit$reference, list(source))
        return(done(list(page = hit$reference)))
    }
    display <- tryCatch({
        manifest <- resolveWorkflowStateValue(context, manifest_reference)
        results <- readWorkflowOmicsDaResults(binding, manager, manifest)
        results$theObject <- readWorkflowOmicsState(binding, manager, manifest$source)
        widget <- binding$volcano_widget(da_results_list = results, selected_contrast = parameters$contrast,
            selected_assay = parameters$assay, da_q_val_thresh = parameters$qvalue)
        if (is.null(widget)) list(page = NULL) else {
            html <- tempfile("omics-volcano-", fileext = ".html")
            on.exit(unlink(html), add = TRUE)
            saveSelfContainedWidget(widget, html)
            reference <- saveWorkflowConsumerEntity(context$store, newSavedWorkflowRenderedAsset(readBin(html, "raw", n = file.size(html)), "text/html", c(1200L, 820L)),
                basename(tempfile("omics-volcano-html-")), identity, "rendered_asset")
            rootWorkflowProduct(manager, reference, list(source))
            list(page = reference)
        }
    }, error = function(error) list(message = substr(conditionMessage(error), 1L, 1024L)))
    done(display)
}

#' The session export's state, serialised here (KIT-082): the controller collected the rest of the session list and publishes
#' the staged file through the owner's existing retry path.
#' @noRd
computeWorkflowOmicsSessionExportAction <- function(manager, request, variant) {
    source <- request$input_references$state
    parameters <- request$parameters
    session_data <- resolveWorkflowStateValue(manager$.context, parameters$session)
    session_data$current_s4_object <- manager$getState(session_data$r6_current_state_name)
    stageWorkflowLegacyExport(session_data, workflowExportPaths(parameters$root, parameters$request_id), parameters$request_id, parameters$timestamp,
        serializeFn = saveRDS, serialize_only = TRUE, prefix = parameters$prefix)
    list(checkpoint = readWorkflowStateMetadata(manager, source)$data, reference = source)
}

#' The summary's integration object, written here (KIT-082). The file replaces any earlier one whole. KIT-095: the
#' binding's `integration_export(manager, parameters, path)` writes it when the binding has one (the transcriptomics
#' RUVIII matrix); lipid and metab write the state's object.
#' @noRd
computeWorkflowOmicsIntegrationExportAction <- function(manager, request, variant) {
    binding <- omicsDeclarationBindings()[[workflowActionSpec(request$operation)$omic_type]]
    source <- request$input_references$state
    parameters <- request$parameters
    temporary <- tempfile("omics-integration-", tmpdir = dirname(parameters$path), fileext = ".RDS")
    on.exit(unlink(temporary), add = TRUE)
    if (is.function(binding$integration_export)) binding$integration_export(manager, parameters, temporary) else
        saveRDS(manager$getState(parameters$state_name), temporary)
    if (!file.rename(temporary, parameters$path)) cli::cli_abort("cannot write the integration object", class = "saved_entity_write")
    # KIT-135: the declared integration input (KIT-061's contract) beside today's file, at `integration/<omic>/`.
    if (is.function(binding$integration_input)) binding$integration_input(binding, manager, file.path(dirname(parameters$path), binding$omic_type))
    list(checkpoint = readWorkflowStateMetadata(manager, source)$data, reference = source)
}

#' The lipidomics and metabolomics Summary's in-process integration writer (a route writer of the shared Summary core,
#' `startWorkflowSummaryIntegrationExport`): today's file, the final object's RDS, unchanged, then the declared
#' integration input (KIT-135) beside it, at `integration/<omic>/`.
#' @param path Today's RDS path; its folder is the project's integration folder.
#' @param object The final object (today's RDS content).
#' @param sourceDir Unused (the route writer's argument).
#' @return TRUE: the save shows today's notice.
#' @noRd
omicsSummaryIntegrationWrite <- function(binding, workflowData, path, object = NULL, sourceDir = NULL, saveRDSFn = saveRDS) {
    saveRDSFn(object, path)
    object <- NULL
    writeOmicsIntegrationInput(binding, workflowData$state_manager, file.path(dirname(path), binding$omic_type))
    TRUE
}

#' The lipidomics Summary's route writer (`omicsSummaryIntegrationWrite`).
#' @noRd
lipidSummaryIntegrationWrite <- function(workflowData, path, object = NULL, sourceDir = NULL, saveRDSFn = saveRDS) {
    omicsSummaryIntegrationWrite(lipidomicsDeclarationBinding, workflowData, path, object = object, sourceDir = sourceDir,
        saveRDSFn = saveRDSFn)
}

#' The metabolomics Summary's route writer (`omicsSummaryIntegrationWrite`).
#' @noRd
metabSummaryIntegrationWrite <- function(workflowData, path, object = NULL, sourceDir = NULL, saveRDSFn = saveRDS) {
    omicsSummaryIntegrationWrite(metabolomicsDeclarationBinding, workflowData, path, object = object, sourceDir = sourceDir,
        saveRDSFn = saveRDSFn)
}

#' Write the lipidomics or metabolomics integration input (KIT-135; KIT-061's contract): one view per assay, named
#' `<omic>_<assay>` (INT-4), from the head's final normalised state; on a session with no such state, or a manager that
#' is not an owner, none is written and an earlier one is removed (K061-3).
#' @param binding The omic's binding.
#' @param manager The session's owner.
#' @param dir The omic's folder, `<project>/integration/<omic>/`.
#' @return TRUE when an input was written, else FALSE, invisibly.
#' @noRd
writeOmicsIntegrationInput <- function(binding, manager, dir) {
    source <- if (inherits(manager, "DIAStateAdapter")) withWorkflowStateScope(manager, function() omicsIntegrationInputSource(binding, manager))
    if (is.null(source)) {
        removeIntegrationInput(dir)
        return(invisible(FALSE))
    }
    writeIntegrationInput(dir, source$manifest, source$samples, source$matrices, source$features, corrected = source$corrected)
    invisible(TRUE)
}

#' The lipidomics or metabolomics integration input (KIT-135; KIT-060 audit 4.3, row C4): the final state nearest the
#' head in the head's own lineage (K061-1), the first of `<prefix>_correlation_filtered`, `<prefix>_norm_complete`,
#' `<prefix>_ruv_corrected` and `<prefix>_normalized`, never a name kept in the history from an earlier run. Its assays
#' are the matrices: the RUV step's corrected values after RUV (decision M10: under a new design the state's table is
#' the corrected matrix, for figures, export and integration, while the DA tests the uncorrected values with W1..Wk in
#' the model), else the normalised values. The role is `corrected` when the lineage's RUV state applied RUV to at least
#' one assay (k above 0), else `normalised`; `correction` names the RUV design and, per view, whether it was applied.
#' Internal standards are kept and marked (INT-3: dropping them is the model's step).
#' @param binding The omic's binding.
#' @param manager The session's owner, inside its state scope.
#' @return `list(manifest, samples, matrices, features)`, or NULL when the lineage has no final state or it is stored
#'   as a whole legacy object.
#' @noRd
omicsIntegrationInputSource <- function(binding, manager) {
    prefix <- binding$prefix
    ruvState <- paste0(prefix, "_ruv_corrected")
    normalised <- paste0(prefix, "_normalized")
    # The whole lineage (`protIntegrationInputLineage`, K061-1): the RUV or normalisation state decides the role, and
    # the internal-standard step's selection marks the standards.
    lineage <- protIntegrationInputLineage(manager, until = paste0(prefix, c("_itsd_norm", "_qc_complete")))
    final <- utils::head(intersect(names(lineage), paste0(prefix, c("_correlation_filtered", "_norm_complete", "_ruv_corrected",
        "_normalized"))), 1L)
    if (!length(final) || is.null(lineage[[final]])) return(NULL)
    decider <- utils::head(intersect(names(lineage), c(ruvState, normalised)), 1L)
    reference <- getWorkflowStateReference(manager, final)
    object <- readWorkflowOmicsState(binding, manager, reference)
    assays <- methods::slot(object, binding$data_slot)
    id_column <- methods::slot(object, binding$id_slot)
    annotation_column <- object@annotation_id_column
    design <- as.data.frame(object@design_matrix)
    sample_column <- object@sample_id
    group_column <- object@group_id
    technical <- object@technical_replicate_id
    # The samples every assay holds, in the design's order.
    held <- Reduce(intersect, lapply(assays, names), as.character(design[[sample_column]]))
    samples <- design[match(held, as.character(design[[sample_column]])), , drop = FALSE]
    rownames(samples) <- NULL
    views <- stats::setNames(names(assays), vapply(names(assays), \(assay) integrationInputViewName(binding$omic_type, assay), ""))
    if (anyDuplicated(names(views))) abortIntegrationInput("two assays give the same view name")
    standards <- omicsIntegrationInputStandards(binding, manager, lineage, object)
    matrices <- list()
    features <- list()
    for (view in names(views)) {
        table <- as.data.frame(assays[[views[[view]]]], check.names = FALSE, stringsAsFactors = FALSE)
        ids <- as.character(table[[id_column]])
        matrices[[view]] <- data.frame(feature_id = ids, table[held], check.names = FALSE, stringsAsFactors = FALSE)
        annotation <- table[setdiff(names(table), c(as.character(design[[sample_column]]), "feature_id", "display_name", "assay",
            "is_internal_standard"))]
        display <- if (annotation_column %in% names(table)) as.character(table[[annotation_column]]) else ids
        features[[view]] <- data.frame(feature_id = ids, display_name = ifelse(is.na(display) | !nzchar(display), ids, display),
            assay = views[[view]], is_internal_standard = standards$mark(views[[view]], table, ids), annotation,
            check.names = FALSE, stringsAsFactors = FALSE)
    }
    object <- assays <- NULL
    correction <- if (identical(decider, ruvState)) omicsIntegrationInputCorrection(binding, manager, ruvState, lineage[[ruvState]],
        views, samples)
    # KIT-153 (INT-23): after RUV the final state's assays are the corrected matrices (of the views RUV was applied to);
    # each data set's own matrix is the normalised one the RUV step read, on the corrected matrix's features and samples.
    corrected <- list()
    normalisation_source <- list(name = final, entity_id = reference$entity_id, sha256 = reference$payload$sha256)
    if (!is.null(correction)) {
        input <- omicsIntegrationInputRuvInput(binding, manager, lineage, ruvState)
        input_assays <- methods::slot(input$object, binding$data_slot)
        input_id <- methods::slot(input$object, binding$id_slot)
        input$object <- NULL
        for (view in names(views)) {
            table <- as.data.frame(input_assays[[views[[view]]]], check.names = FALSE, stringsAsFactors = FALSE)
            uncorrected <- data.frame(feature_id = as.character(table[[input_id]]), table[intersect(held, names(table))],
                check.names = FALSE, stringsAsFactors = FALSE)
            if (isTRUE(correction$views[[view]]$applied)) corrected[[view]] <- matrices[[view]]
            matrices[[view]] <- integrationInputUncorrectedOn(uncorrected, matrices[[view]], view)
        }
        input_assays <- NULL
        normalisation_source <- list(name = input$name, entity_id = input$reference$entity_id,
            sha256 = input$reference$payload$sha256)
    }
    components <- readWorkflowStateMetadata(manager, reference)$data$components[[binding$data_slot]]
    workflow <- manager$workflow_type %||% manager$.context$snapshot$workflow_type %||% NA_character_
    if (!is.character(technical) || length(technical) != 1L || !isTRUE(technical %in% names(samples))) technical <- NA_character_
    manifest <- list(omic_type = binding$omic_type, workflow = workflow,
        matrix_role = "normalised", scale = "log2",
        normalisation = list(method = omicsIntegrationInputNormMethod(manager, lineage, normalised), source = normalisation_source),
        correction = correction,
        source_state = list(name = final, entity_id = reference$entity_id, sha256 = reference$payload$sha256,
            feature_id_column = id_column, annotation_id_column = annotation_column,
            assays = lapply(stats::setNames(views, names(views)), \(assay) list(assay = assay,
                entity_id = components[[assay]]$entity_id, sha256 = components[[assay]]$payload$sha256))),
        sample_id_column = sample_column, group_column = group_column, technical_replicate_column = technical,
        biological_sample_column = integrationInputBiologicalSampleColumn(samples, technical, sample_column),
        internal_standards = standards$rule)
    list(manifest = manifest, samples = samples, matrices = matrices, features = features, corrected = corrected)
}

#' The normalised object a lipidomics or metabolomics RUV step read (KIT-153, INT-23): the between-sample normalisation
#' state (`<prefix>_normalized`) the RUV state was committed on, in the head's lineage (the RUV step runs on the head the
#' normalisation leaves, under today's RUV-III-C and a new design alike). Refused in words when the RUV state does not
#' follow it.
#' @param lineage The head's lineage (`protIntegrationInputLineage`).
#' @param ruvState The RUV state's name.
#' @return `list(name, reference, object)`.
#' @noRd
omicsIntegrationInputRuvInput <- function(binding, manager, lineage, ruvState) {
    normalised <- paste0(binding$prefix, "_normalized")
    if (!identical(lineage[[ruvState]]$previous_state, normalised) || is.null(lineage[[normalised]])) {
        abortIntegrationInput(sprintf(paste("the normalised matrix that entered RUV (%s) is not the state before %s in this",
            "session; run the normalisation again, then save"), normalised, ruvState))
    }
    reference <- getWorkflowStateReference(manager, normalised)
    list(name = normalised, reference = reference, object = readWorkflowOmicsState(binding, manager, reference))
}

#' The between-sample normalisation of a lipidomics or metabolomics input (KIT-153): the method its normalisation step
#' recorded (`normalisation_method` or `norm_method` of the `<prefix>_normalized` state's provenance or configuration),
#' else NA.
#' @noRd
omicsIntegrationInputNormMethod <- function(manager, lineage, normalised) {
    values <- manager$.context$snapshot$step_provenance[[normalised]]$values
    for (method in list(values$normalisation_method, values$norm_method, lineage[[normalised]]$config$normalisation_method,
        lineage[[normalised]]$config$norm_method)) {
        if (is.character(method) && length(method) == 1L && !is.na(method)) return(method)
    }
    NA_character_
}

#' Which features of a lipidomics or metabolomics input are internal standards (KIT-135; INT-3: the writer marks them,
#' the model drops them): a feature whose id or annotation matches the object's `internal_standard_regex`, or one the
#' lineage's internal-standard step was given by hand (`itsd_feature_ids`, per assay, as its provenance recorded it).
#' @param lineage The head's lineage (`protIntegrationInputLineage`).
#' @param object The final state's object.
#' @return `list(mark, rule)`: `mark(assay, table, ids)` the logical marks of one assay's rows; `rule` the manifest's
#'   record of the rule.
#' @noRd
omicsIntegrationInputStandards <- function(binding, manager, lineage, object) {
    regex <- object@internal_standard_regex
    if (!is.character(regex) || length(regex) != 1L || is.na(regex) || !nzchar(regex)) regex <- NA_character_
    step <- paste0(binding$prefix, "_itsd_norm")
    selected <- if (step %in% names(lineage)) manager$.context$snapshot$step_provenance[[step]]$values$itsd_feature_ids
    # A selection too large for the provenance is recorded by class and size only; it is not read back.
    if (!is.list(selected) || identical(names(selected), c("omitted", "bytes"))) selected <- NULL
    columns <- unique(c(methods::slot(object, binding$id_slot), object@annotation_id_column))
    mark <- function(assay, table, ids) {
        marked <- ids %in% as.character(unlist(selected[[assay]]))
        if (!is.na(regex)) {
            for (column in intersect(columns, names(table))) {
                marked <- marked | stringr::str_detect(as.character(table[[column]]), regex) %in% TRUE
            }
        }
        marked
    }
    list(mark = mark, rule = list(regex = regex, columns = I(columns), selected_by_hand = !is.null(selected),
        selected_count = length(unlist(selected))))
}

#' The RUV design behind a corrected lipidomics or metabolomics input (KIT-135; KIT-060 audit row C4), from the RUV
#' state's record as the analysis record reads it (`workflowRecordStepOutcomes`): today's RUV-III-C, or a new design's
#' RUV-III with every field of `RUV_DESIGN_RECORD_FIELDS` (as the proteomics writer names them, KIT-129 seventh pass
#' item 3), and per view whether RUV was applied (k above 0), its k and its control count.
#' @param state The RUV state's name.
#' @param record The RUV state's record.
#' @param views The view names' assays (named by view).
#' @param samples The input's samples (for the grouping variable's levels).
#' @return The manifest's `correction`, or NULL when RUV was applied to no assay (a no-RUV or no-factor state is not
#'   labelled corrected).
#' @noRd
omicsIntegrationInputCorrection <- function(binding, manager, state, record, views, samples) {
    outcomes <- workflowRecordStepOutcomes(manager$.context, state, paste0(binding$prefix, "_ruv"), record, list())
    values <- manager$.context$snapshot$step_provenance[[state]]$values
    declared <- !is.null(outcomes$ruv_design) && !is.na(outcomes$ruv_design)
    k <- if (declared) outcomes$ruv_k_per_assay else outcomes$k_per_assay
    controls <- if (declared) outcomes$ruv_control_count_per_assay else outcomes$negative_control_count_per_assay
    one <- function(value, cast) {
        value <- suppressWarnings(cast(unlist(value)))
        if (length(value) == 1L) value else cast(NA)
    }
    per_view <- lapply(views, function(assay) {
        k_assay <- one(k[[assay]], as.integer)
        applied <- isTRUE(k_assay > 0L) && (!declared || isTRUE(one(outcomes$ruv_applied_per_assay[[assay]], as.logical)))
        list(assay = assay, applied = applied, k = k_assay, control_count = one(controls[[assay]], as.integer),
            not_applied_reason = if (declared && !applied) one(outcomes$ruv_not_applied_reason_per_assay[[assay]], as.character) else NA)
    })
    if (!any(vapply(per_view, `[[`, logical(1), "applied"))) return(NULL)
    grouping <- if (declared) outcomes$ruv_grouping_variable else values$ruv_grouping_variable
    correction <- list(id = "C4", method = if (declared) "RUVIII" else "RUVIII_C_Varying",
        ruv_mode = if (declared) NA else values$ruv_mode,
        replicate_structure = if (declared) outcomes$ruv_replicate_structure else "group", grouping_variable = grouping,
        # Decision M10: the matrix is the RUV step's corrected values; under a new design the DA does not test them.
        matrix = "ruv_corrected_values",
        da_tests = if (declared) "uncorrected_values_with_factors_in_model" else "this_matrix",
        views = per_view)
    if (declared && all(RUV_DESIGN_RECORD_FIELDS %in% names(outcomes))) correction <- c(correction, outcomes[RUV_DESIGN_RECORD_FIELDS])
    correction <- lapply(correction, function(value) if (is.null(value)) NA else value)
    correction$grouping_levels <- I(if (is.character(grouping) && length(grouping) == 1L && grouping %in% names(samples))
        unique(as.character(samples[[grouping]])) else character())
    c(list(state = state), correction)
}

#' What the tab reads from an accepted omics action (KIT-078): the declaration's value (the product the worker saved) and
#' the worker's display handoff. No assay table.
#' @return `list(stateReference, reference, value, display)`.
#' @noRd
readManagedOmicsResult <- function(manager, result, workflow_data) {
    value <- if (!is.null(workflowActionSpec(result$operation)$declaration))
        hydrateSavedWorkflowEntity(manager$.context$store, result$reference)$value
    list(stateReference = result$state_reference, reference = result$reference, value = value, display = result$display)
}

#' Refuse to replay a recorded omics step.
#' @noRd
omicsReplayUnsupported <- function(step) {
    state <- step$state
    operation <- step$operation
    cli::cli_abort("step {.val {state}} cannot be replayed: operation {.val {operation}}", class = "workflow_replay_unsupported")
}

#' The replay of the lipidomics and metabolomics steps (a binding's `replay`, KIT-094), by the name after the prefix:
#' each is `function(binding, manager, step, parameters, config, description)`, the wrapper the module used.
#' @noRd
omicsRouteReplay <- local({
    store_step <- function(binding, manager, step, parameters, config, description) {
        commitWorkflowOmicsStoreStep(manager, step$operation, step$state, config, description)
    }
    list(
        intensity_filter = function(binding, manager, step, parameters, config, description) {
            # KIT-217: the two settings as recorded; a record made before has neither (NULL; its earlier rule reads neither).
            runWorkflowOmicsReferenceIntensityFilter(binding, manager, parameters[[binding$intensity_parameters[[1L]]]],
                parameters[[binding$intensity_parameters[[2L]]]], config, description, state_name = step$state, hydrate = FALSE,
                min_values_per_group = parameters[[binding$intensity_parameters[[3L]]]],
                min_groups = parameters[[binding$intensity_parameters[[4L]]]])
        },
        duplicate_resolution = function(binding, manager, step, parameters, config, description) {
            runWorkflowOmicsReferenceDuplicateResolution(binding, manager, config, description, state_name = step$state, hydrate = FALSE)
        },
        qc_finalize = store_step,
        norm_complete = store_step,
        itsd_normalisation = function(binding, manager, step, parameters, config, description) {
            runWorkflowOmicsReferenceItsdNormalisation(binding, manager, parameters$itsd_aggregation, parameters$itsd_feature_ids,
                config, description, state_name = step$state, hydrate = FALSE)
        },
        log2 = function(binding, manager, step, parameters, config, description) {
            runWorkflowOmicsReferenceLog2(binding, manager, parameters$offset, config, description, state_name = step$state, hydrate = FALSE)
        },
        between_sample_normalisation = function(binding, manager, step, parameters, config, description) {
            runWorkflowOmicsReferenceBetweenSampleNormalisation(binding, manager, parameters$normalisation_method, config, description,
                state_name = step$state, hydrate = FALSE)
        },
        ruv = function(binding, manager, step, parameters, config, description) {
            runWorkflowOmicsReferenceRuv(binding, manager, parameters, config, description, state_name = step$state, hydrate = FALSE)
        },
        pearson_selection = function(binding, manager, step, parameters, config, description) {
            pairs_operation <- paste0(binding$prefix, "_pearson_pairs")
            previous <- Filter(function(p) identical(p$operation, pairs_operation), step$preceded_by %||% list())
            if (!length(previous)) omicsReplayUnsupported(step)
            pairs <- workflowReplayParameters(previous[[1L]]$parameters, workflowOperationEntry(pairs_operation)$parameters)
            runWorkflowOmicsReferenceCorrelation(binding, manager, pairs$correlation_group, parameters$min_pearson_correlation_threshold,
                config, description, state_name = step$state, hydrate = FALSE)
        })
})

#' Replay one recorded omics step (KIT-031/032): the record's operation with its recorded values through the wrapper
#' the module used, committed under the recorded state name with the head's configuration and the recorded
#' description. KIT-094: the step's entry is the binding's `replay` entry for the name after its prefix.
#' @noRd
replayWorkflowOmicsStep <- function(binding, manager, step, parameters) {
    operation <- step$operation
    config <- readWorkflowStateMetadata(manager)$config
    description <- step$description %||% operation
    prefix <- paste0(binding$prefix, "_")
    if (!startsWith(operation, prefix)) omicsReplayUnsupported(step)
    suffix <- substring(operation, nchar(prefix) + 1L)
    if (!suffix %in% names(binding$replay)) omicsReplayUnsupported(step)
    binding$replay[[suffix]](binding, manager, step, parameters, config, description)
    invisible(NULL)
}
