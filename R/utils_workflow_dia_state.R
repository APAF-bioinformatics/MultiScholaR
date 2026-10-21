# DIA-025 experimental adapter. Ordinary WorkflowState remains unchanged.

#' Mark an internal state-publication boundary for deterministic failure checks.
#' @param phase Named boundary, never a scientific callback.
#' @return NULL invisibly.
#' @noRd
diaStateBoundary <- function(phase) invisible(NULL)

#' Observe an active ownership boundary without inspecting scientific payloads.
#' @param context Existing private DIA owner.
#' @param stage Fixed operation boundary label.
#' @return NULL invisibly; the optional observer belongs to this session only.
#' @noRd
observeDIAOwnerBoundary <- function(context, stage) {
    if (is.function(context$on_owner_boundary)) context$on_owner_boundary(stage)
    invisible(NULL)
}

#' Prepare one record using the proven per-section runtime portability policy.
#' @param record Complete legacy history record.
#' @return Portable record with unchanged scientific data and worker counts.
#' @noRd
prepareDIAStateRecord <- function(record) {
    if (!is.list(record) || is.null(names(record)) || anyDuplicated(names(record)) ||
        !"data" %in% names(record)) {
        cli::cli_abort("invalid DIA state record", class = "saved_entity_invalid")
    }
    record <- mapDIAStateConfig(record, function(config) {
        for (section in names(config)) {
            value <- config[[section]]
            if (!is.list(value) || !inherits(value$core_utilisation, "multidplyr_cluster")) next
            config[[section]] <- splitDIAWorkerConfig(stats::setNames(list(value), section))$config[[section]]
        }
        config
    })
    prepareDIAWorkflowHistory(list(state = record))[[1L]]
}

#' Save an immutable record through the existing checkpoint codec and limits.
#' @param context Private adapter owner.
#' @param value Complete record or small reference manifest.
#' @param label Logical identity, not a filesystem path.
#' @param codec Existing codec; scientific records use qs2.
#' @return Accepted immutable descriptor.
#' @noRd
saveDIAStateValue <- function(context, value, label, codec = "qs2") {
    checkWorkflowEntityStore(context$store)
    observeDIAOwnerBoundary(context, paste0("save_begin:", label))
    reference <- saveWorkflowCheckpoint(context$store, value,
        basename(tempfile("dia-state-")),
        list(input = context$store$identity$incarnation, design = "record-owned",
            parameters = label, selectors = "complete-record",
            software = as.character(utils::packageVersion("MultiScholaR"))),
        list(class = class(value), shape = as.integer(dim(value)), factors = list()), codec = codec)
    context$writes <- context$writes + 1L
    releaseDIACheckpointFilePages(context, reference)
    # A write's input is already owned by this action. Reuse it until delivery,
    # rather than decoding a second copy immediately after serializing it.
    if (context$scope_depth > 0L && identical(codec, "qs2")) {
        assign(reference$entity_id, value, context$working)
    }
    observeDIAOwnerBoundary(context, paste0("save_complete:", label))
    reference
}

#' Commit one metadata-only head after all its referenced records are durable.
#' @param context Private adapter owner.
#' @param snapshot References, topology and bounded metadata; never scientific values.
#' @return NULL invisibly; earlier durable heads remain recoverable on failure.
#' @noRd
publishDIAStateHead <- function(context, snapshot) {
    lease <- acquireWorkflowMutationLease(context$store)
    on.exit(if (!is.null(lease)) filelock::unlock(lease), add = TRUE)
    current <- workflowExecutionHead(context$store)
    if (!is.null(context$head) && !identical(current, context$head))
        cli::cli_abort("workflow head changed before mutation", class = "workflow_execution_stale")
    reference <- saveDIAStateValue(context, snapshot, "state-head", codec = "rds")
    temporary <- tempfile("dia-head-", tmpdir = context$store$directory)
    saveRDS(reference, temporary, compress = FALSE, version = 3)
    syncSavedWorkflowPaths(context$store, c(temporary, context$store$directory))
    prepareWorkflowHeadAcceptance(context, reference)
    diaStateBoundary("before_head_publication")
    checkWorkflowEntityStore(context$store)
    if (!file.rename(temporary, context$head_path)) {
        cli::cli_abort("DIA state head publication failed", class = "saved_entity_write")
    }
    syncSavedWorkflowPaths(context$store, c(context$head_path, context$store$directory))
    diaStateBoundary("after_head_publication")
    context$snapshot <- snapshot
    context$head <- reference
    finishWorkflowHeadAcceptance(context)
    if (is.function(context$on_commit)) context$on_commit()
    invisible(NULL)
}

#' Resolve one record, sharing only within an explicitly owned action scope.
#' @param context Private adapter owner.
#' @param reference Immutable descriptor.
#' @return Complete record; missing/corrupt payloads are errors.
#' @noRd
resolveDIAStateValue <- function(context, reference) {
    checkWorkflowEntityStore(context$store)
    key <- reference$entity_id
    if (context$scope_depth > 0L && exists(key, context$working, inherits = FALSE)) {
        return(get(key, context$working, inherits = FALSE))
    }
    value <- hydrateSavedWorkflowCheckpoint(context$store, reference)
    releaseDIACheckpointFilePages(context, reference)
    context$hydrates <- context$hydrates + 1L
    if (context$scope_depth > 0L) assign(key, value, context$working)
    value
}

#' Release pages of one consumed native-session checkpoint, never the file.
#' @param context Current owner with a durably saved or validated/decoded record.
#' @param reference That record's immutable descriptor.
#' @return NULL invisibly; ordinary/legacy routes retain their existing I/O.
#' @noRd
releaseDIACheckpointFilePages <- function(context, reference) {
    if (!identical(context$session_format, "native_v1")) return(invisible(NULL))
    releaseDIANativeFileCache(file.path(context$store$directory, reference$payload$relative_path),
        file.path(context$store$directory, "checkpoint-io.log"))
    invisible(NULL)
}

#' Read the current checkpoint descriptor without loading its scientific data.
#' @param manager Existing native or experimental state manager.
#' @return Current immutable descriptor or NULL, without hydration.
#' @noRd
getDIACurrentReference <- function(manager) {
    if (!inherits(manager, "DIAStateAdapter")) return(NULL)
    manager$.context$snapshot$states[[manager$current_state]]
}

#' Identify the experimental component checkpoint path.
#' @param manager Session-owned state adapter.
#' @return TRUE only after explicit native DIA activation.
#' @noRd
usesDIAReferenceStates <- function(manager) {
    inherits(manager, "DIAStateAdapter") && isTRUE(manager$.context$reference_states)
}

#' Obtain a checkpoint reference without loading its record or assay.
#' @param manager Session owner.
#' @param state_name Existing state name; NULL selects the current state.
#' @return Immutable descriptor or NULL.
#' @noRd
getDIAStateReference <- function(manager, state_name = NULL) {
    if (!inherits(manager, "DIAStateAdapter")) return(NULL)
    if (is.null(state_name)) state_name <- manager$current_state
    manager$.context$snapshot$states[[state_name]]
}

#' Recognize supported reference-backed checkpoint generations.
#' @param value Decoded bounded checkpoint specification.
#' @return TRUE for supported component metadata.
#' @noRd
isDIAComponentCheckpoint <- function(value) {
    is.list(value) && isTRUE(value[["schema", exact = TRUE]] %in%
        c("dia.component_checkpoint.v1", "dia.component_checkpoint.v2"))
}

#' Read only a component checkpoint's bounded record.
#' @param manager Session owner.
#' @param reference Immutable state descriptor.
#' @return Record containing component descriptors, never an assay.
#' @noRd
readDIAStateMetadata <- function(manager, reference = NULL) {
    record <- readComponentMetadata(manager, reference, classes = c("PeptideQuantitativeData", "ProteinQuantitativeData"))
    record
}

#' Read a component checkpoint's bounded record for any declared class (KIT-002).
#' @param manager Session owner.
#' @param reference Immutable state descriptor.
#' @param classes Optional class names the caller accepts; another class is `saved_entity_corrupt`.
#' @return Record containing component descriptors, never an assay.
#' @noRd
readComponentMetadata <- function(manager, reference = NULL, classes = NULL) {
    if (is.null(reference)) reference <- getDIAStateReference(manager)
    if (is.null(reference)) return(NULL)
    if (!identical(reference$payload$codec, "rds") || reference$payload$bytes > 8388608) {
        cli::cli_abort("state requires explicit legacy conversion", class = "saved_entity_unsupported")
    }
    record <- resolveDIAStateValue(manager$.context, reference)
    if (!isDIAComponentCheckpoint(record$data)) {
        cli::cli_abort("state has no component metadata", class = "saved_entity_corrupt")
    }
    spec <- record$data
    accepted <- if (is.null(classes)) names(workflowClassDeclarations()) else classes
    if (!is.character(spec$class) || length(spec$class) != 1L || !spec$class %in% accepted) {
        cli::cli_abort("component checkpoint is incomplete", class = "saved_entity_corrupt")
    }
    declaration <- workflowClassDeclaration(spec$class)
    required <- workflowClassRequiredComponents(declaration)
    complete <- function(field) {
        component <- spec$components[[field]]
        if (field %in% declaration$lists) return(is.list(component) && !isSavedWorkflowReference(component) &&
            all(vapply(component, isSavedWorkflowReference, logical(1))))
        isSavedWorkflowReference(component)
    }
    if (!all(required %in% names(spec$components)) || !all(vapply(required, complete, logical(1)))) {
        cli::cli_abort("component checkpoint is incomplete", class = "saved_entity_corrupt")
    }
    record
}

#' Assemble one explicitly requested compatibility table using bounded reads.
#' @param context Current owner.
#' @param reference Table descriptor.
#' @param max_bytes Retained v2 compatibility allocation, at most 256 MiB;
#'   one existing bounded read batch can overlap this allocation.
#' @return Original data frame; intermediate batches are released on each iteration.
#' @noRd
materializeDIATable <- function(context, reference, selected_columns = NULL,
    scope = NULL, relation = NULL, filter_column = NULL, filter_value = NULL,
    max_bytes = 268435456) {
    if (is.null(scope)) {
        scope <- newWorkflowTableReadScope(context$store, reference,
            reader = if (is.null(selected_columns)) "arrow" else "fastduckdb")
        on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
    }
    metadata <- scope$metadata
    predicate <- NULL
    if (!is.null(relation)) {
        assertWorkflowPhysicalTableScope(scope)
        quoted <- as.character(DBI::dbQuoteIdentifier(scope$connection, relation))
        column <- as.character(DBI::dbQuoteIdentifier(scope$connection, filter_column))
        filter_parameters <- as.list(filter_value[!is.na(filter_value)])
        predicate <- if (length(filter_parameters) == 1L) paste0(column, " = ?") else
            if (length(filter_parameters)) paste0(column, " IN (",
                paste(rep("?", length(filter_parameters)), collapse = ","), ")") else "FALSE"
        if (anyNA(filter_value)) predicate <- paste0("(", predicate, " OR ", column, " IS NULL)")
        metadata$rows <- as.integer(DBI::dbGetQuery(scope$connection,
            paste0("SELECT count(*) AS n FROM ", quoted, " WHERE ", predicate),
            params = filter_parameters)$n)
    }
    if (is.null(selected_columns)) selected_columns <- metadata$columns
    if (!is.character(selected_columns) || anyDuplicated(selected_columns) ||
        !all(selected_columns %in% metadata$columns)) {
        cli::cli_abort("invalid scientific table projection", class = "saved_entity_invalid")
    }
    if (!is.numeric(max_bytes) || length(max_bytes) != 1L || is.na(max_bytes) ||
        !is.finite(max_bytes) || max_bytes < 1 || max_bytes > 268435456)
        cli::cli_abort("invalid selected S4 materialization limit", class = "saved_entity_limit")
    fixed_bytes <- as.double(metadata$rows) * workflowProjectedRowBytes(metadata, selected_columns)
    if (!is.null(scope$view) && fixed_bytes + as.numeric(object.size(metadata)) > max_bytes)
        cli::cli_abort("selected S4 table exceeds its materialization limit",
            class = "saved_entity_limit")
    metadata$columns <- selected_columns
    metadata$column_metadata <- metadata$column_metadata[selected_columns]
    columns <- lapply(metadata$column_metadata, \(spec) {
        if (!is.null(spec$levels)) factor(rep(NA_character_, metadata$rows),
            levels = spec$levels, ordered = isTRUE(spec$ordered)) else vector(spec$typeof, metadata$rows)
    })
    row_names <- if (isTRUE(metadata$automatic_row_names)) NULL else character(metadata$rows)
    source_positions <- if (is.null(relation)) NULL else integer(metadata$rows)
    retained_bytes <- max(fixed_bytes, as.numeric(object.size(columns)) +
        as.numeric(object.size(row_names))) + as.numeric(object.size(metadata))
    if (!is.null(scope$view) && retained_bytes > max_bytes)
        cli::cli_abort("selected S4 table exceeds its materialization limit",
            class = "saved_entity_limit")
    # Charge strings per batch without assuming that strings shared between
    # batches remain interned. Vector slots and restoration metadata are
    # already reserved; this deliberately overestimates shared string storage.
    string_bytes <- function(value) max(0, as.numeric(object.size(value)) - 8 * length(value))
    if (!metadata$rows) {
        for (name in intersect(names(metadata$name_columns), names(columns))) names(columns[[name]]) <- character()
    }
    blocks <- split(metadata$columns, ceiling(seq_along(metadata$columns) / 64L))
    for (block in blocks) {
    query <- if (is.null(relation)) NULL else paste0(
        "SELECT \".multischolar_order\" AS \".source_order\", ",
        paste(DBI::dbQuoteIdentifier(scope$connection,
            workflowTableStorageColumns(scope$metadata, block)), collapse = ", "),
        " FROM ", quoted, " WHERE ", predicate, " ORDER BY \".multischolar_order\"")
    streamSavedWorkflowTable(scope, block, function(data, offset, source_names) {
        if (!nrow(data)) return(invisible(NULL))
        if (!is.null(scope$view)) {
            incoming_bytes <- sum(vapply(data, function(column) {
                if (is.character(column)) string_bytes(column) else
                    if (!is.null(names(column))) string_bytes(names(column)) else 0
            }, numeric(1)))
            if (!is.null(row_names) && identical(block, blocks[[1L]]))
                incoming_bytes <- incoming_bytes + string_bytes(source_names)
            if (retained_bytes + incoming_bytes > max_bytes)
                cli::cli_abort("selected S4 table exceeds its materialization limit",
                    class = "saved_entity_limit")
            retained_bytes <<- retained_bytes + incoming_bytes
        }
        at <- seq.int(offset + 1L, offset + nrow(data))
        for (name in names(data)) {
            columns[[name]][at] <<- data[[name]]
            if (name %in% names(metadata$name_columns)) names(columns[[name]])[at] <<- names(data[[name]])
        }
        if (!is.null(row_names)) row_names[at] <<- source_names
        if (!is.null(source_positions)) source_positions[at] <<- attr(data, "dia_source_order")
        invisible(NULL)
    }, query = query, params = if (is.null(relation)) NULL else filter_parameters,
        expected_rows = metadata$rows)
    if (!is.null(scope$view) && as.numeric(object.size(columns)) > max_bytes)
        cli::cli_abort("selected S4 table exceeds its materialization limit",
            class = "saved_entity_limit")
    }
    for (name in metadata$columns) {
        spec <- metadata$column_metadata[[name]]
        if (is.null(spec$levels) && length(spec$class)) class(columns[[name]]) <- spec$class
        for (attribute in names(spec$attributes)) attr(columns[[name]], attribute) <- spec$attributes[[attribute]]
        columns[[name]] <- restoreWorkflowAttributeOrder(columns[[name]], spec$attribute_order)
    }
    for (attribute in names(metadata$attributes)) attr(columns, attribute) <- metadata$attributes[[attribute]]
    if (!is.null(row_names) && identical(metadata$row_index_type, "integer")) {
        row_names <- as.integer(row_names)
        if (anyNA(row_names)) cli::cli_abort("invalid integer row names", class = "saved_entity_corrupt")
    }
    attr(columns, "row.names") <- row_names %||% .set_row_names(metadata$rows)
    class(columns) <- metadata$class
    columns <- restoreWorkflowAttributeOrder(columns, metadata$attribute_order)
    if (!is.null(scope$view) && as.numeric(object.size(columns)) > max_bytes)
        cli::cli_abort("selected S4 table exceeds its materialization limit",
            class = "saved_entity_limit")
    if (!is.null(source_positions)) attr(columns, "dia_source_order") <- source_positions
    columns
}

#' Open only declared components for an existing R scientific kernel.
#' @param manager Current session owner.
#' @param reference Immutable component checkpoint; defaults to the selected state.
#' @param fields Required scientific tables, never history or support ledgers.
#' @param columns Optional named per-table projections.
#' @return Transient S4 execution view carrying its immutable source descriptor.
#'   Slots not named in `fields` are zero-row prototypes with the stored schema,
#'   listed by the `dia_execution_fields` attribute. Callers may read their column
#'   schema but must request every table a kernel computes on; the checkpoint
#'   writer keeps the stored component for any slot that was not loaded.
#' @noRd
readDIAExecutionComponents <- function(manager, reference = NULL,
    fields = c("protein_quant_table", "design_matrix", "protein_id_table"), columns = list()) {
    if (is.null(reference)) reference <- getDIAStateReference(manager)
    record <- readComponentMetadata(manager, reference)
    spec <- record$data
    object <- rebuildDIACheckpointObject(manager$.context, spec, fields, columns)
    attr(object, "dia_component_source") <- reference
    attr(object, "dia_execution_fields") <- intersect(fields, names(spec$components))
    object
}

#' Rebuild a declared object from a component checkpoint, hydrating only the requested fields.
#' @param context Current owner.
#' @param spec Component checkpoint record (`data` of a state record).
#' @param fields Slot names to hydrate, or `"<slot>/<element>"` for one element of a list slot;
#'   `TRUE` hydrates everything. Unrequested tables become zero-row prototypes, unrequested list
#'   elements zero-row prototypes (tables) or NULL (values), unrequested nested objects are rebuilt
#'   with prototypes, unrequested opaque slots keep the class prototype.
#' @param columns Named list of column selections per table slot.
#' @return The S4 object; support entries are not restored here (see materializeDIAStateRecord).
#' @noRd
rebuildDIACheckpointObject <- function(context, spec, fields = TRUE, columns = list()) {
    declaration <- workflowClassDeclaration(spec$class)
    everything <- isTRUE(fields)
    wanted <- function(field, element = NULL) everything || field %in% fields ||
        (!is.null(element) && paste(field, element, sep = "/") %in% fields)
    prototype <- function(component) {
        scope <- newWorkflowTableReadScope(context$store, component, metadata_only = TRUE)
        tryCatch(readSavedWorkflowTablePrototype(scope), finally = releaseWorkflowTableReadScope(scope))
    }
    read_table <- function(component, field, element = NULL) {
        if (wanted(field, element)) materializeDIATable(context, component, columns[[field]]) else prototype(component)
    }
    slots <- spec$slots
    # Sample projection of a wide assay (KIT-010): the design is hydrated and restricted to the selected
    # sample columns whatever `fields` says, so the rebuilt object's sample sets agree.
    assay_selection <- if (!is.null(declaration$assay)) columns[[declaration$assay]] else NULL
    projected_wide <- identical(declaration$layout, "wide") && !is.null(assay_selection)
    for (field in names(spec$components)) {
        component <- spec$components[[field]]
        if (projected_wide && identical(field, declaration$design) && isWorkflowTableReference(component)) {
            # Only a selection that names sample columns restricts the design; a projection to identifier
            # columns alone keeps the design as `fields` requests it.
            design <- materializeDIATable(context, component)
            sample_slot <- slots[["sample_id"]]
            selected <- if (is.character(sample_slot) && length(sample_slot) == 1L && sample_slot %in% names(design))
                as.character(design[[sample_slot]]) %in% assay_selection else logical()
            if (any(selected)) {
                design <- design[selected, , drop = FALSE]
                rownames(design) <- NULL
                slots[[field]] <- design
                next
            }
            design <- NULL
        }
        if (field %in% declaration$lists) {
            elements <- lapply(names(component), function(name) {
                item <- component[[name]]
                if (isWorkflowTableReference(item)) read_table(item, field, name) else
                    if (wanted(field, name)) resolveDIAStateValue(context, item) else NULL
            })
            names(elements) <- names(component)
            slots[[field]] <- elements
        } else if (field %in% declaration$nested) {
            nested <- resolveDIAStateValue(context, component)$data
            # An unrequested nested object keeps prototypes for its assay but its (small) design
            # table is hydrated, so the nested class's own validity holds.
            nested_fields <- if (wanted(field)) TRUE else
                as.character(workflowClassDeclaration(nested$class)$design %||% character())
            slots[[field]] <- rebuildDIACheckpointObject(context, nested, nested_fields, list())
        } else if (field %in% declaration$opaque) {
            if (wanted(field)) slots[[field]] <- resolveDIAStateValue(context, component)
        } else if (isWorkflowTableReference(component)) {
            slots[[field]] <- read_table(component, field)
        } else if (wanted(field)) {
            slots[[field]] <- resolveDIAStateValue(context, component)
        }
    }
    # KIT-047: a declared owned args entry is reattached whenever the assay is hydrated (a kernel is about to run
    # or the whole object is wanted); a projected read keeps the reference, which the checkpoint writer keeps as is.
    for (section in declaration$args_as_owned %||% character()) {
        reference <- slots$args[[section]]
        if (!isSavedWorkflowReference(reference)) next
        if (!(everything || (!is.null(declaration$assay) && wanted(declaration$assay)))) next
        slots$args[[section]] <- resolveDIAStateValue(context, reference)
    }
    do.call(methods::new, c(list(Class = spec$class), slots))
}

#' Content digest of an owned args entry (KIT-047): a re-checkpoint of a hydrated object reuses the saved value
#' when the entry is unchanged, so the shared protein steps after the DPC quantification write it once.
#' @noRd
workflowOwnedArgsDigest <- function(value) digest::digest(value, algo = "xxhash64")

#' Store the declared owned args entries of an object out of line (KIT-047): each present entry is saved as one
#' owned value and replaced by its reference in `slots$args`; an entry equal to the source state's saved value
#' keeps that reference. The references are returned for the checkpoint's dependency roots.
#' @noRd
checkpointDIAOwnedArgs <- function(context, declaration, slots, previous) {
    owned <- list()
    for (section in declaration$args_as_owned %||% character()) {
        value <- slots$args[[section]]
        if (is.null(value)) next
        if (isSavedWorkflowReference(value)) { owned[[section]] <- value; next }
        before <- previous$slots$args[[section]]
        if (isSavedWorkflowReference(before) &&
            identical(workflowOwnedArgsDigest(resolveDIAStateValue(context, before)), workflowOwnedArgsDigest(value))) {
            owned[[section]] <- before
        } else {
            owned[[section]] <- saveDIAStateValue(context, value, paste0("owned_args:", section), codec = "rds")
        }
        slots$args[[section]] <- owned[[section]]
    }
    list(slots = slots, owned = owned)
}

#' Persist a scientific result as reusable components and bounded metadata.
#' @param manager Current owner.
#' @param object Already-owned scientific result, never obtained via compatibility hydration.
#' @return Component checkpoint specification for the existing state-head commit.
#' @noRd
checkpointDIAExecutionResult <- function(manager, object) {
    # An undeclared class is refused before any component is written.
    declaration <- workflowClassDeclaration(object)
    context <- manager$.context
    source <- attr(object, "dia_component_source", exact = TRUE)
    previous <- if (isSavedWorkflowReference(source)) readDIAStateMetadata(manager, source)$data else NULL
    slots <- stats::setNames(lapply(methods::slotNames(object),
        \(field) methods::slot(object, field)), methods::slotNames(object))
    owned <- list()
    if (length(declaration$args_as_owned) && is.list(slots$args)) {
        stored <- checkpointDIAOwnedArgs(context, declaration, slots, previous)
        slots <- stored$slots
        owned <- stored$owned
    }
    rememberDIAWorkers(context, slots$args)
    slots$args <- prepareDIAStateRecord(list(data = NULL, config = slots$args))$config
    components <- previous$components %||% list()
    components <- components[intersect(names(components), names(slots))]
    support <- previous$support %||% list()
    loaded <- attr(object, "dia_execution_fields", exact = TRUE)
    unchanged <- function(field) !is.null(previous) && field %in% names(previous$components) && !field %in% loaded
    zero_row <- function(value) is.data.frame(value) && !nrow(value)
    for (field in names(slots)) {
        value <- slots[[field]]
        if (field %in% declaration$lists) {
            if (!is.list(value) || is.data.frame(value)) next
            if (unchanged(field) && all(vapply(value, zero_row, logical(1)))) { slots[[field]] <- NULL; next }
            if (is.null(names(value)) && length(value)) cli::cli_abort("list slot {.field {field}} needs element names",
                class = "saved_entity_invalid")
            components[[field]] <- lapply(stats::setNames(names(value), names(value)), function(name) {
                element <- value[[name]]
                if (is.data.frame(element) && ncol(element) > 0L) saveDIATableComponent(context, element) else
                    saveDIAStateValue(context, element, paste0("element:", field, ":", name))
            })
            slots[[field]] <- NULL
            next
        }
        if (field %in% declaration$nested) {
            if (!isS4(value)) next
            if (unchanged(field)) { slots[[field]] <- NULL; next }
            nested <- checkpointDIAExecutionResult(manager, value)
            components[[field]] <- saveDIAStateValue(context, list(schema = "dia.nested_checkpoint.v1", data = nested),
                paste0("nested:", field), codec = "rds")
            slots[[field]] <- NULL
            next
        }
        if (field %in% declaration$opaque) {
            if (is.null(value) || !length(value)) next
            if (unchanged(field)) { slots[[field]] <- NULL; next }
            components[[field]] <- saveDIAStateValue(context, value, paste0("opaque:", field))
            slots[[field]] <- NULL
            next
        }
        table <- is.data.frame(value) && ncol(value) > 0L
        sizeable <- is.matrix(value) && length(value) > 0L
        if (!table && !sizeable) next
        if (unchanged(field) && !nrow(value)) {
            slots[[field]] <- NULL
            next
        }
        components[[field]] <- if (table) saveDIATableComponent(context, value) else
            saveDIAStateValue(context, value, paste0("component:", field))
        slots[[field]] <- NULL
    }
    for (section in names(slots$args)) {
        values <- slots$args[[section]]
        if (section %in% names(owned)) next
        if (is.data.frame(values) && ncol(values) > 0L) {
            support$.arguments[[section]] <- saveDIATableComponent(context, values)
            slots$args[section] <- list(NULL)
            next
        }
        if (!is.list(values) || is.data.frame(values)) next
        if (!identical(values, previous$slots$args[[section]])) support[[section]] <- NULL
        for (field in names(values)) {
            value <- values[[field]]
            if (!(is.data.frame(value) && ncol(value) > 0L) &&
                as.numeric(object.size(value)) < 65536) next
            support[[section]][[field]] <- if (is.data.frame(value) && ncol(value) > 0L)
                saveDIATableComponent(context, value) else
                saveDIAStateValue(context, value, paste0("support:", section, ":", field))
            slots$args[[section]][field] <- list(NULL)
        }
    }
    list(schema = "dia.component_checkpoint.v2", class = class(object)[[1L]],
        declaration = workflowDeclarationField(declaration),
        slots = slots, components = components, support = support,
        products = previous$products %||% list(),
        dependencies = collectDIANativeReferences(list(source, owned, previous$dependencies, previous$components,
            previous$support, previous$products)),
        dimensions = workflowAssayDimensions(object, declaration),
        metadata = if (is.function(declaration$describe)) declaration$describe(object) else NULL)
}

#' Resolve a protein scientific view while leaving public getState compatible.
#' @param manager Current owner.
#' @param state_name Optional state name.
#' @return Transient quantities, annotations and design required by protein kernels.
#' @noRd
readDIAProteinScience <- function(manager, state_name = NULL) {
    if (!usesDIAReferenceStates(manager)) return(if (is.null(state_name)) manager$getState() else manager$getState(state_name))
    assertDIAExportReady(manager)
    readDIAExecutionComponents(manager, getDIAStateReference(manager, state_name))
}

#' Carry an explicit scientific source across native kernels that construct a new S4.
#' @param object Newly calculated protein object.
#' @param manager,source Owning session and immutable input checkpoint.
#' @return The transient result with component provenance, never a hydrated predecessor.
#' @noRd
inheritDIAExecutionSource <- function(object, manager, source) {
    if (!usesDIAReferenceStates(manager)) return(object)
    attr(object, "dia_component_source") <- source
    attr(object, "dia_execution_fields") <- c("protein_quant_table", "design_matrix", "protein_id_table")
    object
}

#' Reconstruct a complete legacy record only at a compatibility boundary.
#' @param context Current owner.
#' @param reference State descriptor.
#' @return Original record and S4 data for public getState/history/export callers.
#' @noRd
materializeDIAStateRecord <- function(context, reference) {
    record <- resolveDIAStateValue(context, reference)
    if (isS4(record$data) || isDIAComponentCheckpoint(record$data)) {
        context$full_s4_reads <- context$full_s4_reads + 1L
    }
    if (!isDIAComponentCheckpoint(record$data)) return(record)
    spec <- record$data
    object <- rebuildDIACheckpointObject(context, spec, TRUE, list())
    if (length(spec$support) && "args" %in% methods::slotNames(object)) {
        args <- methods::slot(object, "args")
        for (section in names(spec$support)) {
            for (field in names(spec$support[[section]])) {
                component <- spec$support[[section]][[field]]
                value <- if (isWorkflowTableReference(component))
                    materializeDIATable(context, component) else resolveDIAStateValue(context, component)
                if (section == ".arguments") args[field] <- list(value) else
                    args[[section]][field] <- list(value)
            }
        }
        methods::slot(object, "args") <- args
    }
    record$data <- object
    record
}

#' Commit a component result using the same durable state head and topology.
#' @param manager Existing state owner.
#' @param state_name,checkpoint,config,description Accepted scientific result and provenance.
#' @param timestamp Original timestamp or the commit time.
#' @return Immutable state reference.
#' @noRd
#' Step provenance on the head snapshot (KIT-025; analysis-record audit AR-02/AR-03)
#'
#' The head records who produced each state and under what software. The two `state_history` writers
#' (`commitDIAReferenceState`, `saveDIAStateRecord`) take the pending provenance the running operation set
#' (`setDIAStepProvenance`: the worker completion from the request, `runSavedOperation` from its declaration,
#' the design builder and the import by name) and write `snapshot$step_provenance[[state_name]]`; a commit
#' with no pending provenance is `legacy`. Nothing is added to a state record, so every fingerprint oracle is
#' unchanged (the record-invariance rule). The software environment of a build is stored once per `build`.
#' @param operation Registry name, or `legacy`, `design_builder`, `import`.
#' @param variant Registry variant or NA.
#' @param parameter_names The registry entry's declared parameter names (values stay in the state's `slots$args`).
#' @param fixed Named list of the declaration's fixed parameter values, when it has any (the IQ rollup).
#' @noRd
#' @param values Named list of the declared parameters' values as the operation received them (KIT-027 audit
#'   PH-01: a step's config does not always record them under the declared names); bounded by
#'   `boundDIAProvenanceValues`.
setDIAStepProvenance <- function(manager, operation, variant = NA_character_, parameter_names = character(), fixed = NULL,
    values = NULL) {
    if (!inherits(manager, "DIAStateAdapter")) return(invisible(NULL))
    previous <- manager$.context$pending_provenance
    pending <- list(operation = operation, variant = variant, parameter_names = as.character(parameter_names), fixed = fixed,
        values = boundDIAProvenanceValues(values, parameter_names))
    # KIT-027 audit PH-02: a second operation before the same commit (correlation's pairs, then its selection) keeps
    # the first as `preceded_by`; the same operation set again (the DA contrasts) replaces it.
    if (!is.null(previous) && !identical(previous$operation, operation))
        pending$preceded_by <- c(previous$preceded_by, list(previous[c("operation", "variant", "parameter_names", "values")]))
    manager$.context$pending_provenance <- pending
    invisible(NULL)
}

#' Keep the declared parameters' values that are small plain values (atomic vectors and lists of them, at most
#' 64 KiB each); a table, a saved reference or a larger value is recorded by class and size only.
#' @noRd
boundDIAProvenanceValues <- function(values, parameter_names) {
    if (!is.list(values) || !length(parameter_names)) return(NULL)
    plain <- function(x) is.null(x) || is.atomic(x) && !is.data.frame(x) ||
        is.list(x) && !is.data.frame(x) && !isSavedWorkflowReference(x) && all(vapply(x, plain, logical(1)))
    out <- list()
    for (name in intersect(as.character(parameter_names), names(values))) {
        value <- values[[name]]
        out[[name]] <- if (plain(value) && as.numeric(object.size(value)) <= 65536) value else
            list(omitted = class(value)[[1L]], bytes = as.numeric(object.size(value)))
    }
    out
}

#' Drop a pending provenance that no commit will consume: a product-only operation (the DA contrasts) sets one
#' through `runSavedOperation` and commits no state (KIT-025 audit PH-01).
#' @noRd
clearDIAStepProvenance <- function(manager) {
    if (inherits(manager, "DIAStateAdapter")) manager$.context$pending_provenance <- NULL
    invisible(NULL)
}

#' Name a registered step's provenance from the section its compute recorded, unless an operation already set
#' it (the worker completion carries the request's variant): the in-process finishing steps (KIT-025 audit PH-02).
#' @noRd
setDIAStepProvenanceForSection <- function(manager, section) {
    if (!inherits(manager, "DIAStateAdapter") || !is.null(manager$.context$pending_provenance)) return(invisible(NULL))
    # KIT-047: a section shared by two routes (the accession cleanup) names the route's own declaration.
    operation <- workflowOperationForSection(section, manager$workflow_type)
    if (!length(operation)) return(invisible(NULL))
    entry <- workflowActionSpec(operation)
    setDIAStepProvenance(manager, operation, entry$default_variant %||% NA_character_, entry$parameters, entry$declaration$fixed)
}

#' @noRd
takeDIAStepProvenance <- function(context) {
    pending <- context$pending_provenance
    context$pending_provenance <- NULL
    if (is.null(pending)) list(operation = "legacy", variant = NA_character_, parameter_names = character(), fixed = NULL) else pending
}

#' @noRd
diaStepProvenanceEntry <- function(snapshot, provenance, timestamp = Sys.time()) {
    software <- workflowSoftwareIdentity()
    environment <- snapshot$software_environments[[software$build]]
    if (is.null(environment)) {
        snapshot$software_environments[[software$build]] <- software[intersect(c("version", "source_commit", "r_version",
            "platform", "packages", "blas", "lapack", "omp_threads", "blas_threads", "duckdb_threads", "os_release", "locale"),
            names(software))]
    } else if (length(setdiff(names(software$packages), names(environment$packages)))) {
        # KIT-047: packages loaded by a later step (the DPC engine) join the build's package list; nothing else changes.
        environment$packages <- software$packages
        snapshot$software_environments[[software$build]] <- environment
    }
    defaults <- system.file("config", "config.ini", package = "MultiScholaR")
    entry <- list(operation = provenance$operation, variant = provenance$variant,
        parameter_names = as.character(provenance$parameter_names), build = software$build,
        source_commit = software$source_commit,
        defaults_sha256 = if (nzchar(defaults)) digest::digest(file = defaults, algo = "sha256") else NA_character_,
        seed = NA_real_, committed_at = format(timestamp, "%Y-%m-%dT%H:%M:%OS3%z"))
    if (!is.null(provenance$fixed)) entry$fixed <- provenance$fixed
    if (length(provenance$values)) entry$values <- provenance$values
    if (length(provenance$preceded_by)) entry$preceded_by <- provenance$preceded_by
    if (isTRUE(provenance$backfilled)) entry$backfilled <- TRUE
    # KIT-042 audit F-01: an entry written at ingest for a state the adopting build did not run.
    if (isTRUE(provenance$adopted)) entry$adopted <- TRUE
    list(snapshot = snapshot, entry = entry)
}

#' @noRd
recordDIAStepProvenance <- function(context, snapshot, state_name, timestamp = Sys.time()) {
    pending <- context$pending_provenance
    # A recommit of an existing state with no operation pending (the DA step commits its arguments onto the head
    # before its contrasts, KIT-022) keeps the entry of the operation that produced the state.
    if (is.null(pending) && !is.null(snapshot$step_provenance[[state_name]])) return(snapshot)
    written <- diaStepProvenanceEntry(snapshot, takeDIAStepProvenance(context), timestamp)
    snapshot <- written$snapshot
    snapshot$step_provenance[[state_name]] <- written$entry
    snapshot
}

#' @noRd
truncateDIAStepProvenance <- function(snapshot) {
    keep <- unlist(snapshot$state_history)
    if (!is.null(snapshot$step_provenance)) snapshot$step_provenance <- snapshot$step_provenance[intersect(names(snapshot$step_provenance), keep)]
    snapshot
}

commitDIAReferenceState <- function(manager, state_name, checkpoint, config, description,
    timestamp = Sys.time()) {
    assertWorkflowController()
    checkSavedEntityString(state_name, "state_name")
    if (!usesDIAReferenceStates(manager) || !isDIAComponentCheckpoint(checkpoint)) {
        cli::cli_abort("invalid component checkpoint", class = "saved_entity_scope")
    }
    context <- manager$.context
    dependencies <- expandWorkflowTableReferences(context$store, collectDIANativeReferences(list(
        checkpoint$components, checkpoint$support, checkpoint$products, checkpoint$dependencies)))
    for (reference in dependencies) validateSavedWorkflowReference(context$store, reference, payload = FALSE)
    native <- context$native
    native$current_state <- context$snapshot$current_state
    native$state_history <- context$snapshot$state_history
    on.exit({ native$states <- native$audit_records <- list() }, add = TRUE)
    native$saveState(state_name, NULL, config, description)
    record <- native$states[[state_name]]
    record$timestamp <- timestamp
    record$data <- checkpoint
    record$s4_class <- checkpoint$class
    if (as.numeric(object.size(record)) > 8388608) {
        cli::cli_abort("component checkpoint metadata exceeds limits", class = "saved_entity_limit")
    }
    reference <- saveDIAStateValue(context, record, state_name, codec = "rds")
    snapshot <- context$snapshot
    for (field in intersect(c("iq_output_file", "peptide_input_file"), names(config))) {
        path <- config[[field]]
        if (is.character(path) && length(path) == 1L && !is.na(path)) {
            dependency <- c(list(role = sub("_file$", "", field), source = path,
                destination = paste("execution", field, basename(path), sep = "/")), describeDIANativeFile(path))
            snapshot$file_dependencies[[field]] <- dependency
            snapshot$file_dependencies[[paste0(field, "_", substr(reference$entity_id, 1L, 16L))]] <- dependency
        }
    }
    snapshot$states[[state_name]] <- reference
    snapshot$component_dependencies[[state_name]] <- dependencies
    snapshot$reference_states <- TRUE
    if (!is.null(checkpoint$metadata)) snapshot$state_metadata[[state_name]] <- checkpoint$metadata
    snapshot$current_state <- state_name
    snapshot$state_history <- native$state_history
    snapshot$sizes[[state_name]] <- list(rows = checkpoint$dimensions$rows,
        columns = length(checkpoint$dimensions$columns),
        note = "component dimensions; no assay hydration or R heap-size estimate")
    snapshot <- recordDIAStepProvenance(context, snapshot, state_name, timestamp)
    publishDIAStateHead(context, snapshot)
    if (is.function(context$store$progress)) context$store$progress(
        paste0("commit:", state_name), 1, 1, reference$entity_id, "checkpoint")
    reference
}

#' Hash the existing version-2 serialized value through a bounded Linux pipe.
#' @param value Exact value previously supplied to digest with serializeVersion=2.
#' @param store Existing entity owner for private scratch output.
#' @return The unchanged SHA256, without retaining a complete serialized raw vector.
#' @noRd
hashWorkflowSerializedValue <- function(value, store) {
    checkWorkflowEntityStore(store)
    commands <- Sys.which(c("bash", "dd", "sha256sum"))
    if (!identical(Sys.info()[["sysname"]], "Linux") || any(!nzchar(commands)))
        cli::cli_abort("bounded serialized hashing is unavailable", class = "saved_entity_unsupported")
    checkSavedEntityPath(store, ".pending")
    # Managed workers use request-owned TMPDIR so controller cleanup also
    # removes a partial hash after an uncatchable worker/process-group stop.
    scratch <- if (identical(Sys.getenv("MULTISCHOLAR_WORKFLOW_WORKER"), "1"))
        tempdir() else file.path(store$directory, ".pending")
    if (!dir.exists(scratch) && !dir.create(scratch))
        cli::cli_abort("cannot create hash scratch directory", class = "saved_entity_write")
    output <- tempfile("column-hash-", tmpdir = scratch)
    on.exit(unlink(output), add = TRUE)
    # digest's binary serializeVersion=2 identity omits the 14-byte R header.
    # pipefail prevents a broken producer from being accepted as a valid digest.
    pipeline <- paste(shQuote(commands[["dd"]]),
        "bs=65536 skip=14 iflag=skip_bytes status=none |",
        shQuote(commands[["sha256sum"]]), ">", shQuote(output))
    connection <- pipe(paste(shQuote(commands[["bash"]]), "-o pipefail -c",
        shQuote(pipeline)), open = "wb")
    on.exit(if (!is.null(connection)) close(connection), add = TRUE)
    serialize(value, connection, ascii = FALSE, xdr = TRUE, version = 2L)
    status <- suppressWarnings(close(connection))
    connection <- NULL
    if (!identical(as.integer(status), 0L) || !file.exists(output) || file.size(output) > 128L)
        cli::cli_abort("bounded serialized hashing failed", class = "saved_entity_write")
    result <- readLines(output, n = 2L, warn = FALSE)
    if (length(result) != 1L || !grepl("^[a-f0-9]{64}  -$", result))
        cli::cli_abort("serialized hash output is invalid", class = "saved_entity_corrupt")
    substr(result, 1L, 64L)
}

#' Persist or reuse an exact DIA table component in the existing registry.
#' @param context Current owner.
#' @param table Complete table already owned by a parser or scientific operation.
#' @param backend Existing writer backend; fast_r is explicitly qualified only.
#' @param consume Optional producer summary callback for each bounded output batch.
#' @return Table descriptor, never a second resident table owner.
#' @noRd
saveDIATableComponent <- function(context, table, backend = "arrow", consume = NULL) {
    table <- prepareDIAParserTable(table)
    metadata <- workflowSavedTableMetadata(table, ".rownames", 3L)
    # Equivalent attribute orders and ALTREP representations have the same
    # scientific values. Hash the codec contract and individual column values.
    column_hashes <- vapply(table, function(value) {
        element_names <- names(value)
        if (!is.null(attributes(value))) attributes(value) <- NULL
        if (identical(Sys.info()[["sysname"]], "Linux"))
            hashWorkflowSerializedValue(list(value = value, names = element_names), context$store) else
            digest::digest(list(value = value, names = element_names), algo = "sha256", serializeVersion = 2L)
    }, character(1))
    content <- digest::digest(list(metadata = metadata, columns = column_hashes,
        row_names = if (isTRUE(metadata$automatic_row_names)) NULL else rownames(table)),
        algo = "sha256", serializeVersion = 2L)
    identity <- newSavedWorkflowScientificIdentity(content,
        "component-owned", "DIA-table-codec-v2", "all-ordered-rows-and-columns", context$software)
    hit <- resolveSavedWorkflowEntity(context$store, "result_table", identity, hydrate = FALSE)
    if (identical(hit$status, "hit")) {
        if (!is.null(consume)) for (start in seq.int(1L, max(1L, nrow(table)), by = 10000L)) {
            at <- if (nrow(table)) seq.int(start, min(nrow(table), start + 9999L)) else integer()
            consume(table[at, , drop = FALSE])
        }
        return(hit$reference)
    }
    reference <- saveWorkflowTable(context$store, table, basename(tempfile("dia-table-")),
        identity, ".rownames", metadata_version = 3L,
        write_fn = function(path, metadata) writeSavedWorkflowParquet(table, ".rownames",
            path, metadata = metadata, backend = backend, consume = consume))
    context$writes <- context$writes + 1L
    reference
}

#' Apply existing design renames, selection and numeric conversion to saved rows.
#' @param manager,reference Current owner and immutable assay table.
#' @param samples Ordered design sample names, or NULL for no selection.
#' @param rename Named original-to-new sample names, applied in their given order.
#' @param numeric_columns Existing DIA quantities requiring numeric conversion.
#' @param strip_parser Match the design builder's existing dplyr reconstruction.
#' @return New table reference, or the unchanged reference when nothing changes.
#' @noRd
transformDIADesignTable <- function(manager, reference, samples = NULL, rename = NULL,
    numeric_columns = character(), strip_parser = FALSE) {
    context <- manager$.context
    identity <- newSavedWorkflowScientificIdentity(reference$entity_id,
        digest::digest(samples, algo = "sha256"),
        digest::digest(list(rename, numeric_columns, strip_parser), algo = "sha256"),
        "DIA-design-ordered-transform", context$software)
    hit <- resolveSavedWorkflowEntity(context$store, "result_table", identity, hydrate = FALSE)
    if (identical(hit$status, "hit")) return(hit$reference)
    if (is.null(samples) && !length(rename) && !length(numeric_columns) && !strip_parser) {
        validateSavedWorkflowReference(context$store, reference, payload = FALSE)
        return(reference)
    }
    scope <- newWorkflowTableReadScope(context$store, reference, reader = "arrow", metadata_only = TRUE)
    on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
    metadata <- scope$metadata
    if (!all(numeric_columns %in% metadata$columns))
        cli::cli_abort("design conversion column is missing", class = "saved_entity_invalid")
    numeric_columns <- numeric_columns[!vapply(numeric_columns, function(name) {
        column <- metadata$column_metadata[[name]]
        identical(column$typeof, "double") && identical(column$class, "numeric") &&
            !length(column$attributes) && !name %in% names(metadata$name_columns)
    }, logical(1))]
    if (length(rename)) {
        same <- !is.na(rename) & !is.na(names(rename)) & unname(rename) == names(rename)
        rename <- rename[!same]
    }
    if (strip_parser && identical(metadata$class, c("spec_tbl_df", "tbl_df", "tbl", "data.frame"))) {
        metadata$class <- c("tbl_df", "tbl", "data.frame")
        metadata$attributes[c("spec", "problems")] <- NULL
        if (!is.null(metadata$attribute_order)) metadata$attribute_order <-
            c("names", "row.names", "class", names(metadata$attributes))
    }
    if (!is.null(samples) && !length(rename)) {
        imported <- getDIAImportMetadata(manager)
        if (identical(imported$reference_id, scope$physical_reference$entity_id)) {
            observed <- as.character(imported$samples$Run)
        } else {
            releaseWorkflowTableReadScope(scope)
            scope <- newWorkflowTableReadScope(context$store, reference, reader = "arrow")
            observed <- character()
            streamSavedWorkflowTable(scope, "Run", function(data, offset, row_names) {
                observed <<- unique(c(observed, as.character(data$Run)))
                if (as.numeric(object.size(observed)) > 8388608)
                    cli::cli_abort("sample catalogue exceeds metadata limits", class = "saved_entity_limit")
            })
        }
        if (all(observed %in% as.character(samples))) samples <- NULL
    }
    if (is.null(samples) && !length(rename) && !length(numeric_columns)) {
        result <- saveWorkflowTableComponent(scope, metadata, basename(tempfile("dia-design-metadata-")), identity)
        return(result)
    }
    releaseWorkflowTableReadScope(scope)
    scope <- newWorkflowTableReadScope(context$store, reference, reader = "arrow")
    if (!metadata$schema %in% c(2L, 3L))
        cli::cli_abort("design streaming requires the component table codec", class = "saved_entity_unsupported")
    metadata$schema <- 3L
    metadata$codec <- "parquet_duckdb_v3"
    metadata$rows <- NULL
    metadata$name_columns <- metadata$name_columns[setdiff(names(metadata$name_columns), numeric_columns)]
    for (column in numeric_columns) metadata$column_metadata[[column]] <- list(
        typeof = "double", class = "numeric", levels = NULL, ordered = FALSE, attributes = list())
    result <- saveWorkflowTable(context$store, NULL, basename(tempfile("dia-design-")),
        identity, metadata$key_column, metadata = metadata, write_fn = function(path, restoration) {
            writer <- openWorkflowTableWriter(path, restoration, context$store$progress)
            on.exit(abortWorkflowTableWriter(writer), add = TRUE)
            streamWorkflowTableRows(scope, scope$metadata$columns, function(data, offset, row_names) {
                for (index in seq_along(rename)) {
                    selected <- !is.na(data$Run) & data$Run == names(rename)[[index]]
                    data$Run[selected] <- unname(rename[[index]])
                }
                for (column in numeric_columns) data[[column]] <- suppressWarnings(as.numeric(data[[column]]))
                keep <- if (is.null(samples)) seq_len(nrow(data)) else which(data$Run %in% samples)
                appendWorkflowTableWriter(writer, data, keep)
            })
            finishWorkflowTableWriter(writer)
        })
    context$writes <- context$writes + 1L
    result
}

#' Stream the required cleaned-data TSV from the accepted table.
#' @param manager,reference Current owner and table.
#' @param path Existing design artifact destination.
#' @return Destination invisibly, after all batches and the database lease close.
#' @noRd
writeDIADesignTable <- function(manager, reference, path, write_fn = NULL) {
    scope <- newWorkflowTableReadScope(manager$.context$store, reference,
        reader = if (reference$metadata$shape[[2L]] > 64L) "duckdb" else "arrow")
    on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
    temporary <- tempfile(".data-cln-", tmpdir = dirname(path))
    on.exit(unlink(temporary), add = TRUE)
    output_owner <- list2env(list(path = temporary, advised_bytes = 0), parent = emptyenv())
    consume <- function(data, offset, row_names) {
        if (is.null(write_fn)) {
            utils::write.table(data, file = temporary, sep = "\t", row.names = FALSE, quote = FALSE,
                col.names = offset == 0L, append = offset != 0L)
        } else write_fn(data, temporary, append = offset > 0L, col_names = offset == 0L)
        releaseWorkflowWrittenPrefix(output_owner)
    }
    streamWorkflowTableRows(scope, scope$metadata$columns, consume)
    syncWorkflowExportPaths(temporary, dirname(path))
    releaseWorkflowWrittenPrefix(output_owner, final = TRUE)
    if (!file.rename(temporary, path)) cli::cli_abort("cannot publish cleaned table", class = "saved_entity_write")
    syncWorkflowExportPaths(c(path, dirname(path)), dirname(path))
    invisible(path)
}

#' Read only the distinct protein identifiers needed by design annotations.
#' @param manager,reference Current owner and accepted table.
#' @param column Declared protein identifier column.
#' @return Small projected identifier table in first-observed order.
#' @noRd
readDIAAnnotationIdentifiers <- function(manager, reference, column) {
    scope <- newWorkflowTableReadScope(manager$.context$store, reference, reader = "arrow")
    on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
    values <- readSavedWorkflowTablePrototype(scope)[[column]]
    streamSavedWorkflowTable(scope, column, function(data, offset, row_names) {
        incoming <- data[[column]]
        at <- which(!duplicated(incoming) & is.na(match(incoming, values)))
        values <<- c(values, incoming[at])
        if (length(values) > 10000L || as.numeric(object.size(values)) > 67108864)
            cli::cli_abort("annotation identifiers exceed read limits", class = "saved_entity_limit")
    })
    result <- stats::setNames(list(values), column)
    attr(result, "row.names") <- .set_row_names(length(values))
    class(result) <- "data.frame"
    metadata <- scope$metadata
    metadata$columns <- column
    metadata$column_metadata <- metadata$column_metadata[column]
    metadata$name_columns <- metadata$name_columns[intersect(names(metadata$name_columns), column)]
    result <- completeWorkflowTableBatch(result, metadata)
    attr(result, "row.names") <- .set_row_names(length(values))
    result
}

#' Create the initial peptide checkpoint from a saved assay and a small design.
#' @param manager Current owner.
#' @param reference Parsed and selected assay table.
#' @param design,args Existing design and effective configuration.
#' @return Component checkpoint specification; no assay payload is reconstructed.
#' @noRd
newDIAPeptideCheckpoint <- function(manager, reference, design, args) {
    context <- manager$.context
    scope <- newWorkflowTableReadScope(context$store, reference, metadata_only = TRUE)
    on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
    empty <- readSavedWorkflowTablePrototype(scope)
    rememberDIAWorkers(context, args)
    args <- prepareDIAStateRecord(list(data = NULL, config = args))$config
    prototype <- PeptideQuantitativeDataDiann(empty, design,
        technical_replicate_id = "tech_rep_group", args = args)
    slots <- stats::setNames(lapply(methods::slotNames(prototype),
        \(name) methods::slot(prototype, name)), methods::slotNames(prototype))
    slots[c("peptide_data", "design_matrix")] <- NULL
    list(schema = "dia.component_checkpoint.v2", class = "PeptideQuantitativeData",
        declaration = workflowDeclarationField(workflowClassDeclaration("PeptideQuantitativeData")),
        slots = slots, components = list(peptide_data = reference,
            design_matrix = saveDIATableComponent(context, design)), support = list(),
        dimensions = list(rows = scope$metadata$rows, columns = scope$metadata$columns))
}

#' Select a saved state without invoking the public materializing Revert result.
#' @param manager Existing owner.
#' @param state_name Existing state name.
#' @return Selected descriptor, or the ordinary legacy result on other routes.
#' @noRd
selectDIAStateReference <- function(manager, state_name) {
    if (!usesDIAReferenceStates(manager)) return(manager$revertToState(state_name))
    reference <- getDIAStateReference(manager, state_name)
    if (is.null(reference)) cli::cli_abort("state does not exist", class = "saved_entity_missing")
    snapshot <- manager$.context$snapshot
    index <- match(state_name, unlist(snapshot$state_history))
    if (is.na(index)) cli::cli_abort("state is outside selected history", class = "saved_entity_scope")
    snapshot$state_history <- snapshot$state_history[seq_len(index)]
    snapshot$current_state <- state_name
    snapshot <- truncateDIAStepProvenance(snapshot)
    publishDIAStateHead(manager$.context, snapshot)
    reference
}

#' Release a consumed predecessor after its successor is durably current.
#' @param manager Existing native or experimental state manager.
#' @param reference Immutable predecessor descriptor captured before the action.
#' @return TRUE invisibly if its working binding was removed; no durable data changes.
#' @noRd
releaseDIAWorkingReference <- function(manager, reference) {
    if (!inherits(manager, "DIAStateAdapter") || is.null(reference)) return(invisible(FALSE))
    context <- manager$.context
    checkWorkflowEntityStore(context$store)
    if (!isSavedWorkflowReference(reference) || !identical(reference$owner, context$store$identity)) {
        cli::cli_abort("consumed DIA reference has a different owner", class = "saved_entity_scope")
    }
    current <- context$snapshot$states[[context$snapshot$current_state]]
    key <- reference$entity_id
    # A failed commit or same-generation call must not evict the active input.
    if (identical(key, current$entity_id) || !exists(key, context$working, inherits = FALSE)) {
        return(invisible(FALSE))
    }
    rm(list = key, envir = context$working)
    observeDIAOwnerBoundary(context, "consumed_reference_released")
    invisible(TRUE)
}

#' Materialize a legacy list explicitly, never for ordinary GUI metadata checks.
#' @param context Private adapter owner.
#' @param share_immutable Reuse exactly equal immutable tables within an export.
#' @return Ordinary complete list with original names and records.
#' @noRd
materializeDIAStates <- function(context, share_immutable = FALSE) {
    context$full_reads <- context$full_reads + 1L
    tables <- list()
    if (share_immutable) context$export_shared_tables <- 0L
    lapply(context$snapshot$states, function(reference) {
        record <- materializeDIAStateRecord(context, reference)
        if (!share_immutable || !isS4(record$data)) return(record)
        for (field in intersect(workflowClassDeclaration(record$data)$tables, methods::slotNames(record$data))) {
            table <- methods::slot(record$data, field)
            if (hasMutableWorkflowEntityOwner(table)) next
            previous <- tables[[field]]
            matches <- which(vapply(previous, \(value)
                identical(dim(value), dim(table)) && identical(value, table), logical(1)))
            if (length(matches)) {
                methods::slot(record$data, field) <- previous[[matches[[1L]]]]
                context$export_shared_tables <- context$export_shared_tables + 1L
            } else tables[[field]] <<- c(previous, list(table))
        }
        # Replace, rather than additionally retain, the export's decoded copy.
        if (context$scope_depth > 0L) assign(reference$entity_id, record, context$working)
        record
    })
}

#' Persist a whole public history replacement before accepting any of it.
#' @param context Private adapter owner.
#' @param records Ordinary complete legacy list, including the empty list.
#' @param topology Optional already validated current/history fields for atomic session restore.
#' @return NULL invisibly; current/history fields retain normal separate-assignment semantics.
#' @noRd
replaceDIAStates <- function(context, records, topology = NULL) {
    if (!is.list(records) || length(records) &&
        (is.null(names(records)) || anyNA(names(records)) || any(!nzchar(names(records))) ||
            anyDuplicated(names(records)))) {
        cli::cli_abort("invalid DIA history replacement", class = "saved_entity_invalid")
    }
    snapshot <- context$snapshot
    references <- lapply(names(records), function(name) {
        saveDIAStateValue(context, prepareDIAStateRecord(records[[name]]), name)
    })
    snapshot$states <- stats::setNames(references, names(records))
    if (!is.null(topology)) snapshot[names(topology)] <- topology
    publishDIAStateHead(context, snapshot)
    invisible(NULL)
}

#' Create the experimental compatible state facade over the original R6 bookkeeping.
#' @param store Existing project-owned store with zero optional payload cache.
#' @param native Ordinary WorkflowState; newly constructed by default.
#' @param reopen Read this store's previously accepted private head when TRUE.
#' @return Internal state facade; complete legacy exports, not this environment, are portable.
#' @noRd
newDIAStateAdapter <- function(store, native = NULL, reopen = FALSE) {
    checkWorkflowEntityStore(store)
    if (store$max_cache_bytes != 0 || !is.logical(reopen) || length(reopen) != 1L || is.na(reopen)) {
        cli::cli_abort("invalid DIA state adapter policy", class = "saved_entity_invalid")
    }
    if (is.null(native)) native <- WorkflowState$new()
    context <- new.env(parent = emptyenv())
    context$store <- store
    context$native <- native
    context$head_path <- file.path(store$directory, "dia-state-head.rds")
    context$head <- NULL
    context$working <- new.env(parent = emptyenv())
    context$scope_depth <- 0L
    context$on_owner_boundary <- NULL
    context$reclaim_enabled <- FALSE
    context$input_file_leases <- list()
    context$input_files_ready <- FALSE
    context$reclaim_pending <- context$reclaim_last <- context$reclaim_result <- NULL
    context$borrowed_workers <- context$owned_workers <- list()
    context$writes <- context$hydrates <- context$full_reads <- 0L
    context$full_s4_reads <- 0L
    binary <- system.file("R", "MultiScholaR.rdb", package = "MultiScholaR")
    context$software <- if (nzchar(binary)) digest::digest(file = binary, algo = "sha256") else
        digest::digest(list(as.character(utils::packageVersion("MultiScholaR")),
            body(newDIAStateAdapter)), algo = "sha256")
    producer_packages <- c("iq", "limma", "ruv", "mixOmics", "qvalue", "matrixStats", "statmod",
        "dplyr", "tidyr", "arrow", "duckdb", "ggplot2", "gridExtra", "ComplexHeatmap", "dynamicTreeCut", "Glimma")
    versions <- stats::setNames(lapply(producer_packages, function(package) {
        tryCatch(as.character(utils::packageVersion(package)), error = function(error) NA_character_)
    }), producer_packages)
    context$software <- digest::digest(list(package_binary = context$software,
        native_libraries = workflowNativeLibraryIdentity(), R = R.version[c("major", "minor", "platform")], producer_packages = versions,
        numerical_libraries = extSoftVersion()[intersect(c("BLAS", "LAPACK"), names(extSoftVersion()))]), algo = "sha256")
    context$snapshot <- list(identity = store$identity, states = list(), audit = list(),
        current_state = native$current_state, state_history = native$state_history,
        workflow_type = native$workflow_type, audit_enabled = native$audit_enabled)
    # KIT-025: the placeholder `initial` state is the import's; its provenance is written with the first head.
    initial <- diaStepProvenanceEntry(context$snapshot, list(operation = "import", variant = NA_character_, parameter_names = character()))
    context$snapshot <- initial$snapshot
    context$snapshot$step_provenance <- list(initial = initial$entry)
    if (isTRUE(reopen)) {
        if (!file.exists(context$head_path) || file.size(context$head_path) > 65536) {
            cli::cli_abort("DIA state head is missing or oversized", class = "saved_entity_missing")
        }
        context$head <- readRDS(context$head_path)
        snapshot <- resolveDIAStateValue(context, context$head)
        if (!identical(snapshot$identity, store$identity) || !is.list(snapshot$states) ||
            !is.list(snapshot$audit) || !(is.list(snapshot$state_history) ||
                is.character(snapshot$state_history))) {
            cli::cli_abort("DIA state head identity differs", class = "saved_entity_scope")
        }
        # KIT-025: a store written before provenance was recorded (or adopted from a legacy session) opens with
        # a `legacy` entry for every state that has none; the entries are written with the next head.
        for (name in setdiff(unlist(snapshot$state_history), names(snapshot$step_provenance))) {
            # KIT-027 audit PH-04: the entry is marked, since its build, defaults and time are the reopening's.
            legacy <- diaStepProvenanceEntry(snapshot, list(operation = "legacy", variant = NA_character_, parameter_names = character(),
                backfilled = TRUE))
            snapshot <- legacy$snapshot
            snapshot$step_provenance[[name]] <- legacy$entry
        }
        # KIT-026: inputs recorded before roles existed open with `role = NA`.
        for (name in names(snapshot$file_dependencies))
            if (is.null(snapshot$file_dependencies[[name]]$role)) snapshot$file_dependencies[[name]]$role <- NA_character_
        context$snapshot <- snapshot
        context$reference_states <- isTRUE(snapshot$reference_states)
        # Rebuild only the bounded index from accepted roots, never from orphan stages.
        for (reference in collectDIANativeReferences(list(context$head, snapshot))) {
            registerSavedWorkflowMetadata(store, reference)
        }
        snapshot <- NULL
    } else {
        if (file.exists(context$head_path)) {
            cli::cli_abort("DIA state head already exists; reopen explicitly", class = "saved_entity_conflict")
        }
        replaceDIAStates(context, native$states)
        if (length(native$audit_records)) {
            snapshot <- context$snapshot
            snapshot$audit <- lapply(native$audit_records, \(record) list(reference =
                saveDIAStateValue(context, record, "audit"), field = NULL))
            publishDIAStateHead(context, snapshot)
            snapshot <- NULL
        }
    }
    native$states <- native$audit_records <- list()
    adapter <- new.env(parent = emptyenv())
    class(adapter) <- c("DIAStateAdapter", "environment")
    adapter$.context <- context
    makeActiveBinding("states", function(value) {
        if (missing(value)) return(materializeDIAStates(context))
        replaceDIAStates(context, value)
    }, adapter)
    for (field in c("current_state", "state_history", "workflow_type", "audit_enabled")) {
        bindDIAStateMetadata(adapter, context, field)
    }
    makeActiveBinding("audit_records", function(value) {
        if (missing(value)) return(adapter$getAuditRecords())
        if (!is.list(value)) cli::cli_abort("invalid DIA audit replacement", class = "saved_entity_invalid")
        snapshot <- context$snapshot
        snapshot$audit <- lapply(value, \(record) list(reference =
            saveDIAStateValue(context, record, "audit"), field = NULL))
        publishDIAStateHead(context, snapshot)
    }, adapter)
    adapter$getState <- function(state_name = NULL) {
        if (is.null(state_name)) state_name <- context$snapshot$current_state
        reference <- context$snapshot$states[[state_name]]
        if (is.null(reference)) return(NULL)
        materializeDIAStateRecord(context, reference)$data
    }
    adapter$saveState <- function(state_name, s4_data_object, config_object, description,
        audit_metadata = NULL) {
        if (usesDIAReferenceStates(adapter) && isS4(s4_data_object)) {
            checkpoint <- checkpointDIAExecutionResult(adapter, s4_data_object)
            return(invisible(commitDIAReferenceState(adapter, state_name, checkpoint,
                config_object, description)))
        }
        saveDIAStateRecord(context, state_name, s4_data_object, config_object, description, audit_metadata)
    }
    adapter$revertToState <- function(state_name) {
        if (!is.null(context$store$execution_owner)) cancelWorkflowAction(context$store$execution_owner, "reverted")
        reference <- context$snapshot$states[[state_name]]
        if (is.null(reference)) stop("State not found: ", state_name)
        record <- materializeDIAStateRecord(context, reference)
        snapshot <- context$snapshot
        index <- which(snapshot$state_history == state_name)
        snapshot$state_history <- snapshot$state_history[1:index]
        snapshot$current_state <- state_name
        snapshot <- truncateDIAStepProvenance(snapshot)
        publishDIAStateHead(context, snapshot)
        record$data
    }
    adapter$getHistory <- function() unlist(context$snapshot$state_history)
    adapter$getStateAudit <- function(state_name = NULL) {
        if (is.null(state_name)) state_name <- context$snapshot$current_state
        reference <- context$snapshot$states[[state_name]]
        metadata <- if (!is.null(reference)) resolveDIAStateValue(context, reference)$audit_metadata
        if (is.null(metadata)) list(provenance_status = "legacy_untracked_state") else metadata
    }
    adapter$getAuditRecords <- function() lapply(context$snapshot$audit, function(entry) {
        value <- resolveDIAStateValue(context, entry$reference)
        if (is.null(entry$field)) value else value[[entry$field]]
    })
    adapter$setWorkflowType <- function(type) {
        # KIT-020: the on-disk owner serves the DIA-NN routes (KIT-046: DIA and DIA_LIMPA) and the protein-level LFQ and TMT routes.
        if (!is.character(type) || length(type) != 1L || !type %in% c("DIA", "DIA_LIMPA", "LFQ", "TMT")) {
            cli::cli_abort("the on-disk session serves the DIA, DIA_LIMPA, LFQ and TMT workflow types", class = "saved_entity_scope")
        }
        adapter$workflow_type <- type
    }
    adapter
}

#' Bind one scalar/topology field without capturing an assay-bearing call frame.
#' @param adapter,context Facade and private owner.
#' @param field Exact existing public metadata field.
#' @return NULL invisibly.
#' @noRd
bindDIAStateMetadata <- function(adapter, context, field) {
    force(field)
    makeActiveBinding(field, function(value) {
        if (missing(value)) return(context$snapshot[[field]])
        if (field == "state_history" && !(is.list(value) || is.character(value)) ||
            field == "audit_enabled" && (!is.logical(value) || length(value) != 1L || is.na(value)) ||
            field %in% c("current_state", "workflow_type") &&
                (!is.character(value) || length(value) != 1L || is.na(value))) {
            cli::cli_abort("invalid DIA state metadata", class = "saved_entity_invalid")
        }
        snapshot <- context$snapshot
        snapshot[[field]] <- value
        publishDIAStateHead(context, snapshot)
    }, adapter)
    invisible(NULL)
}

#' Reuse native state bookkeeping, persisting before the facade advances.
#' @param context Private owner.
#' @param state_name,s4_data_object,config_object,description,audit_metadata Existing saveState inputs.
#' @return NULL invisibly.
#' @noRd
saveDIAStateRecord <- function(context, state_name, s4_data_object, config_object, description,
    audit_metadata = NULL) {
    checkSavedEntityString(state_name, "state_name")
    checkWorkflowEntityStore(context$store)
    native <- context$native
    snapshot <- context$snapshot
    native$current_state <- snapshot$current_state
    native$state_history <- snapshot$state_history
    native$audit_enabled <- snapshot$audit_enabled
    on.exit({ native$states <- native$audit_records <- list() }, add = TRUE)
    native$saveState(state_name, s4_data_object, config_object, description, audit_metadata)
    if (isS4(s4_data_object)) rememberDIAWorkers(context, methods::slot(s4_data_object, "args"))
    rememberDIAWorkers(context, config_object)
    record <- prepareDIAStateRecord(native$states[[state_name]])
    reference <- saveDIAStateValue(context, record, state_name)
    snapshot$states[[state_name]] <- reference
    dependency_fields <- if (identical(context$session_format, "native_v1"))
        intersect(c("iq_output_file", "peptide_input_file"), names(record$config)) else character()
    for (field in dependency_fields) {
        path <- record$config[[field]]
        if (is.character(path) && length(path) == 1L && !is.na(path)) {
            description <- describeDIANativeFile(path)
            snapshot$file_dependencies[[field]] <- c(list(role = sub("_file$", "", field), source = path,
                destination = paste("execution", field, basename(path), sep = "/")), description)
        }
    }
    if (isS4(record$data)) {
        declaration <- workflowClassDeclaration(record$data)
        dimensions <- workflowAssayDimensions(record$data, declaration)
        table <- if (!is.null(declaration$assay)) methods::slot(record$data, declaration$assay) else NULL
        snapshot$sizes[[state_name]] <- list(s4_object_size_bytes = as.numeric(object.size(record$data)),
            table_object_size_bytes = as.numeric(object.size(table)), rows = dimensions$rows,
            columns = length(dimensions$columns),
            note = "object.size of owned values, not unique whole-process heap or RSS")
        if (is.function(declaration$describe)) {
            snapshot$state_metadata[[state_name]] <- declaration$describe(record$data)
        }
    }
    snapshot$current_state <- native$current_state
    snapshot$state_history <- native$state_history
    id <- record$audit_metadata$record_id
    if (!is.null(id)) snapshot$audit[[id]] <- list(reference = reference, field = "audit_metadata")
    snapshot <- recordDIAStepProvenance(context, snapshot, state_name)
    publishDIAStateHead(context, snapshot)
    invisible(NULL)
}

#' Capture export scalars while the accepted protein state is already available.
#' @param object Existing protein S4; no new assay read is performed.
#' @return Plain bounded metadata, never a table, worker or model.
#' @noRd
describeDIAProteinState <- function(object) {
    args <- methods::slot(object, "args")
    global <- args$globalParameters
    table <- methods::slot(object, "protein_quant_table")
    text <- function(value) if (is.character(value) && length(value) == 1L && !is.na(value)) value else NULL
    list(proteins = length(unique(table[[methods::slot(object, "protein_id_column")]])),
        samples = nrow(methods::slot(object, "design_matrix")),
        workflow_type = text(global$workflow_type), report_template = text(global$report_template),
        # KIT-046: the limpa route is a workflow type; neither a flag nor the presence of limpa arguments names it.
        use_limpa = identical(text(global$workflow_type), "DIA_LIMPA"))
}

#' Run one action/export with a scoped decode map, released even on failure.
#' @param manager Native or experimental state manager.
#' @param operation Zero-argument existing action.
#' @return The operation's original result.
#' @noRd
withDIAStateScope <- function(manager, operation) {
    if (!inherits(manager, "DIAStateAdapter")) return(operation())
    context <- manager$.context
    checkWorkflowEntityStore(context$store)
    # KIT-022: the scope names its owner for `getDIAActiveManager()` while it is open.
    depth <- length(.diaActiveScope$stack)
    .diaActiveScope$stack[[depth + 1L]] <- manager
    on.exit(.diaActiveScope$stack <- .diaActiveScope$stack[seq_len(depth)], add = TRUE)
    action <- beginWorkflowTableAction(context$store)
    on.exit(endWorkflowTableAction(context$store, action), add = TRUE)
    context$scope_depth <- context$scope_depth + 1L
    completed <- FALSE
    on.exit({
        if (!completed) context$reclaim_pending <- NULL
        context$scope_depth <- if (isTRUE(context$store$closed)) 0L else context$scope_depth - 1L
        if (context$scope_depth == 0L) {
            tryCatch(observeDIAOwnerBoundary(context, "operation_release"),
                finally = { context$working <- new.env(parent = emptyenv()) })
            observeDIAOwnerBoundary(context, "operation_released")
        }
    }, add = TRUE)
    observeDIAOwnerBoundary(context, "operation_enter")
    result <- withVisible(operation())
    completed <- TRUE
    if (result$visible) result$value else invisible(result$value)
}

#' Read one legacy record without materializing the complete history.
#' @param manager Native or experimental state manager.
#' @param name Requested logical state.
#' @return Complete record or NULL, as ordinary list extraction would return.
#' @noRd
getDIAStateRecord <- function(manager, name) {
    if (!inherits(manager, "DIAStateAdapter")) return(manager$states[[name]])
    reference <- manager$.context$snapshot$states[[name]]
    if (is.null(reference)) return(NULL)
    materializeDIAStateRecord(manager$.context, reference)
}

#' Inspect existence without triggering full-list compatibility hydration.
#' @param manager Native or experimental state manager.
#' @param name Logical state.
#' @return TRUE if a record exists, including an inactive branch record.
#' @noRd
hasDIAState <- function(manager, name) {
    if (!inherits(manager, "DIAStateAdapter")) return(!is.null(manager$states[[name]]))
    !is.null(manager$.context$snapshot$states[[name]])
}

#' Report bounded ownership counters without hydrating history.
#' @param manager Experimental state manager.
#' @return Scalar counts and store information; not a whole-heap census.
#' @noRd
diaStateOwnership <- function(manager) {
    context <- manager$.context
    list(states = length(context$snapshot$states), audits = length(context$snapshot$audit),
        state_sizes = context$snapshot$sizes,
        working_records = length(ls(context$working, all.names = TRUE)),
        native_records = length(context$native$states), native_audits = length(context$native$audit_records),
        writes = context$writes, hydrates = context$hydrates, full_reads = context$full_reads,
        full_s4_reads = context$full_s4_reads, reference_states = isTRUE(context$reference_states),
        component_states = length(context$snapshot$component_dependencies),
        store = workflowEntityCacheInfo(context$store))
}

#' Release the internal adapter without deleting accepted files or other owners.
#' @param manager Experimental state manager.
#' @return NULL invisibly.
#' @noRd
closeDIAStateAdapter <- function(manager) {
    if (!inherits(manager, "DIAStateAdapter")) return(invisible(NULL))
    context <- manager$.context
    for (cancel in context$callbacks) if (is.function(cancel)) cancel()
    context$callbacks <- list()
    context$on_commit <- NULL
    context$on_owner_boundary <- NULL
    context$reclaim_enabled <- FALSE
    context$reclaim_pending <- context$reclaim_last <- context$reclaim_result <- NULL
    context$scope_depth <- 0L
    context$owned_export <- context$export_pending <- NULL
    context$import_metadata_cache <- NULL
    context$working <- new.env(parent = emptyenv())
    context$native$states <- context$native$audit_records <- list()
    if (length(context$owned_workers)) {
        for (cluster in context$owned_workers) stopDIAWorkerCluster(cluster)
        context$owned_workers <- list()
    }
    context$borrowed_workers <- list()
    closeWorkflowEntityStore(context$store)
    invisible(NULL)
}

#' Find only this Shiny session's explicitly enabled experimental manager.
#' @return Internal manager or NULL; never a global live handle.
#' @noRd
getDIAActiveManager <- function() {
    # KIT-022: the innermost `withDIAStateScope` owner first (kit code, workers and scripts), then the Shiny
    # session's owner (the modules), so kit-level readers need no session.
    stack <- .diaActiveScope$stack
    if (length(stack)) return(stack[[length(stack)]])
    session <- shiny::getDefaultReactiveDomain()
    if (is.null(session)) return(NULL)
    manager <- session$userData$dia_capacity_manager
    if (inherits(manager, "DIAStateAdapter")) manager else NULL
}

#' The stack of owners of the open `withDIAStateScope` calls in this process (KIT-022).
#' @noRd
.diaActiveScope <- new.env(parent = emptyenv())
.diaActiveScope$stack <- list()

#' Bind a saved section's exact worker count only at scientific execution.
#' @param manager Native or experimental manager.
#' @param section Existing scientific configuration section.
#' @param state_name Optional historical state; NULL means current.
#' @return Executable S4; generic hydration itself never starts workers.
#' @noRd
getDIAExecutionState <- function(manager, section, state_name = NULL) {
    if (usesDIAReferenceStates(manager)) cli::cli_abort(
        "component actions must request their declared scientific inputs",
        class = "saved_entity_unsupported")
    object <- if (is.null(state_name)) manager$getState() else manager$getState(state_name)
    bindDIAExecutionObject(manager, object, section)
}

#' Bind an already resolved state or roll-up envelope to its existing backend.
#' @param manager Existing state owner.
#' @param object Hydrated S4, without changing its scientific slots.
#' @param section Existing scientific configuration section.
#' @return Executable S4 using the same recorded worker policy.
#' @noRd
bindDIAExecutionObject <- function(manager, object, section) {
    if (!inherits(manager, "DIAStateAdapter") || is.null(object)) return(object)
    args <- methods::slot(object, "args")
    # Rebind only declared execution-file roles on a transient execution copy.
    # Stored historical arguments retain their original provenance paths.
    dependencies <- manager$.context$native_dependencies
    for (entry in dependencies) for (name in names(args)) {
        if (!is.list(args[[name]])) next
        for (field in intersect(c("iq_output_file", "peptide_input_file"), names(args[[name]]))) {
            if (identical(args[[name]][[field]], entry$source_provenance)) {
                if (!identical(workflowTableFileIdentity(entry$path), entry$file_identity)) {
                    cli::cli_abort("restored execution dependency changed", class = "saved_entity_corrupt")
                }
                args[[name]][[field]] <- entry$path
            }
        }
    }
    if (length(dependencies)) methods::slot(object, "args") <- args
    spec <- args[[section]]$core_utilisation
    if (!inherits(spec, "DIAWorkerSpec")) return(object)
    if (!identical(spec, newDIAWorkerSpec(spec$workers))) {
        cli::cli_abort("invalid DIA worker specification", class = "saved_entity_invalid")
    }
    context <- manager$.context
    # These validated scientific helpers already execute in-process. The SQL
    # rollup also does not consume the historical multidplyr pool.
    local_sections <- c("peptideIntensityFiltering", "filterMinNumPeptidesPerProtein",
        "filterMinNumPeptidesPerSample", "removePeptidesWithOnlyOneReplicate",
        "peptideMissingValueImputation", "removeProteinsWithOnlyOneReplicate")
    if (section %in% local_sections ||
        identical(section, "rollUpPrecursorToPeptide") && isTRUE(context$sql_rollup)) {
        args[[section]]$core_utilisation <- NA
        methods::slot(object, "args") <- args
        return(object)
    }
    key <- as.character(spec$workers)
    cluster <- context$borrowed_workers[[key]]
    alive <- function(cluster) !is.null(cluster) &&
        all(vapply(cluster, \(worker) worker$is_alive(), logical(1)))
    if (!alive(cluster)) {
        cluster <- context$owned_workers[[key]]
        if (!alive(cluster)) {
            cluster <- startDIAWorkerCluster(spec)
            context$owned_workers[[key]] <- cluster
        }
    }
    args[[section]] <- bindDIAWorkerConfig(stats::setNames(list(args[[section]]), section), cluster)[[section]]
    methods::slot(object, "args") <- args
    object
}

#' Remember borrowed configuration workers without persisting or taking ownership.
#' @param context Private manager owner.
#' @param config Existing configuration list.
#' @return NULL invisibly.
#' @noRd
rememberDIAWorkers <- function(context, config) {
    for (section in names(config)) {
        value <- config[[section]]
        if (!is.list(value) || !inherits(value$core_utilisation, "multidplyr_cluster")) next
        cluster <- value$core_utilisation
        context$borrowed_workers[[as.character(length(cluster))]] <- cluster
    }
    invisible(NULL)
}

#' Store a named GUI payload using the same whole-record codec, not a new backend.
#' @param value Existing table, S4 or design-builder result.
#' @param label Named current consumer role.
#' @param manager Session manager; NULL preserves the native value outside a trial.
#' @return Private descriptor, or the untouched value on the native path.
#' @noRd
saveDIAOwnedValue <- function(value, label, manager = NULL) {
    if (is.null(manager)) manager <- getDIAActiveManager()
    if (!inherits(manager, "DIAStateAdapter")) return(value)
    if (is.null(value)) {
        if (usesDIAReferenceStates(manager) && !is.null(manager$.context$snapshot$gui_references[[label]])) {
            snapshot <- manager$.context$snapshot
            snapshot$gui_references[[label]] <- NULL
            publishDIAStateHead(manager$.context, snapshot)
        }
        return(NULL)
    }
    if (isSavedWorkflowReference(value)) {
        validateSavedWorkflowReference(manager$.context$store, value, payload = FALSE)
        if (usesDIAReferenceStates(manager) && !identical(manager$.context$snapshot$gui_references[[label]], value)) {
            snapshot <- manager$.context$snapshot
            snapshot$gui_references[[label]] <- value
            publishDIAStateHead(manager$.context, snapshot)
        }
        return(value)
    }
    if (usesDIAReferenceStates(manager) && is.data.frame(value)) {
        reference <- saveDIATableComponent(manager$.context, value)
        snapshot <- manager$.context$snapshot
        snapshot$gui_references[[label]] <- reference
        publishDIAStateHead(manager$.context, snapshot)
        return(reference)
    }
    if (usesDIAReferenceStates(manager) && isS4(value)) {
        checkpoint <- checkpointDIAExecutionResult(manager, value)
        reference <- saveDIAStateValue(manager$.context,
            list(data = checkpoint, config = NULL, description = label), paste0("owner:", label), codec = "rds")
        snapshot <- manager$.context$snapshot
        snapshot$owned_components[[reference$entity_id]] <- collectDIANativeReferences(
            list(reference, checkpoint$components, checkpoint$support, checkpoint$products, checkpoint$dependencies))
        snapshot$gui_references[[label]] <- reference
        publishDIAStateHead(manager$.context, snapshot)
        return(reference)
    }
    if (is.list(value) && !is.data.frame(value) && !is.null(value$config_list)) {
        rememberDIAWorkers(manager$.context, value$config_list)
        value$config_list <- prepareDIAStateRecord(list(data = NULL, config = value$config_list))$config
    }
    if (is.list(value) && !is.data.frame(value)) {
        for (field in intersect(names(value), c("data_cln", "design_matrix"))) {
            value[[field]] <- prepareDIAParserTable(value[[field]])
        }
    }
    if (isS4(value)) rememberDIAWorkers(manager$.context, methods::slot(value, "args"))
    record <- prepareDIAStateRecord(list(data = value, config = NULL, description = label))
    saveDIAStateValue(manager$.context, record, paste0("owner:", label))
}

#' Resolve a private owner reference without altering the caller's stored value.
#' @param value Native value or saved reference.
#' @param manager Session manager; NULL resolves the current session.
#' @return Original scientific value, never a descriptor disguised as an S4.
#' @noRd
resolveDIAOwnedValue <- function(value, manager = NULL) {
    if (!isSavedWorkflowReference(value)) return(value)
    if (is.null(manager)) manager <- getDIAActiveManager()
    if (!inherits(manager, "DIAStateAdapter")) {
        cli::cli_abort("saved GUI payload has no owning DIA session", class = "saved_entity_scope")
    }
    if (isWorkflowTableReference(value)) {
        return(materializeDIATable(manager$.context, value))
    }
    if (usesDIAReferenceStates(manager) && identical(value$payload$codec, "rds")) {
        return(readDIAExecutionComponents(manager, value))
    }
    resolveDIAStateValue(manager$.context, value)$data
}

#' Return a saved GUI reference before a scientific product is looked up.
#' @param owner,name Existing GUI owner and payload field.
#' @return Reference in component mode, original payload on legacy paths.
#' @noRd
getDIAOwnedInput <- function(owner, name) {
    if (usesDIAReferenceStates(getDIAActiveManager())) owner[[name]] else getDIAOwnedPayload(owner, name)
}

#' Read an inventoried GUI payload holder, preserving ordinary callers.
#' @param owner Existing reactiveValues, environment or list.
#' @param name Existing field name.
#' @return Original payload.
#' @noRd
getDIAOwnedPayload <- function(owner, name) resolveDIAOwnedValue(owner[[name]])

#' Update one inventoried GUI payload holder after durable storage succeeds.
#' @param owner Existing mutable GUI owner.
#' @param name Existing field name.
#' @param value Original payload, or NULL to clear the holder.
#' @return Updated owner invisibly, preserving list as well as mutable-owner semantics.
#' @noRd
setDIAOwnedPayload <- function(owner, name, value) {
    manager <- getDIAActiveManager()
    previous <- owner[[name]]
    reference <- saveDIAOwnedValue(value, name, manager)
    owner[[name]] <- reference
    releaseDIAReplacedOwner(manager, previous, reference)
    invisible(owner)
}

#' Evict only a successfully replaced, same-project working binding.
#' @param manager Current DIA manager or NULL on the ordinary path.
#' @param previous,reference Previous and newly installed owner values.
#' @return NULL invisibly; durable files and other reference aliases are untouched.
#' @noRd
releaseDIAReplacedOwner <- function(manager, previous, reference) {
    if (inherits(manager, "DIAStateAdapter") && isSavedWorkflowReference(previous) &&
        identical(previous$owner, manager$.context$store$identity) &&
        !identical(previous, reference)) {
        releaseDIAWorkingReference(manager, previous)
    }
    invisible(NULL)
}

#' Keep a design-builder reactive value cold without changing its get/set contract.
#' @param value Initial native value.
#' @param label Existing consumer role.
#' @param reactiveValFn Existing constructor, injectable by module tests.
#' @param snapshot_field Optional existing-session dependency root for native restore.
#' @return Function whose getter returns the original payload.
#' @noRd
newDIASavedReactiveVal <- function(value, label, reactiveValFn = shiny::reactiveVal,
    snapshot_field = NULL) {
    manager <- getDIAActiveManager()
    reference <- reactiveValFn(saveDIAOwnedValue(value, label, manager))
    value <- NULL
    result <- function(value) {
        manager <- getDIAActiveManager()
        if (missing(value)) return(resolveDIAOwnedValue(reference(), manager))
        previous <- shiny::isolate(reference())
        next_reference <- saveDIAOwnedValue(value, label, manager)
        if (inherits(manager, "DIAStateAdapter") && !is.null(snapshot_field) &&
            !identical(manager$.context$snapshot[[snapshot_field]], next_reference)) {
            snapshot <- manager$.context$snapshot
            snapshot[[snapshot_field]] <- next_reference
            publishDIAStateHead(manager$.context, snapshot)
        }
        reference(next_reference)
        releaseDIAReplacedOwner(manager, previous, next_reference)
        invisible(NULL)
    }
    # Internal reference access avoids hydration for metadata-only consumers.
    attr(result, "dia_reference") <- reference
    result
}

#' Omic-agnostic name of the saved reactive value (KIT-011): a reactive holding a saved reference
#' whose getter resolves the payload through the active state owner. `newDIASavedReactiveVal` stays
#' as the DIA module's name for it.
#' @noRd
savedReactiveVal <- newDIASavedReactiveVal

#' Read an internal reactive reference without hydrating its payload.
#' @param value Existing getter with optional private reference accessor.
#' @return Saved descriptor or unchanged native getter result.
#' @noRd
getDIAReactiveReference <- function(value) {
    reference <- attr(value, "dia_reference", exact = TRUE)
    if (is.function(reference)) reference() else value()
}

#' Attach one adapter lifetime to existing Shiny flush and project boundaries.
#' @param manager Adapter containing descriptors, never a resident history.
#' @param session,workflow_data Existing session and reactive data bus.
#' @param receipt_dir Optional existing harness-only receipt destination.
#' @param scope Existing diagnostic scope label.
#' @return Attached manager.
#' @noRd
attachDIAStateSession <- function(manager, session, workflow_data, receipt_dir = NULL,
    scope = "dia-native-session") {
    manager$.context$on_commit <- function() {
        workflow_data$state_update_trigger <- Sys.time()
    }
    session$userData$dia_capacity_manager <- manager
    manager$.context$store$progress <- if (is.null(receipt_dir)) NULL else
        newWorkflowProgressObserver(receipt_dir, function() manager$.context$snapshot$current_state)
    owner_events <- 0L
    manager$.context$on_owner_boundary <- if (is.null(receipt_dir)) NULL else function(stage) {
        if (owner_events >= 4096L) return(invisible(NULL))
        owner_events <<- owner_events + 1L
        context <- manager$.context
        info <- list(stage = stage, time = as.numeric(Sys.time()),
            current_state = context$snapshot$current_state, scope_depth = context$scope_depth,
            working_records = length(ls(context$working, all.names = TRUE)),
            writes = context$writes, hydrates = context$hydrates,
            open_resources = context$store$open_resources,
            truncated = owner_events == 4096L)
        connection <- file(file.path(receipt_dir, "gui-owner-phases.jsonl"), open = "a")
        on.exit(close(connection), add = TRUE)
        writeLines(jsonlite::toJSON(info, auto_unbox = TRUE, digits = NA), connection)
    }
    write_receipt <- function(stage) {
        if (is.null(receipt_dir)) return(invisible(NULL))
        info <- diaStateOwnership(manager)
        info$stage <- stage
        info$current_state <- manager$current_state
        info$state_names <- names(manager$.context$snapshot$states)
        info$gui_assay_holders <- stats::setNames(lapply(c("data_tbl", "data_cln",
            "ruv_normalised_for_da_analysis_obj"), function(field) {
            value <- shiny::isolate(workflow_data[[field]])
            list(reference = isSavedWorkflowReference(value),
                resident_assay = is.data.frame(value) || isS4(value))
        }), c("data_tbl", "data_cln", "ruv_normalised_for_da_analysis_obj"))
        info$reclamation <- list(enabled = manager$.context$reclaim_enabled,
            pending = manager$.context$reclaim_pending, last = manager$.context$reclaim_result)
        info$scope <- scope
        info$store_directory <- manager$.context$store$directory
        info$time <- as.numeric(Sys.time())
        writeLines(jsonlite::toJSON(info, auto_unbox = TRUE, null = "null", digits = NA),
            file.path(receipt_dir, "gui-state-latest.json"))
        connection <- file(file.path(receipt_dir, "gui-state-phases.jsonl"), open = "a")
        on.exit(close(connection), add = TRUE)
        writeLines(jsonlite::toJSON(info, auto_unbox = TRUE, null = "null", digits = NA), connection)
    }
    # A flush is a delivery boundary, not an ownership scope for heavy payloads.
    callbacks <- list()
    callbacks$before <- session$onFlush(function() {
        observeDIAOwnerBoundary(manager$.context, "flush_begin")
    }, once = FALSE)
    callbacks$after <- session$onFlushed(function() {
        if (!identical(session$userData$dia_capacity_manager, manager) ||
            isTRUE(manager$.context$store$closed)) return(invisible(NULL))
        context <- manager$.context
        if (context$scope_depth != 0L || context$store$open_resources != 0L ||
            !is.null(context$store$table_action)) return(invisible(NULL))
        context$working <- new.env(parent = emptyenv())
        reclaimDIACompletedMemory(manager)
        releaseDIAInputFiles(manager)
        write_receipt("outputs_delivered")
    }, once = FALSE)
    callbacks$end <- session$onSessionEnded(function() {
        write_receipt("session_end")
        closeDIAStateAdapter(manager)
        session$userData$dia_capacity_manager <- NULL
    })
    manager$.context$callbacks <- callbacks
    write_receipt("initialized")
    manager
}

#' Atomically accept validated external session history without retaining its payloads.
#' @param manager Experimental state owner.
#' @param records Complete ordinary legacy history, already externally validated.
#' @param history,current Original active topology.
#' @return NULL invisibly; save failure leaves the previous accepted head intact.
#' @noRd
restoreDIAStateSession <- function(manager, records, history, current) {
    if (!is.list(records) || !is.character(current) || length(current) != 1L ||
        !current %in% names(records) || !all(unlist(history) %in% names(records))) {
        cli::cli_abort("invalid restored DIA topology", class = "saved_entity_invalid")
    }
    replaceDIAStates(manager$.context, records,
        list(state_history = history, current_state = current))
    invisible(NULL)
}

#' Run confidence predicates and global evidence over an accepted saved assay (GUI entry; contract unchanged).
#' KIT-014: the step is the registered declaration `peptide_qvalue_filter`; the GUI runs its compute-only
#' entry point in process with the registry's production variant and commits through the shared finishing step.
#' @param workflow Existing GUI owner.
#' @param qvalue,global_qvalue,proteotypic,global_pg,mode Existing fixed confidence contract.
#' @param now_fn Existing clock seam.
#' @return References and bounded progress summaries, with all database resources released.
#' @noRd
runDIAReferenceQvalue <- function(workflow, qvalue, global_qvalue, proteotypic,
    global_pg, mode, now_fn) {
    manager <- workflow$state_manager
    if (isTRUE(manager$audit_enabled)) {
        cli::cli_abort("reference confidence filtering requires the declared audit-off route",
            class = "saved_entity_unsupported")
    }
    before <- readDIAStateMetadata(manager)
    pending <- resumeDIAQvalueOutput(manager, before)
    if (!is.null(pending)) {
        workflow$qc_params <- pending$qcParams
        return(pending)
    }
    entry <- workflowActionSpec("peptide_qvalue_filter")
    parameters <- list(qvalue_threshold = qvalue, global_qvalue_threshold = global_qvalue,
        global_pg_qvalue_threshold = global_pg, proteotypic_only = isTRUE(proteotypic),
        confidence_provenance_mode = mode)
    validateWorkflowActionParameters(entry$name, parameters, workflowExecutionPolicy(FALSE), resolved = FALSE)
    value <- computeDIAReferenceQvalue(manager, getDIAStateReference(manager), parameters,
        variant = entry$default_variant, now_fn = now_fn)
    finishDIAReferenceQvalueStep(workflow, value, now_fn = now_fn)
}

#' Compute the Q-value and proteotypic filter on a saved state without committing it (KIT-014).
#'
#' The unchanged kernel supplies parameter resolution, validation and the output schema on a zero-row
#' template; the real population runs the exact declared rules: one pass over whole protein groups builds
#' the bounded per-protein identification evidence lookup (the declaration's global statistic), then the
#' rows stream in source order with their global offsets through `produceDIAConfidenceTables` (the
#' `streamed` variant, the reference path as it stood before this ticket).
#' @param manager,reference Owner and the saved input state.
#' @param parameters `qvalue_threshold`, `global_qvalue_threshold`, `global_pg_qvalue_threshold`,
#'   `proteotypic_only` and `confidence_provenance_mode`.
#' @param worker_cap Admitted worker memory in bytes (recorded; the step streams pages and keeps its
#'   evidence lookup under its own 64 MiB bound).
#' @param variant `streamed`.
#' @param now_fn Clock seam for the QC parameter record.
#' @return Compute-only contract: `checkpoint` (components, support, args and dimensions of the output
#'   state, no product), `reference`, `summary` (the filtering progress), `exclusion_summary`,
#'   `dependencies`, `counts`, `execution`, `state_name`, `description`, `config`, `section`, `text`,
#'   `qc` (the step's QC parameter record) and `params` (the resolved kernel parameters).
#' @noRd
diaQvalueWorkingBytes <- 2048

computeDIAReferenceQvalue <- function(manager, reference, parameters, worker_cap = 12 * 1024^3,
    variant = "streamed", now_fn = Sys.time) {
    if (!is.character(variant) || length(variant) != 1L || !identical(variant, "streamed"))
        cli::cli_abort("unsupported Q-value variant", class = "workflow_execution_unsupported")
    if (!is.numeric(worker_cap) || length(worker_cap) != 1L || !is.finite(worker_cap) || worker_cap <= 0)
        cli::cli_abort("invalid admitted worker capacity", class = "workflow_execution_unsupported")
    context <- manager$.context
    before <- readDIAStateMetadata(manager, reference)
    checkpoint <- before$data
    checkpoint$dependencies <- collectDIANativeReferences(list(reference,
        checkpoint$dependencies, checkpoint$components, checkpoint$support, checkpoint$products))
    source <- checkpoint$components$peptide_data
    scope <- newWorkflowTableReadScope(context$store, source, reader = "arrow")
    on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
    columns <- scope$metadata$columns
    empty <- readSavedWorkflowTablePrototype(scope)
    design <- materializeDIATable(context, checkpoint$components$design_matrix)
    slots <- checkpoint$slots
    slots$peptide_data <- empty
    slots$design_matrix <- design
    prototype <- do.call(methods::new, c(list(Class = checkpoint$class), slots))
    local_state <- WorkflowState$new(audit_enabled = FALSE)
    local_state$saveState("raw_data_s4", prototype, before$config, before$description)
    local_workflow <- list2env(list(state_manager = local_state, qc_params = list(peptide_qc = list())),
        parent = emptyenv())
    # The configuration the template keeps in step with the object's args: the session's bus in the GUI,
    # the global copy in a script that has one, else (a worker process or a bare script) a private copy
    # seeded from the saved state's configuration (KIT-009).
    config_env <- workflowConfigEnvironment()
    if (!exists("config_list", envir = config_env, inherits = FALSE)) {
        config_env <- new.env(parent = emptyenv())
        config_env$config_list <- before$config %||% list()
    }
    update_config <- function(...) updateConfigParameter(..., env = config_env)
    # The existing R method supplies parameter resolution, validation and output
    # schema on zero rows. Real populations below use the exact declared rules.
    template <- runPeptideQvalueApplyStep(local_workflow, parameters$qvalue_threshold,
        parameters$global_qvalue_threshold, isTRUE(parameters$proteotypic_only),
        updateConfigParameterFn = update_config, nowFn = now_fn,
        globalPGQvalueThreshold = parameters$global_pg_qvalue_threshold,
        confidenceProvenanceMode = parameters$confidence_provenance_mode)
    result_object <- template$filteredS4
    params <- result_object@args$srlQvalueProteotypicPeptideClean
    protein <- checkpoint$slots$protein_id_column
    peptide <- checkpoint$slots$peptide_sequence_column
    confidence <- unlist(params$confidence_metric_columns, use.names = FALSE)
    retained <- unique(c(params$input_matrix_column_ids, protein, peptide, "Modified.Sequence"))
    retained <- unique(c(retained[retained %in% names(empty) | retained == slots$norm_quantity_column],
        confidence[confidence %in% names(empty)]))
    # KIT-022 (KIT-014 audit PH-02): the admission estimate the family has, before any work: one streamed page
    # (10,000 rows of the retained columns, or fewer within 64 MiB) plus the declared working bytes per row, beside
    # the evidence lookup's own 64 MiB bound; refused above the cap, recorded with the execution.
    row_bytes <- workflowProjectedRowBytes(scope$metadata, retained)
    batch_rows <- min(10000, floor(67108864 / max(1, row_bytes)))
    estimate <- 67108864 + batch_rows * (row_bytes + diaQvalueWorkingBytes)
    if (estimate > worker_cap) cli::cli_abort(
        "Q-value working-memory estimate exceeds admitted worker capacity", class = "workflow_execution_unsupported")
    admission <- list(estimated_bytes = estimate, worker_cap = worker_cap, batch_rows = batch_rows, row_bytes = row_bytes,
        working_bytes_per_row = diaQvalueWorkingBytes, overhead_bytes = 67108864, policy = "qvalue-streamed-page-v1",
        limitation = "admission estimate; kernel cgroup guard remains authoritative")
    produced <- produceDIAConfidenceTables(scope, retained, params, slots,
        result_object@peptide_data, params$biological_exclusion_ledger, context$software)
    result <- produced$reference
    ledger <- produced$ledger
    summary <- produced$summary
    classification_status <- summary$classification_status
    releaseWorkflowTableReadScope(scope)
    accepted_scope <- newWorkflowTableReadScope(context$store, result, reader = "arrow")
    on.exit(releaseWorkflowTableReadScope(accepted_scope), add = TRUE)
    counted <- summarizeSavedDIAFiltering(accepted_scope, prototype, protein, produced$counts)
    checkpoint$dependencies <- collectDIANativeReferences(list(checkpoint$dependencies,
        produced$dependencies, counted$dependencies))
    sample_ids <- counted$runs
    probe <- empty[rep(NA_integer_, length(sample_ids)), , drop = FALSE]
    probe[[slots$sample_id]] <- sample_ids
    prototype@peptide_data <- probe
    cleaned_design <- cleanDesignMatrixPeptide(prototype)@design_matrix
    progress <- counted$summary
    proteins <- counted$proteins
    releaseWorkflowTableReadScope(accepted_scope)
    checkpoint$slots$args <- result_object@args
    checkpoint$slots$args$srlQvalueProteotypicPeptideClean$biological_exclusion_summary <- summary
    checkpoint$slots$args$srlQvalueProteotypicPeptideClean["biological_exclusion_ledger"] <- list(NULL)
    checkpoint$support$srlQvalueProteotypicPeptideClean <- list(biological_exclusion_ledger = ledger)
    checkpoint$components$peptide_data <- result
    checkpoint$components$design_matrix <- saveDIATableComponent(context, cleaned_design)
    checkpoint$dimensions <- list(rows = summary$biological_rows, columns = names(result_object@peptide_data))
    config <- local_state$states$qvalue_filtered$config
    text <- sub("Proteins remaining: 0", paste("Proteins remaining:", proteins), template$resultText, fixed = TRUE)
    text <- sub("Decoy/contaminant rows excluded: 0", paste("Decoy/contaminant rows excluded:",
        summary$excluded_rows), text, fixed = TRUE)
    text <- sub(params$biological_exclusion_summary$classification_status,
        classification_status, text, fixed = TRUE)
    list(checkpoint = checkpoint, reference = result, summary = progress, exclusion_summary = summary,
        dependencies = checkpoint$dependencies, counts = produced$counts,
        execution = list(variant = variant, worker_cap = worker_cap, rows = summary$biological_rows,
            excluded_rows = summary$excluded_rows, classification_status = classification_status, admission = admission),
        state_name = "qvalue_filtered", description = "Applied Q-value and proteotypic peptide filter",
        config = config, section = "srlQvalueProteotypicPeptideClean", text = text,
        qc = local_workflow$qc_params$peptide_qc$qvalue_filter, params = params, proteins = proteins)
}

#' Commit a computed Q-value result with its bounded progress product (KIT-014): the GUI in process and
#' the controller's completion of a worker result share this step.
#' @param workflow Session owner (`state_manager`, `qc_params`, optional `filtering_progress`).
#' @param value Compute-only contract from `computeDIAReferenceQvalue`.
#' @param now_fn Existing clock seam for the commit timestamp.
#' @return Display result containing references and summaries only.
#' @noRd
finishDIAReferenceQvalueStep <- function(workflow, value, now_fn = Sys.time) {
    manager <- workflow$state_manager
    context <- manager$.context
    setDIAStepProvenanceForSection(manager, value$section) # KIT-025: the in-process path names its operation
    checkpoint <- value$checkpoint
    qc_params <- shiny::isolate(workflow$qc_params) %||% list()
    if (is.null(qc_params$peptide_qc)) qc_params$peptide_qc <- list()
    qc_params$peptide_qc$qvalue_filter <- value$qc
    workflow$qc_params <- qc_params
    prior <- shiny::isolate(workflow$filtering_progress)
    prior_slots <- if (methods::is(prior, "FilteringProgress")) stats::setNames(
        lapply(methods::slotNames(prior), \(name) methods::slot(prior, name)), methods::slotNames(prior)) else NULL
    product <- list(filteringSummary = value$summary, resultText = value$text, priorProgress = prior_slots,
        qcParams = qc_params)
    if (as.numeric(object.size(product)) > 67108864) {
        cli::cli_abort("Q-value progress exceeds product limits", class = "saved_entity_limit")
    }
    identity <- newSavedWorkflowScientificIdentity(checkpoint$components$peptide_data$entity_id,
        checkpoint$components$design_matrix$entity_id, digest::digest(value$params, algo = "sha256"),
        "qvalue-progress", context$software)
    checkpoint$products$qvalue_progress <- saveWorkflowConsumerEntity(context$store,
        newSavedWorkflowScientificProduct(product, "dia-qvalue-progress"),
        basename(tempfile("dia-qvalue-progress-")), identity, "scientific_product")
    reference <- commitDIAReferenceState(manager, value$state_name, checkpoint, value$config,
        value$description, timestamp = now_fn())
    c(list(stateReference = reference), product)
}

#' The Q-value and proteotypic filter as a declaration on the feature-group unit (KIT-014): the unchanged
#' `srlQvalueProteotypicPeptideClean` kernel resolves the parameters and the output schema; the population
#' runs the exact declared rules with the per-protein identification evidence computed once over whole
#' protein groups before the rows stream in source order (both global statistics, step map `qvalue`).
#' The kernel runs the compute-only entry point.
#' @return Declaration from `newWorkflowOperation`.
#' @noRd
diaQvalueFilterOperation <- function() {
    newWorkflowOperation(name = "peptide_qvalue_filter", input_class = "PeptideQuantitativeData",
        fields = "design_matrix", unit = "feature_group",
        kernel = function(object, parameters, store) computeDIAReferenceQvalue(store$manager, store$source,
            parameters[c("qvalue_threshold", "global_qvalue_threshold", "global_pg_qvalue_threshold",
                "proteotypic_only", "confidence_provenance_mode")],
            worker_cap = parameters$worker_cap %||% 12 * 1024^3, variant = parameters$variant %||% "streamed"),
        output_class = "PeptideQuantitativeData",
        parameters = c("qvalue_threshold", "global_qvalue_threshold", "global_pg_qvalue_threshold",
            "proteotypic_only", "confidence_provenance_mode"),
        product_type = "dia.qvalue_filter_result",
        global_statistics = c(protein_identification_evidence =
            "distinct peptide and peptidoform counts per protein over all runs, computed once over whole protein groups and joined to every accepted row",
            source_row_offsets = ".source_row_id is a whole-table row ordinal; the rows stream in source order with their global offsets"),
        merge = c(design_matrix = "design", peptide_data = "concatenate", biological_exclusion_ledger = "concatenate",
            counts = "sum"))
}


#' Mark a newly completed expensive publication for one post-delivery release.
#' @param manager,reference Current adapter and accepted scientific reference.
#' @param request_id,label Existing durable output request identity and group.
#' @return TRUE when marked; repeated completed-journal reads must not call this.
#' @noRd
requestDIAReclamation <- function(manager, reference, request_id, label) {
    context <- manager$.context
    if (!isTRUE(context$reclaim_enabled) || isTRUE(context$store$closed) ||
        !label %in% c("import", "design", "qc", "normalization")) return(invisible(FALSE))
    checkSavedEntityString(request_id, "reclamation request")
    token <- list(incarnation = context$store$identity$incarnation,
        state = reference$entity_id, request = request_id, label = label)
    if (identical(token, context$reclaim_last)) return(invisible(FALSE))
    context$reclaim_pending <- token
    invisible(TRUE)
}

#' Consume a completed-action request only when actual resource owners are idle.
#' @param manager Existing adapter; this helper never mutates scientific heads.
#' @return Bounded status; failures preserve scientific/output completion.
#' @noRd
reclaimDIACompletedMemory <- function(manager) {
    context <- manager$.context
    pending <- context$reclaim_pending
    if (!isTRUE(context$reclaim_enabled) || isTRUE(context$store$closed) || is.null(pending))
        return(invisible(list(status = "not_requested")))
    current <- if (identical(pending$label, "import")) context$snapshot$gui_references$data_tbl else
        context$snapshot$states[[context$snapshot$current_state]]
    if (!identical(pending$incarnation, context$store$identity$incarnation) ||
        !identical(pending$state, current$entity_id)) {
        context$reclaim_pending <- NULL
        return(invisible(list(status = "superseded")))
    }
    if (context$scope_depth != 0L || context$store$open_resources != 0L ||
        !is.null(context$store$table_action) || length(ls(context$working, all.names = TRUE)) ||
        length(context$owned_workers) || length(context$borrowed_workers)) {
        return(invisible(list(status = "busy")))
    }
    context$reclaim_pending <- NULL
    context$reclaim_last <- pending
    result <- tryCatch(reclaimWorkflowMemory(), workflow_reclaim_error = function(error) {
        logger::log_warn(paste("Optional memory reclamation failed:", conditionMessage(error)))
        list(status = "failed")
    })
    context$reclaim_result <- c(pending, result)
    invisible(result)
}

#' Record the analysis's external inputs in the head (KIT-026): each file by role with its name, bytes and
#' sha256, each saved input by role with its entity id. These entries carry no `source`, so a session export
#' records them without copying the file (`resolveDIANativeDependencies`).
#' @param manager Current DIA owner; other owners are unchanged.
#' @param paths Named character vector, role to an existing file.
#' @param entities Named list, role to an entity id.
#' @return NULL invisibly.
#' @noRd
recordDIAWorkflowInputs <- function(manager, paths = character(), entities = list()) {
    if (!usesDIAReferenceStates(manager)) return(invisible(NULL))
    snapshot <- manager$.context$snapshot
    for (role in names(paths)) {
        path <- paths[[role]]
        snapshot$file_dependencies[[role]] <- list(role = role, name = basename(path),
            bytes = as.numeric(workflowTableFileIdentity(path)$size), sha256 = hashWorkflowFile(path))
    }
    for (role in names(entities)) snapshot$file_dependencies[[role]] <- list(role = role, entity_id = entities[[role]])
    publishDIAStateHead(manager$.context, snapshot)
    invisible(NULL)
}

#' Record one external service the analysis queried (KIT-026), appended to `snapshot$external_resources`.
#' A worker never publishes the head, so a query made inside a worker is not recorded.
#' @param manager Current DIA owner; other owners are unchanged.
#' @param service Service name (`gprofiler2`, `clusterProfiler`, `string-db`, `uniprot`).
#' @param release Release the service reported or the client queried; NA when unknown.
#' @param release_kind `database_version`, `package_version`, `api_version` or `header`.
#' @param query_date Time of the query, recorded as ISO 8601.
#' @param query_summary Short description of the query.
#' @return NULL invisibly.
#' @noRd
recordWorkflowExternalResource <- function(manager, service, release, release_kind, query_date = Sys.time(),
    query_summary = NA_character_) {
    if (!usesDIAReferenceStates(manager) || identical(Sys.getenv("MULTISCHOLAR_WORKFLOW_WORKER"), "1")) return(invisible(NULL))
    snapshot <- manager$.context$snapshot
    snapshot$external_resources[[length(snapshot$external_resources) + 1L]] <- list(service = service,
        release = if (length(release)) as.character(release[[1L]]) else NA_character_, release_kind = release_kind,
        query_date = format(query_date, "%Y-%m-%dT%H:%M:%OS3%z"), query_summary = query_summary)
    publishDIAStateHead(manager$.context, snapshot)
    invisible(NULL)
}

#' Register immutable source-file identities for one managed import.
#' @param manager Current adapter; ordinary workflows remain unchanged.
#' @param paths Explicit report/FASTA paths consumed by this import.
#' @return NULL invisibly; no file contents are read or copied.
#' @noRd
registerDIAInputFiles <- function(manager, paths) {
    if (!usesDIAReferenceStates(manager)) return(invisible(NULL))
    # KIT-026: the caller names each path by its role (`diann_report`, `fasta`); the names survive to completion.
    paths <- paths[file.exists(paths) & !duplicated(paths)]
    if (length(paths) > 16L) cli::cli_abort("import file leases exceed limits", class = "saved_entity_limit")
    manager$.context$input_file_leases <- lapply(paths, workflowTableFileIdentity)
    manager$.context$input_files_ready <- FALSE
    invisible(NULL)
}

#' Mark consumed input files after successful parsing and durable table handoff.
#' @param manager Current managed import owner.
#' @return NULL invisibly; actual release waits for idle delivery.
#' @noRd
completeDIAInputFiles <- function(manager) {
    if (!usesDIAReferenceStates(manager)) return(invisible(NULL))
    for (identity in manager$.context$input_file_leases) {
        if (!identical(identity, workflowTableFileIdentity(identity$path)))
            cli::cli_abort("source file changed during import", class = "saved_entity_corrupt")
    }
    manager$.context$input_files_ready <- TRUE
    leases <- manager$.context$input_file_leases
    roles <- names(leases) %||% rep("", length(leases))
    if (any(nzchar(roles))) recordDIAWorkflowInputs(manager,
        vapply(leases[nzchar(roles)], `[[`, character(1), "path"))
    reference <- manager$.context$snapshot$gui_references$data_tbl
    if (!is.null(reference)) requestDIAReclamation(manager, reference,
        paste0("import-", reference$entity_id), "import")
    invisible(NULL)
}

#' Release consumed source-file cache after parser/scientific scopes have returned.
#' @param manager Current adapter, checked by its session delivery callback.
#' @return NULL invisibly; originals and accepted checkpoints remain unchanged.
#' @noRd
releaseDIAInputFiles <- function(manager) {
    context <- manager$.context
    if (!isTRUE(context$input_files_ready) || context$scope_depth != 0L ||
        context$store$open_resources != 0L || !is.null(context$store$table_action)) return(invisible(NULL))
    leases <- context$input_file_leases
    context$input_files_ready <- FALSE
    context$input_file_leases <- list()
    for (identity in leases) {
        if (file.exists(identity$path) && identical(identity, workflowTableFileIdentity(identity$path)))
            releaseWorkflowFilePages(identity$path)
    }
    invisible(NULL)
}

#' Apply native confidence predicates and exclusions to one bounded projected batch.
#' @param data,positions Scientific rows and their current source ordinals.
#' @param retained Existing output columns; count-only projections may omit quantities.
#' @param params,slots Resolved native confidence parameters and column roles.
#' @param manifest Already parsed native manifest with unchanged provenance.
#' @return Native exclusion result for confidence-valid rows.
#' @noRd
classifyDIAConfidenceBatch <- function(data, positions, retained, params, slots, manifest) {
    if (!slots$norm_quantity_column %in% names(data) && slots$raw_quantity_column %in% names(data))
        data[[slots$norm_quantity_column]] <- data[[slots$raw_quantity_column]]
    mask <- rep(TRUE, nrow(data))
    for (metric in names(params$confidence_metric_columns)) {
        column <- params$confidence_metric_columns[[metric]]
        if (column %in% names(data)) mask <- mask & if (metric == "proteotypic")
            data[[column]] == 1 else data[[column]] <= 0.01
    }
    data <- dplyr::mutate(data, .source_row_id = .env$positions)
    data <- dplyr::filter(data, .env$mask)
    data <- dplyr::select(data, dplyr::all_of(unique(c(".source_row_id", intersect(retained, names(data))))))
    classifyPeptideBiologicalExclusionsResolved(data, slots$protein_id_column, "Protein.Ids", manifest,
        params$exclude_decoys, params$exclude_contaminants)
}

#' Produce first-QC tables with native complete-protein evidence and ordered streaming.
#' @param scope Current assay lease.
#' @param retained Existing scientific output population/column selection.
#' @param params,slots Resolved native parameters and roles.
#' @param result_prototype,ledger_prototype Native zero-row output schemas.
#' @param software Current numerical software identity.
#' @return Durable accepted/ledger references, exact summary and index dependencies.
#' @noRd
produceDIAConfidenceTables <- function(scope, retained, params, slots,
    result_prototype, ledger_prototype, software) {
    manifest <- params$contaminant_manifest_provenance
    protein <- slots$protein_id_column
    peptide <- slots$peptide_sequence_column
    flags <- .findPeptideExclusionColumns(retained, c("Decoy", "Reverse", "Is.Decoy", "IsDecoy",
        "Contaminant", "Potential.contaminant", "Potential.Contaminant", "Is.Contaminant", "IsContaminant"))
    columns <- intersect(scope$metadata$columns, unique(c(protein, peptide, "Modified.Sequence",
        "Protein.Ids", flags, unlist(params$confidence_metric_columns, use.names = FALSE))))
    cursor <- openWorkflowTableGroups(scope, columns, protein)
    on.exit(cursor$close(), add = TRUE)
    dependencies <- cursor$dependencies
    evidence <- list()
    summary <- NULL
    count_fields <- c("input_rows", "decoy_rows", "contaminant_rows", "partial_manifest_group_rows",
        "excluded_rows", "biological_rows")
    repeat {
        part <- cursor$next_group()
        if (is.null(part)) break
        positions <- attr(part, "workflow_source_order")
        attr(part, "workflow_source_order") <- NULL
        classified <- classifyDIAConfidenceBatch(part, positions, retained, params, slots, manifest)
        if (is.null(summary)) summary <- classified$summary else {
            for (field in count_fields) summary[[field]] <- summary[[field]] + classified$summary[[field]]
            if (classified$summary$classification_status == "classified_from_available_signals")
                summary$classification_status <- "classified_from_available_signals"
        }
        counted <- .calculateProteinIdentificationEvidence(classified$analysis_data,
            protein_id_column = protein, peptide_sequence_column = peptide,
            modified_peptide_sequence_column = "Modified.Sequence")
        # The lookup's key values are matched back to the untouched source column.
        if (is.factor(counted[[protein]])) counted[[protein]] <- as.character(counted[[protein]])
        evidence[[length(evidence) + 1L]] <- counted
        if (as.numeric(object.size(evidence)) > 67108864)
            cli::cli_abort("identification evidence exceeds lookup limits", class = "saved_entity_limit")
        part <- classified <- counted <- NULL
    }
    cursor$close()
    evidence <- dplyr::bind_rows(evidence)
    if (as.numeric(object.size(evidence)) > 67108864 || anyDuplicated(evidence[[protein]]))
        cli::cli_abort("identification evidence is not a bounded unique lookup", class = "saved_entity_limit")
    output_metadata <- workflowSavedTableMetadata(result_prototype, ".rownames", 3L)
    ledger_metadata <- workflowSavedTableMetadata(ledger_prototype, ".rownames", 3L)
    output_metadata$rows <- summary$biological_rows
    ledger_metadata$rows <- summary$excluded_rows
    if (!inherits(ledger_prototype, "tbl_df")) {
        ledger_metadata$automatic_row_names <- FALSE
        ledger_metadata$row_name_column <- ".multischolar_row_key"
        ledger_metadata$row_index_type <- "integer"
    }
    directory <- workflowTableBatchDirectory(scope)
    accepted_writer <- openWorkflowTableWriter(tempfile("accepted-", directory, fileext = ".parquet"), output_metadata, scope$store$progress)
    assign("confidence_data", accepted_writer, scope$writers)
    ledger_writer <- openWorkflowTableWriter(tempfile("excluded-", directory, fileext = ".parquet"), ledger_metadata, scope$store$progress)
    assign("confidence_ledger", ledger_writer, scope$writers)
    counts <- newDIAFilteringAccumulator(result_prototype, protein)
    qvalid_rows <- 0L
    input_columns <- intersect(unique(c(retained, slots$raw_quantity_column)), scope$metadata$columns)
    streamSavedWorkflowTable(scope, input_columns, function(data, offset, row_names) {
        classified <- classifyDIAConfidenceBatch(data, offset + seq_len(nrow(data)), retained, params, slots, manifest)
        keep <- !classified$classified_data$exclude_from_biological_analysis
        accepted <- classified$analysis_data
        at <- match(accepted[[protein]], evidence[[protein]])
        for (field in .proteinIdentificationEvidenceColumns) accepted[[field]] <- evidence[[field]][at]
        accepted <- completeWorkflowTableBatch(accepted, output_metadata)
        ledger <- classified$exclusion_ledger
        if (!isTRUE(ledger_metadata$automatic_row_names))
            attr(ledger, "row.names") <- as.integer(qvalid_rows + which(!keep))
        ledger <- completeWorkflowTableBatch(ledger, ledger_metadata)
        counts$push(accepted)
        appendWorkflowTableWriter(accepted_writer, accepted)
        appendWorkflowTableWriter(ledger_writer, ledger)
        qvalid_rows <<- qvalid_rows + length(keep)
    })
    if (!identical(readPeptideContaminantManifest(params$contaminant_manifest), manifest))
        cli::cli_abort("contaminant manifest changed during filtering", class = "saved_entity_corrupt")
    identity <- newSavedWorkflowScientificIdentity(scope$reference$entity_id, "confidence-population",
        digest::digest(params, algo = "sha256"), "confidence-stream-v1", software)
    accepted <- finishWorkflowTableBatch(scope, "confidence_data", identity)
    ledger <- finishWorkflowTableBatch(scope, "confidence_ledger", identity)
    list(reference = accepted, ledger = ledger, summary = summary, dependencies = dependencies,
        counts = counts$finish())
}
