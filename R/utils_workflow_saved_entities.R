#' Validate a bounded metadata string.
#' @param value String to validate.
#' @param label Field name for errors.
#' @return Value invisibly.
#' @noRd
checkSavedEntityString <- function(value, label) {
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
        !nzchar(value) || nchar(value, type = "bytes") > 256L) {
        cli::cli_abort("invalid saved-entity {.val {label}}", class = "saved_entity_invalid")
    }
    invisible(value)
}

#' Read and validate the small store identity marker.
#' @param marker Marker filename.
#' @return Validated ownership identity.
#' @noRd
readSavedStoreIdentity <- function(marker) {
    link <- Sys.readlink(marker)
    if (!file.exists(marker) || file.info(marker)$size > 8192 ||
        (!is.na(link) && nzchar(link))) {
        cli::cli_abort("invalid or incomplete store marker", class = "saved_entity_invalid")
    }
    identity <- tryCatch(jsonlite::read_json(marker), error = function(error) {
        cli::cli_abort("unreadable store marker", class = "saved_entity_invalid", parent = error)
    })
    if (!is.list(identity) || !identical(sort(names(identity)),
        sort(c("schema", "project_id", "workflow_id", "incarnation"))) ||
        !identical(identity$schema, 1L)) {
        cli::cli_abort("invalid store identity schema", class = "saved_entity_invalid")
    }
    for (name in c("project_id", "workflow_id", "incarnation")) checkSavedEntityString(identity[[name]], name)
    if (!grepl("^[a-f0-9]{64}$", identity$incarnation)) {
        cli::cli_abort("invalid store incarnation", class = "saved_entity_invalid")
    }
    identity
}

#' Internal mutable metadata store, following the package's R6 ownership pattern.
#' @noRd
WorkflowEntityStore <- R6::R6Class("WorkflowEntityStore", cloneable = FALSE,
    public = list(directory = NULL, identity = NULL, records = NULL, count = 0L,
        metadata_bytes = 0, max_entries = NULL, max_metadata_bytes = NULL, sync_command = NULL,
        cache = NULL, cache_order = character(), cache_bytes = 0, max_cache_bytes = 0,
        cache_hits = 0L, payload_reads = 0L, payload_decodes = 0L,
        table_reads = 0L, table_rows = 0L, table_columns = 0L,
        table_bytes = 0, open_resources = 0L, closed = FALSE,
        table_action = NULL, io_counts = NULL, progress = NULL, execution_owner = NULL))

#' Create a metadata-only store handle in a dedicated project subdirectory.
#' @param root Existing project directory.
#' @param project_id,workflow_id Explicit ownership identifiers.
#' @param max_entries Maximum indexed descriptors.
#' @param max_metadata_bytes Conservative metadata size ceiling, at most 64 MiB.
#' @param max_cache_bytes Optional resident payload ceiling, at most 64 MiB; zero disables it.
#' @return Internal environment handle. No scientific payload is read or written.
#' @noRd
newWorkflowEntityStore <- function(root, project_id, workflow_id,
                                   max_entries = 10000L, max_metadata_bytes = 67108864,
                                   max_cache_bytes = 0) {
    checkSavedEntityString(project_id, "project_id")
    checkSavedEntityString(workflow_id, "workflow_id")
    if (!is.character(root) || length(root) != 1L || is.na(root) || !dir.exists(root) ||
        length(max_entries) != 1L || !is.numeric(max_entries) ||
        !is.finite(max_entries) || max_entries < 1 || max_entries != floor(max_entries) ||
        length(max_metadata_bytes) != 1L || !is.numeric(max_metadata_bytes) ||
        !is.finite(max_metadata_bytes) || max_metadata_bytes < 1 || max_metadata_bytes > 67108864 ||
        length(max_cache_bytes) != 1L || !is.numeric(max_cache_bytes) ||
        !is.finite(max_cache_bytes) || max_cache_bytes < 0 || max_cache_bytes > 67108864) {
        cli::cli_abort("invalid saved-entity store limits or root", class = "saved_entity_invalid")
    }
    root <- normalizePath(root, mustWork = TRUE)
    directory <- file.path(root, ".workflow-entities")
    created <- dir.create(directory, showWarnings = FALSE)
    if (!dir.exists(directory) || !startsWith(normalizePath(directory), paste0(root, "/"))) {
        cli::cli_abort("saved-entity directory escapes project", class = "saved_entity_path")
    }
    marker <- file.path(directory, "store.json")
    if (created) {
        identity <- list(schema = 1L, project_id = project_id, workflow_id = workflow_id,
            incarnation = digest::digest(list(Sys.time(), Sys.getpid(), tempfile()), algo = "sha256"))
        jsonlite::write_json(identity, marker, auto_unbox = TRUE)
    }
    identity <- readSavedStoreIdentity(marker)
    if (!identical(identity$schema, 1L) || !identical(identity$project_id, project_id) ||
        !identical(identity$workflow_id, workflow_id)) {
        cli::cli_abort("saved-entity store ownership differs", class = "saved_entity_scope")
    }
    store <- WorkflowEntityStore$new()
    store$directory <- normalizePath(directory)
    store$identity <- identity
    store$records <- new.env(hash = TRUE, parent = emptyenv())
    store$count <- 0L
    store$metadata_bytes <- 0
    store$max_entries <- max_entries
    store$max_metadata_bytes <- max_metadata_bytes
    store$cache <- new.env(hash = TRUE, parent = emptyenv())
    store$cache_order <- character()
    store$cache_bytes <- 0
    store$max_cache_bytes <- max_cache_bytes
    store$io_counts <- list(payload_hash_bytes = 0, payload_hashes = 0L,
        metadata_checks = 0L, recovery_flushes = 0L, staged_bytes = 0,
        output_bytes = 0, batches_written = 0L, query_streams = 0L,
        sequential_streams = 0L, fast_reads = 0L, page_releases = 0L)
    store
}

#' Check that a handle still belongs to the same project incarnation.
#' @param store Internal store handle.
#' @return TRUE invisibly.
#' @noRd
checkWorkflowEntityStore <- function(store) {
    if (!is.environment(store) || !inherits(store, "WorkflowEntityStore")) {
        cli::cli_abort("invalid saved-entity handle", class = "saved_entity_invalid")
    }
    if (isTRUE(store$closed)) {
        cli::cli_abort("saved-entity handle is closed", class = "saved_entity_scope")
    }
    marker <- file.path(store$directory, "store.json")
    if (!dir.exists(store$directory) || !identical(normalizePath(store$directory), store$directory) ||
        !file.exists(marker)) {
        cli::cli_abort("saved-entity project incarnation changed", class = "saved_entity_scope")
    }
    if (!identical(readSavedStoreIdentity(marker), store$identity)) {
        cli::cli_abort("saved-entity project incarnation changed", class = "saved_entity_scope")
    }
    invisible(TRUE)
}

#' Release only optional resident payloads, leaving durable records and stages intact.
#' @param store Open saved-entity store.
#' @return TRUE invisibly.
#' @noRd
releaseWorkflowEntityCache <- function(store) {
    checkWorkflowEntityStore(store)
    store$cache <- new.env(hash = TRUE, parent = emptyenv())
    store$cache_order <- character()
    store$cache_bytes <- 0
    invisible(TRUE)
}

#' Return bounded ownership and actual payload access counters without decoding.
#' @param store Open saved-entity store.
#' @return Cache and resource accounting; nested sharing remains a conservative bound.
#' @noRd
workflowEntityCacheInfo <- function(store) {
    checkWorkflowEntityStore(store)
    ids <- store$cache_order
    entries <- lapply(ids, \(id) get(id, envir = store$cache, inherits = FALSE))
    addresses <- vapply(entries, \(entry) rlang::obj_address(entry$value), character(1))
    unique_top_level_bytes <- sum(vapply(entries[!duplicated(addresses)],
        \(entry) entry$bytes, numeric(1)))
    list(entries = length(ids), bytes = store$cache_bytes,
        unique_top_level_bytes = unique_top_level_bytes,
        top_level_aliases = sum(duplicated(addresses)),
        max_bytes = store$max_cache_bytes, metadata_entries = store$count,
        metadata_bytes = store$metadata_bytes + (store$table_action$receipt_bytes %||% 0), cache_hits = store$cache_hits,
        payload_reads = store$payload_reads, payload_decodes = store$payload_decodes,
        table_reads = store$table_reads, table_rows = store$table_rows,
        table_columns = store$table_columns, table_bytes = store$table_bytes,
        open_resources = store$open_resources,
        io = store$io_counts,
        accounting_note = "bytes is a conservative admission bound; top-level aliases are deduplicated, nested sharing is not")
}

#' Release this handle's memory and executable-path owner without deleting files.
#' @param store Saved-entity store handle.
#' @return TRUE invisibly; repeated close is harmless.
#' @noRd
closeWorkflowEntityStore <- function(store) {
    if (!is.environment(store) || !inherits(store, "WorkflowEntityStore")) {
        cli::cli_abort("invalid saved-entity handle", class = "saved_entity_invalid")
    }
    if (isTRUE(store$closed)) return(invisible(TRUE))
    if (!is.null(store$execution_owner)) closeWorkflowExecutionController(store$execution_owner)
    if (store$open_resources != 0L) {
        cli::cli_abort("saved-entity store still owns open resources",
            class = "saved_entity_scope")
    }
    store$cache <- new.env(hash = TRUE, parent = emptyenv())
    store$cache_order <- character()
    store$cache_bytes <- 0
    store$records <- new.env(hash = TRUE, parent = emptyenv())
    store$count <- 0L
    store$metadata_bytes <- 0
    store$sync_command <- NULL
    store$progress <- NULL
    if (!is.null(store$table_action)) endWorkflowTableAction(store, store$table_action, force = TRUE)
    store$closed <- TRUE
    invisible(TRUE)
}

#' Detect mutable references before optional resident reuse.
#' @param object Validated saved payload subtree.
#' @param seen Addresses already visited in a cyclic structure.
#' @return TRUE if a caller could change a resident value by reference.
#' @noRd
hasMutableWorkflowEntityOwner <- function(object, seen = NULL) {
    if (is.environment(object) || is.function(object) ||
        typeof(object) %in% c("externalptr", "weakref")) return(TRUE)
    if (is.null(seen)) seen <- new.env(hash = TRUE, parent = emptyenv())
    address <- rlang::obj_address(object)
    if (exists(address, envir = seen, inherits = FALSE)) return(FALSE)
    assign(address, TRUE, envir = seen)
    if (isS4(object)) {
        for (slot in methods::slotNames(object)) {
            if (hasMutableWorkflowEntityOwner(methods::slot(object, slot), seen)) return(TRUE)
        }
    } else if (is.list(object) || is.pairlist(object)) {
        for (value in object) if (hasMutableWorkflowEntityOwner(value, seen)) return(TRUE)
    }
    for (value in attributes(object)) {
        if (hasMutableWorkflowEntityOwner(value, seen)) return(TRUE)
    }
    FALSE
}

#' Admit one already validated payload without retaining an oversized value.
#' @param store Open store.
#' @param entity_id Immutable accepted entity identifier.
#' @param object Fully decoded value.
#' @return TRUE invisibly when resident; FALSE for zero/oversized cache.
#' @noRd
cacheWorkflowEntityPayload <- function(store, entity_id, object) {
    if (store$max_cache_bytes == 0 || hasMutableWorkflowEntityOwner(object)) {
        return(invisible(FALSE))
    }
    bytes <- as.numeric(object.size(object)) + 128
    if (bytes > store$max_cache_bytes) return(invisible(FALSE))
    while (length(store$cache_order) >= store$max_entries ||
        store$cache_bytes + bytes > store$max_cache_bytes) {
        oldest <- store$cache_order[[1L]]
        store$cache_bytes <- store$cache_bytes - get(oldest, envir = store$cache)$bytes
        rm(list = oldest, envir = store$cache)
        store$cache_order <- store$cache_order[-1L]
    }
    assign(entity_id, list(value = object, bytes = bytes), envir = store$cache)
    store$cache_order <- c(store$cache_order, entity_id)
    store$cache_bytes <- store$cache_bytes + bytes
    invisible(TRUE)
}

#' Validate a relative payload location without decoding it.
#' @param store Store handle.
#' @param relative_path Store-relative payload path.
#' @return Path invisibly.
#' @noRd
checkSavedEntityPath <- function(store, relative_path) {
    checkSavedEntityString(relative_path, "relative_path")
    parts <- strsplit(relative_path, "/", fixed = TRUE)[[1L]]
    if (grepl("[\\\\:]", relative_path) || startsWith(relative_path, "/") ||
        endsWith(relative_path, "/") || any(parts %in% c("", ".", ".."))) {
        cli::cli_abort("invalid relative payload path", class = "saved_entity_path")
    }
    path <- file.path(store$directory, relative_path)
    cursor <- store$directory
    for (part in parts) {
        cursor <- file.path(cursor, part)
        link <- Sys.readlink(cursor)
        if (!is.na(link) && nzchar(link) && !file.exists(cursor) && !dir.exists(cursor)) {
            cli::cli_abort("unresolved payload symlink", class = "saved_entity_path")
        }
    }
    parent <- path
    while (!file.exists(parent) && !dir.exists(parent)) parent <- dirname(parent)
    resolved <- normalizePath(parent, mustWork = TRUE)
    if (!identical(resolved, store$directory) && !startsWith(resolved, paste0(store$directory, "/"))) {
        cli::cli_abort("payload path escapes store", class = "saved_entity_path")
    }
    invisible(relative_path)
}

#' Validate exact bounded factor metadata without coercion.
#' @param factors Named list of atomic vectors or factors.
#' @return TRUE invisibly.
#' @noRd
checkSavedWorkflowFactors <- function(factors) {
    if (!is.list(factors) || length(factors) > 32L ||
        (length(factors) && (is.null(names(factors)) || anyNA(names(factors)) ||
            any(!nzchar(names(factors))) || anyDuplicated(names(factors))))) {
        cli::cli_abort("invalid factor metadata fields", class = "saved_entity_invalid")
    }
    for (value in factors) {
        if (!typeof(value) %in% c("logical", "integer", "double", "character") ||
            length(value) > 256L || is.object(value) && !is.factor(value) ||
            !all(names(attributes(value)) %in% c("names", "levels", "class"))) {
            cli::cli_abort("unsupported factor metadata", class = "saved_entity_invalid")
        }
    }
    invisible(TRUE)
}

#' Reject executable/reference objects hidden in metadata or its attributes.
#' @param value Metadata subtree.
#' @param depth Recursion bound.
#' @return TRUE invisibly.
#' @noRd
checkSavedMetadataTree <- function(value, depth = 0L) {
    if (depth > 12L || !(is.list(value) || is.atomic(value)) ||
        isS4(value) || is.object(value) && !is.factor(value) ||
        !all(names(attributes(value)) %in% c("names", "class", "levels", "package"))) {
        cli::cli_abort("unsupported metadata value or attributes", class = "saved_entity_invalid")
    }
    if (is.list(value)) for (item in value) checkSavedMetadataTree(item, depth + 1L)
    for (item in attributes(value)) checkSavedMetadataTree(item, depth + 1L)
    invisible(TRUE)
}

#' Build an immutable metadata descriptor, not a saved payload.
#' @param store Store handle.
#' @param kind checkpoint, scientific_product, da_fit, result_table or rendered_asset.
#' @param identity Named immutable input/design/parameters/selectors/software fingerprints.
#'   Input fingerprints must cover every effective component, including Limpa uncertainty.
#' @param payload Named relative_path, bytes, sha256 and codec metadata.
#' @param metadata Named class, shape and factors summaries; bounded, not assay data.
#' @param parent_id Optional immutable parent entity identifier.
#' @param revision Nonnegative integer generation within this lineage.
#' @return Descriptor; DIA-006 owns payload serialization and durable registration.
#' @noRd
newSavedWorkflowDescriptor <- function(store, kind, identity, payload, metadata,
                                      parent_id = NA_character_, revision = 0L) {
    checkWorkflowEntityStore(store)
    checkSavedEntityString(kind, "kind")
    if (!is.list(identity) || !is.list(payload) || !is.list(metadata) ||
        !kind %in% c("checkpoint", "scientific_product", "da_fit", "result_table",
            "rendered_asset") ||
        !identical(sort(names(identity)), sort(c("input", "design", "parameters", "selectors", "software"))) ||
        !identical(sort(names(payload)), sort(c("relative_path", "bytes", "sha256", "codec"))) ||
        !identical(sort(names(metadata)), sort(c("class", "shape", "factors")))) {
        cli::cli_abort("invalid descriptor schema", class = "saved_entity_invalid")
    }
    for (name in names(identity)) checkSavedEntityString(identity[[name]], name)
    checkSavedEntityPath(store, payload$relative_path)
    if (!is.character(payload$codec) || length(payload$codec) != 1L ||
        !payload$codec %in% c("rds", "qs2", "parquet_duckdb") ||
        identical(payload$codec, "qs2") && !identical(kind, "checkpoint") ||
        !is.numeric(payload$bytes) || length(payload$bytes) != 1L ||
        !is.finite(payload$bytes) || payload$bytes < 0 || payload$bytes != floor(payload$bytes) ||
        !is.character(payload$sha256) || length(payload$sha256) != 1L ||
        is.na(payload$sha256) || !grepl("^[a-f0-9]{64}$", payload$sha256) ||
        !is.character(metadata$class) || !length(metadata$class) || length(metadata$class) > 16L || anyNA(metadata$class) ||
        !is.numeric(metadata$shape) || length(metadata$shape) > 4L || anyNA(metadata$shape) ||
        any(!is.finite(metadata$shape) | metadata$shape < 0 | metadata$shape != floor(metadata$shape)) ||
        !is.list(metadata$factors) || length(metadata$factors) > 32L || (length(metadata$factors) &&
            (is.null(names(metadata$factors)) || anyDuplicated(names(metadata$factors)))) ||
        !is.character(parent_id) || length(parent_id) != 1L ||
        (!is.na(parent_id) && !grepl("^[a-f0-9]{64}$", parent_id)) ||
        !is.integer(revision) || length(revision) != 1L || is.na(revision) || revision < 0L) {
        cli::cli_abort("invalid descriptor types", class = "saved_entity_invalid")
    }
    checkSavedWorkflowFactors(metadata$factors)
    record <- list(schema = 1L, owner = store$identity, kind = kind, identity = identity,
        payload = payload, metadata = metadata, parent_id = parent_id, revision = revision)
    checkSavedMetadataTree(record)
    if (as.numeric(object.size(record)) > min(65536, store$max_metadata_bytes)) {
        cli::cli_abort("descriptor exceeds metadata bound", class = "saved_entity_limit")
    }
    record$entity_id <- digest::digest(record, algo = "sha256")
    record
}

#' Index metadata already accepted by the caller; does not publish or save a payload.
#' @param store Store handle.
#' @param record Descriptor returned by newSavedWorkflowDescriptor.
#' @return Immutable entity id; repeated identical registration is idempotent.
#' @noRd
registerSavedWorkflowMetadata <- function(store, record) {
    if (!is.list(record)) cli::cli_abort("invalid descriptor", class = "saved_entity_invalid")
    expected <- newSavedWorkflowDescriptor(store, record$kind, record$identity, record$payload,
        record$metadata, record$parent_id, record$revision)
    if (!identical(record$owner, store$identity)) {
        cli::cli_abort("descriptor belongs to another store", class = "saved_entity_scope")
    }
    if (!identical(record, expected)) {
        cli::cli_abort("immutable descriptor was modified", class = "saved_entity_conflict")
    }
    if (exists(record$entity_id, envir = store$records, inherits = FALSE)) return(record$entity_id)
    # Conservative admission budget, not a retained-heap performance measurement.
    # Account once on insertion instead of repeatedly traversing the whole index.
    record_bytes <- as.numeric(object.size(record)) + 512
    if (store$count >= store$max_entries ||
        store$metadata_bytes + record_bytes > store$max_metadata_bytes) {
        cli::cli_abort("metadata index is full", class = "saved_entity_limit")
    }
    assign(record$entity_id, record, envir = store$records)
    store$count <- store$count + 1L
    store$metadata_bytes <- store$metadata_bytes + record_bytes
    record$entity_id
}

#' Describe a registered entity without opening its payload.
#' @param store Store handle.
#' @param entity_id Immutable entity identifier.
#' @return Descriptor or NULL for a miss.
#' @noRd
describeSavedWorkflowEntity <- function(store, entity_id) {
    checkWorkflowEntityStore(store)
    checkSavedEntityString(entity_id, "entity_id")
    get0(entity_id, envir = store$records, inherits = FALSE)
}

#' Find exact typed metadata matches without choosing a latest generation.
#' @param store Store handle.
#' @param criteria Named subset of kind, identity, factors, parent_id or revision.
#'   identity/factors compare complete typed values; omitted criteria are unconstrained.
#' @param max_matches Positive bound on returned descriptors; exceeding it errors.
#' @return Explicit miss/hit/multiple status and matching descriptors.
#' @noRd
findSavedWorkflowEntities <- function(store, criteria = list(), max_matches = 100L) {
    checkWorkflowEntityStore(store)
    allowed <- c("kind", "identity", "factors", "parent_id", "revision")
    if (!is.list(criteria) || (length(criteria) && (is.null(names(criteria)) ||
        anyDuplicated(names(criteria)) || !all(names(criteria) %in% allowed))) ||
        !is.numeric(max_matches) || length(max_matches) != 1L || !is.finite(max_matches) ||
        max_matches < 1 || max_matches != floor(max_matches)) {
        cli::cli_abort("invalid lookup criteria or bound", class = "saved_entity_invalid")
    }
    for (name in names(criteria)) {
        value <- criteria[[name]]
        valid <- switch(name,
            kind = is.character(value) && length(value) == 1L && !is.na(value),
            parent_id = is.character(value) && length(value) == 1L,
            revision = is.integer(value) && length(value) == 1L && !is.na(value),
            identity = is.list(value), factors = is.list(value))
        if (!isTRUE(valid)) cli::cli_abort("wrong lookup field type", class = "saved_entity_invalid")
        if (name == "factors") checkSavedWorkflowFactors(value)
        if (name == "identity") {
            if (!identical(sort(names(value)), sort(c("input", "design", "parameters", "selectors", "software")))) {
                cli::cli_abort("invalid identity criteria", class = "saved_entity_invalid")
            }
            for (field in names(value)) checkSavedEntityString(value[[field]], field)
        }
    }
    ids <- ls(envir = store$records, all.names = TRUE, sorted = TRUE)
    matches <- ids[vapply(ids, function(id) {
        record <- get(id, envir = store$records, inherits = FALSE)
        all(vapply(names(criteria), function(name) {
            value <- if (name == "factors") record$metadata$factors else record[[name]]
            identical(value, criteria[[name]])
        }, logical(1)))
    }, logical(1))]
    if (length(matches) > max_matches) cli::cli_abort("lookup exceeds result bound", class = "saved_entity_limit")
    list(status = if (!length(matches)) "miss" else if (length(matches) == 1L) "hit" else "multiple",
        matches = lapply(matches, function(id) get(id, envir = store$records, inherits = FALSE)))
}

#' Construct the complete scientific identity used by saved consumer products.
#' @param input,design,parameters,selectors,software Immutable component fingerprints.
#' @return Identity accepted by the single saved-entity store.
#' @noRd
newSavedWorkflowScientificIdentity <- function(input, design, parameters,
                                               selectors, software) {
    identity <- list(input = input, design = design, parameters = parameters,
        selectors = selectors, software = software)
    for (name in names(identity)) checkSavedEntityString(identity[[name]], name)
    identity
}

#' Construct an asset identity separately from its scientific product identity.
#' @param product_reference Accepted scientific product/fit/table descriptor.
#' @param aesthetics Bounded immutable display settings.
#' @param format Asset format.
#' @param dimensions Integer width and height in pixels.
#' @param renderer Renderer software/version identity.
#' @return Standard identity whose input is the immutable product entity ID.
#' @noRd
newSavedWorkflowRenderIdentity <- function(product_reference, aesthetics, format,
                                           dimensions, renderer) {
    if (!is.list(product_reference) || !is.character(product_reference$entity_id) ||
        length(product_reference$entity_id) != 1L ||
        !grepl("^[a-f0-9]{64}$", product_reference$entity_id) ||
        !product_reference$kind %in% c("scientific_product", "da_fit", "result_table") ||
        !is.list(aesthetics) || as.numeric(object.size(aesthetics)) > 65536 ||
        !is.character(format) || length(format) != 1L || is.na(format) ||
        !format %in% c("png", "pdf", "svg", "html") ||
        !is.integer(dimensions) || length(dimensions) != 2L || anyNA(dimensions) ||
        any(dimensions < 1L) || any(dimensions > 100000L)) {
        cli::cli_abort("invalid rendered-asset identity", class = "saved_entity_invalid")
    }
    checkSavedMetadataTree(aesthetics)
    checkSavedEntityString(renderer, "renderer")
    newSavedWorkflowScientificIdentity(product_reference$entity_id,
        "rendered-asset-v1", digest::digest(aesthetics, algo = "sha256"),
        digest::digest(list(format = format, dimensions = dimensions), algo = "sha256"),
        renderer)
}

#' Wrap an existing scientific method output without changing its values.
#' @param value Existing product value, or NULL when unavailable.
#' @param product_type Stable producer/product name.
#' @param provenance Bounded substantive provenance.
#' @param status available or unavailable.
#' @param reason Required reason for unavailable products.
#' @return Explicit scientific-product envelope.
#' @noRd
newSavedWorkflowScientificProduct <- function(value, product_type,
    provenance = list(), status = "available", reason = NA_character_) {
    checkSavedEntityString(product_type, "product_type")
    if (!is.list(provenance) || as.numeric(object.size(provenance)) > 65536 ||
        !is.character(status) || length(status) != 1L || is.na(status) ||
        !status %in% c("available", "unavailable") ||
        !is.character(reason) || length(reason) != 1L ||
        (status == "available" && (is.null(value) || !is.na(reason))) ||
        (status == "unavailable" && (!is.null(value) || is.na(reason) || !nzchar(reason)))) {
        cli::cli_abort("invalid scientific-product envelope", class = "saved_entity_invalid")
    }
    checkSavedMetadataTree(provenance)
    list(schema = "scientific_product_v1", product_type = product_type,
        status = status, value = value, provenance = provenance, reason = reason)
}

#' Wrap a real limma base/contrast fit with its tested populations.
#' @param fit Existing MArrayLM fit.
#' @param feature_ids,sample_ids,design_columns Ordered identities.
#' @return Explicit DA-fit envelope.
#' @noRd
newSavedWorkflowDaFit <- function(fit, feature_ids, sample_ids, design_columns) {
    valid_ids <- function(value) is.character(value) && length(value) &&
        !anyNA(value) && !anyDuplicated(value) && all(nzchar(value))
    coefficients <- if (inherits(fit, "MArrayLM")) fit$coefficients else NULL
    design <- if (inherits(fit, "MArrayLM")) fit$design else NULL
    if (!inherits(fit, "MArrayLM") || !is.matrix(coefficients) || !is.matrix(design) ||
        !valid_ids(feature_ids) || !valid_ids(sample_ids) || !valid_ids(design_columns) ||
        nrow(coefficients) != length(feature_ids) || nrow(design) != length(sample_ids) ||
        !identical(rownames(coefficients), feature_ids) ||
        !identical(rownames(design), sample_ids) ||
        !identical(colnames(design), design_columns)) {
        cli::cli_abort("invalid DA-fit envelope", class = "saved_entity_invalid")
    }
    list(schema = "da_fit_v1", fit = fit, feature_ids = feature_ids,
        sample_ids = sample_ids, design_columns = design_columns)
}

#' Wrap an ordered DA result table with explicit key and contrast metadata.
#' @param table Existing result data frame.
#' @param feature_id Feature column, or .rownames.
#' @param contrast Exact contrast identity.
#' @param qvalue_column Existing adjusted-q-value column.
#' @return Explicit result-table envelope.
#' @noRd
newSavedWorkflowResultTable <- function(table, feature_id, contrast, qvalue_column) {
    for (value in list(feature_id = feature_id, contrast = contrast,
        qvalue_column = qvalue_column)) checkSavedEntityString(value, "result field")
    row_ids <- if (identical(feature_id, ".rownames")) rownames(table) else
        if (is.data.frame(table) && feature_id %in% names(table)) table[[feature_id]] else NULL
    if (!is.data.frame(table) || !qvalue_column %in% names(table) ||
        !is.character(row_ids) || length(row_ids) != nrow(table) ||
        anyNA(row_ids) || anyDuplicated(row_ids) || any(!nzchar(row_ids))) {
        cli::cli_abort("invalid result-table envelope", class = "saved_entity_invalid")
    }
    list(schema = "result_table_v1", table = table, feature_id = feature_id,
        contrast = contrast, qvalue_column = qvalue_column)
}

#' Wrap available bytes or an explicit unavailable rendered panel.
#' @param bytes Raw self-contained asset bytes; empty only when unavailable.
#' @param media_type MIME type, or NA when unavailable.
#' @param dimensions Integer width/height, or two NA integers when unavailable.
#' @param status available or unavailable.
#' @param reason Required unavailable reason.
#' @return Explicit rendered-asset envelope.
#' @noRd
newSavedWorkflowRenderedAsset <- function(bytes, media_type, dimensions,
    status = "available", reason = NA_character_) {
    if (!is.raw(bytes) || !is.character(media_type) || length(media_type) != 1L ||
        !is.integer(dimensions) || length(dimensions) != 2L ||
        !is.character(status) || length(status) != 1L || is.na(status) ||
        !status %in% c("available", "unavailable") ||
        !is.character(reason) || length(reason) != 1L) {
        cli::cli_abort("invalid rendered-asset envelope", class = "saved_entity_invalid")
    }
    available <- status == "available" && length(bytes) > 0L && !is.na(media_type) &&
        nzchar(media_type) && !anyNA(dimensions) && all(dimensions > 0L) && is.na(reason)
    unavailable <- status == "unavailable" && !length(bytes) && is.na(media_type) &&
        all(is.na(dimensions)) && !is.na(reason) && nzchar(reason)
    if (!available && !unavailable) {
        cli::cli_abort("inconsistent rendered-asset availability", class = "saved_entity_invalid")
    }
    list(schema = "rendered_asset_v1", status = status, bytes = bytes,
        media_type = media_type, dimensions = dimensions, reason = reason)
}

#' Validate one explicit consumer envelope without running a producer.
#' @param kind Consumer entity kind.
#' @param object Candidate envelope.
#' @return TRUE invisibly.
#' @noRd
checkSavedWorkflowConsumerPayload <- function(kind, object) {
    expected <- switch(kind,
        scientific_product = if (is.list(object) && identical(sort(names(object)),
            sort(c("schema", "product_type", "status", "value", "provenance", "reason"))))
            newSavedWorkflowScientificProduct(object$value, object$product_type,
                object$provenance, object$status, object$reason) else NULL,
        da_fit = if (is.list(object) && identical(sort(names(object)),
            sort(c("schema", "fit", "feature_ids", "sample_ids", "design_columns"))))
            newSavedWorkflowDaFit(object$fit, object$feature_ids, object$sample_ids,
                object$design_columns) else NULL,
        result_table = if (is.list(object) && identical(sort(names(object)),
            sort(c("schema", "table", "feature_id", "contrast", "qvalue_column"))))
            newSavedWorkflowResultTable(object$table, object$feature_id,
                object$contrast, object$qvalue_column) else NULL,
        rendered_asset = if (is.list(object) && identical(sort(names(object)),
            sort(c("schema", "status", "bytes", "media_type", "dimensions", "reason"))))
            newSavedWorkflowRenderedAsset(object$bytes, object$media_type,
                object$dimensions, object$status, object$reason) else NULL,
        NULL)
    if (is.null(expected) || !identical(object, expected)) {
        cli::cli_abort("consumer payload differs from its schema", class = "saved_entity_invalid")
    }
    invisible(TRUE)
}

#' Save one validated consumer entity through the existing RDS store.
#' @param store,object,request_id,identity Existing store/save inputs.
#' @param kind scientific_product, da_fit, result_table or rendered_asset.
#' @param factors Bounded lookup factors.
#' @param lineage Optional lineage.
#' @return Accepted immutable descriptor.
#' @noRd
saveWorkflowConsumerEntity <- function(store, object, request_id, identity, kind,
                                       factors = list(), lineage = NULL) {
    if (!kind %in% c("scientific_product", "da_fit", "result_table", "rendered_asset")) {
        cli::cli_abort("unsupported consumer entity kind", class = "saved_entity_invalid")
    }
    checkSavedWorkflowConsumerPayload(kind, object)
    checkSavedWorkflowFactors(factors)
    shape <- switch(kind,
        scientific_product = as.integer(dim(object$value)),
        da_fit = as.integer(dim(object$fit$coefficients)),
        result_table = as.integer(dim(object$table)),
        rendered_asset = as.integer(c(length(object$bytes), length(object$dimensions))))
    saveWorkflowEntity(store, object, request_id, identity,
        metadata = list(class = class(object), shape = shape, factors = factors),
        kind = kind, lineage = lineage)
}

#' Report codec capabilities without opening a payload.
#' @param reference Saved descriptor.
#' @param reader Explicit table reader; default capability contract is unchanged.
#' @return Explicit hydration and selective-read capabilities.
#' @noRd
savedWorkflowEntityCapabilities <- function(reference, reader = "duckdb") {
    if (!is.list(reference) || !is.list(reference$payload) ||
        !reference$payload$codec %in% c("rds", "qs2", "parquet_duckdb")) {
        cli::cli_abort("invalid capability reference", class = "saved_entity_invalid")
    }
    selective <- isWorkflowTableReference(reference)
    if (!is.character(reader) || length(reader) != 1L || is.na(reader) ||
        !reader %in% c("duckdb", "fastduckdb")) {
        cli::cli_abort("unsupported table reader", class = "saved_entity_unsupported")
    }
    if (identical(reader, "fastduckdb")) {
        if (!selective) cli::cli_abort("fast reader requires Parquet", class = "saved_entity_unsupported")
        return(list(codec = reference$payload$codec, reader = reader, full_hydration = FALSE,
            column_projection = TRUE, row_ranges = TRUE, predicate_pushdown = FALSE,
            key_lookup = FALSE, ranked_selection = FALSE, lookup_without_computation = TRUE))
    }
    list(codec = reference$payload$codec, full_hydration = !selective,
        column_projection = selective, predicate_pushdown = selective,
        lookup_without_computation = TRUE)
}

#' Resolve file-specific durable flushing for the qualified Linux storage path.
#' @param store Store handle; caches only the validated executable path.
#' @return Absolute GNU sync command, never a filesystem-wide invocation.
#' @noRd
savedWorkflowSyncCommand <- function(store) {
    checkWorkflowEntityStore(store)
    if (!is.null(store$sync_command)) return(store$sync_command)
    command <- unname(Sys.which("sync"))
    if (!identical(Sys.info()[["sysname"]], "Linux") || !nzchar(command)) {
        cli::cli_abort("durable saved-entity storage requires GNU sync on Linux", class = "saved_entity_unsupported")
    }
    version <- system2(command, "--version", stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(version, "status")) || !any(grepl("GNU coreutils", version, fixed = TRUE))) {
        cli::cli_abort("file-specific GNU sync is unavailable", class = "saved_entity_unsupported")
    }
    store$sync_command <- command
    command
}

#' Flush explicit owned files and directory entries to the filesystem.
#' @param store Store handle.
#' @param paths Nonempty paths within the store, plus its project directory.
#' @return TRUE invisibly; errors prevent successful acknowledgement.
#' @noRd
syncSavedWorkflowPaths <- function(store, paths) {
    command <- savedWorkflowSyncCommand(store)
    if (!is.character(paths) || !length(paths) || anyNA(paths) || any(!nzchar(paths)) ||
        !all(file.exists(paths) | dir.exists(paths))) {
        cli::cli_abort("invalid flush targets", class = "saved_entity_path")
    }
    paths <- unique(normalizePath(paths, mustWork = TRUE))
    allowed <- paths == store$directory | paths == dirname(store$directory) |
        startsWith(paths, paste0(store$directory, "/"))
    if (!all(allowed)) cli::cli_abort("flush target escapes project", class = "saved_entity_path")
    status <- system2(command, c("--", shQuote(paths)))
    if (!identical(status, 0L)) cli::cli_abort("durable flush failed", class = "saved_entity_write")
    invisible(TRUE)
}

#' Named storage boundaries for failure tests; never scientific callbacks.
#' @param phase Storage phase name.
#' @return NULL invisibly.
#' @noRd
savedWorkflowStorageBoundary <- function(phase) invisible(NULL)

#' Write the first checkpoint codec using streaming, uncompressed native RDS.
#' @param object Accepted R object.
#' @param path Staged payload file.
#' @return NULL invisibly.
#' @noRd
writeSavedWorkflowPayload <- function(object, path) {
    base::saveRDS(object, path, compress = FALSE, version = 3)
    invisible(NULL)
}

#' Write the explicitly selected experimental qs2 checkpoint codec.
#' @param object Complete accepted checkpoint record.
#' @param path Fresh staged qs2 payload path.
#' @return NULL invisibly; compression policy is fixed for the matched pilot.
#' @noRd
writeSavedWorkflowQsPayload <- function(object, path) {
    if (!requireNamespace("qs2", quietly = TRUE)) {
        cli::cli_abort("qs2 is required for this checkpoint codec",
            class = "saved_entity_unsupported")
    }
    qs2::qs_save(object, path, compress_level = 3L, shuffle = TRUE, nthreads = 2L)
    invisible(NULL)
}

#' Reject nonportable resources without evaluating active or lazy bindings.
#' @param object Value to store.
#' @param seen Environment-identity set used for cycles.
#' @return TRUE invisibly.
#' @noRd
checkSavedWorkflowPayload <- function(object, seen = NULL) {
    # Empty argument bindings in plotting closures are portable R values, not
    # resources. Inspect the sentinel without forcing it as an ordinary symbol.
    if (rlang::is_missing(object)) return(invisible(TRUE))
    if (is.null(seen)) seen <- new.env(parent = emptyenv())
    attributes_to_check <- attributes(object)
    # Native RDS restores readr's parser marker as an unclassed NULL pointer.
    # Preserve that exact historical attribute; it owns no live resource. A live
    # parser pointer (or a process pointer anywhere else) remains unsupported.
    if (inherits(object, "spec_tbl_df") &&
        identical(attr(object, "problems", exact = TRUE), methods::new("externalptr"))) {
        attributes_to_check$problems <- NULL
    }
    if (typeof(object) %in% c("externalptr", "weakref") || inherits(object, "connection")) {
        cli::cli_abort("payload contains a nonportable resource", class = "saved_entity_unsupported")
    }
    if (is.environment(object)) {
        if (isNamespace(object) || identical(object, globalenv()) ||
            identical(object, baseenv()) || identical(object, emptyenv())) return(invisible(TRUE))
        id <- rlang::obj_address(object)
        if (exists(id, envir = seen, inherits = FALSE)) return(invisible(TRUE))
        assign(id, TRUE, envir = seen)
        fields <- ls(object, all.names = TRUE)
        if (any(vapply(fields, bindingIsActive, logical(1), env = object)) ||
            any(rlang::env_binding_are_lazy(object, fields))) {
            cli::cli_abort("payload contains active or lazy bindings", class = "saved_entity_unsupported")
        }
        # Unlike get(), non-dispatching extraction preserves a missing binding.
        # Active/lazy bindings have already been rejected above.
        for (field in fields) checkSavedWorkflowPayload(.subset2(object, field), seen)
        checkSavedWorkflowPayload(parent.env(object), seen)
    } else if (is.function(object) && !is.primitive(object)) {
        checkSavedWorkflowPayload(environment(object), seen)
    } else if (isS4(object)) {
        methods::validObject(object)
        for (field in methods::slotNames(object)) checkSavedWorkflowPayload(methods::slot(object, field), seen)
        attributes_to_check <- attributes_to_check[setdiff(names(attributes_to_check), methods::slotNames(object))]
    } else if (is.list(object) || is.pairlist(object)) {
        for (value in object) checkSavedWorkflowPayload(value, seen)
    }
    for (value in attributes_to_check) checkSavedWorkflowPayload(value, seen)
    invisible(TRUE)
}

#' Transform only durable configuration locations in a DIA checkpoint record.
#' @param state Existing history record with data and config fields.
#' @param transform Configuration transformation; never traverses assay tables.
#' @return Record preserving all non-runtime fields and scientific data.
#' @noRd
mapDIAStateConfig <- function(state, transform) {
    if (!is.list(state) || is.null(names(state)) || anyDuplicated(names(state))) {
        cli::cli_abort("invalid DIA history record", class = "saved_entity_invalid")
    }
    if (isS4(state$data)) {
        if (!methods::is(state$data, "PeptideQuantitativeData") &&
            !methods::is(state$data, "ProteinQuantitativeData")) {
            cli::cli_abort("unsupported DIA checkpoint class", class = "saved_entity_unsupported")
        }
        # Compatibility adapter: these established classes expose args as a slot.
        args <- methods::slot(state$data, "args")
        if (length(args)) methods::slot(state$data, "args") <- transform(args)
    }
    if (length(state$config)) state$config <- transform(state$config)
    state
}

#' Make readr diagnostics portable without changing values or other attributes.
#' @param table A DIA input, annotation or design table.
#' @return Table with live parser diagnostics materialized; original unchanged.
#' @noRd
prepareDIAParserTable <- function(table) {
    if (!inherits(table, "spec_tbl_df")) return(table)
    marker <- attr(table, "problems", exact = TRUE)
    if (typeof(marker) != "externalptr" || identical(marker, methods::new("externalptr"))) {
        return(table)
    }
    diagnostics <- readr::problems(table)
    attr(table, "problems") <- diagnostics
    table
}

#' Prepare historical DIA records without worker handles or live parser pointers.
#' @param history Complete named list of DIA history records.
#' @return Portable history; original data and borrowed workers remain untouched.
#'   Rebind the selected record with bindDIAWorkflowState before execution.
#' @noRd
prepareDIAWorkflowHistory <- function(history) {
    if (!is.list(history) || is.null(names(history)) || anyDuplicated(names(history))) {
        cli::cli_abort("invalid DIA history", class = "saved_entity_invalid")
    }
    for (i in seq_along(history)) {
        state <- mapDIAStateConfig(history[[i]], \(config) splitDIAWorkerConfig(config)$config)
        if (isS4(state$data)) {
            for (slot in intersect(methods::slotNames(state$data),
                c("peptide_data", "protein_quant_table", "protein_id_table", "design_matrix"))) {
                methods::slot(state$data, slot) <- prepareDIAParserTable(methods::slot(state$data, slot))
            }
        } else if (is.data.frame(state$data)) {
            state$data <- prepareDIAParserTable(state$data)
        }
        history[[i]] <- state
    }
    history
}

#' Reconnect a selected hydrated DIA record to caller-owned live workers.
#' @param state Hydrated DIA history record.
#' @param cluster Live cluster matching the saved execution settings.
#' @return Record usable by existing scientific methods; no computation is run.
#' @noRd
bindDIAWorkflowState <- function(state, cluster) {
    mapDIAStateConfig(state, \(config) bindDIAWorkerConfig(config, cluster))
}

#' Resolve a fixed request location, rejecting path escapes.
#' @param store Store handle.
#' @param request_id Caller-owned idempotency token.
#' @return Store-relative committed directory.
#' @noRd
savedWorkflowRequestPath <- function(store, request_id) {
    checkWorkflowEntityStore(store)
    if (!is.character(request_id) || length(request_id) != 1L || is.na(request_id) ||
        !grepl("^[A-Za-z0-9_-]{1,80}$", request_id)) {
        cli::cli_abort("invalid saved request id", class = "saved_entity_invalid")
    }
    relative <- paste0("entities/", request_id)
    checkSavedEntityPath(store, relative)
    relative
}

#' Borrow the existing store's action-local I/O ownership.
#' @param store Open store; no domain-specific state manager is required.
#' @return Shared action with bounded receipts and one optional database connection.
#' @noRd
beginWorkflowTableAction <- function(store) {
    checkWorkflowEntityStore(store)
    if (is.null(store$table_action)) {
        action <- new.env(parent = emptyenv())
        action$depth <- action$scopes <- 0L
        action$connection <- action$spill <- action$batch_directory <- NULL
        action$validated <- new.env(hash = TRUE, parent = emptyenv())
        action$receipt_bytes <- 0
        store$table_action <- action
    }
    store$table_action$depth <- store$table_action$depth + 1L
    store$table_action
}

#' Release shared I/O ownership after its last borrower.
#' @param store,action Existing owner and its exact action.
#' @param force Store-close cleanup after all resource leases have closed.
#' @return NULL invisibly; no accepted payload is deleted.
#' @noRd
endWorkflowTableAction <- function(store, action, force = FALSE) {
    if (!identical(store$table_action, action)) return(invisible(NULL))
    action$depth <- if (force) 0L else action$depth - 1L
    if (action$depth > 0L) return(invisible(NULL))
    if (action$scopes != 0L) cli::cli_abort("table action still owns leases", class = "saved_entity_scope")
    if (!is.null(action$connection)) DBI::dbDisconnect(action$connection, shutdown = TRUE)
    action$connection <- NULL
    if (!is.null(action$spill)) unlink(action$spill, recursive = TRUE)
    action$spill <- NULL
    if (!is.null(action$batch_directory)) unlink(action$batch_directory, recursive = TRUE)
    action$batch_directory <- NULL
    receipts <- as.list(action$validated, all.names = TRUE)
    paths <- vapply(receipts, function(receipt) receipt$file$path, character(1))
    store$table_action <- NULL
    if (length(paths)) releaseWorkflowStoreFilePages(store, paths)
    invisible(NULL)
}

#' Release clean pages for closed, owned files using the established file-specific mechanism.
#' @param store Open store.
#' @param paths Closed files owned by this store.
#' @return NULL invisibly; cache advice never substitutes for a memory limit.
#' @noRd
releaseWorkflowStoreFilePages <- function(store, paths) {
    checkWorkflowEntityStore(store)
    paths <- unique(paths[file.exists(paths)])
    if (!length(paths)) return(invisible(NULL))
    for (path in paths) {
        if (!startsWith(normalizePath(path), paste0(store$directory, "/"))) {
            cli::cli_abort("cache-release file escapes store", class = "saved_entity_path")
        }
    }
    releaseWorkflowFileCache(paths, file.path(store$directory, "checkpoint-io.log"))
    store$io_counts$page_releases <- store$io_counts$page_releases + length(paths)
    invisible(NULL)
}

#' Remember one verified content hash only for the current action and file identity.
#' @param store Open store.
#' @param record Accepted descriptor.
#' @return NULL invisibly; receipts never survive an action or project replacement.
#' @noRd
rememberWorkflowPayloadValidation <- function(store, record) {
    action <- store$table_action
    if (is.null(action)) return(invisible(NULL))
    receipt <- list(sha256 = record$payload$sha256,
        file = workflowTableFileIdentity(file.path(store$directory, record$payload$relative_path)))
    if (!exists(record$entity_id, action$validated, inherits = FALSE)) {
        bytes <- as.numeric(object.size(receipt)) + 512
        if (length(ls(action$validated)) >= store$max_entries ||
            action$receipt_bytes + store$metadata_bytes + bytes > store$max_metadata_bytes) {
            cli::cli_abort("validation receipts exceed metadata limits", class = "saved_entity_limit")
        }
        action$receipt_bytes <- action$receipt_bytes + bytes
    }
    assign(record$entity_id, receipt, action$validated)
    invisible(NULL)
}

#' Hash a consumed payload once per action, refusing mutation before and after hashing.
#' @param store Open store.
#' @param record Descriptor whose metadata and size have been validated.
#' @return NULL invisibly.
#' @noRd
validateWorkflowPayloadOnce <- function(store, record) {
    path <- file.path(store$directory, record$payload$relative_path)
    before <- workflowTableFileIdentity(path)
    action <- store$table_action
    previous <- if (!is.null(action)) get0(record$entity_id, action$validated, inherits = FALSE) else NULL
    if (!is.null(previous) && identical(previous$file, before) &&
        identical(previous$sha256, record$payload$sha256)) return(invisible(NULL))
    actual <- hashWorkflowOwnedFile(store, path)
    store$io_counts$payload_hash_bytes <- store$io_counts$payload_hash_bytes + record$payload$bytes
    store$io_counts$payload_hashes <- store$io_counts$payload_hashes + 1L
    if (!identical(actual, record$payload$sha256) || !identical(before, workflowTableFileIdentity(path))) {
        cli::cli_abort("payload integrity check failed", class = "saved_entity_corrupt")
    }
    rememberWorkflowPayloadValidation(store, record)
    invisible(NULL)
}

#' Recover only a fully published request, validating metadata before decoding it.
#' @param store Store handle, including a freshly reopened one.
#' @param request_id Original idempotency token.
#' @return Validated descriptor or NULL if no committed request exists.
#' @noRd
recoverSavedWorkflowRequest <- function(store, request_id) {
    record <- readCommittedWorkflowRequest(store, request_id)
    if (is.null(record)) return(NULL)
    directory <- file.path(store$directory, dirname(record$payload$relative_path))
    paths <- file.path(directory, c("commit.json", "descriptor.rds", basename(record$payload$relative_path)))
    syncSavedWorkflowPaths(store, c(paths, directory, dirname(directory),
        file.path(store$directory, "store.json"), store$directory, dirname(store$directory)))
    store$io_counts$recovery_flushes <- store$io_counts$recovery_flushes + 1L
    record
}

#' Validate committed metadata, consuming payload bytes only when requested.
#' @param store Existing owner.
#' @param request_id Existing immutable request identifier.
#' @param payload Whether this operation will consume the payload.
#' @return Accepted descriptor or NULL; routine reads never perform recovery flushing.
#' @noRd
readCommittedWorkflowRequest <- function(store, request_id, payload = TRUE) {
    relative <- savedWorkflowRequestPath(store, request_id)
    directory <- file.path(store$directory, relative)
    if (!dir.exists(directory)) return(NULL)
    core_paths <- file.path(directory, c("commit.json", "descriptor.rds"))
    for (name in c("commit.json", "descriptor.rds")) {
        checkSavedEntityPath(store, paste0(relative, "/", name))
    }
    if (!all(file.exists(core_paths)) || file.info(core_paths[[1L]])$size > 8192 ||
        file.info(core_paths[[2L]])$size > 65536) {
        cli::cli_abort("committed request is incomplete", class = "saved_entity_corrupt")
    }
    commit <- tryCatch(jsonlite::read_json(core_paths[[1L]]), error = function(error) {
        cli::cli_abort("invalid commit record", class = "saved_entity_corrupt", parent = error)
    })
    if (!is.list(commit) || !identical(commit$schema, 1L) || !identical(commit$descriptor_sha256,
        digest::digest(file = core_paths[[2L]], algo = "sha256"))) {
        cli::cli_abort("descriptor integrity check failed", class = "saved_entity_corrupt")
    }
    record <- tryCatch(readRDS(core_paths[[2L]]), error = function(error) {
        cli::cli_abort("descriptor cannot be decoded", class = "saved_entity_corrupt", parent = error)
    })
    payload_name <- if (is.list(record) && is.list(record$payload) &&
        identical(record$payload$codec, "parquet_duckdb")) "payload.parquet" else
        if (is.list(record) && is.list(record$payload) &&
            identical(record$payload$codec, "qs2")) "payload.qs2" else "payload.rds"
    payload_relative <- paste0(relative, "/", payload_name)
    checkSavedEntityPath(store, payload_relative)
    payload_path <- file.path(directory, payload_name)
    paths <- c(core_paths, payload_path)
    if (!is.list(record) || !is.list(record$payload) || !identical(record$entity_id, commit$entity_id) ||
        !identical(record$payload$relative_path, payload_relative) ||
        !file.exists(payload_path) ||
        !identical(as.numeric(file.info(payload_path)$size), record$payload$bytes)) {
        cli::cli_abort("payload integrity check failed", class = "saved_entity_corrupt")
    }
    store$io_counts$metadata_checks <- store$io_counts$metadata_checks + 1L
    if (payload) validateWorkflowPayloadOnce(store, record)
    registerSavedWorkflowMetadata(store, record)
    record
}

#' Remove only this call's verified redundant staging copy after successful reuse.
#' @param store Store handle.
#' @param staging Exact directory created by the current save call.
#' @param record Already accepted matching descriptor.
#' @return TRUE invisibly; failed/conflicting stages never call this cleanup.
#' @noRd
discardRedundantWorkflowStage <- function(store, staging, record) {
    checkWorkflowEntityStore(store)
    pending <- normalizePath(file.path(store$directory, ".pending"), mustWork = TRUE)
    target <- normalizePath(staging, mustWork = TRUE)
    files <- list.files(target, all.files = TRUE, no.. = TRUE)
    request_id <- basename(dirname(record$payload$relative_path))
    payload_name <- basename(record$payload$relative_path)
    # The caller already hashed this owned file and matched its full descriptor.
    # Reuse that validation; size is an additional guard, not an integrity proof.
    if (!startsWith(pending, paste0(store$directory, "/")) ||
        !identical(dirname(target), pending) || !startsWith(basename(target), paste0(request_id, "-")) ||
        !payload_name %in% c("payload.rds", "payload.qs2", "payload.parquet") ||
        !all(files %in% c(payload_name, "descriptor.rds", "commit.json")) ||
        any(file.info(file.path(target, files))$isdir) ||
        !identical(as.numeric(file.info(file.path(target, payload_name))$size),
            record$payload$bytes)) {
        cli::cli_abort("refusing unsafe staging cleanup", class = "saved_entity_path")
    }
    if (unlink(target, recursive = TRUE) != 0L || dir.exists(target)) {
        cli::cli_abort("redundant staging cleanup failed", class = "saved_entity_write")
    }
    invisible(TRUE)
}

#' Save a complete object and atomically publish its request-owned directory.
#' @param store Store handle.
#' @param object Complete R/S4/history value; never removed or changed here.
#' @param request_id Caller-owned idempotency token.
#' @param identity Complete effective-input identity fields from DIA-005.
#' @param metadata Bounded class/shape/factor summaries.
#' @param kind Entity capability kind.
#' @param lineage Optional parent_id and integer revision.
#' @param codec Explicit checkpoint codec; rds remains the default.
#' @return Accepted descriptor. Staging survives failures; no GUI state is changed.
#' @noRd
saveWorkflowEntityImpl <- function(store, object, request_id, identity, metadata,
                               kind = "checkpoint", lineage = NULL, codec = "rds") {
    if (!is.character(codec) || length(codec) != 1L || is.na(codec) ||
        !codec %in% c("rds", "qs2") || codec == "qs2" && kind != "checkpoint") {
        cli::cli_abort("unsupported complete-object codec", class = "saved_entity_unsupported")
    }
    if (codec == "qs2" && !requireNamespace("qs2", quietly = TRUE)) {
        cli::cli_abort("qs2 is required for this checkpoint codec", class = "saved_entity_unsupported")
    }
    payload_name <- if (codec == "qs2") "payload.qs2" else "payload.rds"
    relative <- savedWorkflowRequestPath(store, request_id)
    parent_id <- if (is.null(lineage)) NA_character_ else lineage$parent_id
    revision <- if (is.null(lineage)) 0L else lineage$revision
    # Validate metadata before creating any payload or changing the index.
    newSavedWorkflowDescriptor(store, kind, identity,
        list(relative_path = paste0(relative, "/", payload_name), bytes = 0,
            sha256 = strrep("0", 64), codec = codec), metadata, parent_id, revision)
    if (!identical(class(object), metadata$class)) {
        cli::cli_abort("payload class differs from descriptor", class = "saved_entity_invalid")
    }
    checkSavedWorkflowPayload(object)
    savedWorkflowSyncCommand(store)
    for (name in c("entities", ".pending")) {
        checkSavedEntityPath(store, name)
        directory <- file.path(store$directory, name)
        if (!dir.exists(directory) && !dir.create(directory)) {
            cli::cli_abort("cannot create storage directory", class = "saved_entity_write")
        }
    }
    staging <- tempfile(paste0(request_id, "-"), tmpdir = file.path(store$directory, ".pending"))
    if (!dir.create(staging)) cli::cli_abort("cannot stage request", class = "saved_entity_write")
    savedWorkflowStorageBoundary("before_write")
    checkWorkflowEntityStore(store)
    payload <- file.path(staging, payload_name)
    tryCatch({
        if (codec == "rds") writeSavedWorkflowPayload(object, payload) else
            writeSavedWorkflowQsPayload(object, payload)
    }, error = function(error) {
        cli::cli_abort("payload serialization failed", class = "saved_entity_write", parent = error)
    })
    savedWorkflowStorageBoundary("after_payload")
    checkWorkflowEntityStore(store)
    if (codec == "rds") tryCatch(infoRDS(payload), error = function(error) {
        cli::cli_abort("serialized payload header is invalid", class = "saved_entity_write", parent = error)
    })
    store$io_counts$payload_hashes <- store$io_counts$payload_hashes + 1L
    store$io_counts$payload_hash_bytes <- store$io_counts$payload_hash_bytes + as.numeric(file.info(payload)$size)
    record <- newSavedWorkflowDescriptor(store, kind, identity,
        list(relative_path = paste0(relative, "/", payload_name), bytes = as.numeric(file.info(payload)$size),
            sha256 = hashWorkflowOwnedFile(store, payload), codec = codec),
        metadata, parent_id, revision)
    existing <- recoverSavedWorkflowRequest(store, request_id)
    if (!is.null(existing)) {
        if (!identical(existing, record)) cli::cli_abort("request id conflicts with accepted data", class = "saved_entity_conflict")
        discardRedundantWorkflowStage(store, staging, record)
        return(existing)
    }
    if (store$count >= store$max_entries ||
        store$metadata_bytes + as.numeric(object.size(record)) + 512 > store$max_metadata_bytes) {
        cli::cli_abort("metadata index is full", class = "saved_entity_limit")
    }
    descriptor <- file.path(staging, "descriptor.rds")
    tryCatch({
        saveRDS(record, descriptor, compress = FALSE, version = 3)
        jsonlite::write_json(list(schema = 1L, entity_id = record$entity_id,
            descriptor_sha256 = digest::digest(file = descriptor, algo = "sha256")),
            file.path(staging, "commit.json"), auto_unbox = TRUE)
    }, error = function(error) {
        cli::cli_abort("descriptor staging failed", class = "saved_entity_write", parent = error)
    })
    savedWorkflowStorageBoundary("after_record")
    checkWorkflowEntityStore(store)
    syncSavedWorkflowPaths(store, c(payload, descriptor, file.path(staging, "commit.json"),
        file.path(store$directory, "store.json"), staging, dirname(staging),
        file.path(store$directory, "entities"), store$directory, dirname(store$directory)))
    destination <- file.path(store$directory, relative)
    if (!file.rename(staging, destination)) {
        existing <- recoverSavedWorkflowRequest(store, request_id)
        if (!identical(existing, record)) cli::cli_abort("concurrent request conflict", class = "saved_entity_conflict")
        discardRedundantWorkflowStage(store, staging, record)
        return(existing)
    }
    savedWorkflowStorageBoundary("after_rename")
    syncSavedWorkflowPaths(store, c(destination, dirname(destination), dirname(staging)))
    savedWorkflowStorageBoundary("after_commit")
    registerSavedWorkflowMetadata(store, record)
    rememberWorkflowPayloadValidation(store, record)
    store$io_counts$output_bytes <- store$io_counts$output_bytes + record$payload$bytes
    record
}

#' Save a complete entity, reporting interruption without discarding recovery input.
#' @param store,object,request_id,identity,metadata,kind,lineage,codec Inputs to saveWorkflowEntityImpl.
#' @return Accepted descriptor; interruption raises saved_entity_cancelled.
#' @noRd
saveWorkflowEntity <- function(store, object, request_id, identity, metadata,
                               kind = "checkpoint", lineage = NULL, codec = "rds") {
    tryCatch(saveWorkflowEntityImpl(store, object, request_id, identity, metadata, kind, lineage, codec),
        interrupt = function(error) cli::cli_abort("entity save cancelled",
            class = "saved_entity_cancelled", parent = error))
}

#' Validate an accepted descriptor and its owned file without decoding it.
#' @param store Store handle.
#' @param reference Candidate accepted descriptor.
#' @return The matching durable descriptor.
#' @noRd
validateSavedWorkflowReference <- function(store, reference, payload = TRUE) {
    checkWorkflowEntityStore(store)
    if (!is.list(reference)) {
        cli::cli_abort("invalid hydration reference", class = "saved_entity_invalid")
    }
    expected <- newSavedWorkflowDescriptor(store, reference$kind, reference$identity,
        reference$payload, reference$metadata, reference$parent_id, reference$revision)
    if (!identical(reference, expected)) {
        cli::cli_abort("invalid hydration reference", class = "saved_entity_scope")
    }
    relative <- reference$payload$relative_path
    payload_name <- if (identical(reference$payload$codec, "parquet_duckdb"))
        "payload.parquet" else if (identical(reference$payload$codec, "qs2"))
            "payload.qs2" else "payload.rds"
    expected_relative <- paste0("^entities/[A-Za-z0-9_-]{1,80}/",
        gsub("\\.", "\\\\.", payload_name), "$")
    if (!grepl(expected_relative, relative)) {
        cli::cli_abort("reference is not an accepted request", class = "saved_entity_invalid")
    }
    request_id <- strsplit(relative, "/", fixed = TRUE)[[1L]][[2L]]
    accepted <- readCommittedWorkflowRequest(store, request_id, payload)
    if (is.null(accepted)) {
        cli::cli_abort("saved request is missing", class = "saved_entity_missing")
    }
    if (!identical(accepted, reference)) {
        cli::cli_abort("accepted reference differs", class = "saved_entity_conflict")
    }
    accepted
}

#' Hydrate a complete payload only after full metadata/size/hash validation.
#' @param store Store handle.
#' @param reference Descriptor returned by save/recovery, never an arbitrary path.
#' @param validated Internal TRUE only after validateSavedWorkflowReference in the same call.
#' @return Reconstructed R value with native attributes and classes.
#' @noRd
hydrateSavedWorkflowEntity <- function(store, reference, validated = FALSE) {
    if (!is.logical(validated) || length(validated) != 1L || is.na(validated)) {
        cli::cli_abort("invalid hydration validation flag", class = "saved_entity_invalid")
    }
    if (!isTRUE(validated)) validateSavedWorkflowReference(store, reference)
    if (!reference$payload$codec %in% c("rds", "qs2")) {
        cli::cli_abort("selective tables require readSavedWorkflowTable",
            class = "saved_entity_unsupported")
    }
    if (reference$payload$codec == "qs2" && !requireNamespace("qs2", quietly = TRUE)) {
        cli::cli_abort("qs2 is required to hydrate this checkpoint",
            class = "saved_entity_unsupported")
    }
    relative <- reference$payload$relative_path
    cached <- get0(reference$entity_id, envir = store$cache, inherits = FALSE)
    if (!is.null(cached)) {
        store$cache_hits <- store$cache_hits + 1L
        store$cache_order <- c(store$cache_order[store$cache_order != reference$entity_id],
            reference$entity_id)
        return(cached$value)
    }
    store$payload_reads <- store$payload_reads + 1L
    object <- tryCatch({
        path <- file.path(store$directory, relative)
        if (reference$payload$codec == "qs2") {
            qs2::qs_read(path, validate_checksum = TRUE, nthreads = 2L)
        } else readRDS(path)
    }, error = function(error) {
        cli::cli_abort("payload cannot be decoded", class = "saved_entity_corrupt", parent = error)
    })
    store$payload_decodes <- store$payload_decodes + 1L
    if (!identical(class(object), reference$metadata$class)) {
        cli::cli_abort("hydrated class differs from descriptor", class = "saved_entity_corrupt")
    }
    checkSavedWorkflowPayload(object)
    cacheWorkflowEntityPayload(store, reference$entity_id, object)
    object
}

#' Resolve one exact consumer entity without running a producer on a miss.
#' @param store Open store.
#' @param kind Consumer entity kind.
#' @param identity Complete immutable identity.
#' @param factors Exact typed lookup factors.
#' @param hydrate Decode the matched value; FALSE returns only its trusted descriptor.
#' @return Explicit miss or hit with capabilities and optional value.
#' @noRd
resolveSavedWorkflowEntity <- function(store, kind, identity, factors = list(),
                                       hydrate = TRUE) {
    if (!kind %in% c("scientific_product", "da_fit", "result_table", "rendered_asset") ||
        !is.logical(hydrate) || length(hydrate) != 1L || is.na(hydrate) ||
        !is.list(identity) || !identical(sort(names(identity)),
            sort(c("input", "design", "parameters", "selectors", "software")))) {
        cli::cli_abort("invalid consumer resolution request", class = "saved_entity_invalid")
    }
    newSavedWorkflowScientificIdentity(identity$input, identity$design,
        identity$parameters, identity$selectors, identity$software)
    checkSavedWorkflowFactors(factors)
    found <- findSavedWorkflowEntities(store,
        list(kind = kind, identity = identity, factors = factors), max_matches = 2L)
    if (identical(found$status, "miss")) {
        return(list(status = "miss", reference = NULL, capabilities = NULL, value = NULL))
    }
    if (!identical(found$status, "hit")) {
        cli::cli_abort("ambiguous saved consumer identity", class = "saved_entity_conflict")
    }
    reference <- found$matches[[1L]]
    validateSavedWorkflowReference(store, reference, payload = isTRUE(hydrate))
    capabilities <- savedWorkflowEntityCapabilities(reference)
    value <- if (isTRUE(hydrate))
        hydrateSavedWorkflowEntity(store, reference, validated = TRUE) else NULL
    if (isTRUE(hydrate)) checkSavedWorkflowConsumerPayload(kind, value)
    list(status = "hit", reference = reference, capabilities = capabilities, value = value)
}

#' Resolve an asset keyed to a trusted product without hydrating that product.
#' @param store Open store.
#' @param product_reference Trusted source product/fit/table descriptor.
#' @param render_identity Identity returned by newSavedWorkflowRenderIdentity.
#' @param hydrate Decode the asset envelope.
#' @return Explicit asset miss or hit.
#' @noRd
resolveSavedWorkflowAsset <- function(store, product_reference, render_identity,
                                      hydrate = TRUE) {
    validateSavedWorkflowReference(store, product_reference, payload = FALSE)
    if (identical(product_reference$kind, "rendered_asset") ||
        !identical(render_identity$input, product_reference$entity_id) ||
        !identical(render_identity$design, "rendered-asset-v1")) {
        cli::cli_abort("asset identity does not match its product", class = "saved_entity_scope")
    }
    resolveSavedWorkflowEntity(store, "rendered_asset", render_identity,
        factors = list(), hydrate = hydrate)
}

#' Route bounded column/key requests only to a capable saved-table codec.
#' @param store Open store.
#' @param reference Trusted saved descriptor.
#' @param columns Requested columns.
#' @param filters Requested typed filters.
#' @param reader Explicit table reader, or NULL for the supplied scope/default.
#' @param scope Optional action-owned table scope.
#' @return Projected table result, or an explicit unsupported error for RDS.
#' @noRd
hydrateSavedWorkflowProjection <- function(store, reference, columns, filters = list(),
    reader = NULL, scope = NULL) {
    if (!is.character(columns) || !length(columns) || anyNA(columns) ||
        any(!nzchar(columns)) || anyDuplicated(columns) || !is.list(filters)) {
        cli::cli_abort("invalid projection request", class = "saved_entity_invalid")
    }
    if (!is.list(reference) || !is.list(reference$payload)) {
        cli::cli_abort("invalid projection reference", class = "saved_entity_invalid")
    }
    if (isWorkflowTableReference(reference)) {
        allowed <- c("keys", "offset", "limit", "max_rows", "max_columns", "max_bytes")
        if (length(filters) && (is.null(names(filters)) || anyDuplicated(names(filters)) ||
            !all(names(filters) %in% allowed))) {
            cli::cli_abort("unsupported projection filter", class = "saved_entity_invalid")
        }
        defaults <- list(keys = NULL, offset = 0L, limit = 100L,
            max_rows = 10000L, max_columns = 64L, max_bytes = 67108864)
        defaults[names(filters)] <- filters
        return(do.call(readSavedWorkflowTable,
            c(list(store = store, reference = reference, columns = columns,
                reader = reader, scope = scope), defaults)))
    }
    validateSavedWorkflowReference(store, reference)
    cli::cli_abort("complete-object entities do not support projection or predicate pushdown",
        class = "saved_entity_unsupported")
}

#' Save a checkpoint through the single complete-object codec.
#' @param store,object,request_id,identity,metadata,lineage,codec Inputs to saveWorkflowEntity.
#' @return Accepted checkpoint descriptor.
#' @noRd
saveWorkflowCheckpoint <- function(store, object, request_id, identity, metadata,
    lineage = NULL, codec = "rds") {
    saveWorkflowEntity(store, object, request_id, identity, metadata,
        kind = "checkpoint", lineage = lineage, codec = codec)
}

#' Hydrate a checkpoint, refusing descriptors of another entity kind.
#' @param store Store handle.
#' @param reference Accepted checkpoint descriptor.
#' @return Complete checkpoint object.
#' @noRd
hydrateSavedWorkflowCheckpoint <- function(store, reference) {
    if (!is.list(reference) || !identical(reference$kind, "checkpoint")) {
        cli::cli_abort("reference is not a checkpoint", class = "saved_entity_invalid")
    }
    hydrateSavedWorkflowEntity(store, reference)
}


#' Probe or request native reclamation of already-free allocator pages.
#' @param release FALSE probes support; TRUE requests page release.
#' @return Logical result, or NA when unsupported; not a byte count.
#' @useDynLib MultiScholaR, .registration = TRUE, .fixes = "C_"
#' @noRd
workflowAllocatorRelease <- function(release = FALSE) {
    .Call(C_workflow_reclaim_pages, release)
}

#' Reclaim free pages at a caller-verified completed-action boundary.
#' @param native_fn Registered capability/release operation, injectable for checks.
#' @param gc_fn Existing collector, injectable for lifecycle checks.
#' @return Bounded status; no scientific objects or allocation estimates.
#' @noRd
reclaimWorkflowMemory <- function(native_fn = workflowAllocatorRelease, gc_fn = gc) {
    tryCatch({
        if (!isTRUE(native_fn(FALSE))) return(list(status = "unsupported"))
        before <- snapshotWorkflowMemory()
        collected <- gc_fn(full = TRUE)
        after_gc <- snapshotWorkflowMemory()
        released <- native_fn(TRUE)
        after_release <- snapshotWorkflowMemory()
        live <- if (is.matrix(collected) && all(c("Ncells", "Vcells") %in% rownames(collected)))
            sum(collected[c("Ncells", "Vcells"), "used"] * c(if (.Machine$sizeof.pointer == 8L) 56 else 28, 8)) else NULL
        list(status = if (isTRUE(released)) "released" else "no_pages_released",
            before = before, after_gc = after_gc, after_release = after_release,
            live_R_bytes = live, note = "process/allocator observations overlap; not additive ownership totals")
    }, error = function(error) cli::cli_abort("workflow memory reclamation failed",
        class = "workflow_reclaim_error", parent = error))
}

#' Fingerprint this package's installed native libraries, including empty legacy builds.
#' @param package_path Installed MultiScholaR path.
#' @return Named SHA-256 list in stable path order.
#' @noRd
workflowNativeLibraryIdentity <- function(package_path = system.file(package = "MultiScholaR")) {
    root <- file.path(package_path, "libs")
    paths <- sort(list.files(root, pattern = "[.](so|dll|dylib)$", recursive = TRUE, full.names = TRUE))
    stats::setNames(lapply(paths, function(path) digest::digest(file = path, algo = "sha256")),
        substring(paths, nchar(root) + 2L))
}

#' Release clean pages for an explicitly owned file or completed output prefix.
#' @param path Validated regular path; no global cache operations or file deletion.
#' @param offset,length Byte interval; zero length means the closed file's remaining bytes.
#' @param synchronize Synchronize a completed append-only prefix before advising release.
#' @return Logical advice acceptance, or NA where unsupported; never a reclaimed byte count.
#' @noRd
releaseWorkflowFilePages <- function(path, offset = 0, length = 0, synchronize = FALSE) {
    tryCatch(.Call(C_workflow_release_file_pages, path, as.double(offset), as.double(length), synchronize),
        error = function(error) cli::cli_abort("owned file-page release failed",
            class = "saved_entity_write", parent = error))
}

#' Bound the cache of an append-only staged file as complete prefixes accumulate.
#' @param owner Private writer with path and advised_bytes; never a shared mutable input.
#' @param final TRUE after closing the writer, to release the remaining tail.
#' @return NULL invisibly; file contents and publication status are untouched.
#' @noRd
releaseWorkflowWrittenPrefix <- function(owner, final = FALSE) {
    path <- owner$write_path %||% owner$path
    size <- file.size(path)
    if (length(size) != 1L || is.na(size)) cli::cli_abort("staged output disappeared", class = "saved_entity_write")
    previous <- owner$advised_bytes %||% 0
    if (!final && size - previous < 67108864) return(invisible(NULL))
    # 64 KiB is page-aligned on the supported Linux targets; leave the active tail.
    end <- if (final) size else floor(size / 65536) * 65536
    if (end > previous) {
        releaseWorkflowFilePages(path, previous, end - previous, synchronize = TRUE)
        owner$advised_bytes <- end
    }
    invisible(NULL)
}

#' Hash a store-owned payload without retaining a whole large file in page cache.
#' @param store Current immutable store owner.
#' @param path Staged or committed regular file inside that store.
#' @return Complete SHA-256 fingerprint, with unchanged integrity requirements.
#' @noRd
hashWorkflowOwnedFile <- function(store, path) {
    checkWorkflowEntityStore(store)
    if (!startsWith(normalizePath(path, mustWork = TRUE), paste0(store$directory, "/")))
        cli::cli_abort("payload hash escaped store ownership", class = "saved_entity_path")
    hashWorkflowFile(path, log_path = file.path(store$directory, "checkpoint-io.log"))
}

#' Observe bounded process and allocator counters without collection or traversal.
#' @return Scalar counters; missing platform observations are explicit NULL/NA values.
#' @noRd
snapshotWorkflowMemory <- function() {
    allocator <- as.list(.Call(C_workflow_allocator_snapshot))
    process <- NULL
    if (file.exists("/proc/self/status")) {
        lines <- readLines("/proc/self/status", warn = FALSE)
        fields <- c("VmRSS", "RssAnon", "RssFile", "RssShmem")
        process <- stats::setNames(lapply(fields, function(field) {
            line <- grep(paste0("^", field, ":"), lines, value = TRUE)
            if (length(line) != 1L) return(NA_real_)
            as.numeric(strsplit(trimws(sub("^[^:]+:", "", line)), " +")[[1L]][[1L]]) * 1024
        }), fields)
    }
    list(time = as.numeric(Sys.time()), pid = Sys.getpid(), process = process, allocator = allocator)
}

#' Compute a complete file checksum with bounded clean-cache residency.
#' @param path Validated regular input owned by the calling store or publication.
#' @param log_path Optional bounded I/O log; NULL retains ordinary small-file hashing.
#' @return SHA-256; modification during hashing is always an error.
#' @noRd
hashWorkflowFile <- function(path, log_path = NULL) {
    identity <- workflowTableFileIdentity(path)
    hash <- if (is.null(log_path) || identity$size < 4194304) {
        digest::digest(file = path, algo = "sha256")
    } else {
        runDIANativeIO(list(c("sync", "--", path)), log_path)
        output <- runDIANativeIO(list(c("dd", paste0("if=", path), "bs=4M",
            "iflag=nocache", "status=none"), "sha256sum"), log_path)
        if (length(output) != 1L || !grepl("^[a-f0-9]{64}  -$", output))
            cli::cli_abort("streaming fingerprint is invalid", class = "saved_entity_corrupt")
        substr(output, 1L, 64L)
    }
    if (!identical(identity, workflowTableFileIdentity(path)))
        cli::cli_abort("file changed during checksum validation", class = "saved_entity_corrupt")
    hash
}


#' Observe bounded progress without publishing a scientific state or redrawing Shiny.
#' @param directory Existing request-owned evidence directory.
#' @param current_state Zero-argument bounded state-name accessor.
#' @param interval Minimum seconds between updates for an unchanged phase.
#' @return Scalar-only observer; the latest receipt is replaced atomically.
#' @noRd
newWorkflowProgressObserver <- function(directory, current_state, interval = 3) {
    stopifnot(dir.exists(directory), is.function(current_state), is.numeric(interval),
        length(interval) == 1L, is.finite(interval), interval >= 0)
    last_key <- NULL
    last_time <- -Inf
    function(stage, completed, total = NA_real_, entity = "", unit = "rows") {
        stopifnot(is.character(stage), length(stage) == 1L, nchar(stage) <= 256L,
            is.character(entity), length(entity) == 1L, nchar(entity) <= 256L,
            is.character(unit), length(unit) == 1L, nchar(unit) <= 128L,
            is.numeric(completed), length(completed) == 1L, is.finite(completed), completed >= 0,
            is.numeric(total), length(total) == 1L, is.na(total) || is.finite(total) && total >= completed)
        stamp <- unname(proc.time()[["elapsed"]])
        key <- paste(stage, entity, unit, sep = ":")
        complete <- !is.na(total) && completed == total
        if (identical(key, last_key) && stamp - last_time < interval && !complete) return(invisible(NULL))
        state <- current_state()
        checkSavedEntityString(state, "progress state")
        value <- list(time = as.numeric(Sys.time()), current_state = state, stage = stage,
            completed = completed, total = if (is.na(total)) NULL else total,
            unit = unit, entity = entity, phase_complete = complete)
        encoded <- jsonlite::toJSON(value, auto_unbox = TRUE, null = "null", digits = NA)
        if (nchar(encoded, type = "bytes") > 8192L) cli::cli_abort("progress receipt exceeds limit")
        temporary <- tempfile(".progress-", tmpdir = directory)
        on.exit(unlink(temporary), add = TRUE)
        writeLines(encoded, temporary)
        if (!file.rename(temporary, file.path(directory, "gui-progress-latest.json")))
            cli::cli_abort("cannot publish progress receipt", class = "saved_entity_write")
        last_key <<- key
        last_time <<- stamp
        invisible(NULL)
    }
}
