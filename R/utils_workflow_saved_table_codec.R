#' Build bounded R restoration metadata for one supported data frame.
#' @param table Plain data frame.
#' @param key_column Existing unique key column or .rownames.
#' @param metadata_version Explicit restoration format; legacy version 1 is unchanged.
#' @return Safe metadata embedded in the Parquet schema.
#' @noRd
workflowSavedTableMetadata <- function(table, key_column, metadata_version = 1L) {
    if (!metadata_version %in% c(1L, 2L, 3L)) {
        cli::cli_abort("unsupported table metadata version", class = "saved_entity_unsupported")
    }
    supported_class <- identical(class(table), "data.frame") ||
        metadata_version %in% c(2L, 3L) && (identical(class(table), c("tbl_df", "tbl", "data.frame")) ||
            identical(class(table), c("spec_tbl_df", "tbl_df", "tbl", "data.frame")))
    if (!supported_class || !ncol(table) ||
        anyDuplicated(names(table)) || any(grepl("^\\.multischolar_", names(table))) ||
        any(names(table) %in% c(".source_order", ".request_order")) ||
        !is.character(key_column) || length(key_column) != 1L || is.na(key_column) ||
        !key_column %in% c(".rownames", names(table))) {
        cli::cli_abort("unsupported saved table schema", class = "saved_entity_invalid")
    }
    supported <- vapply(table, \(column) is.factor(column) ||
        typeof(column) %in% c("logical", "integer", "double", "character"), logical(1))
    if (!all(supported)) {
        cli::cli_abort("saved table contains an unsupported column type",
            class = "saved_entity_unsupported")
    }
    native_rows <- .row_names_info(table, type = 0L)
    automatic_key <- identical(key_column, ".rownames") && (!nrow(table) ||
        is.integer(native_rows) && length(native_rows) == 2L && is.na(native_rows[[1L]]) && native_rows[[2L]] < 0L)
    key <- if (automatic_key) character() else if (identical(key_column, ".rownames")) rownames(table) else table[[key_column]]
    if (!automatic_key && (is.factor(key) || !(is.character(key) || is.integer(key) || is.double(key)) ||
        length(key) != nrow(table) || anyNA(key) || anyDuplicated(key))) {
        cli::cli_abort("saved table key must be complete and unique",
            class = "saved_entity_invalid")
    }
    column_metadata <- lapply(table, function(column) {
        if (metadata_version %in% c(2L, 3L) && (!is.null(dim(column)) ||
            !is.null(attr(column, "class", exact = TRUE)) &&
            !identical(class(column), "factor") &&
            !identical(class(column), c("ordered", "factor")))) {
            cli::cli_abort("version 2 requires plain atomic or factor columns",
                class = "saved_entity_unsupported")
        }
        extra <- attributes(column)
        extra[c("names", "class", "levels")] <- NULL
        if (is.null(extra)) extra <- list()
        checkSavedMetadataTree(extra)
        spec <- list(typeof = typeof(column), class = class(column),
            levels = if (is.factor(column)) levels(column) else NULL,
            ordered = is.ordered(column), attributes = extra)
        if (metadata_version %in% c(2L, 3L)) spec$attribute_order <- names(attributes(column)) %||% character()
        spec
    })
    table_attributes <- attributes(table)
    table_attributes[c("names", "row.names", "class")] <- NULL
    if (is.null(table_attributes)) table_attributes <- list()
    checkWorkflowTableAttributes(table_attributes)
    metadata <- list(schema = 1L, codec = "parquet_duckdb_v1",
        columns = names(table), rows = nrow(table), key_column = key_column,
        key_type = typeof(key), row_names = if (metadata_version == 1L) rownames(table) else NULL,
        class = class(table), attributes = table_attributes,
        column_metadata = column_metadata)
    if (metadata_version %in% c(2L, 3L)) {
        native_names <- .row_names_info(table, type = 0L)
        automatic <- !nrow(table) || is.integer(native_names) && length(native_names) == 2L &&
            is.na(native_names[[1L]]) && native_names[[2L]] < 0L
        named_columns <- names(table)[vapply(table, \(column) !is.null(names(column)), logical(1))]
        metadata$schema <- metadata_version
        metadata$codec <- paste0("parquet_duckdb_v", metadata_version)
        metadata$row_names <- NULL
        metadata$automatic_row_names <- automatic
        metadata$row_index_type <- typeof(native_names)
        metadata$attribute_order <- names(attributes(table))
        metadata["row_name_column"] <- list(if (automatic) NULL else if (key_column == ".rownames")
            ".multischolar_row_key" else ".multischolar_row_names")
        metadata$name_columns <- stats::setNames(
            paste0(".multischolar_names_", seq_along(named_columns)), named_columns)
        if (!length(named_columns)) metadata$name_columns <- stats::setNames(character(), character())
    }
    if (as.numeric(object.size(metadata)) > 8 * 1024^2) {
        cli::cli_abort("saved table metadata exceeds 8 MiB", class = "saved_entity_limit")
    }
    metadata
}

#' Validate the table attributes used by DIA, including portable readr schemas.
#' @param attributes Table attributes with names, row names and class removed.
#' @return TRUE invisibly; executable values and live parser pointers are rejected.
#' @noRd
checkWorkflowTableAttributes <- function(attributes) {
    spec <- attributes$spec
    if (!is.null(spec)) {
        if (!identical(class(spec), "col_spec") ||
            !all(c("cols", "default") %in% names(spec)) ||
            !all(names(spec) %in% c("cols", "default", "skip", "delim")) || anyDuplicated(names(spec))) {
            cli::cli_abort("unsupported parser schema", class = "saved_entity_unsupported")
        }
        for (collector in c(spec$cols, list(spec$default))) {
            if (length(class(collector)) != 2L || !inherits(collector, "collector") ||
                !startsWith(class(collector)[[1L]], "collector_")) {
                cli::cli_abort("unsupported parser collector", class = "saved_entity_unsupported")
            }
            checkSavedMetadataTree(unclass(collector))
        }
        checkSavedMetadataTree(unclass(spec)[setdiff(names(spec), c("cols", "default"))])
        attributes$spec <- NULL
    }
    problems <- attributes$problems
    if (!is.null(problems)) {
        if (!identical(problems, methods::new("externalptr"))) {
            if (!is.data.frame(problems) || !all(vapply(problems,
                \(column) is.atomic(column) && is.null(attributes(column)), logical(1)))) {
                cli::cli_abort("parser diagnostics must be portable", class = "saved_entity_unsupported")
            }
        }
        attributes$problems <- NULL
    }
    checkSavedMetadataTree(attributes)
    invisible(TRUE)
}

#' Write one table and its bounded R metadata to a Parquet file.
#' @param table Plain data frame.
#' @param key_column Existing unique key column or .rownames.
#' @param path Fresh staged Parquet path.
#' @param metadata_version Explicit restoration format; defaults to version 1.
#' @param metadata Optional metadata already constructed by the save action.
#' @param backend Existing writer backend selected explicitly by the operation.
#' @param consume Optional summary callback on the bounded original rows.
#' @return Restoration metadata invisibly.
#' @noRd
writeSavedWorkflowParquet <- function(table, key_column, path, metadata_version = 1L,
    metadata = NULL, progress = NULL, backend = "arrow", consume = NULL) {
    if (is.null(metadata)) metadata <- workflowSavedTableMetadata(table, key_column, metadata_version)
    writer <- openWorkflowTableWriter(path, metadata, progress, backend = backend)
    on.exit(abortWorkflowTableWriter(writer), add = TRUE)
    width <- 96 + sum(vapply(table, function(column) {
        labels <- if (is.factor(column)) levels(column) else column
        bytes <- if (is.character(labels)) 80 + max(c(0, nchar(labels, type = "bytes")),
            na.rm = TRUE) else 16
        if (!is.null(names(column))) bytes <- bytes + 80 +
            max(c(0, nchar(names(column), type = "bytes")), na.rm = TRUE)
        bytes
    }, numeric(1)))
    if (nrow(table) && width > 67108864)
        cli::cli_abort("one persisted row exceeds its transfer limit", class = "saved_entity_limit")
    batch_rows <- max(1L, as.integer(min(10000, floor(33554432 / max(1, width)))))
    for (start in seq.int(1L, max(1L, nrow(table)), by = batch_rows)) {
        at <- if (nrow(table)) seq.int(start, min(nrow(table), start + batch_rows - 1L)) else integer()
        if (!is.null(consume)) consume(table[at, , drop = FALSE])
        appendWorkflowTableWriter(writer, table, at)
    }
    finishWorkflowTableWriter(writer)
    invisible(metadata)
}

#' Check the installed incremental writer API and qualified native build.
#' @param backend Requested writer mode; only direct R columns are qualified.
#' @return Source-bound capability metadata; release eligibility remains held.
#' @noRd
workflowFastWriterInfo <- function(backend = "fast_r") {
    fail <- function() cli::cli_abort(
        "direct writer requires the qualified incremental fastduckdb build",
        class = "saved_entity_unsupported")
    if (!requireNamespace("fastduckdb", quietly = TRUE) ||
        !"parquet_writer_capabilities" %in% getNamespaceExports("fastduckdb")) fail()
    path <- system.file("duckdb", "linux_amd64", "manifest.rds", package = "MultiScholaR")
    if (!nzchar(path) || file.size(path) > 65536) fail()
    manifest <- readRDS(path)
    actual <- fastduckdb::parquet_writer_capabilities()
    expected <- manifest$fastduckdb
    if (!identical(manifest$schema, 3L) || is.null(expected) ||
        !identical(actual$api, "incremental-parquet-v2") || !isTRUE(actual$supported) ||
        !identical(backend, "fast_r") || !identical(actual$input_modes, "r_columns") ||
        !identical(actual$max_threads, 8L) || !identical(actual$max_batch_rows, 10000L) ||
        !identical(actual$max_batch_bytes, 67108864) ||
        !identical(actual, expected$capability)) fail()
    required <- c("R/fastduckdb.rdb", "R/fastduckdb.rdx", "NAMESPACE")
    if (!identical(names(expected$files), required)) fail()
    paths <- file.path(system.file(package = "fastduckdb"), required)
    if (any(!file.exists(paths))) fail()
    hashes <- stats::setNames(vapply(paths, function(path)
        digest::digest(file = path, algo = "sha256"), character(1)), required)
    if (!identical(hashes, expected$files)) fail()
    binary <- getLoadedDLLs()[["fastduckdb"]][["path"]]
    if (!identical(digest::digest(file = binary, algo = "sha256"), expected$binary_sha256)) fail()
    list(capability = actual, binary_sha256 = expected$binary_sha256,
        runtime_release_status = manifest$runtime_release_status)
}

#' Open one incremental Parquet writer with the existing physical codec.
#' @param path Owned, fresh staged file.
#' @param metadata Exact restoration schema; v3 derives its row count from the footer.
#' @param progress Optional scalar-only observer supplied by the action owner.
#' @param backend Explicit qualified writer; ordinary callers retain Arrow.
#' @param optimization Private measured writer candidate; ordinary callers use none.
#' @return Private writer; callers must finish or abort it.
#' @noRd
openWorkflowTableWriter <- function(path, metadata, progress = NULL, backend = "arrow",
    optimization = "none") {
    if (file.exists(path)) cli::cli_abort("table writer destination exists", class = "saved_entity_conflict")
    if (!backend %in% c("arrow", "fast_arrow", "fast_r")) cli::cli_abort("unsupported table writer")
    if (!optimization %in% c("none", "metadata", "batch") ||
        optimization != "none" && backend != "fast_r") cli::cli_abort("unsupported writer optimization")
    if (backend != "arrow") workflowFastWriterInfo(backend)
    if (backend == "fast_arrow" && !requireNamespace("nanoarrow", quietly = TRUE))
        cli::cli_abort("direct Arrow writer requires nanoarrow", class = "saved_entity_unsupported")
    writer <- new.env(parent = emptyenv())
    writer$path <- path
    writer$backend <- backend
    writer$optimization <- optimization
    writer$schema_metadata <- NULL
    writer$write_path <- NULL
    writer$metadata <- metadata
    writer$rows <- 0L
    writer$written_rows <- 0L
    writer$progress <- progress
    writer$advised_bytes <- 0
    writer$output <- writer$sink <- NULL
    writer$pending <- list()
    writer$pending_rows <- writer$pending_bytes <- 0
    writer$closed <- FALSE
    writer
}

#' Encode a bounded row range without making a whole-table Arrow copy.
#' @param writer Existing incremental writer.
#' @param table Already-owned scientific table.
#' @param at Ordered source rows for this batch.
#' @return NULL invisibly.
#' @noRd
appendWorkflowTableWriter <- function(writer, table, at = seq_len(nrow(table))) {
    if (isTRUE(writer$closed) || length(at) > 10000L) cli::cli_abort("invalid table writer batch")
    metadata <- writer$metadata
    storage <- if (writer$optimization == "batch" && identical(at, seq_len(nrow(table))))
        lapply(table[metadata$columns], identity) else lapply(table[metadata$columns], function(column) {
        if (is.null(dim(column))) column[at] else column[at, , drop = FALSE]
    })
    if (as.numeric(object.size(storage)) > 67108864) cli::cli_abort("table production batch exceeds limits")
    n <- length(at)
    storage$.multischolar_order <- as.integer(writer$rows + seq_len(n))
    source_names <- if (metadata$key_column != ".rownames" && isTRUE(metadata$automatic_row_names)) NULL else
        if (isTRUE(metadata$automatic_row_names)) as.character(storage$.multischolar_order) else
        if (metadata$schema == 1L) metadata$row_names[at] else rownames(table)[at]
    if (metadata$key_column == ".rownames") storage$.multischolar_row_key <- source_names
    if (metadata$schema %in% c(2L, 3L)) {
        if (!is.null(metadata$row_name_column)) storage[[metadata$row_name_column]] <- source_names
        for (column in names(metadata$name_columns)) {
            storage[[metadata$name_columns[[column]]]] <- names(table[[column]])[at]
        }
    }
    attr(storage, "row.names") <- .set_row_names(n)
    class(storage) <- "data.frame"
    schema_metadata <- writer$schema_metadata
    if (is.null(schema_metadata)) {
        encoded_metadata <- metadata
        if (metadata$schema == 3L) encoded_metadata["rows"] <- list(NULL)
        # The workflow codec owns restoration, including row-sized R names.
        schema_metadata <- list(multischolar_table =
            jsonlite::base64_enc(serialize(encoded_metadata, NULL, version = 3)))
        if (writer$optimization == "metadata") writer$schema_metadata <- schema_metadata
    }
    if (writer$backend == "fast_r") {
        if (is.null(writer$output)) {
            writer$output <- fastduckdb::open_parquet_writer_fast(storage[0, , drop = FALSE],
                writer$path, metadata = schema_metadata, n_threads = 8L)
            updateWorkflowFastWriter(writer, fastduckdb::flush_parquet_writer_fast(writer$output))
        }
        if (n) updateWorkflowFastWriter(writer, fastduckdb::append_parquet_writer_fast(writer$output, storage))
        writer$rows <- writer$rows + n
        return(invisible(NULL))
    }
    arrow_table <- arrow::Table$create(storage)
    encoded_table <- arrow_table$ReplaceSchemaMetadata(schema_metadata)
    arrow_table <- encoded_table
    if (writer$backend == "fast_arrow") {
        if (is.null(writer$output)) {
            writer$output <- fastduckdb::open_parquet_writer_fast(
                nanoarrow::as_nanoarrow_array_stream(arrow_table), writer$path, n_threads = 8L)
            updateWorkflowFastWriter(writer, fastduckdb::flush_parquet_writer_fast(writer$output))
        }
        if (n) updateWorkflowFastWriter(writer, fastduckdb::append_parquet_writer_fast(writer$output,
            nanoarrow::as_nanoarrow_array_stream(arrow_table), owner = arrow_table,
            logical_bytes = arrow_table$nbytes()))
        writer$rows <- writer$rows + n
        return(invisible(NULL))
    }
    if (is.null(writer$output)) {
        writer$sink <- arrow::FileOutputStream$create(writer$path)
        writer$output <- arrow::ParquetFileWriter$create(arrow_table$schema, writer$sink,
            properties = arrow::ParquetWriterProperties$create(names(storage), compression = "snappy"))
    }
    # Read batches stay bounded independently of Parquet row groups. Arrow
    # concatenation retains chunks rather than collecting another R table.
    if (n) {
        if (writer$pending_rows && (writer$pending_rows + n > 122880L ||
            writer$pending_bytes + arrow_table$nbytes() > 33554432)) {
            flushWorkflowTableWriter(writer)
        }
        writer$pending[[length(writer$pending) + 1L]] <- arrow_table
        writer$pending_rows <- writer$pending_rows + n
        writer$pending_bytes <- writer$pending_bytes + arrow_table$nbytes()
        if (writer$pending_bytes >= 33554432) flushWorkflowTableWriter(writer)
    }
    writer$rows <- writer$rows + n
    invisible(NULL)
}

#' Apply native writer counters while keeping progress and file-cache ownership unchanged.
#' @param writer Existing private table writer.
#' @param info Scalar native-writer counters and its current physical path.
#' @return NULL invisibly.
#' @noRd
updateWorkflowFastWriter <- function(writer, info) {
    previous <- writer$written_rows
    writer$written_rows <- as.integer(info$written_rows)
    writer$pending_rows <- as.integer(info$pending_rows)
    writer$pending_bytes <- info$pending_bytes
    writer$write_path <- if (isTRUE(info$closed)) writer$path else info$temporary_path
    if (writer$written_rows > previous) {
        if (is.function(writer$progress)) writer$progress("table_written", writer$written_rows,
            if (writer$metadata$schema == 3L) NA_real_ else writer$metadata$rows,
            basename(dirname(writer$path)), "rows")
        releaseWorkflowWrittenPrefix(writer)
    }
    invisible(NULL)
}

#' Flush one bounded Arrow row group without reconstructing an R table.
#' @param writer Private incremental writer.
#' @return NULL invisibly.
#' @noRd
flushWorkflowTableWriter <- function(writer) {
    if (writer$backend %in% c("fast_arrow", "fast_r")) {
        updateWorkflowFastWriter(writer, fastduckdb::flush_parquet_writer_fast(writer$output))
        return(invisible(NULL))
    }
    if (!writer$pending_rows) return(invisible(NULL))
    combined <- do.call(arrow::concat_tables, c(writer$pending, list(unify_schemas = FALSE)))
    writer$output$WriteTable(combined, chunk_size = writer$pending_rows)
    writer$written_rows <- writer$written_rows + writer$pending_rows
    if (is.function(writer$progress)) writer$progress("table_written", writer$written_rows,
        if (writer$metadata$schema == 3L) NA_real_ else writer$metadata$rows,
        basename(dirname(writer$path)), "rows")
    writer$pending <- list()
    writer$pending_rows <- writer$pending_bytes <- 0
    releaseWorkflowWrittenPrefix(writer)
    invisible(NULL)
}

#' Close and check an incremental table without reading its payload.
#' @param writer Existing private writer.
#' @return Footer row count and restoration metadata.
#' @noRd
finishWorkflowTableWriter <- function(writer) {
    if (!writer$closed) {
        if (is.null(writer$output)) cli::cli_abort("table writer has no schema", class = "saved_entity_write")
        flushWorkflowTableWriter(writer)
        if (writer$backend %in% c("fast_arrow", "fast_r")) {
            updateWorkflowFastWriter(writer, fastduckdb::close_parquet_writer_fast(writer$output))
        } else {
            writer$output$Close()
            writer$sink$close()
        }
        writer$output <- writer$sink <- NULL
        writer$closed <- TRUE
        writer$progress <- NULL
        releaseWorkflowWrittenPrefix(writer, final = TRUE)
    }
    metadata <- readSavedWorkflowTableMetadata(writer$path)
    if (metadata$rows != writer$rows) cli::cli_abort("table writer row count differs", class = "saved_entity_corrupt")
    metadata
}

#' Close an interrupted writer and remove only its incomplete staged payload.
#' @param writer Existing private writer.
#' @return NULL invisibly; completed staged payloads remain recoverable.
#' @noRd
abortWorkflowTableWriter <- function(writer) {
    if (isTRUE(writer$closed)) return(invisible(NULL))
    writer$pending <- list()
    writer$pending_rows <- writer$pending_bytes <- 0
    if (!is.null(writer$output)) {
        if (writer$backend %in% c("fast_arrow", "fast_r"))
            try(fastduckdb::abort_parquet_writer_fast(writer$output), silent = TRUE) else
            try(writer$output$Close(), silent = TRUE)
    }
    owned_destination <- identical(writer$backend, "arrow") && !is.null(writer$sink)
    if (!is.null(writer$sink)) try(writer$sink$close(), silent = TRUE)
    writer$output <- writer$sink <- NULL
    writer$closed <- TRUE
    writer$progress <- NULL
    # Native abort owns its private temporary; a final destination may belong to another writer.
    if (owned_destination) unlink(writer$path)
    invisible(NULL)
}

#' Read and validate only the bounded Parquet schema metadata.
#' @param path Validated owned Parquet file.
#' @return Restoration metadata without collecting table rows.
#' @noRd
readSavedWorkflowTableMetadata <- function(path) {
    reader <- tryCatch(arrow::ParquetFileReader$create(path),
        error = function(error) cli::cli_abort("Parquet schema cannot be read",
            class = "saved_entity_corrupt", parent = error))
    on.exit(try(reader$.unsafe_delete(), silent = TRUE), add = TRUE)
    schema <- reader$GetSchema()
    encoded <- schema$metadata$multischolar_table
    if (!is.character(encoded) || length(encoded) != 1L || is.na(encoded) ||
        nchar(encoded, type = "bytes") > 16 * 1024^2) {
        cli::cli_abort("Parquet restoration metadata is missing or oversized",
            class = "saved_entity_corrupt")
    }
    metadata <- tryCatch(unserialize(jsonlite::base64_dec(encoded)),
        error = function(error) cli::cli_abort("Parquet restoration metadata is invalid",
            class = "saved_entity_corrupt", parent = error))
    if (is.list(metadata) && identical(metadata$schema, 3L)) {
        if (!is.null(metadata$rows)) cli::cli_abort("streaming row count must come from the footer",
            class = "saved_entity_corrupt")
        metadata$rows <- as.integer(reader$num_rows)
    }
    required <- c("schema", "codec", "columns", "rows", "key_column", "key_type",
        "row_names", "class", "attributes", "column_metadata")
    version_two <- is.list(metadata) && metadata$schema %in% c(2L, 3L)
    if (version_two) required <- c(setdiff(required, "row_names"),
        "automatic_row_names", "row_name_column", "name_columns")
    if (version_two && "row_index_type" %in% names(metadata)) {
        required <- c(required, "row_index_type")
        if (!is.character(metadata$row_index_type) || length(metadata$row_index_type) != 1L ||
            !metadata$row_index_type %in% c("integer", "character")) {
            cli::cli_abort("invalid row-name restoration type", class = "saved_entity_corrupt")
        }
    }
    if (version_two && "attribute_order" %in% names(metadata)) {
        required <- c(required, "attribute_order")
        order <- metadata$attribute_order
        if (!is.character(order) || anyNA(order) || anyDuplicated(order) ||
            !setequal(order, c("names", "row.names", "class", names(metadata$attributes)))) {
            cli::cli_abort("invalid table attribute ordering", class = "saved_entity_corrupt")
        }
    }
    if (!is.list(metadata) || !identical(sort(names(metadata)), sort(required)) ||
        !metadata$schema %in% c(1L, 2L, 3L) ||
        !identical(metadata$codec, paste0("parquet_duckdb_v", metadata$schema)) ||
        !(identical(metadata$class, "data.frame") || version_two &&
            (identical(metadata$class, c("tbl_df", "tbl", "data.frame")) ||
                identical(metadata$class, c("spec_tbl_df", "tbl_df", "tbl", "data.frame")))) ||
        !is.character(metadata$columns) || anyDuplicated(metadata$columns) ||
        !is.numeric(metadata$rows) || length(metadata$rows) != 1L ||
        metadata$rows < 0 || !version_two && length(metadata$row_names) != metadata$rows ||
        !identical(names(metadata$column_metadata), metadata$columns)) {
        cli::cli_abort("Parquet restoration metadata differs from its schema",
            class = "saved_entity_corrupt")
    }
    if (version_two) {
        extras <- c(metadata$row_name_column, unname(metadata$name_columns))
        if (!is.logical(metadata$automatic_row_names) || length(metadata$automatic_row_names) != 1L ||
            is.na(metadata$automatic_row_names) ||
            !is.character(metadata$name_columns) || anyNA(metadata$name_columns) ||
            anyDuplicated(names(metadata$name_columns)) ||
            !all(names(metadata$name_columns) %in% metadata$columns) ||
            !all(grepl("^\\.multischolar_names_[0-9]+$", metadata$name_columns)) ||
            !is.null(metadata$row_name_column) &&
                !metadata$row_name_column %in% c(".multischolar_row_key", ".multischolar_row_names") ||
            !isTRUE(metadata$automatic_row_names) && is.null(metadata$row_name_column) ||
            anyDuplicated(extras) || !all(extras %in% schema$names) ||
            as.numeric(object.size(metadata)) > 8 * 1024^2) {
            cli::cli_abort("invalid version 2 table metadata", class = "saved_entity_corrupt")
        }
    }
    checkWorkflowTableAttributes(metadata$attributes)
    metadata
}

#' Capture the identity of a validated owned table during one read action.
#' @param path Accepted Parquet path.
#' @return File identity used to detect replacement or mutation within the action.
#' @noRd
workflowTableFileIdentity <- function(path) {
    link <- Sys.readlink(path)
    info <- file.info(path)
    if (!file.exists(path) || !is.na(link) && nzchar(link) || nrow(info) != 1L ||
        is.na(info$size) || isTRUE(info$isdir)) {
        cli::cli_abort("saved table changed during read", class = "saved_entity_corrupt")
    }
    list(path = normalizePath(path, mustWork = TRUE), size = info$size,
        mtime = info$mtime, ctime = info$ctime, mode = info$mode)
}

#' Open one validated Parquet/DuckDB read action.
#' @param store Open saved-entity store.
#' @param reference Accepted Parquet table descriptor.
#' @param reader Explicit duckdb (default) or optional fastduckdb reader.
#' @return Scoped read owner; releaseWorkflowTableReadScope must be called.
#' @noRd
newWorkflowTableReadScope <- function(store, reference, reader = "duckdb", metadata_only = FALSE) {
    checkWorkflowEntityStore(store)
    requested_reference <- reference
    component <- readWorkflowTableComponent(store, reference)
    if (!is.null(component)) reference <- component$source
    if (!is.character(reader) || length(reader) != 1L || is.na(reader) ||
        !reader %in% c("duckdb", "fastduckdb", "arrow")) {
        cli::cli_abort("unsupported table reader", class = "saved_entity_unsupported")
    }
    if (reader == "fastduckdb" && !requireNamespace("fastduckdb", quietly = TRUE)) {
        cli::cli_abort("optional fastduckdb reader is unavailable", class = "saved_entity_unsupported")
    }
    if (!is.list(reference) || !is.list(reference$payload) ||
        !identical(reference$payload$codec, "parquet_duckdb")) {
        cli::cli_abort("table read scope requires a Parquet descriptor",
            class = "saved_entity_unsupported")
    }
    if (!identical(reference$owner, store$identity)) {
        cli::cli_abort("table belongs to another store", class = "saved_entity_scope")
    }
    relative <- reference$payload$relative_path
    checkSavedEntityPath(store, relative)
    if (!grepl("^entities/[A-Za-z0-9_-]{1,80}/payload\\.parquet$", relative)) {
        cli::cli_abort("invalid table payload path", class = "saved_entity_path")
    }
    path <- file.path(store$directory, relative)
    action <- beginWorkflowTableAction(store)
    ready <- FALSE
    on.exit(if (!ready) endWorkflowTableAction(store, action), add = TRUE)
    before <- workflowTableFileIdentity(path)
    validateSavedWorkflowReference(store, reference, payload = !metadata_only)
    if (!identical(before, workflowTableFileIdentity(path))) {
        cli::cli_abort("saved table changed during validation",
            class = "saved_entity_corrupt")
    }
    metadata <- readSavedWorkflowTableMetadata(path)
    physical_metadata <- metadata
    if (!is.null(component)) {
        if (is.null(component$view))
            checkWorkflowTableRestoration(component$metadata, physical_metadata)
        metadata <- component$metadata
    }
    if (!identical(before, workflowTableFileIdentity(path))) {
        cli::cli_abort("saved table changed during metadata read",
            class = "saved_entity_corrupt")
    }
    if (!isTRUE(all.equal(as.numeric(physical_metadata$rows),
        as.numeric(reference$metadata$shape[[1L]])))) {
        cli::cli_abort("table footer differs from descriptor", class = "saved_entity_corrupt")
    }
    if (!metadata_only && reader == "duckdb" && is.null(action$spill)) {
        pending <- file.path(store$directory, ".pending")
        checkSavedEntityPath(store, ".pending")
        if (!dir.exists(pending) && !dir.create(pending)) cli::cli_abort("cannot create table staging directory")
        action$spill <- tempfile("table-spill-", tmpdir = pending)
        if (!dir.create(action$spill)) cli::cli_abort("cannot create table spill directory")
    }
    if (reader == "duckdb" && !metadata_only && is.null(action$connection)) {
        action$connection <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE,
            config = list(temp_directory = "")), dbdir = ":memory:")
        DBI::dbExecute(action$connection, "SET memory_limit = '1GB'")
        DBI::dbExecute(action$connection, "SET threads = 2")
        DBI::dbExecute(action$connection, "SET max_temp_directory_size = '0B'")
        DBI::dbExecute(action$connection, "SET temp_directory = ''")
    }
    connection <- if (reader == "duckdb" && !metadata_only) action$connection else NULL
    spill <- action$spill
    scope <- new.env(parent = emptyenv())
    scope$store <- store
    scope$reference <- requested_reference
    scope$physical_reference <- reference
    scope$physical_metadata <- physical_metadata
    scope$path <- path
    scope$file_identity <- before
    scope$metadata <- metadata
    scope$view <- component$view
    scope$view_id <- if (is.null(component$view)) reference$entity_id else component$view$view_id
    scope$map_path <- if (!is.null(component$view) && !is.null(component$view$rows$map))
        file.path(store$directory, component$view$rows$map$payload$relative_path) else NULL
    scope$map_file_identity <- if (!is.null(scope$map_path)) workflowTableFileIdentity(scope$map_path) else NULL
    scope$connection <- connection
    scope$spill_directory <- spill
    scope$reader <- reader
    scope$metadata_only <- metadata_only
    scope$action <- action
    scope$writers <- new.env(parent = emptyenv())
    scope$readers <- new.env(parent = emptyenv())
    scope$view_scopes <- list()
    action$scopes <- action$scopes + 1L
    scope$released <- FALSE
    store$open_resources <- store$open_resources + 1L
    ready <- TRUE
    scope
}

#' Refuse reuse after project replacement or table mutation.
#' @param scope Active table read owner.
#' @return TRUE invisibly.
#' @noRd
checkWorkflowTableReadScope <- function(scope) {
    if (!is.environment(scope) || isTRUE(scope$released) ||
        !is.environment(scope$store) ||
        identical(scope$reader, "duckdb") && !isTRUE(scope$metadata_only) && is.null(scope$connection)) {
        cli::cli_abort("table read scope is closed", class = "saved_entity_scope")
    }
    checkWorkflowEntityStore(scope$store)
    if (!identical(scope$file_identity, workflowTableFileIdentity(scope$path))) {
        cli::cli_abort("saved table changed during read", class = "saved_entity_corrupt")
    }
    if (!is.null(scope$map_path) &&
        !identical(scope$map_file_identity, workflowTableFileIdentity(scope$map_path)))
        cli::cli_abort("selected row map changed during read", class = "saved_entity_corrupt")
    invisible(TRUE)
}

#' Release one table action and fail if its source changed before delivery.
#' @param scope Table read owner.
#' @return TRUE invisibly; repeated release is harmless.
#' @noRd
releaseWorkflowTableReadScope <- function(scope) {
    if (!is.environment(scope) || is.null(scope$released)) {
        cli::cli_abort("invalid table read scope", class = "saved_entity_invalid")
    }
    if (isTRUE(scope$released)) return(invisible(TRUE))
    # Every owned resource is released even when one close fails; the first
    # failure (source validity first) is raised only after cleanup completes.
    validity <- tryCatch(checkWorkflowTableReadScope(scope), error = \(error) error)
    keep <- function(expr) {
        outcome <- tryCatch({
            force(expr)
            NULL
        }, error = \(error) error)
        if (inherits(outcome, "error") && !inherits(validity, "error")) validity <<- outcome
        invisible(NULL)
    }
    for (child in scope$view_scopes) keep(releaseWorkflowTableReadScope(child))
    scope$view_scopes <- NULL
    keep(observeWorkflowTableScratch(scope))
    for (cursor in as.list(scope$readers, all.names = TRUE)) keep(cursor$close())
    scope$readers <- NULL
    for (writer in as.list(scope$writers, all.names = TRUE)) keep(abortWorkflowTableWriter(writer))
    scope$writers <- NULL
    scope$connection <- NULL
    scope$spill_directory <- NULL
    scope$action$scopes <- scope$action$scopes - 1L
    # Close a borrowed database as soon as the last table scope releases it.
    if (scope$action$scopes == 0L && !is.null(scope$action$connection)) {
        keep(workflowQueryDisconnect(scope$action$connection))
        scope$action$connection <- NULL
    }
    scope$metadata <- NULL
    scope$reference <- NULL
    scope$view <- scope$map_path <- scope$map_file_identity <- NULL
    scope$released <- TRUE
    scope$store$open_resources <- scope$store$open_resources - 1L
    keep(endWorkflowTableAction(scope$store, scope$action))
    scope$action <- NULL
    scope$store <- NULL
    if (inherits(validity, "error")) stop(validity)
    invisible(TRUE)
}

#' Save one selective table through the existing durable request lifecycle.
#' @param store Shared saved-entity store.
#' @param table Plain data frame.
#' @param request_id Idempotent request identity.
#' @param identity Complete scientific identity.
#' @param key_column Existing unique key or .rownames.
#' @param factors Bounded typed lookup factors.
#' @param metadata_version Explicit table restoration format; legacy version 1 remains default.
#' @param metadata,write_fn Internal bounded producer metadata and staged writer.
#' @return Accepted result_table descriptor using parquet_duckdb.
#' @noRd
saveWorkflowTable <- function(store, table, request_id, identity, key_column,
                              factors = list(), metadata_version = 1L,
                              metadata = NULL, write_fn = NULL) {
    checkWorkflowEntityStore(store)
    newSavedWorkflowScientificIdentity(identity$input, identity$design,
        identity$parameters, identity$selectors, identity$software)
    checkSavedWorkflowFactors(factors)
    if (is.null(metadata)) metadata <- workflowSavedTableMetadata(table, key_column, metadata_version)
    relative <- savedWorkflowRequestPath(store, request_id)
    for (name in c("entities", ".pending")) {
        checkSavedEntityPath(store, name)
        path <- file.path(store$directory, name)
        if (!dir.exists(path) && !dir.create(path)) {
            cli::cli_abort("cannot create table storage directory",
                class = "saved_entity_write")
        }
    }
    staging <- tempfile(paste0(request_id, "-"),
        tmpdir = file.path(store$directory, ".pending"))
    if (!dir.create(staging)) {
        cli::cli_abort("cannot stage table request", class = "saved_entity_write")
    }
    payload <- file.path(staging, "payload.parquet")
    savedWorkflowStorageBoundary("before_write")
    if (is.null(write_fn)) {
        writeSavedWorkflowParquet(table, key_column, payload, metadata_version, metadata, store$progress)
    } else {
        write_fn(payload, metadata)
        produced <- readSavedWorkflowTableMetadata(payload)
        if (identical(metadata$schema, 3L) && is.null(metadata$rows)) metadata$rows <- produced$rows
        if (!identical(produced, metadata)) {
            cli::cli_abort("produced table metadata differs", class = "saved_entity_corrupt")
        }
    }
    savedWorkflowStorageBoundary("after_payload")
    payload_info <- file.info(payload)
    store$io_counts$payload_hashes <- store$io_counts$payload_hashes + 1L
    store$io_counts$payload_hash_bytes <- store$io_counts$payload_hash_bytes + as.numeric(payload_info$size)
    record <- newSavedWorkflowDescriptor(store, "result_table", identity,
        list(relative_path = paste0(relative, "/payload.parquet"),
            bytes = as.numeric(payload_info$size),
            sha256 = hashWorkflowOwnedFile(store, payload),
            codec = "parquet_duckdb"),
        list(class = metadata$class, shape = as.integer(c(metadata$rows, length(metadata$columns))),
            factors = factors))
    saveRDS(record, file.path(staging, "descriptor.rds"), version = 3)
    commit <- list(schema = 1L, entity_id = record$entity_id,
        descriptor_sha256 = digest::digest(file = file.path(staging, "descriptor.rds"),
            algo = "sha256"))
    jsonlite::write_json(commit, file.path(staging, "commit.json"), auto_unbox = TRUE)
    savedWorkflowStorageBoundary("after_descriptor")
    existing <- recoverSavedWorkflowRequest(store, request_id)
    if (!is.null(existing)) {
        if (!identical(existing, record)) {
            cli::cli_abort("same table request has different content",
                class = "saved_entity_conflict")
        }
        discardRedundantWorkflowStage(store, staging, record)
        return(existing)
    }
    if (store$count >= store$max_entries ||
        store$metadata_bytes + as.numeric(object.size(record)) + 512 >
            store$max_metadata_bytes) {
        cli::cli_abort("metadata index is full", class = "saved_entity_limit")
    }
    paths <- c(payload, file.path(staging, "descriptor.rds"),
        file.path(staging, "commit.json"))
    syncSavedWorkflowPaths(store, c(paths, staging, dirname(staging),
        file.path(store$directory, "store.json"), store$directory,
        dirname(store$directory), file.path(store$directory, "entities")))
    destination <- file.path(store$directory, relative)
    if (!file.rename(staging, destination)) {
        existing <- recoverSavedWorkflowRequest(store, request_id)
        if (!identical(existing, record)) {
            cli::cli_abort("concurrent table request conflict",
                class = "saved_entity_conflict")
        }
        discardRedundantWorkflowStage(store, staging, record)
        return(existing)
    }
    savedWorkflowStorageBoundary("after_rename")
    syncSavedWorkflowPaths(store, c(destination, dirname(destination), dirname(staging)))
    savedWorkflowStorageBoundary("after_commit")
    registerSavedWorkflowMetadata(store, record)
    rememberWorkflowPayloadValidation(store, record)
    store$io_counts$output_bytes <- store$io_counts$output_bytes + record$payload$bytes
    releaseWorkflowStoreFilePages(store, file.path(store$directory, record$payload$relative_path))
    record
}

#' Persist an action-owned SQL result directly through the existing table commit.
#' @param scope Open table read scope owning the database lease.
#' @param sql Internal query returning physical codec columns in declared row order.
#' @param metadata Restoration metadata for the result, including its row count.
#' @param request_id,identity Existing durable request and scientific identity.
#' @return Immutable table descriptor, with no collected result retained in R.
#' @noRd
saveWorkflowTableQuery <- function(scope, sql, metadata, request_id, identity) {
    checkWorkflowTableReadScope(scope)
    assertWorkflowPhysicalTableScope(scope)
    if (is.null(scope$connection) || !metadata$schema %in% c(2L, 3L)) {
        cli::cli_abort("query production requires the scoped version 2 codec",
            class = "saved_entity_unsupported")
    }
    writer <- function(path, restoration) {
        checkWorkflowTableReadScope(scope)
        encoded_metadata <- restoration
        if (identical(restoration$schema, 3L)) encoded_metadata["rows"] <- list(NULL)
        encoded <- jsonlite::base64_enc(serialize(encoded_metadata, NULL, version = 3))
        DBI::dbExecute(scope$connection, paste0("COPY (", sql, ") TO ",
            DBI::dbQuoteString(scope$connection, path),
            " (FORMAT PARQUET, COMPRESSION SNAPPY, KV_METADATA {multischolar_table: ",
            DBI::dbQuoteString(scope$connection, encoded), "})"))
        size <- readSavedWorkflowTableMetadata(path)$rows
        if (!is.null(restoration$rows) && !isTRUE(size == restoration$rows)) {
            cli::cli_abort("produced table row count differs", class = "saved_entity_corrupt")
        }
    }
    saveWorkflowTable(scope$store, NULL, request_id, identity, metadata$key_column,
        metadata = metadata, write_fn = writer)
}

#' Configure one existing scoped connection for qualified bulk transforms.
#' @param scope Physical saved-table lease owning a DuckDB connection.
#' @return Immutable execution settings attached to the scope.
#' @noRd
configureWorkflowBulkQueryScope <- function(scope) {
    checkWorkflowTableReadScope(scope)
    assertWorkflowPhysicalTableScope(scope)
    if (is.null(scope$connection)) cli::cli_abort(
        "bulk query requires a scoped DuckDB connection",
        class = "saved_entity_unsupported")
    settings <- list(threads = 8L, engine_limit = "12GiB", spill = FALSE)
    workflowQueryExecute(scope$connection, "SET threads = 8")
    workflowQueryExecute(scope$connection, "SET memory_limit = '12GiB'")
    workflowQueryExecute(scope$connection, "SET max_temp_directory_size = '0B'")
    workflowQueryExecute(scope$connection, "SET temp_directory = ''")
    scope$bulk_query_settings <- settings
    settings
}

#' Collect bounded metadata from one bulk query, never an assay result.
#' @param scope Configured physical DuckDB lease.
#' @param sql Fixed internal query without a terminal semicolon.
#' @param params Bound query parameters.
#' @param max_rows,max_bytes Collection ceilings.
#' @return Data frame admitted before R transfer, in query order. The engine
#'   working set remains subject to the separate process guard.
#' @noRd
collectWorkflowBulkMetadata <- function(scope, sql, params = list(),
    max_rows = 10000L, max_bytes = 8388608) {
    checkWorkflowTableReadScope(scope)
    if (is.null(scope$bulk_query_settings) || is.null(scope$connection) ||
        !is.character(sql) || length(sql) != 1L || is.na(sql) || !nzchar(sql) ||
        grepl(";[[:space:]]*$", sql) || !is.list(params) ||
        !is.numeric(max_rows) || length(max_rows) != 1L || is.na(max_rows) ||
        !is.finite(max_rows) || max_rows < 0 || max_rows != floor(max_rows) ||
        max_rows > 10000L || !is.numeric(max_bytes) || length(max_bytes) != 1L ||
        is.na(max_bytes) || !is.finite(max_bytes) || max_bytes < 1 || max_bytes > 67108864)
        cli::cli_abort("invalid bulk metadata query", class = "saved_entity_limit")
    bounded <- paste0("SELECT * FROM (", sql,
        ") AS workflow_bounded_metadata LIMIT ", as.integer(max_rows) + 1L)
    con <- scope$connection
    quote <- \(x) as.character(workflowQueryQuoteIdentifier(con, x))
    fields <- workflowQueryCollect(con, paste("DESCRIBE", bounded), params = params)
    ordinal <- ".workflow_metadata_order"
    supported <- c("BOOLEAN", "TINYINT", "SMALLINT", "INTEGER", "BIGINT", "HUGEINT",
        "UTINYINT", "USMALLINT", "UINTEGER", "UBIGINT", "UHUGEINT", "FLOAT", "DOUBLE", "VARCHAR")
    if (!nrow(fields) || nrow(fields) > 64L || anyDuplicated(fields$column_name) ||
        ordinal %in% fields$column_name || !all(fields$column_type %in% supported))
        cli::cli_abort("unsupported bulk metadata schema", class = "saved_entity_unsupported")
    # Materialize only the bounded metadata, once. A second execution of the
    # scientific aggregation for a size preflight would duplicate its cost.
    table <- quote(basename(tempfile("workflow_metadata_")))
    workflowQueryExecute(con, paste0("CREATE TEMP TABLE ", table,
        " AS SELECT row_number() OVER () AS ", quote(ordinal), ", * FROM (", bounded, ")"),
        params = params)
    on.exit(workflowQueryExecute(con, paste("DROP TABLE", table)), add = TRUE)
    # Include vectors, per-string allocation/alignment and container overhead.
    # Counting repeated strings separately is deliberately conservative.
    widths <- ifelse(fields$column_type == "VARCHAR", paste0(
        "80 + coalesce(octet_length(encode(", quote(fields$column_name), ")), 0)"), "16")
    overhead <- 4096 + 256 * nrow(fields) + sum(nchar(fields$column_name, type = "bytes"))
    footprint <- workflowQueryCollect(con, paste0("SELECT count(*) AS rows, ", overhead,
        " + coalesce(sum(", paste(widths, collapse = " + "), "),0) AS bytes FROM ", table))
    if (footprint$rows > max_rows || footprint$bytes > max_bytes)
        cli::cli_abort("bulk query metadata exceeds collection limits",
            class = "saved_entity_limit")
    value <- workflowQueryCollect(con, paste0("SELECT ",
        paste(quote(fields$column_name), collapse = ","), " FROM ", table,
        " ORDER BY ", quote(ordinal)))
    if (nrow(value) > max_rows || as.numeric(object.size(value)) > max_bytes)
        cli::cli_abort("bulk query metadata exceeds collection limits",
            class = "saved_entity_limit")
    value
}

#' Partition ordered independent groups without splitting a group.
#' @param groups Ordered unique scientific group identifiers.
#' @param batches Positive number of substantial target batches.
#' @param max_groups Maximum supported metadata population.
#' @return Stable list whose concatenation is exactly groups.
#' @noRd
partitionWorkflowCompleteGroups <- function(groups, batches, max_groups = 10000L) {
    scalar <- function(x, minimum) is.numeric(x) && length(x) == 1L &&
        !is.na(x) && is.finite(x) && x >= minimum && x == floor(x)
    if (!is.atomic(groups) || !is.null(dim(groups)) ||
        is.object(groups) && !is.factor(groups) || anyNA(groups) ||
        anyDuplicated(groups) || !scalar(batches, 1L) || batches > .Machine$integer.max ||
        !scalar(max_groups, 1L) ||
        max_groups > 10000L || length(groups) > max_groups)
        cli::cli_abort("invalid complete scientific groups", class = "saved_entity_limit")
    if (!length(groups)) return(list(groups))
    labels <- ceiling(seq_along(groups) * as.double(batches) / length(groups))
    unname(split(groups, labels))
}

#' Inspect declared complete partitions without retaining full grouping state.
#' @param scope Owned bulk query scope.
#' @param relation,partition Trusted relation and declared partition column.
#' @param inspect Callback receiving a parameterized predicate and scalar bindings.
#' @param batches NULL visits each partition value; a count visits that many groups of
#'   complete partitions (one group visits `TRUE`), never splitting a partition.
#' @return NULL invisibly after visiting every complete partition (or empty input).
#' @noRd
visitWorkflowCompletePartitions <- function(scope, relation, partition, inspect, batches = NULL) {
    con <- scope$connection
    q <- function(x) as.character(workflowQueryQuoteIdentifier(con, x))
    stopifnot(is.character(relation), length(relation) == 1L,
        is.character(partition), length(partition) == 1L, is.function(inspect))
    if (is.numeric(batches) && length(batches) == 1L && isTRUE(batches == 1)) {
        inspect("TRUE", list())
        return(invisible(NULL))
    }
    values <- collectWorkflowBulkMetadata(scope, paste0("SELECT DISTINCT ",
        q(partition), " AS partition_value FROM ", q(relation), " ORDER BY ", q(partition)))
    if (!nrow(values)) {
        inspect("FALSE", list())
        return(invisible(NULL))
    }
    if (!is.null(batches)) {
        if (anyNA(values$partition_value)) cli::cli_abort(
            "batched partitions require non-missing partition values", class = "workflow_execution_unsupported")
        for (group in partitionWorkflowCompleteGroups(values$partition_value, batches))
            inspect(paste0(q(partition), " IN (", paste(rep("?", length(group)), collapse = ","), ")"),
                as.list(group))
        return(invisible(NULL))
    }
    for (i in seq_len(nrow(values))) inspect(paste0(q(partition), " IS NOT DISTINCT FROM ?"),
        list(values$partition_value[[i]]))
    invisible(NULL)
}

#' Check operation-declared uniqueness with an identity-preserving partition.
#' @param scope Owned bulk query scope.
#' @param relation,partition Trusted relation and complete partition field.
#' @param keys Declared identity columns, including the partition column.
#' @param batches Passed to `visitWorkflowCompletePartitions`.
#' @return Whether any duplicate identity exists; no additional scientific rules.
#' @noRd
hasWorkflowPartitionDuplicates <- function(scope, relation, partition, keys, batches = NULL) {
    if (!is.character(keys) || !length(keys) || anyNA(keys) || anyDuplicated(keys) ||
        !partition %in% keys) cli::cli_abort(
            "uniqueness partition must be part of the declared identity",
            class = "workflow_execution_unsupported")
    con <- scope$connection
    q <- function(x) as.character(workflowQueryQuoteIdentifier(con, x))
    duplicate <- FALSE
    visitWorkflowCompletePartitions(scope, relation, partition, function(predicate, params) {
        if (!duplicate) duplicate <<- nrow(workflowQueryCollect(con, paste0(
            "SELECT 1 FROM ", q(relation), " WHERE ", predicate, " GROUP BY ",
            paste(q(keys), collapse = ","), " HAVING count(*) > 1 LIMIT 1"), params)) > 0L
    }, batches)
    duplicate
}

#' Stage a bounded R-produced batch using the existing physical table codec.
#' @param scope Action-owned database and spill directory.
#' @param name Trusted internal SQL table name.
#' @param batch Bounded scientific result or support table.
#' @param offset Number of already staged rows.
#' @param key_column Existing unique key, avoiding redundant generated string keys.
#' @return Updated row count; no batch is retained by the owner.
#' @noRd
appendWorkflowTableBatch <- function(scope, name, batch, offset = 0L, key_column = ".rownames") {
    checkWorkflowTableReadScope(scope)
    if (!is.data.frame(batch) || nrow(batch) > 10000L || ncol(batch) > 64L ||
        as.numeric(object.size(batch)) > 67108864) {
        cli::cli_abort("table production batch exceeds limits", class = "saved_entity_limit")
    }
    writer <- get0(name, scope$writers, inherits = FALSE)
    if (is.null(writer)) {
        metadata <- workflowSavedTableMetadata(batch, key_column, 3L)
        writer <- openWorkflowTableWriter(tempfile("producer-", tmpdir = workflowTableBatchDirectory(scope),
            fileext = ".parquet"), metadata, scope$store$progress)
        assign(name, writer, scope$writers)
    }
    if (offset != writer$rows) cli::cli_abort("table producer offset differs", class = "saved_entity_invalid")
    appendWorkflowTableWriter(writer, batch)
    scope$store$io_counts$batches_written <- scope$store$io_counts$batches_written + 1L
    writer$rows
}

#' Commit a produced stream directly, without copying or registering a temporary SQL table.
#' @param scope Existing action with a database lease.
#' @param name Exact name passed to appendWorkflowTableBatch.
#' @param identity Complete input/parameter/software identity of the producing action.
#' @return Accepted immutable table descriptor.
#' @noRd
finishWorkflowTableBatch <- function(scope, name, identity) {
    writer <- get(name, scope$writers, inherits = FALSE)
    metadata <- finishWorkflowTableWriter(writer)
    source <- writer$path
    scope$store$io_counts$staged_bytes <- scope$store$io_counts$staged_bytes + file.size(source)
    identity$selectors <- digest::digest(list(identity$selectors, name, "batch-product-v1"), algo = "sha256")
    reference <- saveWorkflowTable(scope$store, NULL, basename(tempfile("table-batch-")),
        identity, metadata$key_column, metadata = metadata,
        write_fn = function(path, restoration) {
            if (!file.rename(source, path)) cli::cli_abort("cannot publish staged batch table", class = "saved_entity_write")
        })
    rm(list = name, envir = scope$writers)
    reference
}

#' Restore supported R column types and attributes after DuckDB collection.
#' @param data Collected data frame.
#' @param metadata Validated table metadata.
#' @param columns Requested original columns.
#' @return Restored data frame.
#' @noRd
restoreSavedWorkflowTable <- function(data, metadata, columns) {
    result <- unclass(data[columns])
    for (name in columns) {
        spec <- metadata$column_metadata[[name]]
        value <- result[[name]]
        if (!is.null(spec$levels)) {
            value <- factor(value, levels = spec$levels, ordered = isTRUE(spec$ordered))
        } else {
            value <- switch(spec$typeof,
                logical = as.logical(value), integer = as.integer(value),
                double = as.double(value), character = as.character(value), value)
            if (length(spec$class)) class(value) <- spec$class
        }
        for (attribute in names(spec$attributes)) {
            attr(value, attribute) <- spec$attributes[[attribute]]
        }
        result[[name]] <- value
    }
    if (metadata$schema %in% c(2L, 3L)) {
        for (name in intersect(columns, names(metadata$name_columns))) {
            names(result[[name]]) <- data[[metadata$name_columns[[name]]]]
        }
        for (name in columns) result[[name]] <- restoreWorkflowAttributeOrder(result[[name]],
            metadata$column_metadata[[name]]$attribute_order)
    }
    for (attribute in names(metadata$attributes)) {
        attr(result, attribute) <- metadata$attributes[[attribute]]
    }
    class(result) <- metadata$class
    restoreWorkflowAttributeOrder(result, metadata$attribute_order)
}

#' Preserve recorded R attribute order at a supported codec boundary.
#' @param value Restored atomic vector or data frame.
#' @param order Optional bounded attribute names from version-two metadata.
#' @return The same scientific value and attributes in their original order.
#' @noRd
restoreWorkflowAttributeOrder <- function(value, order) {
    if (is.null(order)) return(value)
    if (!is.character(order) || anyNA(order) || anyDuplicated(order) || length(order) > 64L) {
        cli::cli_abort("invalid restoration attribute order", class = "saved_entity_corrupt")
    }
    attributes <- attributes(value)
    # attributes() expands compact automatic row names; keep their signed form.
    if (is.data.frame(value)) attributes$row.names <- .row_names_info(value, type = 0L)
    if (!setequal(names(attributes), order)) {
        cli::cli_abort("restored attributes differ from their schema", class = "saved_entity_corrupt")
    }
    attributes(value) <- if (length(order)) attributes[order] else NULL
    value
}

#' Include row-sized restoration fields only for the requested columns.
#' @param metadata Validated restoration metadata.
#' @param columns Requested logical columns.
#' @return Physical column names; no payload is read.
#' @noRd
workflowTableStorageColumns <- function(metadata, columns) {
    unique(c(columns, metadata$row_name_column,
        unname(metadata$name_columns[intersect(names(metadata$name_columns), columns)])))
}

#' Recover source row names from legacy metadata or selected row data.
#' @param data Collected physical columns.
#' @param metadata Validated restoration metadata.
#' @param source_order Selected original row indices, including missing selectors.
#' @return Character source row names, independent of result-container row names.
#' @noRd
workflowTableSourceRowNames <- function(data, metadata, source_order) {
    if (identical(metadata$schema, 1L)) return(metadata$row_names[as.integer(source_order)])
    if (isTRUE(metadata$automatic_row_names)) return(as.character(source_order))
    as.character(data[[metadata$row_name_column]])
}

#' Restore the declared row-name representation for a projected table.
#' @param data Restored bounded table.
#' @param metadata Existing codec restoration contract.
#' @param source_row_names Valid, unique source row names.
#' @return Table with its original integer or character row-name storage.
#' @noRd
restoreWorkflowTableRowNames <- function(data, metadata, source_row_names) {
    if (inherits(data, "tbl_df")) {
        attr(data, "row.names") <- .set_row_names(nrow(data))
    } else if (!isTRUE(metadata$automatic_row_names) &&
        identical(metadata$row_index_type, "integer")) {
        rows <- suppressWarnings(as.integer(source_row_names))
        if (anyNA(rows)) cli::cli_abort("invalid integer row names", class = "saved_entity_corrupt")
        attr(data, "row.names") <- rows
    } else {
        rownames(data) <- source_row_names
    }
    data
}

#' Query one owned Parquet table with bounded projection and ordered selection.
#' @param store Shared saved-entity store.
#' @param reference Trusted parquet_duckdb descriptor.
#' @param columns Existing original columns to collect.
#' @param keys Optional ordered key selectors; duplicates are preserved.
#' @param offset,limit Bounded page request used when keys are absent.
#' @param max_rows,max_columns,max_bytes Hard collection ceilings.
#' @param scope Optional owner shared by reads in one action.
#' @param reader Explicit reader, or NULL to use the supplied scope/default DuckDB.
#' @return Data, source order, missing selectors and measured collection scope.
#' @noRd
readSavedWorkflowTable <- function(store, reference, columns, keys = NULL,
    offset = 0L, limit = 100L, max_rows = 10000L, max_columns = 64L,
    max_bytes = 67108864, scope = NULL, reader = NULL) {
    if (!isWorkflowTableReference(reference)) {
        cli::cli_abort("table projection is unavailable for this codec",
            class = "saved_entity_unsupported")
    }
    if (is.null(reader)) reader <- if (is.null(scope)) "duckdb" else scope$reader
    if (identical(reader, "fastduckdb") && !is.null(keys)) {
        cli::cli_abort("fast reader does not support key lookup", class = "saved_entity_unsupported")
    }
    if (is.null(scope)) {
        scope <- newWorkflowTableReadScope(store, reference, reader)
        on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
    }
    checkWorkflowTableReadScope(scope)
    if (!identical(scope$store, store) || !identical(scope$reference, reference) ||
        !identical(scope$reader, reader)) {
        cli::cli_abort("table read scope differs from request", class = "saved_entity_scope")
    }
    path <- scope$path
    metadata <- scope$metadata
    scalar_integer <- function(value, minimum) is.numeric(value) &&
        length(value) == 1L && !is.na(value) && is.finite(value) &&
        value >= minimum && value == floor(value)
    if (!is.character(columns) || !length(columns) || anyNA(columns) ||
        anyDuplicated(columns) || !all(columns %in% metadata$columns) ||
        length(columns) > max_columns || !scalar_integer(offset, 0L) ||
        !scalar_integer(limit, 1L) || !scalar_integer(max_rows, 1L) ||
        !scalar_integer(max_columns, 1L) || !is.numeric(max_bytes) ||
        length(max_bytes) != 1L || !is.finite(max_bytes) || max_bytes < 1 ||
        limit > max_rows || !is.null(keys) &&
            (!is.atomic(keys) || anyNA(keys) || length(keys) > max_rows)) {
        cli::cli_abort("invalid or oversized saved table request",
            class = "saved_entity_limit")
    }
    if (!is.null(scope$view)) {
        if (!is.null(keys)) cli::cli_abort("key lookup on a selected table is unsupported",
            class = "saved_entity_unsupported")
        return(readWorkflowViewPage(scope, columns, offset, limit, max_bytes))
    }
    physical_columns <- workflowTableStorageColumns(metadata, columns)
    if (identical(reader, "fastduckdb")) {
        supported <- vapply(metadata$column_metadata[columns], function(spec) {
            plain <- switch(spec$typeof, double = "numeric", integer = "integer",
                logical = "logical", character = "character", NA_character_)
            (identical(spec$class, plain) || identical(spec$class, "factor") ||
                identical(spec$class, c("ordered", "factor"))) && is.null(spec$attributes$dim)
        }, logical(1))
        if (!all(supported)) {
            cli::cli_abort("fast reader does not support the requested column classes",
                class = "saved_entity_unsupported")
        }
        # The pinned native cache is process-wide and entry-, not byte-bounded.
        # Do not retain its footers between calls; action R metadata stays owned.
        collected <- tryCatch(fastduckdb::read_parquet_fast(path,
            columns = c(".multischolar_order", physical_columns),
            n_threads = scope$read_threads %||% 2L, lazy_strings = FALSE,
            filter = list(.multischolar_order = c(offset + 1, offset + limit))),
            finally = fastduckdb::clear_parquet_cache())
        store$io_counts$fast_reads <- store$io_counts$fast_reads + 1L
        if (is.unsorted(collected$.multischolar_order)) {
            collected <- collected[order(collected$.multischolar_order), , drop = FALSE]
        }
        names(collected)[names(collected) == ".multischolar_order"] <- ".source_order"
        requested_keys <- NULL
    } else {
        connection <- scope$connection
        quote_id <- \(name) as.character(DBI::dbQuoteIdentifier(connection, name))
        quoted_columns <- vapply(physical_columns, quote_id, character(1))
        select_sql <- paste(quoted_columns, collapse = ", ")
        joined_select_sql <- paste(paste0("t.", quoted_columns), collapse = ", ")
        order_column <- quote_id(".multischolar_order")
        source_order_alias <- quote_id(".source_order")
        request_order_alias <- quote_id(".request_order")
        key_column <- if (identical(metadata$key_column, ".rownames"))
            ".multischolar_row_key" else metadata$key_column
        quoted_key <- quote_id(key_column)
        if (is.null(keys)) {
            sql <- paste0("SELECT ", order_column, " AS ", source_order_alias, ", ", select_sql,
                " FROM read_parquet(?) WHERE ", order_column, " >= ? AND ", order_column, " <= ?")
            parameters <- list(path, as.double(offset) + 1, as.double(offset) + limit)
            collected <- DBI::dbGetQuery(connection, sql, params = parameters)
            if (is.unsorted(collected[[source_order_alias]]))
                collected <- collected[order(collected[[source_order_alias]]), , drop = FALSE]
            requested_keys <- NULL
        } else {
            placeholders <- paste(rep("(?, ?)", length(keys)), collapse = ", ")
            sql <- paste0("WITH requested(ord, requested_key) AS (VALUES ",
                placeholders, ") SELECT requested.ord AS ", request_order_alias, ", t.",
                order_column, " AS ", source_order_alias, ", ", joined_select_sql,
                " FROM requested LEFT JOIN read_parquet(?) AS t ON t.", quoted_key,
                " = requested.requested_key ORDER BY requested.ord")
            parameters <- unlist(lapply(seq_along(keys), \(index)
                list(as.integer(index), keys[[index]])), recursive = FALSE)
            parameters[[length(parameters) + 1L]] <- path
            collected <- DBI::dbGetQuery(connection, sql, params = parameters)
            requested_keys <- keys
        }
    }
    source_order <- collected$.source_order
    found <- !is.na(source_order)
    data <- restoreSavedWorkflowTable(collected, metadata, columns)
    source_row_names <- workflowTableSourceRowNames(collected, metadata, source_order)
    if (is.null(keys) && !anyDuplicated(source_row_names) && !anyNA(source_row_names)) {
        data <- restoreWorkflowTableRowNames(data, metadata, source_row_names)
    }
    bytes <- as.numeric(object.size(data))
    if (bytes > max_bytes) {
        cli::cli_abort("collected saved table exceeds byte ceiling",
            class = "saved_entity_limit")
    }
    checkWorkflowTableReadScope(scope)
    store$table_reads <- store$table_reads + 1L
    store$table_rows <- store$table_rows + nrow(data)
    store$table_columns <- store$table_columns + ncol(data)
    store$table_bytes <- store$table_bytes + bytes
    list(data = data, requested_keys = requested_keys, found = found,
        source_row_names = source_row_names, source_order = source_order,
        selected_columns = columns, collected_bytes = bytes,
        source_rows = metadata$rows, source_columns = length(metadata$columns),
        codec = "parquet_duckdb", projected = TRUE)
}

#' Read one logical view page through its immutable physical ordinals.
#' @param scope,columns,offset,limit,max_bytes Validated bounded page request.
#' @return The ordinary projected-table page contract.
#' @noRd
readWorkflowViewPage <- function(scope, columns, offset, limit, max_bytes) {
    count <- as.integer(min(limit, max(0, scope$metadata$rows - offset)))
    at <- if (count) as.double(offset + seq_len(count)) else double()
    if (!is.null(scope$view$rows$map) && count) {
        map <- readWorkflowViewOrdinals(workflowViewChildScope(scope, "map"),
            c("logical_ordinal", "source_ordinal"), at)$data
        if (!identical(map$logical_ordinal, at))
            cli::cli_abort("logical row map changed", class = "saved_entity_corrupt")
        at <- map$source_ordinal
    }
    mapping <- scope$view$columns[match(columns, scope$view$columns$output), , drop = FALSE]
    source_columns <- mapping$source
    page <- readWorkflowViewOrdinals(workflowViewChildScope(scope, "source"), source_columns, at)
    physical <- page$data
    source_row_names <- page$source_row_names
    names(physical) <- mapping$output
    for (i in which(mapping$cast == "factor_to_character"))
        physical[[i]] <- as.character(physical[[i]])
    buffer <- newWorkflowTableBuffer(scope$metadata, columns, count)
    assignWorkflowTableBuffer(buffer, physical, seq_len(count), as.integer(at),
        source_row_names)
    data <- finishWorkflowTableBuffer(buffer)
    attr(data, "workflow_source_order") <- NULL
    # A page exposes logical row positions just like a physical table page;
    # whole-table assembly restores compact automatic row names separately.
    if (isTRUE(scope$metadata$automatic_row_names))
        data <- restoreWorkflowTableRowNames(data, scope$metadata,
            as.character(offset + seq_len(count)))
    bytes <- as.numeric(object.size(data))
    if (bytes > max_bytes)
        cli::cli_abort("collected saved table exceeds byte ceiling", class = "saved_entity_limit")
    checkWorkflowTableReadScope(scope)
    scope$store$table_reads <- scope$store$table_reads + 1L
    scope$store$table_rows <- scope$store$table_rows + count
    scope$store$table_columns <- scope$store$table_columns + length(columns)
    scope$store$table_bytes <- scope$store$table_bytes + bytes
    list(data = data, requested_keys = NULL, found = rep(TRUE, count),
        source_row_names = source_row_names, source_order = at,
        selected_columns = columns, collected_bytes = bytes,
        source_rows = scope$metadata$rows,
        source_columns = length(scope$metadata$columns), codec = "parquet_duckdb",
        projected = TRUE)
}

#' Retain a physical source/map lease for the lifetime of a logical read action.
#' @param scope Logical view owner.
#' @param role Physical source or row map.
#' @return Existing child lease using the explicitly requested reader.
#' @noRd
workflowViewChildScope <- function(scope, role) {
    checkWorkflowTableReadScope(scope)
    if (is.null(scope$view) || !role %in% c("source", "map"))
        cli::cli_abort("invalid logical reader owner", class = "saved_entity_scope")
    child <- scope$view_scopes[[role]]
    if (is.null(child)) {
        reference <- if (role == "source") scope$physical_reference else scope$view$rows$map
        child <- newWorkflowTableReadScope(scope$store, reference, reader = scope$reader)
        if (scope$reader == "fastduckdb") child$read_threads <- 8L
        scope$view_scopes[[role]] <- child
    }
    child
}

#' Gather a bounded, unique logical-page selection through existing physical readers.
#' @param scope Physical source or map lease with the requested backend.
#' @param columns Explicit projection, at most 64 fields.
#' @param ordinals Unique physical positions, at most 10,000.
#' @return Restored physical values and immutable source row names.
#' @noRd
readWorkflowViewOrdinals <- function(scope, columns, ordinals) {
    checkWorkflowTableReadScope(scope)
    if (length(ordinals) > 10000L || anyDuplicated(ordinals) ||
        anyNA(ordinals) || any(!is.finite(ordinals) | ordinals < 1 |
            ordinals > scope$metadata$rows | ordinals != floor(ordinals)) ||
        !length(columns) || length(columns) > 64L || anyDuplicated(columns) ||
        !all(columns %in% scope$metadata$columns))
        cli::cli_abort("invalid logical-page ordinal gather", class = "saved_entity_limit")
    if (!length(ordinals)) return(list(data = readSavedWorkflowTablePrototype(scope)[columns],
        source_row_names = character()))
    if (scope$reader == "duckdb") {
        data <- readWorkflowSelectedOrdinals(scope$store, scope$reference, columns, ordinals,
            scope = scope)
        return(list(data = data, source_row_names = row.names(data)))
    }
    buffer <- newWorkflowTableBuffer(scope$metadata, columns, length(ordinals))
    source_names <- character(length(ordinals))
    fill <- function(data, positions, row_names) {
        at <- match(positions, ordinals)
        if (anyNA(at)) cli::cli_abort("logical gather returned unexpected rows",
            class = "saved_entity_corrupt")
        assignWorkflowTableBuffer(buffer, data, at, as.integer(positions), row_names)
        source_names[at] <<- row_names
        if (as.numeric(object.size(buffer$columns)) +
            as.numeric(object.size(source_names)) > 67108864)
            cli::cli_abort("logical-page gather exceeds 64 MiB", class = "saved_entity_limit")
    }
    ordered <- sort(as.integer(ordinals), method = "radix")
    if (scope$reader == "fastduckdb") {
        starts <- c(1L, which(diff(ordered) != 1L) + 1L)
        ends <- c(starts[-1L] - 1L, length(ordered))
        for (i in seq_along(starts)) {
            page <- readSavedWorkflowTable(scope$store, scope$reference, columns,
                offset = ordered[[starts[[i]]]] - 1L,
                limit = ends[[i]] - starts[[i]] + 1L, scope = scope)
            fill(page$data, page$source_order, page$source_row_names)
        }
    } else {
        cursor <- openWorkflowArrowCursor(scope, columns, ordered)
        on.exit(cursor$close(), add = TRUE)
        repeat {
            part <- cursor$fetch()
            if (is.null(part)) break
            restored <- restoreWorkflowCursorBatch(scope, columns, part)
            fill(restored$data, part$.source_order, restored$row_names)
        }
    }
    data <- finishWorkflowTableBuffer(buffer)
    attr(data, "workflow_source_order") <- NULL
    list(data = data, source_row_names = source_names)
}

#' Select DIA heatmap rows from already adjusted full-population statistics.
#' @param store,reference Accepted table owner and descriptor.
#' @param columns Existing feature identifier and normalized sample columns.
#' @param top_n Number of ordered features to return.
#' @param scope Optional shared table read action.
#' @param max_bytes Collection ceiling, unchanged from other table reads.
#' @return Exact ordered table projection without statistical recomputation.
#' @noRd
readSavedWorkflowDiaHeatmapTable <- function(store, reference, columns, top_n,
    scope = NULL, max_bytes = 67108864) {
    if (is.null(scope)) {
        scope <- newWorkflowTableReadScope(store, reference)
        on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
    }
    checkWorkflowTableReadScope(scope)
    assertWorkflowPhysicalTableScope(scope)
    if (identical(scope$reader, "fastduckdb")) {
        cli::cli_abort("fast reader does not support ranked heatmap queries",
            class = "saved_entity_unsupported")
    }
    if (!identical(scope$store, store) || !identical(scope$reference, reference)) {
        cli::cli_abort("table read scope differs from request", class = "saved_entity_scope")
    }
    metadata <- scope$metadata
    if (!is.character(columns) || !length(columns) || anyNA(columns) ||
        anyDuplicated(columns) || length(columns) > 64L ||
        !all(columns %in% metadata$columns) ||
        !all(c("fdr_qvalue", "logFC") %in% metadata$columns) ||
        !identical(metadata$column_metadata$fdr_qvalue$typeof, "double") ||
        !identical(metadata$column_metadata$logFC$typeof, "double") ||
        !is.numeric(top_n) || length(top_n) != 1L || is.na(top_n) ||
        !is.finite(top_n) || top_n < 1 || top_n != floor(top_n) || top_n > 10000L ||
        !is.numeric(max_bytes) || length(max_bytes) != 1L ||
        !is.finite(max_bytes) || max_bytes < 1) {
        cli::cli_abort("invalid DIA heatmap table request", class = "saved_entity_limit")
    }
    connection <- scope$connection
    physical_columns <- workflowTableStorageColumns(metadata, columns)
    quoted <- as.character(DBI::dbQuoteIdentifier(connection, physical_columns))
    order_column <- as.character(DBI::dbQuoteIdentifier(connection, ".multischolar_order"))
    sql <- paste0("SELECT ", order_column, " AS ",
        as.character(DBI::dbQuoteIdentifier(connection, ".source_order")), ", ",
        paste(quoted, collapse = ", "), " FROM read_parquet(?) ORDER BY ",
        as.character(DBI::dbQuoteIdentifier(connection, "fdr_qvalue")),
        " ASC NULLS LAST, ABS(",
        as.character(DBI::dbQuoteIdentifier(connection, "logFC")),
        ") DESC NULLS LAST, ", order_column, " ASC LIMIT ?")
    collected <- DBI::dbGetQuery(connection, sql,
        params = list(scope$path, as.integer(top_n)))
    source_order <- collected$.source_order
    data <- restoreSavedWorkflowTable(collected, metadata, columns)
    source_row_names <- workflowTableSourceRowNames(collected, metadata, source_order)
    data <- restoreWorkflowTableRowNames(data, metadata, source_row_names)
    bytes <- as.numeric(object.size(data))
    if (bytes > max_bytes) {
        cli::cli_abort("collected saved table exceeds byte ceiling",
            class = "saved_entity_limit")
    }
    checkWorkflowTableReadScope(scope)
    store$table_reads <- store$table_reads + 1L
    store$table_rows <- store$table_rows + nrow(data)
    store$table_columns <- store$table_columns + ncol(data)
    store$table_bytes <- store$table_bytes + bytes
    list(data = data, source_row_names = source_row_names,
        source_order = source_order, collected_bytes = bytes,
        full_population_rows = metadata$rows, statistics_recomputed = FALSE)
}

#' Reject abundance-only projection as a Limpa model input.
#' @param store Shared saved-entity store.
#' @param reference Candidate table reference.
#' @return Never returns; DIA-017 must supply a complete aligned EList codec.
#' @noRd
readSavedWorkflowLimpaModelInput <- function(store, reference) {
    validateSavedWorkflowReference(store, reference)
    cli::cli_abort(paste("Limpa model input requires the complete aligned EList",
        "including genes and uncertainty support; table projection is incompatible"),
        class = "saved_entity_unsupported")
}

#' Restore only the declared schema without scanning assay rows.
#' @param scope Existing validated table lease.
#' @return Zero-row table with original types and bounded attributes.
#' @noRd
readSavedWorkflowTablePrototype <- function(scope) {
    checkWorkflowTableReadScope(scope)
    metadata <- scope$metadata
    data <- lapply(metadata$column_metadata, function(spec) vector(spec$typeof, 0L))
    for (name in workflowTableStorageColumns(metadata, character())) data[[name]] <- character()
    for (name in unname(metadata$name_columns)) data[[name]] <- character()
    attr(data, "row.names") <- .set_row_names(0L)
    class(data) <- "data.frame"
    data <- restoreSavedWorkflowTable(data, metadata, metadata$columns)
    attr(data, "row.names") <- .set_row_names(0L)
    data
}

#' Own one Arrow cursor, coalescing sparse selections before R conversion.
#' @param scope Validated table scope, which also owns cursor cleanup.
#' @param columns Logical projected columns.
#' @param source_rows Optional increasing physical ordinals.
#' @return Cursor with fetch and close methods; fetch returns at most 10000 rows.
#' @noRd
openWorkflowArrowCursor <- function(scope, columns, source_rows = NULL) {
    checkWorkflowTableReadScope(scope)
    assertWorkflowPhysicalTableScope(scope)
    dataset <- arrow::open_dataset(scope$path, format = "parquet")
    scanner <- reader <- NULL
    cursor <- new.env(parent = emptyenv())
    cursor$closed <- FALSE
    cursor$tail <- NULL
    cursor$id <- basename(tempfile("arrow-cursor-"))
    cursor$close <- function() {
        if (cursor$closed) return(invisible(NULL))
        cursor$closed <- TRUE
        cursor$tail <- NULL
        if (!is.null(reader)) reader$.unsafe_delete()
        if (!is.null(scanner)) scanner$.unsafe_delete()
        dataset$.unsafe_delete()
        if (!is.null(scope$readers) && exists(cursor$id, scope$readers, inherits = FALSE))
            rm(list = cursor$id, envir = scope$readers)
        invisible(NULL)
    }
    ready <- FALSE
    on.exit(if (!ready) cursor$close(), add = TRUE)
    scanner <- arrow::Scanner$create(dataset,
        projection = unique(c(".multischolar_order", workflowTableStorageColumns(scope$metadata, columns))),
        filter = workflowTableRowFilter(source_rows), batch_size = 10000L, use_threads = FALSE,
        fragment_scan_options = arrow::ParquetFragmentScanOptions$create(pre_buffer = FALSE))
    reader <- scanner$ToRecordBatchReader()
    cursor$fetch <- function() {
        if (cursor$closed) return(NULL)
        checkWorkflowTableReadScope(scope)
        pieces <- vector("list", 1024L)
        count <- rows <- bytes <- 0L
        repeat {
            batch <- cursor$tail
            cursor$tail <- NULL
            if (is.null(batch)) batch <- reader$read_next_batch()
            if (is.null(batch)) break
            if (!batch$num_rows) next
            if (count && bytes + batch$nbytes() > 33554432) {
                cursor$tail <- batch
                break
            }
            take <- min(batch$num_rows, 10000L - rows)
            if (take < batch$num_rows) cursor$tail <- batch$Slice(take, batch$num_rows - take)
            count <- count + 1L
            pieces[[count]] <- if (take == batch$num_rows) batch else batch$Slice(0L, take)
            rows <- rows + take
            bytes <- bytes + pieces[[count]]$nbytes()
            if (is.null(source_rows) || rows == 10000L || count == length(pieces) || bytes >= 33554432) break
        }
        if (!count) { cursor$close(); return(NULL) }
        value <- if (count == 1L) as.data.frame(pieces[[1L]]) else
            as.data.frame(do.call(arrow::Table$create, pieces[seq_len(count)]))
        names(value)[names(value) == ".multischolar_order"] <- ".source_order"
        value
    }
    assign(cursor$id, cursor, scope$readers)
    scope$store$io_counts$sequential_streams <- scope$store$io_counts$sequential_streams + 1L
    ready <- TRUE
    cursor
}

#' Restore one bounded cursor batch and account for its delivered rows.
#' @param scope,columns Validated source and logical projection.
#' @param part Physical bounded batch including source ordinals.
#' @return Logical data and source row names.
#' @noRd
restoreWorkflowCursorBatch <- function(scope, columns, part) {
    checkWorkflowTableReadScope(scope)
    if (nrow(part) > 10000L || as.numeric(object.size(part)) > 67108864)
        cli::cli_abort("table batch exceeds read limits", class = "saved_entity_limit")
    row_names <- workflowTableSourceRowNames(part, scope$metadata, part$.source_order)
    data <- restoreSavedWorkflowTable(part, scope$metadata, columns)
    data <- restoreWorkflowTableRowNames(data, scope$metadata, row_names)
    bytes <- as.numeric(object.size(data))
    if (bytes > 67108864) cli::cli_abort("table batch exceeds read limits", class = "saved_entity_limit")
    scope$store$table_reads <- scope$store$table_reads + 1L
    scope$store$table_rows <- scope$store$table_rows + nrow(data)
    scope$store$table_columns <- scope$store$table_columns + ncol(data)
    scope$store$table_bytes <- scope$store$table_bytes + bytes
    list(data = data, row_names = row_names)
}

#' Consume one ordered table cursor in bounded batches without retaining them.
#' @param scope Existing validated table lease.
#' @param columns Ordered logical columns, within the existing 64-column read limit.
#' @param consume Callback receiving data, zero-based offset and source row names.
#' @param source_rows Optional strictly increasing physical row ordinals; selection is pushed down.
#' @return Row count invisibly; the cursor is closed on completion and failure.
#' @noRd
streamSavedWorkflowTable <- function(scope, columns, consume, query = NULL, params = NULL,
    expected_rows = NULL, source_rows = NULL) {
    checkWorkflowTableReadScope(scope)
    metadata <- scope$metadata
    if (!is.function(consume) || !length(columns) || length(columns) > 64L ||
        anyDuplicated(columns) || !all(columns %in% metadata$columns)) {
        cli::cli_abort("invalid table stream", class = "saved_entity_limit")
    }
    physical <- workflowTableStorageColumns(metadata, columns)
    if (!is.null(source_rows)) {
        if (!is.null(query) || !is.numeric(source_rows) || anyNA(source_rows) ||
            any(source_rows < 1 | source_rows > metadata$rows | source_rows != floor(source_rows)) ||
            is.unsorted(source_rows, strictly = TRUE))
            cli::cli_abort("invalid table row selection", class = "saved_entity_limit")
        expected_rows <- length(source_rows)
    }
    if (is.null(expected_rows)) expected_rows <- metadata$rows
    if (!is.null(scope$view)) {
        if (!is.null(query)) cli::cli_abort("physical SQL cannot stream a selected table",
            class = "saved_entity_unsupported")
        if (!expected_rows) {
            consume(readSavedWorkflowTablePrototype(scope)[columns], 0L, character())
        } else {
            for (start in seq.int(1L, expected_rows, by = 10000L)) {
                at <- seq.int(start, min(expected_rows, start + 9999L))
                logical_rows <- if (is.null(source_rows)) at else source_rows[at]
                pieces <- split(logical_rows, cumsum(c(1L, diff(logical_rows) != 1L)))
                delivered <- 0L
                for (piece in pieces) {
                    page <- readSavedWorkflowTable(scope$store, scope$reference, columns,
                        offset = piece[[1L]] - 1L, limit = length(piece), scope = scope)
                    data <- page$data
                    if (!is.null(source_rows))
                        attr(data, "dia_source_order") <- as.integer(logical_rows[
                            seq.int(delivered + 1L, delivered + nrow(data))])
                    consume(data, start - 1L + delivered, page$source_row_names)
                    delivered <- delivered + nrow(data)
                }
            }
        }
        return(invisible(expected_rows))
    }
    if (is.null(query) && is.null(source_rows) && identical(scope$reader, "fastduckdb")) {
        for (offset in seq.int(0L, max(0L, expected_rows - 1L), by = 10000L)) {
            if (!expected_rows) {
                data <- readSavedWorkflowTablePrototype(scope)[columns]
                consume(data, 0L, character())
            } else {
                page <- readSavedWorkflowTable(scope$store, scope$reference, columns,
                    offset = offset, limit = min(10000L, expected_rows - offset),
                    reader = "fastduckdb", scope = scope)
                consume(page$data, offset, page$source_row_names)
            }
        }
        return(invisible(expected_rows))
    }
    if (is.null(query)) {
        cursor <- openWorkflowArrowCursor(scope, columns, source_rows)
        on.exit(cursor$close(), add = TRUE)
        fetch <- cursor$fetch
    } else {
        con <- scope$connection
        if (is.null(con)) cli::cli_abort("SQL stream requires an explicit query lease")
        result <- DBI::dbSendQuery(con, query)
        on.exit(DBI::dbClearResult(result), add = TRUE)
        if (length(params)) DBI::dbBind(result, params)
        scope$store$io_counts$query_streams <- scope$store$io_counts$query_streams + 1L
        fetch <- function() DBI::dbFetch(result, n = 10000L)
    }
    offset <- 0L
    repeat {
        part <- fetch()
        if (is.null(part) || !is.null(query) && !nrow(part)) {
            if (!offset) consume(readSavedWorkflowTablePrototype(scope)[columns], 0L, character())
            break
        }
        if (!nrow(part)) next
        expected_order <- seq.int(offset + 1L, offset + nrow(part))
        if (!is.null(source_rows)) expected_order <- source_rows[expected_order]
        if (is.null(query) && !identical(as.integer(part$.source_order),
            as.integer(expected_order))) {
            cli::cli_abort("sequential table order differs", class = "saved_entity_corrupt")
        }
        restored <- restoreWorkflowCursorBatch(scope, columns, part)
        data <- restored$data
        row_names <- restored$row_names
        restored <- NULL
        if (!is.null(query) || !is.null(source_rows))
            attr(data, "dia_source_order") <- as.integer(part$.source_order)
        consume(data, offset, row_names)
        offset <- offset + nrow(data)
        if (is.function(scope$store$progress)) scope$store$progress("table_read", offset,
            expected_rows, scope$physical_reference$entity_id, "rows")
        if (!nrow(data)) break
        data <- part <- NULL
    }
    if (offset != expected_rows) cli::cli_abort("table stream is incomplete", class = "saved_entity_corrupt")
    invisible(offset)
}

#' Push ordinal selection into Parquet before materializing scientific columns in R.
#' @param rows Increasing physical ordinals, or NULL for the complete input.
#' @return Arrow filter expression with bounded range-expression complexity.
#' @noRd
workflowTableRowFilter <- function(rows) {
    if (is.null(rows)) return(TRUE)
    if (!length(rows)) return(arrow::Expression$scalar(FALSE))
    field <- arrow::Expression$field_ref(".multischolar_order")
    starts <- c(1L, which(diff(rows) != 1L) + 1L)
    if (length(starts) > 1024L) return(arrow::Expression$create("is_in", field,
        options = list(value_set = arrow::Array$create(rows), skip_nulls = FALSE)))
    ends <- c(starts[-1L] - 1L, length(rows))
    filters <- lapply(seq_along(starts), function(i) arrow::Expression$create("and_kleene",
        arrow::Expression$create("greater_equal", field, rows[[starts[[i]]]]),
        arrow::Expression$create("less_equal", field, rows[[ends[[i]]]])))
    while (length(filters) > 1L) {
        filters <- lapply(seq.int(1L, length(filters), by = 2L), function(i) {
            if (i == length(filters)) filters[[i]] else
                arrow::Expression$create("or_kleene", filters[[i]], filters[[i + 1L]])
        })
    }
    filters[[1L]]
}

#' Allocate one declared scientific group, without retaining fetch batches.
#' @param metadata,columns Restoration schema and explicit projection.
#' @param rows Exact complete-group row count supplied by the grouping query.
#' @return Private assembly buffer.
#' @noRd
newWorkflowTableBuffer <- function(metadata, columns, rows) {
    buffer <- new.env(parent = emptyenv())
    buffer$metadata <- metadata
    buffer$columns <- stats::setNames(lapply(metadata$column_metadata[columns], function(spec) {
        if (!is.null(spec$levels)) factor(rep(NA_character_, rows), levels = spec$levels,
            ordered = isTRUE(spec$ordered)) else vector(spec$typeof, rows)
    }), columns)
    for (name in intersect(columns, names(metadata$name_columns))) names(buffer$columns[[name]]) <- character(rows)
    buffer$row_names <- if (isTRUE(metadata$automatic_row_names)) NULL else character(rows)
    buffer$positions <- integer(rows)
    buffer$filled <- 0L
    buffer$rows <- rows
    buffer
}

#' Copy one bounded physical batch into its final scientific-group allocation.
#' @param buffer Existing exact-size buffer.
#' @param physical Physical rows from the validated cursor.
#' @return NULL invisibly.
#' @noRd
fillWorkflowTableBuffer <- function(buffer, physical) {
    metadata <- buffer$metadata
    columns <- names(buffer$columns)
    data <- restoreSavedWorkflowTable(physical, metadata, columns)
    at <- seq.int(buffer$filled + 1L, buffer$filled + nrow(physical))
    assignWorkflowTableBuffer(buffer, data, at, as.integer(physical$.source_order),
        if (!is.null(buffer$row_names)) workflowTableSourceRowNames(physical, metadata, physical$.source_order))
    invisible(NULL)
}

#' Fill a private scientific buffer without copying its full destination per batch.
#' @param buffer Existing exclusively owned assembly buffer.
#' @param data Restored bounded table batch.
#' @param at Unique destination positions supplied by the validated reader.
#' @param source_order Original physical row ordinals for this batch.
#' @param row_names Original row names when the schema requires them.
#' @return NULL invisibly; the buffer owns the completed assignments.
#' @noRd
assignWorkflowTableBuffer <- function(buffer, data, at, source_order, row_names = NULL) {
    # Transfer each vector out of its container before replacement. Nested
    # environment/list subassignment otherwise duplicates the whole cohort.
    columns <- buffer$columns
    buffer$columns <- NULL
    for (name in names(columns)) {
        column <- columns[[name]]
        columns[name] <- list(NULL)
        column[at] <- data[[name]]
        if (name %in% names(buffer$metadata$name_columns)) {
            column_names <- names(column)
            names(column) <- NULL
            column_names[at] <- names(data[[name]])
            names(column) <- column_names
            column_names <- NULL
        }
        columns[[name]] <- column
        column <- NULL
    }
    buffer$columns <- columns
    columns <- NULL
    positions <- buffer$positions
    buffer$positions <- NULL
    positions[at] <- source_order
    buffer$positions <- positions
    positions <- NULL
    if (!is.null(buffer$row_names)) {
        stored_names <- buffer$row_names
        buffer$row_names <- NULL
        stored_names[at] <- row_names
        buffer$row_names <- stored_names
        stored_names <- NULL
    }
    buffer$filled <- buffer$filled + length(at)
    invisible(NULL)
}

#' Deliver the complete group and release its assembly owner.
#' @param buffer Filled exact-size buffer.
#' @return Original table classes/attributes and an explicit source-row mapping.
#' @noRd
finishWorkflowTableBuffer <- function(buffer) {
    if (buffer$filled != buffer$rows) cli::cli_abort("scientific group is incomplete", class = "saved_entity_corrupt")
    metadata <- buffer$metadata
    columns <- buffer$columns
    for (name in names(columns)) {
        spec <- metadata$column_metadata[[name]]
        if (is.null(spec$levels) && length(spec$class)) class(columns[[name]]) <- spec$class
        for (attribute in names(spec$attributes)) attr(columns[[name]], attribute) <- spec$attributes[[attribute]]
        columns[[name]] <- restoreWorkflowAttributeOrder(columns[[name]], spec$attribute_order)
    }
    for (attribute in names(metadata$attributes)) attr(columns, attribute) <- metadata$attributes[[attribute]]
    row_names <- buffer$row_names
    if (!is.null(row_names) && identical(metadata$row_index_type, "integer")) row_names <- as.integer(row_names)
    attr(columns, "row.names") <- row_names %||% .set_row_names(buffer$rows)
    class(columns) <- metadata$class
    columns <- restoreWorkflowAttributeOrder(columns, metadata$attribute_order)
    attr(columns, "workflow_source_order") <- buffer$positions
    buffer$columns <- buffer$positions <- buffer$row_names <- NULL
    columns
}

#' Deliver physically contiguous scientific groups from one forward cursor.
#' @param scope,columns Existing lease and logical projection, including group.
#' @param group Group identity column.
#' @param directory Validated complete-group boundaries in physical order.
#' @param dependencies Existing index and source references.
#' @return Complete-group iterator retaining at most one read batch and one group output.
#' @noRd
openWorkflowContiguousGroups <- function(scope, columns, group, directory, dependencies) {
    reader <- openWorkflowArrowCursor(scope, columns)
    iterator <- new.env(parent = emptyenv())
    iterator$dependencies <- dependencies
    iterator$at <- 1L
    iterator$loaded <- iterator$delivered <- 0L
    iterator$pending <- iterator$source <- iterator$row_names <- NULL
    iterator$pending_at <- 1L
    iterator$closed <- FALSE
    iterator$close <- function() {
        iterator$closed <- TRUE
        iterator$pending <- iterator$source <- iterator$row_names <- NULL
        reader$close()
        invisible(NULL)
    }
    iterator$next_group <- function() {
        if (iterator$closed) return(NULL)
        checkWorkflowTableReadScope(scope)
        if (iterator$at > nrow(directory)) {
            if (iterator$delivered != scope$metadata$rows)
                cli::cli_abort("contiguous groups are incomplete", class = "saved_entity_corrupt")
            iterator$close()
            return(NULL)
        }
        sizes <- directory$size[seq.int(iterator$at, nrow(directory))]
        count <- max(1L, sum(cumsum(as.double(sizes)) <= 10000L))
        rows <- sum(sizes[seq_len(count)])
        if (rows * workflowProjectedRowBytes(scope$metadata, columns) > 67108864)
            cli::cli_abort("complete scientific group exceeds 64 MiB", class = "saved_entity_limit")
        buffer <- newWorkflowTableBuffer(scope$metadata, columns, rows)
        while (buffer$filled < rows) {
            if (is.null(iterator$pending) || iterator$pending_at > nrow(iterator$pending)) {
                physical <- reader$fetch()
                if (is.null(physical) || !identical(as.integer(physical$.source_order),
                    seq.int(iterator$loaded + 1L, length.out = nrow(physical))))
                    cli::cli_abort("contiguous source is incomplete", class = "saved_entity_corrupt")
                restored <- restoreWorkflowCursorBatch(scope, columns, physical)
                expected <- findInterval(physical$.source_order, directory$first_order)
                if (!identical(match(restored$data[[group]], directory$value), as.integer(expected)))
                    cli::cli_abort("contiguous group identity differs", class = "saved_entity_corrupt")
                iterator$pending <- restored$data
                iterator$source <- as.integer(physical$.source_order)
                iterator$row_names <- restored$row_names
                iterator$loaded <- iterator$loaded + nrow(physical)
                iterator$pending_at <- 1L
                restored <- physical <- NULL
            }
            take <- min(rows - buffer$filled, nrow(iterator$pending) - iterator$pending_at + 1L)
            from <- seq.int(iterator$pending_at, length.out = take)
            to <- seq.int(buffer$filled + 1L, length.out = take)
            for (column in columns) {
                buffer$columns[[column]][to] <- iterator$pending[[column]][from]
                if (column %in% names(scope$metadata$name_columns))
                    names(buffer$columns[[column]])[to] <- names(iterator$pending[[column]])[from]
            }
            buffer$positions[to] <- iterator$source[from]
            if (!is.null(buffer$row_names)) buffer$row_names[to] <- iterator$row_names[from]
            buffer$filled <- buffer$filled + take
            iterator$pending_at <- iterator$pending_at + take
            if (as.numeric(object.size(buffer$columns)) > 67108864)
                cli::cli_abort("scientific group exceeds 64 MiB", class = "saved_entity_limit")
        }
        iterator$at <- iterator$at + count
        iterator$delivered <- iterator$delivered + rows
        if (is.function(scope$store$progress)) scope$store$progress("input_groups", iterator$at - 1L,
            nrow(directory), scope$physical_reference$entity_id, paste0("groups:", group))
        value <- finishWorkflowTableBuffer(buffer)
        if (as.numeric(object.size(value)) > 67108864)
            cli::cli_abort("scientific group exceeds 64 MiB", class = "saved_entity_limit")
        value
    }
    iterator
}

#' Open one projected cursor for complete independence groups, without copying the assay.
#' @param scope Existing SQL action lease.
#' @param columns Required scientific columns, at most 64 per fetch.
#' @param group Existing grouping column; NULL values remain one complete group.
#' @param group_order First-observation order, or native key order for ordered producers.
#' @return Iterator with next_group() and close(); callers must close on failure.
#' @noRd
openWorkflowTableGroups <- function(scope, columns, group, group_order = "first") {
    checkWorkflowTableReadScope(scope)
    if (!is.character(columns) || !length(columns) || anyDuplicated(columns) ||
        length(columns) > 64L || !all(c(columns, group) %in% scope$metadata$columns) ||
        !group_order %in% c("first", "key"))
        cli::cli_abort("invalid scientific group projection", class = "saved_entity_limit")
    index <- saveWorkflowGroupIndex(scope, group)
    directory <- readWorkflowIndexDirectory(scope$store, index$value$groups)
    required <- c("group_id", "value", "size", "start", "first_order")
    if (!is.data.frame(directory) || !identical(names(directory), required) ||
        anyNA(directory[c("group_id", "size", "start", "first_order")]) ||
        !identical(directory$group_id, seq_len(nrow(directory))) || any(directory$size < 1L) ||
        anyDuplicated(directory$value) || sum(as.double(directory$size)) != scope$metadata$rows ||
        !isTRUE(all.equal(directory$start, if (nrow(directory))
            1L + cumsum(c(0L, utils::head(directory$size, -1L))) else integer())))
        cli::cli_abort("group directory is invalid", class = "saved_entity_corrupt")
    if (group_order == "key" && nrow(directory)) directory <- directory[
        order(directory$value, method = "radix", na.last = TRUE), , drop = FALSE]
    contiguous_starts <- if (nrow(directory)) 1 + cumsum(c(0, utils::head(directory$size, -1L))) else numeric()
    if (is.null(scope$view) && nrow(directory) && group %in% columns &&
        identical(as.integer(directory$first_order), as.integer(contiguous_starts))) {
        return(openWorkflowContiguousGroups(scope, columns, group, directory,
            list(index$reference, index$value$groups, index$value$positions, scope$physical_reference)))
    }
    iterator <- new.env(parent = emptyenv())
    iterator$dependencies <- list(index$reference, index$value$groups, index$value$positions,
        scope$physical_reference)
    iterator$at <- 1L
    iterator$delivered_groups <- 0L
    iterator$cohort <- NULL
    iterator$cohort_at <- 1L
    iterator$cohort_groups <- integer()
    iterator$closed <- FALSE
    iterator$empty_delivered <- FALSE
    iterator$close <- function() {
        iterator$cohort <- NULL
        iterator$closed <- TRUE
        invisible(NULL)
    }
    row_bytes <- workflowProjectedRowBytes(scope$metadata, columns)
    iterator$next_group <- function() {
        if (iterator$closed) return(NULL)
        checkWorkflowTableReadScope(scope)
        if (!nrow(directory)) {
            iterator$close()
            value <- readSavedWorkflowTablePrototype(scope)[columns]
            attr(value, "workflow_source_order") <- integer()
            return(value)
        }
        if (is.null(iterator$cohort) || iterator$cohort_at > length(iterator$cohort_groups)) {
            iterator$cohort <- NULL
            if (iterator$at > nrow(directory)) { iterator$close(); return(NULL) }
            remaining <- seq.int(iterator$at, nrow(directory))
            sizes <- directory$size[remaining]
            if (sizes[[1L]] * row_bytes > 67108864)
                cli::cli_abort("complete scientific group exceeds 64 MiB", class = "saved_entity_limit")
            take <- max(1L, sum(cumsum(as.double(sizes)) * row_bytes <= 268435456))
            selected <- remaining[seq_len(take)]
            if (any(directory$size[selected] * row_bytes > 67108864))
                cli::cli_abort("complete scientific group exceeds 64 MiB", class = "saved_entity_limit")
            iterator$cohort <- readWorkflowGroupCohort(scope, columns, index$value$positions,
                directory[selected, , drop = FALSE])
            iterator$cohort_groups <- directory$size[selected]
            iterator$cohort_at <- 1L
            iterator$row_at <- 1L
            iterator$at <- iterator$at + take
        }
        remaining <- seq.int(iterator$cohort_at, length(iterator$cohort_groups))
        count <- max(1L, sum(cumsum(iterator$cohort_groups[remaining]) <= 10000L))
        rows <- sum(iterator$cohort_groups[remaining[seq_len(count)]])
        at <- seq.int(iterator$row_at, length.out = rows)
        positions <- attr(iterator$cohort, "workflow_source_order")[at]
        value <- iterator$cohort[at, , drop = FALSE]
        if (isTRUE(scope$metadata$automatic_row_names)) attr(value, "row.names") <- .set_row_names(rows)
        attr(value, "workflow_source_order") <- positions
        if (as.numeric(object.size(value)) > 67108864)
            cli::cli_abort("scientific group exceeds 64 MiB", class = "saved_entity_limit")
        iterator$cohort_at <- iterator$cohort_at + count
        iterator$row_at <- iterator$row_at + rows
        iterator$delivered_groups <- iterator$delivered_groups + count
        if (is.function(scope$store$progress)) scope$store$progress("input_groups", iterator$delivered_groups,
            nrow(directory), scope$physical_reference$entity_id, paste0("groups:", group))
        value
    }
    iterator
}

#' Complete an output batch from its declared schema, retaining original attributes.
#' @param data Bounded provided columns; omitted columns represent newly generated missing values.
#' @param metadata Final scientific restoration contract.
#' @return One bounded complete row batch.
#' @noRd
completeWorkflowTableBatch <- function(data, metadata) {
    out <- stats::setNames(lapply(metadata$columns, function(name) {
        spec <- metadata$column_metadata[[name]]
        value <- data[[name]]
        if (is.null(value)) value <- if (!is.null(spec$levels)) factor(rep(NA_character_, nrow(data)),
            levels = spec$levels, ordered = isTRUE(spec$ordered)) else rep(switch(spec$typeof,
                integer = NA_integer_, double = NA_real_, logical = NA, character = NA_character_), nrow(data))
        if (is.null(spec$levels) && length(spec$class)) class(value) <- spec$class
        for (attribute in names(spec$attributes)) attr(value, attribute) <- spec$attributes[[attribute]]
        if (name %in% names(metadata$name_columns)) {
            if (is.null(names(value))) names(value) <- rep(NA_character_, length(value))
        } else names(value) <- NULL
        restoreWorkflowAttributeOrder(value, spec$attribute_order)
    }), metadata$columns)
    for (attribute in names(metadata$attributes)) attr(out, attribute) <- metadata$attributes[[attribute]]
    attr(out, "row.names") <- if (isTRUE(metadata$automatic_row_names)) .set_row_names(nrow(data)) else
        .row_names_info(data, type = 0L)
    class(out) <- metadata$class
    restoreWorkflowAttributeOrder(out, metadata$attribute_order)
}

#' Recognize an unchanged table from a complete row map and logical schema.
#' @param scope Existing source lease.
#' @param map_path Saved observed-row map with no generated rows.
#' @param metadata Proposed output restoration contract.
#' @param source_column Stable source ordinal in the map.
#' @return TRUE only when every source row and scientific attribute is unchanged.
#' @noRd
workflowTableChangesAreIdentity <- function(scope, map_path, metadata, source_column) {
    map_scope <- newWorkflowTableReadScope(scope$store, map_path, reader = "arrow", metadata_only = TRUE)
    mapping <- map_scope$metadata
    releaseWorkflowTableReadScope(map_scope)
    if (!identical(mapping$columns, source_column) || mapping$rows != scope$metadata$rows) return(FALSE)
    logical_metadata <- function(value) {
        value[c("schema", "codec", "row_names", "row_name_column")] <- NULL
        value$name_columns <- names(value$name_columns)
        value
    }
    if (!identical(logical_metadata(scope$metadata), logical_metadata(metadata))) return(FALSE)
    cursor <- openWorkflowOrderedRows(scope$store, map_path, source_column, source_column)
    on.exit(cursor$close(), add = TRUE)
    offset <- 0L
    repeat {
        part <- cursor$next_batch()
        if (is.null(part)) break
        expected <- offset + seq_len(nrow(part))
        if (!identical(as.integer(part[[source_column]]), expected)) return(FALSE)
        offset <- offset + nrow(part)
    }
    offset == scope$metadata$rows
}

#' Merge a narrow row/change map with one sequential source pass into final Parquet.
#' @param scope Existing validated source lease.
#' @param map_path,generated_path Saved action-owned tables; generated rows may be absent.
#' @param metadata Final scientific schema and exact row count.
#' @param path Fresh final staging file.
#' @param source_column Stable source-row column in the map.
#' @param generated_order Declared ordering columns for appended rows.
#' @return NULL invisibly; only bounded cursor batches overlap the writer.
#' @noRd
writeWorkflowTableChanges <- function(scope, map_path, generated_path, metadata, path,
    source_column, generated_order = character(), consume = NULL) {
    map_scope <- newWorkflowTableReadScope(scope$store, map_path, reader = "arrow", metadata_only = TRUE)
    map_metadata <- map_scope$metadata
    releaseWorkflowTableReadScope(map_scope)
    map <- openWorkflowOrderedRows(scope$store, map_path, map_metadata$columns, source_column)
    on.exit(map$close(), add = TRUE)
    writer <- openWorkflowTableWriter(path, metadata, scope$store$progress)
    on.exit(abortWorkflowTableWriter(writer), add = TRUE)
    next_map <- map$next_batch
    mapping <- next_map()
    at <- 1L
    previous <- 0L
    changes <- intersect(setdiff(map_metadata$columns,
        c(source_column, generated_order)), metadata$columns)
    streamWorkflowTableRows(scope, scope$metadata$columns, function(data, offset, source_names) {
        if (!nrow(data)) return(invisible(NULL))
        out <- completeWorkflowTableBatch(data, metadata)
        keep <- rep(FALSE, nrow(data))
        repeat {
            if (is.null(mapping) || mapping[[source_column]][[at]] > offset + nrow(data)) break
            stop_at <- min(nrow(mapping), max(which(mapping[[source_column]] <= offset + nrow(data))))
            selected <- seq.int(at, stop_at)
            positions <- as.integer(mapping[[source_column]][selected])
            if (anyNA(positions) || is.unsorted(c(previous, positions), strictly = TRUE) ||
                positions[[1L]] <= offset) cli::cli_abort("invalid source-row change map", class = "saved_entity_corrupt")
            local <- positions - offset
            keep[local] <- TRUE
            for (name in changes) {
                out[[name]][local] <- mapping[[name]][selected]
                if (name %in% names(metadata$name_columns)) names(out[[name]])[local] <- names(mapping[[name]])[selected]
            }
            previous <<- utils::tail(positions, 1L)
            at <<- stop_at + 1L
            if (at > nrow(mapping)) {
                mapping <<- next_map()
                at <<- 1L
            }
        }
        if (!is.null(consume)) consume(out[keep, , drop = FALSE])
        appendWorkflowTableWriter(writer, out, which(keep))
        invisible(NULL)
    })
    if (!is.null(mapping)) cli::cli_abort("change map refers beyond the source", class = "saved_entity_corrupt")
    map$close()
    if (!is.null(generated_path)) {
        generated_scope <- newWorkflowTableReadScope(scope$store, generated_path, reader = "arrow", metadata_only = TRUE)
        generated_metadata <- generated_scope$metadata
        releaseWorkflowTableReadScope(generated_scope)
        generated <- openWorkflowOrderedRows(scope$store, generated_path, generated_metadata$columns,
            c(generated_order, ".multischolar_order"))
        on.exit(generated$close(), add = TRUE)
        repeat {
            data <- generated$next_batch()
            if (is.null(data)) break
            data <- completeWorkflowTableBatch(data, metadata)
            if (!is.null(consume)) consume(data)
            appendWorkflowTableWriter(writer, data)
        }
    }
    if (is.null(writer$output)) appendWorkflowTableWriter(writer,
        completeWorkflowTableBatch(readSavedWorkflowTablePrototype(scope), metadata))
    finishWorkflowTableWriter(writer)
    invisible(NULL)
}

#' Stream wide output rows from bounded column/range reads.
#' @param scope Existing source lease, which may also own a SQL result cursor.
#' @param columns Ordered scientific output columns.
#' @param consume Callback receiving one complete bounded row batch.
#' @return Total consumed rows invisibly.
#' @noRd
streamWorkflowTableRows <- function(scope, columns, consume) {
    if (length(columns) <= 64L) return(streamSavedWorkflowTable(scope, columns, consume))
    projected <- newWorkflowTableReadScope(scope$store, scope$reference, reader = "fastduckdb")
    on.exit(releaseWorkflowTableReadScope(projected), add = TRUE)
    blocks <- split(columns, ceiling(seq_along(columns) / 64L))
    for (offset in seq.int(0L, max(0L, scope$metadata$rows - 1L), by = 10000L)) {
        values <- list()
        for (block in blocks) {
            page <- readSavedWorkflowTable(scope$store, scope$reference, block,
                offset = offset, limit = 10000L, scope = projected, reader = "fastduckdb")
            values[block] <- unclass(page$data)[block]
        }
        attr(values, "row.names") <- .row_names_info(page$data, type = 0L)
        class(values) <- "data.frame"
        metadata <- scope$metadata
        metadata$columns <- columns
        data <- completeWorkflowTableBatch(values, metadata)
        if (as.numeric(object.size(data)) > 67108864) cli::cli_abort("wide output batch exceeds byte limits")
        consume(data, offset, page$source_row_names)
    }
    invisible(scope$metadata$rows)
}

#' Observe bounded scratch-file sizes without reading or capturing their contents.
#' @param scope Current table lease.
#' @return NULL invisibly; this is sampled disk size, not cumulative physical I/O.
#' @noRd
observeWorkflowTableScratch <- function(scope) {
    root <- scope$spill_directory
    if (is.null(root) || !dir.exists(root)) return(invisible(NULL))
    files <- list.files(root, full.names = TRUE, recursive = FALSE, all.files = TRUE, no.. = TRUE)
    if (length(files) > scope$store$max_entries) cli::cli_abort("table scratch inventory exceeds limits")
    sizes <- file.info(files)$size
    scope$store$io_counts$max_scratch_bytes <- max(scope$store$io_counts$max_scratch_bytes %||% 0,
        sum(sizes, na.rm = TRUE))
    invisible(NULL)
}

#' Identify an immutable logical table component stored with the existing RDS codec.
#' @param reference Saved descriptor, or another value.
#' @return TRUE only for the explicit bounded table-component class.
#' @noRd
isWorkflowTableComponent <- function(reference) {
    isSavedWorkflowReference(reference) && identical(reference$kind, "result_table") &&
        identical(reference$payload$codec, "rds") &&
        identical(reference$metadata$class, "WorkflowTableComponent")
}

#' Identify physical and logical tables without hydrating either payload.
#' @param reference Saved descriptor, or another value.
#' @return TRUE for supported table references.
#' @noRd
isWorkflowTableReference <- function(reference) {
    isSavedWorkflowReference(reference) && (identical(reference$payload$codec, "parquet_duckdb") ||
        isWorkflowTableComponent(reference))
}

#' Validate a logical restoration overlay without changing physical column storage.
#' @param metadata Effective restoration metadata.
#' @param physical Immutable Parquet restoration metadata.
#' @return TRUE; value/schema changes require a new physical table.
#' @noRd
checkWorkflowTableRestoration <- function(metadata, physical) {
    mutable <- c("class", "attributes", "attribute_order", "columns", "column_metadata", "name_columns")
    if (!is.list(metadata) || !identical(names(metadata), names(physical)) ||
        !identical(metadata[setdiff(names(metadata), mutable)], physical[setdiff(names(physical), mutable)]) ||
        !is.character(metadata$columns) || !length(metadata$columns) || anyDuplicated(metadata$columns) ||
        !all(metadata$columns %in% physical$columns) ||
        !identical(metadata$column_metadata, physical$column_metadata[metadata$columns]) ||
        !identical(metadata$name_columns, physical$name_columns[intersect(names(physical$name_columns), metadata$columns)]) ||
        !any(vapply(list("data.frame", c("tbl_df", "tbl", "data.frame"),
            c("spec_tbl_df", "tbl_df", "tbl", "data.frame")), identical, logical(1), metadata$class)) ||
        as.numeric(object.size(metadata)) > 8388608) {
        cli::cli_abort("logical table metadata changes physical storage", class = "saved_entity_corrupt")
    }
    checkWorkflowTableAttributes(metadata$attributes)
    if (!is.null(metadata$attribute_order) &&
        !setequal(metadata$attribute_order, c("names", "row.names", "class", names(metadata$attributes))))
        cli::cli_abort("logical table attribute order differs", class = "saved_entity_corrupt")
    invisible(TRUE)
}

#' Read only a bounded logical binding; physical assay bytes remain untouched.
#' @param store Open entity store.
#' @param reference Physical or logical table descriptor.
#' @return Canonical physical reference and effective metadata, or NULL for physical data.
#' @noRd
readWorkflowTableComponent <- function(store, reference) {
    if (!isWorkflowTableComponent(reference)) return(NULL)
    if (reference$payload$bytes > 8388608)
        cli::cli_abort("logical table component exceeds metadata limits", class = "saved_entity_limit")
    value <- hydrateSavedWorkflowEntity(store, reference)
    version <- value$schema
    if (!identical(version, "workflow.table_component.v1") &&
        !identical(version, "workflow.table_component.v2"))
        cli::cli_abort("unknown logical table version", class = "saved_entity_unsupported")
    fields <- if (identical(version, "workflow.table_component.v1"))
        c("schema", "source", "metadata_raw", "indexes") else
        c("schema", "source", "rows", "columns", "metadata_raw", "view_id", "indexes")
    if (!identical(names(value), fields) ||
        !isSavedWorkflowReference(value$source) ||
        !identical(value$source$payload$codec, "parquet_duckdb") ||
        !identical(value$source$owner, store$identity) || !is.raw(value$metadata_raw) ||
        length(value$metadata_raw) > 8388608 || !is.list(value$indexes)) {
        cli::cli_abort("logical table component is invalid", class = "saved_entity_corrupt")
    }
    metadata <- tryCatch(unserialize(value$metadata_raw), error = function(error)
        cli::cli_abort("logical table metadata cannot be decoded", class = "saved_entity_corrupt", parent = error))
    if (!is.list(metadata) || !identical(as.integer(c(metadata$rows, length(metadata$columns))), as.integer(reference$metadata$shape)))
        cli::cli_abort("logical table dimensions differ", class = "saved_entity_corrupt")
    if (identical(version, "workflow.table_component.v1"))
        return(list(source = value$source, metadata = metadata, indexes = value$indexes))
    validateWorkflowTableView(store, reference, value, metadata)
    list(source = value$source, metadata = metadata, indexes = value$indexes, view = value)
}

#' Hash the exact stored logical view tuple without reserializing metadata.
#' @param source Physical table descriptor.
#' @param rows Exact all/mapped row definition.
#' @param columns Ordered output/source/cast mapping.
#' @param metadata_raw Serialized logical restoration metadata bytes.
#' @return Stable SHA256 across supported R readers.
#' @noRd
workflowTableViewId <- function(source, rows, columns, metadata_raw) {
    count <- function(x) is.numeric(x) && length(x) == 1L && !is.na(x) &&
        is.finite(x) && x >= 0 && x == floor(x) && x <= .Machine$integer.max
    if (!isSavedWorkflowReference(source) || !is.list(rows) ||
        !identical(names(rows), c("mode", "map", "count")) ||
        !is.character(rows$mode) || length(rows$mode) != 1L ||
        is.na(rows$mode) || !rows$mode %in% c("all", "mapped") ||
        !is.null(rows$map) && !isSavedWorkflowReference(rows$map) ||
        !count(rows$count) ||
        !is.data.frame(columns) || !identical(names(columns), c("output", "source", "cast")) ||
        !nrow(columns) || !is.raw(metadata_raw))
        cli::cli_abort("invalid logical view identity input", class = "saved_entity_invalid")
    tokens <- c("workflow.table_component.v2", source$entity_id, rows$mode,
        if (is.null(rows$map)) "" else rows$map$entity_id,
        sprintf("%.0f", rows$count), sprintf("%.0f", nrow(columns)),
        as.vector(t(as.matrix(columns[c("output", "source", "cast")]))),
        digest::digest(metadata_raw, algo = "sha256", serialize = FALSE))
    if (!is.character(tokens) || anyNA(tokens) || !all(validUTF8(tokens)) ||
        any(vapply(tokens, \(token) !identical(charToRaw(token), charToRaw(enc2utf8(token))),
            logical(1))))
        cli::cli_abort("logical view tokens must be valid UTF-8", class = "saved_entity_invalid")
    framed <- do.call(c, lapply(tokens, \(token) c(
        charToRaw(as.character(length(charToRaw(token)))), as.raw(58L), charToRaw(token))))
    digest::digest(framed, algo = "sha256", serialize = FALSE)
}

#' Validate one physical row-map table through bounded reads only.
#' @param store Owner of the immutable map.
#' @param reference Physical two-double-column map descriptor.
#' @param rows Expected logical output rows.
#' @param source_rows Physical source rows and exact uniqueness domain.
#' @return TRUE after dense logical and unique source ordinals are proven.
#' @noRd
validateWorkflowTableRowMap <- function(store, reference, rows, source_rows) {
    if (!isWorkflowTableReference(reference) || isWorkflowTableComponent(reference) ||
        !identical(reference$owner, store$identity) ||
        as.double(source_rows) * 4 > 268435456)
        cli::cli_abort("row map exceeds supported ownership or memory", class = "saved_entity_unsupported")
    scope <- newWorkflowTableReadScope(store, reference, reader = "arrow")
    on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
    names <- c("logical_ordinal", "source_ordinal")
    specs <- scope$metadata$column_metadata
    if (!identical(scope$metadata$columns, names) || !identical(scope$metadata$key_column, names[[1L]]) ||
        scope$metadata$rows != rows || rows > source_rows ||
        !all(vapply(specs, \(spec) identical(spec$typeof, "double") &&
            identical(spec$class, "numeric") && is.null(spec$levels) &&
            !length(spec$attributes), logical(1))))
        cli::cli_abort("row map schema differs", class = "saved_entity_corrupt")
    seen <- rep(FALSE, source_rows)
    streamSavedWorkflowTable(scope, names, function(data, offset, row_names) {
        expected <- as.double(offset + seq_len(nrow(data)))
        physical <- data$source_ordinal
        if (!identical(data$logical_ordinal, expected) || anyNA(physical) ||
            any(!is.finite(physical) | physical != floor(physical) |
                physical < 1 | physical > source_rows))
            cli::cli_abort("row map ordinals differ", class = "saved_entity_corrupt")
        at <- as.integer(physical)
        if (anyDuplicated(at) || any(seen[at]))
            cli::cli_abort("row map repeats a source row", class = "saved_entity_corrupt")
        seen[at] <<- TRUE
    })
    invisible(TRUE)
}

#' Validate a v2 view against its immutable physical source and map.
#' @param store,reference,view,metadata Owned descriptor, payload and restoration.
#' @return TRUE with no scientific assay materialization.
#' @noRd
validateWorkflowTableView <- function(store, reference, view, metadata) {
    rows <- view$rows
    columns <- view$columns
    if (!is.list(rows) || !identical(names(rows), c("mode", "map", "count")) ||
        !identical(rows$mode, "all") && !identical(rows$mode, "mapped") ||
        !is.numeric(rows$count) || length(rows$count) != 1L || is.na(rows$count) ||
        !is.finite(rows$count) || rows$count < 0 || rows$count != floor(rows$count) ||
        rows$count > .Machine$integer.max ||
        !is.null(rows$map) && !isSavedWorkflowReference(rows$map) ||
        !is.data.frame(columns) || !identical(names(columns), c("output", "source", "cast")) ||
        !nrow(columns) || !all(vapply(columns, is.character, logical(1))) ||
        anyNA(columns) || any(!nzchar(unlist(columns, use.names = FALSE))) ||
        anyDuplicated(columns$output) || anyDuplicated(columns$source) ||
        !all(columns$cast %in% c("identity", "factor_to_character")) ||
        !is.character(view$view_id) || length(view$view_id) != 1L ||
        !identical(view$view_id, workflowTableViewId(view$source, rows, columns, view$metadata_raw)) ||
        !is.list(view$indexes) || length(view$indexes) > 32L ||
        length(view$indexes) && (is.null(names(view$indexes)) || anyNA(names(view$indexes)) ||
            anyDuplicated(names(view$indexes)) || !all(nzchar(names(view$indexes)))))
        cli::cli_abort("logical view contract differs", class = "saved_entity_corrupt")
    required <- c("schema", "codec", "columns", "rows", "key_column", "key_type",
        "class", "attributes", "column_metadata", "automatic_row_names",
        "row_name_column", "name_columns")
    optional <- c("row_index_type", "attribute_order")
    if (!identical(sort(names(metadata)), sort(c(required,
        intersect(optional, names(metadata))))) ||
        !identical(metadata$codec, paste0("parquet_duckdb_v", metadata$schema)) ||
        !is.character(metadata$key_column) || length(metadata$key_column) != 1L ||
        !metadata$key_column %in% c(".rownames", columns$output) ||
        !is.character(metadata$key_type) || length(metadata$key_type) != 1L ||
        !metadata$key_type %in% c("character", "integer", "double") ||
        !is.list(metadata$attributes) || !is.list(metadata$column_metadata) ||
        !identical(names(metadata$column_metadata), columns$output) ||
        !is.logical(metadata$automatic_row_names) ||
        length(metadata$automatic_row_names) != 1L || is.na(metadata$automatic_row_names) ||
        !is.character(metadata$name_columns) || anyNA(metadata$name_columns) ||
        anyDuplicated(names(metadata$name_columns)) ||
        !all(names(metadata$name_columns) %in% columns$output) ||
        anyDuplicated(unname(metadata$name_columns)) ||
        !all(grepl("^\\.multischolar_names_[0-9]+$", metadata$name_columns)) ||
        isTRUE(metadata$automatic_row_names) && !is.null(metadata$row_name_column) ||
        !isTRUE(metadata$automatic_row_names) &&
            !identical(metadata$row_name_column, if (metadata$key_column == ".rownames")
                ".multischolar_row_key" else ".multischolar_row_names") ||
        !is.null(metadata$row_index_type) &&
            (!is.character(metadata$row_index_type) || length(metadata$row_index_type) != 1L ||
                !metadata$row_index_type %in% c("integer", "character")) ||
        !is.null(metadata$attribute_order) &&
            (!is.character(metadata$attribute_order) || anyNA(metadata$attribute_order) ||
                anyDuplicated(metadata$attribute_order) ||
                !setequal(metadata$attribute_order,
                    c("names", "row.names", "class", names(metadata$attributes)))))
        cli::cli_abort("logical restoration metadata differs", class = "saved_entity_corrupt")
    checkWorkflowTableAttributes(metadata$attributes)
    if (length(view$indexes)) for (index in view$indexes) {
        if (!isSavedWorkflowReference(index) || !identical(index$owner, store$identity))
            cli::cli_abort("logical index belongs to another source", class = "saved_entity_corrupt")
    }
    validateSavedWorkflowReference(store, view$source, payload = TRUE)
    physical <- readSavedWorkflowTableMetadata(file.path(store$directory, view$source$payload$relative_path))
    if (rows$count > physical$rows || !identical(metadata$columns, columns$output) ||
        metadata$rows != rows$count || !metadata$schema %in% c(2L, 3L) ||
        !all(columns$source %in% physical$columns) ||
        !any(vapply(list("data.frame", c("tbl_df", "tbl", "data.frame"),
            c("spec_tbl_df", "tbl_df", "tbl", "data.frame")), identical, logical(1), metadata$class)) ||
        as.numeric(object.size(metadata)) > 8388608)
        cli::cli_abort("logical metadata differs from physical source", class = "saved_entity_corrupt")
    for (at in seq_len(nrow(columns))) {
        source <- columns$source[[at]]
        target <- columns$output[[at]]
        actual <- metadata$column_metadata[[target]]
        original <- physical$column_metadata[[source]]
        if (columns$cast[[at]] == "identity") {
            if (!identical(actual, original) ||
                (target %in% names(metadata$name_columns)) !=
                    (source %in% names(physical$name_columns)))
                cli::cli_abort("logical column metadata differs", class = "saved_entity_corrupt")
        } else {
            if (is.null(original$levels) || source %in% names(physical$name_columns) ||
                !identical(actual$typeof, "character") || !identical(actual$class, "character") ||
                !is.null(actual$levels) || isTRUE(actual$ordered) ||
                length(actual$attributes) || target %in% names(metadata$name_columns))
                cli::cli_abort("factor-to-character metadata differs", class = "saved_entity_corrupt")
        }
    }
    if (rows$mode == "all") {
        if (!is.null(rows$map) || rows$count != physical$rows)
            cli::cli_abort("all-row view count differs", class = "saved_entity_corrupt")
    } else {
        if (is.null(rows$map) || identical(rows$map$entity_id, view$source$entity_id))
            cli::cli_abort("mapped view has no row map", class = "saved_entity_corrupt")
        validateWorkflowTableRowMap(store, rows$map, rows$count, physical$rows)
    }
    invisible(TRUE)
}

#' Save a metadata-only table change without duplicating the physical assay.
#' @param scope Existing metadata/read lease.
#' @param metadata Effective restoration metadata with unchanged physical schema.
#' @param request_id,identity Existing durable request and scientific identity.
#' @return Existing table or a small immutable logical-component descriptor.
#' @noRd
saveWorkflowTableComponent <- function(scope, metadata, request_id, identity) {
    checkWorkflowTableReadScope(scope)
    if (!is.null(scope$view)) {
        checkWorkflowTableRestoration(metadata, scope$metadata)
        columns <- metadata$columns
        return(saveWorkflowTableSelection(scope, NULL,
            data.frame(output = columns, source = columns,
                cast = rep("identity", length(columns))), metadata, request_id, identity))
    }
    checkWorkflowTableRestoration(metadata, scope$physical_metadata)
    if (identical(metadata, scope$metadata)) return(scope$reference)
    if (identical(metadata, scope$physical_metadata)) return(scope$physical_reference)
    identity$input <- scope$physical_reference$entity_id
    identity$design <- "physical-table"
    identity$parameters <- digest::digest(metadata, algo = "sha256")
    identity$selectors <- "logical-table-restoration-v1"
    hit <- resolveSavedWorkflowEntity(scope$store, "result_table", identity, hydrate = FALSE)
    if (identical(hit$status, "hit")) return(hit$reference)
    value <- structure(list(schema = "workflow.table_component.v1", source = scope$physical_reference,
        metadata_raw = serialize(metadata, NULL, version = 3), indexes = list()), class = "WorkflowTableComponent")
    if (as.numeric(object.size(value)) > 8388608)
        cli::cli_abort("logical table component exceeds metadata limits", class = "saved_entity_limit")
    saveWorkflowEntity(scope$store, value, request_id, identity,
        list(class = class(value), shape = as.integer(c(metadata$rows, length(metadata$columns))), factors = list()),
        kind = "result_table", codec = "rds")
}

#' Save a bounded v2 row/column selection over the canonical physical table.
#' @param scope Existing logical or physical table lease.
#' @param rows NULL for all input rows or a physical logical-position map reference.
#' @param columns Ordered data frame with output/source/cast character fields.
#' @param metadata Exact restoration metadata for the logical result.
#' @param request_id,identity Existing durable request and producer identity.
#' @return Small immutable component reference with its source/map dependencies.
#' @noRd
saveWorkflowTableSelection <- function(scope, rows, columns, metadata, request_id, identity) {
    checkWorkflowTableReadScope(scope)
    if (!is.data.frame(columns) || !identical(names(columns), c("output", "source", "cast")) ||
        !nrow(columns) || !all(vapply(columns, is.character, logical(1))) || anyNA(columns) ||
        any(!nzchar(unlist(columns, use.names = FALSE))) || anyDuplicated(columns$output) ||
        anyDuplicated(columns$source) || !all(columns$source %in% scope$metadata$columns) ||
        !all(columns$cast %in% c("identity", "factor_to_character")) ||
        !is.list(metadata) || !metadata$schema %in% c(2L, 3L) ||
        as.numeric(object.size(metadata)) > 8388608 ||
        !is.character(request_id) || length(request_id) != 1L ||
        !is.list(identity))
        cli::cli_abort("invalid logical selection", class = "saved_entity_invalid")
    if (any(vapply(seq_len(nrow(columns)), \(i)
        columns$cast[[i]] == "factor_to_character" &&
            is.null(scope$metadata$column_metadata[[columns$source[[i]]]]$levels), logical(1))))
        cli::cli_abort("factor cast requires a factor input", class = "saved_entity_unsupported")
    if (!is.null(scope$view) && isTRUE(scope$metadata$automatic_row_names) &&
        !isTRUE(metadata$automatic_row_names) &&
        (!identical(scope$view$rows$mode, "all") ||
            !isTRUE(scope$physical_metadata$automatic_row_names)))
        cli::cli_abort(paste("nested selection cannot retain these logical row names;",
            "export the parent with saveWorkflowTablePhysicalExport first"),
            class = "saved_entity_unsupported")
    if (is.null(rows) && identical(metadata, scope$metadata) &&
        identical(columns$output, columns$source) && all(columns$cast == "identity"))
        return(scope$reference)
    count <- if (is.null(rows)) scope$metadata$rows else {
        map <- newWorkflowTableReadScope(scope$store, rows, reader = "arrow", metadata_only = TRUE)
        on.exit(releaseWorkflowTableReadScope(map), add = TRUE)
        n <- map$metadata$rows
        releaseWorkflowTableReadScope(map)
        validateWorkflowTableRowMap(scope$store, rows, n, scope$metadata$rows)
        n
    }
    if (!identical(as.double(metadata$rows), as.double(count)) ||
        !identical(metadata$columns, columns$output))
        cli::cli_abort("logical selection metadata changes rows or columns",
            class = "saved_entity_invalid")
    mapping <- columns
    if (!is.null(scope$view)) {
        parent <- scope$view$columns
        at <- match(columns$source, parent$output)
        mapping$source <- parent$source[at]
        mapping$cast[parent$cast[at] == "factor_to_character"] <- "factor_to_character"
    }
    row_map <- if (is.null(rows)) {
        if (is.null(scope$view)) NULL else scope$view$rows$map
    } else if (is.null(scope$view) || is.null(scope$view$rows$map)) rows else {
        # Compose logical positions through two validated maps in one scoped SQL join.
        sql_scope <- newWorkflowTableReadScope(scope$store, scope$physical_reference,
            reader = "duckdb")
        on.exit(releaseWorkflowTableReadScope(sql_scope), add = TRUE)
        child <- newWorkflowTableReadScope(scope$store, rows, reader = "arrow", metadata_only = TRUE)
        on.exit(releaseWorkflowTableReadScope(child), add = TRUE)
        ancestor <- newWorkflowTableReadScope(scope$store, scope$view$rows$map,
            reader = "arrow", metadata_only = TRUE)
        on.exit(releaseWorkflowTableReadScope(ancestor), add = TRUE)
        con <- sql_scope$connection
        DBI::dbExecute(con, "SET threads = 8")
        DBI::dbExecute(con, "SET memory_limit = '12GiB'")
        DBI::dbExecute(con, "SET max_temp_directory_size = '0B'")
        quote_path <- \(path) as.character(DBI::dbQuoteString(con, path))
        query <- paste0('SELECT CAST(m."logical_ordinal" AS DOUBLE) AS "logical_ordinal", ',
            'CAST(p."source_ordinal" AS DOUBLE) AS "source_ordinal", ',
            'CAST(m."logical_ordinal" AS INTEGER) AS ".multischolar_order" ',
            "FROM read_parquet(", quote_path(child$path), ") m JOIN read_parquet(",
            quote_path(ancestor$path), ") p ON m.\"source_ordinal\" = p.\"logical_ordinal\"",
            ' ORDER BY m."logical_ordinal"')
        map_metadata <- workflowSavedTableMetadata(data.frame(
            logical_ordinal = double(), source_ordinal = double()), "logical_ordinal", 3L)
        map_metadata$rows <- count
        map_identity <- newSavedWorkflowScientificIdentity(scope$view_id,
            rows$entity_id, digest::digest(list(rows$entity_id,
                scope$view$rows$map$entity_id), algo = "sha256"),
            "composed-logical-ordinals-v2", identity$software)
        row_map <- saveWorkflowTableQuery(sql_scope, query, map_metadata,
            basename(tempfile("view-map-")), map_identity)
        releaseWorkflowTableReadScope(child)
        releaseWorkflowTableReadScope(ancestor)
        releaseWorkflowTableReadScope(sql_scope)
        validateWorkflowTableRowMap(scope$store, row_map, count,
            scope$physical_metadata$rows)
        row_map
    }
    if (!is.null(row_map)) validateWorkflowTableRowMap(scope$store, row_map,
        count, scope$physical_metadata$rows)
    row_spec <- list(mode = if (is.null(row_map)) "all" else "mapped",
        map = row_map, count = as.double(count))
    raw_metadata <- serialize(metadata, NULL, version = 3)
    view_id <- workflowTableViewId(scope$physical_reference, row_spec, mapping, raw_metadata)
    value <- structure(list(schema = "workflow.table_component.v2",
        source = scope$physical_reference, rows = row_spec, columns = mapping,
        metadata_raw = raw_metadata, view_id = view_id, indexes = list()),
        class = "WorkflowTableComponent")
    if (as.numeric(object.size(value)) > 8388608)
        cli::cli_abort("logical selection exceeds metadata limits", class = "saved_entity_limit")
    validateWorkflowTableView(scope$store, list(metadata = list(shape = c(count, nrow(mapping)))),
        value, metadata)
    identity$input <- scope$physical_reference$entity_id
    identity$design <- "physical-table"
    identity$parameters <- view_id
    identity$selectors <- "workflow-table-selection-v2"
    hit <- resolveSavedWorkflowEntity(scope$store, "result_table", identity, hydrate = FALSE)
    if (identical(hit$status, "hit")) return(hit$reference)
    saveWorkflowEntity(scope$store, value, request_id, identity,
        list(class = class(value), shape = as.integer(c(count, nrow(mapping))), factors = list()),
        kind = "result_table", codec = "rds")
}

#' Export a selected view as an ordinary physical table for older readers.
#' @param scope Selected view lease.
#' @param request_id,identity Durable request and producer identity.
#' @return Accepted physical Parquet descriptor.
#' @noRd
saveWorkflowTablePhysicalExport <- function(scope, request_id, identity) {
    checkWorkflowTableReadScope(scope)
    if (is.null(scope$view) && identical(scope$metadata, scope$physical_metadata))
        return(scope$physical_reference)
    if (!scope$metadata$schema %in% c(2L, 3L))
        cli::cli_abort("logical physical export requires version 2 or 3 restoration metadata",
            class = "saved_entity_unsupported")
    population <- if (is.null(scope$view)) digest::digest(list(
        scope$physical_reference$entity_id, scope$metadata), algo = "sha256") else scope$view_id
    identity$input <- population
    identity$parameters <- digest::digest(list(population, "physical-export-v2"), algo = "sha256")
    identity$selectors <- "workflow-table-physical-export-v2"
    metadata <- scope$metadata
    if (!metadata$key_column %in% c(".rownames", metadata$columns)) {
        metadata$key_column <- ".rownames"
        metadata$key_type <- "character"
        if (!isTRUE(metadata$automatic_row_names)) metadata$row_name_column <- ".multischolar_row_key"
    }
    saveWorkflowTable(scope$store, NULL, request_id, identity,
        metadata$key_column, metadata = metadata,
        write_fn = function(path, metadata) {
            writer <- openWorkflowTableWriter(path, metadata, scope$store$progress,
                backend = "fast_r")
            on.exit(abortWorkflowTableWriter(writer), add = TRUE)
            streamWorkflowTableRows(scope, metadata$columns, function(data, offset, row_names)
                appendWorkflowTableWriter(writer, data))
            finishWorkflowTableWriter(writer)
        })
}

#' Expand logical table dependencies without hydrating any assay.
#' @param store Open store.
#' @param references Identity-keyed saved descriptor list.
#' @return Descriptor closure including each canonical physical table and indexes.
#' @noRd
expandWorkflowTableReferences <- function(store, references) {
    for (reference in references) {
        component <- readWorkflowTableComponent(store, reference)
        if (is.null(component)) next
        children <- c(list(component$source),
            if (!is.null(component$view) && !is.null(component$view$rows$map))
                list(component$view$rows$map) else list(), component$indexes)
        for (child in children) {
            validateSavedWorkflowReference(store, child, payload = FALSE)
            previous <- references[[child$entity_id]]
            if (!is.null(previous) && !identical(previous, child))
                cli::cli_abort("conflicting table dependency", class = "saved_entity_corrupt")
            references[[child$entity_id]] <- child
        }
    }
    if (length(references) > store$max_entries || as.numeric(object.size(references)) > store$max_metadata_bytes)
        cli::cli_abort("table dependency closure exceeds metadata limits", class = "saved_entity_limit")
    references
}

#' Return the immutable population identity used for derived products.
#' @param scope Validated table lease.
#' @return Physical ID for v1; semantic logical view ID for v2.
#' @noRd
workflowTablePopulationId <- function(scope) {
    if (is.null(scope$view)) scope$physical_reference$entity_id else scope$view_id
}

#' Refuse a direct physical SQL/file access to a selected logical population.
#' @param scope Validated table scope.
#' @return NULL invisibly for physical and legacy metadata-only views.
#' @noRd
assertWorkflowPhysicalTableScope <- function(scope) {
    checkWorkflowTableReadScope(scope)
    if (!is.null(scope$view))
        cli::cli_abort("physical SQL consumer cannot read a selected table",
            class = "saved_entity_unsupported")
    invisible(NULL)
}

#' Estimate fixed per-row group storage before allocating projected buffers.
#' @param metadata Validated restoration metadata.
#' @param columns Required scientific columns.
#' @return Conservative fixed-width bytes; strings are also checked after each fetch.
#' @noRd
workflowProjectedRowBytes <- function(metadata, columns) {
    sizes <- vapply(metadata$column_metadata[columns], function(column)
        if (column$typeof %in% c("integer", "logical")) 4 else 8, numeric(1))
    sum(sizes) + 8 * sum(columns %in% names(metadata$name_columns)) +
        if (isTRUE(metadata$automatic_row_names)) 12 else 20
}

#' Persist a reusable group-to-source-ordinal index without sorting assay values.
#' @param scope Source lease; only the grouping column is scanned.
#' @param group Exact column defining complete independent groups.
#' @return Small index product and its immutable descriptor.
#' @noRd
saveWorkflowGroupIndex <- function(scope, group) {
    checkWorkflowTableReadScope(scope)
    checkSavedEntityString(group, "group column")
    if (!group %in% scope$metadata$columns) cli::cli_abort("group column is missing")
    population_id <- workflowTablePopulationId(scope)
    identity <- newSavedWorkflowScientificIdentity(population_id,
        "physical-row-population", digest::digest(list(group, scope$metadata$column_metadata[[group]]), algo = "sha256"),
        "stable-group-ordinals-v1", digest::digest(list(deparse(body(saveWorkflowGroupIndex),
            control = c("keepInteger", "keepNA"), width.cutoff = 500L),
            R.version[c("major", "minor", "platform")]), algo = "sha256"))
    hit <- resolveSavedWorkflowEntity(scope$store, "scientific_product", identity)
    if (identical(hit$status, "hit")) {
        value <- hit$value$value
        if (!identical(value$source_id, population_id) ||
            !identical(value$group, group)) cli::cli_abort("group index identity differs", class = "saved_entity_corrupt")
        return(list(reference = hit$reference, value = value))
    }
    rows <- scope$metadata$rows
    if (20 * as.double(rows) + 33554432 > 536870912)
        cli::cli_abort("group index workspace exceeds 512 MiB", class = "saved_entity_limit")
    codes <- integer(rows)
    values <- readSavedWorkflowTablePrototype(scope)[[group]]
    names(values) <- NULL
    streamSavedWorkflowTable(scope, group, function(data, offset, row_names) {
        incoming <- data[[group]]
        fresh <- unique(incoming[is.na(match(incoming, values))])
        values <<- c(values, fresh)
        names(values) <<- NULL
        if (as.numeric(object.size(values)) > 33554432)
            cli::cli_abort("group dictionary exceeds 32 MiB", class = "saved_entity_limit")
        if (nrow(data)) codes[seq.int(offset + 1L, length.out = nrow(data))] <<- match(incoming, values)
    })
    if (any(codes == 0L) || anyNA(codes)) cli::cli_abort("group index source coverage differs", class = "saved_entity_corrupt")
    positions <- order(codes, method = "radix")
    counts <- tabulate(codes, nbins = length(values))
    starts <- if (length(counts)) 1L + cumsum(c(0L, utils::head(counts, -1L))) else integer()
    directory <- data.frame(group_id = seq_along(values), value = values, size = counts,
        start = starts, first_order = positions[starts])
    group_identity <- identity
    group_identity$selectors <- "stable-group-directory-v1"
    groups <- saveWorkflowTable(scope$store, directory, basename(tempfile("table-groups-")),
        group_identity, "group_id", metadata_version = 3L)
    prototype <- data.frame(source_order = integer())
    metadata <- workflowSavedTableMetadata(prototype, "source_order", 3L)
    metadata$rows <- rows
    position_identity <- identity
    position_identity$selectors <- "stable-group-positions-v1"
    ordinals <- saveWorkflowTable(scope$store, NULL, basename(tempfile("table-positions-")),
        position_identity, "source_order", metadata = metadata, write_fn = function(path, restoration) {
            writer <- openWorkflowTableWriter(path, restoration, scope$store$progress)
            on.exit(abortWorkflowTableWriter(writer), add = TRUE)
            for (start in seq.int(1L, max(1L, rows), by = 10000L)) {
                at <- if (rows) seq.int(start, min(rows, start + 9999L)) else integer()
                appendWorkflowTableWriter(writer, data.frame(source_order = positions[at]))
            }
            finishWorkflowTableWriter(writer)
        })
    codes <- positions <- directory <- values <- NULL
    value <- list(source_id = population_id, group = group,
        groups = groups, positions = ordinals)
    reference <- saveWorkflowConsumerEntity(scope$store,
        newSavedWorkflowScientificProduct(value, "ordered-table-groups"),
        basename(tempfile("table-group-index-")), identity, "scientific_product")
    list(reference = reference, value = value)
}

#' Read a bounded index directory, never its assay or full ordinal payload.
#' @param store Open store.
#' @param reference Group directory descriptor.
#' @return Small typed directory.
#' @noRd
readWorkflowIndexDirectory <- function(store, reference) {
    scope <- newWorkflowTableReadScope(store, reference, reader = "arrow")
    on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
    if (scope$metadata$rows * 32 > 33554432)
        cli::cli_abort("group directory exceeds 32 MiB", class = "saved_entity_limit")
    parts <- list()
    bytes <- 0
    streamSavedWorkflowTable(scope, scope$metadata$columns, function(data, offset, row_names) {
        bytes <<- bytes + as.numeric(object.size(data))
        if (bytes > 33554432) cli::cli_abort("group directory exceeds 32 MiB", class = "saved_entity_limit")
        parts[[length(parts) + 1L]] <<- data
    })
    do.call(rbind, parts)
}

#' Gather complete groups through pushed-down ordinal selections.
#' @param scope Source lease.
#' @param columns Scientific projection.
#' @param positions Saved group-to-source-ordinal table.
#' @param groups Admitted rows from the small group directory, in execution order.
#' @return Group-contiguous data with stable source identities.
#' @noRd
readWorkflowGroupCohort <- function(scope, columns, positions, groups) {
    rows <- sum(groups$size)
    if (as.double(rows) * workflowProjectedRowBytes(scope$metadata, columns) > 268435456)
        cli::cli_abort("group cohort exceeds 256 MiB", class = "saved_entity_limit")
    selected_rows <- integer(rows)
    index_scope <- newWorkflowTableReadScope(scope$store, positions, reader = "arrow")
    on.exit(releaseWorkflowTableReadScope(index_scope), add = TRUE)
    target <- 0L
    starts <- sizes <- integer(nrow(groups))
    ranges <- 0L
    for (group in seq_len(nrow(groups))) {
        if (ranges && groups$start[[group]] == starts[[ranges]] + sizes[[ranges]]) {
            sizes[[ranges]] <- sizes[[ranges]] + groups$size[[group]]
        } else {
            ranges <- ranges + 1L
            starts[[ranges]] <- groups$start[[group]]
            sizes[[ranges]] <- groups$size[[group]]
        }
    }
    # Scan selected index ranges once in physical order, mapping them back to
    # requested group order. Never reopen/decompress a row group per fetch batch.
    index_rows <- unlist(lapply(seq_len(ranges), function(range)
        seq.int(starts[[range]], length.out = sizes[[range]])), use.names = FALSE)
    index_order <- order(index_rows, method = "radix")
    streamSavedWorkflowTable(index_scope, "source_order", function(data, offset, row_names) {
        page <- data$source_order
        if (anyNA(page) || any(page < 1L | page > scope$metadata$rows) || anyDuplicated(page))
            cli::cli_abort("group ordinal index is invalid", class = "saved_entity_corrupt")
        selected_rows[index_order[offset + seq_along(page)]] <<- page
        target <<- target + length(page)
    }, source_rows = index_rows[index_order])
    index_rows <- index_order <- NULL
    releaseWorkflowTableReadScope(index_scope)
    if (target != rows) cli::cli_abort("group index is incomplete", class = "saved_entity_corrupt")
    if (anyDuplicated(selected_rows))
        cli::cli_abort("group ordinal index repeats source rows", class = "saved_entity_corrupt")
    destinations <- order(selected_rows, method = "radix")
    physical_rows <- selected_rows[destinations]
    selected_rows <- NULL
    buffer <- newWorkflowTableBuffer(scope$metadata, columns, rows)
    unchecked_bytes <- 0
    streamSavedWorkflowTable(scope, columns, function(data, offset, row_names) {
        source <- attr(data, "dia_source_order")
        keep <- seq_len(nrow(data))
        at <- destinations[offset + keep]
        assignWorkflowTableBuffer(buffer, data, at, source[keep], row_names[keep])
        unchecked_bytes <<- unchecked_bytes + as.numeric(object.size(data[keep, , drop = FALSE]))
        if (unchecked_bytes >= 16777216) {
            if (as.numeric(object.size(buffer$columns)) > 268435456)
                cli::cli_abort("group cohort exceeds 256 MiB", class = "saved_entity_limit")
            unchecked_bytes <<- 0
        }
    }, source_rows = physical_rows)
    destinations <- physical_rows <- NULL
    if (as.numeric(object.size(buffer$columns)) > 268435456)
        cli::cli_abort("group cohort exceeds 256 MiB", class = "saved_entity_limit")
    finishWorkflowTableBuffer(buffer)
}

#' Persist a narrow stable ordering index within an explicit memory allowance.
#' @param scope Existing table lease.
#' @param columns Required sort keys; the physical ordinal is the stable tie-breaker.
#' @return Scientific product containing an optional permutation-table reference.
#' @noRd
saveWorkflowOrderIndex <- function(scope, columns) {
    checkWorkflowTableReadScope(scope)
    columns <- setdiff(columns, ".multischolar_order")
    if (!all(columns %in% scope$metadata$columns) || anyDuplicated(columns))
        cli::cli_abort("ordering key is absent", class = "saved_entity_invalid")
    population_id <- workflowTablePopulationId(scope)
    identity <- newSavedWorkflowScientificIdentity(population_id,
        "physical-row-order", digest::digest(list(columns, scope$metadata$column_metadata[columns]), algo = "sha256"),
        "stable-table-order-v1", digest::digest(list(deparse(body(saveWorkflowOrderIndex),
            control = c("keepInteger", "keepNA"), width.cutoff = 500L),
            R.version[c("major", "minor", "platform")]), algo = "sha256"))
    hit <- resolveSavedWorkflowEntity(scope$store, "scientific_product", identity)
    if (identical(hit$status, "hit")) {
        value <- hit$value$value
        if (!identical(value$source_id, population_id) ||
            !identical(value$columns, columns) ||
            !is.null(value$positions) && !isWorkflowTableReference(value$positions))
            cli::cli_abort("ordering index identity differs", class = "saved_entity_corrupt")
        return(list(reference = hit$reference, value = value))
    }
    if (isTRUE(scope$metadata_only) && length(columns)) {
        scope <- newWorkflowTableReadScope(scope$store, scope$reference, reader = "arrow")
        on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
    }
    rows <- scope$metadata$rows
    fixed <- sum(vapply(scope$metadata$column_metadata[columns], function(spec)
        if (spec$typeof %in% c("integer", "logical")) 4 else 8, numeric(1)))
    reserve <- as.double(rows) * 12 + 33554432
    if (as.double(rows) * fixed + reserve > 536870912)
        cli::cli_abort("ordering workspace exceeds 512 MiB", class = "saved_entity_limit")
    keys <- lapply(scope$metadata$column_metadata[columns], function(spec)
        if (!is.null(spec$levels)) factor(rep(NA_character_, rows), levels = spec$levels,
            ordered = isTRUE(spec$ordered)) else vector(spec$typeof, rows))
    unchecked <- 0
    if (length(columns)) streamSavedWorkflowTable(scope, columns, function(data, offset, row_names) {
        at <- offset + seq_len(nrow(data))
        for (column in columns) keys[[column]][at] <<- data[[column]]
        unchecked <<- unchecked + as.numeric(object.size(data))
        if (unchecked >= 16777216) {
            if (as.numeric(object.size(keys)) + reserve > 536870912)
                cli::cli_abort("ordering workspace exceeds 512 MiB", class = "saved_entity_limit")
            unchecked <<- 0
        }
    })
    if (as.numeric(object.size(keys)) + reserve > 536870912)
        cli::cli_abort("ordering workspace exceeds 512 MiB", class = "saved_entity_limit")
    positions <- if (length(keys)) do.call(order, c(unname(keys), list(method = "radix", na.last = TRUE))) else seq_len(rows)
    keys <- NULL
    unchanged <- TRUE
    for (start in seq.int(1L, max(1L, rows), by = 10000L)) {
        at <- if (rows) seq.int(start, min(rows, start + 9999L)) else integer()
        if (!identical(as.integer(positions[at]), at)) { unchanged <- FALSE; break }
    }
    ordinal <- NULL
    if (!unchanged) {
        metadata <- workflowSavedTableMetadata(data.frame(source_order = integer()), "source_order", 3L)
        metadata$rows <- rows
        table_identity <- identity
        table_identity$selectors <- "stable-table-order-positions-v1"
        ordinal <- saveWorkflowTable(scope$store, NULL, basename(tempfile("table-order-")),
            table_identity, "source_order", metadata = metadata, write_fn = function(path, restoration) {
                writer <- openWorkflowTableWriter(path, restoration, scope$store$progress)
                on.exit(abortWorkflowTableWriter(writer), add = TRUE)
                for (start in seq.int(1L, max(1L, rows), by = 10000L)) {
                    at <- if (rows) seq.int(start, min(rows, start + 9999L)) else integer()
                    appendWorkflowTableWriter(writer, data.frame(source_order = positions[at]))
                }
                finishWorkflowTableWriter(writer)
            })
    }
    positions <- NULL
    value <- list(source_id = population_id, columns = columns, positions = ordinal)
    reference <- saveWorkflowConsumerEntity(scope$store,
        newSavedWorkflowScientificProduct(value, "ordered-table-rows"), basename(tempfile("table-order-product-")),
        identity, "scientific_product")
    list(reference = reference, value = value)
}

#' Read an ordered table through bounded cohorts, never a whole-result SQL sort.
#' @param store Existing owner.
#' @param reference Saved input, including a committed change/ledger table.
#' @param columns Projected fields.
#' @param order_columns Native ordering keys.
#' @return Iterator with next_batch(), dependencies and idempotent close().
#' @noRd
openWorkflowOrderedRows <- function(store, reference, columns, order_columns) {
    scope <- newWorkflowTableReadScope(store, reference, reader = "arrow")
    ready <- FALSE
    on.exit(if (!ready) releaseWorkflowTableReadScope(scope), add = TRUE)
    index <- saveWorkflowOrderIndex(scope, order_columns)
    iterator <- new.env(parent = emptyenv())
    iterator$dependencies <- Filter(Negate(is.null), list(reference, index$reference, index$value$positions))
    iterator$cohort <- NULL
    iterator$cohort_at <- 1L
    iterator$offset <- 0L
    iterator$closed <- FALSE
    direct <- if (is.null(index$value$positions) && is.null(scope$view))
        openWorkflowArrowCursor(scope, columns) else NULL
    iterator$close <- function() {
        if (!iterator$closed) {
            iterator$cohort <- NULL
            if (!is.null(direct)) direct$close()
            releaseWorkflowTableReadScope(scope)
            iterator$closed <- TRUE
        }
        invisible(NULL)
    }
    iterator$next_batch <- function() {
        if (iterator$closed) return(NULL)
        if (is.null(index$value$positions)) {
            if (iterator$offset >= scope$metadata$rows) { iterator$close(); return(NULL) }
            if (!is.null(scope$view)) {
                page <- readSavedWorkflowTable(store, reference, columns,
                    offset = iterator$offset,
                    limit = min(10000L, scope$metadata$rows - iterator$offset), scope = scope)
                iterator$offset <- iterator$offset + nrow(page$data)
                return(page$data)
            }
            part <- direct$fetch()
            if (is.null(part) || !identical(as.integer(part$.source_order),
                seq.int(iterator$offset + 1L, length.out = nrow(part))))
                cli::cli_abort("ordered cursor source is incomplete", class = "saved_entity_corrupt")
            iterator$offset <- iterator$offset + nrow(part)
            return(restoreWorkflowCursorBatch(scope, columns, part)$data)
        }
        if (is.null(iterator$cohort) || iterator$cohort_at > nrow(iterator$cohort)) {
            iterator$cohort <- NULL
            if (iterator$offset >= scope$metadata$rows) { iterator$close(); return(NULL) }
            count <- as.integer(min(scope$metadata$rows - iterator$offset,
                floor(268435456 / workflowProjectedRowBytes(scope$metadata, columns))))
            if (count < 1L) cli::cli_abort("ordered row exceeds workspace", class = "saved_entity_limit")
            iterator$cohort <- readWorkflowGroupCohort(scope, columns, index$value$positions,
                data.frame(size = count, start = iterator$offset + 1L))
            iterator$offset <- iterator$offset + count
            iterator$cohort_at <- 1L
        }
        at <- seq.int(iterator$cohort_at, min(nrow(iterator$cohort), iterator$cohort_at + 9999L))
        part <- iterator$cohort[at, , drop = FALSE]
        attr(part, "workflow_source_order") <- NULL
        iterator$cohort_at <- utils::tail(at, 1L) + 1L
        part
    }
    ready <- TRUE
    iterator
}

#' Allocate an owned staging directory for final batch products, separate from engine spill.
#' @param scope Active action lease.
#' @return Directory removed after its writers close or commit.
#' @noRd
workflowTableBatchDirectory <- function(scope) {
    checkWorkflowTableReadScope(scope)
    if (is.null(scope$action$batch_directory)) {
        path <- tempfile("table-batches-", tmpdir = file.path(scope$store$directory, ".pending"))
        if (!dir.create(path)) cli::cli_abort("cannot create table batch staging", class = "saved_entity_write")
        scope$action$batch_directory <- path
    }
    scope$action$batch_directory
}

#' Read a bounded ordered selection of physical row ordinals without a global sort.
#' @param store,reference Existing table owner and immutable descriptor.
#' @param columns At most 64 projected fields.
#' @param ordinals At most 10,000 one-based source ordinals, preserving duplicates.
#' @param scope Optional existing physical DuckDB lease.
#' @return Typed data in exactly requested order.
#' @noRd
readWorkflowSelectedOrdinals <- function(store, reference, columns, ordinals, scope = NULL) {
    if (is.null(scope)) {
        scope <- newWorkflowTableReadScope(store, reference, reader = "duckdb")
        on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
    }
    if (!identical(scope$store, store) || !identical(scope$reference, reference) ||
        !identical(scope$reader, "duckdb"))
        cli::cli_abort("ordinal query scope differs", class = "saved_entity_scope")
    assertWorkflowPhysicalTableScope(scope)
    if (!is.numeric(ordinals) || length(ordinals) > 10000L || anyNA(ordinals) ||
        any(!is.finite(ordinals) | ordinals != floor(ordinals) | ordinals < 1 | ordinals > scope$metadata$rows) ||
        !is.character(columns) || !length(columns) || length(columns) > 64L ||
        !all(columns %in% scope$metadata$columns))
        cli::cli_abort("invalid ordinal projection", class = "saved_entity_limit")
    if (!length(ordinals)) return(readSavedWorkflowTablePrototype(scope)[columns])
    physical <- workflowTableStorageColumns(scope$metadata, columns)
    q <- function(value) as.character(DBI::dbQuoteIdentifier(scope$connection, value))
    selected <- unique(ordinals)
    query <- paste0('SELECT ".multischolar_order" AS ".source_order", ',
        paste(q(physical), collapse = ", "), ' FROM read_parquet(?) WHERE ".multischolar_order" IN (',
        paste(rep("?", length(selected)), collapse = ","), ")")
    data <- DBI::dbGetQuery(scope$connection, query, params = c(list(scope$path), as.list(selected)))
    at <- match(ordinals, data$.source_order)
    if (anyNA(at) || anyDuplicated(data$.source_order))
        cli::cli_abort("ordinal projection differs from source", class = "saved_entity_corrupt")
    data <- data[at, , drop = FALSE]
    source_names <- workflowTableSourceRowNames(data, scope$metadata, data$.source_order)
    data <- restoreSavedWorkflowTable(data, scope$metadata, columns)
    if (!anyDuplicated(source_names) && !anyNA(source_names))
        data <- restoreWorkflowTableRowNames(data, scope$metadata, source_names)
    if (as.numeric(object.size(data)) > 67108864)
        cli::cli_abort("ordinal projection exceeds 64 MiB", class = "saved_entity_limit")
    data
}

#' Execute SQL on the owned native scope or the existing DBI scope.
#' @param con,statement Existing owner and quoted SQL.
#' @param params Optional DBI parameters; native rollup uses quoted fixed statements.
#' @return Existing backend result, invisibly for native execution.
#' @noRd
workflowQueryExecute <- function(con, statement, params = list()) {
    if (!inherits(con, "workflow_fast_query")) return(if (length(params))
        DBI::dbExecute(con, statement, params = params) else DBI::dbExecute(con, statement))
    fastduckdb::execute_query_fast(con, statement, params = params)
}

#' Collect only bounded metadata from the owned native query scope.
#' @param con,statement,params Existing owner, quoted SQL and optional DBI parameters.
#' @return Bounded data frame, with native scientific types restored by the caller's codec.
#' @noRd
workflowQueryCollect <- function(con, statement, params = list()) {
    if (!inherits(con, "workflow_fast_query")) return(if (length(params))
        DBI::dbGetQuery(con, statement, params = params) else DBI::dbGetQuery(con, statement))
    fastduckdb::query_data_fast(con, statement, params = params)
}

#' Quote identifiers through the existing SQL quoting contract.
#' @param con,x Owned connection and identifiers.
#' @return Quoted SQL identifiers.
#' @noRd
workflowQueryQuoteIdentifier <- function(con, x) {
    DBI::dbQuoteIdentifier(if (inherits(con, "workflow_fast_query")) DBI::ANSI() else con, x)
}

#' Quote scalar string values through the existing SQL quoting contract.
#' @param con,x Owned connection and string values.
#' @return Quoted SQL values.
#' @noRd
workflowQueryQuoteString <- function(con, x) {
    DBI::dbQuoteString(if (inherits(con, "workflow_fast_query")) DBI::ANSI() else con, x)
}

#' Release the query owner after its last saved-table scope closes.
#' @param con Owned connection.
#' @return NULL invisibly after all native results close, or the DBI disconnect result.
#' @noRd
workflowQueryDisconnect <- function(con) {
    if (inherits(con, "workflow_fast_query")) return(fastduckdb::close_query_fast(con))
    DBI::dbDisconnect(con, shutdown = TRUE)
}

#' Register only a bounded annotation or sample-rank lookup.
#' @param con,name,value Existing owner, table name and small primitive lookup.
#' @param temporary Must be TRUE for native query registration.
#' @return TRUE invisibly, after bounded inserts finish.
#' @noRd
workflowQueryWriteLookup <- function(con, name, value, temporary = TRUE) {
    if (!inherits(con, "workflow_fast_query")) return(DBI::dbWriteTable(con, name, value,
        temporary = temporary))
    if (!isTRUE(temporary) || !is.data.frame(value) || nrow(value) > 10000L ||
        !length(value) || anyNA(names(value)) || any(!nzchar(names(value))) || anyDuplicated(names(value)))
        cli::cli_abort("invalid native lookup table", class = "saved_entity_limit")
    value <- lapply(value, \(column) if (is.factor(column)) as.character(column) else column)
    types <- vapply(value, typeof, character(1))
    if (any(!types %in% c("double", "integer", "logical", "character")) ||
        any(vapply(value, \(column) !is.null(attributes(column)), logical(1))))
        cli::cli_abort("unsupported native lookup column", class = "saved_entity_unsupported")
    # Admit SQL-escaped text as well as the existing R lookup before constructing statements.
    bytes <- sum(vapply(value, \(column) 32 * (length(column) + 1) +
        if (is.character(column)) 2 * sum(nchar(enc2utf8(column), type = "bytes",
            allowNA = TRUE), na.rm = TRUE) else 0, numeric(1)))
    if (bytes > 67108864) cli::cli_abort("native lookup exceeds 64 MiB", class = "saved_entity_limit")
    quote <- \(x) as.character(workflowQueryQuoteIdentifier(con, x))
    sql_types <- c(double = "DOUBLE", integer = "INTEGER", logical = "BOOLEAN", character = "VARCHAR")
    workflowQueryExecute(con, paste0("CREATE TEMP TABLE ", quote(name), " (",
        paste(paste(quote(names(value)), sql_types[types]), collapse = ","), ")"))
    table <- as.data.frame(value, stringsAsFactors = FALSE, check.names = FALSE)
    fastduckdb::append_query_table_fast(con, name, table)
    invisible(TRUE)
}
