#' The analysis record (KIT-027): a versioned JSON description of one DIA analysis, derived on demand from the
#' state head and each state's metadata record. Nothing is written incrementally and no assay is hydrated; the
#' store is the source and the record is a view of it.
#' @noRd
NULL

#' Encode one argument value for the record without truncation: atomic vectors in full with their type (and
#' names, levels or class when present) so a JSON round trip restores them `identical()`; data frames as a
#' reference-sized summary; lists recursed.
#' @noRd
encodeWorkflowRecordValue <- function(value) {
    if (is.null(value)) return(NULL)
    if (is.data.frame(value)) {
        return(list(type = "data.frame", rows = nrow(value), columns = names(value),
            sha256 = digest::digest(value, algo = "sha256")))
    }
    if (is.atomic(value)) {
        out <- list(type = typeof(value), value = if (is.factor(value)) as.integer(value) else as.vector(unclass(value)))
        if (!is.null(names(value))) out$names <- names(value)
        if (is.factor(value)) out$levels <- levels(value)
        extra <- setdiff(class(value), c(class(as.vector(unclass(value))), "factor"))
        if (!is.factor(value) && length(extra) && !identical(class(value), class(as.vector(unclass(value)))))
            out$class <- class(value)
        out$length <- length(value)
        return(out)
    }
    if (is.list(value)) {
        out <- list(type = "list", value = lapply(value, encodeWorkflowRecordValue), length = length(value))
        if (!is.null(names(value))) out$names <- names(value)
        return(out)
    }
    list(type = "opaque", class = class(value))
}

#' Restore an encoded atomic or list value from a parsed record (the inverse of `encodeWorkflowRecordValue`
#' for atomic vectors and lists; data frames stay summaries).
#' @noRd
decodeWorkflowRecordValue <- function(encoded) {
    if (is.null(encoded)) return(NULL)
    type <- encoded$type
    if (identical(type, "list")) {
        value <- lapply(encoded$value, decodeWorkflowRecordValue)
        if (!length(value)) value <- list()
        if (!is.null(encoded$names)) names(value) <- unlist(encoded$names) else names(value) <- NULL
        return(value)
    }
    if (!type %in% c("logical", "integer", "double", "character", "complex")) return(encoded)
    raw <- encoded$value
    raw <- if (is.list(raw)) lapply(raw, function(x) if (is.null(x)) NA else x) else raw
    value <- unlist(raw, use.names = FALSE)
    if (is.null(value)) value <- vector(type, 0L)
    value <- switch(type, logical = as.logical(value), integer = as.integer(value), double = as.double(value),
        character = as.character(value), complex = as.complex(value))
    if (length(value) != (encoded$length %||% length(value))) value <- rep_len(value, encoded$length)
    if (!is.null(encoded$levels)) value <- factor(unlist(encoded$levels)[value], levels = unlist(encoded$levels))
    if (!is.null(encoded$names)) names(value) <- unlist(encoded$names)
    if (!is.null(encoded$class)) class(value) <- unlist(encoded$class)
    value
}

#' Compare a recorded parameter value with the shipped `config.ini` string for its mapped key.
#' @noRd
workflowParameterMatchesDefault <- function(value, shipped) {
    if (is.null(shipped)) shipped <- ""
    shipped <- trimws(as.character(shipped))
    if (is.null(value) || !length(value) || all(is.na(value))) return(!nzchar(paste(shipped, collapse = "")))
    shipped <- paste(shipped, collapse = ",")
    # KIT-027 audit PH-07: a shipped list is split on commas only for a multi-valued recorded value, so a single
    # value containing a comma (a regex, a formula) is compared whole; numbers compare exactly.
    parts <- if (length(value) > 1L) trimws(strsplit(shipped, ",", fixed = TRUE)[[1L]]) else shipped
    if (length(parts) != length(value)) return(FALSE)
    if (is.logical(value)) return(identical(as.logical(value), as.logical(suppressWarnings(as.numeric(parts)))) ||
        identical(as.logical(value), as.logical(toupper(parts))))
    if (is.numeric(value)) {
        numbers <- suppressWarnings(as.numeric(parts))
        return(!anyNA(numbers) && all(as.numeric(value) == numbers))
    }
    identical(as.character(value), parts)
}

#' One step's declared parameters: each value from the kernel's `slots$args` section, else the values the
#' operation received (head provenance, KIT-027 audit PH-01), else the state's `config`, with its source named.
#' @noRd
workflowRecordDeclaredParameters <- function(names_declared, entry, section, args, config, values, shipped) {
    lapply(names_declared, function(name) {
        in_args <- !is.null(section) && is.list(args[[section]]) && name %in% names(args[[section]])
        in_values <- is.list(values) && name %in% names(values)
        in_config <- is.list(config) && name %in% names(config)
        value <- if (in_args) args[[section]][[name]] else if (in_values) values[[name]] else if (in_config) config[[name]] else NULL
        source <- if (in_args) "args" else if (in_values) "provenance" else if (in_config) "config" else "missing"
        list(name = name, source = source, value = encodeWorkflowRecordValue(value),
            provenance = if (identical(source, "missing")) "unknown" else workflowParameterProvenance(entry$defaults[[name]], value, shipped))
    })
}

#' Provenance of one declared parameter through its registry `defaults` mapping, never by name.
#' @noRd
workflowParameterProvenance <- function(mapping, value, shipped_config) {
    if (is.null(mapping)) return("unknown")
    if (is.character(mapping)) return(switch(mapping, fixed = "fixed", input = "input", no_default = "no_default", "unknown"))
    shipped <- shipped_config[[mapping$section]][[mapping$key]]
    if (is.null(shipped)) return("unknown")
    if (workflowParameterMatchesDefault(value, shipped)) "default" else "user_set"
}

#' Read the shipped `config.ini` as strings, section by key.
#' @noRd
readWorkflowShippedDefaults <- function(path = system.file("config", "config.ini", package = "MultiScholaR")) {
    if (!nzchar(path) || !file.exists(path)) return(list())
    lines <- readLines(path, warn = FALSE)
    out <- list()
    section <- NULL
    for (line in lines) {
        line <- trimws(line)
        if (!nzchar(line) || startsWith(line, ";") || startsWith(line, "#")) next
        if (grepl("^\\[.*\\]$", line)) { section <- sub("^\\[(.*)\\]$", "\\1", line); next }
        if (is.null(section) || !grepl("=", line, fixed = TRUE)) next
        key <- trimws(sub("=.*$", "", line))
        out[[section]][[key]] <- trimws(sub("^[^=]*=", "", line))
    }
    out
}

#' Build the analysis record from head metadata (schema `multischolar.analysis_record.v1`).
#' @param manager On-disk DIA owner (`DIAStateAdapter`).
#' @param now Generation time.
#' @return The record as a list; no state assay is hydrated (`full_s4_reads` unchanged).
#' @noRd
buildWorkflowAnalysisRecord <- function(manager, now = Sys.time()) {
    if (!inherits(manager, "DIAStateAdapter")) cli::cli_abort("the analysis record needs an on-disk DIA owner", class = "saved_entity_scope")
    context <- manager$.context
    snapshot <- context$snapshot
    history <- unlist(snapshot$state_history)
    shipped <- readWorkflowShippedDefaults()
    defaults_path <- system.file("config", "config.ini", package = "MultiScholaR")
    installed_defaults <- if (nzchar(defaults_path)) digest::digest(file = defaults_path, algo = "sha256") else NA_character_
    sections <- stats::setNames(names(workflowSectionOperations), unname(workflowSectionOperations))
    records <- lapply(stats::setNames(history, history), function(state) {
        reference <- snapshot$states[[state]]
        if (is.null(reference)) return(NULL)
        tryCatch(readDIAStateMetadata(manager, reference), saved_entity_unsupported = function(e) {
            # The `initial` placeholder is a small record with no data; a state stored as a whole legacy object is
            # not hydrated here and is listed without its arguments.
            if (identical(state, "initial")) resolveDIAStateValue(context, reference) else list(legacy_object = TRUE)
        })
    })
    steps <- list()
    for (i in seq_along(history)) {
        state <- history[[i]]
        record <- records[[state]]
        provenance <- snapshot$step_provenance[[state]] %||% list(operation = "legacy", variant = NA_character_, parameter_names = character())
        operation <- provenance$operation
        entry <- tryCatch(workflowOperationEntry(operation), error = function(e) NULL)
        args <- record$data$slots$args
        config <- record$config
        section <- unname(sections[operation])
        if (!length(section) || is.na(section)) section <- NULL
        names_declared <- as.character(provenance$parameter_names %||% character())
        parameters <- list()
        if (length(names_declared)) {
            parameters <- workflowRecordDeclaredParameters(names_declared, entry, section, args, config, provenance$values, shipped)
            for (name in names(provenance$fixed)) parameters[[length(parameters) + 1L]] <- list(name = name, source = "fixed",
                value = encodeWorkflowRecordValue(provenance$fixed[[name]]), provenance = "fixed")
        } else if (!is.null(args)) {
            # A legacy step names no parameters: its recorded arguments are listed whole.
            for (name in names(args)) parameters[[length(parameters) + 1L]] <- list(name = name, source = "args",
                value = encodeWorkflowRecordValue(args[[name]]), provenance = "unknown")
        }
        previous <- record$previous_state
        input <- if (is.character(previous) && length(previous) == 1L && !is.na(previous)) snapshot$states[[previous]]$entity_id else NULL
        products <- lapply(unname(record$data$products %||% list()), function(product) {
            if (!isSavedWorkflowReference(product)) return(NULL)
            list(entity_id = product$entity_id, kind = product$kind,
                product_type = product$metadata$product_type %||% product$identity$type %||% NA_character_,
                identity_sha256 = digest::digest(product$identity, algo = "sha256"))
        })
        software <- snapshot$software_environments[[provenance$build %||% ""]]
        sizes <- snapshot$sizes[[state]]
        # KIT-027 audit PH-04: a `legacy` entry written by the reopen backfill carries the reopening build, defaults
        # and time; such a step is dated by its state record and its software is unknown.
        committed_at <- provenance$committed_at %||% NA_character_
        if (identical(operation, "legacy")) {
            stamped <- if (inherits(record$timestamp, "POSIXt")) format(record$timestamp, "%Y-%m-%dT%H:%M:%OS3%z") else NA_character_
            entry_time <- suppressWarnings(as.POSIXct(committed_at, format = "%Y-%m-%dT%H:%M:%OS%z"))
            backfilled <- isTRUE(provenance$backfilled) || inherits(record$timestamp, "POSIXt") && !is.na(entry_time) &&
                abs(as.numeric(difftime(entry_time, record$timestamp, units = "secs"))) > 60
            if (!is.na(stamped)) committed_at <- stamped
            if (backfilled) { software <- NULL; provenance$source_commit <- NA_character_; provenance$defaults_sha256 <- NA_character_ }
        }
        preceded_by <- lapply(provenance$preceded_by %||% list(), function(previous) {
            previous_entry <- tryCatch(workflowOperationEntry(previous$operation), error = function(e) NULL)
            previous_section <- unname(sections[previous$operation])
            if (!length(previous_section) || is.na(previous_section)) previous_section <- NULL
            list(operation = previous$operation, variant = previous$variant %||% NA_character_,
                parameters = workflowRecordDeclaredParameters(as.character(previous$parameter_names %||% character()),
                    previous_entry, previous_section, args, NULL, previous$values, shipped))
        })
        steps[[i]] <- list(order = i, state = state, operation = operation, variant = provenance$variant %||% NA_character_,
            preceded_by = preceded_by,
            step_label = entry$step_label %||% state, committed_at = committed_at,
            seed = provenance$seed %||% NA_real_, result_entity_id = snapshot$states[[state]]$entity_id,
            input_entity_ids = if (is.null(input)) list() else list(input),
            software = list(version = software$version %||% NA_character_, sha = provenance$source_commit %||% NA_character_),
            defaults_sha256 = provenance$defaults_sha256 %||% NA_character_,
            dimensions = list(rows = sizes$rows %||% NA_integer_, columns = sizes$columns %||% NA_integer_),
            description = record$description %||% NA_character_,
            parameters = parameters, products = Filter(Negate(is.null), products),
            stored_as_object = isTRUE(record$legacy_object))
    }
    last <- if (length(history)) records[[history[[length(history)]]]] else NULL
    design <- NULL
    if (!is.null(last$data$components$design_matrix)) {
        rows <- materializeDIATable(context, last$data$components$design_matrix)
        slots <- last$data$slots
        design <- list(entity_id = last$data$components$design_matrix$entity_id, sample_id = slots$sample_id %||% NA_character_,
            group_id = slots$group_id %||% NA_character_, technical_replicate_id = slots$technical_replicate_id %||% NA_character_,
            rows = rows)
    }
    # KIT-027 audit PH-02: the DA step commits no state of its own; its settings are the accepted result manifest's.
    analyses <- list()
    if (isSavedWorkflowReference(snapshot$da_results)) {
        manifest <- resolveDIAStateValue(context, snapshot$da_results)
        da_entry <- tryCatch(workflowOperationEntry("DA-contrast-output"), error = function(e) NULL)
        source_state <- names(snapshot$states)[vapply(snapshot$states, function(reference)
            identical(reference$entity_id, manifest$source$entity_id), logical(1))]
        analyses$differential_abundance <- list(operation = "DA-contrast-output",
            source_state = if (length(source_state)) source_state[[1L]] else NA_character_,
            source_entity_id = manifest$source$entity_id %||% NA_character_,
            parameters = lapply(names(manifest$parameters), function(name) list(name = name, source = "manifest",
                value = encodeWorkflowRecordValue(manifest$parameters[[name]]),
                provenance = workflowParameterProvenance(da_entry$defaults[[name]], manifest$parameters[[name]], shipped))),
            contrasts = lapply(seq_len(nrow(manifest$contrasts %||% data.frame())), function(j) as.list(manifest$contrasts[j, , drop = FALSE])),
            model = manifest$model,
            results = lapply(names(manifest$results), function(name) list(contrast = name,
                entity_id = manifest$results[[name]]$tables$da_proteins_long$entity_id %||% NA_character_,
                rows = manifest$results[[name]]$rows %||% NA_integer_)))
    }
    inputs <- list()
    for (name in names(snapshot$file_dependencies)) {
        dependency <- snapshot$file_dependencies[[name]]
        # The commit also keys each execution input by its state's entity id; the role entry is listed once.
        if (grepl("_[0-9a-f]{16}$", name)) next
        inputs[[length(inputs) + 1L]] <- list(role = dependency$role %||% NA_character_,
            basename = dependency$name %||% if (!is.null(dependency$source)) basename(dependency$source) else NA_character_,
            bytes = dependency$bytes %||% NA_real_, sha256 = dependency$sha256 %||% NA_character_,
            entity_id = dependency$entity_id %||% NA_character_)
    }
    builds <- unique(vapply(snapshot$step_provenance[history], function(p) p$build %||% NA_character_, character(1)))
    builds <- builds[!is.na(builds)]
    current <- snapshot$software_environments[[if (length(builds)) builds[[length(builds)]] else ""]]
    step_defaults <- unique(vapply(steps, function(step) step$defaults_sha256, character(1)))
    list(schema = "multischolar.analysis_record.v1", generated_at = format(now, "%Y-%m-%dT%H:%M:%OS3%z"),
        workflow = snapshot$workflow_type %||% NA_character_,
        store_identity = list(project_id = context$store$identity$project_id %||% NA_character_,
            workflow_id = context$store$identity$workflow_id %||% NA_character_),
        software = list(package = "MultiScholaR", version = current$version %||% NA_character_,
            sha = current$source_commit %||% NA_character_, r_version = current$r_version %||% NA_character_,
            platform = current$platform %||% NA_character_,
            packages = lapply(names(current$packages %||% list()), function(name)
                list(name = name, version = as.character(current$packages[[name]])))),
        defaults_sha256 = installed_defaults,
        defaults_changed = any(!is.na(step_defaults) & step_defaults != installed_defaults),
        inputs = inputs, design = design, steps = steps, analyses = analyses,
        external_resources = lapply(snapshot$external_resources %||% list(), function(resource)
            resource[c("service", "release", "release_kind", "query_date", "query_summary")]),
        limits = list(
            external_resources = "services queried inside a worker process are not recorded (KIT-026); the list covers queries made by the controlling process",
            legacy_steps = "a step committed before step provenance was recorded is `legacy` and lists its arguments without provenance"))
}

#' Write the analysis record as JSON. Doubles are written with 17 significant digits so every value restores
#' exactly.
#' @param manager On-disk DIA owner.
#' @param path Output file.
#' @return The path, invisibly.
#' @noRd
writeWorkflowAnalysisRecord <- function(manager, path, now = Sys.time()) {
    record <- buildWorkflowAnalysisRecord(manager, now)
    jsonlite::write_json(record, path, auto_unbox = TRUE, digits = I(17), null = "null", na = "null", pretty = TRUE)
    invisible(path)
}

#' Check an analysis record's fields and types, and report what is incomplete without failing.
#' @param record Record list, as built or as read with `jsonlite::read_json(simplifyVector = FALSE)`.
#' @return `list(valid, errors, incomplete)`: `valid` is FALSE only for a malformed record; legacy steps, steps
#'   with a missing declared parameter, absent inputs and the external-resources limit are listed in `incomplete`.
#' @noRd
validateWorkflowAnalysisRecord <- function(record) {
    errors <- character()
    incomplete <- character()
    need <- function(condition, message) if (!isTRUE(condition)) errors <<- c(errors, message)
    scalar <- function(x) (is.character(x) || is.numeric(x) || is.logical(x)) && length(x) == 1L
    need(is.list(record), "record is not a list")
    if (!is.list(record)) return(list(valid = FALSE, errors = errors, incomplete = incomplete))
    need(identical(record$schema, "multischolar.analysis_record.v1"), "schema")
    need(scalar(record$generated_at), "generated_at")
    need(is.list(record$store_identity) && scalar(record$store_identity$project_id), "store_identity")
    need(is.list(record$software) && scalar(record$software$package), "software")
    need(is.list(record$inputs), "inputs")
    need(is.list(record$steps), "steps")
    need(is.list(record$external_resources), "external_resources")
    need(is.list(record$limits) && scalar(record$limits$external_resources), "limits")
    for (step in record$steps) {
        need(is.list(step) && is.numeric(step$order) && scalar(step$state) && scalar(step$operation) &&
            scalar(step$step_label) && nzchar(step$step_label) && is.list(step$parameters) && is.list(step$products),
            paste0("step ", step$state %||% "?"))
        if (identical(step$operation, "legacy")) incomplete <- c(incomplete, paste0("legacy step: ", step$state))
        if (isTRUE(step$stored_as_object)) incomplete <- c(incomplete, paste0("state stored as a whole object, arguments not read: ", step$state))
        for (parameter in step$parameters) {
            need(is.list(parameter) && scalar(parameter$name) && scalar(parameter$source) && scalar(parameter$provenance),
                paste0("parameter in ", step$state))
            if (identical(parameter$source, "missing"))
                incomplete <- c(incomplete, paste0("declared parameter not recorded: ", step$state, "/", parameter$name))
        }
    }
    if (!length(record$inputs)) incomplete <- c(incomplete, "no inputs recorded")
    incomplete <- c(incomplete, "external resources: controller-process queries only")
    list(valid = !length(errors), errors = errors, incomplete = incomplete)
}
