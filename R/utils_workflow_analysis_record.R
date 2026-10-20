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
    if (is.character(mapping)) return(switch(mapping, fixed = "fixed", input = "input",
        no_shipped_default = "no_shipped_default", no_default = "no_shipped_default", "unknown"))
    shipped <- shipped_config[[mapping$section]][[mapping$key]]
    if (is.null(shipped)) return("unknown")
    # KIT-038 (IM-02): a requested value that is NULL or NA was not set; it is neither the default nor a user value.
    if (is.null(value) || !length(value) || all(is.na(value))) return("not_set")
    if (workflowParameterMatchesDefault(value, shipped)) "default" else "user_set"
}

#' The closed provenance vocabulary of a record parameter (KIT-038, IM-03). `no_default` is the KIT-027 spelling of
#' `no_shipped_default`, accepted on read and never written.
#' @noRd
workflowRecordProvenanceValues <- c("default", "user_set", "not_set", "no_shipped_default", "fixed", "input", "unknown")
workflowRecordProvenanceAliases <- c(no_default = "no_shipped_default")

#' Realised outcomes of one step (KIT-038, IM-02): the dimensions every step has, whatever the commit put under
#' `config$outcomes` or `config$skipped`, and for the normalisation step the realised RUV result.
#' @noRd
workflowRecordStepOutcomes <- function(context, state, operation, record, sizes) {
    config <- record$config
    outcomes <- list(rows = sizes$rows %||% NA_integer_, columns = sizes$columns %||% NA_integer_)
    if (is.list(config$outcomes)) outcomes[names(config$outcomes)] <- config$outcomes
    if (!is.null(config$skipped)) outcomes$skipped <- isTRUE(config$skipped)
    if (identical(operation, "normalization-ruv")) {
        applied <- identical(state, "ruv_corrected")
        count <- NA_integer_
        product <- record$data$products$normalization_science
        if (applied && isSavedWorkflowReference(product)) {
            hydrated <- tryCatch(hydrateSavedWorkflowEntity(context$store, product), error = function(error) NULL)
            # A product hydrates as an entity whose `value` is the kernel's result; a plain saved value is the list itself.
            saved <- if (is.list(hydrated$value)) hydrated$value else hydrated
            controls <- saved$controls
            if (is.logical(controls)) count <- sum(controls, na.rm = TRUE)
            else if (is.numeric(controls) && !all(is.na(controls))) count <- length(controls[!is.na(controls)])
        }
        outcomes$ruv_applied <- applied
        outcomes$k <- if (applied) config$ruv_k %||% NA_real_ else NA_real_
        outcomes$negative_control_percentage <- if (applied) config$percentage_as_neg_ctrl %||% NA_real_ else NA_real_
        outcomes$negative_control_count <- if (applied) count else NA_integer_
    }
    outcomes
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
        if (length(names_declared))
            parameters <- workflowRecordDeclaredParameters(names_declared, entry, section, args, config, provenance$values, shipped)
        # KIT-038 (IM-01): fixed constants are listed for every registered step, declared parameters or not; only a
        # `legacy` step (no provenance) lists its recorded arguments whole.
        fixed <- provenance$fixed
        fixed_source <- "fixed"
        if (!length(fixed) && !is.null(entry$declaration$fixed)) {
            # The head entry predates the constant's declaration (a store written before KIT-038 declared the P14
            # backend); the kernel passed it then as now, so the registry's constant is listed with its source named.
            fixed <- entry$declaration$fixed
            fixed_source <- "registry"
        }
        for (name in names(fixed)) parameters[[length(parameters) + 1L]] <- list(name = name, source = fixed_source,
            value = encodeWorkflowRecordValue(fixed[[name]]), provenance = "fixed")
        if (identical(operation, "legacy") && !is.null(args)) {
            for (name in names(args)) parameters[[length(parameters) + 1L]] <- list(name = name, source = "args",
                value = encodeWorkflowRecordValue(args[[name]]), provenance = "unknown")
        }
        previous <- record$previous_state
        input <- if (is.character(previous) && length(previous) == 1L && !is.na(previous)) snapshot$states[[previous]]$entity_id else NULL
        products <- lapply(unname(record$data$products %||% list()), function(product) {
            if (!isSavedWorkflowReference(product)) return(NULL)
            list(entity_id = product$entity_id, kind = product$kind,
                product_type = product$identity$selectors %||% NA_character_,
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
            # KIT-042 audit F-01: an `adopted` entry (written at ingest) is treated as a backfill: the step's software is unknown.
            backfilled <- isTRUE(provenance$backfilled) || isTRUE(provenance$adopted) || inherits(record$timestamp, "POSIXt") && !is.na(entry_time) &&
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
            software = list(version = software$version %||% NA_character_, sha = provenance$source_commit %||% NA_character_,
                build = provenance$build %||% NA_character_),
            defaults_sha256 = provenance$defaults_sha256 %||% NA_character_,
            dimensions = list(rows = sizes$rows %||% NA_integer_, columns = sizes$columns %||% NA_integer_),
            outcomes = workflowRecordStepOutcomes(context, state, operation, record, sizes),
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
    # KIT-038 (IM-09): the numerical environment of the build that wrote the last step; an environment entry written
    # by an earlier build carries none of these fields and reads "not recorded".
    environment_field <- function(name) {
        value <- current[[name]]
        if (is.null(value)) "not recorded" else value
    }
    step_defaults <- unique(vapply(steps, function(step) step$defaults_sha256, character(1)))
    list(schema = "multischolar.analysis_record.v1", generated_at = format(now, "%Y-%m-%dT%H:%M:%OS3%z"),
        workflow = snapshot$workflow_type %||% NA_character_,
        store_identity = list(project_id = context$store$identity$project_id %||% NA_character_,
            workflow_id = context$store$identity$workflow_id %||% NA_character_),
        software = list(package = "MultiScholaR", version = current$version %||% NA_character_,
            sha = current$source_commit %||% NA_character_,
            build = if (length(builds)) builds[[length(builds)]] else NA_character_,
            r_version = current$r_version %||% NA_character_,
            platform = current$platform %||% NA_character_,
            blas = environment_field("blas"), lapack = environment_field("lapack"),
            omp_threads = environment_field("omp_threads"), blas_threads = environment_field("blas_threads"),
            duckdb_threads = environment_field("duckdb_threads"),
            os_release = environment_field("os_release"), locale = environment_field("locale"),
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
    # KIT-038 (IM-03): the vocabulary is closed; `no_default` is read as the KIT-027 spelling of `no_shipped_default`.
    provenance_known <- function(x) is.character(x) && length(x) == 1L &&
        x %in% c(workflowRecordProvenanceValues, names(workflowRecordProvenanceAliases))
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
            need(provenance_known(parameter$provenance), paste0("provenance vocabulary in ", step$state, "/", parameter$name %||% "?"))
            if (identical(parameter$source, "missing"))
                incomplete <- c(incomplete, paste0("declared parameter not recorded: ", step$state, "/", parameter$name))
        }
        for (previous in step$preceded_by %||% list()) for (parameter in previous$parameters %||% list())
            need(provenance_known(parameter$provenance), paste0("provenance vocabulary in ", step$state, " (preceded_by)"))
    }
    for (analysis in record$analyses %||% list()) for (parameter in analysis$parameters %||% list())
        need(provenance_known(parameter$provenance), "provenance vocabulary in analyses")
    if (!length(record$inputs)) incomplete <- c(incomplete, "no inputs recorded")
    incomplete <- c(incomplete, "external resources: controller-process queries only")
    list(valid = !length(errors), errors = errors, incomplete = incomplete)
}

#' Decode a record step's parameters into the named list the operation receives: declared names only (a step
#' without declared parameters lists its arguments, which are not request parameters), fixed constants excluded.
#' @noRd
workflowReplayParameters <- function(parameters, declared) {
    out <- list()
    for (parameter in parameters) {
        if (identical(parameter$source, "fixed") || !parameter$name %in% declared) next
        out[parameter$name] <- list(decodeWorkflowRecordValue(parameter$value))
    }
    out
}

#' Replay an analysis record (KIT-028): the record alone, with hash-checked inputs, re-executes the analysis.
#'
#' The import and the GUI-built design are not registry operations: they are preconditions, verified (the
#' store holds exactly those states, with the record's entity ids) and never re-executed. Every later step
#' runs in process in `state_history` order with its recorded operation, variant and parameters: a worker
#' declaration through its registry handler and the controller's completion, a peptide step's progress
#' advanced as the GUI's display does; normalisation and correlation through `runSavedOperation` and the commit
#' their reference wrappers make; DA through `runDIAReferenceDaAnalysis` with the recorded settings and
#' contrasts (publication deferred). The store must hold nothing else, so no step can reuse a saved result.
#' @param record Analysis record (as built, or read with `jsonlite::read_json(simplifyVector = FALSE)`).
#' @param project_root Store root of the prepared project (the directory holding `.workflow-entities`).
#' @param inputs Named character vector, record input role to the file to check (`diann_report`, `fasta`, ...).
#' @param paths Experiment paths of the prepared project.
#' @param scratch Directory for request and publication scratch.
#' @param on_step Optional `function(step, manager)` called after each replayed step.
#' @return `list(states, project_root)` after closing the owner.
#' @noRd
replayWorkflowAnalysisRecord <- function(record, project_root, inputs, paths, scratch, worker_cap = 12 * 1024^3,
    on_step = NULL) {
    check <- validateWorkflowAnalysisRecord(record)
    if (!check$valid) cli::cli_abort(c("the analysis record is malformed", check$errors), class = "workflow_replay_invalid")
    # Inputs first: a record whose files do not match the project is refused before any step runs.
    for (input in record$inputs) {
        if (is.null(input$sha256) || is.na(input$sha256) || !input$role %in% c("diann_report", "maxquant_report",
            "fragpipe_report", "pd_tmt_report", "spectronaut_report", "fasta", "design_file")) next
        path <- if (input$role %in% names(inputs)) inputs[[input$role]] else NULL
        if (is.null(path) || !file.exists(path) || !identical(digest::digest(file = path, algo = "sha256"), input$sha256)) {
            cli::cli_abort("replay input {.val {input$role}} does not match the record's sha256", class = "workflow_replay_inputs")
        }
    }
    identity <- readSavedStoreIdentity(file.path(project_root, ".workflow-entities", "store.json"))
    manager <- newDIAStateAdapter(newWorkflowEntityStore(project_root, identity$project_id, identity$workflow_id), reopen = TRUE)
    on.exit(closeDIAStateAdapter(manager), add = TRUE)
    manager$.context$session_format <- "native_v1"
    manager$.context$reference_states <- TRUE
    manager$.context$experiment_paths <- paths
    steps <- record$steps
    operations <- vapply(steps, function(step) step$operation, character(1))
    preconditions <- seq_len(sum(cumprod(operations %in% c("import", "design_builder"))))
    history <- unlist(manager$.context$snapshot$state_history)
    expected <- vapply(steps[preconditions], function(step) step$state, character(1))
    if (!length(preconditions) || !identical(history, expected) || !all(vapply(steps[preconditions], function(step)
        identical(manager$.context$snapshot$states[[step$state]]$entity_id, step$result_entity_id), logical(1)))) {
        cli::cli_abort("the replay store does not hold exactly the record's import and design states", class = "workflow_replay_precondition")
    }
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    workflow$qc_params <- list(peptide_qc = list())
    progress <- new.env(parent = emptyenv())
    directories <- setupDirectories(file.path(scratch, "progress"), "proteomics")
    directories$proteomics_dia_baseline_smoke <- directories$proteomics
    normalization_final <- NULL
    for (step in steps[-preconditions]) {
        operation <- step$operation
        entry <- if (operation %in% c("legacy", "import", "design_builder")) NULL else
            tryCatch(workflowOperationEntry(operation), error = function(e) NULL)
        if (is.null(entry)) cli::cli_abort("step {.val {step$state}} cannot be replayed: operation {.val {operation}}",
            class = "workflow_replay_unsupported")
        parameters <- workflowReplayParameters(step$parameters, entry$parameters)
        variant <- if (is.null(step$variant) || is.na(step$variant)) entry$default_variant else step$variant
        withDIAStateScope(manager, function() {
            if (entry$compute %in% c("peptide", "protein")) {
                directory <- file.path(scratch, "requests", step$state)
                dir.create(directory, recursive = TRUE, showWarnings = FALSE)
                request <- list(operation = operation, input_references = list(state = getDIAStateReference(manager)),
                    parameters = parameters, worker_cap = worker_cap, directory = directory)
                value <- get(entry$handler, mode = "function")(manager, request, variant)
                value$operation <- operation
                setDIAStepProvenance(manager, operation, variant, entry$parameters, entry$declaration$fixed, values = parameters)
                # The controller's completion (`completeWorkflowAction`) on the replay's one workflow owner.
                workflow$filtering_progress <- progress$filtering_progress
                if (entry$completion == "protein") completeDIAProteinAction(workflow, value)
                else if (identical(value$section, "srlQvalueProteotypicPeptideClean")) finishDIAReferenceQvalueStep(workflow, value)
                else finishDIAReferencePeptideStep(workflow, value$checkpoint, value$state_name, value$description, value$config,
                    value$section, value$counts, observed_inclusion_absence = if (isTRUE(entry$inclusion_absence))
                        value$summary$inclusion_entries_absent_from_data else NULL)
                if (entry$completion == "peptide") {
                    # The GUI's display advances the session's filtering progress with the step's summary; the next
                    # step's progress product records it as prior progress.
                    record_now <- readDIAStateMetadata(manager)
                    product <- hydrateSavedWorkflowEntity(manager$.context$store, workflowPeptideProgressProduct(record_now$data))$value
                    grDevices::pdf(NULL)
                    device <- grDevices::dev.cur()
                    on.exit(if (device %in% grDevices::dev.list()) grDevices::dev.off(device), add = TRUE)
                    updateProteinFiltering(product$filteringSummary, entry$step_label, omic_type = "proteomics",
                        experiment_label = "dia_baseline_smoke", overwrite = TRUE, return_grid = FALSE,
                        project_dirs = directories, progress_owner = progress)
                }
            } else if (identical(operation, "normalization-ruv")) {
                # As `runDIAReferenceNormalization` commits: the operation's final checkpoint with its product.
                run <- runSavedOperation(manager, diaNormalizationOperation(), getDIAStateReference(manager), parameters)
                saved <- run$value
                skipped <- identical(parameters$ruv_mode, "skip")
                final <- readDIAStateMetadata(manager, saved$final)$data
                config <- list(norm_method = parameters$norm_method, ruv_mode = parameters$ruv_mode,
                    ruv_k = saved$best_k, percentage_as_neg_ctrl = saved$percentage)
                if (skipped) config$ruv_applied <- FALSE
                final$products$normalization_science <- run$reference
                commitDIAReferenceState(manager, if (skipped) "normalized" else "ruv_corrected", final, config,
                    if (skipped) "Post-normalization complete: RUV-III skipped by user" else
                        "Post-normalization complete: RUV-III correction and missing value cleanup completed")
                normalization_final <<- saved$final
            } else if (identical(operation, "pearson-selection")) {
                # As `runDIAReferenceCorrelation` commits: pairs, then selection, on the normalisation's final checkpoint.
                previous <- Filter(function(p) identical(p$operation, "pearson-pairs"), step$preceded_by %||% list())
                if (!length(previous) || is.null(normalization_final))
                    cli::cli_abort("step {.val {step$state}} cannot be replayed without its pairs and normalisation", class = "workflow_replay_unsupported")
                pairs_parameters <- workflowReplayParameters(previous[[1L]]$parameters, workflowOperationEntry("pearson-pairs")$parameters)
                pairs <- runSavedOperation(manager, diaPearsonPairsOperation(), normalization_final, pairs_parameters)
                selection <- runSavedOperation(manager, diaPearsonSelectionOperation(), normalization_final,
                    list(pairs = pairs$reference$entity_id, threshold = parameters$threshold, correlations = pairs$value),
                    dependencies = list(pairs$reference))
                spec <- readDIAStateMetadata(manager, selection$value)$data
                commitDIAReferenceState(manager, "correlation_filtered", spec,
                    list(min_pearson_correlation_threshold = parameters$threshold), "Applied final sample correlation filter (chunk 28)")
            } else {
                cli::cli_abort("step {.val {step$state}} cannot be replayed: operation {.val {operation}}", class = "workflow_replay_unsupported")
            }
        })
        if (!identical(manager$current_state, step$state))
            cli::cli_abort("replayed step {.val {step$state}} committed {.val {manager$current_state}}", class = "workflow_replay_invalid")
        if (is.function(on_step)) on_step(step, manager)
    }
    da <- record$analyses$differential_abundance
    if (!is.null(da)) {
        settings <- list()
        for (parameter in da$parameters) settings[parameter$name] <- list(decodeWorkflowRecordValue(parameter$value))
        contrasts <- do.call(rbind, lapply(da$contrasts, function(row) as.data.frame(lapply(row, unlist), stringsAsFactors = FALSE)))
        da_data <- new.env(parent = emptyenv())
        da_data$current_s4_object <- getDIAStateReference(manager)
        da_data$contrasts_available <- contrasts$contrasts
        da_workflow <- new.env(parent = emptyenv())
        da_workflow$state_manager <- manager
        da_workflow$contrasts_tbl <- contrasts
        runDIAReferenceDaAnalysis(list(formula_string = settings$formula, da_q_val_thresh = settings$qvalue,
            treat_lfc_cutoff = settings$lfc), da_data, da_workflow, paths, contrasts = contrasts, publish = FALSE)
        if (is.function(on_step)) on_step(list(state = manager$current_state, operation = "DA-contrast-output"), manager)
    }
    list(states = unlist(manager$.context$snapshot$state_history), project_root = project_root)
}

#' Text of one encoded record value, in full (KIT-029): atomic values comma-joined, lists as compact JSON,
#' tables as their size and hash.
#' @noRd
workflowRecordValueText <- function(encoded) {
    if (is.null(encoded)) return("NULL")
    if (identical(encoded$type, "data.frame"))
        return(sprintf("table: %s rows, columns %s (sha256 %s)", encoded$rows, paste(unlist(encoded$columns), collapse = ", "), encoded$sha256))
    if (identical(encoded$type, "opaque")) return(paste0("<", paste(unlist(encoded$class), collapse = "/"), ">"))
    value <- decodeWorkflowRecordValue(encoded)
    if (is.list(value)) return(as.character(jsonlite::toJSON(value, auto_unbox = TRUE, digits = I(17), null = "null", na = "null")))
    if (!length(value)) return("")
    text <- if (is.double(value)) format(value, digits = 17, trim = TRUE) else as.character(value)
    text[is.na(value)] <- "NA"
    paste(if (is.null(names(value))) text else paste0(names(value), "=", text), collapse = ", ")
}

#' The record as the tables people read (KIT-029): one source for the rendered document and the workbook.
#' @param record Analysis record (built or read back from JSON).
#' @return Named list of data frames: `software`, `packages`, `inputs`, `design`, `steps`, `parameters`,
#'   `analyses`, `external_resources`.
#' @noRd
workflowAnalysisRecordTables <- function(record) {
    chr <- function(x) if (is.null(x) || !length(x)) NA_character_ else as.character(unlist(x)[[1L]])
    frame <- function(rows, columns) {
        if (!length(rows)) return(stats::setNames(as.data.frame(matrix(character(), 0L, length(columns))), columns))
        out <- do.call(rbind, lapply(rows, function(row) as.data.frame(lapply(row[columns], chr), stringsAsFactors = FALSE)))
        names(out) <- columns
        out
    }
    threads <- record$software$duckdb_threads
    software <- data.frame(field = c("package", "version", "source commit", "build", "R", "platform", "BLAS", "LAPACK",
            "OpenMP threads", "BLAS threads", "DuckDB threads", "OS", "locale", "defaults sha256", "defaults changed", "generated"),
        value = c(chr(record$software$package), chr(record$software$version), chr(record$software$sha), chr(record$software$build),
            chr(record$software$r_version), chr(record$software$platform), chr(record$software$blas), chr(record$software$lapack),
            chr(record$software$omp_threads), chr(record$software$blas_threads),
            if (is.list(threads)) paste(names(threads), unlist(threads), sep = " = ", collapse = "; ") else chr(threads),
            chr(record$software$os_release), chr(record$software$locale),
            chr(record$defaults_sha256), chr(record$defaults_changed), chr(record$generated_at)),
        stringsAsFactors = FALSE)
    packages <- frame(record$software$packages %||% list(), c("name", "version"))
    inputs <- frame(record$inputs, c("role", "basename", "bytes", "sha256", "entity_id"))
    design <- if (length(record$design$rows)) as.data.frame(lapply(as.data.frame(do.call(rbind, lapply(record$design$rows, function(row)
        vapply(row, chr, character(1)))), stringsAsFactors = FALSE), as.character), stringsAsFactors = FALSE) else data.frame()
    steps <- frame(lapply(record$steps, function(step) c(step[c("order", "state", "operation", "variant", "step_label", "committed_at")],
        list(description = step$description, software_sha = step$software$sha, rows = step$dimensions$rows, columns = step$dimensions$columns,
            preceded_by = paste(vapply(step$preceded_by %||% list(), function(p) chr(p$operation), ""), collapse = ", ")))),
        c("order", "state", "description", "operation", "variant", "step_label", "committed_at", "software_sha", "rows", "columns", "preceded_by"))
    label_of <- function(operation, name) workflowParameterLabels[[operation]][[name]] %||% name
    parameter_rows <- list()
    for (step in record$steps) {
        add <- function(parameters, operation) for (parameter in parameters) parameter_rows[[length(parameter_rows) + 1L]] <<- list(
            order = step$order, step = step$state, operation = operation, name = parameter$name,
            label = label_of(operation, parameter$name),
            value = workflowRecordValueText(parameter$value), source = parameter$source, provenance = parameter$provenance)
        for (previous in step$preceded_by %||% list()) add(previous$parameters, previous$operation)
        add(step$parameters, step$operation)
    }
    parameters <- frame(parameter_rows, c("order", "step", "operation", "name", "label", "value", "source", "provenance"))
    analyses <- frame(unlist(lapply(names(record$analyses %||% list()), function(name) lapply(record$analyses[[name]]$parameters, function(parameter)
        list(analysis = name, operation = record$analyses[[name]]$operation, source_state = record$analyses[[name]]$source_state,
            name = parameter$name, value = workflowRecordValueText(parameter$value), provenance = parameter$provenance))), recursive = FALSE),
        c("analysis", "operation", "source_state", "name", "value", "provenance"))
    external <- frame(record$external_resources, c("service", "release", "release_kind", "query_date", "query_summary"))
    list(software = software, packages = packages, inputs = inputs, design = design, steps = steps, parameters = parameters,
        analyses = analyses, external_resources = external)
}

#' Write the record's workbook (KIT-029): sheets `software`, `inputs`, `design`, `parameters`, `external_resources`
#' (and `steps`, `analyses`, `packages`), one row per step and parameter in `parameters`.
#' @noRd
writeWorkflowAnalysisRecordWorkbook <- function(record, path) {
    tables <- workflowAnalysisRecordTables(record)
    openxlsx::write.xlsx(tables[c("software", "inputs", "design", "parameters", "external_resources", "steps", "analyses", "packages")],
        path, overwrite = TRUE)
    invisible(path)
}

#' Render the record for people (KIT-029) with the shared template, not `RenderReport` (which reads per-omic
#' templates and takes no parameters).
#' @param record_path Analysis record JSON.
#' @param output_dir Output directory.
#' @param format `docx` (the omic reports' `word_document`) or `html`.
#' @return The rendered file.
#' @noRd
renderWorkflowAnalysisRecord <- function(record_path, output_dir, format = c("docx", "html")) {
    format <- match.arg(format)
    # `rmarkdown::render` changes directory before it reads its arguments, so both paths are made absolute here.
    record_path <- normalizePath(record_path, mustWork = TRUE)
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    output_dir <- normalizePath(output_dir, mustWork = TRUE)
    template <- system.file("reports", "shared", "analysis_record_SI.rmd", package = "MultiScholaR")
    if (!nzchar(template)) cli::cli_abort("the analysis record template is not installed")
    work <- tempfile("analysis-record-render-")
    dir.create(work)
    on.exit(unlink(work, recursive = TRUE), add = TRUE)
    file.copy(template, work)
    rmarkdown::render(file.path(work, basename(template)), output_format = if (format == "docx") "word_document" else "html_document",
        output_file = paste0("analysis_record.", format), output_dir = output_dir, intermediates_dir = work,
        params = list(record_path = record_path), envir = new.env(parent = globalenv()), quiet = TRUE)
}

#' Export the record bundle (KIT-029): `analysis_record.json`, `.docx` and `.xlsx` in one directory.
#' @noRd
exportWorkflowAnalysisRecordBundle <- function(manager, directory) {
    if (!dir.exists(directory)) dir.create(directory, recursive = TRUE)
    json <- file.path(directory, "analysis_record.json")
    writeWorkflowAnalysisRecord(manager, json)
    record <- jsonlite::read_json(json, simplifyVector = FALSE)
    writeWorkflowAnalysisRecordWorkbook(record, file.path(directory, "analysis_record.xlsx"))
    docx <- renderWorkflowAnalysisRecord(json, directory, "docx")
    invisible(c(json = json, docx = docx, xlsx = file.path(directory, "analysis_record.xlsx")))
}

#' Where a report variable comes from (KIT-033): the record's step operation and parameter (`DA` names the
#' differential abundance analysis), and the `study_parameters.txt` section and key the parser reads.
#' @noRd
workflowRecordReportSources <- list(
    qvalue_threshold = list(record = c("peptide_qvalue_filter", "qvalue_threshold"), file = c("srlQvalueProteotypicPeptideClean", "qvalue_threshold")),
    proteotypic_peptide_only = list(record = c("peptide_qvalue_filter", "proteotypic_only"), file = c("srlQvalueProteotypicPeptideClean", "choose_only_proteotypic_peptide")),
    intensity_cutoff_percentile = list(record = c("peptide_intensity_filter", "peptides_intensity_cutoff_percentile"), file = c("peptideIntensityFiltering", "peptides_intensity_cutoff_percentile")),
    proportion_samples_below_cutoff = list(file = c("peptideIntensityFiltering", "peptides_proportion_of_samples_below_cutoff")),
    missing_values_cutoff = list(file = c("removeRowsWithMissingValuesPercent", "groupwise_percentage_cutoff")),
    max_groups_percentage_cutoff = list(file = c("removeRowsWithMissingValuesPercent", "max_groups_percentage_cutoff")),
    min_peptidoforms_per_protein = list(record = c("peptide_protein_evidence_filter", "num_peptidoforms_per_protein_thresh"), file = c("filterMinNumPeptidesPerProtein", "peptidoforms_per_protein_cutoff")),
    min_peptides_per_sample = list(record = c("minimum_peptides_per_sample", "peptides_per_sample_cutoff"), file = c("filterMinNumPeptidesPerSample", "peptides_per_sample_cutoff")),
    normalisation_method = list(record = c("normalization-ruv", "norm_method"), file = c("normaliseBetweenSamples", "method")),
    ruv_number_k = list(record = c("normalization-ruv", "ruv_k"), file = c("ruvIII_C_Varying", "ruv_number_k")),
    s4_de_q_val_thresh = list(record = c("DA", "qvalue"), file = c("deAnalysisParameters", "de_q_val_thresh")),
    s4_de_treat_lfc_cutoff = list(record = c("DA", "lfc"), file = c("deAnalysisParameters", "treat_lfc_cutoff")),
    s4_formula_string = list(record = c("DA", "formula"), file = c("deAnalysisParameters", "formula_string")),
    s4_ebayes_trend = list(file = c("deAnalysisParameters", "eBayes_trend")),
    s4_ebayes_robust = list(file = c("deAnalysisParameters", "eBayes_robust")),
    # KIT-039 (IM-06): the software the record names, and the realised RUV outcome (IM-02) in the words the
    # report's RUV sentence uses (`best_k`, `best_percentage`, `num_neg_ctrl`, `ruv_applied`).
    software_version = list(record = c("software", "version")),
    source_commit = list(record = c("software", "sha")),
    build = list(record = c("software", "build")),
    r_version = list(record = c("software", "r_version")),
    blas = list(record = c("software", "blas")),
    ruv_applied = list(record = c("outcomes", "normalization-ruv", "ruv_applied")),
    best_k = list(record = c("outcomes", "normalization-ruv", "k")),
    best_percentage = list(record = c("outcomes", "normalization-ruv", "negative_control_percentage")),
    num_neg_ctrl = list(record = c("outcomes", "normalization-ruv", "negative_control_count")),
    # KIT-043: the LFQ and TMT templates' variables: the protein intensity filter (the registry step
    # `protein_intensity_filter`; its groupwise and cross-group cutoffs are derived in the state's configuration,
    # not declared, so they stay file-only) and the DA settings the templates print from the UI section.
    proteins_intensity_cutoff_percentile = list(record = c("protein_intensity_filter", "proteins_intensity_cutoff_percentile"),
        file = c("removeRowsWithMissingValuesPercent", "proteins_intensity_cutoff_percentile")),
    groupwise_percentage_cutoff = list(file = c("removeRowsWithMissingValuesPercent", "groupwise_percentage_cutoff")),
    max_groups_percentage_cutoff = list(file = c("removeRowsWithMissingValuesPercent", "max_groups_percentage_cutoff")),
    min_reps_per_group = list(record = c("protein_intensity_filter", "min_reps_per_group")),
    min_groups = list(record = c("protein_intensity_filter", "min_groups")),
    ui_de_q_value_threshold = list(record = c("DA", "qvalue"), file = c("deAnalysisParameters", "de_q_val_thresh")),
    ui_de_log_fold_change_cutoff = list(record = c("DA", "lfc"), file = c("deAnalysisParameters", "treat_lfc_cutoff")))

#' Overlay a report's parameters with the analysis record (KIT-033). Without a record the parameters are
#' returned unchanged. With one, every report variable the record carries takes the record's value, and the
#' variables neither the record nor `study_parameters.txt` carries are named in `record_defaults_note`: the
#' report shows its default for them.
#' @param record_path `analysis_record.json` beside `study_parameters.txt`.
#' @param params The report's extracted parameters (`ExtractStudyParametersNew`).
#' @param study_params The parsed `study_parameters.txt` (`ReadStudyParamsNew`).
#' @noRd
workflowRecordReportParameters <- function(record_path, params, study_params = list()) {
    if (!file.exists(record_path)) return(params)
    record <- jsonlite::read_json(record_path, simplifyVector = FALSE)
    from_record <- function(operation, name) {
        parameters <- if (identical(operation, "DA")) record$analyses$differential_abundance$parameters else {
            steps <- Filter(function(step) identical(step$operation, operation), record$steps)
            if (length(steps)) steps[[length(steps)]]$parameters else list()
        }
        hit <- Filter(function(parameter) identical(parameter$name, name), parameters)
        if (length(hit)) decodeWorkflowRecordValue(hit[[1L]]$value) else NULL
    }
    from_outcomes <- function(operation, name) {
        steps <- Filter(function(step) identical(step$operation, operation), record$steps)
        if (!length(steps)) return(NULL)
        steps[[length(steps)]]$outcomes[[name]]
    }
    file_params <- study_params[["Workflow Parameters"]][["Parameters from Final S4 Object"]] %||% list()
    shown <- character()
    for (variable in names(workflowRecordReportSources)) {
        source <- workflowRecordReportSources[[variable]]
        kind <- if (is.null(source$record)) "file" else if (identical(source$record[[1L]], "software")) "software" else
            if (identical(source$record[[1L]], "outcomes")) "outcomes" else "parameter"
        value <- switch(kind,
            software = record$software[[source$record[[2L]]]],
            outcomes = from_outcomes(source$record[[2L]], source$record[[3L]]),
            parameter = from_record(source$record[[1L]], source$record[[2L]]),
            NULL)
        if (!is.null(value) && is.atomic(value) && length(value) == 1L && !is.na(value)) {
            current <- params[[variable]]
            # A parameter keeps the type the report's parser gave it (KIT-033); a software field or an outcome is
            # taken as the record has it, because the parser's value for those is often a parse artefact (a
            # logical TRUE for "Best k value" on the real journey) that would coerce the number away.
            params[[variable]] <- if (kind %in% c("software", "outcomes")) value else
                if (is.logical(current)) as.logical(value) else if (is.numeric(current)) as.numeric(value) else value
        } else if (kind %in% c("software", "outcomes")) {
            # A record without the field (written before KIT-038) leaves the report's own value, or "not recorded".
            if (is.null(params[[variable]])) params[[variable]] <- "not recorded"
        } else if (!is.null(source$file) && is.null(file_params[[source$file[[1L]]]][[source$file[[2L]]]]) && !is.null(params[[variable]])) {
            shown <- c(shown, sprintf("%s = %s", variable, paste(format(params[[variable]]), collapse = ", ")))
        }
    }
    params$record_present <- TRUE
    params$record_defaults_note <- if (length(shown)) paste0(" Parameter not recorded; default shown: ", paste(shown, collapse = "; "), ".") else ""
    params
}
