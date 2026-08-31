proteomicsPublicationPrivateRoot <- function() {
    publicationPath(
        ".omics-publication-workloads",
        "proteomics",
        "private"
    )
}

proteomicsPublicationSaltPath <- function(root = proteomicsPublicationPrivateRoot()) {
    file.path(root, "fingerprint-salt")
}

proteomicsPublicationCreateSalt <- function(
    path = proteomicsPublicationSaltPath()
) {
    if (file.exists(path)) {
        proteomicsPublicationAbort("private salt already exists", "privacy_error")
    }
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    connection <- file("/dev/urandom", open = "rb", raw = TRUE)
    on.exit(close(connection), add = TRUE)
    bytes <- readBin(connection, what = "raw", n = 32L)
    if (length(bytes) != 32L) {
        proteomicsPublicationAbort("private salt generation failed", "privacy_error")
    }
    salt <- paste(format(bytes), collapse = "")
    writeLines(salt, path, useBytes = TRUE)
    Sys.chmod(path, mode = "0600", use_umask = FALSE)
    if (as.integer(file.info(path)$mode) %% 512L != strtoi("600", base = 8L)) {
        proteomicsPublicationAbort("private salt permissions differ", "privacy_error")
    }
    invisible(path)
}

proteomicsPublicationReadSalt <- function(path = proteomicsPublicationSaltPath()) {
    if (!file.exists(path) || dir.exists(path)) {
        proteomicsPublicationAbort("private salt is unavailable", "privacy_error")
    }
    mode <- as.integer(file.info(path)$mode) %% 512L
    salt <- readLines(path, n = 1L, warn = FALSE)
    if (length(salt) != 1L || !grepl("^[0-9a-f]{64}$", salt) ||
        mode != strtoi("600", base = 8L)) {
        proteomicsPublicationAbort("private salt custody differs", "privacy_error")
    }
    salt
}

proteomicsPublicationPrivateFileState <- function(path) {
    if (!file.exists(path) || dir.exists(path)) {
        proteomicsPublicationAbort("private source is unreadable", "privacy_error")
    }
    list(
        size_bytes = as.numeric(file.info(path)$size),
        modified = as.numeric(file.info(path)$mtime),
        sha256 = publicationFileDigest(path)
    )
}

proteomicsPublicationSaltedFingerprint <- function(state, salt, role) {
    digest::digest(
        paste(salt, role, state$sha256, state$size_bytes, sep = ":"),
        algo = "sha256",
        serialize = FALSE
    )
}

proteomicsPublicationTsvHeaderIndexes <- function(header) {
    required <- c(
        "Protein.Group", "Stripped.Sequence", "Run",
        "Precursor.Quantity", "Precursor.Normalised"
    )
    missing <- setdiff(required, header)
    if (length(missing)) {
        proteomicsPublicationAbort(
            "private DIA-NN source lacks required columns",
            "private_schema_error"
        )
    }
    stats::setNames(match(required, header), required)
}

proteomicsPublicationHashSetAdd <- function(environment, values, salt) {
    for (value in unique(values[nzchar(values)])) {
        key <- digest::digest(
            paste0(salt, ":", value),
            algo = "sha256",
            serialize = FALSE
        )
        environment[[key]] <- TRUE
    }
    invisible(environment)
}

proteomicsPublicationInspectPrivateTsv <- function(
    path,
    salt,
    chunk_size = 5000L
) {
    before <- proteomicsPublicationPrivateFileState(path)
    connection <- file(path, open = "rt", encoding = "UTF-8")
    on.exit(close(connection), add = TRUE)
    header_line <- readLines(connection, n = 1L, warn = FALSE)
    if (length(header_line) != 1L || !nzchar(header_line)) {
        proteomicsPublicationAbort("private DIA-NN header is empty", "private_schema_error")
    }
    header <- strsplit(header_line, "\t", fixed = TRUE)[[1L]]
    indexes <- proteomicsPublicationTsvHeaderIndexes(header)
    runs <- new.env(hash = TRUE, parent = emptyenv())
    proteins <- new.env(hash = TRUE, parent = emptyenv())
    peptides <- new.env(hash = TRUE, parent = emptyenv())
    row_count <- 0L
    nonempty_quantity_count <- 0
    repeat {
        lines <- readLines(connection, n = chunk_size, warn = FALSE)
        if (!length(lines)) break
        fields <- strsplit(lines, "\t", fixed = TRUE)
        widths <- lengths(fields)
        if (any(widths != length(header))) {
            proteomicsPublicationAbort(
                "private DIA-NN source is nonrectangular",
                "private_schema_error"
            )
        }
        matrix <- do.call(rbind, fields)
        proteomicsPublicationHashSetAdd(
            runs,
            matrix[, indexes[["Run"]]],
            salt
        )
        proteomicsPublicationHashSetAdd(
            proteins,
            matrix[, indexes[["Protein.Group"]]],
            salt
        )
        proteomicsPublicationHashSetAdd(
            peptides,
            matrix[, indexes[["Stripped.Sequence"]]],
            salt
        )
        quantity <- matrix[, indexes[["Precursor.Normalised"]]]
        nonempty_quantity_count <- nonempty_quantity_count + sum(nzchar(quantity))
        row_count <- row_count + length(lines)
    }
    close(connection)
    on.exit(NULL, add = FALSE)
    after <- proteomicsPublicationPrivateFileState(path)
    if (!identical(before, after)) {
        proteomicsPublicationAbort("private DIA-NN source changed", "privacy_error")
    }
    list(
        role = "private_diann_report",
        row_count = as.integer(row_count),
        column_count = as.integer(length(header)),
        byte_size = before$size_bytes,
        unique_run_count = as.integer(length(ls(runs, all.names = TRUE))),
        unique_protein_group_count = as.integer(length(ls(
            proteins,
            all.names = TRUE
        ))),
        unique_peptide_count = as.integer(length(ls(peptides, all.names = TRUE))),
        nonempty_quantity_count = as.numeric(nonempty_quantity_count),
        required_schema_complete = TRUE,
        salted_source_fingerprint = proteomicsPublicationSaltedFingerprint(
            before,
            salt,
            "private_diann_report"
        )
    )
}

proteomicsPublicationFastaTerminalNewline <- function(path) {
    size <- as.numeric(file.info(path)$size)
    if (size == 0) return(FALSE)
    connection <- file(path, open = "rb")
    on.exit(close(connection), add = TRUE)
    seek(connection, where = size - 1, origin = "start")
    identical(readBin(connection, what = "raw", n = 1L), charToRaw("\n"))
}

proteomicsPublicationInspectPrivateFasta <- function(path, salt) {
    before <- proteomicsPublicationPrivateFileState(path)
    connection <- file(path, open = "rt", encoding = "UTF-8")
    on.exit(close(connection), add = TRUE)
    header_hashes <- new.env(hash = TRUE, parent = emptyenv())
    record_count <- 0L
    sequence_line_count <- 0L
    invalid_sequence_line_count <- 0L
    empty_record_count <- 0L
    current_has_sequence <- FALSE
    repeat {
        lines <- readLines(connection, n = 10000L, warn = FALSE)
        if (!length(lines)) break
        for (line in lines) {
            if (startsWith(line, ">")) {
                if (record_count > 0L && !current_has_sequence) {
                    empty_record_count <- empty_record_count + 1L
                }
                record_count <- record_count + 1L
                current_has_sequence <- FALSE
                token <- strsplit(sub("^>", "", line), "[[:space:]]+")[[1L]][[1L]]
                proteomicsPublicationHashSetAdd(header_hashes, token, salt)
            } else if (nzchar(line)) {
                sequence_line_count <- sequence_line_count + 1L
                current_has_sequence <- TRUE
                if (!grepl("^[ABCDEFGHIKLMNPQRSTVWYBXZJUO*-]+$", toupper(line))) {
                    invalid_sequence_line_count <- invalid_sequence_line_count + 1L
                }
            }
        }
    }
    if (record_count > 0L && !current_has_sequence) {
        empty_record_count <- empty_record_count + 1L
    }
    close(connection)
    on.exit(NULL, add = FALSE)
    after <- proteomicsPublicationPrivateFileState(path)
    if (!identical(before, after)) {
        proteomicsPublicationAbort("private FASTA source changed", "privacy_error")
    }
    list(
        role = "private_corresponding_fasta",
        byte_size = before$size_bytes,
        record_count = as.integer(record_count),
        unique_header_token_count = as.integer(length(ls(
            header_hashes,
            all.names = TRUE
        ))),
        duplicate_header_token_count = as.integer(
            record_count - length(ls(header_hashes, all.names = TRUE))
        ),
        sequence_line_count = as.integer(sequence_line_count),
        invalid_sequence_line_count = as.integer(invalid_sequence_line_count),
        empty_record_count = as.integer(empty_record_count),
        terminal_newline = proteomicsPublicationFastaTerminalNewline(path),
        syntax_valid = record_count > 0L && sequence_line_count > 0L &&
            invalid_sequence_line_count == 0L && empty_record_count == 0L,
        completeness_status = "unconfirmed_download",
        completeness_reason = "source_filename_has_crdownload_extension",
        salted_source_fingerprint = proteomicsPublicationSaltedFingerprint(
            before,
            salt,
            "private_corresponding_fasta"
        )
    )
}

proteomicsPublicationPrivateEnvelope <- function(tsv_path, fasta_path, salt_path) {
    salt <- proteomicsPublicationReadSalt(salt_path)
    list(
        schema = "multischolar.omics_publication_proteomics_private_envelope",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-065",
        status = "private_report_admissible_fasta_unconfirmed",
        salt_id = digest::digest(salt, algo = "sha256", serialize = FALSE),
        report = proteomicsPublicationInspectPrivateTsv(tsv_path, salt),
        fasta = proteomicsPublicationInspectPrivateFasta(fasta_path, salt),
        private_paths_retained = FALSE,
        identifiers_values_sequences_retained = FALSE,
        distribution_retained = FALSE,
        full_workflow_project_authority = FALSE,
        claim_scope = "private_project_specific_import_scale_non_promotional",
        publication_authority = FALSE
    )
}

proteomicsPublicationForbiddenSourceFields <- function(value) {
    if (!is.list(value)) return(invisible(TRUE))
    forbidden <- c(
        "path", "private_path", "header", "headers", "identifier",
        "identifiers", "sequence", "sequences", "values", "distribution",
        "quantiles", "source_sha256", "unsalted_fingerprint"
    )
    if (length(intersect(names(value), forbidden))) {
        proteomicsPublicationAbort("source authority leaks private fields", "privacy_error")
    }
    lapply(value, proteomicsPublicationForbiddenSourceFields)
    invisible(TRUE)
}

proteomicsPublicationValidatePrivateEnvelope <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "owner_ticket_id", "status", "salt_id",
        "report", "fasta", "private_paths_retained",
        "identifiers_values_sequences_retained", "distribution_retained",
        "full_workflow_project_authority", "claim_scope",
        "publication_authority"
    ), "Proteomics private envelope")
    proteomicsPublicationForbiddenSourceFields(record)
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_proteomics_private_envelope"
    ) && identical(record$schema_version, "1.0.0") &&
        identical(record$owner_ticket_id, "OMICS-ART-065") &&
        identical(
            record$status,
            "private_report_admissible_fasta_unconfirmed"
        ) && grepl("^[0-9a-f]{64}$", record$salt_id) &&
        isTRUE(record$report$required_schema_complete) &&
        isTRUE(record$fasta$syntax_valid) &&
        identical(record$fasta$completeness_status, "unconfirmed_download") &&
        !isTRUE(record$private_paths_retained) &&
        !isTRUE(record$identifiers_values_sequences_retained) &&
        !isTRUE(record$distribution_retained) &&
        !isTRUE(record$full_workflow_project_authority) &&
        identical(
            record$claim_scope,
            "private_project_specific_import_scale_non_promotional"
        ) && !isTRUE(record$publication_authority)
    if (!valid) proteomicsPublicationAbort("private envelope differs", "privacy_error")
    invisible(record)
}

proteomicsPublicationSourceFields <- function() {
    c(
        "source_id", "project_id", "capability_id", "owner_id",
        "independence_owner_id", "license_or_use_authority", "privacy_class",
        "source_kind", "schema_id", "byte_or_aggregate_receipt_digest",
        "selection_rule_id", "exclusion_rule_id", "support_tier",
        "authority_scope", "status", "counts_toward_cross_project"
    )
}

proteomicsPublicationValidateSource <- function(source) {
    publicationRequireNames(
        source,
        proteomicsPublicationSourceFields(),
        "Proteomics source authority"
    )
    text_fields <- setdiff(
        proteomicsPublicationSourceFields(),
        c("counts_toward_cross_project", "byte_or_aggregate_receipt_digest")
    )
    valid <- all(vapply(source[text_fields], publicationScalarString, logical(1))) &&
        source$capability_id %in% vapply(
            proteomicsPublicationCapabilities(),
            `[[`,
            character(1),
            "capability_id"
        ) && source$privacy_class %in% c(
            "public_real", "private_calibrated"
        ) && source$source_kind %in% c("public_real", "private_real") &&
        source$support_tier %in% c(
            "full_workflow", "import_scale_only", "parser_fixture_only"
        ) && source$status %in% c(
            "admissible", "private_report_admissible_fasta_unconfirmed",
            "rejected"
        ) && is.logical(source$counts_toward_cross_project) &&
        length(source$counts_toward_cross_project) == 1L
    proteomicsPublicationRequireDigest(
        source$byte_or_aggregate_receipt_digest,
        "Source receipt"
    )
    if (isTRUE(source$counts_toward_cross_project)) {
        valid <- valid && identical(source$status, "admissible") &&
            identical(source$support_tier, "full_workflow") &&
            identical(source$authority_scope, "independent_real_project")
    }
    if (!valid) proteomicsPublicationAbort("source authority differs")
    invisible(source)
}

proteomicsPublicationValidateSourceDecision <- function(decision, sources) {
    publicationRequireNames(decision, c(
        "capability_id", "required_real_project_count", "source_ids",
        "verified_real_project_count", "cross_project_claim_status",
        "current_claim_scope", "cross_project_source_ready",
        "promotion_eligible"
    ), "Proteomics source decision")
    selected <- Filter(\(source) {
        source$capability_id == decision$capability_id &&
            source$source_id %in% unlist(decision$source_ids, use.names = FALSE)
    }, sources)
    verified <- sum(vapply(
        selected,
        `[[`,
        logical(1),
        "counts_toward_cross_project"
    ))
    ready <- verified >= decision$required_real_project_count
    expected_status <- if (ready) {
        "cross_project_source_ready"
    } else {
        "insufficient_independent_projects"
    }
    expected_scope <- if (ready) {
        "cross_project_source_ready_for_confirmatory_review"
    } else {
        "project_specific_non_promotional"
    }
    valid <- decision$capability_id %in% vapply(
        proteomicsPublicationCapabilities(),
        `[[`,
        character(1),
        "capability_id"
    ) && identical(decision$required_real_project_count, 3L) &&
        identical(decision$verified_real_project_count, as.integer(verified)) &&
        identical(decision$cross_project_claim_status, expected_status) &&
        identical(decision$current_claim_scope, expected_scope) &&
        identical(decision$cross_project_source_ready, ready) &&
        identical(decision$promotion_eligible, ready)
    if (!valid) proteomicsPublicationAbort("source decision differs")
    invisible(decision)
}

proteomicsPublicationValidateSources <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "sources_id", "owner_ticket_id",
        "status", "projects_predecessor", "private_envelope_binding",
        "public_acquisition_binding",
        "minimum_real_projects", "generated_counts_as_real",
        "fixtures_count_as_real", "sources", "capability_decisions",
        "acquisition_outcomes", "unknown_source_policy",
        "publication_authority"
    ), "Proteomics source manifest")
    proteomicsPublicationValidateBinding(
        record$projects_predecessor,
        "Projects predecessor"
    )
    proteomicsPublicationValidateBinding(
        record$private_envelope_binding,
        "Private envelope"
    )
    private <- publicationReadJson(record$private_envelope_binding$path)
    proteomicsPublicationValidatePrivateEnvelope(private)
    proteomicsPublicationValidateBinding(
        record$public_acquisition_binding,
        "Public acquisition"
    )
    acquisition <- publicationReadJson(record$public_acquisition_binding$path)
    proteomicsPublicationValidateAcquisition(acquisition)
    lapply(record$sources, proteomicsPublicationValidateSource)
    ids <- vapply(record$sources, `[[`, character(1), "source_id")
    projects <- vapply(record$sources, `[[`, character(1), "project_id")
    decisions <- vapply(
        record$capability_decisions,
        `[[`,
        character(1),
        "capability_id"
    )
    lapply(
        record$capability_decisions,
        proteomicsPublicationValidateSourceDecision,
        sources = record$sources
    )
    expected_capabilities <- vapply(
        proteomicsPublicationCapabilities(),
        `[[`,
        character(1),
        "capability_id"
    )
    acquisition_outcomes <- acquisition$format_outcomes
    valid_outcomes <- identical(
        publicationObjectDigest(record$acquisition_outcomes),
        publicationObjectDigest(acquisition_outcomes)
    )
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_proteomics_sources"
    ) && identical(record$schema_version, "1.0.0") &&
        identical(record$owner_ticket_id, "OMICS-ART-065") &&
        identical(record$status, "source_inventory_frozen") &&
        identical(record$minimum_real_projects, 3L) &&
        !isTRUE(record$generated_counts_as_real) &&
        !isTRUE(record$fixtures_count_as_real) && !anyDuplicated(ids) &&
        !anyDuplicated(projects) && setequal(decisions, expected_capabilities) &&
        !anyDuplicated(decisions) && valid_outcomes &&
        identical(record$unknown_source_policy, "reject") &&
        !isTRUE(record$publication_authority)
    lapply(record$sources, proteomicsPublicationForbiddenSourceFields)
    lapply(record$acquisition_outcomes, proteomicsPublicationForbiddenSourceFields)
    if (!valid) proteomicsPublicationAbort("source manifest differs")
    invisible(record)
}

proteomicsPublicationPublicColumnNames <- function(path) {
    extension <- tolower(tools::file_ext(path))
    if (identical(extension, "parquet")) {
        if (!requireNamespace("arrow", quietly = TRUE)) {
            proteomicsPublicationAbort("Arrow is required for source schema")
        }
        return(arrow::open_dataset(path)$schema$names)
    }
    if (identical(extension, "xlsx")) {
        if (!requireNamespace("readxl", quietly = TRUE)) {
            proteomicsPublicationAbort("readxl is required for source schema")
        }
        return(names(readxl::read_excel(path, n_max = 0L)))
    }
    line <- readLines(path, n = 1L, warn = FALSE)
    if (length(line) != 1L) {
        proteomicsPublicationAbort("public source header is unavailable")
    }
    strsplit(line, "\t", fixed = TRUE)[[1L]]
}

proteomicsPublicationPublicSchemaDecision <- function(names, format) {
    if (identical(format, "diann")) {
        valid <- all(c(
            "Protein.Group", "Stripped.Sequence", "Run"
        ) %in% names) && any(c(
            "Precursor.Normalised", "Precursor.Quantity"
        ) %in% names)
        reason <- if (valid) NULL else "missing_diann_required_columns"
    } else if (identical(format, "maxquant")) {
        protein <- any(c(
            "Majority.protein.IDs", "Protein.IDs", "Protein.Ids",
            "Protein IDs", "Protein ID"
        ) %in% names)
        valid <- protein && any(startsWith(names, "LFQ.intensity."))
        reason <- if (valid) NULL else "current_parser_rejects_native_lfq_columns"
    } else if (identical(format, "fragpipe")) {
        protein <- any(tolower(names) %in% tolower(c(
            "Protein ID", "Protein.ID", "Protein_ID", "Protein"
        )))
        valid <- protein && any(grepl(
            "MaxLFQ.*Intensity$",
            names,
            ignore.case = TRUE
        ))
        reason <- if (valid) NULL else "maxlfq_columns_unavailable"
    } else if (identical(format, "pd_tmt")) {
        protein <- any(c("Accession", "Master.Protein.Accessions") %in% names)
        valid <- protein && any(grepl("^Abundance: |^Abundance\\.", names))
        reason <- if (valid) NULL else {
            "current_parser_rejects_native_pd_peptide_columns"
        }
    } else {
        proteomicsPublicationAbort("public source format is unsupported")
    }
    list(valid = valid, reason = reason)
}

proteomicsPublicationPrideFileRecord <- function(files, filename) {
    matches <- which(vapply(files, \(record) {
        identical(record$fileName, filename)
    }, logical(1)))
    if (length(matches) != 1L) {
        proteomicsPublicationAbort("PRIDE file metadata is not unique")
    }
    files[[matches]]
}

proteomicsPublicationPublicCandidate <- function(
    project_id,
    format,
    path,
    repository_filename,
    project_metadata_path,
    file_metadata_path
) {
    project <- publicationReadJson(project_metadata_path)
    files <- publicationReadJson(file_metadata_path)
    file <- proteomicsPublicationPrideFileRecord(files, repository_filename)
    columns <- proteomicsPublicationPublicColumnNames(path)
    schema <- proteomicsPublicationPublicSchemaDecision(columns, format)
    owners <- unlist(project$labPIs, use.names = FALSE)
    if (!length(owners)) owners <- unlist(project$submitters, use.names = FALSE)
    owner_digest <- digest::digest(
        paste(sort(owners), collapse = "|"),
        algo = "sha256",
        serialize = FALSE
    )
    valid_bytes <- identical(
        as.numeric(file$fileSizeBytes),
        as.numeric(file.info(path)$size)
    )
    if (!valid_bytes || !identical(project$accession, project_id) ||
        !identical(project$license, "Creative Commons Public Domain (CC0)")) {
        proteomicsPublicationAbort("public source metadata differs")
    }
    list(
        source_id = paste("pride", tolower(project_id), format, sep = "."),
        project_id = project_id,
        format = format,
        source_kind = "public_real",
        privacy_class = "public_real",
        independence_owner_id = paste0("pride-owner-sha256:", owner_digest),
        license_or_use_authority = "Creative Commons Public Domain (CC0)",
        repository = "PRIDE Archive",
        project_api = paste0(
            "https://www.ebi.ac.uk/pride/ws/archive/v2/projects/",
            project_id
        ),
        project_metadata_sha256 = publicationFileDigest(project_metadata_path),
        file_metadata_sha256 = publicationFileDigest(file_metadata_path),
        repository_file_accession = file$accession,
        repository_filename = repository_filename,
        file_size_bytes = as.numeric(file$fileSizeBytes),
        file_sha256 = publicationFileDigest(path),
        column_count = as.integer(length(columns)),
        parser_schema_status = if (schema$valid) {
            "format_schema_admissible"
        } else {
            "rejected_current_parser_schema"
        },
        parser_schema_reason = schema$reason,
        design_authority_status = "not_frozen",
        counts_toward_cross_project = FALSE,
        local_path_retained = FALSE,
        publication_authority = FALSE
    )
}

proteomicsPublicationValidatePublicCandidate <- function(record) {
    publicationRequireNames(record, c(
        "source_id", "project_id", "format", "source_kind", "privacy_class",
        "independence_owner_id", "license_or_use_authority", "repository",
        "project_api", "project_metadata_sha256", "file_metadata_sha256",
        "repository_file_accession", "repository_filename", "file_size_bytes",
        "file_sha256", "column_count", "parser_schema_status",
        "parser_schema_reason", "design_authority_status",
        "counts_toward_cross_project", "local_path_retained",
        "publication_authority"
    ), "Public proteomics acquisition")
    reason_valid <- if (identical(
        record$parser_schema_status,
        "format_schema_admissible"
    )) {
        is.null(record$parser_schema_reason)
    } else {
        publicationScalarString(record$parser_schema_reason)
    }
    valid <- publicationScalarString(record$source_id) &&
        grepl("^PXD[0-9]{6}$", record$project_id) &&
        record$format %in% names(proteomicsPublicationCapabilities()) &&
        identical(record$source_kind, "public_real") &&
        identical(record$privacy_class, "public_real") &&
        startsWith(record$independence_owner_id, "pride-owner-sha256:") &&
        identical(
            record$license_or_use_authority,
            "Creative Commons Public Domain (CC0)"
        ) && identical(record$repository, "PRIDE Archive") &&
        identical(
            record$project_api,
            paste0(
                "https://www.ebi.ac.uk/pride/ws/archive/v2/projects/",
                record$project_id
            )
        ) && grepl("^[0-9a-f]{64}$", record$repository_file_accession) &&
        publicationScalarString(record$repository_filename) &&
        is.numeric(record$file_size_bytes) && record$file_size_bytes > 0 &&
        is.numeric(record$column_count) && record$column_count > 0 &&
        record$parser_schema_status %in% c(
            "format_schema_admissible", "rejected_current_parser_schema"
        ) && reason_valid &&
        identical(record$design_authority_status, "not_frozen") &&
        !isTRUE(record$counts_toward_cross_project) &&
        !isTRUE(record$local_path_retained) &&
        !isTRUE(record$publication_authority)
    for (field in c(
        "project_metadata_sha256", "file_metadata_sha256", "file_sha256"
    )) {
        proteomicsPublicationRequireDigest(record[[field]], field)
    }
    if (!valid) proteomicsPublicationAbort("public acquisition differs")
    invisible(record)
}

proteomicsPublicationValidateAcquisition <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "acquisition_id", "owner_ticket_id",
        "status", "repository_authority", "candidates", "format_outcomes",
        "publication_authority"
    ), "Proteomics public acquisition")
    lapply(record$candidates, proteomicsPublicationValidatePublicCandidate)
    ids <- vapply(record$candidates, `[[`, character(1), "source_id")
    projects <- vapply(record$candidates, `[[`, character(1), "project_id")
    formats <- names(proteomicsPublicationCapabilities())
    outcomes <- vapply(
        record$format_outcomes,
        `[[`,
        character(1),
        "format"
    )
    valid_outcomes <- all(vapply(record$format_outcomes, \(outcome) {
        candidates <- Filter(\(candidate) {
            identical(candidate$format, outcome$format)
        }, record$candidates)
        schema_count <- sum(vapply(candidates, \(candidate) {
            identical(
                candidate$parser_schema_status,
                "format_schema_admissible"
            )
        }, logical(1)))
        setequal(names(outcome), c(
            "format", "reviewed_candidate_count", "schema_admissible_count",
            "full_workflow_design_count", "cross_project_source_count",
            "decision"
        )) && identical(
            outcome$reviewed_candidate_count,
            as.integer(length(candidates))
        ) && identical(outcome$schema_admissible_count, as.integer(schema_count)) &&
            identical(outcome$full_workflow_design_count, 0L) &&
            identical(outcome$cross_project_source_count, 0L) &&
            identical(outcome$decision, "project_specific_non_promotional")
    }, logical(1)))
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_proteomics_public_acquisition"
    ) && identical(record$schema_version, "1.0.0") &&
        identical(record$owner_ticket_id, "OMICS-ART-065") &&
        identical(record$status, "frozen_non_promotional") &&
        identical(record$repository_authority, "PRIDE Archive API v2") &&
        !anyDuplicated(ids) && !anyDuplicated(projects) &&
        setequal(outcomes, formats) && !anyDuplicated(outcomes) &&
        valid_outcomes && !isTRUE(record$publication_authority)
    if (!valid) proteomicsPublicationAbort("public acquisition record differs")
    invisible(record)
}
