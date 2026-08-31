#!/usr/bin/env Rscript

payloadAccessRepoRoot <- function(start = getwd()) {
    path <- normalizePath(start, mustWork = TRUE)
    repeat {
        if (file.exists(file.path(path, "DESCRIPTION")) &&
            dir.exists(file.path(path, "R"))) {
            return(path)
        }
        parent <- dirname(path)
        if (identical(parent, path)) stop("Cannot locate repository root")
        path <- parent
    }
}

payloadAccessSourceFields <- function() {
    c("data_tbl", "data_cln")
}

payloadAccessCanonical <- function(value) {
    if (!is.list(value)) return(value)
    if (!is.null(names(value))) {
        value <- value[order(names(value), method = "radix")]
    }
    lapply(value, payloadAccessCanonical)
}

payloadAccessObjectDigest <- function(value) {
    encoded <- jsonlite::toJSON(
        payloadAccessCanonical(value),
        auto_unbox = TRUE,
        null = "null",
        na = "null",
        digits = 17,
        pretty = FALSE
    )
    digest::digest(encoded, algo = "sha256", serialize = FALSE)
}

payloadAccessSlice <- function(lines, line1, col1, line2, col2) {
    if (line1 == line2) return(substr(lines[[line1]], col1, col2))
    values <- lines[seq.int(line1, line2)]
    values[[1L]] <- substring(values[[1L]], col1)
    values[[length(values)]] <- substr(values[[length(values)]], 1L, col2)
    paste(values, collapse = "\n")
}

payloadAccessNormalizeExpression <- function(value) {
    trimws(gsub("[[:space:]]+", " ", value))
}

payloadAccessDescendants <- function(parse_data, id) {
    found <- integer()
    frontier <- as.integer(id)
    repeat {
        children <- parse_data$id[parse_data$parent %in% frontier]
        children <- setdiff(children, found)
        if (!length(children)) break
        found <- c(found, children)
        frontier <- children
    }
    unique(found)
}

payloadAccessAncestors <- function(parse_data, id) {
    ancestors <- integer()
    current <- as.integer(id)
    repeat {
        row <- parse_data[parse_data$id == current, , drop = FALSE]
        if (!nrow(row) || identical(row$parent[[1L]], 0L)) break
        current <- row$parent[[1L]]
        ancestors <- c(ancestors, current)
    }
    ancestors
}

payloadAccessAssignmentKind <- function(parse_data, access_id) {
    assignment_tokens <- c("LEFT_ASSIGN", "RIGHT_ASSIGN", "EQ_ASSIGN")
    for (ancestor in payloadAccessAncestors(parse_data, access_id)) {
        children <- parse_data[parse_data$parent == ancestor, , drop = FALSE]
        operator <- children[children$token %in% assignment_tokens, , drop = FALSE]
        if (!nrow(operator)) next
        operator <- operator[order(operator$col1), , drop = FALSE][1L, ]
        expressions <- children[children$token == "expr", , drop = FALSE]
        if (!nrow(expressions)) return("read")
        if (identical(operator$token, "RIGHT_ASSIGN")) {
            left <- expressions$line1 > operator$line1 |
                (expressions$line1 == operator$line1 &
                    expressions$col1 > operator$col1)
        } else {
            left <- expressions$line2 < operator$line1 |
                (expressions$line2 == operator$line1 &
                    expressions$col2 < operator$col1)
        }
        left_ids <- expressions$id[left]
        for (left_id in left_ids) {
            if (access_id == left_id ||
                access_id %in% payloadAccessDescendants(parse_data, left_id)) {
                return("write")
            }
        }
        return("read")
    }
    "read"
}

payloadAccessTerminalRows <- function(parse_data, access_id, after_col) {
    ids <- c(access_id, payloadAccessDescendants(parse_data, access_id))
    terminal <- parse_data[
        parse_data$id %in% ids &
            parse_data$terminal &
            parse_data$col1 > after_col,
        ,
        drop = FALSE
    ]
    terminal[order(terminal$line1, terminal$col1), , drop = FALSE]
}

payloadAccessStripString <- function(value) {
    first <- substr(value, 1L, 1L)
    last <- substr(value, nchar(value), nchar(value))
    if (!first %in% c("\"", "'") || !identical(first, last)) return(value)
    substring(value, 2L, nchar(value) - 1L)
}

payloadAccessField <- function(parse_data, operator) {
    terminals <- payloadAccessTerminalRows(
        parse_data,
        operator$parent,
        operator$col1
    )
    if (operator$text %in% c("$", "@")) {
        fields <- terminals[terminals$token %in% c("SYMBOL", "SLOT"), ]
        if (!nrow(fields)) return(list(field = NA_character_, computed = TRUE))
        return(list(field = fields$text[[1L]], computed = FALSE))
    }
    indices <- terminals[terminals$token %in% c(
        "STR_CONST", "SYMBOL", "NUM_CONST"
    ), , drop = FALSE]
    if (!nrow(indices)) return(list(field = NA_character_, computed = TRUE))
    value <- indices[1L, ]
    list(
        field = payloadAccessStripString(value$text),
        computed = !identical(value$token, "STR_CONST")
    )
}

payloadAccessOperatorObject <- function(lines, expression, operator) {
    if (expression$line1 != operator$line1) return("multiline_object")
    value <- payloadAccessSlice(
        lines,
        expression$line1,
        expression$col1,
        operator$line1,
        operator$col1 - 1L
    )
    payloadAccessNormalizeExpression(value)
}

payloadAccessRelevantComputed <- function(object, operator) {
    identical(operator, "[[") && grepl(
        paste0(
            "^(workflow_data|workflowData|workflow_state|workflowState|",
            "state|data_bus)$"
        ),
        object
    )
}

payloadAccessScanFile <- function(path, repo_root, source_fields) {
    relative <- substring(normalizePath(path), nchar(repo_root) + 2L)
    lines <- readLines(path, warn = FALSE)
    parsed <- tryCatch(
        parse(path, keep.source = TRUE),
        error = function(error) error
    )
    if (inherits(parsed, "error")) {
        return(list(
            accesses = list(),
            errors = list(list(
                path = relative,
                condition = conditionMessage(parsed)
            ))
        ))
    }
    parse_data <- utils::getParseData(parsed, includeText = TRUE)
    operators <- parse_data[
        parse_data$token %in% c("'$'", "'@'", "LBB"),
        ,
        drop = FALSE
    ]
    accesses <- lapply(seq_len(nrow(operators)), function(index) {
        operator <- operators[index, ]
        expression <- parse_data[
            parse_data$id == operator$parent,
            ,
            drop = FALSE
        ]
        if (!nrow(expression)) return(NULL)
        field <- payloadAccessField(parse_data, operator)
        object <- payloadAccessOperatorObject(lines, expression, operator)
        relevant <- (!field$computed && field$field %in% source_fields) ||
            (field$computed && payloadAccessRelevantComputed(
                object,
                operator$text
            ))
        if (!relevant) return(NULL)
        expression_text <- payloadAccessSlice(
            lines,
            expression$line1,
            expression$col1,
            expression$line2,
            expression$col2
        )
        normalized <- payloadAccessNormalizeExpression(expression_text)
        access_kind <- payloadAccessAssignmentKind(
            parse_data,
            operator$parent
        )
        stable <- list(
            path = relative,
            expression = normalized,
            field = field$field,
            operator = operator$text,
            access_kind = access_kind,
            computed = field$computed
        )
        list(
            access_key = payloadAccessObjectDigest(stable),
            path = relative,
            expression = normalized,
            object = object,
            field = field$field,
            operator = operator$text,
            access_kind = access_kind,
            computed = field$computed,
            evidence_class = if (startsWith(relative, "R/")) {
                "production"
            } else if (startsWith(relative, "tests/")) {
                "test"
            } else {
                "other"
            },
            location = list(
                line = as.integer(expression$line1),
                column = as.integer(expression$col1)
            )
        )
    })
    list(
        accesses = Filter(Negate(is.null), accesses),
        errors = list()
    )
}

payloadAccessMerge <- function(accesses) {
    if (!length(accesses)) return(list())
    keys <- vapply(accesses, `[[`, character(1), "access_key")
    groups <- split(accesses, keys)
    merged <- lapply(groups, function(values) {
        first <- values[[1L]]
        first$locations <- lapply(values, `[[`, "location")
        first$location <- NULL
        first$occurrence_count <- as.integer(length(values))
        first
    })
    unname(merged[order(vapply(
        merged,
        `[[`,
        character(1),
        "access_key"
    ))])
}

payloadAccessReadAuthority <- function(path) {
    if (is.null(path)) return(NULL)
    jsonlite::read_json(path, simplifyVector = FALSE)
}

payloadAccessValidateAuthority <- function(authority) {
    required <- c(
        "schema", "schema_version", "authority_id", "owner_ticket_id",
        "discovery_digest", "owners", "adapter_obligations",
        "unknown_policy", "candidate_freeze_authority"
    )
    if (!is.list(authority) || !setequal(names(authority), required)) {
        stop("Payload-access owner authority fields differ")
    }
    owner_fields <- c(
        "access_key", "path", "expression", "access_kind", "computed",
        "computed_resolution", "owner_ticket_id", "owner_role",
        "consumer_ticket_ids", "rationale"
    )
    valid_owner <- vapply(authority$owners, function(owner) {
        is.list(owner) && setequal(names(owner), owner_fields) &&
            grepl("^[0-9a-f]{64}$", owner$access_key) &&
            owner$access_kind %in% c("read", "write") &&
            owner$owner_ticket_id %in% sprintf("OMICS-ART-%03d", 69:76) &&
            owner$computed_resolution %in% c(
                "exact_field",
                "contract_source_fields_guarded"
            ) &&
            identical(
                owner$computed_resolution,
                if (isTRUE(owner$computed)) {
                    "contract_source_fields_guarded"
                } else {
                    "exact_field"
                }
            )
    }, logical(1))
    keys <- vapply(authority$owners, `[[`, character(1), "access_key")
    valid <- identical(
        authority$schema,
        "multischolar.omics_payload_access_owners"
    ) && identical(authority$schema_version, "1.0.0") &&
        identical(authority$owner_ticket_id, "OMICS-ART-069") &&
        grepl("^[0-9a-f]{64}$", authority$discovery_digest) &&
        all(valid_owner) && !anyDuplicated(keys) &&
        identical(authority$unknown_policy, "reject") &&
        is.logical(authority$candidate_freeze_authority) &&
        length(authority$candidate_freeze_authority) == 1L &&
        !is.na(authority$candidate_freeze_authority)
    if (!valid) stop("Payload-access owner authority is invalid")
    invisible(authority)
}

payloadAccessReconcile <- function(accesses, authority) {
    production <- Filter(function(value) {
        identical(value$evidence_class, "production")
    }, accesses)
    production_keys <- vapply(production, `[[`, character(1), "access_key")
    if (is.null(authority)) {
        return(list(
            authority_supplied = FALSE,
            unowned_keys = as.list(production_keys),
            stale_owner_keys = list(),
            duplicate_owner_keys = list(),
            mismatched_owner_keys = list(),
            unresolved_computed_keys = as.list(vapply(
                Filter(function(value) isTRUE(value$computed), production),
                `[[`,
                character(1),
                "access_key"
            ))
        ))
    }
    payloadAccessValidateAuthority(authority)
    owners <- authority$owners
    owner_keys <- vapply(owners, `[[`, character(1), "access_key")
    duplicates <- unique(owner_keys[duplicated(owner_keys)])
    access_map <- stats::setNames(production, production_keys)
    owner_map <- stats::setNames(owners, owner_keys)
    shared <- intersect(production_keys, owner_keys)
    mismatched <- shared[!vapply(shared, function(key) {
        access <- access_map[[key]]
        owner <- owner_map[[key]]
        identical(owner$path, access$path) &&
            identical(owner$expression, access$expression) &&
            identical(owner$access_kind, access$access_kind) &&
            identical(isTRUE(owner$computed), isTRUE(access$computed))
    }, logical(1))]
    unresolved <- shared[vapply(shared, function(key) {
        access <- access_map[[key]]
        owner <- owner_map[[key]]
        isTRUE(access$computed) && !identical(
            owner$computed_resolution,
            "contract_source_fields_guarded"
        )
    }, logical(1))]
    list(
        authority_supplied = TRUE,
        unowned_keys = as.list(setdiff(production_keys, owner_keys)),
        stale_owner_keys = as.list(setdiff(owner_keys, production_keys)),
        duplicate_owner_keys = as.list(duplicates),
        mismatched_owner_keys = as.list(mismatched),
        unresolved_computed_keys = as.list(unresolved)
    )
}

payloadAccessAudit <- function(
    repo_root,
    scan_roots = c("R"),
    authority_path = NULL
) {
    source_fields <- payloadAccessSourceFields()
    files <- unlist(lapply(scan_roots, function(root) {
        list.files(
            file.path(repo_root, root),
            pattern = "[.][Rr]$",
            recursive = TRUE,
            full.names = TRUE
        )
    }), use.names = FALSE)
    files <- sort(unique(files), method = "radix")
    scanned <- lapply(files, payloadAccessScanFile, repo_root, source_fields)
    accesses <- payloadAccessMerge(unlist(
        lapply(scanned, `[[`, "accesses"),
        recursive = FALSE
    ))
    errors <- unlist(lapply(scanned, `[[`, "errors"), recursive = FALSE)
    authority <- payloadAccessReadAuthority(authority_path)
    reconciliation <- payloadAccessReconcile(accesses, authority)
    production <- Filter(function(value) {
        identical(value$evidence_class, "production")
    }, accesses)
    unresolved <- reconciliation$unresolved_computed_keys
    status <- if (length(errors)) {
        "parse_failed"
    } else if (length(reconciliation$unowned_keys) ||
        length(reconciliation$stale_owner_keys) ||
        length(reconciliation$duplicate_owner_keys) ||
        length(reconciliation$mismatched_owner_keys) || length(unresolved)) {
        "ownership_incomplete"
    } else if (!is.null(authority) &&
        !isTRUE(authority$candidate_freeze_authority)) {
        "ownership_assigned_pending_adoption"
    } else {
        "passed"
    }
    record <- list(
        schema = "multischolar.omics_payload_access_audit",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-069",
        status = status,
        scan_roots = as.list(scan_roots),
        source_fields = as.list(source_fields),
        file_count = as.integer(length(files)),
        access_count = as.integer(length(accesses)),
        production_access_count = as.integer(length(production)),
        unresolved_computed_count = as.integer(length(unresolved)),
        accesses = accesses,
        parse_errors = errors,
        reconciliation = reconciliation,
        candidate_freeze_allowed = identical(status, "passed"),
        publication_authority = FALSE
    )
    basis <- record
    basis$audit_digest <- NULL
    record$audit_digest <- payloadAccessObjectDigest(basis)
    record
}

payloadAccessWriteJson <- function(value, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    temporary <- tempfile("payload-access-", tmpdir = dirname(path))
    on.exit(unlink(temporary, force = TRUE), add = TRUE)
    jsonlite::write_json(
        value,
        temporary,
        auto_unbox = TRUE,
        pretty = TRUE,
        null = "null",
        na = "null",
        digits = 17
    )
    if (!file.rename(temporary, path)) stop("Could not publish access audit")
    invisible(path)
}

payloadAccessArgs <- function(argv) {
    values <- list(
        repo_root = payloadAccessRepoRoot(),
        scan_roots = "R",
        authority = NULL,
        output = NULL,
        fail_incomplete = "false"
    )
    index <- 1L
    while (index <= length(argv)) {
        token <- argv[[index]]
        key <- gsub("-", "_", sub("^--", "", token), fixed = TRUE)
        if (!startsWith(token, "--") || !key %in% names(values) ||
            index == length(argv)) {
            stop("Payload-access audit arguments are invalid")
        }
        values[[key]] <- argv[[index + 1L]]
        index <- index + 2L
    }
    if (is.null(values$output)) stop("--output is required")
    values$repo_root <- normalizePath(values$repo_root, mustWork = TRUE)
    values$scan_roots <- strsplit(
        values$scan_roots,
        ",",
        fixed = TRUE
    )[[1L]]
    values$fail_incomplete <- identical(
        tolower(values$fail_incomplete),
        "true"
    )
    values
}

payloadAccessMain <- function(argv = commandArgs(trailingOnly = TRUE)) {
    args <- payloadAccessArgs(argv)
    record <- payloadAccessAudit(
        args$repo_root,
        args$scan_roots,
        args$authority
    )
    payloadAccessWriteJson(record, args$output)
    cat(record$audit_digest, "\n")
    if (args$fail_incomplete && !identical(record$status, "passed")) {
        return(invisible(2L))
    }
    invisible(0L)
}

if (identical(environment(), globalenv()) && !interactive()) {
    status <- tryCatch(
        payloadAccessMain(),
        error = function(error) {
            message(conditionMessage(error))
            1L
        }
    )
    quit(status = as.integer(status), save = "no")
}
