publicationCleanupAbort <- function(message) {
    publicationBuildAbort(paste("Comparator cleanup", message))
}

publicationCleanupRelativePath <- function(path, root) {
    normalized_root <- normalizePath(root, mustWork = TRUE)
    normalized_path <- normalizePath(path, mustWork = FALSE)
    prefix <- paste0(normalized_root, .Platform$file.sep)
    if (!startsWith(paste0(normalized_path, .Platform$file.sep), prefix) ||
        identical(normalized_path, normalized_root)) {
        publicationCleanupAbort("path is outside the owned root")
    }
    substring(normalized_path, nchar(prefix) + 1L)
}

publicationCleanupLexicalRelativePath <- function(path, root) {
    normalized_root <- normalizePath(root, mustWork = TRUE)
    prefix <- paste0(normalized_root, .Platform$file.sep)
    if (!startsWith(path, prefix)) {
        publicationCleanupAbort("inventory path is outside the owned root")
    }
    relative <- substring(path, nchar(prefix) + 1L)
    components <- strsplit(relative, .Platform$file.sep, fixed = TRUE)[[1L]]
    if (!nzchar(relative) || any(components %in% c("", ".", ".."))) {
        publicationCleanupAbort("inventory path is not canonical")
    }
    relative
}

publicationCleanupWalk <- function(root) {
    walk <- function(path) {
        entries <- list.files(
            path,
            all.files = TRUE,
            no.. = TRUE,
            full.names = TRUE
        )
        unlist(lapply(entries, \(entry) {
            link <- Sys.readlink(entry)
            if (nzchar(link) || !isTRUE(file.info(entry)$isdir)) return(entry)
            c(entry, walk(entry))
        }), use.names = FALSE)
    }
    walk(root)
}

publicationCleanupInventory <- function(root) {
    paths <- publicationCleanupWalk(root)
    records <- lapply(paths, \(path) {
        link <- Sys.readlink(path)
        info <- file.info(path)
        type <- if (nzchar(link)) {
            "symlink"
        } else if (isTRUE(info$isdir)) {
            "directory"
        } else {
            "file"
        }
        list(
            path = publicationCleanupLexicalRelativePath(path, root),
            type = type,
            size_bytes = if (identical(type, "file")) {
                as.numeric(info$size)
            } else {
                0
            },
            symlink_target = if (identical(type, "symlink")) link else NULL
        )
    })
    list(
        entry_count = as.integer(length(records)),
        file_count = as.integer(sum(vapply(
            records,
            \(record) identical(record$type, "file"),
            logical(1)
        ))),
        symlink_count = as.integer(sum(vapply(
            records,
            \(record) identical(record$type, "symlink"),
            logical(1)
        ))),
        file_bytes = sum(vapply(records, `[[`, numeric(1), "size_bytes")),
        metadata_sha256 = publicationObjectDigest(records),
        records = records
    )
}

publicationCleanupRetainPath <- function(path) {
    startsWith(path, "logs/") || startsWith(path, "scripts/") ||
        identical(path, "project/renv.lock") ||
        (!grepl("/", path, fixed = TRUE) && grepl(
            "\\.(json|R|tar\\.gz)$",
            path
        ))
}

publicationCleanupRetainedFiles <- function(root, inventory) {
    records <- Filter(\(record) {
        identical(record$type, "file") &&
            publicationCleanupRetainPath(record$path)
    }, inventory$records)
    lapply(records, \(record) {
        path <- file.path(root, record$path)
        list(
            path = record$path,
            size_bytes = record$size_bytes,
            sha256 = publicationFileDigest(path)
        )
    })
}

publicationCleanupCopyEvidence <- function(root, archive, retained) {
    if (file.exists(archive) || dir.exists(archive)) {
        publicationCleanupAbort("attempt evidence archive already exists")
    }
    dir.create(archive, recursive = TRUE, showWarnings = FALSE)
    archived <- lapply(retained, \(binding) {
        source <- file.path(root, binding$path)
        destination <- file.path(archive, "retained", binding$path)
        dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
        if (!isTRUE(file.copy(source, destination, copy.mode = TRUE))) {
            publicationCleanupAbort("evidence copy failed")
        }
        if (!identical(publicationFileDigest(destination), binding$sha256)) {
            publicationCleanupAbort("archived evidence bytes differ")
        }
        c(binding, list(archive_path = file.path("retained", binding$path)))
    })
    archived
}

publicationCleanupValidatePlan <- function(
    plan,
    comparator_root,
    attempts_must_exist = TRUE
) {
    publicationRequireNames(plan, c(
        "schema", "schema_version", "owner_ticket_id", "status",
        "comparator_root", "evidence_root", "protected_paths", "attempts",
        "retention_policy", "publication_authority"
    ), "Comparator cleanup plan")
    attempt_paths <- vapply(plan$attempts, `[[`, character(1), "path")
    valid <- identical(
        plan$schema,
        "multischolar.omics_publication_comparator_cleanup_plan"
    ) && identical(plan$schema_version, "1.0.0") &&
        identical(plan$owner_ticket_id, "OMICS-ART-064") &&
        identical(plan$status, "approved_for_execution") &&
        identical(plan$comparator_root, ".omics-publication-comparators") &&
        length(attempt_paths) > 0L && !anyDuplicated(attempt_paths) &&
        all(vapply(plan$attempts, \(attempt) {
            setequal(names(attempt), c("path", "category", "reason")) &&
                publicationScalarString(attempt$path) &&
                publicationScalarString(attempt$category) &&
                publicationScalarString(attempt$reason)
        }, logical(1))) && !isTRUE(plan$publication_authority)
    if (!valid) publicationCleanupAbort("plan differs")
    protected <- vapply(plan$protected_paths, \(path) {
        publicationCleanupRelativePath(file.path(comparator_root, path), comparator_root)
    }, character(1))
    resolved <- vapply(attempt_paths, \(path) {
        full <- file.path(comparator_root, path)
        if (isTRUE(attempts_must_exist) && !dir.exists(full)) {
            publicationCleanupAbort("attempt root is missing")
        }
        publicationCleanupRelativePath(full, comparator_root)
    }, character(1))
    overlap <- vapply(resolved, \(path) {
        any(path == protected) || any(startsWith(path, paste0(protected, "/"))) ||
            any(startsWith(protected, paste0(path, "/")))
    }, logical(1))
    if (any(overlap)) publicationCleanupAbort("plan includes a protected root")
    invisible(plan)
}

publicationValidateCleanupRetained <- function(binding, archive) {
    publicationRequireNames(binding, c(
        "path", "size_bytes", "sha256", "archive_path"
    ), "Retained cleanup evidence")
    path <- file.path(archive, binding$archive_path)
    publicationEvidenceFileCurrent(
        path,
        binding$sha256,
        "retained cleanup evidence",
        archive,
        binding$size_bytes
    )
}

publicationValidateComparatorCleanupResult <- function(record, plan) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "owner_ticket_id", "status", "plan_sha256",
        "attempts", "attempt_count", "removed_count", "removed_file_bytes",
        "publication_authority"
    ), "Comparator cleanup result")
    comparator_root <- normalizePath(plan$comparator_root, mustWork = TRUE)
    publicationCleanupValidatePlan(
        plan,
        comparator_root,
        attempts_must_exist = FALSE
    )
    evidence_root <- file.path(comparator_root, plan$evidence_root)
    attempt_valid <- vapply(seq_along(record$attempts), \(index) {
        attempt <- record$attempts[[index]]
        plan_attempt <- plan$attempts[[index]]
        archive <- file.path(
            evidence_root,
            publicationCleanupArchiveName(plan_attempt)
        )
        archived_record <- publicationReadJson(file.path(
            archive,
            "attempt-record.json"
        ))
        lapply(
            attempt$retained,
            publicationValidateCleanupRetained,
            archive = archive
        )
        identical(
            publicationObjectDigest(attempt),
            publicationObjectDigest(archived_record)
        ) && identical(attempt$path, plan_attempt$path) &&
            identical(attempt$category, plan_attempt$category) &&
            identical(attempt$reason, plan_attempt$reason) &&
            isTRUE(attempt$removed) &&
            !dir.exists(file.path(comparator_root, attempt$path))
    }, logical(1))
    protected <- all(vapply(plan$protected_paths, \(path) {
        dir.exists(file.path(comparator_root, path))
    }, logical(1)))
    removed_bytes <- sum(vapply(record$attempts, \(attempt) {
        attempt$before$file_bytes
    }, numeric(1)))
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_comparator_cleanup_result"
    ) && identical(record$schema_version, "1.0.0") &&
        identical(record$owner_ticket_id, "OMICS-ART-064") &&
        identical(record$status, "passed") &&
        identical(record$plan_sha256, publicationObjectDigest(plan)) &&
        identical(record$attempt_count, as.integer(length(plan$attempts))) &&
        identical(record$removed_count, record$attempt_count) &&
        identical(record$removed_file_bytes, removed_bytes) &&
        all(attempt_valid) && protected && !isTRUE(record$publication_authority)
    if (!valid) publicationCleanupAbort("result differs")
    invisible(record)
}

publicationCleanupArchiveName <- function(attempt) {
    paste(
        gsub("/", "__", attempt$path, fixed = TRUE),
        gsub("[^A-Za-z0-9_.-]", "_", attempt$category),
        sep = "--"
    )
}

publicationCleanupAttempt <- function(
    attempt,
    comparator_root,
    evidence_root,
    execute = FALSE
) {
    root <- file.path(comparator_root, attempt$path)
    inventory <- publicationCleanupInventory(root)
    retained <- publicationCleanupRetainedFiles(root, inventory)
    record <- list(
        schema = "multischolar.omics_publication_cleanup_attempt",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-064",
        path = attempt$path,
        category = attempt$category,
        reason = attempt$reason,
        before = inventory[setdiff(names(inventory), "records")],
        retained = retained,
        execute_requested = isTRUE(execute),
        removed = FALSE,
        publication_authority = FALSE
    )
    if (!isTRUE(execute)) return(record)
    archive <- file.path(evidence_root, publicationCleanupArchiveName(attempt))
    record$retained <- publicationCleanupCopyEvidence(root, archive, retained)
    publicationWriteJson(record, file.path(archive, "attempt-record.json"))
    unlink(root, recursive = TRUE, force = TRUE)
    if (file.exists(root) || dir.exists(root)) {
        publicationCleanupAbort("owned attempt root was not removed")
    }
    record$removed <- TRUE
    publicationWriteJson(record, file.path(archive, "attempt-record.json"))
    record
}

publicationRunComparatorCleanup <- function(plan, execute = FALSE) {
    comparator_root <- normalizePath(plan$comparator_root, mustWork = TRUE)
    publicationCleanupValidatePlan(plan, comparator_root)
    evidence_root <- file.path(comparator_root, plan$evidence_root)
    if (isTRUE(execute) && dir.exists(evidence_root)) {
        publicationCleanupAbort("evidence root already exists")
    }
    attempts <- lapply(plan$attempts, publicationCleanupAttempt,
        comparator_root = comparator_root,
        evidence_root = evidence_root,
        execute = execute
    )
    removed <- vapply(attempts, `[[`, logical(1), "removed")
    list(
        schema = "multischolar.omics_publication_comparator_cleanup_result",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-064",
        status = if (isTRUE(execute) && all(removed)) "passed" else "dry_run",
        plan_sha256 = publicationObjectDigest(plan),
        attempts = attempts,
        attempt_count = as.integer(length(attempts)),
        removed_count = as.integer(sum(removed)),
        removed_file_bytes = sum(vapply(attempts, \(attempt) {
            if (isTRUE(attempt$removed)) attempt$before$file_bytes else 0
        }, numeric(1))),
        publication_authority = FALSE
    )
}
