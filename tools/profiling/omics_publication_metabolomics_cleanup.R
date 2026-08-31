metabCleanupAbort <- function(message) {
    metabPublicationAbort(paste("cleanup", message), "cleanup_error")
}

metabCleanupRoot <- function() {
    publicationPath(".omics-publication-workloads", "metabolomics")
}

metabCleanupLexicalRelative <- function(path, root) {
    prefix <- paste0(normalizePath(root), .Platform$file.sep)
    if (!startsWith(path, prefix)) metabCleanupAbort("path is not owned")
    relative <- substring(path, nchar(prefix) + 1L)
    parts <- strsplit(relative, .Platform$file.sep, fixed = TRUE)[[1L]]
    if (!nzchar(relative) || any(parts %in% c("", ".", ".."))) {
        metabCleanupAbort("path is not canonical")
    }
    relative
}

metabCleanupWalk <- function(root) {
    walk <- function(path) {
        entries <- list.files(
            path,
            all.files = TRUE,
            no.. = TRUE,
            full.names = TRUE
        )
        unlist(lapply(entries, \(entry) {
            if (nzchar(Sys.readlink(entry)) || !isTRUE(file.info(entry)$isdir)) {
                return(entry)
            }
            c(entry, walk(entry))
        }), use.names = FALSE)
    }
    walk(root)
}

metabCleanupInventory <- function(root) {
    paths <- metabCleanupWalk(root)
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
            path = metabCleanupLexicalRelative(path, root),
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
        file_count = as.integer(sum(vapply(records, \(record) {
            identical(record$type, "file")
        }, logical(1)))),
        symlink_count = as.integer(sum(vapply(records, \(record) {
            identical(record$type, "symlink")
        }, logical(1)))),
        file_bytes = sum(vapply(records, `[[`, numeric(1), "size_bytes")),
        metadata_sha256 = publicationObjectDigest(records),
        records = records
    )
}

metabCleanupRetain <- function(record) {
    identical(record$type, "file") && (
        startsWith(record$path, "logs/") ||
            startsWith(record$path, "summaries/") ||
            grepl("(receipt|result|summary|measurement)[.]json$", record$path)
    )
}

metabCleanupCopyRetained <- function(root, archive, inventory) {
    records <- Filter(metabCleanupRetain, inventory$records)
    lapply(records, \(record) {
        source <- file.path(root, record$path)
        target <- file.path(archive, "retained", record$path)
        dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
        if (!isTRUE(file.copy(source, target, copy.mode = TRUE))) {
            metabCleanupAbort("evidence copy failed")
        }
        source_sha256 <- publicationFileDigest(source)
        if (!identical(publicationFileDigest(target), source_sha256)) {
            metabCleanupAbort("retained evidence differs")
        }
        list(
            path = record$path,
            archive_path = file.path("retained", record$path),
            size_bytes = record$size_bytes,
            sha256 = source_sha256
        )
    })
}

metabCleanupValidatePlan <- function(plan, root, require_removals = TRUE) {
    publicationRequireNames(plan, c(
        "schema", "schema_version", "owner_ticket_id", "status", "root",
        "archive_root", "protected_paths", "removals", "publication_authority"
    ), "Metabolomics cleanup plan")
    protected <- unlist(plan$protected_paths, use.names = FALSE)
    removals <- vapply(plan$removals, `[[`, character(1), "path")
    valid <- identical(
        plan$schema,
        "multischolar.omics_publication_metabolomics_cleanup_plan"
    ) && identical(plan$schema_version, "1.0.0") &&
        identical(plan$owner_ticket_id, "OMICS-ART-066") &&
        identical(plan$status, "approved_for_execution") &&
        identical(plan$root, ".omics-publication-workloads/metabolomics") &&
        length(removals) > 0L && !anyDuplicated(removals) &&
        !isTRUE(plan$publication_authority)
    for (relative in c(protected, removals)) {
        full <- file.path(root, relative)
        if (isTRUE(require_removals) && relative %in% removals &&
            !dir.exists(full)) {
            metabCleanupAbort("removal path is missing")
        }
        normalized <- normalizePath(full, mustWork = FALSE)
        if (!startsWith(
            paste0(normalized, .Platform$file.sep),
            paste0(normalizePath(root), .Platform$file.sep)
        )) {
            metabCleanupAbort("plan path is not owned")
        }
    }
    overlap <- vapply(removals, \(removal) {
        any(removal == protected) ||
            any(startsWith(removal, paste0(protected, "/"))) ||
            any(startsWith(protected, paste0(removal, "/")))
    }, logical(1))
    if (!valid || any(overlap)) metabCleanupAbort("plan differs")
    invisible(plan)
}

metabCleanupOne <- function(removal, root, archive_root, execute) {
    publicationRequireNames(
        removal,
        c("path", "reason"),
        "Metabolomics cleanup removal"
    )
    source <- file.path(root, removal$path)
    inventory <- metabCleanupInventory(source)
    record <- list(
        path = removal$path,
        reason = removal$reason,
        before = inventory[setdiff(names(inventory), "records")],
        retained = list(),
        removed = FALSE
    )
    if (!isTRUE(execute)) return(record)
    name <- gsub("/", "__", removal$path, fixed = TRUE)
    archive <- file.path(archive_root, name)
    if (file.exists(archive) || dir.exists(archive)) {
        metabCleanupAbort("archive already exists")
    }
    dir.create(archive, recursive = TRUE, showWarnings = FALSE)
    record$retained <- metabCleanupCopyRetained(source, archive, inventory)
    publicationWriteJson(record, file.path(archive, "cleanup-record.json"))
    unlink(source, recursive = TRUE, force = TRUE)
    if (dir.exists(source)) metabCleanupAbort("removal failed")
    record$removed <- TRUE
    publicationWriteJson(record, file.path(archive, "cleanup-record.json"))
    record
}

metabPublicationRunCleanup <- function(plan, execute = FALSE) {
    root <- normalizePath(plan$root, mustWork = TRUE)
    metabCleanupValidatePlan(plan, root)
    archive_root <- file.path(root, plan$archive_root)
    if (isTRUE(execute) && dir.exists(archive_root)) {
        metabCleanupAbort("archive root already exists")
    }
    records <- lapply(plan$removals, metabCleanupOne,
        root = root,
        archive_root = archive_root,
        execute = execute
    )
    removed <- vapply(records, `[[`, logical(1), "removed")
    list(
        schema = "multischolar.omics_publication_metabolomics_cleanup_result",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-066",
        status = if (isTRUE(execute) && all(removed)) "passed" else "dry_run",
        plan_sha256 = publicationObjectDigest(plan),
        removals = records,
        removed_count = as.integer(sum(removed)),
        removed_file_bytes = sum(vapply(records, \(record) {
            if (isTRUE(record$removed)) record$before$file_bytes else 0
        }, numeric(1))),
        publication_authority = FALSE
    )
}
