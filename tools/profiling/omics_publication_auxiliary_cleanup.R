auxPublicationCleanupInventory <- function(path) {
    resolved <- normalizePath(path, mustWork = TRUE)
    files <- sort(list.files(
        resolved,
        recursive = TRUE,
        full.names = TRUE,
        all.files = TRUE,
        no.. = TRUE
    ), method = "radix")
    files <- files[file.info(files)$isdir %in% FALSE]
    records <- lapply(files, function(file) {
        list(
            relative_path = substring(file, nchar(resolved) + 2L),
            sha256 = publicationFileDigest(file),
            size_bytes = as.numeric(file.info(file)$size)
        )
    })
    list(
        file_count = as.integer(length(records)),
        size_bytes = sum(vapply(records, `[[`, numeric(1), "size_bytes")),
        inventory_sha256 = publicationObjectDigest(records),
        files = records
    )
}

auxPublicationCleanupTargets <- function() {
    root <- ".omics-publication-workloads/auxiliary/execution"
    list(
        file.path(root, "historical", "fixture-campaign"),
        file.path(root, "historical", "fixture-mofa-campaign"),
        file.path(root, "pre-repair", "fixture-campaign")
    )
}

auxPublicationCleanupProtected <- function() {
    root <- ".omics-publication-workloads/auxiliary"
    list(
        file.path(root, "freeze-a"),
        file.path(root, "execution", "historical", "fixture-campaign-v2"),
        file.path(root, "execution", "pre-repair", "fixture-campaign-v2"),
        file.path(
            root,
            "execution",
            "historical",
            "representative-multiomics-v1"
        ),
        file.path(
            root,
            "execution",
            "pre-repair",
            "representative-multiomics-v1"
        ),
        file.path(
            root,
            "execution",
            "historical",
            "representative-phosphosite-v1"
        )
    )
}

auxPublicationBuildCleanupPlan <- function() {
    targets <- Filter(function(path) {
        file.exists(path) || dir.exists(path)
    }, auxPublicationCleanupTargets())
    protected <- auxPublicationCleanupProtected()
    overlap <- intersect(
        vapply(targets, normalizePath, character(1), mustWork = TRUE),
        vapply(protected, normalizePath, character(1), mustWork = TRUE)
    )
    if (length(overlap)) {
        auxPublicationAbort("auxiliary cleanup target is protected")
    }
    records <- lapply(targets, function(path) {
        inventory <- auxPublicationCleanupInventory(path)
        list(
            source_path = path,
            archive_name = gsub("[/]", "__", path),
            file_count = inventory$file_count,
            size_bytes = inventory$size_bytes,
            inventory_sha256 = inventory$inventory_sha256,
            action = "archive",
            deletion_allowed = FALSE
        )
    })
    list(
        schema = "multischolar.omics_publication_auxiliary_cleanup_plan",
        schema_version = .AUX_PUBLICATION_VERSION,
        owner_ticket_id = .AUX_PUBLICATION_OWNER,
        status = "planned",
        archive_root = paste0(
            ".omics-publication-workloads/auxiliary/cleanup-evidence/",
            "omics-art-068-intermediate"
        ),
        targets = records,
        protected_paths = protected,
        deletion_allowed = FALSE,
        publication_authority = FALSE
    )
}

auxPublicationValidateCleanupPlan <- function(plan) {
    publicationRequireNames(plan, c(
        "schema", "schema_version", "owner_ticket_id", "status",
        "archive_root", "targets", "protected_paths", "deletion_allowed",
        "publication_authority"
    ), "Auxiliary cleanup plan")
    target_valid <- all(vapply(plan$targets, function(target) {
        publicationRequireNames(target, c(
            "source_path", "archive_name", "file_count", "size_bytes",
            "inventory_sha256", "action", "deletion_allowed"
        ), "Auxiliary cleanup target")
        current <- auxPublicationCleanupInventory(target$source_path)
        identical(target$action, "archive") &&
            identical(target$file_count, current$file_count) &&
            identical(
                as.numeric(target$size_bytes),
                as.numeric(current$size_bytes)
            ) &&
            identical(target$inventory_sha256, current$inventory_sha256) &&
            !isTRUE(target$deletion_allowed)
    }, logical(1)))
    sources <- vapply(plan$targets, `[[`, character(1), "source_path")
    valid <- identical(
        plan$schema,
        "multischolar.omics_publication_auxiliary_cleanup_plan"
    ) && identical(plan$schema_version, .AUX_PUBLICATION_VERSION) &&
        identical(plan$owner_ticket_id, .AUX_PUBLICATION_OWNER) &&
        identical(plan$status, "planned") && !anyDuplicated(sources) &&
        target_valid && !isTRUE(plan$deletion_allowed) &&
        !isTRUE(plan$publication_authority)
    if (!valid) auxPublicationAbort("auxiliary cleanup plan differs")
    invisible(plan)
}

auxPublicationApplyCleanup <- function(plan, dry_run = TRUE) {
    auxPublicationValidateCleanupPlan(plan)
    if (!dry_run) {
        dir.create(plan$archive_root, recursive = TRUE, showWarnings = FALSE)
    }
    records <- lapply(plan$targets, function(target) {
        destination <- file.path(plan$archive_root, target$archive_name)
        if (file.exists(destination) || dir.exists(destination)) {
            auxPublicationAbort("auxiliary cleanup archive already exists")
        }
        if (!dry_run && !file.rename(target$source_path, destination)) {
            auxPublicationAbort("auxiliary cleanup archive move failed")
        }
        list(
            source_path = target$source_path,
            destination_path = destination,
            action = if (dry_run) "would_archive" else "archived",
            file_count = target$file_count,
            size_bytes = target$size_bytes,
            inventory_sha256 = target$inventory_sha256,
            deletion_performed = FALSE
        )
    })
    list(
        schema = "multischolar.omics_publication_auxiliary_cleanup_result",
        schema_version = .AUX_PUBLICATION_VERSION,
        owner_ticket_id = .AUX_PUBLICATION_OWNER,
        status = if (dry_run) "dry_run_passed" else "passed",
        records = records,
        archived_file_count = as.integer(sum(vapply(
            plan$targets,
            `[[`,
            numeric(1),
            "file_count"
        ))),
        archived_size_bytes = sum(vapply(
            plan$targets,
            `[[`,
            numeric(1),
            "size_bytes"
        )),
        deletion_performed = FALSE,
        publication_authority = FALSE
    )
}
