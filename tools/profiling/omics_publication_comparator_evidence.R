publicationEvidenceAbort <- function(message) {
    publicationBuildAbort(paste("Comparator evidence", message))
}

publicationEvidenceOwnedPath <- function(path, root, label) {
    if (!publicationScalarString(path) || !file.exists(path)) {
        publicationEvidenceAbort(paste(label, "path is missing"))
    }
    normalized_root <- normalizePath(root, mustWork = TRUE)
    normalized_path <- normalizePath(path, mustWork = TRUE)
    prefix <- paste0(normalized_root, .Platform$file.sep)
    if (!startsWith(paste0(normalized_path, .Platform$file.sep), prefix)) {
        publicationEvidenceAbort(paste(label, "path is not owned"))
    }
    normalized_path
}

publicationEvidenceFileCurrent <- function(
    path,
    sha256,
    label,
    root = NULL,
    size_bytes = NULL
) {
    if (!publicationScalarString(path) || !file.exists(path) ||
        !publicationComparatorShaValid(sha256)) {
        publicationEvidenceAbort(paste(label, "binding is incomplete"))
    }
    if (!is.null(root)) publicationEvidenceOwnedPath(path, root, label)
    size_valid <- is.null(size_bytes) || identical(
        as.numeric(file.info(path)$size),
        as.numeric(size_bytes)
    )
    if (!size_valid || !identical(publicationFileDigest(path), sha256)) {
        publicationEvidenceAbort(paste(label, "bytes differ"))
    }
    invisible(path)
}

publicationValidateCommandEvidence <- function(command, root, label) {
    publicationRequireNames(command, c(
        "command", "arguments", "workdir", "environment", "exit_status",
        "timed_out", "elapsed_seconds", "stdout_sha256", "stderr_sha256",
        "stdout_path", "stderr_path"
    ), paste(label, "command"))
    valid <- file.exists(command$command) &&
        identical(normalizePath(command$workdir), normalizePath(root)) &&
        identical(command$exit_status, 0L) && !isTRUE(command$timed_out) &&
        publicationScalarNumber(command$elapsed_seconds) &&
        identical(command$environment$GITHUB_PAT, "") &&
        identical(command$environment$GITHUB_TOKEN, "") &&
        identical(command$environment$OMP_NUM_THREADS, "1") &&
        identical(command$environment$OPENBLAS_NUM_THREADS, "1") &&
        identical(command$environment$MAKEFLAGS, "-j1")
    if (!valid) publicationEvidenceAbort(paste(label, "command differs"))
    publicationEvidenceFileCurrent(
        command$stdout_path,
        command$stdout_sha256,
        paste(label, "stdout"),
        root
    )
    publicationEvidenceFileCurrent(
        command$stderr_path,
        command$stderr_sha256,
        paste(label, "stderr"),
        root
    )
    invisible(command)
}

publicationValidateScriptEvidence <- function(script, root, label) {
    publicationRequireNames(
        script,
        c("path", "sha256", "size_bytes"),
        paste(label, "script")
    )
    publicationEvidenceFileCurrent(
        script$path,
        script$sha256,
        paste(label, "script"),
        root,
        script$size_bytes
    )
}

publicationComparatorBuildReceiptRoot <- function(path) {
    root <- dirname(normalizePath(path, mustWork = TRUE))
    publicationEvidenceOwnedPath(
        root,
        file.path(publicationComparatorRoot(), "comparator-builds"),
        "build root"
    )
    root
}

publicationValidateComparatorSourceEvidence <- function(receipt, root) {
    source <- receipt$source
    publicationRequireNames(source, c(
        "path", "archive_path", "archive_sha256", "description_sha256",
        "source_identity"
    ), "Comparator source evidence")
    publicationEvidenceOwnedPath(source$path, root, "source tree")
    publicationEvidenceFileCurrent(
        source$archive_path,
        source$archive_sha256,
        "source archive",
        root
    )
    description <- file.path(source$path, "DESCRIPTION")
    publicationEvidenceFileCurrent(
        description,
        source$description_sha256,
        "source DESCRIPTION",
        root
    )
    current <- publicationComparatorSourceIdentity(receipt$revision)
    valid <- identical(source$source_identity, current) && identical(
        source$archive_sha256,
        current$archive_sha256
    )
    if (!valid) publicationEvidenceAbort("source identity differs")
    invisible(source)
}

publicationComparatorOverlayEvidence <- function(
    package_library,
    dependency_library,
    lock,
    restore_receipt
) {
    packages <- names(lock$Packages)
    entries <- list.files(package_library, all.files = TRUE, no.. = TRUE)
    if (!setequal(entries, c(packages, "MultiScholaR"))) {
        publicationEvidenceAbort("dependency overlay entries differ")
    }
    records <- lapply(packages, \(package) {
        target <- file.path(package_library, package)
        expected <- normalizePath(file.path(dependency_library, package))
        if (!nzchar(Sys.readlink(target)) ||
            !identical(normalizePath(target), expected)) {
            publicationEvidenceAbort("dependency overlay target differs")
        }
        list(
            package = package,
            version = lock$Packages[[package]]$Version,
            target_sha256 = restore_receipt$library_inventory$packages[[
                package
            ]]$tree$inventory_sha256
        )
    })
    names(records) <- packages
    list(
        package_count = as.integer(length(records)),
        records_sha256 = publicationObjectDigest(records)
    )
}

publicationValidateComparatorEnvironmentEvidence <- function(receipt, root) {
    environment <- receipt$environment
    publicationRequireNames(environment, c(
        "lockfile_sha256", "restore_receipt", "dependency_overlay",
        "dependency_before", "dependency_after"
    ), "Comparator environment evidence")
    lock <- publicationReadJson("renv.lock")
    restore_binding <- publicationRestoreReceiptBinding(
        environment$restore_receipt$path
    )
    restore_receipt <- publicationReadJson(restore_binding$path)
    dependency_library <- file.path(dirname(restore_binding$path), "library")
    package_library <- dirname(receipt$build$installed_inventory$root)
    current <- publicationDependencyFingerprint(dependency_library, lock)
    overlay <- publicationComparatorOverlayEvidence(
        package_library,
        dependency_library,
        lock,
        restore_receipt
    )
    valid <- identical(environment$lockfile_sha256, publicationFileDigest(
        "renv.lock"
    )) && identical(environment$restore_receipt, restore_binding) &&
        identical(environment$dependency_before, current) &&
        identical(environment$dependency_after, current) &&
        identical(environment$dependency_overlay, overlay)
    if (!valid) publicationEvidenceAbort("dependency environment differs")
    publicationEvidenceOwnedPath(package_library, root, "package library")
    invisible(environment)
}

publicationValidateComparatorBuildFiles <- function(receipt, root) {
    build <- receipt$build
    publicationRequireNames(build, c(
        "build_command", "install_command", "archive", "installed_inventory",
        "native_libraries"
    ), "Comparator build evidence")
    publicationValidateCommandEvidence(build$build_command, root, "build")
    publicationValidateCommandEvidence(build$install_command, root, "install")
    publicationEvidenceFileCurrent(
        build$archive$path,
        build$archive$sha256,
        "package archive",
        root,
        build$archive$size_bytes
    )
    installed <- build$installed_inventory$root
    publicationEvidenceOwnedPath(installed, root, "installed package")
    current_inventory <- publicationInstalledInventory(installed)
    current_native <- publicationNativeLibraryInventory(installed)
    if (!identical(
        publicationObjectDigest(build$installed_inventory),
        publicationObjectDigest(current_inventory)
    ) || !identical(
        publicationObjectDigest(build$native_libraries),
        publicationObjectDigest(current_native)
    )) {
        publicationEvidenceAbort("installed package bytes differ")
    }
    invisible(build)
}

publicationValidateComparatorSmokeEvidence <- function(receipt, root) {
    smoke <- receipt$smoke
    publicationRequireNames(smoke, c(
        "status", "script", "command", "output_sha256", "evidence"
    ), "Comparator smoke evidence")
    publicationValidateScriptEvidence(smoke$script, root, "smoke")
    publicationValidateCommandEvidence(smoke$command, root, "smoke")
    output <- file.path(root, "smoke-output.json")
    publicationEvidenceFileCurrent(
        output,
        smoke$output_sha256,
        "smoke output",
        root
    )
    current <- publicationReadJson(output)
    valid <- identical(smoke$status, "passed") &&
        identical(current, smoke$evidence) &&
        identical(current$status, "passed") &&
        all(unlist(current$classes, use.names = FALSE))
    if (!valid) publicationEvidenceAbort("smoke evidence differs")
    invisible(smoke)
}

publicationValidateComparatorBuildReceiptEvidence <- function(path) {
    receipt <- publicationReadJson(path)
    publicationValidateComparatorBuildReceipt(receipt)
    root <- publicationComparatorBuildReceiptRoot(path)
    publicationValidateComparatorSourceEvidence(receipt, root)
    publicationValidateComparatorEnvironmentEvidence(receipt, root)
    publicationValidateComparatorBuildFiles(receipt, root)
    publicationValidateComparatorSmokeEvidence(receipt, root)
    valid <- identical(receipt$session$r_version, as.character(getRversion())) &&
        identical(receipt$session$platform, R.version$platform) &&
        !dir.exists(file.path(publicationComparatorRoot(), "worktrees", receipt$comparator_id))
    if (!valid) publicationEvidenceAbort("build session or cleanup differs")
    list(
        path = path,
        sha256 = publicationFileDigest(path),
        receipt = receipt
    )
}

publicationRestoreReceiptRoot <- function(path) {
    root <- dirname(normalizePath(path, mustWork = TRUE))
    publicationEvidenceOwnedPath(
        root,
        file.path(publicationComparatorRoot(), "restores"),
        "restore root"
    )
    root
}

publicationRestoreCommandEvidence <- function(receipt) {
    commands <- list(
        bootstrap = receipt$bootstrap$command,
        target_renv = receipt$target_renv$command,
        restore = receipt$restore$command,
        audit = receipt$audit$command
    )
    remotes <- lapply(receipt$remote_installs, `[[`, "command")
    c(commands, remotes)
}

publicationValidateRestorePaths <- function(receipt, root) {
    comparator_root <- publicationComparatorRoot()
    expected <- publicationRestorePaths(comparator_root, receipt$restore_id)
    actual <- lapply(receipt$paths, \(path) {
        normalizePath(file.path(comparator_root, path), mustWork = TRUE)
    })
    expected <- lapply(expected, \(path) normalizePath(path, mustWork = TRUE))
    if (!identical(actual, expected)) {
        publicationEvidenceAbort("restore owned paths differ")
    }
    for (path in unlist(actual, use.names = FALSE)) {
        publicationEvidenceOwnedPath(path, root, "restore component")
    }
    invisible(actual)
}

publicationValidateRestoreAuthorityBindings <- function(receipt) {
    bindings <- receipt$environment[c(
        "lockfile", "lock_scope", "repositories", "system_dependencies"
    )]
    valid <- all(vapply(bindings, \(binding) {
        publicationScalarString(binding$path) &&
            publicationComparatorShaValid(binding$sha256) &&
            file.exists(binding$path) &&
            identical(publicationFileDigest(binding$path), binding$sha256)
    }, logical(1)))
    if (!valid) publicationEvidenceAbort("restore authority binding differs")
    invisible(bindings)
}

publicationValidateRestoreReceiptEvidence <- function(path) {
    receipt <- publicationReadJson(path)
    publicationRestoreReceiptBinding(path)
    root <- publicationRestoreReceiptRoot(path)
    publicationValidateRestorePaths(receipt, root)
    publicationValidateRestoreAuthorityBindings(receipt)
    commands <- publicationRestoreCommandEvidence(receipt)
    for (name in names(commands)) {
        publicationValidateCommandEvidence(commands[[name]], root, name)
    }
    publicationValidateScriptEvidence(receipt$restore$script, root, "restore")
    publicationValidateScriptEvidence(receipt$audit$script, root, "audit")
    audit <- publicationReadJson(file.path(root, "restore-audit.json"))
    lock <- publicationReadJson("renv.lock")
    inventory <- publicationRestoredLibraryInventory(
        file.path(root, "library"),
        lock
    )
    valid <- identical(audit, receipt$audit$evidence) && identical(
        publicationObjectDigest(inventory),
        publicationObjectDigest(receipt$library_inventory)
    ) &&
        !isTRUE(receipt$ambient_user_library) &&
        !isTRUE(receipt$cache_shared_with_other_restore) &&
        !isTRUE(receipt$publication_authority)
    if (!valid) publicationEvidenceAbort("restore bytes or isolation differ")
    list(path = path, sha256 = publicationFileDigest(path), receipt = receipt)
}

publicationValidateRestoreReproducibilityEvidence <- function(path) {
    record <- publicationReadJson(path)
    publicationValidateRestoreReproducibilityRecord(record)
    receipt_paths <- vapply(
        record$restore_receipts,
        `[[`,
        character(1),
        "path"
    )
    sweep_paths <- vapply(
        record$namespace_sweeps$bindings,
        `[[`,
        character(1),
        "path"
    )
    lapply(receipt_paths, publicationValidateRestoreReceiptEvidence)
    current <- publicationCreateRestoreReproducibilityRecord(
        dirname(receipt_paths[[1L]]),
        dirname(receipt_paths[[2L]]),
        sweep_paths[[1L]],
        sweep_paths[[2L]]
    )
    if (!identical(record, current)) {
        publicationEvidenceAbort("restore reproducibility record is stale")
    }
    invisible(record)
}

publicationValidateComparatorBuildPairEvidence <- function(path) {
    record <- publicationReadJson(path)
    publicationValidateComparatorBuildPairRecord(record)
    receipt_paths <- vapply(
        record$build_receipts,
        `[[`,
        character(1),
        "path"
    )
    lapply(receipt_paths, publicationValidateComparatorBuildReceiptEvidence)
    current <- publicationCreateComparatorBuildPairRecord(
        receipt_paths[[1L]],
        receipt_paths[[2L]],
        record$dependency_reproducibility$path
    )
    if (!identical(record, current)) {
        publicationEvidenceAbort("build reproducibility record is stale")
    }
    invisible(record)
}
