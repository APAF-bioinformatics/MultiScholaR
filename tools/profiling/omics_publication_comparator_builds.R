publicationComparatorBuildPaths <- function(
    comparator_id,
    repetition,
    attempt = 1L,
    root = publicationComparatorRoot()
) {
    build_root <- file.path(
        root,
        "comparator-builds",
        comparator_id,
        paste0("repetition-", repetition, "-attempt-", attempt)
    )
    list(
        root = build_root,
        staging = file.path(build_root, "staging"),
        package_library = file.path(build_root, "package-library"),
        smoke_output = file.path(build_root, "smoke-output.json"),
        receipt = file.path(build_root, "receipt.json")
    )
}

publicationDependencyFingerprint <- function(library, lock) {
    inventory <- publicationRestoredLibraryInventory(library, lock)
    list(
        package_count = inventory$package_count,
        inventory_sha256 = inventory$inventory_sha256
    )
}

publicationPrepareDependencyOverlay <- function(
    package_library,
    dependency_library,
    lock
) {
    dir.create(package_library, recursive = TRUE, showWarnings = FALSE)
    existing <- list.files(package_library, all.files = TRUE, no.. = TRUE)
    if (length(existing)) {
        publicationBuildAbort("Comparator package overlay is not empty")
    }
    records <- lapply(names(lock$Packages), \(package) {
        source <- normalizePath(
            file.path(dependency_library, package),
            mustWork = TRUE
        )
        target <- file.path(package_library, package)
        if (!isTRUE(file.symlink(source, target))) {
            publicationBuildAbort("Comparator dependency overlay link failed")
        }
        if (!identical(normalizePath(target, mustWork = TRUE), source)) {
            publicationBuildAbort("Comparator dependency overlay target differs")
        }
        list(
            package = package,
            version = lock$Packages[[package]]$Version,
            target_sha256 =
                publicationPackageTreeSummary(source)$inventory_sha256
        )
    })
    names(records) <- names(lock$Packages)
    list(
        package_count = as.integer(length(records)),
        records_sha256 = publicationObjectDigest(records)
    )
}

publicationComparatorSmokeScript <- function(
    package_library,
    dependency_library,
    output
) {
    paste(c(
        paste0(
            "package_library <- ",
            publicationBuildRLiteral(package_library)
        ),
        paste0(
            "dependency_library <- ",
            publicationBuildRLiteral(dependency_library)
        ),
        paste0("output <- ", publicationBuildRLiteral(output)),
        ".libPaths(c(package_library, dependency_library, .Library))",
        paste0(
            "package_path <- find.package(\"MultiScholaR\", ",
            "lib.loc = package_library)"
        ),
        paste0(
            "stopifnot(startsWith(normalizePath(package_path), ",
            "paste0(normalizePath(package_library), .Platform$file.sep)))"
        ),
        "invisible(loadNamespace(\"MultiScholaR\"))",
        "input <- matrix(c(1, 4, 16, 64), nrow = 2L)",
        "transformed <- MultiScholaR::log2Transformation(input)",
        "expected <- log2(input + min(input) / 100)",
        "stopifnot(identical(transformed, expected))",
        paste0(
            "detected <- MultiScholaR::detectProteomicsFormat(",
            "c(\"protein.group\", \"protein.ids\", \"protein.names\", ",
            "\"precursor.id\", \"modified.sequence\", ",
            "\"stripped.sequence\", \"precursor.charge\", ",
            "\"q.value\", \"pg.q.value\", \"run\"), ",
            "\"report.tsv\")"
        ),
        "stopifnot(identical(detected$format, \"diann\"))",
        paste0(
            "classes <- c(\"PeptideQuantitativeData\", ",
            "\"ProteinQuantitativeData\", \"MetaboliteAssayData\", ",
            "\"LipidomicsAssayData\")"
        ),
        "class_status <- vapply(classes, methods::isClass, logical(1))",
        "stopifnot(all(class_status))",
        "exports <- sort(getNamespaceExports(\"MultiScholaR\"))",
        "jsonlite::write_json(list(",
        "    status = \"passed\",",
        "    package_version = as.character(utils::packageVersion(",
        "        \"MultiScholaR\"",
        "    )),",
        "    transformed = as.list(as.numeric(transformed)),",
        "    detected_format = detected$format,",
        "    detected_confidence = detected$confidence,",
        "    classes = as.list(class_status),",
        "    export_count = length(exports),",
        "    exports_sha256 = digest::digest(exports, algo = \"sha256\")",
        "), output, auto_unbox = TRUE, pretty = TRUE, null = \"null\")"
    ), collapse = "\n")
}

publicationRunComparatorSmoke <- function(
    package_library,
    dependency_library,
    output,
    build_root
) {
    script_path <- file.path(build_root, "comparator-smoke.R")
    script <- publicationBuildWriteScript(
        publicationComparatorSmokeScript(
            package_library,
            dependency_library,
            output
        ),
        script_path
    )
    site_library <- file.path(build_root, "empty-site-library")
    environment <- publicationBuildEnvironment(package_library, site_library)
    environment[["HOME"]] <- file.path(build_root, "home")
    environment[["TMPDIR"]] <- file.path(build_root, "tmp")
    command <- publicationBuildRun(
        file.path(R.home("bin"), "Rscript"),
        c("--vanilla", script_path),
        build_root,
        file.path(build_root, "logs"),
        environment,
        timeout_seconds = 600
    )
    publicationBuildRequireSuccess(command, "Comparator smoke")
    evidence <- publicationReadJson(output)
    valid <- identical(evidence$status, "passed") &&
        identical(evidence$package_version, "0.5.0") &&
        identical(evidence$detected_format, "diann") &&
        all(unlist(evidence$classes, use.names = FALSE))
    if (!valid) publicationBuildAbort("Comparator smoke evidence differs")
    list(
        status = "passed",
        script = script,
        command = command,
        output_sha256 = publicationFileDigest(output),
        evidence = evidence
    )
}

publicationComparatorRestoreReceipt <- function(restore_root) {
    path <- file.path(restore_root, "receipt.json")
    binding <- publicationRestoreReceiptBinding(path)
    receipt <- publicationReadJson(path)
    list(binding = binding, receipt = receipt)
}

publicationBuildComparatorOnce <- function(
    comparator_id,
    repetition,
    restore_root,
    attempt = 1L,
    root = publicationComparatorRoot()
) {
    revision <- publicationComparatorRevisions()[[comparator_id]]
    if (is.null(revision)) publicationBuildAbort("Comparator id is not fixed")
    paths <- publicationComparatorBuildPaths(
        comparator_id,
        repetition,
        attempt,
        root
    )
    if (file.exists(paths$root) || dir.exists(paths$root)) {
        publicationBuildAbort("Comparator build root already exists")
    }
    dir.create(paths$root, recursive = TRUE, showWarnings = FALSE)
    restore <- publicationComparatorRestoreReceipt(restore_root)
    lock <- publicationReadJson("renv.lock")
    dependency_library <- file.path(restore_root, "library")
    before <- publicationDependencyFingerprint(dependency_library, lock)
    if (!identical(
        before$inventory_sha256,
        restore$receipt$library_inventory$inventory_sha256
    )) {
        publicationBuildAbort("Dependency library differs before build")
    }
    overlay <- publicationPrepareDependencyOverlay(
        paths$package_library,
        dependency_library,
        lock
    )
    worktree <- publicationCreateComparatorWorktree(
        comparator_id,
        revision,
        root
    )
    worktree_removed <- FALSE
    on.exit({
        if (!worktree_removed && dir.exists(worktree$path)) {
            publicationRemoveComparatorWorktree(worktree, root)
        }
    }, add = TRUE)
    source <- publicationExportComparatorSource(worktree, paths$staging)
    build <- publicationBuildPackage(
        source,
        paths$root,
        paths$package_library,
        dependency_library,
        revision
    )
    smoke <- publicationRunComparatorSmoke(
        paths$package_library,
        dependency_library,
        paths$smoke_output,
        paths$root
    )
    after <- publicationDependencyFingerprint(dependency_library, lock)
    if (!identical(before, after)) {
        publicationBuildAbort("Dependency library changed during build")
    }
    publicationRemoveComparatorWorktree(worktree, root)
    worktree_removed <- TRUE
    receipt <- list(
        schema = "multischolar.omics_publication_comparator_build",
        schema_version = "1.0.0",
        receipt_id = paste(comparator_id, repetition, attempt, sep = "."),
        comparator_id = comparator_id,
        revision = revision,
        source = source,
        environment = list(
            lockfile_sha256 = publicationFileDigest("renv.lock"),
            restore_receipt = restore$binding,
            dependency_overlay = overlay,
            dependency_before = before,
            dependency_after = after
        ),
        restore = list(status = "passed"),
        build = build,
        session = list(
            r_version = as.character(getRversion()),
            platform = R.version$platform
        ),
        smoke = smoke,
        cleanup = list(
            worktree_removed = TRUE,
            dependency_library_retained = TRUE,
            package_library_retained = TRUE
        ),
        status = "verified",
        publication_authority = FALSE
    )
    publicationValidateComparatorBuildReceipt(receipt)
    publicationWriteJson(receipt, paths$receipt)
    receipt
}
