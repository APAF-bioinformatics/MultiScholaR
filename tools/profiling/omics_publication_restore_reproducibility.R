publicationRestoreReproAbort <- function(message) {
    publicationComparatorAbort(message, "restore_reproducibility_error")
}

publicationRestoreFileMap <- function(root) {
    paths <- list.files(
        root,
        recursive = TRUE,
        full.names = TRUE,
        all.files = TRUE,
        no.. = TRUE
    )
    paths <- paths[file.exists(paths) & !dir.exists(paths)]
    relative <- substring(paths, nchar(root) + 2L)
    paths <- paths[order(relative, method = "radix")]
    relative <- relative[order(relative, method = "radix")]
    stats::setNames(paths, relative)
}

publicationNativeKind <- function(path) {
    connection <- file(path, open = "rb")
    on.exit(close(connection), add = TRUE)
    magic <- readBin(connection, what = "raw", n = 8L)
    if (identical(magic[1:4], as.raw(c(0x7f, 0x45, 0x4c, 0x46)))) {
        return("elf")
    }
    if (identical(magic, charToRaw("!<arch>\n")) && endsWith(path, ".a")) {
        return("archive")
    }
    NULL
}

publicationPrintableRuns <- function(bytes) {
    values <- as.integer(bytes)
    printable <- values == 9L | (values >= 32L & values <= 126L)
    encoded <- rle(printable)
    ends <- cumsum(encoded$lengths)
    starts <- ends - encoded$lengths + 1L
    keep <- encoded$values & encoded$lengths >= 4L
    data.frame(start = starts[keep], end = ends[keep])
}

publicationCanonicalBuildString <- function(value) {
    value <- gsub(
        paste0(
            "/home/[^/]+/Projects/APAF/MultiScholaR/",
            "\\.omics-publication-comparators/restores/[^/]+/",
            "tmp/Rtmp[^/]+/R\\.INSTALL[^/]+/"
        ),
        "<BUILD_ROOT>/",
        value
    )
    value <- gsub(
        "/build/multischolar/tmp/Rtmp[^/]+/R\\.INSTALL[^/]+/",
        "<BUILD_ROOT>/",
        value
    )
    sub(
        "^((TBB|TBBmalloc): BUILD_DATE[[:space:]]+).*$",
        "\\1<SOURCE_DATE_EPOCH>",
        value
    )
}

publicationNormalizeVolatileStrings <- function(bytes) {
    runs <- publicationPrintableRuns(bytes)
    normalized <- character()
    for (index in seq_len(nrow(runs))) {
        positions <- runs$start[[index]]:runs$end[[index]]
        value <- rawToChar(bytes[positions])
        canonical <- publicationCanonicalBuildString(value)
        if (!identical(value, canonical)) {
            normalized <- c(normalized, canonical)
            bytes[positions] <- as.raw(0L)
        }
    }
    list(
        bytes = bytes,
        strings = sort(normalized, method = "radix")
    )
}

publicationCommandOutputDigest <- function(command, arguments) {
    result <- processx::run(
        command,
        arguments,
        error_on_status = FALSE,
        echo = FALSE
    )
    if (result$status != 0L) {
        publicationRestoreReproAbort(paste(command, "inspection failed"))
    }
    digest::digest(result$stdout, algo = "sha256", serialize = FALSE)
}

publicationNormalizedElfRecord <- function(path) {
    stripped <- tempfile("publication-elf-")
    on.exit(unlink(stripped), add = TRUE)
    result <- processx::run(
        "objcopy",
        c("--strip-debug", "--remove-section=.note.gnu.build-id", path, stripped),
        error_on_status = FALSE,
        echo = FALSE
    )
    if (result$status != 0L) {
        publicationRestoreReproAbort("ELF normalization failed")
    }
    connection <- file(stripped, open = "rb")
    bytes <- readBin(connection, what = "raw", n = file.info(stripped)$size)
    close(connection)
    normalized <- publicationNormalizeVolatileStrings(bytes)
    list(
        masked_sha256 = digest::digest(
            normalized$bytes,
            algo = "sha256",
            serialize = FALSE
        ),
        volatile_strings = as.list(normalized$strings),
        dynamic_symbols_sha256 = publicationCommandOutputDigest(
            "nm",
            c("-D", "--defined-only", stripped)
        ),
        dynamic_contract_sha256 = publicationCommandOutputDigest(
            "readelf",
            c("-d", stripped)
        )
    )
}

publicationArchiveMembers <- function(path) {
    output <- processx::run(
        "ar",
        c("t", path),
        error_on_status = FALSE,
        echo = FALSE
    )
    if (output$status != 0L) {
        publicationRestoreReproAbort("Static archive listing failed")
    }
    strsplit(trimws(output$stdout), "\n", fixed = TRUE)[[1L]]
}

publicationObjectSectionDigest <- function(path, section) {
    output <- tempfile("publication-section-")
    on.exit(unlink(output), add = TRUE)
    result <- processx::run(
        "objcopy",
        c(paste0("--dump-section=", section, "=", output), path, "/dev/null"),
        error_on_status = FALSE,
        echo = FALSE
    )
    if (result$status != 0L || !file.exists(output)) return(NA_character_)
    connection <- file(output, open = "rb")
    bytes <- readBin(connection, what = "raw", n = file.info(output)$size)
    close(connection)
    normalized <- publicationNormalizeVolatileStrings(bytes)
    publicationObjectDigest(list(
        masked_sha256 = digest::digest(
            normalized$bytes,
            algo = "sha256",
            serialize = FALSE
        ),
        volatile_strings = as.list(normalized$strings)
    ))
}

publicationArchiveMemberRecord <- function(archive, member, index) {
    object <- tempfile("publication-member-")
    on.exit(unlink(object), add = TRUE)
    result <- processx::run(
        "ar",
        c("p", archive, member),
        stdout = object,
        error_on_status = FALSE,
        echo = FALSE
    )
    if (result$status != 0L) {
        publicationRestoreReproAbort("Static archive member extraction failed")
    }
    sections <- c(
        ".text", ".rodata", ".rodata.cst4", ".rodata.cst8",
        ".rodata.cst16", ".rodata.cst32", ".data", ".rela.text",
        ".rela.rodata", ".eh_frame", ".rela.eh_frame"
    )
    list(
        index = as.integer(index),
        member = member,
        runtime_sections = as.list(stats::setNames(vapply(
            sections,
            \(section) publicationObjectSectionDigest(object, section),
            character(1)
        ), sections)),
        global_symbols_sha256 = publicationCommandOutputDigest(
            "nm",
            c("-g", "--defined-only", object)
        )
    )
}

publicationNormalizedArchiveRecord <- function(path) {
    members <- publicationArchiveMembers(path)
    records <- lapply(seq_along(members), \(index) {
        publicationArchiveMemberRecord(path, members[[index]], index)
    })
    list(
        member_count = as.integer(length(records)),
        runtime_contract_sha256 = publicationObjectDigest(records)
    )
}

publicationNativeFileRecords <- function(root) {
    files <- publicationRestoreFileMap(root)
    kinds <- lapply(files, publicationNativeKind)
    keep <- !vapply(kinds, is.null, logical(1))
    files <- files[keep]
    kinds <- unlist(kinds[keep], use.names = FALSE)
    records <- lapply(seq_along(files), \(index) {
        record <- if (identical(kinds[[index]], "elf")) {
            publicationNormalizedElfRecord(files[[index]])
        } else {
            publicationNormalizedArchiveRecord(files[[index]])
        }
        c(list(path = names(files)[[index]], kind = kinds[[index]]), record)
    })
    names(records) <- names(files)
    records
}

publicationCompareNativeFiles <- function(first_root, second_root) {
    records <- list(
        publicationNativeFileRecords(first_root),
        publicationNativeFileRecords(second_root)
    )
    paths_equal <- identical(names(records[[1L]]), names(records[[2L]]))
    mismatches <- if (paths_equal) {
        names(records[[1L]])[vapply(names(records[[1L]]), \(path) {
            !identical(records[[1L]][[path]], records[[2L]][[path]])
        }, logical(1))]
    } else {
        union(names(records[[1L]]), names(records[[2L]]))
    }
    list(
        paths_equal = paths_equal,
        native_file_count = as.integer(length(records[[1L]])),
        mismatches = as.list(mismatches),
        first_sha256 = publicationObjectDigest(records[[1L]]),
        second_sha256 = publicationObjectDigest(records[[2L]]),
        reproducible = paths_equal && !length(mismatches)
    )
}

publicationCompareSourceCaches <- function(first_root, second_root) {
    inventories <- list(
        publicationInstalledInventory(first_root),
        publicationInstalledInventory(second_root)
    )
    list(
        first_file_count = as.integer(inventories[[1L]]$file_count),
        second_file_count = as.integer(inventories[[2L]]$file_count),
        first_sha256 = inventories[[1L]]$inventory_sha256,
        second_sha256 = inventories[[2L]]$inventory_sha256,
        identical = identical(
            inventories[[1L]]$inventory_sha256,
            inventories[[2L]]$inventory_sha256
        )
    )
}

publicationNormalizedDescription <- function(path) {
    description <- read.dcf(path)
    keep <- setdiff(colnames(description), c("Built", "Packaged"))
    values <- as.list(unname(description[1L, keep]))
    names(values) <- keep
    values[order(names(values), method = "radix")]
}

publicationCompareDescriptions <- function(first_root, second_root, packages) {
    mismatches <- packages[vapply(packages, \(package) {
        !identical(
            publicationNormalizedDescription(file.path(
                first_root,
                package,
                "DESCRIPTION"
            )),
            publicationNormalizedDescription(file.path(
                second_root,
                package,
                "DESCRIPTION"
            ))
        )
    }, logical(1))]
    list(
        package_count = as.integer(length(packages)),
        mismatches = as.list(mismatches),
        reproducible = !length(mismatches)
    )
}

publicationInstalledDifferenceCategory <- function(path, kind) {
    if (identical(kind, "elf") || identical(kind, "archive")) return("native")
    if (identical(path, "DESCRIPTION")) return("description")
    if (startsWith(path, "Meta/")) return("generated_metadata")
    if (startsWith(path, "help/")) return("generated_help")
    if (grepl("^R/.*\\.(rdb|rdx)$", path)) return("generated_lazyload")
    if (grepl("(^|/)library/[^/]+/(DESCRIPTION|Meta/|R/|libs/|bin/)", path)) {
        return("nested_tool_runtime")
    }
    if (startsWith(path, "bin/") || grepl(
        "(^|/)lib/(cmake|pkgconfig)/|libhdf5\\.settings$",
        path
    )) {
        return("relocatable_build_text")
    }
    "unexpected"
}

publicationCompareInstalledPaths <- function(first_root, second_root, packages) {
    categories <- character()
    unexpected <- character()
    for (package in packages) {
        maps <- list(
            publicationRestoreFileMap(file.path(first_root, package)),
            publicationRestoreFileMap(file.path(second_root, package))
        )
        paths <- union(names(maps[[1L]]), names(maps[[2L]]))
        changed <- paths[vapply(paths, \(path) {
            if (is.null(maps[[1L]][[path]]) || is.null(maps[[2L]][[path]])) {
                return(TRUE)
            }
            !identical(
                publicationFileDigest(maps[[1L]][[path]]),
                publicationFileDigest(maps[[2L]][[path]])
            )
        }, logical(1))]
        for (path in changed) {
            kind <- if (!is.null(maps[[1L]][[path]])) {
                publicationNativeKind(maps[[1L]][[path]])
            } else {
                NULL
            }
            category <- publicationInstalledDifferenceCategory(path, kind)
            categories <- c(categories, category)
            if (identical(category, "unexpected")) {
                unexpected <- c(unexpected, paste(package, path, sep = "/"))
            }
        }
    }
    counts <- table(factor(categories, levels = c(
        "description", "generated_metadata", "generated_help",
        "generated_lazyload", "native", "nested_tool_runtime",
        "relocatable_build_text", "unexpected"
    )))
    list(
        changed_file_count = as.integer(length(categories)),
        category_counts = as.list(as.integer(counts)),
        category_names = as.list(names(counts)),
        unexpected = as.list(unexpected),
        classified = !length(unexpected)
    )
}

publicationCanonicalRelocatableText <- function(path) {
    lines <- readLines(path, warn = FALSE)
    gsub(
        paste0(
            "/home/[^/]+/Projects/APAF/MultiScholaR/",
            "\\.omics-publication-comparators/restores/[^/]+/",
            "tmp/Rtmp[^/]+/R\\.INSTALL[^/]+/"
        ),
        "<BUILD_ROOT>/",
        lines
    )
}

publicationRelocatableTextPaths <- function(root) {
    files <- publicationRestoreFileMap(root)
    keep <- startsWith(names(files), "bin/") |
        grepl(
            "(^|/)lib/(cmake|pkgconfig)/|libhdf5\\.settings$",
            names(files)
        )
    keep <- keep & vapply(files, \(path) {
        is.null(publicationNativeKind(path))
    }, logical(1))
    files[keep]
}

publicationCompareRelocatableText <- function(first_root, second_root) {
    files <- list(
        publicationRelocatableTextPaths(first_root),
        publicationRelocatableTextPaths(second_root)
    )
    paths_equal <- identical(names(files[[1L]]), names(files[[2L]]))
    mismatches <- if (paths_equal) {
        names(files[[1L]])[vapply(names(files[[1L]]), \(path) {
            !identical(
                publicationCanonicalRelocatableText(files[[1L]][[path]]),
                publicationCanonicalRelocatableText(files[[2L]][[path]])
            )
        }, logical(1))]
    } else {
        union(names(files[[1L]]), names(files[[2L]]))
    }
    list(
        path_count = as.integer(length(files[[1L]])),
        paths_equal = paths_equal,
        mismatches = as.list(mismatches),
        reproducible = paths_equal && !length(mismatches)
    )
}

publicationNamespaceExpression <- function(package, library) {
    paste0(
        ".libPaths(c(", publicationBuildRLiteral(library), ", .Library)); ",
        "path <- find.package(", publicationBuildRLiteral(package),
        ", lib.loc = ", publicationBuildRLiteral(library), "); ",
        "stopifnot(startsWith(normalizePath(path), paste0(normalizePath(",
        publicationBuildRLiteral(library), "), .Platform$file.sep))); ",
        "invisible(loadNamespace(", publicationBuildRLiteral(package), ")); ",
        "cat(unname(read.dcf(file.path(path, \"DESCRIPTION\"))",
        "[1L, \"Version\"]))"
    )
}

publicationNamespaceSweep <- function(library, packages, log_dir) {
    site_library <- file.path(dirname(library), "empty-site-library")
    environment <- publicationBuildEnvironment(library, site_library)
    environment[["HOME"]] <- file.path(dirname(library), "home")
    environment[["TMPDIR"]] <- file.path(dirname(library), "tmp")
    records <- lapply(packages, \(package) {
        started <- proc.time()[["elapsed"]]
        result <- processx::run(
            file.path(R.home("bin"), "Rscript"),
            c("--vanilla", "-e", publicationNamespaceExpression(
                package,
                library
            )),
            env = environment,
            stdout = "|",
            stderr = "|",
            timeout = 300,
            error_on_status = FALSE,
            echo = FALSE
        )
        list(
            package = package,
            status = if (result$status == 0L) "passed" else "failed",
            exit_status = as.integer(result$status),
            timed_out = isTRUE(result$timeout),
            version = trimws(result$stdout),
            stderr_sha256 = digest::digest(
                publicationCanonicalBuildString(result$stderr),
                algo = "sha256",
                serialize = FALSE
            ),
            elapsed_seconds = proc.time()[["elapsed"]] - started
        )
    })
    names(records) <- packages
    dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
    path <- file.path(log_dir, "namespace-sweep.json")
    publicationWriteJson(records, path)
    list(
        package_count = as.integer(length(records)),
        failed = as.list(names(records)[vapply(records, \(record) {
            !identical(record$status, "passed")
        }, logical(1))]),
        records_sha256 = publicationObjectDigest(lapply(records, \(record) {
            record$elapsed_seconds <- NULL
            record
        })),
        evidence_sha256 = publicationFileDigest(path),
        records = records
    )
}

publicationCompareNamespaceSweeps <- function(first, second, lock) {
    packages <- names(lock$Packages)
    versions_valid <- all(vapply(packages, \(package) {
        identical(first$records[[package]]$version, lock$Packages[[package]]$Version) &&
            identical(second$records[[package]]$version, lock$Packages[[package]]$Version)
    }, logical(1)))
    status_equal <- all(vapply(packages, \(package) {
        identical(first$records[[package]]$status, "passed") &&
            identical(second$records[[package]]$status, "passed") &&
            identical(
                first$records[[package]]$stderr_sha256,
                second$records[[package]]$stderr_sha256
            )
    }, logical(1)))
    list(
        bindings = list(first = first$binding, second = second$binding),
        package_count = as.integer(length(packages)),
        first_failed = first$failed,
        second_failed = second$failed,
        versions_valid = versions_valid,
        status_and_stderr_equal = status_equal,
        reproducible = versions_valid && status_equal &&
            !length(first$failed) && !length(second$failed)
    )
}

publicationReadNamespaceSweep <- function(path) {
    records <- publicationReadJson(path)
    list(
        binding = list(
            path = path,
            sha256 = publicationFileDigest(path)
        ),
        records = records,
        failed = as.list(names(records)[vapply(records, \(record) {
            !identical(record$status, "passed")
        }, logical(1))]),
        evidence_sha256 = publicationFileDigest(path)
    )
}

publicationRestoreReceiptBinding <- function(path) {
    receipt <- publicationReadJson(path)
    valid <- identical(receipt$status, "passed") &&
        identical(receipt$library_inventory$package_count, 384L) &&
        isTRUE(receipt$audit$evidence$official_schema_valid) &&
        isTRUE(receipt$audit$evidence$all_packages_owned) &&
        !isTRUE(receipt$ambient_user_library) &&
        !isTRUE(receipt$cache_shared_with_other_restore)
    if (!valid) publicationRestoreReproAbort("Restore receipt is not admissible")
    list(
        path = path,
        sha256 = publicationFileDigest(path),
        restore_id = receipt$restore_id,
        lockfile_sha256 = receipt$environment$lockfile$sha256,
        raw_installed_inventory_sha256 =
            receipt$library_inventory$inventory_sha256
    )
}

publicationRestoreNormalizationPolicy <- function() {
    list(
        raw_tree_identity_required = FALSE,
        exact_source_archive_identity_required = TRUE,
        exact_normalized_description_required = TRUE,
        generated_install_products = as.list(c(
            "Meta", "help", "R_lazyload", "nested_pak_runtime"
        )),
        generated_product_authority =
            "exact_source_archives_plus_clean_namespace_loads",
        elf_removed_sections = as.list(c(
            "debug_sections", ".note.gnu.build-id"
        )),
        elf_volatile_strings = as.list(c(
            "owned_Rtmp_R.INSTALL_prefix", "TBB_BUILD_DATE"
        )),
        elf_executable_bytes_required = TRUE,
        elf_dynamic_symbols_and_dependencies_required = TRUE,
        static_archive_runtime_sections_and_symbols_required = TRUE,
        relocatable_text_path_only_normalization = TRUE,
        unexpected_difference_count_required = 0L
    )
}

publicationCreateRestoreReproducibilityRecord <- function(
    first_root,
    second_root,
    first_namespace_path,
    second_namespace_path
) {
    lock <- publicationReadJson("renv.lock")
    packages <- names(lock$Packages)
    receipt_paths <- c(
        file.path(first_root, "receipt.json"),
        file.path(second_root, "receipt.json")
    )
    bindings <- lapply(receipt_paths, publicationRestoreReceiptBinding)
    same_lock <- identical(
        bindings[[1L]]$lockfile_sha256,
        bindings[[2L]]$lockfile_sha256
    ) && identical(
        bindings[[1L]]$lockfile_sha256,
        publicationFileDigest("renv.lock")
    )
    source_caches <- publicationCompareSourceCaches(
        file.path(first_root, "source-cache"),
        file.path(second_root, "source-cache")
    )
    descriptions <- publicationCompareDescriptions(
        file.path(first_root, "library"),
        file.path(second_root, "library"),
        packages
    )
    installed_paths <- publicationCompareInstalledPaths(
        file.path(first_root, "library"),
        file.path(second_root, "library"),
        packages
    )
    native <- publicationCompareNativeFiles(
        file.path(first_root, "library"),
        file.path(second_root, "library")
    )
    relocatable_text <- publicationCompareRelocatableText(
        file.path(first_root, "library", "Rhdf5lib"),
        file.path(second_root, "library", "Rhdf5lib")
    )
    namespace_sweeps <- publicationCompareNamespaceSweeps(
        publicationReadNamespaceSweep(first_namespace_path),
        publicationReadNamespaceSweep(second_namespace_path),
        lock
    )
    passed <- same_lock && source_caches$identical &&
        descriptions$reproducible && installed_paths$classified &&
        native$reproducible && relocatable_text$reproducible &&
        namespace_sweeps$reproducible
    list(
        schema = "multischolar.omics_publication_restore_reproducibility",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-064",
        status = if (passed) "passed" else "failed",
        lockfile_sha256 = publicationFileDigest("renv.lock"),
        restore_receipts = bindings,
        raw_installed_trees_equal = identical(
            bindings[[1L]]$raw_installed_inventory_sha256,
            bindings[[2L]]$raw_installed_inventory_sha256
        ),
        normalization_policy = publicationRestoreNormalizationPolicy(),
        source_caches = source_caches,
        descriptions = descriptions,
        installed_paths = installed_paths,
        native = native,
        relocatable_text = relocatable_text,
        namespace_sweeps = namespace_sweeps,
        reproducible = passed,
        publication_authority = FALSE
    )
}

publicationValidateRestoreReproducibilityRecord <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "owner_ticket_id", "status",
        "lockfile_sha256", "restore_receipts", "raw_installed_trees_equal",
        "normalization_policy", "source_caches", "descriptions",
        "installed_paths", "native", "relocatable_text", "namespace_sweeps",
        "reproducible", "publication_authority"
    ), "Restore reproducibility record")
    namespace_names <- c(
        "bindings", "package_count", "first_failed", "second_failed",
        "versions_valid", "status_and_stderr_equal", "reproducible"
    )
    namespace_bindings <- record$namespace_sweeps$bindings
    namespace_valid <- is.list(namespace_bindings) &&
        identical(names(namespace_bindings), c("first", "second")) &&
        all(vapply(namespace_bindings, \(binding) {
            publicationScalarString(binding$path) &&
                publicationComparatorShaValid(binding$sha256) &&
                file.exists(binding$path) &&
                identical(publicationFileDigest(binding$path), binding$sha256)
        }, logical(1)))
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_restore_reproducibility"
    ) && identical(record$schema_version, "1.0.0") &&
        identical(record$owner_ticket_id, "OMICS-ART-064") &&
        identical(record$status, "passed") &&
        identical(record$lockfile_sha256, publicationFileDigest("renv.lock")) &&
        !isTRUE(record$raw_installed_trees_equal) &&
        isTRUE(record$source_caches$identical) &&
        isTRUE(record$descriptions$reproducible) &&
        isTRUE(record$installed_paths$classified) &&
        identical(length(record$installed_paths$unexpected), 0L) &&
        isTRUE(record$native$reproducible) &&
        isTRUE(record$relocatable_text$reproducible) &&
        identical(names(record$namespace_sweeps), namespace_names) &&
        namespace_valid &&
        isTRUE(record$namespace_sweeps$reproducible) &&
        isTRUE(record$reproducible) && !isTRUE(record$publication_authority)
    if (!valid) {
        publicationRestoreReproAbort("Restore reproducibility record differs")
    }
    invisible(record)
}
