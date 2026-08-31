publicationExtractPackageArchive <- function(path) {
    root <- tempfile("publication-package-archive-")
    dir.create(root, recursive = TRUE, showWarnings = FALSE)
    result <- processx::run(
        "tar",
        c("-xzf", path, "-C", root),
        error_on_status = FALSE,
        echo = FALSE
    )
    if (result$status != 0L) {
        unlink(root, recursive = TRUE, force = TRUE)
        publicationBuildAbort("Comparator package archive extraction failed")
    }
    root
}

publicationNormalizedArchiveDescription <- function(path) {
    description <- read.dcf(path)
    keep <- setdiff(colnames(description), c("Built", "Packaged"))
    values <- as.list(unname(description[1L, keep]))
    names(values) <- keep
    values[order(names(values), method = "radix")]
}

publicationNormalizedPackageArchive <- function(path) {
    extracted <- publicationExtractPackageArchive(path)
    on.exit(unlink(extracted, recursive = TRUE, force = TRUE), add = TRUE)
    package_root <- file.path(extracted, "MultiScholaR")
    if (!dir.exists(package_root)) {
        publicationBuildAbort("Comparator package archive root differs")
    }
    files <- publicationRestoreFileMap(package_root)
    md5_present <- "MD5" %in% names(files)
    files <- files[names(files) != "MD5"]
    records <- lapply(seq_along(files), \(index) {
        relative <- names(files)[[index]]
        if (identical(relative, "DESCRIPTION")) {
            return(list(
                path = relative,
                kind = "normalized_description",
                sha256 = publicationObjectDigest(
                    publicationNormalizedArchiveDescription(files[[index]])
                )
            ))
        }
        list(
            path = relative,
            kind = "raw_source_file",
            executable = as.integer(file.info(files[[index]])$mode) %% 2L == 1L,
            sha256 = publicationFileDigest(files[[index]])
        )
    })
    list(
        raw_archive_sha256 = publicationFileDigest(path),
        normalized_file_count = as.integer(length(records)),
        md5_manifest_generated = md5_present,
        normalized_content_sha256 = publicationObjectDigest(records)
    )
}

publicationComparatorBuildReceiptBinding <- function(path) {
    receipt <- publicationReadJson(path)
    publicationValidateComparatorBuildReceipt(receipt)
    list(
        path = path,
        sha256 = publicationFileDigest(path),
        receipt_id = receipt$receipt_id,
        comparator_id = receipt$comparator_id,
        revision = receipt$revision,
        raw_archive_sha256 = receipt$build$archive$sha256,
        raw_installed_sha256 =
            receipt$build$installed_inventory$inventory_sha256,
        smoke_sha256 = receipt$smoke$output_sha256
    )
}

publicationCreateComparatorBuildPairRecord <- function(
    first_receipt_path,
    second_receipt_path,
    dependency_reproducibility_path
) {
    receipts <- list(
        publicationReadJson(first_receipt_path),
        publicationReadJson(second_receipt_path)
    )
    bindings <- list(
        publicationComparatorBuildReceiptBinding(first_receipt_path),
        publicationComparatorBuildReceiptBinding(second_receipt_path)
    )
    same_identity <- identical(
        bindings[[1L]]$comparator_id,
        bindings[[2L]]$comparator_id
    ) && identical(bindings[[1L]]$revision, bindings[[2L]]$revision) &&
        identical(
            publicationObjectDigest(receipts[[1L]]$source$source_identity),
            publicationObjectDigest(receipts[[2L]]$source$source_identity)
        )
    archives <- lapply(receipts, \(receipt) {
        publicationNormalizedPackageArchive(receipt$build$archive$path)
    })
    archive_equal <- identical(
        archives[[1L]]$normalized_content_sha256,
        archives[[2L]]$normalized_content_sha256
    )
    package_roots <- lapply(receipts, \(receipt) {
        dirname(receipt$build$installed_inventory$root)
    })
    descriptions <- publicationCompareDescriptions(
        package_roots[[1L]],
        package_roots[[2L]],
        "MultiScholaR"
    )
    installed_paths <- publicationCompareInstalledPaths(
        package_roots[[1L]],
        package_roots[[2L]],
        "MultiScholaR"
    )
    native <- publicationCompareNativeFiles(
        file.path(package_roots[[1L]], "MultiScholaR"),
        file.path(package_roots[[2L]], "MultiScholaR")
    )
    dependency_record <- publicationReadJson(dependency_reproducibility_path)
    publicationValidateRestoreReproducibilityRecord(dependency_record)
    dependency_binding <- list(
        path = dependency_reproducibility_path,
        sha256 = publicationFileDigest(dependency_reproducibility_path)
    )
    smoke_equal <- identical(
        bindings[[1L]]$smoke_sha256,
        bindings[[2L]]$smoke_sha256
    )
    passed <- same_identity && archive_equal && descriptions$reproducible &&
        installed_paths$classified && native$reproducible && smoke_equal
    list(
        schema = "multischolar.omics_publication_comparator_build_pair",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-064",
        comparator_id = bindings[[1L]]$comparator_id,
        revision = bindings[[1L]]$revision,
        status = if (passed) "passed" else "failed",
        build_receipts = bindings,
        dependency_reproducibility = dependency_binding,
        raw_archives_equal = identical(
            bindings[[1L]]$raw_archive_sha256,
            bindings[[2L]]$raw_archive_sha256
        ),
        normalized_archives = archives,
        normalized_archives_equal = archive_equal,
        raw_installed_trees_equal = identical(
            bindings[[1L]]$raw_installed_sha256,
            bindings[[2L]]$raw_installed_sha256
        ),
        descriptions = descriptions,
        installed_paths = installed_paths,
        native = native,
        smoke_equal = smoke_equal,
        smoke_sha256 = bindings[[1L]]$smoke_sha256,
        reproducible = passed,
        publication_authority = FALSE
    )
}

publicationValidateComparatorBuildPairRecord <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "owner_ticket_id", "comparator_id",
        "revision", "status", "build_receipts", "dependency_reproducibility",
        "raw_archives_equal", "normalized_archives",
        "normalized_archives_equal", "raw_installed_trees_equal",
        "descriptions", "installed_paths", "native", "smoke_equal",
        "smoke_sha256", "reproducible", "publication_authority"
    ), "Comparator build pair record")
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_comparator_build_pair"
    ) && identical(record$schema_version, "1.0.0") &&
        identical(record$owner_ticket_id, "OMICS-ART-064") &&
        record$comparator_id %in% names(publicationComparatorRevisions()) &&
        identical(
            record$revision,
            publicationComparatorRevisions()[[record$comparator_id]]
        ) && identical(record$status, "passed") &&
        isTRUE(record$normalized_archives_equal) &&
        isTRUE(record$descriptions$reproducible) &&
        isTRUE(record$installed_paths$classified) &&
        !length(record$installed_paths$unexpected) &&
        isTRUE(record$native$reproducible) && isTRUE(record$smoke_equal) &&
        isTRUE(record$reproducible) && !isTRUE(record$publication_authority)
    if (!valid) publicationBuildAbort("Comparator build pair record differs")
    invisible(record)
}
