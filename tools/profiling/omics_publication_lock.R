publicationLockAbort <- function(message) {
    publicationAbort(message, "multischolar_publication_lock_error")
}

publicationLockRequiredFields <- function() {
    c("Depends", "Imports", "LinkingTo")
}

publicationLockSort <- function(values) {
    sort(unique(values), method = "radix")
}

publicationLockStandardPackages <- function() {
    unique(c(
        tools:::.get_standard_package_names()$base,
        tools:::.get_standard_package_names()$recommended,
        "R"
    ))
}

publicationLockParseField <- function(values) {
    values <- unlist(values, use.names = FALSE)
    values <- values[!is.na(values) & nzchar(values)]
    if (!length(values)) return(character())
    dependencies <- unlist(strsplit(values, ",", fixed = TRUE), use.names = FALSE)
    publicationLockSort(trimws(sub(
        "[[:space:]]*\\(.*$",
        "",
        dependencies
    )))
}

publicationLockGitDescription <- function(revision) {
    publicationComparatorResolveRevision(revision)
    path <- tempfile("multischolar-description-")
    on.exit(unlink(path, force = TRUE), add = TRUE)
    result <- processx::run(
        "git",
        c("show", paste0(revision, ":DESCRIPTION")),
        wd = .PUBLICATION_REPO_ROOT,
        stdout = path,
        stderr = "|",
        error_on_status = FALSE,
        echo = FALSE
    )
    if (result$status != 0L || !file.exists(path)) {
        publicationLockAbort("Comparator DESCRIPTION could not be read")
    }
    description <- read.dcf(path)
    if (nrow(description) != 1L) {
        publicationLockAbort("Comparator DESCRIPTION is not singular")
    }
    required_fields <- publicationLockRequiredFields()
    required <- publicationLockParseField(
        description[1L, intersect(required_fields, colnames(description))]
    )
    suggests <- publicationLockParseField(
        description[1L, intersect("Suggests", colnames(description))]
    )
    remotes <- publicationLockParseField(
        description[1L, intersect("Remotes", colnames(description))]
    )
    list(
        revision = revision,
        description_sha256 = publicationFileDigest(path),
        package = unname(description[1L, "Package"]),
        version = unname(description[1L, "Version"]),
        required = as.list(required),
        suggests = as.list(suggests),
        remotes = as.list(remotes)
    )
}

publicationLockPackageAlias <- function(packages) {
    packages[packages == "GlimmaV2"] <- "Glimma"
    packages
}

publicationLockRootGroups <- function() {
    descriptions <- lapply(
        publicationComparatorRevisions(),
        publicationLockGitDescription
    )
    standard <- publicationLockStandardPackages()
    required <- publicationLockPackageAlias(unique(unlist(
        lapply(descriptions, `[[`, "required"),
        use.names = FALSE
    )))
    suggests <- publicationLockPackageAlias(unique(unlist(
        lapply(descriptions, `[[`, "suggests"),
        use.names = FALSE
    )))
    list(
        package_required = as.list(publicationLockSort(setdiff(
            required,
            standard
        ))),
        package_declared_suggests = as.list(publicationLockSort(setdiff(
            suggests,
            standard
        ))),
        publication_tooling = as.list(c("jsonvalidate", "renv"))
    )
}

publicationLockRoots <- function(root_groups = publicationLockRootGroups()) {
    publicationLockSort(unlist(root_groups, use.names = FALSE))
}

publicationLockRecordDependencies <- function(record) {
    fields <- intersect(publicationLockRequiredFields(), names(record))
    publicationLockParseField(record[fields])
}

publicationLockMetadataFields <- function() {
    c(publicationLockRequiredFields(), "Suggests", "Enhances")
}

publicationLockMetadataValid <- function(lock) {
    all(vapply(lock$Packages, \(record) {
        fields <- intersect(publicationLockMetadataFields(), names(record))
        dependencies_valid <- all(vapply(fields, \(field) {
            values <- record[[field]]
            length(values) > 0L && all(vapply(
                values,
                publicationScalarString,
                logical(1)
            ))
        }, logical(1)))
        date <- record$Date
        date_valid <- is.null(date) || (
            publicationScalarString(date) &&
                grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", date) &&
                !is.na(as.Date(date))
        )
        dependencies_valid && date_valid
    }, logical(1)))
}

publicationLockDependencyClosure <- function(lock, roots) {
    records <- lock$Packages
    standard <- publicationLockStandardPackages()
    queue <- setdiff(unique(roots), standard)
    resolved <- character()
    unresolved <- character()
    while (length(queue)) {
        package <- queue[[1L]]
        queue <- queue[-1L]
        if (package %in% resolved || package %in% unresolved) next
        if (!package %in% names(records)) {
            unresolved <- c(unresolved, package)
            next
        }
        resolved <- c(resolved, package)
        dependencies <- setdiff(
            publicationLockRecordDependencies(records[[package]]),
            standard
        )
        queue <- c(queue, setdiff(dependencies, c(resolved, queue)))
    }
    list(
        packages = publicationLockSort(resolved),
        unresolved = publicationLockSort(unresolved)
    )
}

publicationLockRequireDirectDependency <- function(lock, package, dependency) {
    record <- lock$Packages[[package]]
    valid <- !is.null(record) &&
        dependency %in% publicationLockRecordDependencies(record)
    if (!valid) {
        publicationLockAbort(paste(
            "Lock dependency path differs:",
            paste(package, dependency, sep = " -> ")
        ))
    }
    invisible(TRUE)
}

publicationLockMeasurementIsolation <- function() {
    list(
        forbidden_new_namespaces = as.list(c(
            "chromote", "gt", "jsonvalidate", "juicyjuice", "shinytest2",
            "V8", "webshot2", "websocket"
        )),
        report_dependency_path = as.list(c("gt", "juicyjuice", "V8")),
        schema_dependency_path = as.list(c("jsonvalidate", "V8")),
        measured_worker_outcome = "evidence_invalid",
        installed_but_unloaded_allowed = TRUE
    )
}

publicationLockOfficialSchema <- function() {
    path <- system.file(
        "schema",
        "draft-07.renv.lock.schema.json",
        package = "renv",
        mustWork = TRUE
    )
    list(
        renv_version = as.character(utils::packageVersion("renv")),
        package_relative_path = "schema/draft-07.renv.lock.schema.json",
        sha256 = publicationFileDigest(path),
        validator = "renv::lockfile_validate",
        engine = "ajv",
        draft = "draft-07",
        strict = TRUE
    )
}

publicationCreateLockScopeAuthority <- function(lock_path = "renv.lock") {
    lock <- publicationReadJson(lock_path)
    if (!publicationLockMetadataValid(lock)) {
        publicationLockAbort("Lock package metadata is not schema-compatible")
    }
    root_groups <- publicationLockRootGroups()
    closure <- publicationLockDependencyClosure(
        lock,
        publicationLockRoots(root_groups)
    )
    if (length(closure$unresolved)) {
        publicationLockAbort(paste(
            "Lock has unresolved declared dependencies:",
            paste(closure$unresolved, collapse = ", ")
        ))
    }
    package_names <- names(lock$Packages)
    if (!identical(package_names, publicationLockSort(package_names))) {
        publicationLockAbort("Lock package records are not ordered")
    }
    if (!identical(package_names, closure$packages)) {
        publicationLockAbort("Lock differs from its declared dependency closure")
    }
    revisions <- publicationComparatorRevisions()
    descriptions <- lapply(revisions, publicationLockGitDescription)
    names(descriptions) <- names(revisions)
    list(
        schema = "multischolar.omics_publication_lock_scope",
        schema_version = "1.0.0",
        scope_authority_id =
            "multischolar.omics_publication_lock_scope.2026-08-27.v1",
        owner_ticket_id = "OMICS-ART-064",
        status = "resolved",
        revisions = as.list(revisions),
        descriptions = descriptions,
        root_groups = root_groups,
        resolved = list(
            package_count = as.integer(length(package_names)),
            package_names = as.list(package_names),
            package_names_sha256 = publicationObjectDigest(as.list(package_names)),
            package_records_sha256 = publicationObjectDigest(lock$Packages)
        ),
        lockfile = list(
            path = lock_path,
            sha256 = publicationFileDigest(lock_path)
        ),
        official_schema = publicationLockOfficialSchema(),
        measurement_isolation = publicationLockMeasurementIsolation(),
        policy = list(
            every_declared_suggest_required = TRUE,
            unresolved_root_allowed = FALSE,
            ambient_package_allowed = FALSE,
            undeclared_extra_record_allowed = FALSE,
            package_sources = as.list(c(
                "Repository", "Bioconductor", "GitHub"
            ))
        )
    )
}

publicationValidateLockScopeAuthority <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "scope_authority_id", "owner_ticket_id",
        "status", "revisions", "descriptions", "root_groups", "resolved",
        "lockfile", "official_schema", "measurement_isolation", "policy"
    ), "Lock scope authority")
    expected <- publicationCreateLockScopeAuthority(record$lockfile$path)
    if (!identical(
        publicationObjectDigest(record),
        publicationObjectDigest(expected)
    )) {
        publicationLockAbort("Lock scope authority differs")
    }
    invisible(record)
}

publicationValidateOfficialLockSchema <- function(lock_path, scope) {
    publicationValidateLockScopeAuthority(scope)
    if (!requireNamespace("jsonvalidate", quietly = TRUE)) {
        publicationLockAbort("jsonvalidate is unavailable for official schema validation")
    }
    schema_path <- system.file(
        "schema",
        "draft-07.renv.lock.schema.json",
        package = "renv",
        mustWork = TRUE
    )
    if (!identical(
        publicationFileDigest(schema_path),
        scope$official_schema$sha256
    )) {
        publicationLockAbort("Official renv schema digest differs")
    }
    valid <- renv::lockfile_validate(
        lockfile = lock_path,
        schema = schema_path,
        error = FALSE,
        strict = TRUE
    )
    if (!isTRUE(valid)) {
        publicationLockAbort("renv.lock fails the official strict schema")
    }
    invisible(scope)
}

publicationValidateMeasuredNamespaces <- function(before, after, scope) {
    forbidden <- unlist(
        scope$measurement_isolation$forbidden_new_namespaces,
        use.names = FALSE
    )
    newly_loaded <- setdiff(after, before)
    violations <- publicationLockSort(intersect(newly_loaded, forbidden))
    if (length(violations)) {
        publicationLockAbort(paste(
            "Measured worker loaded verification namespaces:",
            paste(violations, collapse = ", ")
        ))
    }
    list(
        before = as.list(publicationLockSort(before)),
        after = as.list(publicationLockSort(after)),
        newly_loaded = as.list(publicationLockSort(newly_loaded)),
        forbidden_loaded = list(),
        valid = TRUE
    )
}
