publicationBuildAbort <- function(message) {
    publicationComparatorAbort(message, "comparator_build_error")
}

publicationPacmanInventory <- function(arguments) {
    result <- processx::run(
        "pacman",
        arguments,
        error_on_status = FALSE,
        echo = FALSE
    )
    if (result$status != 0L) {
        publicationBuildAbort("Live system package inventory is unavailable")
    }
    lines <- strsplit(trimws(result$stdout), "\n", fixed = TRUE)[[1L]]
    list(
        count = as.integer(length(lines)),
        sha256 = digest::digest(
            result$stdout,
            algo = "sha256",
            serialize = FALSE
        )
    )
}

publicationValidateSystemDependencyAuthority <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "system_authority_id", "owner_ticket_id",
        "captured_at", "status", "host_scope", "package_manager",
        "selected_packages", "missing_optional_packages", "commands",
        "native_libraries", "r_build", "policy"
    ), "System dependency authority")
    packages <- vapply(record$selected_packages, `[[`, character(1), "name")
    commands <- vapply(record$commands, `[[`, character(1), "path")
    libraries <- vapply(record$native_libraries, `[[`, character(1), "path")
    names_only <- publicationPacmanInventory("-Qq")
    names_versions <- publicationPacmanInventory("-Q")
    command_valid <- !anyDuplicated(commands) && all(vapply(
        record$commands,
        \(command) publicationComparatorShaValid(command$sha256) &&
            file.exists(command$path) &&
            identical(publicationFileDigest(command$path), command$sha256),
        logical(1)
    ))
    library_valid <- !anyDuplicated(libraries) && all(vapply(
        record$native_libraries,
        \(library) publicationComparatorShaValid(library$sha256) &&
            file.exists(library$path) &&
            identical(publicationFileDigest(library$path), library$sha256),
        logical(1)
    ))
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_system_dependencies"
    ) && identical(record$schema_version, "1.0.0") &&
        identical(record$owner_ticket_id, "OMICS-ART-064") &&
        identical(record$status, "host_scoped_frozen") &&
        identical(record$host_scope$kernel, unname(Sys.info()[["release"]])) &&
        identical(record$host_scope$machine, unname(Sys.info()[["machine"]])) &&
        identical(record$host_scope$r_version, "4.6.1") &&
        identical(record$host_scope$r_platform, R.version$platform) &&
        identical(
            record$package_manager$installed_package_count,
            names_versions$count
        ) && identical(
            record$package_manager$name_inventory_sha256,
            names_only$sha256
        ) && identical(
            record$package_manager$name_version_inventory_sha256,
            names_versions$sha256
        ) && !anyDuplicated(packages) && command_valid && library_valid &&
        identical(record$policy$host_generalization, FALSE) &&
        isTRUE(record$policy$exact_inventory_required) &&
        isTRUE(record$policy$ambient_system_drift_invalidates_build)
    if (!valid) publicationBuildAbort("System dependency authority differs")
    invisible(record)
}

publicationBuildEnvironment <- function(library = NULL, site_library = NULL) {
    values <- c(
        TZ = "UTC",
        LANG = "C.UTF-8",
        LC_ALL = "C.UTF-8",
        PATH = "/usr/local/sbin:/usr/local/bin:/usr/bin:/bin",
        SHELL = "/bin/bash",
        GITHUB_PAT = "",
        GITHUB_TOKEN = "",
        R_ENVIRON_USER = "",
        R_PROFILE_USER = "",
        OMP_NUM_THREADS = "1",
        OPENBLAS_NUM_THREADS = "1",
        MKL_NUM_THREADS = "1",
        ARROW_NUM_THREADS = "1",
        DUCKDB_THREADS = "1",
        CMAKE_BUILD_PARALLEL_LEVEL = "1",
        MAKEFLAGS = "-j1",
        RENV_CONFIG_SANDBOX_ENABLED = "TRUE",
        RENV_CONFIG_AUTO_SNAPSHOT = "FALSE",
        RENV_CONFIG_SYNCHRONIZED_CHECK = "FALSE"
    )
    if (!is.null(library)) {
        if (is.null(site_library)) {
            publicationBuildAbort("An isolated site library is required")
        }
        values <- c(
            values,
            R_LIBS = library,
            R_LIBS_USER = library,
            R_LIBS_SITE = site_library
        )
    }
    values
}

publicationBuildRun <- function(
    command,
    arguments,
    workdir,
    log_dir,
    environment = publicationBuildEnvironment(),
    timeout_seconds = 3600
) {
    dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
    command_id <- sprintf("command-%03d", length(list.files(log_dir)) / 2 + 1L)
    stdout_path <- file.path(log_dir, paste0(command_id, "-stdout.log"))
    stderr_path <- file.path(log_dir, paste0(command_id, "-stderr.log"))
    started <- proc.time()[["elapsed"]]
    result <- processx::run(
        command,
        arguments,
        wd = workdir,
        env = environment,
        stdout = stdout_path,
        stderr = stderr_path,
        timeout = timeout_seconds,
        error_on_status = FALSE,
        echo = FALSE
    )
    list(
        command = command,
        arguments = as.list(arguments),
        workdir = normalizePath(workdir, mustWork = TRUE),
        environment = as.list(environment),
        exit_status = result$status,
        timed_out = isTRUE(result$timeout),
        elapsed_seconds = proc.time()[["elapsed"]] - started,
        stdout_sha256 = publicationFileDigest(stdout_path),
        stderr_sha256 = publicationFileDigest(stderr_path),
        stdout_path = stdout_path,
        stderr_path = stderr_path
    )
}

publicationBuildRequireSuccess <- function(receipt, label) {
    if (!identical(receipt$exit_status, 0L) || isTRUE(receipt$timed_out)) {
        publicationBuildAbort(paste(label, "failed"))
    }
    invisible(receipt)
}

publicationComparatorRoot <- function() {
    publicationPath(".omics-publication-comparators")
}

publicationComparatorPathInsideRoot <- function(path, root) {
    normalized_root <- normalizePath(root, mustWork = TRUE)
    normalized <- normalizePath(path, mustWork = FALSE)
    startsWith(
        paste0(normalized, .Platform$file.sep),
        paste0(normalized_root, .Platform$file.sep)
    )
}

publicationCreateComparatorWorktree <- function(
    comparator_id,
    revision,
    root = publicationComparatorRoot()
) {
    expected <- publicationComparatorRevisions()[[comparator_id]]
    if (is.null(expected) || !identical(revision, expected)) {
        publicationBuildAbort("Comparator worktree revision is not governed")
    }
    dir.create(root, recursive = TRUE, showWarnings = FALSE)
    path <- file.path(root, "worktrees", comparator_id)
    if (file.exists(path) || dir.exists(path)) {
        publicationBuildAbort("Comparator worktree path already exists")
    }
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    created <- publicationComparatorGit(c(
        "worktree", "add", "--detach", path, revision
    ))
    if (created$status != 0L) {
        publicationBuildAbort("Comparator detached worktree creation failed")
    }
    state <- publicationComparatorWorktreeState(path)
    record <- list(
        comparator_id = comparator_id,
        path = path,
        revision = revision,
        source_identity = publicationComparatorSourceIdentity(revision),
        detached = state$detached,
        clean = state$clean,
        status_sha256 = state$status_sha256
    )
    publicationValidateComparatorWorktree(record)
    record
}

publicationRemoveComparatorWorktree <- function(
    record,
    root = publicationComparatorRoot()
) {
    publicationValidateComparatorWorktree(record)
    if (!publicationComparatorPathInsideRoot(record$path, root)) {
        publicationBuildAbort("Comparator worktree is outside the owned root")
    }
    removed <- publicationComparatorGit(c(
        "worktree", "remove", "--force", record$path
    ))
    if (removed$status != 0L || dir.exists(record$path)) {
        publicationBuildAbort("Comparator worktree cleanup failed")
    }
    invisible(TRUE)
}

publicationExportComparatorSource <- function(record, staging_root) {
    publicationValidateComparatorWorktree(record)
    path <- file.path(staging_root, record$comparator_id)
    archive <- file.path(staging_root, paste0(record$comparator_id, ".tar"))
    if (dir.exists(path) || file.exists(archive)) {
        publicationBuildAbort("Comparator staging path already exists")
    }
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    exported <- publicationComparatorGit(c(
        "archive", "--format=tar", paste0("--output=", archive),
        record$revision
    ))
    if (exported$status != 0L) {
        publicationBuildAbort("Comparator source export failed")
    }
    extracted <- processx::run(
        "tar",
        c("-xf", archive, "-C", path),
        error_on_status = FALSE,
        echo = FALSE
    )
    if (extracted$status != 0L) {
        publicationBuildAbort("Comparator source extraction failed")
    }
    publicationValidateComparatorWorktree(record)
    list(
        path = path,
        archive_path = archive,
        archive_sha256 = publicationFileDigest(archive),
        description_sha256 = publicationFileDigest(file.path(
            path,
            "DESCRIPTION"
        )),
        source_identity = record$source_identity
    )
}

publicationBuildRepositories <- function(record) {
    publicationValidateRepositoryAuthority(record)
    stats::setNames(
        vapply(record$repositories, `[[`, character(1), "url"),
        vapply(record$repositories, `[[`, character(1), "repository_id")
    )
}

publicationBuildRelativePath <- function(path, root) {
    normalized_root <- normalizePath(root, mustWork = TRUE)
    normalized <- normalizePath(path, mustWork = TRUE)
    prefix <- paste0(normalized_root, .Platform$file.sep)
    if (!startsWith(paste0(normalized, .Platform$file.sep), prefix)) {
        publicationBuildAbort("Build path is outside its owned root")
    }
    if (identical(normalized, normalized_root)) return(".")
    substring(normalized, nchar(prefix) + 1L)
}

publicationValidateExternalArchiveEntries <- function(entries) {
    entries <- entries[nzchar(entries)]
    unsafe <- startsWith(entries, "/") | grepl(
        "(^|/)\\.\\.(/|$)",
        entries
    )
    valid <- length(entries) > 0L && !any(unsafe) &&
        all(startsWith(entries, "v8/")) &&
        !anyDuplicated(entries)
    if (!valid) publicationBuildAbort("External archive entries are unsafe")
    invisible(entries)
}

publicationPrepareExternalBuildInput <- function(
    repositories,
    root = publicationComparatorRoot()
) {
    publicationValidateRepositoryAuthority(repositories)
    input <- repositories$external_build_inputs[[1L]]
    dir.create(root, recursive = TRUE, showWarnings = FALSE)
    input_root <- file.path(root, "build-inputs", input$input_id)
    archive <- file.path(input_root, basename(input$url))
    extracted <- file.path(input_root, "extracted")
    logs <- file.path(input_root, "logs")
    if (dir.exists(extracted)) {
        publicationBuildAbort("External build input is already extracted")
    }
    dir.create(input_root, recursive = TRUE, showWarnings = FALSE)
    download <- NULL
    if (!file.exists(archive)) {
        download <- publicationBuildRun(
            "curl",
            c("-fsSL", "-o", archive, input$url),
            input_root,
            logs,
            timeout_seconds = 600
        )
        publicationBuildRequireSuccess(download, "External build input download")
    }
    valid_archive <- identical(publicationFileDigest(archive), input$sha256) &&
        as.numeric(file.info(archive)$size) == as.numeric(input$size_bytes)
    if (!valid_archive) {
        publicationBuildAbort("External build input digest or size differs")
    }
    listing <- publicationBuildRun(
        "tar",
        c("-tzf", archive),
        input_root,
        logs,
        timeout_seconds = 60
    )
    publicationBuildRequireSuccess(listing, "External build input listing")
    entries <- readLines(listing$stdout_path, warn = FALSE)
    publicationValidateExternalArchiveEntries(entries)
    staging <- file.path(input_root, "extracting")
    dir.create(staging, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(staging, recursive = TRUE, force = TRUE), add = TRUE)
    extraction <- publicationBuildRun(
        "tar",
        c("-xzf", archive, "-C", staging),
        input_root,
        logs,
        timeout_seconds = 120
    )
    publicationBuildRequireSuccess(extraction, "External build input extraction")
    source <- file.path(staging, "v8")
    required <- c(
        file.path(source, "include", "v8.h"),
        file.path(source, "lib", "libv8_monolith.a")
    )
    if (!all(file.exists(required)) || !file.rename(source, extracted)) {
        publicationBuildAbort("External build input contents differ")
    }
    inventory <- publicationInstalledInventory(extracted)
    list(
        input_id = input$input_id,
        root = extracted,
        archive = list(
            path = publicationBuildRelativePath(archive, root),
            size_bytes = input$size_bytes,
            sha256 = input$sha256
        ),
        download = download,
        listing = listing,
        extraction = extraction,
        extracted_file_count = inventory$file_count,
        extracted_inventory_sha256 = inventory$inventory_sha256,
        include_dir = file.path(extracted, "include"),
        library_dir = file.path(extracted, "lib"),
        verified_before_extraction = TRUE
    )
}

publicationExternalBuildEnvironment <- function(input) {
    c(
        DISABLE_STATIC_LIBV8 = "1",
        DOWNLOAD_STATIC_LIBV8 = "",
        V8_PKG_CFLAGS = paste0("-I", input$include_dir),
        V8_PKG_LIBS = paste(
            paste0("-L", input$library_dir),
            "-lv8_monolith"
        )
    )
}

publicationValidatePreparedExternalBuildInput <- function(
    input,
    repositories,
    root = publicationComparatorRoot()
) {
    publicationValidateRepositoryAuthority(repositories)
    expected <- repositories$external_build_inputs[[1L]]
    archive <- file.path(root, input$archive$path)
    inventory <- publicationInstalledInventory(input$root)
    valid <- identical(input$input_id, expected$input_id) &&
        file.exists(archive) && dir.exists(input$root) &&
        identical(publicationFileDigest(archive), expected$sha256) &&
        as.numeric(file.info(archive)$size) == as.numeric(expected$size_bytes) &&
        identical(input$archive$sha256, expected$sha256) &&
        identical(
            input$extracted_inventory_sha256,
            inventory$inventory_sha256
        ) && file.exists(file.path(input$include_dir, "v8.h")) &&
        file.exists(file.path(input$library_dir, "libv8_monolith.a")) &&
        isTRUE(input$verified_before_extraction)
    if (!valid) publicationBuildAbort("Prepared external build input differs")
    invisible(input)
}

publicationBuildRLiteral <- function(value) {
    paste(deparse(value, width.cutoff = 500L), collapse = "")
}

publicationBuildWriteScript <- function(lines, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines(lines, path, useBytes = TRUE)
    list(
        path = path,
        sha256 = publicationFileDigest(path),
        size_bytes = as.numeric(file.info(path)$size)
    )
}

publicationRestorePaths <- function(root, restore_id) {
    restore_root <- file.path(root, "restores", restore_id)
    list(
        root = restore_root,
        bootstrap_library = file.path(restore_root, "bootstrap-library"),
        library = file.path(restore_root, "library"),
        site_library = file.path(restore_root, "empty-site-library"),
        cache = file.path(restore_root, "renv-cache"),
        source_cache = file.path(restore_root, "source-cache"),
        state = file.path(restore_root, "renv-state"),
        home = file.path(restore_root, "home"),
        temp = file.path(restore_root, "tmp"),
        project = file.path(restore_root, "project"),
        scripts = file.path(restore_root, "scripts"),
        logs = file.path(restore_root, "logs"),
        audit = file.path(restore_root, "restore-audit.json")
    )
}

publicationInitializeRestorePaths <- function(paths) {
    if (file.exists(paths$root) || dir.exists(paths$root)) {
        publicationBuildAbort("Restore root already exists")
    }
    dir.create(paths$root, recursive = TRUE, showWarnings = FALSE)
    directories <- paths[setdiff(names(paths), c("root", "audit"))]
    invisible(lapply(directories, dir.create, recursive = TRUE))
}

publicationBuildSourceDateEpoch <- function() {
    result <- publicationComparatorGit(c(
        "show", "-s", "--format=%ct",
        publicationComparatorRevisions()[["pre_repair_performance"]]
    ))
    epoch <- trimws(result$stdout)
    if (result$status != 0L || !grepl("^[0-9]+$", epoch)) {
        publicationBuildAbort("Comparator source epoch is unavailable")
    }
    epoch
}

publicationRestoreEnvironment <- function(paths, external_input, arrow_input) {
    arrow_archive <- file.path(arrow_input$root, arrow_input$archive$path)
    c(
        publicationBuildEnvironment(
            paths$bootstrap_library,
            paths$site_library
        ),
        publicationExternalBuildEnvironment(external_input),
        ARROW_DOWNLOADED_BINARIES = arrow_archive,
        ARROW_OFFLINE_BUILD = "true",
        ARROW_R_ENFORCE_CHECKSUM = "true",
        LIBARROW_BINARY = "linux-x86_64",
        LIBARROW_BUILD = "false",
        SOURCE_DATE_EPOCH = publicationBuildSourceDateEpoch(),
        HOME = paths$home,
        TMPDIR = paths$temp,
        RENV_PATHS_ROOT = paths$state,
        RENV_PATHS_CACHE = paths$cache,
        RENV_PATHS_SOURCE = paths$source_cache,
        RENV_CONFIG_CACHE_ENABLED = "FALSE",
        RENV_CONFIG_INSTALL_JOBS = "1",
        RENV_CONFIG_INSTALL_STAGED = "TRUE"
    )
}

publicationBootstrapRenv <- function(
    paths,
    repositories,
    bootstrap_input,
    environment
) {
    publicationValidateRenvBootstrapInput(bootstrap_input, repositories)
    archive <- file.path(bootstrap_input$root, bootstrap_input$archive$path)
    command <- publicationBuildRun(
        file.path(R.home("bin"), "R"),
        c(
            "CMD", "INSTALL",
            paste0("--library=", paths$bootstrap_library),
            archive
        ),
        paths$root,
        paths$logs,
        environment,
        timeout_seconds = 600
    )
    publicationBuildRequireSuccess(command, "renv bootstrap")
    description <- read.dcf(file.path(
        paths$bootstrap_library,
        "renv",
        "DESCRIPTION"
    ))
    version <- unname(description[1L, "Version"])
    if (!identical(version, "1.2.3")) {
        publicationBuildAbort("renv bootstrap version differs")
    }
    list(
        input_id = bootstrap_input$input_id,
        archive_sha256 = bootstrap_input$archive$sha256,
        version = version,
        command = command
    )
}

publicationRestoreCommonLock <- function(
    paths,
    repositories,
    environment,
    lock_path = publicationPath("renv.lock")
) {
    copied <- file.copy(lock_path, file.path(paths$project, "renv.lock"))
    if (!isTRUE(copied)) publicationBuildAbort("Restore lock copy failed")
    repos <- publicationBuildRepositories(repositories)
    script_path <- file.path(paths$scripts, "restore-common-lock.R")
    script <- publicationBuildWriteScript(c(
        paste0("repos <- ", publicationBuildRLiteral(repos)),
        paste0("library_path <- ", publicationBuildRLiteral(paths$library)),
        paste0("project_path <- ", publicationBuildRLiteral(paths$project)),
        paste0(
            "bootstrap <- ",
            publicationBuildRLiteral(paths$bootstrap_library)
        ),
        ".libPaths(c(bootstrap, .Library))",
        "library(renv, lib.loc = bootstrap)",
        paste0(
            "options(repos = repos, pkgType = \"source\", Ncpus = 1L, ",
            "renv.install.timeout = 18000L)"
        ),
        "stopifnot(renv::config$install.jobs() == 1L)",
        "stopifnot(getOption(\"renv.install.timeout\") == 18000L)",
        paste0(
            "renv::restore(project = project_path, library = library_path, ",
            "lockfile = file.path(project_path, \"renv.lock\"), ",
            "rebuild = TRUE, repos = repos, clean = TRUE, strict = TRUE, ",
            "exclude = c(\"renv\", \"RUVIIIC\", \"Glimma\"), ",
            "transactional = FALSE, prompt = FALSE)"
        )
    ), script_path)
    command <- publicationBuildRun(
        file.path(R.home("bin"), "Rscript"),
        c("--vanilla", script_path),
        paths$root,
        paths$logs,
        environment,
        timeout_seconds = 21600
    )
    publicationBuildRequireSuccess(command, "Common lock restore")
    list(script = script, command = command)
}

publicationPackageTreeSummary <- function(path) {
    inventory <- publicationInstalledInventory(path)
    list(
        file_count = inventory$file_count,
        inventory_sha256 = inventory$inventory_sha256
    )
}

publicationRestoredLibraryInventory <- function(library, lock) {
    root <- normalizePath(library, mustWork = TRUE)
    packages <- lapply(names(lock$Packages), \(package) {
        path <- file.path(root, package)
        description_path <- file.path(path, "DESCRIPTION")
        if (!file.exists(description_path)) {
            publicationBuildAbort(paste("Restored package is missing:", package))
        }
        description <- read.dcf(description_path)
        version <- unname(description[1L, "Version"])
        expected <- lock$Packages[[package]]$Version
        if (!identical(version, expected)) {
            publicationBuildAbort(paste("Restored package version differs:", package))
        }
        list(
            package = package,
            version = version,
            source = lock$Packages[[package]]$Source,
            tree = publicationPackageTreeSummary(path),
            native_libraries = publicationNativeLibraryInventory(path)
        )
    })
    names(packages) <- names(lock$Packages)
    list(
        root = root,
        package_count = as.integer(length(packages)),
        packages = packages,
        inventory_sha256 = publicationObjectDigest(packages),
        all_packages_owned = TRUE
    )
}

publicationRestoreAudit <- function(paths, scope, environment) {
    script_path <- file.path(paths$scripts, "audit-restore.R")
    script <- publicationBuildWriteScript(c(
        paste0("library_path <- ", publicationBuildRLiteral(paths$library)),
        paste0("output_path <- ", publicationBuildRLiteral(paths$audit)),
        paste0(
            "repo_root <- ",
            publicationBuildRLiteral(.PUBLICATION_REPO_ROOT)
        ),
        ".libPaths(c(library_path, .Library))",
        paste0(
            "source(file.path(repo_root, ",
            "\"tools/profiling/omics_publication_protocol.R\"))"
        ),
        paste0(
            "source(file.path(repo_root, ",
            "\"tools/profiling/omics_publication_comparators.R\"))"
        ),
        paste0(
            "source(file.path(repo_root, ",
            "\"tools/profiling/omics_publication_lock.R\"))"
        ),
        paste0(
            "scope <- publicationReadJson(file.path(repo_root, ",
            "\"tools/profiling/omics-publication-lock-scope-v1.json\"))"
        ),
        paste0(
            "publicationValidateOfficialLockSchema(",
            "file.path(repo_root, \"renv.lock\"), scope)"
        ),
        paste0(
            "package_paths <- vapply(scope$resolved$package_names, ",
            "\\(package) find.package(package, lib.loc = library_path), ",
            "character(1))"
        ),
        "owned <- startsWith(normalizePath(package_paths), paste0(",
        "    normalizePath(library_path), .Platform$file.sep))",
        "jsonlite::write_json(list(",
        "    status = if (all(owned)) \"passed\" else \"failed\",",
        "    official_schema_valid = TRUE,",
        "    library_paths = as.list(.libPaths()),",
        "    package_paths = as.list(package_paths),",
        "    all_packages_owned = all(owned),",
        "    loaded_namespaces = as.list(sort(loadedNamespaces()))",
        "), output_path, auto_unbox = TRUE, pretty = TRUE, null = \"null\")"
    ), script_path)
    command <- publicationBuildRun(
        file.path(R.home("bin"), "Rscript"),
        c("--vanilla", script_path),
        paths$root,
        paths$logs,
        environment,
        timeout_seconds = 600
    )
    publicationBuildRequireSuccess(command, "Restored library audit")
    audit <- publicationReadJson(paths$audit)
    valid <- identical(audit$status, "passed") &&
        isTRUE(audit$official_schema_valid) &&
        isTRUE(audit$all_packages_owned) &&
        length(audit$package_paths) == scope$resolved$package_count
    if (!valid) publicationBuildAbort("Restored library audit differs")
    list(script = script, command = command, evidence = audit)
}

publicationRestoreComparatorEnvironment <- function(
    restore_id,
    external_input,
    bootstrap_input,
    arrow_input,
    repository_inputs,
    root = publicationComparatorRoot()
) {
    authority <- publicationReadJson(
        "tests/testdata/omics-performance/comparators-v1.json"
    )
    publicationValidateComparatorAuthority(authority)
    repositories <- publicationReadJson(
        authority$common_environment$repositories$path
    )
    publicationValidatePreparedExternalBuildInput(
        external_input,
        repositories,
        root
    )
    publicationValidateRenvBootstrapInput(bootstrap_input, repositories)
    publicationValidateArrowBuildInput(arrow_input, repositories)
    publicationValidateRepositoryInputReceipt(
        repository_inputs,
        repositories,
        root
    )
    scope <- publicationReadJson(authority$common_environment$lock_scope$path)
    lock <- publicationReadJson(authority$common_environment$lockfile$path)
    paths <- publicationRestorePaths(root, restore_id)
    publicationInitializeRestorePaths(paths)
    environment <- publicationRestoreEnvironment(
        paths,
        external_input,
        arrow_input
    )
    bootstrap <- publicationBootstrapRenv(
        paths,
        repositories,
        bootstrap_input,
        environment
    )
    restore <- publicationRestoreCommonLock(
        paths,
        repositories,
        environment
    )
    target_renv <- publicationInstallTargetRenv(
        paths,
        repositories,
        bootstrap_input,
        environment
    )
    remote_installs <- publicationInstallGovernedRemotes(
        paths,
        repositories,
        repository_inputs,
        environment,
        root
    )
    publicationValidateGovernedRemoteInstalls(
        remote_installs,
        repositories,
        lock
    )
    audit <- publicationRestoreAudit(paths, scope, environment)
    library_inventory <- publicationRestoredLibraryInventory(
        paths$library,
        lock
    )
    list(
        schema = "multischolar.omics_publication_restore",
        schema_version = "1.0.0",
        restore_id = restore_id,
        status = "passed",
        environment = authority$common_environment,
        paths = lapply(paths, \(path) {
            if (file.exists(path) || dir.exists(path)) {
                publicationBuildRelativePath(path, root)
            } else {
                path
            }
        }),
        bootstrap = bootstrap,
        restore = restore,
        target_renv = target_renv,
        remote_installs = remote_installs,
        audit = audit,
        library_inventory = library_inventory,
        external_build_input = list(
            input_id = external_input$input_id,
            archive_sha256 = external_input$archive$sha256,
            extracted_inventory_sha256 =
                external_input$extracted_inventory_sha256
        ),
        arrow_build_input = list(
            input_id = arrow_input$input_id,
            archive_sha256 = arrow_input$archive$sha256,
            archive_sha512 = arrow_input$sha512
        ),
        ambient_user_library = FALSE,
        cache_shared_with_other_restore = FALSE,
        publication_authority = FALSE
    )
}

publicationInstalledInventory <- function(path) {
    root <- normalizePath(path, mustWork = TRUE)
    files <- list.files(
        root,
        recursive = TRUE,
        full.names = TRUE,
        all.files = TRUE,
        no.. = TRUE
    )
    files <- files[file.exists(files) & !dir.exists(files)]
    relative <- substring(
        normalizePath(files, mustWork = TRUE),
        nchar(root) + 2L
    )
    ordering <- order(relative, method = "radix")
    files <- files[ordering]
    relative <- relative[ordering]
    info <- file.info(files)
    records <- lapply(seq_along(files), \(index) {
        list(
            path = relative[[index]],
            size_bytes = as.numeric(info$size[[index]]),
            executable = as.integer(info$mode[[index]]) %% 2L == 1L,
            sha256 = publicationFileDigest(files[[index]])
        )
    })
    list(
        root = root,
        file_count = length(records),
        records = records,
        inventory_sha256 = publicationObjectDigest(records)
    )
}

publicationNativeLibraryInventory <- function(installed_package_path) {
    libraries <- list.files(
        file.path(installed_package_path, "libs"),
        pattern = "\\.(so|dll|dylib)$",
        recursive = TRUE,
        full.names = TRUE
    )
    libraries <- sort(libraries, method = "radix")
    lapply(libraries, \(path) {
        list(
            path = basename(path),
            sha256 = publicationFileDigest(path),
            size_bytes = as.numeric(file.info(path)$size)
        )
    })
}

publicationBuildPackage <- function(
    source,
    build_root,
    package_library,
    dependency_library,
    revision
) {
    dir.create(build_root, recursive = TRUE, showWarnings = FALSE)
    dir.create(package_library, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(dependency_library)) {
        publicationBuildAbort("Comparator dependency library is missing")
    }
    dependency_package <- file.path(dependency_library, "MultiScholaR")
    if (dir.exists(dependency_package)) {
        publicationBuildAbort("Dependency library contains MultiScholaR")
    }
    site_library <- file.path(build_root, "empty-site-library")
    dir.create(site_library, recursive = TRUE, showWarnings = FALSE)
    home <- file.path(build_root, "home")
    temp <- file.path(build_root, "tmp")
    dir.create(home, recursive = TRUE, showWarnings = FALSE)
    dir.create(temp, recursive = TRUE, showWarnings = FALSE)
    log_dir <- file.path(build_root, "logs")
    epoch <- publicationComparatorGit(c(
        "show", "-s", "--format=%ct", revision
    ))
    environment <- c(
        publicationBuildEnvironment(package_library, site_library),
        SOURCE_DATE_EPOCH = trimws(epoch$stdout),
        HOME = home,
        TMPDIR = temp
    )
    build <- publicationBuildRun(
        file.path(R.home("bin"), "R"),
        c("CMD", "build", "--no-build-vignettes", "--no-manual", source$path),
        build_root,
        log_dir,
        environment
    )
    publicationBuildRequireSuccess(build, "Comparator package build")
    archives <- list.files(
        build_root,
        pattern = "^MultiScholaR_.*\\.tar\\.gz$",
        full.names = TRUE
    )
    if (length(archives) != 1L) {
        publicationBuildAbort("Comparator build archive is not unique")
    }
    install <- publicationBuildRun(
        file.path(R.home("bin"), "R"),
        c(
            "CMD", "INSTALL",
            paste0("--library=", package_library),
            archives[[1L]]
        ),
        build_root,
        log_dir,
        environment
    )
    publicationBuildRequireSuccess(install, "Comparator package install")
    installed <- file.path(package_library, "MultiScholaR")
    if (!dir.exists(installed) || dir.exists(dependency_package)) {
        publicationBuildAbort("Comparator installed package is missing")
    }
    list(
        build_command = build,
        install_command = install,
        archive = list(
            path = archives[[1L]],
            sha256 = publicationFileDigest(archives[[1L]]),
            size_bytes = as.numeric(file.info(archives[[1L]])$size)
        ),
        installed_inventory = publicationInstalledInventory(installed),
        native_libraries = publicationNativeLibraryInventory(installed)
    )
}

publicationValidateComparatorBuildReceipt <- function(receipt) {
    publicationRequireNames(receipt, c(
        "schema", "schema_version", "receipt_id", "comparator_id",
        "revision", "source", "environment", "restore", "build",
        "session", "smoke", "cleanup", "status", "publication_authority"
    ), "Comparator build receipt")
    revisions <- publicationComparatorRevisions()
    valid <- identical(
        receipt$schema,
        "multischolar.omics_publication_comparator_build"
    ) && identical(receipt$schema_version, "1.0.0") &&
        identical(receipt$revision, revisions[[receipt$comparator_id]]) &&
        identical(receipt$source$source_identity$revision, receipt$revision) &&
        identical(receipt$restore$status, "passed") &&
        identical(receipt$build$build_command$exit_status, 0L) &&
        identical(receipt$build$install_command$exit_status, 0L) &&
        publicationComparatorShaValid(receipt$build$archive$sha256) &&
        publicationComparatorShaValid(
            receipt$build$installed_inventory$inventory_sha256
        ) && identical(receipt$smoke$status, "passed") &&
        isTRUE(receipt$cleanup$worktree_removed) &&
        isTRUE(receipt$cleanup$dependency_library_retained) &&
        isTRUE(receipt$cleanup$package_library_retained) &&
        identical(receipt$status, "verified") &&
        !isTRUE(receipt$publication_authority)
    if (!valid) publicationBuildAbort("Comparator build receipt differs")
    invisible(receipt)
}

publicationCompareBuildReceipts <- function(first, second) {
    publicationValidateComparatorBuildReceipt(first)
    publicationValidateComparatorBuildReceipt(second)
    same_binding <- identical(first$comparator_id, second$comparator_id) &&
        identical(first$revision, second$revision) &&
        identical(
            publicationObjectDigest(first$source$source_identity),
            publicationObjectDigest(second$source$source_identity)
        ) && identical(
            first$environment$lockfile_sha256,
            second$environment$lockfile_sha256
        )
    same_install <- identical(
        first$build$installed_inventory$inventory_sha256,
        second$build$installed_inventory$inventory_sha256
    ) && identical(
        publicationObjectDigest(first$build$native_libraries),
        publicationObjectDigest(second$build$native_libraries)
    )
    same_science <- identical(
        first$smoke$output_sha256,
        second$smoke$output_sha256
    )
    list(
        comparator_id = first$comparator_id,
        revision = first$revision,
        source_and_environment_equal = same_binding,
        installed_and_native_equal = same_install,
        smoke_output_equal = same_science,
        archive_bytes_equal = identical(
            first$build$archive$sha256,
            second$build$archive$sha256
        ),
        reproducible = same_binding && same_install && same_science
    )
}
