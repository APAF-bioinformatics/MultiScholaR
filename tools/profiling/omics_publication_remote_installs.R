publicationRemoteInstalledPackage <- function(remote_id) {
    packages <- c(RUVIIIC = "RUVIIIC", GlimmaV2 = "Glimma")
    package <- unname(packages[[remote_id]])
    if (is.null(package)) {
        publicationBuildAbort("Remote installed package id is unknown")
    }
    package
}

publicationRemoteArchiveReceipt <- function(repository_inputs, remote_id) {
    matches <- vapply(
        repository_inputs$remotes,
        \(remote) identical(remote$input_id, remote_id),
        logical(1)
    )
    if (sum(matches) != 1L) {
        publicationBuildAbort("Remote archive receipt is not unique")
    }
    repository_inputs$remotes[[which(matches)]]
}

publicationTargetLibraryEnvironment <- function(environment, paths) {
    environment[["R_LIBS"]] <- paths$library
    environment[["R_LIBS_USER"]] <- paths$library
    environment[["R_LIBS_SITE"]] <- paths$site_library
    environment
}

publicationInstallTargetRenv <- function(
    paths,
    repositories,
    bootstrap_input,
    environment
) {
    publicationValidateRenvBootstrapInput(bootstrap_input, repositories)
    archive <- file.path(bootstrap_input$root, bootstrap_input$archive$path)
    target_environment <- publicationTargetLibraryEnvironment(
        environment,
        paths
    )
    command <- publicationBuildRun(
        file.path(R.home("bin"), "R"),
        c(
            "CMD", "INSTALL",
            paste0("--library=", paths$library),
            archive
        ),
        paths$root,
        paths$logs,
        target_environment,
        timeout_seconds = 600
    )
    publicationBuildRequireSuccess(command, "Target renv install")
    description <- read.dcf(file.path(
        paths$library,
        "renv",
        "DESCRIPTION"
    ))
    version <- unname(description[1L, "Version"])
    valid <- identical(version, "1.2.3") &&
        identical(publicationFileDigest(archive), bootstrap_input$archive$sha256)
    if (!valid) publicationBuildAbort("Target renv install differs")
    list(
        input_id = bootstrap_input$input_id,
        version = version,
        archive_sha256 = bootstrap_input$archive$sha256,
        command = command
    )
}

publicationInstallGovernedRemotes <- function(
    paths,
    repositories,
    repository_inputs,
    environment,
    root = publicationComparatorRoot()
) {
    publicationValidateRepositoryInputReceipt(
        repository_inputs,
        repositories,
        root
    )
    target_environment <- publicationTargetLibraryEnvironment(
        environment,
        paths
    )
    commands <- lapply(repositories$remotes, \(remote) {
        archive_receipt <- publicationRemoteArchiveReceipt(
            repository_inputs,
            remote$package
        )
        archive <- file.path(
            repository_inputs$verification_root,
            archive_receipt$path
        )
        if (!identical(publicationFileDigest(archive), remote$archive_sha256)) {
            publicationBuildAbort("Governed remote archive digest differs")
        }
        command <- publicationBuildRun(
            file.path(R.home("bin"), "R"),
            c(
                "CMD", "INSTALL",
                paste0("--library=", paths$library),
                archive
            ),
            paths$root,
            paths$logs,
            target_environment,
            timeout_seconds = 1800
        )
        publicationBuildRequireSuccess(
            command,
            paste(remote$package, "local source install")
        )
        package <- publicationRemoteInstalledPackage(remote$package)
        description <- read.dcf(file.path(
            paths$library,
            package,
            "DESCRIPTION"
        ))
        list(
            remote_id = remote$package,
            installed_package = package,
            version = unname(description[1L, "Version"]),
            revision = remote$commit,
            archive_sha256 = remote$archive_sha256,
            command = command
        )
    })
    names(commands) <- vapply(
        repositories$remotes,
        `[[`,
        character(1),
        "package"
    )
    commands
}

publicationValidateGovernedRemoteInstalls <- function(
    installs,
    repositories,
    lock
) {
    valid <- length(installs) == length(repositories$remotes)
    for (remote in repositories$remotes) {
        install <- installs[[remote$package]]
        package <- publicationRemoteInstalledPackage(remote$package)
        valid <- valid && !is.null(install) &&
            identical(install$installed_package, package) &&
            identical(install$version, lock$Packages[[package]]$Version) &&
            identical(install$revision, remote$commit) &&
            identical(install$archive_sha256, remote$archive_sha256) &&
            identical(install$command$exit_status, 0L) &&
            !isTRUE(install$command$timed_out)
    }
    if (!valid) publicationBuildAbort("Governed remote installs differ")
    invisible(installs)
}
