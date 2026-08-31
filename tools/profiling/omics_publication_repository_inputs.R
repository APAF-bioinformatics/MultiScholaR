publicationRepositoryInputPath <- function(root = publicationComparatorRoot()) {
    file.path(root, "build-inputs", "repository-authority")
}

publicationExternalInputById <- function(repositories, input_id) {
    matches <- vapply(
        repositories$external_build_inputs,
        \(input) identical(input$input_id, input_id),
        logical(1)
    )
    if (sum(matches) != 1L) {
        publicationBuildAbort("External build input id is not unique")
    }
    repositories$external_build_inputs[[which(matches)]]
}

publicationPrepareRenvBootstrapInput <- function(
    repositories,
    root = publicationComparatorRoot()
) {
    publicationValidateRepositoryAuthority(repositories)
    input <- publicationExternalInputById(
        repositories,
        "renv-bootstrap-1.2.3"
    )
    input_root <- file.path(root, "build-inputs", input$input_id)
    if (file.exists(input_root) || dir.exists(input_root)) {
        publicationBuildAbort("renv bootstrap input root already exists")
    }
    dir.create(input_root, recursive = TRUE, showWarnings = FALSE)
    downloaded <- publicationDownloadRepositoryInput(
        input$input_id,
        input$url,
        input$sha256,
        input_root,
        file.path(input_root, "logs")
    )
    if (as.numeric(downloaded$size_bytes) != as.numeric(input$size_bytes)) {
        publicationBuildAbort("renv bootstrap input size differs")
    }
    list(
        schema = "multischolar.omics_publication_renv_bootstrap_input",
        schema_version = "1.0.0",
        status = "verified",
        input_id = input$input_id,
        root = input_root,
        archive = downloaded,
        publication_authority = FALSE
    )
}

publicationValidateRenvBootstrapInput <- function(
    receipt,
    repositories
) {
    publicationValidateRepositoryAuthority(repositories)
    expected <- publicationExternalInputById(
        repositories,
        "renv-bootstrap-1.2.3"
    )
    path <- file.path(receipt$root, receipt$archive$path)
    valid <- identical(
        receipt$schema,
        "multischolar.omics_publication_renv_bootstrap_input"
    ) && identical(receipt$schema_version, "1.0.0") &&
        identical(receipt$status, "verified") &&
        identical(receipt$input_id, expected$input_id) &&
        file.exists(path) &&
        identical(publicationFileDigest(path), expected$sha256) &&
        as.numeric(file.info(path)$size) == as.numeric(expected$size_bytes) &&
        !isTRUE(receipt$publication_authority)
    if (!valid) publicationBuildAbort("renv bootstrap input differs")
    invisible(receipt)
}

publicationPrepareArrowBuildInput <- function(
    repositories,
    root = publicationComparatorRoot()
) {
    publicationValidateRepositoryAuthority(repositories)
    input <- publicationExternalInputById(
        repositories,
        "libarrow-linux-x86_64-24.0.0"
    )
    input_root <- file.path(root, "build-inputs", input$input_id)
    if (file.exists(input_root) || dir.exists(input_root)) {
        publicationBuildAbort("Arrow build input root already exists")
    }
    dir.create(input_root, recursive = TRUE, showWarnings = FALSE)
    downloaded <- publicationDownloadRepositoryInput(
        input$input_id,
        input$url,
        input$sha256,
        input_root,
        file.path(input_root, "logs")
    )
    path <- file.path(input_root, downloaded$path)
    sha512 <- digest::digest(file = path, algo = "sha512")
    valid <- as.numeric(downloaded$size_bytes) == as.numeric(input$size_bytes) &&
        identical(sha512, input$sha512)
    if (!valid) publicationBuildAbort("Arrow build input differs")
    list(
        schema = "multischolar.omics_publication_arrow_build_input",
        schema_version = "1.0.0",
        status = "verified",
        input_id = input$input_id,
        root = input_root,
        archive = downloaded,
        sha512 = sha512,
        publication_authority = FALSE
    )
}

publicationValidateArrowBuildInput <- function(receipt, repositories) {
    publicationValidateRepositoryAuthority(repositories)
    expected <- publicationExternalInputById(
        repositories,
        "libarrow-linux-x86_64-24.0.0"
    )
    path <- file.path(receipt$root, receipt$archive$path)
    valid <- identical(
        receipt$schema,
        "multischolar.omics_publication_arrow_build_input"
    ) && identical(receipt$schema_version, "1.0.0") &&
        identical(receipt$status, "verified") &&
        identical(receipt$input_id, expected$input_id) &&
        file.exists(path) &&
        identical(publicationFileDigest(path), expected$sha256) &&
        identical(digest::digest(file = path, algo = "sha512"), expected$sha512) &&
        as.numeric(file.info(path)$size) == as.numeric(expected$size_bytes) &&
        !isTRUE(receipt$publication_authority)
    if (!valid) publicationBuildAbort("Arrow build input receipt differs")
    invisible(receipt)
}

publicationDownloadRepositoryInput <- function(
    id,
    url,
    expected_sha256,
    root,
    logs
) {
    path <- file.path(root, paste0(id, "-", basename(url)))
    command <- publicationBuildRun(
        "curl",
        c("-fsSL", "-o", path, url),
        root,
        logs,
        timeout_seconds = 600
    )
    publicationBuildRequireSuccess(command, paste(id, "download"))
    observed <- publicationFileDigest(path)
    if (!identical(observed, expected_sha256)) {
        publicationBuildAbort(paste(id, "download digest differs"))
    }
    list(
        input_id = id,
        url = url,
        path = basename(path),
        size_bytes = as.numeric(file.info(path)$size),
        sha256 = observed,
        command = command
    )
}

publicationVerifyRepositoryInputs <- function(
    repositories,
    root = publicationComparatorRoot()
) {
    publicationValidateRepositoryAuthority(repositories)
    input_root <- publicationRepositoryInputPath(root)
    if (file.exists(input_root) || dir.exists(input_root)) {
        publicationBuildAbort("Repository input root already exists")
    }
    dir.create(input_root, recursive = TRUE, showWarnings = FALSE)
    logs <- file.path(input_root, "logs")
    indexes <- lapply(repositories$repositories, \(repository) {
        publicationDownloadRepositoryInput(
            repository$repository_id,
            paste0(repository$url, "/", repository$index_path),
            repository$index_sha256,
            input_root,
            logs
        )
    })
    remotes <- lapply(repositories$remotes, \(remote) {
        publicationDownloadRepositoryInput(
            remote$package,
            remote$archive_url,
            remote$archive_sha256,
            input_root,
            logs
        )
    })
    list(
        schema = "multischolar.omics_publication_repository_inputs",
        schema_version = "1.0.0",
        authority_id = repositories$repository_authority_id,
        status = "verified",
        indexes = indexes,
        remotes = remotes,
        verification_root = input_root,
        publication_authority = FALSE
    )
}

publicationValidateRepositoryInputReceipt <- function(
    receipt,
    repositories,
    root = publicationComparatorRoot()
) {
    publicationValidateRepositoryAuthority(repositories)
    expected_root <- publicationRepositoryInputPath(root)
    inputs <- c(receipt$indexes, receipt$remotes)
    valid_files <- all(vapply(inputs, \(input) {
        path <- file.path(expected_root, input$path)
        file.exists(path) &&
            identical(publicationFileDigest(path), input$sha256) &&
            as.numeric(file.info(path)$size) == as.numeric(input$size_bytes)
    }, logical(1)))
    valid <- identical(
        receipt$schema,
        "multischolar.omics_publication_repository_inputs"
    ) && identical(receipt$schema_version, "1.0.0") &&
        identical(receipt$authority_id, repositories$repository_authority_id) &&
        identical(receipt$status, "verified") &&
        identical(normalizePath(receipt$verification_root), expected_root) &&
        length(receipt$indexes) == length(repositories$repositories) &&
        length(receipt$remotes) == length(repositories$remotes) &&
        valid_files && !isTRUE(receipt$publication_authority)
    if (!valid) publicationBuildAbort("Repository input receipt differs")
    invisible(receipt)
}
