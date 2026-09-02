publicationComparatorAbort <- function(message, class = "comparator_error") {
    publicationAbort(
        message,
        paste0("multischolar_publication_", class)
    )
}

publicationComparatorRevisions <- function() {
    c(
        historical_janitor =
            "c7c12851ea4d9e96f91df0e1f9e7b91c2eb51017",
        pre_repair_performance =
            "56725e90c5fda97835775f5f8f57c02703b53120"
    )
}

publicationComparatorShaValid <- function(value, length = 64L) {
    publicationScalarString(value) && grepl(
        paste0("^[0-9a-f]{", length, "}$"),
        value
    )
}

publicationComparatorGit <- function(
    arguments,
    workdir = .PUBLICATION_REPO_ROOT,
    error_on_status = FALSE
) {
    processx::run(
        "git",
        arguments,
        wd = workdir,
        error_on_status = error_on_status,
        echo = FALSE
    )
}

publicationComparatorResolveRevision <- function(revision) {
    if (!publicationComparatorShaValid(revision, 40L)) {
        publicationComparatorAbort("Comparator revision is not a full commit")
    }
    result <- publicationComparatorGit(
        c("rev-parse", "--verify", paste0(revision, "^{commit}"))
    )
    resolved <- trimws(result$stdout)
    if (result$status != 0L || !identical(resolved, revision)) {
        publicationComparatorAbort("Comparator revision does not resolve exactly")
    }
    resolved
}

publicationComparatorSourceIdentity <- function(revision) {
    revision <- publicationComparatorResolveRevision(revision)
    tree <- publicationComparatorGit(c("rev-parse", paste0(revision, "^{tree}")))
    files <- publicationComparatorGit(c(
        "ls-tree", "-r", "--name-only", revision
    ))
    archive <- tempfile("multischolar-comparator-", fileext = ".tar")
    on.exit(unlink(archive, force = TRUE), add = TRUE)
    built <- publicationComparatorGit(c(
        "archive", "--format=tar", paste0("--output=", archive), revision
    ))
    if (tree$status != 0L || files$status != 0L || built$status != 0L) {
        publicationComparatorAbort("Comparator source identity could not be read")
    }
    file_names <- strsplit(trimws(files$stdout), "\n", fixed = TRUE)[[1L]]
    list(
        revision = revision,
        tree_oid = trimws(tree$stdout),
        archive_sha256 = publicationFileDigest(archive),
        tracked_file_count = length(file_names),
        tracked_path_digest = publicationObjectDigest(as.list(file_names))
    )
}

publicationComparatorWorktreeState <- function(path) {
    if (!dir.exists(path)) {
        publicationComparatorAbort("Comparator worktree does not exist")
    }
    resolved <- normalizePath(path, mustWork = TRUE)
    head <- publicationComparatorGit(c("rev-parse", "HEAD"), resolved)
    status <- publicationComparatorGit(
        c("status", "--porcelain=v1", "--untracked-files=all"),
        resolved
    )
    branch <- publicationComparatorGit(
        c("symbolic-ref", "--quiet", "--short", "HEAD"),
        resolved
    )
    list(
        path = resolved,
        revision = trimws(head$stdout),
        detached = branch$status != 0L,
        clean = status$status == 0L && !nzchar(trimws(status$stdout)),
        status_sha256 = digest::digest(
            status$stdout,
            algo = "sha256",
            serialize = FALSE
        )
    )
}

publicationValidateComparatorWorktree <- function(record) {
    publicationRequireNames(record, c(
        "comparator_id", "path", "revision", "source_identity",
        "detached", "clean", "status_sha256"
    ), "Comparator worktree")
    expected_revision <- publicationComparatorRevisions()[[record$comparator_id]]
    observed <- publicationComparatorWorktreeState(record$path)
    valid <- !is.null(expected_revision) &&
        identical(record$revision, expected_revision) &&
        identical(observed$revision, expected_revision) &&
        isTRUE(record$detached) && isTRUE(observed$detached) &&
        isTRUE(record$clean) && isTRUE(observed$clean) &&
        identical(record$status_sha256, observed$status_sha256) &&
        identical(
            publicationObjectDigest(record$source_identity),
            publicationObjectDigest(
                publicationComparatorSourceIdentity(expected_revision)
            )
        )
    if (!valid) {
        publicationComparatorAbort("Comparator worktree binding differs")
    }
    invisible(record)
}

publicationComparatorAuthorityPaths <- function() {
    c(
        lockfile = "renv.lock",
        lock_scope =
            "tools/profiling/omics-publication-lock-scope-v1.json",
        repositories =
            "tools/profiling/omics-publication-repositories-v1.json",
        system_dependencies =
            "tools/profiling/omics-publication-system-dependencies-v1.json"
    )
}

publicationExternalBuildInputAuthority <- function() {
    list(
        list(
            input_id = "libv8-static-12.4.254.21-rocky-8-amd64",
            consumer_package = "V8",
            consumer_version = "8.2.0",
            role = "omics_report_and_lock_schema_verification",
            url = paste0(
                "https://github.com/jeroen/build-v8-static/releases/download/",
                "12.4.254.21/v8-rocky-8-amd64.tar.gz"
            ),
            size_bytes = 19512263L,
            sha256 = paste0(
                "85de4fda7c5dde682792a54572fa0310",
                "8479fe0b84762d2b1d5beef675577493"
            ),
            published_digest_source = "github_release_asset_metadata",
            published_at = "2025-08-23T23:10:27Z",
            verification = "sha256_before_extraction",
            measured_worker_load_allowed = FALSE
        ),
        list(
            input_id = "renv-bootstrap-1.2.3",
            consumer_package = "renv",
            consumer_version = "1.2.3",
            role = "common_environment_bootstrap",
            url = paste0(
                "https://packagemanager.posit.co/cran/2026-08-26/",
                "src/contrib/Archive/renv/renv_1.2.3.tar.gz"
            ),
            size_bytes = 1375285L,
            sha256 = paste0(
                "670719f57c532d292c1dd793f8fa7d06",
                "e2d6f91205dc90dd947e7e8c05a5ae06"
            ),
            published_digest_source = "frozen_repository_archive_capture",
            captured_at = "2026-08-27T16:17:00.000Z",
            verification = "sha256_before_install",
            measured_worker_load_allowed = FALSE
        ),
        list(
            input_id = "libarrow-linux-x86_64-24.0.0",
            consumer_package = "arrow",
            consumer_version = "24.0.0",
            role = "arrow_cpp_runtime",
            url = paste0(
                "https://github.com/apache/arrow/releases/download/",
                "apache-arrow-24.0.0/r-libarrow-linux-x86_64-24.0.0.zip"
            ),
            size_bytes = 43868290L,
            sha256 = paste0(
                "88d227deb516ef6177cf03b6fedfb333",
                "7a9fcea7a7336801d1a184879ad582d1"
            ),
            sha512 = paste0(
                "298dbbfcb34c291ec7b44a207660e647fbe180d8643c0f59f01a48e2c1561ec",
                "999a39ac43ab3e7968885400ec2cdbdd09cdb2a3fcf765454b5cdf61c38bf408e"
            ),
            published_digest_source =
                "github_release_asset_metadata_and_arrow_source_checksum",
            published_at = "2026-04-21T07:39:37Z",
            verification = "sha256_and_sha512_before_install",
            measured_worker_load_allowed = TRUE
        )
    )
}

publicationValidateRepositoryAuthority <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "repository_authority_id",
        "owner_ticket_id", "captured_at", "status", "r_version",
        "bioconductor_version", "cran_snapshot_date", "repositories",
        "remotes", "external_build_inputs", "policy"
    ), "Repository authority")
    repository_ids <- vapply(
        record$repositories,
        `[[`,
        character(1),
        "repository_id"
    )
    expected_ids <- c(
        "CRAN", "BioCsoft", "BioCann", "BioCexp", "BioCworkflows",
        "BioCbooks"
    )
    repositories_valid <- identical(repository_ids, expected_ids) &&
        all(vapply(record$repositories, \(repository) {
            publicationScalarString(repository$url) &&
                !grepl("/latest/?$", repository$url) &&
                identical(repository$index_path, "src/contrib/PACKAGES.gz") &&
                publicationComparatorShaValid(repository$index_sha256)
        }, logical(1)))
    remote_ids <- vapply(record$remotes, `[[`, character(1), "package")
    remotes_valid <- identical(remote_ids, c("RUVIIIC", "GlimmaV2")) &&
        all(vapply(record$remotes, \(remote) {
            publicationComparatorShaValid(remote$commit, 40L) &&
                publicationComparatorShaValid(remote$archive_sha256) &&
                grepl(remote$commit, remote$archive_url, fixed = TRUE)
        }, logical(1)))
    external_valid <- identical(
        publicationObjectDigest(record$external_build_inputs),
        publicationObjectDigest(publicationExternalBuildInputAuthority())
    )
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_repositories"
    ) && identical(record$schema_version, "1.0.0") &&
        identical(record$owner_ticket_id, "OMICS-ART-064") &&
        identical(record$status, "frozen") &&
        identical(record$r_version, "4.6.1") &&
        identical(record$bioconductor_version, "3.23") &&
        identical(record$cran_snapshot_date, "2026-08-26") &&
        repositories_valid && remotes_valid && external_valid &&
        identical(record$policy$latest_urls_allowed, FALSE) &&
        identical(record$policy$branch_remotes_allowed, FALSE) &&
        isTRUE(record$policy$source_packages_required) &&
        isTRUE(record$policy$index_digest_required) &&
        isTRUE(record$policy$external_build_input_digest_required)
    if (!valid) publicationComparatorAbort("Repository authority differs")
    invisible(record)
}

publicationValidateLockAuthority <- function(path, repositories, scope) {
    lock <- publicationReadJson(path)
    publicationRequireNames(
        lock,
        c("R", "Bioconductor", "Packages"),
        "renv lockfile"
    )
    repository_urls <- stats::setNames(
        vapply(
            repositories$repositories,
            `[[`,
            character(1),
            "url"
        ),
        vapply(
            repositories$repositories,
            `[[`,
            character(1),
            "repository_id"
        )
    )
    locked_urls <- stats::setNames(
        vapply(lock$R$Repositories, `[[`, character(1), "URL"),
        vapply(lock$R$Repositories, `[[`, character(1), "Name")
    )
    publicationValidateLockScopeAuthority(scope)
    required <- publicationLockRoots(scope$root_groups)
    sources <- vapply(lock$Packages, `[[`, character(1), "Source")
    package_names_match <- all(vapply(names(lock$Packages), \(package) {
        identical(lock$Packages[[package]]$Package, package)
    }, logical(1)))
    closure <- publicationLockDependencyClosure(lock, required)
    encoded <- jsonlite::toJSON(lock, auto_unbox = TRUE, null = "null")
    remote_valid <- all(vapply(c("RUVIIIC", "Glimma"), \(package) {
        entry <- lock$Packages[[package]]
        identical(entry$Source, "GitHub") &&
            publicationComparatorShaValid(entry$RemoteSha, 40L) &&
            identical(entry$RemoteRef, entry$RemoteSha)
    }, logical(1)))
    valid <- identical(lock$R$Version, "4.6.1") &&
        identical(lock$Bioconductor$Version, "3.23") &&
        identical(locked_urls, repository_urls) &&
        identical(names(lock$Packages), closure$packages) &&
        !length(closure$unresolved) && package_names_match &&
        identical(publicationFileDigest(path), scope$lockfile$sha256) &&
        all(required %in% names(lock$Packages)) &&
        all(sources %in% unlist(scope$policy$package_sources, use.names = FALSE)) &&
        remote_valid &&
        !grepl("/home/|file://|RemoteRef\\\":\\\"HEAD", encoded)
    publicationLockRequireDirectDependency(lock, "gt", "juicyjuice")
    publicationLockRequireDirectDependency(lock, "juicyjuice", "V8")
    publicationLockRequireDirectDependency(lock, "jsonvalidate", "V8")
    if (!valid) publicationComparatorAbort("renv lockfile authority differs")
    invisible(lock)
}

publicationValidateComparatorCommonEnvironment <- function(environment) {
    publicationRequireNames(environment, c(
        "environment_id", "r_version", "platform", "lockfile", "lock_scope",
        "repositories", "system_dependencies", "restore_evidence",
        "restore_policy"
    ), "Comparator common environment")
    paths <- publicationComparatorAuthorityPaths()
    bindings <- c(
        list(environment$lockfile),
        list(environment$lock_scope),
        list(environment$repositories),
        list(environment$system_dependencies)
    )
    bound_paths <- vapply(bindings, `[[`, character(1), "path")
    valid <- publicationScalarString(environment$environment_id) &&
        publicationScalarString(environment$r_version) &&
        publicationScalarString(environment$platform) &&
        identical(unname(bound_paths), unname(paths)) &&
        all(vapply(bindings, \(binding) {
            publicationComparatorShaValid(binding$sha256) &&
                file.exists(publicationPath(binding$path)) &&
                identical(
                    publicationFileDigest(binding$path),
                    binding$sha256
                )
        }, logical(1))) &&
        identical(environment$restore_policy$ambient_user_library, FALSE) &&
        identical(environment$restore_policy$common_lock_required, TRUE) &&
        identical(environment$restore_policy$restore_repetitions, 2L) &&
        identical(environment$restore_evidence$status, "verified") &&
        identical(environment$restore_evidence$restore_repetitions, 2L) &&
        length(environment$restore_evidence$receipt_sha256) == 2L &&
        all(vapply(
            environment$restore_evidence$receipt_sha256,
            publicationComparatorShaValid,
            logical(1)
        )) && publicationComparatorShaValid(
            environment$restore_evidence$reproducibility_sha256
        )
    if (!valid) {
        publicationComparatorAbort("Comparator common environment differs")
    }
    repositories <- publicationReadJson(environment$repositories$path)
    lock_scope <- publicationReadJson(environment$lock_scope$path)
    publicationValidateRepositoryAuthority(repositories)
    publicationValidateLockAuthority(
        environment$lockfile$path,
        repositories,
        lock_scope
    )
    publicationValidateSystemDependencyAuthority(publicationReadJson(
        environment$system_dependencies$path
    ))
    invisible(environment)
}

publicationComparatorRoutes <- function() {
    list(
        list(
            comparison_id = "release_effect.historical_memory_vs_candidate_auto",
            claim_class = "release_effect",
            numerator = "candidate::proposed_auto",
            denominator = "historical_janitor::memory"
        ),
        list(
            comparison_id = "backend_effect.candidate_artifact_vs_memory",
            claim_class = "backend_effect",
            numerator = "candidate::artifact",
            denominator = "candidate::memory"
        ),
        list(
            comparison_id = "repair_effect.candidate_vs_pre_repair_artifact",
            claim_class = "repair_effect",
            numerator = "candidate::artifact",
            denominator = "pre_repair_performance::artifact"
        )
    )
}

publicationValidateComparatorRoutes <- function(routes) {
    expected <- publicationComparatorRoutes()
    ids <- vapply(routes, `[[`, character(1), "comparison_id")
    valid <- !anyDuplicated(ids) && identical(
        publicationObjectDigest(routes),
        publicationObjectDigest(expected)
    ) && length(unique(vapply(
        routes,
        `[[`,
        character(1),
        "claim_class"
    ))) == 3L
    if (!valid) {
        publicationComparatorAbort("Comparator claim routes differ")
    }
    invisible(routes)
}

publicationComparatorEntryContracts <- function() {
    list(
        historical_janitor = list(
            role = "historical_release",
            backend_envelopes = list("memory"),
            claim_classes = list("release_effect")
        ),
        pre_repair_performance = list(
            role = "pre_repair_canary",
            backend_envelopes = list("memory", "artifact"),
            claim_classes = list("repair_effect", "release_effect")
        ),
        candidate = list(
            role = "future_candidate",
            backend_envelopes = list("memory", "artifact", "proposed_auto"),
            claim_classes = list(
                "release_effect", "backend_effect", "repair_effect"
            )
        )
    )
}

publicationValidateComparatorEntry <- function(entry) {
    publicationRequireNames(entry, c(
        "comparator_id", "role", "revision", "status", "source_identity",
        "backend_envelopes", "claim_classes", "build_receipts",
        "execution_authority"
    ), "Comparator entry")
    fixed <- publicationComparatorRevisions()
    contract <- publicationComparatorEntryContracts()[[entry$comparator_id]]
    contract_valid <- !is.null(contract) &&
        identical(entry$role, contract$role) &&
        identical(entry$backend_envelopes, contract$backend_envelopes) &&
        identical(entry$claim_classes, contract$claim_classes)
    if (entry$comparator_id %in% names(fixed)) {
        evidence_valid <- if (identical(entry$status, "verified")) {
            receipt_ids <- vapply(
                entry$build_receipts,
                `[[`,
                character(1),
                "evidence_id"
            )
            evidence_types <- vapply(
                entry$build_receipts,
                `[[`,
                character(1),
                "evidence_type"
            )
            length(entry$build_receipts) == 3L && !anyDuplicated(receipt_ids) &&
                identical(evidence_types, c(
                    "clean_build_receipt", "clean_build_receipt",
                    "build_pair_reproducibility"
                )) && all(vapply(
                entry$build_receipts,
                \(evidence) publicationScalarString(evidence$evidence_id) &&
                    evidence$evidence_type %in% c(
                        "clean_build_receipt",
                        "build_pair_reproducibility"
                    ) && publicationComparatorShaValid(evidence$sha256),
                logical(1)
            ))
        } else {
            !length(entry$build_receipts)
        }
        valid <- contract_valid &&
            identical(entry$revision, fixed[[entry$comparator_id]]) &&
            entry$status %in% c("pending_clean_builds", "verified") &&
            identical(entry$source_identity$revision, entry$revision) &&
            identical(
                publicationObjectDigest(entry$source_identity),
                publicationObjectDigest(publicationComparatorSourceIdentity(
                    entry$revision
                ))
            ) &&
            !isTRUE(entry$execution_authority) &&
            is.list(entry$build_receipts) && evidence_valid
    } else if (identical(entry$comparator_id, "candidate")) {
        valid <- contract_valid && is.null(entry$revision) &&
            identical(entry$status, "pending_omics_art_077") &&
            is.null(entry$source_identity) &&
            !length(entry$build_receipts) &&
            !isTRUE(entry$execution_authority)
    } else {
        valid <- FALSE
    }
    if (!valid) publicationComparatorAbort("Comparator entry differs")
    invisible(entry)
}

publicationValidateComparatorAuthority <- function(authority) {
    publicationRequireNames(authority, c(
        "schema", "schema_version", "authority_id", "owner_ticket_id",
        "status", "common_environment", "comparators", "comparison_routes",
        "execution_authority", "promotion_authority", "immutability"
    ), "Comparator authority")
    ids <- vapply(authority$comparators, `[[`, character(1), "comparator_id")
    valid <- identical(
        authority$schema,
        "multischolar.omics_publication_comparators"
    ) && identical(authority$schema_version, "1.0.0") &&
        identical(authority$owner_ticket_id, "OMICS-ART-064") &&
        identical(
            authority$status,
            "fixed_revisions_verified_candidate_pending"
        ) &&
        identical(ids, c(
            "historical_janitor", "pre_repair_performance", "candidate"
        )) && !isTRUE(authority$execution_authority) &&
        !isTRUE(authority$promotion_authority) &&
        isTRUE(authority$immutability$fixed_worktrees_read_only) &&
        isTRUE(authority$immutability$candidate_successor_required)
    if (!valid) publicationComparatorAbort("Comparator authority differs")
    publicationValidateComparatorCommonEnvironment(authority$common_environment)
    lapply(authority$comparators, publicationValidateComparatorEntry)
    publicationValidateComparatorRoutes(authority$comparison_routes)
    invisible(authority)
}

publicationComparatorFreezeCandidate <- function(
    authority,
    revision,
    source_identity,
    package_sha256,
    bindings,
    readiness
) {
    if (identical(
        authority$schema,
        "multischolar.omics_publication_candidate_freeze_successor"
    )) {
        publicationComparatorAbort("Candidate is already frozen")
    }
    publicationValidateComparatorAuthority(authority)
    candidate <- authority$comparators[[3L]]
    if (!identical(candidate$status, "pending_omics_art_077")) {
        publicationComparatorAbort("Candidate is already frozen")
    }
    publicationComparatorResolveRevision(revision)
    if (!publicationComparatorShaValid(package_sha256)) {
        publicationComparatorAbort("Candidate package digest is invalid")
    }
    publicationValidateCandidateFreezeBindings(bindings)
    publicationValidateCandidateFreezeReadiness(readiness)
    if (!isTRUE(readiness$candidate_freeze_allowed)) {
        publicationComparatorAbort("Candidate freeze readiness is blocked")
    }
    list(
        schema = "multischolar.omics_publication_candidate_freeze_successor",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-077",
        predecessor_authority_id = authority$authority_id,
        predecessor_sha256 = publicationObjectDigest(authority),
        candidate_revision = revision,
        source_identity = source_identity,
        package_sha256 = package_sha256,
        bindings = bindings,
        binding_sha256 = publicationObjectDigest(bindings),
        readiness = readiness,
        frozen_once = TRUE,
        promotion_authority = FALSE
    )
}

publicationCandidateFreezeBindingNames <- function() {
    c(
        "build", "library", "native", "dependency_environment",
        "protocol", "comparators", "workloads", "analysis",
        "scientific_differences", "payload_access", "policy_receipts",
        "handoff", "roles", "projects", "estimands", "splits",
        "threshold_grid", "campaign_budget", "retry_policy", "blind_labels"
    )
}

publicationCandidateFreezeReadinessPaths <- function() {
    c(
        roles = "tests/testdata/omics-performance/roles-v1.json",
        projects = "tests/testdata/omics-performance/projects-v1.json",
        splits = "tests/testdata/omics-performance/splits-v1.json",
        campaign_budget = paste0(
            "tests/testdata/omics-performance/campaign-budget-v1.json"
        ),
        payload_access = "tools/refactor/omics-payload-access-owners-v1.json",
        auxiliary = paste0(
            "tests/testdata/omics-performance/auxiliary/manifest-v1.json"
        ),
        proteomics = paste0(
            "tests/testdata/omics-performance/proteomics/",
            "governance-successor-v1.json"
        ),
        metabolomics = paste0(
            "tests/testdata/omics-performance/metabolomics/",
            "governance-successor-v1.json"
        ),
        lipidomics = paste0(
            "tests/testdata/omics-performance/lipidomics/",
            "governance-successor-v1.json"
        ),
        policy_receipts = paste0(
            "tests/testdata/omics-performance/proposed-policy-receipts-v1.json"
        ),
        handoff = "tests/testdata/omics-performance/handoff-v1.json",
        blind_labels = "tests/testdata/omics-performance/blind-labels-v1.json"
    )
}

publicationCandidateReadinessBinding <- function(path) {
    available <- file.exists(publicationPath(path))
    list(
        path = path,
        sha256 = if (available) publicationFileDigest(path) else NULL,
        available = available
    )
}

publicationCandidateFreezeBlockers <- function(records) {
    blockers <- character()
    roles <- records$roles
    principal_complete <- all(vapply(roles$roles, \(role) {
        publicationScalarString(role$principal_id)
    }, logical(1)))
    if (!principal_complete) blockers <- c(blockers, "role_principals_incomplete")
    if (!isTRUE(roles$readiness$runtime_implementation_authorized)) {
        blockers <- c(blockers, "runtime_implementation_review_missing")
    }
    if (!isTRUE(roles$readiness$campaign_execution_authorized)) {
        blockers <- c(blockers, "campaign_role_handoff_missing")
    }
    projects <- records$projects
    project_counts <- c(
        vapply(projects$full_workflow_claims, \(claim) {
            claim$verified_real_project_count >= claim$required_real_project_count
        }, logical(1)),
        vapply(projects$auxiliary_claims, \(claim) {
            claim$verified_real_project_count >=
                projects$minimum_real_project_authorities_per_cross_project_claim
        }, logical(1))
    )
    if (!all(project_counts)) {
        blockers <- c(blockers, "project_source_authority_incomplete")
    }
    if (!isTRUE(records$splits$readiness$candidate_access_allowed)) {
        blockers <- c(blockers, "split_successor_candidate_access_denied")
    }
    if (!isTRUE(records$campaign_budget$current_status$execution_authorized)) {
        blockers <- c(blockers, "campaign_budget_execution_unauthorized")
    }
    if (!isTRUE(records$payload_access$candidate_freeze_authority)) {
        blockers <- c(blockers, "payload_access_freeze_authority_pending")
    }
    for (name in c("auxiliary", "proteomics", "metabolomics", "lipidomics")) {
        if (!isTRUE(records[[name]]$candidate_access_allowed)) {
            blockers <- c(blockers, paste0(name, "_candidate_access_denied"))
        }
    }
    unavailable <- names(records)[!vapply(records, is.list, logical(1))]
    c(blockers, paste0(unavailable, "_authority_missing"))
}

publicationBuildCandidateFreezeReadiness <- function() {
    paths <- publicationCandidateFreezeReadinessPaths()
    bindings <- lapply(paths, publicationCandidateReadinessBinding)
    records <- lapply(bindings, \(binding) {
        if (!isTRUE(binding$available)) return(FALSE)
        publicationReadJson(binding$path)
    })
    blockers <- unique(publicationCandidateFreezeBlockers(records))
    record <- list(
        schema = "multischolar.omics_publication_candidate_freeze_readiness",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-077",
        status = if (length(blockers)) "blocked" else "ready",
        authority_bindings = bindings,
        blockers = as.list(blockers),
        candidate_freeze_allowed = !length(blockers),
        promotion_authority = FALSE
    )
    record$readiness_digest <- publicationObjectDigest(record)
    record
}

publicationValidateCandidateFreezeReadiness <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "owner_ticket_id", "status",
        "authority_bindings", "blockers", "candidate_freeze_allowed",
        "promotion_authority", "readiness_digest"
    ), "Candidate freeze readiness")
    basis <- record
    basis$readiness_digest <- NULL
    bindings_valid <- length(record$authority_bindings) > 0L && all(vapply(
        record$authority_bindings,
        \(binding) {
            publicationRequireNames(
                binding,
                c("path", "sha256", "available"),
                "Candidate readiness binding"
            )
            publicationScalarString(binding$path) &&
                is.logical(binding$available) && length(binding$available) == 1L &&
                identical(
                    isTRUE(binding$available),
                    publicationComparatorShaValid(binding$sha256)
                ) && if (isTRUE(binding$available)) {
                    file.exists(publicationPath(binding$path)) && identical(
                        binding$sha256,
                        publicationFileDigest(binding$path)
                    )
                } else {
                    !file.exists(publicationPath(binding$path))
                }
        },
        logical(1)
    ))
    allowed <- !length(record$blockers)
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_candidate_freeze_readiness"
    ) && identical(record$schema_version, "1.0.0") &&
        identical(record$owner_ticket_id, "OMICS-ART-077") &&
        identical(record$status, if (allowed) "ready" else "blocked") &&
        identical(isTRUE(record$candidate_freeze_allowed), allowed) &&
        bindings_valid && !isTRUE(record$promotion_authority) &&
        identical(record$readiness_digest, publicationObjectDigest(basis))
    if (!valid) publicationComparatorAbort("Candidate freeze readiness differs")
    invisible(record)
}

publicationValidateCandidateFreezeBindings <- function(bindings) {
    names_valid <- identical(
        names(bindings),
        publicationCandidateFreezeBindingNames()
    )
    digests_valid <- names_valid && all(vapply(
        bindings,
        publicationComparatorShaValid,
        logical(1)
    ))
    if (!digests_valid) {
        publicationComparatorAbort("Candidate freeze bindings differ")
    }
    invisible(bindings)
}

publicationValidateCandidateFreezeSuccessor <- function(
    record,
    authority,
    current_bindings = record$bindings
) {
    publicationValidateComparatorAuthority(authority)
    publicationRequireNames(record, c(
        "schema", "schema_version", "owner_ticket_id",
        "predecessor_authority_id", "predecessor_sha256",
        "candidate_revision", "source_identity", "package_sha256",
        "bindings", "binding_sha256", "readiness", "frozen_once",
        "promotion_authority"
    ), "Candidate freeze successor")
    publicationValidateCandidateFreezeBindings(record$bindings)
    publicationValidateCandidateFreezeBindings(current_bindings)
    publicationValidateCandidateFreezeReadiness(record$readiness)
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_candidate_freeze_successor"
    ) && identical(record$schema_version, "1.0.0") &&
        identical(record$owner_ticket_id, "OMICS-ART-077") &&
        identical(record$predecessor_authority_id, authority$authority_id) &&
        identical(record$predecessor_sha256, publicationObjectDigest(authority)) &&
        publicationComparatorShaValid(record$candidate_revision, 40L) &&
        identical(record$source_identity$revision, record$candidate_revision) &&
        identical(
            publicationObjectDigest(record$source_identity),
            publicationObjectDigest(publicationComparatorSourceIdentity(
                record$candidate_revision
            ))
        ) &&
        publicationComparatorShaValid(record$package_sha256) &&
        identical(record$bindings$build, record$package_sha256) &&
        identical(record$binding_sha256, publicationObjectDigest(record$bindings)) &&
        identical(
            publicationObjectDigest(current_bindings),
            record$binding_sha256
        ) &&
        isTRUE(record$readiness$candidate_freeze_allowed) &&
        isTRUE(record$frozen_once) && !isTRUE(record$promotion_authority)
    if (!valid) publicationComparatorAbort("Candidate freeze successor differs")
    invisible(record)
}

publicationComparatorIncomparableReceipt <- function(reason, evidence) {
    if (!publicationScalarString(reason) || !length(evidence)) {
        publicationComparatorAbort("Incomparable comparator evidence is empty")
    }
    list(
        schema = "multischolar.omics_publication_comparator_outcome",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-064",
        outcome = "incomparable_dependency_environment",
        reason = reason,
        evidence = evidence,
        historical_comparison_authority = FALSE,
        sensitivity_substitution_allowed = FALSE,
        promotion_authority = FALSE
    )
}

publicationValidateComparatorIncomparable <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "owner_ticket_id", "outcome", "reason",
        "evidence", "historical_comparison_authority",
        "sensitivity_substitution_allowed", "promotion_authority"
    ), "Comparator incomparable outcome")
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_comparator_outcome"
    ) && identical(record$schema_version, "1.0.0") &&
        identical(record$owner_ticket_id, "OMICS-ART-064") &&
        identical(record$outcome, "incomparable_dependency_environment") &&
        publicationScalarString(record$reason) && length(record$evidence) > 0L &&
        !isTRUE(record$historical_comparison_authority) &&
        !isTRUE(record$sensitivity_substitution_allowed) &&
        !isTRUE(record$promotion_authority)
    if (!valid) publicationComparatorAbort("Comparator outcome differs")
    invisible(record)
}

publicationScientificDifferenceInventory <- function(base, head) {
    name_status <- publicationComparatorGit(c("diff", "--name-status", base, head))
    commits <- publicationComparatorGit(c(
        "rev-list", "--reverse", paste0(base, "..", head)
    ))
    if (name_status$status != 0L || commits$status != 0L) {
        publicationComparatorAbort("Scientific difference inventory failed")
    }
    paths <- strsplit(name_status$stdout, "\n", fixed = TRUE)[[1L]]
    paths <- paths[nzchar(paths)]
    paths <- sub("^[^\\t]+\\t", "", paths)
    paths <- sub("^.*\\t", "", paths)
    commit_ids <- strsplit(commits$stdout, "\n", fixed = TRUE)[[1L]]
    commit_ids <- commit_ids[nzchar(commit_ids)]
    list(
        range_id = "historical_janitor_to_pre_repair_performance",
        base_revision = base,
        head_revision = head,
        commit_order = "git_rev_list_reverse",
        commit_count = as.integer(length(commit_ids)),
        commit_set_sha256 = digest::digest(
            commits$stdout,
            algo = "sha256",
            serialize = FALSE
        ),
        changed_path_count = as.integer(length(paths)),
        name_status_sha256 = digest::digest(
            name_status$stdout,
            algo = "sha256",
            serialize = FALSE
        ),
        changed_paths = paths,
        commits = commit_ids
    )
}

publicationScientificCommitInRange <- function(commit, inventory) {
    if (!publicationComparatorShaValid(commit, 40L)) return(FALSE)
    commit %in% inventory$commits
}

publicationValidateScientificEvidence <- function(evidence, inventory) {
    valid <- vapply(evidence, \(item) {
        if (startsWith(item, "commit:")) {
            commit <- sub("^commit:", "", item)
            return(publicationScientificCommitInRange(commit, inventory))
        }
        if (startsWith(item, "commit-range:")) {
            range <- strsplit(
                sub("^commit-range:", "", item),
                "..",
                fixed = TRUE
            )[[1L]]
            return(length(range) == 2L && all(vapply(
                range,
                publicationScientificCommitInRange,
                logical(1),
                inventory = inventory
            )))
        }
        if (startsWith(item, "tests:")) {
            test <- sub("^tests:", "", item)
            return(!grepl("/", test, fixed = TRUE) && file.exists(
                publicationPath("tests", "testthat", test)
            ))
        }
        identical(item, paste0("inventory:", inventory$range_id))
    }, logical(1))
    if (!length(evidence) || !all(valid)) {
        publicationComparatorAbort("Scientific difference evidence differs")
    }
    invisible(evidence)
}

publicationValidateScientificFiles <- function(files, inventory) {
    valid <- vapply(files, \(path) {
        if (endsWith(path, "/")) {
            any(startsWith(inventory$changed_paths, path))
        } else {
            path %in% inventory$changed_paths
        }
    }, logical(1))
    if (!length(files) || !all(valid)) {
        publicationComparatorAbort("Scientific difference files differ")
    }
    invisible(files)
}

publicationValidateScientificDifferences <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "difference_authority_id",
        "owner_ticket_id", "status", "revisions", "inventory", "records",
        "candidate_amendment",
        "matched_candidate_rule", "promotion_authority"
    ), "Scientific differences")
    revisions <- publicationComparatorRevisions()
    categories <- c(
        "scientific", "correctness", "storage_runtime", "ui",
        "state_runtime", "documentation", "test_only"
    )
    inventory <- publicationScientificDifferenceInventory(
        revisions[["historical_janitor"]],
        revisions[["pre_repair_performance"]]
    )
    inventory$changed_paths <- NULL
    inventory$commits <- NULL
    current <- publicationScientificDifferenceInventory(
        revisions[["historical_janitor"]],
        revisions[["pre_repair_performance"]]
    )
    publicationRequireNames(record$candidate_amendment, c(
        "owner_ticket_id", "status", "successor_required",
        "predecessor_difference_authority_id", "required_revision_span",
        "freeze_binding", "exact_candidate_memory_artifact_parity_required"
    ), "Scientific difference candidate amendment")
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_scientific_differences"
    ) && identical(record$schema_version, "1.0.0") &&
        identical(record$owner_ticket_id, "OMICS-ART-064") &&
        identical(
            record$status,
            "fixed_revision_inventory_complete_candidate_pending"
        ) &&
        identical(unlist(record$revisions, use.names = TRUE), revisions) &&
        identical(record$inventory, inventory) &&
        identical(
            record$candidate_amendment$owner_ticket_id,
            "OMICS-ART-077"
        ) && identical(record$candidate_amendment$status, "required") &&
        identical(
            record$candidate_amendment$predecessor_difference_authority_id,
            record$difference_authority_id
        ) && identical(
            record$candidate_amendment$required_revision_span,
            "pre_repair_performance_to_candidate"
        ) &&
        identical(
            record$candidate_amendment$freeze_binding,
            "scientific_differences"
        ) && isTRUE(record$candidate_amendment$successor_required) &&
        isTRUE(
            record$candidate_amendment$exact_candidate_memory_artifact_parity_required
        ) &&
        identical(record$matched_candidate_rule, "exact_science_required") &&
        !isTRUE(record$promotion_authority)
    ids <- character()
    for (difference in record$records) {
        publicationRequireNames(difference, c(
            "difference_id", "category", "files", "review_status",
            "file_scope", "expected_scientific_impact", "claim_scope",
            "evidence"
        ), "Scientific difference")
        ids <- c(ids, difference$difference_id)
        publicationValidateScientificFiles(difference$files, current)
        publicationValidateScientificEvidence(difference$evidence, current)
        valid <- valid && difference$category %in% categories &&
            identical(difference$review_status, "reviewed") &&
            identical(difference$file_scope, "focal_owners_not_exhaustive") &&
            publicationScalarString(difference$expected_scientific_impact) &&
            publicationScalarString(difference$claim_scope) &&
            length(difference$evidence) > 0L
    }
    if (!valid || anyDuplicated(ids)) {
        publicationComparatorAbort("Scientific difference authority differs")
    }
    invisible(record)
}
