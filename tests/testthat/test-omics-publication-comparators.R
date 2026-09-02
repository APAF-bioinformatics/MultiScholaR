publicationComparatorAuthorityFixture <- function() {
    publicationReadJson(
        "tests/testdata/omics-performance/comparators-v1.json"
    )
}

publicationComparatorReceiptFixture <- function(
    comparator_id = "historical_janitor",
    receipt_id = "receipt-1"
) {
    authority <- publicationComparatorAuthorityFixture()
    comparator <- authority$comparators[[which(vapply(
        authority$comparators,
        `[[`,
        character(1),
        "comparator_id"
    ) == comparator_id)]]
    digest <- paste(rep("a", 64L), collapse = "")
    command <- list(
        command = "R",
        arguments = list("CMD", "build"),
        workdir = "/tmp/build",
        environment = list(TZ = "UTC"),
        exit_status = 0L,
        timed_out = FALSE,
        elapsed_seconds = 1,
        stdout_sha256 = digest,
        stderr_sha256 = digest,
        stdout_path = "/tmp/stdout",
        stderr_path = "/tmp/stderr"
    )
    list(
        schema = "multischolar.omics_publication_comparator_build",
        schema_version = "1.0.0",
        receipt_id = receipt_id,
        comparator_id = comparator_id,
        revision = comparator$revision,
        source = list(source_identity = comparator$source_identity),
        environment = list(
            lockfile_sha256 = authority$common_environment$lockfile$sha256
        ),
        restore = list(status = "passed"),
        build = list(
            build_command = command,
            install_command = command,
            archive = list(sha256 = digest),
            installed_inventory = list(inventory_sha256 = digest),
            native_libraries = list(list(path = "MultiScholaR.so", sha256 = digest))
        ),
        session = list(r_version = "4.6.1"),
        smoke = list(status = "passed", output_sha256 = digest),
        cleanup = list(
            worktree_removed = TRUE,
            dependency_library_retained = TRUE,
            package_library_retained = TRUE
        ),
        status = "verified",
        publication_authority = FALSE
    )
}

test_that("comparator repository system lock and source authorities validate", {
    authority <- publicationComparatorAuthorityFixture()
    differences <- publicationReadJson(
        "tests/testdata/omics-performance/scientific-differences-v1.json"
    )

    expect_silent(publicationValidateComparatorAuthority(authority))
    expect_silent(publicationValidateScientificDifferences(differences))
    expect_identical(
        vapply(authority$comparators[1:2], `[[`, character(1), "revision"),
        unname(publicationComparatorRevisions())
    )
    expect_null(authority$comparators[[3L]]$revision)
    expect_false(authority$execution_authority)
    expect_false(authority$promotion_authority)

    for (entry in authority$comparators[1:2]) {
        expect_identical(
            publicationObjectDigest(entry$source_identity),
            publicationObjectDigest(
                publicationComparatorSourceIdentity(entry$revision)
            )
        )
    }
})

test_that("comparator authorities reject lock repository system and route drift", {
    authority <- publicationComparatorAuthorityFixture()

    lock_drift <- publicationGovernanceCopy(authority)
    lock_drift$common_environment$lockfile$sha256 <- paste(
        rep("0", 64L),
        collapse = ""
    )
    expect_error(
        publicationValidateComparatorAuthority(lock_drift),
        class = "multischolar_publication_comparator_error"
    )

    source_drift <- publicationGovernanceCopy(authority)
    source_drift$comparators[[1L]]$source_identity$tree_oid <- paste(
        rep("0", 40L),
        collapse = ""
    )
    expect_error(
        publicationValidateComparatorAuthority(source_drift),
        class = "multischolar_publication_comparator_error"
    )

    route_drift <- publicationGovernanceCopy(authority)
    route_drift$comparison_routes[[1L]]$claim_class <- "backend_effect"
    expect_error(
        publicationValidateComparatorAuthority(route_drift),
        class = "multischolar_publication_comparator_error"
    )

    missing_capability <- publicationGovernanceCopy(authority)
    missing_capability$comparators[[2L]]$backend_envelopes <- list("memory")
    expect_error(
        publicationValidateComparatorAuthority(missing_capability),
        class = "multischolar_publication_comparator_error"
    )

    repositories <- publicationReadJson(
        "tools/profiling/omics-publication-repositories-v1.json"
    )
    repositories$repositories[[1L]]$url <-
        "https://packagemanager.posit.co/cran/latest"
    expect_error(
        publicationValidateRepositoryAuthority(repositories),
        class = "multischolar_publication_comparator_error"
    )

    repositories <- publicationReadJson(
        "tools/profiling/omics-publication-repositories-v1.json"
    )
    repositories$external_build_inputs[[1L]]$sha256 <- paste(
        rep("0", 64L),
        collapse = ""
    )
    expect_error(
        publicationValidateRepositoryAuthority(repositories),
        class = "multischolar_publication_comparator_error"
    )

    system <- publicationReadJson(
        "tools/profiling/omics-publication-system-dependencies-v1.json"
    )
    system$commands[[1L]]$sha256 <- paste(rep("0", 64L), collapse = "")
    expect_error(
        publicationValidateSystemDependencyAuthority(system),
        class = "multischolar_publication_comparator_build_error"
    )
})

test_that("lock scope is exact complete and measurement-isolated", {
    scope <- publicationReadJson(
        "tools/profiling/omics-publication-lock-scope-v1.json"
    )
    lock <- publicationReadJson("renv.lock")

    expect_silent(publicationValidateLockScopeAuthority(scope))
    expect_identical(scope$resolved$package_count, 384L)
    expect_true(all(c(
        "dynamicTreeCut", "GSEABase", "MOFA2", "V8", "jsonvalidate"
    ) %in% names(lock$Packages)))
    expect_true(publicationLockMetadataValid(lock))
    expect_silent(publicationLockRequireDirectDependency(
        lock,
        "gt",
        "juicyjuice"
    ))
    expect_silent(publicationLockRequireDirectDependency(
        lock,
        "juicyjuice",
        "V8"
    ))
    expect_silent(publicationLockRequireDirectDependency(
        lock,
        "jsonvalidate",
        "V8"
    ))

    missing_science <- publicationGovernanceCopy(lock)
    missing_science$Packages$MOFA2 <- NULL
    closure <- publicationLockDependencyClosure(
        missing_science,
        publicationLockRoots(scope$root_groups)
    )
    expect_true("MOFA2" %in% closure$unresolved)

    null_metadata <- publicationGovernanceCopy(lock)
    null_metadata$Packages$FNN$Imports <- list(NULL)
    expect_false(publicationLockMetadataValid(null_metadata))

    bad_date <- publicationGovernanceCopy(lock)
    bad_date$Packages$UpSetR$Date <- "2026-5-20"
    expect_false(publicationLockMetadataValid(bad_date))

    changed_scope <- publicationGovernanceCopy(scope)
    changed_scope$resolved$package_count <- 383L
    expect_error(
        publicationValidateLockScopeAuthority(changed_scope),
        class = "multischolar_publication_lock_error"
    )

    expect_silent(publicationValidateMeasuredNamespaces(
        c("base", "methods"),
        c("base", "methods", "MultiScholaR"),
        scope
    ))
    expect_error(
        publicationValidateMeasuredNamespaces(
            c("base", "methods"),
            c("base", "methods", "MultiScholaR", "V8"),
            scope
        ),
        class = "multischolar_publication_lock_error"
    )
})

test_that("detached comparator worktrees reject revision path and dirty drift", {
    root <- withr::local_tempdir(pattern = "publication-comparators-")
    revision <- publicationComparatorRevisions()[["historical_janitor"]]
    record <- publicationCreateComparatorWorktree(
        "historical_janitor",
        revision,
        root
    )
    withr::defer(publicationRemoveComparatorWorktree(record, root))

    expect_silent(publicationValidateComparatorWorktree(record))

    wrong <- publicationGovernanceCopy(record)
    wrong$revision <- publicationComparatorRevisions()[[
        "pre_repair_performance"
    ]]
    expect_error(
        publicationValidateComparatorWorktree(wrong),
        class = "multischolar_publication_comparator_error"
    )

    moved <- publicationGovernanceCopy(record)
    moved$path <- dirname(record$path)
    expect_error(
        publicationValidateComparatorWorktree(moved),
        class = "multischolar_publication_comparator_error"
    )

    dirty_path <- file.path(record$path, "untracked-comparator-drift")
    writeLines("dirty", dirty_path)
    expect_error(
        publicationValidateComparatorWorktree(record),
        class = "multischolar_publication_comparator_error"
    )
    unlink(dirty_path)
    expect_silent(publicationValidateComparatorWorktree(record))
})

test_that("candidate freeze is additive source-bound and non-promotional", {
    authority <- publicationComparatorAuthorityFixture()
    revision <- publicationComparatorRevisions()[[
        "pre_repair_performance"
    ]]
    digest <- paste(rep("a", 64L), collapse = "")
    bindings <- stats::setNames(
        as.list(rep(digest, length(publicationCandidateFreezeBindingNames()))),
        publicationCandidateFreezeBindingNames()
    )
    readiness_path <- "tests/testdata/omics-performance/protocol-v1.json"
    readiness <- list(
        schema = "multischolar.omics_publication_candidate_freeze_readiness",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-077",
        status = "ready",
        authority_bindings = list(list(
            path = readiness_path,
            sha256 = publicationFileDigest(readiness_path),
            available = TRUE
        )),
        blockers = list(),
        candidate_freeze_allowed = TRUE,
        promotion_authority = FALSE
    )
    readiness$readiness_digest <- publicationObjectDigest(readiness)
    successor <- publicationComparatorFreezeCandidate(
        authority,
        revision,
        publicationComparatorSourceIdentity(revision),
        digest,
        bindings,
        readiness
    )

    expect_silent(publicationValidateCandidateFreezeSuccessor(
        successor,
        authority
    ))
    expect_true(successor$frozen_once)
    expect_false(successor$promotion_authority)

    drift <- publicationGovernanceCopy(successor)
    drift$candidate_revision <- paste(rep("0", 40L), collapse = "")
    expect_error(
        publicationValidateCandidateFreezeSuccessor(drift, authority),
        class = "multischolar_publication_comparator_error"
    )

    source_drift <- publicationGovernanceCopy(successor)
    source_drift$source_identity$tree_oid <- paste(rep("0", 40L), collapse = "")
    expect_error(
        publicationValidateCandidateFreezeSuccessor(source_drift, authority),
        class = "multischolar_publication_comparator_error"
    )

    for (binding in publicationCandidateFreezeBindingNames()) {
        changed <- publicationGovernanceCopy(bindings)
        changed[[binding]] <- paste(rep("b", 64L), collapse = "")
        expect_error(
            publicationValidateCandidateFreezeSuccessor(
                successor,
                authority,
                changed
            ),
            class = "multischolar_publication_comparator_error",
            info = binding
        )
    }

    expect_error(
        publicationComparatorFreezeCandidate(
            successor,
            revision,
            publicationComparatorSourceIdentity(revision),
            digest,
            bindings,
            readiness
        ),
        class = "multischolar_publication_comparator_error"
    )
})

test_that("installed inventories and build receipts compare exact identities", {
    first_root <- withr::local_tempdir(pattern = "installed-first-")
    second_root <- withr::local_tempdir(pattern = "installed-second-")
    dir.create(file.path(first_root, "libs"))
    dir.create(file.path(second_root, "libs"))
    writeLines("same", file.path(first_root, "R-code"))
    writeLines("same", file.path(second_root, "R-code"))
    writeBin(as.raw(1:10), file.path(first_root, "libs", "package.so"))
    writeBin(as.raw(1:10), file.path(second_root, "libs", "package.so"))

    first_inventory <- publicationInstalledInventory(first_root)
    second_inventory <- publicationInstalledInventory(second_root)
    expect_identical(
        first_inventory$inventory_sha256,
        second_inventory$inventory_sha256
    )

    first <- publicationComparatorReceiptFixture(receipt_id = "first")
    second <- publicationComparatorReceiptFixture(receipt_id = "second")
    comparison <- publicationCompareBuildReceipts(first, second)
    expect_true(comparison$reproducible)

    second$build$installed_inventory$inventory_sha256 <- paste(
        rep("b", 64L),
        collapse = ""
    )
    expect_false(publicationCompareBuildReceipts(first, second)$reproducible)

    second <- publicationComparatorReceiptFixture(receipt_id = "second-native")
    second$build$native_libraries[[1L]]$sha256 <- paste(
        rep("c", 64L),
        collapse = ""
    )
    expect_false(publicationCompareBuildReceipts(first, second)$reproducible)
})

test_that("restore roots environments and external archives are isolated", {
    root <- withr::local_tempdir(pattern = "publication-restore-")
    first <- publicationRestorePaths(root, "restore-1")
    second <- publicationRestorePaths(root, "restore-2")
    expect_false(identical(first$library, second$library))
    expect_false(identical(first$cache, second$cache))
    expect_false(identical(first$source_cache, second$source_cache))

    publicationInitializeRestorePaths(first)
    expect_true(all(dir.exists(unlist(first[setdiff(
        names(first),
        c("root", "audit")
    )], use.names = FALSE))))
    expect_error(
        publicationInitializeRestorePaths(first),
        class = "multischolar_publication_comparator_build_error"
    )

    environment <- publicationBuildEnvironment(
        first$bootstrap_library,
        first$site_library
    )
    expect_identical(unname(environment[["R_LIBS_USER"]]),
                     first$bootstrap_library)
    expect_identical(unname(environment[["R_LIBS_SITE"]]),
                     first$site_library)
    expect_identical(unname(environment[["GITHUB_PAT"]]), "")
    expect_identical(unname(environment[["MAKEFLAGS"]]), "-j1")
    restore_environment <- publicationRestoreEnvironment(
        first,
        list(
            include_dir = "/tmp/include",
            library_dir = "/tmp/lib"
        ),
        list(
            root = "/tmp",
            archive = list(path = "libarrow.zip")
        )
    )
    expect_identical(
        unname(restore_environment[["RENV_CONFIG_INSTALL_JOBS"]]),
        "1"
    )
    expect_identical(
        unname(restore_environment[["ARROW_OFFLINE_BUILD"]]),
        "true"
    )
    expect_identical(
        unname(restore_environment[["LIBARROW_BUILD"]]),
        "false"
    )
    expect_error(
        publicationBuildEnvironment(first$library),
        class = "multischolar_publication_comparator_build_error"
    )

    expect_silent(publicationValidateExternalArchiveEntries(c(
        "v8/", "v8/include/v8.h", "v8/lib/libv8_monolith.a"
    )))
    expect_error(
        publicationValidateExternalArchiveEntries("../outside"),
        class = "multischolar_publication_comparator_build_error"
    )
    expect_error(
        publicationValidateExternalArchiveEntries("/absolute/path"),
        class = "multischolar_publication_comparator_build_error"
    )

    timeout_root <- file.path(root, "timeout-probe")
    dir.create(timeout_root)
    timeout <- publicationBuildRun(
        "sleep",
        "1",
        timeout_root,
        file.path(timeout_root, "logs"),
        timeout_seconds = 0.05
    )
    expect_true(timeout$timed_out)
    expect_false(identical(timeout$exit_status, 0L))
})

test_that("common environment incompatibility is terminal and cannot substitute", {
    outcome <- publicationComparatorIncomparableReceipt(
        "historical dependency cannot use the common lock",
        list(
            lockfile_sha256 = publicationFileDigest("renv.lock"),
            failed_package = "fixture"
        )
    )
    expect_silent(publicationValidateComparatorIncomparable(outcome))
    expect_false(outcome$historical_comparison_authority)
    expect_false(outcome$sensitivity_substitution_allowed)
    expect_false(outcome$promotion_authority)

    drift <- publicationGovernanceCopy(outcome)
    drift$sensitivity_substitution_allowed <- TRUE
    expect_error(
        publicationValidateComparatorIncomparable(drift),
        class = "multischolar_publication_comparator_error"
    )
})

test_that("tiny comparator envelopes allow only verified frozen builds", {
    authority <- publicationComparatorAuthorityFixture()

    historical <- publicationComparatorEnvelopeDecision(
        "historical_janitor",
        "memory",
        authority
    )
    expect_identical(historical$status, "allowed")
    expect_true(historical$runner_invoked)

    pre_repair <- publicationComparatorEnvelopeDecision(
        "pre_repair_performance",
        "artifact",
        authority
    )
    expect_identical(pre_repair$status, "allowed")

    historical_artifact <- publicationComparatorEnvelopeDecision(
        "historical_janitor",
        "artifact",
        authority
    )
    expect_identical(historical_artifact$status, "rejected")
    expect_false(historical_artifact$runner_invoked)

    for (backend in c("memory", "artifact")) {
        candidate <- publicationComparatorEnvelopeDecision(
            "candidate",
            backend,
            authority
        )
        expect_identical(candidate$status, "rejected")
        expect_identical(
            candidate$reason,
            "candidate_pending_omics_art_077"
        )
        expect_false(candidate$runner_invoked)
    }
})

test_that("comparator envelope summaries require every governed lane", {
    authority <- publicationComparatorAuthorityFixture()
    entries <- lapply(
        publicationComparatorEnvelopeExpectedKeys(),
        \(key) {
            parts <- strsplit(key, "::", fixed = TRUE)[[1L]]
            decision <- publicationComparatorEnvelopeDecision(
                parts[[1L]],
                parts[[2L]],
                authority
            )
            if (identical(decision$status, "allowed")) {
                return(list(
                    comparator_id = parts[[1L]],
                    backend = parts[[2L]],
                    status = "passed"
                ))
            }
            c(
                list(
                    schema = "multischolar.omics_publication_comparator_envelope",
                    schema_version = "1.0.0",
                    comparator_id = parts[[1L]],
                    backend = parts[[2L]]
                ),
                decision,
                list(publication_authority = FALSE)
            )
        }
    )
    expect_identical(length(entries), 7L)
    expect_identical(
        vapply(entries, `[[`, character(1), "comparator_id"),
        c(
            "historical_janitor", "historical_janitor",
            "pre_repair_performance", "pre_repair_performance",
            "candidate", "candidate", "candidate"
        )
    )
    missing <- entries[-length(entries)]
    expect_error(
        publicationCreateComparatorEnvelopeSummary(missing, authority),
        class = "multischolar_publication_comparator_error"
    )
})

test_that("restore normalization is narrow and build attempts are disjoint", {
    first_path <- paste0(
        "/home/user/Projects/APAF/MultiScholaR/",
        ".omics-publication-comparators/restores/restore-1/",
        "tmp/RtmpABC/R.INSTALL123/package/src/file.c"
    )
    second_path <- paste0(
        "/home/user/Projects/APAF/MultiScholaR/",
        ".omics-publication-comparators/restores/restore-2/",
        "tmp/RtmpXYZ/R.INSTALL999/package/src/file.c"
    )
    expect_identical(
        publicationCanonicalBuildString(first_path),
        publicationCanonicalBuildString(second_path)
    )
    expect_identical(
        publicationCanonicalBuildString(
            "TBB: BUILD_DATE\t\tThu Aug 27 12:36:19 UTC 2026"
        ),
        "TBB: BUILD_DATE\t\t<SOURCE_DATE_EPOCH>"
    )
    expect_identical(
        publicationInstalledDifferenceCategory("inst/data/file.csv", NULL),
        "unexpected"
    )

    first <- publicationComparatorBuildPaths("historical_janitor", 1L, 1L)
    retry <- publicationComparatorBuildPaths("historical_janitor", 1L, 2L)
    second <- publicationComparatorBuildPaths("historical_janitor", 2L, 1L)
    expect_false(identical(first$root, retry$root))
    expect_false(identical(first$root, second$root))
})

test_that("dependency overlays contain only owned exact package links", {
    root <- withr::local_tempdir(pattern = "dependency-overlay-")
    dependency <- file.path(root, "dependency")
    overlay <- file.path(root, "overlay")
    dir.create(dependency)
    for (package in c("pkgA", "pkgB")) {
        path <- file.path(dependency, package)
        dir.create(path)
        writeLines(
            c("Package: ", package, "\nVersion: 1.0.0"),
            file.path(path, "DESCRIPTION")
        )
    }
    lock <- list(Packages = list(
        pkgA = list(Version = "1.0.0"),
        pkgB = list(Version = "1.0.0")
    ))
    record <- publicationPrepareDependencyOverlay(overlay, dependency, lock)

    expect_identical(record$package_count, 2L)
    expect_true(all(file.info(file.path(overlay, c("pkgA", "pkgB")))$isdir))
    expect_identical(
        normalizePath(file.path(overlay, "pkgA")),
        normalizePath(file.path(dependency, "pkgA"))
    )
    expect_error(
        publicationPrepareDependencyOverlay(overlay, dependency, lock),
        class = "multischolar_publication_comparator_build_error"
    )
})

test_that("cleanup archives failed roots without following protected links", {
    sandbox <- withr::local_tempdir(pattern = "comparator-cleanup-")
    withr::local_dir(sandbox)
    root <- ".omics-publication-comparators"
    failed <- file.path(root, "restores", "failed-attempt")
    protected <- file.path(root, "restores", "accepted")
    dir.create(file.path(failed, "logs"), recursive = TRUE)
    dir.create(file.path(failed, "library"), recursive = TRUE)
    dir.create(protected, recursive = TRUE)
    writeLines("failed log", file.path(failed, "logs", "command.log"))
    writeLines("partial package", file.path(failed, "library", "payload"))
    writeLines("accepted", file.path(protected, "authority"))
    expect_true(file.symlink(
        normalizePath(file.path(protected, "authority")),
        file.path(failed, "library", "accepted-link")
    ))
    plan <- list(
        schema = "multischolar.omics_publication_comparator_cleanup_plan",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-064",
        status = "approved_for_execution",
        comparator_root = root,
        evidence_root = "failed-attempt-evidence/test",
        protected_paths = list("restores/accepted"),
        attempts = list(list(
            path = "restores/failed-attempt",
            category = "test_failure",
            reason = "Bounded cleanup fixture"
        )),
        retention_policy = list(
            retain_logs = TRUE,
            retain_scripts = TRUE,
            retain_project_lock = TRUE,
            retain_root_json_r_and_package_archives = TRUE,
            discard_partial_libraries_caches_sources_and_temp = TRUE
        ),
        publication_authority = FALSE
    )

    dry_run <- publicationRunComparatorCleanup(plan, execute = FALSE)
    expect_identical(dry_run$status, "dry_run")
    expect_true(dir.exists(failed))
    expect_identical(dry_run$attempts[[1L]]$before$symlink_count, 1L)

    result <- publicationRunComparatorCleanup(plan, execute = TRUE)
    expect_identical(result$status, "passed")
    expect_silent(publicationValidateComparatorCleanupResult(result, plan))
    expect_false(dir.exists(failed))
    expect_true(file.exists(file.path(protected, "authority")))
    archive <- file.path(
        root,
        plan$evidence_root,
        "restores__failed-attempt--test_failure"
    )
    expect_true(file.exists(file.path(
        archive,
        "retained",
        "logs",
        "command.log"
    )))
    expect_false(file.exists(file.path(
        archive,
        "retained",
        "library",
        "accepted-link"
    )))
})

test_that("scientific difference drift cannot excuse candidate backend mismatch", {
    differences <- publicationReadJson(
        "tests/testdata/omics-performance/scientific-differences-v1.json"
    )
    expect_silent(publicationValidateScientificDifferences(differences))
    expect_identical(
        differences$matched_candidate_rule,
        "exact_science_required"
    )

    drift <- publicationGovernanceCopy(differences)
    drift$matched_candidate_rule <- "allow_known_differences"
    expect_error(
        publicationValidateScientificDifferences(drift),
        class = "multischolar_publication_comparator_error"
    )

    inventory_drift <- publicationGovernanceCopy(differences)
    inventory_drift$inventory$commit_count <- 70L
    expect_error(
        publicationValidateScientificDifferences(inventory_drift),
        class = "multischolar_publication_comparator_error"
    )

    amendment_drift <- publicationGovernanceCopy(differences)
    amendment_drift$candidate_amendment$successor_required <- FALSE
    expect_error(
        publicationValidateScientificDifferences(amendment_drift),
        class = "multischolar_publication_comparator_error"
    )

    evidence_drift <- publicationGovernanceCopy(differences)
    evidence_drift$records[[1L]]$evidence[[1L]] <- paste0(
        "commit:",
        paste(rep("0", 40L), collapse = "")
    )
    expect_error(
        publicationValidateScientificDifferences(evidence_drift),
        class = "multischolar_publication_comparator_error"
    )
})
