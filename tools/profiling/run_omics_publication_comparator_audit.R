#!/usr/bin/env Rscript

publicationComparatorAuditRunnerPath <- function() {
    file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(file_arg)) {
        return(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE))
    }
    normalizePath(
        file.path(
            getwd(),
            "tools",
            "profiling",
            "run_omics_publication_comparator_audit.R"
        ),
        mustWork = TRUE
    )
}

.PUBLICATION_COMPARATOR_AUDIT_ROOT <- normalizePath(
    file.path(dirname(publicationComparatorAuditRunnerPath()), "..", ".."),
    mustWork = TRUE
)

for (source_file in c(
    "omics_publication_protocol.R",
    "omics_publication_comparators.R",
    "omics_publication_lock.R",
    "omics_publication_builds.R",
    "omics_publication_repository_inputs.R",
    "omics_publication_remote_installs.R",
    "omics_publication_restore_reproducibility.R",
    "omics_publication_comparator_builds.R",
    "omics_publication_comparator_build_reproducibility.R",
    "omics_publication_comparator_evidence.R",
    "omics_publication_comparator_cleanup.R",
    "omics_publication_comparator_envelopes.R"
)) {
    source(
        file.path(
            .PUBLICATION_COMPARATOR_AUDIT_ROOT,
            "tools",
            "profiling",
            source_file
        ),
        local = FALSE
    )
}

publicationComparatorAuditPaths <- function() {
    list(
        authority = "tests/testdata/omics-performance/comparators-v1.json",
        scientific_differences = paste0(
            "tests/testdata/omics-performance/",
            "scientific-differences-v1.json"
        ),
        restore = "docs/qa/omics-art-064-restore-reproducibility-v1.json",
        historical = paste0(
            "docs/qa/omics-art-064-",
            "historical-build-reproducibility-v1.json"
        ),
        pre_repair = paste0(
            "docs/qa/omics-art-064-",
            "pre-repair-build-reproducibility-v1.json"
        ),
        envelopes = paste0(
            ".omics-publication-comparators/",
            "comparator-envelopes/summary-v1.json"
        ),
        cleanup_plan = "docs/qa/omics-art-064-cleanup-plan-v1.json",
        cleanup_result = "docs/qa/omics-art-064-cleanup-result-v1.json",
        publication_tests = "docs/qa/omics-art-064-publication-tests.xml",
        runtime_tests = "docs/qa/omics-art-064-runtime-tests.xml",
        package_check = "docs/qa/omics-art-064-r-cmd-check.log",
        package_smoke = "docs/qa/omics-art-064-testthat.Rout",
        kernel = "tools/profiling/run_omics_publication_benchmark.R"
    )
}

publicationComparatorAuditImplementationPaths <- function() {
    c(
        "tools/profiling/omics_publication_comparators.R",
        "tools/profiling/omics_publication_lock.R",
        "tools/profiling/omics_publication_builds.R",
        "tools/profiling/omics_publication_repository_inputs.R",
        "tools/profiling/omics_publication_remote_installs.R",
        "tools/profiling/omics_publication_restore_reproducibility.R",
        "tools/profiling/omics_publication_comparator_builds.R",
        paste0(
            "tools/profiling/",
            "omics_publication_comparator_build_reproducibility.R"
        ),
        "tools/profiling/omics_publication_comparator_evidence.R",
        "tools/profiling/omics_publication_comparator_cleanup.R",
        "tools/profiling/omics_publication_comparator_envelopes.R",
        "tools/profiling/run_omics_publication_comparator.R",
        "tools/profiling/run_omics_publication_comparator_audit.R",
        "tools/profiling/run_omics_publication_comparator_cleanup.R",
        "tests/testthat/test-omics-publication-comparators.R"
    )
}

publicationComparatorAuditBinding <- function(path) {
    list(path = path, sha256 = publicationFileDigest(path))
}

publicationComparatorAuditValidateJunit <- function(
    path,
    expected_tests,
    expected_skipped
) {
    document <- xml2::read_xml(path)
    suites <- xml2::xml_find_all(document, "//testsuite")
    total <- \(field) sum(as.numeric(xml2::xml_attr(suites, field)), na.rm = TRUE)
    valid <- identical(total("tests"), as.numeric(expected_tests)) &&
        identical(total("failures"), 0) && identical(total("errors"), 0) &&
        identical(total("skipped"), as.numeric(expected_skipped))
    if (!valid) publicationEvidenceAbort("JUnit evidence differs")
    list(
        tests = as.integer(total("tests")),
        failures = as.integer(total("failures")),
        errors = as.integer(total("errors")),
        skipped = as.integer(total("skipped"))
    )
}

publicationComparatorAuditValidatePackageCheck <- function(path, smoke_path) {
    log <- readLines(path, warn = FALSE)
    smoke <- readLines(smoke_path, warn = FALSE)
    valid <- any(grepl("checking tests \\.\\.\\. OK", log)) &&
        any(grepl("Status: 1 WARNING, 1 NOTE", log, fixed = TRUE)) &&
        !any(grepl("Status:.*ERROR", log)) &&
        any(grepl("FAIL 0.*WARN 0.*SKIP 0.*PASS 3", smoke))
    if (!valid) publicationEvidenceAbort("package check evidence differs")
    list(
        errors = 0L,
        warnings = 1L,
        notes = 1L,
        installed_smoke_passed = TRUE,
        known_legacy_findings = list(
            "duplicate_Rd_aliases",
            "peptide_value_global_binding_note"
        )
    )
}

publicationComparatorAuditWorktrees <- function() {
    result <- publicationComparatorGit(c("worktree", "list", "--porcelain"))
    if (result$status != 0L) {
        publicationEvidenceAbort("git worktree inventory failed")
    }
    paths <- sub("^worktree ", "", grep(
        "^worktree ",
        strsplit(result$stdout, "\n", fixed = TRUE)[[1L]],
        value = TRUE
    ))
    expected <- normalizePath(.PUBLICATION_COMPARATOR_AUDIT_ROOT)
    valid <- identical(normalizePath(paths), expected) &&
        !length(list.files(
            file.path(publicationComparatorRoot(), "worktrees"),
            all.files = TRUE,
            no.. = TRUE
        ))
    if (!valid) publicationEvidenceAbort("temporary worktrees remain")
    list(
        count = as.integer(length(paths)),
        paths = as.list(paths),
        porcelain_sha256 = digest::digest(
            result$stdout,
            algo = "sha256",
            serialize = FALSE
        )
    )
}

publicationRunComparatorEvidenceAudit <- function() {
    paths <- publicationComparatorAuditPaths()
    authority <- publicationReadJson(paths$authority)
    differences <- publicationReadJson(paths$scientific_differences)
    publicationValidateComparatorAuthority(authority)
    publicationValidateScientificDifferences(differences)
    publicationValidateRestoreReproducibilityEvidence(paths$restore)
    publicationValidateComparatorBuildPairEvidence(paths$historical)
    publicationValidateComparatorBuildPairEvidence(paths$pre_repair)
    publicationValidateComparatorEnvelopeSummary(publicationReadJson(
        paths$envelopes
    ))
    cleanup_plan <- publicationReadJson(paths$cleanup_plan)
    publicationValidateComparatorCleanupResult(
        publicationReadJson(paths$cleanup_result),
        cleanup_plan
    )
    verification <- list(
        publication = publicationComparatorAuditValidateJunit(
            paths$publication_tests,
            421L,
            4L
        ),
        runtime = publicationComparatorAuditValidateJunit(
            paths$runtime_tests,
            1244L,
            0L
        ),
        package_check = publicationComparatorAuditValidatePackageCheck(
            paths$package_check,
            paths$package_smoke
        )
    )
    kernel_sha256 <- publicationFileDigest(paths$kernel)
    if (!identical(
        kernel_sha256,
        "823b951e2e53d43d822b80e21bc28e06e1cfc98149812166de71b1f7207657f2"
    )) {
        publicationEvidenceAbort("measurement kernel differs")
    }
    list(
        schema = "multischolar.omics_publication_comparator_evidence_audit",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-064",
        status = "passed",
        validated_at = publicationUtcNow(),
        authorities = lapply(paths[c(
            "authority", "scientific_differences", "restore", "historical",
            "pre_repair", "envelopes", "cleanup_plan", "cleanup_result",
            "publication_tests", "runtime_tests", "package_check",
            "package_smoke", "kernel"
        )], publicationComparatorAuditBinding),
        implementation_sources = lapply(
            as.list(publicationComparatorAuditImplementationPaths()),
            publicationComparatorAuditBinding
        ),
        verification = verification,
        fixed_revisions = as.list(publicationComparatorRevisions()),
        worktrees = publicationComparatorAuditWorktrees(),
        candidate_status = "pending_omics_art_077",
        candidate_execution_authority = FALSE,
        publication_authority = FALSE
    )
}

publicationComparatorAuditParseArgs <- function(argv) {
    if (length(argv) != 2L || !identical(argv[[1L]], "--output")) {
        stop("Usage: --output <json>", call. = FALSE)
    }
    argv[[2L]]
}

publicationComparatorAuditMain <- function(argv = commandArgs(trailingOnly = TRUE)) {
    output <- publicationComparatorAuditParseArgs(argv)
    if (file.exists(output)) {
        stop("Comparator audit output already exists", call. = FALSE)
    }
    record <- publicationRunComparatorEvidenceAudit()
    publicationWriteJson(record, output)
    cat(publicationFileDigest(output), "\n")
    invisible(0L)
}

if (identical(environment(), globalenv()) && !interactive()) {
    status <- tryCatch(
        publicationComparatorAuditMain(),
        error = \(error) {
            message(conditionMessage(error))
            1L
        }
    )
    quit(status = as.integer(status), save = "no")
}
