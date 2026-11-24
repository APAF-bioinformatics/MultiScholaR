# G1(a)'s job pieces (KIT-159 section 2; KIT-199 12.3, D7 / VAL-120 (2); VAL-133 (1)): the
# job-environment recorder, the lockfile check, gates.md, the entry script's refusal, and L1's
# committed DA tables for G13.

envFixtureRoot <- function() {
    root <- withr::local_tempdir(.local_envir = parent.frame())
    writeLines(c("Package: MultiScholaR", "Remotes:",
                 paste0("    fastduckdb=github::o/fast_duckdb/fastduckdb@", strrep("a", 40), ","),
                 paste0("    fastlimpadpc=github::o/fast_limpa_dpc@", strrep("b", 40))),
               file.path(root, "DESCRIPTION"))
    dir.create(file.path(root, "packaging"))
    writeLines("{\"R\": {}}", file.path(root, "packaging", "multischolar.lock.json"))
    writeLines(c("processor\t: 0", "model name\t: AMD EPYC 9V74 80-Core Processor"),
               file.path(root, "cpuinfo"))
    root
}

fakeCommand <- function(command, args, pattern = NULL) {
    switch(basename(command), getconf = "2.35", `dpkg-query` = "3.10.0-2ubuntu1",
           "1.100.0")
}

test_that("T-G1A-ENV-1: the job environment holds the CPU model, the D7 keys and the lockfile digest", {
    root <- envFixtureRoot()
    withr::local_envvar(ImageOS = "ubuntu22", ImageVersion = "20261116.365.1")
    env <- valJobEnvironment(root, cpuinfo = file.path(root, "cpuinfo"), command = fakeCommand)
    expect_identical(names(env), msrval:::G1A_JOB_ENV_KEYS)
    expect_true(all(msrval:::G1A_D7_ENV_KEYS %in% names(env)))
    expect_identical(env[["cpu_model"]], "AMD EPYC 9V74 80-Core Processor")
    expect_identical(env[["image_version"]], "20261116.365.1")
    expect_identical(env[["c_library_version"]], "2.35")
    expect_identical(env[["pin_fastduckdb"]], strrep("a", 40))
    expect_identical(env[["pin_fastlimpadpc"]], strrep("b", 40))
    expect_identical(env[["r_version"]], R.version.string)
    expect_identical(env[["lockfile_sha256"]],
                     sha256File(file.path(root, "packaging", "multischolar.lock.json")))
    path <- file.path(root, "environment.tsv")
    valWriteJobEnvironment(env, path)
    back <- msrval:::valReadEnvironment(path)
    expect_identical(unname(back[names(env)]), unname(env))
    # The recorded values decide D7's mode against a link with the same values.
    link <- list(raw_agreed = TRUE, cpus = env[["cpu_model"]],
                 env = env[msrval:::G1A_D7_ENV_KEYS])
    expect_identical(valG1aDaMode(back, link, lock_unchanged = TRUE)$mode, "exact")
    link$env[["image_version"]] <- "20261201.1.1"
    expect_identical(valG1aDaMode(back, link, lock_unchanged = TRUE)$mode, "swd083")
})

test_that("T-G1A-ENV-2: the lockfile is compared between the link's commit and the job's", {
    repo <- gitRepo()
    a <- gitCommitFiles(repo, list(`packaging/multischolar.lock.json` = "{\"v\": 1}"), "lock")
    b <- gitCommitFiles(repo, list(`R/g.R` = "g <- function() 1"), "code")
    expect_true(valLockfileUnchanged(repo, a, b))
    c1 <- gitCommitFiles(repo, list(`packaging/multischolar.lock.json` = "{\"v\": 2}"), "lock 2")
    expect_false(valLockfileUnchanged(repo, a, c1))
})

test_that("T-G1A-GATES: gates.md carries the environment, part (1), part (2) and the D7 mode", {
    env <- c(cpu_model = "AMD EPYC 7763 64-Core Processor", image_version = "x")
    p1 <- data.frame(lane = c("lipid_canonical", "metab_lc"), status = c("identical", "differs"),
                     first_step = c(NA, "norm"))
    p2 <- data.frame(lane = "lipid_canonical", status = "certified", reason = "", chain = "KIT-172",
                     differences = 0L, classes = "")
    da <- list(link_commit = strrep("c", 40),
               mode = list(mode = "swd083", reasons = "another CPU model"))
    path <- tempfile(fileext = ".md")
    valG1aGatesMd(path, strrep("d", 40), env, p1, p2, da)
    txt <- readLines(path)
    expect_true(any(grepl("AMD EPYC 7763", txt, fixed = TRUE)))
    expect_true(any(grepl("| metab_lc | differs | norm |", txt, fixed = TRUE)))
    expect_true(any(grepl("G1(a) stops at its first differing step", txt, fixed = TRUE)))
    expect_true(any(grepl("| lipid_canonical | certified |", txt, fixed = TRUE)))
    expect_true(any(grepl("mode `swd083` (another CPU model)", txt, fixed = TRUE)))
    valG1aGatesMd(path, strrep("d", 40), env, p1[1L, ])
    txt <- readLines(path)
    expect_true(any(grepl("every compared step identical", txt, fixed = TRUE)))
    expect_true(any(grepl("Not run in this job.", txt, fixed = TRUE)))
})

test_that("T-G1A-ENTRY: val-g1a.R refuses part (2) without a reference table, with its class", {
    script <- file.path(Sys.getenv("VAL_HARNESS_DIR"), "val-g1a.R")
    skip_if(!nzchar(Sys.getenv("VAL_HARNESS_DIR")), "run through val-tests.R")
    expect_true(file.exists(script))
    rt <- childRunnerTemp()
    skip_if(is.na(rt), "msrval is not installed in a val-lib (run through val-tests.R)")
    root <- driverPkgRoot()
    out <- file.path(withr::local_tempdir(), "g1a")
    res <- runEntryScript(script, c(out, "--lanes", "lipid_canonical", "--parts", "2"), root, rt)
    expect_identical(res$status, 3L)
    expect_true(any(grepl("phase=g1a status=refused class=val_spec_invalid",
                          statusLine(res$lines))))
    # The environment is recorded before anything else.
    expect_true(file.exists(file.path(out, "environment.tsv")))
})

test_that("T-PAR-5: L1's committed DA tables are read for G13 and checked against their hash", {
    repo <- gitRepo()
    tab <- "Protein.Ids\tcomparison\tfdr_qvalue\tlog2FC\nP1\tKO-WT\t0.01\t1.5"
    base <- "docs/qa/workflow-execution/kit/kit-159/lifecycle"
    gitCommitFiles(repo, stats::setNames(list(tab), file.path(base, "l1-da-prot_lfq.tsv")), "t")
    sha <- sha256File(file.path(repo, base, "l1-da-prot_lfq.tsv"))
    entry <- list(digest = strrep("d", 64L), fingerprint = list(sha256_10 = strrep("e", 64L)),
                  table = file.path(base, "l1-da-prot_lfq.tsv"), table_sha256 = sha,
                  keys = list("Protein.Ids", "comparison"), q_col = "fdr_qvalue",
                  fc_col = "log2FC")
    path <- file.path(base, "l1-digests.json")
    dir.create(file.path(repo, dirname(path)), recursive = TRUE, showWarnings = FALSE)
    jsonlite::write_json(list(prot_lfq = entry), file.path(repo, path), auto_unbox = TRUE)
    commit <- gitCommitFiles(repo, list(), "l1 evidence")
    req <- list(release_commit = commit, l1_path = path,
                l1 = list(prot_lfq = entry[c("digest", "fingerprint")]))
    got <- valParityL1Tables(req, repo, at = commit)
    expect_identical(got$prot_lfq$keys, c("Protein.Ids", "comparison"))
    expect_match(got$prot_lfq$table, "P1\tKO-WT", fixed = TRUE)
    bad <- entry
    bad$table_sha256 <- strrep("0", 64L)
    jsonlite::write_json(list(prot_lfq = bad), file.path(repo, path), auto_unbox = TRUE)
    c2 <- gitCommitFiles(repo, list(), "wrong hash")
    req$release_commit <- c2
    expect_error(valParityL1Tables(req, repo, at = c2), class = "val_spec_invalid")
    nofc <- entry[setdiff(names(entry), "fc_col")]
    jsonlite::write_json(list(prot_lfq = nofc), file.path(repo, path), auto_unbox = TRUE)
    c3 <- gitCommitFiles(repo, list(), "no fold-change column")
    expect_error(valParityL1Tables(req, repo, at = c3), class = "val_spec_invalid")
})
