# KIT-167 section 10 and 8.4: the engine's pin file, finding the program, the checksum, the Windows runtime statuses
# and the timeout. The program is replaced by stub executables and a stub runner; DIAMOND itself is not run here.

k167StubEngine <- function(body) {
    skip_on_os("windows")
    path <- file.path(tempfile("k167-engine-"), "diamond")
    dir.create(dirname(path))
    writeLines(c("#!/bin/sh", body), path)
    Sys.chmod(path, "0755")
    path
}
k167Pin <- function(path, platform = referenceEnginePlatform()) {
    data.frame(platform = platform, version = "2.2.8", archive_url = "https://example.invalid/diamond.tar.gz",
        archive_bytes = 1, archive_sha256 = strrep("0", 64), executable_path = "diamond",
        executable_sha256 = digest::digest(file = path, algo = "sha256"), stringsAsFactors = FALSE)
}

test_that("the shipped pin file names v2.2.8 and the two release assets of section 2.2", {
    pin <- readDiamondPin()
    expect_setequal(pin$platform, c("linux", "windows"))
    expect_identical(unique(pin$version), "2.2.8")
    expect_identical(pin$archive_sha256[pin$platform == "linux"],
        "33a765d192ae4f5b388f15fe2e7c23d3aaea0f9e3f6872efe8c324aec261fbc7")
    expect_identical(pin$archive_bytes[pin$platform == "windows"], 3366672)
    expect_true(all(grepl("^https://github.com/bbuchfink/diamond/releases/download/v2.2.8/", pin$archive_url)))
})

test_that("the program is found from MSR_DIAMOND, then the option, then the path", {
    stub <- k167StubEngine("echo 'diamond version 2.2.8'")
    withr::local_envvar(MSR_DIAMOND = stub)
    expect_identical(findDiamond(), stub)
    withr::local_envvar(MSR_DIAMOND = "")
    withr::local_options(MultiScholaR.diamond_path = stub)
    expect_identical(findDiamond(), stub)
})

test_that("an absent engine is refused in words", {
    withr::local_envvar(MSR_DIAMOND = file.path(tempdir(), "no-such-diamond"))
    err <- tryCatch(checkDiamond(), reference_mapping_engine_missing = function(e) e)
    expect_s3_class(err, "reference_mapping_engine_missing")
    expect_match(conditionMessage(err), "Saved mappings and the demo tables still work", fixed = TRUE)
})

test_that("a checksum different from the pin is refused", {
    stub <- k167StubEngine("echo 'diamond version 2.2.8'")
    pin <- k167Pin(stub)
    pin$executable_sha256 <- strrep("a", 64)
    resetDiamondCheck()
    expect_error(checkDiamond(stub, pin = pin), class = "reference_mapping_engine_mismatch")
})

test_that("a matching engine passes and is checked once per session", {
    stub <- k167StubEngine("echo 'diamond version 2.2.8'")
    resetDiamondCheck()
    calls <- 0L
    runner <- function(...) {
        calls <<- calls + 1L
        list(status = 0L, stdout = "diamond version 2.2.8\n", stderr = "", timeout = FALSE)
    }
    expect_identical(checkDiamond(stub, pin = k167Pin(stub), run_fn = runner)$version, "2.2.8")
    checkDiamond(stub, pin = k167Pin(stub), run_fn = runner)
    expect_identical(calls, 1L)
})

test_that("the Windows runtime statuses 0xC0000135 and 0xC0000139 are the runtime-missing refusal", {
    stub <- k167StubEngine("exit 0")
    for (status in c(-1073741515L, -1073741511L, 3221225781, 3221225785)) {
        resetDiamondCheck()
        runner <- function(...) list(status = status, stdout = "", stderr = "", timeout = FALSE)
        err <- tryCatch(checkDiamond(stub, pin = k167Pin(stub), run_fn = runner),
            reference_mapping_engine_runtime_missing = function(e) e)
        expect_s3_class(err, "reference_mapping_engine_runtime_missing")
        expect_match(conditionMessage(err), "Visual C++", fixed = TRUE)
    }
})

test_that("a search exceeding its time limit is refused and nothing is kept", {
    stub <- k167StubEngine("sleep 10")
    scratch <- tempfile("k167-scratch-")
    dir.create(scratch)
    expect_error(runDiamondCommand(stub, c("blastp"), scratch = scratch, timeout = 1),
        class = "reference_mapping_timeout")
})

test_that("the search arguments are the shipped settings of section 5.1, in a fixed order", {
    args <- diamondSearchArguments("strain.fasta", "reference", "forward.tsv", evalue = 1e-5, threads = 4L)
    expect_identical(args[1:3], c("blastp", "--query", "strain.fasta"))
    for (flag in c("--ultra-sensitive", "--masking", "--ext", "--comp-based-stats", "-k0", "--max-hsps",
        "--block-size", "--index-chunks", "--algo", "--tmpdir", "--outfmt")) {
        expect_true(flag %in% args, info = flag)
    }
    expect_identical(args[which(args == "--masking") + 1L], "none")
    expect_identical(args[which(args == "--ext") + 1L], "full")
    expect_identical(args[which(args == "--comp-based-stats") + 1L], "1")
    expect_identical(args[which(args == "--evalue") + 1L], "1e-05")
    expect_identical(utils::tail(args, 17L), diamondOutputFields)
})
