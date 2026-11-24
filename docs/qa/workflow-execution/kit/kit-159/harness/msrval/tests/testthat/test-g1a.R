# G1(a) against the real committed control records (KIT-159 section 2, I15 remaining; KIT-199
# 12.3; VAL-120 (2), VAL-133 (2)): the loader of `gui-entry-control.tsv` and `declared-keys.tsv`,
# link (i)'s value check, reference-row adoption by key sets, the release-commit value comparison
# and the part (1) and (2) runners. Fixtures in `g1a/`, copied with their source commit
# (`g1a/SOURCES.tsv`).

g1aFx <- function(name) testthat::test_path("g1a", name)

g1aLaneSteps <- function() {
    t <- utils::read.delim(g1aFx("lane-steps.tsv"), colClasses = "character")
    split(t$step, t$lane)[unique(t$lane)]
}

g1aLink <- function(ticket, lane, keys = TRUE) {
    n <- sub("KIT-", "k", ticket, fixed = TRUE)
    valReadControlLink(g1aFx(paste0(n, "-ctl.tsv")),
                       declared = if (keys) g1aFx(paste0(n, "-keys.tsv")) else NULL,
                       ticket = ticket, lane = lane, lane_steps = g1aLaneSteps()[[lane]])
}

test_that("G1A-LOAD: KIT-199's real record is one certified link of prot_dia_limpa", {
    link <- g1aLink("KIT-199", "prot_dia_limpa")
    expect_identical(link$changed_step, "da")
    expect_identical(nrow(link$base_diffs), 34L)
    expect_identical(sum(link$base_diffs$class == "declared_change (KIT-199)"), 6L)
    expect_identical(nrow(link$before_diffs), 0L)
    expect_identical(nrow(link$after_diffs), 0L)
    expect_identical(nrow(link$declared), 6L)
    expect_true(all(c("base", "live") %in% names(link$base_diffs)))
    ref <- valReadReferenceRows(g1aFx("k199-ref.tsv"))
    expect_identical(valReadChain("prot_dia_limpa", list(link), ref)$status, "certified")
})

test_that("G1A-VAL: link (i) checks each declared row's stored and live digest values", {
    ref <- valReadReferenceRows(g1aFx("k199-ref.tsv"))
    link <- g1aLink("KIT-199", "prot_dia_limpa")
    at <- which(link$base_diffs$class == "declared_change (KIT-199)")[[1L]]
    wrong_live <- link
    wrong_live$base_diffs$live[[at]] <- strrep("e", 64L)
    res <- valReadChain("prot_dia_limpa", list(wrong_live), ref)
    expect_identical(res$status, "not_certified")
    expect_match(res$reason, "(i)", fixed = TRUE)
    wrong_stored <- link
    wrong_stored$base_diffs$base[[at]] <- strrep("e", 64L)
    expect_identical(valReadChain("prot_dia_limpa", list(wrong_stored), ref)$status,
                     "not_certified")
    no_frozen <- link
    no_frozen$declared$frozen <- NULL
    expect_identical(valReadChain("prot_dia_limpa", list(no_frozen), ref)$status,
                     "not_certified")
})

test_that("G1A-LOAD: KIT-172's real record certifies each of its four lanes", {
    ref <- valReadReferenceRows(g1aFx("k172-ref.tsv"))
    for (lane in c("lipid_canonical", "metab_lc", "metab_gc", "metab_combined")) {
        link <- g1aLink("KIT-172", lane)
        expect_identical(link$changed_step, "norm", info = lane)
        expect_gt(nrow(link$before_diffs), 0L)
        expect_true(all(link$before_diffs$class == "declared_change (KIT-172)"), info = lane)
        expect_identical(valReadChain(lane, list(link), ref)$status, "certified", info = lane)
    }
    # a head value before the changed step other than the declared one fails (ii)
    link <- g1aLink("KIT-172", "lipid_canonical")
    link$before_diffs$head_value[[1L]] <- "<ui> 0000"
    res <- valReadChain("lipid_canonical", list(link), ref)
    expect_identical(res$status, "not_certified")
    expect_match(res$reason, "(ii)", fixed = TRUE)
})

test_that("G1A-LOAD: KIT-180 then KIT-177 on the proteomics lanes, in merge order", {
    ref <- valReadReferenceRows(g1aFx("k180-ref.tsv"))
    for (lane in c("prot_tmt", "prot_lfq", "prot_lfq_fragpipe")) {
        l180 <- g1aLink("KIT-180", lane, keys = FALSE)
        l177 <- g1aLink("KIT-177", lane, keys = FALSE)
        expect_identical(l180$changed_step, "design", info = lane)
        expect_identical(nrow(l177$base_diffs), 0L)
        expect_gt(nrow(l177$after_diffs), 0L)
        expect_identical(valReadChain(lane, list(l180, l177), ref)$status, "certified",
                         info = lane)
    }
    l177 <- g1aLink("KIT-177", "prot_tmt", keys = FALSE)
    l177$after_diffs$class[[1L]] <- "unexplained"
    res <- valReadChain("prot_tmt", list(g1aLink("KIT-180", "prot_tmt", keys = FALSE), l177), ref)
    expect_identical(res$status, "not_certified")
    expect_identical(res$failed_link, 2L)
})

test_that("G1A-ADOPT: reference rows are adopted by (lane, step, class) and key sets", {
    earlier <- valReadReferenceRows(g1aFx("k180-ref.tsv"))
    keys <- vapply(earlier$keyset, paste, character(1), collapse = ";")
    own <- data.frame(lane = earlier$lane, step = earlier$step, class = earlier$class,
                      fields = "artifacts", keys = keys, stringsAsFactors = FALSE)
    own <- rbind(own, data.frame(lane = "trx_screen", step = "design",
                                 class = "kit128_technical_replicates", fields = "states",
                                 keys = "states/raw"))
    expect_identical(valAdoptReferenceRows(own, earlier), own)
    other <- own
    other$keys[[1L]] <- sub("samples.tsv", "samples2.tsv", other$keys[[1L]], fixed = TRUE)
    expect_error(valAdoptReferenceRows(other, earlier), class = "val_spec_invalid")
    # KIT-199's two-class steps, its own column form
    k199 <- valReadReferenceRows(g1aFx("k199-ref.tsv"))
    own199 <- data.frame(lane = k199$lane, step = k199$step, class = k199$class,
                         fields = "artifacts",
                         keys = vapply(k199$keyset, paste, character(1), collapse = ";"))
    expect_identical(valAdoptReferenceRows(own199, k199), own199)
    # KIT-172's rows list truncated key lists: its keys must be among G1(r)'s own
    k172 <- valReadReferenceRows(g1aFx("k172-ref.tsv"))
    expect_true(any(k172$truncated))
    own172 <- data.frame(lane = k172$lane, step = k172$step, class = k172$class,
                         fields = "artifacts",
                         keys = vapply(seq_len(nrow(k172)), function(i) {
                             paste(c(k172$keyset[[i]], if (k172$truncated[[i]]) "artifacts/more"),
                                   collapse = ";")
                         }, character(1)))
    expect_identical(valAdoptReferenceRows(own172, k172), own172)
})

test_that("G1A-D7: the release-commit DA values against the link's da-runner table", {
    r1 <- valReadTextTable(g1aFx("k199-da-r1.tsv"))
    r4 <- valReadTextTable(g1aFx("k199-da-r4.tsv"))
    same <- valCompareDaValues(r1, r1)
    expect_identical(same$verdict, "identical")
    near <- valCompareDaValues(r4, r1)
    expect_identical(near$verdict, "agree")
    expect_identical(near$calls_changed, 0L)
    lfc <- near$columns[near$columns$column == "log2FC", ]
    expect_equal(lfc$max_rel, 5.2273398432987143e-11, tolerance = 1e-6)
    expect_identical(valCompareDaValues(r4, r1, mode = "exact")$verdict, "outside")
    called <- r4
    called$fdr_qvalue[called$Protein.Ids == "P00038"] <- "0.06"
    expect_identical(valCompareDaValues(called, r1)$verdict, "outside")
    far <- r4
    far$log2FC[[1L]] <- "3.2"
    expect_identical(valCompareDaValues(far, r1)$verdict, "outside")
    swapped <- r4[c(2L, 1L, 3:5), ]
    expect_identical(valCompareDaValues(swapped, r1)$verdict, "outside")
    missing <- r4
    missing$raw.ko_1.KO[[2L]] <- "NA"
    expect_identical(valCompareDaValues(missing, r1)$verdict, "outside")
})

test_that("G1A-D7: exact only on the link's CPU models with its whole recorded environment", {
    link <- valReadDaLink(g1aFx("k199-da-r1.tsv"),
                          envs = c(g1aFx("k199-env1.tsv"), g1aFx("k199-env2.tsv")),
                          controls = c(g1aFx("k199-da-r1.tsv"), g1aFx("k199-da-r1.tsv")))
    expect_setequal(link$cpus, c("AMD EPYC 9V74 80-Core Processor",
                                 "AMD EPYC 7763 64-Core Processor"))
    job <- stats::setNames(link$env, names(link$env))
    job[["cpu_model"]] <- "AMD EPYC 7763 64-Core Processor"
    expect_identical(valG1aDaMode(job, link, lock_unchanged = TRUE)$mode, "exact")
    expect_identical(valG1aDaMode(job, link, lock_unchanged = FALSE)$mode, "swd083")
    intel <- job
    intel[["cpu_model"]] <- "INTEL(R) XEON(R) PLATINUM 8573C"
    expect_identical(valG1aDaMode(intel, link, lock_unchanged = TRUE)$mode, "swd083")
    image <- job
    image[["image_version"]] <- "20261201.1.1"
    expect_identical(valG1aDaMode(image, link, lock_unchanged = TRUE)$mode, "swd083")
    split_link <- valReadDaLink(g1aFx("k199-da-r1.tsv"),
                                envs = c(g1aFx("k199-env1.tsv"), g1aFx("k199-env2.tsv")),
                                controls = c(g1aFx("k199-da-r1.tsv"), g1aFx("k199-da-r4.tsv")))
    expect_false(split_link$raw_agreed)
    expect_identical(valG1aDaMode(job, split_link, lock_unchanged = TRUE)$mode, "swd083")
})

test_that("G1A-REL: VAL-133 (2): identical captures first, else values; certified with the difference", {
    r1 <- valReadTextTable(g1aFx("k199-da-r1.tsv"))
    r4 <- valReadTextTable(g1aFx("k199-da-r4.tsv"))
    link <- list(table = r1)
    none <- data.frame(lane = character(), step = character(), field = character(),
                       key = character(), base = character(), live = character())
    a <- valG1aReleaseValues(none, r1, link, mode = "exact")
    expect_identical(a$status, "certified")
    expect_true(a$masked_identical)
    da_only <- data.frame(lane = "prot_dia_limpa", step = "da", field = "artifacts",
                          key = msrval:::G1A_DA_KEY, base = "x", live = "y")
    b <- valG1aReleaseValues(da_only, r4, link, mode = "swd083")
    expect_identical(b$status, "certified")
    expect_false(b$masked_identical)
    expect_gt(nrow(b$difference), 0L)
    expect_identical(valG1aReleaseValues(da_only, r4, link, mode = "exact")$status,
                     "not_certified")
    other <- rbind(da_only, data.frame(lane = "prot_dia_limpa", step = "norm", field = "states",
                                       key = "normalized", base = "x", live = "y"))
    expect_identical(valG1aReleaseValues(other, r4, link, mode = "swd083")$status,
                     "not_certified")
    called <- r4
    called$fdr_qvalue[called$Protein.Ids == "P00038"] <- "0.06"
    expect_identical(valG1aReleaseValues(none, called, link, mode = "swd083")$status,
                     "not_certified")
})

test_that("G1A-P1: part (1) compares installed against load_all; prot_dia stops at da", {
    s <- list(import = list(states = list(a = "1")), da = list(states = list(a = "2")),
              enrich = list(artifacts = list(e = "3")))
    t <- s
    t$enrich$artifacts$e <- "4"
    expect_identical(valG1aCompareModes("prot_dia", s, t)$status, "identical")
    expect_identical(valG1aCompareModes("prot_lfq", s, t)$status, "differs")
    expect_identical(valG1aCompareModes("prot_lfq", s, t)$first_step, "enrich")
    runs <- list()
    fake <- function(pkg_root, lane, mode, out_dir, copy_da = NULL) {
        runs[[length(runs) + 1L]] <<- c(lane, mode)
        steps <- if (identical(lane, "metab_lc") && identical(mode, "load_all")) t else s
        path <- file.path(out_dir, paste0(lane, "-", mode, ".json"))
        jsonlite::write_json(list(lane = lane, steps = steps), path, auto_unbox = TRUE)
        path
    }
    out <- withr::local_tempdir()
    res <- valG1aPart1("pkg", c("lipid_canonical", "metab_lc"), out, run_lane = fake)
    expect_identical(res$status, c("identical", "differs"))
    expect_identical(length(runs), 4L)
    expect_true(file.exists(file.path(out, "g1a-part1.tsv")))
})

test_that("G1A-P2: part (2) reads each lane's chain from the history and judges the live run", {
    root <- driverPkgRoot()
    skip_if(is.na(root), "VAL_HARNESS_DIR not set")
    compare_lane <- valLaneComparator(root)
    repo <- gitRepo()
    oracle <- function(steps) list(lane = "lipid_canonical", steps = steps)
    old <- list(import = list(states = list(raw = "a")), norm = list(states = list(n = "b")))
    new <- list(import = list(states = list(raw = "a")), norm = list(states = list(n = "c")))
    path <- "tests/testdata/gui-entry/lipid_canonical.json"
    put <- function(x) jsonlite::toJSON(x, auto_unbox = TRUE, pretty = TRUE)
    since <- gitCommitFiles(repo, stats::setNames(list(put(oracle(old))), path), "oracle")
    ctl <- c("comparison\tlane\tstep\tfield\tkey\tkind\ta\tb\tclass\tchanged_step",
             "base_vs_head\tlipid_canonical\tnorm\tstates\tn\tdigest\tb\tc\tplain_log2 (KIT-900)\tnorm")
    refreeze <- gitCommitFiles(repo, stats::setNames(
        list(put(oracle(new)), ctl),
        c(path, "docs/qa/workflow-execution/kit/kit-900/gui-entry-control.tsv")), "re-freeze")
    ref <- data.frame(lane = character(), step = character(), class = character())
    res <- valG1aPart2(repo, refreeze, "lipid_canonical", list(lipid_canonical = new), ref,
                       compare_lane, since = since)
    expect_identical(res$status, "certified")
    expect_identical(res$chain, "KIT-900")
    live_off <- valG1aPart2(repo, refreeze, "lipid_canonical", list(lipid_canonical = old), ref,
                            compare_lane, since = since)
    expect_identical(live_off$status, "not_certified")
    bare <- gitCommitFiles(repo, stats::setNames(list(put(oracle(old))), path), "no record")
    expect_identical(valG1aPart2(repo, bare, "lipid_canonical", list(lipid_canonical = old), ref,
                                 compare_lane, since = since)$status, "not_certified")
    # no re-freeze since the reference point: the live differences against the reference counts
    expect_identical(valG1aPart2(repo, since, "lipid_canonical", list(lipid_canonical = old),
                                 ref, compare_lane, since = since)$status, "certified")
    expect_identical(valG1aPart2(repo, since, "lipid_canonical", list(lipid_canonical = new),
                                 ref, compare_lane, since = since)$status, "stop")
})
