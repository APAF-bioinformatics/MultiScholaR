# G1(a)'s re-freeze chain reader and the reference rows (KIT-159 G1(a), G1(r); T-CHAIN-1).

refRows <- function() {
    data.frame(lane = "lipid_canonical", step = c("design", "norm", "summary_export"),
               class = c("integration_artifacts", "integration_artifacts", "integration_artifacts"),
               fields = "artifacts", keys = "integration/lipidomics/manifest",
               stringsAsFactors = FALSE)
}

diffRows <- function(steps, class, key = "integration/lipidomics/manifest", field = "artifacts") {
    n <- length(steps)
    data.frame(step = steps, field = rep_len(field, n), key = rep_len(key, n),
               class = rep_len(class, n), stringsAsFactors = FALSE)
}

goodLink <- function(base_diffs) {
    list(ticket = "KIT-172", changed_step = "norm",
         base_diffs = base_diffs,
         before_diffs = diffRows(character(), character()),
         after_diffs = diffRows(c("norm", "da"), "plain_log2 (KIT-172)", key = "log2"),
         declared = data.frame(key = character(), step = character(), expected = character()),
         change_classes = c("plain_log2 (KIT-172)", "declared_change (KIT-172)"))
}

test_that("T-CHAIN-1: a first link absorbing an unintended difference is refused at link one", {
    ref <- refRows()
    bad <- goodLink(diffRows(c("design", "norm"), "integration_artifacts"))
    res <- valReadChain("lipid_canonical", list(bad, goodLink(diffRows(character(), character()))),
                        reference = ref)
    expect_identical(res$status, "not_certified")
    expect_identical(res$failed_link, 1L)
})

test_that("T-CHAIN-1: a correct chain is accepted; a later link's base is compared with zero", {
    ref <- refRows()
    link1 <- goodLink(diffRows(c("design", "norm", "summary_export"), "integration_artifacts"))
    link2 <- goodLink(diffRows(character(), character()))
    expect_identical(valReadChain("lipid_canonical", list(link1, link2), reference = ref)$status,
                     "certified")
    link2_bad <- goodLink(diffRows("norm", "integration_artifacts"))
    res <- valReadChain("lipid_canonical", list(link1, link2_bad), reference = ref)
    expect_identical(res$status, "not_certified")
    expect_identical(res$failed_link, 2L)
})

test_that("T-CHAIN-1: a re-frozen file without a record; a lane not re-frozen against the counts", {
    res <- valReadChain("lipid_canonical", list(NULL), reference = refRows())
    expect_identical(res$status, "not_certified")
    expect_match(res$reason, "control record")
    live <- diffRows(c("design", "norm", "summary_export"), "integration_artifacts")
    expect_identical(valReadChain("lipid_canonical", list(), reference = refRows(),
                                  live_diffs = live)$status, "certified")
    expect_identical(valReadChain("lipid_canonical", list(), reference = refRows(),
                                  live_diffs = live[1:2, ])$status, "stop")
})

test_that("T-CHAIN-1: declared keys before the changed step, and unclassed differences after it", {
    ref <- refRows()
    link <- goodLink(diffRows(c("design", "norm", "summary_export"), "integration_artifacts"))
    link$before_diffs <- data.frame(step = "import_refusal", field = "outputs", key = "offset",
                                    class = "declared_change (KIT-172)", head_value = "0.5")
    link$declared <- data.frame(key = "offset", step = "import_refusal", expected = "0.5")
    expect_identical(valReadChain("lipid_canonical", list(link), reference = ref)$status,
                     "certified")
    wrong <- link
    wrong$before_diffs$head_value <- "0.7"
    expect_identical(valReadChain("lipid_canonical", list(wrong), reference = ref)$status,
                     "not_certified")
    unlisted <- link
    unlisted$declared <- unlisted$declared[0, ]
    expect_identical(valReadChain("lipid_canonical", list(unlisted), reference = ref)$status,
                     "not_certified")
    unclassed <- link
    unclassed$after_diffs$class[2] <- NA
    expect_identical(valReadChain("lipid_canonical", list(unclassed), reference = ref)$status,
                     "not_certified")
})

test_that("T-CHAIN-1: rows a re-freezing ticket committed earlier are adopted when they agree", {
    own <- refRows()
    expect_identical(valAdoptReferenceRows(own, own[own$lane == "lipid_canonical", ]), own)
    other <- own
    other$step[2] <- "da"
    expect_error(valAdoptReferenceRows(own, other), class = "val_spec_invalid")
})
