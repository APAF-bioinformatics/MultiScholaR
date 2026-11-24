# G1(r)'s own reference rows from the committed source blocks (KIT-159 section 2, G1(r); re-audit
# 2, R2.6 item 1): the reader of the four waldo forms into valClassifyLiveDiffs()'s input, the
# classes, the counts and the adoption of KIT-172's, KIT-180's and KIT-199's committed rows. The
# source is the committed `g1a-reference-source.txt` (the 69 GUI-entry blocks of run 44953187101's
# failures.txt); the other fixtures are the copies in `g1a/` (`g1a/SOURCES.tsv`).

g1rSource <- function() {
    harness <- Sys.getenv("VAL_HARNESS_DIR")
    expect_true(nzchar(harness))
    file.path(dirname(harness), "g1a-reference-source.txt")
}

g1rFx <- function(name) testthat::test_path("g1a", name)

g1rLaneSteps <- function() {
    t <- utils::read.delim(g1rFx("lane-steps.tsv"), colClasses = "character")
    split(t$step, t$lane)[unique(t$lane)]
}

# One block in the source's layout: the message lines indented by two spaces, the info line last
g1rBlock <- function(lines, lane = "trx_screen", step = "design") {
    c("file: test-gui-entry-transcriptomics.R", "test: KIT-101 entry lane x", "kind: failure",
      "message:", "  Expected `steps[[id]]` to be identical to `oracle$steps[[id]]`.",
      "  Differences:", paste0("  ", lines), "  ",
      sprintf("  %s at step '%s' (oracle %s)", lane, step, strrep("a", 40L)),
      "trace (first lines):", "   1. guiEntryExpectLane(\"x\")", "")
}

g1rOwnRows <- function() {
    d <- valG1rDiffs(valReadG1rBlocks(g1rSource()))
    cl <- valClassifyLiveDiffs(d, g1rLaneSteps(),
                               declared = valReadTextTable(g1rFx("k199-keys.tsv")))
    list(diffs = cl, rows = valG1rRows(cl, source = "g1a-reference-source.txt",
                                       source_sha256 = strrep("1", 64L)))
}

test_that("G1R-READ: the 69 committed blocks name 69 (lane, step) pairs, 4, 30, 22, 13 per file", {
    blocks <- valReadG1rBlocks(g1rSource())
    expect_identical(length(blocks), 69L)
    files <- vapply(blocks, `[[`, "", "file")
    per_file <- table(files)
    expect_identical(as.integer(per_file[c("test-gui-entry-lipidomics.R",
                                           "test-gui-entry-metabolomics.R",
                                           "test-gui-entry-proteomics.R",
                                           "test-gui-entry-transcriptomics.R")]),
                     c(4L, 30L, 22L, 13L))
    pairs <- paste(vapply(blocks, `[[`, "", "lane"), vapply(blocks, `[[`, "", "step"))
    expect_identical(length(unique(pairs)), 69L)
    d <- valG1rDiffs(blocks)
    expect_identical(names(d)[1:6], c("lane", "step", "field", "key", "base", "live"))
    expect_identical(nrow(d), 524L)
    per_pair <- table(factor(paste(d$lane, d$step), levels = unique(pairs)))
    expect_true(all(per_pair >= 1L))
    expect_identical(sort(unique(d$field)), c("artifacts", "notifications", "outputs", "states"))
    expect_identical(anyDuplicated(paste(d$lane, d$step, d$field, d$key)), 0L)
})

test_that("G1R-FORMS: each waldo form gives (field, key, base = expected, live = actual)", {
    b <- valReadG1rBlocks(g1rBlock(c(
        "`actual$artifacts` is length 3",
        "`expected$artifacts` is length 2",
        "",
        "names(actual$artifacts)[1:3] vs names(expected$artifacts)[1:2]",
        "  \"<store entries>\"",
        "- \"integration/x/a.parquet\"",
        "  \"results/b.tsv\"",
        "",
        "`actual$artifacts$integration/x/a.parquet` is a character vector ('<present>')",
        "`expected$artifacts$integration/x/a.parquet` is absent",
        "",
        "actual$states$s1 vs expected$states$s1",
        "- \"live1\"",
        "+ \"base1\"",
        "",
        "lines(actual$outputs$o-1)[2:4] vs lines(expected$outputs$o-1)[2:4]",
        "  \"same\"",
        "- \"after: 207\"",
        "+ \"after: 188\"",
        "  \"same too\"",
        "",
        "actual$notifications[[1]]$message vs expected$notifications[[1]]$message",
        "- \"m 206\"",
        "+ \"m 187\"",
        "",
        "`actual$notifications[[9]]` is a list",
        "`expected$notifications[[9]]` is absent",
        "",
        "`actual$artifacts$<store entries>` is a character vector ('3')",
        "`expected$artifacts$<store entries>` is a character vector ('2')")))
    expect_identical(length(b), 1L)
    expect_identical(b[[1L]]$lane, "trx_screen")
    expect_identical(b[[1L]]$step, "design")
    expect_identical(b[[1L]]$oracle, strrep("a", 40L))
    d <- valG1rDiffs(b)
    expect_identical(d$field, c("artifacts", "states", "outputs", "notifications",
                                "notifications", "store_entries"))
    expect_identical(d$key, c("integration/x/a.parquet", "s1", "o-1", "[[1]]$message", "[[9]]",
                              "<store entries>"))
    expect_identical(d$base, c(NA, "base1", "after: 188", "m 187", NA, "2"))
    expect_identical(d$live, c("<present>", "live1", "after: 207", "m 206", "<a list>", "3"))
    expect_identical(unique(d$lane), "trx_screen")
    expect_identical(unique(d$step), "design")
})

test_that("G1R-UNREADABLE: a line outside the four forms, or no lane and step, is unreadable", {
    odd <- valReadG1rBlocks(g1rBlock(c("actual$states$s1 vs expected$states$s1", "- \"a\"",
                                       "+ \"b\"", "", "a line in no form")))
    expect_error(valG1rDiffs(odd), class = "val_g1r_unreadable")
    other_key <- valReadG1rBlocks(g1rBlock(c("actual$states$s1 vs expected$states$s2", "- \"a\"",
                                             "+ \"b\"")))
    expect_error(valG1rDiffs(other_key), class = "val_g1r_unreadable")
    no_info <- g1rBlock("actual$states$s1 vs expected$states$s1")
    no_info <- no_info[!grepl(" at step '", no_info, fixed = TRUE)]
    expect_error(valReadG1rBlocks(no_info), class = "val_g1r_unreadable")
    empty <- valReadG1rBlocks(g1rBlock("`actual$artifacts` is length 3"))
    expect_error(valG1rDiffs(empty), class = "val_g1r_unreadable")
})

test_that("G1R-CLASS: every step has a class of G1(a); 73 rows, two per prot_dia_limpa summary step", {
    own <- g1rOwnRows()
    expect_false(anyNA(own$diffs$class))
    rows <- own$rows
    expect_identical(names(rows), c("lane", "step", "class", "fields", "keys",
                                    "names_list_truncated", "oracle_sha", "source",
                                    "source_sha256"))
    expect_identical(nrow(rows), 73L)
    n <- table(rows$class)
    expect_identical(as.integer(n[c("declared_change (KIT-199)", "integration_artifacts",
                                    "kit128_technical_replicates", "metab_fixture_4")]),
                     c(6L, 24L, 13L, 30L))
    expect_identical(length(n), 4L)
    two <- rows[rows$lane == "prot_dia_limpa", ]
    expect_identical(as.integer(table(two$step)[c("da", "summary_refusal", "summary_save",
                                                  "summary_copy", "summary_report",
                                                  "summary_export")]),
                     c(1L, 1L, 2L, 2L, 2L, 2L))
    for (lane in unique(rows$lane)) {
        expect_identical(valStepClassCounts(own$diffs[own$diffs$lane == lane, ]),
                         valClassCounts(rows$class[rows$lane == lane]), info = lane)
    }
    # the declared rows carry KIT-199's stored and live digests (link (i)'s value check)
    dec <- own$diffs[own$diffs$class %in% "declared_change (KIT-199)", ]
    expect_true(valDeclaredValuesOk(dec, valReadTextTable(g1rFx("k199-keys.tsv"))))
    expect_true(all(rows$names_list_truncated == "no"))
    expect_true(all(rows$oracle_sha %in% c("79d6bc726fe1361264ef12a15ca8cbc490e07d52",
                                           "5d38127250d0122ebbe09d52b0477295c870112c")))
})

test_that("G1R-ADOPT: KIT-180's and KIT-199's rows agree; KIT-172's disagree in 18 groups", {
    rows <- g1rOwnRows()$rows
    k180 <- valReadReferenceRows(g1rFx("k180-ref.tsv"))
    expect_identical(nrow(valReferenceDisagreements(rows, k180)), 0L)
    expect_identical(valAdoptReferenceRows(rows, k180), rows)
    k199 <- valReadReferenceRows(g1rFx("k199-ref.tsv"))
    expect_identical(nrow(valReferenceDisagreements(rows, k199)), 0L)
    expect_identical(valAdoptReferenceRows(rows, k199), rows)
    # KIT-172 kept the line window of `lines(actual$outputs$<key>)[a:b]` in one key; without the
    # committed key errata G1(r) stops
    k172 <- valReadReferenceRows(g1rFx("k172-ref.tsv"))
    expect_error(valAdoptReferenceRows(rows, k172, errata = NULL), class = "val_spec_invalid")
    dis <- valReferenceDisagreements(rows, k172, errata = NULL)
    expect_identical(nrow(dis), 18L)
    expect_identical(sort(unique(dis$lane)), c("metab_combined", "metab_gc"))
    expect_true(all(dis$in_own & dis$in_earlier))
    expect_true(all(dis$earlier_only %in%
                        c("outputs/metabolomics_workflow-norm-correlation_filter_summary)[2:8]",
                          "outputs/metabolomics_workflow-norm-correlation_filter_summary)[2:11]")))
    expect_true(all(grepl("outputs/metabolomics_workflow-norm-correlation_filter_summary",
                          dis$own_only, fixed = TRUE)))
})

test_that("G1R-ERRATA: the 18 listed key pairs are accepted, exactly, and nothing else", {
    rows <- g1rOwnRows()$rows
    k172 <- valReadReferenceRows(g1rFx("k172-ref.tsv"))
    errata <- valG1rKeyErrata()
    expect_identical(nrow(errata), 18L)
    expect_identical(sort(unique(errata$lane)), c("metab_combined", "metab_gc"))
    expect_identical(errata$earlier_key,
                     paste0(errata$key, ifelse(errata$lane == "metab_gc", ")[2:8]", ")[2:11]")))
    # with the committed errata (the default) KIT-172's rows are adopted
    expect_identical(nrow(valReferenceDisagreements(rows, k172)), 0L)
    expect_identical(valAdoptReferenceRows(rows, k172), rows)
    # a 19th, unlisted disagreement still stops: G1(r)'s own key changed in an unlisted group
    other <- rows
    at <- which(other$lane == "metab_lc" & other$step == "correlation")
    other$keys[[at]] <- paste0(other$keys[[at]], ";states/unlisted")
    expect_error(valAdoptReferenceRows(other, k172), class = "val_spec_invalid")
    dis <- valReferenceDisagreements(other, k172)
    expect_identical(nrow(dis), 1L)
    expect_identical(c(dis$lane, dis$step), c("metab_lc", "correlation"))
    # a windowed key not on the list (another window) still stops
    k172b <- k172
    at <- which(k172b$lane == "metab_gc" & k172b$step == "correlation")
    k172b$keyset[[at]] <- sub(")[2:8]", ")[2:9]", k172b$keyset[[at]], fixed = TRUE)
    expect_error(valAdoptReferenceRows(rows, k172b), class = "val_spec_invalid")
    # a listed pair applies only to its own (lane, step, class)
    k172c <- k172
    at <- which(k172c$lane == "metab_gc" & k172c$step == "correlation")
    k172c$class[[at]] <- "kit128_technical_replicates"
    expect_error(valAdoptReferenceRows(rows, k172c), class = "val_spec_invalid")
})
