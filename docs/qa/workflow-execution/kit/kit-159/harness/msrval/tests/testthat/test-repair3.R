# Third pass of the pre-release repair (audit I15; D19 answered by VAL-109 (5) and the KIT-199
# plan audit's ruling C1 (a)): G1(a)'s live-difference classifier, bounded by field and key as the
# pass column of G1(a) states, and the DIA-with-limpa link read as `declared_change (KIT-199)`.

LIMPA_KEY <- "results/proteomics/da_proteins/da_proteins_groupKO-groupWT_long_annot.tsv"

g1aSteps <- function() {
    list(trx_screen = c("import", "design", "norm", "da"),
         metab_lc = c("import", "design", "norm", "da"),
         prot_dia = c("import", "design", "qc", "norm", "da", "int"),
         prot_dia_limpa = c("import", "design", "norm", "da", "int"))
}

g1aDiff <- function(lane, step, field, key, base, live) {
    data.frame(lane = lane, step = step, field = field, key = key, base = base, live = live,
               stringsAsFactors = FALSE)
}

test_that("R-I15b: live differences are classed by G1(a)'s bounds, or left outside the classes", {
    steps <- g1aSteps()
    diffs <- rbind(
        g1aDiff("prot_dia", "int", "artifacts", "integration/proteomics/x.parquet", NA, "aa"),
        g1aDiff("prot_dia", "int", "store_entries", "<store entries>", "3", "4"),
        g1aDiff("trx_screen", "design", "states", "raw", "a", "b"),
        g1aDiff("trx_screen", "import", "states", "raw", "a", "b"),
        g1aDiff("metab_lc", "norm", "outputs", "x.tsv", "a", "b"),
        g1aDiff("metab_lc", "design", "outputs", "x.tsv", "a", "b"),
        g1aDiff("prot_dia", "da", "outputs", "results/x.tsv", "a", "b"),
        g1aDiff("prot_dia_limpa", "da", "artifacts", LIMPA_KEY, "80cfedaa11", "8a4c2dd5ff"))
    cl <- valClassifyLiveDiffs(diffs, steps)
    expect_identical(cl$class, c("integration_artifacts", "integration_artifacts",
                                 "kit128_technical_replicates", NA, "metab_fixture_4", NA, NA,
                                 NA))
    # an integration artifact on both sides with another digest is outside the class
    both <- g1aDiff("prot_dia", "int", "artifacts", "integration/proteomics/x.parquet", "bb", "aa")
    expect_true(is.na(valClassifyLiveDiffs(both, steps)$class))
    # a step whose other states digest differs is outside, integration names included
    mixed <- rbind(diffs[1L, ], g1aDiff("prot_dia", "int", "states", "norm", "a", "b"))
    expect_true(all(is.na(valClassifyLiveDiffs(mixed, steps)$class)))
    # a lane-range class takes the step even where it carries integration names (C08)
    trx_int <- g1aDiff("trx_screen", "da", "artifacts", "integration/transcriptomics/y", NA, "c")
    expect_identical(valClassifyLiveDiffs(trx_int, steps)$class, "kit128_technical_replicates")
    # platform_rounding only against a committed reference fingerprint of the same table
    fp <- data.frame(lane = "prot_dia", step = "da", key = "results/x.tsv", rounding_only = TRUE)
    expect_identical(valClassifyLiveDiffs(diffs[7L, ], steps, fingerprints = fp)$class,
                     "platform_rounding")
    expect_error(valClassifyLiveDiffs(g1aDiff("nolane", "da", "states", "k", "a", "b"), steps),
                 class = "val_spec_invalid")
    # the chain reader takes the classed differences
    ref <- data.frame(lane = "metab_lc", step = "norm", class = "metab_fixture_4",
                      fields = "outputs", keys = "x.tsv")
    live <- valClassifyLiveDiffs(diffs[5L, ], steps)
    expect_identical(valReadChain("metab_lc", list(), ref, live_diffs = live)$status, "certified")
    out <- valClassifyLiveDiffs(diffs[6L, ], steps)
    expect_identical(valReadChain("metab_lc", list(), ref, live_diffs = out)$status, "stop")
})

test_that("R-K199: no admission path; the DIA-with-limpa link is a declared change; two classes in a step", {
    expect_false("limpa_da_last_digits" %in% msrval:::G1A_CLASSES)
    expect_false("admit" %in% names(formals(valClassifyLiveDiffs)))
    steps <- g1aSteps()
    # one summary step holding the DA-digest change and the integration-name change
    d <- rbind(g1aDiff("prot_dia_limpa", "int", "artifacts", LIMPA_KEY, "80cfedaa11", "8a4c2dd5ab"),
               g1aDiff("prot_dia_limpa", "int", "artifacts", "integration/proteomics/m.parquet",
                       NA, "cc"))
    expect_true(all(is.na(valClassifyLiveDiffs(d, steps)$class)))
    declared <- data.frame(lane = "prot_dia_limpa", step = "int", key = LIMPA_KEY,
                           ticket = "KIT-199")
    cl <- valClassifyLiveDiffs(d, steps, declared = declared)
    expect_identical(cl$class, c("declared_change (KIT-199)", "integration_artifacts"))
    # KIT-199's link (i) accepts its own declared change; the reference counts compare the rest
    ref <- data.frame(lane = "prot_dia_limpa", step = "int", class = "integration_artifacts",
                      fields = "artifacts", keys = "integration/proteomics/m.parquet")
    none <- data.frame(step = character(), field = character(), key = character(),
                       class = character(), head_value = character())
    link <- list(ticket = "KIT-199", changed_step = "da",
                 base_diffs = cl[, c("step", "field", "key", "class")], before_diffs = none,
                 after_diffs = none[, 1:4],
                 declared = data.frame(key = LIMPA_KEY, step = "int", expected = "8a4c2dd5ab"),
                 change_classes = "declared_change (KIT-199)")
    expect_identical(valReadChain("prot_dia_limpa", list(link), ref)$status, "certified")
    other <- link
    other$base_diffs$class[[1L]] <- "declared_change (KIT-198)"
    expect_identical(valReadChain("prot_dia_limpa", list(other), ref)$status, "not_certified")
    unlisted <- link
    unlisted$declared <- unlisted$declared[0L, ]
    expect_identical(valReadChain("prot_dia_limpa", list(unlisted), ref)$status, "not_certified")
    # with no link, a declared change is not certified by the reference alone
    expect_identical(valReadChain("prot_dia_limpa", list(), ref, live_diffs = cl)$status, "stop")
})
