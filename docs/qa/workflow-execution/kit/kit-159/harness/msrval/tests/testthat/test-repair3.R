# Third pass of the pre-release repair (audit I15; D19 and O10 of VAL-109): G1(a)'s live-difference
# classifier, bounded by field and key as the pass column of G1(a) states, and the admission of
# `limpa_da_last_digits` only after KIT-199's merge is in the commit under test.

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

# A throwaway git repository whose HEAD has the given subject; returns its root and HEAD.
gitRepoWith <- function(subject) {
    skip_if(!nzchar(Sys.which("git")), "no git")
    repo <- withr::local_tempdir(.local_envir = parent.frame())
    git <- function(...) system2("git", c("-C", shQuote(repo), ...), stdout = TRUE, stderr = TRUE)
    git("init", "-q")
    git("-c", "user.email=val@example.invalid", "-c", "user.name=val", "commit", "-q",
        "--allow-empty", "-m", shQuote(subject))
    list(root = repo, head = git("rev-parse", "HEAD")[[1L]])
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

test_that("R-D19: limpa_da_last_digits is admitted only once KIT-199's merge is in the commit under test", {
    expect_error(msrval:::valG1aClasses(list(ruling = "KIT-199 merged", digest = "8a4c2dd5ab")),
                 class = "val_spec_invalid")
    other <- gitRepoWith("Merge KIT-198: something else")
    expect_error(msrval:::valG1aClasses(list(ruling = "KIT-199", digest = "8a4c2dd5ab",
                                            merge_commit = other$head, pkg_root = other$root)),
                 class = "val_spec_invalid")
    merged <- gitRepoWith("Merge KIT-199: re-freeze the DIA-with-limpa reference")
    admit <- list(ruling = "KIT-199", digest = "8a4c2dd5ab", merge_commit = merged$head,
                  pkg_root = merged$root)
    expect_true("limpa_da_last_digits" %in% msrval:::valG1aClasses(admit))
    bad_digest <- admit
    bad_digest$digest <- "80cfedaa"
    expect_error(msrval:::valG1aClasses(bad_digest), class = "val_spec_invalid")
    steps <- g1aSteps()
    d <- g1aDiff("prot_dia_limpa", "da", "artifacts", LIMPA_KEY, "80cfedaa11", "8a4c2dd5ab")
    expect_identical(valClassifyLiveDiffs(d, steps, admit = admit)$class, "limpa_da_last_digits")
    d2 <- d
    d2$live <- "77aa0000"
    expect_true(is.na(valClassifyLiveDiffs(d2, steps, admit = admit)$class))
    d3 <- rbind(d, g1aDiff("prot_dia_limpa", "da", "states", "da", "a", "b"))
    expect_true(all(is.na(valClassifyLiveDiffs(d3, steps, admit = admit)$class)))
    ref <- data.frame(lane = "prot_dia_limpa", step = "da", class = "limpa_da_last_digits",
                      fields = "artifacts", keys = "1")
    live <- data.frame(class = "limpa_da_last_digits", digest = "8a4c2dd5ab")
    expect_identical(valReadChain("prot_dia_limpa", list(), ref, live_diffs = live,
                                  admit = admit)$status, "certified")
    expect_identical(valReadChain("prot_dia_limpa", list(), ref, live_diffs = live)$status,
                     "stop")
})
