# KIT-267 (VAL-310): the eight Enrichment defaults the owner confirmed are no longer presented as awaiting
# confirmation, and the wording table carries the mark for the three organism defaults. Base (4d588ba3): eight of the
# 13 entries are provisional and `mark_organism` is absent. The Enrichment tab's marks are tested in
# test-kit240-integration-default.R (T10); the values are frozen there (T6 (c)).

# A wording-table entry, NA where the table has none (the base commit lacks the entry; the test then fails, not errors).
k267Text <- function(texts, key) if (key %in% names(texts)) texts[[key]] else NA_character_

test_that("KIT-267 T3: no Enrichment default is provisional, and the organism mark is worded", {
    spec <- integrationEnrichmentParameterSpec()
    expect_length(spec, 13L)
    expect_identical(sum(vapply(spec, integrationParameterProvisional, logical(1))), 0L)
    # The wording is the owner's (VAL-41); the entry must exist and say the default is human.
    expect_match(k267Text(integrationGuiTexts(), "mark_organism"), "human", fixed = TRUE)
})
