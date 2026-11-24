# KIT-167 section 10: every status rule of the reciprocal-best-hit resolver (5.2 to 5.4) on hand-built hit tables.

k167Hit <- function(q, s, score, qcov = 100, scov = 100, evalue = 1e-50, pident = 99) {
    data.frame(qseqid = q, sseqid = s, pident = pident, length = 100, mismatch = 1, gapopen = 0, qstart = 1,
        qend = 100, sstart = 1, send = 100, evalue = evalue, bitscore = score / 2, score = score, qlen = 100,
        slen = 100, qcovhsp = qcov, scovhsp = scov, stringsAsFactors = FALSE)
}
k167Hits <- function(...) do.call(rbind, list(...))

k167Resolve <- function(forward, reverse, entities, ...) {
    resolveReciprocalBestHits(forward, reverse, entity_ids = entities, ...)
}
k167Status <- function(resolved, id) resolved$map_status[resolved$entity_id == id]

test_that("a unique reciprocal best hit is one_to_one and only it carries reference_accession_used", {
    forward <- k167Hits(k167Hit("s1", "r1", 100), k167Hit("s1", "r2", 50))
    reverse <- k167Hits(k167Hit("r1", "s1", 100), k167Hit("r2", "s1", 50))
    resolved <- k167Resolve(forward, reverse, "s1")
    expect_identical(k167Status(resolved, "s1"), "one_to_one")
    expect_identical(resolved$reference_accession_used, "r1")
    expect_identical(resolved$top_hit_accession, "r1")
    expect_true(resolved$top_hit_is_candidate)
})

test_that("a forward tie is tie with the tied accessions listed in C-locale order and nothing used", {
    forward <- k167Hits(k167Hit("s2", "r4", 80), k167Hit("s2", "r3", 80))
    reverse <- k167Hits(k167Hit("r3", "s2", 80), k167Hit("r4", "s2", 80))
    resolved <- k167Resolve(forward, reverse, "s2")
    expect_identical(k167Status(resolved, "s2"), "tie")
    expect_identical(resolved$reference_accession, "r3;r4")
    expect_true(is.na(resolved$reference_accession_used))
    expect_identical(resolved$status_detail, "forward_tie")
})

test_that("a reverse tie that includes the entity is tie for both strain entities", {
    forward <- k167Hits(k167Hit("s3", "r5", 90), k167Hit("s4", "r5", 90))
    reverse <- k167Hits(k167Hit("r5", "s3", 90), k167Hit("r5", "s4", 90))
    resolved <- k167Resolve(forward, reverse, c("s3", "s4"))
    expect_identical(resolved$map_status, c("tie", "tie"))
    expect_identical(resolved$status_detail, c("reverse_tie", "reverse_tie"))
    expect_true(all(is.na(resolved$reference_accession_used)))
})

test_that("an in-paralogue is many_to_one and takes no reference identifier", {
    forward <- k167Hits(k167Hit("s5", "r6", 100), k167Hit("s6", "r6", 90))
    reverse <- k167Hits(k167Hit("r6", "s5", 100), k167Hit("r6", "s6", 90))
    resolved <- k167Resolve(forward, reverse, c("s5", "s6"))
    expect_identical(k167Status(resolved, "s5"), "one_to_one")
    expect_identical(k167Status(resolved, "s6"), "many_to_one")
    expect_true(is.na(resolved$reference_accession_used[resolved$entity_id == "s6"]))
    expect_identical(resolved$reference_accession[resolved$entity_id == "s6"], "r6")
})

test_that("a best hit whose reference has no one_to_one partner is one_way", {
    forward <- k167Hits(k167Hit("s7", "r7", 70), k167Hit("s8", "r8", 80))
    reverse <- k167Hits(k167Hit("r7", "s8", 75), k167Hit("r8", "s8", 80))
    resolved <- k167Resolve(forward, reverse, c("s7", "s8"))
    expect_identical(k167Status(resolved, "s7"), "one_way")
    expect_identical(k167Status(resolved, "s8"), "one_to_one")
})

test_that("a fusion (query coverage) and a split (subject coverage) are below_threshold", {
    forward <- k167Hits(k167Hit("s9", "r9", 100, qcov = 45), k167Hit("s10", "r10", 100, scov = 50))
    reverse <- k167Hits(k167Hit("r9", "s9", 100, scov = 45), k167Hit("r10", "s10", 100, qcov = 50))
    resolved <- k167Resolve(forward, reverse, c("s9", "s10"))
    expect_identical(resolved$map_status, c("below_threshold", "below_threshold"))
    expect_identical(resolved$status_detail, c("no_candidate", "no_candidate"))
})

test_that("a top hit failing coverage while a paralogue passes is below_threshold, top_hit_filtered", {
    forward <- k167Hits(k167Hit("s11", "r11", 200, qcov = 50), k167Hit("s11", "r12", 150, qcov = 95, scov = 95))
    reverse <- k167Hits(k167Hit("r12", "s11", 150, qcov = 95, scov = 95), k167Hit("r11", "s11", 200, scov = 50))
    resolved <- k167Resolve(forward, reverse, "s11")
    expect_identical(resolved$map_status, "below_threshold")
    expect_identical(resolved$status_detail, "top_hit_filtered")
    expect_identical(resolved$top_hit_accession, "r11")
    expect_false(resolved$top_hit_is_candidate)
    expect_equal(resolved$top_hit_score, 200)
    expect_true(is.na(resolved$reference_accession_used))
})

test_that("the top-hit rule applies in the reverse direction too", {
    # r13's top hit (s13) fails coverage and outscores its best candidate (s12), so R(r13) is empty.
    forward <- k167Hits(k167Hit("s12", "r13", 120))
    reverse <- k167Hits(k167Hit("r13", "s13", 200, qcov = 40), k167Hit("r13", "s12", 120))
    resolved <- k167Resolve(forward, reverse, c("s12", "s13"))
    expect_identical(k167Status(resolved, "s12"), "one_way")
    expect_identical(k167Status(resolved, "s13"), "no_hit")
})

test_that("hits above the e-value only give no_hit", {
    forward <- k167Hits(k167Hit("s14", "r14", 100, evalue = 1e-3))
    reverse <- k167Hits(k167Hit("r14", "s14", 100, evalue = 1e-3))
    resolved <- k167Resolve(forward, reverse, "s14")
    expect_identical(resolved$map_status, "no_hit")
    expect_true(is.na(resolved$top_hit_accession))
})

test_that("ranking uses the score, never the e-value", {
    forward <- k167Hits(k167Hit("s15", "r15", 100, evalue = 1e-20), k167Hit("s15", "r16", 90, evalue = 1e-30))
    reverse <- k167Hits(k167Hit("r15", "s15", 100, evalue = 1e-20), k167Hit("r16", "s15", 90, evalue = 1e-30))
    resolved <- k167Resolve(forward, reverse, "s15")
    expect_identical(resolved$reference_accession_used, "r15")
})

test_that("the bit score can be the ranking key", {
    forward <- k167Hits(k167Hit("s1", "r1", 100), k167Hit("s1", "r2", 50))
    forward$bitscore <- c(10, 20)
    reverse <- k167Hits(k167Hit("r1", "s1", 100), k167Hit("r2", "s1", 50))
    reverse$bitscore <- c(10, 20)
    resolved <- k167Resolve(forward, reverse, "s1", ranking_key = "bitscore")
    expect_identical(resolved$reference_accession_used, "r2")
})

test_that("min_identity NULL keeps the table; a floor removes candidates in both directions", {
    forward <- k167Hits(k167Hit("s16", "r17", 100, pident = 35))
    reverse <- k167Hits(k167Hit("r17", "s16", 100, pident = 35))
    expect_identical(k167Resolve(forward, reverse, "s16")$map_status, "one_to_one")
    expect_identical(k167Resolve(forward, reverse, "s16", min_identity = NULL)$map_status, "one_to_one")
    floored <- k167Resolve(forward, reverse, "s16", min_identity = 40)
    expect_identical(floored$map_status, "below_threshold")
    # Reverse only: the forward candidate passes, the reverse one fails the floor.
    reverse$pident <- 30
    forward$pident <- 50
    expect_identical(k167Resolve(forward, reverse, "s16", min_identity = 40)$map_status, "one_way")
})

test_that("every entity gets exactly one status and the counts sum to the number of strain proteins", {
    forward <- k167Hits(k167Hit("s1", "r1", 100), k167Hit("s2", "r3", 80), k167Hit("s2", "r4", 80),
        k167Hit("s5", "r6", 100), k167Hit("s6", "r6", 90), k167Hit("s9", "r9", 100, qcov = 45))
    reverse <- k167Hits(k167Hit("r1", "s1", 100), k167Hit("r6", "s5", 100), k167Hit("r6", "s6", 90))
    entities <- c("s1", "s2", "s5", "s6", "s9", "s99")
    resolved <- k167Resolve(forward, reverse, entities)
    expect_identical(nrow(resolved), length(entities))
    expect_setequal(resolved$entity_id, entities)
    expect_true(all(resolved$map_status %in% referenceMappingStatuses))
    members <- data.frame(strain_accession = c("WP_1", "WP_2", "WP_5", "WP_6", "WP_6b", "WP_9", "WP_99"),
        entity_id = c("s1", "s2", "s5", "s6", "s6", "s9", "s99"), stringsAsFactors = FALSE)
    mapping <- assembleReferenceMapping(resolved, members)
    counts <- summariseReferenceMapping(mapping)$counts
    expect_identical(sum(unlist(counts$by_status)), nrow(members))
    expect_identical(counts$strain_proteins, nrow(members))
})

test_that("identical strain sequences are collapsed to one entity and every member inherits its status", {
    entries <- data.frame(accession = c("WP_A", "WP_B", "WP_C"),
        sequence = c("MKTAYIAKQR", "mktayiakqr*", "MSTNPKPQRK"), stringsAsFactors = FALSE)
    collapsed <- collapseIdenticalSequences(entries)
    expect_identical(nrow(collapsed$entities), 2L)
    expect_identical(collapsed$identical_sequence_groups, 1L)
    shared <- collapsed$members$entity_id[collapsed$members$strain_accession == "WP_A"]
    expect_identical(collapsed$members$entity_id[collapsed$members$strain_accession == "WP_B"], shared)
    other <- collapsed$members$entity_id[collapsed$members$strain_accession == "WP_C"]
    forward <- k167Hits(k167Hit(shared, "r1", 100), k167Hit(other, "r2", 100))
    reverse <- k167Hits(k167Hit("r1", shared, 100), k167Hit("r2", other, 100))
    resolved <- k167Resolve(forward, reverse, collapsed$entities$entity_id)
    mapping <- assembleReferenceMapping(resolved, collapsed$members)
    expect_identical(mapping$map_status, rep("one_to_one", 3L))
    expect_identical(mapping$reference_accession_used[mapping$strain_accession %in% c("WP_A", "WP_B")], c("r1", "r1"))
})

test_that("identical reference paralogues give a tie and no reference identifier is guessed", {
    forward <- k167Hits(k167Hit("s1", "rA", 100), k167Hit("s1", "rB", 100))
    reverse <- k167Hits(k167Hit("rA", "s1", 100), k167Hit("rB", "s1", 100))
    resolved <- k167Resolve(forward, reverse, "s1")
    expect_identical(resolved$map_status, "tie")
    expect_true(is.na(resolved$reference_accession_used))
})

test_that("pairs passing the e-value in one direction only are counted", {
    forward <- k167Hits(k167Hit("s1", "r1", 100), k167Hit("s1", "r2", 30, evalue = 1e-6))
    reverse <- k167Hits(k167Hit("r1", "s1", 100))
    resolved <- k167Resolve(forward, reverse, "s1")
    expect_identical(attr(resolved, "one_direction_pairs"), 1L)
})
