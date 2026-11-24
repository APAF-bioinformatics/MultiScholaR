# KIT-167 section 10 and 8.4: strain FASTA header classes, the OX rule, contaminants and the input refusals.

k167Seq <- function(i, length = 60L) {
    set.seed(20261123 + i)
    paste(sample(strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]], length, replace = TRUE), collapse = "")
}
k167WriteFasta <- function(headers, sequences) {
    path <- tempfile(fileext = ".fasta")
    writeLines(paste0(">", headers, "\n", sequences), path)
    path
}
k167RefseqHeaders <- function(n, offset = 0L) sprintf("ref|WP_%09d.1| MULTISPECIES: protein %d [Klebsiella]",
    offset + seq_len(n), seq_len(n))
k167HumanHeaders <- function(n, ox = 9606L) sprintf("sp|P%05d|PROT%d_HUMAN Protein %d OS=Homo sapiens OX=%d GN=G%d PE=1 SV=1",
    seq_len(n), seq_len(n), seq_len(n), ox, seq_len(n))

test_that("each header class gives its accession", {
    entries <- data.frame(header = c(
        "ref|WP_000001.1| MULTISPECIES: thing [Klebsiella]",
        "lcl|NZ_CP000001.1_cds_WP_000002.1_1 [locus_tag=ABC_0001] [protein=x] [protein_id=WP_000002.1] [gbkey=CDS]",
        "WP_000003.1 hypothetical protein",
        "sp|P0A7V8|RS4_ECOLI 30S ribosomal protein S4 OS=Escherichia coli (strain K12) OX=83333 GN=rpsD PE=1 SV=2",
        "tr|A0A0H3|A0A0H3_ECOLI Uncharacterized OS=Escherichia coli OX=83333 PE=4 SV=1",
        "P00761 SWISS-PROT:P00761|TRYP_PIG Trypsin - Sus scrofa (Pig).",
        "CON__P02768"),
        sequence = vapply(1:7, k167Seq, character(1)), stringsAsFactors = FALSE)
    classified <- classifyStrainHeaders(entries)$entries
    expect_identical(classified$class, c("refseq", "ncbi_cds", "plain", "uniprot", "uniprot", "contaminant",
        "contaminant"))
    expect_identical(classified$accession[1:5], c("WP_000001.1", "WP_000002.1", "WP_000003.1", "P0A7V8", "A0A0H3"))
    expect_identical(classified$locus_tag[2], "ABC_0001")
    expect_identical(classified$ox[4:5], c("83333", "83333"))
    expect_false(any(classified$selected[6:7]))
})

test_that("ref|WP| entries plus human entries with no taxon named select no human entry, counted per OX", {
    path <- k167WriteFasta(c(k167RefseqHeaders(60), k167HumanHeaders(20)), vapply(1:80, k167Seq, character(1)))
    entries <- readStrainProteinFasta(path)
    classified <- classifyStrainHeaders(entries)
    expect_identical(sum(classified$entries$selected), 60L)
    expect_false(any(classified$entries$selected[classified$entries$class == "uniprot"]))
    expect_identical(classified$counts$excluded_by_ox, list(`9606` = 20L))
})

test_that("UniProt entries with one OX and no other strain class are strain entries", {
    path <- k167WriteFasta(k167HumanHeaders(55, ox = 83333L), vapply(1:55, k167Seq, character(1)))
    classified <- classifyStrainHeaders(readStrainProteinFasta(path))
    expect_identical(sum(classified$entries$selected), 55L)
})

test_that("UniProt entries only with two OX and no taxon are refused; a named taxon selects its entries", {
    headers <- c(k167HumanHeaders(30, ox = 83333L), sub("P0", "Q0", k167HumanHeaders(30, ox = 9606L)))
    path <- k167WriteFasta(headers, vapply(1:60, k167Seq, character(1)))
    entries <- readStrainProteinFasta(path)
    expect_error(classifyStrainHeaders(entries), class = "reference_mapping_taxon_ambiguous")
    chosen <- classifyStrainHeaders(entries, strain_taxon = 83333)
    expect_identical(sum(chosen$entries$selected), 30L)
    expect_identical(chosen$counts$excluded_by_ox, list(`9606` = 30L))
})

test_that("a duplicated accession with different sequences is refused and listed", {
    headers <- c(k167RefseqHeaders(60), k167RefseqHeaders(1))
    path <- k167WriteFasta(headers, vapply(1:61, k167Seq, character(1)))
    expect_error(prepareStrainEntities(path), class = "reference_mapping_duplicate_id")
    err <- tryCatch(prepareStrainEntities(path), reference_mapping_duplicate_id = function(e) e)
    expect_identical(err$accessions, "WP_000000001.1")
})

test_that("a duplicated accession with the same sequence (one WP_ on two locus tags) is one entry", {
    seqs <- vapply(1:60, k167Seq, character(1))
    headers <- sprintf("lcl|NZ_X_cds_WP_%09d.1_%d [locus_tag=L_%04d] [protein_id=WP_%09d.1]", c(1:60, 1L),
        c(1:60, 61L), c(1:60, 61L), c(1:60, 1L))
    path <- k167WriteFasta(headers, c(seqs, seqs[1]))
    prepared <- prepareStrainEntities(path)
    expect_identical(nrow(prepared$entities), 60L)
    expect_identical(length(unique(prepared$members$strain_accession)), 60L)
})

test_that("the alphabet rules: terminal stop removed, internal stop and nucleotides refused", {
    seqs <- vapply(1:60, k167Seq, character(1))
    stopped <- k167WriteFasta(k167RefseqHeaders(60), paste0(seqs, "*"))
    prepared <- prepareStrainEntities(stopped)
    expect_false(any(grepl("\\*", prepared$entities$sequence)))
    internal <- seqs
    internal[3] <- "MKT*AYIAK"
    expect_error(prepareStrainEntities(k167WriteFasta(k167RefseqHeaders(60), internal)),
        class = "reference_mapping_not_protein")
    dna <- vapply(1:60, function(i) paste(rep("ACGTACGTTG", 6), collapse = ""), character(1))
    expect_error(prepareStrainEntities(k167WriteFasta(k167RefseqHeaders(60), dna)),
        class = "reference_mapping_not_protein")
    letters_ok <- seqs
    letters_ok[1] <- "MKXUOBZJAA"
    expect_silent(prepareStrainEntities(k167WriteFasta(k167RefseqHeaders(60), letters_ok)))
})

test_that("fewer than 50 strain entities are refused", {
    path <- k167WriteFasta(k167RefseqHeaders(49), vapply(1:49, k167Seq, character(1)))
    expect_error(prepareStrainEntities(path), class = "reference_mapping_too_few")
})

test_that("contaminant entries are excluded by class and counted", {
    headers <- c(k167RefseqHeaders(60), "P00761 SWISS-PROT:P00761|TRYP_PIG Trypsin", "CON__P02768",
        "Streptavidin (S.avidinii)")
    path <- k167WriteFasta(headers, vapply(1:63, k167Seq, character(1)))
    prepared <- prepareStrainEntities(path)
    expect_identical(nrow(prepared$members), 60L)
    expect_identical(prepared$counts$by_class$contaminant, 3L)
})
