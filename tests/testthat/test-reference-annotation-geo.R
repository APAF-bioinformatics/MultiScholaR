# KIT-167 ruling (e): the annotation parser on a real cut of the public GEO annotation of strain B36 (GSE152966,
# `GSE152966_GCF_900622635.1_B36_genomic.gff.gz` and `.gtf.gz`; NCBI places no restrictions on GEO data). The cut and
# the expected protein accessions per locus were made on a runner by `docs/qa/workflow-execution/kit/kit-167/
# pass2-inputs.R` (section CUT): the loci of the parser's cases (the first 40 coding loci, multi-row CDS, pseudogenes,
# RNA genes, protein accessions on several loci), the expectation derived with awk from the cut, not by the parser.
# README.md beside the files gives the run, the source files' SHA-256 and the cut's. The GEO GTF files carry no
# `protein_id` (gate F4), so there is no GTF cut: the parser refuses such a file by design.

k167GeoFile <- function(name) testthat::test_path("..", "testdata", "reference-mapping", name)

test_that("the GFF3 cut gives every locus the protein accessions awk finds in its CDS rows", {
    expected <- utils::read.delim(k167GeoFile("CUT-b36-expected.tsv"), colClasses = "character", na.strings = "NA")
    genes <- matchStrainGenesToProteins(k167GeoFile("CUT-b36.gff"), format = "gff3")$genes
    expect_setequal(genes$locus_tag, expected$locus_tag)
    got <- genes$protein_accession[match(expected$locus_tag, genes$locus_tag)]
    expect_identical(got, expected$protein_accession)
    no_protein <- genes$locus_tag[genes$map_status %in% "no_protein"]
    expect_setequal(no_protein, expected$locus_tag[is.na(expected$protein_accession)])
})
