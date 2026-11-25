# KIT-167 section 10 and 5.5: gene to protein from the genome annotation, without sequence search.
# The fixtures are hand-written in the NCBI RefSeq GFF3 and GTF layout of the B36 GEO files (GSE152966); they are not
# cut from those files (no download is made from the laptop), and every identifier in them is made up.

k167Gff <- c(
    "##gff-version 3",
    "NZ_X.1\tRefSeq\tregion\t1\t20000\t.\t+\t.\tID=NZ_X.1:1..20000;Dbxref=taxon:562;gbkey=Src",
    "NZ_X.1\tRefSeq\tgene\t100\t399\t.\t+\t.\tID=gene-T_0001;Name=thrL;gbkey=Gene;gene=thrL;gene_biotype=protein_coding;locus_tag=T_0001",
    "NZ_X.1\tProtein Homology\tCDS\t100\t399\t.\t+\t0\tID=cds-WP_000001.1;Parent=gene-T_0001;gene=thrL;locus_tag=T_0001;product=leader%2C thr;protein_id=WP_000001.1",
    "NZ_X.1\tRefSeq\tpseudogene\t500\t900\t.\t+\t.\tID=gene-T_0002;Name=T_0002;gbkey=Gene;gene_biotype=pseudogene;locus_tag=T_0002;pseudo=true",
    "NZ_X.1\tProtein Homology\tCDS\t500\t900\t.\t+\t0\tID=cds-T_0002;Parent=gene-T_0002;locus_tag=T_0002;pseudo=true;product=IS element",
    "NZ_X.1\tRefSeq\tgene\t1000\t2500\t.\t+\t.\tID=gene-T_0003;Name=T_0003;gbkey=Gene;gene_biotype=rRNA;locus_tag=T_0003",
    "NZ_X.1\tcmsearch\trRNA\t1000\t2500\t.\t+\t.\tID=rna-T_0003;Parent=gene-T_0003;gbkey=rRNA;locus_tag=T_0003;product=16S ribosomal RNA",
    "NZ_X.1\tcmsearch\texon\t1000\t2500\t.\t+\t.\tID=exon-T_0003-1;Parent=rna-T_0003;gbkey=rRNA;locus_tag=T_0003",
    "NZ_X.1\tRefSeq\tgene\t3000\t4100\t.\t+\t.\tID=gene-T_0004;Name=prfB;gbkey=Gene;gene=prfB;gene_biotype=protein_coding;locus_tag=T_0004",
    "NZ_X.1\tProtein Homology\tCDS\t3000\t3074\t.\t+\t0\tID=cds-WP_000004.1;Parent=gene-T_0004;gene=prfB;locus_tag=T_0004;protein_id=WP_000004.1",
    "NZ_X.1\tProtein Homology\tCDS\t3076\t4100\t.\t+\t0\tID=cds-WP_000004.1;Parent=gene-T_0004;gene=prfB;locus_tag=T_0004;protein_id=WP_000004.1",
    "NZ_X.1\tRefSeq\tgene\t5000\t5600\t.\t-\t.\tID=gene-T_0005;gbkey=Gene;gene=yaa%2CB;gene_biotype=protein_coding;locus_tag=T_0005",
    "NZ_X.1\tProtein Homology\tCDS\t5000\t5600\t.\t-\t0\tID=cds-WP_000005.1;Parent=gene-T_0005;gene=yaa%2CB;locus_tag=T_0005;protein_id=WP_000005.1",
    "NZ_X.1\tRefSeq\tgene\t7000\t7600\t.\t+\t.\tID=gene-T_0006;gbkey=Gene;gene_biotype=protein_coding;locus_tag=T_0006",
    "NZ_X.1\tProtein Homology\tCDS\t7000\t7600\t.\t+\t0\tID=cds-WP_000005.1-2;Parent=gene-T_0006;locus_tag=T_0006;protein_id=WP_000005.1",
    "NZ_X.1\tRefSeq\tgene\t8000\t8900\t.\t+\t.\tID=gene-T_0007;gbkey=Gene;gene_biotype=protein_coding;locus_tag=T_0007",
    "NZ_X.1\tProtein Homology\tCDS\t8000\t8900\t.\t+\t0\tID=cds-WP_000007.1;Parent=gene-T_0007;locus_tag=T_0007;protein_id=WP_000007.1",
    "NZ_X.1\tRefSeq\tgene\t9000\t9900\t.\t+\t.\tID=gene-T_0008;gbkey=Gene;gene_biotype=protein_coding;locus_tag=T_0008",
    "NZ_X.1\tProtein Homology\tCDS\t9000\t9900\t.\t+\t0\tID=cds-WP_000081.1;Parent=gene-T_0008;locus_tag=T_0008;protein_id=WP_000081.1",
    "NZ_X.1\tProtein Homology\tCDS\t9030\t9900\t.\t+\t0\tID=cds-WP_000082.1;Parent=gene-T_0008;locus_tag=T_0008;protein_id=WP_000082.1")

k167Mapping <- data.frame(
    strain_accession = c("WP_000001.1", "WP_000004.1", "WP_000005.1", "WP_000081.1", "WP_000082.1"),
    reference_accession_used = c("P1", "P4", NA, "P81", NA),
    map_status = c("one_to_one", "one_to_one", "below_threshold", "one_to_one", "no_hit"),
    stringsAsFactors = FALSE)

k167WriteLines <- function(lines, ext) {
    path <- tempfile(fileext = ext)
    writeLines(lines, path)
    path
}
k167Gene <- function(genes, locus) genes[genes$locus_tag == locus, , drop = FALSE]

test_that("GFF3: every gene gets exactly one status, with the rules of 5.5", {
    path <- k167WriteLines(k167Gff, ".gff")
    result <- suppressWarnings(matchStrainGenesToProteins(path, mapping = k167Mapping))
    genes <- result$genes
    expect_identical(sort(genes$locus_tag), sprintf("T_%04d", 1:8))
    expect_identical(k167Gene(genes, "T_0001")$map_status, "one_to_one")
    expect_identical(k167Gene(genes, "T_0001")$reference_accession_used, "P1")
    expect_identical(k167Gene(genes, "T_0002")$map_status, "no_protein")
    expect_identical(k167Gene(genes, "T_0003")$map_status, "no_protein")
    expect_identical(k167Gene(genes, "T_0003")$biotype, "rRNA")
    expect_identical(k167Gene(genes, "T_0004")$map_status, "one_to_one")
    expect_identical(k167Gene(genes, "T_0005")$map_status, "below_threshold")
    expect_identical(k167Gene(genes, "T_0006")$protein_accession, "WP_000005.1")
    expect_identical(k167Gene(genes, "T_0007")$map_status, "protein_not_mapped")
    expect_identical(k167Gene(genes, "T_0008")$map_status, "ambiguous_gene")
    expect_true(is.na(k167Gene(genes, "T_0008")$reference_accession_used))
    expect_identical(k167Gene(genes, "T_0008")$protein_accession, "WP_000081.1;WP_000082.1")
    expect_identical(sum(unlist(result$counts$by_status)), nrow(genes))
})

test_that("GFF3: percent-encoding is decoded and a multi-row CDS is joined by locus tag", {
    genes <- suppressWarnings(matchStrainGenesToProteins(k167WriteLines(k167Gff, ".gff"), mapping = k167Mapping))$genes
    expect_identical(k167Gene(genes, "T_0005")$gene_name, "yaa,B")
    prfb <- k167Gene(genes, "T_0004")
    expect_identical(prfb$protein_accession, "WP_000004.1")
    expect_identical(prfb$gene_length, 75 + 1025)
    expect_identical(prfb$length_source, "cds_gff3")
})

test_that("gene length takes exon rows, else CDS rows, else the gene row", {
    genes <- suppressWarnings(matchStrainGenesToProteins(k167WriteLines(k167Gff, ".gff"), mapping = k167Mapping))$genes
    expect_identical(k167Gene(genes, "T_0003")$gene_length, 1501)
    expect_identical(k167Gene(genes, "T_0003")$length_source, "exon_gff3")
    expect_identical(k167Gene(genes, "T_0001")$gene_length, 300)
    expect_identical(k167Gene(genes, "T_0002")$length_source, "cds_gff3")
})

test_that("an accession on several genes is counted, and the not-mapped share warns in words", {
    expect_warning(result <- matchStrainGenesToProteins(k167WriteLines(k167Gff, ".gff"), mapping = k167Mapping),
        class = "reference_annotation_not_mapped_share")
    expect_identical(result$counts$accessions_on_several_genes, 1L)
    expect_identical(result$counts$protein_not_mapped, 1L)
    expect_equal(result$counts$protein_not_mapped_share, 1 / 6)
    quiet <- matchStrainGenesToProteins(k167WriteLines(k167Gff, ".gff"), mapping = k167Mapping,
        max_not_mapped_share = 0.5)
    expect_identical(quiet$counts$protein_not_mapped, 1L)
})

test_that("without a mapping the genes carry their proteins and no status from a mapping", {
    result <- matchStrainGenesToProteins(k167WriteLines(k167Gff, ".gff"))
    expect_identical(k167Gene(result$genes, "T_0001")$protein_accession, "WP_000001.1")
    expect_true(is.na(k167Gene(result$genes, "T_0001")$map_status))
    expect_identical(k167Gene(result$genes, "T_0002")$map_status, "no_protein")
})

test_that("a malformed line is refused with its line number", {
    lines <- k167Gff
    lines[5] <- "NZ_X.1\tRefSeq\tgene\t500\t900\t.\t+"
    err <- tryCatch(matchStrainGenesToProteins(k167WriteLines(lines, ".gff")),
        reference_annotation_malformed = function(e) e)
    expect_s3_class(err, "reference_annotation_malformed")
    expect_identical(err$line, 5L)
})

test_that("GTF: CDS rows keyed on locus_tag, CDS length, and a GTF without protein_id refused", {
    gtf <- c(
        "#gtf-version 2.2",
        paste0("NZ_X.1\tRefSeq\tgene\t100\t399\t.\t+\t.\tgene_id \"T_0001\"; transcript_id \"\"; gene \"thrL\";",
            " gene_biotype \"protein_coding\"; locus_tag \"T_0001\";"),
        paste0("NZ_X.1\tProtein Homology\tCDS\t100\t396\t.\t+\t0\tgene_id \"T_0001\"; transcript_id \"unassigned_1\";",
            " locus_tag \"T_0001\"; protein_id \"WP_000001.1\"; gene \"thrL\";"),
        paste0("NZ_X.1\tProtein Homology\tCDS\t3000\t4097\t.\t+\t0\tgene_id \"T_0004\"; transcript_id \"unassigned_2\";",
            " locus_tag \"T_0004\"; protein_id \"WP_000004.1\";"))
    result <- matchStrainGenesToProteins(k167WriteLines(gtf, ".gtf"), mapping = k167Mapping)
    expect_identical(result$provenance$format, "gtf")
    expect_identical(result$provenance$gene_key, "locus_tag")
    expect_identical(k167Gene(result$genes, "T_0001")$gene_length, 297)
    expect_identical(k167Gene(result$genes, "T_0001")$length_source, "cds_gtf")
    expect_identical(k167Gene(result$genes, "T_0004")$reference_accession_used, "P4")
    no_protein <- gsub(" protein_id \"[^\"]*\";", "", gtf)
    expect_error(matchStrainGenesToProteins(k167WriteLines(no_protein, ".gtf")),
        class = "reference_annotation_no_protein_id")
})

test_that("the format is detected from the attribute syntax when not given", {
    expect_identical(detectAnnotationFormat(k167WriteLines(k167Gff, ".txt")), "gff3")
})

test_that("a count table is joined on its gene identifier exactly, with the absent genes counted", {
    genes <- matchStrainGenesToProteins(k167WriteLines(k167Gff, ".gff"))$genes
    joined <- joinGenesToCountTable(genes, c("T_0001", "T_0004", "t_0005", "X_9"))
    expect_identical(joined$in_both, 2L)
    expect_identical(joined$count_only, 2L)
    expect_identical(joined$annotation_only, 6L)
})
