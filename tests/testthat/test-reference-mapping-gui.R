# KIT-167 sections 3 (items 3, 4), 5.6 and 8.3: the import panels, the tab's 2 x 2 counts and download, the taxon
# pre-fill, all behind `MultiScholaR.reference_mapping`. The engine is the stub of helper-reference-mapping.R.

k167Html <- function(tags) as.character(htmltools::renderTags(tags)$html)

test_that("with the gate off no import or enrichment page shows a reference mapping control", {
    withr::local_options(MultiScholaR.reference_mapping = NULL)
    for (omic in c("proteomics", "transcriptomics")) {
        expect_false(grepl("refmap_", k167Html(mod_omic_import_ui("x", omicGuiConfig(omic))), fixed = TRUE), info = omic)
    }
    page <- k167Html(mod_omic_enrich_ui("x", omicGuiConfig("proteomics")$enrichment))
    expect_false(grepl("x-id_source", page, fixed = TRUE))
    expect_false(grepl("reference_mapping_counts", page, fixed = TRUE))
})

test_that("with the gate on the proteomics import offers the mapping and transcriptomics the annotation step", {
    withr::local_options(MultiScholaR.reference_mapping = TRUE)
    prot <- k167Html(mod_omic_import_ui("x", omicGuiConfig("proteomics")))
    for (id in c("x-refmap_strain_fasta", "x-refmap_reference", "x-refmap_coverage", "x-refmap_evalue",
        "x-refmap_min_identity", "x-refmap_strain_taxon", "x-refmap_run", "x-refmap_summary")) {
        expect_true(grepl(id, prot, fixed = TRUE), info = id)
    }
    expect_true(grepl(referenceMappingTexts$applicability, prot, fixed = TRUE))
    trx <- k167Html(mod_omic_import_ui("x", omicGuiConfig("transcriptomics")))
    for (id in c("x-refmap_annotation", "x-refmap_strain_fasta", "x-refmap_mapping_table", "x-refmap_run",
        "x-refmap_summary")) {
        expect_true(grepl(id, trx, fixed = TRUE), info = id)
    }
    enrich <- k167Html(mod_omic_enrich_ui("x", omicGuiConfig("proteomics")$enrichment))
    expect_true(grepl("x-id_source", enrich, fixed = TRUE))
    expect_true(grepl("x-reference_mapping_counts", enrich, fixed = TRUE))
    expect_true(grepl("x-download_unmapped_changed", enrich, fixed = TRUE))
})

test_that("the import's mapping step stores the mapping and reports counts, identity and applicability", {
    manager <- k167Owner()
    strain <- k167StrainFromReference()
    k167UseStubEngine(k167StubDiamond())
    done <- runReferenceMappingImport(manager, strain$path, "UP000000750",
        settings = list(coverage = 80, evalue = 1e-5, min_identity = NULL, strain_taxon = NULL))
    expect_true(done$ok)
    expect_false(is.null(readReferenceMapping(manager)))
    text <- referenceMappingSummaryText(done$result)
    expect_true(any(grepl("one_to_one: 60", text, fixed = TRUE)))
    expect_true(any(grepl("Median identity", text, fixed = TRUE)))
    expect_true(any(grepl(referenceMappingTexts$applicability, text, fixed = TRUE)))
    expect_identical(referenceMappingTaxon(manager), "301447")
})

test_that("a refusal is returned as its plain sentence, nothing is stored", {
    manager <- k167Owner()
    withr::local_envvar(MSR_DIAMOND = file.path(tempdir(), "no-diamond-here"))
    withr::local_options(MultiScholaR.diamond_path = NULL)
    done <- runReferenceMappingImport(manager, k167StrainFromReference()$path, "UP000000750", settings = list())
    expect_false(done$ok)
    expect_identical(done$class, "reference_mapping_engine_missing")
    expect_match(done$message, "Saved mappings and the demo tables still work", fixed = TRUE)
    expect_null(readReferenceMapping(manager))
})

test_that("the identity caution is shown below a median identity of 95 %", {
    result <- list(counts = list(strain_proteins = 2L, by_status = as.list(stats::setNames(c(0L, 0L, 0L, 2L, 0L, 0L),
            referenceMappingStatuses))), identity_summary = list(pairs = 2L, median = 90, p05 = 88, share_below_40 = 0,
            caution = TRUE), provenance = list(reference_name = "R", reference_release = "2026_04"))
    expect_true(any(grepl("below 95%", referenceMappingSummaryText(result), fixed = TRUE)))
})

test_that("the transcriptomics step stores the gene table from the project's mapping or a saved mapping table", {
    manager <- k167Owner()
    gff <- tempfile(fileext = ".gff")
    writeLines(c("##gff-version 3",
        "NZ_X.1\tRefSeq\tgene\t100\t399\t.\t+\t.\tID=gene-T_0001;gbkey=Gene;gene_biotype=protein_coding;locus_tag=T_0001",
        "NZ_X.1\tProtein Homology\tCDS\t100\t399\t.\t+\t0\tID=cds-WP_000000001.1;locus_tag=T_0001;protein_id=WP_000000001.1"),
        gff)
    mapping <- data.frame(strain_accession = "WP_000000001.1", reference_accession_used = "P1", map_status = "one_to_one",
        stringsAsFactors = FALSE)
    table <- tempfile(fileext = ".tsv")
    utils::write.table(mapping, table, sep = "\t", quote = FALSE, row.names = FALSE)
    done <- runReferenceGenesImport(manager, gff, mapping_table = table)
    expect_true(done$ok)
    genes <- readReferenceGenes(manager)
    expect_identical(genes$reference_accession_used, "P1")
    expect_true(any(grepl("Genes: 1", referenceGenesSummaryText(done$result), fixed = TRUE)))
    bad <- tempfile(fileext = ".tsv")
    writeLines("accession\tother", bad)
    refused <- runReferenceGenesImport(manager, gff, mapping_table = bad)
    expect_false(refused$ok)
    expect_identical(refused$class, "reference_mapping_reference_invalid")
})

test_that("the tab's 2 x 2 table gives counts and both shares; the download lists the unmapped changed groups", {
    counts <- list(groups = 5L, mapped = 2L, ambiguous_group = 1L, unmapped = 2L, changed_mapped = 1L,
        changed_unmapped = 2L, unchanged_mapped = 1L, unchanged_unmapped = 1L,
        unmapped_changed_groups = c("WP_4;WP_1", "WP_5"))
    table <- referenceMappingContrastTable(counts)
    expect_identical(table$features, c("mapped", "unmapped"))
    expect_identical(table$changed, c(1L, 2L))
    expect_identical(table$unchanged, c(1L, 1L))
    expect_equal(table$share_of_changed, c(1 / 3, 2 / 3))
    expect_equal(table$share_of_unchanged, c(1 / 2, 1 / 2))
    expect_identical(referenceMappingUnmappedTable(counts)$protein_group, c("WP_4;WP_1", "WP_5"))
})
