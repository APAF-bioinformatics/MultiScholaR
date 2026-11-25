# The generated-input corpus (KIT-159 G1(b), G8, G12, G13, L1a, S0; r1.14): valCorpusSpec,
# valWriteCorpus, valCorpusCell and the two-species input. Every corpus is written by a generator
# the checkout already has, read only (tests/testdata/e2e/generate_fixtures.R for the five
# proteomics formats, its abundances from tools/profiling/kit_protein_corpus.R;
# tools/profiling/kit_omics_corpus.R for transcriptomics, lipidomics and metabolomics).

corpusIds <- function() {
    c("prot_dia", PARITY_TWO_SPECIES, "prot_lfq", "prot_lfq_fragpipe", "prot_tmt", "trx_screen",
      "lipid_canonical", "metab")
}

generatorFiles <- function(root) {
    file.path(root, c("tests/testdata/e2e/generate_fixtures.R",
                      "tools/profiling/kit_protein_corpus.R",
                      "tools/profiling/kit_omics_corpus.R"))
}

test_that("T-CORP-1: each corpus declares its generator, arguments and the routes it serves", {
    spec <- valCorpusSpec()
    expect_setequal(spec$corpus, corpusIds())
    expect_true(all(nzchar(spec$generator)))
    for (a in spec$arguments) expect_type(jsonlite::fromJSON(a), "list")
    expect_setequal(unlist(strsplit(spec$routes, ",", fixed = TRUE)), driverRoutes())
    for (route in driverRoutes()) {
        expect_true(valCorpusFor(route) %in% setdiff(spec$corpus, PARITY_TWO_SPECIES), info = route)
    }
    expect_identical(spec$routes[spec$corpus == PARITY_TWO_SPECIES], "prot_dia")
    expect_error(valWriteCorpus("prot_unknown", withr::local_tempdir(), driverPkgRoot()),
                 class = "val_spec_invalid")
    expect_error(valCorpusFor("prot_unknown"), class = "val_spec_invalid")
})

test_that("T-CORP-2: each corpus is written with its manifest, twice byte-identically; generators untouched", {
    root <- driverPkgRoot()
    gens <- generatorFiles(root)
    before <- vapply(gens, valSha256, character(1))
    watched <- file.path(root, c("tests/testdata/e2e", "tools/profiling"))
    listed <- list.files(watched, recursive = TRUE, all.files = TRUE)
    for (id in corpusIds()) {
        a <- valWriteCorpus(id, withr::local_tempdir(), root)
        b <- valWriteCorpus(id, withr::local_tempdir(), root)
        expect_identical(a$corpus, id)
        expect_gt(nrow(a$manifest), 1L)
        expect_identical(a$manifest$file, b$manifest$file, info = id)
        expect_identical(a$manifest$sha256, b$manifest$sha256, info = id)
        for (i in seq_len(nrow(a$manifest))) {
            p <- file.path(a$dir, a$manifest$file[[i]])
            expect_identical(valSha256(p), a$manifest$sha256[[i]], info = p)
            expect_equal(as.numeric(file.size(p)), as.numeric(a$manifest$bytes[[i]]), info = p)
        }
    }
    # A folder that is not empty (here the generators' own) is refused before anything is written.
    expect_error(valWriteCorpus("prot_dia", dirname(gens[[1L]]), root), class = "val_spec_invalid")
    expect_identical(vapply(gens, valSha256, character(1)), before)
    expect_identical(list.files(watched, recursive = TRUE, all.files = TRUE), listed)
})

test_that("T-CORP-3: a route's cell on its corpus: files, format, design and contrast agree with the data", {
    root <- driverPkgRoot()
    written <- list()
    for (route in driverRoutes()) {
        id <- valCorpusFor(route)
        if (is.null(written[[id]])) written[[id]] <- valWriteCorpus(id, withr::local_tempdir(), root)
        cell <- valCorpusCell(route, written[[id]])
        expect_identical(cell$route, route)
        expect_identical(cell$script$route, route)
        expect_true(all(file.exists(unlist(cell$files))), info = route)
        expect_true(all(c("sample", "group") %in% names(cell$design)), info = route)
        expect_length(cell$contrasts, 1L)
        expect_true(all(cell$contrasts[[1L]] %in% cell$design$group), info = route)
        # default QC: nothing forced on a generated input
        expect_length(cell$forced_inputs, 0L)
        data_file <- cell$files[[1L]]
        present <- if (route %in% c("prot_dia", "prot_dia_limpa")) {
            unique(utils::read.delim(data_file, stringsAsFactors = FALSE)$Run)
        } else {
            names(utils::read.delim(data_file, nrows = 1L, check.names = FALSE))
        }
        if (!is.null(cell$files$assay2)) {
            present <- c(present, names(utils::read.delim(cell$files$assay2, nrows = 1L,
                                                          check.names = FALSE)))
        }
        found <- vapply(cell$design$sample, function(s) any(grepl(s, present, fixed = TRUE)),
                        logical(1))
        expect_true(all(found), info = route)
        if (cell$omic %in% c("lipidomics", "metabolomics")) {
            expect_identical(length(cell$assays), length(Filter(Negate(is.null),
                                                                 cell$files[c("assay1", "assay2")])))
        }
    }
    expect_error(valCorpusCell("trx_screen", written[["prot_dia"]]), class = "val_spec_invalid")
})

test_that("T-CORP-4: a proteomics corpus carries a UniProt cache covering its FASTA; the rule passes", {
    root <- driverPkgRoot()
    corpus <- valWriteCorpus("prot_lfq", withr::local_tempdir(), root)
    cell <- valCorpusCell("prot_lfq", corpus)
    expect_true(file.exists(cell$uniprot_cache))
    cache <- readRDS(cell$uniprot_cache)
    fasta <- grep("^>", readLines(cell$files$fasta), value = TRUE)
    accessions <- sub("^>sp\\|([^|]+)\\|.*$", "\\1", fasta)
    expect_setequal(cache$Entry, accessions)
    expect_true(all(cache$Organism == "Homo sapiens"))
    expect_true(is.na(valAnnotationRule(cache, cell$annotation, character())))
    expect_error(valAnnotationRule(cache[0, ], cell$annotation, character()),
                 class = "val_annotation_empty")
})

test_that("T-CORP-5: the two-species input: the prot_dia report, a FASTA of two organisms, the box forced", {
    root <- driverPkgRoot()
    two <- valWriteCorpus(PARITY_TWO_SPECIES, withr::local_tempdir(), root)
    one <- valWriteCorpus("prot_dia", withr::local_tempdir(), root)
    sha <- function(x, f) x$manifest$sha256[x$manifest$file == f]
    expect_identical(sha(two, "report.tsv"), sha(one, "report.tsv"))
    expect_false(identical(sha(two, "proteins.fasta"), sha(one, "proteins.fasta")))
    cell <- valCorpusCell("prot_dia", two)
    expect_identical(cell$forced_inputs$mixed_species_fasta, TRUE)
    expect_true(nzchar(attr(cell$forced_inputs, "reason")))
    second <- jsonlite::fromJSON(valCorpusSpec()$arguments[valCorpusSpec()$corpus ==
                                                                PARITY_TWO_SPECIES])$second_species
    headers <- grep("^>", readLines(cell$files$fasta), value = TRUE)
    taxa <- sub("^.* OX=([0-9]+) .*$", "\\1", headers)
    expect_setequal(unique(taxa), c("9606", as.character(second$taxon)))
    expect_identical(sum(taxa != "9606"), as.integer(round(length(headers) * second$share)))
    # The app's own organism analysis, as its import runs it on a ticked box
    ns <- asNamespace("MultiScholaR")
    mapping <- ns$extractOrganismsFromFasta(cell$files$fasta)
    report <- utils::read.delim(cell$files$search_results, stringsAsFactors = FALSE)
    dist <- ns$analyzeOrganismDistribution(unique(report$Protein.Group), mapping)
    choices <- ns$prepareProtImportOrganismSelectionChoices(dist)
    expect_setequal(as.character(choices$validOrganisms$taxon_id),
                    c("9606", as.character(second$taxon)))
    expect_identical(choices$selectedTaxon, "9606")
    cache <- readRDS(cell$uniprot_cache)
    expect_setequal(unique(cache$Organism), c("Homo sapiens", second$name))
    expect_true(is.na(valAnnotationRule(cache, cell$annotation, character())))
})

test_that("T-CORP-6: on the two-species input the driver confirms the organism modal at its rendered values", {
    root <- driverPkgRoot()
    corpus <- valWriteCorpus(PARITY_TWO_SPECIES, withr::local_tempdir(), root)
    cell <- valCorpusCell("prot_dia", corpus)
    script <- cell$script
    script$steps <- script$steps[1L]
    expect_identical(script$steps[[1L]]$id, "import")
    res <- valRunWorkflow(script, files = cell$files, design = cell$design,
                          helpers = valGuiHelpers(root), forced_inputs = cell$forced_inputs,
                          uniprot_cache = cell$uniprot_cache, annotation = cell$annotation,
                          pkg_root = root, work_dir = withr::local_tempdir())
    expect_identical(res$exit, "completed",
                     info = paste(res$steps$step, res$steps$outcome, res$steps$reason,
                                  res$steps$message, collapse = "; "))
    led <- res$ledger
    last <- function(id) {
        r <- led[led$input_id == id, , drop = FALSE]
        r[nrow(r), , drop = FALSE]
    }
    mix <- last("setup_import-mixed_species_fasta")
    expect_identical(c(mix$value, mix$source), c("TRUE", "forced"))
    org <- last("setup_import-selected_organism")
    expect_identical(c(org$value, org$source, org$location), c("9606", "ui_default", "modal"))
    filt <- last("setup_import-filter_to_organism")
    expect_identical(c(filt$value, filt$source, filt$location), c("FALSE", "ui_default", "modal"))
    pressed <- res$actions$input_id[res$actions$kind == "press"]
    expect_true("setup_import-confirm_organism" %in% pressed)
})
