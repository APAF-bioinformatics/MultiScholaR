# The generated-input corpus (KIT-159 G1(b), G8, G12, G13, L1a, S0; r1.14).
#
# Every route's input is written by a generator the package checkout already has, read only (never
# edited): the five proteomics formats by the e2e fixture writers of
# tests/testdata/e2e/generate_fixtures.R, whose abundances come from tools/profiling/
# kit_protein_corpus.R (kitProteinRecipe, kitProteinMatrix; FragPipe from the LFQ recipe), and
# transcriptomics, lipidomics and metabolomics by tools/profiling/kit_omics_corpus.R. The writers
# are taken from their files without running the files' top-level code (valExtractHelpers;
# sys.source behind the tools' own `sys.nframe() == 0L` guard). Sizes are declared here for the
# app's default QC. 7.2 places the corpus writer under `lifecycle/`; it lives in msrval instead,
# so that it is tested, linted and pinned with the harness (conventions/R.md: a package, not
# sourced scripts); outcome.md records the divergence.

#' The second organism of the two-species input (A18, G12, G13): its FASTA name, NCBI taxon and
#' the share of the protein universe (the last proteins of the universe) given to it
#' @noRd
CORPUS_SECOND_SPECIES <- list(name = "Escherichia coli (strain K12)", taxon = 83333L, share = 0.3)

#' The declared corpora: id, the routes each serves, the generator and its arguments (JSON)
#'
#' The proteomics corpora share one protein universe (800 synthetic accessions, seed 20261124) and
#' the KIT-021 LFQ or TMT recipe at six (TMT: five, one TMT10 plex) biological samples per group;
#' the two-species corpus is the `prot_dia` corpus with the last 30 % of the universe's FASTA
#' entries moved to a second organism. Transcriptomics, lipidomics and metabolomics: three
#' biological samples per group, each measured twice.
#'
#' @return A data frame: `corpus`, `routes` (comma-separated), `generator`, `arguments`.
#' @export
valCorpusSpec <- function() {
    e2e <- "tests/testdata/e2e/generate_fixtures.R"
    omics <- "tools/profiling/kit_omics_corpus.R"
    prot <- function(writer, seed_file, seed, kit_recipe = "lfq", per_group = 6L, extra = list()) {
        args <- c(list(writer = writer, seed_file = seed_file, kit_recipe = kit_recipe,
                       seed = seed, universe_seed = 20261124L, biological_per_group = per_group,
                       proteins = 800L), extra)
        as.character(jsonlite::toJSON(args, auto_unbox = TRUE))
    }
    omic <- function(kind, features) {
        args <- list(kind = kind, biological_per_group = 3L, features = features, standards = 6L)
        as.character(jsonlite::toJSON(args, auto_unbox = TRUE))
    }
    data.frame(
        corpus = c("prot_dia", PARITY_TWO_SPECIES, "prot_lfq", "prot_lfq_fragpipe", "prot_tmt",
                   "trx_screen", "lipid_canonical", "metab"),
        routes = c("prot_dia,prot_dia_limpa", "prot_dia", "prot_lfq", "prot_lfq_fragpipe",
                   "prot_tmt", "trx_screen", "lipid_canonical", "metab_lc,metab_gc,metab_combined"),
        generator = c(rep(e2e, 5L), rep(omics, 3L)),
        arguments = c(
            prot("generateProtDiannFixture", "report.tsv", 20261125L),
            prot("generateProtDiannFixture", "report.tsv", 20261125L,
                 extra = list(second_species = CORPUS_SECOND_SPECIES)),
            prot("generateProtMaxQuantFixture", "proteinGroups.txt", 20261126L,
                 extra = list(flagged_rows = 10L)),
            prot("generateProtFragPipeFixture", "combined_protein.tsv", 20261127L),
            prot("generateProtPdTmtFixture", "pd_tmt_peptides.tsv", 20261128L,
                 kit_recipe = "tmt", per_group = 5L),
            omic("trx_screen", 2000L),
            omic("lipid_screen", 200L),
            omic("metab_screen", 300L)
        ),
        stringsAsFactors = FALSE
    )
}

#' The corpus a route runs on (the two-species corpus is named, never chosen by route)
#'
#' @param route Route id.
#' @return The corpus id.
#' @export
valCorpusFor <- function(route) {
    spec <- valCorpusSpec()
    spec <- spec[spec$corpus != PARITY_TWO_SPECIES, , drop = FALSE]
    hit <- spec$corpus[vapply(strsplit(spec$routes, ",", fixed = TRUE),
                              function(r) route %in% r, logical(1))]
    if (length(hit) != 1L) valAbort("no generated corpus serves route {route}", "val_spec_invalid")
    hit
}

#' Write one corpus into an empty folder, with its manifest (G1(b): generator, arguments, hashes)
#'
#' @param corpus Corpus id ([valCorpusSpec()]).
#' @param dir An empty (or absent) folder.
#' @param pkg_root The package checkout whose generators are used, read only.
#' @return A list: `corpus`, `dir`, `generator`, `arguments`, `manifest` (corpus, file, bytes,
#'   sha256).
#' @export
valWriteCorpus <- function(corpus, dir, pkg_root) {
    spec <- valCorpusSpec()
    row <- spec[spec$corpus == corpus, , drop = FALSE]
    if (nrow(row) != 1L) valAbort("no generated corpus {corpus}", "val_spec_invalid")
    if (length(list.files(dir, all.files = TRUE, no.. = TRUE))) {
        valAbort("the corpus folder {dir} is not empty", "val_spec_invalid")
    }
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    args <- jsonlite::fromJSON(row$arguments)
    if (identical(row$generator, "tests/testdata/e2e/generate_fixtures.R")) {
        valWriteProtCorpus(args, dir, pkg_root)
    } else {
        valWriteOmicsCorpus(args, dir, pkg_root)
    }
    files <- sort(list.files(dir, recursive = TRUE))
    paths <- file.path(dir, files)
    manifest <- data.frame(corpus = rep(corpus, length(files)), file = files,
                           bytes = as.numeric(file.size(paths)),
                           sha256 = vapply(paths, valSha256, character(1), USE.NAMES = FALSE),
                           stringsAsFactors = FALSE)
    list(corpus = corpus, dir = normalizePath(dir), generator = row$generator,
         arguments = row$arguments, manifest = manifest)
}

#' The four proteomics writers of the e2e fixture generator, extracted from the checkout
#' @noRd
valE2eProteinWriters <- function(pkg_root) {
    parent <- new.env(parent = globalenv())
    # kitProteinTools() reads tools/profiling/kit_protein_corpus.R relative to FIXTURE_ROOT.
    parent$FIXTURE_ROOT <- normalizePath(file.path(pkg_root, "tests", "testdata", "e2e"))
    valExtractHelpers(file.path(parent$FIXTURE_ROOT, "generate_fixtures.R"),
                      c("generateProtDiannFixture", "generateProtMaxQuantFixture",
                        "generateProtFragPipeFixture", "generateProtPdTmtFixture"),
                      parent = parent)
}

#' A proteomics corpus: the e2e writer at the declared size, then (two-species) the second
#' organism in the FASTA, then the UniProt cache from the FASTA
#' @noRd
valWriteProtCorpus <- function(args, dir, pkg_root) {
    writers <- valE2eProteinWriters(pkg_root)
    shape <- args[intersect(c("kit_recipe", "seed", "universe_seed", "biological_per_group",
                              "proteins", "flagged_rows"), names(args))]
    lane <- list(lane_id = "generated", groups = list("WT", "KO"), seed_file = args$seed_file,
                 fasta_file = "proteins.fasta", fixture_recipe = shape)
    # The writers print progress lines; the corpus is their files.
    utils::capture.output(get(args$writer, envir = writers)(lane, dir))
    fasta <- file.path(dir, "proteins.fasta")
    if (!is.null(args$second_species)) valSecondSpeciesFasta(fasta, args$second_species)
    valCorpusUniprotCache(fasta, file.path(dir, "uniprot_annotations.RDS"))
    invisible(dir)
}

#' A transcriptomics, lipidomics or metabolomics corpus from tools/profiling/kit_omics_corpus.R
#' @noRd
valWriteOmicsCorpus <- function(args, dir, pkg_root) {
    tools <- new.env(parent = globalenv())
    sys.source(file.path(pkg_root, "tools", "profiling", "kit_omics_corpus.R"), envir = tools)
    recipe <- tools$kitOmicsRecipe(args$kind, as.integer(args$biological_per_group),
                                   as.integer(args$features),
                                   standards = as.integer(args$standards))
    writer <- switch(args$kind, trx_screen = tools$writeKitTrxCorpus,
                     lipid_screen = tools$writeKitLipidCorpus,
                     metab_screen = tools$writeKitMetabCorpus)
    writer(dir, recipe)
    invisible(dir)
}

#' Move the last share of a FASTA's entries to a second organism (OS= and OX= of the header)
#' @noRd
valSecondSpeciesFasta <- function(path, second) {
    lines <- readLines(path)
    at <- which(startsWith(lines, ">"))
    moved <- at[seq_along(at) > length(at) - round(length(at) * second$share)]
    first <- "OS=Homo sapiens OX=9606"
    if (!all(grepl(first, lines[moved], fixed = TRUE))) {
        valAbort("a FASTA header to move does not name the first organism", "val_spec_invalid")
    }
    lines[moved] <- sub(first, sprintf("OS=%s OX=%d", second$name, as.integer(second$taxon)),
                        lines[moved], fixed = TRUE)
    writeLines(lines, path, useBytes = TRUE)
    invisible(path)
}

#' The offline UniProt cache of a synthetic FASTA, in the columns the design step reads, as
#' tests/testdata/gui-entry/uniprot/build.R builds it for the e2e universe, with each entry's
#' organism taken from its header (no network: KIT-159 ruling 4)
#' @noRd
valCorpusUniprotCache <- function(fasta, path) {
    headers <- grep("^>", readLines(fasta), value = TRUE)
    entry <- sub("^>sp\\|([^|]+)\\|.*$", "\\1", headers)
    cache <- data.frame(
        Entry = entry, From = entry, gene_names = sub("^.* GN=([^ ]+) .*$", "\\1", headers),
        Protein_existence = "Evidence at protein level",
        Protein_names = sub("^>sp\\|[^|]+\\|[^ ]+ (.*) OS=.*$", "\\1", headers),
        annotation_score = 5, Organism = sub("^.* OS=(.*) OX=.*$", "\\1", headers),
        go_id_go_biological_process = "GO:0006915; GO:0007049",
        go_term_go_biological_process = "apoptotic process; cell cycle",
        go_id_go_molecular_function = "GO:0005515",
        go_term_go_molecular_function = "protein binding",
        go_id_go_cellular_compartment = "GO:0005739",
        go_term_go_cellular_compartment = "mitochondrion", stringsAsFactors = FALSE
    )
    saveRDS(cache, path, version = 3)
    invisible(path)
}

#' A route's cell on a written corpus: the accepted (or alternative) script, the files by role,
#' the design, the contrast, the assays and, for proteomics, the UniProt cache and annotation rule;
#' on the two-species corpus the mixed-species box as a declared forced input (7.4, A18)
#'
#' @param route Route id.
#' @param corpus From [valWriteCorpus()].
#' @param alternative `NULL` (the accepted script) or a G1(b) alternative ([valActionScript()]).
#' @return A list: `route`, `omic`, `corpus`, `alternative`, `script`, `files`, `design`,
#'   `contrasts`, `assays`, `forced_inputs`, `uniprot_cache`, `annotation`.
#' @export
valCorpusCell <- function(route, corpus, alternative = NULL) {
    spec <- valCorpusSpec()
    row <- spec[spec$corpus == corpus$corpus, , drop = FALSE]
    if (nrow(row) != 1L || !route %in% strsplit(row$routes, ",", fixed = TRUE)[[1L]]) {
        valAbort("corpus {corpus$corpus} does not serve route {route}", "val_spec_invalid")
    }
    args <- jsonlite::fromJSON(row$arguments)
    r <- valRoute(route)
    at <- function(f) file.path(corpus$dir, f)
    design <- utils::read.delim(at("design.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
    metab <- list(metab_lc = c(LCMS_Pos = "msdial_lcms_pos.txt"),
                  metab_gc = c(GCMS = "msdial_gcms.txt"),
                  metab_combined = c(LCMS_Pos = "msdial_lcms_pos.txt", GCMS = "msdial_gcms.txt"))
    by_assay <- switch(r$omic,
                       lipidomics = c(LCMS_Pos = "lipidsearch_lcms_pos.txt",
                                      LCMS_Neg = "lipidsearch_lcms_neg.txt"),
                       metabolomics = metab[[route]])
    files <- switch(r$omic,
                    proteomics = list(search_results = at(args$seed_file),
                                      fasta = at("proteins.fasta")),
                    transcriptomics = list(counts = at("counts.tsv")),
                    stats::setNames(as.list(at(by_assay)),
                                    paste0("assay", seq_along(by_assay))))
    tool <- if (identical(r$omic, "metabolomics")) "msdial" else r$import_tool
    trx <- identical(r$omic, "transcriptomics")
    contrasts <- list(if (trx) c("Sera", "RPMI") else c("KO", "WT"))
    forced <- list()
    if (!is.null(args$second_species)) {
        forced <- structure(list(mixed_species_fasta = TRUE),
                            reason = paste("two-species generated input (A18): the app's help",
                                           "text asks for the box when the FASTA holds proteins",
                                           "of several species"))
    }
    annotation <- NULL
    cache <- NULL
    if (identical(r$omic, "proteomics")) {
        cache <- at("uniprot_annotations.RDS")
        headers <- grep("^>", readLines(files$fasta), value = TRUE)
        accessions <- sub("^>sp\\|([^|]+)\\|.*$", "\\1", headers)
        human <- accessions[grepl(" OX=9606 ", headers, fixed = TRUE)]
        annotation <- list(expected = "full", min_share = 1, accessions = accessions,
                           in_taxon = function(a) a %in% human)
    }
    script <- valActionScript(route, assays = names(by_assay), import_tool = tool,
                              contrasts = contrasts, alternative = alternative)
    list(route = route, omic = r$omic, corpus = corpus$corpus, alternative = alternative,
         script = script, files = files, design = design, contrasts = contrasts,
         assays = names(by_assay), forced_inputs = forced, uniprot_cache = cache,
         annotation = annotation)
}
