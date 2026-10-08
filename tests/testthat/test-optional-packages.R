# CORE-018: the eleven app-only packages are Suggests, loaded on demand behind classed guards.

optional_packages <- c("Glimma", "cowplot", "ComplexHeatmap", "ggraph", "Biobase", "seqinr",
    "git2r", "devtools", "circlize", "GGally", "ggpubr")

test_that("loading the namespace in a fresh process loads none of the optional packages", {
    script <- tempfile(fileext = ".R")
    on.exit(unlink(script))
    writeLines(c(sprintf(".libPaths(%s)", paste(deparse(.libPaths()), collapse = "")),
        "invisible(suppressPackageStartupMessages(loadNamespace('MultiScholaR')))",
        "cat(loadedNamespaces(), sep = '\\n')"), script)
    loaded <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", script), stdout = TRUE, stderr = FALSE)
    expect_true("MultiScholaR" %in% loaded)
    expect_identical(intersect(loaded, optional_packages), character(0))
})

test_that("requireOptionalPackage refuses with a classed error when the package is absent", {
    testthat::local_mocked_bindings(optionalPackageAvailable = function(package) FALSE, .package = "MultiScholaR")
    expect_error(MultiScholaR:::requireOptionalPackage("ggraph", "a test"), class = "multischolar_package_missing")
    callers <- list(
        ComplexHeatmap = function() MultiScholaR:::renderProtDAHeatmapProduct(),
        Glimma = function() MultiScholaR::getGlimmaVolcanoProteomicsWidget(),
        ggraph = function() MultiScholaR:::cnetplotEdited(),
        GGally = function() MultiScholaR:::compareUmapComponentsPairs(),
        ggpubr = function() MultiScholaR:::rlePcaPlotList(),
        seqinr = function() MultiScholaR:::parseFastaFile(),
        git2r = function() MultiScholaR:::pushProjectToGithub(),
        Biobase = function() MultiScholaR:::getOneRlePlotData())
    for (name in names(callers)) expect_error(callers[[name]](), class = "multischolar_package_missing", label = name)
})

test_that("requireOptionalPackage passes when the package is installed", {
    expect_invisible(MultiScholaR:::requireOptionalPackage("cli", "a test"))
})
