# KIT-116 W2, audit item A103-03: the frozen values of two parity tests of test-omics-design-core.R ("C1 the save
# preflight is today's commit preflight" and "C3 to C5 a builder session through the core gives today's result"),
# which since KIT-103 W4/W5 compare the shared design code with itself (the lipid and metab functions became wrappers
# over it). The values come from the old lipid and metab functions at 616e1147, the last commit where they existed as
# separate code: its child 5d077268 (KIT-103 W4) deleted the seven lipid design helper files and made
# validateLipidDesignDaPreflight and mod_lipid_design_builder_server wrappers over the shared code (metab's followed
# at c56d1c64, W5). The two blocks passed at 616e1147. Run from a checkout of 616e1147:
#
#   Rscript a103-03-design.R <checkout of 616e1147> <output dir>
#
# The script loads the package and test helpers of that checkout and the test file's own helper definitions
# (omicsDesignFixture, omicsDesignBuilderSession, omicsDesignCoreSession, ...), runs each block's old side on the
# block's inputs, and checks each old value against the core of the same commit as the block did. Values are sha256
# of the value the block compared with expect_identical. Output: <output dir>/a103-03-frozen.R and
# a103-03-old-values.rds.

args <- commandArgs(trailingOnly = TRUE)
green <- normalizePath(args[[1]])
out_dir <- args[[2]]
stopifnot(identical(system2("git", c("-C", green, "rev-parse", "--short=8", "HEAD"), stdout = TRUE), "616e1147"))
pkgload::load_all(green, quiet = TRUE)
withr::local_dir(green)
# The helpers defer their temporary owners to testthat's teardown environment, which test_file initialises.
testthat:::local_teardown_env()
testthat::source_test_helpers(file.path(green, "tests", "testthat"), env = globalenv())
library(testthat)

exprs <- parse(file.path(green, "tests", "testthat", "test-omics-design-core.R"), keep.source = FALSE)
for (e in exprs) if (!(is.call(e) && identical(e[[1L]], as.name("test_that")))) eval(e, globalenv())

digestOf <- function(x) digest::digest(x, algo = "sha256")
frozen <- list()
old_values <- list()
unfaithful <- character()
record <- function(group, key, value, core_value) {
    if (is.null(frozen[[group]])) frozen[[group]] <<- list()
    frozen[[group]][key] <<- list(value)
    if (!identical(value, core_value)) unfaithful <<- c(unfaithful, paste(group, key))
}
frozenCode <- function(frozen, prefix) {
    vapply(names(frozen), function(group) {
        values <- frozen[[group]]
        keys <- names(values)
        keys <- ifelse(make.names(keys) == keys, keys, paste0("`", keys, "`"))
        lines <- vapply(seq_along(values), function(i) {
            paste0("    ", keys[[i]], " = ", paste(deparse(values[[i]], width.cutoff = 500L), collapse = " "))
        }, character(1))
        code <- paste0(prefix, group, " <- list(\n", paste(lines, collapse = ",\n"), ")")
        stopifnot(identical(eval(parse(text = code)), values))
        code
    }, character(1))
}

# C1: the save preflight, validateLipidDesignDaPreflight and validateMetabDesignDaPreflight on three results each.
(function() {
    for (omic in c("lipidomics", "metabolomics")) {
        fixture <- omicsDesignFixture(omic, layout = "combined")
        old <- if (identical(omic, "lipidomics")) validateLipidDesignDaPreflight else validateMetabDesignDaPreflight
        drifted <- fixture$payload
        drifted$data_cln[[1]]$WT_1 <- NULL
        no_formula <- fixture$payload
        no_formula$config_list$deAnalysisParameters$formula_string <- ""
        cases <- list(payload = fixture$payload, drifted = drifted, no_formula = no_formula)
        for (name in names(cases)) {
            results <- cases[[name]]
            today <- old(designMatrix = results$design_matrix, assayList = results$data_cln,
                contrastsTbl = results$contrasts_tbl, formulaString = "~ 0 + group",
                columnMapping = results$column_mapping, requireContrasts = FALSE)
            old_values$save_preflight[[paste(omic, name)]] <<- today
            record("save_preflight", paste(omic, name), digestOf(today),
                digestOf(validateOmicsDesignSave(results, results$column_mapping, omic)))
        }
    }
})()

# C3 to C5: one builder session on mod_lipid_design_builder_server (lipid, and trx through its profile) and
# mod_metab_design_builder_server; the saved result.
(function() {
    samples <- c("WT_2", "KO_1", "WT_1", "KO_3", "WT_3", "KO_2")
    lipid_assays <- module_ci_lipid_design_assays("lc_pair", samples = samples)
    lipid_mapping <- module_ci_lipid_design_column_mapping(samples)
    today <- omicsDesignBuilderSession(mod_lipid_design_builder_server, lipid_assays, lipid_mapping)
    core <- omicsDesignCoreSession(lipid_assays, lipid_mapping, omicsDesignTestConfig("lipidomics"))
    old_values$builder$lipidomics <<- today$result
    record("builder", "lipidomics", digestOf(today$result), digestOf(core$payload))

    metab_assays <- module_ci_metab_design_assays("lc_pair", samples = samples)
    metab_mapping <- module_ci_metab_design_column_mapping(samples)
    today <- suppressMessages(omicsDesignBuilderSession(mod_metab_design_builder_server, metab_assays, metab_mapping))
    core <- omicsDesignCoreSession(metab_assays, metab_mapping, omicsDesignTestConfig("metabolomics"))
    old_values$builder$metabolomics <<- today$result
    record("builder", "metabolomics", digestOf(today$result), digestOf(core$payload))

    trx <- omicsDesignTestConfig("transcriptomics")
    today <- omicsDesignBuilderSession(mod_lipid_design_builder_server, lipid_assays, lipid_mapping,
        list(profile = trxOmicProfile()))
    core <- omicsDesignCoreSession(lipid_assays, lipid_mapping, trx)
    core$payload$group_levels <- trx$design$group_levels(core$payload$design_matrix, core$factors)
    core$payload$excluded_runs <- core$removed
    old_values$builder$transcriptomics <<- today$result
    record("builder", "transcriptomics", digestOf(today$result), digestOf(core$payload))
})()

if (length(unfaithful)) stop("the old value is not the green core's for: ", paste(unfaithful, collapse = ", "))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(old_values, file.path(out_dir, "a103-03-old-values.rds"))
writeLines(frozenCode(frozen, "omicsDesignFrozen_"), file.path(out_dir, "a103-03-frozen.R"))
cat("frozen groups:", paste(names(frozen), lengths(frozen), sep = "=", collapse = " "), "\n")
