test_that("checkpoint recipes are explicit and ordinary tests remain small", {
    small <- readWorkflowCheckpointRecipe("small")
    expect_identical(small$id, "small")
    expect_error(readWorkflowCheckpointRecipe("missing"), "resolve uniquely")
    expect_identical(validateWorkflowCheckpointRecipe(small), small)
    for (id in c("representative", "tall", "wide", "annotation_heavy", "cycles", "capacity")) {
        expect_error(validateWorkflowCheckpointRecipe(readWorkflowCheckpointRecipe(id)),
            "explicit opt-in")
    }
    broken <- small
    broken$n_samples <- 9
    expect_error(validateWorkflowCheckpointRecipe(broken), "divisible by 4")
    broken <- small
    broken$n_proteins <- Inf
    expect_error(validateWorkflowCheckpointRecipe(broken), "invalid fixture dimension")
    expect_error(validateWorkflowFixtureRoot(here::here(), here::here()), "outside package")
    if (identical(Sys.info()[["sysname"]], "Linux")) {
        filesystem <- system2("stat", c("-f", "-c", "%T", shQuote(tempdir())), stdout = TRUE)
        if (filesystem %in% c("tmpfs", "ramfs")) {
            expect_error(validateWorkflowFixtureRoot(tempdir(), here::here()), "disk-backed")
        }
    }
})

test_that("precursor chunks are deterministic across chunk boundaries and preserve RNG", {
    recipe <- readWorkflowCheckpointRecipe("small")
    set.seed(77)
    before <- .Random.seed
    whole <- makeWorkflowPrecursorChunk(recipe, 1:4)
    split <- rbind(makeWorkflowPrecursorChunk(recipe, 1:2),
        makeWorkflowPrecursorChunk(recipe, 3:4))
    expect_identical(whole, split)
    expect_identical(.Random.seed, before)
    expect_equal(nrow(whole), 4 * 3 * 2 * 12)
    expect_false(anyDuplicated(whole[c("Protein.Group", "Precursor.Id", "Run")]) > 0L)
    expect_true(all(grepl("^[ACDEFGHIKLMNPQRSTVWY]+$", whole$Stripped.Sequence)))
    expect_true(all(whole$Q.Value <= 0.01))
    expect_true(all(whole$Proteotypic == 1))
    expect_true(anyNA(whole$Precursor.Quantity))
    expect_true(all(whole$Precursor.Normalised ==
        whole$Precursor.Quantity * 1.03, na.rm = TRUE))
    expect_error(makeWorkflowPrecursorChunk(readWorkflowCheckpointRecipe("wide"), 1:1000),
        "allocation boundary")
})

test_that("real peptide and protein fixtures preserve design, scale and metadata", {
    recipe <- readWorkflowCheckpointRecipe("small")
    peptide <- makeWorkflowCheckpointPeptide(recipe)
    protein <- makeWorkflowCheckpointProtein(recipe)
    expect_s4_class(peptide, "PeptideQuantitativeData")
    expect_s4_class(protein, "ProteinQuantitativeData")
    expect_true(methods::validObject(peptide))
    expect_true(methods::validObject(protein))
    expect_identical(peptide@protein_id_column, "Protein.Group")
    expect_false(peptide@is_logged_data)
    expect_identical(setdiff(names(protein@protein_quant_table), "Protein.Group"),
        protein@design_matrix$Run)
    expect_identical(rownames(protein@design_matrix), protein@design_matrix$Run)
    expect_identical(levels(protein@design_matrix$batch),
        c("batch_1", "batch_2", "unused_batch"))
    expect_true(is.ordered(protein@design_matrix$acquisition))
    expect_false(protein@args$fixture$estimator_executed)
    expect_identical(protein, makeWorkflowCheckpointProtein(recipe))
    expect_identical(unserialize(serialize(protein, NULL)), protein)
    expect_identical(unserialize(serialize(peptide, NULL)), peptide)
    raw <- serialize(protein@protein_quant_table, NULL)
    expect_gt(length(memCompress(raw, "gzip")) / length(raw), 0.4)
})

test_that("changed generations do not alias unchanged quantities", {
    recipe <- readWorkflowCheckpointRecipe("small")
    first <- makeWorkflowCheckpointProtein(recipe, generation = 0)
    second <- makeWorkflowCheckpointProtein(recipe, generation = 1)
    shared <- first
    shared@args$fixture$metadata_only <- TRUE
    expect_identical(first@protein_quant_table, shared@protein_quant_table)
    expect_identical(rlang::obj_address(first@protein_quant_table[[2L]]),
        rlang::obj_address(shared@protein_quant_table[[2L]]))
    expect_false(identical(rlang::obj_address(first@protein_quant_table[[2L]]),
        rlang::obj_address(second@protein_quant_table[[2L]])))
    expect_false(identical(first@protein_quant_table, second@protein_quant_table))
    expect_identical(first@protein_id_table, second@protein_id_table)
    expect_identical(first@design_matrix, second@design_matrix)
    expect_false(isTRUE(first@args$fixture$metadata_only))
    # Explicit retained owner: no memory-saving assertion is made while it exists.
    retained_alias <- list(first = first, shared = shared)
    manager <- MultiScholaR::WorkflowState$new()
    manager$saveState("first", first, list(), "synthetic first generation")
    manager$saveState("second", second, list(), "synthetic changed generation")
    expect_identical(manager$revertToState("first"), retained_alias$first)
    expect_true("second" %in% names(manager$states))
    expect_identical(retained_alias$shared@protein_quant_table, first@protein_quant_table)
})

test_that("explicit fixture source root survives an external service working directory", {
    source_root <- normalizePath(testthat::test_path("..", ".."))
    helpers <- new.env(parent = baseenv())
    sys.source(file.path(source_root, "tests", "testthat",
        "helper-workflow-checkpoint-fixtures.R"), envir = helpers)
    observed <- NULL
    helpers$validateWorkflowFixtureRoot <- function(output_root, repo_root) {
        observed <<- list(output_root = output_root, repo_root = repo_root)
        invisible(output_root)
    }
    recipe <- helpers$readWorkflowCheckpointRecipe("small", source_root)
    recipe$annotation_chars <- 257L
    output_root <- withr::local_tempdir()
    object <- withr::with_dir(output_root,
        helpers$makeWorkflowCheckpointProtein(recipe, allow_large = TRUE,
            output_root = output_root, repo_root = source_root))
    expect_identical(observed, list(output_root = output_root, repo_root = source_root))
    expect_true(methods::validObject(object))
})

test_that("codec edge fixtures preserve nonfinite values and selector order exactly", {
    object <- makeWorkflowCheckpointEdges(
        makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small")))
    restored <- unserialize(serialize(object, NULL))
    expect_identical(restored, object)
    values <- restored@protein_quant_table[[2]]
    expect_true(is.na(values[1]) && !is.nan(values[1]))
    expect_true(is.nan(values[2]))
    expect_identical(values[3:6], c(Inf, -Inf, 0, -0))
    expect_true(all(is.na(restored@protein_quant_table[7, -1])))
    expect_true(all(restored@protein_quant_table[8, -1] == 17))
    selected <- restored@protein_quant_table[c(3, 1, 3), c(1, 4, 2), drop = FALSE]
    expect_identical(selected$Protein.Group, object@protein_quant_table$Protein.Group[c(3, 1, 3)])
    expect_identical(names(selected), names(object@protein_quant_table)[c(1, 4, 2)])
    expect_identical(unserialize(serialize(selected, NULL)), selected)
    patched <- restored
    patched@protein_quant_table[9L, 3L] <- 123.456
    expect_identical(restored@protein_quant_table, object@protein_quant_table)
    expect_identical(patched@protein_quant_table[9L, 3L], 123.456)
    expect_identical(patched@design_matrix, restored@design_matrix)
    expect_identical(patched@protein_id_table, restored@protein_id_table)
    expect_false(identical(serialize(patched, NULL), serialize(restored, NULL)))
})

test_that("uncertainty controls differ even when abundance is identical", {
    protein <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    base <- makeWorkflowUncertaintyControl(protein)
    se <- makeWorkflowUncertaintyControl(protein, "standard_error")
    counts <- makeWorkflowUncertaintyControl(protein, "observations")
    reordered <- makeWorkflowUncertaintyControl(protein, "order")
    expect_s4_class(base, "EList")
    expect_identical(base$E, se$E)
    expect_identical(base$E, counts$E)
    expect_false(identical(serialize(base, NULL), serialize(se, NULL)))
    expect_false(identical(serialize(base, NULL), serialize(counts, NULL)))
    expect_identical(dim(base$E), dim(reordered$E))
    expect_identical(reordered$genes$protein.id, rownames(reordered$E))
    expect_identical(dimnames(reordered$E), dimnames(reordered$other$standard.error))
    expect_identical(dimnames(reordered$E), dimnames(reordered$other$n.observations))
    expect_identical(base$E, reordered$E[rev(seq_len(nrow(base$E))), , drop = FALSE])
    expect_true(attr(base, "fixture_only_not_dpcQuant"))
})

test_that("fixtures exercise real precursor aggregation and peptide imputation", {
    recipe <- readWorkflowCheckpointRecipe("small")
    original <- makeWorkflowCheckpointPeptide(recipe)
    # These are actual package methods, not estimator doubles.
    rolled <- MultiScholaR::rollUpPrecursorToPeptide(original, core_utilisation = 1)
    imputed <- MultiScholaR::peptideMissingValueImputation(rolled,
        proportion_missing_values = 0.5, core_utilisation = 1)
    # DIA-003 authorized repair; the original DIA-002 failures remain in evidence.
    expect_true(methods::validObject(rolled))
    expect_true(methods::validObject(imputed))
    expect_identical(rolled@raw_quantity_column, "Peptide.RawQuantity")
    expect_identical(rolled@norm_quantity_column, "Peptide.Normalised")
    expect_equal(nrow(rolled@peptide_data), recipe$n_proteins *
        recipe$peptides_per_protein * recipe$n_samples)
    expect_true(all(c("Peptide.Imputed", "is_imputed") %in% names(imputed@peptide_data)))
    observed <- !is.na(rolled@peptide_data$Peptide.RawQuantity)
    expect_identical(imputed@peptide_data$Peptide.Imputed[observed],
        rolled@peptide_data$Peptide.RawQuantity[observed])
    expect_false(original@is_logged_data)
})

test_that("small fixture is consumed by real IQ without pretending synthetic proteins are fits", {
    recipe <- readWorkflowCheckpointRecipe("small")
    input <- makeWorkflowCheckpointPeptide(recipe)
    rolled <- MultiScholaR::rollUpPrecursorToPeptide(input, core_utilisation = 1)
    imputed <- MultiScholaR::peptideMissingValueImputation(rolled,
        proportion_missing_values = 0.5, core_utilisation = 1)
    iq_input <- imputed@peptide_data
    # These compatibility substitutions are made by the original GUI adapter.
    iq_input$Q.Value <- 0.0009
    iq_input$PG.Q.Value <- 0.009
    iq_input$Peptide.Imputed[is.na(iq_input$Peptide.Imputed)] <- 0
    directory <- withr::local_tempdir()
    input_path <- file.path(directory, "peptide_values_imputed.tsv")
    output_path <- file.path(directory, "iq_output_file.txt")
    vroom::vroom_write(iq_input, input_path, progress = FALSE)
    iq::process_long_format(input_path, output_filename = output_path,
        sample_id = "Run", primary_id = "Protein.Group",
        secondary_id = "Stripped.Sequence", intensity_col = "Peptide.Imputed",
        filter_double_less = c("Q.Value" = "0.01", "PG.Q.Value" = "0.01"),
        normalization = "none", pdf_out = file.path(directory, "qc-plots.pdf"))
    result <- as.data.frame(vroom::vroom(output_path, show_col_types = FALSE,
        progress = FALSE, .name_repair = "minimal"))
    expect_identical(result$Protein.Group, unique(iq_input$Protein.Group))
    expect_identical(names(result)[-1], input@design_matrix$Run)
    expect_true(all(is.finite(as.matrix(result[-1]))))
    actual <- methods::new("ProteinQuantitativeData",
        protein_quant_table = result, protein_id_column = "Protein.Group",
        protein_id_table = unique(iq_input[c("Protein.Group", "Protein.Ids")]),
        design_matrix = input@design_matrix, sample_id = "Run", group_id = "group",
        technical_replicate_id = "tech_rep_group", args = list(
            fixture = list(recipe = recipe$id, estimator = "iq::process_long_format",
                estimator_executed = TRUE, normalization = "none")))
    expect_true(methods::validObject(actual))
    expect_identical(unserialize(serialize(actual, NULL)), actual)
    expect_true(file.size(output_path) > 0)
    expect_true(file.size(file.path(directory, "qc-plots.pdf")) > 0)
})

test_that("small histories distinguish actual shared and independent allocations", {
    recipe <- readWorkflowCheckpointRecipe("small")
    independent <- lapply(seq_len(recipe$history_depth), \(i)
        makeWorkflowCheckpointProtein(recipe, generation = i))
    shared <- lapply(seq_len(recipe$history_depth), \(i) {
        object <- independent[[1L]]
        object@args$fixture$metadata_revision <- i
        object
    })
    independent_addresses <- vapply(independent, \(x)
        rlang::obj_address(x@protein_quant_table[[2L]]), character(1))
    shared_addresses <- vapply(shared, \(x)
        rlang::obj_address(x@protein_quant_table[[2L]]), character(1))
    expect_length(unique(independent_addresses), recipe$history_depth)
    expect_length(unique(shared_addresses), 1L)
    expect_true(all(vapply(independent, \(x)
        isTRUE(methods::validObject(x)), logical(1))))
    # Both histories are intentional caller owners here. Do not claim reclaimed
    # memory from evicting only a manager entry while either list remains alive.
    expect_identical(shared[[1L]]@protein_quant_table, independent[[1L]]@protein_quant_table)
    expect_identical(unserialize(serialize(shared, NULL)), shared)
})

test_that("fixture method and design-factor controls retain exact types", {
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    original_design <- object@design_matrix
    changed_design <- original_design
    changed_design$group <- factor(as.character(changed_design$group),
        levels = rev(levels(changed_design$group)))
    expect_identical(as.character(changed_design$group), as.character(original_design$group))
    expect_false(identical(changed_design, original_design))
    expect_identical(rownames(changed_design), rownames(original_design))
    typed <- list(integer = 1L, double = 1, text = "1",
        unordered = factor("1", levels = c("1", "2")),
        ordered = ordered("1", levels = c("1", "2")),
        missing_integer = NA_integer_, missing_character = NA_character_)
    fingerprints <- vapply(typed, \(x) digest::digest(x, algo = "sha256"), character(1))
    expect_length(unique(fingerprints), length(typed))
    expect_identical(unserialize(serialize(typed, NULL)), typed)
    cases <- list(
        same_input_different_normalization = c("none", "cyclicloess", "quantile", "scale"),
        ruv_modes = c("automatic", "manual", "skip"),
        correlation_actions = c("apply", "skip", "reapply"),
        designs = list(original = original_design, changed = changed_design),
        design_sample_reorder = rev(original_design$Run)
    )
    expect_identical(unserialize(serialize(cases, NULL)), cases)
})

test_that("explicit service preflight preserves success and failure receipts", {
    output_root <- Sys.getenv("DIA_FIXTURE_PREFLIGHT_ROOT")
    skip_if(!nzchar(output_root), "requires explicitly enabled disk-backed cgroup service")
    validateWorkflowFixtureRoot(output_root, here::here())
    recipe <- readWorkflowCheckpointRecipe("small")
    limits <- workflowFixtureCgroupReceipt()
    expect_lte(limits$memory_max, 8 * 1024^3)
    expect_identical(limits$swap_max, 0)
    result <- preflightWorkflowCheckpointRecipe("small", output_root, here::here())
    expect_identical(result$status, "passed_resource_preflight")
    receipt_path <- file.path(output_root, "small", "preflight.json")
    receipt <- jsonlite::read_json(receipt_path, simplifyVector = TRUE)
    expect_identical(receipt$status, "passed_resource_preflight")
    expect_equal(receipt$rows_written, recipe$precursor_rows)
    expect_identical(receipt$resources_at_exit$swap_current, 0L)
    before <- digest::digest(file = receipt_path, algo = "sha256")
    expect_error(preflightWorkflowCheckpointRecipe("small", output_root, here::here()),
        "output already exists")
    expect_identical(digest::digest(file = receipt_path, algo = "sha256"), before)

    failure_root <- file.path(output_root, "injected_failure")
    expect_true(dir.create(failure_root))
    failing <- preflightWorkflowCheckpointRecipe
    environment(failing) <- list2env(list(makeWorkflowPrecursorChunk = function(...) {
        stop("injected fixture construction failure", call. = FALSE)
    }), parent = environment(failing))
    expect_error(failing("small", failure_root, here::here()), "injected fixture")
    failed <- jsonlite::read_json(file.path(failure_root, "small", "preflight.json"),
        simplifyVector = TRUE)
    expect_identical(failed$status, "failed")
    expect_match(failed$error, "injected fixture construction failure")
    expect_equal(failed$report_bytes, 0)
})
