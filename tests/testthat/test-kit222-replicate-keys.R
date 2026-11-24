# KIT-222: the design builder's technical-replicate edit could give two different biological samples one key
# (`tech_rep_group` = group + replicate number). Transcriptomics then summed their counts; the other omics treated them
# as one biological sample (the DA block, RUV-III's replicate sets, the integration pairing). The owner's decisions of
# 2026-11-25 (OD-1 to OD-8): the edit never moves another sample's number, refuses in every mode a number another marked
# sample holds and a re-marking of only some of a sample's runs, then renumbers the group 1..m in its order; Assign
# Metadata clears old marks; saves, imports, opened projects and the consumers refuse a merged design.
# Each E row asserts the fixed behaviour and fails on the base commit at the ticket's base value (section 12); the G rows
# hold on both. Each test also prints what it observed (lines starting "KIT222-OBS") so the base values are recorded.
library(testthat)

kit222Obs <- function(id, value) {
  cat(sprintf("KIT222-OBS %s: %s\n", id, paste(value, collapse = " ")))
  invisible(value)
}

# A fresh builder design: runs in builder-row order, nothing assigned.
kit222Design <- function(runs) {
  data.frame(Run = runs, group = NA_character_, factor1 = NA_character_, factor2 = NA_character_,
    factor3 = NA_character_, batch = NA_character_, replicates = NA_integer_, tech_reps = NA_integer_,
    stringsAsFactors = FALSE)
}

kit222Assign <- function(design, runs, group, start = 1L) {
  applyOmicsDesignMetadataAssignment(design, runs, group, "", "", start)
}

kit222Tech <- function(design, runs, mode, number = NULL) {
  result <- applyOmicsDesignTechRepAssignment(design, runs, mode, number)
  expect_true(isTRUE(result$ok), label = sprintf("the edit on %s (%s) is accepted", paste(runs, collapse = ", "), mode))
  result$designMatrix
}

# The saved design's keys: key -> runs (in design order), keys sorted.
kit222Keys <- function(design) {
  saved <- buildOmicsDesignSaveDesignMatrix(design)
  keys <- split(saved$Run, saved$tech_rep_group)
  keys[order(names(keys))]
}

kit222KeyText <- function(keys) {
  paste(sprintf("%s={%s}", names(keys), vapply(keys, paste, character(1), collapse = ",")), collapse = " ")
}

# The partition of the runs into keys, independent of the labels.
kit222Partition <- function(keys) {
  parts <- lapply(unname(keys), sort)
  parts[order(vapply(parts, `[[`, character(1), 1L))]
}

kit222Assays <- function(runs) {
  values <- as.data.frame(matrix(seq_len(2L * length(runs)), nrow = 2L, dimnames = list(NULL, runs)))
  list(LCMS_Pos = cbind(data.frame(LipidName = c("L1", "L2"), LipidClass = "PC"), values))
}

kit222Save <- function(design) {
  buildOmicsDesignSavePayload(design, character(), kit222Assays(design$Run), list(lipid_id_col = "LipidName",
    annotation_col = "LipidClass"), NULL, list(), "~ 0 + group", "lipid_id_col")
}

# The base commit's technical-replicate edit (cf8af245, R/func_omics_design_shared.R:629-655), frozen for G1 and G5.
kit222BaseTechRep <- function(designMatrix, techRepSamples, assignmentMode, manualReplicateNumber = NULL) {
  selected <- which(designMatrix$Run %in% techRepSamples)
  selectedGroups <- unique(designMatrix$group[selected])
  baseReplicateNumber <- switch(assignmentMode,
    lowest = min(designMatrix$replicates[selected], na.rm = TRUE),
    first = designMatrix$replicates[designMatrix$Run == techRepSamples[1]],
    manual = manualReplicateNumber)
  originalReplicateNumbers <- unique(designMatrix$replicates[selected])
  consolidated <- length(originalReplicateNumbers) - 1
  highest <- max(originalReplicateNumbers, na.rm = TRUE)
  toAdjust <- which(designMatrix$group == selectedGroups[1] & !is.na(designMatrix$group) &
    designMatrix$replicates > highest & !designMatrix$Run %in% techRepSamples)
  if (length(toAdjust) > 0) designMatrix$replicates[toAdjust] <- designMatrix$replicates[toAdjust] - consolidated
  designMatrix$replicates[selected] <- baseReplicateNumber
  techRepMapping <- stats::setNames(seq_along(techRepSamples), techRepSamples)
  designMatrix$tech_reps[selected] <- techRepMapping[designMatrix$Run[selected]]
  designMatrix
}

# The four-run example (ticket section 3): count-file order A_1, B_1, A_2, B_2.
kit222E1Design <- function() {
  design <- kit222Assign(kit222Design(c("A_1", "B_1", "A_2", "B_2")), c("A_1", "B_1", "A_2", "B_2"), "Ctrl")
  design <- kit222Tech(design, c("B_2", "B_1"), "first")
  kit222Tech(design, c("A_2", "A_1"), "first")
}

# The consumers' design (E7's layout): Ctrl holds C (S1, one run), A (S2, S4) and B (S3, S5); Trt holds three single
# runs. "Use first selected sample's replicate", B then A, each clicked higher run first.
kit222ConsumerDesign <- function() {
  runs <- paste0("S", 1:8)
  design <- kit222Assign(kit222Design(runs), runs[1:5], "Ctrl")
  design <- kit222Assign(design, runs[6:8], "Trt")
  design <- kit222Tech(design, c("S5", "S3"), "first")
  kit222Tech(design, c("S4", "S2"), "first")
}

kit222ConsumerPartition <- list(c("S1"), c("S2", "S4"), c("S3", "S5"), c("S6"), c("S7"), c("S8"))

# E1's base design written directly (one key Ctrl_3 for the runs of A and B; tech_reps 2, 2, 1, 1), and a Trt group.
kit222MergedRows <- function() {
  data.frame(Run = c("A_1", "B_1", "A_2", "B_2", "T_1", "T_2"), group = rep(c("Ctrl", "Trt"), c(4L, 2L)),
    replicates = c(3L, 3L, 3L, 3L, 1L, 2L), tech_reps = c(2L, 2L, 1L, 1L, NA, NA),
    tech_rep_group = c(rep("Ctrl_3", 4L), "Trt_1", "Trt_2"), stringsAsFactors = FALSE)
}

kit222Data <- function(runs, features = 30L, seed = 222L) {
  withr::with_seed(seed, matrix(stats::rnorm(features * length(runs), 10, 1), nrow = features,
    dimnames = list(sprintf("F%02d", seq_len(features)), runs)))
}

# ------------------------------------------------------------------------------------------------------------------
# The defect on base: each numbering mode merges two biological samples (E1 to E4, E3b).
# ------------------------------------------------------------------------------------------------------------------

test_that("KIT-222 E1: 'Use first selected sample's replicate' on A_1, B_1, A_2, B_2 keeps A and B apart", {
  keys <- kit222Keys(kit222E1Design())
  kit222Obs("E1", kit222KeyText(keys))
  expect_identical(keys, list(Ctrl_1 = c("A_1", "A_2"), Ctrl_2 = c("B_1", "B_2")))
})

test_that("KIT-222 E2: 'Specify replicate number' B=2, C=3, A=1 on sorted rows keeps three samples", {
  runs <- c("A_1", "A_2", "B_1", "B_2", "C_1", "C_2")
  design <- kit222Assign(kit222Design(runs), runs, "Ctrl")
  design <- kit222Tech(design, c("B_1", "B_2"), "manual", 2)
  design <- kit222Tech(design, c("C_1", "C_2"), "manual", 3)
  design <- kit222Tech(design, c("A_1", "A_2"), "manual", 1)
  keys <- kit222Keys(design)
  kit222Obs("E2", kit222KeyText(keys))
  expect_identical(keys, list(Ctrl_1 = c("A_1", "A_2"), Ctrl_2 = c("B_1", "B_2"), Ctrl_3 = c("C_1", "C_2")))
})

# KIT-209's U8 corpus (trx_60k_200; origin/run-kit209-memory tools/profiling/kit_omics_corpus.R:53-63 and :265): two
# groups (RPMI, Sera), n samples per group, two runs each, design order by replicate then group; the count columns, and
# so the transcriptomics builder's rows, in the permutation `withr::with_seed(20261108L, sample.int(4n))`. The journey
# (.github/k209/k209-journey.R:279-302) assigns each group from replicate 1, then marks each sample's two runs with
# "Specify replicate number" b, b = 1..n, group by group.
kit222Corpus <- function(per_group) {
  biological <- expand.grid(replicate = seq_len(per_group), group = c("RPMI", "Sera"), stringsAsFactors = FALSE)
  biological <- biological[order(biological$replicate, biological$group), , drop = FALSE]
  biological$sample <- sprintf("%s_%02d", biological$group, biological$replicate)
  design <- biological[rep(seq_len(nrow(biological)), each = 2L), , drop = FALSE]
  design$run <- sprintf("%s_r%d", design$sample, rep(1:2, times = nrow(biological)))
  rownames(design) <- NULL
  order <- withr::with_seed(20261108L, sample.int(nrow(design)))
  list(design = design, rows = design$run[order])
}

kit222CorpusMarked <- function(per_group) {
  corpus <- kit222Corpus(per_group)
  design <- kit222Design(corpus$rows)
  for (group in c("RPMI", "Sera")) design <- kit222Assign(design, corpus$design$run[corpus$design$group == group], group)
  for (group in c("RPMI", "Sera")) for (b in seq_len(per_group)) {
    runs <- corpus$design$run[corpus$design$group == group & corpus$design$replicate == b]
    design <- kit222Tech(design, runs, "manual", b)
  }
  list(corpus = corpus, design = design)
}

kit222CorpusExpected <- function(corpus) {
  samples <- unique(corpus$design$sample)
  keys <- sub("_0?([0-9]+)$", "_\\1", samples)
  expected <- stats::setNames(lapply(samples, function(s) sort(corpus$design$run[corpus$design$sample == s])), keys)
  expected[order(names(expected))]
}

kit222MergedText <- function(keys, corpus) {
  merged <- keys[lengths(keys) > 2L]
  sample_of <- stats::setNames(corpus$design$sample, corpus$design$run)
  vapply(names(merged), function(key) sprintf("%s=%s", key, paste(unique(sample_of[merged[[key]]]), collapse = "+")),
    character(1))
}

test_that("KIT-222 E3: the 100-sample corpus marked with 'Specify replicate number' gives 100 samples", {
  marked <- kit222CorpusMarked(50L)
  keys <- kit222Keys(marked$design)
  kit222Obs("E3 keys", length(keys))
  kit222Obs("E3 merged", kit222MergedText(keys, marked$corpus))
  expect_identical(length(keys), 100L)
  expect_identical(lapply(keys, sort), kit222CorpusExpected(marked$corpus))
})

test_that("KIT-222 E3b: summing the corpus's technical replicates adds only each sample's own two runs", {
  skip_if_not_installed("edgeR")
  marked <- kit222CorpusMarked(50L)
  saved <- buildOmicsDesignSaveDesignMatrix(marked$design)
  runs <- saved$Run
  # One gene whose count is the run's column index.
  table <- cbind(data.frame(GeneID = "G1", stringsAsFactors = FALSE),
    as.data.frame(matrix(seq_along(runs), nrow = 1L, dimnames = list(NULL, runs))))
  summed <- trxSumTechnicalReplicates(list(counts = table), saved)
  units <- setdiff(names(summed$count_data$counts), "GeneID")
  sample_of <- stats::setNames(marked$corpus$design$sample, marked$corpus$design$run)
  contributors <- lapply(summed$summed, function(r) unique(sample_of[r]))
  several <- contributors[lengths(contributors) > 1L]
  kit222Obs("E3b units", length(units))
  kit222Obs("E3b summed from two samples", vapply(names(several), function(u) sprintf("%s=%s (count %s = %s)", u,
    paste(several[[u]], collapse = "+"), summed$count_data$counts[[u]],
    paste(match(summed$summed[[u]], runs), collapse = "+")), character(1)))
  expect_identical(length(units), 100L)
  expect_true(all(lengths(summed$summed) == 2L))
  expect_true(all(lengths(contributors) == 1L))
  for (unit in units) expect_identical(summed$count_data$counts[[unit]], sum(match(summed$summed[[unit]], runs)))
})

test_that("KIT-222 E3c (G): the 24-run corpus gives 12 samples on base and on the fix", {
  marked <- kit222CorpusMarked(6L)
  keys <- kit222Keys(marked$design)
  kit222Obs("E3c keys", length(keys))
  expect_identical(lapply(keys, sort), kit222CorpusExpected(marked$corpus))
})

test_that("KIT-222 E4: 'Use lowest replicate number' with three interleaved runs keeps B and C apart", {
  runs <- c("A_1", "A_2", "B_1", "A_3", "B_2", "C_1", "C_2")
  design <- kit222Assign(kit222Design(runs), runs, "Ctrl")
  design <- kit222Tech(design, c("B_1", "B_2"), "lowest")
  design <- kit222Tech(design, c("A_1", "A_2", "A_3"), "lowest")
  design <- kit222Tech(design, c("C_1", "C_2"), "lowest")
  keys <- kit222Keys(design)
  kit222Obs("E4", kit222KeyText(keys))
  expect_identical(keys, list(Ctrl_1 = c("A_1", "A_2", "A_3"), Ctrl_2 = c("B_1", "B_2"), Ctrl_3 = c("C_1", "C_2")))
})

# ------------------------------------------------------------------------------------------------------------------
# E7: each omic's builder and the proteomics workbook applet, through their own inputs.
# ------------------------------------------------------------------------------------------------------------------

# The clicks of the consumers' design through a builder's inputs (`set` writes inputs by name), then save.
kit222BuilderClicks <- function(set) {
  set(new_factor = "Ctrl"); set(add_factor = 1L)
  set(new_factor = "Trt"); set(add_factor = 2L)
  set(selected_runs = paste0("S", 1:5), factor1_select = "Ctrl", factor2_select = "", factor3_select = "",
    replicate_start = 1L)
  set(assign_metadata = 1L)
  set(selected_runs = paste0("S", 6:8), factor1_select = "Trt", replicate_start = 1L)
  set(assign_metadata = 2L)
  set(tech_rep_samples = c("S5", "S3"), tech_rep_assignment_mode = "first")
  set(assign_tech_reps = 1L)
  set(tech_rep_samples = c("S4", "S2"))
  set(assign_tech_reps = 2L)
  set(contrast_group1 = "Trt", contrast_group2 = "Ctrl")
  set(add_contrast = 1L)
  set(formula_string = "~ 0 + group")
  set(save_results = 1L)
}

kit222SharedBuilderSave <- function(omic) {
  runs <- paste0("S", 1:8)
  input <- switch(omic,
    lipidomics = list(assays = module_ci_lipid_design_assays(layout = "single", samples = runs),
      mapping = module_ci_lipid_design_column_mapping(runs)),
    metabolomics = list(assays = module_ci_metab_design_assays(layout = "single", samples = runs),
      mapping = module_ci_metab_design_column_mapping(runs)),
    transcriptomics = list(assays = list(counts = cbind(data.frame(GeneID = c("G1", "G2"), GeneName = c("N1", "N2")),
      as.data.frame(matrix(1:16, nrow = 2L, dimnames = list(NULL, runs))))),
      mapping = list(lipid_id_col = "GeneID", annotation_col = "GeneName", sample_columns = runs)))
  saved <- NULL
  suppressMessages(shiny::testServer(mod_omic_design_builder_server, args = list(data_tbl = function() input$assays,
    config_list = shiny::reactive(list(deAnalysisParameters = list(formula_string = "~ 0 + group"))),
    column_mapping = shiny::reactive(input$mapping), config = omicGuiConfig(omic)), {
    session$flushReact()
    kit222BuilderClicks(function(...) session$setInputs(...))
    saved <<- shiny::isolate(session$returned())
  }))
  saved
}

kit222SavedPartition <- function(saved) {
  kit222Partition(split(saved$design_matrix$Run, saved$design_matrix$tech_rep_group))
}

test_that("KIT-222 E7: the lipidomics, metabolomics and transcriptomics builders keep A and B apart", {
  local_mocked_bindings(showNotification = function(...) invisible(NULL), .package = "shiny")
  for (omic in c("lipidomics", "metabolomics", "transcriptomics")) {
    saved <- kit222SharedBuilderSave(omic)
    expect_type(saved, "list")
    kit222Obs(paste("E7", omic), kit222KeyText(split(saved$design_matrix$Run, saved$design_matrix$tech_rep_group)))
    expect_identical(kit222SavedPartition(saved), kit222ConsumerPartition, label = paste(omic, "saved partition"))
  }
})

# A DIA owner holding a freshly imported peptide table (test-prot-04b-design-server-shared.R's owner).
kit222ProtOwner <- function(label, runs, envir = parent.frame()) {
  owner <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(.local_envir = envir), label, "DIA"),
    WorkflowState$new(audit_enabled = FALSE))
  owner$.context$reference_states <- TRUE
  withr::defer(closeDIAStateAdapter(owner), envir = envir)
  rows <- sprintf("%s\tP1\t%d", runs, seq_along(runs))
  table <- prepareDIAParserTable(readr::read_tsv(I(paste0("Run\tProtein.Group\tPrecursor.Normalised\n",
    paste(rows, collapse = "\n"), "\n")), show_col_types = FALSE))
  mapping <- list(run_col = "Run", protein_col = "Protein.Group", quantity_col = "Precursor.Normalised")
  settings <- list(deAnalysisParameters = list(formula_string = "~ 0 + group"))
  references <- prepareDIAGUIReferences(owner, list(data_tbl = table, data_cln = table, data_type = "peptide",
    column_mapping = mapping, config_list = settings))
  owner$.context$design_adoption_pending <- NULL
  list(manager = owner, reference = references$data_tbl, mapping = mapping, settings = settings)
}

# The proteomics builder over a DIA reference: `clicks(set)` drives it; returns the saved result and the notifications.
kit222ProtBuilder <- function(runs, clicks) {
  owner <- kit222ProtOwner("kit222-prot", runs)
  reference <- shiny::reactiveVal(owner$reference)
  provider <- function() resolveDIAOwnedValue(reference())
  attr(provider, "dia_reference") <- reference
  seen <- new.env(parent = emptyenv())
  seen$notes <- character()
  local_mocked_bindings(showNotification = function(ui, ..., type = "default") {
    seen$notes <- c(seen$notes, sprintf("%s: %s", type, paste(as.character(ui), collapse = " ")))
    invisible(NULL)
  }, .package = "shiny")
  suppressMessages(shiny::testServer(mod_prot_design_builder_server, args = list(data_tbl = provider,
    config_list = shiny::reactiveVal(owner$settings), column_mapping = shiny::reactiveVal(owner$mapping)), {
    session$userData$dia_capacity_manager <- owner$manager
    session$flushReact()
    clicks(function(...) session$setInputs(...))
    seen$saved <- shiny::isolate(session$returned())
  }))
  list(saved = seen$saved, notes = seen$notes)
}

test_that("KIT-222 E7: the proteomics builder keeps A and B apart", {
  result <- kit222ProtBuilder(paste0("S", 1:8), kit222BuilderClicks)
  expect_type(result$saved, "list")
  kit222Obs("E7 proteomics", kit222KeyText(split(result$saved$design_matrix$Run, result$saved$design_matrix$tech_rep_group)))
  expect_identical(kit222SavedPartition(result$saved), kit222ConsumerPartition)
})

test_that("KIT-222 E7: the proteomics builder refuses a number another sample holds with the shared text", {
  result <- kit222ProtBuilder(paste0("S", 1:4), function(set) {
    set(new_factor = "Ctrl"); set(add_factor = 1L)
    set(selected_runs = paste0("S", 1:4), factor1_select = "Ctrl", factor2_select = "", factor3_select = "",
      replicate_start = 1L)
    set(assign_metadata = 1L)
    set(tech_rep_samples = c("S1", "S2"), tech_rep_assignment_mode = "manual", manual_replicate_number = 1)
    set(assign_tech_reps = 1L)
    set(tech_rep_samples = c("S3", "S4"))
    set(assign_tech_reps = 2L)
  })
  kit222Obs("E7 proteomics number held", result$notes)
  expect_length(result$notes, 2L)
  expect_match(result$notes[[2L]], "^warning: ")
  expect_match(result$notes[[2L]], "S1, S2", fixed = TRUE)
  expect_no_match(result$notes[[2L]], "same group", fixed = TRUE)
})

test_that("KIT-222 E7: the proteomics workbook applet keeps A and B apart", {
  skip_if_not(workflowOwnerCapable())
  source_dir <- withr::local_tempdir("kit222-applet-source-")
  output_dir <- withr::local_tempdir("kit222-applet-output-")
  report_path <- file.path(withr::local_tempdir("kit222-applet-data-"), "report.tsv")
  runs <- paste0("S", 1:8)
  utils::write.table(data.frame(Run = rep(runs, each = 2L), Protein.Group = rep(c("P1", "P2"), 8L),
    Stripped.Sequence = rep(c("AAK", "CCK"), 8L), Precursor.Quantity = seq_len(16L) * 10, Q.Value = 0.001),
    report_path, sep = "\t", row.names = FALSE, quote = FALSE)
  applet <- new.env(parent = emptyenv())
  local_mocked_bindings(
    runApp = function(appDir, ...) {
      shiny::testServer(appDir, {
        session$flushReact()
        kit222BuilderClicks(function(...) {
          values <- list(...)
          names(values) <- paste0("designBuilder-", names(values))
          do.call(session$setInputs, values)
        })
      })
      applet$stopped
    },
    stopApp = function(returnValue = invisible()) {
      applet$stopped <- returnValue
      invisible(NULL)
    },
    updateSelectizeInput = function(...) invisible(NULL),
    showNotification = function(...) invisible(NULL),
    .package = "shiny")
  workbook <- local({
    project_dirs <- list(proteomics_exp = list(source_dir = source_dir, output_dir = output_dir))
    data_tbl <- vroom::vroom(report_path, show_col_types = FALSE)
    config_list <- list(deAnalysisParameters = list(formula_string = "~ 0 + group"))
    returned <- suppressMessages(RunApplet("designMatrix", "proteomics", "exp", force = TRUE))
    environment()
  })
  returned <- workbook$returned
  expect_type(returned, "list")
  kit222Obs("E7 applet", kit222KeyText(split(returned$design_matrix$Run, returned$design_matrix$tech_rep_group)))
  expect_identical(kit222SavedPartition(returned), kit222ConsumerPartition)
  written <- utils::read.delim(file.path(source_dir, "design_matrix.tab"), stringsAsFactors = FALSE)
  expect_identical(kit222Partition(split(written$Run, written$tech_rep_group)), kit222ConsumerPartition)
})

# ------------------------------------------------------------------------------------------------------------------
# What the consumers do with the builder's design (the consumers' design through the shared edits).
# ------------------------------------------------------------------------------------------------------------------

test_that("KIT-222 consumers: transcriptomics sums only each sample's own runs", {
  skip_if_not_installed("edgeR")
  saved <- buildOmicsDesignSaveDesignMatrix(kit222ConsumerDesign())
  table <- cbind(data.frame(GeneID = "G1", stringsAsFactors = FALSE),
    as.data.frame(matrix(seq_len(8L), nrow = 1L, dimnames = list(NULL, saved$Run))))
  summed <- trxSumTechnicalReplicates(list(counts = table), saved)
  kit222Obs("trx sum", vapply(names(summed$summed), function(u) sprintf("%s={%s}=%s", u,
    paste(summed$summed[[u]], collapse = ","), summed$count_data$counts[[u]]), character(1)))
  expect_identical(kit222Partition(summed$summed), kit222ConsumerPartition)
  # A's unit is S2 + S4 (2 + 4) and B's is S3 + S5 (3 + 5); on base one unit holds S2 + S3 + S4 + S5 (14).
  sums <- vapply(summed$summed, function(r) sum(match(r, saved$Run)), numeric(1))
  expect_setequal(unname(sums[lengths(summed$summed) == 2L]), c(6, 8))
})

test_that("KIT-222 consumers: the DA model's unit (the lipid fit's block) keeps A and B apart", {
  skip_if_not_installed("limma")
  saved <- as.data.frame(buildOmicsDesignSaveDesignMatrix(kit222ConsumerDesign()))
  blocks <- new.env(parent = emptyenv())
  original <- resolveProtTechnicalReplicateBlock
  local_mocked_bindings(resolveProtTechnicalReplicateBlock = function(...) {
    result <- original(...)
    blocks$block <- result$block
    result
  })
  fit <- suppressWarnings(fitLipidDaAssayModel(kit222Data(saved$Run), saved, "~ 0 + group", "Run", "group",
    "tech_rep_group"))
  kit222Obs("DA block", paste(saved$Run, blocks$block, sep = "="))
  expect_true(isTRUE(fit$uses_duplicate_correlation))
  expect_identical(kit222Partition(split(saved$Run, blocks$block)), kit222ConsumerPartition)
})

test_that("KIT-222 consumers: RUV-III's replicate sets keep A and B apart", {
  saved <- as.data.frame(buildOmicsDesignSaveDesignMatrix(kit222ConsumerDesign()))
  sets <- resolveProtRuvReplicateBlocks(saved, "Run", "group", "tech_rep_group", "group")
  kit222Obs("RUV-III sets", paste(sets$design$Run, sets$design$.ruv_replicate_block, sep = "="))
  expect_identical(sets$structure, "technical_replicates")
  expect_identical(kit222Partition(split(sets$design$Run, sets$design$.ruv_replicate_block)), kit222ConsumerPartition)
})

test_that("KIT-222 consumers: the integration averages each sample's own runs and pairs each name with one sample", {
  saved <- as.data.frame(buildOmicsDesignSaveDesignMatrix(kit222ConsumerDesign()))
  identity <- c(S1 = "C", S2 = "A", S3 = "B", S4 = "A", S5 = "B", S6 = "D", S7 = "E", S8 = "F")
  column <- integrationInputBiologicalSampleColumn(saved, "tech_rep_group", "Run")
  expect_identical(column, "tech_rep_group")
  values <- matrix(seq_len(8L), nrow = 1L, dimnames = list("F1", saved$Run))
  averaged <- integrationAverageRuns(values, saved[[column]])
  kit222Obs("integration averages", paste(colnames(averaged), averaged[1L, ], sep = "="))
  expect_identical(ncol(averaged), 6L)
  # Transcriptomics of the same six samples, one run each, numbered in its own row order C, A, B.
  trx <- data.frame(sample = c("C", "A", "B", "D", "E", "F"), group = rep(c("Ctrl", "Trt"), each = 3L),
    biological_sample = c("Ctrl_1", "Ctrl_2", "Ctrl_3", "Trt_1", "Trt_2", "Trt_3"), stringsAsFactors = FALSE)
  design <- rbind(data.frame(omic = "proteomics", biological_sample = saved[[column]], group = saved$group),
    data.frame(omic = "transcriptomics", biological_sample = trx$biological_sample, group = trx$group))
  aligned <- integrationAlignSamples(design, c(prot = "proteomics", trx = "transcriptomics"), c("prot", "trx"),
    list(prot = unique(saved[[column]]), trx = trx$biological_sample), "drop_incomplete")
  pairs <- vapply(aligned$shared, function(name) paste(unique(identity[saved$Run[saved[[column]] == name]]), "with",
    trx$sample[trx$biological_sample == name]), character(1))
  kit222Obs("integration pairs", paste(aligned$shared, pairs, sep = ": "))
  expect_identical(length(aligned$shared), 6L)
  expect_identical(unname(pairs), c("C with C", "A with A", "B with B", "D with D", "E with E", "F with F"))
})

# ------------------------------------------------------------------------------------------------------------------
# The refusals and the renumbering (fail on base: base has neither).
# ------------------------------------------------------------------------------------------------------------------

test_that("KIT-222 E5: a typed number held by another marked sample is refused in 'Specify replicate number'", {
  runs <- c("A_1", "A_2", "B_1", "B_2")
  design <- kit222Tech(kit222Assign(kit222Design(runs), runs, "Ctrl"), c("A_1", "A_2"), "manual", 1)
  second <- applyOmicsDesignTechRepAssignment(design, c("B_1", "B_2"), "manual", 1)
  kit222Obs("E5", c(second$ok, kit222KeyText(kit222Keys(second$designMatrix))))
  expect_false(isTRUE(second$ok))
  expect_identical(second$code, "number_held")
  expect_identical(second$designMatrix, design)
  expect_match(second$message, "A_1, A_2", fixed = TRUE)
})

test_that("KIT-222 E6: re-marking only some of a sample's runs is refused", {
  runs <- c("A_1", "A_2", "A_3")
  design <- kit222Tech(kit222Assign(kit222Design(runs), runs, "Ctrl"), c("A_1", "A_2"), "lowest")
  second <- applyOmicsDesignTechRepAssignment(design, c("A_1", "A_3"), "lowest")
  kit222Obs("E6", c(second$ok, "tech_reps", second$designMatrix$tech_reps))
  expect_false(isTRUE(second$ok))
  expect_identical(second$code, "partial_sample")
  expect_identical(second$designMatrix, design)
  for (run in runs) expect_match(second$message, run, fixed = TRUE)
})

test_that("KIT-222 E14: 'Specify replicate number' A=2, then 'Use lowest replicate number' on B is refused", {
  runs <- c("A_1", "B_1", "A_2", "B_2")
  design <- kit222Tech(kit222Assign(kit222Design(runs), runs, "Ctrl"), c("A_1", "A_2"), "manual", 2)
  second <- applyOmicsDesignTechRepAssignment(design, c("B_1", "B_2"), "lowest")
  kit222Obs("E14", c(second$ok, kit222KeyText(kit222Keys(second$designMatrix)), "tech_reps", second$designMatrix$tech_reps))
  expect_false(isTRUE(second$ok))
  expect_identical(second$code, "number_held")
  expect_match(second$message, "A_1, A_2", fixed = TRUE)
  # B_1 shares A's number unmarked (1 after the renumbering): the save refuses it (K128-2).
  expect_identical(design$replicates[design$Run == "B_1"], design$replicates[design$Run == "A_1"])
  refused <- kit222Save(design)
  expect_type(refused, "character")
  expect_match(refused, "B_1", fixed = TRUE)
})

test_that("KIT-222 E15: 'Specify replicate number' A=2, then B=2: two samples, A shown as 1", {
  runs <- c("A_1", "B_1", "A_2", "B_2")
  design <- kit222Tech(kit222Assign(kit222Design(runs), runs, "Ctrl"), c("A_1", "A_2"), "manual", 2)
  design <- kit222Tech(design, c("B_1", "B_2"), "manual", 2)
  keys <- kit222Keys(design)
  kit222Obs("E15", kit222KeyText(keys))
  expect_identical(keys, list(Ctrl_1 = c("A_1", "A_2"), Ctrl_2 = c("B_1", "B_2")))
})

test_that("KIT-222 E16: re-assigning a marked run to another group clears its mark; the save refuses the shared number", {
  runs <- c("A_1", "A_2", "A_3", "C_1", "C_2")
  design <- kit222Assign(kit222Design(runs), runs[1:3], "Ctl")
  design <- kit222Assign(design, runs[4:5], "Trt")
  design <- kit222Tech(design, runs[1:3], "lowest")
  design <- kit222Tech(design, runs[4:5], "lowest")
  design <- kit222Assign(design, "A_3", "Trt")
  kit222Obs("E16 tech_reps", paste(design$Run, design$tech_reps, sep = "="))
  expect_true(is.na(design$tech_reps[design$Run == "A_3"]))
  refused <- kit222Save(design)
  kit222Obs("E16 save", if (is.character(refused)) refused else kit222KeyText(split(refused$design_matrix$Run,
    refused$design_matrix$tech_rep_group)))
  expect_type(refused, "character")
  expect_match(refused, "A_3, C_1, C_2", fixed = TRUE)
})

test_that("KIT-222 E16b: the proteomics builder's re-assignment clears the mark; the proteomics save refuses", {
  result <- kit222ProtBuilder(paste0("S", 1:5), function(set) {
    set(new_factor = "Ctl"); set(add_factor = 1L)
    set(new_factor = "Trt"); set(add_factor = 2L)
    set(selected_runs = paste0("S", 1:3), factor1_select = "Ctl", factor2_select = "", factor3_select = "",
      replicate_start = 1L)
    set(assign_metadata = 1L)
    set(selected_runs = paste0("S", 4:5), factor1_select = "Trt", replicate_start = 1L)
    set(assign_metadata = 2L)
    set(tech_rep_samples = paste0("S", 1:3), tech_rep_assignment_mode = "lowest")
    set(assign_tech_reps = 1L)
    set(tech_rep_samples = paste0("S", 4:5))
    set(assign_tech_reps = 2L)
    set(selected_runs = "S3", factor1_select = "Trt", replicate_start = 1L)
    set(assign_metadata = 3L)
    set(formula_string = "~ 0 + group")
    set(save_results = 1L)
  })
  kit222Obs("E16b", c(if (is.null(result$saved)) "no result" else kit222KeyText(split(result$saved$design_matrix$Run,
    result$saved$design_matrix$tech_rep_group)), result$notes))
  expect_null(result$saved)
  expect_true(any(grepl("^error: .*S3, S4, S5", result$notes)))
})

test_that("KIT-222 E8: the shared save refuses a design whose key holds two marking actions", {
  rows <- kit222MergedRows()
  refused <- kit222Save(rows[, setdiff(names(rows), "tech_rep_group")])
  kit222Obs("E8", if (is.character(refused)) refused else class(refused))
  expect_type(refused, "character")
  expect_match(refused, "Ctrl_3", fixed = TRUE)
  expect_match(refused, "A_1, B_1, A_2, B_2", fixed = TRUE)
})

kit222D2Rows <- function() {
  # Sequence 8 on base: A_1, B (one run), A_2, C (one run); 'first' with A_2 clicked first moved C into A's key.
  data.frame(Run = c("A_1", "B", "A_2", "C", "T_1", "T_2"), group = rep(c("Ctrl", "Trt"), c(4L, 2L)),
    replicates = c(3L, 2L, 3L, 3L, 1L, 2L), tech_reps = c(2L, NA, 1L, NA, NA, NA), stringsAsFactors = FALSE)
}

test_that("KIT-222 E9: the proteomics save refuses D1, D2 and an unmarked shared number", {
  d1 <- kit222MergedRows()[, c("Run", "group", "replicates", "tech_reps")]
  # OD-5: two Assign clicks to one group from replicate 1 (A_1 and A_3 share 1, A_2 and A_4 share 2), unmarked.
  od5 <- data.frame(Run = c("A_1", "A_2", "A_3", "A_4", "T_1", "T_2"), group = rep(c("Ctrl", "Trt"), c(4L, 2L)),
    replicates = c(1L, 2L, 1L, 2L, 1L, 2L), tech_reps = NA_integer_, stringsAsFactors = FALSE)
  for (case in list(list(id = "D1", design = d1, runs = "A_1"), list(id = "D2", design = kit222D2Rows(), runs = "C"),
    list(id = "OD-5", design = od5, runs = "A_1"))) {
    owner <- kit222ProtOwner(paste0("kit222-e9-", case$id), case$design$Run)
    outcome <- tryCatch(suppressMessages(withDIAStateScope(owner$manager, function() {
      buildProtDesignSaveResultsPayload(case$design, character(), resolveDIAOwnedValue(owner$reference, owner$manager),
        NULL, owner$settings, "~ 0 + group")
    })), error = identity)
    kit222Obs(paste("E9 save", case$id), if (is.character(outcome)) outcome else if (inherits(outcome, "condition"))
      paste(class(outcome)[[1L]], conditionMessage(outcome)) else kit222KeyText(split(outcome$design_matrix$Run,
      outcome$design_matrix$tech_rep_group)))
    expect_type(outcome, "character")
    if (is.character(outcome)) expect_match(outcome, case$runs, fixed = TRUE)
  }
})

test_that("KIT-222 E9: the proteomics checkpoint (builder commit and import) refuses D1 and D2", {
  d1 <- kit222MergedRows()[, c("Run", "group", "replicates", "tech_reps")]
  for (case in list(list(id = "D1", design = d1), list(id = "D2", design = kit222D2Rows()))) {
    workflow <- new.env(parent = emptyenv())
    workflow$design_matrix <- case$design
    outcome <- tryCatch(suppressMessages(buildProtDesignStateCheckpoint(workflow, "DIA", "kit222")), error = identity)
    kit222Obs(paste("E9 checkpoint", case$id), class(outcome))
    expect_s3_class(outcome, "omics_design_replicate_key_collision")
    if (inherits(outcome, "condition")) expect_match(conditionMessage(outcome), "A_1", fixed = TRUE)
  }
})

test_that("KIT-222 E10: the lipid import refuses a design file carrying D1 or D2", {
  for (case in list(list(id = "D1", tech = c(1L, 1L)), list(id = "D2", tech = c(1L, NA)))) {
    payload <- module_ci_lipid_design_payload(layout = "lc_pair")
    payload$design_matrix$replicates[payload$design_matrix$Run %in% c("WT_1", "WT_2")] <- 1L
    payload$design_matrix$tech_reps[payload$design_matrix$Run %in% c("WT_1", "WT_2")] <- case$tech
    root <- module_ci_lipid_design_write_import_pack(withr::local_tempdir("kit222-import-"), payload)
    workflow <- module_ci_lipid_design_workflow(payload)
    config <- omicGuiConfig("lipidomics")
    outcome <- tryCatch(shiny::withReactiveDomain(omicsOwnerSession(), {
      imported <- readOmicsDesignImport(root, config)
      commitOmicsDesignImport(workflow, imported, root, config)
    }), error = identity)
    kit222Obs(paste("E10", case$id), class(outcome))
    expect_s3_class(outcome, "omics_design_refused")
    if (inherits(outcome, "condition")) expect_match(conditionMessage(outcome), "WT_1, WT_2", fixed = TRUE)
  }
})

test_that("KIT-222 E11: the consumers refuse a merged design; the display opts out", {
  rows <- kit222MergedRows()
  class <- "omics_design_replicate_key_collision"
  outcome <- function(expr) tryCatch(expr, error = identity)
  checks <- list(
    resolver = outcome(resolveProtTechnicalReplicateBlock(rows, "group", "tech_rep_group")),
    ruv_iii = outcome(resolveProtRuvReplicateBlocks(rows, "Run", "group", "tech_rep_group", "group")),
    lipid_fit = outcome(suppressWarnings(fitLipidDaAssayModel(kit222Data(rows$Run), rows, "~ 0 + group", "Run",
      "group", "tech_rep_group"))),
    metab_fit = outcome(suppressWarnings(fitMetabDaAssayModel(kit222Data(rows$Run), rows, "~ 0 + group", "Run",
      "group", "tech_rep_group"))),
    integration = outcome(integrationInputBiologicalSampleColumn(rows, "tech_rep_group", "Run")))
  if (requireNamespace("edgeR", quietly = TRUE)) {
    table <- cbind(data.frame(GeneID = "G1"), as.data.frame(matrix(1:6, nrow = 1L, dimnames = list(NULL, rows$Run))))
    checks$trx_sum <- outcome(trxSumTechnicalReplicates(list(counts = table), rows))
  }
  for (name in names(checks)) {
    kit222Obs(paste("E11", name), class(checks[[name]]))
    expect_s3_class(checks[[name]], class)
  }
  if (!is.null(checks$trx_sum)) expect_s3_class(checks$trx_sum, "trx_error")
  display <- outcome(resolveProtTechnicalReplicateBlock(rows, "group", "tech_rep_group", check = FALSE))
  expect_type(display, "list")
})

test_that("KIT-222 E11: peptide imputation refuses a merged design", {
  rows <- kit222MergedRows()
  peptides <- data.frame(Protein.Ids = "P1", Stripped.Sequence = rep(c("SEQ1", "SEQ2"), each = 6L),
    Proteotypic.Sequence = rep(c("SEQ1", "SEQ2"), each = 6L), Run = rep(rows$Run, 2L),
    Precursor.Quantity = c(10, NA, 30, 40, 50, 60, 11, 21, 31, 41, 51, 61), Precursor.Normalised = 1,
    Q.Value = 0.001, Global.Q.Value = 0.001, Global.PG.Q.Value = 0.001, stringsAsFactors = FALSE)
  object <- methods::new("PeptideQuantitativeData", peptide_data = peptides,
    peptide_matrix = matrix(1, nrow = 2L, ncol = 6L, dimnames = list(c("SEQ1", "SEQ2"), rows$Run)),
    protein_id_column = "Protein.Ids", peptide_sequence_column = "Stripped.Sequence", q_value_column = "Q.Value",
    global_q_value_column = "Global.Q.Value", proteotypic_peptide_sequence_column = "Proteotypic.Sequence",
    raw_quantity_column = "Precursor.Quantity", norm_quantity_column = "Precursor.Normalised", is_logged_data = FALSE,
    design_matrix = rows, sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group",
    args = list())
  result <- tryCatch(suppressMessages(peptideMissingValueImputation(object, proportion_missing_values = 0.5)),
    error = identity)
  kit222Obs("E11 peptide imputation", class(result))
  expect_s3_class(result, "omics_design_replicate_key_collision")
})

# ------------------------------------------------------------------------------------------------------------------
# Guards: unchanged where base did not merge (hold on base and on the fix).
# ------------------------------------------------------------------------------------------------------------------

kit222Permutations <- function(x) {
  if (length(x) <= 1L) return(list(x))
  out <- list()
  for (i in seq_along(x)) for (rest in kit222Permutations(x[-i])) out[[length(out) + 1L]] <- c(x[i], rest)
  out
}

# Marks every sample (two runs each) of one group in `order`, through both edits; TRUE when the designs are identical.
kit222SameAsBase <- function(rows, samples, order, mode, typed = NULL) {
  base <- head <- kit222Assign(kit222Design(rows), rows, "G")
  for (i in seq_along(order)) {
    runs <- rows[startsWith(rows, paste0(order[[i]], "_"))]
    if (identical(mode, "first")) runs <- runs[order(base$replicates[match(runs, base$Run)])]
    number <- if (identical(mode, "manual")) typed[[i]] else NULL
    base <- kit222BaseTechRep(base, runs, mode, number)
    head <- applyOmicsDesignTechRepAssignment(head, runs, mode, number)$designMatrix
  }
  identical(base, head)
}

test_that("KIT-222 G1: merge-free designs are byte-identical to base (exhaustive at 2 and 3 samples, 200 random at 50)", {
  failures <- character()
  for (k in 2:3) {
    samples <- LETTERS[seq_len(k)]
    runs <- as.vector(t(outer(samples, 1:2, paste, sep = "_")))
    for (rows in kit222Permutations(runs)) for (order in kit222Permutations(samples)) {
      for (mode in c("lowest", "first")) if (!kit222SameAsBase(rows, samples, order, mode)) {
        failures <- c(failures, sprintf("%s %s %s", mode, paste(rows, collapse = ","), paste(order, collapse = ",")))
      }
      # 'Specify replicate number' typed 1..k in increasing order: sample order[[i]] typed i.
      if (!kit222SameAsBase(rows, samples, order, "manual", seq_len(k))) {
        failures <- c(failures, sprintf("manual %s %s", paste(rows, collapse = ","), paste(order, collapse = ",")))
      }
    }
  }
  samples <- sprintf("S%02d", 1:50)
  runs <- as.vector(t(outer(samples, 1:2, paste, sep = "_")))
  withr::with_seed(222L, for (i in seq_len(200L)) {
    rows <- runs[sample.int(length(runs))]
    mode <- c("lowest", "first", "manual")[[(i %% 3L) + 1L]]
    if (!kit222SameAsBase(rows, samples, samples, mode, seq_along(samples))) {
      failures <- c(failures, sprintf("random %d %s", i, mode))
    }
  })
  kit222Obs("G1 failures", length(failures))
  expect_identical(failures, character())
})

test_that("KIT-222 G4: duplicate-injection proteomics against single-run transcriptomics pairs six of six", {
  runs <- as.vector(t(outer(sprintf("S%d", 1:6), 1:2, paste, sep = "_")))
  prot <- kit222Assign(kit222Design(runs), runs, "Ctrl")
  for (s in sprintf("S%d", 1:6)) prot <- kit222Tech(prot, paste0(s, c("_1", "_2")), "lowest")
  prot <- as.data.frame(buildOmicsDesignSaveDesignMatrix(prot))
  trx <- as.data.frame(buildOmicsDesignSaveDesignMatrix(kit222Assign(kit222Design(sprintf("S%d", 1:6)),
    sprintf("S%d", 1:6), "Ctrl")))
  design <- rbind(data.frame(omic = "proteomics", biological_sample = prot$tech_rep_group, group = prot$group),
    data.frame(omic = "transcriptomics", biological_sample = trx$tech_rep_group, group = trx$group))
  aligned <- integrationAlignSamples(design, c(prot = "proteomics", trx = "transcriptomics"), c("prot", "trx"),
    list(prot = unique(prot$tech_rep_group), trx = trx$tech_rep_group), "drop_incomplete")
  expect_identical(length(aligned$shared), 6L)
  for (name in aligned$shared) {
    expect_identical(unique(sub("_[12]$", "", prot$Run[prot$tech_rep_group == name])), trx$Run[trx$tech_rep_group == name])
  }
})

test_that("KIT-222 G5: the peptide replicate filter's grouping is unchanged on a merge-free design", {
  runs <- c(as.vector(t(outer(sprintf("A%d", 1:2), 1:3, paste, sep = "_"))),
    as.vector(t(outer(sprintf("B%d", 1:2), 1:2, paste, sep = "_"))))
  base <- head <- kit222Assign(kit222Assign(kit222Design(runs), runs[1:6], "GA"), runs[7:10], "GB")
  for (s in c("A1", "A2", "B1", "B2")) {
    selection <- runs[startsWith(runs, paste0(s, "_"))]
    base <- kit222BaseTechRep(base, selection, "lowest")
    head <- applyOmicsDesignTechRepAssignment(head, selection, "lowest")$designMatrix
  }
  expect_identical(head, base)
  peptides <- data.frame(Run = rep(runs, each = 2L), Protein.Ids = "P1", Stripped.Sequence = rep(c("AAK", "CCK"), 10L),
    stringsAsFactors = FALSE)
  peptides <- peptides[!(peptides$Stripped.Sequence == "CCK" & peptides$Run %in% c("A1_1", "A1_2", "B1_1")), ]
  filter <- function(design) removePeptidesWithOnlyOneReplicateHelper(peptides, design, Run, Run, replicates,
    Protein.Ids, Stripped.Sequence)
  expect_identical(filter(head), filter(base))
})
