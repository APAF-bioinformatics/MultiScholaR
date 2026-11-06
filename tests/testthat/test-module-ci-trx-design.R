# KIT-096 MCI-027: the transcriptomics design, the lipid design builder and its commit under the transcriptomics profile.
# Runs are matched to count columns by name only and keep the counts file's order; KIT-093's design validation runs at
# save, with at least two runs per group; the level order entered is recorded; technical replicates are committed.

trxDesignStateObject <- function(owner) withWorkflowStateScope(owner, function() readTrxWorkflowState(owner))

#' The object without the attributes that name its store (run-specific).
trxDesignObjectValue <- function(object) {
  attr(object, "dia_component_source") <- NULL
  attr(object, "dia_execution_fields") <- NULL
  object
}

trxDesignRefusal <- function(commit, class) {
  shown <- vapply(Filter(function(note) identical(note$type, "error"), commit$notes), `[[`, "", "message")
  testthat::expect_length(shown, 1L)
  testthat::expect_match(shown, class, fixed = TRUE)
}

test_that("MCI-027.1 the builder keeps the counts file's run order and names; the commit is the design state", {
  skip_if_not(workflowOwnerCapable())
  truth <- module_ci_trx_truth()
  fixture <- module_ci_trx_fixture_owner()
  module_ci_trx_import(fixture, file.path(module_ci_trx_corpus(), "counts.tsv"))
  design <- module_ci_trx_corpus_design()
  result <- module_ci_trx_build_design(fixture, design, c("RPMI", "Sera"))
  # The fresh design is the counts file's runs, in its order (not mixedsort), unsanitised.
  expect_identical(attr(result, "initial_runs"), truth$column_order)
  expect_identical(result$design_matrix$Run, truth$column_order)
  expect_identical(result$group_levels, c("RPMI", "Sera"))
  commit <- module_ci_trx_commit(fixture, result)
  owner <- fixture$workflow$state_manager
  expect_identical(owner$current_state, "trx_raw_data_s4")
  expect_identical(owner$.context$store$identity$workflow_id, "TRANSCRIPTOMICS")
  expect_identical(fixture$workflow$tab_status$design_matrix, "complete")
  expect_true(isTRUE(commit$qc))
  messages <- vapply(commit$notes, `[[`, "", "message")
  expect_true("Design matrix and contrasts saved successfully!" %in% messages)
  object <- trxDesignStateObject(owner)
  expect_s4_class(object, "TranscriptomicsAssayData")
  expect_identical(object@args$group_levels, c("RPMI", "Sera"))
  expect_identical(object@gene_id_column, "GeneID")
  expect_identical(object@annotation_id_column, "GeneName")
  expect_identical(object@technical_replicate_id, "tech_rep_group")
  # Names and order reach the kernels unchanged: the class's inputs are the corpus's runs, by name.
  inputs <- trxAssayInputs(object)
  expect_identical(colnames(inputs$counts), truth$column_order)
  expect_identical(as.character(inputs$design$Run), truth$column_order)
  expect_identical(levels(inputs$group), c("RPMI", "Sera"))
  expect_identical(as.character(inputs$group), design$group[match(truth$column_order, design$Run)])
  # The technical replicates partition the runs as the corpus's biological samples do.
  partition <- function(runs, sets) {
    parts <- lapply(split(as.character(runs), sets), sort)
    unname(parts[order(vapply(parts, `[`, "", 1L))])
  }
  committed <- inputs$design
  expect_identical(partition(committed$Run, committed$tech_rep_group), partition(design$Run, design$tech_rep_group))
  expect_true(validateWorkflowAnalysisRecord(buildWorkflowAnalysisRecord(owner))$valid)
  for (file in c("design_matrix.tab", "contrast_strings.tab", "data_cln_counts.tab")) {
    expect_true(file.exists(file.path(fixture$paths$source_dir, file)), info = file)
  }
})

test_that("MCI-027.1 the level order is the order the user entered", {
  skip_if_not(workflowOwnerCapable())
  fixture <- module_ci_trx_fixture_owner()
  module_ci_trx_import(fixture, file.path(module_ci_trx_corpus(), "counts.tsv"))
  result <- module_ci_trx_build_design(fixture, module_ci_trx_corpus_design(), c("Sera", "RPMI"), technical = FALSE)
  expect_identical(result$group_levels, c("Sera", "RPMI"))
  module_ci_trx_commit(fixture, result)
  expect_identical(trxDesignStateObject(fixture$workflow$state_manager)@args$group_levels, c("Sera", "RPMI"))
})

test_that("MCI-027.2 the design commit's state and provenance are those of KIT-095's journey setup", {
  skip_if_not(workflowOwnerCapable())
  fixture <- module_ci_trx_fixture_owner()
  module_ci_trx_designed(fixture)
  gui <- fixture$workflow$state_manager
  committed <- trxDesignStateObject(gui)@design_matrix
  # The journey's setup (helper-workflow-trx-declarations.R) on the same corpus and the design the GUI committed.
  journey <- trxDeclarationsOwner(withr::local_tempdir(), module_ci_trx_corpus(), design = committed,
    technical_replicate_id = "tech_rep_group")
  withr::defer(closeDIAStateAdapter(journey$owner))
  owner <- journey$owner
  expect_identical(gui$getHistory(), owner$getHistory())
  expect_identical(gui$current_state, owner$current_state)
  expect_identical(gui$.context$store$identity$workflow_id, owner$.context$store$identity$workflow_id)
  expect_identical(trxDesignObjectValue(trxDesignStateObject(gui)), trxDesignObjectValue(trxDesignStateObject(owner)))
  metadata <- function(manager) {
    meta <- readWorkflowOmicsMetadata(transcriptomicsDeclarationBinding, manager)
    meta[c("class", "slots", "assay_names", "assay_rows", "design_matrix")]
  }
  expect_identical(metadata(gui), metadata(owner))
  provenance <- function(manager) manager$.context$snapshot$step_provenance$trx_raw_data_s4[c("operation", "parameter_names", "values")]
  expect_identical(provenance(gui), provenance(owner))
  inputs <- function(manager) lapply(manager$.context$snapshot$file_dependencies, `[`, c("role", "bytes", "sha256"))
  expect_identical(inputs(gui), inputs(owner))
  expect_identical(names(inputs(gui)), c("trx_counts_1", "design_file"))
})

test_that("MCI-027.3 mismatched, unassigned, duplicate, unknown-level, empty and single-run groups are refused by class", {
  skip_if_not(workflowOwnerCapable())
  design <- module_ci_trx_corpus_design()
  fresh <- function() {
    fixture <- module_ci_trx_fixture_owner(env = parent.frame())
    module_ci_trx_import(fixture, file.path(module_ci_trx_corpus(), "counts.tsv"))
    fixture
  }
  refused <- function(fixture, result, class) {
    commit <- module_ci_trx_commit(fixture, result)
    trxDesignRefusal(commit, class)
    expect_identical(fixture$workflow$state_manager$current_state %||% "initial", "initial", info = class)
    expect_false(identical(fixture$workflow$tab_status$design_matrix, "complete"), info = class)
  }
  fixture <- fresh()
  good <- module_ci_trx_build_design(fixture, design, c("RPMI", "Sera"), technical = FALSE)
  # A design whose run is not a count column: refused by name, with no positional mapping.
  mismatched <- good
  extra <- good$design_matrix[1L, ]
  extra$Run <- "not_a_count_column"
  mismatched$design_matrix <- rbind(good$design_matrix, extra)
  refused(fixture, mismatched, "trx_design_runs_missing")
  duplicated <- good
  duplicated$design_matrix <- rbind(good$design_matrix, good$design_matrix[1L, ])
  refused(fixture, duplicated, "trx_design_runs_duplicated")
  unknown <- good
  unknown$group_levels <- c("RPMI", "Other")
  refused(fixture, unknown, "trx_design_level_unknown")
  empty <- good
  empty$group_levels <- c("RPMI", "Sera", "Unused")
  refused(fixture, empty, "trx_design_group_empty")
  # An unassigned run (not removed by the user) is a count column outside the design.
  partial <- design[design$Run != design$Run[[1L]], ]
  fixture <- fresh()
  refused(fixture, module_ci_trx_build_design(fixture, partial, c("RPMI", "Sera"), technical = FALSE), "trx_design_runs_extra")
  # A group of one run.
  single <- design
  single$group[single$Run == single$Run[single$group == "RPMI"][[1L]]] <- "Solo"
  fixture <- fresh()
  refused(fixture, module_ci_trx_build_design(fixture, single, c("RPMI", "Sera", "Solo"), technical = FALSE), "trx_design_group_singleton")
})

test_that("MCI-027.3 removed runs leave the design and the counts; three groups are accepted with their order", {
  skip_if_not(workflowOwnerCapable())
  design <- module_ci_trx_corpus_design()
  fixture <- module_ci_trx_fixture_owner()
  module_ci_trx_import(fixture, file.path(module_ci_trx_corpus(), "counts.tsv"))
  removed <- design$Run[design$tech_rep_group == design$tech_rep_group[[1L]]]
  kept <- design[!design$Run %in% removed, ]
  module_ci_trx_commit(fixture, module_ci_trx_build_design(fixture, kept, c("RPMI", "Sera"), technical = FALSE, remove = removed))
  inputs <- trxAssayInputs(trxDesignStateObject(fixture$workflow$state_manager))
  expect_false(any(removed %in% colnames(inputs$counts)))
  expect_identical(sort(colnames(inputs$counts)), sort(kept$Run))
  three <- design
  three$group[three$group == "RPMI" & grepl("_0[23]_", three$Run)] <- "RPMI_late"
  fixture <- module_ci_trx_fixture_owner()
  module_ci_trx_import(fixture, file.path(module_ci_trx_corpus(), "counts.tsv"))
  module_ci_trx_commit(fixture, module_ci_trx_build_design(fixture, three, c("Sera", "RPMI", "RPMI_late"), technical = FALSE))
  object <- trxDesignStateObject(fixture$workflow$state_manager)
  expect_identical(object@args$group_levels, c("Sera", "RPMI", "RPMI_late"))
})

test_that("MCI-027.4 the design tab under the transcriptomics profile: the shared builder ids, matching by name stated", {
  html <- htmltools::renderTags(mod_lipid_design_ui("design", profile = trxOmicProfile()))$html
  expect_match(html, "by sample name", fixed = TRUE)
  expect_false(grepl("design-show_import_modal", html, fixed = TRUE))
  expect_false(grepl("Multi-Assay Workflow", html, fixed = TRUE))
  for (id in c("main_tabs", "new_factor", "add_factor", "selected_runs", "factor1_select", "assign_metadata", "save_results")) {
    expect_match(html, sprintf("id=\"design-builder-%s\"", id), fixed = TRUE)
  }
  lipid <- htmltools::renderTags(mod_lipid_design_ui("design"))$html
  expect_match(lipid, "design-show_import_modal", fixed = TRUE)
  expect_match(lipid, "Multi-Assay Workflow", fixed = TRUE)
})
