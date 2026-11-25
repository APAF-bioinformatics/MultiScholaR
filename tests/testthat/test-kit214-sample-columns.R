# KIT-214 (VAL-173 (1); KIT-202 D5 (a), VAL-130): lipidomics and metabolomics duplicate resolution (C1) and the QC
# displays (C2 to C5) count only the design's samples: D = setdiff(intersect(table columns, design samples), id column),
# in the table's order (KIT-202 7.1). Tests E1 to E10, T6 (i) to (v), G1, G2, G3 and G7 of the ticket's section 10.
# G5 and G6 are other files (R2). Base values are G0's (docs/qa/workflow-execution/kit/kit-214/g0-base-values.R) and
# the ticket's hand computations (section 2.2, section 10).
# KIT-258 (VAL-254, VAL-284) declared edits (its section 6, items 1 and 2): `k214Earlier` also sets the rank option to
# KIT-214's rule ("mean"), so it means "a build before KIT-214" again; `k214OracleKept` states KIT-258's rule.

# ---- Helper section: the application's import path ------------------------------------------------------------------
# Copied from KIT-202's tests (branch kit-202-intensity-filter-samples at f9d7d910,
# tests/testthat/test-kit202-intensity-filter-samples.R, `k202RunInput`), itself from KIT-198's G0 script
# (docs/qa/workflow-execution/kit/kit-198/g0-numeric-columns.R at b0e913b8, `run_input`): the configuration's preview,
# the panel's column resolvers, `importOmicsTables`, `buildOmicsDesignSaveDataList` with the input's design, then the
# configuration's `design$build_object`. Returns the design object.
k214RunInput <- function(omic, vendor_format, assay1_file, assay1_name, assay2_file = NULL, assay2_name = NULL,
                         custom = NULL, sanitize = TRUE, design = NULL) {
  config <- omicGuiConfig(omic)
  import <- config$import
  preview <- config$import$preview(assay1_file, import = import)
  headers <- names(preview$assay_data)
  if (identical(vendor_format, "custom")) {
    id_col <- resolveOmicsImportColumn(headers, vendor_format, NULL, custom$id)
    annotation_col <- resolveOmicsImportColumn(headers, vendor_format, NULL, custom$annotation)
    sample_pattern <- custom$samples
    is_pattern <- custom$is_pattern
  } else {
    id_col <- preview$updates$id$selected
    annotation_col <- preview$updates$annotation$selected
    sample_pattern <- ""
    is_pattern <- preview$updates$is_pattern %||% ""
  }
  exclude <- if (identical(omic, "lipidomics")) OMICS_IMPORT_NORMALISED_COLUMNS
  sample_columns <- resolveOmicsImportSampleColumns(preview$assay_data, vendor_format, sample_pattern,
    preview$import_result, exclude)
  source_dir <- file.path(tempfile("kit214-import-"), "source")
  dir.create(source_dir, recursive = TRUE)
  request <- list(assay1_name = assay1_name, assay1_data = preview$assay_data, assay1_file = assay1_file,
    assay2_file = assay2_file, assay2_name = assay2_name, vendor_format = vendor_format,
    detected_format = preview$detected_format, id_col = id_col, annotation_col = annotation_col,
    sample_columns = sample_columns, sanitize_names = sanitize, is_pattern = is_pattern,
    experiment_paths = list(source_dir = source_dir))
  imported <- importOmicsTables(request, config)
  if (is.null(design)) {
    id_and_annotation <- c(imported$column_mapping[[import$id_mapping_key]], imported$column_mapping$annotation_col)
    samples <- setdiff(imported$sample_columns, id_and_annotation)
    design <- data.frame(sample = samples, group = sub("_.*$", "", samples), stringsAsFactors = FALSE)
  }
  data_cln <- buildOmicsDesignSaveDataList(imported$assays, design$sample, imported$column_mapping, import$id_mapping_key)
  design_matrix <- data.frame(Run = design$sample, group = design$group,
    tech_rep_group = paste(design$group, ave(design$sample, design$group, FUN = seq_along), sep = "_"),
    stringsAsFactors = FALSE)
  config$design$build_object(list(data_cln = data_cln, design_matrix = design_matrix, config_list = list()),
    imported$column_mapping)
}

k214E2e <- function(...) normalizePath(testthat::test_path("..", "testdata", "e2e", ...), mustWork = TRUE)

# A table written to a temporary file of its own and imported through MS-DIAL (panel defaults; every resolved sample
# assigned).
k214Msdial <- function(omic, table, file_name) {
  path <- file.path(tempfile("kit214-msdial-"), file_name)
  dir.create(dirname(path), recursive = TRUE)
  utils::write.csv(table, path, row.names = FALSE)
  k214RunInput(omic, "msdial", path, "LCMS_Pos")
}

k214Seed <- function() {
  utils::read.csv(k214E2e("metab_lc", "seed_lcms_pos.csv"), check.names = FALSE, stringsAsFactors = FALSE)
}

# The GUI-entry lanes with their manifest import settings and design (KIT-202's `k202Lane`: custom lanes id Feature.Name,
# no annotation, samples "^(WT|KO)_", internal standards "^ITSD|_d[0-9]+$"; sanitising off).
k214Lane <- function(lane_id) {
  lane <- guiEntryLane(lane_id)
  dir <- k214E2e(lane$fixture_dir)
  assays <- unlist(lane$assays)
  design <- utils::read.delim(file.path(dir, "design.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
  second <- !is.null(lane$assay2_file) && length(assays) >= 2L
  custom <- list(id = "Feature.Name", annotation = "", samples = "^(WT|KO)_", is_pattern = "^ITSD|_d[0-9]+$")
  k214RunInput(lane$omic_type, lane$import_tool, file.path(dir, lane$seed_file), assays[[1L]],
    assay2_file = if (second) file.path(dir, lane$assay2_file), assay2_name = if (second) assays[[2L]],
    custom = if (identical(lane$import_tool, "custom")) custom, sanitize = FALSE, design = design)
}

# ---- Helper section: the synthetic tables of section 2.2 -------------------------------------------------------------

k214Design <- function() {
  data.frame(Run = c("S1", "S2", "S3", "S4"), group = c("A", "A", "B", "B"), tech_rep_group = c("A_1", "A_2", "B_1", "B_2"),
    stringsAsFactors = FALSE)
}

k214Build <- function(omic, table, id_column, annotation_column) {
  if (identical(omic, "lipidomics")) {
    createLipidomicsAssayData(lipid_data = list(LCMS_Pos = table), design_matrix = k214Design(),
      lipid_id_column = id_column, annotation_id_column = annotation_column, sample_id = "Run", group_id = "group",
      technical_replicate_id = "tech_rep_group")
  } else {
    createMetaboliteAssayData(metabolite_data = list(LCMS_Pos = table), design_matrix = k214Design(),
      metabolite_id_column = id_column, annotation_id_column = annotation_column, sample_id = "Run", group_id = "group",
      technical_replicate_id = "tech_rep_group")
  }
}

# 2.2 example 1: a numeric (integer) id duplicated. Row A (50, NA, NA, NA), row B (60, 60, 60, 60) of id 1000; id 2000
# once. Over every numeric column A's mean is (1000 + 50) / 2 = 525 and B's (1000 + 240) / 5 = 248, so A is kept;
# over the design's samples A's is 50 and B's 60, so B is kept. The kept row is read from S1 (A 50, B 60).
k214E1Table <- function() {
  data.frame(feature_id = c(1000L, 1000L, 2000L), name = c("A", "B", "C"), S1 = c(50, 60, 100), S2 = c(NA, 60, 110),
    S3 = c(NA, 60, 120), S4 = c(NA, 60, 130), stringsAsFactors = FALSE)
}

# 2.2 example 2: a text id duplicated with a numeric annotation column. Row A (10 x 4; mz 100), row B (12 x 4; mz 10)
# of id M1. Every numeric column: A (40 + 100) / 5 = 28, B (48 + 10) / 5 = 11.6, A kept; the design's samples: B kept.
# The kept row is read from S1 (A 10, B 12).
k214E2Table <- function() {
  data.frame(feature_id = c("M1", "M1", "M2"), mz = c(100, 10, 50), S1 = c(10, 12, 30), S2 = c(10, 12, 31),
    S3 = c(10, 12, 32), S4 = c(10, 12, 33), stringsAsFactors = FALSE)
}

k214E1 <- function(omic) k214Build(omic, k214E1Table(), "feature_id", "name")
k214E2 <- function(omic) k214Build(omic, k214E2Table(), "feature_id", "mz")

# The kept rows of section 2.2: list(id, S1 value of the row each rule keeps).
k214Dup <- list(e1 = list(id = 1000L, head = 60, base = 50), e2 = list(id = "M1", head = 12, base = 10))

# ---- Helper section: fixtures, cached -----------------------------------------------------------------------------

k214Cache <- new.env(parent = emptyenv())
k214Object <- function(name) {
  if (is.null(k214Cache[[name]])) {
    k214Cache[[name]] <- suppressWarnings(suppressMessages(switch(name,
      # G0's import path: the metabolite seed (`Alignment ID` 1 to 8) through the metabolomics MS-DIAL import.
      metab_msdial = k214RunInput("metabolomics", "msdial", k214E2e("metab_lc", "seed_lcms_pos.csv"), "LCMS_Pos"),
      # The seed with `Alignment ID` 0 to 7 (the KIT harness corpora's shape).
      metab_msdial_id0 = {
        seed <- k214Seed()
        seed[["Alignment ID"]] <- seed[["Alignment ID"]] - 1L
        k214Msdial("metabolomics", seed, "seed_lcms_pos_id0.csv")
      },
      # module_ci_lipid_msdial_assay() with `Peak ID` = 1:3 through the lipid MS-DIAL import.
      lipid_msdial = {
        assay <- module_ci_lipid_msdial_assay()
        assay[["Peak ID"]] <- seq_len(nrow(assay))
        k214Msdial("lipidomics", assay, "msdial_lipid_pos.csv")
      },
      module_ci_lipid = module_ci_lipid_qc_object(),
      module_ci_metab = module_ci_metab_qc_object(),
      e1_lipid = k214E1("lipidomics"), e1_metab = k214E1("metabolomics"),
      e2_lipid = k214E2("lipidomics"), e2_metab = k214E2("metabolomics"),
      k214Lane(name))))
  }
  k214Cache[[name]]
}

# ---- Helper section: reading an object, running the step -----------------------------------------------------------

k214Lipid <- function(object) methods::is(object, "LipidomicsAssayData")
k214Omic <- function(object) if (k214Lipid(object)) "lipidomics" else "metabolomics"
k214Binding <- function(object) if (k214Lipid(object)) lipidomicsDeclarationBinding else metabolomicsDeclarationBinding
k214Assays <- function(object) methods::slot(object, k214Binding(object)$data_slot)
k214IdColumn <- function(object) methods::slot(object, k214Binding(object)$id_slot)
k214Samples <- function(object) as.character(object@design_matrix[[object@sample_id]])
k214Config <- function(object) omicGuiConfig(k214Omic(object))

# D of one assay (the definition, ticket 2.3), and D restricted to the columns that hold numbers (the display sites).
k214D <- function(table, object) setdiff(intersect(names(table), k214Samples(object)), k214IdColumn(object))
k214DNumeric <- function(table, object) {
  columns <- k214D(table, object)
  columns[vapply(columns, function(column) is.numeric(table[[column]]), logical(1))]
}

k214SetAssays <- function(object, assays) {
  methods::slot(object, k214Binding(object)$data_slot, check = FALSE) <- assays
  object
}

# The S1 value of the row kept for `id` in the first assay of a resolved object.
k214KeptS1 <- function(object, id) {
  table <- k214Assays(object)[[1L]]
  rows <- table[as.character(table[[k214IdColumn(object)]]) == as.character(id), , drop = FALSE]
  if (nrow(rows) != 1L) return(sprintf("%d rows of id %s", nrow(rows), id))
  rows$S1
}

# The declared step on its own owner in its own store root (a second run on one owner would return the first run's saved
# product: KIT-202 F2). Returns the resolved object, the variant the provenance names and the analysis record.
k214Step <- function(object, variant = c("whole", "partitioned")) {
  variant <- match.arg(variant)
  binding <- k214Binding(object)
  root <- withr::local_tempdir()
  fixture <- if (k214Lipid(object)) lipidDeclarationsOwner(root, object) else metabDeclarationsOwner(root, object)
  withr::defer(closeDIAStateAdapter(fixture$owner))
  step <- omicsDuplicateResolutionStep(binding)
  kit172Quiet(runWorkflowOmicsReferenceOperation(fixture$owner, workflowOperationEntry(step$operation)$declaration,
    c(step$parameters, list(variant = variant)), step$state_name, fixture$workflow$config_list, step$description))
  list(object = suppressWarnings(fixture$owner$getState(step$state_name)),
    variant = fixture$owner$.context$snapshot$step_provenance[[step$state_name]]$variant,
    record = buildWorkflowAnalysisRecord(fixture$owner), state = step$state_name)
}

k214Marker <- function(record, state) {
  step <- Filter(function(s) identical(s$state, state), record$steps)
  if (!length(step)) return(NULL)
  rows <- Filter(function(p) identical(p$name, "duplicate_resolution_samples"), step[[1L]]$parameters)
  if (!length(rows)) return(NULL)
  list(source = rows[[1L]]$source, value = decodeWorkflowRecordValue(rows[[1L]]$value))
}

k214Earlier <- list(MultiScholaR.duplicate_resolution_samples = "numeric_columns",
  MultiScholaR.duplicate_resolution_rank = "mean")

# The C2 value: the number the QC data summary shows after "Number of Samples:".
k214SummarySamples <- function(object) {
  html <- paste(as.character(omicQcDataSummaryUi(object, k214Config(object))), collapse = "")
  found <- regmatches(html, regexec("Number of Samples:</strong>\\s*</td>\\s*<td>([^<]*)</td>", html))[[1L]]
  if (length(found) < 2L) return(NA_character_)
  found[[2L]]
}

# An owner holding `object` as its head (one per call; closed when the calling test ends).
k214Owner <- function(object, env = parent.frame()) {
  omicsTestOwner(object, "raw", config = list(qc = list(mode = "kit214")), env = env)$state_manager
}

# C3 both ways: the whole object's statistics and the store's (the shell with its streamed missingness).
k214Stats <- function(object, env = parent.frame()) {
  config <- k214Config(object)
  shell <- readWorkflowOmicsShell(k214Binding(object), k214Owner(object, env), measures = c("features", "missing"))
  list(frame = omicQcAssayStatsFrame(object, config), shell = omicQcAssayStatsFrame(shell, config))
}

# ---- Oracles: base R only, from the rule (ticket 2.3, 10 S1), not from package code ------------------------------------

# Per id with more than one row (KIT-258's rule, VAL-254 and VAL-284; KIT-214's was the mean alone): a value is measured
# when neither NA nor 0; the rows with the most measured values over `columns`; among them the largest mean of a row's
# measured values (-Inf when none); among those the first row in input order is kept. Returns the sorted "id|.row"
# pairs of every kept row.
k214OracleKept <- function(table, id_column, columns) {
  values <- do.call(cbind, lapply(columns, function(column) as.numeric(table[[column]])))
  counts <- apply(values, 1L, function(x) sum(!is.na(x) & x != 0))
  means <- apply(values, 1L, function(x) {
    x <- x[!is.na(x) & x != 0]
    if (length(x)) mean(x) else -Inf
  })
  keys <- as.character(table[[id_column]])
  kept <- vapply(unique(keys), function(key) {
    rows <- which(keys == key)
    rows <- rows[counts[rows] == max(counts[rows])]
    rows[[which.max(means[rows])]]
  }, integer(1))
  sort(paste(keys[kept], table$.row[kept], sep = "|"))
}

k214KeptPairs <- function(table, id_column) sort(paste(as.character(table[[id_column]]), table$.row, sep = "|"))

# Percent of the cells of `columns` that are NA or zero, rounded to 1 decimal (the display's rule).
k214OracleMissing <- function(table, columns) {
  if (!length(columns)) return(NA_real_)
  cells <- unlist(lapply(columns, function(column) table[[column]]), use.names = FALSE)
  round(sum(is.na(cells) | cells == 0) / (nrow(table) * length(columns)) * 100, 1)
}

# Per feature row, the mean (missing removed) and the CV in percent (sd over mean when the mean is positive) over
# `columns`; sorted by feature id, then mean.
k214OracleRowStats <- function(table, id_column, columns) {
  values <- do.call(cbind, lapply(columns, function(column) as.numeric(table[[column]])))
  means <- apply(values, 1L, function(x) if (all(is.na(x))) NaN else mean(x[!is.na(x)]))
  sds <- apply(values, 1L, function(x) stats::sd(x, na.rm = TRUE))
  out <- data.frame(feature_id = as.character(table[[id_column]]), mean_intensity = means,
    cv_percent = ifelse(!is.na(means) & means > 0, sds / means * 100, NA_real_), stringsAsFactors = FALSE)
  out <- out[order(out$feature_id, out$mean_intensity), , drop = FALSE]
  rownames(out) <- NULL
  out
}

k214SortedTable <- function(table) {
  out <- data.frame(feature_id = as.character(table$feature_id), mean_intensity = table$mean_intensity,
    cv_percent = table$cv_percent, stringsAsFactors = FALSE)
  out <- out[order(out$feature_id, out$mean_intensity), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# The internal-standard metrics the QC tab shows, from the rule: the first of the id, annotation and name columns in
# which `pattern` matches; per matched value the pooled values of its rows over `columns`, missing removed: mean, and
# sd / mean x 100 (two values or more).
k214OracleItsd <- function(table, object, pattern, columns) {
  itsd <- k214Config(object)$qc$itsd
  search <- unique(c(k214IdColumn(object), object@annotation_id_column, itsd$name_columns))
  search <- search[search %in% names(table)]
  used <- Filter(function(column) any(grepl(pattern, as.character(table[[column]]))), search)
  if (!length(used)) return(NULL)
  keys <- as.character(table[[used[[1L]]]])
  ids <- unique(keys[grepl(pattern, keys)])
  values <- do.call(cbind, lapply(columns, function(column) as.numeric(table[[column]])))
  stats <- lapply(ids, function(id) {
    x <- as.vector(values[keys == id, , drop = FALSE])
    x <- x[!is.na(x)]
    data.frame(is_id = id, mean_intensity = mean(x), cv = if (length(x) >= 2L) stats::sd(x) / mean(x) * 100 else NA_real_,
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, stats)
  out[order(out$is_id), , drop = FALSE]
}

# The annotation column the selection table searches; NULL (the builder's default search) when the object names none
# that the table holds (the custom lanes have no annotation column). The sample columns are what G1 and G7 check.
k214Annotation <- function(object, table) {
  annotation <- object@annotation_id_column
  if (length(annotation) == 1L && !is.na(annotation) && annotation %in% names(table)) annotation
}

k214Itsd <- function(object, owner, pattern) {
  kit172Quiet(analyzeWorkflowOmicsItsd(k214Binding(object), owner, pattern, k214Config(object)$qc$itsd))
}

# ---- E1, E2: C1 on the declared step, both omics, both variants ------------------------------------------------------

for (example in c("e1", "e2")) {
  for (omic in c("lipid", "metab")) {
    for (variant in c("whole", "partitioned")) {
      local({
        name <- paste(example, omic, sep = "_")
        variant <- variant
        expected <- k214Dup[[example]]
        test_that(sprintf("KIT-214 %s (%s, %s): duplicate resolution keeps the row with the highest mean over the design's samples",
          toupper(example), name, variant), {
          skip_if_not(workflowOwnerCapable())
          run <- k214Step(k214Object(name), variant)
          expect_identical(run$variant, variant)
          expect_identical(k214KeptS1(run$object, expected$id), expected$head)
        })
      })
    }
  }
}

# ---- E3: the metabolite assay helper takes the design's samples -------------------------------------------------------

test_that("KIT-214 E3: resolveMetabDuplicateAssayData averages the design's samples it is given", {
  result <- tryCatch(kit172Quiet(resolveMetabDuplicateAssayData(assayList = list(LCMS_Pos = k214E1Table()),
    metaboliteIdCol = "feature_id", designSamples = c("S1", "S2", "S3", "S4"), logWarnFn = function(...) NULL)),
    error = function(e) NULL)
  expect_false(is.null(result))
  if (!is.null(result)) {
    table <- result$resolvedAssayList$LCMS_Pos
    expect_identical(table$S1[table$feature_id == 1000], 60)
    expect_identical(result$statsList$LCMS_Pos$removed, 1L)
  }
})

# ---- E4 to E7: the QC data summary and the per-assay statistics --------------------------------------------------------

test_that("KIT-214 E4: the QC summary's sample count is the design's samples (metabolite MS-DIAL seed)", {
  expect_identical(k214SummarySamples(k214Object("metab_msdial")), "4")
})

# Samples and Missingness per assay, the whole object's and the store's (base: the id column counted).
k214StatsExpected <- list(
  E5 = list(name = "metab_msdial", samples = 4L, missing = 0),
  E6 = list(name = "metab_msdial_id0", samples = 4L, missing = 0),
  E7 = list(name = "lipid_msdial", samples = 4L, missing = 16.7))

for (row in names(k214StatsExpected)) {
  local({
    row <- row
    expected <- k214StatsExpected[[row]]
    test_that(sprintf("KIT-214 %s (%s): the per-assay Samples and Missingness count the design's samples", row,
      expected$name), {
      skip_if_not(workflowOwnerCapable())
      stats <- k214Stats(k214Object(expected$name))
      for (path in c("frame", "shell")) {
        expect_identical(as.integer(stats[[path]]$Samples), expected$samples, info = path)
        expect_equal(as.numeric(stats[[path]]$Missingness), expected$missing, info = path)
      }
    })
  })
}

# ---- E8: the internal-standard metrics of the QC tab ------------------------------------------------------------------

test_that("KIT-214 E8: the ITSD metrics, the CV > 30 % count and the long table count the design's samples", {
  skip_if_not(workflowOwnerCapable())
  object <- k214Object("metab_msdial")
  analysis <- k214Itsd(object, k214Owner(object), "ITSD")
  samples <- c("ctrl_1", "ctrl_2", "treat_1", "treat_2")
  seed <- k214Seed()
  for (standard in c("ITSD_Caffeine_d9", "ITSD_Leucine_d10")) {
    values <- unlist(seed[seed$Name == standard, samples], use.names = FALSE)
    metrics <- as.data.frame(analysis$metrics)
    found <- metrics[metrics$is_id == standard, , drop = FALSE]
    expect_identical(nrow(found), 1L, info = standard)
    expect_equal(found$mean_intensity, mean(values), info = standard)
    expect_equal(found$cv, stats::sd(values) / mean(values) * 100, info = standard)
  }
  expect_match(analysis$result_text, "IS with CV > 30%: 0", fixed = TRUE)
  expect_true(all(unique(analysis$long_data$Sample) %in% samples))
})

# ---- E9: the normalisation's ITSD selection table ---------------------------------------------------------------------

test_that("KIT-214 E9: the ITSD selection table averages the sample columns it is given (numeric annotation column)", {
  seed <- k214Seed()
  samples <- c("ctrl_1", "ctrl_2", "treat_1", "treat_2")
  table <- seed[, c("Name", "Average Mz", samples)]
  build <- omicGuiConfig("metabolomics")$norm$itsd_table
  result <- tryCatch(build(table, "Name", "Average Mz", sample_columns = samples), error = function(e) NULL)
  expect_false(is.null(result))
  if (!is.null(result)) {
    expect_equal(k214SortedTable(result)[c("feature_id", "mean_intensity")],
      k214OracleRowStats(table, "Name", samples)[c("feature_id", "mean_intensity")])
  }
})

# ---- E10: the refusal under the current rule only --------------------------------------------------------------------

# A module-CI object (duplicated text ids) whose first design sample is text in its first assay; validity bypassed.
k214TextSample <- function(name) {
  object <- k214Object(name)
  assays <- k214Assays(object)
  sample <- k214Samples(object)[[1L]]
  assays[[1L]][[sample]] <- as.character(assays[[1L]][[sample]])
  k214SetAssays(object, assays)
}

for (name in c("module_ci_lipid", "module_ci_metab")) {
  local({
    name <- name
    test_that(sprintf("KIT-214 E10 (%s): a text design sample is refused before any row is chosen", name), {
      object <- k214TextSample(name)
      binding <- k214Binding(object)
      expect_error(kit172Quiet(binding$resolve_duplicates(object)), class = "omics_duplicate_resolution_samples_invalid")
      withr::with_options(k214Earlier, expect_no_error(kit172Quiet(binding$resolve_duplicates(object))))
    })
  })
}

# ---- T6: the replay marker --------------------------------------------------------------------------------------------

test_that("KIT-214 T6 (i): the replay rule of duplicate resolution follows the step's marker", {
  markers <- get0("workflowReplayRuleMarkers", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  expect_true(!is.null(markers[["duplicate_resolution_samples"]]))
  if (!is.null(markers[["duplicate_resolution_samples"]])) {
    expect_identical(markers[["duplicate_resolution_samples"]], list(option = "MultiScholaR.duplicate_resolution_samples",
      current = "design", earlier = "numeric_columns"))
    entry <- function(source) list(name = "duplicate_resolution_samples", source = source,
      value = encodeWorkflowRecordValue("design"), provenance = "fixed")
    for (operation in c("lipid_duplicate_resolution", "metab_duplicate_resolution")) {
      step <- function(parameters) list(state = "s", operation = operation, parameters = parameters)
      expect_identical(workflowReplayRule(step(list(entry("fixed"))), "duplicate_resolution_samples"), "design",
        info = operation)
      expect_identical(workflowReplayRule(step(list()), "duplicate_resolution_samples"), "numeric_columns",
        info = operation)
      expect_identical(workflowReplayRule(step(list(entry("registry"))), "duplicate_resolution_samples"),
        "numeric_columns", info = operation)
    }
  }
})

for (example in c("e1", "e2")) {
  for (omic in c("lipid", "metab")) {
    local({
      name <- paste(example, omic, sep = "_")
      expected <- k214Dup[[example]]
      test_that(sprintf("KIT-214 T6 (ii) (%s): under the earlier rule the step keeps the base's row", name), {
        skip_if_not(workflowOwnerCapable())
        withr::local_options(k214Earlier)
        for (variant in c("whole", "partitioned")) {
          # New owners in new store roots (k214Step), never one holding a product of the other rule.
          run <- k214Step(k214Object(name), variant)
          expect_identical(run$variant, variant)
          expect_identical(k214KeptS1(run$object, expected$id), expected$base, info = variant)
        }
      })
    })
  }
}

for (omic in c("lipid", "metab")) {
  local({
    name <- paste("e1", omic, sep = "_")
    test_that(sprintf("KIT-214 T6 (iii) (%s): a duplicate step made now records the marker", name), {
      skip_if_not(workflowOwnerCapable())
      run <- k214Step(k214Object(name))
      expect_identical(k214Marker(run$record, run$state), list(source = "fixed", value = "design"))
    })
  })
}

# T6 (iv): as KIT-202's T6 (iv) (branch kit-202-intensity-filter-samples at f9d7d910, `k202Replay`), itself after
# KIT-190's T7 (iv) and KIT-187's T8 (iv): the owner's record written, the owner closed, a copy cut to the design state
# and the record replayed on it; the replay compared with the origin.
k214Paths <- function(root, lipid) {
  qc <- if (lipid) "lipid_qc_dir" else "metabolite_qc_dir"
  stats::setNames(list(file.path(root, "source"), file.path(root, "qc"), file.path(root, "da"), file.path(root, "pub")),
    c("source_dir", qc, "da_output_dir", "publication_graphs_dir"))
}

# One owner per rule in its own store root; `earlier = TRUE` computes the step under the earlier rule, so its provenance
# has no marker, as a build before this change wrote it. The step runs as the module runs it (its default variant).
k214ReplayOwner <- function(root, object, earlier) {
  fixture <- if (k214Lipid(object)) lipidDeclarationsOwner(root, object) else metabDeclarationsOwner(root, object)
  run <- function() kit172Quiet(runWorkflowOmicsReferenceDuplicateResolution(k214Binding(object), fixture$owner,
    fixture$workflow$config_list))
  if (earlier) withr::with_options(k214Earlier, run()) else run()
  fixture
}

k214Replay <- function(fixture, root, object, id) {
  prefix <- k214Binding(object)$prefix
  state <- paste0(prefix, "_duplicates_resolved")
  kept <- k214KeptS1(suppressWarnings(fixture$owner$getState(state)), id)
  json <- file.path(root, "analysis_record.json")
  writeWorkflowAnalysisRecord(fixture$owner, json)
  record <- jsonlite::read_json(json, simplifyVector = FALSE)
  closeDIAStateAdapter(fixture$owner)
  source_dir <- file.path(root, "source")
  store_root <- dirname(dirname(list.files(source_dir, "^store[.]json$", recursive = TRUE, all.files = TRUE,
    full.names = TRUE)[[1L]]))
  copy <- tempfile("k214-copy-")
  dir.create(copy)
  stopifnot(system2("cp", c("-a", file.path(source_dir, "."), copy)) == 0L)
  harness <- new.env(parent = globalenv())
  sys.source(file.path(testthat::test_path("..", ".."), "tools", "profiling", "workflow_kit.R"), harness)
  api <- function(name) get(name, asNamespace("MultiScholaR"))
  copy_root <- file.path(copy, substring(store_root, nchar(source_dir) + 2L))
  prepared <- harness$kitPrepareReplayProject(copy_root, c("initial", paste0(prefix, "_raw_data_s4")), api)
  inputs <- c(stats::setNames(file.path(source_dir, "import_sources", "LCMS_Pos_export.txt"), paste0(prefix, "_assay_1")),
    design_file = file.path(source_dir, "design_matrix.tab"))
  replay <- kit172Quiet(harness$kitReplayJourney(record, prepared, store_root, inputs, k214Paths(root, k214Lipid(object)),
    tempfile("k214-scratch-"), api))
  identity <- readSavedStoreIdentity(file.path(prepared, ".workflow-entities", "store.json"))
  replayed <- newDIAStateAdapter(newWorkflowEntityStore(prepared, identity$project_id, identity$workflow_id), reopen = TRUE)
  replayed_kept <- k214KeptS1(suppressWarnings(replayed$getState(state)), id)
  replayed_record <- buildWorkflowAnalysisRecord(replayed)
  closeDIAStateAdapter(replayed)
  list(record = record, replayed_record = replayed_record, mismatches = as.character(replay$mismatches), kept = kept,
    replayed_kept = replayed_kept, state = state)
}

for (omic in c("lipid", "metab")) {
  local({
    name <- paste("e1", omic, sep = "_")
    test_that(sprintf("KIT-214 T6 (iv) (%s): records made before and after the change replay, each with its own rule", name), {
      skip_if_not(workflowOwnerCapable())
      object <- k214Object(name)
      id <- k214Dup$e1$id
      root_a <- withr::local_tempdir()
      root_b <- withr::local_tempdir()
      a <- k214Replay(k214ReplayOwner(root_a, object, earlier = TRUE), root_a, object, id)
      b <- k214Replay(k214ReplayOwner(root_b, object, earlier = FALSE), root_b, object, id)
      # Non-vacuity: the two rules keep different rows on this table.
      expect_false(identical(a$kept, b$kept))
      expect_identical(a$mismatches, character())
      expect_identical(b$mismatches, character())
      expect_identical(a$replayed_kept, a$kept)
      expect_identical(b$replayed_kept, b$kept)
      expect_null(k214Marker(a$record, a$state))
      expect_null(k214Marker(a$replayed_record, a$state))
      expect_identical(k214Marker(b$record, b$state), list(source = "fixed", value = "design"))
      expect_identical(k214Marker(b$replayed_record, b$state), list(source = "fixed", value = "design"))
    })
  })
}

# T6 (v): a managed step, as KIT-202's T6 (v): the worker's entry point (omicsWorkerRun, as
# tests/testthat/test-kit129-run-path.R:189), completed with completeWorkflowAction.
test_that("KIT-214 T6 (v): a managed duplicate step records the marker", {
  skip_if_not(workflowOwnerCapable())
  testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64),
    assertWorkflowController = function(...) invisible(TRUE), .package = "MultiScholaR")
  for (lipid in c(TRUE, FALSE)) {
    binding <- if (lipid) lipidomicsDeclarationBinding else metabolomicsDeclarationBinding
    root <- withr::local_tempdir()
    fixture <- if (lipid) lipidDeclarationsOwner(root) else metabDeclarationsOwner(root)
    withr::defer(closeDIAStateAdapter(fixture$owner))
    step <- omicsDuplicateResolutionStep(binding)
    worker <- omicsWorkerRun(fixture, step$operation, step$parameters)
    expect_identical(worker$status, "complete", info = paste(binding$prefix, worker$error))
    kit172Quiet(completeWorkflowAction(list(manager = fixture$owner, accepting = TRUE), step$operation, worker$result,
      context = list(state_name = step$state_name, config = list(), description = step$description),
      values = step$parameters))
    expect_identical(k214Marker(buildWorkflowAnalysisRecord(fixture$owner), step$state_name),
      list(source = "fixed", value = "design"), info = binding$prefix)
  }
})

# ---- G1: unaffected objects (passes on base and head) ----------------------------------------------------------------

k214Unaffected <- c("lipid_canonical", "metab_lc", "metab_gc", "metab_combined", "module_ci_lipid", "module_ci_metab")

# The object with a text column `.row` (the row's position) in each assay: neither numeric nor a design sample.
k214WithRow <- function(object) {
  assays <- lapply(k214Assays(object), function(table) {
    table$.row <- as.character(seq_len(nrow(table)))
    table
  })
  k214SetAssays(object, assays)
}

for (name in k214Unaffected) {
  local({
    name <- name
    test_that(sprintf("KIT-214 G1 (%s): an object whose numeric columns are its samples gives the oracle's values", name), {
      object <- k214Object(name)
      id <- k214IdColumn(object)
      assays <- k214Assays(object)
      config <- k214Config(object)
      # C1: the kept rows, against the oracle over D.
      marked <- k214WithRow(object)
      resolved <- kit172Quiet(k214Binding(object)$resolve_duplicates(marked)$object)
      for (assay in names(assays)) {
        table <- k214Assays(marked)[[assay]]
        expect_identical(k214KeptPairs(k214Assays(resolved)[[assay]], id), k214OracleKept(table, id, k214D(table, object)),
          info = paste("C1", assay))
      }
      # C2 and C3 (the whole object's statistics).
      expect_identical(k214SummarySamples(object), as.character(length(k214D(assays[[1L]], object))))
      stats <- omicQcAssayStatsFrame(object, config)
      for (assay in names(assays)) {
        table <- assays[[assay]]
        row <- stats[stats$Assay == assay, , drop = FALSE]
        expect_identical(as.integer(row$Samples), length(k214D(table, object)), info = paste("C3", assay))
        expect_equal(as.numeric(row$Missingness), k214OracleMissing(table, k214DNumeric(table, object)),
          info = paste("C3", assay))
      }
      # C5: the selection table as the panel's builder makes it.
      for (assay in names(assays)) {
        table <- assays[[assay]]
        built <- config$norm$itsd_table(table, id, k214Annotation(object, table))
        expect_equal(k214SortedTable(built), k214OracleRowStats(table, id, k214DNumeric(table, object)),
          info = paste("C5", assay))
      }
      # C4 on the module-CI objects (their pattern "^IS_" finds standards in every assay).
      if (startsWith(name, "module_ci")) {
        skip_if_not(workflowOwnerCapable())
        analysis <- k214Itsd(object, k214Owner(object), "^IS_")
        for (assay in names(assays)) {
          table <- assays[[assay]]
          oracle <- k214OracleItsd(table, object, "^IS_", k214DNumeric(table, object))
          found <- as.data.frame(analysis$metrics_by_assay[[assay]])
          found <- found[order(found$is_id), c("is_id", "mean_intensity", "cv"), drop = FALSE]
          rownames(found) <- NULL
          rownames(oracle) <- NULL
          expect_equal(found, oracle, info = paste("C4", assay))
        }
      }
    })
  })
}

# ---- G2, G3: C1's other contracts (pass on base and head) ------------------------------------------------------------

for (omic in c("lipidomics", "metabolomics")) {
  local({
    omic <- omic
    test_that(sprintf("KIT-214 G2 (%s): a table with a numeric id and no duplicated id is returned unchanged", omic), {
      object <- k214Build(omic, k214E1Table()[c(1L, 3L), , drop = FALSE], "feature_id", "name")
      resolved <- kit172Quiet(k214Binding(object)$resolve_duplicates(object)$object)
      table <- k214Assays(object)[[1L]]
      out <- k214Assays(resolved)[[1L]]
      cat(sprintf("KIT-214 G2 %s: input %s output %s\n", omic, digest::digest(table, algo = "sha256"),
        digest::digest(out, algo = "sha256")))
      expect_identical(out, table)
    })

    test_that(sprintf("KIT-214 G3 (%s): C1 keeps the columns, their order and every numeric column as double", omic), {
      object <- k214E1(omic)
      table <- k214Assays(object)[[1L]]
      out <- k214Assays(kit172Quiet(k214Binding(object)$resolve_duplicates(object)$object))[[1L]]
      expected <- vapply(table, function(column) if (is.numeric(column)) "numeric" else class(column)[[1L]], character(1))
      expect_identical(names(out), names(table))
      expect_identical(vapply(out, function(column) class(column)[[1L]], character(1)), expected)
      expect_identical(nrow(out), 2L)
    })
  })
}

# ---- G7: a text design sample at the display sites -------------------------------------------------------------------
# (a) The per-assay Samples counts it (it is a design sample: E-8, RC1-02). On base Samples counts the numeric columns,
# so (a) fails on base at one fewer (step-0 correction of the ticket's G7 base column). (b) It is left out of the
# missingness, the ITSD metrics and the selection table; an error fails the test (passes on base and head).

for (name in c("module_ci_lipid", "module_ci_metab")) {
  local({
    name <- name
    test_that(sprintf("KIT-214 G7 (a) (%s): the per-assay Samples counts a text design sample", name), {
      skip_if_not(workflowOwnerCapable())
      object <- k214TextSample(name)
      stats <- k214Stats(object)
      first <- names(k214Assays(object))[[1L]]
      for (path in c("frame", "shell")) {
        row <- stats[[path]][stats[[path]]$Assay == first, , drop = FALSE]
        expect_identical(as.integer(row$Samples), length(k214Samples(object)), info = path)
      }
    })

    test_that(sprintf("KIT-214 G7 (b) (%s): a text design sample is left out of the statistics, without error", name), {
      skip_if_not(workflowOwnerCapable())
      object <- k214TextSample(name)
      id <- k214IdColumn(object)
      assays <- k214Assays(object)
      config <- k214Config(object)
      stats <- k214Stats(object)
      for (assay in names(assays)) {
        table <- assays[[assay]]
        columns <- k214DNumeric(table, object)
        for (path in c("frame", "shell")) {
          row <- stats[[path]][stats[[path]]$Assay == assay, , drop = FALSE]
          expect_equal(as.numeric(row$Missingness), k214OracleMissing(table, columns), info = paste(path, assay))
        }
        built <- config$norm$itsd_table(table, id, k214Annotation(object, table))
        expect_equal(k214SortedTable(built), k214OracleRowStats(table, id, columns), info = paste("C5", assay))
      }
      analysis <- k214Itsd(object, k214Owner(object), "^IS_")
      for (assay in names(assays)) {
        table <- assays[[assay]]
        oracle <- k214OracleItsd(table, object, "^IS_", k214DNumeric(table, object))
        found <- as.data.frame(analysis$metrics_by_assay[[assay]])
        found <- found[order(found$is_id), c("is_id", "mean_intensity", "cv"), drop = FALSE]
        rownames(found) <- NULL
        rownames(oracle) <- NULL
        expect_equal(found, oracle, info = paste("C4", assay))
      }
    })
  })
}
