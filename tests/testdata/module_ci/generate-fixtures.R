#!/usr/bin/env Rscript

fixture_classes <- c(
  "happy_path",
  "missingness",
  "duplicates",
  "invalid_design",
  "multi_factor",
  "bad_names",
  "nonfinite",
  "small_n",
  "large_enough_for_plots",
  "multi_assay_mixed"
)

import_formats <- list(
  proteomics = c("dia", "dda", "tmt", "dia_limpa"),
  metabolomics = c("lc", "gc"),
  lipidomics = c("lc", "gc")
)

fixture_class_contracts <- list(
  list(
    fixture_class = "happy_path",
    case_category = "normal",
    ci_lane = "push",
    push_safe = TRUE,
    expected_warnings = character(),
    expected_validation_errors = character()
  ),
  list(
    fixture_class = "missingness",
    case_category = "malformed",
    ci_lane = "push",
    push_safe = TRUE,
    expected_warnings = c("missing_values_detected", "row_missingness_filter_candidate"),
    expected_validation_errors = character()
  ),
  list(
    fixture_class = "duplicates",
    case_category = "malformed",
    ci_lane = "push",
    push_safe = TRUE,
    expected_warnings = c("duplicate_feature_ids_detected"),
    expected_validation_errors = c("duplicate_feature_ids")
  ),
  list(
    fixture_class = "invalid_design",
    case_category = "malformed",
    ci_lane = "push",
    push_safe = TRUE,
    expected_warnings = character(),
    expected_validation_errors = c("design_sample_mismatch", "duplicate_design_samples")
  ),
  list(
    fixture_class = "multi_factor",
    case_category = "normal",
    ci_lane = "push",
    push_safe = TRUE,
    expected_warnings = character(),
    expected_validation_errors = character()
  ),
  list(
    fixture_class = "bad_names",
    case_category = "malformed",
    ci_lane = "push",
    push_safe = TRUE,
    expected_warnings = c("unsafe_sample_names_detected"),
    expected_validation_errors = c("unsafe_sample_names")
  ),
  list(
    fixture_class = "nonfinite",
    case_category = "malformed",
    ci_lane = "push",
    push_safe = TRUE,
    expected_warnings = c("nonfinite_intensity_values_detected"),
    expected_validation_errors = c("nonfinite_intensity_values")
  ),
  list(
    fixture_class = "small_n",
    case_category = "boundary",
    ci_lane = "push",
    push_safe = TRUE,
    expected_warnings = c("single_replicate_guard"),
    expected_validation_errors = c("insufficient_replicates_for_model")
  ),
  list(
    fixture_class = "large_enough_for_plots",
    case_category = "boundary",
    ci_lane = "nightly",
    push_safe = FALSE,
    expected_warnings = character(),
    expected_validation_errors = character()
  ),
  list(
    fixture_class = "multi_assay_mixed",
    case_category = "multi_assay",
    ci_lane = "push",
    push_safe = TRUE,
    expected_warnings = c("multi_assay_input_detected"),
    expected_validation_errors = character()
  )
)

feature_column_for <- function(omic, import_format) {
  if (identical(omic, "proteomics")) {
    if (identical(import_format, "tmt")) {
      return("Annotated.Sequence")
    }
    if (identical(import_format, "dda")) {
      return("Protein")
    }
    return("Protein.Group")
  }
  if (identical(omic, "lipidomics")) {
    return("LipidName")
  }
  "Feature.Name"
}

assay_names_for <- function(omic, import_format, fixture_class) {
  if (identical(fixture_class, "multi_assay_mixed")) {
    if (identical(omic, "proteomics")) {
      return(c("DIA", "TMT"))
    }
    if (identical(omic, "metabolomics")) {
      return(c("LCMS_Pos", "GCMS"))
    }
    return(c("LCMS_Pos", "GCMS"))
  }
  if (identical(omic, "proteomics")) {
    switch(
      import_format,
      dia = "DIA",
      dda = "LFQ",
      tmt = "TMT",
      dia_limpa = c("DIA", "LIMPA_IMPUTED")
    )
  } else if (identical(import_format, "gc")) {
    "GCMS"
  } else {
    "LCMS_Pos"
  }
}

case_contract <- function(fixture_class) {
  fixture_class_contracts[[match(fixture_class, fixture_classes)]]
}

sample_columns_for <- function(fixture_class) {
  switch(
    fixture_class,
    bad_names = c("WT 1", "WT-2", "KO/1", "KO.2"),
    small_n = c("WT_1", "KO_1"),
    large_enough_for_plots = c(sprintf("WT_%d", 1:6), sprintf("KO_%d", 1:6)),
    c("WT_1", "WT_2", "KO_1", "KO_2")
  )
}

feature_count_for <- function(fixture_class) {
  switch(
    fixture_class,
    small_n = 4L,
    large_enough_for_plots = 16L,
    6L
  )
}

feature_prefix_for <- function(omic) {
  switch(
    omic,
    proteomics = "P",
    metabolomics = "M",
    lipidomics = "L"
  )
}

format_offset <- function(import_format) {
  match(import_format, c("dia", "dda", "tmt", "dia_limpa", "lc", "gc")) * 7L
}

feature_ids_for <- function(omic, import_format, fixture_class, n_features) {
  prefix <- paste0(feature_prefix_for(omic), toupper(gsub("[^A-Za-z0-9]", "", substr(import_format, 1L, 3L))))
  ids <- sprintf("%s_%03d", prefix, seq_len(n_features))
  if (identical(fixture_class, "duplicates") && length(ids) >= 2L) {
    ids[2L] <- ids[1L]
  }
  if (identical(fixture_class, "bad_names") && length(ids) >= 2L) {
    ids[1L] <- paste(ids[1L], "ambiguous")
    ids[2L] <- paste0(ids[2L], "/isoform")
  }
  ids
}

make_design <- function(sample_columns, fixture_class) {
  group <- ifelse(grepl("^WT|^WT\\b|^WT[-_ .]", sample_columns), "WT", "KO")
  design <- data.frame(
    sample = sample_columns,
    group = group,
    batch = rep(c("B1", "B2"), length.out = length(sample_columns)),
    sex = rep(c("F", "M"), length.out = length(sample_columns)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (identical(fixture_class, "multi_factor")) {
    design$site <- rep(c("site_a", "site_b"), length.out = nrow(design))
    design$age <- seq(41L, by = 3L, length.out = nrow(design))
  }
  if (identical(fixture_class, "invalid_design") && nrow(design) >= 3L) {
    design <- design[-nrow(design), , drop = FALSE]
    design[nrow(design), "sample"] <- design[1L, "sample"]
  }
  design
}

make_data <- function(omic, import_format, fixture_class) {
  feature_col <- feature_column_for(omic, import_format)
  n_features <- feature_count_for(fixture_class)
  sample_columns <- sample_columns_for(fixture_class)
  feature_ids <- feature_ids_for(omic, import_format, fixture_class, n_features)
  assay_names <- assay_names_for(omic, import_format, fixture_class)
  assay <- rep(assay_names, length.out = n_features)

  values <- outer(
    seq_len(n_features),
    seq_along(sample_columns),
    function(feature_idx, sample_idx) {
      1000L + feature_idx * 101L + sample_idx * 17L + format_offset(import_format)
    }
  )
  colnames(values) <- sample_columns

  if (identical(fixture_class, "missingness")) {
    values[n_features, seq_len(max(1L, floor(length(sample_columns) / 2L)))] <- NA_real_
    values[2L, 1L] <- NA_real_
  }
  if (identical(fixture_class, "nonfinite")) {
    values[1L, 1L] <- Inf
    if (length(sample_columns) >= 2L) {
      values[2L, 2L] <- NaN
    }
  }

  metadata <- switch(
    omic,
    proteomics = {
      if (identical(import_format, "tmt")) {
        data.frame(
          Annotated.Sequence = feature_ids,
          Master.Protein.Accessions = paste0("ACC", seq_len(n_features)),
          Gene = paste0("GENE", seq_len(n_features)),
          Assay = assay,
          ReporterIon.QC = ifelse(seq_len(n_features) %% 2L == 0L, "pass", "review"),
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      } else if (identical(import_format, "dda")) {
        data.frame(
          Protein = feature_ids,
          `Protein ID` = paste0("sp|", feature_ids, "|", feature_ids),
          Gene = paste0("GENE", seq_len(n_features)),
          Assay = assay,
          Peptide.Count = seq_len(n_features) + 1L,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      } else {
        data.frame(
          Protein.Group = feature_ids,
          Protein.Ids = paste0("sp|", feature_ids, "|", feature_ids),
          Protein.Names = paste0("Protein ", seq_len(n_features)),
          Genes = paste0("GENE", seq_len(n_features)),
          Assay = assay,
          Q.Value = seq(0.001, by = 0.002, length.out = n_features),
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      }
    },
    metabolomics = data.frame(
      Feature.Name = feature_ids,
      Annotation = paste0(ifelse(import_format == "gc", "GC", "LC"), "_compound_", seq_len(n_features)),
      m.z = round(100 + seq_len(n_features) * 1.013, 4L),
      Retention.Time = round(2 + seq_len(n_features) * 0.37, 3L),
      Assay = assay,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    lipidomics = data.frame(
      LipidName = feature_ids,
      LipidClass = rep(c("PC", "TG", "SM", "PE"), length.out = n_features),
      FattyAcid = rep(c("16:0/18:1", "18:1/18:2", "18:0/20:4"), length.out = n_features),
      Assay = assay,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  )

  cbind(metadata, as.data.frame(values, check.names = FALSE))
}

canonical_table_checksum <- function(data, columns) {
  temp <- tempfile(fileext = ".tsv")
  on.exit(unlink(temp), add = TRUE)
  utils::write.table(
    data[, columns, drop = FALSE],
    file = temp,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
  unname(tools::md5sum(temp))
}

write_tsv <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    data,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
}

make_oracle <- function(omic, import_format, fixture_class, data, design, data_path, design_path) {
  feature_col <- feature_column_for(omic, import_format)
  sample_columns <- sample_columns_for(fixture_class)
  contract <- case_contract(fixture_class)
  data_checksum_columns <- c(feature_col, sample_columns)
  design_checksum_columns <- intersect(c("sample", "group", "batch", "sex", "site", "age"), names(design))
  feature_ids <- as.character(data[[feature_col]])
  duplicate_feature_ids <- unique(feature_ids[duplicated(feature_ids)])
  nonfinite_columns <- sample_columns[vapply(
    sample_columns,
    function(column) any(!is.finite(data[[column]]) & !is.na(data[[column]])),
    logical(1)
  )]

  expected_na_by_column <- stats::setNames(
    as.list(vapply(sample_columns, function(column) sum(is.na(data[[column]])), integer(1))),
    sample_columns
  )
  expected_dropped_rows <- character()
  if (identical(fixture_class, "missingness")) {
    expected_dropped_rows <- unique(feature_ids[rowSums(is.na(data[, sample_columns, drop = FALSE])) >= 2L])
  }
  if (identical(fixture_class, "duplicates")) {
    expected_dropped_rows <- duplicate_feature_ids
  }

  list(
    schema_version = "1.0.0",
    ticket_id = "MCI-002",
    fixture_class = fixture_class,
    case_category = contract$case_category,
    omic = omic,
    import_format = import_format,
    ci_lane = contract$ci_lane,
    push_safe = contract$push_safe,
    data_path = data_path,
    design_path = design_path,
    schema = list(
      feature_id_col = feature_col,
      sample_columns = as.list(sample_columns),
      required_columns = as.list(unique(c(feature_col, "Assay", sample_columns)))
    ),
    assay_names = as.list(unique(as.character(data$Assay))),
    sample_order = as.list(sample_columns),
    feature_ids = as.list(feature_ids),
    duplicate_feature_ids = as.list(duplicate_feature_ids),
    expected_dropped_rows = as.list(expected_dropped_rows),
    expected_warnings = as.list(contract$expected_warnings),
    expected_validation_errors = as.list(contract$expected_validation_errors),
    expected_missingness_by_column = expected_na_by_column,
    expected_nonfinite_columns = as.list(nonfinite_columns),
    expected_exported_schema = list(
      long_table_columns = as.list(c(
        "omic",
        "import_format",
        "assay",
        "feature_id",
        "sample",
        "group",
        "batch",
        "intensity",
        "qc_flag"
      )),
      report_parameter_columns = as.list(c(
        "fixture_class",
        "case_category",
        "ci_lane",
        "feature_id_col",
        "sample_count",
        "feature_count"
      ))
    ),
    expected_transformed_columns = as.list(c(
      "feature_id",
      "assay",
      "sample",
      "intensity",
      "group",
      "batch",
      "source_fixture",
      "import_format"
    )),
    input_checksums = list(
      data_columns = as.list(data_checksum_columns),
      data_selected_columns_md5 = canonical_table_checksum(data, data_checksum_columns),
      design_columns = as.list(design_checksum_columns),
      design_selected_columns_md5 = canonical_table_checksum(design, design_checksum_columns)
    )
  )
}

pack_id_for <- function(omic, import_format, fixture_class) {
  paste("MCI-002", omic, gsub("_", "-", import_format), gsub("_", "-", fixture_class), sep = "-")
}

generate_module_ci_fixtures <- function(root) {
  root <- normalizePath(root, mustWork = FALSE)
  packs <- list()

  for (omic in names(import_formats)) {
    for (import_format in import_formats[[omic]]) {
      for (fixture_class in fixture_classes) {
        pack_id <- pack_id_for(omic, import_format, fixture_class)
        file_stem <- sub("^MCI-002-", "", pack_id)
        data_rel <- file.path("fixtures", omic, paste0(file_stem, "_data.tsv"))
        design_rel <- file.path("fixtures", omic, paste0(file_stem, "_design.tsv"))
        oracle_rel <- file.path("oracles", omic, paste0(file_stem, ".json"))
        data <- make_data(omic, import_format, fixture_class)
        design <- make_design(sample_columns_for(fixture_class), fixture_class)
        oracle <- make_oracle(omic, import_format, fixture_class, data, design, data_rel, design_rel)
        contract <- case_contract(fixture_class)

        write_tsv(data, file.path(root, data_rel))
        write_tsv(design, file.path(root, design_rel))
        dir.create(dirname(file.path(root, oracle_rel)), recursive = TRUE, showWarnings = FALSE)
        jsonlite::write_json(
          oracle,
          path = file.path(root, oracle_rel),
          auto_unbox = TRUE,
          pretty = TRUE,
          null = "null"
        )

        packs[[length(packs) + 1L]] <- list(
          pack_id = pack_id,
          ticket_id = "MCI-002",
          item_ids = as.list(sprintf("MCI-002.%d", 1:6)),
          omic = omic,
          import_format = import_format,
          fixture_class = fixture_class,
          case_category = contract$case_category,
          ci_lane = contract$ci_lane,
          push_safe = contract$push_safe,
          data_path = data_rel,
          design_path = design_rel,
          oracle_path = oracle_rel,
          generator = "tests/testdata/module_ci/generate-fixtures.R",
          estimated_rows = nrow(data),
          estimated_columns = ncol(data)
        )
      }
    }
  }

  manifest <- list(
    schema_version = "1.0.0",
    ticket_id = "MCI-002",
    generated_by = "tests/testdata/module_ci/generate-fixtures.R",
    deterministic_seed = 4202L,
    description = "Deterministic module-CI fixture pack catalog for omic module and transformation tests.",
    fixture_classes = fixture_class_contracts,
    import_formats = import_formats,
    packs = packs
  )
  jsonlite::write_json(
    manifest,
    path = file.path(root, "fixture_packs.json"),
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  )
  invisible(manifest)
}

args <- commandArgs(trailingOnly = TRUE)
root <- file.path("tests", "testdata", "module_ci")
if (length(args) >= 2L && identical(args[[1L]], "--root")) {
  root <- args[[2L]]
} else if (length(args) == 1L) {
  root <- args[[1L]]
}

set.seed(4202L)
generate_module_ci_fixtures(root)
