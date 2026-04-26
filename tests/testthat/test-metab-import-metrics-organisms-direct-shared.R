# fidelity-coverage-compare: shared
library(testthat)

makeFunctionWithOverrides <- function(fun, replacements) {
  fun_override <- fun
  environment(fun_override) <- list2env(replacements, parent = environment(fun))
  fun_override
}

createMetaboliteAssayData <- get(
  "createMetaboliteAssayData",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
getMetaboliteQuantData <- get(
  "getMetaboliteQuantData",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
importMetabMSDIALData <- get(
  "importMetabMSDIALData",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
countUniqueMetabolites <- get(
  "countUniqueMetabolites",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
countMetabolitesPerSample <- get(
  "countMetabolitesPerSample",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
calculateMissingness <- get(
  "calculateMissingness",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
calculateSumIntensityPerSample <- get(
  "calculateSumIntensityPerSample",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
calculateTotalUniqueMetabolitesAcrossAssays <- get(
  "calculateTotalUniqueMetabolitesAcrossAssays",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
createMetabDaResultsLongFormat <- get(
  "createMetabDaResultsLongFormat",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
extractOrganismsFromFasta <- get(
  "extractOrganismsFromFasta",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
analyzeOrganismDistribution <- get(
  "analyzeOrganismDistribution",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)

test_that("metabolomics import helpers preserve constructor, sample extraction, and MS-DIAL parsing", {
  skip_if_not_installed("vroom")
  skip_if_not_installed("seqinr")

  assay_data <- data.frame(
    metabolite_id = c("M1", "M2", NA_character_),
    annotation = c("alpha", "beta", "gamma"),
    mz = c(101.1, 202.2, 303.3),
    S1 = c(10, 0, NA_real_),
    S2 = c(20, 5, 7),
    note = c("a", "b", "c"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  design_matrix <- data.frame(
    Sample_ID = c("S1", "S2"),
    group = c("A", "B"),
    stringsAsFactors = FALSE
  )

  assay_object <- createMetaboliteAssayData(
    metabolite_data = list(LCMS_Pos = assay_data),
    design_matrix = design_matrix,
    metabolite_id_column = "metabolite_id",
    annotation_id_column = "annotation",
    sample_id = "Sample_ID",
    group_id = "group",
    database_identifier_type = "HMDB",
    internal_standard_regex = "^ISTD",
    args = list(source = "direct-shared")
  )

  expect_s4_class(assay_object, "MetaboliteAssayData")
  expect_identical(assay_object@sample_id, "Sample_ID")
  expect_identical(assay_object@metabolite_id_column, "metabolite_id")
  expect_identical(assay_object@args$source, "direct-shared")

  explicit_quant <- getMetaboliteQuantData(assay_data, sample_columns = c("S1", "S2"))
  expect_identical(colnames(explicit_quant$quant_data), c("S1", "S2"))
  expect_true(all(c("metabolite_id", "annotation", "mz", "note") %in% colnames(explicit_quant$annotation_data)))

  expect_warning(
    fallback_quant <- getMetaboliteQuantData(assay_data, sample_columns = "missing_sample"),
    "None of the provided sample_columns exist in assay_data",
    fixed = TRUE
  )
  expect_identical(fallback_quant$sample_names, c("mz", "S1", "S2"))

  guessed_quant <- getMetaboliteQuantData(assay_data)
  expect_identical(guessed_quant$sample_names, c("mz", "S1", "S2"))

  import_msdial <- makeFunctionWithOverrides(
    importMetabMSDIALData,
    list(
      getMetabolomicsColumnDefaults = function(format_name) {
        expect_identical(format_name, "msdial")
        list(
          metabolite_id = c("Peak ID"),
          annotation = c("Metabolite name"),
          rt = c("RT (min)"),
          mz = c("Precursor m/z"),
          adduct = c("Adduct"),
          is_pattern = FALSE
        )
      },
      findMetabMatchingColumn = function(headers, candidates) {
        matches <- intersect(headers, candidates)
        if (length(matches) == 0) {
          return(NULL)
        }
        matches[[1]]
      }
    )
  )

  expect_error(
    import_msdial(tempfile("missing-msdial-")),
    "File not found",
    fixed = TRUE
  )

  csv_path <- tempfile(fileext = ".csv")
  writeLines(
    c(
      "Peak ID,Metabolite name,RT (min),Precursor m/z,Adduct,Height,Sample_A,Sample_B,Comment",
      "1,Alanine,2.1,100.1,[M+H]+,999,100,200,first",
      "2,Glycine,3.2,150.2,[M+Na]+,888,110,210,second"
    ),
    csv_path
  )
  withr::defer(unlink(csv_path, force = TRUE))

  imported <- import_msdial(csv_path)
  expect_identical(imported$detected_columns$metabolite_id, "Peak ID")
  expect_identical(imported$detected_columns$annotation, "Metabolite name")
  expect_identical(imported$sample_columns, c("Sample_A", "Sample_B"))
  expect_true(all(c("Peak ID", "Height", "Comment") %in% imported$annotation_columns))
  expect_identical(imported$format, "msdial")
  expect_false(imported$is_pattern)
})

test_that("metabolomics QC metrics and DA long-format helpers preserve counting and intensity expansion behavior", {
  assay_data <- data.frame(
    metabolite_id = c("M1", "M2", NA_character_, "M1"),
    annotation = c("alpha", "beta", "gamma", "delta"),
    mz = c(101.1, 202.2, 303.3, 404.4),
    S1 = c(10, 0, NA_real_, 4),
    S2 = c(20, 5, 7, 0),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  expect_identical(countUniqueMetabolites(assay_data, "metabolite_id"), 2L)
  expect_warning(
    missing_count <- countUniqueMetabolites(assay_data, "missing_id"),
    "Metabolite ID column 'missing_id' not found in assay data.",
    fixed = TRUE
  )
  expect_identical(missing_count, 0)

  per_sample <- countMetabolitesPerSample(
    assay_data,
    sample_id_col = "Run",
    metabolite_id_col = "metabolite_id",
    sample_columns = c("S1", "S2")
  )
  expect_identical(per_sample$Run, c("S1", "S2"))
  expect_identical(per_sample$n_detected, c(2L, 3L))

  missingness <- calculateMissingness(
    assay_data,
    sample_id_col = "Run",
    sample_columns = c("S1", "S2")
  )
  expect_equal(missingness, 37.5)

  expect_warning(
    empty_missingness <- calculateMissingness(
      assay_data["annotation"],
      sample_id_col = "Run",
      sample_columns = "missing_sample"
    ),
    "No valid data for missingness calculation",
    fixed = TRUE
  )
  expect_true(is.na(empty_missingness))

  per_sample_sum <- calculateSumIntensityPerSample(
    assay_data,
    sample_id_col = "Run",
    sample_columns = c("S1", "S2")
  )
  expect_identical(per_sample_sum$Run, c("S1", "S2"))
  expect_equal(per_sample_sum$sum_intensity, c(14, 32))

  total_unique <- calculateTotalUniqueMetabolitesAcrossAssays(
    list(
      assay_data,
      data.frame(
        metabolite_id = c("M3", "M2"),
        stringsAsFactors = FALSE
      )
    ),
    metabolite_id_col = "metabolite_id"
  )
  expect_identical(total_unique, 3L)
  expect_identical(calculateTotalUniqueMetabolitesAcrossAssays(list(), "metabolite_id"), 0)

  expr_matrix <- matrix(
    c(
      10, 20, 30, 40,
      15, 25, 35, 45
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("M1", "M2"), c("S1", "S2", "S3", "S4"))
  )
  design_matrix <- data.frame(
    sample_id = c("S1", "S2", "S3", "S4"),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )
  lfc_qval_tbl <- data.frame(
    metabolite_id = c("M1", "M2"),
    comparison = c("A-B", "A-B"),
    logFC = c(1.5, -0.5),
    fdr_qvalue = c(0.01, 0.2),
    stringsAsFactors = FALSE
  )

  da_long <- createMetabDaResultsLongFormat(
    lfc_qval_tbl = lfc_qval_tbl,
    expr_matrix = expr_matrix,
    design_matrix = design_matrix,
    sample_id_col = "sample_id",
    group_id_col = "group",
    metabolite_id_col = "metabolite_id"
  )

  expect_true(all(c("numerator", "denominator", "intensity.S1.A", "intensity.S4.B") %in% colnames(da_long)))
  expect_identical(da_long$numerator, c("A", "A"))
  expect_identical(da_long$denominator, c("B", "B"))
  expect_equal(da_long$intensity.S1.A, c(10, 15))
  expect_equal(da_long$intensity.S4.B, c(40, 45))
})

test_that("organism import helpers preserve FASTA parsing and unmatched distribution reporting", {
  skip_if_not_installed("seqinr")

  fasta_path <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">sp|P12345|PROT_HUMAN OS=Homo sapiens OX=9606 GN=GENE1",
      "MAAA",
      ">tr|Q99999|PROT_MOUSE OS=Mus musculus OX=10090",
      "MBBB",
      ">CUSTOM_HEADER OS=Unknown species",
      "MCCC"
    ),
    fasta_path
  )
  withr::defer(unlink(fasta_path, force = TRUE))

  expect_error(
    extractOrganismsFromFasta(tempfile("missing-fasta-")),
    "FASTA file not found",
    fixed = TRUE
  )

  organism_tbl <- extractOrganismsFromFasta(fasta_path)
  expect_identical(organism_tbl$uniprot_acc, c("P12345", "Q99999", "CUSTOM_HEADER"))
  expect_identical(organism_tbl$organism_name, c("Homo sapiens", "Mus musculus", "Unknown species"))
  expect_identical(organism_tbl$taxon_id, c(9606L, 10090L, NA_integer_))

  analyze_distribution <- makeFunctionWithOverrides(
    analyzeOrganismDistribution,
    list(
      normalizeUniprotAccession = function(acc, remove_isoform = TRUE) {
        sub("-\\d+$", "", acc)
      }
    )
  )

  distribution <- analyze_distribution(
    protein_ids = c("P12345-2;Q99999", "CUSTOM_HEADER", "UNMATCHED"),
    organism_mapping = organism_tbl
  )

  expect_true(all(c("organism_name", "taxon_id", "protein_count", "percentage", "rank") %in% colnames(distribution)))
  expect_equal(distribution$protein_count[distribution$organism_name == "Homo sapiens"], 1)
  expect_equal(distribution$protein_count[distribution$organism_name == "Mus musculus"], 1)
  expect_equal(distribution$protein_count[distribution$organism_name == "Unknown species"], 1)
  expect_equal(distribution$protein_count[distribution$organism_name == "[Unmatched/Unknown]"], 1)
  expect_equal(sum(distribution$protein_count), 4)
  expect_true(all(diff(distribution$rank) >= 0))
})
