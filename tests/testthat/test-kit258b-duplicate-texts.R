# KIT-258 follow-up (VAL-306): the four in-app texts state the duplicate rule as merged (VAL-254, VAL-284, VAL-337):
# among the rows of one duplicated identifier, the row measured (neither missing nor 0) in the most samples is kept; among
# rows measured in equally many, the one with the highest mean of its measured values; if still tied, the first row in the
# table. The texts are: (1) the Duplicates tab's help, (2) the step's result text, (3) the kernels' and S4 methods' console
# messages, (4) the exported kernels' help titles (roxygen and man/). None may still say "highest average intensity"
# (KIT-214's rule). The wordings are the audit's drafts (AUDIT-KIT-258-impl-20261126 I01) and the owner's to set, so each
# expectation matches the rule's parts, not one exact sentence. Tests only: no value stored anywhere is read here.

k258bFlat <- function(text) gsub("[[:space:]]+", " ", paste(text, collapse = " "))

# The rule's parts as a text can state them: measured in the most samples; then the highest average (mean) of the
# measured values; then the table's order. The earlier rule's words.
k258bMost <- "(measured in the most|most completely measured|measured in more)"
k258bThenMean <- "highest (average|mean) of (its|the) measured values"
k258bTableOrder <- "first (row )?in the table"
k258bOld <- "highest average intensity|lower-intensity|highest intensity"

# ---- (1) The Duplicates tab's help ------------------------------------------------------------------------------------

for (omic in c("lipidomics", "metabolomics")) {
  local({
    omic <- omic
    test_that(sprintf("KIT-258b (1) %s: the Duplicates help states the rule as merged", omic), {
      help <- k258bFlat(omicGuiConfig(omic)$qc$texts[["duplicates_help"]])
      expect_match(help, k258bMost, ignore.case = TRUE)
      expect_match(help, k258bThenMean, ignore.case = TRUE)
      expect_match(help, k258bTableOrder, ignore.case = TRUE)
      expect_false(grepl(k258bOld, help, ignore.case = TRUE))
    })
  })
}

# ---- (2) The step's result text ---------------------------------------------------------------------------------------

for (prefix in c("lipid", "metab")) {
  local({
    prefix <- prefix
    test_that(sprintf("KIT-258b (2) %s: the step's result text states the rule as merged", prefix), {
      resolved <- list(LCMS_Pos = list(original = 5L, resolved = 4L, removed = 1L))
      text <- k258bFlat(omicQcDuplicatesText(resolved, prefix))
      expect_match(text, k258bMost, ignore.case = TRUE)
      expect_match(text, k258bThenMean, ignore.case = TRUE)
      expect_false(grepl(k258bOld, text, ignore.case = TRUE))
    })
  })
}

# ---- (3) The console messages -----------------------------------------------------------------------------------------

# Six design samples, two groups of three; id D twice (row A measured in 2 of 6, row B in 6 of 6), id E once.
k258bSamples <- sprintf("S%d", 1:6)
k258bTable <- function() {
  values <- rbind(c(200000, 220000, NA, NA, NA, NA), c(140000, 145000, 150000, 150000, 155000, 160000),
    c(1000, 1100, 1200, 1300, 1400, 1500))
  colnames(values) <- k258bSamples
  cbind(data.frame(feature_id = c("D", "D", "E"), name = c("A", "B", "C"), stringsAsFactors = FALSE),
    as.data.frame(values))
}
k258bDesign <- function() {
  data.frame(Run = k258bSamples, group = rep(c("A", "B"), each = 3L),
    tech_rep_group = paste(rep(c("A", "B"), each = 3L), rep(1:3, 2L), sep = "_"), stringsAsFactors = FALSE)
}

# The messages must say the rule once, say what was removed without the earlier rule's words, and none may carry them.
k258bExpectMessages <- function(messages, removed_pattern, info) {
  flat <- k258bFlat(messages)
  expect_true(any(grepl(k258bMost, messages, ignore.case = TRUE)), info = info)
  expect_true(any(grepl(removed_pattern, messages)), info = info)
  expect_false(grepl(k258bOld, flat, ignore.case = TRUE), info = info)
}

for (kernel in c("resolveDuplicateFeaturesByIntensity", "resolveLipidDuplicateFeaturesByIntensity")) {
  local({
    kernel <- kernel
    test_that(sprintf("KIT-258b (3) %s: its console messages state the rule as merged", kernel), {
      messages <- testthat::capture_messages(
        get(kernel, asNamespace("MultiScholaR"))(k258bTable(), "feature_id", k258bSamples))
      k258bExpectMessages(messages, "Removed 1 duplicate feature row", kernel)
    })
  })
}

for (omic in c("lipidomics", "metabolomics")) {
  local({
    omic <- omic
    test_that(sprintf("KIT-258b (3) %s: the S4 method resolveDuplicateFeatures's console messages state the rule", omic), {
      table <- k258bTable()
      # The S4 methods require their descriptive-name column (`lipid`, `metabolite`); no internal-standard regex.
      table[[if (identical(omic, "lipidomics")) "lipid" else "metabolite"]] <- table$name
      object <- if (identical(omic, "lipidomics")) {
        createLipidomicsAssayData(lipid_data = list(LCMS_Pos = table), design_matrix = k258bDesign(),
          lipid_id_column = "feature_id", annotation_id_column = "name", sample_id = "Run", group_id = "group",
          technical_replicate_id = "tech_rep_group")
      } else {
        createMetaboliteAssayData(metabolite_data = list(LCMS_Pos = table), design_matrix = k258bDesign(),
          metabolite_id_column = "feature_id", annotation_id_column = "name", sample_id = "Run", group_id = "group",
          technical_replicate_id = "tech_rep_group")
      }
      messages <- testthat::capture_messages(suppressWarnings(resolveDuplicateFeatures(object)))
      k258bExpectMessages(messages, "Removed 1 duplicate non-ITSD feature row", omic)
    })
  })
}

# ---- (4) The exported kernels' help titles: roxygen comment and the hand-edited man/ page -------------------------------

k258bTitles <- list(
  resolveDuplicateFeaturesByIntensity = "func_omics_qc_shared.R",
  resolveLipidDuplicateFeaturesByIntensity = "func_lipid_qc_support_helpers.R")

for (name in names(k258bTitles)) {
  local({
    name <- name
    source_file <- k258bTitles[[name]]
    test_that(sprintf("KIT-258b (4) %s: the help title (man/ and roxygen) states the rule as merged", name), {
      root <- normalizePath(testthat::test_path("..", ".."), mustWork = FALSE)
      rd <- file.path(root, "man", paste0(name, ".Rd"))
      code <- file.path(root, "R", source_file)
      testthat::skip_if_not(file.exists(rd) && file.exists(code), "the package source tree is not present")
      rd_lines <- readLines(rd, warn = FALSE)
      rd_title <- sub("^\\\\title\\{(.*)\\}$", "\\1", grep("^\\\\title\\{", rd_lines, value = TRUE))
      code_lines <- readLines(code, warn = FALSE)
      at <- grep(paste0("^", name, " <- function"), code_lines)
      # The roxygen title is the first "#' " line of the comment block that ends at the function.
      block_start <- max(grep("^# -{20,}", code_lines[seq_len(at[[1L]])]))
      title_line <- code_lines[block_start - 1L + which(startsWith(code_lines[block_start:at[[1L]]], "#' "))[[1L]]]
      roxygen_title <- sub("^#' ", "", title_line)
      expect_match(rd_title, "most completely measured", ignore.case = TRUE)
      expect_false(grepl("highest average intensity", rd_title, ignore.case = TRUE))
      expect_match(roxygen_title, "most completely measured", ignore.case = TRUE)
      expect_false(grepl("highest average intensity", roxygen_title, ignore.case = TRUE))
    })
  })
}
