# KIT-055: the proteomics report templates state only what the record and the produced outputs support. These are
# static checks on the template sources (deterministic, no render); the rendered-text checks run in
# docs/qa/workflow-execution/kit/kit-055/render-check.R on copies of finished kit projects.
library(testthat)

template_dir <- system.file("reports", "proteomics", package = "MultiScholaR")
if (!nzchar(template_dir)) template_dir <- file.path(rprojroot::find_package_root_file(), "inst", "reports", "proteomics")
templates <- c("LFQ_report.rmd", "TMT_report.rmd", "DIANN_report.rmd", "DIANN_limpa_report.rmd")

# A glue string the template would print: the text inside the Format* functions, with the R code comments removed.
printed_text <- function(lines) paste(lines[!grepl("^\\s*#", lines)], collapse = "\n")

test_that("no proteomics report template prints a canned finding or an unfilled placeholder", {
  canned <- c("All samples passed initial QC criteria", "meets industry standards", "batches merged effectively",
    "Proceed with downstream analysis using normalized dataset", "demonstrates robust technical quality",
    "Batch effects were substantially mitigated", "showed notable batch effects")
  placeholders <- c("SOP xxx", "Your_fancy_self", "Salesforce", "SOP number(s) and name(s)")
  for (template in templates) {
    text <- printed_text(readLines(file.path(template_dir, template), warn = FALSE))
    for (sentence in c(canned, placeholders)) expect_false(grepl(sentence, text, fixed = TRUE), info = paste(template, sentence))
    expect_true(grepl("[To be completed by the analyst", text, fixed = TRUE), info = template)
    expect_true(grepl("FormatAnalystInterpretation", text, fixed = TRUE), info = template)
    expect_false(grepl("FormatComments|FormatOpinions", text), info = template)
  }
})

test_that("the templates name the outputs as the publication copy writes them and keep warnings out of the report body", {
  for (template in templates) {
    lines <- readLines(file.path(template_dir, template), warn = FALSE)
    text <- printed_text(lines)
    expect_false(grepl("NumSigDeMolecules", text, fixed = TRUE), info = template)
    expect_false(grepl("DE_proteins_results.xlsx", text, fixed = TRUE), info = template)
    expect_true(grepl("DA_results_{params$omic_type}.xlsx", text, fixed = TRUE), info = template)
    expect_false(grepl("consensus clusters|pair of time points", text), info = template)
    expect_false(any(grepl('cat(glue::glue("Warning:', lines, fixed = TRUE)), info = template)
    if (!grepl("limpa", template)) {
      expect_true(grepl("NumSigDaMolecules", text, fixed = TRUE), info = template)
      expect_true(grepl('paste0("da_", params$omic_type, "_num_sig_da_molecules.tab")', text, fixed = TRUE), info = template)
    }
  }
})

test_that("the enrichment section is conditional and the software block is the last section", {
  for (template in templates) {
    lines <- readLines(file.path(template_dir, template), warn = FALSE)
    expect_true(any(grepl("enrichment_ran <- ", lines, fixed = TRUE)), info = template)
    expect_true(any(grepl("GO Enrichment Plots, echo=FALSE, eval=enrichment_ran", lines, fixed = TRUE)), info = template)
    expect_true(any(grepl("No enrichment analysis was run for this project", lines, fixed = TRUE)), info = template)
    # Markdown headings only: lines inside code chunks (R comments) are not headings.
    in_chunk <- cumsum(grepl("^```", lines)) %% 2 == 1
    headings <- lines[grepl("^# ", lines) & !in_chunk]
    expect_identical(headings[[length(headings)]], "# Software and build", info = template)
    expect_false(any(headings == "# Version Information"), info = template)
    expect_true(any(headings == "# Analyst interpretation"), info = template)
  }
})

test_that("the LFQ template names the import format it is rendered for", {
  lines <- readLines(file.path(template_dir, "LFQ_report.rmd"), warn = FALSE)
  text <- paste(lines, collapse = "\n")
  expect_true(grepl("import_format <- tolower(as.character(params$workflow_name", text, fixed = TRUE))
  for (format in c("maxquant", "fragpipe", "spectronaut")) expect_true(grepl(paste0("    ", format, " = paste0("), text, fixed = TRUE), info = format)
  expect_false(grepl("imported from FragPipe or MaxQuant", text, fixed = TRUE))
})
