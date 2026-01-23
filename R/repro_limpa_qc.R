# Set working directory to where the R files are
setwd("/Users/ignatiuspang/Workings/2026/MultiScholaR/R")

# Load necessary libraries (suppress warnings)
suppressPackageStartupMessages({
    library(methods)
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(tibble)
    library(tidyr)
})

# Source the files
# We need to suppress the output of source if possible, or just run it
source("helper_functions.R")
source("proteinVsSamplesS4Objects.R")
# We need to mock 'project_dirs' if it's used globally in generateLimpaQCPlots?
# generateLimpaQCPlots takes 'save_dir' argument.
source("limpa_functions.R")

# Create a mock ProteinQuantitativeData object
# args is where the issue lies
args_list <- list(
    # Simulating inheritance: This object has BOTH results!
    limpa_dpc_results = list(
        dpc_parameters = c(1, 0.5), # intercept, slope
        slope_interpretation = "test mechanism"
    ),
    limpa_dpc_quant_results = list(
        dpc_parameters_used = c(1, 0.5),
        slope_interpretation = "test mechanism",
        y_peptide_for_dpc = matrix(rnorm(100), nrow = 10, ncol = 10), # Pseudo peptide data
        dpc_method = "limpa_dpc_quant"
    )
)

# Mock protein data
protein_data <- data.frame(
    Protein.Ids = paste0("Prot", 1:10),
    Sample1 = rnorm(10),
    Sample2 = rnorm(10)
)
colnames(protein_data)[2:3] <- c("S1", "S2")

design_matrix <- data.frame(
    Sample_id = c("S1", "S2"),
    group = c("G1", "G1")
)

# Initialize object
obj <- new("ProteinQuantitativeData",
    protein_quant_table = protein_data,
    protein_id_column = "Protein.Ids",
    design_matrix = design_matrix,
    sample_id = "Sample_id",
    group_id = "group",
    technical_replicate_id = "replicates", # Optional?
    args = args_list
)

# Attempt to run generateLimpaQCPlots
# We expect this to FAIL to produce meaningful plots because it will pick 'limpa_dpc_results'
# and try to access @peptide_matrix
message("Running generateLimpaQCPlots...")
plots <- generateLimpaQCPlots(obj, save_plots = FALSE, verbose = TRUE)

message("\nChecking plots...")
# Check if 'missing_comparison' is the "Unavailable" plot
# We can check the labels in the plot
check_plot <- function(p, name) {
    if (inherits(p, "ggplot")) {
        # Extract data or layers to see if it's a dummy plot
        # Dummy plots usually have geom_text with specific label
        # or empty data
        is_dummy <- any(grepl("unavailable", p$layers[[1]]$data$label))
        if (is_dummy) {
            message(paste("FAIL:", name, "is unavailable dummy plot."))
        } else {
            # Also check if it's the "Comparison requires before/after" dummy
            # In the failure mode I suspect, it enters "peptide_imputation" mode
            # checks for peptide_matrix, fails, and returns dummy.
            message(paste("PASS:", name, "seems to be generated."))
        }
    } else {
        message(paste("FAIL:", name, "is not a ggplot."))
    }
}

# In the current buggy version:
# It should pick "peptide_imputation"
# It tries to access obj@peptide_matrix -> FAILS (returns NULL)
# missing_comparison -> NULL -> generates "Unavailable" plot

if (!is.null(plots$missing_comparison)) {
    # Inspect the plot title or content
    print(plots$missing_comparison$labels$title)
}
# We expect "Missing Values Comparison" title, but the content will be "Missing value comparison unavailable"
