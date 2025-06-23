\name{findBestNegCtrlPercentage}
\alias{findBestNegCtrlPercentage}
\title{Find the Best Negative Control Percentage for RUV-III Analysis}
\usage{
findBestNegCtrlPercentage(
  normalised_protein_matrix_obj,
  percentage_range = seq(1, 20, by = 1),
  num_components_to_impute = 5,
  ruv_grouping_variable = "group",
  ruv_qval_cutoff = 0.05,
  ruv_fdr_method = "qvalue",
  separation_metric = "max_difference",
  k_penalty_weight = 0.5,
  max_acceptable_k = 3,
  adaptive_k_penalty = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{normalised_protein_matrix_obj}{A ProteinQuantitativeData object containing the normalized protein quantification data}

\item{percentage_range}{A numeric vector specifying the range of percentages to test. Default is seq(1, 20, by = 1) for testing 1\% to 20\% in 1\% increments}

\item{num_components_to_impute}{Number of components to use for imputation in ruvCancor. Default is 5}

\item{ruv_grouping_variable}{The grouping variable to use for RUV analysis. Default is "group"}

\item{ruv_qval_cutoff}{The FDR threshold for negative control selection. Default is 0.05}

\item{ruv_fdr_method}{The FDR calculation method. Default is "qvalue"}

\item{separation_metric}{The metric to use for evaluating separation quality. Options: "max_difference" (default), "mean_difference", "auc", "weighted_difference"}

\item{k_penalty_weight}{Weight for penalizing high k values in composite score. Default is 0.5. Higher values penalize high k more strongly}

\item{max_acceptable_k}{Maximum acceptable k value. k values above this get heavy penalty. Default is 3 (ignored when adaptive_k_penalty = TRUE)}

\item{adaptive_k_penalty}{Whether to automatically adjust max_acceptable_k based on sample size. Default is TRUE (recommended). Set to FALSE only if you need exact reproducibility with previous results}

\item{verbose}{Whether to print progress messages. Default is TRUE}
}
\value{
A list containing:
\itemize{
  \item best_percentage: The optimal percentage as a numeric value
  \item best_k: The optimal k value from findBestK() for the best percentage
  \item best_control_genes_index: The control genes index for the best percentage
  \item best_separation_score: The separation score for the best percentage
  \item best_composite_score: The composite score (separation penalized by k value)
  \item optimization_results: A data frame with all tested percentages and their scores
  \item best_cancor_plot: The canonical correlation plot for the best percentage
  \item separation_metric_used: The separation metric that was used
  \item k_penalty_weight: The k penalty weight that was used
  \item max_acceptable_k: The maximum acceptable k value that was used
}
}
\description{
This function automatically determines the optimal percentage of proteins to use
as negative controls for RUV-III analysis by testing different percentages and
evaluating both the separation quality between "All" and "Control" groups and
the resulting k value from findBestK(). This prevents over-optimization towards
percentages that give good separation but unreasonably high k values that would
remove biological signal. This removes the subjectivity from manually selecting
the percentage_as_neg_ctrl parameter.
}
\details{
The function works by:
\enumerate{
  \item Iterating through a range of percentage values (default 1\% to 20\%)
  \item For each percentage, generating negative control genes using \code{getNegCtrlProtAnova}
  \item Creating canonical correlation plots using \code{ruvCancor}
  \item Calculating the optimal k value using \code{findBestK}
  \item Evaluating separation quality using one of four metrics:
    \itemize{
      \item "max_difference": Maximum difference between All and Control curves
      \item "mean_difference": Average difference across all k values
      \item "auc": Area under the curve of differences
      \item "weighted_difference": K-weighted average of differences
    }
  \item Computing a composite score that balances separation quality with k value penalty
  \item Returning the percentage that gives the best composite score
}

The composite scoring mechanism prevents over-correction by penalizing high k values:
\itemize{
  \item k values <= max_acceptable_k get linear penalty (0 to k_penalty_weight)
  \item k values > max_acceptable_k get exponential penalty (approaching 1)
  \item Final score = separation_score * (1 - k_penalty)
}

The function requires at least 5 control genes for each percentage tested to ensure
reliable RUV-III analysis. Percentages that result in fewer control genes are skipped.
}
\examples{
\dontrun{
# Basic usage with default parameters
result <- findBestNegCtrlPercentage(normalised_frozen_protein_matrix_obj)

# Use the optimal percentage, k value, and control genes
best_percentage <- result$best_percentage
best_k <- result$best_k  # No need to call findBestK() separately!
control_genes_index <- result$best_control_genes_index

# View the optimization results for all tested percentages
print(result$optimization_results)

# Display the best canonical correlation plot
print(result$best_cancor_plot)

# Custom percentage range, separation metric, and k penalty settings
result_custom <- findBestNegCtrlPercentage(
  normalised_frozen_protein_matrix_obj,
  percentage_range = seq(2, 15, by = 0.5),
  separation_metric = "auc",
  k_penalty_weight = 0.7,  # Higher penalty for high k values
  max_acceptable_k = 2,    # More conservative k threshold
  verbose = TRUE
)

# Integration with existing workflow
config_list$ruvParameters$percentage_as_neg_ctrl <- result$best_percentage
control_genes_index <- result$best_control_genes_index
best_k <- result$best_k

# Continue with existing RUV-III analysis using optimal parameters
}
} 