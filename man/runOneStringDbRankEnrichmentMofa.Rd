% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multiomics_enrichment_functions.R
\name{runOneStringDbRankEnrichmentMofa}
\alias{runOneStringDbRankEnrichmentMofa}
\title{Run STRING DB Rank Enrichment Analysis for MOFA Factors or Similar Data}
\usage{
runOneStringDbRankEnrichmentMofa(
  input_table,
  identifier_column_name = "protein_id",
  value_column_name = "score",
  result_label,
  results_dir,
  api_key = NULL,
  species = "9606",
  ge_fdr = 0.05,
  ge_enrichment_rank_direction = -1,
  polling_interval_seconds = 10,
  max_polling_attempts = 30
)
}
\arguments{
\item{input_table}{A data frame containing the data to be submitted for enrichment.
Must include columns specified by \code{identifier_column_name} and \code{value_column_name}.}

\item{identifier_column_name}{Character string. The name of the column in
\code{input_table} containing the identifiers (e.g., protein IDs, gene symbols).
Default: \code{"protein_id"}.}

\item{value_column_name}{Character string. The name of the column in \code{input_table}
containing the numerical values (e.g., scores, weights, logFC) associated
with each identifier. This column must be numeric. Default: \code{"score"}.}

\item{result_label}{Character string. A label used for naming output files
(e.g., "MOFA_Factor1_Proteomics").}

\item{results_dir}{Character string. The path to the directory where enrichment
results (table, URL file, graph image) will be saved. The directory will be
created if it doesn't exist.}

\item{api_key}{Character string or NULL. Your personal STRING API key.
If NULL, the API might work with limitations or require it for full access.
Default: \code{NULL}.}

\item{species}{Character string. NCBI/STRING species identifier.
Default: \code{"9606"} (Homo sapiens).}

\item{ge_fdr}{Numeric. FDR threshold for gene expression enrichment.
Default: \code{0.05}.}

\item{ge_enrichment_rank_direction}{Integer. Direction for enrichment rank
(-1, 0, or 1). Default: \code{-1}.}

\item{polling_interval_seconds}{Numeric. Seconds to wait between polling attempts
for job status. Default: \code{10}.}

\item{max_polling_attempts}{Numeric. Maximum polling attempts before timing out.
Default: \code{30}.}
}
\value{
A data frame containing the enrichment results from STRING DB
(specifically \code{output_tbl$enrichment_data} from the retrieval step).
Returns \code{NULL} or an empty structure if the process fails.
}
\description{
This function automates the process of submitting data to the STRING database
for values/ranks enrichment analysis, polling for results, and saving the
output. It is a wrapper around \code{submitStringDBEnrichment} and
\code{retrieveStringDBEnrichmentResults}, tailored for analyses where inputs
might come from MOFA factor weights or similar ranked lists.

Results, including the enrichment table, a URL list, and a graph image,
are saved to the specified \code{results_dir}.
}
\seealso{
\code{\link{submitStringDBEnrichment}}, \code{\link{retrieveStringDBEnrichmentResults}}
}
