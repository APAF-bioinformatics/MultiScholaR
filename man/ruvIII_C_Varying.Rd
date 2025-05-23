% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metaboliteVsSamplesS4Objects.R
\name{ruvIII_C_Varying,MetaboliteAssayData-method}
\alias{ruvIII_C_Varying,MetaboliteAssayData-method}
\title{Apply RUV-III Correction with Varying K}
\usage{
\S4method{ruvIII_C_Varying}{MetaboliteAssayData}(
  theObject,
  ruv_grouping_variable = NULL,
  ruv_number_k = NULL,
  ctrl = NULL
)
}
\arguments{
\item{theObject}{A \code{MetaboliteAssayData} object.}

\item{ruv_grouping_variable}{Character string. The column name in the
\code{design_matrix} that defines the replicate structure for RUV-III (e.g.,
biological groups where variation \emph{within} the group is considered noise).
Defaults are looked up via \code{checkParamsObjectFunctionSimplify} using the
key \code{"ruv_grouping_variable"}. Must be provided.}

\item{ruv_number_k}{An integer or a named list/vector.
\itemize{
\item If an integer, this number of factors (\code{k}) is removed from all assays.
\item If a named list/vector, the names must correspond to the assay names
in \code{metabolite_data}. The value associated with each name specifies the
\code{k} for that assay. Assays not named will use a default \code{k}.
Defaults are looked up via \code{checkParamsObjectFunctionSimplify} using the key
\code{"ruv_number_k"}.
}}

\item{ctrl}{A logical vector, numeric vector, character vector, or a named list.
\itemize{
\item If a vector, it specifies the control features used for all assays.
Can be logical (matching features), numeric indices, or character IDs.
\item If a named list, the names must correspond to assay names. Each element
should be a vector specifying controls for that specific assay.
Defaults are looked up via \code{checkParamsObjectFunctionSimplify} using the key
\code{"ctrl"}. Must be provided.
}}

\item{num_components_to_impute}{Integer. The number of principal components
to use for NIPALS imputation if missing values are present before RUV.
Defaults are looked up via \code{checkParamsObjectFunctionSimplify} using the key
\code{"num_components_to_impute"}.}
}
\value{
A modified \code{MetaboliteAssayData} object where the \code{metabolite_data}
slot contains the RUV-corrected assay data. Features or samples with only
NA/NaN values after correction are removed. The \code{design_matrix} is updated
via \code{cleanDesignMatrix}.
}
\description{
Applies the RUV-III correction method to each assay within a
\code{MetaboliteAssayData} object. This method accounts for unwanted variation
using control features and a replicate structure matrix. It allows for
potentially different numbers of factors (\code{k}) to be removed for each assay.
}
\section{Functions}{
\itemize{
\item \code{ruvIII_C_Varying(MetaboliteAssayData)}: Method for MetaboliteAssayData

}}
