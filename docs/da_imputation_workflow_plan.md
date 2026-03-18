# Differential Abundance (DA) Imputation Workflow Plan

## 1. Objective
To construct a robust Differential Abundance (DA) testing logic that dynamically adapts to the missing-value imputation strategy chosen by the user. Most notably, this workflow will incorporate specialized handling for features imputed via `limpa`, ensuring standard error uncertainty is propagated into the DA statistics.

## 2. Core Branching Logic

Before running DA testing, the workflow must define the imputation context (e.g., checking user-selected analytical arguments or metadata).

### Branch 1: `limpa` Specific Pathway
If the data was imputed using `limpa::dpcImpute()` (or `dpcQuant()`), the resultant quantitative object will be a specialized `EList` containing both the imputed expression values (`$E`) and standard errors (`$other$standard.error`). 

*   **Function to Use**: `limpa::dpcDE(data, design)`
*   **Reason**: `dpcDE()` leverages the matrix of standard errors calculated during imputation to correctly penalize the standard deviations of features heavily reliant on missing value imputation. Standard `limma` functions ignore these standard errors, which would result in artificially inflated confidence levels.

### Branch 2: Standard `limma` Pathway
If the data was imputed using any other method (e.g., `knn`, `missForest`), or if no imputation was performed (`none`), the quantitative data remains a standard data frame or matrix of intensities where all non-missing values are assumed equally certain.

*   **Function to Use**: `limma::lmFit(data, design)`
*   **Reason**: `lmFit` fits multiple linear models (one per feature) to the expression matrix given the defined design matrix. This is the canonical start of the `limma` pipeline for normally prepared continuous data.

## 3. Downstream Processing (Common Pipeline)

Regardless of the branching logic described above, both `dpcDE()` and `lmFit()` output compatible `fit` objects that funnel back into the standard `limma` pipeline.

1.  **Empirical Bayes Moderation**:
    ```R
    fit <- limma::eBayes(fit)
    ```
    This applies Moderated t-statistics to the fitted models, allowing information to be borrowed across features to stabilize variance estimates.

2.  **Results Extraction**:
    ```R
    results <- limma::topTable(fit, coef = "your_contrast", number = Inf)
    ```
    This extracts the final DA results (fold changes, t-statistics, p-values, and adjusted p-values) for the desired experimental contrasts.
