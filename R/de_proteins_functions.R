#' @title Remove Rows with All-Missing Values
#' @description This function filters a data frame to remove rows where all specified
#' numeric columns (matching a pattern) contain either `NA` or `0`.
#'
#' @details The function identifies columns based on `col_pattern`. For each row, it
#' checks if all values in these columns are `NA` or `0`. Rows that meet this
#' condition are removed from the data frame.
#'
#' @param input_table A data frame or tibble with columns containing protein/feature abundances.
#' @param col_pattern A regular expression string that identifies the columns to check
#'   for missing values (e.g., "Reporter.intensity.corrected").
#' @param row_id The column name that serves as the unique identifier for each row.
#'   This should be an unquoted column name (e.g., `Protein.Ids`).
#'
#' @return A data frame identical in structure to `input_table` but with rows
#'   containing only missing or zero values removed.
#'
#' @importFrom dplyr mutate if_all inner_join select
#' @importFrom rlang sym enquo as_string as_name
#' @export
#'
#' @examples
#' df <- data.frame(
#'   Protein.ID = c("A", "B", "C", "D"),
#'   Sample1 = c(10, 0, 5, 0),
#'   Sample2 = c(20, NA, 8, 0),
#'   Sample3 = c(30, 0, 0, NA)
#' )
#' # Row B and D should be removed
#' removeEmptyRows(df, col_pattern = "Sample", row_id = Protein.ID)
removeEmptyRows <- function(input_table, col_pattern, row_id) {

  temp_col_name <- paste0("temp_", as_string(as_name(enquo(row_id))))

  temp_input_table <- input_table |>
    dplyr::mutate(!!rlang::sym(temp_col_name) := row_number())

  sites_to_accept <- temp_input_table |>
    mutate(across(matches(col_pattern, perl = TRUE), \(x){ (is.na(x) | x == 0) })) |>
    dplyr::filter(!if_all(matches(col_pattern, perl = TRUE), \(x){ x == TRUE })) |>
    dplyr::select({ { temp_col_name } })

  ## Removing entries where all the "Reporter intensity corrected" rows are zero
  filtered_table <- temp_input_table |>
    inner_join(sites_to_accept, by = temp_col_name) |>
    dplyr::select(-temp_col_name)

  return(filtered_table)
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' @title Plot Number of Missing Values per Sample
#' @description This function takes a data matrix, calculates the number of missing
#' values (`NA` or non-finite after log2 transformation) for each sample (column),
#' and generates a bar plot visualizing these counts.
#'
#' @details The function first applies a `log2` transformation to the input data.
#' It then counts values that are not finite (i.e., `NA`, `NaN`, `Inf`, `-Inf`),
#' which is a common way to define missingness after log transformation.
#'
#' @param input_table A data frame or matrix where rows are features (e.g., proteins)
#'   and columns are samples. The function expects numeric data suitable for
#'   log2 transformation.
#'
#' @return A ggplot2 bar plot object showing the number of missing values for each
#'   sample. The x-axis represents sample IDs (original column names) and the
#'   y-axis represents the count of missing values.
#'
#' @importFrom ggplot2 ggplot aes geom_bar theme element_text
#' @importFrom tibble rownames_to_column
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' data_matrix <- data.frame(
#'   SampleA = c(10, 20, 0, 40),
#'   SampleB = c(5, NaN, 15, 25),
#'   SampleC = c(NA, 10, 30, 50)
#' )
#' # After log2, 0 becomes -Inf. So SampleA has 1 missing, SampleB has 1, SampleC has 1.
#' plotNumMissingValues(data_matrix)
plotNumMissingValues <- function(input_table) {

  plot_num_missing_values <- apply(data.matrix(log2(input_table)), 2,
                                   function(x) { length(which(!is.finite(x))) }) |>
    t() |>
    t() |>
    set_colnames("No. of Missing Values") |>
    as.data.frame() |>
    rownames_to_column("Samples ID") |>
    ggplot(aes(x = `Samples ID`, y = `No. of Missing Values`)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90))

  plot_num_missing_values
}




## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' @title Plot Number of Valid Values per Sample (Log-Transformed)
#' @description This function counts the number of valid (non-NA) values for each
#' sample (column) after a log2 transformation and creates a bar plot of these counts.
#'
#' @param input_table A data frame or matrix where rows are features and columns are
#'   samples. Data should be numeric.
#'
#' @return A ggplot2 bar plot showing the number of valid data points for each sample
#'   after log2 transformation.
#' @export
plotNumOfValues <- function(input_table) {

  plot_num_missing_values <- apply(data.matrix(log2(input_table)), 2,
                                   function(x) { length(which(!is.na(x))) }) |>
    t() |>
    t() |>
    set_colnames("No. of Values") |>
    as.data.frame() |>
    rownames_to_column("Samples ID") |>
    ggplot(aes(x = `Samples ID`, y = `No. of Values`)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90))

  plot_num_missing_values
}

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @title Plot Number of Valid Values per Sample (No Transformation)
#' @description This function counts the number of valid (non-NA) values for each
#' sample (column) from the raw input data and creates a bar plot of these counts.
#'
#' @param input_table A data frame or matrix where rows are features and columns are
#'   samples.
#'
#' @return A ggplot2 bar plot showing the number of valid data points for each sample.
#' @export
plotNumOfValuesNoLog <- function(input_table) {

  plot_num_missing_values <- apply(data.matrix(input_table), 2,
                                   function(x) { length(which(!is.na(x))) }) |>
    t() |>
    t() |>
    set_colnames("No. of Values") |>
    as.data.frame() |>
    rownames_to_column("Samples ID") |>
    ggplot(aes(x = `Samples ID`, y = `No. of Values`)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90))

  plot_num_missing_values
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @title Filter Rows Based on Maximum Missing Values per Group
#' @description This function filters a data table by removing rows (features) that
#' have more than a specified number of missing values within any experimental group.
#'
#' @details The function first reshapes the data into a long format. It then counts
#' the number of missing values (defined as `NA` or below an `abundance_threshold`)
#' for each feature within each experimental group. If any group exceeds the
#' `max_num_samples_miss_per_group` for a given feature, that feature (row) is
#' removed from the original table.
#'
#' @param input_table A data frame with a feature ID column and abundance columns for samples.
#' @param cols A tidyselect helper (e.g., `starts_with("Abundance")`) to specify the
#'   abundance columns.
#' @param design_matrix A data frame mapping sample IDs to experimental groups.
#' @param sample_id The unquoted column name in `design_matrix` for sample identifiers.
#' @param row_id The unquoted column name in `input_table` for feature identifiers.
#' @param grouping_variable The unquoted column name in `design_matrix` for experimental groups.
#' @param max_num_samples_miss_per_group An integer. The maximum number of samples
#'   allowed to have missing values within a single group for a feature to be retained.
#' @param abundance_threshold A numeric threshold below which values are considered missing.
#' @param temporary_abundance_column A string for the temporary column name used during
#'   pivoting. Defaults to "Abundance".
#'
#' @return A filtered data frame with the same structure as `input_table`.
#' @export
removeRowsWithMissingValues <- function(input_table, cols, design_matrix, sample_id, row_id, grouping_variable, max_num_samples_miss_per_group, abundance_threshold
                                        , temporary_abundance_column = "Abundance") {

  abundance_long <- input_table |>
    pivot_longer(cols = { { cols } },
                 names_to =  as_string(as_name(enquo(sample_id))) ,
                 values_to = temporary_abundance_column  ) |>
    mutate( {{sample_id}} := purrr::map_chr(   {{sample_id}}  , as.character)   ) |>
    left_join(design_matrix |>
                mutate(  {{sample_id}} := purrr::map_chr(    {{sample_id}} , as.character)   )
              , by = as_string(as_name(enquo(sample_id))))

  count_missing_values_per_group <- abundance_long |>
    mutate(is_missing = ifelse(!is.na( !!sym(temporary_abundance_column)) & !!sym(temporary_abundance_column) > abundance_threshold, 0, 1)) |>
    group_by( {{ row_id }}, {{ grouping_variable }} ) |>
    summarise(num_missing_values = sum(is_missing)) |>
    ungroup()

  remove_rows_temp <- count_missing_values_per_group |>
    dplyr::filter(max_num_samples_miss_per_group < num_missing_values) |>
    dplyr::select(-num_missing_values, -{ { grouping_variable } }) |>
    distinct({ { row_id } })

  filtered_tbl <- input_table |>
    dplyr::anti_join(remove_rows_temp, by = as_string(as_name(enquo(row_id))))

  return(filtered_tbl)

}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @title Advanced Filtering Based on Missing Value Percentages
#' @description This function provides advanced, two-tiered filtering of features
#' (rows) based on the percentage of missing values within and across experimental groups.
#'
#' @details The function operates in several steps:
#' 1.  It calculates a dynamic intensity cutoff based on the specified `proteins_intensity_cutoff_percentile` of all observed abundance values. Any value below this is treated as missing.
#' 2.  For each feature and each experimental group, it calculates the percentage of missing values.
#' 3.  It checks if this percentage exceeds the `groupwise_percentage_cutoff`.
#' 4.  It then calculates the percentage of *groups* that failed the check in step 3.
#' 5.  If this percentage of failing groups exceeds `max_groups_percentage_cutoff`, the entire feature is removed from the dataset.
#'
#' This method is useful for removing features that are consistently missing in a significant number of experimental groups.
#'
#' @param input_table A data frame with a feature ID column and abundance columns.
#' @param cols A character vector specifying the abundance column names. This is a key difference from the other helper and requires an explicit vector of names.
#' @param design_matrix A data frame mapping sample IDs to experimental groups.
#' @param sample_id The unquoted column name in `design_matrix` for sample identifiers.
#' @param row_id The unquoted column name in `input_table` for feature identifiers.
#' @param grouping_variable The unquoted column name in `design_matrix` for experimental groups.
#' @param groupwise_percentage_cutoff The maximum allowed percentage of missing values within a single group (0-100). Defaults to 1.
#' @param max_groups_percentage_cutoff The maximum allowed percentage of groups that can fail the `groupwise_percentage_cutoff` (0-100). Defaults to 50.
#' @param proteins_intensity_cutoff_percentile The percentile (0-100) of all intensities used to set a minimum abundance threshold. Defaults to 1.
#' @param temporary_abundance_column A string for the temporary column name used during pivoting. Defaults to "Abundance".
#'
#' @return A filtered data frame with the same structure as `input_table`.
#'
#' @importFrom rlang ensym as_string sym
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate case_when left_join filter pull distinct group_by summarise ungroup anti_join
#' @export
removeRowsWithMissingValuesPercentHelper <- function(input_table
                                                     , cols
                                                     , design_matrix
                                                     , sample_id # symbol, e.g. Run
                                                     , row_id    # symbol
                                                     , grouping_variable # symbol
                                                     , groupwise_percentage_cutoff = 1
                                                     , max_groups_percentage_cutoff = 50
                                                     , proteins_intensity_cutoff_percentile = 1
                                                     , temporary_abundance_column = "Abundance") {

  # Ensure the sample ID column name is a string for robust use
  sample_id_col_name_string <- rlang::as_string(rlang::ensym(sample_id))

  abundance_long <- input_table |>
    tidyr::pivot_longer(cols = !all_of(cols) # Use !all_of() with the cols string directly
                 , names_to = sample_id_col_name_string # Use the explicit string name
                 , values_to = temporary_abundance_column  ) |>
    # Convert the newly created sample ID column to character, using its string name
    dplyr::mutate( !!rlang::sym(sample_id_col_name_string) := as.character(!!rlang::sym(sample_id_col_name_string)) ) |>
    dplyr::mutate( !!sym(temporary_abundance_column) := dplyr::case_when (is.nan(!!sym(temporary_abundance_column)) ~ NA_real_
                                                            , TRUE ~ !!sym(temporary_abundance_column) ) ) |>
    dplyr::left_join(design_matrix |>
                # Convert the sample ID column in design_matrix to character
                # {{sample_id}} correctly refers to the column (e.g., 'Run') in design_matrix here
                dplyr::mutate( {{sample_id}} := as.character( {{sample_id}} ) )
              , by = sample_id_col_name_string ) # Join using the explicit string name

  min_protein_intensity_threshold <- ceiling( quantile( abundance_long |>
                                                          dplyr::filter( !is.nan(!!sym(temporary_abundance_column)) & !is.infinite(!!sym(temporary_abundance_column))) |>
                                                          dplyr::pull(!!sym(temporary_abundance_column))
                                                        , na.rm=TRUE
                                                        , probs = c(proteins_intensity_cutoff_percentile/100) ))[1]

  count_values_per_group <- abundance_long |>
    distinct( {{ sample_id }}, {{ grouping_variable }} ) |>
    group_by(  {{ grouping_variable }} ) |>
    summarise(  num_per_group = n()) |>
    ungroup()

  count_values_missing_per_group <- abundance_long |>
    mutate(is_missing = ifelse( !is.na( !!sym(temporary_abundance_column))
                                & !!sym(temporary_abundance_column) > min_protein_intensity_threshold
                                , 0, 1)) |>
    group_by( {{ row_id }}, {{ grouping_variable }} ) |>
    summarise( num_missing_per_group = sum(is_missing)) |>
    ungroup()

  count_percent_missing_per_group <- count_values_missing_per_group |>
    full_join( count_values_per_group,
               by = join_by( {{ grouping_variable }} )) |>
    mutate(  perc_missing_per_group = num_missing_per_group / num_per_group * 100 )

  total_num_of_groups <- count_values_per_group |> nrow()

  remove_rows_temp <- count_percent_missing_per_group |>
    dplyr::filter(groupwise_percentage_cutoff <  perc_missing_per_group) |>
    group_by( { { row_id } }) |>
    summarise( percent  = n()/total_num_of_groups*100 ) |>
    ungroup() |>
    dplyr::filter(percent > max_groups_percentage_cutoff)

  filtered_tbl <- input_table |>
    dplyr::anti_join(remove_rows_temp, by = join_by({{row_id}}))

  return(filtered_tbl)

}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @title Identify Features to Keep Based on Minimum Valid Samples Per Group
#' @description This function evaluates each feature (row) within each experimental group
#' to determine if it meets a minimum threshold of valid samples. It returns a list
#' of features to keep for each group.
#'
#' @details A value is considered "valid" if it is not `NA` and is above the specified
#' `abundance_threshold`. The function is useful for downstream processing where
#' different sets of features might be used for analyses involving different groups.
#'
#' @param input_table A data frame with a feature ID column and abundance columns.
#' @param cols A tidyselect helper (e.g., `starts_with("Abundance")`) to specify the abundance columns.
#' @param design_matrix A data frame mapping sample IDs to experimental groups.
#' @param sample_id The unquoted column name for sample identifiers.
#' @param row_id The unquoted column name for feature identifiers.
#' @param grouping_variable The unquoted column name for experimental groups.
#' @param min_num_samples_per_group An integer specifying the minimum number of valid samples required per group for a feature to be kept for that group.
#' @param abundance_threshold A numeric threshold below which values are considered missing.
#'
#' @return A named list. Each name corresponds to an experimental group, and each
#'   element is a character vector of feature IDs that passed the filter for that group.
#'
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr left_join mutate group_by summarise ungroup filter select nest
#' @importFrom purrr map
#' @importFrom rlang enquo as_name
#' @export
getRowsToKeepList <- function(input_table, cols, design_matrix, sample_id, row_id, grouping_variable, min_num_samples_per_group, abundance_threshold) {

  abundance_long <- input_table |>
    pivot_longer(cols = { { cols } },
                 names_to = as_string(as_name(enquo(sample_id))),
                 values_to = "Abundance") |>
    left_join(design_matrix, by = as_string(as_name(enquo(sample_id))))


  count_values_per_group <- abundance_long |>
    mutate(has_value = ifelse(!is.na(Abundance) & Abundance > abundance_threshold, 1, 0)) |>
    group_by({ { row_id } }, { { grouping_variable } }) |>
    summarise(num_values = sum(has_value)) |>
    ungroup()


  kept_rows_temp <- count_values_per_group |>
    dplyr::filter(num_values >= min_num_samples_per_group) |>
    dplyr::select(-num_values) |>
    group_by({ { grouping_variable } }) |>
    nest(data = c({ { row_id } })) |>
    ungroup() |>
    mutate(data = purrr::map(data, \(x){ x[, as_name(enquo(row_id))][[1]] }))


  sample_rows_lists <- kept_rows_temp$data
  names(sample_rows_lists) <- kept_rows_temp[, as_name(enquo(grouping_variable))][[1]]

  return(sample_rows_lists)

}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @title Impute Missing Values Using a Down-Shifted Normal Distribution
#' @description This function imputes missing values (NA or non-finite) in a numeric
#' vector by drawing random values from a normal distribution. The distribution is
#' defined relative to the observed data, with a shifted mean and scaled standard deviation.
#' This method is designed to simulate values that are missing-not-at-random (MNAR),
#' such as those below the limit of detection.
#'
#' @param temp A numeric vector containing values to be imputed.
#' @param width A numeric factor to scale the standard deviation of the imputed distribution. A smaller value (e.g., 0.3) creates a narrower distribution, implying less uncertainty. Defaults to 0.3.
#' @param downshift A numeric factor determining how far to shift the mean of the imputed distribution downwards from the observed mean, in units of the observed standard deviation. Defaults to 1.8.
#'
#' @return A numeric vector of the same length as the input, with missing values replaced by imputed values.
#'
#' @importFrom stats sd rnorm
#' @export
imputePerCol <- function(temp, width = 0.3, downshift = 1.8) {

  temp[!is.finite(temp)] <- NA

  temp.sd <- width * sd(temp, na.rm = TRUE)   # shrink sd width
  temp.mean <- mean(temp, na.rm = TRUE) -
    downshift * sd(temp, na.rm = TRUE)   # shift mean of imputed values

  n.missing <- sum(is.na(temp))
  temp[is.na(temp)] <- rnorm(n.missing, mean = temp.mean, sd = temp.sd)
  return(temp)
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @title Create a Replicate Matrix for RUV-III Analysis
#' @description This function transforms a standard design matrix into the specific
#' binary replicate matrix format required by the `RUVIII` function from the `ruv` package.
#'
#' @param design_matrix A data frame containing sample and group information.
#' @param sample_id_column The unquoted column name for unique sample identifiers.
#' @param grouping_variable The unquoted column name for the experimental groups that define the replicates.
#' @param temp_column The name of a temporary column used internally. Defaults to `is_replicate_temp`.
#'
#' @return A numeric matrix where rows are samples, columns are experimental groups,
#'   and a `1` indicates that a sample belongs to the group represented by the column.
#'   This matrix is suitable for the `M` argument of `RUVIII`.
#'
#' @importFrom dplyr select mutate
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames
#' @importFrom rlang enquo as_name as_string
#' @export
getRuvIIIReplicateMatrixHelper <- function(design_matrix, sample_id_column, grouping_variable, temp_column = is_replicate_temp) {

  ruvIII_replicates_matrix <- design_matrix |>
    dplyr::select({ { sample_id_column } }, { { grouping_variable } }) |>
    mutate({ { temp_column } } := 1) |>
    pivot_wider(id_cols = as_string( as_name(enquo(sample_id_column ))),
                names_from = { { grouping_variable } },
                values_from = { { temp_column } },
                values_fill = 0) |>
    column_to_rownames(as_string(as_name(enquo(sample_id_column)))) |>
    as.matrix()

  ruvIII_replicates_matrix
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @title Helper Function for PCA Plotting
#' @description A flexible helper function to perform Principal Component Analysis (PCA)
#' using `mixOmics` and generate a scatter plot using `ggplot2`.
#'
#' @details This function calculates principal components, merges them with sample
#' metadata from the design matrix, and creates a ggplot object. It supports
#' coloring by one variable, shaping by another, and adding text labels.
#'
#' @param data A numeric matrix or data frame where columns are samples and rows are features.
#' @param design_matrix A data frame containing sample metadata.
#' @param sample_id_column The column name in `design_matrix` that contains sample IDs,
#'   matching the column names of `data`. Defaults to "Sample_ID".
#' @param grouping_variable The column name in `design_matrix` to use for coloring points.
#'   Defaults to "group".
#' @param shape_variable Optional. The column name in `design_matrix` to use for point shapes.
#' @param label_column Optional. The column name in `design_matrix` to use for text labels.
#' @param title The title for the plot.
#' @param geom.text.size The font size for labels if `label_column` is provided.
#' @param ncomp The number of principal components to compute. Defaults to 2.
#' @param ... Additional arguments (not currently used).
#'
#' @return A ggplot object representing the PCA plot.
#'
#' @importFrom mixOmics pca
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggplot2 ggplot aes geom_point xlab ylab labs theme coord_cartesian
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join
#' @export
plotPcaHelper <- function(data,
                    design_matrix,
                    sample_id_column = "Sample_ID",
                    grouping_variable = "group",
                    shape_variable = NULL,
                    label_column = NULL,
                    title, geom.text.size = 11, ncomp = 2,
                    ...) {
    
  # Ensure design_matrix is a data frame
  design_matrix <- as.data.frame(design_matrix)
  
  pca.res <- mixOmics::pca(t(as.matrix(data)), ncomp = ncomp)
  proportion_explained <- pca.res$prop_expl_var

  temp_tbl <- pca.res$variates$X |>
    as.data.frame() |>
    rownames_to_column(var = sample_id_column) |>
    left_join(design_matrix, by = sample_id_column)

  # More defensive check for grouping variables
  if (!grouping_variable %in% colnames(temp_tbl)) {
    stop(sprintf("Grouping variable '%s' not found in the data", grouping_variable))
  }
  
  if (!is.null(shape_variable) && !shape_variable %in% colnames(temp_tbl)) {
    stop(sprintf("Shape variable '%s' not found in the data", shape_variable))
  }

  # Create base plot with appropriate aesthetics based on whether shape_variable is NULL
  if (is.null(label_column) || label_column == "") {
    if (is.null(shape_variable)) {
      # No shape variation, only color
      base_plot <- temp_tbl |>
        ggplot(aes(PC1, PC2, color = !!sym(grouping_variable)))
    } else {
      # Both color and shape
      base_plot <- temp_tbl |>
        ggplot(aes(PC1, PC2, color = !!sym(grouping_variable), shape = !!sym(shape_variable)))
    }
  } else {
    if (!label_column %in% colnames(temp_tbl)) {
      stop(sprintf("Label column '%s' not found in the data", label_column))
    }
    
    if (is.null(shape_variable)) {
      # No shape variation, only color, with labels
      base_plot <- temp_tbl |>
        ggplot(aes(PC1, PC2, color = !!sym(grouping_variable), label = !!sym(label_column)))
    } else {
      # Both color and shape, with labels
      base_plot <- temp_tbl |>
        ggplot(aes(PC1, PC2, color = !!sym(grouping_variable), shape = !!sym(shape_variable), 
                   label = !!sym(label_column)))
    }
  }

  output <- base_plot +
    geom_point(size = 3) +
    xlab(paste("PC1 (", round(proportion_explained$X[["PC1"]] * 100, 0), "%)", sep = "")) +
    ylab(paste("PC2 (", round(proportion_explained$X[["PC2"]] * 100, 0), "%)", sep = "")) +
    labs(title = title) +
    theme(legend.title = element_blank())

  # Calculate axis limits based on the full range of data
  pc1_range <- range(temp_tbl$PC1, na.rm = TRUE)
  pc2_range <- range(temp_tbl$PC2, na.rm = TRUE)
  buffer_pc1 <- (pc1_range[2] - pc1_range[1]) * 0.05 # 5% buffer
  buffer_pc2 <- (pc2_range[2] - pc2_range[1]) * 0.05 # 5% buffer

  output <- output + coord_cartesian(
    xlim = c(pc1_range[1] - buffer_pc1, pc1_range[2] + buffer_pc1),
    ylim = c(pc2_range[1] - buffer_pc2, pc2_range[2] + buffer_pc2)
  )
  
  if (!is.null(label_column) && label_column != "") {
    output <- output + geom_text_repel(size = geom.text.size, show.legend = FALSE)
  }

  output
}



#' @title Plot a List of PCA Plots
#' @description This function takes a dataset and generates a list of PCA plots,
#' where each plot is colored by a different grouping variable specified in a list.
#'
#' @param data A numeric matrix or data frame (features x samples).
#' @param design_matrix A data frame with sample metadata.
#' @param sample_id_column The column name for sample IDs.
#' @param grouping_variables_list A character vector of column names in `design_matrix`
#'   to use for coloring the points in each respective plot.
#' @param label_column Optional. The column name for text labels on the plots.
#' @param title The base title for the plots.
#' @param geom.text.size The font size for labels.
#' @param ncomp The number of principal components to compute.
#' @param ... Additional arguments (not currently used).
#'
#' @return A named list of ggplot objects, where each element is a PCA plot colored
#'   by one of the `grouping_variables_list`.
#'
#' @importFrom purrr map
#' @export
plotPcaListHelper <- function(data,
                              design_matrix,
                              sample_id_column = "Sample_ID",
                              grouping_variables_list = c("group"),
                              label_column = NULL,
                              title, geom.text.size = 11, ncomp = 2,
                              ...) {

  pca.res <- mixOmics::pca(t(as.matrix(data)), ncomp = ncomp)
  proportion_explained <- pca.res$prop_expl_var

  temp_tbl <- pca.res$variates$X |>
    as.data.frame() |>
    rownames_to_column(var = sample_id_column) |>
    left_join(design_matrix, by = sample_id_column)


  plotOneGgplotPca <- function( grouping_variable ) {
    unique_groups <- temp_tbl |> distinct(!!sym(grouping_variable)) |> dplyr::pull(!!sym(grouping_variable))

    if (is.null(label_column) || label_column == "") {
      output <- temp_tbl |>
        ggplot(aes(PC1, PC2, col = !!sym(grouping_variable))) +
        geom_point() +
        xlab(paste("PC1 (", round(proportion_explained$X[["PC1"]] * 100, 0), "%)", sep = "")) +
        ylab(paste("PC2 (", round(proportion_explained$X[["PC2"]] * 100, 0), "%)", sep = "")) +
        labs(title = title) +
        theme(legend.title = element_blank())
    } else {
      output <- temp_tbl |>
        ggplot(aes(PC1, PC2, col = !!sym(grouping_variable), label = !!sym(label_column))) +
        geom_point() +
        geom_text_repel(size = geom.text.size, show.legend = FALSE) +
        xlab(paste("PC1 (", round(proportion_explained$X[["PC1"]] * 100, 0), "%)", sep = "")) +
        ylab(paste("PC2 (", round(proportion_explained$X[["PC2"]] * 100, 0), "%)", sep = "")) +
        labs(title = title) +
        theme(legend.title = element_blank())
    }

    return(output)
  }

  output_list <- purrr::map( grouping_variables_list, plotOneGgplotPca)

  output_list
}





#' @title Create a PCA-based GGPAIRS Plot
#' @description Performs PCA on a data matrix and visualizes the relationships between
#' the principal components using a `GGally::ggpairs` plot. This creates a matrix
#' of plots, including scatter plots, density plots, and correlations for the PCs.
#'
#' @param data_matrix A numeric matrix or data frame (features x samples).
#' @param design_matrix A data frame with sample metadata.
#' @param grouping_variable The unquoted column name in `design_matrix` to use for coloring points.
#' @param sample_id_column The unquoted column name for sample identifiers.
#' @param ncomp The number of principal components to compute and plot.
#'
#' @return A `ggpairs` object.
#'
#' @importFrom mixOmics pca
#' @importFrom purrr map_chr
#' @importFrom tibble as.data.frame rownames_to_column
#' @importFrom dplyr left_join
#' @importFrom GGally ggpairs
#' @export
plotPcaGgpairs <- function( data_matrix
                            , design_matrix
                            , grouping_variable
                            , sample_id_column
                            , ncomp=2 ) {

  pca.res <- mixOmics::pca(t(as.matrix(data_matrix)), ncomp=ncomp)


  pca_prop_explained_helper <- function( pca_obj, comp_idx ) {
    proportion_explained <- pca.res$prop_expl_var

    pc_label <- paste0("PC", comp_idx)

    perc_label <- paste( paste0(pc_label, " ("), round(proportion_explained$X[[pc_label]]*100, 0),"%)", sep="")

    perc_label
  }

  pc_list <- purrr::map_chr( seq_len(ncomp)
                             , \(comp_idx){ pca_prop_explained_helper(pca.res, comp_idx)})

  pca_variates_x <- pca.res$variates$X

  colnames(pca_variates_x) <- pc_list

  pca_plot_ggpairs <- pca_variates_x |>
    as.data.frame() |>
    rownames_to_column( sample_id_column ) |>
    left_join( design_matrix
               , by = join_by( !!sym(sample_id_column) ==  !!sym(sample_id_column)) ) |>
    ggpairs( columns=pc_list, aes( colour = !!sym(grouping_variable), fill= !!sym(grouping_variable), alpha=0.4)
             , legend = 1)

  pca_plot_ggpairs
}

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


#' @title Plot Relative Log Expression (RLE)
#' @description Generates a Relative Log Expression (RLE) plot for a data matrix.
#' RLE plots are useful for visualizing unwanted variation in high-throughput data.
#'
#' @details This function calculates the RLE for each sample by subtracting the
#' median abundance of each feature from all samples. It then creates a boxplot
#' for each sample, showing the distribution of these relative log expressions.
#' The `rowinfo` parameter can be used to color the boxplots by experimental group.
#'
#' @param Y A numeric matrix where rows are samples and columns are features. Note
#'   that the input is typically a transposed abundance matrix.
#' @param rowinfo An optional vector or single-column data frame containing metadata
#'   (e.g., group assignments) for each sample (row) in `Y`.
#' @param probs A numeric vector of probabilities for calculating quantiles for the boxplot.
#' @param yaxis_limit A numeric vector of length 2 specifying the y-axis limits.
#'
#' @return A ggplot object representing the RLE plot.
#'
#' @import ggplot2
#' @importFrom dplyr arrange distinct pull mutate
#' @export
plotRleHelper <- function(Y, rowinfo = NULL, probs = c(0.05, 0.25, 0.5, 0.75,
                                                       0.95), yaxis_limit = c(-0.5, 0.5))
{
  #  checks = check.ggplot()
  # if (checks) {
  rle <- t(apply(t(Y) - apply(Y, 2, function(x){median(x, na.rm=TRUE)}), 2, function(x){quantile(x, probs = probs, na.rm=TRUE)}))
  colnames(rle) <- c("min", "lower", "middle", "upper",
                     "max")
  df <- cbind(data.frame(rle.x.factor = rownames(rle)), data.frame(rle))

  if (!is.null(rowinfo)) {
    rowinfo <- data.frame(rowinfo = rowinfo)
    df_temp <- cbind(df, rowinfo)

    my.x.factor.levels <- df_temp |>
      arrange(rowinfo) |>
      distinct(rle.x.factor) |>
      dplyr::pull(rle.x.factor)

    df <- df_temp |>
      mutate(rle.x.factor = factor(rle.x.factor,
                                   levels = my.x.factor.levels)) |>
      arrange(rowinfo)
  }

  rleplot <- ggplot(df, aes(x = .data[["rle.x.factor"]])) +
    geom_boxplot(aes(lower = .data[["lower"]]
                     , middle = .data[["middle"]]
                     , upper = .data[["upper"]]
                     , max = .data[["max"]]
                     , min = .data[["min"]]),
                 stat = "identity") +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90) #, axis.ticks.x = element_blank()
    ) +
    theme(axis.title.y = element_blank(), axis.text.y = element_text(size = rel(1.5))) +
    geom_hline(yintercept = 0)


  if( length( yaxis_limit ) ==2 ) {

    rleplot <- rleplot +
      coord_cartesian(ylim = yaxis_limit)

  }


  if (!is.null(rowinfo)) {
    if (ncol(rowinfo) == 1) {
      rleplot <- rleplot + aes(fill = rowinfo) + labs(fill = "")
    }
  }

  return(rleplot)
  # }
  # else return(FALSE)
}


#' @title Get Min/Max Y-Axis Limits from a Boxplot
#' @description A helper function to extract the minimum and maximum data points
#' from a ggplot boxplot object and apply an adjustment factor. This is useful for
#' setting coordinated y-axis limits across multiple plots.
#'
#' @param input_boxplot A ggplot object containing a `geom_boxplot` layer where the
#'   data has `min` and `max` columns (as created by `plotRleHelper`).
#' @param adjust_factor A numeric value to expand the range by. For example, 0.05
#'   will expand the range by 5% at both ends.
#'
#' @return A numeric vector of length 2 containing the adjusted minimum and maximum values.
#' @export
getMaxMinBoxplot <- function( input_boxplot, adjust_factor = 0.05) {

  df_min <- min( input_boxplot$data$min, na.rm=TRUE)

  df_max <- max( input_boxplot$data$max, na.rm=TRUE )

  if( df_min > 0 ) {
    df_min <- df_min*(1-adjust_factor)
  } else {
    df_min <- df_min*(1+adjust_factor)

  }

  if( df_max > 0 ) {
    df_max <- df_max*(1+adjust_factor)
  } else {
    df_max <- df_max*(1-adjust_factor)

  }

  return( c(df_min, df_max))
}

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' @title Create a Grid of RLE and PCA Plots
#' @description This function generates and arranges a series of RLE and PCA plots
#' for different datasets or analysis steps (e.g., before and after normalization).
#'
#' @param list_of_data_matrix A list of data matrices to be plotted.
#' @param list_of_design_matrix A list of corresponding design matrices.
#' @param sample_id_column The unquoted column name for sample identifiers.
#' @param grouping_variable The unquoted column name for the grouping variable used for coloring.
#' @param list_of_descriptions A character vector of titles for each corresponding plot pair.
#'
#' @return A `ggarrange` object containing the grid of plots.
#'
#' @importFrom purrr pmap
#' @importFrom ggpubr ggarrange
#' @importFrom rlang enquo as_name
#' @export
rlePcaPlotList <- function(list_of_data_matrix, list_of_design_matrix,
                           sample_id_column = Sample_ID, grouping_variable = group, list_of_descriptions) {

  rle_list <- purrr::pmap( list( data_matrix=list_of_data_matrix, description=list_of_descriptions, design_matrix=list_of_design_matrix),
                           function( data_matrix, description, design_matrix) { plotRleHelper(t(as.matrix(data_matrix)),
                                                                                              rowinfo = design_matrix[colnames(data_matrix), as_name(enquo(grouping_variable))]  )  +
                               labs(title = description)} )

  pca_list <- purrr::pmap(list( data_matrix=list_of_data_matrix, description=list_of_descriptions, design_matrix=list_of_design_matrix),
                          function( data_matrix, description, design_matrix) { plotPcaHelper(data_matrix,
                                                                                             design_matrix = design_matrix,
                                                                                             sample_id_column = sample_id_column ,
                                                                                             grouping_variable =  grouping_variable ,
                                                                                             title = description, cex = 7) })

  list_of_plots <- c(rle_list, pca_list)

  rle_pca_plots_arranged <- ggarrange(plotlist = list_of_plots, nrow = 2, ncol = length(list_of_descriptions),
                                      common.legend = FALSE, legend = "bottom", widths = 10, heights = 10)

  rle_pca_plots_arranged
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' @title Count Significantly Differentially Expressed Genes/Proteins
#' @description This function takes a data frame of differential expression results
#' and counts the number of features that are "Significant and Up", "Significant and Down",
#' or "Not significant" based on log-fold-change and q-value thresholds.
#'
#' @param data A data frame containing differential expression results.
#' @param lfc_thresh The absolute log-fold-change threshold for significance. Defaults to 0.
#' @param q_val_thresh The q-value (adjusted p-value) threshold for significance. Defaults to 0.05.
#' @param log_fc_column The unquoted column name for the log-fold-change values.
#' @param q_value_column The unquoted column name for the q-values.
#'
#' @return A data frame with two columns: `status` and `counts`.
#'
#' @importFrom dplyr mutate case_when group_by summarise ungroup left_join
#' @export
countStatDeGenes <- function(data,
                             lfc_thresh = 0,
                             q_val_thresh = 0.05,
                             log_fc_column = log2FC,
                             q_value_column = fdr_qvalue) {

  # comparison <- as.data.frame(data) |>
  #   distinct(comparison) |>
  #   pull(comparison)

  selected_data <- data |>
    dplyr::mutate(status = case_when( { { q_value_column } } >= q_val_thresh ~ "Not significant"
                                      , { { log_fc_column } } >= lfc_thresh & { { q_value_column } } < q_val_thresh ~ "Significant and Up"
                                      , { { log_fc_column } } < lfc_thresh & { { q_value_column } } < q_val_thresh ~ "Significant and Down"
                                      , TRUE ~ "Not significant"))

  counts <- selected_data |>
    group_by(status) |>
    summarise(counts = n()) |>
    ungroup()

  all_possible_status <- data.frame( status = c( "Not significant", "Significant and Up", "Significant and Down"))

  results <- all_possible_status |>
    left_join( counts, by=c("status" = "status")) |>
    mutate( counts = ifelse(is.na(counts), 0, counts))

  return(results)
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @title Count DE Genes/Proteins for a List of Analyses (Helper)
#' @description A helper function that applies `countStatDeGenes` to a list of
#' differential expression result tables. It then combines the counts into a single
#' data frame, adding columns for the analysis description and comparison details.
#'
#' @param de_table A named list of data frames, where each data frame is a DE result table.
#' @param description A string describing the analysis type (e.g., "RUV-corrected").
#' @param facet_column The unquoted column name to store the `description`.
#' @param comparison_column The unquoted column name to store the comparison name.
#' @param expression_column The unquoted column name to store the expression part of the contrast.
#'
#' @return A single data frame summarizing the counts of significant features across all comparisons.
#'
#' @importFrom purrr map map2
#' @importFrom dplyr bind_rows mutate
#' @importFrom tidyr separate_wider_delim
#' @export
countStatDeGenesHelper <- function(de_table
                                   , description
                                   , facet_column = analysis_type
                                   , comparison_column = "comparison"
                                   , expression_column = "expression") {

  # print(head(de_table))

  de_table_updated <- purrr::map(de_table, \(x){  countStatDeGenes(x,
                                                                   lfc_thresh = 0,
                                                                   q_val_thresh = 0.05,
                                                                   log_fc_column = logFC,
                                                                   q_value_column = fdr_qvalue)})

  list_of_tables <- purrr::map2(de_table_updated
                                ,names(de_table_updated)
                                ,\(.x, .y){ .x |>
                                    mutate(!!sym(comparison_column) := .y) })

  # print(head(temp[[1]]))

  merged_tables <- list_of_tables |>
    bind_rows() |>
    mutate({ { facet_column } } := description) |>
    separate_wider_delim( !!sym(comparison_column ),
                          delim = "=",
                          names = c(comparison_column,
                                    expression_column))

  merged_tables
}

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @title Generate a Bar Plot of DE Gene/Protein Counts
#' @description This function processes a list of differential expression results from
#' multiple analysis types (e.g., "raw", "normalized"), counts the number of
#' significant features for each, and generates a faceted bar plot summarizing these counts.
#'
#' @param list_of_de_tables A named list where each element is a list of DE result tables
#'   (e.g., `list(raw = list_of_raw_de_tables, normalized = ...)`).
#' @param list_of_descriptions A character vector of descriptions corresponding to `list_of_de_tables`.
#' @param formula_string A formula string (e.g., `"analysis_type ~ comparison"`) for faceting the plot.
#' @param facet_column The unquoted column name for the analysis type facet.
#' @param comparison_column The unquoted column name for the comparison facet.
#' @param expression_column The unquoted column name for the expression part of the contrast.
#'
#' @return A list containing two elements: `plot` (the ggplot object) and `table` (the summary data frame).
#'
#' @importFrom purrr map2
#' @importFrom dplyr bind_rows filter
#' @importFrom ggplot2 ggplot aes geom_bar geom_text theme element_text facet_grid
#' @export
printCountDeGenesTable <- function(  list_of_de_tables
                                     , list_of_descriptions
                                     , formula_string = "analysis_type ~ comparison"
                                     , facet_column = analysis_type
                                     , comparison_column = "comparison"
                                     , expression_column = "expression") {




  num_significant_de_genes_all <- purrr::map2(list_of_de_tables,
                                              list_of_descriptions,
                                              function(a, b) { countStatDeGenesHelper( de_table = a
                                                                                       , description = b
                                                                                       , facet_column = {{facet_column}}
                                                                                       , comparison_column = comparison_column
                                                                                       , expression_column = expression_column) }) |>
    bind_rows()

  num_sig_de_genes_barplot <- num_significant_de_genes_all |>
    dplyr::filter(status != "Not significant") |>
    ggplot(aes(x = status, y = counts)) +
    geom_bar(stat = "identity") +
    geom_text(stat = 'identity', aes(label = counts), vjust = -0.5) +
    theme(axis.text.x = element_text(angle = 90))

  # print(head(num_sig_de_genes_barplot))

  if (!is.na(formula_string)) {

    num_sig_de_genes_barplot <- num_sig_de_genes_barplot +
      facet_grid(as.formula(formula_string))
  }


  return(list(plot = num_sig_de_genes_barplot, table = num_significant_de_genes_all))
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @title Prepare Data for Volcano Plots
#' @description This function takes a list of differential expression results, combines
#' them into a single long-format data frame, and prepares them for plotting with
#' `plotVolcano`. It calculates -log10(q-value) and assigns a color category
#' based on significance thresholds.
#'
#' @param list_of_de_tables A list of DE result lists.
#' @param list_of_descriptions A character vector of descriptions for each list in `list_of_de_tables`.
#' @param row_id The unquoted column name for feature identifiers.
#' @param p_value_column The unquoted column name for raw p-values.
#' @param q_value_column The unquoted column name for q-values.
#' @param fdr_value_column The unquoted column name for another FDR-adjusted p-value column.
#' @param log_q_value_column The unquoted column name to store the -log10(q-value).
#' @param log_fc_column The unquoted column name for log-fold-change values.
#' @param comparison_column The unquoted column name for the comparison.
#' @param expression_column The unquoted column name for the expression part of the contrast.
#' @param facet_column The unquoted column name for the analysis type facet.
#' @param q_val_thresh The q-value threshold for determining significance.
#'
#' @return A single data frame ready for volcano plotting, with columns for plotting
#'   aesthetics (`lqm`, `colour`, etc.).
#'
#' @importFrom purrr map map2
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr select mutate bind_rows
#' @importFrom tidyr separate_wider_delim
#' @importFrom rlang enquo as_name as_string
#' @export
getSignificantData <- function( list_of_de_tables
                                , list_of_descriptions
                                , row_id = uniprot_acc
                                , p_value_column = raw_pvalue
                                , q_value_column = fdr_qvalue
                                , fdr_value_column = fdr_value_bh_adjustment
                                , log_q_value_column = lqm
                                , log_fc_column = logFC
                               , comparison_column = "comparison"
                               , expression_column = "log_intensity"
                               , facet_column = analysis_type
                               , q_val_thresh = 0.05) {

  get_row_binded_table <- function(de_table_list, description) {
    output <- purrr::map(de_table_list,
                         function(tbl) { tbl |>
                           rownames_to_column(as_string(as_name(enquo(row_id)))) |>
                           dplyr::select({ { row_id } },
                                         { { p_value_column } },
                                         { { q_value_column } },
                                         { { fdr_value_column } },
                                         { { log_fc_column } }) }) |>
      purrr::map2(names(de_table_list), \(.x, .y){ .x |>
        mutate({ { comparison_column } } := .y) }) |>
      bind_rows() |>
      mutate({ { facet_column } } := description) |>
      separate_wider_delim({ { comparison_column } },
               delim = "=",
               names = c( comparison_column,
                          expression_column) )

  }

  logfc_tbl_all <- purrr::map2(list_of_de_tables, list_of_descriptions,
                               function(a, b) { get_row_binded_table(de_table_list = a, description = b) }) |>
    bind_rows()

  selected_data <- logfc_tbl_all |>
    mutate({ { log_q_value_column } } := -log10(fdr_qvalue)) |>
    dplyr::select({ { row_id } }, { { log_q_value_column } }, { { q_value_column } }, { { p_value_column } }, { { log_fc_column } },
                  { { comparison_column } }, { { expression_column } },
                  { { facet_column } }) |>
    dplyr::mutate(colour = case_when(abs({ { log_fc_column } }) >= 1 & { { q_value_column } } >= q_val_thresh ~ "orange",
                                     abs({ { log_fc_column } }) >= 1 & { { q_value_column } } < q_val_thresh ~ "purple",
                                     abs({ { log_fc_column } }) < 1 & { { q_value_column } } < q_val_thresh ~ "blue",
                                     TRUE ~ "black")) |>
    dplyr::mutate(colour = factor(colour, levels = c("black", "orange", "blue", "purple")))

  selected_data

}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' @title Draw a Faceted Volcano Plot
#' @description Generates a volcano plot from pre-processed differential expression
#' data. The plot displays -log10(q-value) vs. log-fold-change and can be faceted
#' to show multiple comparisons or analysis types.
#'
#' @param selected_data A data frame prepared by `getSignificantData`, containing
#'   columns for plotting aesthetics (`log_q_value_column`, `log_fc_column`, `colour`).
#' @param log_q_value_column The unquoted column name for the -log10(q-value) y-axis.
#' @param log_fc_column The unquoted column name for the log-fold-change x-axis.
#' @param q_val_thresh A numeric threshold for q-value significance, used to draw a horizontal line.
#' @param formula_string A formula string (e.g., `"analysis_type ~ comparison"`) for faceting with `facet_grid`.
#'
#' @return A ggplot object representing the volcano plot.
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_colour_manual geom_vline geom_hline theme_bw xlab ylab theme facet_grid
#' @export
plotVolcano <- function( selected_data
                         , log_q_value_column = lqm
                         , log_fc_column = logFC
                         , q_val_thresh = 0.05
                         , formula_string = "analysis_type ~ comparison" ) {

  volplot_gg.all <- selected_data |>
    ggplot(aes(y = { { log_q_value_column } }, x = { { log_fc_column } })) +
    geom_point(aes(col = colour)) +
    scale_colour_manual(values = c(levels(selected_data$colour)),
                        labels = c(paste0("Not significant, logFC > ",
                                          1),
                                   paste0("Significant, logFC >= ",
                                          1),
                                   paste0("Significant, logFC <",
                                          1),
                                   "Not Significant")) +
    geom_vline(xintercept = 1, colour = "black", linewidth = 0.2) +
    geom_vline(xintercept = -1, colour = "black", linewidth = 0.2) +
    geom_hline(yintercept = -log10(q_val_thresh)) +
    theme_bw() +
    xlab("Log fold changes") +
    ylab("-log10 q-value") +
    theme(legend.position = "none")

  volplot_gg.plot <- volplot_gg.all
  if( !is.na(formula_string) | formula_string != "" ) {
    volplot_gg.plot <- volplot_gg.all +
      facet_grid( as.formula(formula_string),
                  labeller = labeller(facet_category = label_wrap_gen(width = 10)))
  }

  volplot_gg.plot
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @title Draw a Single Volcano Plot for Publication
#' @description Creates a single, publication-quality volcano plot with more control
#' over aesthetics.
#'
#' @param input_data A data frame containing the data for the plot.
#' @param input_title The title of the plot.
#' @param log_q_value_column The unquoted column name for the y-axis (-log10 q-value).
#' @param log_fc_column The unquoted column name for the x-axis (log2 fold-change).
#' @param points_type_label The unquoted column name containing labels for the color legend.
#' @param points_color The unquoted column name containing hex codes or color names for points.
#' @param q_val_thresh The q-value threshold for the significance line.
#' @param log2FC_thresh The log2-fold-change threshold for vertical lines.
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_colour_manual geom_hline geom_vline theme_bw xlab ylab labs theme element_text
#' @importFrom dplyr distinct pull
#' @export
plotOneVolcano <- function( input_data, input_title,
                            log_q_value_column = lqm,
                            log_fc_column = logFC,
                            points_type_label = label,
                            points_color = colour,
                            q_val_thresh=0.05 ,
                            log2FC_thresh = 1) {

  colour_tbl <- input_data |>
    distinct( {{points_type_label}}, {{points_color}} )

  # print(colour_tbl)

  colour_map <- colour_tbl |>
    dplyr::pull({{points_color}} ) |>
    as.vector()

  names( colour_map ) <- colour_tbl |>
    dplyr::pull({{points_type_label}} )

  avail_labels <- input_data |>
    distinct({{points_type_label}}) |>
    dplyr::pull({{points_type_label}})

  avail_colours <- colour_map[avail_labels]

  # print(avail_labels)
  # print(avail_colours)

  volcano_plot <-  input_data |>
    ggplot(aes(y = {{log_q_value_column}},
               x = {{log_fc_column}} )) +
    geom_point(aes(col = label)) +
    scale_colour_manual(values = avail_colours)  +
    geom_hline(yintercept = -log10(q_val_thresh)) +
    theme_bw() +
    xlab(expression(Log[2](`fold-change`))) +
    ylab(expression(-log[10](`q-value`))) +
    labs(title = input_title)+  # Remove legend title
    theme(legend.title = element_blank()) +
    # theme(legend.position = "none")  +
    theme(axis.text.x = element_text(size = 13))   +
    theme(axis.text.y = element_text(size = 13))  +
    theme(axis.title.x = element_text(size = 12))  +
    theme(axis.title.y = element_text(size = 12))  +
    theme(plot.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 12)) # +
  # theme(legend.title = element_text(size = 12))

  if( !is.na(log2FC_thresh) ) {
  volcano_plot <- volcano_plot+
    geom_vline(xintercept = log2FC_thresh, colour = "black", linewidth = 0.2) +
    geom_vline(xintercept = -log2FC_thresh, colour = "black", linewidth = 0.2)

  }

  volcano_plot
}

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' @title Prepare Data for a Labeled Volcano Plot
#' @description Processes a differential expression results table to prepare it for
#' creating a volcano plot with labels for the most significant points. It joins
#' with UniProt data to get gene names and ranks features to select top candidates for labeling.
#'
#' @param input_table The DE results data frame.
#' @param protein_id_column The unquoted column name for protein identifiers.
#' @param uniprot_table A data frame with UniProt annotations.
#' @param uniprot_protein_id_column The unquoted column name for protein IDs in `uniprot_table`.
#' @param gene_name_column The unquoted column name for gene names in `uniprot_table`.
#' @param number_of_genes The number of top up- and down-regulated genes to label.
#' @param fdr_threshold The FDR/q-value significance threshold.
#' @param fdr_column The unquoted column name for FDR/q-values.
#' @param log2FC_column The unquoted column name for log2 fold changes.
#'
#' @return A data frame with added columns (`lqm`, `colour`, `label`, `gene_name_significant`, etc.)
#'   ready for plotting.
#'
#' @importFrom dplyr mutate relocate select left_join case_when arrange
#' @importFrom purrr map_chr
#' @importFrom stringr str_split
#' @importFrom rlang enquo as_name as_string
#' @export
prepareDataForVolcanoPlot <- function(input_table
                                      , protein_id_column = uniprot_acc
                                      , uniprot_table
                                      , uniprot_protein_id_column = uniprot_acc_first
                                      , gene_name_column = gene_name
                                      , number_of_genes = 3000
                                      , fdr_threshold = 0.05
                                      , fdr_column = q.mod
                                      , log2FC_column = log2FC){

  temp_col_name <-  as_string(as_name(enquo(protein_id_column)))

  proteomics_volcano_tbl <- input_table |>
    dplyr::mutate( uniprot_acc_first  = purrr::map_chr( {{protein_id_column}}, \(x) { str_split( x, ":")[[1]][1] } ) ) |>
    dplyr::relocate( uniprot_acc_first, .after=temp_col_name) |>
    dplyr::select (uniprot_acc_first, {{fdr_column}}, {{log2FC_column}}) |>
    left_join(uniprot_table
              , by = join_by( uniprot_acc_first == {{uniprot_protein_id_column}})) |>
    mutate( colour = case_when ( {{fdr_column}} < fdr_threshold & {{log2FC_column}} > 0 ~ "red"
                                 , {{fdr_column}} < fdr_threshold & {{log2FC_column}} < 0 ~ "blue"
                                 , TRUE ~ "grey" )) |>
    mutate( lqm = -log10({{fdr_column}})) |>
    mutate( label = case_when (  {{fdr_column}} < fdr_threshold & {{log2FC_column}} > 0 ~ "Significant Increase"
                                 , {{fdr_column}} < fdr_threshold & {{log2FC_column}} < 0 ~ "Significant Decrease"
                                 , TRUE ~ "Not significant" )) |>
    mutate( label = factor( label, levels = c( "Significant Increase"
                                               , "Significant Decrease"
                                               , "Not significant" ))) |>
    mutate ( rank_positive = case_when( {{log2FC_column}} > 0 ~ {{fdr_column}}
                                        , TRUE ~ NA_real_) |> rank() ) |>
    mutate ( rank_negative = case_when( {{log2FC_column}} < 0 ~ {{fdr_column}}
                                        , TRUE ~ NA_real_) |> rank() ) |>
    mutate ( {{gene_name_column}} := purrr::map_chr( {{gene_name_column}}, \(x) { str_split( x, " ")[[1]][1] } ) ) |>
    mutate( gene_name_significant = case_when( {{fdr_column}} < fdr_threshold &
                                                 ( rank_positive <= number_of_genes |
                                                   rank_negative <= number_of_genes ) ~  {{gene_name_column}}
                                               , TRUE ~ NA ) )

  proteomics_volcano_tbl
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @title Draw a Volcano Plot Without Vertical Fold-Change Lines
#' @description A variation of `plotOneVolcano` that creates a volcano plot without
#' the vertical lines indicating a fold-change threshold.
#'
#' @param input_data A data frame prepared for plotting.
#' @param input_title The title for the plot.
#' @param log_q_value_column The unquoted column name for the y-axis.
#' @param log_fc_column The unquoted column name for the x-axis.
#' @param points_type_label The unquoted column name for legend labels.
#' @param points_color The unquoted column name for point colors.
#' @param q_val_thresh The q-value threshold for the significance line.
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_colour_manual geom_hline theme_bw xlab ylab labs theme element_text
#' @importFrom dplyr distinct pull
#' @export
plotOneVolcanoNoVerticalLines <- function( input_data, input_title,
                            log_q_value_column = lqm,
                            log_fc_column = logFC,
                            points_type_label = label,
                            points_color = colour,
                            q_val_thresh=0.05) {

  colour_tbl <- input_data |>
    distinct( {{points_type_label}}, {{points_color}} )

  colour_map <- colour_tbl |>
    dplyr::pull({{points_color}} ) |>
    as.vector()

  names( colour_map ) <- colour_tbl |>
    dplyr::pull({{points_type_label}} )

  avail_labels <- input_data |>
    distinct({{points_type_label}}) |>
    dplyr::pull({{points_type_label}})

  avail_colours <- colour_map[avail_labels]

  volcano_plot <-  input_data |>
    ggplot(aes(y = {{log_q_value_column}},
               x = {{log_fc_column}},
               col={{points_type_label}})) +
    geom_point()

  volcano_plot <-   volcano_plot +
    scale_colour_manual(values = avail_colours) +
    # geom_vline(xintercept = 1, colour = "black", size = 0.2) +
    # geom_vline(xintercept = -1, colour = "black", size = 0.2) +
    geom_hline(yintercept = -log10(q_val_thresh)) +
    theme_bw() +
    xlab(expression(Log[2](`fold-change`))) +
    ylab(expression(-log[10](FDR))) +
    labs(title = input_title)+  # Remove legend title
    theme(legend.title = element_blank()) +
    # theme(legend.position = "none")  +
    theme(axis.text.x = element_text(size = 13))   +
    theme(axis.text.y = element_text(size = 13))  +
    theme(axis.title.x = element_text(size = 12))  +
    theme(axis.title.y = element_text(size = 12))  +
    theme(plot.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 12)) # +
  # theme(legend.title = element_text(size = 12))

  volcano_plot
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' @title Create a Labeled Volcano Plot
#' @description A wrapper function that combines data preparation and plotting to
#' create a volcano plot where the most significant proteins are labeled with gene names.
#'
#' @param input_table The raw DE results table.
#' @param uniprot_table A data frame with UniProt annotations for mapping gene names.
#' @param protein_id_column Unquoted column name for protein IDs in `input_table`.
#' @param uniprot_protein_id_column Unquoted column name for protein IDs in `uniprot_table`.
#' @param gene_name_column Unquoted column name for gene names in `uniprot_table`.
#' @param number_of_genes The number of top features to label.
#' @param fdr_threshold FDR/q-value significance threshold.
#' @param fdr_column Unquoted column name for FDR/q-values.
#' @param log2FC_column Unquoted column name for log2 fold changes.
#' @param input_title The plot title.
#' @param include_protein_label A logical flag to enable/disable labeling.
#' @param max.overlaps The `max.overlaps` parameter passed to `ggrepel::geom_text_repel`.
#'
#' @return A ggplot object of the labeled volcano plot.
#'
#' @importFrom ggrepel geom_text_repel
#' @export
printOneVolcanoPlotWithProteinLabel <- function( input_table
                                                 , uniprot_table
                                                 , protein_id_column = Protein.Ids
                                                 , uniprot_protein_id_column = Entry
                                                 , gene_name_column = gene_name
                                                 , number_of_genes = 100
                                                 , fdr_threshold = 0.05
                                                 , fdr_column = fdr_qvalue
                                                 , log2FC_column = log2FC
                                                 , input_title = "Proteomics"
                                                 , include_protein_label = TRUE
                                                 , max.overlaps = 20) {
  proteomics_volcano_tbl <- prepareDataForVolcanoPlot( input_table
                                                       , protein_id_column = {{protein_id_column}}
                                                       , uniprot_table = uniprot_table
                                                       , uniprot_protein_id_column = {{uniprot_protein_id_column}}
                                                       , gene_name_column = {{gene_name_column}}
                                                       , number_of_genes = number_of_genes
                                                       , fdr_threshold = fdr_threshold
                                                       , fdr_column = {{fdr_column}}
                                                       , log2FC_column = {{log2FC_column}})

  proteomics_volcano_plot <- plotOneVolcanoNoVerticalLines ( proteomics_volcano_tbl,
                                                                                 input_title = input_title,
                                                                                 log_q_value_column = lqm,
                                                                                 log_fc_column = log2FC,
                                                                                 points_type_label = label,
                                                                                 points_color = colour,
                                                                                 q_val_thresh= fdr_threshold)

  proteomics_volcano_plot_with_proteins_label <- proteomics_volcano_plot

  if( include_protein_label == TRUE) {

    proteomics_volcano_plot_with_proteins_label <- proteomics_volcano_plot +
      ggrepel::geom_text_repel( aes( label = gene_name_significant)
                                , show.legend = FALSE
                                , max.overlaps = max.overlaps)

  }

  proteomics_volcano_plot_with_proteins_label

}



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @title Create an Interactive Glimma Volcano Plot for Proteomics
#' @description Generates an interactive volcano plot using the `Glimma` package and
#' saves it as a self-contained HTML file. This plot allows for interactive
#' exploration of differential expression results.
#'
#' @param r_obj A `limma` `MArrayLM` object (the output of `eBayes`).
#' @param coef An integer specifying which coefficient (contrast) from the fit object to plot.
#' @param volcano_plot_tab A data frame for annotation, mapping UniProt accessions to gene names.
#' @param uniprot_column The unquoted column name for UniProt accessions in `volcano_plot_tab`.
#' @param gene_name_column The unquoted column name for gene names in `volcano_plot_tab`.
#' @param display_columns A character vector of additional columns from `volcano_plot_tab` to display in the plot.
#' @param additional_annotations An optional data frame with more annotations to join.
#' @param additional_annotations_join_column The unquoted column name to join `additional_annotations` by.
#' @param counts_tbl An optional matrix of expression counts to display alongside the plot.
#' @param groups A vector specifying the group for each sample in `counts_tbl`.
#' @param output_dir The directory where the output HTML file will be saved.
#'
#' @return This function is called for its side effect of writing an HTML file and does not return a value.
#'
#' @importFrom Glimma glimmaVolcano
#' @importFrom qvalue qvalue
#' @importFrom htmlwidgets saveWidget
#' @importFrom dplyr select distinct left_join mutate pull
#' @importFrom purrr map_chr
#' @importFrom stringr str_split
#' @export
getGlimmaVolcanoProteomics <- function( r_obj
                                        , coef
                                        , volcano_plot_tab
                                        , uniprot_column = best_uniprot_acc
                                        , gene_name_column = gene_name
                                        , display_columns = c(  "PROTEIN_NAMES"   )
                                        , additional_annotations = NULL
                                        , additional_annotations_join_column = NULL
                                        , counts_tbl = NULL
                                        , groups = NULL
                                        , output_dir) {

  if( coef <= ncol(r_obj$coefficients )) {

    best_uniprot_acc <- str_split(rownames(r_obj@.Data[[1]]), " |:" ) |>
      purrr::map_chr(1)

    volcano_plot_tab_cln <- volcano_plot_tab |>
      # dplyr::rename( best_uniprot_acc =  {{uniprot_column}}
      #                , gene_name = {{gene_name_column}}   ) |>
      dplyr::select( {{uniprot_column}}
                     , {{gene_name_column}}, any_of( display_columns) ) |>
      distinct()


    if( !is.null( additional_annotations )
        & !is.null( additional_annotations_join_column ) ) {

      volcano_plot_tab_cln <- volcano_plot_tab_cln |>
        left_join( additional_annotations
                   , by = join_by( {{uniprot_column}} == {{additional_annotations_join_column}} ) ) |>
        dplyr::select( {{uniprot_column}}
                       , {{gene_name_column}}
                       ,any_of( display_columns))

    }

    anno_tbl <- data.frame( uniprot_acc = rownames(r_obj@.Data[[1]])
                            , temp_column = best_uniprot_acc ) |>
      dplyr::rename( {{uniprot_column}} := temp_column)

      anno_tbl <- anno_tbl |>
        left_join( volcano_plot_tab_cln
                   , by = join_by({{uniprot_column}} == {{uniprot_column}}) )  |>
        mutate( gene_name = case_when( is.na( gene_name) ~ {{uniprot_column}},
                                       TRUE ~ gene_name) )

      gene_names <- anno_tbl |>
        dplyr::pull(gene_name)

      rownames( r_obj@.Data[[1]] ) <- gene_names

      r_obj$p.value[,coef] <- qvalue( r_obj$p.value[,coef])$qvalues

      htmlwidgets::saveWidget( widget = Glimma::glimmaVolcano(r_obj
                                                      , coef=coef
                                                      , anno=anno_tbl
                                                      , counts = counts_tbl
                                                      , groups = groups
                                                      , display.columns = colnames(anno_tbl )
                                                      , status=decideTests(r_obj, adjust.method="none")
                                                      , p.adj.method = "none"
                                                      , transform.counts='none'
                                                      ) #the plotly object
                               , file = file.path( output_dir
                                                   , paste0(colnames(r_obj$coefficients)[coef], ".html"))  #the path & file name
                               , selfcontained = TRUE #creates a single html file
      )

    }

}


#' @title Create an Interactive Glimma Volcano Plot Widget for Proteomics
#' @description Generates an interactive volcano plot widget using the `Glimma` package.
#' Unlike `getGlimmaVolcanoProteomics`, this function returns the widget object directly,
#' making it suitable for embedding in R Markdown documents or Shiny apps.
#'
#' @param r_obj A `limma` `MArrayLM` object.
#' @param coef An integer specifying which coefficient (contrast) to plot.
#' @param volcano_plot_tab A data frame for annotation.
#' @param uniprot_column The unquoted column name for UniProt accessions.
#' @param gene_name_column The unquoted column name for gene names.
#' @param display_columns Additional columns to display.
#' @param additional_annotations Optional additional annotation data frame.
#' @param additional_annotations_join_column Unquoted column name for joining additional annotations.
#' @param counts_tbl Optional matrix of expression counts.
#' @param groups A vector of group assignments for samples.
#'
#' @return An `htmlwidget` object.
#'
#' @importFrom Glimma glimmaVolcano
#' @export
getGlimmaVolcanoProteomicsWidget <- function( r_obj
                                        , coef
                                        , volcano_plot_tab
                                        , uniprot_column = best_uniprot_acc
                                        , gene_name_column = gene_name
                                        , display_columns = c(  "PROTEIN_NAMES"   )
                                        , additional_annotations = NULL
                                        , additional_annotations_join_column = NULL
                                        , counts_tbl = NULL
                                        , groups = NULL) {

  if( coef <= ncol(r_obj$coefficients )) {

    best_uniprot_acc <- str_split(rownames(r_obj@.Data[[1]]), " |:" ) |>
      purrr::map_chr(1)

    # print(paste("nrow = ", nrow(r_obj@.Data[[1]])))
    # print(head(best_uniprot_acc))

    volcano_plot_tab_cln <- volcano_plot_tab  |>
      dplyr::select ( {{uniprot_column}}
                       , {{gene_name_column}}, any_of( display_columns)) |>
      distinct()

    # print (head( volcano_plot_tab_cln))

    if( !is.null( additional_annotations )
        & !is.null( additional_annotations_join_column ) ) {

      volcano_plot_tab_cln <- volcano_plot_tab_cln |>
        left_join( additional_annotations
                   , by = join_by( {{uniprot_column}} == {{additional_annotations_join_column}} ) ) |>
        dplyr::select( {{uniprot_column}}
                       , {{gene_name_column}}
                       , any_of( display_columns))
    }

    anno_tbl <- data.frame( uniprot_acc = rownames(r_obj@.Data[[1]]) # This uniprot_acc does not matter, only shows in glimma Volcano table
                            , temp_column = best_uniprot_acc ) |>
      dplyr::rename( {{uniprot_column}} := temp_column) |>
      left_join( volcano_plot_tab_cln
                 , by = join_by({{uniprot_column}} == {{uniprot_column}}) )  |>
      mutate( gene_name = case_when( is.na( gene_name) ~ {{uniprot_column}},
                                     TRUE ~ gene_name) )

    gene_names <- anno_tbl |>
      dplyr::pull(gene_name)

    rownames( r_obj@.Data[[1]] ) <- gene_names

    r_obj$p.value[,coef] <- qvalue( r_obj$p.value[,coef])$qvalues

     glimmaVolcano(r_obj
                   , coef=coef
                   , counts = counts_tbl
                   , groups = groups
                   , anno=anno_tbl
                   , display.columns = display_columns
                   , status=decideTests(r_obj, adjust.method="none")
                   , p.adj.method="none"
                   , transform.counts='none') #the plotly object

  }

}

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @title Create an Interactive Glimma Volcano Plot for Phosphoproteomics
#' @description Generates an interactive volcano plot using `Glimma` specifically
#' for phosphoproteomics data, using phosphosite IDs for annotation. Saves the
#' plot as a self-contained HTML file.
#'
#' @param r_obj A `limma` `MArrayLM` object.
#' @param coef An integer specifying which coefficient (contrast) to plot.
#' @param volcano_plot_tab A data frame for annotation, mapping site IDs to other information.
#' @param sites_id_column The unquoted column name for phosphosite identifiers.
#' @param sites_id_display_column The unquoted column name for the labels to display in the plot.
#' @param display_columns Additional columns from `volcano_plot_tab` to show.
#' @param additional_annotations Optional additional annotation data frame.
#' @param additional_annotations_join_column Unquoted column name for joining additional annotations.
#' @param counts_tbl Optional matrix of expression counts.
#' @param output_dir The directory to save the output HTML file.
#'
#' @return This function is called for its side effect of writing an HTML file.
#'
#' @importFrom Glimma glimmaVolcano
#' @export
getGlimmaVolcanoPhosphoproteomics <- function( r_obj
                                        , coef
                                        , volcano_plot_tab
                                        , sites_id_column = sites_id
                                        , sites_id_display_column = sites_id_short
                                        , display_columns = c(  "sequence", "PROTEIN_NAMES"   )
                                        , additional_annotations = NULL
                                        , additional_annotations_join_column = NULL
                                        , counts_tbl = NULL
                                        , output_dir) {

  if( coef <= ncol(r_obj$coefficients )) {

    volcano_plot_tab_cln <- volcano_plot_tab |>
      dplyr::distinct( {{sites_id_column}}
                       , {{sites_id_display_column}}
                       , any_of( display_columns) )

    if( !is.null( additional_annotations )
        & !is.null( additional_annotations_join_column ) ) {

      volcano_plot_tab_cln <- volcano_plot_tab_cln |>
        left_join( additional_annotations
                   , by = join_by( {{sites_id_column}} == {{additional_annotations_join_column}} ) ) |>
        dplyr::select( any_of( display_columns))
    }

    anno_tbl <-  data.frame(  sites_id = rownames(r_obj@.Data[[1]])) |>
      left_join( volcano_plot_tab_cln
                 , by = join_by(sites_id == {{sites_id_column}} ) )

    sites_id_short_list <- anno_tbl |>
                   dplyr::pull(sites_id_short)

    rownames( r_obj@.Data[[1]] ) <- sites_id_short_list

    #coef <- seq_len( ncol(r_obj$coefficients))[1]

    r_obj$p.value[,coef] <- qvalue( r_obj$p.value[,coef])$qvalues


    htmlwidgets::saveWidget( widget = glimmaVolcano(r_obj, coef=coef
                                                    , counts = counts_tbl
                                                    , anno=anno_tbl
                                                    , display.columns=display_columns
                                                    , p.adj.method = "none"
                                                    , transform.counts='none' ) #the plotly object
                             , file = file.path( output_dir
                                                 , paste0(colnames(r_obj$coefficients)[coef], ".html"))  #the path & file name
                             , selfcontained = TRUE #creates a single html file
    )
  }

}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' @title Plot P-Value Distribution
#' @description Generates a histogram of raw p-values to visualize their distribution.
#' This is a common diagnostic plot to assess the validity of the statistical model.
#' A uniform distribution is often expected under the null hypothesis.
#'
#' @param selected_data A data frame containing a column of p-values.
#' @param p_value_column The unquoted column name for the raw p-values.
#' @param formula_string A formula string (e.g., `"is_ruv_applied ~ comparison"`) for faceting with `facet_grid`.
#'
#' @return A ggplot object showing the p-value distribution histogram.
#'
#' @importFrom ggplot2 ggplot aes theme xlab geom_histogram facet_grid
#' @export
printPValuesDistribution <- function(selected_data, p_value_column = raw_pvalue, formula_string = "is_ruv_applied ~ comparison") {

  breaks <- c(0, 0.001, 0.01, 0.05,
              seq(0.1, 1, by = 0.1))

# after_stat(density)
  pvalhist <- ggplot(selected_data, aes({ { p_value_column } })) +
    theme(axis.title.y = element_blank()) +
    xlab("P-value") +
    geom_histogram(aes(y = after_stat(density)),
                   breaks = breaks,
                   position = "identity",
                   color = "black") +
    geom_histogram(aes(y = after_stat(density)),
                   breaks = breaks,
                   position = "identity")

  if (!is.na(formula_string)) {

    pvalhist <- pvalhist +
      facet_grid(as.formula(formula_string))
  }


  pvalhist

}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' @title Fit a Linear Model and Compute Empirical Bayes Statistics
#' @description A wrapper for the core `limma` workflow. It fits a linear model to
#' expression data, computes contrasts, and applies Empirical Bayes moderation to
#' compute moderated t-statistics, p-values, and q-values.
#'
#' @param data A numeric matrix or data frame of expression data (features x samples).
#' @param design A design matrix created by `model.matrix`.
#' @param contr.matrix A contrast matrix created by `makeContrasts`.
#'
#' @return A list containing two elements:
#'   - `table`: A data frame with detailed statistics for the first contrast, including
#'     logFC, t-statistics, p-values, and q-values.
#'   - `fit.eb`: The `MArrayLM` fit object returned by `eBayes`.
#'
#' @importFrom limma lmFit contrasts.fit eBayes
#' @importFrom qvalue qvalue
#' @export
ebFit <- function(data, design, contr.matrix)
{
  fit <- lmFit(data, design)
  fit.c <- contrasts.fit(fit, contrasts = contr.matrix)

  fit.eb <- suppressWarnings(eBayes(fit.c))

  logFC <- fit.eb$coefficients[, 1]
  df.r <- fit.eb$df.residual
  df.0 <- rep(fit.eb$df.prior, dim(data)[1])
  s2.0 <- rep(fit.eb$s2.prior, dim(data)[1])
  s2 <- (fit.eb$sigma)^2
  s2.post <- fit.eb$s2.post
  t.ord <- fit.eb$coefficients[, 1] /
    fit.eb$sigma /
    fit.eb$stdev.unscaled[, 1]
  t.mod <- fit.eb$t[, 1]
  p.ord <- 2 * pt(-abs(t.ord), fit.eb$df.residual)
  raw_pvalue <- fit.eb$p.value[, 1]
  q.ord <- qvalue(p.ord)$q
  fdr_qvalue <- qvalue(raw_pvalue)$q

  return(list(table = data.frame(logFC, t.ord, t.mod, p.ord, raw_pvalue, q.ord, fdr_qvalue, df.r, df.0, s2.0, s2, s2.post),
              fit.eb = fit.eb))
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @title Run a Differential Expression Test for a Single Contrast
#' @description This function performs a differential expression test for a single
#' pairwise comparison between two experimental groups using the `limma` package.
#'
#' @param ID A vector of feature identifiers (e.g., protein accessions).
#' @param A A numeric matrix of expression data for group A (features x samples).
#' @param B A numeric matrix of expression data for group B (features x samples).
#' @param group_A A string naming group A.
#' @param group_B A string naming group B.
#' @param design_matrix A design matrix for the samples in groups A and B.
#' @param formula_string A formula string for `model.matrix` (e.g., `"~ 0 + group"`).
#' @param contrast_variable The name of the contrast variable in the design matrix.
#' @param weights An optional numeric matrix of weights for the linear model fit.
#'
#' @return A list containing two elements:
#'   - `table`: A data frame summarizing the results (logFC, p-values, q-values, etc.).
#'   - `fit.eb`: The `MArrayLM` fit object from `eBayes`.
#'
#' @importFrom limma makeContrasts
#' @export
runTest <- function(ID, A, B, group_A, group_B, design_matrix, formula_string,
                    contrast_variable = "group",
                    weights = NA) {


  ff <- as.formula(formula_string)
  mod_frame <- model.frame(ff, design_matrix)
  design_m <- model.matrix(ff, mod_frame)


  # print("My design matrix")
  # print(design_m)
  # print( paste( "nrow(weights)", nrow(weights), "nrow(design_m)", nrow(design_m)))

  if (!is.na(weights)) {
    if (nrow(weights) == nrow(design_m)) {
      design_m <- cbind(design_m, weights)
    } else {
      stop("Stop: nrow(weights) should be equal to nrow(design_m)")
    }

  }

  # print(paste("group_A = ", group_A))
  # print(paste("group_B = ", group_B))

  contr.matrix <- makeContrasts(contrasts = paste0(group_B, "vs", group_A, "=", contrast_variable, group_B, "-", contrast_variable, group_A),
                                levels = colnames(design_m))

  eb_fit_list <- ebFit(cbind(A, B), design_m, contr.matrix = contr.matrix)

  r <- eb_fit_list$table
  fit.eb <- eb_fit_list$fit.eb

  return(list(table = data.frame(row.names = row.names(r),
                                 comparison = paste("log(", group_B, ") minus log(", group_A, ")", sep = ""),
                                 meanA = rowMeans(A),
                                 meanB = rowMeans(B),
                                 logFC = r$logFC,
                                 tstats = r$t.ord,
                                 tmod = r$t.mod,
                                 pval = r$p.ord,
                                 raw_pvalue = r$raw_pvalue,
                                 qval = r$q.ord,
                                 fdr_qvalue = r$fdr_qvalue),
              fit.eb = fit.eb))
}

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' @title Get Grouping Information from Design Matrix
#' @description This function processes a design matrix to create a named list
#' where each name is an experimental group and each element is a character vector
#' of the sample IDs belonging to that group.
#'
#' @param design_matrix A data frame containing the experimental design.
#' @param group_id The unquoted column name for the grouping variable.
#' @param sample_id The unquoted column name for the sample identifiers.
#'
#' @return A named list mapping group names to sample ID vectors.
#'
#' @importFrom dplyr select group_by summarise ungroup pull
#' @importFrom rlang sym
#' @export
getTypeOfGrouping <- function(design_matrix, group_id, sample_id) {
  temp_type_of_grouping <- design_matrix |>
    dplyr::select(!!rlang::sym(group_id), !!rlang::sym(sample_id)) |>
    group_by(!!rlang::sym(group_id)) |>
    summarise(!!rlang::sym(sample_id) := list(!!rlang::sym(sample_id))) |>
    ungroup()

  type_of_grouping <- temp_type_of_grouping |> dplyr::pull(!!rlang::sym(sample_id))
  names(type_of_grouping) <- temp_type_of_grouping |> dplyr::pull(!!rlang::sym(group_id))

  return(type_of_grouping)

}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' @title Run Pairwise Differential Expression Tests for Multiple Comparisons
#' @description This function iterates through a set of pairwise comparisons defined
#' in `test_pairs` and runs a differential expression test for each using the `runTest` helper function.
#'
#' @details
#' For each pair of groups to be compared, this function subsets the data and design
#' matrix. It can also handle pre-filtered feature lists (`sample_rows_list`) to
#' test only a subset of features relevant to the specific comparison.
#'
#' @param ID A vector of all feature identifiers in the `data` matrix.
#' @param data The full expression data matrix (features x samples).
#' @param test_pairs A data frame with two columns ("A" and "B") defining the pairwise comparisons.
#' @param sample_columns A character vector of all sample column names to be included in the analysis.
#' @param sample_rows_list An optional named list (output of `getRowsToKeepList`) where each
#'   element contains the feature IDs to be tested for a specific group.
#' @param type_of_grouping A named list (output of `getTypeOfGrouping`) mapping group names to sample IDs.
#' @param design_matrix The full design matrix for the experiment.
#' @param formula_string A formula string for the linear model.
#' @param contrast_variable The name of the contrast variable in the design matrix.
#' @param weights An optional numeric matrix of weights for the linear model.
#'
#' @return A named list of results. Each element corresponds to a comparison and contains:
#'   - `results`: A data frame of DE statistics from `runTest`.
#'   - `counts`: The subset of the data matrix used for the test.
#'   - `fit.eb`: The `MArrayLM` fit object from `eBayes`.
#' @export
runTests <- function(ID, data, test_pairs, sample_columns, sample_rows_list = NA, type_of_grouping, design_matrix, formula_string, contrast_variable = "group", weights = NA) {
  r <- list()
  for (i in 1:nrow(test_pairs)) {

    rows_to_keep <- rownames(data)


    if (length(sample_rows_list) > 0) {
      if (!is.na(sample_rows_list) &
        #  Check that sample group exists as names inside sample_rows_list
        length(which(c(test_pairs[i, "A"], test_pairs[i, "B"]) %in% names(sample_rows_list))) > 0) {

        rows_to_keep <- unique(sample_rows_list[[ test_pairs[[i, "A"]]]],
                               sample_rows_list[[ test_pairs[[i, "B"]]]])
      }
    }

    tmp <- data[rows_to_keep, sample_columns]
    rep <- colnames(tmp)

    # print( paste( test_pairs[i,]$A, test_pairs[i,]$B) )
    A <- tmp[, type_of_grouping[test_pairs[i,]$A][[1]]]
    B <- tmp[, type_of_grouping[test_pairs[i,]$B][[1]]]

    subset_weights <- NA

    if (!is.na(weights)) {
      subset_weights <- weights[c(colnames(A), colnames(B)),]
    }

    # print(colnames(A))
    # print(colnames(B))
    tmp <- unname(cbind(A, B))
    Aname <- paste(test_pairs[i,]$A, 1:max(1, ncol(A)), sep = "_")
    Bname <- paste(test_pairs[i,]$B, 1:max(1, ncol(B)), sep = "_")
    colnames(tmp) <- c(Aname, Bname)

    selected_sample_ids <- c(type_of_grouping[test_pairs[i,]$A][[1]], type_of_grouping[test_pairs[i,]$B][[1]])
    design_matrix_subset <- design_matrix[selected_sample_ids, , drop = FALSE]

    # print("My design matrix 1")
    # print( selected_sample_ids)
    # print(design_matrix)
    # print( dim(design_matrix))

    group_A <- test_pairs[i,]$A
    group_B <- test_pairs[i,]$B

    x <- runTest(ID, A, B, group_A, group_B, design_matrix = design_matrix_subset,
                 formula_string = formula_string, contrast_variable = contrast_variable,
                 weights = subset_weights)

    comparison <- paste(group_B, " vs ", group_A, sep = "")

    r[[comparison]] <- list(results = x$table, counts = t(cbind(A, B)), fit.eb = x$fit.eb)
  }
  r
}

## -----------

#' @title Run Differential Expression Analysis Using `limma` Contrasts
#' @description A flexible function to perform differential expression analysis for
#' a set of specified contrasts using the `limma` package. It supports `eBayes`
#' moderation and the `treat` method for testing against a log-fold-change threshold.
#'
#' @param data The expression data matrix (features x samples).
#' @param contrast_strings A character vector of contrast formulas (e.g., "GroupA-GroupB").
#' @param design_matrix The full design matrix for the experiment.
#' @param formula_string A formula string for the linear model (e.g., `"~ 0 + group"`).
#' @param p_value_column The unquoted column name to store the raw p-values in the output tables.
#' @param q_value_column The unquoted column name to store q-values (from the `qvalue` package).
#' @param fdr_value_column The unquoted column name to store BH-adjusted p-values.
#' @param weights An optional numeric matrix of weights.
#' @param treat_lfc_cutoff An optional log-fold-change threshold for `limma::treat`.
#'   If provided, tests if the fold change is significantly greater than this value.
#' @param eBayes_trend A logical value passed to `limma::eBayes`, indicating whether
#'   to allow for a trend in the prior variance.
#' @param eBayes_robust A logical value passed to `limma::eBayes`, indicating whether
#'   to use robust estimation for the prior variance.
#'
#' @return A list containing two elements:
#'   - `results`: A named list of data frames, where each data frame contains the DE
#'     results for one of the specified contrasts.
#'   - `fit.eb`: The final `MArrayLM` fit object from `eBayes` or `treat`.
#'
#' @importFrom limma makeContrasts lmFit contrasts.fit eBayes treat topTable topTreat
#' @importFrom qvalue qvalue
#' @importFrom purrr map
#' @export
runTestsContrasts <- function(data,
                              contrast_strings,
                              design_matrix,
                              formula_string,
                              p_value_column = raw_pvalue,
                              q_value_column = fdr_qvalue,
                              fdr_value_column = fdr_value_bh_adjustment,
                              weights = NA,
                              treat_lfc_cutoff = NA,
                              eBayes_trend = FALSE,
                              eBayes_robust = FALSE) {

  ff <- as.formula(formula_string)
  mod_frame <- model.frame(ff, design_matrix)
  design_m <- model.matrix(ff, mod_frame)

  data_subset <- data[, rownames( design_m)]

  ## Make contrasts
  contr.matrix <- makeContrasts(contrasts = contrast_strings,
                                levels = colnames(design_m))

  ## Attach weights
  if (!is.na(weights)) {
    if (nrow(weights) == nrow(design_m)) {
      design_m <- cbind(design_m, weights)
    } else {
      stop("Stop: nrow(weights) should be equal to nrow(design_m)")
    }

  }

  fit <- lmFit(data_subset, design = design_m)

  cfit <- contrasts.fit(fit, contrasts = contr.matrix)

  eb.fit <- eBayes( cfit, trend = eBayes_trend, robust = eBayes_robust )

  ## Run treat over here
  t.fit <- NA
  result_tables <- NA
  if( !is.na( treat_lfc_cutoff)) {
    t.fit <- treat(eb.fit, lfc=as.double(treat_lfc_cutoff)) ## assign log fold change threshold below which is scientifically not relevant

    result_tables <- purrr::map(contrast_strings,
                                function(contrast) {
                                  de_tbl <- topTreat(t.fit, coef = contrast, n = Inf) |>
                                    mutate({ { q_value_column } } := qvalue(P.Value)$q) |>
                                    mutate({ { fdr_value_column } } := p.adjust(P.Value, method="BH")) |>

                                    dplyr::rename({ { p_value_column } } := P.Value)
                                }
    )
  } else {
    t.fit <- eb.fit

    result_tables <- purrr::map(contrast_strings,
                                function(contrast) {
                                  de_tbl <- topTable(t.fit, coef = contrast, n = Inf) |>
                                    mutate({ { q_value_column } } := qvalue(P.Value)$q) |>
                                    mutate({ { fdr_value_column } } := p.adjust(P.Value, method="BH")) |>

                                    dplyr::rename({ { p_value_column } } := P.Value)
                                }
    )
  }



  names(result_tables) <- contrast_strings


  return(list(results = result_tables,
              fit.eb = t.fit))
}



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @title Extract RUV-Corrected Results
#' @description A simple helper function to extract the `results` data frame from
#' each element of a list returned by `runTests` or a similar function.
#'
#' @param results_list A named list where each element is a list containing a `results` data frame.
#'
#' @return A named list containing only the `results` data frames.
#'
#' @importFrom purrr map
#' @export
extractRuvResults <- function(results_list) {

  extracted <- purrr::map(results_list, \(x){ x$results })

  names(extracted) <- names(results_list)

  return(extracted)
}


#' @title Extract DE Results Tables
#' @description A simple helper function to extract the `results` data frame from
#' each element of a list returned by `runTests` or a similar function. It is
#' functionally identical to `extractRuvResults` but may be used for semantic clarity.
#'
#' @param results_list A named list where each element is a list containing a `results` data frame.
#'
#' @return A named list containing only the `results` data frames.
#'
#' @importFrom purrr map
#' @export
extractResults <- function(results_list) {

  extracted <- purrr::map(results_list, \(x){ x$results })

  names(extracted) <- names(results_list)

  return(extracted)
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' @title Save a List of DE Result Tables to Files
#' @description Iterates over a named list of differential expression result tables
#' and saves each table to a separate file.
#'
#' @param list_of_de_tables A named list of data frames. The names of the list
#'   elements are used to construct the output filenames.
#' @param row_id The unquoted column name to use for the feature identifiers when converting
#'   rownames to a column before saving.
#' @param sort_by_column The unquoted column name to sort each table by before saving.
#' @param results_dir The directory where the files will be saved.
#' @param file_suffix The suffix (including extension, e.g., "_results.tsv") to
#'   append to each filename.
#'
#' @return This function is called for its side effect of writing files and does not return a value.
#'
#' @importFrom vroom vroom_write
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr arrange
#' @importFrom purrr walk2
#' @export
saveDeProteinList <- function(list_of_de_tables, row_id, sort_by_column = fdr_qvalue, results_dir, file_suffix) {

  purrr::walk2(list_of_de_tables, names(list_of_de_tables),

               \(.x) {vroom::vroom_write(.x |>
                                     rownames_to_column(row_id) |>
                                     arrange({ { sort_by_column } }),
                                   path = file.path(results_dir, paste0(.y, file_suffix)))} )

}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' @title Identify Negative Control Proteins Using ANOVA
#' @description This function identifies a set of negative control features (e.g., proteins)
#' for use in normalization methods like RUV. It performs an ANOVA test for each
#' feature across experimental groups and selects features that are least likely
#' to be differentially expressed (i.e., have the highest p-values/q-values).
#'
#' @details
#' The function calculates p-values using ANOVA (`stats::aov`). It then adjusts
#' these p-values using the specified `ruv_fdr_method` (`qvalue` or `BH`).
#' Features with adjusted p-values above `ruv_qval_cutoff` are considered candidates.
#' The top `num_neg_ctrl` features from this candidate list (those with the highest
#' adjusted p-values) are selected as negative controls.
#'
#' @param data_matrix A numeric matrix of expression data (features x samples).
#' @param design_matrix A data frame with sample metadata.
#' @param grouping_variable The column name in `design_matrix` for the experimental groups.
#' @param percentage_as_neg_ctrl The percentage of total features to select as negative controls.
#'   This is used to calculate `num_neg_ctrl` if `num_neg_ctrl` is not provided. Defaults to 10.
#' @param num_neg_ctrl The number of negative control features to select. Overrides `percentage_as_neg_ctrl`.
#' @param ruv_qval_cutoff The q-value threshold. Only features with a q-value *above*
#'   this threshold are considered for the control set.
#' @param ruv_fdr_method The method for FDR calculation ("qvalue" or "BH").
#'
#' @return A logical vector with the same length as the number of rows in `data_matrix`.
#'   `TRUE` indicates that the feature is selected as a negative control. The vector is
#'   named with the feature IDs from the `data_matrix` rownames.
#'
#' @importFrom stats aov
#' @importFrom qvalue qvalue
#' @export
getNegCtrlProtAnovaHelper <- function(data_matrix
                                , design_matrix
                                , grouping_variable = "group"
                                , percentage_as_neg_ctrl = 10
                                , num_neg_ctrl = round( nrow(data_matrix)*percentage_as_neg_ctrl/100, 0)
                                , ruv_qval_cutoff = 0.05
                                , ruv_fdr_method = "qvalue") {

  ## Both percentage_as_neg_ctrl and num_neg_ctrl is missing, and number of proteins >= 50 use only 10 percent of the proteins as negative control by default
  if((is.null(percentage_as_neg_ctrl) ||
     is.na(percentage_as_neg_ctrl) ) &&
     (is.null(num_neg_ctrl) ||
      is.na(num_neg_ctrl) ) &&
     nrow(data_matrix) >= 50 ) {
    num_neg_ctrl <- round( nrow(data_matrix)*10/100, 0)
    warnings( paste0( getFunctionName(), ": Using 10% of proteins from the input matrix as negative controls by default.\n"))
  } else if (!is.null(percentage_as_neg_ctrl) &
             !is.na(percentage_as_neg_ctrl)) {
    num_neg_ctrl <- round( nrow(data_matrix)*percentage_as_neg_ctrl/100, 0)
  } else if(!is.null(num_neg_ctrl) &
       !is.na(num_neg_ctrl)) {
    num_neg_ctrl <- as.integer(num_neg_ctrl)
  } else {
    stop(paste0( getFunctionName(), ": Please provide either percentage_as_neg_ctrl or num_neg_ctrl.\n"))
  }

  ## Inspired by matANOVA function from PhosR package: http://www.bioconductor.org/packages/release/bioc/html/PhosR.html

  grps <- design_matrix[colnames(data_matrix), grouping_variable]

  ps <- apply(data_matrix, 1, function(x) {
       if( length( unique( grps[!is.na(x)] )  ) > 1 ) {
         summary(stats::aov(as.numeric(x) ~ grps))[[1]][["Pr(>F)"]][1]
       } else {
          return(NA_real_)
       }
    })

  ps[is.na(ps)] <- 1

  aov <- c()

  if ( ruv_fdr_method == "qvalue") {
    aov <- qvalue(unlist(ps))$qvalues
  } else if ( ruv_fdr_method == "BH") {
    aov <- qvalue(unlist(ps), pi0=1)$qvalues
  } else {
    error( paste( "Input FDR method", ruv_fdr_method, "not valid") )
  }

  filtered_list <- aov[aov > ruv_qval_cutoff]

  list_size <- ifelse(num_neg_ctrl > length(filtered_list), length(filtered_list), num_neg_ctrl)

  control_genes <- names(sort(filtered_list, decreasing = TRUE)[1:list_size])

  #nrow(data_matrix) - length(control_genes)
  control_genes_index <- rownames(data_matrix) %in% control_genes
  names(control_genes_index) <- rownames(data_matrix)

  return(control_genes_index)

}



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @title Find the Best K for RUV Normalization
#' @description Analyzes a canonical correlation plot from `ruv::RUV_cancor` to
#' determine the optimal number of factors of unwanted variation (k).
#'
#' @details The function works by finding the value of k that maximizes the
#' difference between the canonical correlation of all features and the canonical
#' correlation of control features. This point is assumed to represent the best
#' trade-off in removing unwanted variation without removing the signal of interest.
#'
#' @param cancorplot_r1 A ggplot object produced by `ruv::RUV_cancor`. The plot's
#'   data layer must contain `featureset`, `cc` (canonical correlation), and `K` columns.
#'
#' @return An integer representing the determined optimal value of k.
#' @export
findBestK <- function( cancorplot_r1) {
  controls_idx <- which(cancorplot_r1$data$featureset == "Control")
  all_idx <- which( cancorplot_r1$data$featureset == "All")
  difference_between_all_ctrl <- cancorplot_r1$data$cc[all_idx] - cancorplot_r1$data$cc[controls_idx]
  max_difference <- max(difference_between_all_ctrl, na.rm=TRUE)
  best_idx <- which( difference_between_all_ctrl == max_difference)
  best_k <- (cancorplot_r1$data$K[controls_idx] )[best_idx]
  return( best_k)
}

#' @title Find the Best K Value for RUV Across Multiple Assays
#' @description This function iterates over a list of canonical correlation plots
#' (e.g., from a multi-assay experiment) and applies `findBestK` to each to
#' determine the optimal number of unwanted variation factors (k) for each assay.
#'
#' @param cancor_plots_list A named list where each element is a ggplot object
#'   returned by `ruv_cancorplot` (the output of `ruvCancor` for a multi-assay object).
#'
#' @return A named list where keys are the assay names (from the input list)
#'   and values are the determined best K for each assay. Returns `NA_integer_`
#'   for assays where `findBestK` fails or returns an invalid result.
#'
#' @importFrom purrr map
#' @importFrom logger log_warn
#' @export
findBestKForAssayList <- function(cancor_plots_list) {

    # --- Input Validation ---
    if (!is.list(cancor_plots_list)) {
        stop("`cancor_plots_list` must be a list.")
    }
    if (length(cancor_plots_list) == 0) {
        log_warn("`cancor_plots_list` is empty. Returning an empty list.")
        return(list())
    }
    if (is.null(names(cancor_plots_list)) || any(names(cancor_plots_list) == "")) {
        log_warn("Input list is not fully named. Assigning default names (Plot_1, Plot_2, ...)")
        names(cancor_plots_list) <- paste0("Plot_", seq_along(cancor_plots_list))
    }

    # --- Iterate and Apply findBestK ---
    best_k_list <- purrr::map(seq_along(cancor_plots_list), function(i) {
        assay_name <- names(cancor_plots_list)[i]
        current_plot <- cancor_plots_list[[i]]

        # Check if it's a ggplot object
        if (!inherits(current_plot, "ggplot")) {
            log_warn("Element '{assay_name}' in the list is not a ggplot object. Skipping.", .logr = TRUE)
            return(NA_integer_)
        }

        # Call the original findBestK function
        best_k_assay <- tryCatch({
            k <- findBestK(current_plot)
            # Add a check for Inf or non-numeric result which might cause issues later
            if (is.infinite(k) || !is.numeric(k) || length(k) != 1) {
                 log_warn("Assay '{assay_name}': findBestK returned an invalid value ({k}). Returning NA.", .logr = TRUE)
                 NA_integer_
            } else {
                 # Ensure it's an integer
                 as.integer(k)
            }
        }, warning = function(w) {
            # Catch the specific warning from findBestK if max returns -Inf
            if (grepl("no non-missing arguments to max", w$message, ignore.case = TRUE)) {
                log_warn("Assay '{assay_name}': Warning in findBestK (likely no valid data in plot): {w$message}. Returning NA.", .logr = TRUE)
            } else {
                # Log other warnings but still try to return NA
                log_warn("Assay '{assay_name}': Warning during findBestK: {w$message}. Returning NA.", .logr = TRUE)
            }
            NA_integer_
        }, error = function(e) {
            log_warn("Assay '{assay_name}': Error calling findBestK: {e$message}. Returning NA.", .logr = TRUE)
            NA_integer_
        })

        return(best_k_assay)
    })

    # Set names for the result list
    names(best_k_list) <- names(cancor_plots_list)

    return(best_k_list)
}

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' @title Average Values from Technical Replicates
#' @description This function takes a wide-format data table and averages the values
#' of technical replicates, returning a new matrix with one column per biological sample.
#'
#' @param input_table A data frame or matrix with features in rows and samples in columns.
#' @param design_matrix A data frame mapping sample IDs to biological replicate IDs.
#' @param group_pattern A regular expression to identify the sample columns in `input_table`.
#' @param row_id The unquoted column name for feature identifiers.
#' @param sample_id The unquoted column name for the sample identifiers (technical replicates).
#' @param average_replicates_id The unquoted column name in `design_matrix` that identifies
#'   the biological samples to average over.
#'
#' @return A matrix with features in rows and averaged biological samples in columns.
#'
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom dplyr left_join group_by summarise ungroup
#' @importFrom rlang sym
#' @export
averageValuesFromReplicates <- function(input_table, design_matrix, group_pattern, row_id, sample_id, average_replicates_id) {

  output_table <- input_table |>
    as.data.frame() |>
    rownames_to_column(  row_id   ) |>
    pivot_longer( cols=matches(group_pattern),
                  names_to = sample_id,
                  values_to = "value") |>
    left_join( design_matrix, by = sample_id ) |>
    group_by( !!rlang::sym(average_replicates_id) ,  !!rlang::sym(row_id) ) |>
    summarise( value = mean(value, na.rm = TRUE)) |>
    ungroup() |>
    pivot_wider( names_from = !!rlang::sym(average_replicates_id),
                 values_from = "value") |>
    column_to_rownames(row_id) |>
    as.matrix()

  return( output_table )
}



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# up <- UniProt.ws(taxId=10090)
# keytypes(up)
# columns(up)
# test <- batch_query_evidence(subset_tbl, Proteins)

#' @title Clean Isoform Numbers from Protein Accessions
#' @description A simple string manipulation function that removes isoform suffixes
#' (e.g., "-1", "-2") from UniProt protein accession numbers.
#'
#' @param string A character vector of protein accession numbers.
#'
#' @return A character vector with isoform numbers removed.
#'
#' @importFrom stringr str_replace
#' @export
cleanIsoformNumber <- function(string) {
  # "Q8K4R4-2"
  str_replace(string, "-\\d+$", "")

}


# Filter for a batch and run analysis on that batch of uniprot accession keys only.
subsetQuery <- function(data, subset, accessions_col_name, uniprot_handle, uniprot_columns = c("EXISTENCE", "SCORE", "REVIEWED", "GENENAME", "PROTEIN-NAMES", "LENGTH"),
                        uniprot_keytype = "UniProtKB") {


  # print(subset)
  my_keys <- data |>
    dplyr::filter(round == subset) |>
    dplyr::pull({ { accessions_col_name } })

   print(head( my_keys) )

  # print(uniprot_keytype)
  # Learning in progress

  output <- UniProt.ws::select(uniprot_handle,
                     keys = my_keys,
                     columns = uniprot_columns,
                     keytype = uniprot_keytype)

  print(head(output))
}


# The UniProt.ws::select function limits the number of keys queried to 100. This gives a batch number for it to be queried in batches.
batchQueryEvidenceHelper <- function(uniprot_acc_tbl, uniprot_acc_column) {

  all_uniprot_acc <- uniprot_acc_tbl |>
    dplyr::select({ { uniprot_acc_column } }) |>
    mutate(Proteins = str_split({ { uniprot_acc_column } }, ";")) |>
    unnest(Proteins) |>
    distinct() |>
    arrange(Proteins) |>
    mutate(Proteins = cleanIsoformNumber(Proteins)) |>
    dplyr::mutate(round = ceiling(row_number() / 100))  ## 100 is the maximum number of queries at one time
}

#' @title Batch Query UniProt for Protein Evidence
#' @description This function queries the UniProt database in batches to retrieve
#' annotation information for a list of protein accessions. It is designed to
#' handle a large number of queries by splitting them into smaller chunks, respecting
#' UniProt's API limits.
#'
#' @param uniprot_acc_tbl A data frame containing the protein accessions to query.
#' @param uniprot_acc_column The unquoted column name in `uniprot_acc_tbl` that contains the protein accessions.
#' @param uniprot_handle An active `UniProt.ws` handle.
#' @param uniprot_columns A character vector of the UniProt annotation columns to retrieve.
#' @param uniprot_keytype The keytype to use for the query (e.g., "UniProtKB").
#'
#' @return A single data frame containing the combined annotation results from all batches.
#'
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' @export
batchQueryEvidence <- function(uniprot_acc_tbl, uniprot_acc_column, uniprot_handle,
                               uniprot_columns = c("EXISTENCE", "SCORE", "REVIEWED", "GENENAME", "PROTEIN-NAMES", "LENGTH"),
                               uniprot_keytype = "UniProtKB") {

  all_uniprot_acc <- batchQueryEvidenceHelper(uniprot_acc_tbl,
                                              { { uniprot_acc_column } })

  partial_subset_query <- partial(subsetQuery,
                                  data = all_uniprot_acc,
                                  accessions_col_name = { { uniprot_acc_column } },
                                  uniprot_handle = uniprot_handle,
                                  uniprot_columns = uniprot_columns,
                                  uniprot_keytype = uniprot_keytype)

  rounds_list <- all_uniprot_acc |>
    dplyr::distinct(round) |>
    dplyr::arrange(round) |>
    dplyr::pull(round)

  all_uniprot_evidence <- purrr::map(rounds_list, \(x){ partial_subset_query(subset = x) }) |>
    dplyr::bind_rows()

  return(all_uniprot_evidence)
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# The UniProt.ws::select function limits the number of keys queried to 100. This gives a batch number for it to be queried in batches.
batchQueryEvidenceHelperGeneId <- function(input_tbl, gene_id_column, delim =":") {

  all_uniprot_acc <- input_tbl |>
    dplyr::select({ { gene_id_column } }) |>
    separate_longer_delim({ { gene_id_column } }, delim= delim ) |>
    dplyr::arrange({ { gene_id_column } }) |>
    dplyr::mutate(round = ceiling(dplyr::row_number() / 25))  ## 25 is the maximum number of queries at one time

  return(all_uniprot_acc)
}

#' @title Batch Query UniProt by Gene ID
#' @description This function queries the UniProt database in batches using gene
#' identifiers. It splits a list of gene IDs into manageable chunks for querying.
#'
#' @param input_tbl A data frame containing the gene IDs.
#' @param gene_id_column The unquoted column name for the gene identifiers.
#' @param uniprot_handle An active `UniProt.ws` handle.
#' @param uniprot_keytype The keytype for the query.
#' @param uniprot_columns A character vector of UniProt columns to retrieve.
#'
#' @return A single data frame containing the combined annotation results from all batches.
#'
#' @importFrom purrr map
#' @importFrom dplyr bind_rows distinct arrange pull
#' @export
batchQueryEvidenceGeneId <- function(input_tbl, gene_id_column, uniprot_handle, uniprot_keytype = "UniProtKB",
                                     uniprot_columns = c("EXISTENCE", "SCORE", "REVIEWED", "GENENAME", "PROTEIN-NAMES", "LENGTH")) {


  all_gene_id <- batchQueryEvidenceHelperGeneId(input_tbl, {{ gene_id_column }})

  partial_subset_query <- partial(subsetQuery,
                                  data = all_gene_id,
                                  accessions_col_name = {{ gene_id_column }},
                                  uniprot_handle = uniprot_handle,
                                  uniprot_columns = uniprot_columns,
                                  uniprot_keytype = uniprot_keytype)

  rounds_list <- all_gene_id |>
    dplyr::distinct(round) |>
    dplyr::arrange(round) |>
    dplyr::pull(round)

  all_uniprot_evidence <- purrr::map(rounds_list, \(x) { partial_subset_query(subset = x) }) |>
    dplyr::bind_rows()

  return(all_uniprot_evidence)
}



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @title Convert GO IDs to GO Terms
#' @description This function takes a string of semicolon-separated Gene Ontology (GO)
#' IDs, splits them, and maps each ID to its corresponding term name and ontology
#' (BP, CC, MF). It then reshapes the output into a wide format with one column per ontology.
#'
#' @param go_string A single character string containing one or more GO IDs separated by `sep`.
#' @param sep The separator used in `go_string`. Defaults to "; ".
#' @param goterms A named vector or list from `AnnotationDbi::Term(GO.db::GOTERM)`.
#' @param gotypes A named vector or list from `AnnotationDbi::Ontology(GO.db::GOTERM)`.
#'
#' @return A tibble with one row and columns for each ontology (`go_biological_process`,
#'   `go_cellular_compartment`, `go_molecular_function`), containing semicolon-separated
#'   GO term names. Returns `NA` if the input string is `NA`.
#'
#' @importFrom tidyr separate_rows pivot_wider
#' @importFrom dplyr mutate group_by summarise ungroup case_when
#' @importFrom purrr map_chr
#' @export
#' @examples
#' \dontrun{
#' goterms <- AnnotationDbi::Term(GO.db::GOTERM)
#' gotypes <- AnnotationDbi::Ontology(GO.db::GOTERM)
#' go_string <- "GO:0016021;GO:0005575"
#' goIdToTerm(go_string, sep = ";", goterms, gotypes)
#' }
goIdToTerm <- function(go_string, sep = "; ", goterms, gotypes) {

  if (!is.na(go_string)) {
    go_string_tbl <- data.frame(go_id = go_string) |>
      separate_rows(go_id, sep = sep) |>
      mutate(go_term = purrr::map_chr(go_id, function(x) { if (x %in% names(goterms)) { return(goterms[[x]]) }; return(NA) })) |>
      mutate(go_type = purrr::map_chr(go_id, function(x) { if (x %in% names(gotypes)) { return(gotypes[[x]]) }; return(NA) })) |>
      group_by(go_type) |>
      summarise(go_term = paste(go_term, collapse = "; ")) |>
      ungroup() |>
      mutate(go_type = case_when(go_type == "BP" ~ "go_biological_process",
                                 go_type == "CC" ~ "go_cellular_compartment",
                                 go_type == "MF" ~ "go_molecular_function")) |>
      pivot_wider(names_from = go_type,
                  values_from = go_term)

    return(go_string_tbl)
  }
  return(NA)
}

# go_string <- "GO:0016021; GO:0030659; GO:0031410; GO:0035915; GO:0042742; GO:0045087; GO:0045335; GO:0050829; GO:0050830"
# go_id_to_term(go_string)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' @title Annotate UniProt Data with GO Terms
#' @description This function takes a data frame of UniProt data containing a column of
#' GO IDs and enriches it with GO term names and IDs, organized by ontology type
#' (BP, CC, MF) in separate columns.
#'
#' @param uniprot_dat A data frame containing UniProt query results.
#' @param uniprot_id_column The unquoted column name for UniProt accessions.
#' @param go_id_column The unquoted column name for the GO IDs (often semicolon-separated).
#' @param sep The separator for the GO IDs.
#' @param goterms A named vector/list from `AnnotationDbi::Term(GO.db::GOTERM)`.
#' @param gotypes A named vector/list from `AnnotationDbi::Ontology(GO.db::GOTERM)`.
#'
#' @return The input `uniprot_dat` data frame with new columns added for each GO
#'   ontology, containing the corresponding term names and IDs (e.g., `go_term_go_biological_process`,
#'   `go_id_go_biological_process`). The original `go_id_column` is removed.
#'
#' @importFrom tidyr separate_rows pivot_wider
#' @importFrom dplyr distinct filter mutate group_by summarise ungroup left_join select
#' @importFrom purrr map_chr
#' @importFrom GO.db GOTERM
#' @importFrom AnnotationDbi Term Ontology
#' @export
uniprotGoIdToTerm <- function(uniprot_dat, uniprot_id_column = UNIPROTKB
                              , go_id_column = `GO-IDs`,  sep = "; "
                              , goterms = AnnotationDbi::Term(GO.db::GOTERM)
                              , gotypes = AnnotationDbi::Ontology(GO.db::GOTERM)) {

  uniprot_acc_to_go_id <- uniprot_dat |>
    dplyr::distinct({{uniprot_id_column}}, {{go_id_column}}) |>
    separate_rows({{go_id_column}}, sep = sep) |>
    dplyr::distinct({{uniprot_id_column}}, {{go_id_column}}) |>
    dplyr::filter(!is.na({{go_id_column}}))

  go_term_temp <- uniprot_acc_to_go_id |>
    dplyr::distinct({{go_id_column}}) |>
    mutate(go_term = purrr::map_chr({{go_id_column}}, function(x) { if (x %in% names(goterms)) { return(goterms[[x]]) }; return(NA) })) |>
    mutate(go_type = purrr::map_chr({{go_id_column}}, function(x) { if (x %in% names(gotypes)) { return(gotypes[[x]]) }; return(NA) })) |>
    mutate(go_type = case_when(go_type == "BP" ~ "go_biological_process",
                               go_type == "CC" ~ "go_cellular_compartment",
                               go_type == "MF" ~ "go_molecular_function"))

  uniprot_acc_to_go_term <- uniprot_acc_to_go_id |>
    left_join(go_term_temp, by = join_by({{go_id_column}} == {{go_id_column}})) |>
    dplyr::filter(!is.na(go_term)) |>
    group_by({{uniprot_id_column}}, go_type) |>
    summarise(go_id = paste({{go_id_column}}, collapse = "; ")
              , go_term = paste(go_term, collapse = "; ")
              , .groups = 'drop' ) |>
    pivot_wider(id_cols = {{uniprot_id_column}},
                names_from = go_type,
                values_from = c(go_term, go_id))


  output_uniprot_dat <- uniprot_dat |>
    left_join(uniprot_acc_to_go_term, by = join_by({{uniprot_id_column}} == {{uniprot_id_column}})) |>
    dplyr::select( -{{go_id_column}})

  return(output_uniprot_dat)

}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @title Create a Long-Format Differential Expression Results Table
#' @description This function combines differential expression statistics with both
#' normalized and raw abundance data for each replicate, creating a comprehensive,
#' long-format table.
#'
#' @param lfc_qval_tbl A data frame containing DE results (logFC, q-value, etc.).
#' @param norm_counts_input_tbl A wide data frame or matrix of normalized abundance data.
#' @param raw_counts_input_tbl A wide data frame or matrix of raw abundance data.
#' @param row_id The unquoted column name for feature identifiers.
#' @param sample_id The unquoted column name for sample identifiers.
#' @param group_id The unquoted column name for group identifiers.
#' @param group_pattern A regular expression to identify the abundance columns.
#' @param design_matrix_norm The design matrix corresponding to the normalized data.
#' @param design_matrix_raw The design matrix corresponding to the raw data.
#' @param expression_column The unquoted column name for the expression/contrast.
#' @param protein_id_table A data frame for joining additional protein information.
#'
#' @return A long-format data frame combining DE stats with replicate-level abundance data.
#'
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer pivot_wider separate_wider_delim
#' @importFrom dplyr select mutate left_join group_by arrange ungroup
#' @importFrom purrr map_chr
#' @importFrom rlang set_names sym
#' @importFrom stringr str_replace_all
#' @export
createDeResultsLongFormat <- function( lfc_qval_tbl,
                                       norm_counts_input_tbl,
                                       raw_counts_input_tbl,
                                       row_id,
                                       sample_id,
                                       group_id,
                                       group_pattern,
                                       design_matrix_norm,
                                       design_matrix_raw,
                                       expression_column = log_intensity,
                                       protein_id_table
) {

  norm_counts <- norm_counts_input_tbl |>
    as.data.frame() |>
    rownames_to_column(row_id) |>
    pivot_longer(cols = matches(group_pattern),
                 names_to = sample_id,
                 values_to = "log2norm") |>
    left_join(design_matrix_norm, by = sample_id) |>
    group_by(!!sym(row_id), !!sym(group_id)) |>
    arrange(!!sym(row_id), !!sym(group_id), !!sym(sample_id)) |>
    mutate(replicate_number = paste0("log2norm.", row_number())) |>
    ungroup() |>
    pivot_wider(id_cols = c(!!sym(row_id), !!sym(group_id)),
                names_from = replicate_number,
                values_from = log2norm) |>
    mutate( {{group_id}} := purrr::map_chr( !!sym(group_id), as.character))


  # print(head(norm_counts))

  raw_counts <- raw_counts_input_tbl |>
    as.data.frame() |>
    rownames_to_column(row_id) |>
    pivot_longer(cols = matches(group_pattern),
                 names_to = sample_id,
                 values_to = "raw") |>
    left_join(design_matrix_raw, by = sample_id) |>
    group_by(!!sym(row_id), !!sym(group_id)) |>
    arrange(!!sym(row_id), !!sym(group_id), !!sym(sample_id)) |>
    mutate(replicate_number = paste0("raw.", row_number())) |>
    ungroup() |>
    pivot_wider(id_cols = c(!!sym(row_id), !!sym(group_id)),
                names_from = replicate_number,
                values_from = raw)  |>
    mutate( {{group_id}} := purrr::map_chr( !!sym(group_id), as.character))

  # print(head(raw_counts))

  left_join_columns <- rlang::set_names(c(row_id, group_id ),
                                        c(row_id, "left_group"))

  right_join_columns <- rlang::set_names(c(row_id, group_id ),
                                         c(row_id, "right_group"))

   # print(head(lfc_qval_tbl))

  print( row_id)
  print(colnames( protein_id_table)[1])

  de_proteins_long <- lfc_qval_tbl |>
    dplyr::select(-lqm, -colour, -analysis_type) |>
    dplyr::mutate( {{expression_column}} := str_replace_all({{expression_column}}, group_id, "")) |>
    separate_wider_delim( {{expression_column}}, delim = "-", names = c("left_group", "right_group"))  |>
   left_join(norm_counts, by = left_join_columns) |>
    left_join(norm_counts, by = right_join_columns,
              suffix = c(".left", ".right")) |>
    left_join(raw_counts, by = left_join_columns) |>
    left_join(raw_counts, by = right_join_columns,
              suffix = c(".left", ".right")) |>
  left_join( protein_id_table
               , by = join_by( !!sym(row_id) == !!sym( colnames( protein_id_table)[1]))) |>
    arrange( comparison, fdr_qvalue, log2FC) |>
    distinct()

  de_proteins_long
}



#' @title Save a ggplot Object with Logging
#' @description A wrapper for `ggplot2::ggsave` that captures any messages or
#' warnings generated during the saving process and logs them using `logger::log_debug`.
#' It can save a plot in multiple formats.
#'
#' @param input_plot The ggplot object to save.
#' @param file_name_part The base path and filename for the output file, without the extension.
#' @param plots_format A character vector of file extensions to save the plot as (e.g., `c(".png", ".pdf")`).
#' @param width The width of the plot.
#' @param height The height of the plot.
#'
#' @return This function is called for its side effect of writing files.
#'
#' @importFrom ggplot2 ggsave
#' @importFrom logger log_debug
#' @export
gg_save_logging <- function( input_plot
                             , file_name_part
                             , plots_format
                             , width=7
                             , height=7) {
  for( format_ext in plots_format) {
    file_name <- paste0(file_name_part, format_ext)
    captured_output<-capture.output(
      ggsave(plot=input_plot
             , filename = file_name
             , width=width
             , height=height )
      ,type = "message"
    )
    logdebug(captured_output)
  }
}




## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @title Calculate Correlation Between Technical Replicates
#' @description This function calculates the Pearson and Spearman correlation coefficients
#' between pairs of technical replicates for each feature.
#'
#' @param design_matrix_tech_rep A design matrix that maps technical replicate sample
#'   IDs to biological sample IDs.
#' @param data_matrix The abundance data matrix (features x samples).
#' @param protein_id_column The column name for protein/feature IDs.
#' @param sample_id_column The column name for sample IDs (technical replicates).
#' @param tech_rep_column The column name identifying the biological sample.
#' @param tech_rep_num_column The column name identifying the replicate number (e.g., 1 or 2).
#' @param tech_rep_remove_regex A regular expression to filter out samples that are not
#'   part of the replicate analysis (e.g., pools).
#'
#' @return A data frame with one row per feature, containing columns for Pearson
#'   and Spearman correlation coefficients.
#'
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom tidyr pivot_longer pivot_wider nest
#' @importFrom dplyr left_join filter select mutate
#' @importFrom purrr map_dbl
#' @importFrom stats cor
#' @export
proteinTechRepCorrelationHelper <- function( design_matrix_tech_rep, data_matrix
                                             , protein_id_column = "Protein.Ids"
                                             , sample_id_column="Sample_ID", tech_rep_column = "replicates", tech_rep_num_column = "tech_rep_num", tech_rep_remove_regex = "pool" ) {

  tech_reps_list <- design_matrix_tech_rep |> dplyr::pull( !!sym(tech_rep_num_column )) |> unique()

  frozen_protein_matrix_tech_rep <- data_matrix  |>
    as.data.frame() |>
    rownames_to_column(protein_id_column) |>
    pivot_longer( cols=!matches( protein_id_column)
                  , values_to = "log2_intensity"
                  , names_to = sample_id_column) |>
    left_join( design_matrix_tech_rep
               , by = join_by( !!sym(sample_id_column) == !!sym(sample_id_column))) |>
    dplyr::filter( !str_detect(  !!sym(tech_rep_column) , tech_rep_remove_regex ) ) |>
    dplyr::select(!!sym( protein_id_column), !!sym(tech_rep_column), log2_intensity, !!sym(tech_rep_num_column)) |>
    dplyr::filter( !!sym(tech_rep_num_column ) %in% tech_reps_list ) |>
    pivot_wider( id_cols = c(!!sym( protein_id_column), !!sym(tech_rep_column))
                 , names_from = !!sym(tech_rep_num_column)
                 , values_from = log2_intensity) |>
    nest( data=!matches(protein_id_column)) |>
    mutate( data = purrr::map( data, \(x){ x |> column_to_rownames(tech_rep_column)} ) ) |>
    mutate( pearson = purrr::map_dbl( data, \(x){  if( length(which(!is.na(x[,1]))) > 0 & length(which(!is.na(x[,2]))) > 0) { cor(x, use="pairwise.complete.obs")[1,2] } else { NA_real_ }   })) |>
    mutate( spearman = purrr::map_dbl( data, \(x){  if( length(which(!is.na(x[,1]))) > 0 & length(which(!is.na(x[,2]))) > 0) { cor(x, use="pairwise.complete.obs", method="spearman")[1,2] } else { NA_real_ }  }))

  frozen_protein_matrix_tech_rep
}

#' @title Download and Process UniProt Annotations
#' @description Downloads protein information from UniProt for a list of protein IDs,
#' processes the results including Gene Ontology annotations, and caches
#' the result for future use.
#'
#' @param input_tbl Data frame containing protein IDs in a column named 'Protein.Ids'.
#' @param cache_dir Directory path for caching the results.
#' @param taxon_id Taxonomic identifier for the organism (e.g., 9606 for human).
#' @param force_download Logical; if TRUE, forces new download even if cache exists.
#' @param batch_size Number of protein IDs to query in each batch.
#' @param timeout Timeout in seconds for the download operation.
#' @param api_delay Sleep time in seconds between API calls.
#'
#' @return A data frame containing UniProt annotations and GO terms.
#'
#' @export
getUniprotAnnotations <- function(input_tbl, 
                                 cache_dir, 
                                 taxon_id,
                                 force_download = FALSE,
                                 batch_size = 25,
                                 timeout = 600,
                                 api_delay = 1) {
  
  # Ensure cache directory exists
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }
  
  # Define cache file paths
  cache_file <- file.path(cache_dir, "uniprot_annotations.RDS")
  raw_results_file <- file.path(cache_dir, "uniprot_results.tsv")
  
  # Check if cache exists and should be used
  if (!force_download && file.exists(cache_file)) {
    message("Loading cached UniProt annotations...")
    return(readRDS(cache_file))
  }
  
  # Download annotations if needed
  message("Fetching UniProt annotations...")
  annotations <- directUniprotDownload(
    input_tbl = input_tbl,
    output_path = raw_results_file,
    taxon_id = taxon_id,
    batch_size = batch_size,
    timeout = timeout,
    api_delay = api_delay
  )
  
  # Process annotations or create empty table if download failed
  if (!is.null(annotations) && nrow(annotations) > 0) {
    message("Processing GO terms...")
    processed_annotations <- annotations |>
      uniprotGoIdToTerm(
        uniprot_id_column = Entry,
        go_id_column = Gene.Ontology.IDs,
        sep = "; "
      )
    
    # Standardize column names
    uniprot_dat_cln <- standardizeUniprotColumns(processed_annotations)
    
    # Save to cache
    saveRDS(uniprot_dat_cln, cache_file)
    message("UniProt annotations saved to cache.")
  } else {
    warning("Failed to retrieve UniProt annotations. Using empty table.")
    uniprot_dat_cln <- createEmptyUniprotTable()
    saveRDS(uniprot_dat_cln, cache_file)
  }
  
  return(uniprot_dat_cln)
}

#' @title Download Protein Data Directly from UniProt REST API
#' @description Downloads protein information from the UniProt REST API for a list of
#' protein IDs. It processes proteins in batches to avoid overwhelming the API.
#'
#' @param input_tbl Data frame containing protein IDs in a column named 'Protein.Ids'.
#' @param output_path File path to save the raw results as a TSV file.
#' @param taxon_id Taxonomic identifier for the organism.
#' @param batch_size Number of protein IDs to query in each batch.
#' @param timeout Timeout in seconds for the download operation.
#' @param api_delay Sleep time in seconds between API calls.
#'
#' @return A data frame containing the raw UniProt results, or `NULL` if the download fails.
#'
#' @importFrom httr GET status_code content timeout
#' @importFrom purrr imap compact reduce
#' @export
directUniprotDownload <- function(input_tbl, 
                                 output_path, 
                                 taxon_id, 
                                 batch_size = 25,
                                 timeout = 600, 
                                 api_delay = 1) {
  # Set a longer timeout
  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout))
  options(timeout = timeout)
  
  # Extract unique protein IDs
  protein_ids <- unique(input_tbl$Protein.Ids)
  message(paste("Found", length(protein_ids), "unique protein IDs to query"))
  
  # Split into batches
  chunks <- split(protein_ids, ceiling(seq_along(protein_ids)/batch_size))
  message(paste("Split into", length(chunks), "chunks for processing"))
  
  # Function to process one chunk
  process_chunk <- function(chunk, chunk_idx, total_chunks) {
    message(paste("Processing chunk", chunk_idx, "of", total_chunks, "with", length(chunk), "IDs"))
    
    # Create query for this batch
    query <- paste0("(", paste(chunk, collapse=" OR "), ") AND organism_id:", taxon_id)
    
    # Use httr to download
    response <- httr::GET(
      url = "https://rest.uniprot.org/uniprotkb/search",
      query = list(
        query = query,
        format = "tsv",
        fields = "accession,id,protein_name,gene_names,organism_name,length,go_id,reviewed"
      ),
      httr::timeout(30)
    )
    
    # Be nice to the API
    Sys.sleep(api_delay)
    
    # Check if successful
    if (httr::status_code(response) == 200) {
      content <- httr::content(response, "text", encoding = "UTF-8")
      temp_file <- tempfile(fileext = ".tsv")
      writeLines(content, temp_file)
      
      chunk_result <- suppressWarnings(
        read.delim(temp_file, sep="\t", quote="", stringsAsFactors=FALSE)
      )
      
      if (nrow(chunk_result) > 0) {
        message(paste("  Found", nrow(chunk_result), "results"))
        return(chunk_result)
      }
    } else {
      message(paste("  Request failed with status", httr::status_code(response)))
    }
    
    return(NULL)
  }
  
  # Process all chunks using imap (provides both value and index)
  total_chunks <- length(chunks)
  results <- purrr::imap(chunks, ~ process_chunk(.x, .y, total_chunks)) |>
    purrr::compact() # Remove NULL results
  
  # Combine results
  if (length(results) > 0) {
    all_results <- purrr::reduce(results, rbind)
    
    # Standardize column names for downstream processing
    names(all_results) <- gsub(" ", ".", names(all_results))
    
    # Add From column needed for downstream processing
    all_results$From <- all_results$Entry
    
    # Write to file
    write.table(all_results, output_path, sep="\t", quote=FALSE, row.names=FALSE)
    message(paste("Successfully retrieved", nrow(all_results), "entries from UniProt"))
    return(all_results)
  } else {
    message("Failed to retrieve any data from UniProt")
    return(NULL)
  }
}

#' @title Standardize UniProt Column Names
#' @description An internal helper to standardize column names from UniProt results
#' for downstream processing. It handles missing columns gracefully by adding them
#' with `NA` values.
#'
#' @param df Data frame with UniProt results.
#'
#' @return Data frame with standardized column names.
#'
#' @keywords internal
standardizeUniprotColumns <- function(df) {
  # Handle Protein existence column
  if ("Protein.existence" %in% colnames(df)) {
    df <- df |> dplyr::rename(Protein_existence = "Protein.existence")
  } else {
    df$Protein_existence <- NA_character_
  }
  
  # Handle Protein names column
  if ("Protein.names" %in% colnames(df)) {
    df <- df |> dplyr::rename(Protein_names = "Protein.names")
  } else {
    df$Protein_names <- NA_character_
  }
  
  # Handle Gene Names column (different versions of API may use different names)
  gene_names_col <- grep("Gene\\.Names", colnames(df), value = TRUE)
  if (length(gene_names_col) > 0) {
    df <- df |> dplyr::rename_with(~"gene_names", .cols = matches("Gene.Names"))
  } else {
    df$gene_names <- NA_character_
  }
  
  return(df)
}

#' @title Create Empty UniProt Table
#' @description An internal helper to create an empty tibble with a standard set of
#' UniProt columns. This is used as a fallback when a UniProt download fails,
#' ensuring downstream functions receive a valid data frame.
#'
#' @return An empty data frame with standard UniProt columns.
#'
#' @keywords internal
createEmptyUniprotTable <- function() {
  data.frame(
    Entry = character(0),
    From = character(0),
    "gene_names" = character(0),
    "Protein_existence" = character(0),
    "Protein_names" = character(0),
    stringsAsFactors = FALSE
  )
}