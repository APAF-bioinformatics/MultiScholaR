# Author(s): Ignatius Pang, Pablo Galaviz
# Email: cmri-bioinformatics@cmri.org.au
# Children’s Medical Research Institute, finding cures for childhood genetic diseases

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'Remove rows in the table where the columns specified by the column regular expression pattern are all zero or NA value.

#'@param input_table Input table with columns recording protein abundances for each sample. The name of these columns matches a regular expression pattern, defined by 'col_pattern'. Remove rows with all samples having no protein abundance.
#'@param col_pattern String representing regular expression pattern that matches the name of columns containing the protein abundance values.
#'@param row_id The column name with the row_id, tidyverse style name.
#'@return A data frame with the rows without abundance values removed.
#' @export
removeEmptyRows <- function(input_table, col_pattern, row_id) {

  temp_col_name <- paste0("temp_", quo_name(enquo(row_id)))

  temp_input_table <- input_table %>%
    dplyr::mutate(!!rlang::sym(temp_col_name) := row_number())

  sites_to_accept <- temp_input_table %>%
    mutate(across(matches(col_pattern, perl = TRUE), ~{ (is.na(.) | . == 0) })) %>%
    dplyr::filter(!if_all(matches(col_pattern, perl = TRUE), ~. == TRUE)) %>%
    dplyr::select({ { temp_col_name } })

  ## Removing entries where all the "Reporter intensity corrected" rows are zero
  filtered_table <- temp_input_table %>%
    inner_join(sites_to_accept, by = temp_col_name) %>%
    dplyr::select(-temp_col_name)

  return(filtered_table)
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' Plot the number of missing values in each sample
#'@param input_table  Data matrix with each row as a protein and each column a sample.
#'@return A ggplot2 bar plot showing the number of missing values per column.
#'@export
plotNumMissingValues <- function(input_table) {

  plot_num_missing_values <- apply(data.matrix(log2(input_table)), 2,
                                   function(x) { length(which(is.infinite(x))) }) %>%
    t %>%
    t %>%
    set_colnames("No. of Missing Values") %>%
    as.data.frame %>%
    rownames_to_column("Samples ID") %>%
    ggplot(aes(x = `Samples ID`, y = `No. of Missing Values`)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90))

  plot_num_missing_values
}




## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' Plot the number of values in each sample
#'@param input_table  Data matrix with each row as a protein and each column a sample.
#'@return A ggplot2 bar plot showing the number of missing values per column.
#'@export
plotNumOfValues <- function(input_table) {

  plot_num_missing_values <- apply(data.matrix(log2(input_table)), 2,
                                   function(x) { length(which(!is.infinite(x))) }) %>%
    t %>%
    t %>%
    set_colnames("No. of Missing Values") %>%
    as.data.frame %>%
    rownames_to_column("Samples ID") %>%
    ggplot(aes(x = `Samples ID`, y = `No. of Missing Values`)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90))

  plot_num_missing_values
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' For each experimental group, identify proteins that have more than accepted number of missing values per group.
#'@param input_table An input table with a column containing the row ID and the rest of the columns representing abundance values for each sample.
#'@param cols A tidyselect command to select the columns. This includes the functions starts_with(), ends_with(), contains(), matches(), and num_range()
#'@param design_matrix A data frame with a column containing the sample ID (as per the sample_id param) and the experimental group (as per the group param). Each row as the sample ID as row name in the data frame.
#'@param sample_id The name of the column in design_matrix table that has the sample ID.
#'@param row_id A unique ID for each row of the 'input_table' variable.
#'@param group_column The name of the column in design_matrix table that has the experimental group.
#'@param max_num_samples_miss_per_group An integer representing the maximum number of samples with missing values per group.
#'@param abundance_threshold Abundance threshold in which the protein in the sample must be above for it to be considered for inclusion into data analysis.
#'@return A list, the name of each element is the sample ID and each element is a vector containing the protein accessions (e.g. row_id) with enough number of values.
#'@export
removeRowsWithMissingValues <- function(input_table, cols, design_matrix, sample_id, row_id, group_column, max_num_samples_miss_per_group, abundance_threshold) {

  abundance_long <- input_table %>%
    pivot_longer(cols = { { cols } },
                 names_to = quo_name(enquo(sample_id)),
                 values_to = "Abundance") %>%
    left_join(design_matrix, by = quo_name(enquo(sample_id)))


  count_values_per_group <- abundance_long %>%
    mutate(has_value = ifelse(!is.na(Abundance) & Abundance > abundance_threshold, 0, 1)) %>%
    group_by({ { row_id } }, { { group_column } }) %>%
    summarise(num_values = sum(has_value)) %>%
    ungroup()


  remove_rows_temp <- count_values_per_group %>%
    dplyr::filter(max_num_samples_miss_per_group < num_values) %>%
    dplyr::select(-num_values, -{ { group_column } }) %>%
    distinct({ { row_id } })

  filtered_tbl <- input_table %>%
    dplyr::anti_join(remove_rows_temp, by = quo_name(enquo(row_id)))

  return(filtered_tbl)

}

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' For each experimental group, identify proteins that does have enough number of samples with abundance values.
#'@param input_table An input table with a column containing the row ID and the rest of the columns representing abundance values for each sample.
#'@param cols A tidyselect command to select the columns. This includes the functions starts_with(), ends_with(), contains(), matches(), and num_range()
#'@param design_matrix A data frame with a column containing the sample ID (as per the sample_id param) and the experimental group (as per the group param). Each row as the sample ID as row name in the data frame.
#'@param sample_id The name of the column in design_matrix table that has the sample ID.
#'@param row_id A unique ID for each row of the 'input_table' variable.
#'@param group_column The name of the column in design_matrix table that has the experimental group.
#'@param min_num_samples_per_group An integer representing the minimum number of samples per group.
#'@param abundance_threshold Abundance threshold in which the protein in the sample must be above for it to be considered for inclusion into data analysis.
#'@return A list, the name of each element is the sample ID and each element is a vector containing the protein accessions (e.g. row_id) with enough number of values.
#'@export
getRowsToKeepList <- function(input_table, cols, design_matrix, sample_id, row_id, group_column, min_num_samples_per_group, abundance_threshold) {

  abundance_long <- input_table %>%
    pivot_longer(cols = { { cols } },
                 names_to = quo_name(enquo(sample_id)),
                 values_to = "Abundance") %>%
    left_join(design_matrix, by = quo_name(enquo(sample_id)))


  count_values_per_group <- abundance_long %>%
    mutate(has_value = ifelse(!is.na(Abundance) & Abundance > abundance_threshold, 1, 0)) %>%
    group_by({ { row_id } }, { { group_column } }) %>%
    summarise(num_values = sum(has_value)) %>%
    ungroup()


  kept_rows_temp <- count_values_per_group %>%
    dplyr::filter(num_values >= min_num_samples_per_group) %>%
    dplyr::select(-num_values) %>%
    group_by({ { group_column } }) %>%
    nest(data = c({ { row_id } })) %>%
    ungroup() %>%
    mutate(data = purrr::map(data, ~{ .[, quo_name(enquo(row_id))][[1]] }))


  sample_rows_lists <- kept_rows_temp$data
  names(sample_rows_lists) <- kept_rows_temp[, quo_name(enquo(group_column))][[1]]

  return(sample_rows_lists)

}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Data imputation function
#'@param df Data matrix
#'@param width Adjustment factor to the observed standard deviation
#'@param downshift Downshift the mean value by this downshift factor multiplied by the observed standard deviation.
#'@return Data matrix with the missing values from each column replaced with a value randomly sampled from the normal distribution with adjusted mean and standard deviation. The normal distribution parameters are based on the observed distribution of the same column.
#'@export
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
#'Converts a design matrix to a biological replicate matrix for use with ruvIII.
#'@param design_matrix The design matrix with the sample ID in one column and the experimental group in another column
#'@param sample_id_column The name of the column with the sample ID, tidyverse style input.
#'@param group_column The name of the column with the experimental group, tidyverse style input.
#'@param temp_column The name of the temporary column that indicates which samples are biological replicates of the same experimental group.
#'@return A numeric matrix with rows as samples, columns as experimental group, and a value of 1 for samples within the same experimental group represented by the same column, and a value of zero otherwise.
#'@export
getRuvIIIReplicateMatrix <- function(design_matrix, sample_id_column, group_column, temp_column = is_replicate_temp) {

  ruvIII_replicates_matrix <- design_mat_cln %>%
    dplyr::select({ { sample_id_column } }, { { group_column } }) %>%
    mutate({ { temp_column } } := 1) %>%
    pivot_wider(id_cols = { { sample_id_column } },
                names_from = { { group_column } },
                values_from = { { temp_column } },
                values_fill = 0) %>%
    column_to_rownames(quo_name(enquo(sample_id_column))) %>%
    as.matrix

  ruvIII_replicates_matrix
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
plotPca <- function(data,
                    design_matrix,
                    sample_id_column = Sample_ID,
                    group_column = group,
                    title,  geom.text.size=11,
                   ...) {

  pca.res <- pca(t(as.matrix(data)))

  temp_tbl <- pca.res$variates$X %>%
    as.data.frame %>%
    rownames_to_column(quo_name(enquo(sample_id_column))) %>%
    left_join(design_matrix, by = quo_name(enquo(sample_id_column)))

  unique_groups <- temp_tbl %>% distinct( {{group_column}}) %>% pull( {{group_column}})

  output <- temp_tbl %>%
    ggplot(aes(PC1, PC2, col = {{group_column}}, label = {{sample_id_column}})) +
    geom_point() +
    geom_text_repel(size  = geom.text.size, show.legend=FALSE) +
    labs(title = title) +
    theme(legend.title = element_blank())

  output
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
plotRle <- function(Y, rowinfo = NULL, probs = c(0.05, 0.25, 0.5, 0.75,
                                                 0.95), ylim = c(-0.5, 0.5))
{
  #  checks = check.ggplot()
  # if (checks) {
  rle = t(apply(t(Y) - apply(Y, 2, median), 2, quantile,
                probs = probs))
  colnames(rle) = c("min", "lower", "middle", "upper",
                    "max")
  df = cbind(data.frame(rle.x.factor = rownames(rle)), data.frame(rle))

  if (!is.null(rowinfo)) {
    rowinfo = data.frame(rowinfo = rowinfo)
    df_temp = cbind(df, rowinfo)

    my.x.factor.levels <- df_temp %>%
      arrange(rowinfo) %>%
      distinct(rle.x.factor) %>%
      pull(rle.x.factor)

    df <- df_temp %>%
      mutate(rle.x.factor = factor(rle.x.factor,
                                   levels = my.x.factor.levels)) %>%
      arrange(rowinfo)
  }

  rleplot = ggplot(df, aes_string(x = "rle.x.factor")) +
    geom_boxplot(aes_string(lower = "lower", middle = "middle",
                            upper = "upper", max = "max", min = "min"),
                 stat = "identity") +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90) #, axis.ticks.x = element_blank()
    ) +
    theme(axis.title.y = element_blank(), axis.text.y = element_text(size = rel(1.5))) +
    geom_hline(yintercept = 0) +
    coord_cartesian(ylim = ylim)
  if (!is.null(rowinfo))
    if (ncol(rowinfo) == 1)
      rleplot = rleplot + aes(fill = rowinfo) + labs(fill = "")
  return(rleplot)
  # }
  # else return(FALSE)
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export
rlePcaPlotList <- function(list_of_data_matrix, list_of_design_matrix,
                           sample_id_column = Sample_ID, group_column = group, list_of_descriptions) {

  rle_list <- purrr::pmap( list( data_matrix=list_of_data_matrix, description=list_of_descriptions, design_matrix=list_of_design_matrix),
                          function( data_matrix, description, design_matrix) { plotRle(t(as.matrix(data_matrix)),
                                   rowinfo = design_matrix[colnames(data_matrix), quo_name(enquo(group_column))]  )  +
                            labs(title = description)} )

  pca_list <- purrr::pmap(list( data_matrix=list_of_data_matrix, description=list_of_descriptions, design_matrix=list_of_design_matrix),
                          function( data_matrix, description, design_matrix) { plotPca(data_matrix,
                                   design_matrix = design_matrix,
                                   sample_id_column = { { sample_id_column } },
                                   group_column = { { group_column } },
                                   title = description, cex = 7) })

  list_of_plots <- c(rle_list, pca_list)

  rle_pca_plots_arranged <- ggarrange(plotlist = list_of_plots, nrow = 2, ncol = length(list_of_descriptions),
                                      common.legend = FALSE, legend = "bottom", widths = 10, heights = 10)

  rle_pca_plots_arranged
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' Count the number of statistically significant differentially expressed proteins (according to user-defined threshold)
#' @param lfc_thresh A numerical value specifying the log fold-change threhold (absolute value) for calling statistically significant proteins.
#' @param q_val_thresh A numerical value specifying the q-value threshold for statistically significant proteins.
#' @param log_fc_column The name of the log fold-change column (tidyverse style).
#' @param q_value_column The name of the q-value column (tidyverse style).
#' @return A table with the following columns:
#' status: The status could be Significant and Up, Significant and Down or Not significant
#' counts: The number of proteins wit this status
#'@export
countStatDeGenes <- function(data,
                             lfc_thresh = 0,
                             q_val_thresh = 0.05,
                             log_fc_column = logFC,
                             q_value_column = q.mod) {

  comparison <- unique(data$comparison)

  selected_data <- data %>%
    dplyr::mutate(status = case_when({ { q_value_column } } >= q_val_thresh ~ "Not significant",
                                     { { log_fc_column } } >= lfc_thresh & { { q_value_column } } < q_val_thresh ~ "Significant and Up",
                                     { { log_fc_column } } < lfc_thresh & { { q_value_column } } < q_val_thresh ~ "Significant and Down",
                                     TRUE ~ "Not significant"))

  counts <- selected_data %>%
    group_by(status) %>%
    summarise(counts = n()) %>%
    ungroup()

  all_possible_status <- data.frame( status = c( "Not significant", "Significant and Up", "Significant and Down"))

  results <- all_possible_status %>%
    left_join( counts, by=c("status" = "status")) %>%
    mutate( counts = ifelse(is.na(counts), 0, counts))

  return(results)
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Format results table for use in volcano plots, counting number of significant proteins, p-values distribution histogram.
#' @param list_of_de_tables A list with each element being a results table with log fold-change and q-value per protein.
#' @param list_of_descriptions  A list of strings describing the parameters used to generate the result table.
#' @param formula_string The formula string used in the facet_grid command for the ggplot scatter plot.
#' @param facet_column The name of the column describing the type of analysis or parameters used to generate the result table (tidyverse style). This is related to the \code{list_of_descriptions} parameter above.
#' @param comparison_column The name of the column describing the contrasts or comparison between groups (tidyverse style).
#' @param expression_column The name of the column that will contain the formula expressions of the contrasts.
#'@export
printCountDeGenesTable <- function(list_of_de_tables,
                                   list_of_descriptions,
                                   formula_string = "analysis_type ~ comparison",
                                   facet_column = analysis_type,
                                   comparison_column = comparison,
                                   expression_column = expression) {

  count_stat_de_genes_helper <- function(de_table, description) {
    purrr::map(de_table, ~countStatDeGenes(.,
                                           lfc_thresh = 0,
                                           q_val_thresh = 0.05,
                                           log_fc_column = logFC,
                                           q_value_column = q.mod)) %>%
      purrr::map2(names(de_table),
                  ~{ .x %>%
                    mutate({ { comparison_column } } := .y) }) %>%
      bind_rows %>%
      mutate({ { facet_column } } := description) %>%
      separate({ { comparison_column } },
               sep = "=",
               into = c(quo_name(enquo(comparison_column)),
                        quo_name(enquo(expression_column))))
  }


  num_significant_de_genes_all <- purrr::map2(list_of_de_tables,
                                              list_of_descriptions,
                                              function(a, b) { count_stat_de_genes_helper(de_table = a,
                                                                                          description = b) }) %>%
    bind_rows()

  num_sig_de_genes_barplot <- num_significant_de_genes_all %>%
    dplyr::filter(status != "Not significant") %>%
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
#' Format results table for use in volcano plots, counting number of significant proteins, p-values distribution histogram.
#' @param list_of_de_tables A list with each element being a results table with log fold-change and q-value per protein.
#' @param list_of_descriptions  A list of strings describing the parameters used to generate the result table.
#' @param row_id The name of the row ID column (tidyverse style).
#' @param p_value_column The name of the raw p-value column (tidyverse style).
#' @param q_value_column The name of the q-value column (tidyverse style).
#' @param log_q_value_column The name of the log q-value column (tidyverse style).
#' @param log_fc_column The name of the log fold-change column (tidyverse style).
#' @param comparison_column The name of the column describing the contrasts or comparison between groups (tidyverse style).
#' @param expression_column The name of the column that will contain the formula expressions of the contrasts.
#' @param facet_column The name of the column describing the type of analysis or parameters used to generate the result table (tidyverse style). This is related to the \code{list_of_descriptions} parameter above.
#' @param q_val_thresh A numerical value specifying the q-value threshold for statistically significant proteins.
#' @return A table with the following columns:
#' row_id:  The protein ID, this column is derived from the input to the row_id column.
#' log_q_value_column: The log (base 10) q-value, this column name is derived from the input to the log_q_value_column.
#' q_value_column:  The q-value, this column name is derived from the input to the q_value_column.
#'   p_value_column: The p-value, this column name is derived from the input to p_value_column.
#'   log_fc_column: The log fold-change, this column name is derived from the input to log_fc_column.
#'   comparison_column: The comparison, this column name is derived from the input to comparison_column.
#'   expression_column: The formula expression for the contrasts, this column name is derived from the input to expression_column.
#'   facet_column: The analysis type, this column name is derived from the input to facet_column.
#'   colour The colour of the dots used in the volcano plot.
#'   orange = Absolute Log fold-change >= 1 and q-value >= threshold
#'   purple = Absolute Log fold-change >= 1 and q-value < threshold
#'   blue = Absolute Log fold-change < 1 and q-value < threshold
#'   black = all other values
#' @export
getSignificantData <- function(list_of_de_tables,
                               list_of_descriptions,
                               row_id = uniprot_acc,
                               p_value_column = p.mod,
                               q_value_column = q.mod,
                               fdr_value_column = fdr.mod,
                               log_q_value_column = lqm,
                               log_fc_column = logFC,
                               comparison_column = comparison,
                               expression_column = expression,
                               facet_column = analysis_type,
                               q_val_thresh = 0.05) {

  get_row_binded_table <- function(de_table_list, description) {
    output <- purrr::map(de_table_list,
                         function(tbl) { tbl %>%
                           rownames_to_column(quo_name(enquo(row_id))) %>%
                           dplyr::select({ { row_id } },
                                         { { p_value_column } },
                                         { { q_value_column } },
                                         { { fdr_value_column } },
                                         { { log_fc_column } }) }) %>%
      purrr::map2(names(de_table_list), ~{ .x %>%
        mutate({ { comparison_column } } := .y) }) %>%
      bind_rows() %>%
      mutate({ { facet_column } } := description) %>%
      separate({ { comparison_column } },
               sep = "=",
               into = c(quo_name(enquo(comparison_column)),
                        quo_name(enquo(expression_column))))

  }

  logfc_tbl_all <- purrr::map2(list_of_de_tables, list_of_descriptions,
                               function(a, b) { get_row_binded_table(de_table_list = a, description = b) }) %>%
    bind_rows()

  selected_data <- logfc_tbl_all %>%
    mutate({ { log_q_value_column } } := -log10(q.mod)) %>%
    dplyr::select({ { row_id } }, { { log_q_value_column } }, { { q_value_column } }, { { p_value_column } }, { { log_fc_column } },
                  { { comparison_column } }, { { expression_column } },
                  { { facet_column } }) %>%
    dplyr::mutate(colour = case_when(abs({ { log_fc_column } }) >= 1 & { { q_value_column } } >= q_val_thresh ~ "orange",
                                     abs({ { log_fc_column } }) >= 1 & { { q_value_column } } < q_val_thresh ~ "purple",
                                     abs({ { log_fc_column } }) < 1 & { { q_value_column } } < q_val_thresh ~ "blue",
                                     TRUE ~ "black")) %>%
    dplyr::mutate(colour = factor(colour, levels = c("black", "orange", "blue", "purple")))

  selected_data

}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' Draw the volcano plot.
#' @param selected_data A table that is generated by running the function \code{\link{get_significant_data}}.
#' @param log_q_value_column The name of the column representing the log q-value.
#' @param log_fc_column The name of the column representing the log fold-change.
#' @param q_val_thresh A numerical value specifying the q-value threshold for statistically significant proteins.
#' @param formula_string The formula string used in the facet_grid command for the ggplot scatter plot.
#'@export
plotVolcano <- function(selected_data,
                        log_q_value_column = lqm,
                        log_fc_column = logFC,
                        q_val_thresh = 0.05,
                        formula_string = "analysis_type ~ comparison") {

  volplot_gg.all <- selected_data %>%
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
    geom_vline(xintercept = 1, colour = "black", size = 0.2) +
    geom_vline(xintercept = -1, colour = "black", size = 0.2) +
    geom_hline(yintercept = -log10(q_val_thresh)) +
    theme_bw() +
    xlab("Log fold changes") +
    ylab("-log10 q-value") +
    theme(legend.position = "none") +
    facet_grid( as.formula(formula_string),
               labeller = labeller(facet_category = label_wrap_gen(width = 10)))

  volplot_gg.all
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' Draw the p-values distribution plot.
#' @param selected_data A table that is generated by running the function \code{\link{get_significant_data}}.
#' @param log_p_value_column The name of the column representing the p-value.
#' @param formula_string The formula string used in the facet_grid command for the ggplot scatter plot.
#'@export
printPValuesDistribution <- function(selected_data, p_value_column = p.mod, formula_string = "is_ruv_applied ~ comparison") {

  breaks <- c(0, 0.001, 0.01, 0.05,
              seq(0.1, 1, by = 0.1))


  pvalhist <- ggplot(selected_data, aes({ { p_value_column } })) +
    theme(axis.title.y = element_blank()) +
    xlab("P-value") +
    geom_histogram(aes_string(y = "..density.."),
                   breaks = breaks,
                   position = "identity",
                   color = "black") +
    geom_histogram(aes_string(y = "..density.."),
                   breaks = breaks,
                   position = "identity")

  if (!is.na(formula_string)) {

    pvalhist <- pvalhist +
      facet_grid(as.formula(formula_string))
  }


  pvalhist

}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' Run the Empircal Bayes Statistics for Differential Expression in the limma package
#' @param ID List of protein accessions / row names.
#' @param design Output from running the function \code{\link{model.matrix}}.
#' @param contr.matrix Output from the function \code{\link{makeContrasts}}.
#' @seealso \code{\link{model.matrix}}
#' @seealso \code{\link{makeContrasts}}
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
  p.mod <- fit.eb$p.value[, 1]
  q.ord <- qvalue(p.ord)$q
  q.mod <- qvalue(p.mod)$q

  return(list(table = data.frame(logFC, t.ord, t.mod, p.ord, p.mod, q.ord, q.mod, df.r, df.0, s2.0, s2, s2.post),
              fit.eb = fit.eb))
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Analyse one contrast (e.g. compare a pair of experimental groups) and output the q-values per protein.
#'@param ID List of protein accessions / row names.
#'@param A String representing the name of experimental group A for pairwise comparison of B - A.
#'@param B String representing the name of experimental group B for pairwise comparison of B - A.
#'@param group_A Names of all the columns / samples that are in experimental group A.
#'@param group_B Names of all the columns / samples that are in experimental group B.
#'@param design_matrix A data frame with a column containing the sample ID (as per the sample_id param) and the experimental group (as per the group param). Each row as the sample ID as row name in the data frame.
#'@param formula_string A formula string representing the experimental design. e.g. ("~ 0 + group")
#'@param contrast_variable String representing the contrast variable, which is also used in the formula string. (e.g. "group")
#'@param weights Numeric matrix for adjusting each sample and gene.
#'@return A data frame with the following columns:
#' row.names = the protein accessions
#' comparison A string showing log({group B's name}) minus log({group A's name})
#' meanA     mean of the normalized log abundance value of the gene across samples from experimental group A
#' meanB     mean of the normalized log abundance value of the gene across samples from experimental group B
#' logFC     log fold-change
#' tstats    t-test statistics
#' tmod      moderated t-test statistics
#' pval      t-test p-value
#' pmod      moderated t-test p-value
#' qval      t-test q-value
#' q.mod     moderated t-test q-value
#'@export
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
                                 pmod = r$p.mod,
                                 qval = r$q.ord,
                                 q.mod = r$q.mod),
              fit.eb = fit.eb))
}

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'Assign experimental group list
#'@param design_matrix A data frame representing the design matrix.
#'@param group_id A string representing the name of the group ID column used in the design matrix.
#'@param sample_id A string representing the name of the sample ID column used in the design matrix.
#'@return A list where each element name is the name of a treatment group and each element is a vector containing the sample IDs within the treatment group.
#'@export
getTypeOfGrouping <- function(design_matrix, group_id, sample_id) {
  temp_type_of_grouping <- design_matrix %>%
    dplyr::select(!!rlang::sym(group_id), !!rlang::sym(sample_id)) %>%
    group_by(!!rlang::sym(group_id)) %>%
    summarise(!!rlang::sym(sample_id) := list(!!rlang::sym(sample_id))) %>%
    ungroup

  type_of_grouping <- temp_type_of_grouping %>% pull(!!rlang::sym(sample_id))
  names(type_of_grouping) <- temp_type_of_grouping %>% pull(!!rlang::sym(group_id))

  return(type_of_grouping)

}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' Compare a pair of experimental groups and output the log fold-change and q-values per protein.
#'@param ID List of protein accessions / row names.
#'@param data Data frame containing the log (base 2) protein abundance values where each column represents a sample and each row represents a protein group, and proteins as rows. The data is preferably median-scaled, with missing values imputed, and batch-effects removed.
#'@param test_pairs Input file with a table listing all the pairs of experimental groups to compare. First column represents group A and second column represents group B. Linear model comparisons (e.g. Contrasts) would be group B minus group A.
#'@param sample_columns A vector of column names (e.g. strings) representing samples which would be used in the statistical tests. Each column contains protein abundance values.
#'@param sample_rows_list A list, the name of each element is the sample ID and each element is a vector containing the protein accessions (e.g. row_id) with enough number of values. It is usually the output from the function \code{get_rows_to_keep_list}.
#'@param type_of_grouping A list where each element name is the name of a treatment group and each element is a vector containing the sample IDs within the treatment group. It is usually the output from the function \code{get_type_of_grouping}.
#'@param design_matrix A data frame with a column containing the sample ID (as per the sample_id param) and the experimental group (as per the group param). Each row as the sample ID as row name in the data frame.
#'@param formula_string A formula string representing the experimental design. e.g. ("~ 0 + group")
#'@param contrast_variable String representing the contrast variable, which is also used in the formula string. (e.g. "group")
#'@param weights Numeric matrix for adjusting each sample and gene.
#'@return A list of data frames, the name of each element represents each pairwise comparison. Each data frame has the following columns:
#' row.names = the protein accessions
#' comparison A string showing log({group B's name}) minus log({group A's name})
#' meanA     mean of the normalized log abundance value of the gene across samples from experimental group A
#' meanB     mean of the normalized log abundance value of the gene across samples from experimental group B
#' logFC     log fold-change
#' tstats    t-test statistics
#' tmod      moderated t-test statistics
#' pval      t-test p-value
#' pmod      moderated t-test p-value
#' qval      t-test q-value
#' q.mod     moderated t-test q-value
#' @seealso \code{\link{get_rows_to_keep_list}}
#' @seealso \code{\link{get_type_of_grouping}}
#'@export
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

#' Run the linear model fitting and statistical tests for a set of contrasts, then adjust with Empirical Bayes function
#'@param data Data frame containing the log (base 2) protein abundance values where each column represents a sample and each row represents a protein group, and proteins as rows. The data is preferably median-scaled, with missing values imputed, and batch-effects removed.
#'@param contrast_strings Input file with a table listing all the experimental contrasts to analyse. It will be in the format required for the function \code{makeContrasts} in the limma package.
#'The contrast string consists of variable that each consist of concatenating the column name (e.g. group) and the string representing the group type (e.g. A) in the design matrix.
#'@param design_matrix A data frame with a column containing the sample ID (as per the sample_id param) and the experimental group (as per the group param). Each row as the sample ID as row name in the data frame.
#'@param formula_string A formula string representing the experimental design. e.g. ("~ 0 + group")
#' @param p_value_column The name of the raw p-value column (tidyverse style).
#' @param q_value_column The name of the q-value column (tidyverse style).
#' @param fdr_value_column The name of the fdr-value column (tidyverse style).
#'@return A list containing two elements. $results returns a list of tables containing logFC and q-values. $fit.eb returns the Empiracle Bayes output object.
#'@export
runTestsContrasts <- function(data,
                              contrast_strings,
                              design_matrix,
                              formula_string,
                              p_value_column = p.mod,
                              q_value_column = q.mod,
                              fdr_value_column = fdr.mod,
                              weights = NA,
                              treat_lfc_cutoff = NA,
                              eBayes_trend = FALSE,
                              eBayes_robust = FALSE) {

  ff <- as.formula(formula_string)
  mod_frame <- model.frame(ff, design_matrix)
  design_m <- model.matrix(ff, mod_frame)

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

  fit <- lmFit(data, design = design_m)

  cfit <- contrasts.fit(fit, contrasts = contr.matrix)

  eb.fit <- eBayes( cfit, trend = eBayes_trend, robust = eBayes_robust )

  ## Run treat over here
  t.fit <- NA
  if( !is.na( treat_lfc_cutoff)) {
    t.fit <- treat(eb.fit, lfc=as.double(treat_lfc_cutoff)) ## assign log fold change threshold below which is scientifically not relevant
  } else {
    t.fit <- eb.fit
  }

  result_tables <- purrr::map(contrast_strings,
                              function(contrast) {
                                de_tbl <- topTreat(t.fit, coef = contrast, n = Inf) %>%
                                  mutate({ { q_value_column } } := qvalue(P.Value)$q) %>%
                                  mutate({ { fdr_value_column } } := p.adjust(P.Value, method="BH")) %>%

                                  dplyr::rename({ { p_value_column } } := P.Value)
                              }
  )

  names(result_tables) <- contrast_strings


  return(list(results = result_tables,
              fit.eb = t.fit))
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## edited from missMethyl source code https://rdrr.io/bioc/missMethyl/src/R/RUVfunctions.R
#'@export
cmriRUVfit <- function(Y, X, ctl, Z = 1, k = NULL, method = c("inv", "rinv",
                                                              "ruv4", "ruv3", "ruv2"), M = NULL, ...)
{
  method <- match.arg(method)
  if ((method %in% c("ruv4", "ruv3", "ruv2")) & is.null(k))
    stop("'k' cannot be NULL if method is 'ruv4', 'ruv3' or 'ruv2'.")
  if (mode(ctl) != "logical")
    stop("'ctl' must be a logical vector.")
  if (is.data.frame(Y))
    Y <- data.matrix(Y)
  if (mode(Y) != "numeric")
    stop("'Y' must be a numeric matrix.")
  if (method == "ruv3" & is.null(M))
    stop("'M' cannot be NULL if method is 'ruv3'")
  Y <- t(Y)
  fit <- switch(method, inv = ruv::RUVinv(Y = Y, X = X,
                                          ctl = ctl, Z = Z, ...),
                rinv = ruv::RUVrinv(Y = Y, X = X, ctl = ctl,
                                    Z = Z, k = k, ...),
                ruv4 = ruv::RUV4(Y = Y, X = X, ctl = ctl,
                                 k = k, Z = Z, ...),
                ruv2 = ruv::RUV2(Y = Y, X = X, ctl = ctl,
                                 k = k, Z = Z, ...),
                ruv3 = ruv::RUVIII(Y = Y, M = M, ctl = ctl, k = k,
                                   ...))
  return(fit)
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
extractRuvResults <- function(results_list) {

  extracted <- purrr::map(results_list, ~{ .$results })

  names(extracted) <- names(results_list)

  return(extracted)
}


#'@export
extractResults <- function(results_list) {

  extracted <- purrr::map(results_list, ~{ .$results })

  names(extracted) <- names(results_list)

  return(extracted)
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'Save the list of output tables from differential expression analysis of proteins or phosphopeptides into a file and in a specific directory.
#'@param list_of_de_tables A list, each element is a table of log fold-change and q-values from differential expression analysis of proteins / phosphopeptides. Each element in the list has a name, usually the name of the pairwise comparison.
#'@param row_id Add row ID to the output table based on the name (protein or phosphopeptid ID) of each row
#'@param sort_by_column Each table in the list_of_de_tables is sorted in ascending order
#'@param results_dir The results directory to store the output file
#'@param file_suffix The file suffix string to aadd to the name of each comparison from the list_of_de_tables.
#'@export
saveDeProteinList <- function(list_of_de_tables, row_id, sort_by_column = q.mod, results_dir, file_suffix) {

  purrr::walk2(list_of_de_tables, names(list_of_de_tables),

               ~vroom::vroom_write(.x %>%
                                     rownames_to_column(row_id) %>%
                                     arrange({ { sort_by_column } }),
                                   path = file.path(results_dir, paste0(.y, file_suffix))))

}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
analyseRanking <- function(data, uniprot_acc_column = uniprot_acc) {

  results_tbl <- data %>%
    as.data.frame %>%
    rownames_to_column(quo_name(enquo(uniprot_acc_column))) %>%
    dplyr::select(one_of(c(quo_name(enquo(uniprot_acc_column)), "q.mod", "logFC"))) %>%
    arrange(desc(q.mod)) %>%
    mutate(ctrl_gene_rank = row_number())

  return(results_tbl)
}

#'@export
getControlGenes <- function(data,
                            q.value = 0.05,
                            logFC_threshold = 1,
                            uniprot_acc_column = uniprot_acc) {

  temp <- purrr::map(data, ~analyseRanking(., uniprot_acc_column = { { uniprot_acc_column } }))

  ctrl_genes_list <- temp %>%
    bind_rows(.id = "Test") %>%
    dplyr::filter(q.mod >= q.value &
                    abs(logFC) <= logFC_threshold) %>%
    group_by({ { uniprot_acc_column } }) %>%
    summarise(total_num_set = n()) %>%
    ungroup() %>%
    dplyr::filter(total_num_set == length(data)) %>%
    dplyr::select({ { uniprot_acc_column } })

  return(ctrl_genes_list)

}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' Identify negative control proteins for use in removal of unwanted variation, using an ANOVA test.
#' @param data_matrix A matrix containing the log (base 2) protein abundance values where each column represents a sample and each row represents a protein group, and proteins as rows. The row ID are the protein accessions. The data is preferably median-scaled with missing values imputed.
#' @param group_column The name of the column with the experimental group, as a string.
#' @param num_neg_ctrl The number of negative control genes to select. Typically the number of genes with the highest q-value (e.g. least statistically significant).
#' @param q_val_thresh The q-value threshold. No proteins with q-values lower than this value are included in the list of negative control proteins. This means the number of negative control proteins could be less than the number specified in \code{num_neg_ctrl} when some were excluded by this threshold.
#' @return A boolean vector which indicates which row in the input data matrix is a control gene. The row is included if the value is TRUE. The names of each element is the row ID / protein accessions of the input data matrix.
#'@export
getNegCtrlProtAnova <- function(data_matrix, design_matrix, group_column = "group", num_neg_ctrl = 500, q_val_thresh = 0.05) {

  ## Inspired by matANOVA function from PhosR package: http://www.bioconductor.org/packages/release/bioc/html/PhosR.html

  grps <- design_matrix[colnames(data_matrix), group_column]

  ps <- apply(data_matrix, 1, function(x) {
    summary(stats::aov(as.numeric(x) ~ grps))[[1]][["Pr(>F)"]][1]
  })

  aov <- qvalue(ps)$qvalues

  filtered_list <- aov[aov > q_val_thresh]

  list_size <- ifelse(num_neg_ctrl > length(filtered_list), length(filtered_list), num_neg_ctrl)

  control_genes <- names(sort(filtered_list, decreasing = TRUE)[1:list_size])

  #nrow(data_matrix) - length(control_genes)
  control_genes_index <- rownames(data_matrix) %in% control_genes
  names(control_genes_index) <- rownames(data_matrix)

  return(control_genes_index)

}

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@param design_matrix Contains the sample_id column and the average_replicates_id column
#'@export
averageValuesFromReplicates <- function(input_table, design_matrix, row_id, sample_id, average_replicates_id) {

  output_table <- input_table %>%
    as.data.frame %>%
    rownames_to_column(  row_id   ) %>%
    pivot_longer( cols=colnames(input_table),
                  names_to = sample_id,
                  values_to = "value") %>%
    left_join( design_matrix, by = sample_id ) %>%
    group_by( !!rlang::sym(average_replicates_id) ,  !!rlang::sym(row_id) ) %>%
    summarise( value = mean(value, na.rm = TRUE)) %>%
    ungroup %>%
    pivot_wider( names_from = !!rlang::sym(average_replicates_id),
                 values_from = "value") %>%
    column_to_rownames(row_id) %>%
    as.matrix

  return( output_table )
}



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# up <- UniProt.ws(taxId=10090)
# keytypes(up)
# columns(up)
# test <- batch_query_evidence(subset_tbl, Proteins)

#'@export
cleanIsoformNumber <- function(string) {
  # "Q8K4R4-2"
  str_replace(string, "-\\d+$", "")

}


# Filter for a batch and run analysis on that batch of uniprot accession keys only.
subsetQuery <- function(data, subset, accessions_col_name, uniprot_handle, uniprot_columns = c("EXISTENCE", "SCORE", "REVIEWED", "GENENAME", "PROTEIN-NAMES", "LENGTH"),
                        uniprot_keytype = "UNIPROTKB") {


  # print(subset)
  my_keys <- data %>%
    dplyr::filter(round == subset) %>%
    pull({ { accessions_col_name } })

   # print(my_keys)

  # print(uniprot_keytype)

  UniProt.ws::select(up,
                     keys = my_keys,
                     columns = uniprot_columns,
                     keytype = uniprot_keytype)
}


# The UniProt.ws::select function limits the number of keys queried to 100. This gives a batch number for it to be queried in batches.
batchQueryEvidenceHelper <- function(uniprot_acc_tbl, uniprot_acc_column) {

  all_uniprot_acc <- uniprot_acc_tbl %>%
    dplyr::select({ { uniprot_acc_column } }) %>%
    mutate(Proteins = str_split({ { uniprot_acc_column } }, ";")) %>%
    unnest(Proteins) %>%
    distinct %>%
    arrange(Proteins) %>%
    mutate(Proteins = cleanIsoformNumber(Proteins)) %>%
    dplyr::mutate(round = ceiling(row_number() / 100))  ## 100 is the maximum number of queries at one time
}

## Run evidence collection online, giving a table of keys (uniprot_acc_tbl) and the column name (uniprot_acc_column)
#'@export
batchQueryEvidence <- function(uniprot_acc_tbl, uniprot_acc_column, uniprot_handle,
                               uniprot_columns = c("EXISTENCE", "SCORE", "REVIEWED", "GENENAME", "PROTEIN-NAMES", "LENGTH"),
                               uniprot_keytype = "UNIPROTKB") {

  # uniprot_evidence_levels <- c("Evidence at protein level",
  #                              "Evidence at transcript level",
  #                              "Inferred from homology",
  #                              "Predicted",
  #                              "Uncertain",
  #                              NA)

  all_uniprot_acc <- batchQueryEvidenceHelper(uniprot_acc_tbl,
                                              { { uniprot_acc_column } })

  partial_subset_query <- partial(subsetQuery,
                                  data = all_uniprot_acc,
                                  accessions_col_name = { { uniprot_acc_column } },
                                  uniprot_handle = uniprot_handle,
                                  uniprot_columns = uniprot_columns,
                                  uniprot_keytype = uniprot_keytype)

  rounds_list <- all_uniprot_acc %>%
    distinct(round) %>%
    arrange(round) %>%
    pull(round)

  all_uniprot_evidence <- purrr::map(rounds_list, ~{ partial_subset_query(subset = .) }) %>%
    bind_rows

  return(all_uniprot_evidence)
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# The UniProt.ws::select function limits the number of keys queried to 100. This gives a batch number for it to be queried in batches.
batchQueryEvidenceHelperGeneId <- function(input_tbl, gene_id_column) {

  all_uniprot_acc <- input_tbl %>%
    dplyr::select({ { gene_id_column } }) %>%
    arrange({ { gene_id_column } }) %>%
    dplyr::mutate(round = ceiling(row_number() / 100))  ## 100 is the maximum number of queries at one time
}

batchQueryEvidenceGeneId <- function(input_tbl, gene_id_column, uniprot_handle,
                               uniprot_columns = c("EXISTENCE", "SCORE", "REVIEWED", "GENENAME", "PROTEIN-NAMES", "LENGTH")) {


  all_gene_id <- batchQueryEvidenceHelperGeneId(input_tbl,
                                              { { gene_id_column } })


  partial_subset_query <- partial(subsetQuery,
                                  data = all_gene_id,
                                  accessions_col_name = { { gene_id_column } },
                                  uniprot_handle = uniprot_handle,
                                  uniprot_columns = uniprot_columns,
                                  uniprot_keytype = "GENENAME")

  rounds_list <- all_gene_id %>%
    distinct(round) %>%
    arrange(round) %>%
    pull(round)

  all_uniprot_evidence <- purrr::map(rounds_list, ~{ partial_subset_query(subset = .) }) %>%
    bind_rows

  return(all_uniprot_evidence)
}



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Convert a list of Gene Ontology IDs to their respective human readable name (e.g. GO term).
#' @param go_string A string consisting of a list of Gene Ontology ID, separated by a delimiter
#' @param goterms Output from running \code{goterms <- Term(GOTERM)} from the GO.db library.
#' @param gotypes Output from running \code{gotypes <- Ontology(GOTERM)} from the GO.db library.
#' @return A table with three columns. go_biological_process, go_celluar_compartment, and go_molecular_function. Each column is a list of gene ontology terms, separated by '; '.
#' @export
#' @examples
#' go_string <- "GO:0016021; GO:0030659; GO:0031410; GO:0035915; GO:0042742; GO:0045087; GO:0045335; GO:0050829; GO:0050830"
#' go_id_to_term(go_string)
goIdToTerm <- function(go_string, sep = "; ", goterms, gotypes) {

  if (!is.na(go_string)) {
    go_string_tbl <- data.frame(go_id = go_string) %>%
      separate_rows(go_id, sep = sep) %>%
      mutate(go_term = purrr::map_chr(go_id, function(x) { if (x %in% names(goterms)) { return(goterms[[x]]) }; return(NA) })) %>%
      mutate(go_type = purrr::map_chr(go_id, function(x) { if (x %in% names(gotypes)) { return(gotypes[[x]]) }; return(NA) })) %>%
      group_by(go_type) %>%
      summarise(go_term = paste(go_term, collapse = "; ")) %>%
      ungroup() %>%
      mutate(go_type = case_when(go_type == "BP" ~ "go_biological_process",
                                 go_type == "CC" ~ "go_cellular_compartment",
                                 go_type == "MF" ~ "go_molecular_function")) %>%
      pivot_wider(names_from = go_type,
                  values_from = go_term)

    return(go_string_tbl)
  }
  return(NA)
}

# go_string <- "GO:0016021; GO:0030659; GO:0031410; GO:0035915; GO:0042742; GO:0045087; GO:0045335; GO:0050829; GO:0050830"
# go_id_to_term(go_string)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' @param uniprot_dat  a table with uniprot accessions and a column with GO-ID
#' @param goterms Output from running \code{goterms <- Term(GOTERM)} from the GO.db library.
#' @param gotypes Output from running \code{gotypes <- Ontology(GOTERM)} from the GO.db library.
#' @return A table with three columns. go_biological_process, go_celluar_compartment, and go_molecular_function. Each column is a list of gene ontology terms, separated by '; '.
#' @export
uniprotGoIdToTerm <- function(uniprot_dat, sep = "; ", goterms, gotypes) {


  uniprot_acc_to_go_id <- uniprot_dat %>%
    dplyr::distinct(UNIPROTKB, `GO-ID`) %>%
    separate_rows(`GO-ID`, sep = sep) %>%
    dplyr::distinct(UNIPROTKB, `GO-ID`) %>%
    dplyr::filter(!is.na(`GO-ID`))

  go_term_temp <- uniprot_acc_to_go_id %>%
    dplyr::distinct(`GO-ID`) %>%
    mutate(go_term = purrr::map_chr(`GO-ID`, function(x) { if (x %in% names(goterms)) { return(goterms[[x]]) }; return(NA) })) %>%
    mutate(go_type = purrr::map_chr(`GO-ID`, function(x) { if (x %in% names(gotypes)) { return(gotypes[[x]]) }; return(NA) })) %>%
    mutate(go_type = case_when(go_type == "BP" ~ "go_biological_process",
                               go_type == "CC" ~ "go_cellular_compartment",
                               go_type == "MF" ~ "go_molecular_function"))

  uniprot_acc_to_go_term <- uniprot_acc_to_go_id %>%
    left_join(go_term_temp, by = c("GO-ID" = "GO-ID")) %>%
    dplyr::filter(!is.na(go_term)) %>%
    group_by(UNIPROTKB, go_type) %>%
    summarise(go_term = paste(go_term, collapse = "; ")) %>%
    ungroup() %>%
    pivot_wider(id_cols = "UNIPROTKB",
                names_from = go_type,
                values_from = go_term)


  output_uniprot_dat <- uniprot_dat %>%
    left_join(uniprot_acc_to_go_term, by = c("UNIPROTKB" = "UNIPROTKB")) %>%
    relocate(KEYWORDS, .before = "GO-ID")

  return(output_uniprot_dat)

}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'Choose the best accession
#'@param input_tbl Contain the following columns, 'group_id' which is the Id for each protein group, 'accessions_column' which is the column with the accession of the protein
#'@param acc_detail_tabl The out table from running the function 'parse_fasta_file'
#'@param accessions_column The name of the column with the list of protein accessions, separated by ';' semi-colon. No need to quote the name as we are using tidyverse programming quosure.
#'@param group_id The name of the column with the group ID for each protein group. No need to quote the name as we are using tidyverse programming quosure.
#'@export
chooseBestProteinAccession <- function(input_tbl, acc_detail_tab, accessions_column, row_id_column = uniprot_acc, group_id) {


  join_condition <- rlang::set_names(c(quo_name(enquo(row_id_column)), "cleaned_acc"),
                                     c(quo_name(enquo(row_id_column)), "cleaned_acc"))

  resolve_acc_helper <- input_tbl %>%
    dplyr::select({ { group_id } }, { { accessions_column } }) %>%
    mutate({ { row_id_column } } := str_split({ { accessions_column } }, ";")) %>%
    unnest({ { row_id_column } }) %>%
    mutate(cleaned_acc = cleanIsoformNumber({ { accessions_column } })) %>%
    left_join(acc_detail_tab,
              by = join_condition) %>%
    dplyr::select({ { group_id } }, one_of(c(quo_name(enquo(row_id_column)), "gene_name", "cleaned_acc",
                                             "protein_evidence", "status", "is_isoform", "isoform_num", "seq_length"))) %>%
    distinct %>%
    arrange({ { group_id } }, protein_evidence, status, is_isoform, desc(seq_length), isoform_num)

  score_isoforms <- resolve_acc_helper %>%
    mutate(gene_name = ifelse(is.na(gene_name) | gene_name == "", "NA", gene_name)) %>%
    group_by({ { group_id } }, gene_name) %>%
    arrange({ { group_id } }, protein_evidence,
            status, is_isoform, desc(seq_length), isoform_num, cleaned_acc) %>%
    mutate(ranking = row_number()) %>%
    ungroup


  ## For each gene name find the uniprot_acc with the lowest ranking

  my_group_id <- enquo(group_id)

  join_names <- rlang::set_names(c(quo_name(my_group_id), "ranking", "gene_name"),
                                 c(quo_name(my_group_id), "ranking", "gene_name"))

  group_gene_names_and_uniprot_accs <- score_isoforms %>%
    distinct({ { group_id } }, gene_name, ranking) %>%
    dplyr::filter(ranking == 1) %>%
    left_join(score_isoforms %>%
                dplyr::select({ { group_id } }, ranking, gene_name, uniprot_acc),
              by = join_names) %>%
    dplyr::select(-ranking) %>%
    group_by({ { group_id } }) %>%
    summarise(num_gene_names = n(),
              gene_names = paste(gene_name, collapse = ":"),
              uniprot_acc = paste(uniprot_acc, collapse = ":")) %>%
    ungroup() %>%
    mutate(is_unique = case_when(num_gene_names == 1 ~ "Unique",
                                 TRUE ~ "Multimapped"))


  return(group_gene_names_and_uniprot_accs)

}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
cmriCamera <- function(contrast_name, index_name, abundance_mat, replicates_mat, lists_of_contrasts, list_of_gene_sets, min_set_size = 4) {

  print(paste("contrast_name =", contrast_name))
  print(paste("index_name =", index_name))

  this_contrast <- lists_of_contrasts[[contrast_name]]

  groupA <- str_replace( names( this_contrast[this_contrast == 1] )[1], paste0("^", quo_name(enquo(group_id))), "" )
  groupB <- str_replace( names( this_contrast[this_contrast == -1] )[1], paste0("^", quo_name(enquo(group_id))), "" )

  # print(paste(groupA, groupB))

  design_choose_column <- replicates_mat[, c(groupA, groupB)]

  design_trimmed <- design_choose_column[rowSums(design_choose_column) > 0,]

  # print(design_trimmed)
  #
  # print(head( abundance_mat[[contrast_name]][ , rownames(design_trimmed)]))

  abundance_mat_trimmed <- abundance_mat %>%
    dplyr::filter( comparison == contrast_name) %>%
    dplyr::pull( data) %>%
    .[[1]] %>%
    .[,rownames(design_trimmed)]

  contrast_mat_trimmed <-this_contrast[colnames(replicates_mat) %in% colnames(design_trimmed)]

  index <- list_of_gene_sets[[index_name]]

  msigdb_ids <- geneIds(index)

  #convert gene sets into a list of gene indices
  camera_indices <- ids2indices(msigdb_ids,
                                rownames(abundance_mat_trimmed))

  ## At least two genes in the gene set
  camera_indices_filt <- camera_indices[purrr::map(camera_indices, length) >= min_set_size]


  camera_result <- NA
  if (length(camera_indices_filt) > 0) {
    camera_result <- camera(y = abundance_mat_trimmed, design = design_trimmed, index = camera_indices_filt, contrast = contrast_mat_trimmed)
  }

  info_list <- list(camera = camera_result,
                    y = abundance_mat_trimmed,
                    design = design_trimmed,

                    index_name = index_name,
                    index = camera_indices_filt,

                    contrast_name = contrast_name,
                    contrast = contrast_mat_trimmed)

  return(info_list)

}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
runGsea <- function(index_name, contrast_name, list_of_de_proteins, list_of_gene_sets, min_set_size = 4) {

  gene_list <- list_of_de_proteins[[contrast_name]]

  msigdb_gene_set <- geneIds(list_of_gene_sets[[index_name]])

  query_gene_list <- data.frame(gene = names(gene_list))

  term_to_gene_tab <- tibble(term = names(msigdb_gene_set), gene = msigdb_gene_set) %>%
    unnest(gene) %>%
    dplyr::inner_join(query_gene_list, by = c("gene"))

  terms_to_keep <- term_to_gene_tab %>%
    group_by(term) %>%
    summarise(counts = n()) %>%
    ungroup() %>%
    dplyr::filter(counts >= min_set_size) %>%
    dplyr::select(-counts)

  term_to_gene_tab_filt <- term_to_gene_tab %>%
    inner_join(terms_to_keep, by = "term") %>%
    mutate(gene = as.character(gene))

  ## Check that there is overlap
  # intersect( names( gene_list_final) ,  unique( term_to_gene_tab_filt$gene )) %>% length


  gsea_results <- GSEA(geneList = gene_list, TERM2GENE = as.data.frame(term_to_gene_tab_filt))

  return(gsea_results)

}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
runEnricher <- function(index_name, contrast_name, list_of_de_proteins, list_of_gene_sets, min_set_size = 4) {

  gene_list <- list_of_de_proteins[[contrast_name]]

  msigdb_gene_set <- geneIds(list_of_gene_sets[[index_name]])

  query_gene_list <- data.frame(gene = gene_list)

  term_to_gene_tab <- tibble(term = names(msigdb_gene_set), gene = msigdb_gene_set) %>%
    unnest(gene) %>%
    dplyr::inner_join(query_gene_list, by = c("gene"))

  terms_to_keep <- term_to_gene_tab %>%
    group_by(term) %>%
    summarise(counts = n()) %>%
    ungroup() %>%
    dplyr::filter(counts >= min_set_size) %>%
    dplyr::select(-counts)

  term_to_gene_tab_filt <- term_to_gene_tab %>%
    inner_join(terms_to_keep, by = "term") %>%
    mutate(gene = as.character(gene))

  ## Check that there is overlap
  # intersect( names( gene_list_final) ,  unique( term_to_gene_tab_filt$gene )) %>% length

  print(intersect(gene_list, unique(term_to_gene_tab_filt$gene)) %>% length)


  gsea_results <- enricher(gene = gene_list, TERM2GENE = as.data.frame(term_to_gene_tab_filt))

  return(gsea_results)

}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Convert uniprot_acc to entrez_id
## gene with two proteins, choose the one with best q-value
#'@export
mergeWithEntrezId <- function(input_table, lookup_table) {

  de_prot_for_camera_helper <- input_table %>%
    mutate(protein_id = row_number()) %>%
    inner_join(lookup_table %>% dplyr::filter(!is.na(ENTREZ_GENE)), by = c("uniprot_acc" = "UNIPROTKB"))

  without_entrez_id <- input_table %>%
    left_join(lookup_table %>% dplyr::filter(!is.na(ENTREZ_GENE)), by = c("uniprot_acc" = "UNIPROTKB")) %>%
    dplyr::filter(is.na(ENTREZ_GENE))

  ## Find rows with duplicated ids
  duplicated_row_id <- de_prot_for_camera_helper %>%
    dplyr::group_by(ENTREZ_GENE, comparison) %>%
    summarise(counts = n()) %>%
    ungroup() %>%
    dplyr::filter(counts > 1) %>%
    dplyr::select(-counts)

  to_be_selected <- de_prot_for_camera_helper %>%
    inner_join(duplicated_row_id, by = c("ENTREZ_GENE", "comparison"))

  ## Duplicates with best q-value
  selected_one_helper <- to_be_selected %>%
    inner_join(to_be_selected %>%
                 group_by(comparison, ENTREZ_GENE) %>%
                 summarise(q.mod = min(q.mod)) %>%
                 ungroup %>%
                 dplyr::select( comparison, ENTREZ_GENE, q.mod),
               by = c("comparison" = "comparison",
                      "ENTREZ_GENE" = "ENTREZ_GENE",
                      "q.mod" = "q.mod"))

  selected_one  <- selected_one_helper  %>%
    inner_join(selected_one_helper %>%
                 group_by(comparison, ENTREZ_GENE) %>%
                 mutate( gene_rank  = row_number()) %>%
                 ungroup %>%
                 dplyr::filter( gene_rank == 1) %>%
                 dplyr::select(-gene_rank) %>%
                 dplyr::select( comparison, ENTREZ_GENE, uniprot_acc, protein_id),
               by = c("comparison" = "comparison",
                      "uniprot_acc" = "uniprot_acc",
                      "ENTREZ_GENE" = "ENTREZ_GENE",
                      "protein_id" = "protein_id")  )

  ## List of rows that were discarded as they are duplicates
  excluded_duplicates <- to_be_selected %>%
    anti_join(selected_one, by = c("protein_id" = "protein_id",
                                   "comparison" = "comparison"))

  ## Combine those that are not duplicates
  de_prot_for_camera_no_dup <- de_prot_for_camera_helper %>%
    anti_join(duplicated_row_id, by = c("ENTREZ_GENE", "comparison"))

  de_prot_for_camera_cleaned <- selected_one %>%
    bind_rows(de_prot_for_camera_no_dup) %>%
    dplyr::arrange(comparison, q.mod)

  de_prot_for_camera_tab <- de_prot_for_camera_cleaned %>%
    dplyr::select(-protein_id) %>%
    relocate(ENTREZ_GENE, .before = comparison)

  return(list(de_proteins = de_prot_for_camera_tab,
              without_entrez_id = without_entrez_id,
              excluded_duplicates = excluded_duplicates))
}

