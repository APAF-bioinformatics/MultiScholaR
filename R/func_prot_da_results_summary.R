# countStatDaGenes
# ----------------------------------------------------------------------------
#' Count the number of statistically significant differentially abundant proteins (according to user-defined threshold)
#' @param lfc_thresh A numerical value specifying the log fold-change threhold (absolute value) for calling statistically significant proteins.
#' @param q_val_thresh A numerical value specifying the q-value threshold for statistically significant proteins.
#' @param log_fc_column The name of the log fold-change column (tidyverse style).
#' @param q_value_column The name of the q-value column (tidyverse style).
#' @return A table with the following columns:
#' status: The status could be Significant and Up, Significant and Down or Not significant
#' counts: The number of proteins wit this status
#' @export
countStatDaGenes <- function(data,
                             lfc_thresh = 0,
                             q_val_thresh = 0.05,
                             log_fc_column = log2FC,
                             q_value_column = fdr_qvalue) {
  # comparison <- as.data.frame(data) |>
  #   distinct(comparison) |>
  #   pull(comparison)

  selected_data <- data |>
    dplyr::mutate(status = case_when(
      {{ q_value_column }} >= q_val_thresh ~ "Not significant",
      {{ log_fc_column }} >= lfc_thresh & {{ q_value_column }} < q_val_thresh ~ "Significant and Up",
      {{ log_fc_column }} < lfc_thresh & {{ q_value_column }} < q_val_thresh ~ "Significant and Down",
      TRUE ~ "Not significant"
    ))

  counts <- selected_data |>
    group_by(status) |>
    summarise(counts = n()) |>
    ungroup()

  all_possible_status <- data.frame(status = c("Not significant", "Significant and Up", "Significant and Down"))

  results <- all_possible_status |>
    left_join(counts, by = c("status" = "status")) |>
    mutate(counts = ifelse(is.na(counts), 0, counts))

  return(results)
}

# ----------------------------------------------------------------------------
# countStatDaGenesHelper
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @export
countStatDaGenesHelper <- function(
  da_table,
  description,
  facet_column = analysis_type,
  comparison_column = "comparison",
  expression_column = "expression",
  lfc_thresh = 0,
  q_val_thresh = 0.05,
  log_fc_column = logFC,
  q_value_column = fdr_qvalue
) {
  message("--- Entering countStatDaGenesHelper (DEBUG66) ---")
  message(paste("   countStatDaGenesHelper: da_table class =", class(da_table)))
  message(paste("   countStatDaGenesHelper: da_table is list =", is.list(da_table)))
  message(paste("   countStatDaGenesHelper: da_table length =", length(da_table)))
  message("   countStatDaGenesHelper: da_table names:")
  print(names(da_table))

  # CRITICAL FIX 1: If da_table is a list of tables from limma, we need to process each one
  da_table_updated <- if (is.list(da_table) && !is.data.frame(da_table)) {
    purrr::imap(da_table, \(x, n) {
      if (is.data.frame(x)) {
        # Ensure comparison column is set if missing
        if (!"comparison" %in% colnames(x)) {
          x <- x |> dplyr::mutate(comparison = n)
        }
      }
      countStatDaGenes(x,
                       lfc_thresh = lfc_thresh,
                       q_val_thresh = q_val_thresh,
                       log_fc_column = {{ log_fc_column }},
                       q_value_column = {{ q_value_column }})
    })
  } else {
    list(countStatDaGenes(da_table,
                          lfc_thresh = lfc_thresh,
                       q_val_thresh = q_val_thresh,
                       log_fc_column = {{ log_fc_column }},
                       q_value_column = {{ q_value_column }}))
  }

  message("   countStatDaGenesHelper: da_table_updated created")
  message(paste("   countStatDaGenesHelper: da_table_updated length =", length(da_table_updated)))
  message("   countStatDaGenesHelper: da_table_updated names:")
  print(names(da_table_updated))

  list_of_tables <- purrr::map2(
    da_table_updated,
    names(da_table_updated),
    \(.x, .y){
      message(paste("      [map2] Processing element with name:", .y))
      .x |> mutate(!!sym(comparison_column) := .y)
    }
  )

  message("   countStatDaGenesHelper: list_of_tables created")
  message("   countStatDaGenesHelper: About to bind_rows...")

  bound_tables <- list_of_tables |> bind_rows()
  message(paste("   countStatDaGenesHelper: bound_tables dims =", nrow(bound_tables), "x", ncol(bound_tables)))

  faceted_tables <- bound_tables |> mutate({{ facet_column }} := description)
  message("   countStatDaGenesHelper: facet column added")
  message("   countStatDaGenesHelper: Checking comparison column values:")
  print(unique(faceted_tables[[comparison_column]]))

  message(paste("   countStatDaGenesHelper: About to separate_wider_delim on column:", comparison_column))
  message(paste("   countStatDaGenesHelper: Looking for delimiter: ="))

  # Check if any values contain "="
  has_delimiter <- any(grepl("=", faceted_tables[[comparison_column]]))
  message(paste("   countStatDaGenesHelper: Any values contain '=' ?", has_delimiter))

  if (has_delimiter) {
    merged_tables <- faceted_tables |>
      separate_wider_delim(!!sym(comparison_column),
        delim = "=",
        names = c(
          comparison_column,
          expression_column
        ),
        too_few = "align_start"
      )
    message("   countStatDaGenesHelper: Separation complete")
  } else {
    message("   countStatDaGenesHelper: No '=' found, skipping separation and using original comparison names")
    # Ensure expression column exists even if no separation
    merged_tables <- faceted_tables
    if (!expression_column %in% colnames(merged_tables)) {
      merged_tables[[expression_column]] <- merged_tables[[comparison_column]]
    }
  }

  message("--- Exiting countStatDaGenesHelper (DEBUG66) ---")
  merged_tables
}

# ----------------------------------------------------------------------------
# printCountDaGenesTable
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Format results table for use in volcano plots, counting number of significant proteins, p-values distribution histogram.
#' @param list_of_da_tables A list with each element being a results table with log fold-change and q-value per protein.
#' @param list_of_descriptions  A list of strings describing the parameters used to generate the result table.
#' @param formula_string The formula string used in the facet_grid command for the ggplot scatter plot.
#' @param facet_column The name of the column describing the type of analysis or parameters used to generate the result table (tidyverse style). This is related to the \code{list_of_descriptions} parameter above.
#' @param comparison_column The name of the column describing the contrasts or comparison between groups (tidyverse style).
#' @param expression_column The name of the column that will contain the formula expressions of the contrasts.
#' @export
printCountDaGenesTable <- function(
  list_of_da_tables,
  list_of_descriptions,
  formula_string = "analysis_type ~ comparison",
  facet_column = analysis_type,
  comparison_column = "comparison",
  expression_column = "expression"
) {
  num_significant_da_genes_all <- purrr::map2(
    list_of_da_tables,
    list_of_descriptions,
    function(a, b) {
      countStatDaGenesHelper(
        da_table = a,
        description = b,
        facet_column = {{ facet_column }},
        comparison_column = comparison_column,
        expression_column = expression_column
      )
    }
  ) |>
    bind_rows()

  num_sig_da_genes_barplot <- num_significant_da_genes_all |>
    dplyr::filter(status != "Not significant") |>
    ggplot(aes(x = status, y = counts)) +
    geom_bar(stat = "identity") +
    geom_text(stat = "identity", aes(label = counts), vjust = -0.5) +
    theme(axis.text.x = element_text(angle = 90))

  # print(head(num_sig_da_genes_barplot))

  if (!is.na(formula_string)) {
    num_sig_da_genes_barplot <- num_sig_da_genes_barplot +
      facet_grid(as.formula(formula_string))
  }


  return(list(plot = num_sig_da_genes_barplot, table = num_significant_da_genes_all))
}

# ----------------------------------------------------------------------------
# getSignificantData
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Format results table for use in volcano plots, counting number of significant proteins, p-values distribution histogram.
#' @param list_of_da_tables A list with each element being a results table with log fold-change and q-value per protein.
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
getSignificantData <- function(
  list_of_da_tables,
  list_of_descriptions,
  row_id = uniprot_acc,
  p_value_column = raw_pvalue,
  q_value_column = fdr_qvalue,
  fdr_value_column = fdr_value_bh_adjustment,
  log_q_value_column = lqm,
  log_fc_column = logFC,
  comparison_column = "comparison",
  expression_column = "log_intensity",
  facet_column = analysis_type,
  q_val_thresh = 0.05
) {
  message("--- Entering getSignificantData (DEBUG66) ---")
  message(paste("   getSignificantData: list_of_da_tables class =", class(list_of_da_tables)))
  message(paste("   getSignificantData: list_of_da_tables length =", length(list_of_da_tables)))
  message("   getSignificantData: list_of_da_tables structure:")
  str(list_of_da_tables, max.level = 2)

  get_row_binded_table <- function(de_table_list, description) {
    message("   --- Entering get_row_binded_table (DEBUG66) ---")
    message(paste("      get_row_binded_table: de_table_list class =", class(de_table_list)))
    message(paste("      get_row_binded_table: de_table_list length =", length(de_table_list)))
    message("      get_row_binded_table: de_table_list names:")
    print(names(de_table_list))

    # This internal helper is now more robust. It checks if the table
    # already has the row_id as a column. If not, it converts rownames.
    # This handles both old and new data structures.

    processed_list <- purrr::map(de_table_list, function(tbl) {
      row_id_col <- as_string(as_name(enquo(row_id)))
      message(paste("         [map] Processing table, row_id_col =", row_id_col))
      message(paste("         [map] Table class =", class(tbl)))
      message(paste("         [map] Table is data.frame =", is.data.frame(tbl)))

      if (!row_id_col %in% colnames(tbl)) {
        # If row_id is not a column, convert rownames
        message(paste("         [map] row_id not in columns, converting rownames"))
        tbl <- tbl |> rownames_to_column(var = row_id_col)
      } else {
        message(paste("         [map] row_id already in columns"))
      }

      return(tbl)
    })

    message("      get_row_binded_table: processed_list created")
    message(paste("      get_row_binded_table: processed_list length =", length(processed_list)))

    output <- processed_list |>
      purrr::map2(names(processed_list), \(.x, .y){
        message(paste("         [map2] Processing element with name:", .y))
        message(paste("         [map2] comparison_column already exists:", comparison_column %in% colnames(.x)))
        # If the 'comparison' column doesn't already exist, create it from the list name
        if (!comparison_column %in% colnames(.x)) {
          message(paste("         [map2] Adding comparison column with value:", .y))
          .x <- .x |> mutate({{ comparison_column }} := .y)
        }
        .x
      }) |>
      bind_rows() |>
      mutate({{ facet_column }} := description)

    message("      get_row_binded_table: output table created")
    message(paste("      get_row_binded_table: output dims =", nrow(output), "x", ncol(output)))
    message("      get_row_binded_table: Checking comparison column values:")
    print(unique(output[[comparison_column]]))

    # This separator logic assumes a specific format like 'ContrastName=ExpressionType'
    # in the comparison column. We need to make this conditional as well.
    # Check if any values in the comparison column contain '=' before trying to separate.
    has_delimiter <- any(grepl("=", output[[comparison_column]]))
    message(paste("      get_row_binded_table: Any values contain '=' ?", has_delimiter))

    if (has_delimiter) {
      message("      get_row_binded_table: Separating on '=' delimiter")
      output <- output |>
        separate_wider_delim(
          {{ comparison_column }},
          delim = "=",
          names = c(comparison_column, expression_column)
        )
      message("      get_row_binded_table: Separation complete")
    } else {
      message("      get_row_binded_table: No '=' delimiter found, skipping separation")
    }

    message("   --- Exiting get_row_binded_table (DEBUG66) ---")
    return(output)
  }

  message("   getSignificantData: About to call get_row_binded_table for each list element...")

  logfc_tbl_all <- purrr::map2(
    list_of_da_tables, list_of_descriptions,
    function(a, b) {
      get_row_binded_table(de_table_list = a, description = b)
    }
  ) |>
    bind_rows()

  message("   getSignificantData: All tables bound")
  message(paste("   getSignificantData: logfc_tbl_all dims =", nrow(logfc_tbl_all), "x", ncol(logfc_tbl_all)))

  # CRITICAL FIX: Ensure expression_column exists if it doesn't
  expr_col_name <- expression_column
  if (!expr_col_name %in% colnames(logfc_tbl_all)) {
    message(sprintf("   getSignificantData: expression_column '%s' not found, creating dummy", expr_col_name))
    logfc_tbl_all[[expr_col_name]] <- NA_real_
  }

  selected_data <- logfc_tbl_all |>
    mutate({{ log_q_value_column }} := -log10(fdr_qvalue)) |>
    dplyr::select(
      {{ row_id }}, {{ log_q_value_column }}, {{ q_value_column }}, {{ p_value_column }}, {{ log_fc_column }},
      {{ comparison_column }}, {{ expression_column }},
      {{ facet_column }}
    ) |>
    dplyr::mutate(colour = case_when(
      abs({{ log_fc_column }}) >= 1 & {{ q_value_column }} >= q_val_thresh ~ "orange",
      abs({{ log_fc_column }}) >= 1 & {{ q_value_column }} < q_val_thresh ~ "purple",
      abs({{ log_fc_column }}) < 1 & {{ q_value_column }} < q_val_thresh ~ "blue",
      TRUE ~ "black"
    )) |>
    dplyr::mutate(colour = factor(colour, levels = c("black", "orange", "blue", "purple")))

  message("--- Exiting getSignificantData (DEBUG66) ---")
  selected_data
}
