#' @title Choose Best Protein Accession for ProteinQuantitativeData
#' @name chooseBestProteinAccession,ProteinQuantitativeData-method
#' @export
#' @param theObject The object of class ProteinQuantitativeData
#' @param delim The delimiter used to split the protein accessions
#' @param seqinr_obj The object of class Seqinr::seqinr
#' @param seqinr_accession_column The column in the seqinr object that contains the protein accessions
#' @param replace_zero_with_na Replace zero values with NA
#' @param aggregation_method Method to aggregate protein values: "sum", "mean", or "median" (default: "sum")
setMethod(
  f = "chooseBestProteinAccession",
  signature = "ProteinQuantitativeData",
  definition = function(
    theObject, delim = NULL, seqinr_obj = NULL,
    seqinr_accession_column = NULL,
    replace_zero_with_na = NULL,
    aggregation_method = NULL
  ) {
    cat("\n\n+===========================================================================+\n")
    cat("|           ENTERING chooseBestProteinAccession (DEBUG66)               |\n")
    cat("+===========================================================================+\n\n")

    protein_quant_table <- theObject@protein_quant_table
    protein_id_column <- theObject@protein_id_column

    cat(">>> STEP 1: PARAMETER INPUTS (before checkParamsObjectFunctionSimplify) <<<\n")
    cat(sprintf("   Input parameter 'delim' (user provided): %s\n", ifelse(is.null(delim), "NULL", delim)))
    cat(sprintf("   Input parameter 'seqinr_obj' is NULL: %s\n", is.null(seqinr_obj)))
    cat(sprintf("   Input parameter 'seqinr_accession_column': %s\n", ifelse(is.null(seqinr_accession_column), "NULL", seqinr_accession_column)))
    cat(sprintf("   Input parameter 'replace_zero_with_na': %s\n", ifelse(is.null(replace_zero_with_na), "NULL", replace_zero_with_na)))
    cat(sprintf("   Input parameter 'aggregation_method': %s\n", ifelse(is.null(aggregation_method), "NULL", aggregation_method)))
    cat(sprintf("   Protein ID column from S4: %s\n", protein_id_column))
    cat(sprintf("   Protein table dimensions: %d rows x %d cols\n", nrow(protein_quant_table), ncol(protein_quant_table)))
    cat(sprintf("   First 5 Protein IDs: %s\n", paste(head(protein_quant_table[[protein_id_column]], 5), collapse = ", ")))
    cat("\n")

    delim <- checkParamsObjectFunctionSimplify(theObject, "delim", default_value = " |;|:|\\|")
    seqinr_obj <- checkParamsObjectFunctionSimplify(theObject, "seqinr_obj", default_value = NULL)
    seqinr_accession_column <- checkParamsObjectFunctionSimplify(theObject,
      "seqinr_accession_column",
      default_value = NULL
    )
    replace_zero_with_na <- checkParamsObjectFunctionSimplify(theObject,
      "replace_zero_with_na",
      default_value = FALSE
    )
    aggregation_method <- checkParamsObjectFunctionSimplify(theObject,
      "aggregation_method",
      default_value = "sum"
    )

    cat(">>> STEP 2: PARAMETERS AFTER checkParamsObjectFunctionSimplify <<<\n")
    cat(sprintf("   Resolved 'delim': '%s'\n", delim))
    cat(sprintf("   Resolved 'seqinr_obj' is NULL: %s\n", is.null(seqinr_obj)))
    cat(sprintf("   Resolved 'seqinr_accession_column': %s\n", ifelse(is.null(seqinr_accession_column), "NULL", seqinr_accession_column)))
    cat(sprintf("   Resolved 'replace_zero_with_na': %s\n", replace_zero_with_na))
    cat(sprintf("   Resolved 'aggregation_method': %s\n", aggregation_method))
    cat("\n")

    if (!aggregation_method %in% c("sum", "mean", "median")) {
      stop("aggregation_method must be one of: 'sum', 'mean', 'median'")
    }

    # --- Robustly handle list input for seqinr_obj ---
    if (is.list(seqinr_obj) && !is.data.frame(seqinr_obj) && "aa_seq_tbl_final" %in% names(seqinr_obj)) {
      cat("   [DEBUG66] Extracting data frame from seqinr_obj list...\n")
      seqinr_obj <- seqinr_obj$aa_seq_tbl_final
    }


    theObject <- updateParamInObject(theObject, "delim")
    theObject <- updateParamInObject(theObject, "seqinr_obj")
    theObject <- updateParamInObject(theObject, "seqinr_accession_column")
    theObject <- updateParamInObject(theObject, "replace_zero_with_na")
    theObject <- updateParamInObject(theObject, "aggregation_method")

    evidence_tbl_cleaned <- protein_quant_table |>
      distinct() |>
      mutate(row_id = row_number() - 1)

    cat(">>> STEP 3: CALLING chooseBestProteinAccessionHelper <<<\n")
    cat(sprintf("   Passing 'delim' to helper: '%s'\n", delim))
    cat(sprintf("   First 5 Protein IDs to process: %s\n", paste(head(evidence_tbl_cleaned[[protein_id_column]], 5), collapse = ", ")))
    cat(sprintf("   Number of rows to process: %d\n", nrow(evidence_tbl_cleaned)))
    cat("\n")

    accession_gene_name_tbl <- chooseBestProteinAccessionHelper(
      input_tbl = evidence_tbl_cleaned,
      acc_detail_tab = seqinr_obj,
      accessions_column = !!sym(protein_id_column),
      row_id_column = seqinr_accession_column,
      group_id = row_id,
      delim = delim
    )

    cat(">>> STEP 4: RETURNED FROM chooseBestProteinAccessionHelper <<<\n")
    cat(sprintf("   Number of rows returned: %d\n", nrow(accession_gene_name_tbl)))
    if (nrow(accession_gene_name_tbl) > 0) {
      cat(sprintf("   First 5 accessions: %s\n", paste(head(accession_gene_name_tbl[[seqinr_accession_column]], 5), collapse = ", ")))
    }
    cat("\n")

    protein_log2_quant_cln <- evidence_tbl_cleaned |>
      left_join(
        accession_gene_name_tbl |>
          dplyr::distinct(row_id, !!sym(as.character(seqinr_accession_column))),
        by = join_by(row_id)
      ) |>
      mutate(!!sym(theObject@protein_id_column) := !!sym(as.character(seqinr_accession_column))) |>
      dplyr::select(-row_id, -!!sym(as.character(seqinr_accession_column)))

    protein_id_table <- evidence_tbl_cleaned |>
      left_join(
        accession_gene_name_tbl |>
          dplyr::distinct(row_id, !!sym(as.character(seqinr_accession_column))),
        by = join_by(row_id)
      ) |>
      distinct(uniprot_acc, !!sym(protein_id_column)) |>
      mutate(!!sym(paste0(protein_id_column, "_list")) := !!sym(protein_id_column)) |>
      mutate(!!sym(protein_id_column) := !!sym("uniprot_acc")) |>
      distinct(!!sym(protein_id_column), !!sym(paste0(protein_id_column, "_list"))) |>
      group_by(!!sym(protein_id_column)) |>
      summarise(!!sym(paste0(protein_id_column, "_list")) := paste(!!sym(paste0(protein_id_column, "_list")), collapse = ";")) |>
      ungroup() |>
      mutate(!!sym(paste0(protein_id_column, "_list")) := purrr::map_chr(
        !!sym(paste0(protein_id_column, "_list")),
        \(x){
          paste(unique(sort(str_split(x, ";")[[1]])), collapse = ";")
        }
      ))

    summed_data <- protein_log2_quant_cln |>
      mutate(!!sym(protein_id_column) := purrr::map_chr(!!sym(protein_id_column), \(x){
        str_split(x, delim)[[1]][1]
      })) |>
      pivot_longer(
        cols = !matches(protein_id_column),
        names_to = "sample_id",
        values_to = "temporary_values_choose_accession"
      ) |>
      group_by(!!sym(protein_id_column), sample_id) |>
      summarise(
        is_na = sum(is.na(temporary_values_choose_accession)),
        temporary_values_choose_accession = case_when(
          all(is.na(temporary_values_choose_accession)) ~ NA_real_,
          aggregation_method == "sum" ~ sum(temporary_values_choose_accession, na.rm = TRUE),
          aggregation_method == "mean" ~ mean(temporary_values_choose_accession, na.rm = TRUE),
          aggregation_method == "median" ~ median(temporary_values_choose_accession, na.rm = TRUE)
        ),
        num_values = n()
      ) |>
      ungroup() |>
      pivot_wider(
        id_cols = !!sym(protein_id_column),
        names_from = sample_id,
        values_from = temporary_values_choose_accession,
        values_fill = NA_real_
      )

    if (replace_zero_with_na == TRUE) {
      summed_data[is.na(summed_data)] <- NA
    }

    cat(">>> STEP 5: CALLING rankProteinAccessionHelper <<<\n")
    cat(sprintf("   Passing 'delim' to helper: '%s'\n", delim))
    cat(sprintf("   Number of protein_id_table rows to rank: %d\n", nrow(protein_id_table)))
    if (nrow(protein_id_table) > 0) {
      list_col_name <- paste0(protein_id_column, "_list")
      cat(sprintf("   First 5 protein ID lists: %s\n", paste(head(protein_id_table[[list_col_name]], 5), collapse = ", ")))
    }
    cat("\n")

    protein_id_table <- rankProteinAccessionHelper(
      input_tbl = protein_id_table,
      acc_detail_tab = seqinr_obj,
      accessions_column = !!sym(paste0(protein_id_column, "_list")),
      row_id_column = seqinr_accession_column,
      group_id = !!sym(protein_id_column),
      delim = delim
    ) |>
      dplyr::rename(!!sym(paste0(protein_id_column, "_list")) := seqinr_accession_column) |>
      dplyr::select(-num_gene_names, -gene_names, -is_unique)

    cat(">>> STEP 6: RETURNED FROM rankProteinAccessionHelper <<<\n")
    cat(sprintf("   Number of rows returned: %d\n", nrow(protein_id_table)))
    cat("\n")

    theObject@protein_id_table <- protein_id_table
    theObject@protein_quant_table <- summed_data[, colnames(protein_quant_table)]

    cat(">>> STEP 7: FINAL OUTPUT <<<\n")
    cat(sprintf("   Final protein_quant_table dimensions: %d rows x %d cols\n", nrow(theObject@protein_quant_table), ncol(theObject@protein_quant_table)))
    cat(sprintf("   Final protein_id_table dimensions: %d rows x %d cols\n", nrow(theObject@protein_id_table), ncol(theObject@protein_id_table)))
    cat(sprintf("   First 5 final Protein IDs: %s\n", paste(head(theObject@protein_quant_table[[protein_id_column]], 5), collapse = ", ")))
    cat("\n")
    cat("+===========================================================================+\n")
    cat("|           EXITING chooseBestProteinAccession (DEBUG66)                |\n")
    cat("+===========================================================================+\n\n")

    return(theObject)
  }
)

#' @title Choose Best Protein Accession Sum Duplicates
#' @name chooseBestProteinAccessionSumDuplicates,ProteinQuantitativeData-method
#' @export
setMethod(
  f = "chooseBestProteinAccessionSumDuplicates",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, delim = ";", quant_columns_pattern = "\\d+", islogged = TRUE) {
    protein_quant_table <- theObject@protein_quant_table
    protein_id_column <- theObject@protein_id_column

    protein_log2_quant_cln <- protein_quant_table |>
      mutate(!!sym(protein_id_column) := str_split_i(!!sym(protein_id_column), delim, 1)) |>
      group_by(!!sym(protein_id_column)) |>
      summarise(across(
        matches(quant_columns_pattern),
        \(x){
          if (islogged == TRUE) {
            log2(sum(2^x, na.rm = TRUE))
          } else {
            sum(x, na.rm = TRUE)
          }
        }
      )) |>
      ungroup()

    theObject@protein_quant_table <- protein_log2_quant_cln

    theObject
  }
)
