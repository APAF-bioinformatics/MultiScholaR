# enrichProteinsPathwaysHelper
# ----------------------------------------------------------------------------
#' Helper function to perform protein pathway enrichment analysis for a single contrast
#'
#' @description
#' This function performs pathway enrichment analysis on protein data using GO terms and gene symbols
#' downloaded from UniProt. The data is cached for future use to improve performance.
#'
#' @param da_analysis_results Output from deAnalysisWrapperFunction containing differential expression results
#' @param organism_taxid NCBI taxonomy ID for the organism (e.g., "9606" for human)
#' @param min_gene_set_size Minimum number of genes in a gene set (default: 4)
#' @param max_gene_set_size Maximum number of genes in a gene set (default: 200)
#' @param p_val_thresh P-value threshold for enrichment significance (default: 0.05)
#' @param protein_p_val_thresh P-value threshold for protein significance (default: 0.05)
#' @param cache_dir Directory to store cached UniProt data (default: "cache")
#' @param output_dir Directory for output files (default: "proteins_pathways_enricher")
#' @param use_cached Whether to use cached data if available (default: TRUE)
#' @param protein_id_delimiter Delimiter used in protein IDs (default: ":")
#' @param protein_id_column Name of the protein ID column (default: "Protein.Ids")
#'
#' @return A data frame containing enrichment results
#'
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @rawNamespace import(plotly, except = last_plot)
#' @importFrom purrr map map_chr walk
#' @importFrom stringr str_split str_replace_all
#'
#' @export
#' @param cache_file Runtime inputs used by this function; see the usage section for accepted values.
enrichProteinsPathwaysHelper <- function(da_analysis_results,
                                  organism_taxid,
                                  min_gene_set_size = 4,
                                  max_gene_set_size = 200,
                                  p_val_thresh = 0.05,
                                  protein_p_val_thresh = 0.05,
                                  cache_dir = "cache",
                                  cache_file = "uniprot_annotations.RDS",
                                  output_dir = "proteins_pathways_enricher",
                                  use_cached = TRUE,
                                  protein_id_delimiter = ":",
                                  protein_id_column = "Protein.Ids"
                                  ) {

  # Create directories if they don't exist
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Initialize UniProt.ws handle
  up <- UniProt.ws::UniProt.ws(taxId = organism_taxid)

  # Cache file paths
  go_cache_file <- file.path(cache_dir, cache_file)

  # Extract protein IDs from all rows and clean them
  # Handle both da_proteins_wide and de_proteins_wide for backward compatibility
  wide_table <- da_analysis_results$da_proteins_wide
  if (is.null(wide_table)) {
    wide_table <- da_analysis_results$de_proteins_wide
  }
  
  if (is.null(wide_table)) {
    stop("da_analysis_results must contain either da_proteins_wide or de_proteins_wide")
  }

  protein_data <- wide_table |>
    dplyr::select(!!sym(protein_id_column)) |>
    distinct() |>
    dplyr::mutate(!!sym(protein_id_column) := purrr::map_chr(!!sym(protein_id_column), \(x) {
      str_split(x, protein_id_delimiter)[[1]][1]
    })) |>
    dplyr::mutate(uniprot_acc = .cleanProteinIds(!!sym(protein_id_column)))

  # Get cached or fresh data
  uniprot_data <-  download_uniprot_data(protein_ids = protein_data
                                , cache_file = go_cache_file
                                , uniprot_handle = up
                                , protein_id_delimiter = protein_id_delimiter)
  # Process protein data for each comparison
  enrichment_results <- list()

    # Get significant proteins for this comparison and clean their IDs
    sig_data <- da_analysis_results$significant_rows |>
      dplyr::filter(analysis_type == "RUV applied") |>
      dplyr::mutate(!!sym(protein_id_column) := purrr::map_chr(!!sym(protein_id_column), \(x) {
        str_split(x, protein_id_delimiter)[[1]][1]
      })) |>
      dplyr::mutate(uniprot_acc = .cleanProteinIds(!!sym(protein_id_column)))

    # Get positive and negative proteins
    pos_prots <- sig_data |>
      dplyr::filter(fdr_qvalue < protein_p_val_thresh,
                   log2FC > 0) |>
      dplyr::pull(uniprot_acc)

    neg_prots <- sig_data |>
      dplyr::filter(fdr_qvalue < protein_p_val_thresh,
                   log2FC < 0) |>
      dplyr::pull(uniprot_acc)

    # Background proteins - use all proteins from the dataset
    background_proteins <- protein_data |>
      dplyr::pull(uniprot_acc)


    # Create GO annotation data frame from UniProt data
    uniprot_data_updated <- uniprot_data
    if( "Entry" %in% colnames(uniprot_data) & !("uniprot_acc" %in%  colnames(uniprot_data)) ) {
      uniprot_data_updated <- uniprot_data|>
        dplyr::rename( uniprot_acc = "Entry")
    }

   goterms <- AnnotationDbi::Term(GO.db::GOTERM)
   gotypes <- AnnotationDbi::Ontology(GO.db::GOTERM)

    go_annot_temp <- uniprot_data_updated |>
      pivot_longer( cols=matches("go_id"), names_to = "go_type", values_to = "go_id") |>
      dplyr::select(-go_type) |>
      separate_rows(go_id, sep = "; ")

    go_id_to_term <- go_annot_temp |>
      distinct(go_id)    |>
      dplyr::mutate(
        go_term = purrr::map_chr(go_id, \(x) {
          if (x %in% names(goterms)) {
            return(goterms[[x]])
          }
          return(NA_character_)
        }),
        go_type = purrr::map_chr(go_id, \(x) {
          if (x %in% names(gotypes)) {
            type <- gotypes[[x]]
            return(dplyr::case_when(
              type == "BP" ~ "Biological Process",
              type == "CC" ~ "Cellular Compartment",
              type == "MF" ~ "Molecular Function",
              TRUE ~ NA_character_
            ))
          }
          return(NA_character_)
        }) )

    go_annot <- go_annot_temp |>
      left_join( go_id_to_term, by = "go_id") |>
      dplyr::distinct(uniprot_acc, go_id, go_term, go_type)  |>
      dplyr::filter(!is.na(go_id))

  # Create empty result tibble
  empty_result <- tibble::tibble(
    ID = character(0),
    Description = character(0),
    GeneRatio = character(0),
    BgRatio = character(0),
    pvalue = numeric(0),
    p.adjust = numeric(0),
    qvalue = numeric(0),
    geneID = character(0),
    Count = integer(0),
    term = character(0),
    go_type = character(0)
  )

  # Perform enrichment analysis
  pos_enrich <- runOneGoEnrichmentInOutFunction(
    significant_proteins = pos_prots,
    background_proteins = background_proteins,
    go_annotations = go_annot,
    uniprot_data = uniprot_data,
    p_val_thresh = p_val_thresh,
    min_gene_set_size = min_gene_set_size,
    max_gene_set_size = max_gene_set_size
  )

  neg_enrich <- runOneGoEnrichmentInOutFunction(
    significant_proteins = neg_prots,
    background_proteins = background_proteins,
    go_annotations = go_annot,
    uniprot_data = uniprot_data,
    p_val_thresh = p_val_thresh,
    min_gene_set_size = min_gene_set_size,
    max_gene_set_size = max_gene_set_size
  )

  # Handle NULL results
  pos_enrich <- if(is.null(pos_enrich)) empty_result else pos_enrich
  neg_enrich <- if(is.null(neg_enrich)) empty_result else neg_enrich

  # Add directionality
  enrichment_results <- bind_rows(
    pos_enrich |> mutate(directionality = "positive"),
    neg_enrich |> mutate(directionality = "negative")
  )

  return(enrichment_results)
}


# ----------------------------------------------------------------------------
# enrichProteinsPathways
# ----------------------------------------------------------------------------
#' Perform protein pathway enrichment analysis across multiple contrasts
#'
#' @description
#' This function performs pathway enrichment analysis across multiple contrasts in a proteomics dataset.
#' It processes each contrast separately and combines the results into a single data frame.
#'
#' @param da_analysis_results_list List of differential expression results for each contrast
#' @param taxon_id NCBI taxonomy ID for the organism (e.g., "9606" for human)
#' @param protein_id_delimiter Delimiter used in protein IDs (default: ":")
#' @param protein_p_val_thresh P-value threshold for protein significance (default: 0.05)
#' @param min_gene_set_size Minimum number of genes in a gene set (default: 4)
#' @param max_gene_set_size Maximum number of genes in a gene set (default: 200)
#' @param p_val_thresh P-value threshold for enrichment significance (default: 0.05)
#' @param cache_dir Directory to store cached UniProt data (default: "cache")
#' @param cache_file Name of the cache file for UniProt annotations (default: "uniprot_annotations.RDS")
#' @param use_cached Whether to use cached data if available (default: TRUE)
#'
#' @return A data frame containing combined enrichment results across all contrasts
#'
#' @import dplyr
#' @importFrom purrr map set_names
#'
#' @export
#' @param de_analysis_results_list Runtime inputs used by this function; see the usage section for accepted values.
enrichProteinsPathways <- function(da_analysis_results_list = NULL,
                                 taxon_id,
                                 protein_id_delimiter = ":",
                                 protein_p_val_thresh = 0.05,
                                 min_gene_set_size = 4,
                                 max_gene_set_size = 200,
                                 p_val_thresh = 0.05,
                                 cache_dir = "cache",
                                 cache_file = "uniprot_annotations.RDS",
                                 use_cached = TRUE,
                                 de_analysis_results_list = NULL) {
  
  # Handle alias for backward compatibility
  if (is.null(da_analysis_results_list)) {
    if (!is.null(de_analysis_results_list)) {
      da_analysis_results_list <- de_analysis_results_list
    } else {
      stop("Either da_analysis_results_list or de_analysis_results_list must be provided")
    }
  }

  # Create a list to store all enrichment results
  all_enrichment_results_by_group <- names(da_analysis_results_list) |>
    purrr::set_names() |>  # Keep the contrast names
    purrr::map(\(contrast_name) {
      message(paste("Processing enrichment for contrast:", contrast_name))

      # Get the DE results for this contrast
      da_results <- da_analysis_results_list[[contrast_name]]

      # Run enrichment analysis
      enrichment_result <- enrichProteinsPathwaysHelper(
        da_analysis_results = da_results,
        organism_taxid = as.character(taxon_id),
        protein_p_val_thresh = protein_p_val_thresh,
        min_gene_set_size = min_gene_set_size,
        max_gene_set_size = max_gene_set_size,
        p_val_thresh = p_val_thresh,
        cache_dir = cache_dir,
        cache_file = cache_file,
        use_cached = use_cached,
        protein_id_delimiter = protein_id_delimiter )

      return(enrichment_result)
    })

  # Combine results from all contrasts
  go_results_table_by_group <- bind_rows(all_enrichment_results_by_group, .id = "comparison")

  return(go_results_table_by_group)
}


# ----------------------------------------------------------------------------
# download_uniprot_data
# ----------------------------------------------------------------------------
#'@export
download_uniprot_data <- function(protein_ids, cache_file, uniprot_handle, protein_id_delimiter = ":") {
  # Ensure cache directory exists
  cache_dir <- dirname(cache_file)
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

  # Check if cache exists
  if (file.exists(cache_file)) {
    tryCatch({
      cached_data <- readRDS(cache_file)

      # Find which proteins are not in cache
      missing_proteins <- protein_ids |>
        dplyr::filter(!uniprot_acc %in% cached_data$Entry)

      if (nrow(missing_proteins) == 0) {
        return(cached_data)
      }

      # Download only missing proteins
      message("Downloading data for ", nrow(missing_proteins), " new protein IDs")
      new_data <- batchQueryEvidenceGeneId(
        data.frame(uniprot_acc = unique(missing_proteins$uniprot_acc)  ),
        gene_id_column = "uniprot_acc",
        uniprot_handle = uniprot_handle,
        uniprot_columns = c(
          "protein_existence",
          "annotation_score",
          "reviewed",
          "gene_names",
          "protein_name",
          "length",
          "xref_ensembl",
          "go_id",  # This is the correct column name from UniProt.ws
          "keyword"
        )
      )



      # Process new data
      if (!is.null(new_data) && nrow(new_data) > 0) {

        new_data_processed <- new_data |>
          uniprotGoIdToTerm(
            uniprot_id_column = Entry,
            go_id_column = Gene.Ontology.IDs,
            sep = "; "
          ) |>
          dplyr::rename(
            Protein_existence = "Protein.existence",
            Protein_names = "Protein.names"
          ) |>
          mutate( Protein_existence = purrr::map_dbl(Protein_existence, as.double ))


        # Combine with cached data
        combined_data <- dplyr::bind_rows(cached_data, new_data_processed)

        # Try to update cache
        tryCatch({
          saveRDS(combined_data, cache_file)
        }, error = function(e) {
          warning("Could not update cache file: ", e$message)
        })

        return(combined_data)
      }
      return(cached_data)
    }, error = function(e) {
      warning("Error reading cache: ", e$message, ". Downloading all data fresh.")
      # Fall through to download all data
    })
  }

  # If no cache exists or cache read failed, download all data
  message("Downloading data for all protein IDs")

  all_data <- batchQueryEvidenceGeneId(
    data.frame(uniprot_acc = unique(protein_ids$uniprot_acc) ),
    gene_id_column = "uniprot_acc",
    uniprot_handle = uniprot_handle,
    uniprot_columns = c(
      "protein_existence",
      "annotation_score",
      "reviewed",
      "gene_names",
      "protein_name",
      "length",
      "xref_ensembl",
      "go_id",  # This is the correct column name from UniProt.ws
      "keyword"
    )
  )

  if (!is.null(all_data) && nrow(all_data) > 0) {
    processed_data <- all_data |>
      uniprotGoIdToTerm(
        uniprot_id_column = Entry,
        go_id_column = Gene.Ontology.IDs,
        sep = "; "
      ) |>
      dplyr::rename(
        Protein_existence = "Protein.existence",
        Protein_names = "Protein.names"
      )

    # Try to save to cache
    tryCatch({
      saveRDS(processed_data, cache_file)
    }, error = function(e) {
      warning("Could not save to cache file: ", e$message)
    })

    return(processed_data)
  }

  return(NULL)
}


# ----------------------------------------------------------------------------
# uniprotGoIdToTermSimple
# ----------------------------------------------------------------------------
#' Convert UniProt GO IDs to terms without grouping or pivoting
#' @param uniprot_dat  a table with uniprot accessions and a column with GO-ID
#' @param uniprot_id_column The name of the column with the uniprot accession, as a tidyverse header format, not a string
#' @param go_id_column The name of the column with the GO-ID, as a tidyverse header format, not a string
#' @param goterms Output from running \code{goterms <- Term(GOTERM)} from the GO.db library.
#' @param gotypes Output from running \code{gotypes <- Ontology(GOTERM)} from the GO.db library.
#' @return A table with columns for uniprot_id, go_id, go_term, and go_type
#' @export
#' @param gene_name_column,sep Runtime inputs used by this function; see the usage section for accepted values.
uniprotGoIdToTermSimple <- function(uniprot_dat
                                   , uniprot_id_column = UNIPROTKB
                                   , go_id_column = `GO-IDs`
                                   , gene_name_column = Gene.Names
                                   , sep = "; "
                                   , goterms = AnnotationDbi::Term(GO.db::GOTERM)
                                   , gotypes = AnnotationDbi::Ontology(GO.db::GOTERM)) {

  print("Run uniprotGoIdToTermSimple")

  uniprot_acc_to_go_term <- uniprot_dat |>
    dplyr::distinct({{uniprot_id_column}}, {{gene_name_column}}, {{go_id_column}}) |>
    tidyr::separate_rows({{go_id_column}}, sep = sep) |>
    dplyr::filter(!is.na({{go_id_column}})) |>
    dplyr::mutate(
      go_term = purrr::map_chr({{go_id_column}}, \(x) {
        if (x %in% names(goterms)) {
          return(goterms[[x]])
        }
        return(NA_character_)
      }),
      go_type = purrr::map_chr({{go_id_column}}, \(x) {
        if (x %in% names(gotypes)) {
          type <- gotypes[[x]]
          return(dplyr::case_when(
            type == "BP" ~ "Biological Process",
            type == "CC" ~ "Cellular Compartment",
            type == "MF" ~ "Molecular Function",
            TRUE ~ NA_character_
          ))
        }
        return(NA_character_)
      })
    ) |>
    dplyr::filter(!is.na(go_term)) |>
    dplyr::rename(
      uniprot_acc = {{uniprot_id_column}},
      go_id = {{go_id_column}}
    )

  return(uniprot_acc_to_go_term)
}



# ----------------------------------------------------------------------------
