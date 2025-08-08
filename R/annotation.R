#' @title Retrieve Protein Annotations from UniProt
#' @description This function fetches detailed protein annotations from the UniProt
#' database for a given list of protein identifiers. It handles protein accessions
#' that may contain multiple IDs (e.g., from protein grouping), queries UniProt,
#' processes the results, and caches them locally to an RDS file to speed up
#' subsequent runs.
#'
#' @details
#' The function performs the following steps:
#' 1. Extracts and cleans protein accession numbers from the specified column.
#'    It can handle delimited strings of multiple accessions in a single row.
#' 2. Checks for a local cache file (`uniprot_dat.rds`) in `output_dir`. If found,
#'    it loads and returns the cached data.
#' 3. If no cache is found, it connects to UniProt via `UniProt.ws` using the
#'    provided taxonomy ID.
#' 4. It queries UniProt for a predefined set of annotation columns, including
#'    gene names, protein names, GO IDs, and keywords.
#' 5. It processes the Gene Ontology (GO) IDs to retrieve the full GO terms.
#' 6. It aggregates the annotations back to match the original input protein IDs,
#'    concatenating results with a colon where one ID mapped to multiple entries.
#' 7. The final annotation table is saved to the cache file for future use.
#'
#' @param input_table A data frame or tibble containing the protein identifiers.
#' @param taxonomy_id A numeric value specifying the NCBI taxonomy ID of the organism.
#'   Defaults to 9606 (Homo sapiens).
#' @param protein_id_column A character string naming the column in `input_table`
#'   that contains the protein identifiers. Defaults to "Protein.Ids".
#' @param protein_id_delim A character string used as a delimiter for splitting
#'   entries in `protein_id_column` that contain multiple accessions. Defaults to ":".
#' @param output_dir A character string specifying the directory path where the
#'   UniProt annotation data should be cached as "uniprot_dat.rds". Defaults to the
#'   current working directory (".").
#'
#' @return A data frame where each row corresponds to an entry in the original
#'   `protein_id_column`, and columns contain the fetched and aggregated annotations
#'   from UniProt.
#'
#' @importFrom UniProt.ws UniProt.ws columns keytypes
#' @importFrom AnnotationDbi select
#' @importFrom GO.db GOTERM
#' @importFrom dplyr mutate distinct group_by ungroup left_join arrange summarise rename
#' @importFrom tidyr separate_rows
#' @importFrom magrittr %>%
#' @importFrom rlang !! sym
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create a dummy data frame with human protein IDs
#' my_proteins <- data.frame(
#'   Protein.Ids = c("P04637", "P02768;Q02880"),
#'   other_data = c(10, 20)
#' )
#'
#' # Create a temporary output directory for caching
#' temp_dir <- tempdir()
#'
#' # Get annotations
#' annotations <- getUniProtAnnotation(
#'   input_table = my_proteins,
#'   taxonomy_id = 9606,
#'   protein_id_column = "Protein.Ids",
#'   protein_id_delim = ";",
#'   output_dir = temp_dir
#' )
#'
#' # View the first few columns of the result
#' print(head(annotations))
#' }
getUniProtAnnotation <- function(   input_table, taxonomy_id  =9606, protein_id_column = "Protein.Ids",  protein_id_delim=":", output_dir = "." ) {


  uniprot_dat <- NA
  uniprot_acc_tbl <- input_table |>
    mutate( uniprot_acc_copy = !!sym(protein_id_column) )  |>
    separate_rows(uniprot_acc_copy, sep=protein_id_delim )  |>
    mutate( join_uniprot_acc = cleanIsoformNumber(uniprot_acc_copy))  |>
    dplyr::distinct( !!sym(protein_id_column) , join_uniprot_acc)  |>
    group_by( !!sym(protein_id_column) )  |>
    mutate( acc_order_id = row_number())  |>
    ungroup()

  uniprot_file<-file.path( output_dir, "uniprot_dat.rds")
  if( ! file.exists( uniprot_file )) {

    up <- UniProt.ws(taxId= taxonomy_id)
    list_of_sp_columns <- c("EXISTENCE"
                            , "SCORE"
                            , "REVIEWED"
                            , "GENENAME"
                            , "PROTEIN-NAMES"
                            , "LENGTH"
                            , "ENSEMBL"
                            , "GO-ID"
                            , "KEYWORDS"

                            ,"protein_existence"
                            ,"annotation_score"#?
                            ,"reviewed"
                            ,"gene_names"
                            ,"protein_name"
                            ,"length"
                            ,"xref_ensembl"
                            , "go_id"
                            , "keyword"
    )
    up_cls<-unlist(columns(up))
    list_intersect<-intersect(list_of_sp_columns,up_cls)
    if(length(setdiff( list_of_sp_columns,list_intersect)) > 0)
    {
      print(paste("UniProt fields not found:", setdiff( list_of_sp_columns,list_intersect),sep=", "))
    }

    my_keytype <- "UniProtKB"
    if( "UNIPROTKB" %in% keytypes(up) ) {
      my_keytype <- "UNIPROTKB"
    }

    uniprot_dat <- batchQueryEvidence(uniprot_acc_tbl, join_uniprot_acc, uniprot_handle=up,
                                      uniprot_columns = list_intersect, uniprot_keytype=my_keytype)

    if( my_keytype == "UniProtKB") {
      uniprot_dat <- uniprot_dat %>%
        dplyr::select(-From) %>%
        dplyr::rename( UNIPROTKB = "Entry",
                       EXISTENCE = "Protein.existence",
                       SCORE = "Annotation",
                       REVIEWED = "Reviewed",
                       GENENAME = "Gene.Names",
                       `PROTEIN-NAMES` = "Protein.names",
                       LENGTH = "Length",
                       ENSEMBL = "Ensembl",
                       `GO-ID` = "Gene.Ontology.IDs",
                       KEYWORDS   = "Keywords")


    }


    ## Merge with Gene Ontology terms.
    goterms <- Term(GOTERM)
    gotypes <- Ontology(GOTERM)


    uniprot_dat_cln <- uniprotGoIdToTerm(uniprot_dat, sep="; ", goterms, gotypes  )


    uniprot_dat_multiple_acc <- uniprot_acc_tbl %>%
      left_join( uniprot_dat_cln, by=c("join_uniprot_acc" = "UNIPROTKB") ) %>%
      arrange( !!sym(protein_id_column) , acc_order_id) %>%
      group_by(!!sym(protein_id_column)  ) %>%
      summarise( across( .cols=setdiff( colnames( uniprot_dat_cln), "UNIPROTKB")   , ~paste(., collapse=":"))   ) %>%
      ungroup() %>%
      dplyr::rename( UNIPROT_GENENAME = "GENENAME")

    saveRDS( uniprot_dat_multiple_acc, uniprot_file)

    return( uniprot_dat_multiple_acc)

  } else {
    uniprot_dat <- readRDS(uniprot_file)

    return( uniprot_dat)

  }

}

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




