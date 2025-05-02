

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export

#'@export
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




