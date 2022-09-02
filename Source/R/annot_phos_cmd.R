#!/usr/bin/env Rscript

# Author(s): Ignatius Pang, Pablo Galaviz
# Email: cmri-bioinformatics@cmri.org.au
# Children’s Medical Research Institute, finding cures for childhood genetic diseases


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Test if BioManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
   BiocManager::install(version = "3.13")
}

# load pacman package manager
if(!require(pacman)){
    install.packages("pacman")
    library(pacman)
}

p_load(tidyverse)
p_load(vroom)
p_load(magrittr)
p_load(rlang)
p_load(ProteomeRiver)

p_load(UniProt.ws)
p_load(biomaRt)
p_load(GO.db)

p_load(optparse)
p_load(tictoc)
p_load(configr)
p_load(logging)

tic()


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

command_line_options <- commandArgs(trailingOnly = TRUE)
parser <- OptionParser(add_help_option = TRUE)
#Note: options with default values are ignored in the configuration file parsing.
parser <- add_option(parser, c("-d", "--debug"), action = "store_true", default = FALSE,
                     help = "Print debugging output")

parser <- add_option(parser, c("-s", "--silent"), action = "store_true", default = FALSE,
                     help = "Only print critical information to the console.")

parser <- add_option(parser, c("-n", "--no_backup"), action = "store_true", default = FALSE,
                     help = "Deactivate backup of previous run.")

parser <- add_option(parser, c("-c","--config"), type = "character", default = "config_phos.ini", dest = "config",
                     help = "Configuration file.",
                     metavar = "string")

parser <- add_option(parser, c("-o","--output_dir"), type = "character", dest = "output_dir",
                     help = "Directory path for all results files.",
                     metavar = "string")

parser <- add_option(parser, c("-t","--tmp_dir"), type = "character", default = "cache", dest = "tmp_dir",
                     help = "Directory path for temporary files.",
                     metavar = "string")

parser <- add_option(parser, c("-l","--log_file"), type = "character", default = "output.log", dest = "log_file",
                     help = "Name of the logging file.",
                     metavar = "string")

#Options without a default value have the following priority: configuration file < command line argument
parser <- add_option(parser,  "--taxonomy_id", type="integer", dest = "taxonomy_id",
                     help="The NCBI taxonomy ID of the organism being investigated (e.g. M. musculus=10090, H. sapien=9606).",
                     metavar="integer")

parser <- add_option(parser, "--input_wide_file", type="character", dest = "input_wide_file",
                       help="Results table with values in wider format.",
                       metavar="string")

parser <- add_option(parser, "--input_long_file", type="character",  dest = "input_long_file",
                     help="Results table with values in longer format.",
                     metavar="string")

parser <- add_option(parser, "--raw_counts_file", type="character", dest = "raw_counts_file",
                     help="Input file with the protein abundance data.",
                     metavar="string")

parser <- add_option(parser,  "--output_wide_file", type="character",  dest = "output_wide_file",
                     help="Results table with values in wider format.",
                     metavar="string")

parser <- add_option(parser,  "--output_long_file", type="character", dest = "output_long_file",
                     help="Results table with values in longer format.",
                     metavar="string")

parser <- add_option(parser,  "--phosphosite_db_dir", type="character", dest = "phosphosite_db_dir",
                     help="File path for location of PhosphoSitePlus DB data files (https://www.phosphosite.org/staticDownloads).",
                     metavar="string")


parser <- add_option(parser, c( "--near_ptm_num_residues"), type="integer", dest = "near_ptm_num_residues",
                     help="Find other PTM near within +/- this many number of residues from the phosphosites.",
                     metavar="integer")

parser <- add_option(parser,  "--reactome_file", type="character",  dest = "reactome_file",
                     help="Name of the reactome data (Download and save if it does not exists). ",
                     metavar="string")

parser <- add_option(parser,  "--uniprot_file", type="character",  dest = "uniprot_file",
                     help="Name of the uniprot data (Download and save if it does not exists). ",
                     metavar="string")

#parse command line arguments first.
args <- parse_args(parser)

#parse and merge the configuration file options.
if (args$config != "") {
   args <- config.list.merge(eval.config(file = args$config, config = "annot_phos"), args)
}

args <- setArgsDefault(args, "output_dir", as_func=as.character, default_val="annot_phos" )

createOutputDir(args$output_dir, args$no_backup)
createDirectoryIfNotExists(args$tmp_dir)

## Logger configuration
logReset()
addHandler(writeToConsole , formatter = cmriFormatter)
addHandler(writeToFile, file = file.path(args$output_dir, args$log_file), formatter = cmriFormatter)

level <- ifelse(args$debug, loglevels["DEBUG"], loglevels["INFO"])
setLevel(level = ifelse(args$silent, loglevels["ERROR"], level))

#parse and merge the configuration file options.
if (args$config != "") {
  args <- config.list.merge(eval.config(file = args$config, config = "annot_phos"), args)
}

cmriWelcome("ProteomeRiver", c("Ignatius Pang", "Pablo Galaviz"))
loginfo("Reading configuration file %s", args$config)
loginfo("Argument: Value")
loginfo("----------------------------------------------------")
for (v in names(args))
{
  loginfo("%s : %s", v, args[v])
}
loginfo("----------------------------------------------------")

testRequiredArguments(args, c(
  "input_long_file"
  ,"raw_counts_file"
  ,"output_wide_file"
  ,"output_long_file"
  ,"taxonomy_id"
  ,"phosphosite_db_dir"
  ,"near_ptm_num_residues"
  ,"reactome_file"
  ,"uniprot_file"
))

ks_file <-  file.path( args$phosphosite_db_dir, "Kinase_Substrate_Dataset")
reg_sites_file <- file.path( args$phosphosite_db_dir, "Regulatory_sites" )
disease_file <- file.path( args$phosphosite_db_dir, "Disease-associated_sites")

testRequiredFiles(c(
  args$input_long_file
  , args$raw_counts_file
  ,args$ids_file
  ,ks_file
  ,reg_sites_file
  ,disease_file
))

list_of_ptm_file_names <- c("Acetylation_site_dataset",
"Methylation_site_dataset",
"O-GalNAc_site_dataset",
"O-GlcNAc_site_dataset",
"Phosphorylation_site_dataset",
"Sumoylation_site_dataset",
"Ubiquitination_site_dataset")

list_of_ptm_files <- purrr::map_chr( list_of_ptm_file_names, ~file.path(args$phosphosite_db_dir, .))
testRequiredFiles(list_of_ptm_files)


#unused
#args<-parseType(args,
#  c("some_double_par")
#  ,as.double)

args<-parseType(args,
  c("taxonomy_id"
  ,"near_ptm_num_residues")
                ,as.integer)



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Read abundance table.")
captured_output<-capture.output(
  abundance_tbl <- vroom::vroom(args$raw_counts_file)
  ,type = "message"
)
logdebug(captured_output)

loginfo("Read differentially abundant phosphopeptides table in long format")
captured_output<-capture.output(
  de_phos_long <- vroom::vroom( args$input_long_file) %>%
  mutate( sites_id_copy = sites_id, .after="sites_id") %>%
  separate( sites_id_copy, sep="!", into=c("uniprot_acc", "gene_name", "position", "sequence")) %>%
  mutate( residue= purrr::map_chr( sequence, ~{ str_replace_all( ., "[A-Z_]{7}(.)[A-Z_]{7}([\\:;\\|]*)", "\\1\\2"    )  }  )) %>%
  dplyr::relocate(residue, .before="position")
  ,type = "message"
)
logdebug(captured_output)

loginfo("Read differentially abundant phosphopeptides table in wide format")
captured_output<-capture.output(
  de_phos_wide <- vroom::vroom( args$input_wide_file) %>%
  mutate( sites_id_copy = sites_id, .after="sites_id") %>%
  separate( sites_id_copy, sep="!", into=c("uniprot_acc", "gene_name", "position", "sequence")) %>%
  mutate( residue= purrr::map_chr( sequence, ~{ str_replace_all( ., "[A-Z_]{7}(.)[A-Z_]{7}([\\:;\\|]*)", "\\1\\2"    )  }  )) %>%
  dplyr::relocate(residue, .before="position")
  ,type = "message"
)
logdebug(captured_output)



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Read PhosphoSitePlus (PSP) kinase-substrate table.")
captured_output<-capture.output(
  ks_tbl <- vroom::vroom( ks_file, skip=3 ) %>%
  mutate( SUB_MOD_RSD_CLN = str_replace_all(SUB_MOD_RSD, "([A-Z])(\\d+)", "\\1 \\2")) %>%
  separate( SUB_MOD_RSD_CLN, into=c("residue", "position"))
  ,type = "message"
)
logdebug(captured_output)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Read PSP regulatory sites table.")
captured_output<-capture.output(
  reg_sites_tbl <- vroom::vroom( reg_sites_file, skip=3)  %>%
      mutate( MOD_RSD_CLN=  str_replace_all(  MOD_RSD, "([A-Z])(\\d+)-(.*)", "\\1 \\2 \\3") ) %>%
      separate( MOD_RSD_CLN, sep=" ", into=c("residue", "position", "ptm_type")) %>%
      relocate (ACC_ID, residue, position, ptm_type, .before="GENE" ) %>%
      dplyr::select(-`...21`) %>%
      dplyr::rename( REG_SITES_PMIDs = "PMIDs") %>%
      dplyr::rename( REG_SITES_NOTES = "NOTES")
  ,type = "message"
)
logdebug(captured_output)



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo( "Read PSP disease table")
captured_output<-capture.output(
  disease_tbl <- vroom::vroom( disease_file, skip=3 ) %>%
  dplyr::select( ACC_ID, MOD_RSD, DISEASE, ALTERATION, NOTES ) %>%
  mutate( MOD_RSD_CLN=  str_replace_all(  MOD_RSD, "([A-Z])(\\d+)-(.*)", "\\1 \\2 \\3") ) %>%
  separate( MOD_RSD_CLN, sep=" ", into=c("residue", "position", "ptm_type")) %>%
  dplyr::rename( DISEASE_NOTES = "NOTES") %>%
  dplyr::select(-MOD_RSD, -ptm_type)
  ,type = "message"
)
logdebug(captured_output)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

loginfo("Read PSP post-translational modification (PTM) tables")
logdebug(list_of_ptm_files)
captured_output<-capture.output(
  ptm_tbl <- vroom::vroom( list_of_ptm_files, skip=3, id="ptm")  %>%
  dplyr::select( ACC_ID, MOD_RSD, ptm) %>%
  mutate( MOD_RSD_CLN=  str_replace_all(  MOD_RSD, "([A-Z])(\\d+)-(.*)", "\\1 \\2 \\3") ) %>%
  separate( MOD_RSD_CLN, sep=" ", into=c("residue", "position", "ptm_type")) %>%
  mutate( uniprot_acc = ACC_ID) %>%
  dplyr::mutate( position = as.integer(position)) %>%
  mutate( ptm= purrr::map_chr( ptm, ~{ temp_vec <- str_split(., "/")[[1]]
  temp_vec[length(temp_vec)] } )) %>%
  dplyr::mutate( ptm = str_replace_all( ptm, "_site_dataset", ""))
  ,type = "message"
)
logdebug(captured_output)



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo( "Find other PTM nearby +/- %s amino acid wihtin the phosphorylation sites.", args$near_ptm_num_residue)
get_nearby_ptm <- de_phos_long  %>%
  dplyr::distinct( sites_id, uniprot_acc,  position) %>%
  separate_rows( uniprot_acc, position, sep="\\:") %>%
  separate_rows( uniprot_acc, position, sep="\\|") %>%
  separate_rows( uniprot_acc, position, sep=";") %>%
  mutate( position = str_replace_all( position, "\\(|\\)", "") %>% purrr::map_int(as.integer)) %>%
  mutate( residue_window =  purrr::map(position,  ~seq( from=as.integer( .)-args$near_ptm_num_residues, to= as.integer( .)+args$near_ptm_num_residues, by=1 )) ) %>%
  unnest( residue_window) %>%
  left_join(ptm_tbl %>%
              mutate(ptm_count=1) %>%
              dplyr::select(-residue), by=c( "uniprot_acc" = "uniprot_acc",
                                 "residue_window" = "position")) %>%
  distinct %>%
  dplyr::filter( !( residue_window == position & ptm_type == "p") ) %>% ## Avoid counting the same phosphorylation site as the query itself
  dplyr::select(-MOD_RSD, -uniprot_acc, -position, -ptm_type)

nearby_ptm_count <- get_nearby_ptm %>%
  dplyr::select(-ACC_ID) %>%
  dplyr::filter( !is.na(ptm_count)) %>%
  group_by( across(.cols=setdiff(colnames(get_nearby_ptm), c("ptm_count", "residue_window", "ACC_ID") ))) %>%
  distinct() %>%
  summarise( ptm_count = sum(ptm_count) ) %>%
  ungroup  %>%
  distinct %>%
  pivot_wider( id_cols = sites_id,
               values_from = ptm_count,
               names_from="ptm",
               names_prefix = paste0("nearby_+/-", args$near_ptm_num_residues, "_") )




## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
num_phos_sites <- de_phos_long  %>%
  dplyr::distinct( sites_id, position) %>%
  dplyr::mutate( num_sites =   str_split(position, ":") %>%
                   purrr::map_chr(1) %>%
                   str_split( "\\|")  %>%
                 purrr::map_chr(1) %>%
                   str_split(";") %>%
                   purrr::map_int(length) ) %>%
  dplyr::select(-position)




## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# de_phos_long %>% dplyr::filter( str_detect( uniprot_acc, "[^[:alnum:]]+" ) )

phosphosite_plus_tbl <- de_phos_long  %>%
  dplyr::distinct( sites_id, uniprot_acc, residue, position, sequence) %>%
  separate_rows( uniprot_acc, position, sequence, residue, sep="\\:") %>%
  separate_rows( uniprot_acc, position, sequence, residue, sep="\\|") %>%
  separate_rows( uniprot_acc, position, sequence, residue, sep=";") %>%
  mutate( position = str_replace_all( position, "\\(|\\)", "") %>% purrr::map_int(as.integer)) %>%
  left_join (reg_sites_tbl %>%
               dplyr::filter(ptm_type == "p") %>%
               dplyr::mutate( position = purrr::map_int(position, as.integer)),
             by=c("uniprot_acc" = "ACC_ID",
                                 "residue" = "residue",
                                 "position" = "position"
                                 )) %>%
  left_join( ks_tbl %>%
               dplyr::rename(KINASE_GENE = "GENE") %>%
               dplyr::select(-DOMAIN, - `SITE_+/-7_AA`) %>%
               dplyr::mutate( position = purrr::map_int(position, as.integer)),
             by=c("uniprot_acc" = "SUB_ACC_ID",
                          "residue"="residue",
                          "position" = "position",
                          "SITE_GRP_ID" = "SITE_GRP_ID"))   %>%
  left_join( disease_tbl %>%
               dplyr::mutate( position = purrr::map_int(position, as.integer)),
             by=c( "uniprot_acc" = "ACC_ID",
                                 "residue" = "residue",
                                 "position" = "position")) %>%
  group_by(sites_id, uniprot_acc ) %>%
  summarise( across( .cols=everything()  , ~paste(unique(.), collapse="//"))   ) %>%
  ungroup() %>%
  group_by(sites_id ) %>%
  summarise( across( .cols=everything()  , ~paste(unique(.), collapse=":"))   ) %>%
  ungroup()

# de_phos_long_annot %>% dplyr::filter( str_detect( KINASE, "//" ) )
#
# # colnames( de_phos_long_annot)[ which( !is.na( str_match( colnames( de_phos_long_annot), "\\.y" ) )  ) ]
#
# de_phos_long_annot %>%
#   dplyr::filter( !is.na(ptm_type) | !is.na( KINASE))



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

reactome_file<-file.path(args$tmp_dir,args$reactome_file)
if(!file.exists(reactome_file))
{
  logwarn("Download Reactome UniProt to pathways file.")
  status <- download.file(url="https://reactome.org/download/current/UniProt2Reactome.txt", destfile=reactome_file)
  loginfo(status)
}
loginfo("Reading Reactome UniProt to pathways file.")
captured_output<-capture.output(
  reactome_map <- vroom::vroom( reactome_file ,
                              col_names = c("uniprot_acc", "reactome_id", "url", "reactome_term", "evidence", "organism") )
  ,type = "message"
)
logdebug(captured_output)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

loginfo("Get the best UniProt accession per row.")

uniprot_acc_tbl <- de_phos_long %>%
  mutate( uniprot_acc_copy = uniprot_acc ) %>%
  separate_rows(uniprot_acc_copy, sep=":" ) %>%
  mutate( join_uniprot_acc = cleanIsoformNumber(uniprot_acc_copy)) %>%
  dplyr::distinct( uniprot_acc, join_uniprot_acc) %>%
  group_by( uniprot_acc) %>%
  mutate( acc_order_id = row_number()) %>%
  ungroup




## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Download information from UniProt.")
uniprot_file<-file.path(args$tmp_dir,args$uniprot_file)
if( ! file.exists( uniprot_file )) {

  up <- UniProt.ws(taxId=args$taxonomy_id )
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
    logwarn("UniProt fields not found: %s",paste(setdiff( list_of_sp_columns,list_intersect),sep=", "))
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

  saveRDS( uniprot_dat, uniprot_file)

}

loginfo("Reading UniProt data from %s.",uniprot_file)
uniprot_dat <- readRDS( uniprot_file )



## ----eval=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##  # uniprot_dat <- batch_query_evidence(uniprot_acc_tbl %>% head(100), best_uniprot_acc, uniprot_handle=up,
##  #                                      uniprot_columns = list_of_sp_columns)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

loginfo("Merge with Gene Ontology terms.")
goterms <- Term(GOTERM)
gotypes <- Ontology(GOTERM)

uniprot_dat_cln <- uniprotGoIdToTerm(uniprot_dat, sep="; ", goterms, gotypes  )

uniprot_dat_multiple_acc <- uniprot_acc_tbl %>%
  left_join( uniprot_dat_cln, by=c("join_uniprot_acc" = "UNIPROTKB") )   %>%
  arrange( uniprot_acc, acc_order_id) %>%
  group_by(uniprot_acc ) %>%
  summarise( across( .cols=setdiff( colnames( uniprot_dat_cln), "UNIPROTKB")   , ~paste(., collapse=":"))   ) %>%
  ungroup()   %>%
  dplyr::rename( UNIPROT_GENENAME = "GENENAME")



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Add reactome pathways annotation.")
reactome_term_tbl <- uniprot_acc_tbl %>%
  left_join( reactome_map, by=c("join_uniprot_acc" = "uniprot_acc") )   %>%
  dplyr::filter(reactome_term != "NA" ) %>%
  group_by(uniprot_acc, join_uniprot_acc) %>%
  summarise( reactome_term = paste(reactome_term, collapse="; ") ) %>%
  ungroup()     %>%
  mutate(reactome_term = str_replace_all( reactome_term , ":", "-")) %>%
  group_by(uniprot_acc ) %>%
  summarise( reactome_term = paste(reactome_term, collapse=":") ) %>%
  ungroup()



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# print("Output wider format results table with protein annotation.")
# de_phos_wide_annot <- de_phos_wide %>%
#   left_join( num_phos_sites, by =c("sites_id" = "sites_id")) %>%
#   left_join( uniprot_dat_multiple_acc, by = c("uniprot_acc" = "uniprot_acc") ) %>%
#   left_join( reactome_term_tbl, by = c("uniprot_acc" = "uniprot_acc"))  %>%
#   left_join( phosphosite_plus_tbl %>%
#                dplyr::select(-uniprot_acc, -position, -residue, -sequence),
#              by=c("sites_id" = "sites_id")) %>%
#   left_join(  abundance_tbl %>%
#    dplyr::select( sites_id, maxquant_row_ids ),
#    by=c("sites_id" = "sites_id")) %>%
#   relocate(maxquant_row_ids, .after="sites_id") %>%
#   left_join( get_nearby_ptm,
#              by=c("sites_id" = "sites_id")) %>%
#   arrange( comparison, q.mod, log2FC) %>%
#   distinct()
#
# head( de_phos_wide_annot )
#
# vroom::vroom_write(de_phos_wide_annot, path=file.path(args$output_dir,args$output_wide_file )



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Output longer format results table with protein annotation.")
de_phos_long_annot <- de_phos_long %>%
  left_join( num_phos_sites, by =c("sites_id" = "sites_id")) %>%
  left_join( uniprot_dat_multiple_acc, by = c("uniprot_acc" = "uniprot_acc") ) %>%
  left_join( reactome_term_tbl, by = c("uniprot_acc" = "uniprot_acc"))  %>%
  left_join( phosphosite_plus_tbl %>%
               dplyr::select(-uniprot_acc, -position, -residue, -sequence),
             by=c("sites_id" = "sites_id"))  %>%
  left_join(  abundance_tbl %>%
   dplyr::select( sites_id, maxquant_row_ids ),
   by=c("sites_id" = "sites_id")) %>%
  relocate(maxquant_row_ids, .after="sites_id")  %>%
  left_join( nearby_ptm_count,
             by=c("sites_id" = "sites_id")) %>%
  arrange( comparison, q.mod, log2FC) %>%
  distinct()

vroom::vroom_write(de_phos_long_annot, file.path(args$output_dir,args$output_long_file ))


list_of_long_columns <- intersect(colnames(de_phos_long_annot), c("protein_names",
                                                                  "ENSEMBL",
                                                                  "PROTEIN-NAMES",
                                                                  "KEYWORDS",
                                                                  "GO-ID",
                                                                  "go_biological_process",
                                                                  "go_cellular_compartment",
                                                                  "go_molecular_function",
                                                                  "reactome_term",
                                                                  "majority_protein_ids") )

writexl::write_xlsx(de_phos_long_annot %>%
                      mutate_at( list_of_long_columns, ~substr(., 1, 32760) ),
                     file.path(args$output_dir,  str_replace(args$output_long_file, "\\..*", ".xlsx")  ))


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
te<-toc(quiet = TRUE)
loginfo("%f sec elapsed",te$toc-te$tic)
writeLines(capture.output(sessionInfo()), file.path(args$output_dir,"sessionInfo.txt"))

