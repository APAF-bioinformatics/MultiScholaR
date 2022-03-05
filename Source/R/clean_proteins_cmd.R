#!/usr/bin/env Rscript

# Author(s): Ignatius Pang, Pablo Galaviz
# Email: cmri-bioinformatics@cmri.org.au
# Children’s Medical Research Institute, finding cures for childhood genetic diseases


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Test if BioManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install()
}

# load pacman package manager
if (!require(pacman)) {
  install.packages("pacman")
  library(pacman)
}

p_load(tidyverse)
p_load(vroom)
p_load(magrittr)
p_load(knitr)
p_load(rlang)
p_load(optparse)

p_load(seqinr)
p_load(ProteomeRiver)
p_load(janitor)
p_load(tictoc)
p_load(configr)
p_load(logging)


tic()


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

command_line_options <- commandArgs(trailingOnly = TRUE)
parser <- OptionParser(add_help_option = TRUE)
#Note: options with default values are ignored in the configuration file parsing.
parser <- add_option(parser, c("-d", "--debug"), action = "store_true", default = FALSE,
                     help = "Print debugging output")

parser <- add_option(parser, c("-s", "--silent"), action = "store_true", default = FALSE,
                     help = "Only print critical information to the console.")

parser <- add_option(parser, c("-n", "--no_backup"), action = "store_true", default = FALSE,
                     help = "Deactivate backup of previous run.")

parser <- add_option(parser, c("-c","--config"), type = "character", default = "config.ini", dest = "config",
                     help = "Configuration file.",
                     metavar = "string")

parser <- add_option(parser, c("-o","--output_dir"), type = "character",  dest = "output_dir",
                     help = "Directory path for all results files.",
                     metavar = "string")

parser <- add_option(parser, c("-t","--tmp_dir"), type = "character", default = "cache", dest = "tmp_dir",
                     help = "Directory path for temporary files.",
                     metavar = "string")

parser <- add_option(parser, c("-l","--log_file"), type = "character", default = "output.log", dest = "log_file",
                     help = "Name of the logging file.",
                     metavar = "string")

#Options without a default value have the following priority: configuration file < command line argument
parser <- add_option(parser, "--accession_record_file", type = "character", dest = "accession_record_file",
                     help = "File to link the cleaned list of accessions to the original list of protein groups from MaxQuant file. Also, contains the MaxQuant output ID.",
                     metavar = "string")

parser <- add_option(parser, "--fasta_file", type = "character", dest = "fasta_file",
                     help = "Input protein sequence FASTA file with UniProt FASTA header format.",
                     metavar = "string")

parser <- add_option(parser, "--raw_counts_file", type = "character", dest = "raw_counts_file",
                     help = "Input file with the protein abundance data.",
                     metavar = "string")

parser <- add_option(parser, "--output_counts_file", type = "character", dest = "output_counts_file",
                     help = "String representing the name of the output counts data file which will be saved in the directory specified with the --output-dir flag.",
                     metavar = "string")

parser <- add_option(parser, "--column_pattern", type = "character", dest = "column_pattern",
                     help = "String pattern, together with the experimental group pattern, that matches the abundance data columns.",
                     metavar = "string")

parser <- add_option(parser, "--group_pattern", type = "character", dest = "group_pattern",
                     help = "Regular expression pattern to identify columns with abundance values belonging to the experiment. [default %default]",
                     metavar = "string")

parser <- add_option(parser, "--pattern_suffix", type = "character", dest = "pattern_suffix",
                     help = "Regular expression pattern to identify columns with abundance values belonging to the experiment. [default %default]",
                     metavar = "string")

parser <- add_option(parser, "--extract_patt_suffix", type = "character", dest = "extract_patt_suffix",
                     help = "Regular expression pattern to identify columns with abundance values belonging to the experiment. [default %default]",
                     metavar = "string")


parser <- add_option(parser, "--razor_unique_peptides_group_thresh", type = "integer", dest = "razor_unique_peptides_group_thresh",
                     help = "Number of razor and unique peptides for the specified experiemtal group needs to be higher than this threshold for the protein to be included for the analysis.",
                     metavar = "integer")

parser <- add_option(parser, "--unique_peptides_group_thresh", type = "integer", dest = "unique_peptides_group_thresh",
                     help = "Number of unique peptides for the specified experiemtal group needs to be higher than this threshold for the protein to be included for the analysis.\nThe file will be saved in the results directory",
                     metavar = "integer")

parser <- add_option(parser, "--fasta_meta_file", type = "character", dest = "fasta_meta_file",
                     help = "R object storage file that records all the sequence and header information in the FASTA file.",
                     metavar = "string")

#parse comand line arguments first.
args <- parse_args(parser)

#parse and merge the configuration file options.
if (args$config != "") {
  args <- config.list.merge(eval.config(file = args$config, config = "clean_proteins"), args)
}

args <- setArgsDefault(args, "output_dir", as_func=as.character, default_val="clean_proteins" )


createOutputDir(args$output_dir, args$no_backup)
createDirectoryIfNotExists(args$tmp_dir)

## Logger configuration
logReset()
addHandler(writeToConsole , formatter = cmriFormatter)
addHandler(writeToFile, file = file.path(args$output_dir, args$log_file), formatter = cmriFormatter)

level <- ifelse(args$debug, loglevels["DEBUG"], loglevels["INFO"])
setLevel(level = ifelse(args$silent, loglevels["ERROR"], level))




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
  "output_counts_file"
  ,"razor_unique_peptides_group_thresh"
  ,"unique_peptides_group_thresh"
  ,"fasta_meta_file"
  ,"group_pattern"
  ,"accession_record_file"
  ,"fasta_file"
  ,"raw_counts_file"
))

testRequiredFiles(c(
  args$fasta_file
  ,args$raw_counts_file))

args <- setArgsDefault(args, "pattern_suffix", as_func=as.character, default_val="_\\d+" )
args <- setArgsDefault(args, "extract_patt_suffix", as_func=as.character, default_val="_(\\d+)" )


args<-parseType(args,
  c("razor_unique_peptides_group_thresh"
    ,"unique_peptides_group_thresh"
  ),as.integer)

args<-parseString(args,
  c("group_pattern"
    ,"column_pattern"
    ,"pattern_suffix"
    ,"extract_patt_suffix" ))



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Reading the counts file.")

captured_output<-capture.output(
  dat_tbl <- vroom::vroom(args$raw_counts_file)
  ,type = "message"
)
logdebug(captured_output)

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Clean counts table header name.")
dat_cln <- janitor::clean_names(dat_tbl)

colnames(dat_cln) <- str_replace(colnames(dat_cln), "_i_ds", "_ids")


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Prepare regular expressions to select counts data columns.")

pattern_suffix <- args$pattern_suffix
if (args$group_pattern != "") {
  pattern_suffix <- paste(args$pattern_suffix, tolower(args$group_pattern), sep = "_")
}

extract_patt_suffix <- args$extract_patt_suffix
if (args$group_pattern != "") {
  extract_patt_suffix <- paste0(args$extract_patt_suffix, "_(", tolower(args$group_pattern), ")")
}

column_pattern <- tolower(paste0(make_clean_names(args$column_pattern), pattern_suffix))

extract_replicate_group <- tolower(paste0(make_clean_names(args$column_pattern), extract_patt_suffix))


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Prepare column names associated with peptide counts.")

razor_unique_peptides_group_col <- "razor_unique_peptides"
unique_peptides_group_col <- "unique_peptides"

if("group_pattern" %in% args) {
  if (args$group_pattern != "") {
    razor_unique_peptides_group_col <- paste0("razor_unique_peptides_", tolower(args$group_pattern))
    unique_peptides_group_col <- paste0("unique_peptides_", tolower(args$group_pattern))
  }
}

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Read fasta file
# colnames(aa_seq_tbl)

fasta_meta_file<-file.path(args$tmp_dir,args$fasta_meta_file)
loginfo("Reading the FASTA file %s and saving the meta-data file.",fasta_meta_file)
if (file.exists(fasta_meta_file))
{
  loginfo("Mata-data exists, reading RDS.")
  aa_seq_tbl <- readRDS(fasta_meta_file)
}else{
  loginfo("Mata-data does not exists, parsing fasta and saving RDS.")
  aa_seq_tbl <- parseFastaFile(args$fasta_file)
  saveRDS(aa_seq_tbl, fasta_meta_file)
}

## Add the row id column and create a column containing the cleaned  peptide
loginfo("Get row ID and get cleaned peptide sequence.")
evidence_tbl <- dat_cln %>%
  mutate(maxquant_row_id = id)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Filtering counts table by number of peptides available, remove decoy proteins and protein contaminants.")

peptides_count_helper <- evidence_tbl %>%
  dplyr::select(maxquant_row_id,
                protein_ids,
                !!rlang::sym(razor_unique_peptides_group_col),
                !!rlang::sym(unique_peptides_group_col),
                reverse,
                potential_contaminant,
                matches(column_pattern)) %>%
  dplyr::filter(is.na(reverse) &
                  is.na(potential_contaminant)) %>%
  dplyr::filter(!str_detect(protein_ids, "CON__") &
                  !str_detect(protein_ids, "REV__")) %>%
  dplyr::mutate(protein_ids = str_split(protein_ids, ";")) %>%
  dplyr::mutate(!!rlang::sym(razor_unique_peptides_group_col) := str_split(!!rlang::sym(razor_unique_peptides_group_col), ";")) %>%
  dplyr::mutate(!!rlang::sym(unique_peptides_group_col) := str_split(!!rlang::sym(unique_peptides_group_col), ";")) %>%
  unnest(cols = c(protein_ids,
                  !!rlang::sym(razor_unique_peptides_group_col),
                  !!rlang::sym(unique_peptides_group_col)))


evidence_tbl_cleaned <- NA

evidence_tbl_cleaned <- peptides_count_helper %>%
  dplyr::filter(!!rlang::sym(razor_unique_peptides_group_col) >= args$razor_unique_peptides_group_thresh &
                  !!rlang::sym(unique_peptides_group_col) >= args$unique_peptides_group_thresh)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loginfo("Identify best UniProt accession per entry, extract sample number and simplify column header")


accession_gene_name_tbl <- chooseBestProteinAccession(input_tbl = evidence_tbl_cleaned,
                                                      acc_detail_tab = aa_seq_tbl,
                                                      accessions_column = protein_ids,
                                                      row_id_column = uniprot_acc,
                                                      group_id = maxquant_row_id)


accession_gene_name_tbl_record <- accession_gene_name_tbl %>%
  left_join(evidence_tbl %>% dplyr::select(maxquant_row_id, protein_ids), by = c("maxquant_row_id"))


evidence_tbl_filt <- NA

evidence_tbl_filt <- evidence_tbl_cleaned %>%
  inner_join(accession_gene_name_tbl %>%
               dplyr::select(maxquant_row_id, uniprot_acc), by = "maxquant_row_id") %>%
  dplyr::select(uniprot_acc, matches(column_pattern), -contains(c("razor", "unique"))) %>%
  distinct

#TODO: This part need improvement. There is potential for bugs.
extraction_pattern <- "\\1"
if (args$group_pattern != "") {
  extraction_pattern <- "\\1_\\2"
}

captured_output <- capture.output({
  print("Cleaned table headers:")
  print(toupper(str_replace_all(colnames(evidence_tbl_filt), tolower(extract_replicate_group), extraction_pattern))%>%
          str_replace_all( "UNIPROT_ACC", "uniprot_acc")) }
  ,type = "message"
)
logdebug(captured_output)

colnames(evidence_tbl_filt) <- str_replace_all(colnames(evidence_tbl_filt), tolower(extract_replicate_group), extraction_pattern) %>%
  toupper( ) %>%
  str_replace_all( "UNIPROT_ACC", "uniprot_acc")


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cleaned_counts<-file.path(args$output_dir, args$output_counts_file)
loginfo("Save the cleaned counts data into %s",cleaned_counts)
vroom::vroom_write(evidence_tbl_filt, cleaned_counts)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
accession_record_file<-file.path(args$output_dir, args$accession_record_file)
loginfo("Save the  data into accession_records to %s",accession_record_file)
vroom::vroom_write(accession_gene_name_tbl_record, accession_record_file)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
te<-toc(quiet = TRUE)
loginfo("%f sec elapsed",te$toc-te$tic)

writeLines(capture.output(sessionInfo()), file.path(args$output_dir,"sessionInfo.txt"))
