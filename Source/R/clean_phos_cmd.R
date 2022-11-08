#!/usr/bin/env Rscript

# Author(s): Ignatius Pang, Pablo Galaviz
# Email: cmri-bioinformatics@cmri.org.au
# Children’s Medical Research Institute, finding cures for childhood genetic diseases


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
p_load(ggplot2)
p_load(knitr)
p_load(furrr)
p_load(seqinr)
p_load(seqinr)
p_load(here)
p_load(rlang)
p_load(glue)
p_load(ProteomeRiver)
p_load(optparse)
p_load(tictoc)
p_load(configr)
p_load(logging)
p_load(janitor)

tic()


command_line_options <- commandArgs(trailingOnly = TRUE)
parser <- OptionParser(add_help_option = TRUE)
#Note: options with default values are ignored in the configuration file parsing.
parser <- add_option(parser, c("-d", "--debug"), action = "store_true", default = FALSE,
                     help = "Print debugging output")

parser <- add_option(parser, c("-s", "--silent"), action = "store_true", default = FALSE,
                     help = "Only print critical information to the console.")

parser <- add_option(parser, c("-n", "--no_backup"), action = "store_true", default = FALSE,
                     help = "Deactivate backup of previous run.")


parser <- add_option(parser, c("-c","--config"), type = "character", default = "/home/ubuntu/Workings/2022/Neuropsych_RussellDale_BMP_17_20220530/Source/config_phos.ini", dest = "config",
                     help = "Configuration file.",
                     metavar = "string")

parser <- add_option(parser, c("-o","--output_dir"), type = "character", dest = "output_dir",
                     help = "Directory path for all results files.",
                     metavar = "string")

parser <- add_option(parser, c("-t","--tmp_dir"), type = "character", default = "../Results/cache", dest = "tmp_dir",
                     help = "Directory path for temporary files.",
                     metavar = "string")

parser <- add_option(parser, c("-l","--log_file"), type = "character", default = "output.log", dest = "log_file",
                     help = "Name of the logging file.",
                     metavar = "string")

#Options without a default value have the following priority: configuration file < command line argument
parser <- add_option(parser, "--fasta_file", type="character", dest = "fasta_file",
                     help="Input protein sequence FASTA file with UniProt FASTA header format.",
                     metavar="string")

parser <- add_option(parser, "--fasta_meta_file", type = "character", dest = "fasta_meta_file",
                   help = "R object storage file that records all the sequence and header information in the FASTA file.",
                   metavar = "string")

parser <- add_option(parser, "--raw_counts_file", type = "character", dest = "raw_counts_file",
                     help = "Input file with the protein abundance data.",
                     metavar = "string")

 parser <- add_option(parser, "--site_prob_threshold", type="numeric", dest = "site_prob_threshold",
                     help="Numeric value representing the probability threshold for accepting a primary modification site.",
                     metavar="string")

parser <- add_option(parser, "--recover_site_prob_thresh", type="numeric", dest = "recover_site_prob_thresh",
                     help="Numeric value representing the probability threshold for recovering the abundance data for a lower quality modification site, if the same site has already been found as high probability site.",
                     metavar="string")

parser <- add_option(parser, "--col_pattern_string", type="character", dest = "col_pattern_string",
                     help="String pattern that matches the abundance data columns. It also facilitate the ID of the sample to be extracted.",
                     metavar="string")

parser <- add_option(parser, "--pattern_suffix", type = "character", dest = "pattern_suffix",
                     help = "Regular expression pattern to identify columns with abundance values belonging to the experiment. [default %default]",
                     metavar = "string")

parser <- add_option(parser, "--extract_patt_suffix", type = "character", dest = "extract_patt_suffix",
                     help = "Regular expression pattern to identify columns with abundance values belonging to the experiment. [default %default]",
                     metavar = "string")

parser <- add_option(parser, "--add_cols_string", type="character", dest = "add_cols_string",
                     help="A string listing all the additional columns to be included (e.g. experiment group column). Each column is separated by a comma.",
                     metavar="string")

#parse comand line arguments first.
args <- parse_args(parser)

#parse and merge the configuration file options.
if (args$config != "") {
  args <- config.list.merge(eval.config(file = args$config, config = "clean_phos"), args)
}

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

args <- setArgsDefault(args, "pattern_suffix", as_func=as.character, default_val="_\\d+" )
args <- setArgsDefault(args, "extract_patt_suffix", as_func=as.character, default_val="_(\\d+)" )
args <- setArgsDefault(args, "output_dir", as_func=as.character, default_val="clean_phos" )


testRequiredArguments(args, c(
  "fasta_file"
  ,"raw_counts_file"
  ,"fasta_meta_file"
  ,"site_prob_threshold"
  ,"recover_site_prob_thresh"
  ,"col_pattern_string"
  ,"add_cols_string"

))

testRequiredFiles(c(
  args$fasta_file
  ,args$raw_counts_file))

args<-parseType(args,
  c("site_prob_threshold"
    ,"recover_site_prob_thresh"
  ),as.double)

args<-parseString(args,
  c("add_cols_string"
    ,"col_pattern_string"
    ,"pattern_suffix"
    ,"extract_patt_suffix"
  ))


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

additional_cols <- str_split(  args$add_cols_string, ",")[[1]]
col_pattern <-   janitor::make_clean_names( args$col_pattern_string) # tolower(args$col_pattern_string)  #

captured_output<-capture.output(
  evidence_tbl <- vroom::vroom( args$raw_counts_file)
  ,type = "message"
)
logdebug(captured_output)

evidence_janitor <- janitor::clean_names(evidence_tbl)
colnames(evidence_janitor) <- str_replace(colnames(evidence_janitor), "_i_ds", "_ids" )



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Read fasta file
fasta_meta_file<-file.path(args$tmp_dir,args$fasta_meta_file)
loginfo( "Reading the fasta file: %s",fasta_meta_file)
aa_seq_tbl <-NA
if(!file.exists( file.path( fasta_meta_file))) {
  aa_seq_tbl <- parseFastaFile(args$fasta_file)
  saveRDS( aa_seq_tbl,  file.path( fasta_meta_file))
} else {
  aa_seq_tbl <- readRDS( file.path( fasta_meta_file))
}

## Add the row id column and create a column containing the cleaned  peptide
loginfo("Get row ID and get cleaned peptide sequence.")
evidence_tbl_cleaned <- addColumnsToEvidenceTbl(evidence_janitor )


## Get best accession per entry, work out peptides mapped to multiple genes
loginfo("Use decision tree to get best accession per phosphosite evidence entry")
#TODO:leading_proteins and evidence_id are hard coded, remove or make args
captured_output<-capture.output(
  accession_gene_name_tbl <- chooseBestPhosphositeAccession(evidence_tbl_cleaned,
                                               aa_seq_tbl,
                                               leading_proteins ,
                                               evidence_id)
)
logdebug(captured_output)

## Remove peptides without abundance values at all
loginfo("Remove peptides without abundance values at all")
evidence_tbl_filt <- removePeptidesWithoutAbundances(evidence_tbl_cleaned, col_pattern)

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## For all the multi-phosphosites peptide extract their intensity, filter peptide with no intensity across all samples, extract site probabilities
loginfo("Filter peptides with no intensity across all samples, extract intensity data, extract sites")
sites_probability_tbl <- filterPeptideAndExtractProbabilities (evidence_tbl_filt ,
                                                               accession_gene_name_tbl,
                                                               col_pattern,
                                                               accession_col = leading_proteins,
                                                               phospho_site_prob_col = phospho_sty_probabilities,
                                                               num_phospho_site_col = phospho_sty )

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Get the peptide start and end position for each peptide
loginfo("Add peptide start and end position")
peptide_start_and_end <- addPeptideStartAndEnd(sites_probability_tbl , aa_seq_tbl )


## Get the phosphosites position string
loginfo("Add string listing the positions of phosphosites")
phosphosite_pos_string_tbl <- addPhosphositesPositionsString(peptide_start_and_end )

## Get the string listing all the 15-mer sequences, each sequence has the phosphorylation site in the middle
loginfo("Add string listing all 15-mer sequences, each sequence has phosphosite in the center")
get_15_mer_tbl <- addXMerStrings(phosphosite_pos_string_tbl, padding_length=7)

## Get peptides with at least one phosphosite over threshold. Find all peptides with same sites as another peptides that contained at least one phosphosies >= threshold.
loginfo("Get high conf. peptides (e.g. phosphosites >= threshold). Get peptide W/ same sites as high conf. peptide.")
get_15_mer_tbl_filt <- filterByScoreAndGetSimilarPeptides(get_15_mer_tbl,
                                                            args$site_prob_threshold,
                                                            args$recover_site_prob_thresh)

## Pivot the phosphosites to a longer table
loginfo("Pivot phosphosites/phosphopeptide table to long format")
all_phos_sites_long_tbl <- allPhosphositesPivotLonger(get_15_mer_tbl_filt,
                                                      additional_cols ,
                                                      col_pattern,
                                                      pattern_suffix = args$pattern_suffix,
                                                      extract_patt_suffix=args$extract_patt_suffix)

## Group peptides from paralog proteins
loginfo("Group peptides from paralog proteins ")
paralog_sites_long <- groupParalogPeptides (all_phos_sites_long_tbl, additional_cols )

## Pivot the phosphosites data to a wide format
loginfo("Pivot phosphosites/phosphopeptide table to wide format")
all_phos_sites_wide_tbl <- allPhosphositesPivotWider(paralog_sites_long,
                                                     additional_cols )

## Summarise the abundance values for each unique phosphosites (mean, median, sum), return table in long format
loginfo("Summarise abundance values for each unique phosphosites, long format")
summarised_long_tbl_list <- uniquePhosphositesSummariseLongList(paralog_sites_long,
                                                                additional_cols )

## Summarise the abundance values for each unique phosphosites (mean, median, sum), return table in wide format
loginfo("Summarise abundance values for each unique phosphosites, wide format")
summarised_wide_tbl_list <- uniquePhosphositesSummariseWideList(summarised_long_tbl_list,
                                                                additional_cols)

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
results_list  <- list( summarised_wide_list = summarised_wide_tbl_list,
                summarised_long_list = summarised_long_tbl_list,
                all_phos_sites_wide  = all_phos_sites_wide_tbl,
                all_phos_sites_long  = all_phos_sites_long_tbl )

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# output_files <- map_chr( list(  "mean", "median", "sum"),
#                          ~file.path(args$output_dir, paste0(., "_phosphosites.tsv" )  ) )
#
# walk2(  results_list$summarised_wide_list,
#        output_files,
#        ~vroom::vroom_write( .x, .y))


output_files <- map_chr( list(  "mean", "median", "sum"),
                         ~file.path(args$output_dir, paste0(., "_phoshpsites.tsv" )  ) )

walk2(  results_list$summarised_wide_list,
        output_files,
        ~vroom::vroom_write( .x, .y))

vroom::vroom_write( results_list$all_phos_sites_wide,
         file.path(args$output_dir, "all_phos_sites_wide_tbl.tsv" ))

saveRDS( results_list$summarised_long_list,
         file.path(args$output_dir, "summarised_long_tbl_list.RDS" ))

saveRDS( results_list$summarised_wide_list,
         file.path(args$output_dir, "summarised_wide_tbl_list.RDS" ))

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
sample_names_file<-file.path(args$output_dir, "sample_names.tab")
loginfo("Save the sample names to %s", "sample_names.tab")
sample_names <- colnames(results_list$summarised_wide_list[[1]])[c(-1,-2)]
vroom::vroom_write(data.frame( sample_names=t(t(sample_names) ) ), sample_names_file)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
te<-toc(quiet = TRUE)
loginfo("%f sec elapsed",te$toc-te$tic)

writeLines(capture.output(sessionInfo()), file.path(args$output_dir,"sessionInfo.txt"))

