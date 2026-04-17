# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# ============================================================================
# func_prot_annotation.R
# ============================================================================
# Purpose: Annotation and protein ID management functions
# 
# This file contains functions for protein annotation, UniProt data retrieval,
# protein ID cleaning and selection, FASTA file processing, and phosphoproteomics
# annotation. Functions in this file are used across proteomics workflows.
#
# Functions to extract here:
# - UniProt annotation functions (getUniprotAnnotations, etc.)
# - Protein ID cleaning functions (.cleanProteinIds, cleanIsoformNumber, etc.)
# - Best accession selection functions (chooseBestProteinAccession, etc.)
# - FASTA processing functions (parseFastaFile, processFastaFile, etc.)
# - Phosphoproteomics annotation functions (addPhosphositesPositionsString, etc.)
# - ID conversion functions (convertEnsemblToUniprot, etc.)
# - Additional annotation helper functions
#
# Dependencies:
# - UniProt.ws, seqinr, GO.db
# - func_general_helpers.R (for utility functions)
# ============================================================================

# TODO: Extract the following functions from their current locations:

# === UniProt Annotation Functions ===

# Function 1: getUniprotAnnotations()
# Current location: R/annotation.R
# Description: Gets UniProt annotations
# getUniprotAnnotations <- function(...) {
#   # Extract from R/annotation.R
# }

# Function 2: getUniprotAnnotationsFull()
# Current location: R/annotation.R
# Description: Gets full UniProt annotations
# getUniprotAnnotationsFull <- function(...) {
#   # Extract from R/annotation.R
# }

# Function 3: getUniProtAnnotation()
# Current location: R/annotation.R
# Description: Gets UniProt annotation for input table
# getUniProtAnnotation <- function(...) {
#   # Extract from R/annotation.R
# }

# Function 4: download_uniprot_data()
# Current location: R/annotation.R
# Description: Downloads UniProt data
# download_uniprot_data <- function(...) {
#   # Extract from R/annotation.R
# }

# Function 5: directUniprotDownload()
# Current location: R/annotation.R
# Description: Direct UniProt download
# directUniprotDownload <- function(...) {
#   # Extract from R/annotation.R
# }

# Function 6: buildAnnotationIdToAnnotationNameDictionary()
# Current location: R/enrichment_functions.R
# Description: Builds annotation ID to name dictionary
# buildAnnotationIdToAnnotationNameDictionary <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 7: buildOneProteinToAnnotationList()
# Current location: R/enrichment_functions.R
# Description: Builds one protein to annotation list
# buildOneProteinToAnnotationList <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 8: convertIdToAnnotation()
# Current location: R/enrichment_functions.R
# Description: Converts ID to annotation
# convertIdToAnnotation <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 9: convertKeyToAttribute()
# Current location: R/helper_functions.R
# Description: Converts key to attribute (also in helpers)
# convertKeyToAttribute <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 10: createIdToAttributeHash()
# Current location: R/helper_functions.R
# Description: Creates ID to attribute hash (also in helpers)
# createIdToAttributeHash <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 11: convertProteinAccToGeneSymbol()
# Current location: R/enrichment_functions.R
# Description: Converts protein accession to gene symbol
# convertProteinAccToGeneSymbol <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 12: getUniprotAccToGeneSymbolDictionary()
# Current location: R/enrichment_functions.R
# Description: Gets UniProt accession to gene symbol dictionary
# getUniprotAccToGeneSymbolDictionary <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 13: convertEnsemblToUniprot()
# Current location: R/annotation.R
# Description: Converts Ensembl IDs to UniProt
# convertEnsemblToUniprot <- function(...) {
#   # Extract from R/annotation.R
# }

# Function 14: detectEnsemblIds()
# Current location: R/annotation.R
# Description: Detects Ensembl IDs
# detectEnsemblIds <- function(...) {
#   # Extract from R/annotation.R
# }

# Function 15: matchAnnotations()
# Current location: R/annotation.R
# Description: Matches annotations
# matchAnnotations <- function(...) {
#   # Extract from R/annotation.R
# }

# Function 16: convertDpcDAToStandardFormat()
# Current location: R/limpa_functions.R
# Description: Converts DPC DE to standard format
# convertDpcDAToStandardFormat <- function(...) {
#   # Extract from R/limpa_functions.R
# }

# === FASTA Processing Functions ===

# Function 17: getFastaFields()
# Current location: R/get_best_accession_helper.R
# Description: Gets FASTA fields
# getFastaFields <- function(...) {
#   # Extract from R/get_best_accession_helper.R
# }

# Function 18: parseFastaFile()
# Current location: R/get_best_accession_helper.R
# Description: Parses FASTA file
# parseFastaFile <- function(...) {
#   # Extract from R/get_best_accession_helper.R
# }

# Function 19: parseFastaObject()
# Current location: R/get_best_accession_helper.R
# Description: Parses FASTA object
# parseFastaObject <- function(...) {
#   # Extract from R/get_best_accession_helper.R
# }

# Function 20: processFastaFile()
# Current location: R/get_best_accession_helper.R
# Description: Processes FASTA file
# processFastaFile <- function(...) {
#   # Extract from R/get_best_accession_helper.R
# }

# === Protein ID Selection Functions ===

# Function 21: chooseBestProteinAccession()
# Current location: R/proteinVsSamplesS4Objects.R, R/get_best_accession_helper.R
# Type: S4 method (exportMethods) and helper
# Description: Chooses best protein accession
# setMethod(f = "chooseBestProteinAccession", ...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 22: chooseBestProteinAccessionHelper()
# Current location: R/get_best_accession_helper.R
# Description: Helper for choosing best protein accession
# chooseBestProteinAccessionHelper <- function(...) {
#   # Extract from R/get_best_accession_helper.R
# }

# Function 23: chooseBestProteinAccession_s3()
# Current location: R/helper_functions.R
# Description: S3 version of choose best protein accession
# chooseBestProteinAccession_s3 <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 24: rankProteinAccessionHelper()
# Current location: R/get_best_accession_helper.R
# Description: Ranks protein accessions
# rankProteinAccessionHelper <- function(...) {
#   # Extract from R/get_best_accession_helper.R
# }

# === Protein ID Cleaning Functions ===

# Function 25: .cleanProteinIds()
# Current location: R/get_best_accession_helper.R
# Description: Cleans protein IDs (internal function)
# .cleanProteinIds <- function(...) {
#   # Extract from R/get_best_accession_helper.R
# }

# Function 26: cleanIsoformNumber()
# Current location: R/get_best_accession_helper.R
# Description: Cleans isoform number from accession
# cleanIsoformNumber <- function(...) {
#   # Extract from R/get_best_accession_helper.R
# }

# Function 27: cleanMaxQuantProteins()
# Current location: R/get_best_accession_helper.R
# Description: Cleans MaxQuant protein IDs
# cleanMaxQuantProteins <- function(...) {
#   # Extract from R/get_best_accession_helper.R
# }

# Function 28: updateProteinIDs()
# Current location: R/get_best_accession_helper.R
# Description: Updates protein IDs
# updateProteinIDs <- function(...) {
#   # Extract from R/get_best_accession_helper.R
# }

# === Phosphoproteomics Annotation Functions ===

# Function 29: getBestPosition()
# Current location: R/phosphoproteomics_helper.R
# Description: Gets best position for phosphosite
# getBestPosition <- function(...) {
#   # Extract from R/phosphoproteomics_helper.R
# }

# Function 30: getBestPositionFutureMap()
# Current location: R/phosphoproteomics_helper.R
# Description: Gets best position using future map
# getBestPositionFutureMap <- function(...) {
#   # Extract from R/phosphoproteomics_helper.R
# }

# Function 31: getMaxProb()
# Current location: R/phosphoproteomics_helper.R
# Description: Gets maximum probability
# getMaxProb <- function(...) {
#   # Extract from R/phosphoproteomics_helper.R
# }

# Function 32: getMaxProbFutureMap()
# Current location: R/phosphoproteomics_helper.R
# Description: Gets max probability using future map
# getMaxProbFutureMap <- function(...) {
#   # Extract from R/phosphoproteomics_helper.R
# }

# Function 33: getPosString()
# Current location: R/phosphoproteomics_helper.R
# Description: Gets position string
# getPosString <- function(...) {
#   # Extract from R/phosphoproteomics_helper.R
# }

# Function 34: getXMerString()
# Current location: R/phosphoproteomics_helper.R
# Description: Gets X-mer string
# getXMerString <- function(...) {
#   # Extract from R/phosphoproteomics_helper.R
# }

# Function 35: getXMersList()
# Current location: R/phosphoproteomics_helper.R
# Description: Gets list of X-mers
# getXMersList <- function(...) {
#   # Extract from R/phosphoproteomics_helper.R
# }

# Function 36: addXMerStrings()
# Current location: R/phosphoproteomics_helper.R
# Description: Adds X-mer strings
# addXMerStrings <- function(...) {
#   # Extract from R/phosphoproteomics_helper.R
# }

# Function 37: addPeptideStartAndEnd()
# Current location: R/phosphoproteomics_helper.R
# Description: Adds peptide start and end positions
# addPeptideStartAndEnd <- function(...) {
#   # Extract from R/phosphoproteomics_helper.R
# }

# Function 38: addPhosphositesPositionsString()
# Current location: R/phosphoproteomics_helper.R
# Description: Adds phosphosite position strings
# addPhosphositesPositionsString <- function(...) {
#   # Extract from R/phosphoproteomics_helper.R
# }

# Function 39: chooseBestPhosphositeAccession()
# Current location: R/get_best_accession_helper.R
# Description: Chooses best phosphosite accession
# chooseBestPhosphositeAccession <- function(...) {
#   # Extract from R/get_best_accession_helper.R
# }

# Function 40: getUniprotAccRankFromSitesId()
# Current location: R/phosphoproteomics_helper.R
# Description: Gets UniProt accession rank from sites ID
# getUniprotAccRankFromSitesId <- function(...) {
#   # Extract from R/phosphoproteomics_helper.R
# }

# Function 41: formatPhosphositePosition()
# Current location: R/phosphoproteomics_helper.R
# Description: Formats phosphosite position
# formatPhosphositePosition <- function(...) {
#   # Extract from R/phosphoproteomics_helper.R
# }

# Function 42: allPhosphositesPivotLonger()
# Current location: R/phosphoproteomics_helper.R
# Description: Pivots all phosphosites to long format
# allPhosphositesPivotLonger <- function(...) {
#   # Extract from R/phosphoproteomics_helper.R
# }

# Function 43: allPhosphositesPivotWider()
# Current location: R/phosphoproteomics_helper.R
# Description: Pivots all phosphosites to wide format
# allPhosphositesPivotWider <- function(...) {
#   # Extract from R/phosphoproteomics_helper.R
# }

# Function 44: uniquePhosphositesSummariseLongList()
# Current location: R/phosphoproteomics_helper.R
# Description: Summarizes unique phosphosites in long list
# uniquePhosphositesSummariseLongList <- function(...) {
#   # Extract from R/phosphoproteomics_helper.R
# }

# Function 45: uniquePhosphositesSummariseWideList()
# Current location: R/phosphoproteomics_helper.R
# Description: Summarizes unique phosphosites in wide list
# uniquePhosphositesSummariseWideList <- function(...) {
#   # Extract from R/phosphoproteomics_helper.R
# }

# Function 46: processMultisiteEvidence()
# Current location: R/phosphoproteomics_helper.R
# Description: Processes multisite evidence
# processMultisiteEvidence <- function(...) {
#   # Extract from R/phosphoproteomics_helper.R
# }

# Function 47: addColumnsToEvidenceTbl()
# Current location: R/phosphoproteomics_helper.R
# Description: Adds columns to evidence table
# addColumnsToEvidenceTbl <- function(...) {
#   # Extract from R/phosphoproteomics_helper.R
# }

# Function 48: batchQueryEvidence()
# Current location: R/annotation.R
# Description: Batch queries UniProt evidence
# batchQueryEvidence <- function(...) {
#   # Extract from R/annotation.R
# }

# Function 49: batchQueryEvidenceGeneId()
# Current location: R/annotation.R
# Description: Batch queries evidence with gene ID
# batchQueryEvidenceGeneId <- function(...) {
#   # Extract from R/annotation.R
# }
































































