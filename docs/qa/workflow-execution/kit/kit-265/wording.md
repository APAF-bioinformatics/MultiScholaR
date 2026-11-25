# KIT-265: user-facing texts (drafts; the owner sets the final wording, VAL-41)

All texts below are as implemented. The code constants and builders are in `R/mod_prot_import_organism_helpers.R` and
`R/mod_prot_import_orchestration_helpers.R`. Messages run 40 to 50 words (ticket section 3).

## In the application

| Id | Where | Text |
|---|---|---|
| N1 | Warning notice, shown until closed, under the fixed id `prot_import_organism_cancelled`; route (i): no organism with a taxonomy ID could be matched (FASTA headers without `OS=`/`OX=`); the import continues. The form's organism and taxon are filled in | No organism with a taxonomy ID could be matched to the proteins in the search results (the FASTA headers name no organism or taxonomy ID for them, or the proteins are not in the FASTA file), so no organism window is shown. This analysis will use the organism entered in the form: `<name>` (taxon `<id>`). To use another organism, untick 'Mixed species FASTA', enter it and press 'Process Imported Data' again. |
| N1b | The second sentence of N1 when the form's Taxonomy ID is empty (the taxon alone sets the UniProt query and the enrichment organism; AUDIT-KIT-265-impl I01, correction pass; enrichment's behaviour with no taxon is not checked, so the sentence does not mention it). For the wording pass: with a blank form the third sentence ("To use another organism") reads oddly | The form has no Taxonomy ID, so UniProt annotation will have no organism to use. |
| N1c | The second sentence of N1 when the form's organism name is empty and its Taxonomy ID is not (I01) | This analysis will use the Taxonomy ID entered in the form: `<id>` (the form gives no organism name). |
| R-ii | Error (shown as "Error: ..." for 10 s, the processing window closed): counting the organisms failed. `<text>` is the error's text | The import was stopped: the organisms in the FASTA file could not be counted (`<text>`). The organism is chosen from that count, so nothing was imported. Check the FASTA file, or untick 'Mixed species FASTA', enter the organism yourself and press 'Process Imported Data' again. |
| R-iii | Error: the search results have no protein column. `<protein_col>` is the column's name; when the import has no name for it the text reads "no protein column" | The import was stopped: the search results have no '`<protein_col>`' column, so their proteins cannot be matched to the organisms in the FASTA file. Nothing was imported. Export the search results with that column and import them again. |
| R-iv | Error: FASTA processing failed (neutral: the same handler also catches a missing `seqinr` and failed writes of the cache copies; its reason goes only to the log) | The import was stopped: the FASTA file could not be read or processed, so the organisms in it cannot be counted. Nothing was imported. Check that the file is a FASTA file of protein sequences, or untick 'Mixed species FASTA', enter the organism yourself and press 'Process Imported Data' again. |
| R-c | Error: route (i)(c), the FASTA names organisms and none of the data's proteins is in it (owner's decision VAL-304: refuse). "Could be found", because a FASTA whose accessions are written differently from the search engine's is also reported as not found | The import was stopped: none of the proteins in the search results could be found in the FASTA file, so their organisms cannot be counted. Nothing was imported. Check that this is the FASTA file the search used, or untick 'Mixed species FASTA', enter the organism yourself and press 'Process Imported Data' again. |
| W-D2 | The organism window, a help text directly after the radio buttons (VAL-291 D-2). It deliberately does not say that the other organisms' proteins "get no gene names": whether the differential-analysis results fill `gene_names` for them from the identifier table (`R/func_prot_da_results_methods.R:69-79`) was not traced | Only the organism you choose here is looked up in UniProt and tested for enrichment; the other organisms' proteins are not, and they stay in the normalisation, RUV and the differential analysis (and in its correction for multiple testing) unless you tick 'Filter data to keep only proteins from selected organism'. |

## ALTERNATIVE draft for R-c (AUDIT-KIT-265-impl I04, option (c), recommended by the auditor; NOT implemented)

The shipped R-c is unchanged; the owner decides. The register's draft (VAL-41), conditional on why the proteins were not
found (it also covers I03: identifiers written differently from the FASTA's):

> The import was stopped: none of the proteins in the search results could be found in the FASTA file. Check that this
> is the FASTA file the search used. If it is, the search results may name the proteins differently from the FASTA
> file: untick 'Mixed species FASTA', enter the organism yourself and press 'Process Imported Data' again.

(48 words. Option (a) keeps the shipped text; option (b) drops the untick clause: "Supply the FASTA file the search
used and import again.")

## README

`README.md`, Proteomics row, "Key Features" cell. Before: "Peptide rollup (IQ), Auto-RUV-IIIc, Mixed Species support".
After (draft): "Peptide rollup (IQ), Auto-RUV-IIIc, Mixed-species FASTA files (one organism annotated and enriched)".

## Owed elsewhere (not written here)

- The user guide paragraph: `docs/qa/workflow-execution/kit/kit-250/wording.md` (O3), to be placed when the guide is
  written.
- The paper's sentence (the manuscripts repository, in the owner's voice): one organism per project is annotated and
  tested for enrichment; the other organisms' proteins remain in normalisation and testing unless filtered out at import.

## Left to the wording pass (pre-existing, not changed)

The window's heading "Multiple organisms detected" shown for one organism, and the unmatched row's percentage
denominator (AUDIT-KIT-250-impl I06).
