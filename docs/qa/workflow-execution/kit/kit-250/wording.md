# KIT-250: user-facing texts (drafts; the owner sets the final wording, VAL-41)

## In the application (as implemented)

| Id | Where | Text |
|---|---|---|
| W1 | Persistent warning after Cancel in the "Select Primary Organism" window; the import page's status panel after Cancel | Import not finished: no organism was chosen. Press 'Process Imported Data' again and choose the organism in the window that follows, or untick 'Mixed species FASTA' and enter the organism yourself. |
| W2 | Import page's status panel while the window is open | Choose the organism in the window to finish the import. |
| W3 | Data Summary's organism line while the window is open | Organism: not confirmed |
| — | Unchanged | "Primary organism set to: <name> (Taxon: <id>)", "Data import successful!" |

Constants: `PROT_IMPORT_ORGANISM_CANCELLED_TEXT` (W1) and `PROT_IMPORT_ORGANISM_WAITING_TEXT` (W2), `R/mod_prot_import_ui_helpers.R`.

## O3: the single-organism limit for mixed-species data (the owner's decision (a), VAL-240)

To be stated in the user guide. If the paper mentions support for mixed-species FASTA files, it should be stated there too, in the owner's wording. Draft:

> **Mixed-species FASTA files.** When 'Mixed species FASTA' is ticked, the import counts the proteins of each organism in the FASTA file, and you choose one of them. Protein annotation (gene names from UniProt) and enrichment analysis use that organism only. Proteins of the other organisms get no gene names and are not tested for enrichment. Unless you tick 'Filter data to keep only proteins from selected organism' in the same window, those proteins stay in the data: they take part in normalisation, in RUV's choice of control proteins, and in the differential analysis, and they count in its correction for multiple testing. For data from two organisms studied together (for example host and pathogen, or a xenograft), decide before the import whether the other organism's proteins should stay in the analysis.

Basis: `R/func_prot_import_organisms.R:202-224` (one organism chosen, the most frequent preselected); `R/func_prot_annotation_uniprot_helpers.R:254` (one `organism_id` in the UniProt query); `R/mod_prot_import_organism_helpers.R:31-35` (the filter is off by default). KIT-250 does not change any of these.
