# KIT-250: user-facing texts (drafts; the owner sets the final wording, VAL-41)

## In the application (as implemented)

| Id | Where | Text |
|---|---|---|
| W1 | Warning after Cancel in the "Select Primary Organism" window, under the fixed id `prot_import_organism_cancelled` (one copy; it stays until closed or until the next import starts, which removes it); the import page's status panel after Cancel | Import not finished: no organism was chosen. Press 'Process Imported Data' again and choose the organism in the window that follows, or untick 'Mixed species FASTA' and enter the organism yourself. |
| W2 | Import page's status panel while the window is open | Choose the organism in the window to finish the import. |
| W3 | Data Summary's organism line while the window is open | Organism: not confirmed |
| — | Unchanged | "Primary organism set to: <name> (Taxon: <id>)", "Data import successful!" |

Constants: `PROT_IMPORT_ORGANISM_CANCELLED_TEXT` (W1) and `PROT_IMPORT_ORGANISM_WAITING_TEXT` (W2), `R/mod_prot_import_ui_helpers.R`.

## O3: the single-organism limit for mixed-species data (the owner's decision (a), VAL-240)

**Open question for the owner (AUDIT-KIT-250-impl I06, D-2): where should this be stated?** The repository has no user guide: `docs/manual/` holds only the design-builder how-to. `README.md:95` advertises "Mixed Species support" with no limit, and the window itself says only "Please select the primary organism for this analysis."

The audit recommends, in place and wording the owner's (VAL-41):
- one sentence in the organism window, under the radio buttons, where the choice and the filter are made;
- README's cell qualified ("Mixed-species FASTA: one organism annotated and enriched");
- the full paragraph below in the user guide when it is written;
- the paper, if it mentions mixed-species support.

Nothing is placed by KIT-250: O3 is not met until the owner chooses the place.

Draft (corrected per AUDIT-KIT-250-impl I06):

> **Mixed-species FASTA files.** When 'Mixed species FASTA' is ticked, the import counts, for each organism named in the FASTA file's headers, how many of your data's protein accessions belong to it (the accessions of a protein group are counted separately), and you choose one organism. Protein annotation (gene names from UniProt) and enrichment analysis use that organism only. Proteins of the other organisms get no gene names and are not tested for enrichment. Unless you tick 'Filter data to keep only proteins from selected organism' in the same window, those proteins stay in the data: they take part in normalisation, in RUV's choice of control proteins, and in the differential analysis, and they count in its correction for multiple testing. For data from two organisms studied together (for example host and pathogen, or a xenograft), decide in that window whether the other organism's proteins stay in the analysis: keeping them changes the normalisation, the control proteins RUV chooses and the number of tests the correction for multiple testing counts; removing them analyses one organism only.

Basis:
- `R/func_prot_import_organisms.R:130-199`: the per-organism count of matched accessions.
- `R/func_prot_import_organisms.R:202-224`: one organism chosen, the most frequent preselected.
- `R/func_prot_annotation_uniprot_helpers.R:254`: one `organism_id` in the UniProt query.
- `R/mod_prot_import_organism_helpers.R:31-35`: the filter is off by default and is chosen in the window.

KIT-250 changes none of these.

Two pre-existing inaccuracies in the window, for the same wording pass (I06):
- its heading says "Multiple organisms detected" when only one organism is found (G0: Homo sapiens, 400, 100 %);
- the unmatched row's percentage uses a different denominator from the organisms' rows (`R/func_prot_import_organisms.R:178-191`).
