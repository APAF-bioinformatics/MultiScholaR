# KIT-247 wording (DRAFTS, VAL-41: the owner sets the final wording)

| Where | Function | Draft |
|---|---|---|
| Confirmation before a save that changes the design table (title "Save the changed design?", buttons Cancel / Save) | `omicsDesignResaveConfirmText` | "This design changes the samples, groups, technical replicates, factors or batch. Saving it sets aside the Quality Control, normalisation and later results made with the previous design; you will apply them again." |
| Notice after the later steps were set aside (warning, shown until closed) | `omicsDesignCutText` | "The design was saved again. The Quality Control, normalisation and later results made with the previous design are set aside: apply Quality Control again, then run the normalisation and the later steps." With the integration input removed: "... The earlier integration input of this omic was removed." |
| Comparison-only re-save | `omicsDesignUnchangedText(TRUE)` | "Comparisons saved. Quality Control and normalisation are kept; run the differential analysis again for the new comparisons." |
| Unchanged re-save | `omicsDesignUnchangedText(FALSE)` | "Design saved; nothing changed." |
| Lipid/metab/trx Reset with nothing current (warning, shown until closed) | `omicNormResetNothingText` | "Normalisation reset: the table the normalisation started from is no longer current (for example, the design was saved again), so nothing was brought back. Apply Quality Control again, then run the normalisation." |
| KIT-222's collision text, amended sentence | `omicsDesignReplicateKeyCollisionText` | "... Re-open the design, mark the runs of each biological sample together, each sample with its own number, save it, then apply Quality Control again and run the normalisation and the later steps. ..." |

## From the implementation audit (AUDIT-KIT-247-impl-20261126 I09; the owner sets the wording, VAL-41)

- (a) The question says "This design changes the samples, groups, technical replicates, factors or batch" but a change
  of the data table alone also asks (e.g. after the data were imported again); it does not say that the omic's
  integration input will be removed (only the notice after does); "you will apply them again" reads as a promise ("you
  will need to apply Quality Control and run the later steps again" is the instruction).
- (b) "Design saved; nothing changed." is shown when something did change but is not compared: a formula-only change
  (written to the configuration and `config.ini`; the DA is not set to pending; the DA tabs take their formula from
  their own input or the analysed table, not from the builder); a proteomics Import Existing Design whose organism,
  taxonomy id or FASTA differ (the unchanged path applies none of them); a transcriptomics change of the group level
  order only (`results$group_levels`, not compared; no number moves with explicit comparisons, INFERRED).
- (c) The builder's own toast "Design saved successfully!" / "Design saved successfully. You can close this builder."
  shows on every Save press, before the question; after Cancel the user has been told "saved".
- (d) "integration input of this omic" is internal vocabulary (a plain form: "the file this omic had saved for the
  multi-omics integration was removed; press 'Save Workflow Arguments' again to write it").

Owner's questions:
1. Should a formula-only change set the differential analysis to pending (as a change of comparisons does)?
2. Should an organism, taxonomy or FASTA change on Import Existing Design count as a change?

Answered: VAL-288, then withdrawn from this ticket by VAL-331 and decided for KIT-277 by VAL-345.

## After the withdrawal of VAL-288 (VAL-331)

- A formula-only Save, and on proteomics Import Existing Design an organism, taxonomy id or FASTA-only Save, show
  "Design saved; nothing changed." (`omicsDesignUnchangedText(FALSE)`): no question, nothing reset. The limited question
  and its notices (the closure check's C04) do not exist on this branch; VAL-288 is its own ticket, where its wording is
  set.
- A confirmed Import Existing Design that is refused shows "Error during import: <reason>" (VAL-277; I01).
- The sentence is not literally true in these cases (the formula is written to the configuration; the import dialog's
  organism and the folder's FASTA table are ignored): KIT-277's wording. That holds for a formula-only, organism-only or
  FASTA-only Save until KIT-277.
- I09 (a), (c), (d) stay with the owner.
