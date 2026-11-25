# KIT-247 wording (DRAFTS, VAL-41: the owner sets the final wording)

| Where | Function | Draft |
|---|---|---|
| Confirmation before a save that changes the design table (title "Save the changed design?", buttons Cancel / Save) | `omicsDesignResaveConfirmText` | "This design changes the samples, groups, technical replicates, factors or batch. Saving it sets aside the Quality Control, normalisation and later results made with the previous design; you will apply them again." |
| Notice after the later steps were set aside (warning, shown until closed) | `omicsDesignCutText` | "The design was saved again. The Quality Control, normalisation and later results made with the previous design are set aside: apply Quality Control again, then run the normalisation and the later steps." With the integration input removed: "... The earlier integration input of this omic was removed." |
| Comparison-only re-save | `omicsDesignUnchangedText(TRUE)` | "Comparisons saved. Quality Control and normalisation are kept; run the differential analysis again for the new comparisons." |
| Unchanged re-save | `omicsDesignUnchangedText(FALSE)` | "Design saved; nothing changed." |
| Lipid/metab/trx Reset with nothing current (warning, shown until closed) | `omicNormResetNothingText` | "Normalisation reset: the table the normalisation started from is no longer current (for example, the design was saved again), so nothing was brought back. Apply Quality Control again, then run the normalisation." |
| KIT-222's collision text, amended sentence | `omicsDesignReplicateKeyCollisionText` | "... Re-open the design, mark the runs of each biological sample together, each sample with its own number, save it, then apply Quality Control again and run the normalisation and the later steps. ..." |
