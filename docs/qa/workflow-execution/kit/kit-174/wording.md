# KIT-174 wording for the owner's pass (VAL-41)

Implemented and tested as proposed in ticket section 5.1; the owner edits before release. Part (a)'s strings are
listed here; part (b)'s are added when part (b) is implemented.

| Where | Text (exact, as implemented) | Pinned by |
|---|---|---|
| Proteomics normalisation tab, help under "RUV Grouping Variable:" (`R/utils_omic_gui_config.R`, `omicGuiProteomicsNorm`, `ruv_iii$help[["grouping"]]`) | Choose the column that holds your experimental groups (normally 'group'). Samples that share a value in this column are treated as replicates of each other: RUV-III takes what differs between them to be unwanted variation and builds its correction from it. In the Automatic and Manual modes, and in both RUV modes of DIA with limpa, the same column is also used to pick the control proteins (those that differ least between its values), and, except in Manual mode, the number of factors; the sample-correlation filter on this tab compares samples within it. Do not choose 'batch', or any column in which a value contains samples of more than one group: the correction would then be built from the differences between your groups, and the differences between batches would stay in the data. If you know each sample's batch, add it to the Model Formula on the Differential Expression tab instead, for example ~ 0 + group + batch. If the design declares technical replicates, those are the replicates RUV-III uses; this column still picks the control proteins. | T1 |
| Run-time warning (notification of type warning, the log; on lipidomics and metabolomics also the tab's log, prefixed "WARNING:"), `R/mod_omic_norm.R`, `omicNormRuvGroupingWarningText` | The RUV grouping column '{column}' puts samples of different groups together ({n} of its {m} values contain samples of more than one group). RUV-III treats samples with the same value as replicates, so the correction would be built from the differences between your groups. Choose 'group' unless each value of the column lies inside one group. | T10 |
| NEWS, part (a) | The help for the RUV grouping variable in proteomics advised selecting 'batch' for TMT. That advice was wrong: with 'batch' the correction is built from the differences between your groups and the differences between batches stay. To check an earlier analysis, look for 'RUV grouping variable' in its study_parameters.txt (or the normalisation step of its analysis record); if it is not 'group', run the normalisation again with 'group' and, if you know the batches, add batch to the model formula. The normalisation tab now warns when the chosen column puts samples of different groups together (KIT-174). | none |

T2's demonstration (R1a, run 45494891087; the same figures on head in R2a, run 45496986018): one simulated data set (seed 174) of 2 groups x 2 batches x 3 runs and 300 proteins, complete data, noise SD 0.2 log2; a batch effect whose size per protein is drawn from N(0, 1); a group difference of 1 log2 on every protein, with random sign, so that whichever proteins are chosen as controls carry it; RUV-III with k fixed at 1 and 30 controls (10 %) chosen by the ANOVA on the same column. Median over the 300 proteins of the absolute coefficient in a per-protein fit of `~ group + batch`: before correction group 0.978, batch 0.689; with 'group' as replicate sets group 0.977, batch 0.0015; with 'batch' as replicate sets group 0.0013, batch 0.689. This shows both halves of the help's sentence on a design built to show them. It does not measure how much of a real study's group difference would be lost: that depends on how many control proteins carry group differences, on their size against the other variation within batches, and on k. The four thresholds were fixed in the ticket (revision 1.1, 2026-11-23) before the run and were met on the first run.

Note for the owner: the help is about 160 words (the ticket's estimate); the owner may want it shorter.

## Open for the owner (implementation audit AUDIT-KIT-174a-impl-20261125, section 6; not applied)

- **O1 ("built from the differences between your groups").** More exactly, the correction is built from the variation
  within each value of the column, which includes the group differences, and it removes them to the extent the
  controls carry them. The auditor's suggested texts, if the owner wants the precise form:
  - Help: replace "the correction would then be built from the differences between your groups, and the differences
    between batches would stay in the data." with "the correction would then be built partly from the differences
    between your groups, so it can remove them, and the differences between that column's values (for example between
    batches) would stay in the data."
  - Warning: replace "so the correction would be built from the differences between your groups." with "so the
    correction would be built partly from the differences between your groups and can remove them."
  - NEWS: replace "with 'batch' the correction is built from the differences between your groups and the differences
    between batches stay" with "with 'batch' the correction is built partly from the differences between your groups,
    so it can remove them, and the differences between batches stay".
- **O2 (lipidomics and metabolomics users).** The NEWS bullet tells only proteomics users how to check earlier
  analyses; lipidomics and metabolomics offered the same choice and now warn too. Suggested sentence to append: "The
  same applies to lipidomics and metabolomics: their normalisation tab offered the same choice (with a neutral help
  text) and now warns too; to check an earlier analysis, look for 'RUV grouping variable' in the normalisation step of
  its analysis record."

If either is taken, T1's and T10's strings and this list change in the same commit.
