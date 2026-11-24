# KIT-174 wording for the owner's pass (VAL-41)

Implemented and tested as proposed in ticket section 5.1; the owner edits before release. Part (a)'s strings are
listed here; part (b)'s are added when part (b) is implemented.

| Where | Text (exact, as implemented) | Pinned by |
|---|---|---|
| Proteomics normalisation tab, help under "RUV Grouping Variable:" (`R/utils_omic_gui_config.R`, `omicGuiProteomicsNorm`, `ruv_iii$help[["grouping"]]`) | Choose the column that holds your experimental groups (normally 'group'). Samples that share a value in this column are treated as replicates of each other: RUV-III takes what differs between them to be unwanted variation and builds its correction from it. In the Automatic and Manual modes, and in both RUV modes of DIA with limpa, the same column is also used to pick the control proteins (those that differ least between its values), and, except in Manual mode, the number of factors; the sample-correlation filter on this tab compares samples within it. Do not choose 'batch', or any column in which a value contains samples of more than one group: the correction would then be built from the differences between your groups, and the differences between batches would stay in the data. If you know each sample's batch, add it to the Model Formula on the Differential Expression tab instead, for example ~ 0 + group + batch. If the design declares technical replicates, those are the replicates RUV-III uses; this column still picks the control proteins. | T1 |
| Run-time warning (notification of type warning, the log; on lipidomics and metabolomics also the tab's log, prefixed "WARNING:"), `R/mod_omic_norm.R`, `omicNormRuvGroupingWarningText` | The RUV grouping column '{column}' puts samples of different groups together ({n} of its {m} values contain samples of more than one group). RUV-III treats samples with the same value as replicates, so the correction would be built from the differences between your groups. Choose 'group' unless each value of the column lies inside one group. | T10 |
| NEWS, part (a) | The help for the RUV grouping variable in proteomics advised selecting 'batch' for TMT. That advice was wrong: with 'batch' the correction is built from the differences between your groups and the differences between batches stay. To check an earlier analysis, look for 'RUV grouping variable' in its study_parameters.txt (or the normalisation step of its analysis record); if it is not 'group', run the normalisation again with 'group' and, if you know the batches, add batch to the model formula. The normalisation tab now warns when the chosen column puts samples of different groups together (KIT-174). | none |

T2's demonstration (R1a, run 45494891087; seed 174, 2 groups x 2 batches x 3 runs, 300 proteins, k = 1, controls
from the ANOVA at 10 %, 30 controls each), median absolute coefficients of `~ group + batch`, beside the help: before
correction group 0.978, batch 0.689; with 'group' as replicate sets group 0.977, batch 0.0015; with 'batch' as
replicate sets group 0.0013, batch 0.689. The thresholds of the ticket were met on the first run (not reset).

Note for the owner: the help is about 160 words (the ticket's estimate); the owner may want it shorter.
