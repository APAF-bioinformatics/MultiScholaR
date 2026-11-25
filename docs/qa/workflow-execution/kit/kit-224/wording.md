# KIT-224 user-facing texts (draft; the owner sets the final wording)

No message, notice or screen text is added or changed. The only new user-facing text is the NEWS entry (ticket
section 12, VAL-41), under `### Scientific changes` of the unreleased heading:

> Lipidomics and metabolomics RUV-III correction (Automatic and Manual, the default): each feature now keeps its mean
> abundance after the correction, as in proteomics. Before, every feature's corrected values were shifted by a
> constant that depended on the control features' levels and on the unit the intensities were expressed in. Fold
> changes are unchanged where the sample-correlation filter keeps the same samples. The moderated p- and q-values of
> the differential analysis can change, because limma's variance trend is fitted against each feature's mean
> abundance, and the sample-correlation filter can keep different samples. The abundances in the exported and
> displayed tables, the sample-correlation plot, the session file and the RUV-corrected matrices offered to the
> integration change; PCA, RLE and the PC boxplots do not. Analyses recorded before this change replay as they were
> made. Transcriptomics: the RUVIII export file and the RUV-corrected integration input keep each gene's logCPM
> level; the differential analysis is not affected (KIT-224).

In the older entry (NEWS.md, KIT-127's lipidomics and metabolomics RUV entry, a historical section), the true sentence is
kept and a pointer added (implementation audit N1):

> Without `ruv_factors` every result is as before: `ruvIII_C_Varying()` for lipid and metabolite data is not changed and
> still shifts each feature by a constant (KIT-224 later removed the shift; see the unreleased section).

The manifest field `correction$levels = "restored_to_tmm_logcpm"` (transcriptomics integration input) is a machine
token, not shown on a screen.
