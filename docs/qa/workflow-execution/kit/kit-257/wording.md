# KIT-257 wording for the owner's pass (VAL-41)

Drafts only; the owner sets the final wording.

| Where | Draft (as on the branch) | Note |
|---|---|---|
| The normalisation tab, every omic: the run's refusal (`omicNormGuardTexts()$norm_unremoved`, `R/mod_omic_norm.R`; shown by the run rule's notice in transcriptomics, lipidomics and metabolomics, and by proteomics' error handler as "Error in normalization: ...") | An earlier file of the normalisation could not be removed, so this normalisation was not run: <paths>. Close any program that has it open, then run the normalisation again. | The ticket's section 7 text. Shown only when a file of the stage cannot be removed (in practice a file open in another program on Windows, or a folder of that name). |
| Proteomics, the correlation filter's refusal (`omicNormGuardTexts()$correlation_unremoved`; shown as "Error in correlation filtering: ...") | An earlier file of the normalisation could not be removed, so this correlation filter was not run: <paths>. Close any program that has it open, then apply or skip the correlation filter again. | Added at step 0 (r1.2, section 7): the same text naming its own step. |
| Proteomics, the error notice of a normalisation or correlation filter that stops after its start removed files (`omicNormGuardTexts()$removed_then_stopped`, appended to the first line of "Error in normalization: ..." / "Error in correlation filtering: ...") | The earlier RUV figures and normalised table of this analysis were removed when this run started; apply or skip the correlation filter again before Copy. | AUDIT-KIT-257-impl I05, the coordinator's ruling that the user is told; the register's draft. |
| The normalisation tab's log (transcriptomics, lipidomics, metabolomics) | Removed the earlier run's files this run replaces: <file names> | One line, only when a file was removed (R1-02). Proteomics has no tab log: the logger line only. |
| The logger (every omic) | Removed the earlier run's files: <paths> | Diagnostic; not user-facing. |
| `NEWS.md`, unreleased section, `### Fixes` | Normalisation, every omic: each run now first removes the earlier run's files that depend on RUV (and, in proteomics, the earlier normalised table), and writes again those it makes. Before, Copy and the report could deliver files from an earlier, different normalisation: (three cases). A file that cannot be removed (for example, open in another program) stops the run before the normalisation is computed. If a proteomics normalisation or correlation filter stops after it has started (for example, a manual k larger than the data allow), the analysis is unchanged but its RUV figures and normalised table are no longer in the project, and the error notice says so: apply or skip the correlation filter again to write the table, and run the normalisation again for the figures, before Copy. On DIA with limpa, pressing "Run Normalization & RUV" again on the correlation's or the RUV covariates' result keeps that result's files. Projects normalised before this release keep such files until the normalisation is run again. (KIT-257) | The ticket's section 12 text with Re-check 1's R1-02 replacements. It names the three cases only (F1, the transcriptomics RUVIII export, is not fixed here). |

## Open wording defects of the I05 notice (closure check C1, C2; for the owner's pass)

The notice always says "The earlier RUV figures and normalised table of this analysis were removed", although the
removed paths vary (E10 itself removed only the table; a correlation start never removes figures), and it is joined to
the kernel's first line by a single space. Drafts (the register's):

| Case | Draft |
|---|---|
| An `ruv_corrected_` file was removed | The earlier RUV figures and normalised table of this analysis were removed when this run started; apply or skip the correlation filter again, and run the normalisation again for the figures, before Copy. |
| Otherwise | The earlier normalised table of this analysis was removed when this run started; apply or skip the correlation filter again before Copy. |
| Separator (C2) | join with ". " when the kernel's first line does not end in punctuation |

The small code change (choosing the string by the removed paths, which the condition carries as `removed_own_files`,
and the separator) follows the owner's wording pass; KIT-257's correction pass is spent.
