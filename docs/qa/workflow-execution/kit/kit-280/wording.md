# KIT-280 user-facing texts: drafts for the owner's wording pass (VAL-41)

Each text is implemented and tested as written; the owner edits the wording. Where a text sits, so that an edit is found:

| Text | File | Tested by |
|---|---|---|
| Label "Low-intensity marking percentile (%)" | `R/mod_prot_qc_protein_intensity.R`; step label "Low-intensity marking percentile" in `R/utils_workflow_operations.R` | T10b |
| Help: "Values at or below this percentile of all values are marked in the results table; marking does not remove proteins. The figure quoted with the flag was measured at 1 %." | `R/mod_prot_qc_protein_intensity.R` | T10b |
| Help: "A value at or below the percentile below counts as present: only an absent value counts as missing." | same | T10b |
| DIA sentence: "On the DIA route the peptide intensity filter runs first, defaults to 2 groups and counts values below its own threshold as missing, so at the defaults this protein-level change removes almost no additional protein; the flag still marks proteins that hold a low value." | `PROTEIN_DIA_PEPTIDE_FILTER_SENTENCE` in `R/func_prot_qc_filtering_helpers.R` (tab help, DIA report, NEWS) | T10, T10b |
| Report sentence, marker present: "A value at or below the Nth percentile of all intensities in the table counted as present: only an absent value counted as missing. Proteins kept with such a value in one of the two groups compared are marked in the results table (column low_intensity_value). In simulated data (KIT-211; ten settings pooled), unchanged proteins kept with such a value in one of the two groups compared reached p < 0.05 1.9 times as often as other kept proteins (two-sided 95 % interval 1.83 to 1.96; 1.5 to 3.3 times in the individual simulated settings; measured at the 1st percentile only). Treat such calls with caution." | `proteinThresholdValuesSentence` in `R/func_prot_qc_filtering_helpers.R` | T10 |
| Report sentence, no record: "Whether intensities at or below the Nth percentile of all intensities counted as missing in this filter is not recorded for this analysis." | same; the substitution in `test-kit227-report-ruv-sentences.R` G1 | T10, KIT-227 G1 |
| Apply text line: "Proteins kept with a value at or below the threshold: N (marked in the results)" | `computeDIAReferenceProteinIntensityFilter` | T8 |
| NEWS bullet (Scientific changes), the replaced clause of KIT-187's bullet, the ledger rows `VAL-<merge>` and `VAL-<merge+1>`, the KIT-197 correction notes | `NEWS.md`, `docs/methods-decision-ledger.md`, `docs/qa/workflow-execution/kit/kit-197/` | diff |

The figure 1.89 [1.83, 1.96] is KIT-211's R-c under rule B: unchanged proteins kept with a low value in one of the two groups compared, against the other kept proteins, simulated data, ten settings pooled, measured at the 1st percentile only. Every text that quotes it carries "simulated data" and "1st percentile".
