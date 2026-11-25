# KIT-248 wording: drafts for the owner (VAL-41)

Every text below is a draft. The owner sets the final wording. No on-screen text changes in this ticket: the Skip's
notice ("Correlation filtering SKIPPED. All samples retained.") and the supplement's step heading ("Skipped final sample
correlation filter", the state's description) are unchanged.

## 1. What the study parameters now print for a skipped correlation filter

Written by `createWorkflowArgsFromConfig` from the analysis record (the step's order number varies):

```
[4. correlation_filtered: pearson-skip]
  skipped [fixed] = TRUE
```

Before: after Apply then Skip, the block named `pearson-selection` with the pairs parameters and `threshold = <t>`;
after a Skip with no Apply, no block at all.

## 2. NEWS entry (### Fixes), as committed

> Proteomics: skipping the correlation filter is now recorded as a skip in the analysis record and in the study
> parameters ('pearson-skip', skipped = TRUE), and a replay of the record keeps every sample, as the run did. Before, a
> skip with no filter applied before it made the record impossible to replay, and a skip after an applied filter was
> recorded, and listed in the study parameters, as the filter at its threshold (a replay then removed the samples the
> filter had removed, so its differential abundance results differed from the reported ones). Results of the analysis
> itself did not change. Records made before this version are replayed as the skips they were (KIT-248).

## 3. Open for the owner

- The operation name `pearson-skip` appears to readers of the study parameters and the record. A plainer reader's label
  (for example "Correlation filter skipped") would need a step label, which this ticket does not add.
- Parity (R1-04): lipidomics and metabolomics skip steps still print nothing in the study parameters (KIT-254).
- The fixed constant prints as 'skipped' (no reader's label), in the study parameters and in the supplement, beside the
  internal codes `correlation_filtered` and `pearson-skip` (AUDIT-KIT-248-impl I09). The six reports (before and after
  KIT-236) print nothing for the step, as for a proteomics Apply.
- Read this NEWS bullet together with KIT-254's (AUDIT-KIT-248-impl I13): proteomics replays an older record whose
  outcome says skipped as the skip and says nothing (VAL-238); lipidomics and metabolomics replay as recorded with a
  statement of the samples left out (VAL-247). Consider "Analysis records made by the development builds" in place of
  "Records made before this version" in the last sentence ('before this version' can be read as covering `main`,
  which has no record).
