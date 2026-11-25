# KIT-236 gate G0 (ticket r1.2, section 4)

Status: G0 passed as predicted; the ticket stops here. The tests-only commit waits for the owner's answer to OD-1. D8 and D19 wait for KIT-242 (the owner's decision, then its merge; R1-N08).

## What was run

- Run 45577606908. Package commit 10198a0d on branch `kit-236-g0`: KIT-235's head dc2c638f (staging 3eeaf100, KIT-227 5e7efe5f, the KIT-235 tests and the parser fix F1) plus `tools/kit236/g0.R`. No package change; the branch is never merged.
- Preconditions (`g0-45577606908/meta.txt`): R 4.6.1, pandoc 2.9.2.1, C.UTF-8.
- G0a: `guiEntryRunLane("lipid_canonical")` and `guiEntryRunLane("metab_lc")` (70 s and 67 s, 20 captures each). A tracer on `guiEntryArtifacts` copied the lane's project (store excluded) at each capture. The report docx was read as text (`</w:p>` to a line break, other tags removed, curly quotes made straight, whitespace collapsed), and each section-1 sentence was matched by a fixed phrase.
- G0b: the copied project, one condition changed, rendered with `RenderReport` as the Summary calls it (`project_dirs` re-pointed to the copy). V0 is an unchanged copy, added as an instrument control. V1 to V7 are as section 4 lists them.
- G0c: the templates' function definitions (each chunk's code parsed, function definitions evaluated in a child of the namespace) called directly, and the `QC_metrics` chunk evaluated. A plotting-scale fixture was added. Part B covered the four proteomics templates.

## Observed

- **G0a, the lanes:** all 54 lipid and 55 metab sentence rows are as predicted (`g0a/<lane>-sentences.tsv`, with the paragraph each one prints in).
  - Present: S11 to S29/S30, S30/S31, M1 to M4b, M6, M7 to M15, F1, F2, F3.
  - Present for metab only: F4.
  - Absent: S1 to S10, M5, F5; F4 is absent for lipid.
- **Record (`g0a/<lane>-record.txt`, `-steps.tsv`), both lanes:**
  - Normalisation is `cyclicloess`, the log2 offset is 0, and there is no RUV step. The ITSD step's outcomes include `itsd_applied: true`.
  - DA parameters: formula `~ 0 + group`, qvalue 0.05, lfc 0. Contrasts: one, `KO_vs_WT=groupKO-groupWT`.
  - DA `model` verbatim: `{"model":"limma","formula_used":"~ 0 + group","eBayes_trend":true,"eBayes_robust":true,"treat_applied":false,"block_reason":"no_duplicates","uses_duplicate_correlation":false}`. `block_reason` is one string, as R1-N04 read.
  - DA `results` hold `{contrast, assay, entity_id, rows}` per assay, with no `fdr_method`.
  - Software: version 0.9.0; build recorded; sha null.
  - Design: `technical_replicate_id = tech_rep_group`. There are 6 record rows with 6 distinct ids, so no technical replicates are declared. `design_matrix.tab` has 6 rows, so n = 6 and k = 6.
  - The `<prefix>_pearson_selection` step is present with threshold 0.5.
  - Step dimensions, total rows across the two assays:
    - lipid: 20 after design, 17 after the intensity filter and at every later step;
    - metab: 6 throughout.
  - Columns are 8 for lipid and 7 for metab, which is samples plus annotation columns, not samples.
- **File (`g0a/<lane>-file.txt`):**
  - The ITSD section reads Applied Yes, average_centered, median, removed Yes. There is **no** `ITSD Features Per Assay` block.
  - There is no DA UI block under either name, and no `User Interface Parameters` section.
  - The log section reads `* Transform: log2` (no Offset line).
  - Version Information lines are written (Repository, Branch N/A, Commit, Commit Timestamp).
- **Parser (`g0a/<lane>-parsed.txt`):**
  - The extractor has 41 fields, and AnalyzeSamples returns 84 numbers with the record overlay.
  - `ui_de_treat_enabled` is FALSE (the block is absent), so M14's cutoff rule is never reached. `da_block_reason` is `no_duplicates`. `software_version` is 0.9.0.
- **G0b (`g0b/summary.tsv`):** every variant rendered on both lanes, and every check of section 4's prediction is as predicted.
  - V0 reproduces all 54/55 sentence rows.
  - V1 prints "using the 'none' method from the 'limma' R package" with "Cyclic loess normalization (rightmost column".
  - V2 prints "Lipid intensities were transformed as log2(x + 1). Missing values (zeros or non-detects) were coded as NA …".
  - V3 still prints "6 biological samples with no technical replicates."
  - V4: the volcano and count images are drawn (drawings 2→3 for lipid, 1→3 for metab) and no "not found" line prints.
  - V5 still names `DE_<…>_results.xlsx`.
  - V6 prints "…applied using the average_centered method with median aggregation across ITSD features." with no word on the other assay.
  - V7 prints "The standard empirical Bayes method was used to test for any deviation from zero log-fold-change." It has no fold-change clause, although the overlay gives `ui_de_log_fold_change_cutoff` 0.5.
- **G0c (`g0c/summary.tsv`):** every case with a prediction is as predicted.
  - With RUV applied, S1 to S10 and M5 print.
  - With three groups, "comparing each pair" prints.
  - With `block_reason` (single or joined), no block clause prints.
  - With the composite absent, "Composite QC figure not found" prints.
  - Part B, all four templates:
    - B1 prints "no technical replicates" for a design with shared `tech_rep_group`;
    - B2 prints "'none' method from the 'limma'";
    - B3 prints no "Benjamini";
    - B4 (lfc 0.5) prints the standard empirical Bayes sentence;
    - B5 prints no duplicateCorrelation;
    - B6 (limpa, non-limpa branch) prints "Technical replicate averaging", with and without `block_reason`.
  - R1-N03 (limpa non-limpa branch, file cutoff 0.5, no record): both forms print the standard empirical Bayes sentence, not TREAT. Form (b) adds the writer's `[Differential Expression UI Parameters]` block.
  - The proteomics default texts are recorded (`g0c/<T>-default.txt`).
- **Plotting scale (`g0c/plotting-scale.txt`):**
  - No line of `plotPca`, `plotRle`, `plotDensity`, `plotPearson` or their helpers names a log function (lipid and metab).
  - On a fixture of raw intensities, plotRle's per-sample medians are −250, −100, 100 and 500. These are the raw-scale values; log2 values would lie within ±2. PCA scores are in the thousands.

## Inferred

- KIT-235's parser returns a written `* ITSD Features Per Assay:` block as its **own top-level section** (`'ITSD Features Per Assay'` appears in V6's section list). So the extractor's `itsd_features_per_assay` stays "NA", and the template's "The ITSD features used per assay were:" clause never prints. D19's per-assay source must read that section.
- The composite's left column is drawn on the post-filter intensities as they are, with no transformation. That agrees with D1's left-column wording.
- The Version Information block prints values once Main runs. The P16 hand-off's "print empty" described the pre-KIT-235 state, when Main failed.

## Acceptance criteria

- G0 "observed as predicted": **met**, by run 45577606908 and the tables above.
- Section 5's rows: not started. The tests-only commit waits for OD-1.
- R2 and R CMD check: not started.

## Not checked

- Whether `qvalue()` fell back to BH on the lanes. The record holds no `fdr_method`, and the log was not searched for `qvalue_warnings`.
- The parsed contents of the top-level `ITSD Features Per Assay` section; only its presence among the section names was recorded.
- Scripted (non-GUI) lipid/metab routes.
- `plotPearson` and `plotDensity` were not run on the fixture; only their source was checked.
- V0 to V7 render seconds read 0 in `summary.tsv`. This is likely a rounding or timing artefact of the measurement; the texts differ per variant as the edits require.
