# KIT-198 outcome: stopped at gate G0 (D8)

Ticket `tickets/kit/KIT-198.md` r1.2 (owner decisions D1 to D8 as recommended, methods ledger VAL-118). Branch
`kit-198-intensity-filter-wording` from `v0.9-staging` 19997323. By the coordinator's order G0 ran first and alone,
on a gate commit that adds only `g0-numeric-columns.R` (nothing under `R/`, `man/` or `tests/`).

## What was run

| Run | Commit | What |
|---|---|---|
| 45314460795 (`kit198-g0-b0e913b8`) | b0e913b8 | G0: `Rscript docs/qa/workflow-execution/kit/kit-198/g0-numeric-columns.R`, one scripts job, R 4.6.1, reference BLAS, ubuntu-22.04; exit 0 after 7 s of script time; 10.8 job-minutes in all (0.72 core-hours) |

Returned: `g0/g0-numeric-columns.tsv`, `g0/g0-threshold-effect.tsv`, `g0/g0-decision.txt`, `g0/console.txt`.

## Observed

Six inputs, all run (none failed). Each went through the configuration's preview, the import panel's column
resolvers, `importOmicsTables()`, `buildOmicsDesignSaveDataList()` with the input's design, and the configuration's
`design$build_object()`; the design object's assay tables count the same numeric columns as the design save's
(`object_counts_same` TRUE on every row).

| Input | Assay | Id column (class) | Columns the filter counts | Non-sample numeric |
|---|---|---|---|---|
| `lipid_canonical` (LipidSearch) | LCMS_Pos, LCMS_Neg | LipidName (character) | the 6 samples | none |
| `metab_lc` (custom) | LCMS_Pos, LCMS_Neg | Feature.Name (character) | the 6 samples | none |
| `metab_gc` (custom) | GCMS | Feature.Name (character) | the 6 samples | none |
| `metab_combined` (custom) | LCMS_Pos, GCMS | Feature.Name (character) | the 6 samples | none |
| `module_ci_lipid_msdial_assay()`, `Peak ID` = 1..3, MS-DIAL lipid import | LCMS_Pos | Peak ID (numeric) | Peak ID + 4 samples (5) | **Peak ID** |
| `tests/testdata/e2e/metab_lc/seed_lcms_pos.csv`, MS-DIAL metab import | LCMS_Pos | Alignment ID (numeric) | Alignment ID + 4 samples (5) | **Alignment ID** |

The MS-DIAL default id was proposed by the preview, kept numeric through the import and the design save, and counted
by the filter as a fifth sample. The other numeric MS-DIAL columns (retention time, m/z, score) were dropped by the
design save, as the audit read.

Threshold effect (reported only; `g0-threshold-effect.tsv`). Thresholds over the counted columns against over the
design's samples only; features removed at proportion 0.5 and the 10th percentile under each:

| Input | Features | p1 counted / samples only | p10 counted / samples only | Removed at 0.5, p10: counted / samples only |
|---|---|---|---|---|
| `lipid_canonical` LCMS_Pos | 10 | 20568.17 / 20568.17 | 30111.8 / 30111.8 | 1 / 1 |
| `lipid_canonical` LCMS_Neg | 10 | 18526.8 / 18526.8 | 25782 / 25782 | 2 / 2 |
| `metab_lc` LCMS_Pos, LCMS_Neg | 3 each | 2439.358 / same | 2583.5 / same | 0 / 0 |
| `metab_gc` GCMS | 3 | 4927.253 / same | 5100.24 / same | 0 / 0 |
| `metab_combined` LCMS_Pos, GCMS | 3 each | 6499.901 / same | 6697.865 / same | 0 / 0 |
| MS-DIAL lipid (module-CI, 4 samples) | 3 | 1.12 / 100.09 | 2.2 / 100.9 | 0 / 0 |
| MS-DIAL metab seed (4 samples) | 8 | 1.39 / 9731 | 4.9 / 16410 | **0 / 2** |

## Inferred

- On an MS-DIAL import with its default integer id, the filter's denominator is n + 1 and its threshold is taken
  partly from the id numbers: at the 1st and 10th percentiles the threshold falls among the ids (1.39 and 4.9 against
  9731 and 16410 over the samples), so no intensity is below it and the filter removes nothing that a
  samples-only threshold would remove (2 of 8 features on the seed). This is F2 and U1 of the ticket, now shown.
- "Samples" in W1 to W4 is not exact for every application input. By the ticket's G0 decision (section 9) and the
  owner's D8 (VAL-118 (1)), the ticket stops before the fix commit; the correction (the filter and its percentile on
  the design's samples only) is its own ticket, before the release commit and before KIT-165, KIT-166 and KIT-170 if
  their inputs carry such an id.

## Acceptance criteria

| # | State | Evidence |
|---|---|---|
| A0 | G0 ran and found a numeric non-sample column on two inputs; the ticket stopped (the owner's D8 answer is already (i) with (ii)) | `g0/g0-numeric-columns.tsv`, run 45314460795 |
| A1 to A7 | not started (no tests commit, no fix, no R2, C, R3) | – |

## Deviation from the ticket's text

G0 ran on a gate commit (b0e913b8) that holds only the G0 script, not on the tests-only commit, by the coordinator's
order (G0 first and alone). The tree under `R/`, `man/` and `tests/` is `v0.9-staging` 19997323's.

## Not checked

- Real MS-DIAL exports (only the repository's seed and the module-CI fixture with its id made integer); whether the
  lipid and metabolite inputs of KIT-165, KIT-166 and KIT-170 carry such an id.
- A user who picks a text column as the id, or the full GUI-entry driver on an MS-DIAL input (the panel's functions
  were called directly with the preview's proposals, which the browser double applies).
- The partitioned variant's own column read (`partition$numeric_columns()`); it was not run. The design object's tables
  were checked, not the state store's copy.
- The threshold effect at other proportions or percentiles.
