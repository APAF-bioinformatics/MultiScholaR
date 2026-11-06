# KIT-099 sizes, fixed by the stepped trial (declaration sections 2 to 4), 2026-11-06

Recorded before the large oracle (declaration section 5, item 5). Trial records: `size-trial.txt` (one line per point), `power-state.txt` (host per point). All points ran in process on `kit_098`, one attempt each, in ascending order. The trial stopped at 800, the first point outside the limits. A trial point is not a measurement.

| Size | Corpus key | Samples | Raw genes | Retained genes | Raw cells | Retained cells | Trial point | Finished | Elapsed | Peak RSS | Manifest sha256 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| small | `trx_60k_024` | 24 | 60,000 | 21,079 | 1,440,000 | 505,896 | 24 | yes | 58 s | 1.68 GiB | `651f0a4b15f7…` |
| medium | `trx_60k_100` | 100 | 60,000 | 21,133 | 6,000,000 | 2,113,300 | 100 | yes | 132 s | 2.48 GiB | `6d1de529daaf…` |
| large | `trx_60k_200` | 200 | 60,000 | 21,247 | 12,000,000 | 4,249,400 | 200 | yes | 259 s | 2.83 GiB | `c233aaa00786…` |
| (between) | `trx_60k_400` | 400 | 60,000 | 21,068 | 24,000,000 | 8,427,200 | 400 | yes | 780 s | 4.14 GiB | `c21a8df54edf…` |
| **ultra** | `trx_60k_600` | 600 | 60,000 | 21,280 | 36,000,000 | 12,768,000 | 600 | yes | 2,056 s | 5.75 GiB | `eb5e8e7746c3…` |
| limit | `trx_60k_800` | 800 | 60,000 | 21,273 | 48,000,000 | 17,018,400 | 800 | **no**: stopped at the 1,800 s step limit in normalisation (TMM/logCPM, grid, RUVs) | 2,160 s | 5.30 GiB | `d767df400708…` |

**Recipe (every size):**
- `kit_omics_corpus.R` `trx_screen`, seed 20261106, 60,000 genes, low-count fraction 0.6;
- two groups (RPMI, Sera), two technical replicates per biological sample, three batches;
- 10% of the genes differential, 30% with a batch effect;
- 50 all-zero genes, 1 constant gene;
- biological samples per group: 6, 25, 50, 100, 150, 200.

**Retention:** 35.1 to 35.5% at every size, as declared (about 35%).

**Per-step seconds in process (`size-trial.txt`):**
- normalisation: 26, 57, 123, 499 and 1,635 s, and at 800 more than 1,810 s (stopped);
- DA: 11, 29, 51, 106 and 139 s;
- design: 5 to 181 s.

Normalisation grows about as n^2.5 between 200 and 600 samples. It is the step that ends the curve.

**Against the two-assay omics on KIT-088's axis (matrix cells):**
- KIT-084's lipid and metab large: 2 assays x 20,000 features x 200 samples = 8.0 M cells;
- KIT-068's screens: 9,600 (lipid) and 14,400 (metab).

Transcriptomics large is 12.0 M raw cells and 4.25 M retained cells at the same 200 samples. Its small is already 1.44 M raw cells: the gene table is held at a realistic size and only the samples step.

**Not in the trial:** the QC figures, the report and the session export are GUI publications that the in-process recorder does not make. The large journey measures them.

**Proxy caveat:** the in-process step limit is conservative. The in-process normalisation step holds three worker steps (TMM/logCPM, grid, RUVs) that each have their own 1,800 s deadline in a journey. A managed run at 800 samples might finish where this proxy stopped, or might not. The point's peak RSS (5.30 GiB) was far from 15.5 GiB, so time ended it, not memory.

The limit is a transcriptomics RUV or grid ceiling candidate for KIT-091 (the call profile below names the call).
