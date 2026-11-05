# KIT-086 outcome: the complete-controls rule, and missing values at size

**With the KIT-087 fix, the large journeys with missing values now finish and match their records exactly on both omics. Metabolomics peaks at 3.94 GiB of anonymous memory (5.17 GiB whole slice) in 1,220 s, the same as complete data. The same journey on the build before the fix was stopped at the worker deadline at 2,124 s, having already reached 10.42 GiB. The complete-controls rule is in place but never fired on these corpora, so only its unit test exercises it.**

User direction, 2026-11-05: "i think do all these before we move on no?", then the KIT-087 sidestep: "when thats done we get to dive back into here, run it again, close it out". Synthetic corpora: no real-data claim. No default changed.

## The large journey, before and after the fix

Same corpus: 200 samples, 20,000 features per assay, 2 assays, the declared missing values. `kit-journey-086m` ran on `kit_086` = aKld9N, `kit-journey-086p` on `kit_087` = 6pAjP4.

| Step | Before: time | Before: peak | After: time | After: peak |
|---|---|---|---|---|
| Import | 19 s | 1.15 GiB | 10 s | 1.15 GiB |
| Design | 80 s | 1.23 GiB | 76 s | 1.22 GiB |
| Intensity filter | 113 s | 3.25 GiB | 54 s | 3.23 GiB |
| Duplicates | 60 s | 2.28 GiB | 26 s | 2.26 GiB |
| QC finalise | 32 s | 2.45 GiB | 19 s | 2.45 GiB |
| Normalisation | **stopped at 1,804 s** | **10.42 GiB** | 781 s | 3.94 GiB |
| Correlation | did not run | | 27 s | 2.12 GiB |
| Session export | did not run | | 17 s | 2.19 GiB |
| DA | did not run | | 81 s | 3.13 GiB |
| Report | did not run | | 55 s | 3.62 GiB |
| **Whole journey** | **stopped, 2,124 s** | **10.42 GiB; 10.52 whole slice** | **complete, 1,220 s, exact** | **3.94 GiB; 5.17 whole slice** |

Peak is the slice's anonymous memory.

**Three limits on reading this table:**
- **No finite ratio.** The "before" run did not finish, so there is no speedup figure, and its peak is only what it had reached. It stopped in the first QC figures, before RUV.
- **The early steps aren't a time comparison.** Load was 5-8 during the "before" run and about 2 during the "after" run, and the code for those steps is unchanged.
- **Which step stopped it.** The first `metab_qc_figures` step ran 1,800.4 s before the worker deadline. On the fixed build the same step takes 29.7 s. `plotDensity` runs PCA with missing values on the full matrix, per assay, which is why PCA joined KIT-087's scope.

## All members on the fixed build

| Member | Corpus | Decision | Exact against | Files | Time | Peak | Whole slice | Application |
|---|---|---|---|---|---|---|---|---|
| `086p` | Metab large | accepted | yes, `metab86` | 199 | 1,220 s | 3.94 GiB | 5.17 GiB | 3.26 GiB |
| `086q` | Lipid large | accepted | yes, `lipid86` | 208 | 1,246 s | 4.13 GiB | 5.24 GiB | 2.92 GiB |
| `086u` | Lipid screen | accepted | yes, `lipid86s` | 208 | 186 s | 1.80 GiB | 1.94 GiB | 1.01 GiB |
| `086v` | Metab screen | accepted | yes, `metab86s` | 199 | 180 s | 1.85 GiB | 1.95 GiB | 0.86 GiB |

**What "exact" proves differs by member:**
- **`086p`, `086u`, `086v`: old code against new.** Their oracles were recorded on `kit_086`, with `mixOmics`. So on these corpora the new imputation and PCA give the old scientific result through a whole journey.
- **`086q`: journey against in-process, both new.** There was no `kit_086` oracle for lipid at large, because the declared order ran it only after metabolomics completed, and by then the fix was in. Old-against-new for large lipid rests on KIT-087's matrix comparison, which covers imputation and PCA on both assays.

No large lipid journey with missing values ran on the old build, so there's no before-and-after for it here. KIT-088 is due to measure that.

## Missing values against complete data

| Journey | Whole slice | Peak | Time |
|---|---|---|---|
| Metab large, complete (`084c`, `kit_084`) | 5.36 GiB | 3.97 GiB | 1,235 s |
| Metab large, missing values (`086p`) | 5.17 GiB | 3.94 GiB | 1,220 s |
| Lipid large, complete (`084d`, `kit_084`) | 5.63 GiB | 3.79 GiB | 2,136 s |
| Lipid large, missing values (`086q`) | 5.24 GiB | 4.13 GiB | 1,246 s |

On this pattern of missing values, imputation no longer shows in the peak or the time. The lipid times don't show a speedup: `084d` ran with other work on the host, and its RUV took 1,065 s against 617 s here on the same third-party code. On a quiet host RUV takes about the same on both omics (603 s and 617 s), so lipid isn't inherently slower, which answers the question I'd left open in KIT-084.

## What the ticket asked to report

- **Where the peak goes.** Before: full-matrix NIPALS in the figures, 5.7 GiB per call on top of the worker's other data. After: the highest point is the worker in normalisation. The application's own peak stays at 3.26 GiB (metabolomics) and 2.92 GiB (lipid) at the report step, unexplained and still open with the user.
- **The 1,800 s deadline.** Every step now fits. The largest are RUV at 603 s (metabolomics) and 617 s (lipid), and the former QC-figure stop is 30 s.
- **RUV's time.** RUV is still 77% of normalisation, almost all in the third-party `RUVIIIC`. This ticket didn't change that.
- **How often the rule skips (`rule-outcome.txt`).** Never: 0 of the 240 percentages across 12 assay sweeps (6 journeys, 2 assays each). The chosen k is 1 or 2, and these corpora have many complete controls. The screen corpora's only skips are 1-2 percentages per assay from the existing `skipped_insufficient_controls`.

## Corrections to the record

- **"At least 13 hours" and "at least 40 minutes" (KIT-085's outcome)** came from uniform missingness. Measured on the declared clustered pattern, old imputation took about 5 minutes per assay. Corrected in that outcome.
- **"10.3 GiB" for the large metabolomics oracle** is 9.83 GiB (10,306,124 KiB). The oracle is one process holding every step, so it isn't comparable with a journey's slice.

## Deviations from the frozen declaration

- **One build, one attempt each.** The small and large runs on `kit_086` stand as recorded (`086t`, `086s`, `086m`). The rule held, and the deadline stop was a valid result. Under the user's direction, the fix went in through KIT-087 and the journeys reran on `kit_087` under new labels.
- **The lipid large oracle** was recorded on `kit_087`, not `kit_086`, as above.

## Checks

- **The rule's test:** `test-workflow-omics-ruv-complete-controls.R`, 16 pass. Since the rule never fired in a journey, this test is its only evidence.
- **Corpora:** regenerated `metab_screen` is byte-identical at 0 missing values.
- **Complete data:** `086l` and `086n` accepted and exact on `kit_086`. Not rerun on `kit_087`: complete data returns unchanged from the imputation and goes to `mixOmics::pca` as before.
- **Dependents:** 52 files, 8,774 pass, 0 fail (`dependents.txt`, `kit_086`). On the final build, KIT-087's 61 files, 8,932 pass, 0 fail.
- **`R CMD check`:** at baseline, 3 WARNINGs and 3 NOTEs, on both builds.
- **Full suite:** not run since KIT-084. KIT-085, 086 and 087 are three tickets; under the test rule the next checkpoint is due when KIT-040 closes.

## Handoff

KIT-088 reruns all paths against `kit_086`, where a member stopped by a limit is plotted as not finishing. Still open with the user: whether any default moves to `partitioned`, and the report step's application peak. Next by the selector: KIT-040.
