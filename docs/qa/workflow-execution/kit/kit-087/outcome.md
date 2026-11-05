# KIT-087 outcome: NIPALS without the dense matrix, same result

**Imputing missing values at 200 × 18,285 now peaks at 0.63 GiB instead of 5.7 GiB and takes about 4 s instead of about 5 minutes, and the result is bit-identical to `mixOmics` wherever `mixOmics` finished. The one gap: `mixOmics` did not finish imputing LFQ and TMT at full size inside the 3,600 s ceiling, so those two are compared at 4,000 features.**

User direction, 2026-11-05: "bring that number down but return exactly the same scientific result for all omics. DONT RUN THE WHOLE CHAIN." No journey ran for this proof. Build: `kit_087` = 6pAjP4, 248261df.

## The change

`mixOmics::nipals` builds a feature-by-feature matrix on every iteration and reads only its diagonal: about 2.5 GiB at 18,285 features. `R/func_general_nipals.R` computes that diagonal directly from the mask of observed values. Nothing else about the method changes.

- `imputeNipalsMatrix` replaces `mixOmics::impute.nipals` at the four RUV imputation sites (`func_prot_s4_norm_methods.R`, `func_metab_norm.R`, `func_metab_s4_ruv_methods.R`, `func_lipid_s4_ruv_methods.R`).
- `pcaWithMissing` replaces `mixOmics::pca` at the six PCA sites that may see missing values (`func_metab_s4_plotting_methods.R`, `func_lipid_s4_plotting_methods.R`, `func_prot_s4_qc_plot_methods.R`, three in `func_general_plotting_pca_rle_helpers.R`). Complete data still goes to `mixOmics::pca`.

**PCA was added after the plan.** The plan left PCA alone, and I expected it not to matter. `kit-journey-086m` then stopped at the worker deadline inside `plotDensity`, which runs PCA on the full matrix. PCA with missing values uses the same dense matrix, so it now uses the same core.

## Same result

Each comparison ran each version in its own process on the stored state, with no stored file changed.

| Comparison | Assays | Result |
|---|---|---|
| Imputation, full size, stored or declared missing values | 10: DIA, DIA-Limpa, lipid large ×2, lipid screen ×2, metab large ×2, metab screen ×2 | all `identical` TRUE, largest difference 0 |
| Imputation, 4,000 features | LFQ, TMT | `identical` TRUE |
| PCA, full size | 12: the 10 above, LFQ, TMT | all `identical` TRUE |
| PCA, 4,000 features | lipid large ×2, metab large ×2 | all `identical` TRUE |
| Real sweep and correction, new code, compared with stored corrected tables | metab screen (k 1,1, 3% and 7%), lipid screen (k 2,2, 3% and 7%), metab large (k 1,1, 1% and 5%) | all TRUE |
| Final build's imputation digests against the first build's | 12 | all equal |

The KC real data was excluded. Only aggregates and hashes were recorded.

## The reduction

Peak memory is the process's maximum resident size in GiB (`ru_maxrss`, in KiB, ÷ 1,048,576), including about 0.38 GiB for R and the package. Each figure is one run on a shared host, so memory is firm and time is indicative.

**Imputation (2 components on lipid and metabolomics, 5 on the proteomics paths):**

| Assay | Samples × features | Old peak | New peak | Old time | New time |
|---|---|---|---|---|---|
| Metab large, GCMS | 200 × 18,285 | 5.71 GiB | 0.63 GiB | 287 s | 4.4 s |
| Metab large, LCMS pos | 200 × 18,285 | 5.73 GiB | 0.63 GiB | 324 s | 4.1 s |
| Lipid large, neg | 200 × 18,301 | 5.77 GiB | 0.63 GiB | 309 s | 3.2 s |
| Lipid large, pos | 200 × 18,208 | 5.73 GiB | 0.62 GiB | 272 s | 3.9 s |
| DIA | 297 × 2,500 | 0.66 GiB | 0.45 GiB | 1,755 s | 10.3 s |
| DIA-Limpa | 300 × 2,500 | 0.66 GiB | 0.46 GiB | 1,552 s | 13.9 s |
| LFQ | 300 × 8,000 | 2.04 GiB at the ceiling | 0.59 GiB | over 3,600 s | 46.5 s |
| TMT | 300 × 7,987 | 2.10 GiB at the ceiling | 0.60 GiB | over 3,600 s | 13.4 s |
| Screen corpora | 24 × about 200 | 0.38-0.39 GiB | 0.38 GiB | 0.2 s | 0.0 s |

The DIA and Limpa times are long because those assays take many iterations; the saving there is mostly time. At 4,000 features, LFQ took 1,317 s old and 3.2 s new, TMT 980 s and 2.4 s.

**PCA (2 components, full size):**

| Assay | Old peak | New peak | Old time | New time |
|---|---|---|---|---|
| Metab large, GCMS | 5.70 GiB | 0.59 GiB | 492 s | 2.5 s |
| Metab large, LCMS pos | 5.70 GiB | 0.59 GiB | 477 s | 2.8 s |
| Lipid large, neg | 5.71 GiB | 0.59 GiB | 611 s | 2.4 s |
| Lipid large, pos | 5.64 GiB | 0.59 GiB | 632 s | 2.6 s |
| LFQ | 2.01 GiB | 0.56 GiB | 2,547 s | 7.0 s |
| TMT | 2.08 GiB | 0.54 GiB | 217 s | 2.1 s |
| DIA | 0.66 GiB | 0.42 GiB | 185 s | 0.7 s |
| DIA-Limpa | 0.68 GiB | 0.43 GiB | 781 s | 3.8 s |

These are per call. What a whole journey saves is for KIT-086 and KIT-088 to show.

**Corrected figures.** Earlier reports and the plan gave these peaks as "6 GB" and the large metabolomics oracle as "10.3 GiB". Both divided KiB by 10⁶. The oracle's peak is 10,306,124 KiB = **9.83 GiB**, and the per-call peak is 5.7 GiB. The evidence files hold the raw KiB and are unchanged.

## Limits

- **LFQ and TMT imputation at full size were not compared.** `mixOmics` passed the declared 3,600 s ceiling on both, so `nipals-compare.txt` holds two `exit 124` lines and two read errors, annotated there. The evidence for these two is imputation at 4,000 features plus PCA at full size on the same core. The ceiling was not raised.
- **Bit-identity is shown on this R and BLAS only** (R 4.6.1, reference BLAS 3.12.0). On another BLAS `mixOmics` may itself round differently, so results may differ from it in the last digits there.
- **Only the sweep and correction on lipid and metabolomics were rerun from stored state.** Proteomics stores were recorded on older builds, so a mismatch there could not be attributed to this change. Proteomics relies on the matrix comparison.
- **One component is refused for imputation** (`ncomp < 2`), because `mixOmics::impute.nipals` crashes there. `pcaWithMissing` accepts one component, as `mixOmics::pca` does; identical on one 30 × 20 check, not a test in the suite.
- **Licence, for the user's review before release:** written from the published NIPALS method and checked against `mixOmics` as a black box. It is not a strict clean-room, because the `mixOmics` source was read to find the cause. MultiScholaR stays LGPL.

## Checks

- **Test first (`test-first.txt`):** `test-workflow-nipals-identity.R` on `kit_086`: all 8 tests error, as the function doesn't exist there. On `kit_087`: 77 expectations pass, 0 fail, 0 error. This was run at closure on the two installed builds, not recorded when the test was first written.
- **Dependents (`dependents.txt`):** 61 files, 8,932 pass, 0 fail, 0 error, 59 skip.
- **`R CMD check` (`check.log`):** at baseline, 3 WARNINGs and 3 NOTEs.
- **Test pins moved to the new functions:** `test-workflow-omics-ruv-cancor-once.R`, `test-metab-02z-ruv-cancor-characterization.R`, `test-prot-s4-norm-methods-direct-shared.R` (quarantined), `helper-module-ci-prot-norm.R`, `test-general-plotting-pca-helper-characterization.R`, `test-metab-01i-s4-plot-density-characterization.R`.
- **Attempts (`power-state.txt`):** one correction pass used. The first extraction named the wrong directory for LFQ and TMT, and no comparison had run.

## Handoff

`kit_087` goes to KIT-086's journeys. KIT-088 reruns all paths against `kit_086` for the exact reduction. Transcriptomics has no NIPALS call today; KIT-040 has a note to use these functions if it needs them.
