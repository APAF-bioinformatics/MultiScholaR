# KIT-090 outcome: a clean-room NIPALS, same results

**`R/func_general_nipals.R` is now code written in a clean room: an implementer that never saw `mixOmics`' source or the previous file wrote it from a behaviour-only specification. On this R and BLAS, it gives bit-identical results to `mixOmics` everywhere KIT-087 checked. MultiScholaR stays LGPL.**

User decision, 2026-11-06: stay LGPL by a clean-room rewrite, rather than moving to GPL. Build: `kit_090` = y0OrcU, bbe1fae9.

## How it was written

- **Specifier (me, who had read `mixOmics`):** [`nipals-spec.md`](nipals-spec.md), prose and equations only. It covers the method (Wold, 1966), the interface, the refusals and messages, the memory rule, and the output conventions to match. It leaves out the floating-point order of operations.
- **Implementer:** a fresh agent working in a separate folder holding only the spec, the identity test and a loader that swaps its functions in for the installed ones ([`load-candidate.R`](load-candidate.R), [`run-tests.R`](run-tests.R)).
  - **Forbidden:** `mixOmics`' source, the installed package's functions, R library files, the web and this repository.
  - **Allowed:** calling `mixOmics` as a black box only.
- **Its log** ([`implementer-log.md`](implementer-log.md), 10 entries): it read only the spec, the test and the helpers. It compared 32 orderings of the sums against `mixOmics`' outputs ([`implementer-explore-impute.R`](implementer-explore-impute.R)), and only all-BLAS sums inside the iteration matched bit for bit.
- **The result:** the file builds nothing larger than the input's own size. It keeps no transposed mask, so it uses slightly less memory than KIT-087's version.

## Same results

| Check | Result |
|---|---|
| Identity test on the installed build | 77 of 77 |
| Imputation, 10 stored assays, against `mixOmics`' saved results | all `identical`, largest difference 0 |
| Imputation digests, all 12 assays, against KIT-087's | all equal |
| LFQ and TMT imputation at 4,000 features, against a fresh `mixOmics` run | `identical` (old 1,391 s and 1,006 s, new 2.7 s and 2.0 s) |
| PCA, 16 comparisons against `mixOmics`' saved results | all `identical` |
| Real RUV sweep and correction from stored state | corrected tables identical on metab screen, lipid screen and metab large, with the same k and percentage |

## Checks

- **Tests:** the identity test and the five RUV and PCA test files that use it pass on the source (162 expectations).
- **Dependents:** 61 files, 8,932 pass, 0 fail.
- **`R CMD check`:** at baseline, 3 WARNINGs and 3 NOTEs.
- **Full suite checkpoint** (KIT-085 to 090, one process per file, tree 86a32799): 363 of 363 files, 31,484 pass, 0 fail. The 2 errors are the same two known diagnostic tests that error at every checkpoint (`test-workflow-perf101-diagnostic.R`, `test-workflow-ruv-pearson-diagnostic.R`).
- **Attempts:** one correction pass, recorded in `power-state.txt`. My runner didn't set the working-folder variable for the three RUV-output comparisons, so they were rerun.

## Limits

- **A clean room run by an AI model.** The implementer is the same kind of model, and has probably seen `mixOmics`' public source in training. This is a documented independent derivation, not a legal guarantee. The institution's IP contact should make the final call before release.
- **Bit-identity is shown on this R and BLAS only.**
- **Two spec points the implementer read literally:**
  - `imputeNipalsMatrix` doesn't round a fractional `ncomp`; it only range-checks it.
  - The start column is chosen once, which makes no difference while the missing-value pattern is fixed.
