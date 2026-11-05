# Behaviour specification: NIPALS with missing values (KIT-090)

Written by the specifier, who has read `mixOmics`' source. This document gives behaviour only: what each function takes,
returns and refuses, and the mathematics. It contains no code. The implementer writes `func_general_nipals.R` from this
spec and the published method, and checks the result against `mixOmics` as a black box only.

Reference method: H. Wold (1966), "Estimation of principal components and related models by iterative least squares",
in P. R. Krishnaiah (ed.), *Multivariate Analysis*, Academic Press, pp. 391-420. NIPALS: Non-linear Iterative PArtial
Least Squares.

## What must hold

1. **Bit-identical results.** On this machine (R 4.6.1, reference BLAS), the outputs below must be `identical()` to the
   reference functions named for each, not merely close. Floating-point results depend on the order and the kind of each
   arithmetic operation. This spec deliberately gives the mathematics and not that order: find it by comparing against
   the reference. One general R fact that matters: `sum()`, `colSums()` and `rowSums()` accumulate in extended
   precision, and BLAS matrix products (`%*%`, `crossprod()`) do not. The two give different last digits for the same
   mathematical sum.
2. **Memory.** No intermediate may be a features-by-features or samples-by-samples matrix. At 200 samples by 18,000
   features, a features-by-features double matrix is about 2.4 GiB; the whole point of this file is to avoid that. Work
   in vectors and in matrices of the input's own size (samples by features), or smaller.
3. **Names and signatures** exactly as below. They are internal (`@noRd`) functions of the package MultiScholaR.

## Notation

`X` is a numeric matrix, samples in rows (`n` rows), features in columns (`p` columns), with `NA` for missing values.
A cell is *observed* when it is not `NA`. The missing pattern of the input never changes during the computation.

## 1. `nipalsComponents(X, ncomp, max_iter, tol)`

The components of `X` by NIPALS over observed cells only.

**Refusal.** If any column or any row of `X` has no observed cell, stop with the message
`Some rows or columns are entirely missing: remove them first.` (an error without the call, as `stop(..., call. = FALSE)`).

**Starting column.** One column index `s` is used to start every component: the column of `X` with the fewest missing
cells; on a tie, the lowest index.

**For each component `h` = 1, ..., `ncomp`**, with `R` the current residual matrix (`R` = `X` for `h` = 1; missing cells
of `R` stay missing):

1. Score vector `t` (length `n`): column `s` of `R`, with its missing entries set to 0.
2. `R0`: `R` with every missing cell set to 0.
3. Previous loading `q` (length `p`): every entry equal to 1 / sqrt(`p`).
4. Iteration counter `k` = 1 and change `c` = 1. While `c` > `tol` and `k` <= `max_iter`:
   - Loading `w` (length `p`): for each column `j`, w_j = (sum over all rows i of R0_ij * t_i) divided by (sum over the
     rows i where X_ij is observed of t_i squared).
   - Normalise `w` to unit Euclidean length: divide it by the square root of the sum of its squared entries.
   - Score `t`: for each row `i`, t_i = (sum over all columns j of R0_ij * w_j) divided by (sum over the columns j where
     X_ij is observed of w_j squared).
   - `c` = sum over j of (w_j - q_j) squared, ignoring missing terms; then `q` = `w`; `k` = `k` + 1.
5. If the loop ended because `k` exceeded `max_iter`, emit a message (not a warning):
   `NIPALS reached the iteration limit on component <h>.` with `<h>` the component number.
6. Deflate: `R` = `R` minus the outer product of `t` and `w` (missing cells stay missing).
7. Record: loading column `h` = `w`; singular value d_h = square root of (sum of t_i squared, ignoring missing terms);
   score column `h` = `t` / d_h.

**Returns** `list(scores = <n x ncomp matrix>, loadings = <p x ncomp matrix>, singular_values = <numeric, length ncomp>)`,
with no dimnames required on scores and loadings.

## 2. `imputeNipalsMatrix(X, ncomp = 2L, max_iter = 500L, tol = 1e-06)`

Missing values of `X` filled from the rank-`ncomp` NIPALS reconstruction.

- **No missing value:** return `X` unchanged (the identical object).
- **Refusal:** if `ncomp` < 2 or `ncomp` > min(n, p), stop (no call) with
  `` `ncomp` must be between 2 and the smaller dimension of `X`. `` (backticks included).
- Otherwise run `nipalsComponents(X, ncomp, max_iter, tol)`. The reconstruction is scores times diag(singular values)
  times the transposed loadings, an `n x p` matrix. Return `X` with each missing cell replaced by the reconstruction's
  value in that cell, and every observed cell, the dimnames and the attributes unchanged.

**Reference:** `mixOmics::impute.nipals(X, ncomp = ncomp)` (and with `max.iter`, `tol` where those are passed). The
result must be `identical()` to it.

## 3. `pcaWithMissing(X, ncomp = 2L)`

Principal components for the QC plots.

- **No missing value:** return `mixOmics::pca(X, ncomp = ncomp)` itself (one call; the result is that object).
- **Refusals,** in this order, each stop without the call:
  - any infinite value: `` `X` must hold finite values or missing values. ``
  - after rounding `ncomp` to the nearest integer, `ncomp` < 1 or `ncomp` > min(n, p): `Use a smaller `ncomp`.`
    (the text `Use a smaller ` followed by `ncomp` in backticks and a full stop).
- Otherwise:
  - Centre: subtract from each column its mean over its observed cells, as base R's `scale(X, center = TRUE, scale =
    FALSE)` does, keeping the attribute that function sets. Call the result `C`.
  - Run `nipalsComponents(C, ncomp, 500L, 1e-09)`.
  - Denominator `m` = the larger of 1 and n - 1.
  - Component names `PC1`, ..., `PC<ncomp>`.
  - Loadings: the components' loadings, rows named by `colnames(C)`, columns by the component names.
  - Variates: `C` with missing cells set to 0, times the loadings; rows named by `rownames(C)`, columns by the
    component names.
  - `sdev`: singular values divided by sqrt(`m`), named by the component names.
  - `var.tot`: the sum of the squared observed cells of `C`, divided by `m`.
  - Proportion explained: `sdev` squared divided by `var.tot`.
  - Return `list(X = C, ncomp = ncomp, sdev = sdev, loadings = list(X = <loadings>), variates = list(X = <variates>),
    prop_expl_var = list(X = <proportion>), var.tot = <var.tot>, cum.var = <cumulative sum of the proportion>)`.

**Reference:** `mixOmics::pca(X, ncomp = ncomp)` on data with missing values. These must be `identical()` to the
reference's fields of the same names: `variates$X`, `prop_expl_var$X`, `sdev`, `loadings$X`, `var.tot`, `cum.var`.

## Documentation in the file

Roxygen comments with `@noRd` on each function. The file header cites Wold (1966) and says the file was written in a
clean room from this specification (KIT-090).
