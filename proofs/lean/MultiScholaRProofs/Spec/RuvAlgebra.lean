module

public import Mathlib.LinearAlgebra.Matrix.NonsingularInverse
public import Mathlib.LinearAlgebra.Matrix.Rank
public import Mathlib.Basic.Real.Basic

/-!
# Group R, the RUV-III algebra (KIT-193 sections 7.1, 7.2)

Notation of section 7: `Y` runs × features (as `ruv` takes it); `M` the runs × replicate-sets matrix;
`R_M = I - M (M'M)⁻¹ M'`; `ctl` the control columns; `k` the number of factors. RUV-III (`ruv::RUVIII`, ruv 0.9.7.1
(the locked version), `R/RUVIII.R:35-42`, read at the statements gate: the lines are the model's): `Y0 = R_M Y`; `U` the first `k` left singular
vectors of `Y0` (the solver's output, abstracted: a field of the structure, constrained only by the hypotheses each
statement names); `α = U'Y`; `a_c = α[, ctl]`; `W = Y[, ctl] a_c'(a_c a_c')⁻¹`; `Y* = Y - Wα`. The package's centred
form (`computeRuvCovariates`, `R/func_prot_norm_optimization_helpers.R:1172-1183`; `computeOmicsRuvAssayFactors`,
`R/func_omics_ruv_factors.R:170-173`): `W_c = (Y[, ctl] - 1 m_ctl') a_c'(a_c a_c')⁻¹`, adjusted `= Y - W_c α`.
`restoreRuvFeatureLevels` (`R/func_prot_norm_optimization_helpers.R:1133-1140`): each feature shifted so its mean
over its finite corrected runs equals its mean before; on complete data the mean over all runs.

Every matrix is real and complete (no missing value): RUV-III-C with missing values is `Spec/RuvIIIC.lean`.

Hypotheses used below, each with its reason (section 5.1: no junk values):
* `IsUnit (Mᵀ * M).det`: every replicate set is non-empty (`ruv::replicate.matrix` builds no empty set); makes
  `residop M` the real projection.
* `IsUnit (ac * acᵀ).det`: mirrors `solve(ac %*% t(ac))` failing on a singular matrix (the code's refusal).
* `residop M * U = U`: the solver's `U` lies in the replicate-residual space (true of the left singular vectors of
  `R_M Y` whenever `R_M Y` has rank at least `k`; `svd` returns arbitrary directions for zero singular values).
* `residop M *ᵥ 1 = 0`: every run is in exactly one replicate set (true of an indicator `M` whose sets are non-empty).
-/

@[expose] public section

namespace MultiScholaRProofs.Ruv

open Matrix

/-- `R_X = I - X (X'X)⁻¹ X'`, the projection off the columns of `X`. Every statement that uses it at a variable `X`
assumes `IsUnit (Xᵀ * X).det` (otherwise Lean's inverse is the junk `0` and this is `1`), except where a docstring
says the conclusion is meant for that case too. -/
noncomputable def residop {n : ℕ} {m : Type} [Fintype m] [DecidableEq m] (X : Matrix (Fin n) m ℝ) :
    Matrix (Fin n) (Fin n) ℝ :=
  1 - X * (Xᵀ * X)⁻¹ * Xᵀ

/-- The column means of a matrix, `1'A / n`. At `n = 0` the matrix is empty and Lean's `0 / 0 = 0` is never read. -/
noncomputable def colMean {n : ℕ} {m : Type} (A : Matrix (Fin n) m ℝ) : m → ℝ :=
  fun j => (∑ i, A i j) / n

/-- Each column centred on its mean, `A - 1 m'` (R's `sweep(A, 2, colMeans(A))`). -/
noncomputable def centre {n : ℕ} {m : Type} (A : Matrix (Fin n) m ℝ) : Matrix (Fin n) m ℝ :=
  A - Matrix.of fun _ j => colMean A j

/-- `restoreRuvFeatureLevels(corrected, original)` on complete data (`:1133-1140`, transposed to runs × features):
each feature (column) of `Z` shifted by `mean(Y_j) - mean(Z_j)`. -/
noncomputable def restore {n p : ℕ} (Z Y : Matrix (Fin n) (Fin p) ℝ) : Matrix (Fin n) (Fin p) ℝ :=
  Z + Matrix.of fun _ j => colMean Y j - colMean Z j

/-- One RUV-III computation: data `Y`, replicate matrix `M`, controls `ctl`, and the solver's `U` (abstract). -/
structure RuvIII (n p s k : ℕ) where
  Y : Matrix (Fin n) (Fin p) ℝ
  M : Matrix (Fin n) (Fin s) ℝ
  ctl : Finset (Fin p)
  U : Matrix (Fin n) (Fin k) ℝ

namespace RuvIII

variable {n p s k : ℕ}

/-- `α = U'Y` (`fullalpha[1:k, ]`). -/
def alpha (r : RuvIII n p s k) : Matrix (Fin k) (Fin p) ℝ := r.Uᵀ * r.Y

/-- The controls' columns of `Y`, `Y[, ctl]`. -/
def Yc (r : RuvIII n p s k) : Matrix (Fin n) r.ctl ℝ := r.Y.submatrix id fun c : r.ctl => (c : Fin p)

/-- `a_c = α[, ctl]`. -/
def ac (r : RuvIII n p s k) : Matrix (Fin k) r.ctl ℝ := r.alpha.submatrix id fun c : r.ctl => (c : Fin p)

/-- `W = Y[, ctl] a_c'(a_c a_c')⁻¹` (`ruv::RUVIII`). -/
noncomputable def W (r : RuvIII n p s k) : Matrix (Fin n) (Fin k) ℝ := r.Yc * r.acᵀ * (r.ac * r.acᵀ)⁻¹

/-- `Y* = Y - Wα` (`newY`). -/
noncomputable def Ystar (r : RuvIII n p s k) : Matrix (Fin n) (Fin p) ℝ := r.Y - r.W * r.alpha

/-- The centred factors `W_c = (Y[, ctl] - 1 m_ctl') a_c'(a_c a_c')⁻¹` (`computeRuvCovariates`, `:1179-1180`). -/
noncomputable def Wc (r : RuvIII n p s k) : Matrix (Fin n) (Fin k) ℝ := centre r.Yc * r.acᵀ * (r.ac * r.acᵀ)⁻¹

/-- `adjusted = Y - W_c α` (`computeRuvCovariates`, `:1182`). -/
noncomputable def adjusted (r : RuvIII n p s k) : Matrix (Fin n) (Fin p) ℝ := r.Y - r.Wc * r.alpha

/-- The per-feature offset amount of R-E1, `1 - 1_ctl' A α_j` with `A = a_c'(a_c a_c')⁻¹`. -/
noncomputable def offsetAmount (r : RuvIII n p s k) (j : Fin p) : ℝ :=
  1 - ∑ l : r.ctl, (r.acᵀ * (r.ac * r.acᵀ)⁻¹ * r.alpha) l j

/-- Every replicate set non-empty. -/
def Hset (r : RuvIII n p s k) : Prop := IsUnit (r.Mᵀ * r.M).det

/-- `a_c a_c'` invertible (mirrors `solve()`'s refusal). -/
def Hac (r : RuvIII n p s k) : Prop := IsUnit (r.ac * r.acᵀ).det

/-- `R_M U = U`. -/
def HU (r : RuvIII n p s k) : Prop := residop r.M * r.U = r.U

/-- `R_M 1 = 0`: every run in exactly one replicate set. -/
def Hone (r : RuvIII n p s k) : Prop := residop r.M *ᵥ (fun _ => (1 : ℝ)) = 0

end RuvIII

open RuvIII

/-- R-A1 (kind E). Correction form and subspace: `Y - Y* = Wα`; every feature's removed vector lies in `col(W)`;
`rank(Y - Y*) ≤ k`. Covers every RUV-III and RUV-III-C route (per target). No hypothesis: the conclusion is meant
for the junk inverse too (then `W = 0` and nothing is removed). -/
def A1_stmt : Prop :=
  ∀ (n p s k : ℕ) (r : RuvIII n p s k),
    r.Y - r.Ystar = r.W * r.alpha ∧
    (∀ j : Fin p, (fun i => (r.Y - r.Ystar) i j) ∈ LinearMap.range (Matrix.mulVecLin r.W)) ∧
    (r.Y - r.Ystar).rank ≤ k

/-- R-A2 (kind H, rank `a_c = k`). Projection identities: `U'W = I_k`; `Π = WU'` is idempotent; `U'Y* = 0` (the
corrected values are orthogonal to the replicate-residual directions `U`, not to `W`); the controls' corrected values
are `Y_ctl (I - a_c'(a_c a_c')⁻¹ a_c)`. Hypothesis `Hac`. Covers all routes. -/
def A2_stmt : Prop :=
  ∀ (n p s k : ℕ) (r : RuvIII n p s k), Hac r →
    r.Uᵀ * r.W = 1 ∧ (r.W * r.Uᵀ) * (r.W * r.Uᵀ) = r.W * r.Uᵀ ∧ r.Uᵀ * r.Ystar = 0 ∧
    r.Ystar.submatrix id (fun c : r.ctl => (c : Fin p)) = r.Yc * (1 - r.acᵀ * (r.ac * r.acᵀ)⁻¹ * r.ac)

/-- R-A3 (kind E). `k = 0` is the identity (empty `W`). No hypothesis. The "S" half of the row (the code's `k = 0`
paths write no correction, P4 and L2 `none`) is a property of the R code, checked by the paired test, not stated
here. -/
def A3_stmt : Prop :=
  ∀ (n p s : ℕ) (r : RuvIII n p s 0), r.Ystar = r.Y

/-- R-B1 (kind H, `R_M U = U`). Per-feature level shifts `Y ↦ Y + 1c'`: `R_M Y` unchanged (hence `U`, a function of
`R_M Y`, which the model keeps fixed), `α` unchanged; `Y*` changes by the constant `c_j - v'α_j` per feature; every
group contrast `λ` (`λ'1 = 0`) of `Y*` is unchanged. Hypotheses `Hset`, `Hac`, `HU`, `Hone`. Covers all RUV-III;
RUV-III-C only when no feature has a missing value (R-J2). -/
def B1_stmt : Prop :=
  ∀ (n p s k : ℕ) (r : RuvIII n p s k) (c : Fin p → ℝ), Hset r → Hac r → HU r → Hone r →
    let r' : RuvIII n p s k := { r with Y := r.Y + Matrix.of fun _ j => c j }
    residop r'.M * r'.Y = residop r.M * r.Y ∧ r'.alpha = r.alpha ∧
    (∀ i j, r'.Ystar i j - r.Ystar i j = c j - ∑ l : r.ctl, c l * (r.acᵀ * (r.ac * r.acᵀ)⁻¹ * r.alpha) l j) ∧
    (∀ lam : Fin n → ℝ, ∑ i, lam i = 0 → lam ᵥ* r'.Ystar = lam ᵥ* r.Ystar)

/-- R-B2 (kind E given B1's hypotheses). After restoring levels (complete data), the output shifts by exactly `c`.
Hypotheses as `B1_stmt`. Covers P1-P3, L2 (restored). -/
def B2_stmt : Prop :=
  ∀ (n p s k : ℕ) (r : RuvIII n p s k) (c : Fin p → ℝ), Hset r → Hac r → HU r → Hone r →
    let r' : RuvIII n p s k := { r with Y := r.Y + Matrix.of fun _ j => c j }
    restore r'.Ystar r'.Y = restore r.Ystar r.Y + Matrix.of fun _ j => c j

/-- R-C1 (kind E). Restoring is the centred formula: on complete data `restoreRuvFeatureLevels(Y - Wα, Y) =
Y - W_c α`; and the second centring of `R/func_omics_ruv_factors.R:173` is the identity (`1'W_c = 0`). Hypothesis
`Hac` (the code's refusal). Covers P1-P3, L2, P5. -/
def C1_stmt : Prop :=
  ∀ (n p s k : ℕ) (r : RuvIII n p s k), Hac r →
    restore r.Ystar r.Y = r.adjusted ∧ centre r.Wc = r.Wc

/-- R-D1 (kind E, `Q` invertible). Basis invariance: replacing `U` by `UQ` gives the same `Y*` (sign flips and
rotations of the eigen-solver do not matter). Hypotheses `IsUnit Q.det`, `Hac`. Covers all. -/
def D1_stmt : Prop :=
  ∀ (n p s k : ℕ) (r : RuvIII n p s k) (Q : Matrix (Fin k) (Fin k) ℝ), IsUnit Q.det → Hac r →
    ({ r with U := r.U * Q } : RuvIII n p s k).Ystar = r.Ystar

/-- R-D2 (kind H, the eigen gap stated as subspace equality). Scaling `Y ↦ sY`: if the solver returns, for `sY`, a
basis `UQ` of the same subspace (true when the `k`-th and `(k+1)`-th eigenvalues differ; with R's `svd`, `Q = I`),
then `Y*(sY) = s Y*(Y)` and `W` is unchanged up to that basis change (`W' Q' = W`; with `Q = I`, `W' = W`).
Hypotheses `s ≠ 0`, `IsUnit Q.det`, `Hac`. Covers all. -/
def D2_stmt : Prop :=
  ∀ (n p s k : ℕ) (r : RuvIII n p s k) (t : ℝ) (Q : Matrix (Fin k) (Fin k) ℝ), t ≠ 0 → IsUnit Q.det → Hac r →
    let r' : RuvIII n p s k := { r with Y := t • r.Y, U := r.U * Q }
    r'.Ystar = t • r.Ystar ∧ r'.W * Qᵀ = r.W

/-- R-E1 (kind H, `R_M U = U`). Offset `Y ↦ Y + c 11'` (every value shifted, a change of log unit): `R_M Y` and `α`
unchanged; `Y*` changes for feature `j` by `c (1 - 1_ctl' A α_j)`, an amount that differs per feature (KIT-188 2.4
item 4 (iii)); after restoring, by exactly `c`; group contrasts unchanged either way. Hypotheses as `B1_stmt`.
Covers the per-feature amount: L1 and T1's export; exactly `c`: P1 to P3, L2, P5, P6. -/
def E1_stmt : Prop :=
  ∀ (n p s k : ℕ) (r : RuvIII n p s k) (c : ℝ), Hset r → Hac r → HU r → Hone r →
    let r' : RuvIII n p s k := { r with Y := r.Y + Matrix.of fun _ _ => c }
    residop r'.M * r'.Y = residop r.M * r.Y ∧ r'.alpha = r.alpha ∧
    (∀ i j, r'.Ystar i j - r.Ystar i j = c * r.offsetAmount j) ∧
    (∀ i j, restore r'.Ystar r'.Y i j - restore r.Ystar r.Y i j = c) ∧
    (∀ lam : Fin n → ℝ, ∑ i, lam i = 0 → lam ᵥ* r'.Ystar = lam ᵥ* r.Ystar)

/-- R-E2 (kind E, an instance). The amounts differ: an instance with 4 runs, 3 features and `k = 1` that satisfies
R-E1's hypotheses has two features with different offset amounts (the proof gives rational entries). -/
def E2_stmt : Prop :=
  ∃ (s : ℕ) (r : RuvIII 4 3 s 1), Hset r ∧ Hac r ∧ HU r ∧ Hone r ∧
    ∃ j j' : Fin 3, r.offsetAmount j ≠ r.offsetAmount j'

/-- R-K1 (i) (kind H, `k = n - sets`). If `U U' = R_M`: `R_M Y* = 0`, the corrected values are constant within every
replicate set (a test whose groups are the replicate sets has residual sum of squares 0). Hypotheses `Hset`, `Hac`,
`U U' = R_M`. Covers all. -/
def K1_stmt : Prop :=
  ∀ (n p s k : ℕ) (r : RuvIII n p s k), Hset r → Hac r → r.U * r.Uᵀ = residop r.M →
    residop r.M * r.Ystar = 0

/-- R-K1 (ii). If also `S = Y_c Y_c'` is invertible (at least as many controls as runs), `Y* = M(M'S⁻¹M)⁻¹M'S⁻¹Y`,
the closed form RUVIIIC computes in its `currentResidualDimensions == k` branch, read at the locked commit c534151f,
`src/RUVIIIC_Varying.cpp:173-184`: `productControls = Y_C Y_C'` (`:182`), `.inverse()` with no rank check (`:183`),
`M_r (M_r' S⁻¹ M_r)⁻¹ M_r' S⁻¹ y` (`:184`), `M_r` the rows of R and the replicate-set columns with a run in R
(`:151-165`); and `ruv::RUVIII` at `k = NULL`. Hypotheses `Hset`; `U U' = R_M`; `IsUnit S.det` (the branch inverts `S` and does not check it, audit
U3); `IsUnit (M'S⁻¹M).det` and `Hac`, which follow from the first three (kept explicit unless a lemma in `Audit.lean`
derives them). -/
def K1_gls_stmt : Prop :=
  ∀ (n p s k : ℕ) (r : RuvIII n p s k), Hset r → r.U * r.Uᵀ = residop r.M →
    IsUnit (r.Yc * r.Ycᵀ).det → IsUnit (r.Mᵀ * (r.Yc * r.Ycᵀ)⁻¹ * r.M).det → Hac r →
      r.Ystar = r.M * (r.Mᵀ * (r.Yc * r.Ycᵀ)⁻¹ * r.M)⁻¹ * (r.Mᵀ * (r.Yc * r.Ycᵀ)⁻¹) * r.Y

/-! ## Witnesses (section 5.1)

Each witness asks for a non-degenerate instance (`0 < k`, at least one more run than sets) that satisfies all the
statement's hypotheses together. -/

/-- Witness for `A2_stmt`, `C1_stmt`. -/
def A2_witness_stmt : Prop :=
  ∃ (n p s k : ℕ) (r : RuvIII n p s k), 0 < k ∧ s < n ∧ Hac r

/-- Witness for `B1_stmt`, `B2_stmt`, `E1_stmt` (the four hypotheses together). -/
def B1_witness_stmt : Prop :=
  ∃ (n p s k : ℕ) (r : RuvIII n p s k), 0 < k ∧ s < n ∧ Hset r ∧ Hac r ∧ HU r ∧ Hone r

/-- Witness for `D1_stmt` and `D2_stmt`: a non-identity `Q` and `t ∉ {0, 1}`. -/
def D_witness_stmt : Prop :=
  ∃ (n p s k : ℕ) (r : RuvIII n p s k) (Q : Matrix (Fin k) (Fin k) ℝ) (t : ℝ),
    0 < k ∧ Hac r ∧ IsUnit Q.det ∧ Q ≠ 1 ∧ t ≠ 0 ∧ t ≠ 1

/-- Witness for `K1_stmt`. -/
def K1_witness_stmt : Prop :=
  ∃ (n p s k : ℕ) (r : RuvIII n p s k), 0 < k ∧ 0 < s ∧ Hset r ∧ Hac r ∧ r.U * r.Uᵀ = residop r.M

/-- Witness for `K1_gls_stmt`. -/
def K1_gls_witness_stmt : Prop :=
  ∃ (n p s k : ℕ) (r : RuvIII n p s k), 0 < k ∧ 0 < s ∧ Hset r ∧ r.U * r.Uᵀ = residop r.M ∧
    IsUnit (r.Yc * r.Ycᵀ).det ∧ IsUnit (r.Mᵀ * (r.Yc * r.Ycᵀ)⁻¹ * r.M).det ∧ Hac r

end MultiScholaRProofs.Ruv
