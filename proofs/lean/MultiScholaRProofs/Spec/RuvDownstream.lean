module

public import MultiScholaRProofs.Spec.RuvAlgebra
public import Mathlib.Data.Matrix.ColumnRowPartitioned
public import Mathlib.LinearAlgebra.Matrix.PosDef
public import Mathlib.Analysis.Real.Sqrt

/-!
# Group R, what consumes the corrected values or the factors (KIT-193 section 7.2)

Least squares from the normal equations, over any finite column index (so `[X, C]` is `Matrix.fromCols X C`, indexed
by `Fin g ⊕ Fin q`). The fits the package calls (limma's `lmFit`, `lm.fit` through `contrasts.fit`, limpa's `dpcDE`
weighted least squares, `duplicateCorrelation`'s generalised least squares with a fixed block correlation) enter
only through these exact formulas; their empirical Bayes and dispersion steps are functions of per-feature inputs
and are not modelled (section 7.3). Every matrix is real and complete.

Sources: KIT-184 r1.3 5.2 (the span argument; its audit `:185-193`), 7.1's inventory rows P1-P6, L1, L2, T1, T2;
`ruvExactContrastsFit` (`R/func_general_ruv_factor_count.R:383-397`, `limma::contrastAsCoef`;
`R/func_prot_da_model_stats.R:424-471`); `dpcDE` weights (`R/func_prot_da_model_methods.R:316-380`);
`duplicateCorrelation` (`R/func_prot_da_model_stats.R:259-261`); edgeR quasi-likelihood (`R/func_trx_da.R:66-77`).
-/

@[expose] public section

namespace MultiScholaRProofs.Ruv

open Matrix RuvIII

variable {n : ℕ} {m : Type} [Fintype m] [DecidableEq m]

/-- Least-squares coefficients `(X'X)⁻¹ X'y`. Used with `IsUnit (Xᵀ * X).det` (full column rank, which the code
requires: limma reports non-estimable coefficients and the RUV designs are refused when rank-deficient). -/
noncomputable def betaHat (X : Matrix (Fin n) m ℝ) (y : Fin n → ℝ) : m → ℝ := (Xᵀ * X)⁻¹ *ᵥ (Xᵀ *ᵥ y)

/-- Residuals `y - X β̂`. -/
noncomputable def resid (X : Matrix (Fin n) m ℝ) (y : Fin n → ℝ) : Fin n → ℝ := y - X *ᵥ betaHat X y

/-- Residual sum of squares. -/
noncomputable def rss (X : Matrix (Fin n) m ℝ) (y : Fin n → ℝ) : ℝ := resid X y ⬝ᵥ resid X y

/-- The ordinary t statistic of contrast `c`: `c'β̂ / √(σ̂² c'(X'X)⁻¹c)`, `σ̂² = RSS / (n - #columns)`. The square
root's argument is `≥ 0` whenever `n > #columns` and `X'X` is positive definite; a zero argument gives Lean's
`x / 0 = 0` (a zero t), and statements that use `tStat` say whether they are meant for that case. -/
noncomputable def tStat (X : Matrix (Fin n) m ℝ) (y : Fin n → ℝ) (c : m → ℝ) : ℝ :=
  (c ⬝ᵥ betaHat X y) /
    Real.sqrt (rss X y / ((n : ℝ) - Fintype.card m) * (c ⬝ᵥ ((Xᵀ * X)⁻¹ *ᵥ c)))

/-- Weighted (generalised) least squares with a fixed weight matrix `Ω`: `(X'ΩX)⁻¹ X'Ωy`. -/
noncomputable def betaHatW (Ω : Matrix (Fin n) (Fin n) ℝ) (X : Matrix (Fin n) m ℝ) (y : Fin n → ℝ) : m → ℝ :=
  (Xᵀ * Ω * X)⁻¹ *ᵥ (Xᵀ *ᵥ (Ω *ᵥ y))

/-- Weighted residuals `y - X β̂_Ω`. -/
noncomputable def residW (Ω : Matrix (Fin n) (Fin n) ℝ) (X : Matrix (Fin n) m ℝ) (y : Fin n → ℝ) : Fin n → ℝ :=
  y - X *ᵥ betaHatW Ω X y

/-- R-F1 (kind E; the briefed target 6's span lemma, KIT-184 5.2). For `X` of full column rank and any `δ`, fitting
`y - Xδ` instead of `y` shifts the coefficients by exactly `δ` and leaves the residual vector (hence the residual df
`n - rank X`) unchanged; so for `δ` supported on the columns not tested (`[1, f, C]` or `[groups, W]`, `δ` on `C` or
`W`) the tested coefficient is the same. Hypothesis `IsUnit (Xᵀ * X).det` (full rank, which the code requires).
Covers P3, L2, P5, I1, and T1/T2's `removeBatchEffect` export. -/
def F1_stmt : Prop :=
  ∀ (n : ℕ) (m : Type) [Fintype m] [DecidableEq m] (X : Matrix (Fin n) m ℝ) (y : Fin n → ℝ) (δ : m → ℝ),
    IsUnit (Xᵀ * X).det →
      betaHat X (y - X *ᵥ δ) = betaHat X y - δ ∧ resid X (y - X *ᵥ δ) = resid X y

/-- R-F2 (kind E; the briefed target 6, Frisch-Waugh-Lovell). For `[X1, X2]` of full column rank: the coefficients
on `X1` of `y` on `[X1, X2]` equal those of `M_2 y` on `M_2 X1` (`M_2 = residop X2`), and the residuals are equal.
Hypothesis: `IsUnit` of the determinant of `[X1, X2]'[X1, X2]` only; that `X2'X2` and `(M_2 X1)'(M_2 X1)` are then
invertible (so `residop X2` and the second fit are the real ones) is to be proved as a lemma, not assumed. -/
def F2_stmt : Prop :=
  ∀ (n a b : ℕ) (X1 : Matrix (Fin n) (Fin a) ℝ) (X2 : Matrix (Fin n) (Fin b) ℝ) (y : Fin n → ℝ),
    IsUnit ((fromCols X1 X2)ᵀ * fromCols X1 X2).det →
      (fun i => betaHat (fromCols X1 X2) y (Sum.inl i)) = betaHat (residop X2 * X1) (residop X2 *ᵥ y) ∧
      resid (fromCols X1 X2) y = resid (residop X2 * X1) (residop X2 *ᵥ y)

/-- R-F3 (kind E, `Ω` fixed), the weighted span lemma: R-F1 with a fixed weight matrix `Ω` (diagonal precision
weights of `dpcDE`; the inverse of a fixed block correlation of `duplicateCorrelation`, treated as known).
Hypothesis `IsUnit (Xᵀ * Ω * X).det`; positive definiteness of `Ω` is not needed here. Covers P5; P1 to P3 and L1, L2
with declared technical replicates. -/
def F3_span_stmt : Prop :=
  ∀ (n : ℕ) (m : Type) [Fintype m] [DecidableEq m] (Ω : Matrix (Fin n) (Fin n) ℝ) (X : Matrix (Fin n) m ℝ)
    (y : Fin n → ℝ) (δ : m → ℝ),
    IsUnit (Xᵀ * Ω * X).det →
      betaHatW Ω X (y - X *ᵥ δ) = betaHatW Ω X y - δ ∧ residW Ω X (y - X *ᵥ δ) = residW Ω X y

/-- R-F3, the weighted Frisch-Waugh-Lovell: with `M_2 = I - X2 (X2'ΩX2)⁻¹ X2'Ω`, the `X1` coefficients of the
weighted fit of `y` on `[X1, X2]` equal those of `M_2 y` on `M_2 X1`, and the weighted residuals are equal.
Hypotheses: `Ω` positive definite (so that the full-rank hypothesis gives `X2'ΩX2` invertible, a lemma);
`IsUnit` of `[X1, X2]'Ω[X1, X2]`'s determinant. -/
def F3_fwl_stmt : Prop :=
  ∀ (n a b : ℕ) (Ω : Matrix (Fin n) (Fin n) ℝ) (X1 : Matrix (Fin n) (Fin a) ℝ) (X2 : Matrix (Fin n) (Fin b) ℝ)
    (y : Fin n → ℝ),
    Ω.PosDef → IsUnit ((fromCols X1 X2)ᵀ * Ω * fromCols X1 X2).det →
      let M2 : Matrix (Fin n) (Fin n) ℝ := 1 - X2 * (X2ᵀ * Ω * X2)⁻¹ * X2ᵀ * Ω
      (fun i => betaHatW Ω (fromCols X1 X2) y (Sum.inl i)) = betaHatW Ω (M2 * X1) (M2 *ᵥ y) ∧
      residW Ω (fromCols X1 X2) y = residW Ω (M2 * X1) (M2 *ᵥ y)

/-- R-G0 (kind H, `R_M U = U`). Residual space of corrected values. If every replicate set is non-empty, `R_M U = U`,
`a_c a_c'` is invertible and the tested design `X` lies in the span of the replicate sets (`R_M X = 0`: each set
inside one group, as with the default formula `~ 0 + group`): `U'(Y* - XB) = 0` for every `B`, so the residuals of
the corrected values from any fit whose fitted values lie in `col(X)` (ordinary, weighted, or with a block
correlation) are orthogonal to `col(U)`; with `X` of full column rank the least-squares residuals have rank at most
`n - g - k`. The same for `Y - W_c α` (= the restored values, R-C1) when every run is in exactly one set.
Hypotheses `Hset`, `Hac`, `HU`, `R_M X = 0`; `IsUnit (Xᵀ * X).det` for the rank bound; `Hone` for the adjusted
values (`U'1 = 0`). Covers P1, P2, P6, L1 with the default formula, on complete data. -/
def G0_stmt : Prop :=
  ∀ (n p s k g : ℕ) (r : RuvIII n p s k) (X : Matrix (Fin n) (Fin g) ℝ),
    Hset r → Hac r → HU r → residop r.M * X = 0 →
      (∀ B : Matrix (Fin g) (Fin p) ℝ, r.Uᵀ * (r.Ystar - X * B) = 0) ∧
      (IsUnit (Xᵀ * X).det → (residop X * r.Ystar).rank + g + k ≤ n) ∧
      (Hone r →
        (∀ B : Matrix (Fin g) (Fin p) ℝ, r.Uᵀ * (r.adjusted - X * B) = 0) ∧
        (IsUnit (Xᵀ * X).det → (residop X * r.adjusted).rank + g + k ≤ n))

/-- R-G1 (kind E; an illustration, not a shipped route). If `C'X = 0` and `y* = (I - C(C'C)⁻¹C')y`, the fit of `y*`
on `X` has the same coefficients and residual sum of squares as `y` on `[X, C]`, and the ordinary t of any contrast
from testing `y*` on `X` is `√((n - g)/(n - g - q))` times the t of the full model (the size of the overstatement in
the orthogonal case only). Hypotheses: `g + q < n` (both residual dfs positive; the divisions are real);
`X'C = 0`; `X'X` and `C'C` invertible. The t identity is meant also where a t is Lean's junk `0` (a zero contrast or
a zero residual sum of squares): it then holds as `0 = 0`. -/
def G1_stmt : Prop :=
  ∀ (n g q : ℕ) (X : Matrix (Fin n) (Fin g) ℝ) (C : Matrix (Fin n) (Fin q) ℝ) (y : Fin n → ℝ),
    g + q < n → Xᵀ * C = 0 → IsUnit (Xᵀ * X).det → IsUnit (Cᵀ * C).det →
      let ystar := y - C *ᵥ ((Cᵀ * C)⁻¹ *ᵥ (Cᵀ *ᵥ y))
      betaHat X ystar = (fun i => betaHat (fromCols X C) y (Sum.inl i)) ∧
      rss X ystar = rss (fromCols X C) y ∧
      ∀ c : Fin g → ℝ,
        tStat X ystar c = Real.sqrt (((n : ℝ) - g) / ((n : ℝ) - g - q)) * tStat (fromCols X C) y (Sum.elim c 0)

/-- R-G2 (kind E). General case: whenever `y - y* ∈ col(C)`, `RSS(y* on X) ≥ RSS(y on [X, C]) = RSS(y* on [X, C])`.
Hypotheses: `[X, C]` of full column rank; `X'X` invertible (follows from the first; kept explicit); `y - y* = Cγ`.
Covers P1 to P3, L1, L2, P5, P6. -/
def G2_stmt : Prop :=
  ∀ (n g q : ℕ) (X : Matrix (Fin n) (Fin g) ℝ) (C : Matrix (Fin n) (Fin q) ℝ) (y ystar : Fin n → ℝ),
    IsUnit ((fromCols X C)ᵀ * fromCols X C).det → IsUnit (Xᵀ * X).det → (∃ γ, y - ystar = C *ᵥ γ) →
      rss (fromCols X C) y = rss (fromCols X C) ystar ∧ rss (fromCols X C) y ≤ rss X ystar

/-- The ratio of the two variance estimates of R-G2: `(RSS(y* on X)/(n - g)) / (RSS(y on [X, C])/(n - g - q))`. -/
noncomputable def varRatio {n g q : ℕ} (X : Matrix (Fin n) (Fin g) ℝ) (C : Matrix (Fin n) (Fin q) ℝ)
    (y ystar : Fin n → ℝ) : ℝ :=
  (rss X ystar / ((n : ℝ) - g)) / (rss (fromCols X C) y / ((n : ℝ) - g - q))

/-- R-G2 (an exact two-case instance). No general direction for the variance estimate's ratio: with the same `X`, `C`
and `y` (4 runs, one tested column, one factor), two corrections `y*` in `y + col(C)` give ratios above and below 1,
with a positive full-model residual sum of squares. -/
def G2_sign_stmt : Prop :=
  ∃ (X : Matrix (Fin 4) (Fin 1) ℝ) (C : Matrix (Fin 4) (Fin 1) ℝ) (y y1 y2 : Fin 4 → ℝ),
    IsUnit ((fromCols X C)ᵀ * fromCols X C).det ∧ IsUnit (Xᵀ * X).det ∧
    (∃ γ, y - y1 = C *ᵥ γ) ∧ (∃ γ, y - y2 = C *ᵥ γ) ∧ 0 < rss (fromCols X C) y ∧
    1 < varRatio X C y y1 ∧ varRatio X C y y2 < 1

/-- R-H1 (kind E). Exact per-comparison standard errors. For `X` of full column rank and contrast `c`, the
reparametrised design `X T⁻¹` (`T` invertible with first row `c'`, as `limma::contrastAsCoef` builds) has first
coefficient `c'β̂` and first unscaled variance `c'(X'X)⁻¹c`. `X` is any full-rank matrix, so on a feature's observed
rows `X_o` (of full rank) the refit gives `c'β̂_o` and `c'(X_o'X_o)⁻¹c` exactly, which `contrasts.fit` approximates
when rows are missing (decision M29). Hypotheses: `IsUnit (Xᵀ * X).det`; `IsUnit T.det`; `T`'s first row is `c`.
Covers P3, P5, L2. -/
def H1_stmt : Prop :=
  ∀ (n g : ℕ) (X : Matrix (Fin n) (Fin (g + 1)) ℝ) (T : Matrix (Fin (g + 1)) (Fin (g + 1)) ℝ)
    (c : Fin (g + 1) → ℝ) (y : Fin n → ℝ),
    IsUnit (Xᵀ * X).det → IsUnit T.det → (∀ j, T 0 j = c j) →
      betaHat (X * T⁻¹) y 0 = c ⬝ᵥ betaHat X y ∧
      ((X * T⁻¹)ᵀ * (X * T⁻¹))⁻¹ 0 0 = c ⬝ᵥ ((Xᵀ * X)⁻¹ *ᵥ c)

/-- R-H1 with fixed weights `w` on the rows (limpa's DPC precision weights, P5): the same with `X'X` replaced by
`X' diag(w) X`. Hypotheses: `IsUnit (Xᵀ * diagonal w * X).det` (positivity of `w` is not needed for the identity);
`IsUnit T.det`; `T`'s first row is `c`. -/
def H1_weighted_stmt : Prop :=
  ∀ (n g : ℕ) (X : Matrix (Fin n) (Fin (g + 1)) ℝ) (T : Matrix (Fin (g + 1)) (Fin (g + 1)) ℝ)
    (c : Fin (g + 1) → ℝ) (w y : Fin n → ℝ),
    IsUnit (Xᵀ * diagonal w * X).det → IsUnit T.det → (∀ j, T 0 j = c j) →
      betaHatW (diagonal w) (X * T⁻¹) y 0 = c ⬝ᵥ betaHatW (diagonal w) X y ∧
      ((X * T⁻¹)ᵀ * diagonal w * (X * T⁻¹))⁻¹ 0 0 = c ⬝ᵥ ((Xᵀ * diagonal w * X)⁻¹ *ᵥ c)

/-- R-M1 (kind E, abstract). edgeR's design enters through its span: `[X, W]` and `[X, WQ]` (`Q` invertible) have the
same column space and rank (so the residual df `n - rank` is the same), so any fit defined as a function of the set
of linear predictors `{Xβ}` (an optimum of any criterion over it) is the same; and every predictor of the first is
the predictor of the second with the `X` coefficients unchanged (a general `Q` changes what a `W` coefficient
means, not an `X` one). The estimator is abstract; RUVs' and RUVg's own algebra is stated after the gate reads
RUVSeq (2.3). Hypothesis `IsUnit Q.det`. Covers T1, T2. -/
def M1_stmt : Prop :=
  ∀ (n a k : ℕ) (X : Matrix (Fin n) (Fin a) ℝ) (W : Matrix (Fin n) (Fin k) ℝ) (Q : Matrix (Fin k) (Fin k) ℝ),
    IsUnit Q.det →
      LinearMap.range (fromCols X W).mulVecLin = LinearMap.range (fromCols X (W * Q)).mulVecLin ∧
      (fromCols X W).rank = (fromCols X (W * Q)).rank ∧
      ∀ (βX : Fin a → ℝ) (βW : Fin k → ℝ),
        fromCols X W *ᵥ Sum.elim βX βW = fromCols X (W * Q) *ᵥ Sum.elim βX (Q⁻¹ *ᵥ βW)

/-! ## Witnesses (section 5.1) -/

/-- Witness for `F1_stmt`, `F2_stmt`, `F3_span_stmt` (with `Ω = 1`), `F3_fwl_stmt`, `G2_stmt`: a full-rank
`[X1, X2]` with both blocks non-empty and more runs than columns. -/
def F_witness_stmt : Prop :=
  ∃ (n a b : ℕ) (X1 : Matrix (Fin n) (Fin a) ℝ) (X2 : Matrix (Fin n) (Fin b) ℝ),
    0 < a ∧ 0 < b ∧ a + b < n ∧ IsUnit ((fromCols X1 X2)ᵀ * fromCols X1 X2).det ∧ IsUnit (X1ᵀ * X1).det

/-- Witness for `F3_fwl_stmt`: a positive definite, non-identity `Ω` with a full-rank weighted Gram matrix. -/
def F3_fwl_witness_stmt : Prop :=
  ∃ (n a b : ℕ) (Ω : Matrix (Fin n) (Fin n) ℝ) (X1 : Matrix (Fin n) (Fin a) ℝ) (X2 : Matrix (Fin n) (Fin b) ℝ),
    0 < a ∧ 0 < b ∧ a + b < n ∧ Ω.PosDef ∧ Ω ≠ 1 ∧ IsUnit ((fromCols X1 X2)ᵀ * Ω * fromCols X1 X2).det

/-- Witness for `G0_stmt`: all its hypotheses together (including the optional ones), `0 < k`, `0 < g`. -/
def G0_witness_stmt : Prop :=
  ∃ (n p s k g : ℕ) (r : RuvIII n p s k) (X : Matrix (Fin n) (Fin g) ℝ),
    0 < k ∧ 0 < g ∧ Hset r ∧ Hac r ∧ HU r ∧ Hone r ∧ residop r.M * X = 0 ∧ IsUnit (Xᵀ * X).det

/-- Witness for `G1_stmt`. -/
def G1_witness_stmt : Prop :=
  ∃ (n g q : ℕ) (X : Matrix (Fin n) (Fin g) ℝ) (C : Matrix (Fin n) (Fin q) ℝ),
    0 < g ∧ 0 < q ∧ g + q < n ∧ Xᵀ * C = 0 ∧ IsUnit (Xᵀ * X).det ∧ IsUnit (Cᵀ * C).det

/-- Witness for `H1_stmt` and `H1_weighted_stmt`: a full-rank `X`, a contrast that is not a coordinate, an invertible
`T` with that first row, and weights other than 1. -/
def H1_witness_stmt : Prop :=
  ∃ (n g : ℕ) (X : Matrix (Fin n) (Fin (g + 1)) ℝ) (T : Matrix (Fin (g + 1)) (Fin (g + 1)) ℝ)
    (c : Fin (g + 1) → ℝ) (w : Fin n → ℝ),
    0 < g ∧ IsUnit (Xᵀ * X).det ∧ IsUnit (Xᵀ * diagonal w * X).det ∧ IsUnit T.det ∧ (∀ j, T 0 j = c j) ∧
    (∀ j, c j ≠ 0) ∧ (∃ i, w i ≠ 1)

/-- Witness for `M1_stmt`: a non-identity invertible `Q` with `0 < k`. -/
def M1_witness_stmt : Prop :=
  ∃ (k : ℕ) (Q : Matrix (Fin k) (Fin k) ℝ), 0 < k ∧ IsUnit Q.det ∧ Q ≠ 1

end MultiScholaRProofs.Ruv
