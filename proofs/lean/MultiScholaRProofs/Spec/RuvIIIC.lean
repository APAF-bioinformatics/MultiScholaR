module

public import MultiScholaRProofs.Spec.RuvAlgebra
public import Mathlib.Analysis.Real.Sqrt
public import Mathlib.Data.Fin.Tuple.Basic

/-!
# Group R, RUV-III-C with missing values (KIT-193 section 7.2: R-J1, R-J2, R-L1, R-L2)

RUV-III-C (RUVIIIC 1.0.19, KIT-185 2.1, read there): per target feature `j`, rows `R_j` (its observed runs), `P_j`
the projection off the replicate sets on `R_j`, `U_j` the `k` leading eigenvectors of `P_j Y0[R_j,] Y0[R_j,]' P_j`
with `Y0` the matrix with every missing value set to 0, `a = U_j' Y0[R_j,]`, `C_j` the controls observed on all of
`R_j`, `A_j = a_C'(a_C a_C')⁻¹`, `W_j = Y[R_j, C_j] A_j`, `y*_j = y_j - W_j a_j` on `R_j`; missing runs stay missing.
Called at `R/func_prot_s4_norm_methods.R:248-254` (P1-P3) and `R/func_lipid_s4_ruv_methods.R:691-697` (L1).

Model. Row restriction to `R` is the diagonal 0/1 mask `mask R` on full-size matrices. The projection `P_R` and the
eigen-solver are abstract fields (how RUVIIIC builds `P_R` when a set has no run in `R` was not read; the statements
below need `P_R` only through the solver's input); the solver is a function of its input matrix only. The
closed-form branch (residual dimension on `R_j` equal to `k`; R-K1 (ii)) and the all-missing branch (fewer than `k`
dimensions or controls; R-J3) are not in this model.
-/

@[expose] public section

namespace MultiScholaRProofs.Ruv

open Matrix Finset

/-- One RUV-III-C computation: values `Y` (entries at missing positions are ignored: zero-filled below), the observed
pattern `obs`, the replicate matrix `M`, the controls, the code's projection off the replicate sets on a row set
(abstract), and the eigen-solver (abstract, a function of its input). -/
structure RuvIIIC (n p s k : ℕ) where
  Y : Matrix (Fin n) (Fin p) ℝ
  obs : Fin n → Fin p → Bool
  M : Matrix (Fin n) (Fin s) ℝ
  ctl : Finset (Fin p)
  proj : Finset (Fin n) → Matrix (Fin n) (Fin n) ℝ
  solver : Matrix (Fin n) (Fin n) ℝ → Matrix (Fin n) (Fin k) ℝ

/-- The diagonal 0/1 matrix that keeps the rows in `R`. -/
def mask {n : ℕ} (R : Finset (Fin n)) : Matrix (Fin n) (Fin n) ℝ := diagonal fun i => if i ∈ R then 1 else 0

namespace RuvIIIC

variable {n p s k : ℕ}

/-- `Y0`: every missing value set to 0. -/
def Y0 (r : RuvIIIC n p s k) : Matrix (Fin n) (Fin p) ℝ := Matrix.of fun i j => if r.obs i j then r.Y i j else 0

/-- `R_j`, the observed runs of feature `j`. -/
def rows (r : RuvIIIC n p s k) (j : Fin p) : Finset (Fin n) := univ.filter fun i => r.obs i j = true

/-- `U_R`: the solver applied to `P_R Y0[R,] Y0[R,]' P_R`. -/
def Uon (r : RuvIIIC n p s k) (R : Finset (Fin n)) : Matrix (Fin n) (Fin k) ℝ :=
  r.solver (r.proj R * mask R * r.Y0 * r.Y0ᵀ * mask R * r.proj R)

/-- `C_R`, the controls observed on every run of `R`. -/
def ctlOn (r : RuvIIIC n p s k) (R : Finset (Fin n)) : Finset (Fin p) :=
  r.ctl.filter fun c => ∀ i ∈ R, r.obs i c = true

/-- `a = U_R' Y0[R,]` (all features). -/
def aOn (r : RuvIIIC n p s k) (R : Finset (Fin n)) : Matrix (Fin k) (Fin p) ℝ := (r.Uon R)ᵀ * mask R * r.Y0

/-- `a_C`, the columns of `a` for `C_R`. -/
def acOn (r : RuvIIIC n p s k) (R : Finset (Fin n)) : Matrix (Fin k) (r.ctlOn R) ℝ :=
  (r.aOn R).submatrix id fun c : r.ctlOn R => (c : Fin p)

/-- `W_R = Y[R, C_R] a_C'(a_C a_C')⁻¹` (rows outside `R` zero). -/
noncomputable def WOn (r : RuvIIIC n p s k) (R : Finset (Fin n)) : Matrix (Fin n) (Fin k) ℝ :=
  mask R * r.Y0.submatrix id (fun c : r.ctlOn R => (c : Fin p)) * (r.acOn R)ᵀ * (r.acOn R * (r.acOn R)ᵀ)⁻¹

/-- `y*_j`: on `R_j`, `y_j - W_j a_j`; missing runs stay missing. -/
noncomputable def ystar (r : RuvIIIC n p s k) (j : Fin p) : Fin n → Option ℝ :=
  fun i => if r.obs i j then some (r.Y i j - (r.WOn (r.rows j) * r.aOn (r.rows j)) i j) else none

/-- The RUV-III computation the complete features share: `Y0`, `M`, the controls complete in every run, and the
solver's output on all runs. -/
def complete (r : RuvIIIC n p s k) : RuvIII n p s k := ⟨r.Y0, r.M, r.ctlOn univ, r.Uon univ⟩

end RuvIIIC

/-- The mean of a vector (used with a non-empty index). -/
noncomputable def mean {m : ℕ} (v : Fin m → ℝ) : ℝ := (∑ i, v i) / m

/-- R-J1 (kinds S + E). RUV-III-C on complete features: for every two features observed in every run, `W` is the
same (so all such features share one `W`), and each one's corrected values are those of the single RUV-III
computation `complete` (so R-A to R-E apply to them with that `W`); a feature with a missing run has its own
`W_j` on `R_j` (by the definition). Hypotheses: features `j`, `j'` observed in every run. Covers P1-P3, L1, L2. -/
def J1_stmt : Prop :=
  ∀ (n p s k : ℕ) (r : RuvIIIC n p s k) (j j' : Fin p), (∀ i, r.obs i j = true) → (∀ i, r.obs i j' = true) →
    r.WOn (r.rows j) = r.WOn (r.rows j') ∧
    ∀ i, r.ystar j i = some (r.complete.Ystar i j)

/-- R-J2 (kind E). What a zero-filled value contributes: on a replicate set of `m + 1` runs where a feature has `m`
observed values `x` and one zero-filled run, the sum of squared deviations from the set's mean is the observed
spread plus `(m/(m+1)) x̄²`, a term of the size of the squared level (KIT-185 2.1's `μ²(1 - 1/m)` is the case of
constant `x`, with `m` there the set size). Hypothesis `0 < m` (the observed mean is a real quotient). -/
def J2_stmt : Prop :=
  ∀ (m : ℕ) (x : Fin m → ℝ), 0 < m →
    ∑ i, ((Fin.snoc x 0 : Fin (m + 1) → ℝ) i - mean (Fin.snoc x 0 : Fin (m + 1) → ℝ)) ^ 2 =
      ∑ i, (x i - mean x) ^ 2 + ((m : ℝ) / (m + 1)) * mean x ^ 2

/-- R-J2, the link to the direction matrix: on one replicate set (an all-ones `M` column), a feature's contribution
to `tr(P Y0 Y0' P)`, `‖P v‖²`, is its sum of squared deviations from the set's mean. No hypothesis (`m + 1 > 0`). -/
def J2_contribution_stmt : Prop :=
  ∀ (m : ℕ) (v : Fin (m + 1) → ℝ),
    (residop (Matrix.of fun (_ : Fin (m + 1)) (_ : Fin 1) => (1 : ℝ)) *ᵥ v) ⬝ᵥ
      (residop (Matrix.of fun (_ : Fin (m + 1)) (_ : Fin 1) => (1 : ℝ)) *ᵥ v) =
      ∑ i, (v i - mean v) ^ 2

/-- The corrected target of R-L1, over its rows `R` and controls `C`: `y - Y_C a_C'(a_C a_C')⁻¹ a`. -/
noncomputable def ystarT {R C : Type} [Fintype R] [Fintype C] {k : ℕ} (y : R → ℝ) (YC : Matrix R C ℝ)
    (aC : Matrix (Fin k) C ℝ) (a : Fin k → ℝ) : R → ℝ :=
  y - (YC * aCᵀ * (aC * aCᵀ)⁻¹) *ᵥ a

/-- KIT-188's fold-change error of a contrast `λ` on target `y` with direction matrix `U` (`a = U'y`,
`a_C = U'Y_C`): `e = -d_C' A a` with `d_C = Y_C' λ`, `A = a_C'(a_C a_C')⁻¹`. -/
noncomputable def foldErr {R C : Type} [Fintype R] [Fintype C] {k : ℕ} (U : Matrix R (Fin k) ℝ) (y : R → ℝ)
    (YC : Matrix R C ℝ) (lam : R → ℝ) : ℝ :=
  -((YCᵀ *ᵥ lam) ⬝ᵥ ((((Uᵀ * YC)ᵀ * ((Uᵀ * YC) * (Uᵀ * YC)ᵀ)⁻¹) *ᵥ (Uᵀ *ᵥ y))))

/-- R-L1 (kind E, identity). For a contrast `λ` on `R_j`: `λ'y*_j = λ'y_j - d_C' A_j a_j`, `d_C = Y[R_j, C_j]'λ`; so
the fold-change error is `e_j = -d_C' A_j a_j`. No hypothesis: the identity holds for any `λ` and is meant for the
junk inverse too (then `A_j = 0` and `e_j = 0`). Covers P1, P2, L1 and every RUV-III route. -/
def L1_stmt : Prop :=
  ∀ (R C : Type) [Fintype R] [Fintype C] (k : ℕ) (y : R → ℝ) (YC : Matrix R C ℝ) (aC : Matrix (Fin k) C ℝ)
    (a : Fin k → ℝ) (lam : R → ℝ),
    lam ⬝ᵥ ystarT y YC aC a = lam ⬝ᵥ y - (YCᵀ *ᵥ lam) ⬝ᵥ ((aCᵀ * (aC * aCᵀ)⁻¹) *ᵥ a)

/-- R-L1 (kind H, `U_j'1 = 0`), invariance: for a fixed `U_j` with `U_j'1 = 0` and a contrast (`λ'1 = 0`), `e_j` is
unchanged when a constant is added to feature `j` or to one control column. (`U_j` itself is unchanged by such
constants only when no feature has a missing value on `R_j`, R-J2; not claimed here.) Hypotheses `U'1 = 0`,
`∑ λ = 0`; meant for the junk inverse too (`a_C` does not change). -/
def L1_invariance_stmt : Prop :=
  ∀ (R C : Type) [Fintype R] [Fintype C] [DecidableEq C] (k : ℕ) (U : Matrix R (Fin k) ℝ) (y : R → ℝ)
    (YC : Matrix R C ℝ) (lam : R → ℝ) (t : ℝ) (l : C),
    Uᵀ *ᵥ (fun _ => (1 : ℝ)) = 0 → ∑ i, lam i = 0 →
      foldErr U (y + fun _ => t) YC lam = foldErr U y YC lam ∧
      foldErr U y (YC + Matrix.of fun _ c => if c = l then t else 0) lam = foldErr U y YC lam

/-- R-L2 (kind E). Its bound: if `κ > 0` and `x'(a_C a_C')x ≥ κ‖x‖²` for all `x`, then
`|d_C' A a| ≤ ‖d_C‖ ‖a‖ / √κ` (`√κ` is `σ_min(a_C)` at the best `κ`; stated without singular values). Hypotheses
`0 < κ` and the quadratic bound (which make `a_C a_C'` invertible and the square roots real). Priority 3. -/
def L2_stmt : Prop :=
  ∀ (C : Type) [Fintype C] (k : ℕ) (aC : Matrix (Fin k) C ℝ) (dC : C → ℝ) (a : Fin k → ℝ) (κ : ℝ),
    0 < κ → (∀ x : Fin k → ℝ, κ * (x ⬝ᵥ x) ≤ x ⬝ᵥ ((aC * aCᵀ) *ᵥ x)) →
      |dC ⬝ᵥ ((aCᵀ * (aC * aCᵀ)⁻¹) *ᵥ a)| ≤ Real.sqrt (dC ⬝ᵥ dC) * Real.sqrt (a ⬝ᵥ a) / Real.sqrt κ

/-! ## Witnesses (section 5.1) -/

/-- Witness for `J1_stmt`: two distinct complete features and a third with a missing run. -/
def J1_witness_stmt : Prop :=
  ∃ (n p s k : ℕ) (r : RuvIIIC n p s k) (j j' j'' : Fin p),
    j ≠ j' ∧ (∀ i, r.obs i j = true) ∧ (∀ i, r.obs i j' = true) ∧ ∃ i, r.obs i j'' = false

/-- Witness for `J2_stmt`: a non-zero observed mean, so the extra term is positive. -/
def J2_witness_stmt : Prop := ∃ (m : ℕ) (x : Fin m → ℝ), 0 < m ∧ mean x ≠ 0

/-- Witness for `L1_invariance_stmt`: `U'1 = 0` with `U ≠ 0`, a non-zero contrast. -/
def L1_invariance_witness_stmt : Prop :=
  ∃ (U : Matrix (Fin 4) (Fin 1) ℝ) (lam : Fin 4 → ℝ),
    Uᵀ *ᵥ (fun _ => (1 : ℝ)) = 0 ∧ U ≠ 0 ∧ ∑ i, lam i = 0 ∧ lam ≠ 0

/-- Witness for `L2_stmt`: an `a_C` with `k = 2` and a positive `κ` that satisfies the bound. -/
def L2_witness_stmt : Prop :=
  ∃ (aC : Matrix (Fin 2) (Fin 3) ℝ) (κ : ℝ), 0 < κ ∧ ∀ x : Fin 2 → ℝ, κ * (x ⬝ᵥ x) ≤ x ⬝ᵥ ((aC * aCᵀ) *ᵥ x)

end MultiScholaRProofs.Ruv
