module

public import MultiScholaRProofs.Spec.RuvDownstream
public import MultiScholaRProofs.Proof.RuvAlgebra
import Mathlib.Tactic

/-!
# Proofs for group R, what consumes the corrected values or the factors (KIT-193 section 7.2), priority 1

Theorems whose types are the frozen statements of `Spec/RuvDownstream.lean`: R-F1 (span lemma) with `F_witness`;
R-G0 (residual space of corrected values) with `G0_witness`; R-G2 (`G2` and the two-case instance `G2_sign`); R-H1
(`H1`, `H1_weighted`) with `H1_witness`.
-/

public section

namespace MultiScholaRProofs.Ruv

open Matrix RuvIII

/-! ## Least squares -/

section LeastSquares

variable {n : ℕ} {m : Type} [Fintype m] [DecidableEq m]

/-- Shifting the response by `Xδ` shifts the coefficients by `δ`. -/
theorem betaHat_shift (X : Matrix (Fin n) m ℝ) (y : Fin n → ℝ) (δ : m → ℝ) (h : IsUnit (Xᵀ * X).det) :
    betaHat X (y - X *ᵥ δ) = betaHat X y - δ := by
  unfold betaHat
  simp only [Matrix.mulVec_sub, Matrix.mulVec_mulVec, Matrix.nonsing_inv_mul _ h, Matrix.one_mulVec]

theorem resid_shift (X : Matrix (Fin n) m ℝ) (y : Fin n → ℝ) (δ : m → ℝ) (h : IsUnit (Xᵀ * X).det) :
    resid X (y - X *ᵥ δ) = resid X y := by
  unfold resid
  rw [betaHat_shift X y δ h, Matrix.mulVec_sub]
  abel

/-- The normal equations: `X'(y - Xβ̂) = 0`. -/
theorem normal_eq (X : Matrix (Fin n) m ℝ) (y : Fin n → ℝ) (h : IsUnit (Xᵀ * X).det) :
    Xᵀ *ᵥ resid X y = 0 := by
  unfold resid betaHat
  simp only [Matrix.mulVec_sub, Matrix.mulVec_mulVec, ← Matrix.mul_assoc, Matrix.mul_nonsing_inv _ h,
    Matrix.one_mul, sub_self]

theorem dot_self_nonneg (v : Fin n → ℝ) : 0 ≤ v ⬝ᵥ v := by
  unfold dotProduct
  exact Finset.sum_nonneg fun i _ => mul_self_nonneg _

/-- Least squares is least: `RSS ≤ ‖y - Xb‖²` for every `b`. -/
theorem rss_le (X : Matrix (Fin n) m ℝ) (y : Fin n → ℝ) (b : m → ℝ) (h : IsUnit (Xᵀ * X).det) :
    rss X y ≤ (y - X *ᵥ b) ⬝ᵥ (y - X *ᵥ b) := by
  set r := resid X y with hr
  set d := betaHat X y - b with hd
  have he : y - X *ᵥ b = r + X *ᵥ d := by
    rw [hr, hd, resid, Matrix.mulVec_sub]
    abel
  have hcross : r ⬝ᵥ (X *ᵥ d) = 0 := by
    rw [Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose, hr, normal_eq X y h, zero_dotProduct]
  have hcross' : (X *ᵥ d) ⬝ᵥ r = 0 := by rw [dotProduct_comm, hcross]
  rw [he, add_dotProduct, dotProduct_add, dotProduct_add, hcross, hcross']
  unfold rss
  rw [← hr]
  linarith [dot_self_nonneg (X *ᵥ d)]

end LeastSquares

/-! ## R-F1 -/

/-- R-F1: the span lemma. -/
theorem F1 : F1_stmt := by
  intro n m _ _ X y δ h
  exact ⟨betaHat_shift X y δ h, resid_shift X y δ h⟩

/-- Witness for `F1_stmt`, `F2_stmt`, `F3_span_stmt`, `F3_fwl_stmt`, `G2_stmt`: `X1 = e_1`, `X2 = e_2` in `ℝ^3`. -/
theorem F_witness : F_witness_stmt := by
  let X1 : Matrix (Fin 3) (Fin 1) ℝ := Matrix.of fun i _ => if i = 0 then 1 else 0
  let X2 : Matrix (Fin 3) (Fin 1) ℝ := Matrix.of fun i _ => if i = 1 then 1 else 0
  have hG : (fromCols X1 X2)ᵀ * fromCols X1 X2 = 1 := by
    ext (a | a) (b | b) <;> fin_cases a <;> fin_cases b <;>
      simp [X1, X2, Matrix.mul_apply]
  have hG1 : X1ᵀ * X1 = 1 := by
    ext a b
    fin_cases a
    fin_cases b
    simp [X1, Matrix.mul_apply]
  refine ⟨3, 1, 1, X1, X2, by norm_num, by norm_num, by norm_num, ?_, ?_⟩
  · rw [hG, Matrix.det_one]
    exact isUnit_one
  · rw [hG1, Matrix.det_one]
    exact isUnit_one

/-! ## R-G0 -/

section G0

variable {n g k : ℕ}

theorem Ut_X_zero {s : ℕ} (M : Matrix (Fin n) (Fin s) ℝ) (U : Matrix (Fin n) (Fin k) ℝ)
    (X : Matrix (Fin n) (Fin g) ℝ) (hU : residop M * U = U) (hX : residop M * X = 0) : Uᵀ * X = 0 := by
  rw [← hU, Matrix.transpose_mul, residop_transpose, Matrix.mul_assoc, hX, Matrix.mul_zero]

theorem Xt_residop (X : Matrix (Fin n) (Fin g) ℝ) (hX : IsUnit (Xᵀ * X).det) : Xᵀ * residop X = 0 := by
  unfold residop
  rw [Matrix.mul_sub, Matrix.mul_one, ← Matrix.mul_assoc, ← Matrix.mul_assoc, Matrix.mul_nonsing_inv _ hX,
    Matrix.one_mul, sub_self]

theorem residop_X (X : Matrix (Fin n) (Fin g) ℝ) (hX : IsUnit (Xᵀ * X).det) : residop X * X = 0 := by
  unfold residop
  rw [Matrix.sub_mul, Matrix.one_mul, Matrix.mul_assoc, Matrix.mul_assoc, Matrix.nonsing_inv_mul _ hX, Matrix.mul_one,
    sub_self]

theorem Ut_residop (X : Matrix (Fin n) (Fin g) ℝ) (U : Matrix (Fin n) (Fin k) ℝ) (hUX : Uᵀ * X = 0) :
    Uᵀ * residop X = Uᵀ := by
  unfold residop
  rw [Matrix.mul_sub, Matrix.mul_one, ← Matrix.mul_assoc, ← Matrix.mul_assoc, hUX, Matrix.zero_mul,
    Matrix.zero_mul, sub_zero]

theorem residop_U (X : Matrix (Fin n) (Fin g) ℝ) (U : Matrix (Fin n) (Fin k) ℝ) (hUX : Uᵀ * X = 0) :
    residop X * U = U := by
  have hXU : Xᵀ * U = 0 := by
    rw [← transpose_transpose (Xᵀ * U), transpose_mul, transpose_transpose, hUX, transpose_zero]
  unfold residop
  rw [Matrix.sub_mul, Matrix.one_mul, Matrix.mul_assoc, hXU, Matrix.mul_zero, sub_zero]

/-- The rank bound: if `U'V = 0`, `U'X = 0`, `U'W = I` and `X'X` is invertible, then the residuals of `V` on `X`
have rank at most `n - g - k` (they lie in the orthogonal complement of `[X, U]`, which has rank `g + k`). -/
theorem rank_bound {p : ℕ} (X : Matrix (Fin n) (Fin g) ℝ) (U W : Matrix (Fin n) (Fin k) ℝ)
    (V : Matrix (Fin n) (Fin p) ℝ) (hX : IsUnit (Xᵀ * X).det) (hUX : Uᵀ * X = 0) (hUW : Uᵀ * W = 1)
    (hUV : Uᵀ * V = 0) : (residop X * V).rank + g + k ≤ n := by
  have hXU : Xᵀ * U = 0 := by
    rw [← transpose_transpose (Xᵀ * U), transpose_mul, transpose_transpose, hUX, transpose_zero]
  have hWU : Wᵀ * U = 1 := by
    rw [← transpose_transpose (Wᵀ * U), transpose_mul, transpose_transpose, hUW, transpose_one]
  have hQZ : (fromCols X U)ᵀ * (residop X * V) = 0 := by
    rw [transpose_fromCols, fromRows_mul, ← Matrix.mul_assoc, ← Matrix.mul_assoc, Xt_residop X hX,
      Ut_residop X U hUX, hUV, Matrix.zero_mul, fromRows_zero]
  have hLQ : fromRows ((Xᵀ * X)⁻¹ * Xᵀ) (Wᵀ * residop X) * fromCols X U = 1 := by
    rw [fromRows_mul_fromCols, ← fromBlocks_one, fromBlocks_inj]
    refine ⟨?_, ?_, ?_, ?_⟩
    · rw [Matrix.mul_assoc, Matrix.nonsing_inv_mul _ hX]
    · rw [Matrix.mul_assoc, hXU, Matrix.mul_zero]
    · rw [Matrix.mul_assoc, residop_X X hX, Matrix.mul_zero]
    · rw [Matrix.mul_assoc, residop_U X U hUX, hWU]
  have hrank1 : g + k ≤ (fromCols X U).rank := by
    have h := Matrix.rank_mul_le_right (fromRows ((Xᵀ * X)⁻¹ * Xᵀ) (Wᵀ * residop X)) (fromCols X U)
    rw [hLQ, Matrix.rank_one, Fintype.card_sum, Fintype.card_fin, Fintype.card_fin] at h
    exact h
  have hsyl := Matrix.rank_mul_ge (fromCols X U)ᵀ (residop X * V)
  rw [hQZ, Matrix.rank_zero, Matrix.rank_transpose, Fintype.card_fin, zero_add] at hsyl
  omega

end G0

/-- R-G0: the residual space of corrected values. -/
theorem G0 : G0_stmt := by
  intro n p s k g r X _ hac hU hX
  have hUX := Ut_X_zero r.M r.U X hU hX
  have hUW := UtW r hac
  have hUY := Ut_Ystar r hac
  refine ⟨fun B => ?_, fun hXX => rank_bound X r.U r.W r.Ystar hXX hUX hUW hUY, fun hone => ?_⟩
  · rw [Matrix.mul_sub, hUY, ← Matrix.mul_assoc, hUX, Matrix.zero_mul, sub_zero]
  · have hU1 : r.Uᵀ *ᵥ (fun _ => (1 : ℝ)) = 0 := by
      have h := hone
      unfold Hone at h
      rw [← hU, transpose_mul, residop_transpose, ← Matrix.mulVec_mulVec, h, Matrix.mulVec_zero]
    have hUcen : r.Uᵀ * centre r.Yc = r.ac := by
      unfold centre
      rw [Matrix.mul_sub, Ut_Yc]
      have h0 : r.Uᵀ * (Matrix.of (fun (_ : Fin n) (l : r.ctl) => colMean r.Yc l) : Matrix (Fin n) r.ctl ℝ) = 0 := by
        ext a l
        have h1 := congrFun hU1 a
        simp only [Matrix.mulVec, dotProduct, mul_one, Pi.zero_apply, Matrix.transpose_apply] at h1
        simp only [Matrix.mul_apply, Matrix.transpose_apply, Matrix.of_apply, Matrix.zero_apply]
        rw [← Finset.sum_mul, h1, zero_mul]
      rw [h0, sub_zero]
    have hUWc : r.Uᵀ * r.Wc = 1 := by
      unfold RuvIII.Wc
      rw [← Matrix.mul_assoc, ← Matrix.mul_assoc, hUcen]
      exact Matrix.mul_nonsing_inv _ hac
    have hUadj : r.Uᵀ * r.adjusted = 0 := by
      unfold RuvIII.adjusted
      rw [Matrix.mul_sub, ← Matrix.mul_assoc, hUWc, Matrix.one_mul]
      unfold RuvIII.alpha
      exact sub_self _
    refine ⟨fun B => ?_, fun hXX => rank_bound X r.U r.W r.adjusted hXX hUX hUW hUadj⟩
    rw [Matrix.mul_sub, hUadj, ← Matrix.mul_assoc, hUX, Matrix.zero_mul, sub_zero]

/-- The witness for R-G0: two runs in one set, one feature (a control), `U = (1, -1)`, `X = M = (1, 1)`. -/
def witnessG0 : RuvIII 2 1 1 1 where
  Y := Matrix.of fun i _ => if i = 0 then 1 else 0
  M := Matrix.of fun _ _ => 1
  ctl := Finset.univ
  U := Matrix.of fun i _ => if i = 0 then 1 else -1

theorem witnessG0_MtM_inv : (witnessG0.Mᵀ * witnessG0.M)⁻¹ = Matrix.of fun _ _ => (1 / 2 : ℝ) := by
  apply Matrix.inv_eq_left_inv
  ext a b
  fin_cases a
  fin_cases b
  simp [witnessG0, Matrix.mul_apply]

theorem witnessG0_residop : residop witnessG0.M =
    Matrix.of fun i j => if i = j then (1 / 2 : ℝ) else -1 / 2 := by
  unfold residop
  rw [witnessG0_MtM_inv]
  ext i j
  fin_cases i <;> fin_cases j <;> simp [witnessG0, Matrix.mul_apply] <;> norm_num

theorem G0_witness : G0_witness_stmt := by
  refine ⟨2, 1, 1, 1, 1, witnessG0, witnessG0.M, by norm_num, by norm_num, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · have h : witnessG0.Mᵀ * witnessG0.M = Matrix.of fun _ _ => (2 : ℝ) := by
      ext a b
      fin_cases a
      fin_cases b
      simp [witnessG0, Matrix.mul_apply]
    unfold Hset
    rw [h, Matrix.det_unique]
    simp
  · have h : witnessG0.ac * witnessG0.acᵀ = 1 := by
      ext a b
      fin_cases a
      fin_cases b
      rw [Matrix.mul_apply]
      simp only [witnessG0, RuvIII.ac, RuvIII.alpha, Matrix.submatrix_apply, Matrix.transpose_apply,
        Matrix.mul_apply, Matrix.of_apply, Matrix.one_apply, id]
      simp
      rfl
    unfold Hac
    rw [h, Matrix.det_one]
    exact isUnit_one
  · unfold HU
    rw [witnessG0_residop]
    ext i j
    fin_cases i <;> fin_cases j <;> simp [witnessG0, Matrix.mul_apply, Fin.sum_univ_two] <;> norm_num
  · unfold Hone
    rw [witnessG0_residop]
    ext i
    fin_cases i <;> simp [Matrix.mulVec, dotProduct, Fin.sum_univ_two] <;> norm_num
  · rw [witnessG0_residop]
    ext i j
    fin_cases i <;> fin_cases j <;> simp [witnessG0, Matrix.mul_apply, Fin.sum_univ_two] <;> norm_num
  · have h : witnessG0.Mᵀ * witnessG0.M = Matrix.of fun _ _ => (2 : ℝ) := by
      ext a b
      fin_cases a
      fin_cases b
      simp [witnessG0, Matrix.mul_apply]
    rw [h, Matrix.det_unique]
    simp

/-! ## R-G2 -/

/-- R-G2: `RSS(y on [X, C]) = RSS(y* on [X, C]) ≤ RSS(y* on X)` when `y - y* ∈ col(C)`. -/
theorem G2 : G2_stmt := by
  intro n g q X C y ystar hQ _ ⟨γ, hγ⟩
  have hys : ystar = y - fromCols X C *ᵥ Sum.elim 0 γ := by
    rw [fromCols_mulVec_sumElim, Matrix.mulVec_zero, zero_add, ← hγ]
    abel
  have heq : rss (fromCols X C) y = rss (fromCols X C) ystar := by
    unfold rss
    rw [hys, resid_shift _ _ _ hQ]
  refine ⟨heq, ?_⟩
  rw [heq]
  have h := rss_le (fromCols X C) ystar (Sum.elim (betaHat X ystar) 0) hQ
  rw [fromCols_mulVec_sumElim, Matrix.mulVec_zero, add_zero] at h
  exact h

set_option linter.unusedSimpArgs false in
-- A style linter only (a simp argument not needed in every one of the `fin_cases` branches); no bearing on the proof.
/-- The instance of R-G2: `X = e_1`, `C = e_2` in `ℝ^4`, `y = e_3`; `y1 = y - C` and `y2 = y`. -/
theorem G2_sign : G2_sign_stmt := by
  let X : Matrix (Fin 4) (Fin 1) ℝ := Matrix.of fun i _ => if i = 0 then 1 else 0
  let C : Matrix (Fin 4) (Fin 1) ℝ := Matrix.of fun i _ => if i = 1 then 1 else 0
  let y : Fin 4 → ℝ := fun i => if i = 2 then 1 else 0
  let y1 : Fin 4 → ℝ := fun i => if i = 2 then 1 else if i = 1 then -1 else 0
  have hG : (fromCols X C)ᵀ * fromCols X C = 1 := by
    ext (a | a) (b | b) <;> fin_cases a <;> fin_cases b <;>
      simp [X, C, Matrix.mul_apply, Fin.sum_univ_four]
  have hG1 : Xᵀ * X = 1 := by
    ext a b
    fin_cases a
    fin_cases b
    simp [X, Matrix.mul_apply, Fin.sum_univ_four]
  have hrssQ : ∀ v : Fin 4 → ℝ, rss (fromCols X C) v = v 2 ^ 2 + v 3 ^ 2 := by
    intro v
    have hp : fromCols X C *ᵥ ((fromCols X C)ᵀ *ᵥ v) =
        fun i => if i = 0 then v 0 else if i = 1 then v 1 else 0 := by
      ext i
      fin_cases i <;>
        simp only [Matrix.mulVec, dotProduct, Fintype.sum_sum_type, Matrix.fromCols_apply_inl,
          Matrix.fromCols_apply_inr, Matrix.transpose_apply, Fin.sum_univ_four, Fin.sum_univ_one, X, C,
          Matrix.of_apply] <;> norm_num
    rw [rss, resid, betaHat, hG, inv_one, Matrix.one_mulVec, hp]
    simp only [dotProduct, Fin.sum_univ_four, Pi.sub_apply]
    norm_num
    ring
  have hrssX : ∀ v : Fin 4 → ℝ, rss X v = v 1 ^ 2 + v 2 ^ 2 + v 3 ^ 2 := by
    intro v
    have hp : X *ᵥ (Xᵀ *ᵥ v) = fun i => if i = 0 then v 0 else 0 := by
      ext i
      fin_cases i <;>
        simp only [Matrix.mulVec, dotProduct, Matrix.transpose_apply, Fin.sum_univ_four, Fin.sum_univ_one, X,
          Matrix.of_apply] <;> norm_num
    rw [rss, resid, betaHat, hG1, inv_one, Matrix.one_mulVec, hp]
    simp only [dotProduct, Fin.sum_univ_four, Pi.sub_apply]
    norm_num
    ring
  have hdet : IsUnit ((fromCols X C)ᵀ * fromCols X C).det := by
    rw [hG, Matrix.det_one]
    exact isUnit_one
  have hdet1 : IsUnit (Xᵀ * X).det := by
    rw [hG1, Matrix.det_one]
    exact isUnit_one
  refine ⟨X, C, y, y1, y, hdet, hdet1, ⟨fun _ => 1, ?_⟩, ⟨0, ?_⟩, ?_, ?_, ?_⟩
  · ext i
    fin_cases i <;> simp [y, y1, C, Matrix.mulVec, dotProduct]
  · ext i
    fin_cases i <;> simp [y, Matrix.mulVec, dotProduct]
  · rw [hrssQ]
    simp [y]
  · unfold varRatio
    rw [hrssX, hrssQ]
    simp [y, y1]
    norm_num
  · unfold varRatio
    rw [hrssX, hrssQ]
    simp [y]
    norm_num

/-! ## R-H1 -/

section H1

variable {g : ℕ}

theorem inv_conj (G T : Matrix (Fin (g + 1)) (Fin (g + 1)) ℝ) (hT : IsUnit T.det) :
    ((T⁻¹)ᵀ * G * T⁻¹)⁻¹ = T * G⁻¹ * Tᵀ := by
  have hTt : IsUnit Tᵀ.det := by rwa [Matrix.det_transpose]
  rw [Matrix.mul_inv_rev, Matrix.mul_inv_rev, Matrix.nonsing_inv_nonsing_inv _ hT, transpose_nonsing_inv,
    Matrix.nonsing_inv_nonsing_inv _ hTt, Matrix.mul_assoc]

theorem TtTinvt (T : Matrix (Fin (g + 1)) (Fin (g + 1)) ℝ) (hT : IsUnit T.det) : Tᵀ * (T⁻¹)ᵀ = 1 := by
  rw [← transpose_mul, Matrix.nonsing_inv_mul _ hT, transpose_one]

theorem quad_entry (T H : Matrix (Fin (g + 1)) (Fin (g + 1)) ℝ) (c : Fin (g + 1) → ℝ) (hc : ∀ j, T 0 j = c j) :
    (T * H * Tᵀ) 0 0 = c ⬝ᵥ (H *ᵥ c) := by
  simp only [Matrix.mul_apply, Matrix.transpose_apply, hc, dotProduct, Matrix.mulVec, Finset.sum_mul,
    Finset.mul_sum]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun j _ => Finset.sum_congr rfl fun l _ => ?_
  ring

theorem first_entry (T : Matrix (Fin (g + 1)) (Fin (g + 1)) ℝ) (c b : Fin (g + 1) → ℝ) (hc : ∀ j, T 0 j = c j) :
    (T *ᵥ b) 0 = c ⬝ᵥ b := by
  simp [Matrix.mulVec, dotProduct, hc]

end H1

/-- R-H1: exact per-comparison standard errors through the reparametrised design. -/
theorem H1 : H1_stmt := by
  intro n g X T c y _ hT hc
  have hG : (X * T⁻¹)ᵀ * (X * T⁻¹) = (T⁻¹)ᵀ * (Xᵀ * X) * T⁻¹ := by
    rw [transpose_mul]
    simp only [Matrix.mul_assoc]
  have hinv : ((X * T⁻¹)ᵀ * (X * T⁻¹))⁻¹ = T * (Xᵀ * X)⁻¹ * Tᵀ := by rw [hG, inv_conj _ _ hT]
  refine ⟨?_, ?_⟩
  · have hb : betaHat (X * T⁻¹) y = T *ᵥ betaHat X y := by
      unfold betaHat
      rw [hinv, transpose_mul]
      simp only [Matrix.mulVec_mulVec]
      congr 1
      simp only [Matrix.mul_assoc]
      rw [← Matrix.mul_assoc Tᵀ, TtTinvt T hT, Matrix.one_mul]
    rw [hb, first_entry T c _ hc]
  · rw [hinv, quad_entry T _ c hc]

/-- R-H1 with fixed weights. -/
theorem H1_weighted : H1_weighted_stmt := by
  intro n g X T c w y _ hT hc
  have hG : (X * T⁻¹)ᵀ * diagonal w * (X * T⁻¹) = (T⁻¹)ᵀ * (Xᵀ * diagonal w * X) * T⁻¹ := by
    rw [transpose_mul]
    simp only [Matrix.mul_assoc]
  have hinv : ((X * T⁻¹)ᵀ * diagonal w * (X * T⁻¹))⁻¹ = T * (Xᵀ * diagonal w * X)⁻¹ * Tᵀ := by
    rw [hG, inv_conj _ _ hT]
  refine ⟨?_, ?_⟩
  · have hb : betaHatW (diagonal w) (X * T⁻¹) y = T *ᵥ betaHatW (diagonal w) X y := by
      unfold betaHatW
      rw [hinv, transpose_mul]
      simp only [Matrix.mulVec_mulVec]
      congr 1
      simp only [Matrix.mul_assoc]
      rw [← Matrix.mul_assoc Tᵀ, TtTinvt T hT, Matrix.one_mul]
    rw [hb, first_entry T c _ hc]
  · rw [hinv, quad_entry T _ c hc]

/-- Witness for `H1_stmt` and `H1_weighted_stmt`: `X = I_2`, `T` with first row `(1, -1)`, weights `(2, 1)`. -/
theorem H1_witness : H1_witness_stmt := by
  let T : Matrix (Fin 2) (Fin 2) ℝ := !![1, -1; 0, 1]
  let w : Fin 2 → ℝ := ![2, 1]
  refine ⟨2, 1, 1, T, ![1, -1], w, by norm_num, ?_, ?_, ?_, ?_, ?_, ⟨0, ?_⟩⟩
  · simp
  · simp [w]
  · simp [T, Matrix.det_fin_two]
  · intro j
    fin_cases j <;> simp [T]
  · intro j
    fin_cases j <;> simp
  · simp [w]

/-! ## R-F2 (Frisch-Waugh-Lovell) -/

section F2

variable {n : ℕ} {m : Type} [Fintype m] [DecidableEq m]

/-- A Gram matrix of a matrix with injective `mulVec` is invertible. -/
theorem gram_unit (A : Matrix (Fin n) m ℝ) (h : Function.Injective A.mulVec) : IsUnit (Aᵀ * A).det := by
  rw [← Matrix.isUnit_iff_isUnit_det, ← Matrix.mulVec_injective_iff_isUnit]
  intro v w hvw
  have h0 : (Aᵀ * A) *ᵥ (v - w) = 0 := by
    rw [Matrix.mulVec_sub, hvw, sub_self]
  have h1 : (A *ᵥ (v - w)) ⬝ᵥ (A *ᵥ (v - w)) = 0 := by
    rw [Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose, Matrix.mulVec_mulVec, h0, zero_dotProduct]
  have h2 : A *ᵥ (v - w) = 0 := dotProduct_self_eq_zero.1 h1
  rw [Matrix.mulVec_sub, sub_eq_zero] at h2
  exact h h2

/-- An invertible Gram matrix makes `mulVec` injective. -/
theorem mulVec_inj_of_gram (A : Matrix (Fin n) m ℝ) (h : IsUnit (Aᵀ * A).det) : Function.Injective A.mulVec := by
  intro v w hvw
  have h1 : (Aᵀ * A) *ᵥ v = (Aᵀ * A) *ᵥ w := by
    rw [← Matrix.mulVec_mulVec, ← Matrix.mulVec_mulVec]
    exact congrArg (Aᵀ *ᵥ ·) hvw
  have h2 := congrArg ((Aᵀ * A)⁻¹ *ᵥ ·) h1
  simp only [Matrix.mulVec_mulVec, Matrix.nonsing_inv_mul _ h, Matrix.one_mulVec] at h2
  exact h2

/-- The least-squares coefficient is the unique solution of the normal equations. -/
theorem betaHat_of_normal (A : Matrix (Fin n) m ℝ) (z : Fin n → ℝ) (b : m → ℝ) (h : IsUnit (Aᵀ * A).det)
    (hn : Aᵀ *ᵥ (z - A *ᵥ b) = 0) : betaHat A z = b := by
  have h1 : Aᵀ *ᵥ z = (Aᵀ * A) *ᵥ b := by
    rw [Matrix.mulVec_sub, sub_eq_zero] at hn
    rw [hn, Matrix.mulVec_mulVec]
  unfold betaHat
  rw [h1, Matrix.mulVec_mulVec, Matrix.nonsing_inv_mul _ h, Matrix.one_mulVec]

end F2

/-- R-F2: Frisch-Waugh-Lovell. -/
theorem F2 : F2_stmt := by
  intro n a b X1 X2 y hQ
  set Q := fromCols X1 X2 with hQdef
  have hQinj := mulVec_inj_of_gram Q hQ
  have h2inj : Function.Injective X2.mulVec := by
    intro v w hvw
    have : Q *ᵥ Sum.elim 0 v = Q *ᵥ Sum.elim 0 w := by
      rw [hQdef, fromCols_mulVec_sumElim, fromCols_mulVec_sumElim, hvw]
    have h := hQinj this
    funext i
    exact congrFun h (Sum.inr i)
  have h2 : IsUnit (X2ᵀ * X2).det := gram_unit X2 h2inj
  set R := residop X2 with hR
  have hRX2 : R * X2 = 0 := residop_X X2 h2
  have hRt : Rᵀ = R := residop_transpose X2
  have h1inj : Function.Injective (R * X1).mulVec := by
    intro v w hvw
    have hd : (R * X1) *ᵥ (v - w) = 0 := by rw [Matrix.mulVec_sub, hvw, sub_self]
    set u := (X2ᵀ * X2)⁻¹ *ᵥ (X2ᵀ *ᵥ (X1 *ᵥ (v - w))) with hu
    have hX1 : X1 *ᵥ (v - w) = X2 *ᵥ u := by
      rw [← Matrix.mulVec_mulVec, hR, residop, Matrix.sub_mulVec, Matrix.one_mulVec, sub_eq_zero] at hd
      rw [hd, hu]
      simp only [Matrix.mulVec_mulVec, Matrix.mul_assoc]
    have hz : Q *ᵥ Sum.elim (v - w) (-u) = Q *ᵥ Sum.elim 0 0 := by
      rw [hQdef, fromCols_mulVec_sumElim, fromCols_mulVec_sumElim, hX1, Matrix.mulVec_neg, Matrix.mulVec_zero,
        Matrix.mulVec_zero, add_zero, add_neg_cancel]
    have h := hQinj hz
    funext i
    have := congrFun h (Sum.inl i)
    simp only [Sum.elim_inl, Pi.sub_apply, Pi.zero_apply, sub_eq_zero] at this
    exact this
  have h1 : IsUnit ((R * X1)ᵀ * (R * X1)).det := gram_unit (R * X1) h1inj
  set β := betaHat Q y with hβ
  set r := resid Q y with hr
  have hnorm : Qᵀ *ᵥ r = 0 := normal_eq Q y hQ
  rw [hQdef, transpose_fromCols, fromRows_mulVec] at hnorm
  have hX1r : X1ᵀ *ᵥ r = 0 := by
    funext i
    exact congrFun hnorm (Sum.inl i)
  have hX2r : X2ᵀ *ᵥ r = 0 := by
    funext i
    exact congrFun hnorm (Sum.inr i)
  have hRr : R *ᵥ r = r := by
    rw [hR, residop, Matrix.sub_mulVec, Matrix.one_mulVec, ← Matrix.mulVec_mulVec, hX2r, Matrix.mulVec_zero,
      sub_zero]
  have hy : y = X1 *ᵥ (β ∘ Sum.inl) + X2 *ᵥ (β ∘ Sum.inr) + r := by
    rw [hr, resid, ← hβ, hQdef, fromCols_mulVec]
    abel
  have hRy : R *ᵥ y - (R * X1) *ᵥ (β ∘ Sum.inl) = r := by
    rw [hy, Matrix.mulVec_add, Matrix.mulVec_add]
    simp only [Matrix.mulVec_mulVec, hRX2, Matrix.zero_mulVec, hRr, add_zero]
    abel
  have hn1 : (R * X1)ᵀ *ᵥ (R *ᵥ y - (R * X1) *ᵥ (β ∘ Sum.inl)) = 0 := by
    rw [hRy, transpose_mul, hRt, ← Matrix.mulVec_mulVec, hRr, hX1r]
  have hb := betaHat_of_normal (R * X1) (R *ᵥ y) (β ∘ Sum.inl) h1 hn1
  refine ⟨?_, ?_⟩
  · rw [hb]
    rfl
  · unfold resid
    rw [hb, hRy]

end MultiScholaRProofs.Ruv
