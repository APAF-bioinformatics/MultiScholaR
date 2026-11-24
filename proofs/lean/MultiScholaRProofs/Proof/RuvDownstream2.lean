module

public import MultiScholaRProofs.Spec.RuvDownstream
public import MultiScholaRProofs.Proof.RuvDownstream
import Mathlib.Tactic
import Mathlib.Algebra.Order.Star.Real

/-!
# Proofs for group R, what consumes the corrected values or the factors (KIT-193 section 7.2), priorities 2 to 4

Theorems whose types are the frozen statements of `Spec/RuvDownstream.lean`: R-G1 (the orthogonal illustration) with
`G1_witness`; R-F3 (`F3_span`, `F3_fwl`) with `F3_fwl_witness`; R-M1 with `M1_witness`.
-/

public section

namespace MultiScholaRProofs.Ruv

open Matrix

/-! ## R-G1 -/

/-- The real identity behind R-G1's t ratio, for any `N`, `ρ`, `v` and positive `a`, `b` (a non-positive square-root
argument gives Lean's `x / 0 = 0` on both sides). -/
theorem t_ratio (N ρ v a b : ℝ) (ha : 0 < a) (hb : 0 < b) :
    N / Real.sqrt (ρ / a * v) = Real.sqrt (a / b) * (N / Real.sqrt (ρ / b * v)) := by
  have h1 : ρ / a * v = (ρ / b * v) * (b / a) := by
    field_simp
  rw [h1, Real.sqrt_mul' _ (div_pos hb ha).le, div_mul_eq_div_div, div_eq_inv_mul (N / _), ← Real.sqrt_inv,
    inv_div]

/-- R-G1: the orthogonal illustration. -/
theorem G1 : G1_stmt := by
  intro n g q X C y hdf hXC hX hC ystar
  have hCX : Cᵀ * X = 0 := by
    rw [← transpose_transpose (Cᵀ * X), transpose_mul, transpose_transpose, hXC, transpose_zero]
  have hgram : (fromCols X C)ᵀ * fromCols X C = fromBlocks (Xᵀ * X) 0 0 (Cᵀ * C) := by
    rw [transpose_fromCols, fromRows_mul_fromCols, hXC, hCX]
  have hQ : IsUnit ((fromCols X C)ᵀ * fromCols X C).det := by
    rw [hgram, Matrix.det_fromBlocks_zero₂₁]
    exact hX.mul hC
  have hginv : ((fromCols X C)ᵀ * fromCols X C)⁻¹ = fromBlocks (Xᵀ * X)⁻¹ 0 0 (Cᵀ * C)⁻¹ := by
    rw [hgram]
    apply Matrix.inv_eq_left_inv
    rw [fromBlocks_multiply]
    simp only [Matrix.mul_zero, Matrix.zero_mul, add_zero, zero_add, Matrix.nonsing_inv_mul _ hX,
      Matrix.nonsing_inv_mul _ hC, fromBlocks_one]
  have hys : ystar = residop C *ᵥ y := by
    show y - C *ᵥ ((Cᵀ * C)⁻¹ *ᵥ (Cᵀ *ᵥ y)) = residop C *ᵥ y
    rw [residop, Matrix.sub_mulVec, Matrix.one_mulVec]
    simp only [Matrix.mulVec_mulVec, Matrix.mul_assoc]
  have hRX : residop C * X = X := by
    rw [residop, Matrix.sub_mul, Matrix.one_mul, Matrix.mul_assoc, hCX, Matrix.mul_zero, sub_zero]
  obtain ⟨hb, hr⟩ := F2 n g q X C y hQ
  rw [hRX, ← hys] at hb hr
  have hrss : rss X ystar = rss (fromCols X C) y := by
    unfold rss
    rw [hr]
  refine ⟨hb.symm, hrss, fun c => ?_⟩
  have hnum : Sum.elim c 0 ⬝ᵥ betaHat (fromCols X C) y = c ⬝ᵥ betaHat X ystar := by
    rw [← Sum.elim_comp_inl_inr (betaHat (fromCols X C) y), sumElim_dotProduct_sumElim, zero_dotProduct,
      add_zero]
    congr 1
  have hquad : Sum.elim c 0 ⬝ᵥ (((fromCols X C)ᵀ * fromCols X C)⁻¹ *ᵥ Sum.elim c 0) =
      c ⬝ᵥ ((Xᵀ * X)⁻¹ *ᵥ c) := by
    rw [hginv, fromBlocks_mulVec]
    simp only [Sum.elim_comp_inl, Sum.elim_comp_inr, Matrix.mulVec_zero, Matrix.zero_mulVec, add_zero,
      sumElim_dotProduct_sumElim, zero_dotProduct]
  have hcard : ((n : ℝ) - (Fintype.card (Fin g ⊕ Fin q) : ℕ)) = (n : ℝ) - g - q := by
    rw [Fintype.card_sum, Fintype.card_fin, Fintype.card_fin]
    push_cast
    ring
  unfold tStat
  rw [hnum, hquad, hcard, Fintype.card_fin, hrss]
  have ha : 0 < (n : ℝ) - g := by
    have : g < n := by omega
    have : (g : ℝ) < n := by exact_mod_cast this
    linarith
  have hb' : 0 < (n : ℝ) - g - q := by
    have : ((g + q : ℕ) : ℝ) < n := by exact_mod_cast hdf
    push_cast at this
    linarith
  exact t_ratio _ _ _ _ _ ha hb'

/-- Witness for `G1_stmt`: `X = e_1`, `C = e_2` in `ℝ^3`. -/
theorem G1_witness : G1_witness_stmt := by
  let X : Matrix (Fin 3) (Fin 1) ℝ := Matrix.of fun i _ => if i = 0 then 1 else 0
  let C : Matrix (Fin 3) (Fin 1) ℝ := Matrix.of fun i _ => if i = 1 then 1 else 0
  have hXX : Xᵀ * X = 1 := by
    ext a b
    fin_cases a
    fin_cases b
    simp [X, Matrix.mul_apply]
  have hCC : Cᵀ * C = 1 := by
    ext a b
    fin_cases a
    fin_cases b
    simp [C, Matrix.mul_apply]
  have hXC : Xᵀ * C = 0 := by
    ext a b
    fin_cases a
    fin_cases b
    simp [X, C, Matrix.mul_apply]
  refine ⟨3, 1, 1, X, C, by norm_num, by norm_num, by norm_num, hXC, ?_, ?_⟩
  · rw [hXX, Matrix.det_one]
    exact isUnit_one
  · rw [hCC, Matrix.det_one]
    exact isUnit_one

/-! ## R-F3 (weighted least squares with a fixed weight matrix) -/

section F3

variable {n : ℕ} {m : Type} [Fintype m] [DecidableEq m]

/-- R-F3, the weighted span lemma. -/
theorem F3_span : F3_span_stmt := by
  intro n m _ _ Ω X y δ h
  have hb : betaHatW Ω X (y - X *ᵥ δ) = betaHatW Ω X y - δ := by
    unfold betaHatW
    rw [Matrix.mulVec_sub, Matrix.mulVec_sub, Matrix.mulVec_sub]
    simp only [Matrix.mulVec_mulVec, ← Matrix.mul_assoc, Matrix.nonsing_inv_mul _ h, Matrix.one_mulVec]
  refine ⟨hb, ?_⟩
  unfold residW
  rw [hb, Matrix.mulVec_sub]
  abel

/-- With `Ω` positive definite, an injective `A` has an invertible weighted Gram matrix `A'ΩA`. -/
theorem gramW_unit (Ω : Matrix (Fin n) (Fin n) ℝ) (hΩ : Ω.PosDef) (A : Matrix (Fin n) m ℝ)
    (h : Function.Injective A.mulVec) : IsUnit (Aᵀ * Ω * A).det := by
  have h1 := hΩ.conjTranspose_mul_mul_same h
  rw [Matrix.conjTranspose_eq_transpose_of_trivial] at h1
  exact (Matrix.isUnit_iff_isUnit_det _).1 h1.isUnit

/-- An invertible weighted Gram matrix makes `mulVec` injective. -/
theorem mulVec_inj_of_gramW (Ω : Matrix (Fin n) (Fin n) ℝ) (A : Matrix (Fin n) m ℝ) (h : IsUnit (Aᵀ * Ω * A).det) :
    Function.Injective A.mulVec := by
  intro v w hvw
  have h1 := congrArg ((Aᵀ * Ω) *ᵥ ·) hvw
  simp only [Matrix.mulVec_mulVec] at h1
  have h2 := congrArg ((Aᵀ * Ω * A)⁻¹ *ᵥ ·) h1
  simp only [Matrix.mulVec_mulVec, Matrix.nonsing_inv_mul _ h, Matrix.one_mulVec] at h2
  exact h2

/-- The weighted least-squares coefficient is the unique solution of the weighted normal equations. -/
theorem betaHatW_of_normal (Ω : Matrix (Fin n) (Fin n) ℝ) (A : Matrix (Fin n) m ℝ) (z : Fin n → ℝ) (b : m → ℝ)
    (h : IsUnit (Aᵀ * Ω * A).det) (hn : Aᵀ *ᵥ (Ω *ᵥ (z - A *ᵥ b)) = 0) : betaHatW Ω A z = b := by
  have h1 : Aᵀ *ᵥ (Ω *ᵥ z) = (Aᵀ * Ω * A) *ᵥ b := by
    rw [Matrix.mulVec_sub, Matrix.mulVec_sub, sub_eq_zero] at hn
    rw [hn, Matrix.mulVec_mulVec, Matrix.mulVec_mulVec]
  unfold betaHatW
  rw [h1, Matrix.mulVec_mulVec, Matrix.nonsing_inv_mul _ h, Matrix.one_mulVec]

/-- The weighted normal equations: `X'Ω(y - Xβ̂_Ω) = 0`. -/
theorem normalW_eq (Ω : Matrix (Fin n) (Fin n) ℝ) (X : Matrix (Fin n) m ℝ) (y : Fin n → ℝ)
    (h : IsUnit (Xᵀ * Ω * X).det) : Xᵀ *ᵥ (Ω *ᵥ residW Ω X y) = 0 := by
  unfold residW betaHatW
  rw [Matrix.mulVec_sub, Matrix.mulVec_sub]
  simp only [Matrix.mulVec_mulVec, ← Matrix.mul_assoc, Matrix.mul_nonsing_inv _ h, Matrix.one_mul, sub_self]

end F3

/-- R-F3, the weighted Frisch-Waugh-Lovell. -/
theorem F3_fwl : F3_fwl_stmt := by
  intro n a b Ω X1 X2 y hΩ hQ M2
  set Q := fromCols X1 X2 with hQdef
  have hQinj := mulVec_inj_of_gramW Ω Q hQ
  have h2inj : Function.Injective X2.mulVec := by
    intro v w hvw
    have : Q *ᵥ Sum.elim 0 v = Q *ᵥ Sum.elim 0 w := by
      rw [hQdef, fromCols_mulVec_sumElim, fromCols_mulVec_sumElim, hvw]
    have h := hQinj this
    funext i
    exact congrFun h (Sum.inr i)
  have h2 : IsUnit (X2ᵀ * Ω * X2).det := gramW_unit Ω hΩ X2 h2inj
  set G2 := X2ᵀ * Ω * X2 with hG2
  have hM2 : M2 = 1 - X2 * G2⁻¹ * X2ᵀ * Ω := rfl
  have hM2X2 : M2 * X2 = 0 := by
    rw [hM2, Matrix.sub_mul, Matrix.one_mul]
    have : X2 * G2⁻¹ * X2ᵀ * Ω * X2 = X2 * (G2⁻¹ * (X2ᵀ * Ω * X2)) := by simp only [Matrix.mul_assoc]
    rw [this, ← hG2, Matrix.nonsing_inv_mul _ h2, Matrix.mul_one, sub_self]
  have h1inj : Function.Injective (M2 * X1).mulVec := by
    intro v w hvw
    have hd : (M2 * X1) *ᵥ (v - w) = 0 := by rw [Matrix.mulVec_sub, hvw, sub_self]
    set u := G2⁻¹ *ᵥ (X2ᵀ *ᵥ (Ω *ᵥ (X1 *ᵥ (v - w)))) with hu
    have hX1 : X1 *ᵥ (v - w) = X2 *ᵥ u := by
      rw [← Matrix.mulVec_mulVec, hM2, Matrix.sub_mulVec, Matrix.one_mulVec, sub_eq_zero] at hd
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
  have h1 : IsUnit ((M2 * X1)ᵀ * Ω * (M2 * X1)).det := gramW_unit Ω hΩ (M2 * X1) h1inj
  set β := betaHatW Ω Q y with hβ
  set r := residW Ω Q y with hr
  have hnorm : Qᵀ *ᵥ (Ω *ᵥ r) = 0 := normalW_eq Ω Q y hQ
  rw [hQdef, transpose_fromCols, fromRows_mulVec] at hnorm
  have hX1r : X1ᵀ *ᵥ (Ω *ᵥ r) = 0 := by
    funext i
    exact congrFun hnorm (Sum.inl i)
  have hX2r : X2ᵀ *ᵥ (Ω *ᵥ r) = 0 := by
    funext i
    exact congrFun hnorm (Sum.inr i)
  have hM2r : M2 *ᵥ r = r := by
    rw [hM2, Matrix.sub_mulVec, Matrix.one_mulVec]
    have : (X2 * G2⁻¹ * X2ᵀ * Ω) *ᵥ r = (X2 * G2⁻¹) *ᵥ (X2ᵀ *ᵥ (Ω *ᵥ r)) := by
      simp only [Matrix.mulVec_mulVec, Matrix.mul_assoc]
    rw [this, hX2r, Matrix.mulVec_zero, sub_zero]
  have hy : y = X1 *ᵥ (β ∘ Sum.inl) + X2 *ᵥ (β ∘ Sum.inr) + r := by
    rw [hr, residW, ← hβ, hQdef, fromCols_mulVec]
    abel
  have hM2y : M2 *ᵥ y - (M2 * X1) *ᵥ (β ∘ Sum.inl) = r := by
    rw [hy, Matrix.mulVec_add, Matrix.mulVec_add]
    simp only [Matrix.mulVec_mulVec, hM2X2, Matrix.zero_mulVec, hM2r, add_zero]
    abel
  have hM2tΩr : M2ᵀ *ᵥ (Ω *ᵥ r) = Ω *ᵥ r := by
    rw [hM2, transpose_sub, transpose_one, Matrix.sub_mulVec, Matrix.one_mulVec]
    have : (X2 * G2⁻¹ * X2ᵀ * Ω)ᵀ *ᵥ (Ω *ᵥ r) = (Ωᵀ * X2 * G2⁻¹ᵀ) *ᵥ (X2ᵀ *ᵥ (Ω *ᵥ r)) := by
      simp only [transpose_mul, transpose_transpose, Matrix.mulVec_mulVec, Matrix.mul_assoc]
    rw [this, hX2r, Matrix.mulVec_zero, sub_zero]
  have hn1 : (M2 * X1)ᵀ *ᵥ (Ω *ᵥ (M2 *ᵥ y - (M2 * X1) *ᵥ (β ∘ Sum.inl))) = 0 := by
    rw [hM2y, transpose_mul, ← Matrix.mulVec_mulVec, hM2tΩr, hX1r]
  have hb := betaHatW_of_normal Ω (M2 * X1) (M2 *ᵥ y) (β ∘ Sum.inl) h1 hn1
  refine ⟨?_, ?_⟩
  · rw [hb]
    rfl
  · unfold residW
    rw [hb, hM2y]

/-- Witness for `F3_fwl_stmt`: `Ω = diag(2, 1, 1)`, `X1 = e_1`, `X2 = e_2` in `ℝ^3`. -/
theorem F3_fwl_witness : F3_fwl_witness_stmt := by
  let X1 : Matrix (Fin 3) (Fin 1) ℝ := Matrix.of fun i _ => if i = 0 then 1 else 0
  let X2 : Matrix (Fin 3) (Fin 1) ℝ := Matrix.of fun i _ => if i = 1 then 1 else 0
  let Ω : Matrix (Fin 3) (Fin 3) ℝ := diagonal fun i => if i = 0 then 2 else 1
  have hG : (fromCols X1 X2)ᵀ * Ω * fromCols X1 X2 = diagonal fun i => Sum.elim (fun _ => (2 : ℝ)) (fun _ => 1) i := by
    ext (a | a) (b | b) <;> fin_cases a <;> fin_cases b <;>
      simp [X1, X2, Ω, Matrix.mul_apply, Matrix.diagonal_apply]
  refine ⟨3, 1, 1, Ω, X1, X2, by norm_num, by norm_num, by norm_num, ?_, ?_, ?_⟩
  · rw [Matrix.posDef_diagonal_iff]
    intro i
    split_ifs <;> norm_num
  · intro h
    have := congrFun (congrFun h 0) 0
    simp [Ω] at this
  · rw [hG, Matrix.det_diagonal]
    simp [Fintype.prod_sum_type]

/-! ## R-M1 -/

/-- R-M1: the design enters through its span. -/
theorem M1 : M1_stmt := by
  intro n a k X W Q hQ
  have hid : ∀ (βX : Fin a → ℝ) (βW : Fin k → ℝ),
      fromCols X W *ᵥ Sum.elim βX βW = fromCols X (W * Q) *ᵥ Sum.elim βX (Q⁻¹ *ᵥ βW) := by
    intro βX βW
    rw [fromCols_mulVec_sumElim, fromCols_mulVec_sumElim, Matrix.mulVec_mulVec, Matrix.mul_assoc,
      Matrix.mul_nonsing_inv _ hQ, Matrix.mul_one]
  have hid' : ∀ (βX : Fin a → ℝ) (βW : Fin k → ℝ),
      fromCols X (W * Q) *ᵥ Sum.elim βX βW = fromCols X W *ᵥ Sum.elim βX (Q *ᵥ βW) := by
    intro βX βW
    rw [fromCols_mulVec_sumElim, fromCols_mulVec_sumElim, Matrix.mulVec_mulVec]
  have hrange : LinearMap.range (fromCols X W).mulVecLin = LinearMap.range (fromCols X (W * Q)).mulVecLin := by
    ext v
    simp only [LinearMap.mem_range, Matrix.mulVecLin_apply]
    constructor
    · rintro ⟨x, rfl⟩
      refine ⟨Sum.elim (x ∘ Sum.inl) (Q⁻¹ *ᵥ (x ∘ Sum.inr)), ?_⟩
      rw [← hid, Sum.elim_comp_inl_inr]
    · rintro ⟨x, rfl⟩
      refine ⟨Sum.elim (x ∘ Sum.inl) (Q *ᵥ (x ∘ Sum.inr)), ?_⟩
      rw [← hid', Sum.elim_comp_inl_inr]
  refine ⟨hrange, ?_, hid⟩
  unfold Matrix.rank
  rw [hrange]

/-- Witness for `M1_stmt`: `k = 1`, `Q = (2)`. -/
theorem M1_witness : M1_witness_stmt := by
  refine ⟨1, Matrix.of fun _ _ => 2, by norm_num, ?_, ?_⟩
  · rw [Matrix.det_unique]
    simp
  · intro h
    have := congrFun (congrFun h 0) 0
    simp at this

end MultiScholaRProofs.Ruv
