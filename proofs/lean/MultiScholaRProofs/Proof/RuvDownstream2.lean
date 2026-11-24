module

public import MultiScholaRProofs.Spec.RuvDownstream
public import MultiScholaRProofs.Proof.RuvDownstream
import Mathlib.Tactic

/-!
# Proofs for group R, what consumes the corrected values or the factors (KIT-193 section 7.2), priorities 2 to 4

Theorems whose types are the frozen statements of `Spec/RuvDownstream.lean`: R-G1 (the orthogonal illustration) with
`G1_witness`.
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

end MultiScholaRProofs.Ruv
