module

public import MultiScholaRProofs.Spec.RuvIIIC
public import MultiScholaRProofs.Proof.RuvAlgebra
import Mathlib.Tactic

/-!
# Proofs for group R, RUV-III-C with missing values (KIT-193 section 7.2)

Theorems whose types are the frozen statements of `Spec/RuvIIIC.lean`: R-J1 with `J1_witness`; R-J2 (`J2`,
`J2_contribution`) with `J2_witness`; R-L1 (`L1`, `L1_invariance`) with `L1_invariance_witness`; R-L2 with `L2_witness`.
-/

public section

namespace MultiScholaRProofs.Ruv

open Matrix Finset

/-! ## R-J1 -/

theorem mask_univ {n : ℕ} : mask (univ : Finset (Fin n)) = 1 := by
  unfold mask
  simp

theorem rows_complete {n p s k : ℕ} (r : RuvIIIC n p s k) (j : Fin p) (hj : ∀ i, r.obs i j = true) :
    r.rows j = univ := by
  unfold RuvIIIC.rows
  exact Finset.filter_true_of_mem fun i _ => hj i

/-- R-J1: complete features share one `W`, and their corrected values are those of one RUV-III computation. -/
theorem J1 : J1_stmt := by
  intro n p s k r j j' hj hj'
  refine ⟨by rw [rows_complete r j hj, rows_complete r j' hj'], fun i => ?_⟩
  have hY0 : r.Y0 i j = r.Y i j := by simp [RuvIIIC.Y0, hj i]
  have haOn : r.aOn univ = r.complete.alpha := by
    unfold RuvIIIC.aOn
    rw [mask_univ, Matrix.mul_one]
    rfl
  have hWa : r.WOn univ * r.aOn univ = r.complete.W * r.complete.alpha := by
    unfold RuvIIIC.WOn RuvIIIC.acOn
    rw [haOn, mask_univ, Matrix.one_mul]
    rfl
  unfold RuvIIIC.ystar
  simp only [hj i, ↓reduceIte]
  rw [rows_complete r j hj, hWa]
  show some (r.Y i j - _) = some (r.Y0 i j - _)
  rw [hY0]

/-- Witness for `J1_stmt`: one run, three features, the third missing. -/
theorem J1_witness : J1_witness_stmt := by
  refine ⟨1, 3, 0, 0, ⟨0, fun _ j => decide (j ≠ 2), 0, ∅, fun _ => 0, fun _ => 0⟩, 0, 1, 2, by decide,
    fun _ => rfl, fun _ => rfl, ⟨0, rfl⟩⟩

/-! ## R-J2 -/

/-- `Σ (x_i - c)² = Σ x_i² - 2c Σ x_i + m c²`. -/
theorem sum_sq_dev {m : ℕ} (x : Fin m → ℝ) (c : ℝ) :
    ∑ i, (x i - c) ^ 2 = ∑ i, x i ^ 2 - 2 * c * ∑ i, x i + m * c ^ 2 := by
  have h : ∀ i, (x i - c) ^ 2 = x i ^ 2 - 2 * c * x i + c ^ 2 := fun i => by ring
  simp only [h, Finset.sum_add_distrib, Finset.sum_sub_distrib, ← Finset.mul_sum, Finset.sum_const,
    Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]

/-- R-J2: the zero-filled run adds `(m/(m+1)) x̄²` to the observed spread. -/
theorem J2 : J2_stmt := by
  intro m x hm
  have hm' : (m : ℝ) ≠ 0 := by exact_mod_cast hm.ne'
  have hm1 : (m : ℝ) + 1 ≠ 0 := by positivity
  set v : Fin (m + 1) → ℝ := Fin.snoc x 0 with hv
  have hsum : ∑ i, v i = ∑ i, x i := by
    rw [Fin.sum_univ_castSucc]
    simp [hv]
  have hsq : ∑ i, v i ^ 2 = ∑ i, x i ^ 2 := by
    rw [Fin.sum_univ_castSucc]
    simp [hv]
  rw [sum_sq_dev, sum_sq_dev, hsum, hsq]
  unfold mean
  rw [hsum]
  push_cast
  field_simp
  ring

/-- R-J2, the link: on one replicate set, `‖P v‖²` is the sum of squared deviations from the set's mean. -/
theorem J2_contribution : J2_contribution_stmt := by
  intro m v
  set M : Matrix (Fin (m + 1)) (Fin 1) ℝ := Matrix.of fun _ _ => 1 with hM
  have hm1 : ((m : ℝ) + 1) ≠ 0 := by positivity
  have hMtM : Mᵀ * M = Matrix.of fun _ _ => ((m : ℝ) + 1) := by
    ext a b
    simp [hM, Matrix.mul_apply]
  have hinv : (Mᵀ * M)⁻¹ = Matrix.of fun _ _ => 1 / ((m : ℝ) + 1) := by
    apply Matrix.inv_eq_left_inv
    rw [hMtM]
    ext a b
    fin_cases a
    fin_cases b
    simp [Matrix.mul_apply]
    field_simp
  have hP : residop M *ᵥ v = fun i => v i - mean v := by
    ext i
    rw [residop, hinv, Matrix.sub_mulVec, Matrix.one_mulVec]
    simp only [Pi.sub_apply, Matrix.mulVec, dotProduct, Matrix.mul_apply, Matrix.transpose_apply, hM,
      Matrix.of_apply, mean]
    simp only [Finset.sum_const, Finset.card_univ, Fintype.card_fin, mul_one, one_mul]
    push_cast
    rw [Finset.sum_div]
    congr 1
    refine Finset.sum_congr rfl fun j _ => ?_
    simp only [one_smul]
    ring
  rw [hP]
  simp only [dotProduct, sq]

/-- Witness for `J2_stmt`: `m = 1`, `x = 1`. -/
theorem J2_witness : J2_witness_stmt :=
  ⟨1, fun _ => 1, by norm_num, by simp [mean]⟩

/-! ## R-L1 -/

/-- R-L1: the contrast of the corrected target. -/
theorem L1 : L1_stmt := by
  intro R C _ _ k y YC aC a lam
  unfold ystarT
  rw [dotProduct_sub, Matrix.mul_assoc, ← Matrix.mulVec_mulVec, Matrix.dotProduct_mulVec,
    ← Matrix.mulVec_transpose]

/-- `A'(1c') = 0` style: columns that are constant multiples of `1` are killed by `U'` when `U'1 = 0`. -/
theorem Ut_const_rows {R : Type} [Fintype R] {k : ℕ} {q : Type} (U : Matrix R (Fin k) ℝ) (e : q → ℝ)
    (h : Uᵀ *ᵥ (fun _ => (1 : ℝ)) = 0) : Uᵀ * (Matrix.of fun (_ : R) c => e c) = 0 := by
  ext a c
  have h1 := congrFun h a
  simp only [Matrix.mulVec, dotProduct, mul_one, Pi.zero_apply, Matrix.transpose_apply] at h1
  simp only [Matrix.mul_apply, Matrix.transpose_apply, Matrix.of_apply, Matrix.zero_apply]
  rw [← Finset.sum_mul, h1, zero_mul]

/-- R-L1, invariance of the fold-change error under constants added to the target or to one control. -/
theorem L1_invariance : L1_invariance_stmt := by
  intro R C _ _ _ k U y YC lam t l hU hlam
  refine ⟨?_, ?_⟩
  · unfold foldErr
    have h : Uᵀ *ᵥ (y + fun _ => t) = Uᵀ *ᵥ y := by
      rw [Matrix.mulVec_add]
      have h0 : Uᵀ *ᵥ (fun _ : R => t) = 0 := by
        have : (fun _ : R => t) = t • (fun _ : R => (1 : ℝ)) := by
          funext i
          simp
        rw [this, Matrix.mulVec_smul, hU, smul_zero]
      rw [h0, add_zero]
    rw [h]
  · unfold foldErr
    set E : Matrix R C ℝ := Matrix.of fun _ c => if c = l then t else 0 with hE
    have hUE : Uᵀ * (YC + E) = Uᵀ * YC := by
      rw [Matrix.mul_add, hE, Ut_const_rows U _ hU, add_zero]
    have hEl : (YC + E)ᵀ *ᵥ lam = YCᵀ *ᵥ lam := by
      rw [transpose_add, Matrix.add_mulVec]
      have h0 : Eᵀ *ᵥ lam = 0 := by
        funext c
        simp only [hE, Matrix.mulVec, dotProduct, Matrix.transpose_apply, Matrix.of_apply, Pi.zero_apply]
        rw [← Finset.mul_sum, hlam, mul_zero]
      rw [h0, add_zero]
    rw [hUE, hEl]

/-- Witness for `L1_invariance_stmt`: `U = λ = (1, -1, 0, 0)'`. -/
theorem L1_invariance_witness : L1_invariance_witness_stmt := by
  refine ⟨Matrix.of fun i _ => if i = 0 then 1 else if i = 1 then -1 else 0,
    fun i => if i = 0 then 1 else if i = 1 then -1 else 0, ?_, ?_, ?_, ?_⟩
  · ext a
    fin_cases a
    simp [Matrix.mulVec, dotProduct, Fin.sum_univ_four]
  · intro h
    have := congrFun (congrFun h 0) 0
    simp at this
  · simp [Fin.sum_univ_four]
  · intro h
    have := congrFun h 0
    simp at this

/-! ## R-L2 -/

/-- Cauchy-Schwarz for the dot product: `(u'v)² ≤ (u'u)(v'v)`. -/
theorem dot_sq_le {ι : Type} [Fintype ι] (u v : ι → ℝ) : (u ⬝ᵥ v) ^ 2 ≤ (u ⬝ᵥ u) * (v ⬝ᵥ v) := by
  have h := Finset.sum_mul_sq_le_sq_mul_sq Finset.univ u v
  simp only [dotProduct, ← sq]
  exact h

theorem dot_self_nonneg' {ι : Type} [Fintype ι] (u : ι → ℝ) : 0 ≤ u ⬝ᵥ u :=
  Finset.sum_nonneg fun _ _ => mul_self_nonneg _

/-- R-L2: the bound on the fold-change error. -/
theorem L2 : L2_stmt := by
  intro C _ k aC dC a κ hκ hbound
  set G := aC * aCᵀ with hG
  have hGinj : Function.Injective G.mulVec := by
    intro v w hvw
    have h0 : G *ᵥ (v - w) = 0 := by rw [Matrix.mulVec_sub, hvw, sub_self]
    have h1 := hbound (v - w)
    rw [h0, dotProduct_zero] at h1
    have h2 : (v - w) ⬝ᵥ (v - w) = 0 :=
      le_antisymm (by nlinarith [dot_self_nonneg' (v - w)]) (dot_self_nonneg' _)
    exact sub_eq_zero.1 (dotProduct_self_eq_zero.1 h2)
  have hGu : IsUnit G := Matrix.mulVec_injective_iff_isUnit.1 hGinj
  have hGd : IsUnit G.det := (Matrix.isUnit_iff_isUnit_det G).1 hGu
  set z := G⁻¹ *ᵥ a with hz
  have hGz : G *ᵥ z = a := by rw [hz, Matrix.mulVec_mulVec, Matrix.mul_nonsing_inv _ hGd, Matrix.one_mulVec]
  set u := (aCᵀ * G⁻¹) *ᵥ a with hu
  have hu' : u = aCᵀ *ᵥ z := by rw [hu, hz, Matrix.mulVec_mulVec]
  -- `u'u = z'a`
  have huu : u ⬝ᵥ u = z ⬝ᵥ a := by
    rw [hu', Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose, transpose_transpose, Matrix.mulVec_mulVec,
      ← hGz, hG, dotProduct_comm]
  -- `κ z'z ≤ z'a` and `(z'a)² ≤ (z'z)(a'a)` give `z'a ≤ a'a / κ`
  have hza : z ⬝ᵥ a ≤ (a ⬝ᵥ a) / κ := by
    have h1 : κ * (z ⬝ᵥ z) ≤ z ⬝ᵥ a := by
      have := hbound z
      rwa [hGz] at this
    have h2 := dot_sq_le z a
    have hzz := dot_self_nonneg' z
    have haa := dot_self_nonneg' a
    rw [le_div_iff₀ hκ]
    have hw : 0 ≤ z ⬝ᵥ a := le_trans (mul_nonneg hκ.le hzz) h1
    have h3 : κ * (z ⬝ᵥ a) ^ 2 ≤ (z ⬝ᵥ a) * (a ⬝ᵥ a) := by
      nlinarith [mul_le_mul_of_nonneg_right h1 haa, mul_le_mul_of_nonneg_left h2 hκ.le]
    rcases hw.eq_or_lt with h0 | hpos
    · rw [← h0]
      simpa using haa
    · by_contra hc
      rw [not_le] at hc
      nlinarith [mul_lt_mul_of_pos_left hc hpos]
  have hdu : (dC ⬝ᵥ u) ^ 2 ≤ (dC ⬝ᵥ dC) * ((a ⬝ᵥ a) / κ) := by
    calc (dC ⬝ᵥ u) ^ 2 ≤ (dC ⬝ᵥ dC) * (u ⬝ᵥ u) := dot_sq_le dC u
      _ ≤ (dC ⬝ᵥ dC) * ((a ⬝ᵥ a) / κ) := by
          rw [huu]
          exact mul_le_mul_of_nonneg_left hza (dot_self_nonneg' dC)
  have h := Real.abs_le_sqrt hdu
  rw [Real.sqrt_mul (dot_self_nonneg' dC), Real.sqrt_div' _ hκ.le, ← mul_div_assoc] at h
  exact h

/-- Witness for `L2_stmt`: `a_C = [I_2, 0]`, `κ = 1`. -/
theorem L2_witness : L2_witness_stmt := by
  refine ⟨Matrix.of fun i j => if (i : ℕ) = (j : ℕ) then 1 else 0, 1, by norm_num, fun x => ?_⟩
  have h : (Matrix.of fun (i : Fin 2) (j : Fin 3) => if (i : ℕ) = (j : ℕ) then (1 : ℝ) else 0) *
      (Matrix.of fun (i : Fin 2) (j : Fin 3) => if (i : ℕ) = (j : ℕ) then (1 : ℝ) else 0)ᵀ = 1 := by
    ext a b
    fin_cases a <;> fin_cases b <;> simp [Matrix.mul_apply, Fin.sum_univ_three]
  rw [h, Matrix.one_mulVec, one_mul]

end MultiScholaRProofs.Ruv
