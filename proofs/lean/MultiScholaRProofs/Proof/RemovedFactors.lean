module

public import MultiScholaRProofs.Spec.RemovedFactors
public import MultiScholaRProofs.Proof.RuvAlgebra
public import MultiScholaRProofs.Proof.RuvAlgebra2
import Mathlib.Tactic

/-!
# Proofs for group R, R-I1: KIT-184's removed factors (KIT-193 section 7.2)

Theorems whose types are the frozen statements of `Spec/RemovedFactors.lean` (frozen at 2ab4c189): `I1_UtWc`,
`I1_centred_diff`, `I1_centred_diff_rows`, `I1_span`, `I1_rank`, `I1_ctl_rank`, `I1_average`, and the witnesses
`I1_witness` (on `witnessE2`, `F = {0}`) and `I1_average_witness`. The key facts: centring commutes with
multiplication on the right, so `W_c = centre W` and `W - W_c` has constant rows; `U'1 = 0` gives `U'W_c = I`.
-/

public section

namespace MultiScholaRProofs.Ruv

open Matrix RuvIII

variable {n p s k : ℕ}

/-! ## Centring -/

/-- Centring commutes with multiplication on the right. -/
theorem centre_mul {m : ℕ} {a o : Type} [Fintype a] (A : Matrix (Fin m) a ℝ) (B : Matrix a o ℝ) :
    centre (A * B) = centre A * B := by
  unfold centre
  rw [Matrix.sub_mul, colMean_mul]
  congr 1

theorem centre_centre {m : ℕ} {a : Type} (A : Matrix (Fin m) a ℝ) : centre (centre A) = centre A :=
  centre_eq_self _ (colMean_centre A)

/-- Removing a constant row from every row does not change the centred matrix. -/
theorem centre_sub_rows {m : ℕ} {a : Type} (A : Matrix (Fin m) a ℝ) (c : a → ℝ) :
    centre (A - Matrix.of fun _ j => c j) = centre A := by
  ext i j
  have hn : (m : ℝ) ≠ 0 := by exact_mod_cast i.pos.ne'
  simp only [centre, colMean, Matrix.sub_apply, Matrix.of_apply, Finset.sum_sub_distrib, Finset.sum_const,
    Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
  field_simp
  ring

/-- `W_c = centre W`. -/
theorem Wc_eq_centre (r : RuvIII n p s k) : r.Wc = centre r.W := by
  unfold RuvIII.Wc RuvIII.W
  rw [Matrix.mul_assoc, Matrix.mul_assoc, centre_mul]

theorem Y_sub_Ystar (r : RuvIII n p s k) : r.Y - r.Ystar = r.W * r.alpha := by
  unfold RuvIII.Ystar
  exact sub_sub_cancel _ _

theorem Y_sub_adjusted (r : RuvIII n p s k) : r.Y - r.adjusted = r.Wc * r.alpha := by
  unfold RuvIII.adjusted
  exact sub_sub_cancel _ _

theorem submatrix_rows_mul {m a o : Type} [Fintype a] {m' : Type} (A : Matrix m a ℝ) (B : Matrix a o ℝ)
    (e : m' → m) : (A * B).submatrix e id = A.submatrix e id * B := by
  ext i j
  simp [Matrix.mul_apply]

/-! ## R-I1 (a), (b), (b') -/

theorem UtWc (r : RuvIII n p s k) (hac : Hac r) (hU : HU r) (hone : Hone r) : r.Uᵀ * r.Wc = 1 := by
  have hc : r.Uᵀ * centre r.Yc = r.ac := by
    unfold centre
    rw [Matrix.mul_sub, Ut_Yc, mul_const_rows _ _ (Ut_one r hU hone), sub_zero]
  unfold RuvIII.Wc
  rw [← Matrix.mul_assoc, ← Matrix.mul_assoc, hc]
  exact Matrix.mul_nonsing_inv _ hac

/-- R-I1 (a). -/
theorem I1_UtWc : I1_UtWc_stmt := by
  intro n p s k r hac hU hone
  have h1 := UtWc r hac hU hone
  refine ⟨h1, le_antisymm ?_ ?_⟩
  · simpa using Matrix.rank_le_width r.Wc
  · have h2 : (r.Uᵀ * r.Wc).rank ≤ r.Wc.rank := Matrix.rank_mul_le_right _ _
    rw [h1, Matrix.rank_one] at h2
    simpa using h2

/-- R-I1 (b). -/
theorem I1_centred_diff : I1_centred_diff_stmt := by
  intro n p s k r
  refine ⟨?_, ?_⟩
  · rw [Y_sub_Ystar, centre_mul, ← Wc_eq_centre]
  · rw [Y_sub_adjusted, centre_mul, Wc_eq_centre, centre_centre]

/-- R-I1 (b'): runs removed after the correction. -/
theorem I1_centred_diff_rows : I1_centred_diff_rows_stmt := by
  intro n p s k m r e
  refine ⟨?_, ?_⟩
  · rw [Y_sub_Ystar, submatrix_rows_mul, centre_mul]
  · rw [Y_sub_adjusted, submatrix_rows_mul, centre_mul, Wc_eq_centre]
    congr 1
    have h : (centre r.W).submatrix e id = r.W.submatrix e id - Matrix.of fun _ j => colMean r.W j := by
      ext i j
      simp [centre]
    rw [h, centre_sub_rows]

/-! ## R-I1 (c), (c'), (c'') -/

/-- R-I1 (c). -/
theorem I1_span : I1_span_stmt := by
  intro n p s k r F
  refine ⟨?_, fun hrank => ?_⟩
  · rw [Matrix.mulVecLin_mul]
    exact LinearMap.range_comp_le_range _ _
  · rw [Matrix.mulVecLin_mul]
    apply LinearMap.range_comp_of_range_eq_top
    apply Submodule.eq_top_of_finrank_eq
    rw [Module.finrank_fin_fun]
    exact hrank

/-- R-I1 (c'). -/
theorem I1_rank : I1_rank_stmt := by
  intro n p s k r F hac hU hone hrank
  set B := r.alpha.submatrix id (fun j : F => (j : Fin p))
  refine le_antisymm ((Matrix.rank_mul_le_right _ _).trans hrank.le) ?_
  have h : r.Uᵀ * (r.Wc * B) = B := by
    rw [← Matrix.mul_assoc, UtWc r hac hU hone, Matrix.one_mul]
  calc k = B.rank := hrank.symm
    _ = (r.Uᵀ * (r.Wc * B)).rank := by rw [h]
    _ ≤ (r.Wc * B).rank := Matrix.rank_mul_le_right _ _

/-- R-I1 (c''). -/
theorem I1_ctl_rank : I1_ctl_rank_stmt := by
  intro n p s k r F hac hF
  have hsel : r.alpha.submatrix id (fun j : F => (j : Fin p)) *
      (Matrix.of fun (f : F) (c : r.ctl) => if (f : Fin p) = (c : Fin p) then (1 : ℝ) else 0) = r.ac := by
    ext a c
    rw [Matrix.mul_apply]
    rw [Finset.sum_eq_single (⟨(c : Fin p), hF c.2⟩ : F)]
    · show r.alpha a (c : Fin p) * (if (c : Fin p) = (c : Fin p) then (1 : ℝ) else 0) = r.alpha a (c : Fin p)
      simp
    · intro f _ hf
      have : (f : Fin p) ≠ (c : Fin p) := fun h => hf (Subtype.ext h)
      show r.alpha a (f : Fin p) * (if (f : Fin p) = (c : Fin p) then (1 : ℝ) else 0) = 0
      simp [this]
    · intro h
      exact absurd (Finset.mem_univ _) h
  have hunit : IsUnit (r.ac * r.acᵀ) := (Matrix.isUnit_iff_isUnit_det _).mpr hac
  have h1 : (r.ac * r.acᵀ).rank = k := by
    rw [Matrix.rank_of_isUnit _ hunit]
    simp
  refine le_antisymm ?_ ?_
  · exact (Matrix.rank_le_card_height _).trans (Fintype.card_fin k).le
  · calc k = (r.ac * r.acᵀ).rank := h1.symm
      _ ≤ r.ac.rank := Matrix.rank_mul_le_left _ _
      _ = (r.alpha.submatrix id (fun j : F => (j : Fin p)) * _).rank := by rw [hsel]
      _ ≤ (r.alpha.submatrix id (fun j : F => (j : Fin p))).rank := Matrix.rank_mul_le_left _ _

/-! ## R-I1 (d) -/

/-- An averaging map sends `W_c a + t 1` into `span[1, A W_c]`. -/
theorem average_mem {d : ℕ} (A : Matrix (Fin d) (Fin n) ℝ) (hA : A *ᵥ (fun _ => (1 : ℝ)) = fun _ => 1)
    (V : Matrix (Fin n) (Fin k) ℝ) (a : Fin k → ℝ) (t : ℝ) :
    A *ᵥ (fun i => (V *ᵥ a) i + t) ∈
      LinearMap.range (fromCols (Matrix.of fun (_ : Fin d) (_ : Fin 1) => (1 : ℝ)) (A * V)).mulVecLin := by
  refine ⟨Sum.elim (fun _ => t) a, ?_⟩
  rw [Matrix.mulVecLin_apply, fromCols_mulVec_sumElim]
  have h1 : (fun i => (V *ᵥ a) i + t) = V *ᵥ a + t • (fun _ => (1 : ℝ)) := by
    funext i
    simp
  have h2 : (Matrix.of fun (_ : Fin d) (_ : Fin 1) => (1 : ℝ)) *ᵥ (fun _ => t) = t • fun _ => (1 : ℝ) := by
    funext x
    simp [Matrix.mulVec, dotProduct]
  rw [h1, Matrix.mulVec_add, Matrix.mulVec_mulVec, Matrix.mulVec_smul, hA, h2]
  exact add_comm _ _

/-- R-I1 (d). -/
theorem I1_average : I1_average_stmt := by
  intro n p s k d r A hA j
  have hrow : ∀ i l, r.W i l = r.Wc i l + colMean r.W l := by
    intro i l
    rw [Wc_eq_centre]
    simp [centre]
  refine ⟨?_, ?_⟩
  · have h : (fun i => (r.Y - r.Ystar) i j) =
        fun i => (r.Wc *ᵥ fun l => r.alpha l j) i + ∑ l, colMean r.W l * r.alpha l j := by
      funext i
      rw [Y_sub_Ystar]
      simp only [Matrix.mulVec, dotProduct, Matrix.mul_apply, hrow, add_mul, Finset.sum_add_distrib]
    rw [h]
    exact average_mem A hA r.Wc _ _
  · have h : (fun i => (r.Y - r.adjusted) i j) = fun i => (r.Wc *ᵥ fun l => r.alpha l j) i + 0 := by
      funext i
      rw [Y_sub_adjusted]
      simp [Matrix.mulVec, dotProduct, Matrix.mul_apply]
    rw [h]
    exact average_mem A hA r.Wc _ _

/-! ## Witnesses -/

theorem witnessE2_W : witnessE2.W = Matrix.of fun i _ => if i = 0 then (1 : ℝ) else 0 := by
  unfold RuvIII.W
  rw [witnessE2_acac, inv_one, Matrix.mul_one]
  ext i a
  rw [Matrix.mul_apply, witnessE2_sum_ctl]
  fin_cases a
  simp only [RuvIII.Yc, RuvIII.ac, Matrix.submatrix_apply, Matrix.transpose_apply, id, witnessE2_alpha]
  simp [witnessE2]

/-- Witness for `I1_UtWc_stmt`, `I1_rank_stmt` and the equality of `I1_span_stmt`: `witnessE2`, `F = {0}`. -/
theorem I1_witness : I1_witness_stmt := by
  refine ⟨4, 3, 1, 1, witnessE2, {0}, by norm_num, by norm_num, ?_, ?_, witnessE2_hac, witnessE2_hU,
    witnessE2_hone, I1_ctl_rank 4 3 1 1 witnessE2 {0} witnessE2_hac (by simp [witnessE2])⟩
  · intro h
    have h1 : (1 : Fin 3) ∈ witnessE2.ctl := h ▸ Finset.mem_univ _
    simp [witnessE2] at h1
  · intro h
    have h1 := congrFun (congrFun h 1) 0
    rw [Wc_eq_centre, witnessE2_W] at h1
    simp [centre, colMean] at h1

/-- Witness for `I1_average_stmt`: two runs averaged into one sample. -/
theorem I1_average_witness : I1_average_witness_stmt := by
  refine ⟨1, 2, Matrix.of fun _ _ => (1 / 2 : ℝ), by norm_num, by norm_num, ?_⟩
  funext x
  simp [Matrix.mulVec, dotProduct]

end MultiScholaRProofs.Ruv
