module

public import MultiScholaRProofs.Spec.SharperWitnesses
public import MultiScholaRProofs.Proof.RuvAlgebra
import Mathlib.Tactic

/-!
# Proofs of the two sharper witness statements (KIT-193 proofs audit P03)

`perm_measure_moving_witness`: the uniform law on the four unit vectors of `ℝ^4`, `T f = f 0`, the atom `e_0` and the
swap of the first two samples. `G0_interior_witness`: three runs in one replicate set, `X = M = (1, 1, 1)'`,
`U = (1, -1, 0)'`, two features of which the first is the only control; the corrected second feature is `(0, 0, 1)`,
whose residual on `X` is `(-1/3, -1/3, 2/3)`.
-/

public section

-- A style linter only (a simp argument not needed in every `fin_cases` branch); no bearing on what is proved.
set_option linter.unusedSimpArgs false

namespace MultiScholaRProofs.Permutation

open MeasureTheory ENNReal

/-- The unit vectors of `ℝ^4`. -/
noncomputable def unitVec (i : Fin 4) : Fin 4 → ℝ := Pi.single i 1

theorem unitVec_comp (i : Fin 4) (σ : Equiv.Perm (Fin 4)) : unitVec i ∘ ⇑σ = unitVec (σ.symm i) := by
  funext j
  simp only [unitVec, Function.comp_apply, Pi.single_apply, Equiv.eq_symm_apply]

/-- The uniform law on the four unit vectors. -/
noncomputable def unitLaw : Measure (Fin 4 → ℝ) := ∑ i : Fin 4, (4⁻¹ : ℝ≥0∞) • Measure.dirac (unitVec i)

theorem unitLaw_apply (S : Set (Fin 4 → ℝ)) :
    unitLaw S = ∑ i : Fin 4, (4⁻¹ : ℝ≥0∞) * S.indicator 1 (unitVec i) := by
  simp [unitLaw, Measure.dirac_apply]

theorem perm_measure_moving_witness : perm_measure_moving_witness_stmt := by
  refine ⟨unitLaw, fun f => f 0, ⟨?_⟩, measurable_pi_apply 0, fun σ => ?_, unitVec 0, Equiv.swap 0 1, ?_, ?_, ?_⟩
  · rw [unitLaw_apply]
    simp only [Set.indicator_univ, Pi.one_apply, mul_one, Finset.sum_const, Finset.card_univ, Fintype.card_fin,
      nsmul_eq_mul]
    exact ENNReal.mul_inv_cancel (by norm_num) (by norm_num)
  · have hc : Measurable fun f : Fin 4 → ℝ => f ∘ ⇑σ := Measurable.of_eval fun i => measurable_pi_apply (σ i)
    ext S hS
    rw [Measure.map_apply hc hS, unitLaw_apply, unitLaw_apply]
    have h : ∀ i : Fin 4, ((fun f : Fin 4 → ℝ => f ∘ ⇑σ) ⁻¹' S).indicator 1 (unitVec i) =
        S.indicator (1 : (Fin 4 → ℝ) → ℝ≥0∞) (unitVec (σ.symm i)) := by
      intro i
      have hiff : unitVec i ∈ (fun f : Fin 4 → ℝ => f ∘ ⇑σ) ⁻¹' S ↔ unitVec (σ.symm i) ∈ S := by
        rw [Set.mem_preimage, unitVec_comp]
      by_cases h : unitVec (σ.symm i) ∈ S
      · rw [Set.indicator_of_mem (hiff.2 h), Set.indicator_of_mem h]
        rfl
      · rw [Set.indicator_of_notMem (fun h' => h (hiff.1 h')), Set.indicator_of_notMem h]
    simp only [h]
    exact Equiv.sum_comp σ.symm (fun j => (4⁻¹ : ℝ≥0∞) * S.indicator 1 (unitVec j))
  · rw [unitLaw_apply]
    refine lt_of_lt_of_le ?_ (Finset.single_le_sum (fun _ _ => bot_le) (Finset.mem_univ (0 : Fin 4)))
    simp
  · intro h
    have := congrFun h 0
    simp [unitVec] at this
  · simp [unitVec]

end MultiScholaRProofs.Permutation

namespace MultiScholaRProofs.Ruv

open Matrix RuvIII

/-- The interior instance: three runs in one set, two features, the first the only control. -/
def witnessG0i : RuvIII 3 2 1 1 where
  Y := Matrix.of fun i j => if j = 0 then (if i = 0 then 1 else 0) else (if i = 2 then 1 else 0)
  M := Matrix.of fun _ _ => 1
  ctl := {0}
  U := Matrix.of fun i _ => if i = 0 then 1 else if i = 1 then -1 else 0

theorem witnessG0i_mem0 : (0 : Fin 2) ∈ witnessG0i.ctl := by simp [witnessG0i]

theorem witnessG0i_sum_ctl (f : witnessG0i.ctl → ℝ) : ∑ c, f c = f ⟨0, witnessG0i_mem0⟩ := by
  have hu : ∀ c : witnessG0i.ctl, c = ⟨0, witnessG0i_mem0⟩ := fun c =>
    Subtype.ext (by simpa [witnessG0i] using c.2)
  rw [Finset.sum_eq_single (⟨0, witnessG0i_mem0⟩ : witnessG0i.ctl) (fun c _ hc => absurd (hu c) hc)
    (fun h => absurd (Finset.mem_univ _) h)]

theorem witnessG0i_MtM : witnessG0i.Mᵀ * witnessG0i.M = Matrix.of fun _ _ => (3 : ℝ) := by
  ext a b
  fin_cases a
  fin_cases b
  simp [witnessG0i, Matrix.mul_apply, Fin.sum_univ_three]

theorem witnessG0i_MtM_inv : (witnessG0i.Mᵀ * witnessG0i.M)⁻¹ = Matrix.of fun _ _ => (1 / 3 : ℝ) := by
  apply Matrix.inv_eq_left_inv
  rw [witnessG0i_MtM]
  ext a b
  fin_cases a
  fin_cases b
  simp [Matrix.mul_apply]

theorem witnessG0i_residop : residop witnessG0i.M =
    Matrix.of fun i j => if i = j then (2 / 3 : ℝ) else -1 / 3 := by
  unfold residop
  rw [witnessG0i_MtM_inv]
  ext i j
  fin_cases i <;> fin_cases j <;> simp [witnessG0i, Matrix.mul_apply] <;> norm_num

theorem witnessG0i_alpha : witnessG0i.alpha = Matrix.of fun _ j => if j = 0 then (1 : ℝ) else 0 := by
  ext a j
  fin_cases a
  fin_cases j <;> simp [witnessG0i, RuvIII.alpha, Matrix.mul_apply, Fin.sum_univ_three, Function.swap]

theorem witnessG0i_acac : witnessG0i.ac * witnessG0i.acᵀ = 1 := by
  ext a b
  fin_cases a
  fin_cases b
  rw [Matrix.mul_apply]
  simp only [RuvIII.ac, Matrix.submatrix_apply, Matrix.transpose_apply, id]
  rw [witnessG0i_sum_ctl]
  simp [witnessG0i_alpha]

theorem witnessG0i_W : witnessG0i.W = Matrix.of fun i _ => if i = 0 then (1 : ℝ) else 0 := by
  unfold RuvIII.W
  rw [witnessG0i_acac, inv_one, Matrix.mul_one]
  ext i a
  fin_cases a
  rw [Matrix.mul_apply]
  simp only [RuvIII.Yc, RuvIII.ac, Matrix.submatrix_apply, Matrix.transpose_apply, id]
  rw [witnessG0i_sum_ctl, witnessG0i_alpha]
  fin_cases i <;> simp [witnessG0i]

theorem witnessG0i_Ystar : witnessG0i.Ystar = Matrix.of fun i j => if j = 1 ∧ i = 2 then (1 : ℝ) else 0 := by
  unfold RuvIII.Ystar
  rw [witnessG0i_W, witnessG0i_alpha]
  ext i j
  fin_cases i <;> fin_cases j <;> simp [witnessG0i, Matrix.mul_apply]

theorem G0_interior_witness : G0_interior_witness_stmt := by
  have hdet : IsUnit (witnessG0i.Mᵀ * witnessG0i.M).det := by
    rw [witnessG0i_MtM, Matrix.det_unique]
    simp
  refine ⟨3, 2, 1, 1, 1, witnessG0i, witnessG0i.M, by norm_num, by norm_num, by norm_num, ?_, hdet, ?_, ?_, ?_,
    ?_, hdet, ?_⟩
  · intro h
    have h1 : (1 : Fin 2) ∈ witnessG0i.ctl := by rw [h]; exact Finset.mem_univ _
    simp [witnessG0i] at h1
  · unfold Hac
    rw [witnessG0i_acac, Matrix.det_one]
    exact isUnit_one
  · unfold HU
    rw [witnessG0i_residop]
    ext i j
    fin_cases i <;> fin_cases j <;> simp [witnessG0i, Matrix.mul_apply, Fin.sum_univ_three] <;> norm_num
  · unfold Hone
    rw [witnessG0i_residop]
    ext i
    fin_cases i <;> simp [Matrix.mulVec, dotProduct, Fin.sum_univ_three] <;> norm_num
  · rw [witnessG0i_residop]
    ext i j
    fin_cases i <;> fin_cases j <;> simp [witnessG0i, Matrix.mul_apply, Fin.sum_univ_three] <;> norm_num
  · intro h
    have h1 := congrFun (congrFun h 2) 1
    rw [witnessG0i_residop, witnessG0i_Ystar, Matrix.mul_apply] at h1
    simp [Fin.sum_univ_three] at h1

end MultiScholaRProofs.Ruv
