module

public import MultiScholaRProofs.Spec.RuvAlgebra
import Mathlib.Tactic

/-!
# Proofs for group R, the RUV-III algebra (KIT-193 section 7.2), priority 1

Theorems whose types are the frozen statements of `Spec/RuvAlgebra.lean`: R-A1 (correction form and subspace), R-A2
(projection identities), R-A3 (`k = 0`), R-C1 (restoring is the centred formula), with the witness `A2_witness`
(shared by R-A2 and R-C1). The lemmas `residop_transpose`, `Ut_Yc`, `UtW`, `colMean_mul` are used again by
`Proof/RuvDownstream.lean`.
-/

public section

namespace MultiScholaRProofs.Ruv

open Matrix RuvIII

/-! ## Shared lemmas -/

/-- `R_X` is symmetric (the inverse of the symmetric `X'X` is symmetric, junk or not). -/
theorem residop_transpose {n : ℕ} {m : Type} [Fintype m] [DecidableEq m] (X : Matrix (Fin n) m ℝ) :
    (residop X)ᵀ = residop X := by
  unfold residop
  rw [transpose_sub, transpose_one, transpose_mul, transpose_mul, transpose_transpose, transpose_nonsing_inv,
    transpose_mul, transpose_transpose, Matrix.mul_assoc]

/-- `U'Y[, ctl] = a_c`. -/
theorem Ut_Yc {n p s k : ℕ} (r : RuvIII n p s k) : r.Uᵀ * r.Yc = r.ac := by
  ext i c
  simp [RuvIII.Yc, RuvIII.ac, RuvIII.alpha, Matrix.mul_apply]

/-- `U'W = I` when `a_c a_c'` is invertible. -/
theorem UtW {n p s k : ℕ} (r : RuvIII n p s k) (h : Hac r) : r.Uᵀ * r.W = 1 := by
  unfold RuvIII.W
  rw [← Matrix.mul_assoc, ← Matrix.mul_assoc, Ut_Yc]
  exact Matrix.mul_nonsing_inv _ h

/-- `U'Y* = 0` when `a_c a_c'` is invertible. -/
theorem Ut_Ystar {n p s k : ℕ} (r : RuvIII n p s k) (h : Hac r) : r.Uᵀ * r.Ystar = 0 := by
  unfold RuvIII.Ystar
  rw [Matrix.mul_sub, ← Matrix.mul_assoc, UtW r h, Matrix.one_mul]
  unfold RuvIII.alpha
  exact sub_self _

/-- Column means are linear: the mean of a difference. -/
theorem colMean_sub {n : ℕ} {m : Type} (A B : Matrix (Fin n) m ℝ) :
    colMean (A - B) = fun j => colMean A j - colMean B j := by
  funext j
  simp [colMean, Finset.sum_sub_distrib, sub_div]

/-- Column means of a product: `1'(AB)/n = (1'A/n) B`. -/
theorem colMean_mul {n : ℕ} {m o : Type} [Fintype m] (A : Matrix (Fin n) m ℝ) (B : Matrix m o ℝ) :
    colMean (A * B) = fun j => ∑ l, colMean A l * B l j := by
  funext j
  simp only [colMean, Matrix.mul_apply]
  rw [Finset.sum_comm, Finset.sum_div]
  refine Finset.sum_congr rfl fun l _ => ?_
  rw [div_mul_eq_mul_div, Finset.sum_mul]

/-- A centred matrix has zero column means (at `n = 0` both are Lean's `0 / 0 = 0` over an empty sum). -/
theorem colMean_centre {n : ℕ} {m : Type} (A : Matrix (Fin n) m ℝ) : colMean (centre A) = 0 := by
  funext j
  rcases Nat.eq_zero_or_pos n with hn | hn
  · subst hn
    simp [colMean]
  · have hn' : (n : ℝ) ≠ 0 := by exact_mod_cast hn.ne'
    simp only [centre, colMean_sub, Pi.zero_apply]
    simp only [colMean, Matrix.of_apply, Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
    field_simp
    ring

theorem centre_eq_self {n : ℕ} {m : Type} (A : Matrix (Fin n) m ℝ) (h : colMean A = 0) : centre A = A := by
  ext i j
  simp [centre, h]

/-! ## R-A1, R-A2, R-A3 -/

/-- R-A1: correction form and subspace. -/
theorem A1 : A1_stmt := by
  intro n p s k r
  have h1 : r.Y - r.Ystar = r.W * r.alpha := by
    unfold RuvIII.Ystar
    exact sub_sub_cancel _ _
  refine ⟨h1, fun j => ⟨fun l => r.alpha l j, ?_⟩, ?_⟩
  · rw [h1]
    funext i
    simp [Matrix.mulVec, dotProduct, Matrix.mul_apply]
  · rw [h1]
    exact (Matrix.rank_mul_le_left _ _).trans (Matrix.rank_le_width _)

/-- R-A2: projection identities. -/
theorem A2 : A2_stmt := by
  intro n p s k r h
  have hW := UtW r h
  refine ⟨hW, ?_, Ut_Ystar r h, ?_⟩
  · rw [Matrix.mul_assoc, ← Matrix.mul_assoc r.Uᵀ, hW, Matrix.one_mul]
  · have hsub : r.Ystar.submatrix id (fun c : r.ctl => (c : Fin p)) = r.Yc - r.W * r.ac := by
      ext i c
      simp [RuvIII.Ystar, RuvIII.Yc, RuvIII.ac, Matrix.mul_apply]
    rw [hsub, Matrix.mul_sub, Matrix.mul_one]
    unfold RuvIII.W
    simp only [Matrix.mul_assoc]

/-- R-A3: `k = 0` is the identity. -/
theorem A3 : A3_stmt := by
  intro n p s r
  ext i j
  simp [RuvIII.Ystar]

/-! ## R-C1 -/

/-- R-C1: restoring levels of `Y - Wα` gives `Y - W_c α`; `W_c` is already centred. -/
theorem C1 : C1_stmt := by
  intro n p s k r _
  set D := r.acᵀ * (r.ac * r.acᵀ)⁻¹ * r.alpha with hD
  have hWα : r.W * r.alpha = r.Yc * D := by
    unfold RuvIII.W
    simp only [hD, Matrix.mul_assoc]
  have hWcα : r.Wc * r.alpha = centre r.Yc * D := by
    unfold RuvIII.Wc
    simp only [hD, Matrix.mul_assoc]
  have hcm := colMean_mul r.Yc D
  refine ⟨?_, ?_⟩
  · ext i j
    have hm : (Matrix.of (fun _ l => colMean r.Yc l) * D) i j = colMean (r.Yc * D) j := by
      rw [hcm]
      simp [Matrix.mul_apply]
    simp only [restore, RuvIII.adjusted, RuvIII.Ystar, hWα, hWcα, centre, Matrix.sub_mul, colMean_sub,
      Matrix.add_apply, Matrix.sub_apply, Matrix.of_apply]
    rw [hm]
    ring
  · apply centre_eq_self
    unfold RuvIII.Wc
    rw [Matrix.mul_assoc, colMean_mul, colMean_centre]
    funext j
    simp

/-! ## Witness -/

/-- The witness instance: two runs, one feature (a control), one set, `k = 1`, `U = e_1`. -/
def witnessA2 : RuvIII 2 1 1 1 where
  Y := Matrix.of fun _ _ => 1
  M := Matrix.of fun _ _ => 1
  ctl := Finset.univ
  U := Matrix.of fun i _ => if i = 0 then 1 else 0

theorem witnessA2_hac : Hac witnessA2 := by
  have h : witnessA2.ac * witnessA2.acᵀ = 1 := by
    ext a b
    fin_cases a
    fin_cases b
    rw [Matrix.mul_apply]
    simp only [witnessA2, RuvIII.ac, RuvIII.alpha, Matrix.submatrix_apply, Matrix.transpose_apply, Matrix.mul_apply,
      Matrix.of_apply, Matrix.one_apply, id]
    simp
    rfl
  unfold Hac
  rw [h, Matrix.det_one]
  exact isUnit_one

/-- Witness for `A2_stmt` and `C1_stmt`. -/
theorem A2_witness : A2_witness_stmt :=
  ⟨2, 1, 1, 1, witnessA2, by norm_num, by norm_num, witnessA2_hac⟩

end MultiScholaRProofs.Ruv
