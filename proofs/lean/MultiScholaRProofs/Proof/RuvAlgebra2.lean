module

public import MultiScholaRProofs.Spec.RuvAlgebra
public import MultiScholaRProofs.Proof.RuvAlgebra
import Mathlib.Tactic

/-!
# Proofs for group R, the RUV-III algebra (KIT-193 section 7.2), priority 2

Theorems whose types are the frozen statements of `Spec/RuvAlgebra.lean`: R-B1 (per-feature level shifts), R-B2
(after restoring levels), R-D1 (basis invariance), R-D2 (scaling), R-E1 (offset), R-E2 (the amounts differ, an
instance), R-K1 (`K1`, `K1_gls`), with the witnesses `B1_witness`, `D_witness`, `K1_witness`, `K1_gls_witness`.
-/

public section

namespace MultiScholaRProofs.Ruv

open Matrix RuvIII

variable {n p s k : ℕ}

/-! ## Level shifts and offsets (R-B1, R-B2, R-E1) -/

/-- The data shifted by `1c'` (every run of feature `j` shifted by `c j`). -/
@[expose] def shifted (r : RuvIII n p s k) (c : Fin p → ℝ) : RuvIII n p s k :=
  { r with Y := r.Y + Matrix.of fun _ j => c j }

/-- `U'1 = 0` from `R_M U = U` and `R_M 1 = 0`. -/
theorem Ut_one (r : RuvIII n p s k) (hU : HU r) (hone : Hone r) : r.Uᵀ *ᵥ (fun _ => (1 : ℝ)) = 0 := by
  unfold Hone at hone
  unfold HU at hU
  rw [← hU, transpose_mul, residop_transpose, ← Matrix.mulVec_mulVec, hone, Matrix.mulVec_zero]

/-- A matrix whose every row is `c'` is killed on the left by any `A` with `A 1 = 0`. -/
theorem mul_const_rows {m : Type} [Fintype m] {q : Type} (A : Matrix m (Fin n) ℝ) (c : q → ℝ)
    (h : A *ᵥ (fun _ => (1 : ℝ)) = 0) : A * (Matrix.of fun (_ : Fin n) j => c j) = 0 := by
  ext a j
  have h1 := congrFun h a
  simp only [Matrix.mulVec, dotProduct, mul_one, Pi.zero_apply] at h1
  simp only [Matrix.mul_apply, Matrix.of_apply, Matrix.zero_apply]
  rw [← Finset.sum_mul, h1, zero_mul]

theorem shifted_residop (r : RuvIII n p s k) (c : Fin p → ℝ) (hone : Hone r) :
    residop r.M * (shifted r c).Y = residop r.M * r.Y := by
  show residop r.M * (r.Y + _) = _
  rw [Matrix.mul_add, mul_const_rows _ _ hone, add_zero]

theorem shifted_alpha (r : RuvIII n p s k) (c : Fin p → ℝ) (hU : HU r) (hone : Hone r) :
    (shifted r c).alpha = r.alpha := by
  show r.Uᵀ * (r.Y + _) = r.Uᵀ * r.Y
  rw [Matrix.mul_add, mul_const_rows _ _ (Ut_one r hU hone), add_zero]

theorem shifted_ac (r : RuvIII n p s k) (c : Fin p → ℝ) (hU : HU r) (hone : Hone r) :
    (shifted r c).ac = r.ac := by
  ext a l
  show (shifted r c).alpha a l = r.alpha a l
  rw [shifted_alpha r c hU hone]

theorem shifted_Yc (r : RuvIII n p s k) (c : Fin p → ℝ) :
    (shifted r c).Yc = (r.Yc + Matrix.of fun _ (l : r.ctl) => c l : Matrix (Fin n) r.ctl ℝ) := by
  ext i l
  rfl

/-- The shifted corrected values: `Y* + 1 d'` with `d_j = c_j - Σ_{l ∈ ctl} c_l D_{lj}`. -/
theorem shifted_Ystar (r : RuvIII n p s k) (c : Fin p → ℝ) (hU : HU r) (hone : Hone r) :
    (shifted r c).Ystar = r.Ystar + Matrix.of fun _ j =>
      c j - ∑ l : r.ctl, c l * (r.acᵀ * (r.ac * r.acᵀ)⁻¹ * r.alpha) l j := by
  have hW : (shifted r c).W * (shifted r c).alpha =
      r.W * r.alpha + (Matrix.of fun _ (l : r.ctl) => c l) * (r.acᵀ * (r.ac * r.acᵀ)⁻¹ * r.alpha) := by
    rw [shifted_alpha r c hU hone]
    have h1 : (shifted r c).W = r.W + (Matrix.of fun _ (l : r.ctl) => c l) * r.acᵀ * (r.ac * r.acᵀ)⁻¹ := by
      show (shifted r c).Yc * (shifted r c).acᵀ * ((shifted r c).ac * (shifted r c).acᵀ)⁻¹ = _
      rw [shifted_ac r c hU hone]
      show (r.Yc + Matrix.of fun _ (l : r.ctl) => c l) * r.acᵀ * (r.ac * r.acᵀ)⁻¹ = _
      rw [Matrix.add_mul, Matrix.add_mul]
      rfl
    rw [h1, Matrix.add_mul]
    simp only [Matrix.mul_assoc]
  unfold RuvIII.Ystar
  rw [hW]
  ext i j
  simp only [shifted, Matrix.sub_apply, Matrix.add_apply, Matrix.of_apply, Matrix.mul_apply]
  ring

/-- Column means of `A + 1e'`. -/
theorem colMean_add_rows {m : Type} (A : Matrix (Fin n) m ℝ) (e : m → ℝ) (hn : (n : ℝ) ≠ 0) (j : m) :
    colMean (A + Matrix.of fun _ l => e l) j = colMean A j + e j := by
  simp only [colMean, Matrix.add_apply, Matrix.of_apply, Finset.sum_add_distrib, Finset.sum_const,
    Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
  field_simp

/-- After restoring levels, a shift `1c'` of the data shifts the output by exactly `1c'`. -/
theorem shifted_restore (r : RuvIII n p s k) (c : Fin p → ℝ) (hU : HU r) (hone : Hone r) :
    restore (shifted r c).Ystar (shifted r c).Y = restore r.Ystar r.Y + Matrix.of fun _ j => c j := by
  ext i j
  have hn : (n : ℝ) ≠ 0 := by
    have := i.pos
    exact_mod_cast this.ne'
  rw [shifted_Ystar r c hU hone]
  show (restore (r.Ystar + _) (r.Y + _)) i j = _
  simp only [restore, Matrix.add_apply, Matrix.of_apply, colMean_add_rows _ _ hn]
  ring

/-- Contrasts (`λ'1 = 0`) are unchanged by adding `1d'`. -/
theorem vecMul_add_rows (A : Matrix (Fin n) (Fin p) ℝ) (d : Fin p → ℝ) (lam : Fin n → ℝ) (h : ∑ i, lam i = 0) :
    lam ᵥ* (A + Matrix.of fun _ j => d j) = lam ᵥ* A := by
  rw [Matrix.vecMul_add]
  have h0 : lam ᵥ* (Matrix.of fun (_ : Fin n) j => d j) = 0 := by
    funext j
    simp only [Matrix.vecMul, dotProduct, Matrix.of_apply, Pi.zero_apply]
    rw [← Finset.sum_mul, h, zero_mul]
  rw [h0, add_zero]

/-- R-B1: per-feature level shifts. -/
theorem B1 : B1_stmt := by
  intro n p s k r c _ _ hU hone
  refine ⟨shifted_residop r c hone, shifted_alpha r c hU hone, fun i j => ?_, fun lam hlam => ?_⟩
  · show (shifted r c).Ystar i j - r.Ystar i j = _
    rw [shifted_Ystar r c hU hone]
    simp only [Matrix.add_apply, Matrix.of_apply]
    ring
  · show lam ᵥ* (shifted r c).Ystar = lam ᵥ* r.Ystar
    rw [shifted_Ystar r c hU hone]
    exact vecMul_add_rows _ _ lam hlam

/-- R-B2: after restoring levels, fully equivariant. -/
theorem B2 : B2_stmt := by
  intro n p s k r c _ _ hU hone
  exact shifted_restore r c hU hone

/-- R-E1: an offset of every value. -/
theorem E1 : E1_stmt := by
  intro n p s k r c _ _ hU hone
  refine ⟨shifted_residop r (fun _ => c) hone, shifted_alpha r (fun _ => c) hU hone, fun i j => ?_, fun i j => ?_,
    fun lam hlam => ?_⟩
  · show (shifted r (fun _ => c)).Ystar i j - r.Ystar i j = _
    rw [shifted_Ystar r _ hU hone]
    simp only [Matrix.add_apply, Matrix.of_apply, offsetAmount, mul_sub, mul_one, Finset.mul_sum]
    ring
  · show restore (shifted r (fun _ => c)).Ystar (shifted r (fun _ => c)).Y i j - restore r.Ystar r.Y i j = c
    rw [shifted_restore r _ hU hone]
    simp
  · show lam ᵥ* (shifted r (fun _ => c)).Ystar = lam ᵥ* r.Ystar
    rw [shifted_Ystar r _ hU hone]
    exact vecMul_add_rows _ _ lam hlam

/-! ## Basis invariance and scaling (R-D1, R-D2) -/

/-- Multiplying on the left commutes with taking the control columns. -/
theorem mul_submatrix_cols {a b q : Type} [Fintype b] {q' : Type} (A : Matrix a b ℝ) (B : Matrix b q ℝ)
    (f : q' → q) : (A * B).submatrix id f = A * B.submatrix id f := by
  ext i j
  simp [Matrix.mul_apply]

/-- R-D1: replacing `U` by `UQ` leaves `Y*` unchanged. -/
theorem D1 : D1_stmt := by
  intro n p s k r Q hQ _
  have hQt : IsUnit Qᵀ.det := by rwa [Matrix.det_transpose]
  set r' : RuvIII n p s k := { r with U := r.U * Q } with hr'
  have hal : r'.alpha = Qᵀ * r.alpha := by
    show (r.U * Q)ᵀ * r.Y = Qᵀ * (r.Uᵀ * r.Y)
    rw [transpose_mul, Matrix.mul_assoc]
  have hac : r'.ac = Qᵀ * r.ac := by
    show r'.alpha.submatrix id _ = Qᵀ * r.alpha.submatrix id _
    rw [hal, mul_submatrix_cols]
  have hYc : r'.Yc = r.Yc := rfl
  have hW : r'.W = r.W * (Qᵀ)⁻¹ := by
    unfold RuvIII.W
    rw [hac, hYc, transpose_mul, transpose_transpose]
    have : Qᵀ * r.ac * (r.acᵀ * Q) = Qᵀ * (r.ac * r.acᵀ) * Q := by simp only [Matrix.mul_assoc]
    rw [this, Matrix.mul_inv_rev, Matrix.mul_inv_rev]
    simp only [Matrix.mul_assoc]
    rw [← Matrix.mul_assoc Q, Matrix.mul_nonsing_inv _ hQ, Matrix.one_mul]
  unfold RuvIII.Ystar
  rw [hW, hal]
  show r.Y - r.W * (Qᵀ)⁻¹ * (Qᵀ * r.alpha) = r.Y - r.W * r.alpha
  rw [Matrix.mul_assoc, ← Matrix.mul_assoc (Qᵀ)⁻¹, Matrix.nonsing_inv_mul _ hQt, Matrix.one_mul]

/-- R-D2: scaling the data by `t ≠ 0`, with the solver's basis `UQ` of the same subspace. -/
theorem D2 : D2_stmt := by
  intro n p s k r t Q ht hQ hac0
  unfold Hac at hac0
  have hQt : IsUnit Qᵀ.det := by rwa [Matrix.det_transpose]
  have htt : (t * t)⁻¹ * (t * t) = 1 := inv_mul_cancel₀ (mul_ne_zero ht ht)
  intro r'
  have hal : r'.alpha = t • (Qᵀ * r.alpha) := by
    show (r.U * Q)ᵀ * (t • r.Y) = t • (Qᵀ * (r.Uᵀ * r.Y))
    rw [transpose_mul, Matrix.mul_smul, Matrix.mul_assoc]
  have hac : r'.ac = t • (Qᵀ * r.ac) := by
    show r'.alpha.submatrix id _ = t • (Qᵀ * r.alpha.submatrix id _)
    rw [hal, ← mul_submatrix_cols]
    rfl
  have hYc : r'.Yc = t • r.Yc := rfl
  have hact : r'.acᵀ = t • (r.acᵀ * Q) := by
    rw [hac, transpose_smul, transpose_mul, transpose_transpose]
  have hG : r'.ac * r'.acᵀ = (t * t) • (Qᵀ * (r.ac * r.acᵀ) * Q) := by
    rw [hact, hac, Matrix.smul_mul, Matrix.mul_smul, smul_smul]
    simp only [Matrix.mul_assoc]
  have hGinv : (r'.ac * r'.acᵀ)⁻¹ = (t * t)⁻¹ • (Q⁻¹ * (r.ac * r.acᵀ)⁻¹ * (Qᵀ)⁻¹) := by
    rw [hG]
    apply Matrix.inv_eq_left_inv
    rw [Matrix.smul_mul, Matrix.mul_smul, smul_smul, htt, one_smul]
    calc Q⁻¹ * (r.ac * r.acᵀ)⁻¹ * (Qᵀ)⁻¹ * (Qᵀ * (r.ac * r.acᵀ) * Q)
        = Q⁻¹ * ((r.ac * r.acᵀ)⁻¹ * ((Qᵀ)⁻¹ * Qᵀ) * (r.ac * r.acᵀ)) * Q := by
          simp only [Matrix.mul_assoc]
      _ = 1 := by
          rw [Matrix.nonsing_inv_mul _ hQt, Matrix.mul_one, Matrix.nonsing_inv_mul _ hac0, Matrix.mul_one,
            Matrix.nonsing_inv_mul _ hQ]
  have hW : r'.W = r.W * (Qᵀ)⁻¹ := by
    unfold RuvIII.W
    rw [hGinv, hact, hYc]
    simp only [Matrix.smul_mul, Matrix.mul_smul, smul_smul]
    rw [mul_comm t t, htt, one_smul]
    calc r.Yc * (r.acᵀ * Q) * (Q⁻¹ * (r.ac * r.acᵀ)⁻¹ * (Qᵀ)⁻¹)
        = r.Yc * r.acᵀ * (Q * Q⁻¹) * (r.ac * r.acᵀ)⁻¹ * (Qᵀ)⁻¹ := by simp only [Matrix.mul_assoc]
      _ = r.Yc * r.acᵀ * (r.ac * r.acᵀ)⁻¹ * (Qᵀ)⁻¹ := by rw [Matrix.mul_nonsing_inv _ hQ, Matrix.mul_one]
  refine ⟨?_, ?_⟩
  · unfold RuvIII.Ystar
    rw [hW, hal]
    show t • r.Y - r.W * (Qᵀ)⁻¹ * (t • (Qᵀ * r.alpha)) = t • (r.Y - r.W * r.alpha)
    rw [Matrix.mul_smul, Matrix.mul_assoc, ← Matrix.mul_assoc (Qᵀ)⁻¹, Matrix.nonsing_inv_mul _ hQt,
      Matrix.one_mul, smul_sub]
  · rw [hW, Matrix.mul_assoc, Matrix.nonsing_inv_mul _ hQt, Matrix.mul_one]

/-! ## k equal to the residual dimension (R-K1) -/

/-- R-K1 (i): `R_M Y* = 0` when `U U' = R_M`. -/
theorem K1 : K1_stmt := by
  intro n p s k r _ hac h
  rw [← h, Matrix.mul_assoc, Ut_Ystar r hac, Matrix.mul_zero]

/-- A real matrix with `A'A = 0` is zero. -/
theorem eq_zero_of_gram_zero {a b : Type} [Fintype a] (A : Matrix a b ℝ) (h : Aᵀ * A = 0) : A = 0 := by
  ext i j
  have h2 : (fun l => A l j) ⬝ᵥ (fun l => A l j) = 0 := by
    have := congrFun (congrFun h j) j
    simpa [Matrix.mul_apply, dotProduct] using this
  exact congrFun (dotProduct_self_eq_zero.1 h2) i

/-- R-K1 (ii): the closed form. -/
theorem K1_gls : K1_gls_stmt := by
  intro n p s k r hset hUU hS hMSM hac
  unfold Hac at hac
  unfold Hset at hset
  have hact : r.acᵀ = r.Ycᵀ * r.U := by rw [← Ut_Yc r, transpose_mul, transpose_transpose]
  -- `Y* = (1 - S U G⁻¹ U') Y`, `S = Y_c Y_c'`, `G = a_c a_c'`
  have hYstar : r.Ystar = (1 - r.Yc * r.Ycᵀ * r.U * (r.ac * r.acᵀ)⁻¹ * r.Uᵀ) * r.Y := by
    have hW : r.W = r.Yc * r.Ycᵀ * r.U * (r.ac * r.acᵀ)⁻¹ := by
      unfold RuvIII.W
      rw [show r.Yc * r.acᵀ = r.Yc * r.Ycᵀ * r.U by rw [hact, Matrix.mul_assoc]]
    unfold RuvIII.Ystar
    rw [hW, Matrix.sub_mul, Matrix.one_mul]
    unfold RuvIII.alpha
    simp only [Matrix.mul_assoc]
  -- `U'SU = G`
  have hUSU : r.Uᵀ * (r.Yc * r.Ycᵀ) * r.U = r.ac * r.acᵀ := by
    rw [hact, ← Ut_Yc r]
    simp only [Matrix.mul_assoc]
  set S := r.Yc * r.Ycᵀ with hSdef
  set G := r.ac * r.acᵀ with hGdef
  set P := r.M * (r.Mᵀ * r.M)⁻¹ * r.Mᵀ with hP
  have hUU' : r.U * r.Uᵀ = 1 - P := hUU
  have hRM : (1 - P) * r.M = 0 := by
    rw [hP, Matrix.sub_mul, Matrix.one_mul, Matrix.mul_assoc, Matrix.mul_assoc, Matrix.nonsing_inv_mul _ hset,
      Matrix.mul_one, sub_self]
  have hUM : r.Uᵀ * r.M = 0 := by
    apply eq_zero_of_gram_zero
    rw [transpose_mul, transpose_transpose, Matrix.mul_assoc, ← Matrix.mul_assoc r.U, hUU', hRM,
      Matrix.mul_zero]
  have hMU : r.Mᵀ * r.U = 0 := by
    rw [← transpose_transpose (r.Mᵀ * r.U), transpose_mul, transpose_transpose, hUM, transpose_zero]
  set H := r.M * (r.Mᵀ * S⁻¹ * r.M)⁻¹ * (r.Mᵀ * S⁻¹) with hH
  set K := 1 - S * r.U * G⁻¹ * r.Uᵀ - H with hK
  have hKM : K * r.M = 0 := by
    rw [hK, hH, Matrix.sub_mul, Matrix.sub_mul, Matrix.one_mul]
    have h1 : S * r.U * G⁻¹ * r.Uᵀ * r.M = 0 := by rw [Matrix.mul_assoc, hUM, Matrix.mul_zero]
    have h2 : r.M * (r.Mᵀ * S⁻¹ * r.M)⁻¹ * (r.Mᵀ * S⁻¹) * r.M = r.M := by
      rw [Matrix.mul_assoc, Matrix.mul_assoc, Matrix.nonsing_inv_mul _ hMSM, Matrix.mul_one]
    rw [h1, h2, sub_zero, sub_self]
  have hKSU : K * (S * r.U) = 0 := by
    rw [hK, hH, Matrix.sub_mul, Matrix.sub_mul, Matrix.one_mul]
    have h1 : S * r.U * G⁻¹ * r.Uᵀ * (S * r.U) = S * r.U := by
      have : r.Uᵀ * (S * r.U) = G := by rw [← Matrix.mul_assoc, hUSU]
      rw [Matrix.mul_assoc, this, Matrix.mul_assoc (S * r.U), Matrix.nonsing_inv_mul _ hac, Matrix.mul_one]
    have h2 : r.M * (r.Mᵀ * S⁻¹ * r.M)⁻¹ * (r.Mᵀ * S⁻¹) * (S * r.U) = 0 := by
      have : r.Mᵀ * S⁻¹ * (S * r.U) = 0 := by
        rw [Matrix.mul_assoc, ← Matrix.mul_assoc S⁻¹, Matrix.nonsing_inv_mul _ hS, Matrix.one_mul, hMU]
      rw [Matrix.mul_assoc, this, Matrix.mul_zero]
    rw [h1, h2, sub_self, sub_zero]
  -- `1 = M A + S U B`
  have hdec : (1 : Matrix (Fin n) (Fin n) ℝ) =
      r.M * ((r.Mᵀ * r.M)⁻¹ * r.Mᵀ * (1 - S * r.U * (G⁻¹ * r.Uᵀ))) + S * r.U * (G⁻¹ * r.Uᵀ) := by
    have hP' : P = 1 - r.U * r.Uᵀ := by rw [hUU']; abel
    have hMA : r.M * ((r.Mᵀ * r.M)⁻¹ * r.Mᵀ * (1 - S * r.U * (G⁻¹ * r.Uᵀ))) =
        P * (1 - S * r.U * (G⁻¹ * r.Uᵀ)) := by
      rw [hP]
      simp only [Matrix.mul_assoc]
    have h3 : r.Uᵀ * (S * r.U * (G⁻¹ * r.Uᵀ)) = r.Uᵀ := by
      calc r.Uᵀ * (S * r.U * (G⁻¹ * r.Uᵀ)) = (r.Uᵀ * S * r.U * G⁻¹) * r.Uᵀ := by simp only [Matrix.mul_assoc]
        _ = r.Uᵀ := by rw [hUSU, Matrix.mul_nonsing_inv _ hac, Matrix.one_mul]
    have h4 : r.U * r.Uᵀ * (S * r.U * (G⁻¹ * r.Uᵀ)) = r.U * r.Uᵀ := by
      rw [Matrix.mul_assoc, h3]
    rw [hMA, hP', Matrix.sub_mul, Matrix.one_mul, Matrix.mul_sub, Matrix.mul_one, h4]
    abel
  have hK0 : K = 0 := by
    rw [← Matrix.mul_one K, hdec, Matrix.mul_add, ← Matrix.mul_assoc, hKM, Matrix.zero_mul, ← Matrix.mul_assoc,
      hKSU, Matrix.zero_mul, add_zero]
  have hEq : 1 - S * r.U * G⁻¹ * r.Uᵀ = H := by
    rw [← sub_eq_zero]
    exact hK0
  rw [hYstar, hEq, hH]

/-! ## The instance of R-E2 and the witnesses -/

/-- Four runs in one set, three features, the first the only control, `k = 1`, `U = (1, -1, 0, 0)'`. -/
@[expose] def witnessE2 : RuvIII 4 3 1 1 where
  Y := Matrix.of fun i j => if i = 0 ∧ j = 0 then 1 else 0
  M := Matrix.of fun _ _ => 1
  ctl := {0}
  U := Matrix.of fun i _ => if i = 0 then 1 else if i = 1 then -1 else 0

theorem witnessE2_mem0 : (0 : Fin 3) ∈ witnessE2.ctl := by simp [witnessE2]

theorem witnessE2_sum_ctl (f : witnessE2.ctl → ℝ) : ∑ c, f c = f ⟨0, witnessE2_mem0⟩ := by
  have hu : ∀ c : witnessE2.ctl, c = ⟨0, witnessE2_mem0⟩ := fun c =>
    Subtype.ext (by simpa [witnessE2] using c.2)
  rw [Finset.sum_eq_single (⟨0, witnessE2_mem0⟩ : witnessE2.ctl) (fun c _ hc => absurd (hu c) hc)
    (fun h => absurd (Finset.mem_univ _) h)]

theorem witnessE2_MtM : witnessE2.Mᵀ * witnessE2.M = Matrix.of fun _ _ => (4 : ℝ) := by
  ext a b
  fin_cases a
  fin_cases b
  simp [witnessE2, Matrix.mul_apply]

theorem witnessE2_MtM_inv : (witnessE2.Mᵀ * witnessE2.M)⁻¹ = Matrix.of fun _ _ => (1 / 4 : ℝ) := by
  apply Matrix.inv_eq_left_inv
  rw [witnessE2_MtM]
  ext a b
  fin_cases a
  fin_cases b
  simp [Matrix.mul_apply]

theorem witnessE2_residop : residop witnessE2.M =
    Matrix.of fun i j => if i = j then (3 / 4 : ℝ) else -1 / 4 := by
  unfold residop
  rw [witnessE2_MtM_inv]
  ext i j
  fin_cases i <;> fin_cases j <;> simp [witnessE2, Matrix.mul_apply] <;> norm_num

theorem witnessE2_alpha : witnessE2.alpha = Matrix.of fun _ j => if j = 0 then (1 : ℝ) else 0 := by
  ext a j
  fin_cases a
  fin_cases j <;> simp [witnessE2, RuvIII.alpha, Matrix.mul_apply]

theorem witnessE2_acac : witnessE2.ac * witnessE2.acᵀ = 1 := by
  ext a b
  fin_cases a
  fin_cases b
  rw [Matrix.mul_apply]
  simp only [RuvIII.ac, Matrix.submatrix_apply, Matrix.transpose_apply, id]
  rw [witnessE2_sum_ctl]
  simp [witnessE2_alpha]

theorem witnessE2_hset : Hset witnessE2 := by
  unfold Hset
  rw [witnessE2_MtM, Matrix.det_unique]
  simp

theorem witnessE2_hac : Hac witnessE2 := by
  unfold Hac
  rw [witnessE2_acac, Matrix.det_one]
  exact isUnit_one

theorem witnessE2_hU : HU witnessE2 := by
  unfold HU
  rw [witnessE2_residop]
  ext i j
  fin_cases i <;> fin_cases j <;> simp [witnessE2, Matrix.mul_apply, Fin.sum_univ_four] <;> norm_num

theorem witnessE2_hone : Hone witnessE2 := by
  unfold Hone
  rw [witnessE2_residop]
  ext i
  fin_cases i <;> simp [Matrix.mulVec, dotProduct, Fin.sum_univ_four] <;> norm_num

theorem witnessE2_offset (j : Fin 3) : witnessE2.offsetAmount j = if j = 0 then 0 else 1 := by
  unfold offsetAmount
  rw [witnessE2_acac, inv_one, Matrix.mul_one, witnessE2_sum_ctl]
  rw [Matrix.mul_apply]
  simp only [RuvIII.ac, Matrix.submatrix_apply, Matrix.transpose_apply, id]
  rw [witnessE2_alpha]
  fin_cases j <;> simp

/-- R-E2: the offset amounts differ (features 0 and 1 of `witnessE2`: 0 and 1). -/
theorem E2 : E2_stmt := by
  refine ⟨1, witnessE2, witnessE2_hset, witnessE2_hac, witnessE2_hU, witnessE2_hone, 0, 1, ?_⟩
  rw [witnessE2_offset, witnessE2_offset]
  simp

/-- Witness for `B1_stmt`, `B2_stmt`, `E1_stmt`: `witnessE2` (`k = 1`, one set, four runs). -/
theorem B1_witness : B1_witness_stmt :=
  ⟨4, 3, 1, 1, witnessE2, by norm_num, by norm_num, witnessE2_hset, witnessE2_hac, witnessE2_hU, witnessE2_hone⟩

/-- Witness for `D1_stmt`, `D2_stmt`: `witnessE2`, `Q = (2)`, `t = 2`. -/
theorem D_witness : D_witness_stmt := by
  refine ⟨4, 3, 1, 1, witnessE2, Matrix.of fun _ _ => 2, 2, by norm_num, witnessE2_hac, ?_, ?_, by norm_num,
    by norm_num⟩
  · rw [Matrix.det_unique]
    simp
  · intro h
    have h1 := congrFun (congrFun h 0) 0
    simp at h1

/-! ## The instance of R-K1: four runs in one set, `k = 3`, four controls (`Y = I`), `U` from a Hadamard matrix -/

/-- `U`'s columns `(1, -1, 1, -1)/2`, `(1, 1, -1, -1)/2`, `(1, -1, -1, 1)/2`: an orthonormal basis of the
contrasts of four runs, so `U U' = I - 11'/4 = R_M`. -/
@[expose] noncomputable def hadamardU : Matrix (Fin 4) (Fin 3) ℝ :=
  !![1 / 2, 1 / 2, 1 / 2; -1 / 2, 1 / 2, -1 / 2; 1 / 2, -1 / 2, -1 / 2; -1 / 2, -1 / 2, 1 / 2]

@[expose] noncomputable def witnessK1 : RuvIII 4 4 1 3 where
  Y := 1
  M := Matrix.of fun _ _ => 1
  ctl := Finset.univ
  U := hadamardU

theorem witnessK1_sum_ctl (f : witnessK1.ctl → ℝ) :
    ∑ c, f c = ∑ c : Fin 4, f ⟨c, Finset.mem_univ c⟩ :=
  Fintype.sum_equiv (Equiv.subtypeUnivEquiv fun c => Finset.mem_univ c) _ _ fun _ => rfl

theorem witnessK1_MtM : witnessK1.Mᵀ * witnessK1.M = Matrix.of fun _ _ => (4 : ℝ) := by
  ext a b
  fin_cases a
  fin_cases b
  simp [witnessK1, Matrix.mul_apply]

theorem witnessK1_MtM_inv : (witnessK1.Mᵀ * witnessK1.M)⁻¹ = Matrix.of fun _ _ => (1 / 4 : ℝ) := by
  apply Matrix.inv_eq_left_inv
  rw [witnessK1_MtM]
  ext a b
  fin_cases a
  fin_cases b
  simp [Matrix.mul_apply]

theorem witnessK1_residop : residop witnessK1.M =
    Matrix.of fun i j => if i = j then (3 / 4 : ℝ) else -1 / 4 := by
  unfold residop
  rw [witnessK1_MtM_inv]
  ext i j
  fin_cases i <;> fin_cases j <;> simp [witnessK1, Matrix.mul_apply] <;> norm_num

theorem witnessK1_hset : Hset witnessK1 := by
  unfold Hset
  rw [witnessK1_MtM, Matrix.det_unique]
  simp

theorem witnessK1_UU : witnessK1.U * witnessK1.Uᵀ = residop witnessK1.M := by
  rw [witnessK1_residop]
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [witnessK1, hadamardU, Matrix.mul_apply, Fin.sum_univ_three] <;> norm_num

theorem witnessK1_YcYc : witnessK1.Yc * witnessK1.Ycᵀ = 1 := by
  ext i j
  rw [Matrix.mul_apply, witnessK1_sum_ctl]
  simp only [RuvIII.Yc, Matrix.submatrix_apply, Matrix.transpose_apply, witnessK1]
  fin_cases i <;> fin_cases j <;> simp [Matrix.one_apply]

theorem witnessK1_acac : witnessK1.ac * witnessK1.acᵀ = 1 := by
  ext a b
  rw [Matrix.mul_apply, witnessK1_sum_ctl]
  simp only [RuvIII.ac, RuvIII.alpha, Matrix.submatrix_apply, Matrix.transpose_apply, id, witnessK1,
    Matrix.mul_one]
  fin_cases a <;> fin_cases b <;>
    simp [hadamardU, Fin.sum_univ_four] <;> norm_num

theorem witnessK1_hac : Hac witnessK1 := by
  unfold Hac
  rw [witnessK1_acac, Matrix.det_one]
  exact isUnit_one

theorem K1_witness : K1_witness_stmt :=
  ⟨4, 4, 1, 3, witnessK1, by norm_num, by norm_num, witnessK1_hset, witnessK1_hac, witnessK1_UU⟩

theorem K1_gls_witness : K1_gls_witness_stmt := by
  refine ⟨4, 4, 1, 3, witnessK1, by norm_num, by norm_num, witnessK1_hset, witnessK1_UU, ?_, ?_, witnessK1_hac⟩
  · rw [witnessK1_YcYc, Matrix.det_one]
    exact isUnit_one
  · rw [witnessK1_YcYc, inv_one, Matrix.mul_one, witnessK1_MtM, Matrix.det_unique]
    simp

end MultiScholaRProofs.Ruv
