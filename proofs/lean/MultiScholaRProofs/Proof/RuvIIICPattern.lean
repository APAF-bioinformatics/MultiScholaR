module

public import MultiScholaRProofs.Spec.RuvIIICPattern
import Mathlib.Tactic

/-!
# Proofs for group R, R-J3: RUV-III-C's missing-value pattern (KIT-193 section 7.2)

Theorems whose types are the frozen statements of `Spec/RuvIIICPattern.lean` (frozen at 2ab4c189): `J3_pattern`,
`J3_counts` (specification lemmas: the content is the definition `out`, read back), and the witness `J3_witness`
(four runs in one set, `k = 1`, a feature observed on three runs with its control complete there, a feature
observed on one run).
-/

public section

namespace MultiScholaRProofs.Ruv

open Matrix Finset

/-- R-J3: the pattern. -/
theorem J3_pattern : J3_pattern_stmt := by
  intro n p s k r j
  refine ⟨fun h i => ?_, fun h i => ?_⟩
  · cases hobs : r.obs i j <;> simp [RuvIIICBr.out, h, hobs]
  · simp [RuvIIICBr.out, h]

/-- R-J3: the missing values the clean-up meets. -/
theorem J3_counts : J3_counts_stmt := by
  intro n p s k r j g ⟨i0, hi0⟩
  have hc : r.correctable j = true := by
    by_contra hc
    simp only [Bool.not_eq_true] at hc
    simp [RuvIIICBr.out, hc] at hi0
  congr 1
  apply Finset.filter_congr
  intro i _
  cases hobs : r.obs i j <;> simp [RuvIIICBr.out, hc, hobs]

/-- The witness instance: four runs in one set, `k = 1`, two features, the first the only control. Feature 0 is
observed on runs 0 to 2 (residual dimension 3 - 1 = 2 > k, the general branch, the eigen step succeeding); feature 1
on run 0 only (residual dimension 0 < k: not corrected). The residual dimension is the code's (`:166`): observed runs
minus the replicate sets with an observed run (one set here). -/
def witnessJ3 : RuvIIICBr 4 2 1 1 where
  Y := 0
  obs := fun i j => if j = 0 then decide (i ≠ 3) else decide (i = 0)
  M := Matrix.of fun _ _ => 1
  ctl := {0}
  proj := fun _ => 1
  solver := fun _ => 0
  resDim := fun R => R.card - 1
  eigenOK := fun _ => true
  closedForm := fun _ _ => 0

/-- Witness for `J3_counts_stmt`. -/
theorem J3_witness : J3_witness_stmt :=
  ⟨4, 2, 1, 1, witnessJ3, 0, 1, by norm_num, by decide, by decide, by decide, ⟨3, by decide⟩, ⟨0, by decide⟩⟩

end MultiScholaRProofs.Ruv
