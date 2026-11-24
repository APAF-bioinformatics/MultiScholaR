module

public import MultiScholaRProofs.Spec.Permutation
public import MultiScholaRProofs.Spec.RuvAlgebra

/-!
# Two sharper witness statements (KIT-193 proofs audit P03)

Source: `tickets/kit/audits/AUDIT-KIT-193-proofs-20261124.md`, P03, adopted by the coordinator on 2026-11-24 (an
engineering decision). The witnesses `perm_measure_witness_stmt` (`Spec/Permutation.lean`) and `G0_witness_stmt`
(`Spec/RuvDownstream.lean`) stay frozen and proved; these two ask for instances away from the corner cases those
allow: a law in which a relabelling moves an atom and changes the statistic, and an RUV-III instance with residual
degrees of freedom left, a feature that is not a control, and corrected values whose residuals are not zero. The two
definitions are written exactly as the register gives them. This file is new and has its own line in `SPEC.sha256`;
no other Spec file changes. It is reviewed at the next statements review.
-/

@[expose] public section

namespace MultiScholaRProofs.Permutation

open MeasureTheory

/-- A witness for `perm_measure_stmt` in which a relabelling moves an atom and changes the statistic. -/
def perm_measure_moving_witness_stmt : Prop :=
  ∃ (μ : Measure (Fin 4 → ℝ)) (T : (Fin 4 → ℝ) → ℝ),
    IsProbabilityMeasure μ ∧ Measurable T ∧ (∀ σ : Equiv.Perm (Fin 4), μ.map (fun f => f ∘ σ) = μ) ∧
    ∃ (f : Fin 4 → ℝ) (σ : Equiv.Perm (Fin 4)), 0 < μ {f} ∧ f ∘ σ ≠ f ∧ T (f ∘ σ) ≠ T f

end MultiScholaRProofs.Permutation

namespace MultiScholaRProofs.Ruv

open Matrix RuvIII

/-- A witness for `G0_stmt` away from the boundary: residual degrees of freedom left, a feature that is not a
control, and corrected values whose residuals are not zero. -/
def G0_interior_witness_stmt : Prop :=
  ∃ (n p s k g : ℕ) (r : RuvIII n p s k) (X : Matrix (Fin n) (Fin g) ℝ),
    0 < k ∧ 0 < g ∧ g + k < n ∧ r.ctl ≠ Finset.univ ∧ Hset r ∧ Hac r ∧ HU r ∧ Hone r ∧
    residop r.M * X = 0 ∧ IsUnit (Xᵀ * X).det ∧ residop X * r.Ystar ≠ 0

end MultiScholaRProofs.Ruv
