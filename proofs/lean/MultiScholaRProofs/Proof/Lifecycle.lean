module

public import MultiScholaRProofs.Spec.FilterCount

/-! KIT-193 lifecycle L0 (branch `kit-193-lean-lifecycle` only): a trivial theorem the record job must pass. -/

public section

namespace MultiScholaRProofs.Lifecycle

theorem l0_trivial : (1 : Nat) + 1 = 2 := rfl

/-- L0 push 3: a planted new axiom; the axiom audit must fail. -/
axiom l0_planted_axiom : (1 : Nat) = 2

theorem l0_uses_planted_axiom : (1 : Nat) = 2 := l0_planted_axiom

end MultiScholaRProofs.Lifecycle
