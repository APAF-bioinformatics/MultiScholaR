module

public import MultiScholaRProofs.Spec.FilterCount

/-! KIT-193 lifecycle L0 (branch `kit-193-lean-lifecycle` only): a trivial theorem the record job must pass. -/

public section

namespace MultiScholaRProofs.Lifecycle

theorem l0_trivial : (1 : Nat) + 1 = 2 := rfl

/-- L0 push 2: a planted gap; the build must fail under `--wfail`. -/
theorem l0_planted_gap : (1 : Nat) = 2 := sorry

end MultiScholaRProofs.Lifecycle
