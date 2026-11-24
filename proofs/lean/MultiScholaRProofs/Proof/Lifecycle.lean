module

public import MultiScholaRProofs.Spec.FilterCount

/-! KIT-193 lifecycle L0 (branch `kit-193-lean-lifecycle` only): a trivial theorem the record job must pass. -/

public section

namespace MultiScholaRProofs.Lifecycle

theorem l0_trivial : (1 : Nat) + 1 = 2 := rfl

end MultiScholaRProofs.Lifecycle
