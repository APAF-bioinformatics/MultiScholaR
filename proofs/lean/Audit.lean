import MultiScholaRProofs

/-!
# KIT-193 audit (sections 5.1, 5.2)

For every proved item: `example : <name>_stmt := <theorem>` (a theorem whose type is not its frozen statement does not
compile) and `#print axioms <theorem>` (only `propext`, `Classical.choice` and `Quot.sound` are allowed). Witness
theorems are listed the same way. Run from `proofs/lean` with `lake env lean Audit.lean`; the record job checks
the printed axioms (`.github/workflows/kit-193-lean.yml`, step "axioms").
-/

open MultiScholaRProofs
