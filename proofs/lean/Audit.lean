import MultiScholaRProofs

/-!
# KIT-193 audit (sections 5.1, 5.2)

For every proved item: `example : <name>_stmt := <theorem>` (a theorem whose type is not its frozen statement does not
compile) and `#print axioms <theorem>` (only `propext`, `Classical.choice` and `Quot.sound` are allowed). Witness
theorems are listed the same way. Run from `proofs/lean` with `lake env lean Audit.lean`; the record job checks
the printed axioms (`.github/workflows/kit-193-lean.yml`, step "axioms").
-/

open MultiScholaRProofs

/-! ## Target 3: Gene Ontology propagation (`Spec/GoClosure.lean`) -/

example : GoClosure.bfs_terminates_stmt := GoClosure.bfs_terminates
example : GoClosure.bfs_correct_stmt := GoClosure.bfs_correct
example : GoClosure.propagate_fin_stmt := GoClosure.propagate_fin
example : GoClosure.idempotent_stmt := GoClosure.idempotent
example : GoClosure.monotone_stmt := GoClosure.monotone
example : GoClosure.extensive_kept_stmt := GoClosure.extensive_kept
example : GoClosure.not_extensive_stmt := GoClosure.not_extensive
example : GoClosure.labels_mono_stmt := GoClosure.labels_mono
example : GoClosure.root_removed_stmt := GoClosure.root_removed
example : GoClosure.all_absent_stmt := GoClosure.all_absent
example : GoClosure.bfs_all_stmt := GoClosure.bfs_all
example : GoClosure.aspect_closed_stmt := GoClosure.aspect_closed
example : GoClosure.monotone_witness_stmt := GoClosure.monotone_witness
example : GoClosure.labels_mono_witness_stmt := GoClosure.labels_mono_witness
example : GoClosure.all_absent_witness_stmt := GoClosure.all_absent_witness
example : GoClosure.aspect_closed_witness_stmt := GoClosure.aspect_closed_witness

#print axioms GoClosure.bfs_terminates
#print axioms GoClosure.bfs_correct
#print axioms GoClosure.propagate_fin
#print axioms GoClosure.idempotent
#print axioms GoClosure.monotone
#print axioms GoClosure.extensive_kept
#print axioms GoClosure.not_extensive
#print axioms GoClosure.labels_mono
#print axioms GoClosure.root_removed
#print axioms GoClosure.all_absent
#print axioms GoClosure.bfs_all
#print axioms GoClosure.aspect_closed
#print axioms GoClosure.monotone_witness
#print axioms GoClosure.labels_mono_witness
#print axioms GoClosure.all_absent_witness
#print axioms GoClosure.aspect_closed_witness
