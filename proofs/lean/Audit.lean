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

/-! ## Target 2: replay-rule selection (`Spec/ReplayRule.lean`) -/

example : ReplayRule.rule_first_match_stmt := ReplayRule.rule_first_match
example : ReplayRule.totality_stmt := ReplayRule.totality
example : ReplayRule.faithful_stmt := ReplayRule.faithful
example : ReplayRule.without_marker_stmt := ReplayRule.without_marker
example : ReplayRule.managed_new_design_stmt := ReplayRule.managed_new_design
example : ReplayRule.round_trip_stmt := ReplayRule.round_trip
example : ReplayRule.locality_stmt := ReplayRule.locality
example : ReplayRule.unfixed_default_stmt := ReplayRule.unfixed_default
example : ReplayRule.faithful_witness_stmt := ReplayRule.faithful_witness
example : ReplayRule.without_marker_witness_stmt := ReplayRule.without_marker_witness
example : ReplayRule.managed_new_design_witness_stmt := ReplayRule.managed_new_design_witness
example : ReplayRule.unfixed_default_witness_stmt := ReplayRule.unfixed_default_witness

#print axioms ReplayRule.rule_first_match
#print axioms ReplayRule.totality
#print axioms ReplayRule.faithful
#print axioms ReplayRule.without_marker
#print axioms ReplayRule.managed_new_design
#print axioms ReplayRule.round_trip
#print axioms ReplayRule.locality
#print axioms ReplayRule.unfixed_default
#print axioms ReplayRule.faithful_witness
#print axioms ReplayRule.without_marker_witness
#print axioms ReplayRule.managed_new_design_witness
#print axioms ReplayRule.unfixed_default_witness

/-! ## Target 1: the missing-value filter's integer rule (`Spec/FilterCount.lean`) -/

example : FilterCount.never_stricter_stmt := FilterCount.never_stricter
example : FilterCount.exact_lcm_stmt := FilterCount.exact_lcm
example : FilterCount.exact_316_stmt := FilterCount.exact_316
example : FilterCount.balanced_stmt := FilterCount.balanced
example : FilterCount.share_mono_stmt := FilterCount.share_mono
example : FilterCount.peptide_floor_stmt := FilterCount.peptide_floor
example : FilterCount.counterexample_stmt := FilterCount.counterexample
example : FilterCount.counterexample_any_rounding_stmt := FilterCount.counterexample_any_rounding
example : FilterCount.never_stricter_witness_stmt := FilterCount.never_stricter_witness
example : FilterCount.exact_lcm_witness_stmt := FilterCount.exact_lcm_witness
example : FilterCount.exact_316_witness_stmt := FilterCount.exact_316_witness
example : FilterCount.balanced_witness_stmt := FilterCount.balanced_witness
example : FilterCount.share_mono_witness_stmt := FilterCount.share_mono_witness
example : FilterCount.counterexample_any_rounding_witness_stmt := FilterCount.counterexample_any_rounding_witness

#print axioms FilterCount.never_stricter
#print axioms FilterCount.exact_lcm
#print axioms FilterCount.exact_316
#print axioms FilterCount.balanced
#print axioms FilterCount.share_mono
#print axioms FilterCount.peptide_floor
#print axioms FilterCount.counterexample
#print axioms FilterCount.counterexample_any_rounding
#print axioms FilterCount.never_stricter_witness
#print axioms FilterCount.exact_lcm_witness
#print axioms FilterCount.exact_316_witness
#print axioms FilterCount.balanced_witness
#print axioms FilterCount.share_mono_witness
#print axioms FilterCount.counterexample_any_rounding_witness

/-! ## Target 4: q-value against BH (`Spec/QvalueBh.lean`) -/

example : QvalueBh.bh_sorted_neg_ties_counterexample_stmt := QvalueBh.bh_sorted_neg_ties_counterexample

#print axioms QvalueBh.bh_sorted_neg_ties_counterexample
