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

example : QvalueBh.bh_sorted_eq_stmt := QvalueBh.bh_sorted_eq
example : QvalueBh.bh_sorted_eq_witness_stmt := QvalueBh.bh_sorted_eq_witness

#print axioms QvalueBh.bh_sorted_eq
#print axioms QvalueBh.bh_sorted_eq_witness

example : QvalueBh.qflag_contains_bh_stmt := QvalueBh.qflag_contains_bh
example : QvalueBh.flag_none_stmt := QvalueBh.flag_none
example : QvalueBh.qflag_extra_stmt := QvalueBh.qflag_extra
example : QvalueBh.bh_ge_p_stmt := QvalueBh.bh_ge_p
example : QvalueBh.bh_mono_stmt := QvalueBh.bh_mono
example : QvalueBh.bh_ties_stmt := QvalueBh.bh_ties
example : QvalueBh.bh_stepup_stmt := QvalueBh.bh_stepup
example : QvalueBh.grid_mass_stmt := QvalueBh.grid_mass
example : QvalueBh.grid_superuniform_stmt := QvalueBh.grid_superuniform
example : QvalueBh.fdr_n20_stmt := QvalueBh.fdr_n20
example : QvalueBh.fdr_n20_strict_stmt := QvalueBh.fdr_n20_strict
example : QvalueBh.fdr_n25_stmt := QvalueBh.fdr_n25
example : QvalueBh.qflag_contains_bh_witness_stmt := QvalueBh.qflag_contains_bh_witness
example : QvalueBh.qflag_extra_witness_stmt := QvalueBh.qflag_extra_witness
example : QvalueBh.bh_unit_witness_stmt := QvalueBh.bh_unit_witness
example : QvalueBh.grid_witness_stmt := QvalueBh.grid_witness

#print axioms QvalueBh.qflag_contains_bh
#print axioms QvalueBh.flag_none
#print axioms QvalueBh.qflag_extra
#print axioms QvalueBh.bh_ge_p
#print axioms QvalueBh.bh_mono
#print axioms QvalueBh.bh_ties
#print axioms QvalueBh.bh_stepup
#print axioms QvalueBh.grid_mass
#print axioms QvalueBh.grid_superuniform
#print axioms QvalueBh.fdr_n20
#print axioms QvalueBh.fdr_n20_strict
#print axioms QvalueBh.fdr_n25
#print axioms QvalueBh.qflag_contains_bh_witness
#print axioms QvalueBh.qflag_extra_witness
#print axioms QvalueBh.bh_unit_witness
#print axioms QvalueBh.grid_witness

/-! ## Target 5: the permutation p-value (`Spec/Permutation.lean`) -/

example : Permutation.identity_clause_stmt := Permutation.identity_clause
example : Permutation.perm_tolerance_mono_stmt := Permutation.perm_tolerance_mono
example : Permutation.perm_floor_stmt := Permutation.perm_floor
example : Permutation.perm_floor_reflect_stmt := Permutation.perm_floor_reflect
example : Permutation.bh_floor_stmt := Permutation.bh_floor
example : Permutation.bh_floor_strict_stmt := Permutation.bh_floor_strict
example : Permutation.bh_floor_sentence_counterexample_stmt := Permutation.bh_floor_sentence_counterexample
example : Permutation.perm_floor_reflect_witness_stmt := Permutation.perm_floor_reflect_witness
example : Permutation.bh_floor_witness_stmt := Permutation.bh_floor_witness
example : Permutation.perm_finite_core_stmt := Permutation.perm_finite_core
example : Permutation.perm_witness_stmt := Permutation.perm_witness
example : Permutation.perm_measure_stmt := Permutation.perm_measure
example : Permutation.perm_measure_witness_stmt := Permutation.perm_measure_witness

#print axioms Permutation.identity_clause
#print axioms Permutation.perm_tolerance_mono
#print axioms Permutation.perm_floor
#print axioms Permutation.perm_floor_reflect
#print axioms Permutation.bh_floor
#print axioms Permutation.bh_floor_strict
#print axioms Permutation.bh_floor_sentence_counterexample
#print axioms Permutation.perm_floor_reflect_witness
#print axioms Permutation.bh_floor_witness
#print axioms Permutation.perm_finite_core
#print axioms Permutation.perm_witness
#print axioms Permutation.perm_measure
#print axioms Permutation.perm_measure_witness

/-! ## Target 5's oracle against the Spec (`Proof/PermutationOracle.lean`; ticket 5.5)

`Oracle/Permutation.lean` counts over lists of rationals; these theorems state, with their types written out here,
that its counts are the Spec's `unchangedDraws` and `reflectedDraws`, and that its 24 orderings are exactly the
permutations of `Fin 4`, each once. -/

example : ∀ {n B : ℕ} (f : Fin n → ℚ) (π : Fin B → Equiv.Perm (Fin n)),
    Oracle.Permutation.unchanged (List.ofFn f) (Permutation.encDraws π) = Permutation.unchangedDraws f π :=
  Permutation.unchanged_eq_unchangedDraws
example : ∀ {n B : ℕ} (f : Fin n → ℚ) (π : Fin B → Equiv.Perm (Fin n)),
    Oracle.Permutation.reflected (List.ofFn f) (Permutation.encDraws π) =
      Permutation.reflectedDraws (fun i => (f i : ℝ)) π :=
  Permutation.reflected_eq_reflectedDraws
example : ∀ l ∈ Oracle.Permutation.orderings, ∃ σ : Equiv.Perm (Fin 4), Permutation.encPerm σ = l :=
  Permutation.orderings_sound
example : ∀ σ : Equiv.Perm (Fin 4), Permutation.encPerm σ ∈ Oracle.Permutation.orderings :=
  Permutation.orderings_complete
example : Oracle.Permutation.orderings.Nodup ∧ Oracle.Permutation.orderings.length = 24 :=
  Permutation.orderings_nodup

#print axioms Permutation.unchanged_eq_unchangedDraws
#print axioms Permutation.reflected_eq_reflectedDraws
#print axioms Permutation.row_counts
#print axioms Permutation.orderings_sound
#print axioms Permutation.orderings_complete
#print axioms Permutation.orderings_nodup

/-! ## Group R priority 1 (`Spec/RuvAlgebra.lean`, `Spec/RuvDownstream.lean`) -/

example : Ruv.A1_stmt := Ruv.A1
example : Ruv.A2_stmt := Ruv.A2
example : Ruv.A3_stmt := Ruv.A3
example : Ruv.C1_stmt := Ruv.C1
example : Ruv.A2_witness_stmt := Ruv.A2_witness
example : Ruv.F1_stmt := Ruv.F1
example : Ruv.F_witness_stmt := Ruv.F_witness
example : Ruv.G0_stmt := Ruv.G0
example : Ruv.G0_witness_stmt := Ruv.G0_witness
example : Ruv.G2_stmt := Ruv.G2
example : Ruv.G2_sign_stmt := Ruv.G2_sign
example : Ruv.H1_stmt := Ruv.H1
example : Ruv.H1_weighted_stmt := Ruv.H1_weighted
example : Ruv.H1_witness_stmt := Ruv.H1_witness

#print axioms Ruv.A1
#print axioms Ruv.A2
#print axioms Ruv.A3
#print axioms Ruv.C1
#print axioms Ruv.A2_witness
#print axioms Ruv.F1
#print axioms Ruv.F_witness
#print axioms Ruv.G0
#print axioms Ruv.G0_witness
#print axioms Ruv.G2
#print axioms Ruv.G2_sign
#print axioms Ruv.H1
#print axioms Ruv.H1_weighted
#print axioms Ruv.H1_witness

/-! ## Sharper witnesses (proofs audit P03; `Spec/SharperWitnesses.lean`) -/

example : Permutation.perm_measure_moving_witness_stmt := Permutation.perm_measure_moving_witness
example : Ruv.G0_interior_witness_stmt := Ruv.G0_interior_witness

#print axioms Permutation.perm_measure_moving_witness
#print axioms Ruv.G0_interior_witness

/-! ## Group R priority 2 -/

example : Ruv.F2_stmt := Ruv.F2

#print axioms Ruv.F2

example : Ruv.B1_stmt := Ruv.B1
example : Ruv.B2_stmt := Ruv.B2
example : Ruv.D1_stmt := Ruv.D1
example : Ruv.D2_stmt := Ruv.D2
example : Ruv.E1_stmt := Ruv.E1
example : Ruv.E2_stmt := Ruv.E2
example : Ruv.K1_stmt := Ruv.K1
example : Ruv.K1_gls_stmt := Ruv.K1_gls
example : Ruv.B1_witness_stmt := Ruv.B1_witness
example : Ruv.D_witness_stmt := Ruv.D_witness
example : Ruv.K1_witness_stmt := Ruv.K1_witness
example : Ruv.K1_gls_witness_stmt := Ruv.K1_gls_witness

#print axioms Ruv.B1
#print axioms Ruv.B2
#print axioms Ruv.D1
#print axioms Ruv.D2
#print axioms Ruv.E1
#print axioms Ruv.E2
#print axioms Ruv.K1
#print axioms Ruv.K1_gls
#print axioms Ruv.B1_witness
#print axioms Ruv.D_witness
#print axioms Ruv.K1_witness
#print axioms Ruv.K1_gls_witness

example : Ruv.G1_stmt := Ruv.G1
example : Ruv.G1_witness_stmt := Ruv.G1_witness
example : Ruv.J1_stmt := Ruv.J1
example : Ruv.J1_witness_stmt := Ruv.J1_witness
example : Ruv.J2_stmt := Ruv.J2
example : Ruv.J2_contribution_stmt := Ruv.J2_contribution
example : Ruv.J2_witness_stmt := Ruv.J2_witness
example : Ruv.L1_stmt := Ruv.L1
example : Ruv.L1_invariance_stmt := Ruv.L1_invariance
example : Ruv.L1_invariance_witness_stmt := Ruv.L1_invariance_witness

#print axioms Ruv.G1
#print axioms Ruv.G1_witness
#print axioms Ruv.J1
#print axioms Ruv.J1_witness
#print axioms Ruv.J2
#print axioms Ruv.J2_contribution
#print axioms Ruv.J2_witness
#print axioms Ruv.L1
#print axioms Ruv.L1_invariance
#print axioms Ruv.L1_invariance_witness

/-! ## Group R priority 3 -/

example : Ruv.F3_span_stmt := Ruv.F3_span
example : Ruv.F3_fwl_stmt := Ruv.F3_fwl
example : Ruv.F3_fwl_witness_stmt := Ruv.F3_fwl_witness
example : Ruv.L2_stmt := Ruv.L2
example : Ruv.L2_witness_stmt := Ruv.L2_witness

#print axioms Ruv.F3_span
#print axioms Ruv.F3_fwl
#print axioms Ruv.F3_fwl_witness
#print axioms Ruv.L2
#print axioms Ruv.L2_witness

/-! ## Group R priority 4 -/

example : Ruv.M1_stmt := Ruv.M1
example : Ruv.M1_witness_stmt := Ruv.M1_witness

#print axioms Ruv.M1
#print axioms Ruv.M1_witness

/-! ## Target 4d-ii (`Proof/KeptFamily.lean`) -/

example : QvalueBh.kept_family_unit_a_counterexample_stmt := QvalueBh.kept_family_unit_a_counterexample
example : QvalueBh.kept_family_formula_stmt := QvalueBh.kept_family_formula
example : QvalueBh.kept_family_instance_stmt := QvalueBh.kept_family_instance
example : QvalueBh.kept_family_witness_stmt := QvalueBh.kept_family_witness

#print axioms QvalueBh.kept_family_unit_a_counterexample
#print axioms QvalueBh.kept_family_formula
#print axioms QvalueBh.kept_family_instance
#print axioms QvalueBh.kept_family_witness

/-! ## Group R, R-I1 (`Proof/RemovedFactors.lean`) -/

example : Ruv.I1_UtWc_stmt := Ruv.I1_UtWc
example : Ruv.I1_centred_diff_stmt := Ruv.I1_centred_diff
example : Ruv.I1_centred_diff_rows_stmt := Ruv.I1_centred_diff_rows
example : Ruv.I1_span_stmt := Ruv.I1_span
example : Ruv.I1_rank_stmt := Ruv.I1_rank
example : Ruv.I1_ctl_rank_stmt := Ruv.I1_ctl_rank
example : Ruv.I1_average_stmt := Ruv.I1_average
example : Ruv.I1_witness_stmt := Ruv.I1_witness
example : Ruv.I1_average_witness_stmt := Ruv.I1_average_witness

#print axioms Ruv.I1_UtWc
#print axioms Ruv.I1_centred_diff
#print axioms Ruv.I1_centred_diff_rows
#print axioms Ruv.I1_span
#print axioms Ruv.I1_rank
#print axioms Ruv.I1_ctl_rank
#print axioms Ruv.I1_average
#print axioms Ruv.I1_witness
#print axioms Ruv.I1_average_witness

/-! ## Group R, R-J3 (`Proof/RuvIIICPattern.lean`) -/

example : Ruv.J3_pattern_stmt := Ruv.J3_pattern
example : Ruv.J3_counts_stmt := Ruv.J3_counts
example : Ruv.J3_witness_stmt := Ruv.J3_witness

#print axioms Ruv.J3_pattern
#print axioms Ruv.J3_counts
#print axioms Ruv.J3_witness

/-! ## Group R, R-N1 (`Proof/RuvControlRules.lean`) -/

example : RuvControls.N1_stable_stmt := RuvControls.N1_stable
example : RuvControls.N1_ties_stmt := RuvControls.N1_ties
example : RuvControls.N1_first_pass_stmt := RuvControls.N1_first_pass
example : RuvControls.N1_trx_stmt := RuvControls.N1_trx
example : RuvControls.N1_rules_differ_stmt := RuvControls.N1_rules_differ
example : RuvControls.N1_stable_witness_stmt := RuvControls.N1_stable_witness

#print axioms RuvControls.N1_stable
#print axioms RuvControls.N1_ties
#print axioms RuvControls.N1_first_pass
#print axioms RuvControls.N1_trx
#print axioms RuvControls.N1_rules_differ
#print axioms RuvControls.N1_stable_witness
