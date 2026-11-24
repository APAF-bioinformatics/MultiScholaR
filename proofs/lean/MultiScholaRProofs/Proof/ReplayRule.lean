module

public import MultiScholaRProofs.Spec.ReplayRule
import Mathlib.Tactic

/-!
# Proofs for target 2 (KIT-193 section 6.2)

Theorems whose types are the frozen statements of `Spec/ReplayRule.lean`. The round trip (2c) is faithfulness (2b)
under the option state the replay sets, whose value for a fixed marker is that marker's recorded rule.
-/

public section

namespace MultiScholaRProofs.ReplayRule

section lemmas

-- A style linter only (an instance argument a lemma does not use); it has no bearing on what is proved.
set_option linter.unusedSectionVars false

variable {N O V : Type} [DecidableEq N] [DecidableEq O] [DecidableEq V]

theorem recordsCurrent_iff (μ : Marker N O V) (e : Entry N V) :
    recordsCurrent μ e = true ↔ e.name = μ.name ∧ e.source = Source.fixed ∧ e.value = μ.current := by
  simp [recordsCurrent, and_assoc]

theorem rule_eq_current_or (step : Step N V) (μ : Marker N O V) :
    rule step μ = μ.current ∨ rule step μ = μ.earlier := by
  unfold rule
  split_ifs
  · exact Or.inl rfl
  · exact Or.inr rfl

theorem rule_of_exists (step : Step N V) (μ : Marker N O V) (h : ∃ e ∈ step, recordsCurrent μ e = true) :
    rule step μ = μ.current := by
  unfold rule
  simp only [h, ↓reduceIte]

theorem rule_of_not_exists (step : Step N V) (μ : Marker N O V) (h : ¬ ∃ e ∈ step, recordsCurrent μ e = true) :
    rule step μ = μ.earlier := by
  unfold rule
  simp only [h, ↓reduceIte]

/-- A built step has no current-rule entry for `μ` when the provenance has no entry named `μ` and no other entry has
source `fixed`. -/
theorem rule_build_without (tbl : List (Marker N O V)) (μ : Marker N O V) (others : Step N V)
    (prov decl : Fixed N V) (hprov : μ.name ∉ prov.map Prod.fst)
    (hothers : ∀ e ∈ others, e.source ≠ Source.fixed) :
    rule (build tbl others prov decl) μ = μ.earlier := by
  apply rule_of_not_exists
  rintro ⟨e, he, hrc⟩
  rw [recordsCurrent_iff] at hrc
  obtain ⟨hname, hsrc, -⟩ := hrc
  simp only [build, List.mem_append] at he
  rcases he with he | he
  · exact hothers e he hsrc
  · split_ifs at he
    · simp only [List.mem_map] at he
      obtain ⟨x, hx, rfl⟩ := he
      exact hprov (List.mem_map.2 ⟨x, hx, hname⟩)
    · simp only [List.mem_map] at he
      obtain ⟨x, _, rfl⟩ := he
      exact absurd hsrc (by simp)

theorem eq_of_name (tbl : List (Marker N O V)) (htbl : (tbl.map Marker.name).Nodup) {μ ν : Marker N O V}
    (hμ : μ ∈ tbl) (hν : ν ∈ tbl) (h : μ.name = ν.name) : μ = ν :=
  List.inj_on_of_nodup_map htbl hμ hν h

theorem eq_of_option (tbl : List (Marker N O V)) (htbl : (tbl.map Marker.option).Nodup) {μ ν : Marker N O V}
    (hμ : μ ∈ tbl) (hν : ν ∈ tbl) (h : μ.option = ν.option) : μ = ν :=
  List.inj_on_of_nodup_map htbl hμ hν h

theorem find_name (tbl : List (Marker N O V)) (htbl : TableOK tbl) (μ : Marker N O V) (hμ : μ ∈ tbl) :
    tbl.find? (fun ν => decide (ν.name = μ.name)) = some μ := by
  obtain ⟨ν, hν⟩ := Option.isSome_iff_exists.1
    (List.find?_isSome.2 ⟨μ, hμ, by simp⟩ : (tbl.find? fun ν => decide (ν.name = μ.name)).isSome)
  have hνmem : ν ∈ tbl := List.mem_of_find?_eq_some hν
  have hνname : ν.name = μ.name := by simpa using List.find?_some hν
  rw [hν, eq_of_name tbl htbl.1 hνmem hμ hνname]

theorem mem_markersFixed (tbl : List (Marker N O V)) (htbl : TableOK tbl) (decl : Fixed N V)
    (ν : Marker N O V) : ν ∈ markersFixed tbl decl ↔ ν ∈ tbl ∧ ν.name ∈ decl.map Prod.fst := by
  simp only [markersFixed, List.mem_filterMap, List.mem_dedup]
  constructor
  · rintro ⟨nm, hnm, hfind⟩
    have hmem : ν ∈ tbl := List.mem_of_find?_eq_some hfind
    have hname : ν.name = nm := by simpa using List.find?_some hfind
    exact ⟨hmem, hname ▸ hnm⟩
  · rintro ⟨hν, hname⟩
    exact ⟨ν.name, hname, find_name tbl htbl ν hν⟩

theorem lookup_options (tbl : List (Marker N O V)) (htbl : TableOK tbl) (f : Marker N O V → V) :
    ∀ (L : List (Marker N O V)), (∀ ν ∈ L, ν ∈ tbl) → ∀ μ ∈ L,
      (L.map fun ν => (ν.option, f ν)).lookup μ.option = some (f μ)
  | [], _, μ, hμ => absurd hμ (by simp)
  | ν :: L, hL, μ, hμ => by
    simp only [List.map_cons, List.lookup]
    by_cases hopt : μ.option = ν.option
    · have : μ = ν := eq_of_option tbl htbl.2 (hL μ hμ) (hL ν (by simp)) hopt
      subst this
      simp
    · have hne : μ ≠ ν := fun h => hopt (h ▸ rfl)
      have hμL : μ ∈ L := by
        rcases List.mem_cons.1 hμ with h | h
        · exact absurd h hne
        · exact h
      have hbeq : (μ.option == ν.option) = false := by simpa using hopt
      rw [hbeq]
      exact lookup_options tbl htbl f L (fun x hx => hL x (by simp [hx])) μ hμL

theorem lookup_options_none (f : Marker N O V → V) (o : O) :
    ∀ (L : List (Marker N O V)), (∀ ν ∈ L, ν.option ≠ o) → (L.map fun ν => (ν.option, f ν)).lookup o = none
  | [], _ => rfl
  | ν :: L, hL => by
    simp only [List.map_cons, List.lookup]
    have hne : (o == ν.option) = false := by simpa using (hL ν (by simp)).symm
    rw [hne]
    exact lookup_options_none f o L (fun x hx => hL x (by simp [hx]))

/-- Faithfulness for every option state, as a lemma. -/
theorem faithful_lemma (tbl : List (Marker N O V)) (μ : Marker N O V) (others : Step N V) (decl : Fixed N V)
    (o : Opts O V) (htbl : TableOK tbl) (hμ : μ ∈ tbl) (hdecl : (μ.name, μ.current) ∈ decl)
    (hothers : ∀ e ∈ others, e.source ≠ Source.fixed) :
    rule (build tbl others (compute tbl decl o) decl) μ = helperRule o μ := by
  unfold helperRule
  split_ifs with hopt
  · apply rule_build_without tbl μ others _ decl _ hothers
    simp only [compute, inForce, List.mem_map, List.mem_filter, decide_eq_true_eq, not_exists, not_and]
    rintro ⟨nm, v⟩ ⟨_, hkeep⟩ hnm
    simp only at hnm
    exact hkeep μ hμ hnm.symm hopt
  · apply rule_of_exists
    have hkept : (μ.name, μ.current) ∈ compute tbl decl o := by
      simp only [compute, inForce, List.mem_filter, decide_eq_true_eq, not_exists, not_and]
      refine ⟨hdecl, fun ν hν hname hν' => hopt ?_⟩
      rw [← eq_of_name tbl htbl.1 hν hμ hname] at *
      exact hν'
    have hne : compute tbl decl o ≠ [] := List.ne_nil_of_mem hkept
    refine ⟨⟨μ.name, Source.fixed, μ.current⟩, ?_, by simp [recordsCurrent]⟩
    simp only [build, List.mem_append, hne, ne_eq, not_false_eq_true, ↓reduceIte, List.mem_map]
    exact Or.inr ⟨_, hkept, rfl⟩

end lemmas

/-- `ruleLoop` (first match) equals `rule`. -/
theorem rule_first_match : rule_first_match_stmt := by
  intro N O V _ _ _ step μ
  unfold ruleLoop rule
  cases h : step.find? (recordsCurrent μ) with
  | none =>
    rw [List.find?_eq_none] at h
    have hn : ¬ ∃ e ∈ step, recordsCurrent μ e = true := fun ⟨e, he, hrc⟩ => h e he hrc
    simp only [hn, ↓reduceIte]
  | some e =>
    have hy : ∃ e ∈ step, recordsCurrent μ e = true := ⟨e, List.mem_of_find?_eq_some h, List.find?_some h⟩
    simp only [hy, ↓reduceIte]

/-- (2a) totality. -/
theorem totality : totality_stmt := by
  intro N O V _ _ _ tbl step decl htbl
  refine ⟨fun μ => rule_eq_current_or step μ, fun o => ?_⟩
  simp only [replayOptions, List.map_map, List.mem_map, Function.comp_def, mem_markersFixed tbl htbl]
  constructor
  · rintro ⟨ν, ⟨hν, hname⟩, rfl⟩
    exact ⟨ν, hν, rfl, hname⟩
  · rintro ⟨ν, hν, rfl, hname⟩
    exact ⟨ν, ⟨hν, hname⟩, rfl⟩

/-- (2b) faithfulness. -/
theorem faithful : faithful_stmt := by
  intro N O V _ _ _ tbl μ others decl o htbl hμ hdecl hothers
  exact faithful_lemma tbl μ others decl o htbl hμ hdecl (fun e he => (hothers e he).1)

/-- (2b) legacy steps, (2e) new-design steps. -/
theorem without_marker : without_marker_stmt := by
  intro N O V _ _ _ tbl μ others prov decl hprov hothers
  exact rule_build_without tbl μ others prov decl hprov (fun e he => (hothers e he).1)

/-- (2e) the managed new-design writer. -/
theorem managed_new_design : managed_new_design_stmt := by
  intro N O V _ _ _ tbl μ dropped others decl o hdrop hothers
  apply rule_build_without tbl μ others _ decl _ (fun e he => (hothers e he).1)
  simp only [computeManagedNewDesign, inForce, List.mem_map, List.mem_filter, decide_eq_true_eq, not_exists,
    not_and]
  rintro ⟨nm, v⟩ ⟨⟨_, hnd⟩, _⟩ hnm
  simp only at hnm hnd
  exact hnd (hnm ▸ hdrop)

/-- (2c) round trip. -/
theorem round_trip : round_trip_stmt := by
  intro N O V _ _ _ tbl μ step others decl ambient htbl hμ hdecl hothers
  rw [faithful_lemma tbl μ others decl _ htbl hμ hdecl (fun e he => (hothers e he).1)]
  have hmem : μ ∈ markersFixed tbl decl :=
    (mem_markersFixed tbl htbl decl μ).2 ⟨hμ, List.mem_map.2 ⟨_, hdecl, rfl⟩⟩
  have hopt : withOpts ambient (replayOptions tbl step decl) μ.option = some (rule step μ) := by
    unfold withOpts replayOptions
    rw [lookup_options tbl htbl (fun ν => rule step ν) _
      (fun ν hν => ((mem_markersFixed tbl htbl decl ν).1 hν).1) μ hmem]
  unfold helperRule
  rw [hopt]
  rcases rule_eq_current_or step μ with h | h
  · rw [h]
    split_ifs with h'
    · exact (Option.some.inj h').symm
    · rfl
  · rw [h]
    simp

/-- (2d) locality. -/
theorem locality : locality_stmt := by
  intro N O V _ _ _ step μ
  unfold rule
  congr 1
  apply propext
  constructor
  · rintro ⟨e, he, hrc⟩
    refine ⟨e, List.mem_filter.2 ⟨he, ?_⟩, hrc⟩
    simpa using ((recordsCurrent_iff μ e).1 hrc).1
  · rintro ⟨e, he, hrc⟩
    exact ⟨e, (List.mem_filter.1 he).1, hrc⟩

/-- (2d) an unfixed marker sees the helper's default. -/
theorem unfixed_default : unfixed_default_stmt := by
  intro N O V _ _ _ tbl ν step decl ambient htbl hν hname hamb
  have hnone : withOpts ambient (replayOptions tbl step decl) ν.option = none := by
    unfold withOpts replayOptions
    rw [lookup_options_none (fun x => rule step x) ν.option _ ?_]
    · exact hamb
    · intro x hx hopt
      obtain ⟨hxt, hxn⟩ := (mem_markersFixed tbl htbl decl x).1 hx
      have := eq_of_option tbl htbl.2 hxt hν hopt
      exact hname (this ▸ hxn)
  unfold helperRule
  rw [hnone]
  simp

/-! Witnesses: the marker table after KIT-190. -/

/-- KIT-187's row. -/
def markerMissing : Marker String String String :=
  ⟨"missing_value_comparison", "MultiScholaR.missing_filter_rule", "counts", "rounded_percentage"⟩

/-- KIT-190's row. -/
def markerRuv : Marker String String String :=
  ⟨"ruv_control_grouping", "MultiScholaR.ruv_control_grouping", "kept", "dropped"⟩

theorem twoRow : TwoRowTable [markerMissing, markerRuv] := by
  refine ⟨⟨?_, ?_⟩, rfl, ?_⟩
  · simp [markerMissing, markerRuv]
  · simp [markerMissing, markerRuv]
  · intro μ hμ
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hμ
    rcases hμ with rfl | rfl <;> simp [markerMissing, markerRuv]

theorem faithful_witness : faithful_witness_stmt :=
  ⟨[markerMissing, markerRuv], markerMissing, [⟨"x", Source.args, "1"⟩],
    [(markerMissing.name, markerMissing.current)], twoRow, by simp, by simp, by simp, by simp⟩

theorem without_marker_witness : without_marker_witness_stmt :=
  ⟨[markerMissing, markerRuv], markerMissing, [], [("other", "v")], twoRow, by simp, by simp,
    by simp [markerMissing], by simp⟩

theorem managed_new_design_witness : managed_new_design_witness_stmt :=
  ⟨[markerMissing, markerRuv], markerRuv, ["ruv_control_grouping"], twoRow, by simp, by simp [markerRuv]⟩

theorem unfixed_default_witness : unfixed_default_witness_stmt :=
  ⟨[markerMissing, markerRuv], markerRuv, [(markerMissing.name, markerMissing.current)], twoRow, by simp,
    by simp, by simp [markerMissing, markerRuv]⟩

end MultiScholaRProofs.ReplayRule
