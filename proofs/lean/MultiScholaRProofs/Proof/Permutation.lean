module

public import MultiScholaRProofs.Spec.Permutation
public import MultiScholaRProofs.Proof.QvalueBh
import Mathlib.Tactic

/-!
# Proofs for target 5 (KIT-193 section 6.5)

Theorems whose types are the frozen statements of `Spec/Permutation.lean`: the identity clause, the tolerance, the
floors, and the BH floor with its counterexample sentence.
-/

public section

namespace MultiScholaRProofs.Permutation

open Finset MultiScholaRProofs.QvalueBh

/-- (5a') the identity clause is the exact count. -/
theorem identity_clause : identity_clause_stmt := by
  intro n X T f σ ε hε
  constructor
  · rintro (h | h)
    · exact h
    · rw [h]
      linarith
  · exact Or.inl

/-- A filter's count is monotone in its predicate (the instances are the goal's own). -/
theorem card_filter_mono {α : Type} [Fintype α] (P R : α → Prop) {hP : DecidablePred P} {hR : DecidablePred R}
    (h : ∀ x, P x → R x) : (@Finset.filter α P hP univ).card ≤ (@Finset.filter α R hR univ).card :=
  card_le_card fun x hx => mem_filter.2 ⟨mem_univ _, h x (mem_filter.1 hx).2⟩

/-- Two disjoint filters together count at most a filter that contains both. -/
theorem card_add_le_of {α : Type} [Fintype α] (P Q R : α → Prop) {hP : DecidablePred P} {hQ : DecidablePred Q}
    {hR : DecidablePred R} (hd : ∀ x, P x → ¬ Q x) (hs : ∀ x, P x ∨ Q x → R x) :
    (@Finset.filter α P hP univ).card + (@Finset.filter α Q hQ univ).card ≤ (@Finset.filter α R hR univ).card := by
  classical
  rw [← card_union_of_disjoint (disjoint_filter.2 fun x _ h1 h2 => hd x h1 h2)]
  exact card_le_card fun x hx => by
    rcases mem_union.1 hx with h | h
    · exact mem_filter.2 ⟨mem_univ _, hs x (Or.inl (mem_filter.1 h).2)⟩
    · exact mem_filter.2 ⟨mem_univ _, hs x (Or.inr (mem_filter.1 h).2)⟩

/-- The tolerance only raises the p-value. -/
theorem perm_tolerance_mono : perm_tolerance_mono_stmt := by
  intro n B X T f π ε hε
  unfold pval
  apply div_le_div_of_nonneg_right _ (by positivity)
  have := card_filter_mono (fun b => T (f ∘ π b) ≥ T f - 0 ∨ f ∘ π b = f)
    (fun b => T (f ∘ π b) ≥ T f - ε ∨ f ∘ π b = f) (hP := by classical infer_instance)
    (hR := by classical infer_instance) (fun b hb => by
      rcases hb with h | h
      · exact Or.inl (by linarith)
      · exact Or.inr h)
  gcongr

/-- (5c) the floor `(1 + d)/(1 + B)`. -/
theorem perm_floor : perm_floor_stmt := by
  intro n B X T ε f π
  unfold pval unchangedDraws
  apply div_le_div_of_nonneg_right _ (by positivity)
  gcongr
  exact Or.inr

/-- (5c, fixture C) reflected draws are counted too. -/
theorem perm_floor_reflect : perm_floor_reflect_stmt := by
  intro n B T ε f π _ hT hε
  unfold pval unchangedDraws reflectedDraws
  apply div_le_div_of_nonneg_right _ (by positivity)
  rw [add_assoc]
  gcongr
  rw [← Nat.cast_add]
  refine Nat.cast_le.2 (card_add_le_of _ _ _ (fun b h1 h2 => ?_) fun b hb => ?_)
  · exact h2.2 h1
  · rcases hb with h | ⟨h, _⟩
    · exact Or.inr h
    · left
      rw [h, hT f]
      linarith

/-- The share `m p_min / q` is at most `m` when `p_min ≤ 1 ≤ q`. -/
theorem floor_le_m (m : ℕ) (q pmin : ℝ) (hq : 1 ≤ q) (hp1 : pmin ≤ 1) :
    (Nat.ceil (m * pmin / q) : ℝ) ≤ m := by
  have hm : (0 : ℝ) ≤ m := Nat.cast_nonneg m
  have : (m : ℝ) * pmin / q ≤ m := by
    rw [div_le_iff₀ (by linarith)]
    nlinarith
  exact_mod_cast Nat.ceil_le.2 (by exact_mod_cast this)

/-- (5c) the BH floor, "adjusted p at most q". -/
theorem bh_floor : bh_floor_stmt := by
  intro m p q pmin hq hpmin hlo hhi ⟨i0, hi0⟩
  have hunit : UnitP p := fun j => ⟨le_trans hpmin.le (hlo j), hhi j⟩
  rcases (bh_le_iff p i0 q).1 hi0 with h1 | ⟨j, _, hj⟩
  · have hall : (univ.filter fun i => bh p i ≤ q) = univ :=
      filter_true_of_mem fun i _ => le_trans (min_le_left _ _) h1
    rw [hall, card_univ, Fintype.card_fin]
    exact floor_le_m m q pmin h1 (le_trans (hlo i0) (hhi i0))
  · set c := cntLe p (p j) with hc
    have hcpos : (0 : ℝ) < c := by exact_mod_cast cntLe_pos p j
    have hsub : c ≤ (univ.filter fun i => bh p i ≤ q).card :=
      card_le_card fun l hl => mem_filter.2 ⟨mem_univ _,
        (bh_le_iff p l q).2 (Or.inr ⟨j, (mem_filter.1 hl).2, hj⟩)⟩
    have hmc : (m : ℝ) * pmin ≤ q * c := by
      unfold bhTerm at hj
      rw [← hc, div_le_iff₀ hcpos] at hj
      have : (m : ℝ) * pmin ≤ m * p j := mul_le_mul_of_nonneg_left (hlo j) (Nat.cast_nonneg m)
      linarith
    have hceil : Nat.ceil (m * pmin / q) ≤ c := Nat.ceil_le.2 (by rw [div_le_iff₀ hq]; linarith)
    exact_mod_cast le_trans hceil hsub

/-- (5c) the BH floor, "adjusted p below q". -/
theorem bh_floor_strict : bh_floor_strict_stmt := by
  intro m p q pmin hq hpmin hlo hhi ⟨i0, hi0⟩
  rcases (bh_lt_iff p i0 q).1 hi0 with h1 | ⟨j, _, hj⟩
  · have hall : (univ.filter fun i => bh p i < q) = univ :=
      filter_true_of_mem fun i _ => lt_of_le_of_lt (min_le_left _ _) h1
    rw [hall, card_univ, Fintype.card_fin]
    exact floor_le_m m q pmin h1.le (le_trans (hlo i0) (hhi i0))
  · set c := cntLe p (p j) with hc
    have hcpos : (0 : ℝ) < c := by exact_mod_cast cntLe_pos p j
    have hsub : c ≤ (univ.filter fun i => bh p i < q).card :=
      card_le_card fun l hl => mem_filter.2 ⟨mem_univ _,
        (bh_lt_iff p l q).2 (Or.inr ⟨j, (mem_filter.1 hl).2, hj⟩)⟩
    have hmc : (m : ℝ) * pmin ≤ q * c := by
      unfold bhTerm at hj
      rw [← hc, div_lt_iff₀ hcpos] at hj
      have : (m : ℝ) * pmin ≤ m * p j := mul_le_mul_of_nonneg_left (hlo j) (Nat.cast_nonneg m)
      linarith
    have hceil : Nat.ceil (m * pmin / q) ≤ c := Nat.ceil_le.2 (by rw [div_le_iff₀ hq]; linarith)
    exact_mod_cast le_trans hceil hsub

/-- (5c') the floor sentence of KIT-191 r1.2 was false. -/
theorem bh_floor_sentence_counterexample : bh_floor_sentence_counterexample_stmt := by
  show (∀ i, (1 : ℚ) / 1001 ≤ (fun i : Fin 100 => if (i : ℕ) < 4 then (2 : ℚ) / 1001 else 1) i) ∧
    (∃ i, bh (fun i : Fin 100 => if (i : ℕ) < 4 then (2 : ℚ) / 1001 else 1) i < 1 / 20) ∧
    Nat.ceil ((100 : ℚ) * (1 / 1001) / (1 / 20)) = 2 ∧
    (univ.filter fun i : Fin 100 => (if (i : ℕ) < 4 then (2 : ℚ) / 1001 else 1) ≤ (1 / 20 : ℚ) * 2 / 100).card = 0
  set p : Fin 100 → ℚ := fun i => if (i : ℕ) < 4 then 2 / 1001 else 1 with hp
  refine ⟨fun i => ?_, ⟨0, ?_⟩, ?_, ?_⟩
  · simp only [hp]
    split_ifs <;> norm_num
  · have hcnt : cntLe p (p 0) = 4 := by
      have hfil : (univ.filter fun l : Fin 100 => p l ≤ p 0) = univ.filter fun l : Fin 100 => (l : ℕ) < 4 := by
        apply filter_congr
        intro l _
        simp only [hp]
        split_ifs <;> norm_num <;> omega
      unfold cntLe
      rw [hfil]
      decide
    rw [bh_lt_iff]
    right
    refine ⟨0, le_rfl, ?_⟩
    unfold bhTerm
    rw [hcnt]
    simp only [hp]
    norm_num
  · rw [Nat.ceil_eq_iff (by norm_num)]
    norm_num
  · rw [card_eq_zero, filter_eq_empty_iff]
    intro i _
    split_ifs <;> norm_num

/-! Witnesses -/

theorem perm_floor_reflect_witness : perm_floor_reflect_witness_stmt := by
  refine ⟨fun h => |h 0 - h 1|, fun h => ?_, fun _ => 0, fun i => if (i : ℕ) = 0 then 1 else 0, ?_⟩
  · simp only [reflect]
    have key : ∀ c : ℝ, |c - h 0 - (c - h 1)| = |h 0 - h 1| := fun c => by
      rw [show c - h 0 - (c - h 1) = -(h 0 - h 1) by ring, abs_neg]
    exact key _
  · norm_num

theorem bh_floor_witness : bh_floor_witness_stmt := by
  refine ⟨2, fun _ => 1 / 2, 9 / 10, 1 / 2, by norm_num, by norm_num, fun _ => le_rfl, fun _ => by norm_num,
    ⟨0, ?_⟩, ?_⟩
  · rw [bh_const]
    norm_num
  · exact Nat.lt_ceil.2 (by norm_num)

end MultiScholaRProofs.Permutation
