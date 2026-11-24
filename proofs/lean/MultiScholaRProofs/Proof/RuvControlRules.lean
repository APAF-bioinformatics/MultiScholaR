module

public import MultiScholaRProofs.Spec.RuvControlRules
import Mathlib.Tactic

/-!
# Proofs for group R, R-N1: the control-selection specifications (KIT-193 section 7.2)

Theorems whose types are the frozen statements of `Spec/RuvControlRules.lean` (frozen at 2ab4c189): `N1_stable`,
`N1_ties`, `N1_first_pass`, `N1_trx`, `N1_rules_differ`, and the witness `N1_stable_witness`.

How: (a) and (a') from the sortedness and stability of `List.mergeSort` (Lean core: `pairwise_mergeSort`,
`mergeSort_of_pairwise`, `pair_sublist_mergeSort`). (b) and (c) from one counting lemma (`card_rank_gt`): under a
strict total order on a finite set `S` of `N` elements, the rank `1 + #{j before i}` is a bijection from `S` onto
`1..N`, so exactly `N - m` elements have rank above `m`; and an element before another has the smaller rank.
-/

public section

namespace MultiScholaRProofs.RuvControls

open Finset

/-! ## Ranks under a strict total order -/

section Rank

variable {α : Type} [Fintype α] (lt : α → α → Prop) [DecidableRel lt]

theorem rank_lt_rank (htrans : ∀ a b c, lt a b → lt b c → lt a c) (hirr : ∀ a, ¬ lt a a) {a b : α}
    (h : lt a b) : (univ.filter fun j => lt j a).card < (univ.filter fun j => lt j b).card := by
  apply Finset.card_lt_card
  rw [Finset.ssubset_iff_of_subset]
  · exact ⟨a, by simp [h], by simp [hirr]⟩
  · intro j
    simp only [Finset.mem_filter, Finset.mem_univ, true_and]
    exact fun hj => htrans _ _ _ hj h

/-- Exactly `N - m` elements of `S` have rank above `m`, for `m ≤ N = #S`. -/
theorem card_rank_gt [DecidableEq α] (S : Finset α) (hS : ∀ a b, lt a b → a ∈ S)
    (htrans : ∀ a b c, lt a b → lt b c → lt a c) (hirr : ∀ a, ¬ lt a a)
    (htot : ∀ a ∈ S, ∀ b ∈ S, a ≠ b → lt a b ∨ lt b a) (m : ℕ) (hm : m ≤ S.card) :
    (S.filter fun i => m < 1 + (univ.filter fun j => lt j i).card).card = S.card - m := by
  set f : α → ℕ := fun i => 1 + (univ.filter fun j => lt j i).card with hf
  have hinj : Set.InjOn f S := by
    intro a ha b hb hab
    by_contra hne
    rcases htot a ha b hb hne with h | h
    · have := rank_lt_rank lt htrans hirr h
      simp only [hf] at hab
      omega
    · have := rank_lt_rank lt htrans hirr h
      simp only [hf] at hab
      omega
  have hmaps : ∀ i ∈ S, f i ∈ Icc 1 S.card := by
    intro i hi
    have hsub : (univ.filter fun j => lt j i) ⊆ S.erase i := by
      intro j hj
      simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hj
      exact Finset.mem_erase.mpr ⟨fun h => hirr i (h ▸ hj), hS _ _ hj⟩
    have h1 := Finset.card_le_card hsub
    rw [Finset.card_erase_of_mem hi] at h1
    have h2 : 0 < S.card := Finset.card_pos.mpr ⟨i, hi⟩
    simp only [Finset.mem_Icc, hf]
    omega
  have himage : S.image f = Icc 1 S.card := by
    apply Finset.eq_of_subset_of_card_le
    · intro x hx
      obtain ⟨i, hi, rfl⟩ := Finset.mem_image.mp hx
      exact hmaps i hi
    · rw [Nat.card_Icc, Finset.card_image_of_injOn hinj]
      omega
  have hsub : (S.filter fun i => m < f i).image f = (Icc 1 S.card).filter fun x => m < x := by
    rw [← himage, Finset.filter_image]
  have hcard : ((S.filter fun i => m < f i).image f).card = (S.filter fun i => m < f i).card :=
    Finset.card_image_of_injOn (hinj.mono (by intro x hx; exact (Finset.mem_filter.mp hx).1))
  have hIcc : (Icc 1 S.card).filter (fun x => m < x) = Icc (m + 1) S.card := by
    ext x
    simp only [Finset.mem_filter, Finset.mem_Icc]
    omega
  show (S.filter fun i => m < f i).card = S.card - m
  rw [← hcard, hsub, hIcc, Nat.card_Icc]
  omega

end Rank

/-! ## The order "smaller p, then earlier row" -/

theorem before_iff {F : ℕ} (pv : Fin F → Option ℚ) (j i : Fin F) :
    before pv j i ↔ ∃ a b, pv j = some a ∧ pv i = some b ∧ (a < b ∨ (a = b ∧ j < i)) := by
  cases hj : pv j <;> cases hi : pv i <;> simp [before, hj, hi]

theorem lex_trans {F : ℕ} {x y z : ℚ} {a b c : Fin F} :
    (x < y ∨ (x = y ∧ a < b)) → (y < z ∨ (y = z ∧ b < c)) → (x < z ∨ (x = z ∧ a < c)) := by
  rintro (h1 | ⟨rfl, h1⟩) (h2 | ⟨rfl, h2⟩)
  · exact Or.inl (h1.trans h2)
  · exact Or.inl h1
  · exact Or.inl h2
  · exact Or.inr ⟨rfl, h1.trans h2⟩

theorem before_trans {F : ℕ} (pv : Fin F → Option ℚ) (a b c : Fin F) :
    before pv a b → before pv b c → before pv a c := by
  rw [before_iff, before_iff, before_iff]
  rintro ⟨x, y, hx, hy, h1⟩ ⟨y', z, hy', hz, h2⟩
  rw [hy] at hy'
  cases hy'
  exact ⟨x, z, hx, hz, lex_trans h1 h2⟩

theorem before_irrefl {F : ℕ} (pv : Fin F → Option ℚ) (a : Fin F) : ¬ before pv a a := by
  rw [before_iff]
  rintro ⟨x, y, hx, hy, h⟩
  rw [hx] at hy
  cases hy
  rcases h with h | ⟨_, h⟩
  · exact lt_irrefl _ h
  · exact lt_irrefl _ h

theorem before_present {F : ℕ} (pv : Fin F → Option ℚ) (j i : Fin F) (h : before pv j i) :
    (pv j).isSome = true ∧ (pv i).isSome = true := by
  rw [before_iff] at h
  obtain ⟨x, y, hx, hy, _⟩ := h
  simp [hx, hy]

theorem before_total {F : ℕ} (pv : Fin F → Option ℚ) (a b : Fin F) (ha : (pv a).isSome = true)
    (hb : (pv b).isSome = true) (hne : a ≠ b) : before pv a b ∨ before pv b a := by
  obtain ⟨x, hx⟩ := Option.isSome_iff_exists.mp ha
  obtain ⟨y, hy⟩ := Option.isSome_iff_exists.mp hb
  rw [before_iff, before_iff]
  rcases lt_trichotomy x y with h | rfl | h
  · exact Or.inl ⟨x, y, hx, hy, Or.inl h⟩
  · rcases lt_or_gt_of_ne hne with h | h
    · exact Or.inl ⟨x, x, hx, hy, Or.inr ⟨rfl, h⟩⟩
    · exact Or.inr ⟨x, x, hy, hx, Or.inr ⟨rfl, h⟩⟩
  · exact Or.inr ⟨y, x, hy, hx, Or.inl h⟩

/-! ## R's `order()`: missing last, then "smaller p, then earlier row" -/

theorem orderBefore_trans {G : ℕ} (pv : Fin G → Option ℚ) (a b c : Fin G) :
    orderBefore pv a b → orderBefore pv b c → orderBefore pv a c := by
  unfold orderBefore
  rintro (⟨ha, hb⟩ | h1 | ⟨ha, hb, h1⟩) (⟨hb', hc⟩ | h2 | ⟨hb', hc, h2⟩)
  · simp [hb] at hb'
  · have := (before_present pv b c h2).1
    simp [hb] at this
  · exact Or.inl ⟨ha, hc⟩
  · exact Or.inl ⟨(before_present pv a b h1).1, hc⟩
  · exact Or.inr (Or.inl (before_trans pv a b c h1 h2))
  · have := (before_present pv a b h1).2
    simp [hb'] at this
  · simp [hb] at hb'
  · have := (before_present pv b c h2).1
    simp [hb] at this
  · exact Or.inr (Or.inr ⟨ha, hc, h1.trans h2⟩)

theorem orderBefore_irrefl {G : ℕ} (pv : Fin G → Option ℚ) (a : Fin G) : ¬ orderBefore pv a a := by
  unfold orderBefore
  rintro (⟨ha, ha'⟩ | h | ⟨_, _, h⟩)
  · simp [ha'] at ha
  · exact before_irrefl pv a h
  · exact lt_irrefl _ h

theorem orderBefore_total {G : ℕ} (pv : Fin G → Option ℚ) (a b : Fin G) (hne : a ≠ b) :
    orderBefore pv a b ∨ orderBefore pv b a := by
  unfold orderBefore
  cases ha : pv a <;> cases hb : pv b
  · rcases lt_or_gt_of_ne hne with h | h
    · exact Or.inl (Or.inr (Or.inr ⟨rfl, rfl, h⟩))
    · exact Or.inr (Or.inr (Or.inr ⟨rfl, rfl, h⟩))
  · exact Or.inr (Or.inl ⟨by simp, rfl⟩)
  · exact Or.inl (Or.inl ⟨by simp, rfl⟩)
  · rcases before_total pv a b (by simp [ha]) (by simp [hb]) hne with h | h
    · exact Or.inl (Or.inr (Or.inl h))
    · exact Or.inr (Or.inr (Or.inl h))

/-! ## R-N1 (b) and (c) -/

theorem quarter_lt_iff (N r : ℕ) : (1 / 4 : ℚ) * (N : ℚ) < (r : ℚ) ↔ N / 4 < r := by
  rw [Nat.div_lt_iff_lt_mul (by norm_num)]
  constructor
  · intro h
    have h' : (N : ℚ) < (r : ℚ) * 4 := by linarith
    exact_mod_cast h'
  · intro h
    have h' : (N : ℚ) < (r : ℚ) * 4 := by exact_mod_cast h
    linarith

/-- R-N1 (b). -/
theorem N1_first_pass : N1_first_pass_stmt := by
  intro F pv
  have hC : firstPassControls pv = (univ.filter fun i => (pv i).isSome = true).filter
      fun i => nPresent pv / 4 < rankFirst pv i := by
    ext i
    simp only [firstPassControls, Finset.mem_filter, Finset.mem_univ, true_and]
    rw [quarter_lt_iff]
  refine ⟨?_, ?_⟩
  · rw [hC]
    exact card_rank_gt (before pv) _ (fun a b h => by simpa using (before_present pv a b h).1)
      (before_trans pv) (before_irrefl pv)
      (fun a ha b hb hne => before_total pv a b (by simpa using ha) (by simpa using hb) hne) _
      (Nat.div_le_self _ _)
  · intro i j hi hiC hjC
    have hi' : ¬ nPresent pv / 4 < rankFirst pv i := by
      intro h
      apply hiC
      rw [hC]
      simp only [Finset.mem_filter, Finset.mem_univ, true_and]
      exact ⟨hi, h⟩
    rw [hC] at hjC
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hjC
    have hne : i ≠ j := by
      rintro rfl
      exact hi' hjC.2
    rcases before_total pv i j hi hjC.1 hne with h | h
    · exact h
    · have := rank_lt_rank (before pv) (before_trans pv) (before_irrefl pv) h
      have ei : rankFirst pv i = 1 + (univ.filter fun k => before pv k i).card := rfl
      have ej : rankFirst pv j = 1 + (univ.filter fun k => before pv k j).card := rfl
      exfalso
      apply hi'
      omega

/-- R-N1 (c). -/
theorem N1_trx : N1_trx_stmt := by
  intro G pv
  have hm : ⌈(1 / 4 : ℚ) * G⌉₊ ≤ G := by
    rw [Nat.ceil_le]
    have : (0 : ℚ) ≤ G := by positivity
    linarith
  refine ⟨?_, ?_⟩
  · have h : (trxFirstPassControls pv).card = (univ : Finset (Fin G)).card - ⌈(1 / 4 : ℚ) * G⌉₊ :=
      card_rank_gt (orderBefore pv) univ (fun a b _ => Finset.mem_univ a) (orderBefore_trans pv)
        (orderBefore_irrefl pv) (fun a _ b _ hne => orderBefore_total pv a b hne) ⌈(1 / 4 : ℚ) * G⌉₊
        (by rw [Finset.card_univ, Fintype.card_fin]; exact hm)
    rw [h, Finset.card_univ, Fintype.card_fin]
  · intro i j hiC hjC
    simp only [trxFirstPassControls, Finset.mem_filter, Finset.mem_univ, true_and] at hiC hjC
    have hne : i ≠ j := by
      rintro rfl
      exact hiC hjC
    rcases orderBefore_total pv i j hne with h | h
    · exact h
    · have := rank_lt_rank (orderBefore pv) (orderBefore_trans pv) (orderBefore_irrefl pv) h
      exfalso
      apply hiC
      unfold orderPos at hjC ⊢
      omega

/-! ## R-N1 (a) and (a') -/

/-- R-N1 (a). -/
theorem N1_stable : N1_stable_stmt := by
  intro p q cutoff num hq
  unfold sweepControls
  rw [List.mergeSort_of_pairwise (List.pairwise_of_forall fun i j => by simp [hq i j])]

/-- In a list in strictly increasing order, two members in that order form a sublist. -/
theorem pair_sublist_of_lt {p : ℕ} {l : List (Fin p)} (h : l.Pairwise (· < ·)) {a b : Fin p} (ha : a ∈ l)
    (hb : b ∈ l) (hab : a < b) : List.Sublist [a, b] l := by
  induction l with
  | nil => simp at ha
  | cons x t ih =>
    rw [List.pairwise_cons] at h
    rcases List.mem_cons.mp ha with rfl | ha'
    · have hb' : b ∈ t := by
        rcases List.mem_cons.mp hb with rfl | hb'
        · exact absurd hab (lt_irrefl _)
        · exact hb'
      exact (List.singleton_sublist.mpr hb').cons_cons _
    · rcases List.mem_cons.mp hb with rfl | hb'
      · exact absurd (h.1 a ha') (lt_asymm hab)
      · exact (ih h.2 ha' hb').cons _

/-- In a list without repeats, two distinct values cannot appear in both orders. -/
theorem not_both_orders {p : ℕ} {l : List (Fin p)} (hl : l.Nodup) {a b : Fin p} (hne : a ≠ b)
    (h1 : List.Sublist [a, b] l) (h2 : List.Sublist [b, a] l) : False := by
  induction l with
  | nil => simp at h1
  | cons x t ih =>
    rw [List.nodup_cons] at hl
    rcases List.sublist_cons_iff.mp h1 with h1' | ⟨r, hr, h1'⟩ <;>
      rcases List.sublist_cons_iff.mp h2 with h2' | ⟨r', hr', h2'⟩
    · exact ih hl.2 h1' h2'
    · simp only [List.cons.injEq] at hr'
      obtain ⟨rfl, rfl⟩ := hr'
      exact hl.1 (h1'.subset (by simp))
    · simp only [List.cons.injEq] at hr
      obtain ⟨rfl, rfl⟩ := hr
      exact hl.1 (h2'.subset (by simp))
    · simp only [List.cons.injEq] at hr hr'
      exact hne (hr.1.trans hr'.1.symm)

/-- R-N1 (a'). -/
theorem N1_ties : N1_ties_stmt := by
  intro p q cutoff
  set L := (List.finRange p).filter fun i => decide (cutoff < q i) with hL
  have htrans : ∀ a b c : Fin p, decide (q b ≤ q a) = true → decide (q c ≤ q b) = true →
      decide (q c ≤ q a) = true := by
    intro a b c h1 h2
    simp only [decide_eq_true_eq] at h1 h2 ⊢
    exact h2.trans h1
  have htot : ∀ a b : Fin p, (decide (q b ≤ q a) || decide (q a ≤ q b)) = true := by
    intro a b
    simp only [Bool.or_eq_true, decide_eq_true_eq]
    exact le_total _ _
  have hsorted := List.pairwise_mergeSort htrans htot L
  have hnodupL : L.Nodup := (List.nodup_finRange p).sublist List.filter_sublist
  have hnodupR : (L.mergeSort fun i j => decide (q j ≤ q i)).Nodup :=
    (List.mergeSort_perm L _).symm.nodup hnodupL
  have hLlt : L.Pairwise (· < ·) := (List.pairwise_lt_finRange p).sublist List.filter_sublist
  apply List.pairwise_of_forall_sublist
  intro a b hab
  have hle : q b ≤ q a := by
    have := List.pairwise_iff_forall_sublist.mp hsorted hab
    simpa using this
  rcases hle.lt_or_eq with hlt | heq
  · exact Or.inl hlt
  · refine Or.inr ⟨heq, ?_⟩
    have hne : a ≠ b := by
      have := hnodupR.sublist hab
      simpa using this
    rcases lt_or_gt_of_ne hne with h | h
    · exact h
    · exfalso
      have haR : a ∈ L := List.mem_mergeSort.mp (hab.subset (by simp))
      have hbR : b ∈ L := List.mem_mergeSort.mp (hab.subset (by simp))
      have hba := List.pair_sublist_mergeSort (a := b) (b := a) htrans htot
        (by rw [decide_eq_true_eq, heq]) (pair_sublist_of_lt hLlt hbR haR h)
      exact not_both_orders hnodupR hne hab hba

/-! ## R-N1 (d) and the witness -/

/-- R-N1 (d). -/
theorem N1_rules_differ : N1_rules_differ_stmt := by
  refine ⟨?_, ?_⟩
  · intro pv h
    have hn : nPresent pv = 5 := by
      simp [nPresent, pv]
    have hceil : ⌈(1 / 4 : ℚ) * ((5 : ℕ) : ℚ)⌉₊ = 2 := by
      rw [Nat.ceil_eq_iff (by norm_num)]
      norm_num
    have h1 := (N1_first_pass 5 pv).1
    have h2 := (N1_trx 5 pv).1
    rw [h, h2, hn, hceil] at h1
    norm_num at h1
  · intro pv
    have hceil : ⌈(1 / 4 : ℚ) * ((4 : ℕ) : ℚ)⌉₊ = 1 := by
      norm_num
    refine ⟨?_, ?_⟩
    · simp only [trxFirstPassControls, Finset.mem_filter, Finset.mem_univ, true_and]
      rw [hceil]
      unfold orderPos
      have : 0 < (univ.filter fun j => orderBefore pv j 0).card :=
        Finset.card_pos.mpr ⟨1, by
          simp only [Finset.mem_filter, Finset.mem_univ, true_and]
          exact Or.inl ⟨by simp [pv], by simp [pv]⟩⟩
      omega
    · intro hmem
      rw [firstPassControls, Finset.mem_filter] at hmem
      have h0 := hmem.2.1
      simp [pv] at h0

/-- Witness for `N1_stable_stmt`: two rows, every `q` equal to 1, cutoff 0, one control asked for. -/
theorem N1_stable_witness : N1_stable_witness_stmt :=
  ⟨2, fun _ => 1, 0, 1, fun _ _ => rfl, fun _ => by norm_num, by norm_num, by norm_num⟩

end MultiScholaRProofs.RuvControls
