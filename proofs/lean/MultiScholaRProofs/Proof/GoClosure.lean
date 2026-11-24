module

public import MultiScholaRProofs.Spec.GoClosure
import Mathlib.Tactic

/-!
# Proofs for target 3 (KIT-193 section 6.3)

Theorems whose types are the frozen statements of `Spec/GoClosure.lean`. The search `bfs` is proved correct by a loop
invariant (every seen pair is a (term, ancestor) pair; the start pairs are seen; every seen pair that has left the
frontier has all its successors seen) and a potential that falls by at least one on each pass that does not stop.
-/

public section

namespace MultiScholaRProofs.GoClosure

open Finset Relation

section bfs

variable {T L : Type} [DecidableEq T] [DecidableEq L]

theorem mem_keptEdges (Γ : GoGraph T L) (Ls : Finset L) (x y : T) :
    (x, y) ∈ keptEdges Γ Ls ↔ parentRel Γ Ls x y := by
  simp only [keptEdges, List.mem_map, List.mem_filter, Bool.and_eq_true, decide_eq_true_eq, parentRel]
  constructor
  · rintro ⟨⟨c, p, l⟩, ⟨hmem, hl, hp⟩, he⟩
    simp only [Prod.mk.injEq] at he
    obtain ⟨rfl, rfl⟩ := he
    exact ⟨l, hmem, hl, hp⟩
  · rintro ⟨l, hmem, hl, hp⟩
    exact ⟨⟨x, y, l⟩, ⟨hmem, hl, hp⟩, rfl⟩

theorem mem_bfsStep_fst (E : List (T × T)) (st : Finset (T × T) × Finset (T × T)) (z : T × T) :
    z ∈ (bfsStep E st).1 ↔ z ∈ st.1 ∨ ∃ x ∈ st.2, (x.2, z.2) ∈ E ∧ z.1 = x.1 := by
  obtain ⟨z1, z2⟩ := z
  simp only [bfsStep, mem_union, mem_sdiff, mem_biUnion, List.mem_toFinset, List.mem_map, List.mem_filter,
    decide_eq_true_eq]
  constructor
  · rintro (h | ⟨⟨x, hx, e, ⟨he, hex⟩, hz⟩, _⟩)
    · exact Or.inl h
    · right
      refine ⟨x, hx, ?_, ?_⟩
      · simp only [Prod.mk.injEq] at hz
        rw [← hz.2, ← hex]
        exact he
      · simp only [Prod.mk.injEq] at hz
        exact hz.1.symm
  · rintro (h | ⟨x, hx, he, hz⟩)
    · exact Or.inl h
    · by_cases hs : (z1, z2) ∈ st.1
      · exact Or.inl hs
      · right
        refine ⟨⟨x, hx, (x.2, z2), ⟨he, rfl⟩, ?_⟩, hs⟩
        simp [hz]

theorem mem_bfsStep_snd (E : List (T × T)) (st : Finset (T × T) × Finset (T × T)) (z : T × T) :
    z ∈ (bfsStep E st).2 ↔ z ∉ st.1 ∧ ∃ x ∈ st.2, (x.2, z.2) ∈ E ∧ z.1 = x.1 := by
  obtain ⟨z1, z2⟩ := z
  simp only [bfsStep, mem_sdiff, mem_biUnion, List.mem_toFinset, List.mem_map, List.mem_filter, decide_eq_true_eq]
  constructor
  · rintro ⟨⟨x, hx, e, ⟨he, hex⟩, hz⟩, hs⟩
    refine ⟨hs, x, hx, ?_, ?_⟩
    · simp only [Prod.mk.injEq] at hz
      rw [← hz.2, ← hex]
      exact he
    · simp only [Prod.mk.injEq] at hz
      exact hz.1.symm
  · rintro ⟨hs, x, hx, he, hz⟩
    refine ⟨⟨x, hx, (x.2, z2), ⟨he, rfl⟩, ?_⟩, hs⟩
    simp [hz]

/-- The loop invariant, for a successor relation given by the list `E` and start terms `terms`. -/
structure Inv (E : List (T × T)) (terms : Finset T) (st : Finset (T × T) × Finset (T × T)) : Prop where
  front_sub : st.2 ⊆ st.1
  sound : ∀ z ∈ st.1, z.1 ∈ terms ∧ ReflTransGen (fun x y => (x, y) ∈ E) z.1 z.2
  start : ∀ t ∈ terms, (t, t) ∈ st.1
  closed : ∀ z ∈ st.1, z ∉ st.2 → ∀ y, (z.2, y) ∈ E → (z.1, y) ∈ st.1
  within : st.1 ⊆ terms ×ˢ nodes E terms

theorem inv_start (E : List (T × T)) (terms : Finset T) : Inv E terms (bfsStart terms) := by
  refine ⟨subset_refl _, ?_, ?_, ?_, ?_⟩
  · intro z hz
    simp only [bfsStart, mem_image] at hz
    obtain ⟨t, ht, rfl⟩ := hz
    exact ⟨ht, ReflTransGen.refl⟩
  · intro t ht
    simp only [bfsStart, mem_image]
    exact ⟨t, ht, rfl⟩
  · intro z hz hz'
    exact absurd hz hz'
  · intro z hz
    simp only [bfsStart, mem_image] at hz
    obtain ⟨t, ht, rfl⟩ := hz
    simp only [mem_product, nodes, mem_union]
    exact ⟨ht, Or.inl ht⟩

theorem inv_step (E : List (T × T)) (terms : Finset T) (st : Finset (T × T) × Finset (T × T))
    (h : Inv E terms st) : Inv E terms (bfsStep E st) := by
  refine ⟨?_, ?_, ?_, ?_, ?_⟩
  · intro z hz
    rw [mem_bfsStep_snd] at hz
    rw [mem_bfsStep_fst]
    exact Or.inr hz.2
  · intro z hz
    rw [mem_bfsStep_fst] at hz
    rcases hz with hz | ⟨x, hx, he, hz1⟩
    · exact h.sound z hz
    · obtain ⟨hterm, hanc⟩ := h.sound x (h.front_sub hx)
      rw [hz1]
      exact ⟨hterm, hanc.tail he⟩
  · intro t ht
    rw [mem_bfsStep_fst]
    exact Or.inl (h.start t ht)
  · intro z hz hzf y hy
    rw [mem_bfsStep_snd] at hzf
    have hz' : z ∈ st.1 := by
      rw [mem_bfsStep_fst] at hz
      rcases hz with hz | ⟨x, hx, he, hz1⟩
      · exact hz
      · by_contra hs
        exact hzf ⟨hs, x, hx, he, hz1⟩
    rw [mem_bfsStep_fst]
    by_cases hfront : z ∈ st.2
    · by_cases hseen : (z.1, y) ∈ st.1
      · exact Or.inl hseen
      · exact Or.inr ⟨z, hfront, hy, rfl⟩
    · exact Or.inl (h.closed z hz' hfront y hy)
  · intro z hz
    rw [mem_bfsStep_fst] at hz
    rcases hz with hz | ⟨x, hx, he, hz1⟩
    · exact h.within hz
    · have hxw := h.within (h.front_sub hx)
      simp only [mem_product, nodes, mem_union, List.mem_toFinset, List.mem_map] at hxw ⊢
      refine ⟨hz1 ▸ hxw.1, Or.inr ⟨(x.2, z.2), he, rfl⟩⟩

theorem inv_loop (E : List (T × T)) (terms : Finset T) :
    ∀ (f : ℕ) (st : Finset (T × T) × Finset (T × T)), Inv E terms st → Inv E terms (bfsLoop E f st)
  | 0, st, h => by simpa [bfsLoop] using h
  | f + 1, st, h => by
    simp only [bfsLoop]
    split_ifs
    · exact h
    · exact inv_loop E terms f _ (inv_step E terms st h)

/-- The potential: pairs not yet seen, plus one while the frontier is not empty. -/
def potential (E : List (T × T)) (terms : Finset T) (st : Finset (T × T) × Finset (T × T)) : ℕ :=
  (terms ×ˢ nodes E terms).card - st.1.card + (if st.2 = ∅ then 0 else 1)

theorem potential_step (E : List (T × T)) (terms : Finset T) (st : Finset (T × T) × Finset (T × T))
    (h : Inv E terms st) (hne : st.2 ≠ ∅) :
    potential E terms (bfsStep E st) + 1 ≤ potential E terms st := by
  have hstep := inv_step E terms st h
  have hsub : st.1 ⊆ (bfsStep E st).1 := fun z hz => (mem_bfsStep_fst E st z).2 (Or.inl hz)
  have hcardU : (bfsStep E st).1.card ≤ (terms ×ˢ nodes E terms).card := card_le_card hstep.within
  simp only [potential, hne, ↓reduceIte]
  by_cases hfresh : (bfsStep E st).2 = ∅
  · simp only [hfresh, ↓reduceIte]
    have := card_le_card hsub
    omega
  · simp only [hfresh, ↓reduceIte]
    obtain ⟨z, hz⟩ := nonempty_iff_ne_empty.2 hfresh
    have hzn : z ∉ st.1 := ((mem_bfsStep_snd E st z).1 hz).1
    have hz1 : z ∈ (bfsStep E st).1 := hstep.front_sub hz
    have hlt : st.1.card < (bfsStep E st).1.card := card_lt_card ⟨hsub, fun hss => hzn (hss hz1)⟩
    omega

theorem loop_terminates (E : List (T × T)) (terms : Finset T) :
    ∀ (f : ℕ) (st : Finset (T × T) × Finset (T × T)), Inv E terms st → potential E terms st ≤ f →
      (bfsLoop E f st).2 = ∅
  | 0, st, h, hp => by
    simp only [bfsLoop]
    by_contra hne
    simp only [potential, hne, ↓reduceIte] at hp
    omega
  | f + 1, st, h, hp => by
    simp only [bfsLoop]
    split_ifs with hne
    · exact hne
    · apply loop_terminates E terms f _ (inv_step E terms st h)
      have := potential_step E terms st h hne
      omega

theorem bfsState_snd (Γ : GoGraph T L) (Ls : Finset L) (terms : Finset T) : (bfsState Γ Ls terms).2 = ∅ := by
  apply loop_terminates _ terms _ _ (inv_start _ terms)
  simp only [potential, fuel, card_product]
  split_ifs <;> omega

theorem anc_iff (Γ : GoGraph T L) (Ls : Finset L) (t a : T) :
    anc Γ Ls t a ↔ ReflTransGen (fun x y => (x, y) ∈ keptEdges Γ Ls) t a := by
  have : (fun x y => (x, y) ∈ keptEdges Γ Ls) = parentRel Γ Ls := by
    funext x y
    exact propext (mem_keptEdges Γ Ls x y)
  rw [this]
  rfl

theorem mem_bfs (Γ : GoGraph T L) (Ls : Finset L) (terms : Finset T) (t a : T) :
    (t, a) ∈ bfs Γ Ls terms ↔ t ∈ terms ∧ anc Γ Ls t a := by
  have hinv := inv_loop (keptEdges Γ Ls) terms (fuel (keptEdges Γ Ls) terms) _ (inv_start _ terms)
  have hempty := bfsState_snd Γ Ls terms
  rw [anc_iff]
  constructor
  · intro h
    exact hinv.sound (t, a) h
  · rintro ⟨ht, hanc⟩
    induction hanc with
    | refl => exact hinv.start t ht
    | tail _ hyz ih =>
      refine hinv.closed _ ih ?_ _ hyz
      unfold bfsState at hempty
      rw [hempty]
      exact notMem_empty _

end bfs

/-- (3a) termination. -/
theorem bfs_terminates : bfs_terminates_stmt := by
  intro T L _ _ Γ Ls terms
  exact bfsState_snd Γ Ls terms

/-- (3a) correctness. -/
theorem bfs_correct : bfs_correct_stmt := by
  intro T L _ _ Γ Ls terms t a
  exact mem_bfs Γ Ls terms t a

theorem anc_all {T L : Type} (Γ : GoGraph T L) (Ls : Finset L) (t : T) (h : anc Γ Ls t Γ.all) : t = Γ.all := by
  rcases ReflTransGen.cases_tail h with h | ⟨c, _, l, _, _, hp⟩
  · exact h.symm
  · exact absurd rfl hp

/-- (3b) the code's propagation is the closure image. -/
theorem propagate_fin : propagate_fin_stmt := by
  intro T G L Asp _ _ _ _ Γ Ls aspectOf o A
  ext ⟨a, g⟩
  simp only [propagateFin, coe_image, coe_filter, Set.mem_image, Set.mem_ofPred_eq, mem_product, mem_filter,
    propagate, kept, mem_coe, Prod.exists, Prod.mk.injEq, decide_eq_true_eq]
  constructor
  · rintro ⟨t, g', t', a', ⟨⟨⟨hA, hasp⟩, hpair⟩, rfl⟩, rfl, rfl⟩
    rw [mem_bfs] at hpair
    exact ⟨t, ⟨hA, hasp⟩, hpair.2⟩
  · rintro ⟨t, ⟨hA, hasp⟩, hanc⟩
    refine ⟨t, g, t, a, ⟨⟨⟨hA, hasp⟩, ?_⟩, rfl⟩, rfl, rfl⟩
    rw [mem_bfs]
    refine ⟨?_, hanc⟩
    simp only [mem_image, mem_filter, Prod.exists, exists_and_right, exists_eq_right]
    exact ⟨⟨g, hA⟩, hasp⟩

/-- (3b) idempotent. -/
theorem idempotent : idempotent_stmt := by
  intro T G L Asp _ _ _ _ Γ Ls aspectOf o A
  ext ⟨a, g⟩
  simp only [propagate, kept, Set.mem_ofPred_eq]
  constructor
  · rintro ⟨t', ⟨⟨t, ⟨hA, hasp⟩, h1⟩, _⟩, h2⟩
    exact ⟨t, ⟨hA, hasp⟩, h1.trans h2⟩
  · rintro ⟨t, ⟨hA, hasp⟩, h⟩
    exact ⟨t, ⟨⟨t, ⟨hA, hasp⟩, ReflTransGen.refl⟩, hasp⟩, h⟩

/-- (3b) monotone. -/
theorem monotone : monotone_stmt := by
  intro T G L Asp _ _ _ _ Γ Ls aspectOf o A B hAB x hx
  obtain ⟨t, ⟨hA, hasp⟩, h⟩ := hx
  exact ⟨t, ⟨hAB hA, hasp⟩, h⟩

/-- (3b) extensive on what is kept. -/
theorem extensive_kept : extensive_kept_stmt := by
  intro T G L Asp _ _ _ _ Γ Ls aspectOf o A x hx
  exact ⟨x.1, hx, ReflTransGen.refl⟩

/-- (3b) not extensive on the input: an annotation of a term GO.db does not hold is dropped. -/
theorem not_extensive : not_extensive_stmt := by
  refine ⟨⟨[], 0⟩, ∅, fun _ => none, 0, {(1, 1)}, ?_⟩
  intro h
  obtain ⟨t, ⟨_, hasp⟩, _⟩ := h (Set.mem_singleton (1, 1))
  exact absurd hasp (by simp)

/-- (3c) more labels, more pairs. -/
theorem labels_mono : labels_mono_stmt := by
  intro T G L Asp _ _ _ _ Γ Ls Ls' aspectOf o A hL x hx
  obtain ⟨t, hk, h⟩ := hx
  refine ⟨t, hk, ReflTransGen.mono (fun c p hcp => ?_) _ _ h⟩
  obtain ⟨l, hmem, hl, hp⟩ := hcp
  exact ⟨l, hmem, hL hl, hp⟩

/-- (3d) roots removed. -/
theorem root_removed : root_removed_stmt := by
  intro T G L Asp _ _ _ _ Γ Ls aspectOf o root A
  refine ⟨fun x hx => hx.2, fun t g ht => ?_⟩
  exact ⟨fun h => h.1, fun h => ⟨h, ht⟩⟩

/-- (3d) `all` is never a term of a propagated pair. -/
theorem all_absent : all_absent_stmt := by
  intro T G L Asp _ _ _ _ Γ Ls aspectOf o A hall g h
  obtain ⟨t, ⟨_, hasp⟩, hanc⟩ := h
  have ht := anc_all Γ Ls t hanc
  rw [ht, hall] at hasp
  exact absurd hasp (by simp)

/-- (3d) of `goAncestorPairs` alone. -/
theorem bfs_all : bfs_all_stmt := by
  intro T L _ _ Γ Ls terms t h
  rw [mem_bfs] at h
  have ht := anc_all Γ Ls t h.2
  exact ⟨ht, ht ▸ h.1⟩

/-- (3e) ancestors stay in the aspect when no followed edge leaves it. -/
theorem aspect_closed : aspect_closed_stmt := by
  intro T L Asp _ _ _ Γ Ls aspectOf hedge t a h
  induction h with
  | refl => rfl
  | tail _ hbc ih =>
    obtain ⟨l, hmem, hl, hp⟩ := hbc
    rw [hedge _ _ l hmem hl hp, ih]

/-! Witnesses -/

theorem monotone_witness : monotone_witness_stmt := by
  refine ⟨⟨[], 0⟩, ∅, fun _ => some 0, 0, ∅, {(1, 1)}, Set.empty_subset _, ?_⟩
  intro h
  have hx : ((1 : ℕ), (1 : ℕ)) ∈ propagate (⟨[], 0⟩ : GoGraph ℕ ℕ) ∅ (fun _ => some 0) 0 {(1, 1)} :=
    ⟨1, ⟨Set.mem_singleton _, rfl⟩, ReflTransGen.refl⟩
  rw [← h] at hx
  obtain ⟨t, ⟨hA, _⟩, _⟩ := hx
  exact hA

theorem labels_mono_witness : labels_mono_witness_stmt := by
  refine ⟨⟨[(0, 1, 7)], 99⟩, ∅, {7}, fun _ => some 0, 0, {(0, 5)}, Finset.empty_subset _, ?_⟩
  intro h
  have hx : ((1 : ℕ), (5 : ℕ)) ∈ propagate (⟨[(0, 1, 7)], 99⟩ : GoGraph ℕ ℕ) {7} (fun _ => some 0) 0 {(0, 5)} :=
    ⟨0, ⟨Set.mem_singleton _, rfl⟩, ReflTransGen.single ⟨7, by simp, by simp, by decide⟩⟩
  rw [← h] at hx
  obtain ⟨t, ⟨hA, _⟩, hanc⟩ := hx
  simp only [Set.mem_singleton_iff, Prod.mk.injEq] at hA
  obtain ⟨rfl, -⟩ := hA
  rcases ReflTransGen.cases_tail hanc with h1 | ⟨c, _, l, _, hl, _⟩
  · exact absurd h1 (by decide)
  · exact absurd hl (by simp)

theorem all_absent_witness : all_absent_witness_stmt := by
  refine ⟨⟨[(0, 99, 7)], 99⟩, {7}, fun t => if t = 99 then none else some 0, 0, {(0, 5)}, by simp,
    ⟨0, 7, by simp, by simp⟩, ⟨(0, 5), 0, ⟨Set.mem_singleton _, by simp⟩, ReflTransGen.refl⟩⟩

theorem aspect_closed_witness : aspect_closed_witness_stmt := by
  refine ⟨⟨[(0, 1, 7)], 99⟩, {7}, fun t => if t ≤ 1 then some 0 else some 1, ?_, ⟨0, 1, 7, by simp, by simp,
    by decide, by decide⟩, ⟨0, 2, by simp⟩⟩
  intro c p l hmem _ _
  simp only [List.mem_singleton, Prod.mk.injEq] at hmem
  obtain ⟨rfl, rfl, -⟩ := hmem
  simp

end MultiScholaRProofs.GoClosure
