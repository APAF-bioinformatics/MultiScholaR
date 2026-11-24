module

public import MultiScholaRProofs.Spec.QvalueBh
import Mathlib.Tactic

/-!
# Target 4 (KIT-193 section 6.4): the counterexample that amended the freeze

`bh_sorted_eq_stmt` as frozen at 4183eb10 (R's transcription `bhSorted` equals the closed form `bh` for every `p`) was
false: with negative tied values the later member of a tie group gets a smaller cumulative minimum. The statements
gate's re-review 1 added the hypothesis `∀ j, 0 ≤ p j` and froze the counter-instance as
`bh_sorted_neg_ties_counterexample_stmt` (`K = ℚ`, `m = 2`, `p = (-1, -1)`), proved here.
-/

public section

namespace MultiScholaRProofs.QvalueBh

open Finset

section general

variable {K : Type} [Field K] [LinearOrder K] [IsStrictOrderedRing K]

-- A style linter only (an instance argument a lemma does not use); it has no bearing on what is proved.
set_option linter.unusedSectionVars false

theorem bh_le_iff {m : ℕ} (p : Fin m → K) (i : Fin m) (α : K) :
    bh p i ≤ α ↔ 1 ≤ α ∨ ∃ j, p i ≤ p j ∧ bhTerm p j ≤ α := by
  unfold bh
  rw [min_le_iff, Finset.fold_min_le]
  constructor
  · rintro (h | h | ⟨j, hj, hle⟩)
    · exact Or.inl h
    · exact Or.inr ⟨i, le_rfl, h⟩
    · exact Or.inr ⟨j, (Finset.mem_filter.1 hj).2, hle⟩
  · rintro (h | ⟨j, hij, hle⟩)
    · exact Or.inl h
    · exact Or.inr (Or.inr ⟨j, Finset.mem_filter.2 ⟨Finset.mem_univ _, hij⟩, hle⟩)

theorem bh_lt_iff {m : ℕ} (p : Fin m → K) (i : Fin m) (α : K) :
    bh p i < α ↔ 1 < α ∨ ∃ j, p i ≤ p j ∧ bhTerm p j < α := by
  unfold bh
  rw [min_lt_iff, Finset.fold_min_lt]
  constructor
  · rintro (h | h | ⟨j, hj, hle⟩)
    · exact Or.inl h
    · exact Or.inr ⟨i, le_rfl, h⟩
    · exact Or.inr ⟨j, (Finset.mem_filter.1 hj).2, hle⟩
  · rintro (h | ⟨j, hij, hle⟩)
    · exact Or.inl h
    · exact Or.inr (Or.inr ⟨j, Finset.mem_filter.2 ⟨Finset.mem_univ _, hij⟩, hle⟩)

theorem le_bh_iff {m : ℕ} (p : Fin m → K) (i : Fin m) (α : K) :
    α ≤ bh p i ↔ α ≤ 1 ∧ ∀ j, p i ≤ p j → α ≤ bhTerm p j := by
  unfold bh
  rw [le_min_iff, Finset.le_fold_min]
  constructor
  · rintro ⟨h1, hi, hall⟩
    exact ⟨h1, fun j hj => hall j (Finset.mem_filter.2 ⟨Finset.mem_univ _, hj⟩)⟩
  · rintro ⟨h1, hall⟩
    exact ⟨h1, hall i le_rfl, fun j hj => hall j (Finset.mem_filter.1 hj).2⟩

theorem cntLe_pos {m : ℕ} (p : Fin m → K) (j : Fin m) : 0 < cntLe p (p j) :=
  Finset.card_pos.2 ⟨j, Finset.mem_filter.2 ⟨Finset.mem_univ _, le_rfl⟩⟩

theorem cntLe_le_m {m : ℕ} (p : Fin m → K) (v : K) : cntLe p v ≤ m := by
  unfold cntLe
  calc (Finset.univ.filter fun l => p l ≤ v).card ≤ Finset.univ.card := Finset.card_filter_le _ _
    _ = m := by simp

theorem cntLe_mono {m : ℕ} (p : Fin m → K) {v w : K} (h : v ≤ w) : cntLe p v ≤ cntLe p w :=
  Finset.card_le_card fun _ hl => Finset.mem_filter.2 ⟨Finset.mem_univ _, le_trans (Finset.mem_filter.1 hl).2 h⟩

theorem bhTerm_nonneg {m : ℕ} (p : Fin m → K) (j : Fin m) (h : 0 ≤ p j) : 0 ≤ bhTerm p j := by
  unfold bhTerm
  positivity

theorem bh_nonneg {m : ℕ} (p : Fin m → K) (i : Fin m) (h : ∀ j, 0 ≤ p j) : 0 ≤ bh p i :=
  (le_bh_iff p i 0).2 ⟨zero_le_one, fun j _ => bhTerm_nonneg p j (h j)⟩

theorem bh_stepup_gen {m : ℕ} (p : Fin m → K) (i : Fin m) (α : K) (hp : UnitP p) :
    bh p i ≤ α ↔ ∃ k : ℕ, 1 ≤ k ∧ p i ≤ k * α / m ∧ k ≤ cntLe p (k * α / m) := by
  have hm : 0 < m := Fin.pos i
  have hm' : (0 : K) < m := by exact_mod_cast hm
  rw [bh_le_iff]
  constructor
  · rintro (h1 | ⟨j, hij, hle⟩)
    · refine ⟨m, hm, ?_, ?_⟩
      · rw [mul_div_cancel_left₀ _ hm'.ne']
        exact le_trans (hp i).2 h1
      · rw [mul_div_cancel_left₀ _ hm'.ne']
        apply le_of_eq
        unfold cntLe
        rw [Finset.filter_true_of_mem (fun l _ => le_trans (hp l).2 h1)]
        simp
    · set c := cntLe p (p j) with hc
      have hcpos : 0 < c := cntLe_pos p j
      have hc' : (0 : K) < c := by exact_mod_cast hcpos
      have hpj : p j ≤ c * α / m := by
        unfold bhTerm at hle
        rw [← hc, div_le_iff₀ hc'] at hle
        rw [le_div_iff₀ hm']
        linarith
      exact ⟨c, hcpos, le_trans hij hpj, cntLe_mono p hpj⟩
  · rintro ⟨k, hk, hpi, hkc⟩
    right
    set t := (k : K) * α / m with ht
    obtain ⟨j, hjT, hjmax⟩ := Finset.exists_max_image (Finset.univ.filter fun l => p l ≤ t) p
      ⟨i, Finset.mem_filter.2 ⟨Finset.mem_univ _, hpi⟩⟩
    have hjt : p j ≤ t := (Finset.mem_filter.1 hjT).2
    refine ⟨j, hjmax i (Finset.mem_filter.2 ⟨Finset.mem_univ _, hpi⟩), ?_⟩
    have hsub : cntLe p t ≤ cntLe p (p j) :=
      Finset.card_le_card fun l hl => Finset.mem_filter.2 ⟨Finset.mem_univ _, hjmax l hl⟩
    have hk' : (0 : K) < k := by exact_mod_cast hk
    have hkc' : (k : K) ≤ cntLe p (p j) := by exact_mod_cast le_trans hkc hsub
    have hpj0 : 0 ≤ p j := (hp j).1
    unfold bhTerm
    calc (m : K) * p j / cntLe p (p j) ≤ m * p j / k :=
          div_le_div_of_nonneg_left (by positivity) hk' hkc'
      _ ≤ m * t / k := by
          apply div_le_div_of_nonneg_right _ hk'.le
          exact mul_le_mul_of_nonneg_left hjt hm'.le
      _ = α := by
          rw [ht]
          field_simp

/-! ### R's sort-and-cumulative-minimum form equals the closed form on non-negative values

At 0-based position `r` of the decreasing order `L` the cumulative minimum is the minimum over `s ≤ r` of
`m / (m - s) * p L[s]`. Each such `L[s]` has `p ≥ p i`, and at least `m - s` values lie at or below it (the positions
from `s` on), so its closed-form term `m p / #{≤ p}` is at most the sorted term (this needs `p ≥ 0`). Conversely, for
`j` with `p j ≥ p i` and `c = #{l : p l ≤ p j}`, the position `s = m - c` is at most `r` and carries a value at most
`p j`, so its sorted term `m / c * p L[s]` is at most `bhTerm p j`. -/

theorem length_cumMin : ∀ l : List K, (cumMin l).length = l.length
  | [] => rfl
  | a :: l => by simp [cumMin, length_cumMin l]

/-- The prefix minimum is the minimum of the prefix: `α ≤ cm[r] ↔ α ≤ l[s]` for every `s ≤ r`. -/
theorem le_cumMin_iff : ∀ (l : List K) (r : ℕ) (hr : r < (cumMin l).length) (α : K),
    α ≤ (cumMin l)[r] ↔ ∀ s (hs : s < l.length), s ≤ r → α ≤ l[s]
  | [], r, hr, α => by simp [cumMin] at hr
  | a :: l, 0, hr, α => by
      simp only [cumMin, List.getElem_cons_zero]
      constructor
      · intro h s hs hsr
        have hs0 : s = 0 := by omega
        subst hs0
        simpa using h
      · intro h
        have h0 := h 0 (by simp) le_rfl
        rwa [List.getElem_cons_zero] at h0
  | a :: l, r + 1, hr, α => by
      have hr' : r < (cumMin l).length := by
        simp only [cumMin, List.length_cons, List.length_map] at hr
        omega
      simp only [cumMin, List.getElem_cons_succ, List.getElem_map, le_min_iff]
      rw [le_cumMin_iff l r hr' α]
      constructor
      · rintro ⟨ha, hl⟩ s hs hsr
        cases s with
        | zero => simpa using ha
        | succ s => simpa using hl s (by simpa using hs) (by omega)
      · intro h
        refine ⟨?_, fun s hs hsr => ?_⟩
        · have h0 := h 0 (by simp) (by omega)
          rwa [List.getElem_cons_zero] at h0
        · have h1 := h (s + 1) (by simpa using hs) (by omega)
          rwa [List.getElem_cons_succ] at h1

/-- `bhSorted = bh` on non-negative values, whatever order the sort gives tied values. -/
theorem bhSorted_eq_bh {m : ℕ} (p : Fin m → K) (hp : ∀ j, 0 ≤ p j) : bhSorted p = bh p := by
  funext i
  unfold bhSorted
  set L := (List.finRange m).mergeSort fun a b => decide (p b ≤ p a) with hL
  have hperm : L.Perm (List.finRange m) := List.mergeSort_perm _ _
  have hlen : L.length = m := by rw [hperm.length_eq, List.length_finRange]
  have hnd : L.Nodup := hperm.nodup_iff.2 (List.nodup_finRange m)
  have hmem : ∀ l : Fin m, l ∈ L := fun l => hperm.mem_iff.2 (List.mem_finRange l)
  have hsort : ∀ (a b : ℕ) (ha : a < L.length) (hb : b < L.length), a < b → p L[b] ≤ p L[a] := by
    have hpw := List.pairwise_mergeSort (le := fun a b => decide (p b ≤ p a))
      (fun a b c hab hbc => by
        simp only [decide_eq_true_eq] at hab hbc ⊢
        exact le_trans hbc hab)
      (fun a b => by simp [le_total]) (List.finRange m)
    rw [List.pairwise_iff_getElem] at hpw
    intro a b ha hb hab
    simpa using hpw a b ha hb hab
  have hsort' : ∀ (a b : ℕ) (ha : a < L.length) (hb : b < L.length), a ≤ b → p L[b] ≤ p L[a] := by
    intro a b ha hb hab
    rcases hab.lt_or_eq with h | h
    · exact hsort a b ha hb h
    · subst h
      exact le_rfl
  have hget : ∀ (u : ℕ) (hu : u < L.length), L.getD u i = L[u] := by
    intro u hu
    rw [List.getD_eq_getElem?_getD, List.getElem?_eq_getElem hu, Option.getD_some]
  -- (A) a value at position `t` has at least `m - t` values at or below it
  have hA : ∀ (t : ℕ) (ht : t < L.length) (v : K), p L[t] ≤ v → m - t ≤ cntLe p v := by
    intro t ht v hv
    have key := card_le_card_of_injOn (s := Finset.Ico t m) (t := univ.filter fun l => p l ≤ v)
      (fun u => L.getD u i) (fun u hu => ?_) (fun u hu w hw huw => ?_)
    · simpa [cntLe] using key
    · rw [Finset.coe_Ico, Set.mem_Ico] at hu
      have hu' : u < L.length := by omega
      rw [Finset.coe_filter]
      simp only [mem_univ, true_and, Set.mem_ofPred_eq]
      rw [hget u hu']
      exact le_trans (hsort' t u ht hu' hu.1) hv
    · rw [Finset.coe_Ico, Set.mem_Ico] at hu hw
      have hu' : u < L.length := by omega
      have hw' : w < L.length := by omega
      simp only [hget u hu', hget w hw'] at huw
      exact (hnd.getElem_inj_iff).1 huw
  -- (B) values strictly below the value at position `s` sit after `s`
  have hB : ∀ (s : ℕ) (hs : s < L.length) (v : K), v < p L[s] → cntLe p v + s + 1 ≤ m := by
    intro s hs v hv
    have key := card_le_card_of_injOn (s := univ.filter fun l => p l ≤ v) (t := Finset.Ico (s + 1) m)
      (fun l => L.idxOf l) (fun l hl => ?_) (fun a _ b _ hab => (List.idxOf_inj (hmem a)).1 hab)
    · rw [Nat.card_Ico] at key
      unfold cntLe
      omega
    · rw [Finset.coe_filter] at hl
      simp only [mem_univ, true_and, Set.mem_ofPred_eq] at hl
      have hlt : L.idxOf l < L.length := List.idxOf_lt_length_of_mem (hmem l)
      rw [Finset.coe_Ico, Set.mem_Ico]
      beta_reduce
      refine ⟨?_, by omega⟩
      by_contra hc
      have h1 := hsort' (L.idxOf l) s hlt hs (by omega)
      rw [List.getElem_idxOf hlt] at h1
      exact absurd (lt_of_lt_of_le hv h1) (not_lt.2 hl)
  have hrlt : L.idxOf i < L.length := List.idxOf_lt_length_of_mem (hmem i)
  have hLr : L[L.idxOf i] = i := List.getElem_idxOf hrlt
  set r := L.idxOf i with hr
  set sc := L.zipIdx.map fun x => (m : K) / ((m - x.2 : ℕ) : K) * p x.1 with hsc
  have hsclen : sc.length = L.length := by simp [hsc]
  have hsc_get : ∀ (s : ℕ) (hs : s < sc.length),
      sc[s] = (m : K) / ((m - s : ℕ) : K) * p (L[s]'(by omega)) := by
    intro s hs
    simp [hsc, List.getElem_zipIdx]
  have hcm : r < (cumMin sc).length := by rw [length_cumMin]; omega
  rw [List.getD_eq_getElem?_getD, List.getElem?_eq_getElem hcm, Option.getD_some]
  apply eq_of_forall_le_iff
  intro α
  rw [le_min_iff, le_cumMin_iff sc r hcm α, le_bh_iff]
  refine and_congr_right fun _ => ⟨fun h j hij => ?_, fun h s hs hsr => ?_⟩
  · set c := cntLe p (p j) with hc
    have hc1 : 1 ≤ c := cntLe_pos p j
    have hcm' : c ≤ m := cntLe_le_m p _
    have hsr : m - c ≤ r := by
      have := hA r hrlt (p j) (by rw [hLr]; exact hij)
      omega
    have hslt : m - c < L.length := by omega
    have hval : p L[m - c] ≤ p j := by
      by_contra hcon
      have := hB (m - c) hslt (p j) (lt_of_not_ge hcon)
      omega
    have h1 := h (m - c) (by omega) hsr
    rw [hsc_get] at h1
    have hmc : ((m - (m - c) : ℕ) : K) = (c : K) := by
      congr 1
      omega
    rw [hmc] at h1
    have hcpos : (0 : K) < c := by exact_mod_cast hc1
    unfold bhTerm
    rw [← hc]
    calc α ≤ (m : K) / c * p L[m - c] := h1
      _ ≤ (m : K) / c * p j := mul_le_mul_of_nonneg_left hval (by positivity)
      _ = m * p j / c := by ring
  · have hsL : s < L.length := by omega
    have hij : p i ≤ p L[s] := by
      have := hsort' s r hsL hrlt hsr
      rw [hLr] at this
      exact this
    have h1 := h L[s] hij
    have hcnt : m - s ≤ cntLe p (p L[s]) := hA s hsL _ le_rfl
    rw [hsc_get]
    have hms : (0 : K) < ((m - s : ℕ) : K) := by
      have : 0 < m - s := by omega
      exact_mod_cast this
    have hcK : ((m - s : ℕ) : K) ≤ (cntLe p (p L[s]) : K) := by exact_mod_cast hcnt
    have hp0 : 0 ≤ p L[s] := hp _
    unfold bhTerm at h1
    calc α ≤ m * p L[s] / cntLe p (p L[s]) := h1
      _ ≤ m * p L[s] / ((m - s : ℕ) : K) := div_le_div_of_nonneg_left (by positivity) hms hcK
      _ = m / ((m - s : ℕ) : K) * p L[s] := by ring

end general

/-- (4) R's `p.adjust` form equals the closed form on non-negative values (the corrected statement). -/
theorem bh_sorted_eq : bh_sorted_eq_stmt := fun _ _ _ _ _ p hp => bhSorted_eq_bh p hp

/-- Witness for `bh_sorted_eq_stmt`: `p = (1, 1, 2)`. -/
theorem bh_sorted_eq_witness : bh_sorted_eq_witness_stmt :=
  ⟨![1, 1, 2], fun j => by fin_cases j <;> norm_num, by norm_num, by norm_num, by norm_num⟩

/-- The counter-instance: two tied p-values of `-1`. -/
def tiedNeg : Fin 2 → ℚ := fun _ => -1

theorem bh_tiedNeg (i : Fin 2) : bh tiedNeg i = -1 := by
  have hcnt : cntLe tiedNeg (-1) = 2 := by
    simp [cntLe, tiedNeg]
  have hterm : ∀ j, bhTerm tiedNeg j = -1 := by
    intro j
    simp only [bhTerm, tiedNeg, hcnt]
    norm_num
  have hfold : ((univ.filter fun j => tiedNeg i ≤ tiedNeg j).fold min (bhTerm tiedNeg i) (bhTerm tiedNeg)) = -1 := by
    apply le_antisymm
    · exact (Finset.fold_min_le _).2 (Or.inl (le_of_eq (hterm i)))
    · exact (Finset.le_fold_min _).2 ⟨le_of_eq (hterm i).symm, fun j _ => le_of_eq (hterm j).symm⟩
  unfold bh
  rw [hfold]
  norm_num

/-- The counter-instance (ticket 5.4; statements gate re-review 1). -/
theorem bh_sorted_neg_ties_counterexample : bh_sorted_neg_ties_counterexample_stmt := by
  show (∃ i, bhSorted tiedNeg i = -2) ∧ ∀ i, bh tiedNeg i = -1
  refine ⟨?_, bh_tiedNeg⟩
  have hperm := List.mergeSort_perm (List.finRange 2) (fun a b => decide (tiedNeg b ≤ tiedNeg a))
  have hlen : ((List.finRange 2).mergeSort fun a b => decide (tiedNeg b ≤ tiedNeg a)).length = 2 := by
    rw [hperm.length_eq]
    simp
  obtain ⟨a, b, hab_eq⟩ :
      ∃ a b : Fin 2, ((List.finRange 2).mergeSort fun a b => decide (tiedNeg b ≤ tiedNeg a)) = [a, b] := by
    match hm : ((List.finRange 2).mergeSort fun a b => decide (tiedNeg b ≤ tiedNeg a)), hlen with
    | [a, b], _ => exact ⟨a, b, rfl⟩
  have hnd : [a, b].Nodup := by
    rw [← hab_eq]
    exact hperm.nodup_iff.2 (List.nodup_finRange 2)
  have hab : a ≠ b := by
    intro hab
    subst hab
    simp at hnd
  refine ⟨b, ?_⟩
  unfold bhSorted
  rw [hab_eq]
  have hidx : [a, b].idxOf b = 1 := by
    simp [List.idxOf_cons_ne _ hab]
  rw [hidx]
  simp [cumMin, tiedNeg]
  norm_num

/-- (4a) the q-value flag contains the BH flag. -/
theorem qflag_contains_bh : qflag_contains_bh_stmt := by
  intro m p π0 hπ hp
  have hle : ∀ i, qv π0 p i ≤ bh p i := by
    intro i
    have := bh_nonneg p i hp
    unfold qv
    nlinarith
  exact ⟨hle, fun i α h => ⟨qv π0 p i, rfl, lt_of_le_of_lt (hle i) h⟩⟩

/-- (4a) a missing q-value flags nothing. -/
theorem flag_none : flag_none_stmt := by
  intro α h
  obtain ⟨v, hv, _⟩ := h
  exact absurd hv (by simp)

/-- (4b) the terms the q-value flags and BH does not. -/
theorem qflag_extra : qflag_extra_stmt := by
  intro m p π0 α i hπ
  have hflag : flag (some (qv π0 p i)) α ↔ bh p i < α / π0 := by
    rw [lt_div_iff₀ hπ]
    unfold flag qv
    constructor
    · rintro ⟨v, hv, hlt⟩
      rw [← Option.some.inj hv] at hlt
      linarith
    · intro h
      exact ⟨_, rfl, by linarith⟩
  rw [hflag, not_lt]
  exact ⟨fun h => ⟨h.2, h.1⟩, fun h => ⟨h.2, h.1⟩⟩

/-- (4c) `p_i ≤ bh p i`. -/
theorem bh_ge_p : bh_ge_p_stmt := by
  intro m p i hp
  rw [le_bh_iff]
  refine ⟨(hp i).2, fun j hij => le_trans hij ?_⟩
  have hc : (0 : ℝ) < cntLe p (p j) := by exact_mod_cast cntLe_pos p j
  have hcm : (cntLe p (p j) : ℝ) ≤ m := by exact_mod_cast cntLe_le_m p (p j)
  unfold bhTerm
  rw [le_div_iff₀ hc]
  nlinarith [(hp j).1]

/-- (4c) ties. -/
theorem bh_ties : bh_ties_stmt := by
  intro m p i j h
  have key : ∀ k, bh p k = min 1 ((Finset.univ.filter fun l => p k ≤ p l).fold min
      ((m : ℝ) * p k / (cntLe p (p k) : ℝ)) (bhTerm p)) := fun k => rfl
  rw [key i, key j, h]

/-- (4c) the step-up characterisation. -/
theorem bh_stepup : bh_stepup_stmt := by
  intro m p i α hp
  exact bh_stepup_gen p i α hp

/-- (4c) monotone. -/
theorem bh_mono : bh_mono_stmt := by
  intro m p p' i hp hp' hle
  obtain ⟨k, hk, hpi, hkc⟩ := (bh_stepup_gen p' i (bh p' i) hp').1 le_rfl
  refine (bh_stepup_gen p i (bh p' i) hp).2 ⟨k, hk, le_trans (hle i) hpi, le_trans hkc ?_⟩
  exact Finset.card_le_card fun l hl =>
    Finset.mem_filter.2 ⟨Finset.mem_univ _, le_trans (hle l) (Finset.mem_filter.1 hl).2⟩

section witnesses4

theorem bh_const {K : Type} [Field K] [LinearOrder K] [IsStrictOrderedRing K] {m : ℕ} (v : K) (i : Fin m) :
    bh (fun _ : Fin m => v) i = min 1 v := by
  have hm : (0 : K) < m := by exact_mod_cast Fin.pos i
  have hcnt : cntLe (fun _ : Fin m => v) v = m := by simp [cntLe]
  have hterm : ∀ j, bhTerm (fun _ : Fin m => v) j = v := by
    intro j
    simp only [bhTerm, hcnt]
    field_simp
  apply le_antisymm
  · rw [le_min_iff]
    refine ⟨min_le_left _ _, (bh_le_iff _ i v).2 (Or.inr ⟨i, le_rfl, le_of_eq (hterm i)⟩)⟩
  · rw [le_bh_iff]
    exact ⟨min_le_left _ _, fun j _ => le_trans (min_le_right _ _) (le_of_eq (hterm j).symm)⟩

theorem qflag_contains_bh_witness : qflag_contains_bh_witness_stmt :=
  ⟨fun _ => 0, 1 / 2, by norm_num, by norm_num, fun _ => le_rfl⟩

theorem qflag_extra_witness : qflag_extra_witness_stmt := by
  refine ⟨fun _ => 1 / 2, 1 / 2, 1 / 2, 0, by norm_num, by norm_num, ?_, ?_⟩ <;>
    rw [bh_const] <;> norm_num

theorem bh_unit_witness : bh_unit_witness_stmt := by
  refine ⟨![1 / 4, 1 / 2], ![1 / 2, 3 / 4], ?_, ?_, ?_, ?_, ?_, ?_⟩ <;>
    simp [UnitP, Fin.forall_fin_two] <;> norm_num

theorem grid_witness : grid_witness_stmt := by
  unfold grid_witness_stmt
  norm_num

end witnesses4

section grid

/-- The call condition of `calls`. -/
def callCond (strict : Bool) (q v : ℚ) : Prop := if strict then v < q else v ≤ q

theorem calls_pos_iff (strict : Bool) (q : ℚ) {m : ℕ} (p : Fin m → ℚ) :
    0 < calls strict q p ↔ ∃ i, callCond strict q (bh p i) := by
  unfold calls callCond
  rw [Finset.card_pos]
  constructor
  · rintro ⟨i, hi⟩
    exact ⟨i, (Finset.mem_filter.1 hi).2⟩
  · rintro ⟨i, hi⟩
    exact ⟨i, Finset.mem_filter.2 ⟨Finset.mem_univ _, hi⟩⟩

/-- The two-term call condition: some value at most (below) `q/2`, or both at most (below) `q`. -/
def twoCond (strict : Bool) (q x y : ℚ) : Prop :=
  callCond strict q (2 * x) ∨ callCond strict q (2 * y) ∨ (callCond strict q x ∧ callCond strict q y)

theorem cnt_two (p : Fin 2 → ℚ) (v : ℚ) :
    (cntLe p v : ℚ) = (if p 0 ≤ v then 1 else 0) + (if p 1 ≤ v then 1 else 0) := by
  simp only [cntLe, Finset.card_filter, Fin.sum_univ_two]
  push_cast
  split_ifs <;> norm_num

/-- Closes one ordering case of `calls_two`: both directions of the iff, every disjunct. -/
macro "two_close" : tactic => `(tactic| (constructor; (rintro (h | h) <;> first | (left; linarith) | (right; left; linarith) | (right; right; exact ⟨by linarith, by linarith⟩)); (rintro (h | h | ⟨h, h'⟩) <;> first | (left; linarith) | (right; linarith))))

theorem calls_two (strict : Bool) (q : ℚ) (p : Fin 2 → ℚ) (h1 : 0 ≤ p 1) (hq0 : 0 < q)
    (hq : q < 1) : 0 < calls strict q p ↔ twoCond strict q (p 0) (p 1) := by
  rw [calls_pos_iff]
  have hex : (∃ i, callCond strict q (bh p i)) ↔ ∃ j, callCond strict q (bhTerm p j) := by
    unfold callCond
    cases strict
    · simp only [Bool.false_eq_true, ↓reduceIte]
      constructor
      · rintro ⟨i, hi⟩
        rcases (bh_le_iff p i q).1 hi with h | ⟨j, _, hj⟩
        · exact absurd h (not_le.2 hq)
        · exact ⟨j, hj⟩
      · rintro ⟨j, hj⟩
        exact ⟨j, (bh_le_iff p j q).2 (Or.inr ⟨j, le_rfl, hj⟩)⟩
    · simp only [↓reduceIte]
      constructor
      · rintro ⟨i, hi⟩
        rcases (bh_lt_iff p i q).1 hi with h | ⟨j, _, hj⟩
        · exact absurd h (not_lt.2 hq.le)
        · exact ⟨j, hj⟩
      · rintro ⟨j, hj⟩
        exact ⟨j, (bh_lt_iff p j q).2 (Or.inr ⟨j, le_rfl, hj⟩)⟩
  rw [hex, Fin.exists_fin_two]
  have t0 : bhTerm p 0 = 2 * p 0 / ((if p 1 ≤ p 0 then 1 else 0) + 1) := by
    simp only [bhTerm, cnt_two, le_refl, ↓reduceIte]
    push_cast
    ring_nf
  have t1 : bhTerm p 1 = 2 * p 1 / ((if p 0 ≤ p 1 then 1 else 0) + 1) := by
    simp only [bhTerm, cnt_two, le_refl, ↓reduceIte]
    push_cast
    ring_nf
  rw [t0, t1]
  unfold twoCond callCond
  have d2 : ∀ x : ℚ, 2 * x / ((1 : ℚ) + 1) = x := fun x => by ring
  have d1 : ∀ x : ℚ, 2 * x / ((0 : ℚ) + 1) = 2 * x := fun x => by ring
  rcases le_total (p 0) (p 1) with hxy | hyx
  · by_cases hyx : p 1 ≤ p 0
    · have he : p 0 = p 1 := le_antisymm hxy hyx
      simp only [hxy, hyx, ↓reduceIte, d2]
      rw [he]
      cases strict <;> simp only [Bool.false_eq_true, ↓reduceIte] <;> two_close
    · have hlt : p 0 < p 1 := lt_of_not_ge hyx
      simp only [hxy, hyx, ↓reduceIte, d1, d2]
      cases strict <;> simp only [Bool.false_eq_true, ↓reduceIte] <;> two_close
  · by_cases hxy : p 0 ≤ p 1
    · have he : p 0 = p 1 := le_antisymm hxy hyx
      simp only [hxy, hyx, ↓reduceIte, d2]
      rw [he]
      cases strict <;> simp only [Bool.false_eq_true, ↓reduceIte] <;> two_close
    · have hlt : p 1 < p 0 := lt_of_not_ge hxy
      simp only [hxy, hyx, ↓reduceIte, d1, d2]
      cases strict <;> simp only [Bool.false_eq_true, ↓reduceIte] <;> two_close

/-- BH's calls on a list (the form of `callsKept` and `callsWhole`). -/
def callsL (strict : Bool) (q : ℚ) (l : List ℚ) : ℕ := calls strict q fun i : Fin l.length => l.get i

theorem callsL_nil (strict : Bool) (q : ℚ) : callsL strict q [] = 0 := by
  simp [callsL, calls]

theorem callsL_one (strict : Bool) (q : ℚ) (x : ℚ) (hq : q < 1) :
    0 < callsL strict q [x] ↔ callCond strict q x := by
  unfold callsL
  rw [calls_pos_iff]
  have hp : (fun i : Fin [x].length => [x].get i) = fun _ : Fin 1 => x := by
    funext i
    fin_cases i
    rfl
  rw [hp]
  change (∃ i : Fin 1, callCond strict q (bh (fun _ : Fin 1 => x) i)) ↔ _
  simp only [bh_const, Fin.exists_fin_one]
  unfold callCond
  cases strict
  · simp only [Bool.false_eq_true, ↓reduceIte, min_le_iff, not_le.2 hq, false_or]
  · simp only [↓reduceIte, min_lt_iff, not_lt.2 hq.le, false_or]

theorem callsL_two (strict : Bool) (q : ℚ) (x y : ℚ) (hy : 0 ≤ y) (hq0 : 0 < q) (hq : q < 1) :
    0 < callsL strict q [x, y] ↔ twoCond strict q x y :=
  calls_two strict q (fun i : Fin 2 => [x, y].get i) hy hq0 hq

theorem gridVal_zero (N : ℕ) : gridVal N 0 = 1 := by simp [gridVal]

theorem gridVal_succ (N : ℕ) (b : Fin N) : gridVal N b.succ = ((b : ℕ) + 1 : ℚ) / (2 * N) := by
  simp [gridVal, Fin.val_succ]

theorem gridVal_succ_nonneg (N : ℕ) (b : Fin N) : 0 ≤ gridVal N b.succ := by
  rw [gridVal_succ]
  positivity

theorem gridVal_succ_lt (N : ℕ) (b : Fin N) : gridVal N b.succ < 1 := by
  rw [gridVal_succ]
  have hN : (0 : ℚ) < N := by exact_mod_cast (Fin.pos b)
  rw [div_lt_one (by positivity)]
  have : ((b : ℕ) : ℚ) + 1 ≤ N := by exact_mod_cast b.isLt
  linarith

theorem gridWt_zero (N : ℕ) : gridWt N 0 = 1 / 2 := by simp [gridWt]

theorem gridWt_succ (N : ℕ) (b : Fin N) : gridWt N b.succ = 1 / (2 * N) := by simp [gridWt, Fin.val_succ]

/-- The outcome's call indicator. -/
def ind (N : ℕ) (strict : Bool) (q : ℚ) (keptOnly : Bool) (a b : Fin (N + 1)) : ℚ :=
  if 0 < (if keptOnly then callsKept strict q [gridVal N a, gridVal N b]
    else callsWhole strict q [gridVal N a, gridVal N b]) then 1 else 0

theorem fdr2_split (N : ℕ) (strict : Bool) (q : ℚ) (keptOnly : Bool) :
    fdr2 N strict q keptOnly = 1 / 2 * (1 / 2) * ind N strict q keptOnly 0 0 +
      ∑ b : Fin N, 1 / 2 * (1 / (2 * N)) * ind N strict q keptOnly 0 b.succ +
      ∑ a : Fin N, (1 / (2 * N) * (1 / 2) * ind N strict q keptOnly a.succ 0 +
        ∑ b : Fin N, 1 / (2 * N) * (1 / (2 * N)) * ind N strict q keptOnly a.succ b.succ) := by
  unfold fdr2
  rw [Fintype.sum_prod_type, Fin.sum_univ_succ]
  simp only [Fin.sum_univ_succ, gridWt_zero, gridWt_succ]
  rfl

instance (strict : Bool) (q v : ℚ) : Decidable (callCond strict q v) :=
  inferInstanceAs (Decidable (if strict then v < q else v ≤ q))

instance (strict : Bool) (q x y : ℚ) : Decidable (twoCond strict q x y) :=
  inferInstanceAs (Decidable (callCond strict q (2 * x) ∨ callCond strict q (2 * y) ∨
    (callCond strict q x ∧ callCond strict q y)))

theorem ind_kept_00 (N : ℕ) (strict : Bool) (q : ℚ) : ind N strict q true 0 0 = 0 := by
  have h : callsKept strict q [gridVal N 0, gridVal N 0] = callsL strict q [] := by
    show callsL strict q _ = _
    simp [gridVal_zero]
  simp [ind, h, callsL_nil]

theorem ind_kept_0s (N : ℕ) (strict : Bool) (q : ℚ) (hq : q < 1) (b : Fin N) :
    ind N strict q true 0 b.succ = if callCond strict q (gridVal N b.succ) then 1 else 0 := by
  have h : callsKept strict q [gridVal N 0, gridVal N b.succ] = callsL strict q [gridVal N b.succ] := by
    show callsL strict q _ = _
    simp [gridVal_zero, gridVal_succ_lt]
  simp only [ind, ↓reduceIte, h, callsL_one strict q _ hq]

theorem ind_kept_s0 (N : ℕ) (strict : Bool) (q : ℚ) (hq : q < 1) (a : Fin N) :
    ind N strict q true a.succ 0 = if callCond strict q (gridVal N a.succ) then 1 else 0 := by
  have h : callsKept strict q [gridVal N a.succ, gridVal N 0] = callsL strict q [gridVal N a.succ] := by
    show callsL strict q _ = _
    simp [gridVal_zero, gridVal_succ_lt]
  simp only [ind, ↓reduceIte, h, callsL_one strict q _ hq]

theorem ind_kept_ss (N : ℕ) (strict : Bool) (q : ℚ) (hq0 : 0 < q) (hq : q < 1) (a b : Fin N) :
    ind N strict q true a.succ b.succ =
      if twoCond strict q (gridVal N a.succ) (gridVal N b.succ) then 1 else 0 := by
  have h : callsKept strict q [gridVal N a.succ, gridVal N b.succ] =
      callsL strict q [gridVal N a.succ, gridVal N b.succ] := by
    show callsL strict q _ = _
    simp [gridVal_succ_lt]
  simp only [ind, ↓reduceIte, h, callsL_two strict q _ _ (gridVal_succ_nonneg N b) hq0 hq]

theorem twoCond_one_left (strict : Bool) (q y : ℚ) (hq : q < 1) :
    twoCond strict q 1 y ↔ callCond strict q (2 * y) := by
  unfold twoCond callCond
  cases strict <;> simp only [Bool.false_eq_true, ↓reduceIte] <;> constructor <;>
    (try rintro (h | h | ⟨h, h'⟩)) <;> (try intro h) <;> first | linarith | (right; left; linarith)

theorem twoCond_one_right (strict : Bool) (q x : ℚ) (hq : q < 1) :
    twoCond strict q x 1 ↔ callCond strict q (2 * x) := by
  unfold twoCond callCond
  cases strict <;> simp only [Bool.false_eq_true, ↓reduceIte] <;> constructor <;>
    (try rintro (h | h | ⟨h, h'⟩)) <;> (try intro h) <;> first | linarith | (left; linarith)

theorem ind_whole_00 (N : ℕ) (strict : Bool) (q : ℚ) (hq0 : 0 < q) (hq : q < 1) :
    ind N strict q false 0 0 = 0 := by
  have h := callsL_two strict q 1 1 zero_le_one hq0 hq
  have hn : ¬ twoCond strict q 1 1 := by
    rw [twoCond_one_left strict q 1 hq]
    unfold callCond
    cases strict <;> simp only [Bool.false_eq_true, ↓reduceIte] <;> intro h' <;> linarith
  simp only [ind, Bool.false_eq_true, ↓reduceIte, gridVal_zero]
  show (if 0 < callsL strict q [1, 1] then (1 : ℚ) else 0) = 0
  have hf : ¬ 0 < callsL strict q [1, 1] := fun h' => hn (h.1 h')
  simp only [hf, ↓reduceIte]

theorem ind_whole_0s (N : ℕ) (strict : Bool) (q : ℚ) (hq0 : 0 < q) (hq : q < 1) (b : Fin N) :
    ind N strict q false 0 b.succ = if callCond strict q (2 * gridVal N b.succ) then 1 else 0 := by
  simp only [ind, Bool.false_eq_true, ↓reduceIte, gridVal_zero]
  show (if 0 < callsL strict q [1, gridVal N b.succ] then (1 : ℚ) else 0) = _
  simp only [callsL_two strict q 1 _ (gridVal_succ_nonneg N b) hq0 hq,
    twoCond_one_left strict q _ hq]

theorem ind_whole_s0 (N : ℕ) (strict : Bool) (q : ℚ) (hq0 : 0 < q) (hq : q < 1) (a : Fin N) :
    ind N strict q false a.succ 0 = if callCond strict q (2 * gridVal N a.succ) then 1 else 0 := by
  simp only [ind, Bool.false_eq_true, ↓reduceIte, gridVal_zero]
  show (if 0 < callsL strict q [gridVal N a.succ, 1] then (1 : ℚ) else 0) = _
  simp only [callsL_two strict q _ 1 zero_le_one hq0 hq,
    twoCond_one_right strict q _ hq]

theorem ind_whole_ss (N : ℕ) (strict : Bool) (q : ℚ) (hq0 : 0 < q) (hq : q < 1) (a b : Fin N) :
    ind N strict q false a.succ b.succ =
      if twoCond strict q (gridVal N a.succ) (gridVal N b.succ) then 1 else 0 := by
  simp only [ind, Bool.false_eq_true, ↓reduceIte]
  show (if 0 < callsL strict q [gridVal N a.succ, gridVal N b.succ] then (1 : ℚ) else 0) = _
  simp only [callsL_two strict q _ _ (gridVal_succ_nonneg N b) hq0 hq]

theorem natle_iff (n T : ℕ) (r : ℚ) (h1 : (T : ℚ) ≤ r) (h2 : r < T + 1) : (n : ℚ) ≤ r ↔ n ≤ T := by
  constructor
  · intro h
    have : (n : ℚ) < T + 1 := lt_of_le_of_lt h h2
    have : n < T + 1 := by exact_mod_cast this
    omega
  · intro h
    exact le_trans (by exact_mod_cast h) h1

theorem natlt_iff (n T : ℕ) (r : ℚ) (h1 : (T : ℚ) < r) (h2 : r ≤ T + 1) : (n : ℚ) < r ↔ n ≤ T := by
  constructor
  · intro h
    have : (n : ℚ) < T + 1 := lt_of_lt_of_le h h2
    have : n < T + 1 := by exact_mod_cast this
    omega
  · intro h
    exact lt_of_le_of_lt (by exact_mod_cast h) h1

/-- A grid value's call condition as a whole-number condition on its index. -/
theorem cond_grid (N : ℕ) (strict : Bool) (q c : ℚ) (hc : 0 < c) (b : Fin N) (T : ℕ)
    (hT : if strict then (T : ℚ) < q * (2 * N) / c ∧ q * (2 * N) / c ≤ T + 1
      else (T : ℚ) ≤ q * (2 * N) / c ∧ q * (2 * N) / c < T + 1) :
    callCond strict q (c * gridVal N b.succ) ↔ (b : ℕ) + 1 ≤ T := by
  have hN : (0 : ℚ) < 2 * N := by
    have : (0 : ℚ) < N := by exact_mod_cast Fin.pos b
    linarith
  rw [gridVal_succ]
  have hx : ((b : ℕ) : ℚ) + 1 = (((b : ℕ) + 1 : ℕ) : ℚ) := by push_cast; ring
  rw [hx]
  unfold callCond
  cases strict
  · simp only [Bool.false_eq_true, ↓reduceIte] at hT ⊢
    rw [← natle_iff _ T _ hT.1 hT.2, show c * ((((b : ℕ) + 1 : ℕ) : ℚ) / (2 * N)) =
      (((b : ℕ) + 1 : ℕ) : ℚ) * c / (2 * N) by ring, div_le_iff₀ hN, le_div_iff₀ hc]
  · simp only [↓reduceIte] at hT ⊢
    rw [← natlt_iff _ T _ hT.1 hT.2, show c * ((((b : ℕ) + 1 : ℕ) : ℚ) / (2 * N)) =
      (((b : ℕ) + 1 : ℕ) : ℚ) * c / (2 * N) by ring, div_lt_iff₀ hN, lt_div_iff₀ hc]

theorem cond_grid1 (N : ℕ) (strict : Bool) (q : ℚ) (b : Fin N) (T : ℕ)
    (hT : if strict then (T : ℚ) < q * (2 * N) / 1 ∧ q * (2 * N) / 1 ≤ T + 1
      else (T : ℚ) ≤ q * (2 * N) / 1 ∧ q * (2 * N) / 1 < T + 1) :
    callCond strict q (gridVal N b.succ) ↔ (b : ℕ) + 1 ≤ T := by
  rw [← cond_grid N strict q 1 one_pos b T hT, one_mul]

end grid



theorem fdr_n20_kept : fdr2 20 false (1 / 20) true = 3 / 40 := by
  have hq0 : (0 : ℚ) < 1 / 20 := by norm_num
  have hq : (1 / 20 : ℚ) < 1 := by norm_num
  have c1 : ∀ b : Fin 20, callCond false (1 / 20) (gridVal 20 b.succ) ↔ (b : ℕ) + 1 ≤ 2 :=
    fun b => cond_grid1 20 false _ b 2 (by norm_num)
  have c2 : ∀ b : Fin 20, callCond false (1 / 20) (2 * gridVal 20 b.succ) ↔ (b : ℕ) + 1 ≤ 1 :=
    fun b => cond_grid 20 false _ 2 (by norm_num) b 1 (by norm_num)
  rw [fdr2_split]
  simp only [ind_kept_00, ind_kept_0s _ _ _ hq, ind_kept_s0 _ _ _ hq, ind_kept_ss _ _ _ hq0 hq, twoCond, c1, c2]
  simp [Fin.sum_univ_succ]
  norm_num

theorem fdr_n20_whole : fdr2 20 false (1 / 20) false = 1 / 20 := by
  have hq0 : (0 : ℚ) < 1 / 20 := by norm_num
  have hq : (1 / 20 : ℚ) < 1 := by norm_num
  have c1 : ∀ b : Fin 20, callCond false (1 / 20) (gridVal 20 b.succ) ↔ (b : ℕ) + 1 ≤ 2 :=
    fun b => cond_grid1 20 false _ b 2 (by norm_num)
  have c2 : ∀ b : Fin 20, callCond false (1 / 20) (2 * gridVal 20 b.succ) ↔ (b : ℕ) + 1 ≤ 1 :=
    fun b => cond_grid 20 false _ 2 (by norm_num) b 1 (by norm_num)
  rw [fdr2_split]
  simp only [ind_whole_00 _ _ _ hq0 hq, ind_whole_0s _ _ _ hq0 hq, ind_whole_s0 _ _ _ hq0 hq, ind_whole_ss _ _ _ hq0 hq, twoCond, c1, c2]
  simp [Fin.sum_univ_succ]
  norm_num

theorem fdr_n20_kept_strict : fdr2 20 true (1 / 20) true = 41 / 1600 := by
  have hq0 : (0 : ℚ) < 1 / 20 := by norm_num
  have hq : (1 / 20 : ℚ) < 1 := by norm_num
  have c1 : ∀ b : Fin 20, callCond true (1 / 20) (gridVal 20 b.succ) ↔ (b : ℕ) + 1 ≤ 1 :=
    fun b => cond_grid1 20 true _ b 1 (by norm_num)
  have c2 : ∀ b : Fin 20, callCond true (1 / 20) (2 * gridVal 20 b.succ) ↔ (b : ℕ) + 1 ≤ 0 :=
    fun b => cond_grid 20 true _ 2 (by norm_num) b 0 (by norm_num)
  rw [fdr2_split]
  simp only [ind_kept_00, ind_kept_0s _ _ _ hq, ind_kept_s0 _ _ _ hq, ind_kept_ss _ _ _ hq0 hq, twoCond, c1, c2]
  simp [Fin.sum_univ_succ]
  norm_num

theorem fdr_n25_kept : fdr2 25 false (1 / 20) true = 3 / 50 := by
  have hq0 : (0 : ℚ) < 1 / 20 := by norm_num
  have hq : (1 / 20 : ℚ) < 1 := by norm_num
  have c1 : ∀ b : Fin 25, callCond false (1 / 20) (gridVal 25 b.succ) ↔ (b : ℕ) + 1 ≤ 2 :=
    fun b => cond_grid1 25 false _ b 2 (by norm_num)
  have c2 : ∀ b : Fin 25, callCond false (1 / 20) (2 * gridVal 25 b.succ) ↔ (b : ℕ) + 1 ≤ 1 :=
    fun b => cond_grid 25 false _ 2 (by norm_num) b 1 (by norm_num)
  rw [fdr2_split]
  simp only [ind_kept_00, ind_kept_0s _ _ _ hq, ind_kept_s0 _ _ _ hq, ind_kept_ss _ _ _ hq0 hq, twoCond, c1, c2]
  simp [Fin.sum_univ_succ]
  norm_num

theorem fdr_n25_whole : fdr2 25 false (1 / 20) false = 1 / 25 := by
  have hq0 : (0 : ℚ) < 1 / 20 := by norm_num
  have hq : (1 / 20 : ℚ) < 1 := by norm_num
  have c1 : ∀ b : Fin 25, callCond false (1 / 20) (gridVal 25 b.succ) ↔ (b : ℕ) + 1 ≤ 2 :=
    fun b => cond_grid1 25 false _ b 2 (by norm_num)
  have c2 : ∀ b : Fin 25, callCond false (1 / 20) (2 * gridVal 25 b.succ) ↔ (b : ℕ) + 1 ≤ 1 :=
    fun b => cond_grid 25 false _ 2 (by norm_num) b 1 (by norm_num)
  rw [fdr2_split]
  simp only [ind_whole_00 _ _ _ hq0 hq, ind_whole_0s _ _ _ hq0 hq, ind_whole_s0 _ _ _ hq0 hq, ind_whole_ss _ _ _ hq0 hq, twoCond, c1, c2]
  simp [Fin.sum_univ_succ]
  norm_num

theorem fdr_n25_kept_strict : fdr2 25 true (1 / 20) true = 3 / 50 := by
  have hq0 : (0 : ℚ) < 1 / 20 := by norm_num
  have hq : (1 / 20 : ℚ) < 1 := by norm_num
  have c1 : ∀ b : Fin 25, callCond true (1 / 20) (gridVal 25 b.succ) ↔ (b : ℕ) + 1 ≤ 2 :=
    fun b => cond_grid1 25 true _ b 2 (by norm_num)
  have c2 : ∀ b : Fin 25, callCond true (1 / 20) (2 * gridVal 25 b.succ) ↔ (b : ℕ) + 1 ≤ 1 :=
    fun b => cond_grid 25 true _ 2 (by norm_num) b 1 (by norm_num)
  rw [fdr2_split]
  simp only [ind_kept_00, ind_kept_0s _ _ _ hq, ind_kept_s0 _ _ _ hq, ind_kept_ss _ _ _ hq0 hq, twoCond, c1, c2]
  simp [Fin.sum_univ_succ]
  norm_num

theorem fdr_n25_whole_strict : fdr2 25 true (1 / 20) false = 1 / 25 := by
  have hq0 : (0 : ℚ) < 1 / 20 := by norm_num
  have hq : (1 / 20 : ℚ) < 1 := by norm_num
  have c1 : ∀ b : Fin 25, callCond true (1 / 20) (gridVal 25 b.succ) ↔ (b : ℕ) + 1 ≤ 2 :=
    fun b => cond_grid1 25 true _ b 2 (by norm_num)
  have c2 : ∀ b : Fin 25, callCond true (1 / 20) (2 * gridVal 25 b.succ) ↔ (b : ℕ) + 1 ≤ 1 :=
    fun b => cond_grid 25 true _ 2 (by norm_num) b 1 (by norm_num)
  rw [fdr2_split]
  simp only [ind_whole_00 _ _ _ hq0 hq, ind_whole_0s _ _ _ hq0 hq, ind_whole_s0 _ _ _ hq0 hq, ind_whole_ss _ _ _ hq0 hq, twoCond, c1, c2]
  simp [Fin.sum_univ_succ]
  norm_num

/-- (4d-i) `N = 20`, standard rule. -/
theorem fdr_n20 : fdr_n20_stmt := ⟨fdr_n20_kept, fdr_n20_whole⟩

/-- (4d-i) `N = 20`, strict rule over the kept family. -/
theorem fdr_n20_strict : fdr_n20_strict_stmt := fdr_n20_kept_strict

/-- (4d-i) `N = 25`, either comparison. -/
theorem fdr_n25 : fdr_n25_stmt := by
  intro strict
  cases strict
  · exact ⟨fdr_n25_kept, fdr_n25_whole⟩
  · exact ⟨fdr_n25_kept_strict, fdr_n25_whole_strict⟩

/-- The grid law's weights sum to 1. -/
theorem grid_mass : grid_mass_stmt := by
  intro N hN
  have hN' : (0 : ℚ) < N := by exact_mod_cast hN
  rw [Fin.sum_univ_succ]
  simp only [gridWt_zero, gridWt_succ, Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
  field_simp
  ring

/-- Each grid p-value is super-uniform. -/
theorem grid_superuniform : grid_superuniform_stmt := by
  intro N t hN ht0 ht1
  have hN' : (0 : ℚ) < N := by exact_mod_cast hN
  rcases lt_or_eq_of_le ht1 with htl | hteq
  · rw [Finset.sum_filter, Fin.sum_univ_succ]
    have h0 : ¬ gridVal N 0 ≤ t := by rw [gridVal_zero]; exact not_le.2 htl
    simp only [h0, ↓reduceIte, zero_add, gridWt_succ]
    rw [Finset.sum_ite, Finset.sum_const_zero, add_zero, Finset.sum_const, nsmul_eq_mul]
    set S := Finset.univ.filter fun b : Fin N => gridVal N b.succ ≤ t with hS
    have hcard : S.card ≤ ⌊2 * N * t⌋₊ := by
      have : S.card ≤ (Finset.Icc 1 ⌊2 * N * t⌋₊).card := by
        apply Finset.card_le_card_of_injOn (fun b : Fin N => (b : ℕ) + 1)
        · intro b hb
          simp only [hS, Finset.coe_filter, Finset.mem_univ, true_and, Set.mem_ofPred_eq] at hb
          rw [gridVal_succ, div_le_iff₀ (by positivity)] at hb
          simp only [Finset.coe_Icc, Set.mem_Icc]
          refine ⟨by omega, Nat.le_floor ?_⟩
          push_cast
          linarith
        · intro a _ b _ hab
          exact Fin.ext (by simpa using hab)
      simpa using this
    have hfl : (⌊2 * N * t⌋₊ : ℚ) ≤ 2 * N * t := Nat.floor_le (by positivity)
    have hc : (S.card : ℚ) ≤ 2 * N * t := le_trans (by exact_mod_cast hcard) hfl
    rw [← div_eq_mul_one_div, div_le_iff₀ (by positivity)]
    linarith
  · subst hteq
    calc ∑ j ∈ Finset.univ.filter (fun j => gridVal N j ≤ 1), gridWt N j ≤ ∑ j : Fin (N + 1), gridWt N j :=
          Finset.sum_le_sum_of_subset_of_nonneg (Finset.filter_subset _ _) (fun j _ _ => by
            unfold gridWt; split_ifs <;> positivity)
      _ = 1 := grid_mass N hN

end MultiScholaRProofs.QvalueBh
