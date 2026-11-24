module

public import MultiScholaRProofs.Spec.QvalueBh
public import MultiScholaRProofs.Proof.QvalueBh
import Mathlib.Tactic

/-!
# Target 4d-ii (KIT-193 section 6.4): the kept family on a fine grid

Theorems whose types are the frozen statements of `Spec/QvalueBh.lean`: `kept_family_unit_a_counterexample` (the
formula fails at `a = 1`), `kept_family_formula`, `kept_family_instance`, `kept_family_witness`.

How. (1) `kept_event_iff`: BH over the terms with `p < 1` calls at least one term exactly when some `r ≥ 1` has at
least `r` kept p-values at most `r q / K` (`K` the number kept), by the step-up characterisation `bh_stepup_gen`.
(2) `simes_grid`, Simes' equality on a grid: for `n` values uniform on `1..G` and thresholds `λ r` that are whole
numbers wherever they are at most `G`, the number of outcomes with some `r` such that at least `r` values are at
most `λ r` is `min(λ n, G) G^(n-1)`. By strong induction on `n`: with `g = λ n < G`, split the outcomes by the set
`S` of values above `g` (frozen out by `freeze_set`); the rest is the same problem on `T \ S` with `G = g`.
(3) The false discovery rate splits by the set of kept terms; each set of size `k ≥ 1` contributes
`a^k (1 - a)^(m - k) min(q/a, 1)`, and the binomial theorem sums these to `min(q/a, 1) (1 - (1 - a)^m)`.
-/

public section

namespace MultiScholaRProofs.QvalueBh

open Finset

/-! ## Counting on lists -/

theorem card_filter_get (l : List ℚ) (P : ℚ → Prop) [DecidablePred P] :
    (univ.filter fun i : Fin l.length => P (l.get i)).card = l.countP (fun x => decide (P x)) := by
  induction l with
  | nil => simp
  | cons x l ih =>
    rw [Finset.card_filter, List.countP_cons, ← ih, Finset.card_filter]
    change (∑ i : Fin (l.length + 1), (if P ((x :: l).get i) then 1 else 0)) = _
    rw [Fin.sum_univ_succ]
    simp only [List.get_eq_getElem, Fin.val_zero, List.getElem_cons_zero, Fin.val_succ,
      List.getElem_cons_succ]
    by_cases h : P x <;> simp [h, add_comm]

theorem countP_ofFn {m : ℕ} (v : Fin m → ℚ) (P : ℚ → Bool) :
    (List.ofFn v).countP P = (univ.filter fun b => P (v b) = true).card := by
  induction m with
  | zero => simp
  | succ m ih =>
    rw [List.ofFn_succ, List.countP_cons, ih, Finset.card_filter, Finset.card_filter, Fin.sum_univ_succ]
    by_cases h : P (v 0) <;> simp [h, add_comm]

/-- The kept list: the values below 1, in order. -/
abbrev keptL {m : ℕ} (v : Fin m → ℚ) : List ℚ := (List.ofFn v).filter fun x => decide (x < 1)

theorem keptL_length {m : ℕ} (v : Fin m → ℚ) : (keptL v).length = (univ.filter fun b => v b < 1).card := by
  rw [← List.countP_eq_length_filter, countP_ofFn]
  simp

theorem keptL_cnt {m : ℕ} (v : Fin m → ℚ) (t : ℚ) :
    cntLe (fun i : Fin (keptL v).length => (keptL v).get i) t = (univ.filter fun b => v b < 1 ∧ v b ≤ t).card := by
  unfold cntLe
  change (univ.filter fun i : Fin (keptL v).length => (keptL v).get i ≤ t).card = _
  rw [card_filter_get (keptL v) (fun x => x ≤ t), List.countP_filter, countP_ofFn]
  congr 1
  ext b
  simp [and_comm]

theorem keptL_get {m : ℕ} (v : Fin m → ℚ) (i : Fin (keptL v).length) : ∃ b, v b = (keptL v).get i ∧ v b < 1 := by
  have hmem : (keptL v).get i ∈ keptL v := List.get_mem _ _
  rw [List.mem_filter, List.mem_ofFn] at hmem
  obtain ⟨⟨b, hb⟩, hlt⟩ := hmem
  exact ⟨b, hb, by rw [hb]; simpa using hlt⟩

/-- BH over the terms with `p < 1` calls at least one term exactly when some `r ≥ 1` has at least `r` kept values
at most `r q / K`, `K` the number kept. -/
theorem kept_event_iff {m : ℕ} (v : Fin m → ℚ) (hv : ∀ b, 0 ≤ v b) (q : ℚ) :
    0 < callsKept false q (List.ofFn v) ↔
      ∃ r : ℕ, 1 ≤ r ∧
        r ≤ (univ.filter fun b => v b < 1 ∧ v b ≤ r * q / ((univ.filter fun b => v b < 1).card : ℚ)).card := by
  set p : Fin (keptL v).length → ℚ := fun i => (keptL v).get i with hp
  have hU : UnitP p := by
    intro i
    obtain ⟨b, hb, hlt⟩ := keptL_get v i
    rw [hp]
    dsimp only
    rw [← hb]
    exact ⟨hv b, hlt.le⟩
  have hcalls : callsKept false q (List.ofFn v) = (univ.filter fun i => bh p i ≤ q).card := by
    unfold callsKept calls
    simp only [Bool.false_eq_true, ↓reduceIte]
    rfl
  rw [hcalls, Finset.card_pos]
  constructor
  · rintro ⟨i, hi⟩
    rw [Finset.mem_filter] at hi
    obtain ⟨k, hk, -, hkc⟩ := (bh_stepup_gen p i q hU).1 hi.2
    refine ⟨k, hk, ?_⟩
    rw [hp, keptL_cnt, keptL_length] at hkc
    exact hkc
  · rintro ⟨k, hk, hkc⟩
    have hkc' : k ≤ cntLe p (k * q / ((keptL v).length : ℚ)) := by
      rw [hp, keptL_cnt, keptL_length]
      exact hkc
    have hpos : 0 < cntLe p (k * q / ((keptL v).length : ℚ)) := lt_of_lt_of_le hk hkc'
    unfold cntLe at hpos
    obtain ⟨i, hi⟩ := Finset.card_pos.1 hpos
    refine ⟨i, Finset.mem_filter.2 ⟨Finset.mem_univ _, ?_⟩⟩
    exact (bh_stepup_gen p i q hU).2 ⟨k, hk, (Finset.mem_filter.1 hi).2, hkc'⟩

/-! ## Freezing coordinates a predicate does not read -/

section Freeze

variable {ι : Type} [Fintype ι] [DecidableEq ι]

theorem freeze_one (A : ι → Finset ℕ) (j : ι) (P : (ι → ℕ) → Prop) [DecidablePred P]
    (hP : ∀ u x, P (Function.update u j x) ↔ P u) :
    ((Fintype.piFinset A).filter P).card =
      (A j).card * ((Fintype.piFinset (Function.update A j {0})).filter P).card := by
  rw [← Finset.card_product]
  apply Finset.card_nbij' (fun u => (u j, Function.update u j 0)) (fun xv => Function.update xv.2 j xv.1)
  · intro u hu
    simp only [Finset.coe_filter, Set.mem_ofPred_eq, Fintype.mem_piFinset] at hu
    simp only [Finset.coe_product, Set.mem_prod, Finset.mem_coe, Finset.mem_filter, Fintype.mem_piFinset]
    refine ⟨hu.1 j, fun i => ?_, (hP u 0).2 hu.2⟩
    by_cases h : i = j
    · subst h
      simp
    · simp only [Function.update_of_ne h]
      exact hu.1 i
  · intro xv hxv
    simp only [Finset.coe_product, Set.mem_prod, Finset.mem_coe, Finset.mem_filter, Fintype.mem_piFinset] at hxv
    simp only [Finset.coe_filter, Set.mem_ofPred_eq, Fintype.mem_piFinset]
    refine ⟨fun i => ?_, (hP _ _).2 hxv.2.2⟩
    by_cases h : i = j
    · subst h
      simpa using hxv.1
    · have := hxv.2.1 i
      simp only [Function.update_of_ne h] at this ⊢
      exact this
  · intro u _
    simp
  · intro xv hxv
    simp only [Finset.coe_product, Set.mem_prod, Finset.mem_coe, Finset.mem_filter, Fintype.mem_piFinset] at hxv
    have h0 : xv.2 j = 0 := by simpa using hxv.2.1 j
    refine Prod.ext (by simp) ?_
    funext i
    by_cases h : i = j
    · subst h
      simp [h0]
    · simp [Function.update_of_ne h]

theorem freeze_set (P : (ι → ℕ) → Prop) [DecidablePred P] (S : Finset ι) :
    ∀ A : ι → Finset ℕ, (∀ u v : ι → ℕ, (∀ i, i ∉ S → u i = v i) → (P u ↔ P v)) →
      ((Fintype.piFinset A).filter P).card =
        (∏ i ∈ S, (A i).card) * ((Fintype.piFinset fun i => if i ∈ S then {0} else A i).filter P).card := by
  induction S using Finset.induction_on with
  | empty =>
    intro A _
    simp
  | insert j S hj ih =>
    intro A hP
    have hPj : ∀ u x, P (Function.update u j x) ↔ P u := fun u x =>
      hP _ _ fun i hi => by
        have hij : i ≠ j := fun h => hi (h ▸ Finset.mem_insert_self _ _)
        simp [Function.update_of_ne hij]
    rw [freeze_one A j P hPj,
      ih (Function.update A j {0}) fun u v h => hP u v fun i hi => h i fun hS => hi (Finset.mem_insert_of_mem hS),
      Finset.prod_insert hj]
    have h1 : ∏ i ∈ S, (Function.update A j {0} i).card = ∏ i ∈ S, (A i).card :=
      Finset.prod_congr rfl fun i hi => by
        have hij : i ≠ j := fun h => hj (h ▸ hi)
        rw [Function.update_of_ne hij]
    have h2 : (fun i => if i ∈ S then ({0} : Finset ℕ) else Function.update A j {0} i) =
        fun i => if i ∈ insert j S then {0} else A i := by
      funext i
      by_cases hi : i = j
      · subst hi
        simp [hj]
      · simp [Finset.mem_insert, hi]
    rw [h1, h2, mul_assoc]

end Freeze

/-! ## Simes' equality on a grid -/

section Simes

variable {ι : Type} [Fintype ι] [DecidableEq ι]

/-- The outcomes: coordinates in `T` range over `1..G`, the others are `0`. -/
def dom (T : Finset ι) (G : ℕ) : Finset (ι → ℕ) := Fintype.piFinset fun i => if i ∈ T then Icc 1 G else {0}

/-- The number of coordinates of `T` whose value is at most `t`. -/
def cntT (T : Finset ι) (u : ι → ℕ) (t : ℚ) : ℕ := (T.filter fun i => (u i : ℚ) ≤ t).card

/-- Simes' event on `T` with thresholds `λ r`. -/
abbrev simesEv (T : Finset ι) (lam : ℚ) (u : ι → ℕ) : Prop := ∃ r ∈ Icc 1 T.card, r ≤ cntT T u (lam * r)

theorem mem_dom (T : Finset ι) (G : ℕ) (u : ι → ℕ) :
    u ∈ dom T G ↔ ∀ i, (i ∈ T → 1 ≤ u i ∧ u i ≤ G) ∧ (i ∉ T → u i = 0) := by
  unfold dom
  rw [Fintype.mem_piFinset]
  refine forall_congr' fun i => ?_
  by_cases h : i ∈ T <;> simp [h]

theorem card_dom (T : Finset ι) (G : ℕ) : (dom T G).card = G ^ T.card := by
  unfold dom
  rw [Fintype.card_piFinset]
  rw [show (∏ i, (if i ∈ T then Icc 1 G else {0} : Finset ℕ).card) = ∏ i, (if i ∈ T then G else 1) from
    Finset.prod_congr rfl fun i _ => by split_ifs <;> simp]
  rw [Finset.prod_ite_mem, Finset.univ_inter, Finset.prod_const]

omit [Fintype ι] [DecidableEq ι] in
theorem cntT_le (T : Finset ι) (u : ι → ℕ) (t : ℚ) : cntT T u t ≤ T.card := Finset.card_filter_le _ _

/-- When every value is at most the top threshold, every outcome is in the event. -/
theorem simes_full (T : Finset ι) (G : ℕ) (lam : ℚ) (hn : 1 ≤ T.card) (hG : (G : ℚ) ≤ lam * T.card) :
    ((dom T G).filter (simesEv T lam)).card = G ^ T.card := by
  rw [Finset.filter_true_of_mem, card_dom]
  intro u hu
  refine ⟨T.card, Finset.mem_Icc.2 ⟨hn, le_rfl⟩, le_of_eq ?_⟩
  unfold cntT
  rw [Finset.filter_true_of_mem]
  intro i hi
  have := ((mem_dom T G u).1 hu i).1 hi
  exact le_trans (by exact_mod_cast this.2) hG

/-- `Σ_k C(n,k) x^k (n - k) y^(n-k-1) = n (x + y)^(n-1)`. -/
theorem binom_deriv (n : ℕ) (x y : ℚ) :
    ∑ k ∈ range (n + 1 + 1), ((n + 1).choose k : ℚ) * (x ^ k * (((n + 1 - k : ℕ) : ℚ) * y ^ (n + 1 - k - 1))) =
      (n + 1) * (x + y) ^ n := by
  rw [Finset.sum_range_succ, Nat.sub_self, Nat.cast_zero, zero_mul, mul_zero, mul_zero, add_zero, add_pow,
    Finset.mul_sum]
  refine Finset.sum_congr rfl fun k hk => ?_
  have hk' : k ≤ n := Nat.lt_succ_iff.1 (Finset.mem_range.1 hk)
  have hc : ((n + 1).choose k : ℚ) * ((n + 1 - k : ℕ) : ℚ) = (n.choose k : ℚ) * (n + 1) := by
    have := Nat.choose_mul_succ_eq n k
    exact_mod_cast this.symm
  have he : n + 1 - k - 1 = n - k := by omega
  rw [he]
  linear_combination (x ^ k * y ^ (n - k)) * hc

theorem simes_grid (n : ℕ) : ∀ (T : Finset ι) (G : ℕ) (lam : ℚ), T.card = n → 0 ≤ lam →
    (∀ r : ℕ, r ≤ n → lam * r ≤ G → ∃ j : ℕ, lam * r = j) →
    (((dom T G).filter (simesEv T lam)).card : ℚ) = min (lam * n) G * (G : ℚ) ^ (n - 1) := by
  induction n using Nat.strong_induction_on with
  | _ n ih =>
  intro T G lam hT hlam hgrid
  rcases Nat.eq_zero_or_pos n with hn | hn
  · subst hn
    rw [Finset.filter_false_of_mem]
    · simp
    · intro u _
      simp only [simesEv, hT]
      simp
  · by_cases hG : (G : ℚ) ≤ lam * n
    · rw [simes_full T G lam (hT ▸ hn) (hT ▸ hG), min_eq_right hG, hT]
      push_cast
      rw [← pow_succ', Nat.sub_add_cancel hn]
    · rw [not_le] at hG
      obtain ⟨g, hg⟩ := hgrid n le_rfl hG.le
      have hgG : g < G := by
        have : (g : ℚ) < G := hg ▸ hG
        exact_mod_cast this
      set F := (dom T G).filter (simesEv T lam) with hF
      set high : (ι → ℕ) → Finset ι := fun u => T.filter fun i => g < u i with hhigh
      have hmaps : (F : Set (ι → ℕ)).MapsTo high (T.powerset : Finset (Finset ι)) := by
        intro u _
        simp only [Finset.mem_coe, Finset.mem_powerset, hhigh]
        exact Finset.filter_subset _ _
      rw [Finset.card_eq_sum_card_fiberwise hmaps]
      push_cast
      -- each fibre
      have hfib : ∀ S ∈ T.powerset, ((F.filter fun u => high u = S).card : ℚ) =
          ((G - g : ℕ) : ℚ) ^ S.card * (lam * ((n - S.card : ℕ) : ℚ) * (g : ℚ) ^ (n - S.card - 1)) := by
        intro S hS
        rw [Finset.mem_powerset] at hS
        set A : ι → Finset ℕ := fun i =>
          if i ∈ T \ S then Icc 1 g else if i ∈ S then Icc (g + 1) G else {0} with hA
        have hset : F.filter (fun u => high u = S) = (Fintype.piFinset A).filter (simesEv (T \ S) lam) := by
          ext u
          simp only [hF, Finset.mem_filter, Fintype.mem_piFinset, hA, hhigh]
          constructor
          · rintro ⟨⟨hu, hev⟩, hhs⟩
            have hu' := (mem_dom T G u).1 hu
            have hiS : ∀ i, i ∈ S ↔ i ∈ T ∧ g < u i := fun i => by
              rw [← hhs, Finset.mem_filter]
            refine ⟨fun i => ?_, ?_⟩
            · by_cases hiT : i ∈ T
              · by_cases his : i ∈ S
                · have := (hiS i).1 his
                  simp [his, Finset.mem_sdiff]
                  exact ⟨by omega, (hu' i).1 hiT |>.2⟩
                · have hle : u i ≤ g := by
                    by_contra hc
                    exact his ((hiS i).2 ⟨hiT, by omega⟩)
                  simp [his, hiT, Finset.mem_sdiff]
                  exact ⟨(hu' i).1 hiT |>.1, hle⟩
              · have his : i ∉ S := fun h => hiT (hS h)
                simp [his, hiT, Finset.mem_sdiff, (hu' i).2 hiT]
            · obtain ⟨r, hr, hrc⟩ := hev
              have hrn : r ≤ n := hT ▸ (Finset.mem_Icc.1 hr).2
              have hcnt : cntT T u (lam * r) = cntT (T \ S) u (lam * r) := by
                unfold cntT
                congr 1
                ext i
                simp only [Finset.mem_filter, Finset.mem_sdiff]
                constructor
                · rintro ⟨hiT, hle⟩
                  refine ⟨⟨hiT, fun his => ?_⟩, hle⟩
                  have := ((hiS i).1 his).2
                  have h1 : (g : ℚ) < u i := by exact_mod_cast this
                  have h2 : lam * r ≤ lam * n := mul_le_mul_of_nonneg_left (by exact_mod_cast hrn) hlam
                  linarith
                · rintro ⟨⟨hiT, _⟩, hle⟩
                  exact ⟨hiT, hle⟩
              refine ⟨r, Finset.mem_Icc.2 ⟨(Finset.mem_Icc.1 hr).1, ?_⟩, hcnt ▸ hrc⟩
              exact le_trans (hcnt ▸ hrc) (cntT_le _ _ _)
          · rintro ⟨hu, hev⟩
            have hmem : ∀ i, (i ∈ T → 1 ≤ u i ∧ u i ≤ G) ∧ (i ∉ T → u i = 0) ∧ (i ∈ T → (g < u i ↔ i ∈ S)) := by
              intro i
              have := hu i
              by_cases hiT : i ∈ T
              · by_cases his : i ∈ S
                · simp [his, Finset.mem_sdiff] at this
                  refine ⟨fun _ => ⟨by omega, this.2⟩, fun h => absurd hiT h, fun _ => ⟨fun _ => his, fun _ => by omega⟩⟩
                · simp [his, hiT, Finset.mem_sdiff] at this
                  refine ⟨fun _ => ⟨this.1, by omega⟩, fun h => absurd hiT h, fun _ => ⟨fun h => by omega, fun h => absurd h his⟩⟩
              · have his : i ∉ S := fun h => hiT (hS h)
                simp [his, hiT, Finset.mem_sdiff] at this
                exact ⟨fun h => absurd h hiT, fun _ => this, fun h => absurd h hiT⟩
            refine ⟨⟨(mem_dom T G u).2 fun i => ⟨(hmem i).1, (hmem i).2.1⟩, ?_⟩, ?_⟩
            · obtain ⟨r, hr, hrc⟩ := hev
              have hTS : (T \ S).card ≤ T.card := Finset.card_le_card Finset.sdiff_subset
              refine ⟨r, Finset.mem_Icc.2 ⟨(Finset.mem_Icc.1 hr).1, le_trans (Finset.mem_Icc.1 hr).2 hTS⟩,
                le_trans hrc (Finset.card_le_card ?_)⟩
              intro i hi
              simp only [Finset.mem_filter, Finset.mem_sdiff] at hi ⊢
              exact ⟨hi.1.1, hi.2⟩
            · ext i
              simp only [Finset.mem_filter]
              constructor
              · rintro ⟨hiT, hgi⟩
                exact ((hmem i).2.2 hiT).1 hgi
              · intro his
                exact ⟨hS his, ((hmem i).2.2 (hS his)).2 his⟩
        rw [hset]
        have hP : ∀ u v : ι → ℕ, (∀ i, i ∉ S → u i = v i) → (simesEv (T \ S) lam u ↔ simesEv (T \ S) lam v) := by
          intro u v h
          have hc : ∀ t, cntT (T \ S) u t = cntT (T \ S) v t := by
            intro t
            unfold cntT
            congr 1
            ext i
            simp only [Finset.mem_filter, Finset.mem_sdiff]
            constructor
            · rintro ⟨⟨hiT, his⟩, hle⟩
              exact ⟨⟨hiT, his⟩, h i his ▸ hle⟩
            · rintro ⟨⟨hiT, his⟩, hle⟩
              exact ⟨⟨hiT, his⟩, (h i his).symm ▸ hle⟩
          simp only [simesEv, hc]
        rw [freeze_set (simesEv (T \ S) lam) S A hP]
        have hprod : ∏ i ∈ S, (A i).card = (G - g) ^ S.card := by
          rw [Finset.prod_congr rfl fun i hi => (show (A i).card = G - g by
            have : i ∉ T \ S := fun h => (Finset.mem_sdiff.1 h).2 hi
            simp only [hA, this, hi, ↓reduceIte, Nat.card_Icc]
            omega)]
          rw [Finset.prod_const]
        have hdom : (Fintype.piFinset fun i => if i ∈ S then ({0} : Finset ℕ) else A i) = dom (T \ S) g := by
          unfold dom
          congr 1
          funext i
          by_cases his : i ∈ S
          · have : i ∉ T \ S := fun h => (Finset.mem_sdiff.1 h).2 his
            simp [his, this]
          · simp [his, hA]
        rw [hprod, hdom]
        push_cast
        have hcard : (T \ S).card = n - S.card := by
          rw [Finset.card_sdiff_of_subset hS, hT]
        -- the inner count
        have hinner : (((dom (T \ S) g).filter (simesEv (T \ S) lam)).card : ℚ) =
            lam * ((n - S.card : ℕ) : ℚ) * (g : ℚ) ^ (n - S.card - 1) := by
          by_cases hS0 : S.card = 0
          · have hSe : S = ∅ := Finset.card_eq_zero.1 hS0
            have hTS : T \ S = T := by rw [hSe, Finset.sdiff_empty]
            rw [hTS, simes_full T g lam (hT ▸ hn) (by rw [hT, hg])]
            rw [hS0, Nat.sub_zero, hT, ← hg]
            push_cast
            rw [mul_comm, ← pow_succ, Nat.sub_add_cancel hn, hg]
          · rw [ih (n - S.card) (by omega) (T \ S) g lam hcard hlam (fun r hr hrg =>
              hgrid r (by omega) (le_trans hrg (by exact_mod_cast hgG.le)))]
            rw [min_eq_left]
            rw [hg.symm]
            exact mul_le_mul_of_nonneg_left (by exact_mod_cast Nat.sub_le n S.card) hlam
        rw [hinner]
      have hsum := Finset.sum_powerset_apply_card
        (fun k => ((G - g : ℕ) : ℚ) ^ k * (lam * ((n - k : ℕ) : ℚ) * (g : ℚ) ^ (n - k - 1))) (x := T)
      rw [Finset.sum_congr rfl hfib, hsum, hT]
      obtain ⟨n', rfl⟩ : ∃ n', n = n' + 1 := ⟨n - 1, by omega⟩
      have hx : ((G - g : ℕ) : ℚ) + g = G := by
        rw [Nat.cast_sub hgG.le]
        ring
      have hb := binom_deriv n' ((G - g : ℕ) : ℚ) g
      rw [hx] at hb
      rw [min_eq_left hG.le, hg, Nat.add_sub_cancel]
      have hfac : ∑ k ∈ range (n' + 1 + 1), ((n' + 1).choose k) • (((G - g : ℕ) : ℚ) ^ k *
          (lam * ((n' + 1 - k : ℕ) : ℚ) * (g : ℚ) ^ (n' + 1 - k - 1))) =
          lam * ∑ k ∈ range (n' + 1 + 1), ((n' + 1).choose k : ℚ) * (((G - g : ℕ) : ℚ) ^ k *
            (((n' + 1 - k : ℕ) : ℚ) * (g : ℚ) ^ (n' + 1 - k - 1))) := by
        rw [Finset.mul_sum]
        refine Finset.sum_congr rfl fun k _ => ?_
        rw [nsmul_eq_mul]
        ring
      rw [hfac, hb, ← hg]
      push_cast
      ring

end Simes

/-! ## The false discovery rate of the kept family -/

section Formula

/-- An outcome as a natural number: `0` is `p = 1`, `j ≥ 1` is `p = a j / N`. -/
def valN (a : ℚ) (N j : ℕ) : ℚ := if j = 0 then 1 else a * j / N

/-- Its probability. -/
def wtN (a : ℚ) (N j : ℕ) : ℚ := if j = 0 then 1 - a else a / N

/-- The kept terms of an outcome vector. -/
def keptSet {m : ℕ} (v : Fin m → ℕ) : Finset (Fin m) := univ.filter fun b => v b ≠ 0

/-- The thresholds' slope on the outcome scale: `a j / N ≤ r q / k ↔ j ≤ λ r`. -/
def lamOf (N : ℕ) (a q : ℚ) (k : ℕ) : ℚ := N * q / (a * k)

theorem valN_nonneg (a : ℚ) (N j : ℕ) (ha : 0 < a) : 0 ≤ valN a N j := by
  unfold valN
  split_ifs
  · exact zero_le_one
  · positivity

theorem valN_lt_one_iff (a : ℚ) (N j : ℕ) (hN : 0 < N) (ha : 0 < a) (ha1 : a < 1) (hj : j ≤ N) :
    valN a N j < 1 ↔ j ≠ 0 := by
  unfold valN
  by_cases h : j = 0
  · simp [h]
  · simp only [h, ↓reduceIte, ne_eq, not_false_eq_true, iff_true]
    have hN' : (0 : ℚ) < N := by exact_mod_cast hN
    rw [div_lt_one hN']
    have : (j : ℚ) ≤ N := by exact_mod_cast hj
    nlinarith

theorem event_eq {m N : ℕ} (a q : ℚ) (hN : 0 < N) (ha : 0 < a) (ha1 : a < 1) (v : Fin m → ℕ)
    (hv : ∀ b, v b ≤ N) :
    0 < callsKept false q (List.ofFn fun b => valN a N (v b)) ↔
      simesEv (keptSet v) (lamOf N a q (keptSet v).card) v := by
  rw [kept_event_iff _ (fun b => valN_nonneg a N (v b) ha)]
  have hS : (univ.filter fun b => valN a N (v b) < 1) = keptSet v := by
    unfold keptSet
    congr 1
    funext b
    rw [valN_lt_one_iff a N (v b) hN ha ha1 (hv b)]
  have hN' : (0 : ℚ) < N := by exact_mod_cast hN
  have hcnt : ∀ r : ℕ, (univ.filter fun b => valN a N (v b) < 1 ∧
      valN a N (v b) ≤ r * q / ((univ.filter fun b => valN a N (v b) < 1).card : ℚ)).card =
      cntT (keptSet v) v (lamOf N a q (keptSet v).card * r) := by
    intro r
    rw [hS]
    unfold cntT keptSet
    rw [Finset.filter_filter]
    congr 1
    apply Finset.filter_congr
    intro b _
    rw [valN_lt_one_iff a N (v b) hN ha ha1 (hv b)]
    apply and_congr_right
    intro hb
    simp only [valN, hb, ↓reduceIte, lamOf]
    rw [div_le_iff₀ hN']
    constructor
    · intro h
      have : (v b : ℚ) * a ≤ r * q / ((univ.filter fun b => v b ≠ 0).card : ℚ) * N := by linarith
      calc (v b : ℚ) = (v b : ℚ) * a / a := by field_simp
        _ ≤ r * q / ((univ.filter fun b => v b ≠ 0).card : ℚ) * N / a := by
            exact div_le_div_of_nonneg_right this ha.le
        _ = _ := by ring
    · intro h
      have h1 : (v b : ℚ) * a ≤ N * q / (a * ((univ.filter fun b => v b ≠ 0).card : ℚ)) * r * a :=
        mul_le_mul_of_nonneg_right h ha.le
      have h2 : N * q / (a * ((univ.filter fun b => v b ≠ 0).card : ℚ)) * r * a =
          r * q / ((univ.filter fun b => v b ≠ 0).card : ℚ) * N := by
        field_simp
      linarith
  simp only [hcnt]
  constructor
  · rintro ⟨r, hr, hrc⟩
    exact ⟨r, Finset.mem_Icc.2 ⟨hr, le_trans hrc (cntT_le _ _ _)⟩, hrc⟩
  · rintro ⟨r, hr, hrc⟩
    exact ⟨r, (Finset.mem_Icc.1 hr).1, hrc⟩

theorem mixVal_eq (a : ℚ) (N : ℕ) (j : Fin (N + 1)) : mixVal a N j = valN a N j := by
  simp [mixVal, valN]

theorem mixWt_eq (a : ℚ) (N : ℕ) (j : Fin (N + 1)) : mixWt a N j = wtN a N j := by
  simp [mixWt, wtN]

/-- The outcomes as natural-number vectors. -/
theorem fdrKeptMix_eq (m N : ℕ) (a q : ℚ) :
    fdrKeptMix m N a q false = ∑ v ∈ Fintype.piFinset (fun _ : Fin m => range (N + 1)),
      (∏ b, wtN a N (v b)) * (if 0 < callsKept false q (List.ofFn fun b => valN a N (v b)) then 1 else 0) := by
  unfold fdrKeptMix
  apply Finset.sum_nbij' (fun ω b => ((ω b : Fin (N + 1)) : ℕ))
    (fun v b => (⟨v b % (N + 1), Nat.mod_lt _ (Nat.succ_pos N)⟩ : Fin (N + 1)))
  · intro ω _
    rw [Fintype.mem_piFinset]
    intro b
    exact Finset.mem_range.2 (ω b).isLt
  · intro v _
    exact Finset.mem_univ _
  · intro ω _
    funext b
    exact Fin.ext (Nat.mod_eq_of_lt (ω b).isLt)
  · intro v hv
    rw [Fintype.mem_piFinset] at hv
    funext b
    exact Nat.mod_eq_of_lt (Finset.mem_range.1 (hv b))
  · intro ω _
    simp only [mixWt_eq, mixVal_eq]

theorem fiber_eq_dom {m N : ℕ} (S : Finset (Fin m)) :
    (Fintype.piFinset fun _ : Fin m => range (N + 1)).filter (fun v => keptSet v = S) = dom S N := by
  ext v
  rw [Finset.mem_filter, Fintype.mem_piFinset, mem_dom]
  constructor
  · rintro ⟨hv, hS⟩
    intro b
    have hb : b ∈ S ↔ v b ≠ 0 := by rw [← hS]; simp [keptSet]
    have := Finset.mem_range.1 (hv b)
    refine ⟨fun h => ⟨by have := hb.1 h; omega, by omega⟩, fun h => ?_⟩
    by_contra hc
    exact h (hb.2 hc)
  · intro h
    refine ⟨fun b => ?_, ?_⟩
    · rw [Finset.mem_range]
      by_cases hb : b ∈ S
      · have := ((h b).1 hb).2
        omega
      · rw [(h b).2 hb]
        omega
    · ext b
      simp only [keptSet, Finset.mem_filter, Finset.mem_univ, true_and]
      constructor
      · intro hb
        by_contra hc
        exact hb ((h b).2 hc)
      · intro hb
        have := ((h b).1 hb).1
        omega

theorem prod_wt_dom {m N : ℕ} (a : ℚ) (S : Finset (Fin m)) (v : Fin m → ℕ) (hv : v ∈ dom S N) :
    ∏ b, wtN a N (v b) = (a / N) ^ S.card * (1 - a) ^ (m - S.card) := by
  rw [mem_dom] at hv
  have h : ∀ b, wtN a N (v b) = if b ∈ S then a / N else 1 - a := by
    intro b
    unfold wtN
    by_cases hb : b ∈ S
    · have := ((hv b).1 hb).1
      simp [hb, show v b ≠ 0 by omega]
    · simp [hb, (hv b).2 hb]
  rw [Finset.prod_congr rfl fun b _ => h b, Finset.prod_ite, Finset.prod_const, Finset.prod_const]
  congr 2
  · congr 1
    ext b
    simp
  · rw [Finset.filter_not, Finset.card_sdiff_of_subset (Finset.filter_subset _ _)]
    simp

/-- (4d-ii) The kept-family formula. -/
theorem kept_family_formula : kept_family_formula_stmt := by
  intro m N a q hm hN ha ha1 hq hgrid
  have hN' : (0 : ℚ) < N := by exact_mod_cast hN
  rw [fdrKeptMix_eq]
  have hmaps : ∀ v ∈ Fintype.piFinset (fun _ : Fin m => range (N + 1)),
      keptSet v ∈ (univ : Finset (Fin m)).powerset := fun v _ => Finset.mem_powerset.2 (Finset.subset_univ _)
  rw [← Finset.sum_fiberwise_of_maps_to hmaps]
  -- each fibre
  set h : ℕ → ℚ := fun k => (a / N) ^ k * (1 - a) ^ (m - k) * (min (lamOf N a q k * k) N * (N : ℚ) ^ (k - 1))
    with hh
  have hfib : ∀ S ∈ (univ : Finset (Fin m)).powerset,
      (∑ v ∈ (Fintype.piFinset fun _ : Fin m => range (N + 1)).filter (fun v => keptSet v = S),
        (∏ b, wtN a N (v b)) * (if 0 < callsKept false q (List.ofFn fun b => valN a N (v b)) then 1 else 0)) =
        h S.card := by
    intro S _
    rw [fiber_eq_dom]
    have hterm : ∀ v ∈ dom S N, (∏ b, wtN a N (v b)) *
        (if 0 < callsKept false q (List.ofFn fun b => valN a N (v b)) then (1 : ℚ) else 0) =
        (a / N) ^ S.card * (1 - a) ^ (m - S.card) *
          (if simesEv S (lamOf N a q S.card) v then 1 else 0) := by
      intro v hv
      have hvN : ∀ b, v b ≤ N := by
        intro b
        have := (mem_dom S N v).1 hv b
        by_cases hb : b ∈ S
        · exact (this.1 hb).2
        · rw [this.2 hb]
          exact Nat.zero_le _
      have hkS : keptSet v = S := by
        have := Finset.mem_filter.1 ((fiber_eq_dom (N := N) S).symm ▸ hv)
        exact this.2
      rw [prod_wt_dom a S v hv]
      congr 1
      have he := event_eq a q hN ha ha1 v hvN
      rw [hkS] at he
      simp only [he]
    rw [Finset.sum_congr rfl hterm, ← Finset.mul_sum, Finset.sum_boole]
    rw [simes_grid S.card S N (lamOf N a q S.card) rfl (by unfold lamOf; positivity)]
    intro r hr hrN
    rcases Nat.eq_zero_or_pos r with h0 | hpos
    · exact ⟨0, by rw [h0]; simp⟩
    · have hk : 1 ≤ S.card := le_trans hpos hr
      have hk' : (0 : ℚ) < S.card := by exact_mod_cast hk
      have hthr : q / a * ((r : ℚ) / S.card) ≤ 1 := by
        have : lamOf N a q S.card * r = N * (q / a * ((r : ℚ) / S.card)) := by
          unfold lamOf
          field_simp
        rw [this] at hrN
        have := (mul_le_iff_le_one_right hN').1 (by simpa using hrN)
        exact this
      obtain ⟨j, hj⟩ := hgrid S.card r hpos hr (le_trans (Finset.card_le_univ S) (by simp)) hthr
      refine ⟨j, ?_⟩
      unfold lamOf
      rw [show N * q / (a * S.card) * r = N * (q / a * ((r : ℚ) / S.card)) by field_simp, hj]
      field_simp
  rw [Finset.sum_congr rfl hfib, Finset.sum_powerset_apply_card h, Finset.card_univ, Fintype.card_fin]
  -- the binomial sum
  rw [Finset.sum_range_succ']
  have h0 : h 0 = 0 := by
    rw [hh]
    simp
  have hk : ∀ k : ℕ, h (k + 1) = a ^ (k + 1) * (1 - a) ^ (m - (k + 1)) * min (q / a) 1 := by
    intro k
    rw [hh]
    simp only
    have hl : lamOf N a q (k + 1) * ((k + 1 : ℕ) : ℚ) = N * (q / a) := by
      unfold lamOf
      push_cast
      field_simp
    have hmin : min ((N : ℚ) * (q / a)) N = N * min (q / a) 1 := by
      rw [mul_min_of_nonneg _ _ hN'.le, mul_one]
    rw [hl, hmin, Nat.add_sub_cancel, div_pow]
    field_simp
    ring
  simp only [h0, smul_zero, add_zero, hk, nsmul_eq_mul]
  have hbin := add_pow a (1 - a) m
  rw [add_sub_cancel, one_pow, Finset.sum_range_succ'] at hbin
  simp only [pow_zero, one_mul, Nat.sub_zero, Nat.choose_zero_right, Nat.cast_one, mul_one] at hbin
  have : ∑ k ∈ range m, ((m.choose (k + 1) : ℕ) : ℚ) * (a ^ (k + 1) * (1 - a) ^ (m - (k + 1)) * min (q / a) 1) =
      min (q / a) 1 * ∑ k ∈ range m, a ^ (k + 1) * (1 - a) ^ (m - (k + 1)) * (m.choose (k + 1) : ℚ) := by
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun k _ => ?_
    ring
  rw [this]
  linear_combination (-(min (q / a) 1)) * hbin

end Formula

/-! ## The formula fails at `a = 1` -/

theorem ex_r_le_two (c : ℕ → ℕ) (hc : ∀ r, c r ≤ 2) : (∃ r : ℕ, 1 ≤ r ∧ r ≤ c r) ↔ (1 ≤ c 1 ∨ 2 ≤ c 2) := by
  constructor
  · rintro ⟨r, h1, h2⟩
    have h3 : r ≤ 2 := le_trans h2 (hc r)
    interval_cases r
    · exact Or.inl h2
    · exact Or.inr h2
  · rintro (h | h)
    · exact ⟨1, le_rfl, h⟩
    · exact ⟨2, by norm_num, h⟩

/-- BH over the kept terms of two p-values. -/
theorem event_two (v : Fin 2 → ℚ) (hv : ∀ b, 0 ≤ v b) (q : ℚ) :
    0 < callsKept false q (List.ofFn v) ↔
      (1 ≤ (univ.filter fun b => v b < 1 ∧ v b ≤ (1 : ℕ) * q / ((univ.filter fun b => v b < 1).card : ℚ)).card ∨
        2 ≤ (univ.filter fun b => v b < 1 ∧ v b ≤ (2 : ℕ) * q / ((univ.filter fun b => v b < 1).card : ℚ)).card) := by
  rw [kept_event_iff v hv q]
  exact ex_r_le_two (fun r => (univ.filter fun b => v b < 1 ∧
    v b ≤ r * q / ((univ.filter fun b => v b < 1).card : ℚ)).card)
    fun r => le_trans (Finset.card_filter_le _ _) (by simp)

/-- (4d-ii) At `a = 1` the formula's `1/2` is not the kept family's `5/8`. -/
theorem kept_family_unit_a_counterexample : kept_family_unit_a_counterexample_stmt := by
  refine ⟨?_, ?_, by norm_num⟩
  · intro K i hi hiK hK _
    interval_cases K <;> interval_cases i
    · exact ⟨2, by norm_num⟩
    · exact ⟨1, by norm_num⟩
    · exact ⟨2, by norm_num⟩
  · unfold fdrKeptMix
    rw [← Equiv.sum_comp (piFinTwoEquiv fun _ => Fin 5).symm, Fintype.sum_prod_type]
    have hnn : ∀ (x y : Fin 5) (b : Fin 2),
        0 ≤ mixVal 1 4 ((piFinTwoEquiv fun _ => Fin 5).symm (x, y) b) := by
      intro x y b
      unfold mixVal
      split_ifs
      · exact zero_le_one
      · positivity
    simp only [fun x y => event_two _ (hnn x y)]
    simp only [Fin.sum_univ_five, Fin.prod_univ_two, piFinTwoEquiv_symm_apply, Fin.cons_zero, Fin.cons_one,
      Finset.card_filter, Fin.sum_univ_two, mixVal, mixWt]
    norm_num

/-! ## KIT-189's instance -/

theorem dvd_foldl_lcm_acc : ∀ (l : List ℕ) (acc : ℕ), acc ∣ l.foldl Nat.lcm acc
  | [], _ => dvd_refl _
  | x :: l, acc => (Nat.dvd_lcm_left acc x).trans (dvd_foldl_lcm_acc l (Nat.lcm acc x))

theorem dvd_foldl_lcm_mem : ∀ (l : List ℕ) (acc x : ℕ), x ∈ l → x ∣ l.foldl Nat.lcm acc
  | [], _, _, h => absurd h List.not_mem_nil
  | y :: l, acc, x, h => by
    rcases List.mem_cons.1 h with rfl | h
    · exact (Nat.dvd_lcm_right acc x).trans (dvd_foldl_lcm_acc l _)
    · exact dvd_foldl_lcm_mem l _ x h

theorem foldl_lcm_pos : ∀ (l : List ℕ) (acc : ℕ), 0 < acc → (∀ x ∈ l, 0 < x) → 0 < l.foldl Nat.lcm acc
  | [], _, h, _ => h
  | x :: l, _, h, hl =>
    foldl_lcm_pos l _ (Nat.lcm_pos h (hl x List.mem_cons_self)) fun y hy => hl y (List.mem_cons_of_mem _ hy)

theorem lcmUpTo_pos (n : ℕ) : 0 < lcmUpTo n :=
  foldl_lcm_pos _ 1 one_pos fun x hx => by
    rw [List.mem_range'_1] at hx
    omega

theorem dvd_lcmUpTo (n K : ℕ) (h1 : 1 ≤ K) (h2 : K ≤ n) : K ∣ lcmUpTo n :=
  dvd_foldl_lcm_mem _ 1 K (by rw [List.mem_range'_1]; omega)

/-- Witness for `kept_family_formula_stmt`: KIT-189's instance satisfies every hypothesis. -/
theorem kept_family_witness : kept_family_witness_stmt := by
  have hL := lcmUpTo_pos 150
  refine ⟨by norm_num, by omega, by norm_num, by norm_num, by norm_num, ?_⟩
  intro K i hi hiK hK _
  obtain ⟨t, ht⟩ := dvd_lcmUpTo 150 K (le_trans hi hiK) hK
  have hK0 : (K : ℚ) ≠ 0 := by
    have : 0 < K := lt_of_lt_of_le hi hiK
    exact_mod_cast this.ne'
  have hL0 : (lcmUpTo 150 : ℚ) ≠ 0 := by exact_mod_cast hL.ne'
  refine ⟨5 * i * t, ?_⟩
  have hLt : (lcmUpTo 150 : ℚ) = K * t := by exact_mod_cast ht
  push_cast
  rw [hLt] at hL0 ⊢
  have ht0 : (t : ℚ) ≠ 0 := right_ne_zero_of_mul hL0
  field_simp
  ring

/-- (4d-ii) KIT-189's figure: `(5/27)(1 - (73/100)^150) > 0.185`. -/
theorem kept_family_instance : kept_family_instance_stmt := by
  obtain ⟨hm, hN, ha, ha1, hq, hgrid⟩ := kept_family_witness
  have hf := kept_family_formula 150 (27 * lcmUpTo 150) (27 / 100) (1 / 20) hm hN ha ha1 hq hgrid
  have hmin : min ((1 / 20 : ℚ) / (27 / 100)) 1 = 5 / 27 := by norm_num
  rw [hmin, show (1 : ℚ) - 27 / 100 = 73 / 100 by norm_num] at hf
  refine ⟨hf, ?_⟩
  have hp : ((73 : ℚ) / 100) ^ 150 ≤ (73 / 100) ^ 30 :=
    pow_le_pow_of_le_one (by norm_num) (by norm_num) (by norm_num)
  have h30 : ((73 : ℚ) / 100) ^ 30 < 1 / 1000 := by norm_num
  linarith

end MultiScholaRProofs.QvalueBh
