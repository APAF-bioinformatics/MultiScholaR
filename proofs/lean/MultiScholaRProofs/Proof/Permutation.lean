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

/-! ## (5a) The finite core

The plan recorded in `outcome.md`: stack the relabelling and the draws as `σ = (τ, τ π_1, …, τ π_B)` (an injection of
`Perm × Perm^B` into `Perm^(B+1)`); the count against index 0 of `σ` is at most the p-value's count; the `B + 1`
coordinate swaps show every index has as many tuples with a small count as index 0; and within one tuple at most `A`
indices have `1 + count ≤ A` (the one of least statistic among them has every other counted against it). -/

/-- `cntGe s i`: the indices other than `i` whose value is at least `s i`. -/
noncomputable def cntGe {m : ℕ} (s : Fin m → ℝ) (i : Fin m) : ℕ :=
  (univ.filter fun j => j ≠ i ∧ s i ≤ s j).card

/-- Reindexing the values by a permutation `e` moves the count from `i` to `e i`. -/
theorem cntGe_comp {m : ℕ} (s : Fin m → ℝ) (e : Equiv.Perm (Fin m)) (i : Fin m) :
    cntGe (s ∘ e) i = cntGe s (e i) := by
  unfold cntGe
  refine card_nbij' e e.symm (fun j hj => ?_) (fun j hj => ?_) (fun j _ => by simp) (fun j _ => by simp)
  · simp only [coe_filter, mem_univ, true_and, Set.mem_ofPred_eq, Function.comp_apply] at hj ⊢
    exact ⟨fun h => hj.1 (e.injective h), hj.2⟩
  · simp only [coe_filter, mem_univ, true_and, Set.mem_ofPred_eq, Function.comp_apply,
      Equiv.apply_symm_apply] at hj ⊢
    exact ⟨fun h => hj.1 (by rw [← h, Equiv.apply_symm_apply]), hj.2⟩

/-- Within one tuple of values, at most `A` indices have `1 + count ≤ A` (ties counted against). -/
theorem card_small_count {m : ℕ} (s : Fin m → ℝ) (A : ℝ) (hA : 0 ≤ A) :
    ((univ.filter fun i => 1 + (cntGe s i : ℝ) ≤ A).card : ℝ) ≤ A := by
  set I := univ.filter fun i => 1 + (cntGe s i : ℝ) ≤ A with hI
  rcases I.eq_empty_or_nonempty with h | h
  · rw [h, card_empty, Nat.cast_zero]
    exact hA
  obtain ⟨i0, hi0, hmin⟩ := I.exists_min_image s h
  have hsub : I.erase i0 ⊆ univ.filter fun j => j ≠ i0 ∧ s i0 ≤ s j := by
    intro j hj
    rw [mem_erase] at hj
    exact mem_filter.2 ⟨mem_univ _, hj.1, hmin j hj.2⟩
  have h1 : I.card - 1 ≤ cntGe s i0 := by
    rw [← card_erase_of_mem hi0]
    exact card_le_card hsub
  have h2 : 1 + (cntGe s i0 : ℝ) ≤ A := (mem_filter.1 hi0).2
  have h3 : 1 ≤ I.card := card_pos.2 h
  have h4 : I.card ≤ 1 + cntGe s i0 := by omega
  have h5 : (I.card : ℝ) ≤ 1 + cntGe s i0 := by exact_mod_cast h4
  linarith

/-- The statistics of a stacked tuple of `B + 1` permutations applied to `g`. -/
noncomputable def stats {n B : ℕ} {X : Type} (g : Fin n → X) (T : (Fin n → X) → ℝ)
    (σ : Fin (B + 1) → Equiv.Perm (Fin n)) : Fin (B + 1) → ℝ :=
  fun j => T (g ∘ ⇑(σ j))

/-- The tuples whose index `i` has `1 + count ≤ A`. -/
noncomputable def smallAt {n B : ℕ} {X : Type} (g : Fin n → X) (T : (Fin n → X) → ℝ) (A : ℝ) (i : Fin (B + 1)) :
    Finset (Fin (B + 1) → Equiv.Perm (Fin n)) :=
  univ.filter fun σ => 1 + (cntGe (stats g T σ) i : ℝ) ≤ A

/-- The stacking map `(τ, π) ↦ (τ, τ π_1, …, τ π_B)`. -/
def stack {n B : ℕ} (x : Equiv.Perm (Fin n) × (Fin B → Equiv.Perm (Fin n))) : Fin (B + 1) → Equiv.Perm (Fin n) :=
  Fin.cons x.1 fun b => x.1 * x.2 b

theorem stack_injective {n B : ℕ} : Function.Injective (stack (n := n) (B := B)) := by
  intro x y hxy
  have h0 : x.1 = y.1 := by simpa [stack] using congrFun hxy 0
  refine Prod.ext h0 (funext fun b => ?_)
  have hb := congrFun hxy b.succ
  simp only [stack, Fin.cons_succ] at hb
  rw [h0] at hb
  exact mul_left_cancel hb

/-- Every index has as many tuples with a small count as index 0 (swap coordinates `0` and `i`). -/
theorem card_smallAt {n B : ℕ} {X : Type} (g : Fin n → X) (T : (Fin n → X) → ℝ) (A : ℝ) (i : Fin (B + 1)) :
    (smallAt g T A i).card = (smallAt (B := B) g T A 0).card := by
  have hst : ∀ σ : Fin (B + 1) → Equiv.Perm (Fin n),
      stats g T (σ ∘ ⇑(Equiv.swap 0 i)) = stats g T σ ∘ ⇑(Equiv.swap 0 i) := fun σ => rfl
  have hinv : ∀ σ : Fin (B + 1) → Equiv.Perm (Fin n), (σ ∘ ⇑(Equiv.swap 0 i)) ∘ ⇑(Equiv.swap 0 i) = σ :=
    fun σ => funext fun j => by simp
  unfold smallAt
  refine card_nbij' (fun σ => σ ∘ ⇑(Equiv.swap 0 i)) (fun σ => σ ∘ ⇑(Equiv.swap 0 i))
    (fun σ hσ => ?_) (fun σ hσ => ?_) (fun σ _ => hinv σ) (fun σ _ => hinv σ)
  · simp only [coe_filter, mem_univ, true_and, Set.mem_ofPred_eq] at hσ ⊢
    rw [hst, cntGe_comp, Equiv.swap_apply_left]
    exact hσ
  · simp only [coe_filter, mem_univ, true_and, Set.mem_ofPred_eq] at hσ ⊢
    rw [hst, cntGe_comp, Equiv.swap_apply_right]
    exact hσ

/-- The count against index 0 of the stacked tuple is at most the p-value's count. -/
theorem cntGe_stack_le {n B : ℕ} {X : Type} (g : Fin n → X) (T : (Fin n → X) → ℝ) (ε : ℝ) (hε : 0 ≤ ε)
    (x : Equiv.Perm (Fin n) × (Fin B → Equiv.Perm (Fin n))) (P : Fin B → Prop) {hP : DecidablePred P}
    (hPdef : ∀ b, P b ↔ (T ((g ∘ ⇑x.1) ∘ ⇑(x.2 b)) ≥ T (g ∘ ⇑x.1) - ε ∨ (g ∘ ⇑x.1) ∘ ⇑(x.2 b) = g ∘ ⇑x.1)) :
    cntGe (stats g T (stack x)) 0 ≤ (@Finset.filter _ P hP univ).card := by
  classical
  unfold cntGe
  refine le_trans (card_le_card (t := (@Finset.filter _ P hP univ).image Fin.succ) fun j hj => ?_) card_image_le
  rw [mem_filter] at hj
  obtain ⟨b, rfl⟩ := Fin.exists_succ_eq.2 hj.2.1
  refine mem_image.2 ⟨b, mem_filter.2 ⟨mem_univ _, (hPdef b).2 (Or.inl ?_)⟩, rfl⟩
  have h := hj.2.2
  simp only [stats, stack, Fin.cons_zero, Fin.cons_succ, Equiv.Perm.coe_mul] at h
  rw [Function.comp_assoc]
  linarith

/-- (5a) the finite core. -/
theorem perm_finite_core : perm_finite_core_stmt := by
  intro n B X g T α ε hα hε
  set A : ℝ := α * (1 + B) with hA
  have hA0 : 0 ≤ A := by positivity
  have h1 : (univ.filter fun x : Equiv.Perm (Fin n) × (Fin B → Equiv.Perm (Fin n)) =>
      pval T ε (g ∘ x.1) x.2 ≤ α).card ≤ (smallAt (B := B) g T A 0).card := by
    refine card_le_card_of_injOn stack (fun x hx => ?_) (fun x _ y _ h => stack_injective h)
    simp only [coe_filter, mem_univ, true_and, Set.mem_ofPred_eq] at hx
    simp only [smallAt, coe_filter, mem_univ, true_and, Set.mem_ofPred_eq]
    unfold pval at hx
    rw [div_le_iff₀ (by positivity)] at hx
    exact le_trans ((add_le_add_iff_left 1).2 (Nat.cast_le.2 (cntGe_stack_le g T ε hε x _ (fun b => Iff.rfl)))) hx
  have h3 : ∑ i : Fin (B + 1), (smallAt g T A i).card =
      ∑ σ : Fin (B + 1) → Equiv.Perm (Fin n), (univ.filter fun i => 1 + (cntGe (stats g T σ) i : ℝ) ≤ A).card := by
    simp only [smallAt, card_filter]
    exact sum_comm
  have hcard : (Fintype.card (Fin (B + 1) → Equiv.Perm (Fin n)) : ℝ) = (n.factorial : ℝ) ^ (B + 1) := by
    simp [Fintype.card_perm]
  have key : ((B : ℝ) + 1) * (smallAt (B := B) g T A 0).card ≤ ((B : ℝ) + 1) * (α * (n.factorial : ℝ) ^ (B + 1)) := by
    calc ((B : ℝ) + 1) * (smallAt (B := B) g T A 0).card = ∑ i : Fin (B + 1), ((smallAt g T A i).card : ℝ) := by
          simp [card_smallAt]
      _ = ∑ σ : Fin (B + 1) → Equiv.Perm (Fin n),
            ((univ.filter fun i => 1 + (cntGe (stats g T σ) i : ℝ) ≤ A).card : ℝ) := by
          exact_mod_cast h3
      _ ≤ ∑ _σ : Fin (B + 1) → Equiv.Perm (Fin n), A := sum_le_sum fun σ _ => card_small_count _ A hA0
      _ = (Fintype.card (Fin (B + 1) → Equiv.Perm (Fin n)) : ℝ) * A := by simp
      _ = ((B : ℝ) + 1) * (α * (n.factorial : ℝ) ^ (B + 1)) := by rw [hcard, hA]; ring
  have h5 : ((smallAt (B := B) g T A 0).card : ℝ) ≤ α * (n.factorial : ℝ) ^ (B + 1) :=
    le_of_mul_le_mul_left key (by positivity)
  exact le_trans (by exact_mod_cast h1) h5

/-! ## (5b) The measure-theoretic form

`Measure.prod_apply` writes the probability as the `μ`-integral of the draws' probability given the factor; the
invariance of `μ` under each relabelling `τ` moves that integrand to `f ∘ τ`; summing over the `n!` relabellings, the
integrand becomes the count of 5a over `(τ, π)` (with `g = f`) divided by `(n!)^B`, at most `α n!` at every `f`. -/

section MeasureForm

open MeasureTheory ENNReal

/-- Every set of permutations is measurable in the discrete structure, so singletons are. -/
local instance permMeasurableSingleton (n : ℕ) : MeasurableSingletonClass (Equiv.Perm (Fin n)) :=
  ⟨fun _ => MeasurableSpace.measurableSet_top⟩

/-- The event `{p_ε ≤ α}` as a set of (factor, draws). -/
def permEvent (n B : ℕ) (T : (Fin n → ℝ) → ℝ) (ε α : ℝ) : Set ((Fin n → ℝ) × (Fin B → Equiv.Perm (Fin n))) :=
  {x | pval T ε x.1 x.2 ≤ α}

theorem measurable_comp_perm {n : ℕ} (σ : Equiv.Perm (Fin n)) : Measurable fun f : Fin n → ℝ => f ∘ ⇑σ :=
  Measurable.of_eval fun i => measurable_pi_apply (σ i)

/-- For fixed draws, `p_ε` is a measurable function of the factor. -/
theorem measurable_pval_left {n B : ℕ} (T : (Fin n → ℝ) → ℝ) (hT : Measurable T) (ε : ℝ)
    (π : Fin B → Equiv.Perm (Fin n)) : Measurable fun f : Fin n → ℝ => pval T ε f π := by
  have hset : ∀ b, MeasurableSet {f : Fin n → ℝ | T (f ∘ ⇑(π b)) ≥ T f - ε ∨ f ∘ ⇑(π b) = f} := fun b => by
    rw [Set.ofPred_or]
    exact (measurableSet_le (hT.sub_const ε) (hT.comp (measurable_comp_perm (π b)))).union
      (measurableSet_eq_fun (measurable_comp_perm (π b)) measurable_id)
  unfold pval
  refine Measurable.div_const (measurable_const.add ?_) _
  simp only [card_filter, Nat.cast_sum, Nat.cast_ite, Nat.cast_one, Nat.cast_zero]
  exact Finset.measurable_sum _ fun b _ => Measurable.ite (hset b) measurable_const measurable_const

theorem measurableSet_permEvent {n B : ℕ} (T : (Fin n → ℝ) → ℝ) (hT : Measurable T) (ε α : ℝ) :
    MeasurableSet (permEvent n B T ε α) :=
  measurableSet_le (measurable_from_prod_countable_left fun π => measurable_pval_left T hT ε π) measurable_const

local instance uniformDraws_finite (n B : ℕ) : IsFiniteMeasure (uniformDraws n B) :=
  ⟨by
    rw [uniformDraws, Measure.uniformOfFinset_apply]
    exact ENNReal.div_lt_top (by simp) (by simp)⟩

/-- At every factor `f`, the draws' probabilities of the event summed over the relabellings are at most `α n!`. -/
theorem sum_relabel_le {n B : ℕ} (T : (Fin n → ℝ) → ℝ) (ε α : ℝ) (hα : 0 ≤ α) (hε : 0 ≤ ε) (f : Fin n → ℝ) :
    ∑ τ : Equiv.Perm (Fin n), uniformDraws n B (Prod.mk (f ∘ ⇑τ) ⁻¹' permEvent n B T ε α) ≤
      ENNReal.ofReal α * (n.factorial : ℝ≥0∞) := by
  classical
  have h5a := perm_finite_core n B ℝ f T α ε hα hε
  set C := (univ.filter fun x : Equiv.Perm (Fin n) × (Fin B → Equiv.Perm (Fin n)) =>
    pval T ε (f ∘ ⇑x.1) x.2 ≤ α).card with hCdef
  have hC : ∑ τ : Equiv.Perm (Fin n), (univ.filter fun π : Fin B → Equiv.Perm (Fin n) =>
      π ∈ Prod.mk (f ∘ ⇑τ) ⁻¹' permEvent n B T ε α).card = C := by
    rw [hCdef, card_eq_sum_ones (univ.filter _), sum_filter, Fintype.sum_prod_type]
    refine sum_congr rfl fun τ _ => ?_
    rw [card_eq_sum_ones, sum_filter]
    refine sum_congr rfl fun π _ => ?_
    by_cases hp : pval T ε (f ∘ ⇑τ) π ≤ α
    · rw [ite_eq_left hp, ite_eq_left (by simpa [permEvent] using hp)]
    · rw [ite_eq_right hp, ite_eq_right (by simpa [permEvent] using hp)]
  have hM : ((univ : Finset (Fin B → Equiv.Perm (Fin n))).card : ℝ≥0∞) = (n.factorial : ℝ≥0∞) ^ B := by
    simp [Fintype.card_perm]
  have hCle : (C : ℝ≥0∞) ≤ ENNReal.ofReal α * (n.factorial : ℝ≥0∞) ^ (B + 1) := by
    rw [← ENNReal.ofReal_natCast C, ← ENNReal.ofReal_natCast n.factorial, ← ENNReal.ofReal_pow (by positivity),
      ← ENNReal.ofReal_mul hα]
    exact ENNReal.ofReal_le_ofReal h5a
  have h0 : ((n.factorial : ℝ≥0∞) ^ B) ≠ 0 := by simp [Nat.factorial_ne_zero]
  have htop : ((n.factorial : ℝ≥0∞) ^ B) ≠ ∞ := by simp
  simp only [uniformDraws, Measure.uniformOfFinset_apply, div_eq_mul_inv]
  rw [← sum_mul, ← Nat.cast_sum, hC, hM]
  calc (C : ℝ≥0∞) * ((n.factorial : ℝ≥0∞) ^ B)⁻¹
      ≤ ENNReal.ofReal α * (n.factorial : ℝ≥0∞) ^ (B + 1) * ((n.factorial : ℝ≥0∞) ^ B)⁻¹ := by gcongr
    _ = ENNReal.ofReal α * (n.factorial : ℝ≥0∞) * (n.factorial : ℝ≥0∞) ^ B * ((n.factorial : ℝ≥0∞) ^ B)⁻¹ := by
        ring
    _ = ENNReal.ofReal α * (n.factorial : ℝ≥0∞) := ENNReal.mul_inv_cancel_right h0 htop

/-- (5b) the measure-theoretic form. -/
theorem perm_measure : perm_measure_stmt := by
  intro n B μ T α ε hμ hT hinv hα hε
  change (μ.prod (uniformDraws n B)) (permEvent n B T ε α) ≤ ENNReal.ofReal α
  have hE := measurableSet_permEvent (B := B) T hT ε α
  rw [Measure.prod_apply hE]
  have hhm : Measurable fun f => uniformDraws n B (Prod.mk f ⁻¹' permEvent n B T ε α) :=
    measurable_measure_prodMk_left hE
  have hshift : ∀ τ : Equiv.Perm (Fin n),
      ∫⁻ f, uniformDraws n B (Prod.mk f ⁻¹' permEvent n B T ε α) ∂μ =
        ∫⁻ f, uniformDraws n B (Prod.mk (f ∘ ⇑τ) ⁻¹' permEvent n B T ε α) ∂μ := fun τ => by
    rw [← lintegral_map hhm (measurable_comp_perm τ), hinv τ]
  have hN0 : (n.factorial : ℝ≥0∞) ≠ 0 := by simp [Nat.factorial_ne_zero]
  have hNtop : (n.factorial : ℝ≥0∞) ≠ ∞ := by simp
  have hsum : (n.factorial : ℝ≥0∞) * ∫⁻ f, uniformDraws n B (Prod.mk f ⁻¹' permEvent n B T ε α) ∂μ ≤
      (n.factorial : ℝ≥0∞) * ENNReal.ofReal α := by
    calc (n.factorial : ℝ≥0∞) * ∫⁻ f, uniformDraws n B (Prod.mk f ⁻¹' permEvent n B T ε α) ∂μ
        = ∑ _τ : Equiv.Perm (Fin n), ∫⁻ f, uniformDraws n B (Prod.mk f ⁻¹' permEvent n B T ε α) ∂μ := by
          simp [Fintype.card_perm]
      _ = ∑ τ : Equiv.Perm (Fin n), ∫⁻ f, uniformDraws n B (Prod.mk (f ∘ ⇑τ) ⁻¹' permEvent n B T ε α) ∂μ :=
          sum_congr rfl fun τ _ => hshift τ
      _ = ∫⁻ f, ∑ τ : Equiv.Perm (Fin n), uniformDraws n B (Prod.mk (f ∘ ⇑τ) ⁻¹' permEvent n B T ε α) ∂μ :=
          (lintegral_finsetSum _ fun τ _ => hhm.comp (measurable_comp_perm τ)).symm
      _ ≤ ∫⁻ _f, ENNReal.ofReal α * (n.factorial : ℝ≥0∞) ∂μ :=
          lintegral_mono fun f => sum_relabel_le T ε α hα hε f
      _ = (n.factorial : ℝ≥0∞) * ENNReal.ofReal α := by
          rw [lintegral_const, measure_univ, mul_one, mul_comm]
  exact (ENNReal.mul_le_mul_iff_right hN0 hNtop).1 hsum

/-- Witness for 5b: the law putting 1/2 on each of two constant vectors (invariant under every relabelling). -/
theorem perm_measure_witness : perm_measure_witness_stmt := by
  let μ : Measure (Fin 4 → ℝ) :=
    (2⁻¹ : ℝ≥0∞) • Measure.dirac (fun _ => 0) + (2⁻¹ : ℝ≥0∞) • Measure.dirac (fun _ => 1)
  refine ⟨μ, fun f => f 0, ⟨?_⟩, measurable_pi_apply 0, fun σ => ?_, fun _ => 0, fun _ => 1, ?_, ?_, ?_⟩
  · simp [μ, ENNReal.inv_two_add_inv_two]
  · have hc := measurable_comp_perm σ
    simp only [μ, Measure.map_add _ _ hc, Measure.map_smul _ hc.aemeasurable, Measure.map_dirac' hc]
    rfl
  · intro h
    have := congrFun h 0
    norm_num at this
  · simp [μ]
  · simp [μ]

end MeasureForm

/-! Witnesses -/

/-- Witness for 5a: `g = (0, 1)`, `T h = -h 0`, the swap as the one draw: it is not counted, so `p = 1/2 ≤ α`. -/
theorem perm_witness : perm_witness_stmt := by
  refine ⟨fun i => (i : ℝ), fun h => -h 0, 1 / 2, 1 / 10 ^ 9, by norm_num, by norm_num, rfl, ?_⟩
  refine card_pos.2 ⟨(1, fun _ => Equiv.swap 0 1), mem_filter.2 ⟨mem_univ _, ?_⟩⟩
  have hempty : ∀ b : Fin 1, ¬ ((-(((fun i : Fin 2 => (i : ℝ)) ∘ ⇑(1 : Equiv.Perm (Fin 2))) ∘
      ⇑((fun _ : Fin 1 => Equiv.swap (0 : Fin 2) 1) b)) 0 ≥
        -((fun i : Fin 2 => (i : ℝ)) ∘ ⇑(1 : Equiv.Perm (Fin 2))) 0 - 1 / 10 ^ 9) ∨
      ((fun i : Fin 2 => (i : ℝ)) ∘ ⇑(1 : Equiv.Perm (Fin 2))) ∘ ⇑((fun _ : Fin 1 => Equiv.swap (0 : Fin 2) 1) b) =
        (fun i : Fin 2 => (i : ℝ)) ∘ ⇑(1 : Equiv.Perm (Fin 2))) := by
    intro b h
    rcases h with h | h
    · simp at h
      norm_num at h
    · have := congrFun h 0
      simp at this
  unfold pval
  rw [filter_false_of_mem fun b _ => hempty b]
  norm_num

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
