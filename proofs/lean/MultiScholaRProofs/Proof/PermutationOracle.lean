module

public import MultiScholaRProofs.Spec.Permutation
public import MultiScholaRProofs.Oracle.Permutation
import Mathlib.Tactic

/-!
# The target 5 oracle equals the Spec's counts (KIT-193 section 5.5)

`Oracle/Permutation.lean` writes `perm-floor.tsv` with its own counts over lists of rationals (`unchanged`,
`reflected`), because the Spec's `unchangedDraws` and `reflectedDraws` count with classical decidability and cannot
be evaluated. Ticket 5.5 allows a separate executable oracle only with a theorem that it equals the Spec definition.
This file proves it:

- `unchanged_eq_unchangedDraws`: for every factor `f : Fin n → ℚ` and draws `π : Fin B → Perm (Fin n)`, the
  oracle's `unchanged`, on `f` and the draws written as lists, is the Spec's `unchangedDraws f π`.
- `reflected_eq_reflectedDraws`: the same for `reflected` and the Spec's `reflectedDraws` at the real factor
  `i ↦ (f i : ℝ)` (the Spec's `reflect` is stated over `ℝ`; the oracle's mean is the rational one, and the cast
  carries one to the other exactly).
- `row_counts`: one table row, `f = g ∘ τ` and the two draws `π₁, π₂`, as the oracle's `row` builds it: `d`, `r` and
  `(1 + d) / 3` are the Spec's `unchangedDraws`, `reflectedDraws` and the floor of `perm_floor_stmt` at `B = 2`.
- `orderings_sound`, `orderings_complete`, `orderings_nodup`: the oracle's 24 orderings are exactly the encodings
  of the permutations of `Fin 4`, each once, so the table's `24³` rows per fixture are the cases `(τ, π₁, π₂)`.

The encoding (`encPerm`, `encDraws`) is the oracle's: a permutation `σ` is the list `σ 0, …, σ (n - 1)` of
0-based values, the draws are the list of their encodings, a factor is `List.ofFn f`. Not covered: that a fixture's
literal list is `List.ofFn` of its values (`List.ofFn_get`, by definition), the one-based labels and the text of
the rows (string formatting), and Lean's interpreter that evaluates the oracle (trusted for the tables, 9.1).
-/

public section

namespace MultiScholaRProofs.Permutation

open Finset MultiScholaRProofs.Oracle.Permutation

/-- A permutation as the oracle writes an ordering: the list `σ 0, …, σ (n - 1)` of 0-based values. -/
def encPerm {n : ℕ} (σ : Equiv.Perm (Fin n)) : List ℕ := List.ofFn fun i => (σ i : ℕ)

/-- The draws as the oracle takes them: the list of their orderings. -/
def encDraws {n B : ℕ} (π : Fin B → Equiv.Perm (Fin n)) : List (List ℕ) := List.ofFn fun b => encPerm (π b)

/-- The oracle's `compose` is composition: `f[σ] = f ∘ σ`. -/
theorem compose_ofFn {n : ℕ} (f : Fin n → ℚ) (σ : Equiv.Perm (Fin n)) :
    compose (List.ofFn f) (encPerm σ) = List.ofFn (f ∘ σ) := by
  simp [compose, encPerm, List.map_ofFn, Function.comp_def, List.getD_eq_getElem?_getD]

/-- The oracle's `reflectL` is the reflection about the rational mean. -/
theorem reflectL_ofFn {n : ℕ} (f : Fin n → ℚ) :
    reflectL (List.ofFn f) = List.ofFn fun i => 2 * ((∑ j, f j) / n) - f i := by
  simp [reflectL, List.map_ofFn, List.sum_ofFn, Function.comp_def]

/-- A list filter's length over `List.ofFn` is a `Finset` filter's card (any decidability instance). -/
theorem length_filter_ofFn {α : Type} {B : ℕ} (g : Fin B → α) (p : α → Bool) (P : Fin B → Prop)
    [DecidablePred P] (h : ∀ b, p (g b) = true ↔ P b) :
    ((List.ofFn g).filter p).length = (univ.filter P).card := by
  rw [List.ofFn_eq_map, List.filter_map, List.length_map]
  have : (univ.filter P) = ((List.finRange B).filter (p ∘ g)).toFinset := by
    ext b
    simp [h]
  rw [this, List.toFinset_card_of_nodup ((List.nodup_finRange B).filter _)]

/-- The reflection of the real factor is the real image of the rational reflection. -/
theorem reflect_cast {n : ℕ} (f : Fin n → ℚ) :
    reflect (fun i => (f i : ℝ)) = fun i => ((2 * ((∑ j, f j) / n) - f i : ℚ) : ℝ) := by
  funext i
  simp [reflect]

/-- `d`: the oracle's `unchanged` is the Spec's `unchangedDraws`. -/
theorem unchanged_eq_unchangedDraws {n B : ℕ} (f : Fin n → ℚ) (π : Fin B → Equiv.Perm (Fin n)) :
    unchanged (List.ofFn f) (encDraws π) = unchangedDraws f π := by
  have h := length_filter_ofFn (fun b => encPerm (π b))
    (fun σ => decide (compose (List.ofFn f) σ = List.ofFn f)) (fun b => f ∘ π b = f)
    (fun b => by simp [compose_ofFn, List.ofFn_inj])
  unfold unchanged encDraws unchangedDraws
  rw [h]
  exact congrArg Finset.card (Finset.filter_congr_decidable _ _ _).symm

/-- `r`: the oracle's `reflected` is the Spec's `reflectedDraws` at the real factor. -/
theorem reflected_eq_reflectedDraws {n B : ℕ} (f : Fin n → ℚ) (π : Fin B → Equiv.Perm (Fin n)) :
    reflected (List.ofFn f) (encDraws π) = reflectedDraws (fun i => (f i : ℝ)) π := by
  have hinj : Function.Injective (fun g : Fin n → ℚ => fun i => (g i : ℝ)) :=
    Rat.cast_injective.comp_left
  have h := length_filter_ofFn (fun b => encPerm (π b))
    (fun σ => decide (compose (List.ofFn f) σ = reflectL (List.ofFn f) ∧ compose (List.ofFn f) σ ≠ List.ofFn f))
    (fun b => (fun i => (f i : ℝ)) ∘ π b = reflect (fun i => (f i : ℝ)) ∧
      (fun i => (f i : ℝ)) ∘ π b ≠ fun i => (f i : ℝ))
    (fun b => by
      have e1 : ((fun i => (f i : ℝ)) ∘ π b = reflect (fun i => (f i : ℝ))) ↔
          (f ∘ π b = fun i => 2 * ((∑ j, f j) / n) - f i) := by
        rw [reflect_cast]
        exact hinj.eq_iff
      have e2 : ((fun i => (f i : ℝ)) ∘ π b = fun i => (f i : ℝ)) ↔ f ∘ π b = f := hinj.eq_iff
      simp only [e1, ne_eq, e2]
      simp [compose_ofFn, reflectL_ofFn, List.ofFn_inj])
  unfold reflected encDraws reflectedDraws
  -- the Spec's instance and the synthesized one agree here (equality on `ℝ` is decided classically in both)
  rw [h]

/-- One row of `perm-floor.tsv` as `row` computes it, `f = g ∘ τ` and the draws `π₁, π₂`: its `d`, its `r` and its
`(1 + d) / 3` are the Spec's `unchangedDraws`, `reflectedDraws` and the floor of `perm_floor_stmt` at `B = 2`. -/
theorem row_counts {n : ℕ} (g : Fin n → ℚ) (τ π₁ π₂ : Equiv.Perm (Fin n)) :
    unchanged (compose (List.ofFn g) (encPerm τ)) [encPerm π₁, encPerm π₂] = unchangedDraws (g ∘ τ) ![π₁, π₂] ∧
    reflected (compose (List.ofFn g) (encPerm τ)) [encPerm π₁, encPerm π₂] =
      reflectedDraws (fun i => ((g ∘ τ) i : ℝ)) ![π₁, π₂] ∧
    (((1 + unchanged (compose (List.ofFn g) (encPerm τ)) [encPerm π₁, encPerm π₂] : ℚ) / 3 : ℚ) : ℝ) =
      (1 + (unchangedDraws (g ∘ τ) ![π₁, π₂] : ℝ)) / (1 + (2 : ℕ)) := by
  have hd : [encPerm π₁, encPerm π₂] = encDraws ![π₁, π₂] := by
    simp [encDraws, List.ofFn_succ]
  rw [hd, compose_ofFn, unchanged_eq_unchangedDraws, reflected_eq_reflectedDraws]
  refine ⟨rfl, rfl, ?_⟩
  push_cast
  norm_num

/-- A list of length `n` without repeats and with values below `n` is the encoding of a permutation. -/
theorem exists_encPerm {n : ℕ} (l : List ℕ) (hlen : l.length = n) (hnd : l.Nodup) (hlt : ∀ x ∈ l, x < n) :
    ∃ σ : Equiv.Perm (Fin n), encPerm σ = l := by
  subst hlen
  let f : Fin l.length → Fin l.length := fun i => ⟨l.get i, hlt _ (List.get_mem l i)⟩
  have hf : Function.Injective f := by
    intro i j hij
    exact (List.Nodup.get_inj_iff hnd).1 (congrArg Fin.val hij)
  refine ⟨Equiv.ofBijective f (Finite.injective_iff_bijective.1 hf), ?_⟩
  unfold encPerm
  conv_rhs => rw [← List.ofFn_get l]
  rfl

/-- Every ordering the oracle lists is the encoding of a permutation of `Fin 4`. -/
theorem orderings_sound : ∀ l ∈ orderings, ∃ σ : Equiv.Perm (Fin 4), encPerm σ = l := by
  intro l hl
  have hp : l.Perm [0, 1, 2, 3] := List.mem_permutations.1 hl
  refine exists_encPerm l hp.length_eq (hp.nodup_iff.2 (by decide)) fun x hx => ?_
  have : x ∈ [0, 1, 2, 3] := hp.subset hx
  simp at this
  omega

/-- Every permutation of `Fin 4` is among the oracle's orderings. -/
theorem orderings_complete (σ : Equiv.Perm (Fin 4)) : encPerm σ ∈ orderings := by
  unfold orderings encPerm
  rw [List.mem_permutations]
  refine (List.perm_ext_iff_of_nodup ?_ (by decide)).2 fun x => ?_
  · exact List.nodup_ofFn.2 (Fin.val_injective.comp σ.injective)
  · simp only [List.mem_ofFn, List.mem_cons, List.not_mem_nil, or_false]
    constructor
    · rintro ⟨i, rfl⟩
      have := (σ i).isLt
      omega
    · intro hx
      refine ⟨σ.symm ⟨x, by omega⟩, ?_⟩
      simp

/-- The oracle's orderings are listed once each (so, with the two lemmas above, exactly the 24 permutations). -/
theorem orderings_nodup : orderings.Nodup ∧ orderings.length = 24 :=
  ⟨List.nodup_permutations _ (by decide), by simp [orderings, List.length_permutations, Nat.factorial]⟩

end MultiScholaRProofs.Permutation
