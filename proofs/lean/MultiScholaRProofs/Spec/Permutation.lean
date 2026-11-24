module

public import MultiScholaRProofs.Spec.QvalueBh
public import Mathlib.Data.Fintype.Perm
public import Mathlib.Algebra.Order.Floor.Semiring
public import Mathlib.MeasureTheory.Measure.Prod
public import Mathlib.Probability.Distributions.Uniform

/-!
# Target 5: the permutation p-value (KIT-193 section 6.5)

Source: KIT-191 r1.3 section 5.2 steps 4, 5, 7, 8 (step 5's tie tolerance and step 8's corrected floor, from this
ticket's first audit, U2 and U1); the implementation on `kit-191-enrichment-test-options` at 2f54469e,
`integrationEnrichmentPermutation()`, `R/utils_integration_enrichment.R:1590` onward; `:1631-1656`: the draws
`perms[, b]`, the permuted factor `f[perms[, b]]` (that is `f ∘ π_b`), the count
`abs(t*) >= abs(t_obs) - INTEGRATION_ENRICHMENT_PERMUTATION_TIE` (the constant 1e-9, `:857`), draws that leave `f`
unchanged counted without comparing (`identical(f[perms[, b]], f)`), `pvalue = (1 + exceed) / (1 + permutations)`;
VAL-97 (1), (2): B = 1,000 by default or 10,000; no covariates. The floor sentence, `R/utils_integration_gui.R:820-822`.

Model. The statistic `T` is any function (the two-sided `|t|` of Fisher-transformed correlations is one); the
observed and permuted statistics are in exact arithmetic; the tolerance `ε ≥ 0` is a parameter (the code's is
`10⁻⁹`, `ε = 0` is Phipson and Smyth's count); draws are with replacement, `B` of them, independent and uniform on
permutations and independent of `f` (the random-number generator's uniformity is an idealisation). Exchangeability:
the law of `f` is invariant under permuting the samples; it holds for a factor that does not depend on the data `y`;
a factor computed from `y` is, given `y`, a fixed vector, and for it the theorem describes the reference distribution
the test uses (its error rates are measured, KIT-183, KIT-191). That the rounding error of `|t*| - |t|` stays below
`ε` is assumed, not proved (9.1).
-/

@[expose] public section

namespace MultiScholaRProofs.Permutation

open Finset MeasureTheory MultiScholaRProofs.QvalueBh

open scoped Classical in
/-- The permutation p-value with tolerance `ε`: `(1 + #{b : T(f ∘ π_b) ≥ T(f) - ε or f ∘ π_b = f}) / (1 + B)`. -/
noncomputable def pval {n B : ℕ} {X : Type} (T : (Fin n → X) → ℝ) (ε : ℝ) (f : Fin n → X)
    (π : Fin B → Equiv.Perm (Fin n)) : ℝ :=
  (1 + ((univ.filter fun b => T (f ∘ π b) ≥ T f - ε ∨ f ∘ π b = f).card : ℝ)) / (1 + B)

open scoped Classical in
/-- `d`, the number of draws that leave `f` unchanged (`identity_draws`). -/
noncomputable def unchangedDraws {n B : ℕ} {X : Type} (f : Fin n → X) (π : Fin B → Equiv.Perm (Fin n)) : ℕ :=
  (univ.filter fun b => f ∘ π b = f).card

/-- The reflection of a real factor about its mean, `2 mean(f) - f` (fixture C of the paired test). At `n = 0`
the mean is Lean's `0 / 0 = 0` and the vector is empty; the definition is used with `0 < n`. -/
noncomputable def reflect {n : ℕ} (f : Fin n → ℝ) : Fin n → ℝ :=
  fun i => 2 * ((∑ j, f j) / n) - f i

open scoped Classical in
/-- `r`, the number of draws that map `f` to its reflection about its mean. -/
noncomputable def reflectedDraws {n B : ℕ} (f : Fin n → ℝ) (π : Fin B → Equiv.Perm (Fin n)) : ℕ :=
  (univ.filter fun b => f ∘ π b = reflect f ∧ f ∘ π b ≠ f).card

/-- The discrete measurable structure on permutations of `Fin n` (every set measurable). -/
scoped instance permMeasurableSpace (n : ℕ) : MeasurableSpace (Equiv.Perm (Fin n)) := ⊤

/-- The law of the `B` draws: uniform on `(Fin B → Perm (Fin n))`, i.e. independent and uniform. -/
noncomputable def uniformDraws (n B : ℕ) : Measure (Fin B → Equiv.Perm (Fin n)) :=
  Measure.uniformOfFinset univ

/-- (5a) Finite core: for every `g`, `T`, `α ≥ 0` and `ε ≥ 0`, among the `(n!)^(B+1)` pairs (relabelling `τ`, draws
`π`), at most `α (n!)^(B+1)` give `p_ε(g ∘ τ, π) ≤ α`. Hypotheses: `0 ≤ α` (at `α < 0` the left side is 0 and the
right side negative); `0 ≤ ε` (the code's tolerance). -/
def perm_finite_core_stmt : Prop :=
  ∀ (n B : ℕ) (X : Type) (g : Fin n → X) (T : (Fin n → X) → ℝ) (α ε : ℝ), 0 ≤ α → 0 ≤ ε →
    ((univ.filter fun x : Equiv.Perm (Fin n) × (Fin B → Equiv.Perm (Fin n)) =>
        pval T ε (g ∘ x.1) x.2 ≤ α).card : ℝ) ≤ α * (Nat.factorial n : ℝ) ^ (B + 1)

/-- (5a') The identity clause is the exact count: for every `ε ≥ 0`, a draw that leaves `f` unchanged is counted by
the comparison itself in exact arithmetic. Hypothesis `0 ≤ ε`. -/
def identity_clause_stmt : Prop :=
  ∀ (n : ℕ) (X : Type) (T : (Fin n → X) → ℝ) (f : Fin n → X) (σ : Equiv.Perm (Fin n)) (ε : ℝ), 0 ≤ ε →
    ((T (f ∘ σ) ≥ T f - ε ∨ f ∘ σ = f) ↔ T (f ∘ σ) ≥ T f - ε)

/-- The tolerance can only raise the p-value: `p_0 ≤ p_ε` for `ε ≥ 0` (C01). Hypothesis `0 ≤ ε`. -/
def perm_tolerance_mono_stmt : Prop :=
  ∀ (n B : ℕ) (X : Type) (T : (Fin n → X) → ℝ) (f : Fin n → X) (π : Fin B → Equiv.Perm (Fin n)) (ε : ℝ),
    0 ≤ ε → pval T 0 f π ≤ pval T ε f π

/-- (5b) Measure-theoretic form (D3): for a probability law `μ` of the factor invariant under every permutation of
the samples, a measurable `T`, and draws independent and uniform, `P(p_ε ≤ α) ≤ α`. Hypotheses: `μ` a probability
measure; `T` measurable; `μ.map (· ∘ σ) = μ` for every `σ` (exchangeability); `0 ≤ α`; `0 ≤ ε`. -/
def perm_measure_stmt : Prop :=
  ∀ (n B : ℕ) (μ : Measure (Fin n → ℝ)) (T : (Fin n → ℝ) → ℝ) (α ε : ℝ),
    IsProbabilityMeasure μ → Measurable T → (∀ σ : Equiv.Perm (Fin n), μ.map (fun f => f ∘ σ) = μ) →
    0 ≤ α → 0 ≤ ε →
      (μ.prod (uniformDraws n B)) {x | pval T ε x.1 x.2 ≤ α} ≤ ENNReal.ofReal α

/-- (5c) The floor: `p ≥ (1 + d) / (1 + B)`, `d` the draws that leave `f` unchanged (`min_p`). No hypothesis. -/
def perm_floor_stmt : Prop :=
  ∀ (n B : ℕ) (X : Type) (T : (Fin n → X) → ℝ) (ε : ℝ) (f : Fin n → X) (π : Fin B → Equiv.Perm (Fin n)),
    (1 + (unchangedDraws f π : ℝ)) / (1 + B) ≤ pval T ε f π

/-- (5c, fixture C) For a statistic invariant under reflecting the factor about its mean (`|t|` is: the correlations
change sign), every draw that reflects `f` is counted too: `p ≥ (1 + d + r) / (1 + B)`. Hypotheses: `0 < n` (the
mean is a real quotient); `T (reflect f) = T f` for every `f`; `0 ≤ ε`. -/
def perm_floor_reflect_stmt : Prop :=
  ∀ (n B : ℕ) (T : (Fin n → ℝ) → ℝ) (ε : ℝ) (f : Fin n → ℝ) (π : Fin B → Equiv.Perm (Fin n)),
    0 < n → (∀ h : Fin n → ℝ, T (reflect h) = T h) → 0 ≤ ε →
      (1 + (unchangedDraws f π : ℝ) + (reflectedDraws f π : ℝ)) / (1 + B) ≤ pval T ε f π

/-- (5c) The BH floor, "adjusted p at most q": if every p-value is at least `p_min > 0` and at most 1, and BH at
level `q > 0` calls at least one term, it calls at least `⌈m p_min / q⌉`. Hypotheses: `0 < q`, `0 < p_min` (the
division is real; the floor is positive); `p_min ≤ p_i ≤ 1`; one call. -/
def bh_floor_stmt : Prop :=
  ∀ (m : ℕ) (p : Fin m → ℝ) (q pmin : ℝ), 0 < q → 0 < pmin →
    (∀ i, pmin ≤ p i) → (∀ i, p i ≤ 1) → (∃ i, bh p i ≤ q) →
      (Nat.ceil (m * pmin / q) : ℝ) ≤ ((univ.filter fun i => bh p i ≤ q).card : ℝ)

/-- (5c) The BH floor for the route's comparison, "adjusted p below q". Hypotheses as `bh_floor_stmt`. -/
def bh_floor_strict_stmt : Prop :=
  ∀ (m : ℕ) (p : Fin m → ℝ) (q pmin : ℝ), 0 < q → 0 < pmin →
    (∀ i, pmin ≤ p i) → (∀ i, p i ≤ 1) → (∃ i, bh p i < q) →
      (Nat.ceil (m * pmin / q) : ℝ) ≤ ((univ.filter fun i => bh p i < q).card : ℝ)

/-- (5c') The floor sentence of KIT-191 r1.2 (5.2 step 8), "at least k_min sets must have p ≤ q k_min / m", was
false: `m = 100`, `p_min = 1/1001`, `q = 1/20`; four terms at `2/1001` and 96 at 1 are called (adjusted p
`50/1001`), `k_min = 2`, and no term has `p ≤ 1/1000`. Exact rationals; closed statement. -/
def bh_floor_sentence_counterexample_stmt : Prop :=
  let p : Fin 100 → ℚ := fun i => if (i : ℕ) < 4 then 2 / 1001 else 1
  (∀ i, (1 : ℚ) / 1001 ≤ p i) ∧ (∃ i, bh p i < 1 / 20) ∧
    Nat.ceil ((100 : ℚ) * (1 / 1001) / (1 / 20)) = 2 ∧
    (univ.filter fun i => p i ≤ (1 / 20 : ℚ) * 2 / 100).card = 0

/-! ## Witnesses (section 5.1) -/

open scoped Classical in
/-- Witness for `perm_finite_core_stmt`, `identity_clause_stmt`, `perm_tolerance_mono_stmt`: two samples, one draw,
the code's tolerance, `0 < α < 1`, and a pair `(τ, π)` that the bound counts. -/
def perm_witness_stmt : Prop :=
  ∃ (g : Fin 2 → ℝ) (T : (Fin 2 → ℝ) → ℝ) (α ε : ℝ), 0 < α ∧ α < 1 ∧ ε = 1 / 10 ^ 9 ∧
    0 < (univ.filter fun x : Equiv.Perm (Fin 2) × (Fin 1 → Equiv.Perm (Fin 2)) =>
      pval T ε (g ∘ x.1) x.2 ≤ α).card

/-- Witness for `perm_measure_stmt`: a permutation-invariant probability law on `Fin 4 → ℝ` that is not a point mass,
with a measurable statistic. -/
def perm_measure_witness_stmt : Prop :=
  ∃ (μ : Measure (Fin 4 → ℝ)) (T : (Fin 4 → ℝ) → ℝ),
    IsProbabilityMeasure μ ∧ Measurable T ∧ (∀ σ : Equiv.Perm (Fin 4), μ.map (fun f => f ∘ σ) = μ) ∧
    ∃ f f' : Fin 4 → ℝ, f ≠ f' ∧ 0 < μ {f} ∧ 0 < μ {f'}

/-- Witness for `perm_floor_reflect_stmt`: a reflection-invariant statistic (`|f 0 - f 1|`-like) on 4 samples. -/
def perm_floor_reflect_witness_stmt : Prop :=
  ∃ T : (Fin 4 → ℝ) → ℝ, (∀ h : Fin 4 → ℝ, T (reflect h) = T h) ∧ ∃ h h' : Fin 4 → ℝ, T h ≠ T h'

/-- Witness for `bh_floor_stmt` and `bh_floor_strict_stmt`: their hypotheses hold with a floor of at least 2. -/
def bh_floor_witness_stmt : Prop :=
  ∃ (m : ℕ) (p : Fin m → ℝ) (q pmin : ℝ), 0 < q ∧ 0 < pmin ∧ (∀ i, pmin ≤ p i) ∧ (∀ i, p i ≤ 1) ∧
    (∃ i, bh p i < q) ∧ 2 ≤ Nat.ceil (m * pmin / q)

end MultiScholaRProofs.Permutation
