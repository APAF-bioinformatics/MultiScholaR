module

public import Mathlib.Algebra.Order.Field.Rat
public import Mathlib.Basic.Real.Basic
public import Mathlib.Algebra.BigOperators.Fin
public import Mathlib.Data.Fintype.Pi
public import Mathlib.Data.List.Sort

/-!
# Target 4: q-value against Benjamini-Hochberg (KIT-193 section 6.4)

Source: KIT-189 r1.2 section 2.1 (confirmed items 2 and 4, as replaced by AUDIT-KIT-189 A02); `qvalue` as KIT-189
2.1 reads it (`qvalue/R/qvalue.R:111-119`: `pi0 * pmin(1, cummin(p[o] * m / i))[ro]`); `stats::p.adjust(method =
"BH")` (`i <- lp:1L; o <- order(p, decreasing = TRUE); ro <- order(o); pmin(1, cummin(n/i * p[o]))[ro]`); the
route's flag `significant = qvalue < q_cutoff`, `R/func_general_enrichment_orchestration.R:544`.

Model. `bh` is stated once over any linearly ordered field `K` (ℝ in the theorems 4a-4c, ℚ in 4d and in the
oracle). `bhSorted` transcribes `p.adjust`'s code (decreasing order, cumulative minimum of `m/i * p`, capped
at 1, back to the original order); `bh` is its closed form, `min(1, min over j with p_j ≥ p_i of m p_j / #{l : p_l ≤
p_j})`, which the theorems use; `bh_sorted_eq_stmt` states that they are equal. Real arithmetic: `q = pi0 * BH`
holds in R only up to rounding (`p[o] * m / i` against `n / i * p[o]`); `pi0` is a given number (how `qvalue`
estimates it, and its fallbacks, are outside the model; KIT-189 measures them); in 4d the p-values' law is the stated
mixture, not the hypergeometric p-values of a real list. A missing q-value is `none`; R's `NA < alpha` is `NA`, which
flags nothing.
-/

@[expose] public section

namespace MultiScholaRProofs.QvalueBh

open Finset

variable {K : Type} [Field K] [LinearOrder K] [IsStrictOrderedRing K]

/-- The prefix minimum of a list (R's `cummin`). -/
def cumMin : List K → List K
  | [] => []
  | a :: l => a :: (cumMin l).map (min a)

/-- `p.adjust(p, "BH")` as R computes it, for `m` p-values without `NA`: order decreasing (`List.mergeSort` with
"≥"; the order R gives tied values is not modelled and does not matter: `bh_sorted_eq_stmt` equates this definition,
whatever order ties take here, with the closed form, which does not depend on it), scale the value at 0-based position `r` by `m / (m - r)` (R's `n / i`,
`i = m - r ≥ 1`), take the cumulative minimum, cap at 1, return to the original order. The `getD` default is never
used (every index is in the order). -/
def bhSorted {m : ℕ} (p : Fin m → K) : Fin m → K :=
  let o : List (Fin m) := (List.finRange m).mergeSort fun a b => decide (p b ≤ p a)
  let scaled : List K := o.zipIdx.map fun x => (m : K) / ((m - x.2 : ℕ) : K) * p x.1
  let cm : List K := cumMin scaled
  fun i => min 1 (cm.getD (o.idxOf i) 0)

/-- `#{l : p_l ≤ v}`. -/
def cntLe {m : ℕ} (p : Fin m → K) (v : K) : ℕ := (univ.filter fun l => p l ≤ v).card

/-- `m p_j / #{l : p_l ≤ p_j}`; the count is at least 1 (it counts `j`), so the division is the real one. -/
def bhTerm {m : ℕ} (p : Fin m → K) (j : Fin m) : K := (m : K) * p j / (cntLe p (p j) : K)

/-- The BH adjusted p-value in closed form: `min(1, min over j with p_i ≤ p_j of bhTerm p j)`. The fold starts from
`bhTerm p i`, a member of the set (it contains `i`), so it is the minimum over a non-empty set. -/
def bh {m : ℕ} (p : Fin m → K) (i : Fin m) : K :=
  min 1 ((univ.filter fun j => p i ≤ p j).fold min (bhTerm p i) (bhTerm p))

/-- The q-value with a given `pi0` (`qvalue::qvalue(p, pi0 = pi0)`, KIT-189 2.1). -/
def qv {m : ℕ} (π0 : K) (p : Fin m → K) (i : Fin m) : K := π0 * bh p i

/-- The route's flag on an optional value: `TRUE` exactly when the value is present and below `α` (`NA < α` flags
nothing). -/
def flag (q : Option K) (α : K) : Prop := ∃ v, q = some v ∧ v < α

/-- `bhSorted` (R's code) equals the closed form `bh`. No hypothesis. Specification lemma. -/
def bh_sorted_eq_stmt : Prop :=
  ∀ (K : Type) [Field K] [LinearOrder K] [IsStrictOrderedRing K] (m : ℕ) (p : Fin m → K), bhSorted p = bh p

/-- (4a) For `pi0 ≤ 1` and p-values `≥ 0`: `q ≤ bh` at every term, hence the q-value flag contains the BH flag.
Hypotheses: `π0 ≤ 1` (KIT-189's `pi0` is at most 1); `0 ≤ p j` for all `j` (p-values). The ticket's `0 ≤ π0` is not
needed for this conclusion and is not assumed (a weaker hypothesis, so a stronger statement). -/
def qflag_contains_bh_stmt : Prop :=
  ∀ (m : ℕ) (p : Fin m → ℝ) (π0 : ℝ), π0 ≤ 1 → (∀ j, 0 ≤ p j) →
    (∀ i, qv π0 p i ≤ bh p i) ∧ (∀ (i : Fin m) (α : ℝ), bh p i < α → flag (some (qv π0 p i)) α)

/-- (4a) A missing q-value flags nothing. No hypothesis. -/
def flag_none_stmt : Prop := ∀ α : ℝ, ¬ flag (none : Option ℝ) α

/-- (4b) For `pi0 > 0`, the q-value flags a term BH does not exactly when `α ≤ bh < α / pi0`. Hypothesis
`0 < π0` (the division is the real one). The ticket's `π0 ≤ 1` is not needed for this equivalence and is not
assumed. -/
def qflag_extra_stmt : Prop :=
  ∀ (m : ℕ) (p : Fin m → ℝ) (π0 α : ℝ) (i : Fin m), 0 < π0 →
    (flag (some (qv π0 p i)) α ∧ ¬ bh p i < α ↔ α ≤ bh p i ∧ bh p i < α / π0)

/-- The p-values lie in `[0, 1]`. -/
def UnitP {m : ℕ} (p : Fin m → K) : Prop := ∀ j, 0 ≤ p j ∧ p j ≤ 1

/-- (4c) For p-values in `[0, 1]`: `p_i ≤ bh p i`. Hypothesis `UnitP p` (`p_i ≤ 1` is needed: `bh` is capped at 1). -/
def bh_ge_p_stmt : Prop :=
  ∀ (m : ℕ) (p : Fin m → ℝ) (i : Fin m), UnitP p → p i ≤ bh p i

/-- (4c) Monotone in the p-values. Hypotheses `UnitP p`, `UnitP p'`, `p ≤ p'` pointwise. -/
def bh_mono_stmt : Prop :=
  ∀ (m : ℕ) (p p' : Fin m → ℝ) (i : Fin m), UnitP p → UnitP p' → (∀ j, p j ≤ p' j) → bh p i ≤ bh p' i

/-- (4c) Ties get equal adjusted values. Hypothesis `UnitP p` (the range of the ticket's 4c; not needed for this
item, kept as the ticket states 4c's range once for all four facts). -/
def bh_ties_stmt : Prop :=
  ∀ (m : ℕ) (p : Fin m → ℝ) (i j : Fin m), UnitP p → p i = p j → bh p i = bh p j

/-- (4c) The step-up characterisation: `bh p i ≤ α` iff for some `k ≥ 1`, `p_i ≤ k α / m` and at least `k` p-values
are `≤ k α / m`. Hypothesis `UnitP p`. The division by `m` is real: the term `i` exists, so `m ≥ 1`. -/
def bh_stepup_stmt : Prop :=
  ∀ (m : ℕ) (p : Fin m → ℝ) (i : Fin m) (α : ℝ), UnitP p →
    (bh p i ≤ α ↔ ∃ k : ℕ, 1 ≤ k ∧ p i ≤ k * α / m ∧ k ≤ cntLe p (k * α / m))

/-! ## 4d-i: exact finite counterexamples

Two independent p-values, each `1` with probability 1/2 and otherwise uniform on `{1/(2N), …, N/(2N)}`. Outcome
`j : Fin (N + 1)`: `j = 0` is `p = 1`, `j ≥ 1` is `p = j / (2N)`. Under the complete null every call is false, so the
false discovery rate is the probability of at least one call. -/

/-- The p-value of outcome `j`. -/
def gridVal (N : ℕ) (j : Fin (N + 1)) : ℚ := if (j : ℕ) = 0 then 1 else (j : ℚ) / (2 * N)

/-- The probability of outcome `j` (used with `0 < N`). -/
def gridWt (N : ℕ) (j : Fin (N + 1)) : ℚ := if (j : ℕ) = 0 then 1 / 2 else 1 / (2 * N)

/-- The number of terms BH calls at level `q`: adjusted p `≤ q` (BH's standard rule) or, with `strict`, `< q` (the
route's flag, `:544`). -/
def calls (strict : Bool) (q : ℚ) {m : ℕ} (p : Fin m → ℚ) : ℕ :=
  (univ.filter fun i => if strict then bh p i < q else bh p i ≤ q).card

/-- BH over the terms with `p < 1` (the kept family). -/
def callsKept (strict : Bool) (q : ℚ) (l : List ℚ) : ℕ :=
  calls strict q fun i : Fin (l.filter fun x => decide (x < 1)).length => (l.filter fun x => decide (x < 1)).get i

/-- BH over the whole family. -/
def callsWhole (strict : Bool) (q : ℚ) (l : List ℚ) : ℕ :=
  calls strict q fun i : Fin l.length => l.get i

/-- The false discovery rate under the complete null of two independent grid p-values: the probability of at least
one call, BH over the kept family (`kept = true`) or over both terms. -/
def fdr2 (N : ℕ) (strict : Bool) (q : ℚ) (keptOnly : Bool) : ℚ :=
  ∑ ω : Fin (N + 1) × Fin (N + 1), gridWt N ω.1 * gridWt N ω.2 *
    (let l := [gridVal N ω.1, gridVal N ω.2]
     if 0 < (if keptOnly then callsKept strict q l else callsWhole strict q l) then 1 else 0)

/-- The grid law is a probability: its weights sum to 1. Hypothesis `0 < N`. -/
def grid_mass_stmt : Prop := ∀ N : ℕ, 0 < N → ∑ j : Fin (N + 1), gridWt N j = 1

/-- Each grid p-value is super-uniform: `P(p ≤ t) ≤ t` for every `t` in `[0, 1]`. Hypotheses `0 < N`,
`0 ≤ t ≤ 1`. -/
def grid_superuniform_stmt : Prop :=
  ∀ (N : ℕ) (t : ℚ), 0 < N → 0 ≤ t → t ≤ 1 → ∑ j ∈ univ.filter (fun j => gridVal N j ≤ t), gridWt N j ≤ t

/-- (4d-i) BH's standard rule (`≤ q`), `N = 20`, `q = 1/20`: false discovery rate exactly `3/40` over the terms with
`p < 1`, exactly `1/20` over both. Closed statement. -/
def fdr_n20_stmt : Prop := fdr2 20 false (1 / 20) true = 3 / 40 ∧ fdr2 20 false (1 / 20) false = 1 / 20

/-- (4d-i) The parenthesis of 6.4: on the `N = 20` grid the strict rule over the kept family gives `41/1600`.
Closed statement. -/
def fdr_n20_strict_stmt : Prop := fdr2 20 true (1 / 20) true = 41 / 1600

/-- (4d-i) `N = 25`, `q = 1/20` (neither `q` nor `q/2` is a grid point): exactly `3/50` over the terms with `p < 1`
against `1/25` over both, for either comparison. Closed statement. -/
def fdr_n25_stmt : Prop :=
  ∀ strict : Bool, fdr2 25 strict (1 / 20) true = 3 / 50 ∧ fdr2 25 strict (1 / 20) false = 1 / 25

/-! ## 4d-ii: the kept-family formula on a fine grid (stretch, D2) -/

/-- The p-value of outcome `j` under the law "1 with probability `1 - a`, else uniform on `{a j / N}`". -/
def mixVal (a : ℚ) (N : ℕ) (j : Fin (N + 1)) : ℚ := if (j : ℕ) = 0 then 1 else a * j / N

/-- The probability of outcome `j` (used with `0 < N`, `0 < a ≤ 1`). -/
def mixWt (a : ℚ) (N : ℕ) (j : Fin (N + 1)) : ℚ := if (j : ℕ) = 0 then 1 - a else a / N

/-- The false discovery rate under the complete null of `m` independent p-values with the mixture law, BH (standard
rule unless `strict`) over the terms with `p < 1`: the probability of at least one call. Independence is the product
of the weights. -/
def fdrKeptMix (m N : ℕ) (a q : ℚ) (strict : Bool) : ℚ :=
  ∑ ω : Fin m → Fin (N + 1), (∏ b, mixWt a N (ω b)) *
    (if 0 < callsKept strict q (List.ofFn fun b => mixVal a N (ω b)) then 1 else 0)

/-- Every BH threshold, on the scale of `p / a`, that is at most 1 is a grid point: `(q/a)(i/K) = j/N` for some
whole `j`, for all `1 ≤ i ≤ K ≤ m`. -/
def GridFine (m N : ℕ) (a q : ℚ) : Prop :=
  ∀ K i : ℕ, 1 ≤ i → i ≤ K → K ≤ m → q / a * ((i : ℚ) / K) ≤ 1 → ∃ j : ℕ, q / a * ((i : ℚ) / K) = (j : ℚ) / N

/-- The least common multiple of `1, …, n`. -/
def lcmUpTo (n : ℕ) : ℕ := (List.range' 1 n).foldl Nat.lcm 1

/-- (4d-ii) On a grid fine enough that every threshold is a grid point, BH (at most `q`) over the terms with `p < 1`
has false discovery rate exactly `min(q/a, 1) (1 - (1 - a)^m)` under the complete null. Hypotheses: `0 < m`, `0 < N`,
`0 < a ≤ 1` (a probability that a p-value is below 1), `0 < q`, `GridFine`. -/
def kept_family_formula_stmt : Prop :=
  ∀ (m N : ℕ) (a q : ℚ), 0 < m → 0 < N → 0 < a → a ≤ 1 → 0 < q → GridFine m N a q →
    fdrKeptMix m N a q false = min (q / a) 1 * (1 - (1 - a) ^ m)

/-- (4d-ii) KIT-189's figure: at `a = 27/100`, `q = 1/20`, `m = 150`, `N = 27 lcm(1..150)`, the kept family's false
discovery rate is `(5/27)(1 - (73/100)^150)`, which exceeds 0.185. Closed statement. -/
def kept_family_instance_stmt : Prop :=
  fdrKeptMix 150 (27 * lcmUpTo 150) (27 / 100) (1 / 20) false = 5 / 27 * (1 - (73 / 100) ^ 150) ∧
    (185 / 1000 : ℚ) < 5 / 27 * (1 - (73 / 100) ^ 150)

/-! ## Witnesses (section 5.1) -/

/-- Witness for `qflag_contains_bh_stmt`: `0 < π0 < 1`, two p-values in `[0, 1]`. -/
def qflag_contains_bh_witness_stmt : Prop :=
  ∃ (p : Fin 2 → ℝ) (π0 : ℝ), 0 < π0 ∧ π0 < 1 ∧ ∀ j, 0 ≤ p j

/-- Witness for `qflag_extra_stmt`: a term the q-value flags and BH does not. -/
def qflag_extra_witness_stmt : Prop :=
  ∃ (p : Fin 2 → ℝ) (π0 α : ℝ) (i : Fin 2), 0 < π0 ∧ π0 < 1 ∧ α ≤ bh p i ∧ bh p i < α / π0

/-- Witness for the 4c statements: distinct p-values strictly inside `(0, 1)`, and a second vector above it. -/
def bh_unit_witness_stmt : Prop :=
  ∃ (p p' : Fin 2 → ℝ), UnitP p ∧ UnitP p' ∧ (∀ j, p j < p' j) ∧ p 0 ≠ p 1 ∧ 0 < p 0 ∧ 0 < p 1

/-- Witness for `grid_superuniform_stmt`, `grid_mass_stmt`: `N = 20` and `N = 25` satisfy `0 < N`. -/
def grid_witness_stmt : Prop := 0 < (20 : ℕ) ∧ 0 < (25 : ℕ)

/-- Witness for `kept_family_formula_stmt`: KIT-189's instance satisfies every hypothesis. -/
def kept_family_witness_stmt : Prop :=
  0 < (150 : ℕ) ∧ 0 < 27 * lcmUpTo 150 ∧ (0 : ℚ) < 27 / 100 ∧ (27 / 100 : ℚ) ≤ 1 ∧ (0 : ℚ) < 1 / 20 ∧
    GridFine 150 (27 * lcmUpTo 150) (27 / 100) (1 / 20)

end MultiScholaRProofs.QvalueBh
