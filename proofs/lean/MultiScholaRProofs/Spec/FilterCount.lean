module

public import Mathlib.Algebra.Order.Field.Rat

/-!
# Target 1; the missing-value filter's integer rule (KIT-193 section 6.1)

Source: KIT-187 r1.3 sections 5.1 and 7.1, AUDIT-KIT-187 A09. Code (KIT-187's branch
`kit-187-filter-by-count` at e6c42d1d):

* `exceedsRoundedPercent(count, size, cutoff)` = `200000 * count > size * roundedPercentNumerator(cutoff)`,
  `R/func_prot_qc_filtering_helpers.R:304-306`;
* `roundedPercentNumerator(cutoff)` = `2 * round(1000 * cutoff) + 1`, `:319`;
* `allowedByRoundedPercent(size, cutoff)` = `floor(size * roundedPercentNumerator(cutoff) / 200000)`, `:313-315`
  (the peptide path's floor form);
* the cutoffs, `R/func_general_parameter_helpers.R:254-267`: per group
  `missing_percent = round(((n_reps - pmin(n_reps, min_reps)) / n_reps) * 100, 3)`, the `max` over groups
  (`:263`), and at the group level `round(((G - g) / G) * 100, 3)` with `g <= G` (stop at `:249-251`).

Model. Integers and rationals only. `P` is the cutoff in thousandths of a percent, an integer within 1/2 of
`100000 s` (`RoundedOK`), so R's `round(x, 3)`, `round(1000 * cutoff)` and the floating products are replaced by
the one hypothesis they must satisfy; whether the shipped arithmetic satisfies it is the paired test's job (6.1 (i)).
Not modelled: IEEE doubles (the products are exact in doubles below 2^53: stated, not proved); NA (the helper counts a
non-finite or not-above-threshold value as missing before this rule; `k` is that count); the percentile threshold;
the peptide path's double division and `floor` (exact while the product is below 2^53 and the quotient below 2^31:
stated, not proved).

Every statement here concerns `(n*, m)` for a group (`n_reps`, `min_reps_per_group`) and equally `(G, g)` for the
group level (section 6.1 (e)): the same definitions with the group counts substituted.
-/

@[expose] public section

namespace MultiScholaRProofs.FilterCount

/-- The intended allowed share of missing values for a group of `nstar` samples that needs `m` present:
`(n* - min(n*, m)) / n*` (`R/func_general_parameter_helpers.R:254-258`, `max_missing / n_reps`). The subtraction is
in `ℕ` and never truncates (`min nstar m ≤ nstar`). At `nstar = 0` the division is Lean's `x / 0 = 0`; every
statement that uses this share at a variable `nstar` assumes `0 < nstar` (no group has zero samples). -/
def allowedShare (nstar m : ℕ) : ℚ := ((nstar - min nstar m : ℕ) : ℚ) / nstar

/-- The shipped rule: a group of `n` samples with `k` missing fails when `200000 k > n (2P + 1)`
(`exceedsRoundedPercent`, `R/func_prot_qc_filtering_helpers.R:304-306`, with `roundedPercentNumerator` = `2P + 1`,
`:319`). Exact integer arithmetic. -/
def failsRule (P : ℤ) (n k : ℕ) : Prop := (200000 : ℤ) * k > (n : ℤ) * (2 * P + 1)

/-- The intended rule: fail when the missing share `k / n` exceeds `s` (KIT-187 5.1). Used with `0 < n` only. -/
def failsIntended (s : ℚ) (n k : ℕ) : Prop := (k : ℚ) / n > s

/-- `P` is a rounding of `100000 s` to an integer, ties either way: `|P - 100000 s| ≤ 1/2`. This replaces R's
`round(x, 3)` of the percentage and `round(1000 * cutoff)` (section 6.1, "Which version is proved"). -/
def RoundedOK (P : ℤ) (s : ℚ) : Prop := |(P : ℚ) - 100000 * s| ≤ 1 / 2

/-- (6.1 (a)) The rule is never stricter than the intended rule: whatever the share `s` and its admissible
rounding `P`, a group the rule fails has a missing share above `s`.
Hypotheses: `0 < n` (a group has at least one sample; `k / n` is the real quotient); `RoundedOK P s` (the cutoff is
a rounding of the share). -/
def never_stricter_stmt : Prop :=
  ∀ (s : ℚ) (P : ℤ) (n k : ℕ), 0 < n → RoundedOK P s → failsRule P n k → failsIntended s n k

/-- (6.1 (b)) Exactness within the lcm bound: when `lcm(n, n*) < 100000`, the rule and the intended rule agree for
every `k`, `n*` being the size of the group whose share the cutoff rounds (the largest group, by (d)).
Hypotheses: `0 < nstar`, `0 < n` (groups are non-empty; the share and `k / n` are real quotients);
`RoundedOK P (allowedShare nstar m)` (the shipped cutoff rounds that share; checked by the paired test (i));
`Nat.lcm n nstar < 100000` (the bound; sufficient, not sharp). -/
def exact_lcm_stmt : Prop :=
  ∀ (nstar m n k : ℕ) (P : ℤ), 0 < nstar → 0 < n →
    RoundedOK P (allowedShare nstar m) → Nat.lcm n nstar < 100000 →
    (failsRule P n k ↔ failsIntended (allowedShare nstar m) n k)

/-- (6.1 (b)) The same with the bound every design meets whose groups all have at most 316 samples
(`lcm(a, b) ≤ a b ≤ 316^2 = 99856 < 100000` for `a, b ≤ 316`).
Hypotheses as `exact_lcm_stmt`, with `n ≤ 316 ∧ nstar ≤ 316` in place of the lcm bound. -/
def exact_316_stmt : Prop :=
  ∀ (nstar m n k : ℕ) (P : ℤ), 0 < nstar → 0 < n →
    RoundedOK P (allowedShare nstar m) → n ≤ 316 → nstar ≤ 316 →
    (failsRule P n k ↔ failsIntended (allowedShare nstar m) n k)

/-- (6.1 (c)) Balanced designs below 100,000 per group: the rule is "at least `min(n, m)` present", that is it fails
exactly when more than `n - min(n, m)` values are missing. The subtraction in `ℕ` never truncates.
Hypotheses: `0 < n` (non-empty group); `n < 100000` (the lcm bound for `lcm(n, n) = n`); `RoundedOK` as above. -/
def balanced_stmt : Prop :=
  ∀ (n m k : ℕ) (P : ℤ), 0 < n → n < 100000 →
    RoundedOK P (allowedShare n m) → (failsRule P n k ↔ n - min n m < k)

/-- (6.1 (d)) The allowed share grows with group size, so the `max` over groups at
`R/func_general_parameter_helpers.R:263` is the largest group's share (the code takes the `max` of the rounded
percentages; that R's `round` is monotone is not modelled).
Hypotheses: `0 < a` (non-empty group; the smaller size's share is a real quotient); `a ≤ b`. -/
def share_mono_stmt : Prop :=
  ∀ (m a b : ℕ), 0 < a → a ≤ b → allowedShare a m ≤ allowedShare b m

/-- (6.1 (f)) The peptide path's floor form is the same rule: `k` passes exactly when
`k ≤ ⌊n (2P + 1) / 200000⌋` (`allowedByRoundedPercent`, `:313-315`). Lean's `/` on `ℤ` rounds down for a positive
divisor whatever the sign of the dividend, so no hypothesis on `P` is needed. -/
def peptide_floor_stmt : Prop :=
  ∀ (P : ℤ) (n k : ℕ), ((k : ℤ) ≤ ((n : ℤ) * (2 * P + 1)) / 200000 ↔ ¬ failsRule P n k)

/-- (6.1 (g)) Beyond the bound the rule can pass a group the intended rule fails, with one admissible rounding:
groups of 313 and 320, `m = 137`, `k = 179`, `P = 57188` (`100000 * 183/320 = 57187.5`, a tie; with `P = 57187`
the rule is exact at this pair). Closed statement, no hypothesis. -/
def counterexample_stmt : Prop :=
  RoundedOK 57188 (allowedShare 320 137) ∧ failsIntended (allowedShare 320 137) 313 179 ∧
    ¬ failsRule 57188 313 179

/-- (6.1 (g)) The same whatever the rounding: groups of 321 and 320, `m = 320`, `k = 1`; the only admissible `P` is
312 (`100000 / 321 = 311.53…`). Hypothesis: `RoundedOK P (allowedShare 321 320)`, any admissible rounding. -/
def counterexample_any_rounding_stmt : Prop :=
  ∀ P : ℤ, RoundedOK P (allowedShare 321 320) →
    failsIntended (allowedShare 321 320) 320 1 ∧ ¬ failsRule P 320 1

/-! ## Witnesses (section 5.1, no vacuous statement)

Each statement with hypotheses has a witness statement: its hypotheses hold together on an explicit instance with a
non-degenerate group (`0 < m < n*`, so the allowed share is strictly between 0 and 1). The witness theorems are
proved in `Proof/FilterCount.lean`. -/

/-- Witness for `never_stricter_stmt`. -/
def never_stricter_witness_stmt : Prop :=
  ∃ (s : ℚ) (P : ℤ) (n k : ℕ), 0 < s ∧ s < 1 ∧ 0 < n ∧ RoundedOK P s ∧ failsRule P n k

/-- Witness for `exact_lcm_stmt`. -/
def exact_lcm_witness_stmt : Prop :=
  ∃ (nstar m n : ℕ) (P : ℤ), 0 < m ∧ m < nstar ∧ 0 < n ∧
    RoundedOK P (allowedShare nstar m) ∧ Nat.lcm n nstar < 100000

/-- Witness for `exact_316_stmt`. -/
def exact_316_witness_stmt : Prop :=
  ∃ (nstar m n : ℕ) (P : ℤ), 0 < m ∧ m < nstar ∧ 0 < n ∧
    RoundedOK P (allowedShare nstar m) ∧ n ≤ 316 ∧ nstar ≤ 316

/-- Witness for `balanced_stmt`. -/
def balanced_witness_stmt : Prop :=
  ∃ (n m : ℕ) (P : ℤ), 0 < m ∧ m < n ∧ n < 100000 ∧ RoundedOK P (allowedShare n m)

/-- Witness for `share_mono_stmt`. -/
def share_mono_witness_stmt : Prop :=
  ∃ (m a b : ℕ), 0 < m ∧ m < a ∧ a < b

/-- Witness for `counterexample_any_rounding_stmt`: an admissible rounding exists. -/
def counterexample_any_rounding_witness_stmt : Prop :=
  ∃ P : ℤ, RoundedOK P (allowedShare 321 320)

end MultiScholaRProofs.FilterCount
