module

public import MultiScholaRProofs.Spec.FilterCount
public import Mathlib.Algebra.Order.Floor.Ring
public import Mathlib.Data.Rat.Floor

/-!
# Target 1's oracle table (`filter-thresholds.tsv`, KIT-193 sections 5.5 and 6.1 (ii))

One row per `(n*, m, n)` with `1 ≤ m ≤ n* ≤ 60` and `1 ≤ n ≤ n*` (73,810 rows): `P`, the integer nearest
`100000 · allowedShare n* m` (`⌊100000 s + 1/2⌋`; below 64 there is no tie, so it is the only integer with
`RoundedOK`), `RoundedOK P s` decided from the Spec's definition (`1` when it holds), and `kmax`, the largest `k` that
does not fail, by the proved floor form `⌊n (2P + 1) / 200000⌋` (`peptide_floor_stmt`).
-/

@[expose] public section

namespace MultiScholaRProofs.Oracle.FilterCount

open MultiScholaRProofs.FilterCount

instance (P : ℤ) (s : ℚ) : Decidable (RoundedOK P s) := by
  unfold RoundedOK
  infer_instance

/-- The nearest integer to `100000 s`. -/
def roundP (s : ℚ) : ℤ := ⌊100000 * s + 1 / 2⌋

/-- The largest `k` that passes, by the floor form. -/
def kmax (P : ℤ) (n : ℕ) : ℤ := ((n : ℤ) * (2 * P + 1)) / 200000

/-- One row. -/
def row (nstar m n : ℕ) : String :=
  let s := allowedShare nstar m
  let P := roundP s
  String.intercalate "\t"
    [toString nstar, toString m, toString n, toString P, if decide (RoundedOK P s) then "1" else "0",
      toString (kmax P n)]

/-- The whole table, with its header. -/
def table : List String :=
  ["nstar\tm\tn\tP\trounded_ok\tkmax"] ++
    (List.range' 1 60).flatMap fun nstar =>
      (List.range' 1 nstar).flatMap fun m => (List.range' 1 nstar).map fun n => row nstar m n

end MultiScholaRProofs.Oracle.FilterCount
