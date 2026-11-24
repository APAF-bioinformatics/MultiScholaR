module

public import MultiScholaRProofs.Spec.QvalueBh

/-!
# Target 4 d-i's oracle table (`bh-counterexample.tsv`, KIT-193 sections 5.5 and 6.4 (iii))

For the two grids of `fdr_n20_stmt`, `fdr_n20_strict_stmt` and `fdr_n25_stmt` (`N = 20` and `N = 25`, `q = 1/20`),
one row per outcome `(j₁, j₂)` (441 and 676 rows): the number of terms BH calls over the terms with `p < 1` and over
both terms, at "adjusted p at most q" and at "below q". The Spec's own definitions evaluated (`callsKept`,
`callsWhole` on `[gridVal N j₁, gridVal N j₂]`); `j = 0` is `p = 1`, `j ≥ 1` is `p = j / (2N)`.
-/

@[expose] public section

namespace MultiScholaRProofs.Oracle.QvalueBh

open MultiScholaRProofs.QvalueBh

/-- One row. -/
def row (N : ℕ) (j₁ j₂ : Fin (N + 1)) : String :=
  let l := [gridVal N j₁, gridVal N j₂]
  let q : ℚ := 1 / 20
  String.intercalate "\t"
    [toString N, toString (j₁ : ℕ), toString (j₂ : ℕ), toString (callsKept false q l), toString (callsKept true q l),
      toString (callsWhole false q l), toString (callsWhole true q l)]

/-- The rows of one grid. -/
def grid (N : ℕ) : List String :=
  (List.finRange (N + 1)).flatMap fun j₁ => (List.finRange (N + 1)).map fun j₂ => row N j₁ j₂

/-- The whole table, with its header. -/
def table : List String :=
  ["N\tj1\tj2\tkept_le\tkept_lt\twhole_le\twhole_lt"] ++ grid 20 ++ grid 25

end MultiScholaRProofs.Oracle.QvalueBh
