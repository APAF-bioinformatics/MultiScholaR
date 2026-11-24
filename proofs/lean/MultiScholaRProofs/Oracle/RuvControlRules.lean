module

public import MultiScholaRProofs.Spec.RuvControlRules

/-!
# R-N1's oracle table for the first-pass rules (`first-pass.tsv`, KIT-193 sections 5.5 and 7.4)

One row per feature of each case: the case, the feature's 1-based row, its p-value in eighths (`p8`, the p-value is
`p8 / 8`; `NA` for a missing one), whether it is a first-pass control of proteomics, lipidomics and metabolomics
(`firstPassControls`) and whether it is one of transcriptomics (`trxFirstPassControls`): the Spec's own definitions
evaluated. Cases: hand-made ones (every p equal, tied blocks, a missing first row, the five genes of
`N1_rules_differ_stmt`, every p missing, one feature, ties at 0); every vector of four over
`{NA, 1/8, 4/8, 1}` (256 cases); 100 vectors of nine over `{NA, 0, 1/8, …, 1}` from a fixed linear congruential
generator. Dyadic values, so R's doubles hold them exactly and ties are exact. The sweep's rule (`sweepControls`) is
not in this table: its pair waits for KIT-190.
-/

@[expose] public section

namespace MultiScholaRProofs.Oracle.RuvControlRules

open MultiScholaRProofs.RuvControls

/-- A p-value in eighths. -/
def eighth (k : ℕ) : Option ℚ := some ((k : ℚ) / 8)

/-- The rows of one case. -/
def rows (name : String) (l : List (Option ℕ)) : List String :=
  let pv : Fin l.length → Option ℚ := fun i => (l.get i).bind eighth
  let fp := firstPassControls pv
  let trx := trxFirstPassControls pv
  (List.finRange l.length).map fun (i : Fin l.length) =>
    String.intercalate "\t"
      [name, toString ((i : ℕ) + 1), (match l.get i with | some k => toString k | none => "NA"),
        if i ∈ fp then "1" else "0", if i ∈ trx then "1" else "0"]

/-- Hand-made cases. -/
def handCases : List (String × List (Option ℕ)) :=
  [("all_equal", List.replicate 8 (some 4)),
    ("tied_blocks", [some 1, some 1, some 1, some 4, some 4, some 4, some 8, some 8]),
    ("missing_first", [none, some 2, some 1, some 3]),
    ("five_genes", [some 1, some 2, some 3, some 4, some 5]),
    ("missing_one_of_four", [none, some 1, some 2, some 3]),
    ("all_missing", [none, none, none]),
    ("one", [some 3]),
    ("zeros_tied", [some 0, some 0, some 0, some 0])]

/-- Every vector of four over `{NA, 1/8, 4/8, 1}`. -/
def fourCases : List (String × List (Option ℕ)) :=
  let vals : List (Option ℕ) := [none, some 1, some 4, some 8]
  let all := vals.flatMap fun a => vals.flatMap fun b => vals.flatMap fun c => vals.map fun d => [a, b, c, d]
  (List.range all.length).zip all |>.map fun x => ("four_" ++ toString (x.1 + 1), x.2)

/-- The generator `x ↦ (1103515245 x + 12345) mod 2^31`, as the target 3 oracle's. -/
def lcg (x : ℕ) : ℕ := (1103515245 * x + 12345) % 2147483648

/-- Nine values from state `x`: `NA` or `0 … 8` eighths, and the next state. -/
def nine : ℕ → ℕ → List (Option ℕ) × ℕ
  | 0, x => ([], x)
  | c + 1, x =>
    let y := lcg x
    let v := (y / 65536) % 10
    let (rest, z) := nine c y
    ((if v = 9 then none else some v) :: rest, z)

/-- 100 generated vectors of nine. -/
def genCases : List (String × List (Option ℕ)) :=
  ((List.range 100).foldl (fun (acc : List (String × List (Option ℕ)) × ℕ) i =>
      let (l, x') := nine 9 acc.2
      (acc.1 ++ [("gen_" ++ toString (i + 1), l)], x')) ([], 20261124)).1

/-- The whole table, with its header. -/
def table : List String :=
  ["case\trow\tp8\tfirst_pass\ttrx"] ++ (handCases ++ fourCases ++ genCases).flatMap fun c => rows c.1 c.2

end MultiScholaRProofs.Oracle.RuvControlRules
