import MultiScholaRProofs

/-!
# KIT-193 oracle tables (section 5.5)

Run with Lean's interpreter, from `proofs/lean`:

    lake env lean --run Kit193Oracles.lean <output directory>

It writes each table the Oracle modules define, one file per table. The interpreter runs the Spec definitions'
compiled intermediate code from the build (Mathlib's included), so no native code of Mathlib is compiled; Lean's
compiler and runtime are trusted for the tables (section 9.1), not for the theorems.
-/

open MultiScholaRProofs

/-- Write one table. -/
def writeTable (dir name : String) (lines : List String) : IO Unit :=
  IO.FS.writeFile (dir ++ "/" ++ name) (String.intercalate "\n" lines ++ "\n")

def main (args : List String) : IO UInt32 := do
  match args with
  | [dir] =>
    IO.FS.createDirAll dir
    writeTable dir "go-closure.tsv" Oracle.GoClosure.table
    writeTable dir "filter-thresholds.tsv" Oracle.FilterCount.table
    writeTable dir "replay-rule.tsv" Oracle.ReplayRule.ruleTable
    writeTable dir "replay-inforce.tsv" Oracle.ReplayRule.inforceTable
    writeTable dir "replay-build.tsv" Oracle.ReplayRule.buildTable
    writeTable dir "replay-options.tsv" Oracle.ReplayRule.optionsTable
    return 0
  | _ =>
    IO.eprintln "usage: lake env lean --run Kit193Oracles.lean <output directory>"
    return 2
