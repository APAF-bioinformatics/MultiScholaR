module

public import MultiScholaRProofs.Spec.ReplayRule

/-!
# Target 2's oracle tables (KIT-193 sections 5.5 and 6.2)

The Spec's `rule`, `inForce`, `build` and `replayOptions` evaluated over the marker table after KIT-190
(`missing_value_comparison`: option `MultiScholaR.missing_filter_rule`, `counts` / `rounded_percentage`;
`ruv_control_grouping`: option `MultiScholaR.ruv_control_grouping`, `kept` / `dropped`). An entry kind is a name in
{the two markers, `other`} by a source in {`fixed`, `registry`, `args`} by a value in {`counts`, `kept`, `other`}:
27 kinds. A step is written `name|source|value` per entry, `;` between entries, `-` when empty.

* `replay-rule.tsv`: every step of at most two entries (1 + 27 + 729 = 757) and 300 three-entry steps from a fixed
  linear congruential generator; the rule per marker.
* `replay-inforce.tsv`: every subset of the fixed names {two markers, `other`} (each with its current token, `other`
  with `x`) by every option state per marker in {unset, current, earlier, `another`} (8 × 16 = 128 rows); the names
  `workflowFixedInForce` keeps.
* `replay-build.tsv`: every such provenance (8) by each declaration shape (the subsets of the two markers, with
  `other`: 4), 32 rows; the built fixed entries and the rule per marker.
* `replay-options.tsv`: the steps of `replay-rule.tsv` by the four declaration shapes; the options set.
-/

@[expose] public section

namespace MultiScholaRProofs.Oracle.ReplayRule

open MultiScholaRProofs.ReplayRule

/-- KIT-187's row. -/
def mMissing : Marker String String String :=
  ⟨"missing_value_comparison", "MultiScholaR.missing_filter_rule", "counts", "rounded_percentage"⟩

/-- KIT-190's row. -/
def mRuv : Marker String String String :=
  ⟨"ruv_control_grouping", "MultiScholaR.ruv_control_grouping", "kept", "dropped"⟩

/-- The table after KIT-190. -/
def tbl : List (Marker String String String) := [mMissing, mRuv]

def sourceText : Source → String
  | .fixed => "fixed"
  | .registry => "registry"
  | .args => "args"
  | .declared => "declared"

def entryText (e : Entry String String) : String := e.name ++ "|" ++ sourceText e.source ++ "|" ++ e.value

def stepText (s : Step String String) : String :=
  if s = [] then "-" else String.intercalate ";" (s.map entryText)

def fixedText (f : Fixed String String) : String :=
  if f = [] then "-" else String.intercalate ";" (f.map fun x => x.1 ++ "|" ++ x.2)

/-- The 27 entry kinds. -/
def kinds : List (Entry String String) :=
  [mMissing.name, mRuv.name, "other"].flatMap fun nm =>
    [Source.fixed, Source.registry, Source.args].flatMap fun src =>
      ["counts", "kept", "other"].map fun v => ⟨nm, src, v⟩

/-- The generator `x ↦ (1103515245 x + 12345) mod 2^31`. -/
def lcg (x : ℕ) : ℕ := (1103515245 * x + 12345) % 2147483648

/-- 300 three-entry steps from seed 1930. -/
def threeSteps : List (Step String String) :=
  ((List.range 300).foldl (fun (acc : List (Step String String) × ℕ) _ =>
    let x1 := lcg acc.2
    let x2 := lcg x1
    let x3 := lcg x2
    let pick := fun x => kinds.getD ((x / 65536) % 27) ⟨"", Source.args, ""⟩
    (acc.1 ++ [[pick x1, pick x2, pick x3]], x3)) ([], 1930)).1

/-- Every step of at most two entries, then the 300 three-entry steps. -/
def steps : List (Step String String) :=
  [[]] ++ kinds.map (fun e => [e]) ++ kinds.flatMap (fun e => kinds.map fun e' => [e, e']) ++ threeSteps

def ruleTable : List String :=
  ["step_id\tstep\trule_missing_value_comparison\trule_ruv_control_grouping"] ++
    (steps.zipIdx.map fun x =>
      String.intercalate "\t" [toString (x.2 + 1), stepText x.1, rule x.1 mMissing, rule x.1 mRuv])

/-- An option state for one marker: unset, current, earlier or another token. -/
def optStates (μ : Marker String String String) : List (String × Option String) :=
  [("unset", none), ("current", some μ.current), ("earlier", some μ.earlier), ("another", some "another")]

/-- The option state with the two markers' options set as given. -/
def optsOf (a b : Option String) : Opts String String :=
  fun o => if o = mMissing.option then a else if o = mRuv.option then b else none

/-- The subsets of the fixed names, each with its value (markers their current token, `other` `x`). -/
def fixedSubsets : List (Fixed String String) :=
  [[], [(mMissing.name, "counts")], [(mRuv.name, "kept")], [("other", "x")],
    [(mMissing.name, "counts"), (mRuv.name, "kept")], [(mMissing.name, "counts"), ("other", "x")],
    [(mRuv.name, "kept"), ("other", "x")], [(mMissing.name, "counts"), (mRuv.name, "kept"), ("other", "x")]]

def inforceTable : List String :=
  ["fixed\topt_missing_value_comparison\topt_ruv_control_grouping\tin_force"] ++
    fixedSubsets.flatMap fun f =>
      (optStates mMissing).flatMap fun a =>
        (optStates mRuv).map fun b =>
          String.intercalate "\t" [fixedText f, a.1, b.1, fixedText (inForce tbl f (optsOf a.2 b.2))]

/-- The declaration shapes: the subsets of the two markers, each with `other`. -/
def declShapes : List (String × Fixed String String) :=
  [("none", [("other", "x")]), ("missing", [(mMissing.name, "counts"), ("other", "x")]),
    ("ruv", [(mRuv.name, "kept"), ("other", "x")]),
    ("both", [(mMissing.name, "counts"), (mRuv.name, "kept"), ("other", "x")])]

def buildTable : List String :=
  ["provenance_fixed\tdeclaration\tbuilt\trule_missing_value_comparison\trule_ruv_control_grouping"] ++
    fixedSubsets.flatMap fun prov =>
      declShapes.map fun d =>
        let b := build tbl [] prov d.2
        String.intercalate "\t" [fixedText prov, d.1, stepText b, rule b mMissing, rule b mRuv]

def optionsText (l : List (String × String)) : String :=
  if l = [] then "-" else String.intercalate ";" (l.map fun x => x.1 ++ "=" ++ x.2)

def optionsTable : List String :=
  ["step_id\tdeclaration\toptions"] ++
    (steps.zipIdx.flatMap fun x =>
      declShapes.map fun d =>
        String.intercalate "\t" [toString (x.2 + 1), d.1, optionsText (replayOptions tbl x.1 d.2)])

end MultiScholaRProofs.Oracle.ReplayRule
