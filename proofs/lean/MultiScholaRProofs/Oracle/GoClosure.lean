module

public import MultiScholaRProofs.Spec.GoClosure
public import Mathlib.Data.String.Basic
public import Mathlib.Data.Finset.Sort

/-!
# Target 3's oracle table (`go-closure.tsv`, KIT-193 sections 5.5 and 6.3)

The table is the Spec's own definitions evaluated: `bfs` (the transcription of `goAncestorPairs`),
`propagateFin` (of `propagateGoAnnotations`) and the root removal of `buildUniprotGoGeneSets`, on hand-made graphs
and on 200 graphs of 6 terms from a fixed linear congruential generator. Every case is in the aspect `BP`, whose
root is `GO:0008150` (`goAspectRoots`). One long table, `case`, `record`, `x`, `y`, `z`:

* inputs: `labels` (the relation set's name, `is_a_part_of` or `is_a_part_of_regulates`), `edge` (child, parent,
  relation), `aspect` (term, aspect: the terms `goDbTermAspects` holds), `annot` (term, gene);
* outputs: `pair` (term, node: `goAncestorPairs` from the distinct kept terms), `prop` (term, gene:
  `propagateGoAnnotations`), `build` (term, gene: the propagated pairs without the root).
-/

@[expose] public section

namespace MultiScholaRProofs.Oracle.GoClosure

open MultiScholaRProofs.GoClosure

/-- One oracle case: its edges, aspect map (as a list), annotations and relation set. -/
structure Case where
  name : String
  labels : String
  edges : List (String × String × String)
  aspects : List (String × String)
  annots : List (String × String)

/-- The two relation sets of `goPropagationRelations`. -/
def relationSet (labels : String) : Finset String :=
  if labels = "is_a_part_of" then {"isa", "part of"}
  else {"isa", "part of", "regulates", "positively regulates", "negatively regulates"}

/-- The aspect map as a function (`none` for a term GO.db does not hold). -/
def aspectOf (c : Case) (t : String) : Option String := (c.aspects.lookup t)

/-- The aspect root of `BP`. -/
def rootBP : String := "GO:0008150"

/-- The case's graph. -/
def graph (c : Case) : GoGraph String String := ⟨c.edges, "all"⟩

/-- The kept annotations' distinct terms, the search's start terms. -/
def startTerms (c : Case) : Finset String :=
  ((c.annots.toFinset.filter fun x => decide (aspectOf c x.1 = some "BP")).image Prod.fst)

/-- Sorted text lines of a set of pairs. -/
def pairLines (s : Finset (String × String)) : List String :=
  (s.image fun x => x.1 ++ "\t" ++ x.2).sort (· ≤ ·)

/-- The case's rows of the table. -/
def rows (c : Case) : List String :=
  let pre := c.name ++ "\t"
  let Ls := relationSet c.labels
  let pairs := bfs (graph c) Ls (startTerms c)
  let prop := propagateFin (graph c) Ls (aspectOf c) "BP" c.annots.toFinset
  let built := prop.filter fun x => decide (x.1 ≠ rootBP)
  [pre ++ "labels\t" ++ c.labels ++ "\t\t"] ++
  c.edges.map (fun e => pre ++ "edge\t" ++ e.1 ++ "\t" ++ e.2.1 ++ "\t" ++ e.2.2) ++
  c.aspects.map (fun a => pre ++ "aspect\t" ++ a.1 ++ "\t" ++ a.2 ++ "\t") ++
  c.annots.map (fun a => pre ++ "annot\t" ++ a.1 ++ "\t" ++ a.2 ++ "\t") ++
  (pairLines pairs).map (fun l => pre ++ "pair\t" ++ l ++ "\t") ++
  (pairLines prop).map (fun l => pre ++ "prop\t" ++ l ++ "\t") ++
  (pairLines built).map (fun l => pre ++ "build\t" ++ l ++ "\t")

/-- Term names for the hand-made cases. -/
def t (i : ℕ) : String := "GO:91000" ++ (if i < 10 then "0" else "") ++ toString i

/-- Hand-made cases (6.3's list): chain, diamond, cycle, self-loop, an edge to `all`, a cross-aspect edge, a
`regulates` edge under both relation sets, the aspect root, an input term the aspect map does not hold, an input
term of another aspect. -/
def handCases : List Case :=
  let bp := fun (ts : List String) => ts.map fun x => (x, "BP")
  [ ⟨"chain", "is_a_part_of", [(t 1, t 2, "isa"), (t 2, t 3, "part of"), (t 3, t 4, "isa")],
      bp [t 1, t 2, t 3, t 4], [(t 1, "g1"), (t 3, "g2")]⟩,
    ⟨"diamond", "is_a_part_of", [(t 1, t 2, "isa"), (t 1, t 3, "isa"), (t 2, t 4, "part of"), (t 3, t 4, "isa")],
      bp [t 1, t 2, t 3, t 4], [(t 1, "g1"), (t 2, "g2")]⟩,
    ⟨"cycle", "is_a_part_of", [(t 1, t 2, "isa"), (t 2, t 3, "isa"), (t 3, t 1, "isa")],
      bp [t 1, t 2, t 3], [(t 1, "g1")]⟩,
    ⟨"self_loop", "is_a_part_of", [(t 1, t 1, "isa"), (t 1, t 2, "isa")], bp [t 1, t 2], [(t 1, "g1")]⟩,
    ⟨"edge_to_all", "is_a_part_of", [(t 1, t 2, "isa"), (t 2, "all", "isa")], bp [t 1, t 2], [(t 1, "g1")]⟩,
    ⟨"cross_aspect", "is_a_part_of", [(t 1, t 2, "isa"), (t 2, t 3, "part of")],
      [(t 1, "BP"), (t 2, "MF"), (t 3, "BP")], [(t 1, "g1"), (t 2, "g2")]⟩,
    ⟨"regulates_two", "is_a_part_of", [(t 1, t 2, "regulates"), (t 1, t 3, "isa")], bp [t 1, t 2, t 3],
      [(t 1, "g1")]⟩,
    ⟨"regulates_five", "is_a_part_of_regulates", [(t 1, t 2, "regulates"), (t 1, t 3, "isa")], bp [t 1, t 2, t 3],
      [(t 1, "g1")]⟩,
    ⟨"aspect_root", "is_a_part_of", [(t 1, t 2, "isa"), (t 2, rootBP, "isa"), (rootBP, "all", "isa")],
      bp [t 1, t 2, rootBP], [(t 1, "g1"), (rootBP, "g2")]⟩,
    ⟨"not_in_go_db", "is_a_part_of", [(t 1, t 2, "isa")], bp [t 2], [(t 1, "g1"), (t 2, "g2")]⟩,
    ⟨"other_aspect_input", "is_a_part_of", [(t 1, t 2, "isa"), (t 3, t 2, "isa")],
      [(t 1, "BP"), (t 2, "BP"), (t 3, "CC")], [(t 1, "g1"), (t 3, "g3")]⟩ ]

/-- The generator `x ↦ (1103515245 x + 12345) mod 2^31` (fixed, so the table is the same on every machine). -/
def lcg (x : ℕ) : ℕ := (1103515245 * x + 12345) % 2147483648

/-- Draw a number below `k` (`k > 0`) and the next state, from the generator's high bits. -/
def draw (x k : ℕ) : ℕ × ℕ := let y := lcg x; ((y / 65536) % k, y)

/-- `count` draws, each turned into a value by `f`. -/
def draws {α : Type} (f : ℕ → ℕ → α × ℕ) : ℕ → ℕ → List α × ℕ
  | 0, x => ([], x)
  | count + 1, x =>
    let (a, x') := f x 0
    let (rest, x'') := draws f count x'
    (a :: rest, x'')

/-- The six terms of a generated graph: the `BP` root and five others. -/
def genTerms : List String := [rootBP, t 51, t 52, t 53, t 54, t 55]

/-- The five relation labels GO.db's parent maps return. -/
def genLabels : List String := ["isa", "part of", "regulates", "positively regulates", "negatively regulates"]

/-- One generated case from state `x`: 4 to 10 edges (child among the six terms, parent among them and `all`,
self-loops allowed), each term `BP` (4 in 6), `MF` (1 in 6) or not held (1 in 6) with the root always `BP`, and 3 to
6 annotations over genes `g1` to `g4`; the relation set alternates with the case number. -/
def genCase (i x : ℕ) : Case × ℕ :=
  let (ne, x) := draw x 7
  let (edges, x) := draws (fun x _ =>
      let (c, x) := draw x 6
      let (p, x) := draw x 7
      let (l, x) := draw x 5
      ((genTerms.getD c "", (genTerms ++ ["all"]).getD p "", genLabels.getD l ""), x)) (ne + 4) x
  let (asp, x) := draws (fun x _ => draw x 6) 6 x
  let aspects := (genTerms.zip asp).filterMap fun ta =>
    if ta.1 = rootBP then some (ta.1, "BP")
    else if ta.2 < 4 then some (ta.1, "BP") else if ta.2 = 4 then some (ta.1, "MF") else none
  let (na, x) := draw x 4
  let (annots, x) := draws (fun x _ =>
      let (tt, x) := draw x 6
      let (g, x) := draw x 4
      ((genTerms.getD tt "", "g" ++ toString (g + 1)), x)) (na + 3) x
  (⟨"gen" ++ toString i, if i % 2 = 0 then "is_a_part_of" else "is_a_part_of_regulates", edges, aspects, annots⟩, x)

/-- The 200 generated cases, from seed 193. -/
def genCases : List Case :=
  ((List.range 200).foldl (fun (acc : List Case × ℕ) i =>
    let (c, x) := genCase (i + 1) acc.2
    (acc.1 ++ [c], x)) ([], 193)).1

/-- The whole table, with its header. -/
def table : List String :=
  ["case\trecord\tx\ty\tz"] ++ (handCases ++ genCases).flatMap rows

end MultiScholaRProofs.Oracle.GoClosure
