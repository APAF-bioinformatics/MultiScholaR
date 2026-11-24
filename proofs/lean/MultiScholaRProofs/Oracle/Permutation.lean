module

public import MultiScholaRProofs.Spec.Permutation

/-!
# Target 5's oracle table (`perm-floor.tsv`, KIT-193 sections 5.5 and 6.5)

For each fixture of 6.5's paired test (A: `g = (1, 2, 4, 8)/8`; B: `(1, 1, 4, 8)/8`, ties; C: `(1, 2, 3, 4)/4`,
symmetric about its mean) and each of the `24³ = 13,824` cases `(τ, π₁, π₂)` of four samples and two draws: the
factor `f = g ∘ τ`, `d` = the number of draws with `f ∘ π_b = f` (the Spec's `unchangedDraws f π`), `r` = the number
of draws with `f ∘ π_b = reflect f` and `f ∘ π_b ≠ f` (the Spec's `reflectedDraws f π`), and the floor `(1 + d) / 3`
(`perm_floor_stmt` with `B = 2`). Compositions of orderings and equality of rational vectors only; no statistic
(the model's `T` is abstract).

The Spec's `unchangedDraws` and `reflectedDraws` count with classical decidability (they are stated for any type of
values), so they cannot be evaluated; the counts below are the same filters written over lists of rationals with
their decidable equality. That the two agree (a filter does not depend on the decidability instance,
`Finset.filter_congr_decidable`, and a function on `Fin 4` is its list of values) is the bridge lemma ticket 5.5
asks of a separate executable function; it is not in the Proof files yet (recorded in `outcome.md`).

Orderings are written one-based as R indexes them: `τ = "2143"` is `f = g[c(2, 1, 4, 3)]`, `π = "4321"` is the draw
`f[c(4, 3, 2, 1)]` (`perms[, b]` of the shipped code). They are listed in the order of `List.permutations`.
-/

@[expose] public section

namespace MultiScholaRProofs.Oracle.Permutation

/-- The 24 orderings of four samples, 0-based, as lists of values `σ(0), …, σ(3)`. -/
def orderings : List (List ℕ) := List.permutations [0, 1, 2, 3]

/-- `f ∘ σ` for a factor and an ordering, both as lists (`f[σ]`). -/
def compose (f : List ℚ) (σ : List ℕ) : List ℚ := σ.map fun i => f.getD i 0

/-- The reflection about the mean, `2 mean(f) - f` (the Spec's `reflect`). -/
def reflectL (f : List ℚ) : List ℚ := f.map fun x => 2 * (f.sum / f.length) - x

/-- `d`: draws that leave `f` unchanged. -/
def unchanged (f : List ℚ) (draws : List (List ℕ)) : ℕ := (draws.filter fun π => decide (compose f π = f)).length

/-- `r`: draws that map `f` to its reflection and change it. -/
def reflected (f : List ℚ) (draws : List (List ℕ)) : ℕ :=
  (draws.filter fun π => decide (compose f π = reflectL f ∧ compose f π ≠ f)).length

/-- One-based digits of an ordering. -/
def label (σ : List ℕ) : String := String.join (σ.map fun i => toString (i + 1))

/-- The three fixtures. -/
def fixtures : List (String × List ℚ) :=
  [("A", [1 / 8, 2 / 8, 4 / 8, 8 / 8]), ("B", [1 / 8, 1 / 8, 4 / 8, 8 / 8]), ("C", [1 / 4, 2 / 4, 3 / 4, 4 / 4])]

/-- One row: fixture, τ, π₁, π₂, d, r, (1 + d)/3. -/
def row (name : String) (g : List ℚ) (τ π₁ π₂ : List ℕ) : String :=
  let f := compose g τ
  let d := unchanged f [π₁, π₂]
  String.intercalate "\t"
    [name, label τ, label π₁, label π₂, toString d, toString (reflected f [π₁, π₂]), toString ((1 + d : ℚ) / 3)]

/-- The whole table, with its header. -/
def table : List String :=
  ["fixture\ttau\tpi1\tpi2\td\tr\tmin_p"] ++
    fixtures.flatMap fun x => orderings.flatMap fun τ => orderings.flatMap fun π₁ =>
      orderings.map fun π₂ => row x.1 x.2 τ π₁ π₂

end MultiScholaRProofs.Oracle.Permutation
