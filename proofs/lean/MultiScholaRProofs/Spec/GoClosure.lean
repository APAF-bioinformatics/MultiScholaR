module

public import Mathlib.Logic.Relation
public import Mathlib.Data.Finset.Prod
public import Mathlib.Data.Finset.Card
public import Mathlib.Data.Set.Basic

/-!
# Target 3: Gene Ontology propagation (KIT-193 section 6.3)

Source: KIT-173 sections 2.2 and 7.2; the code merged at d12f04ae, read on `v0.9-staging` at 51959495:
`R/func_general_geneset_kernel.R:72-92` (`goAncestorPairs`: edges filtered by relation label and
`parent != "all"`, `:79`; a breadth-first search from each term with a seen set, starting from (term, term), `:80`),
`:109-148` (`propagateGoAnnotations`: input terms GO.db does not hold, or of another aspect, dropped, `:128-131`;
(term, gene) joined to (term, node); distinct; sorted), `:191-236` (`buildUniprotGoGeneSets`: per aspect, the
aspect root removed, `:204-208`, then pooled, distinct, sorted). `goDbTermAspects` (`:52-56`) keeps BP, MF, CC, so
`all` has no aspect.

Model. Terms `T`, genes `G`, relation labels `L` and aspects `Asp` are abstract types with decidable equality
(identifiers). Data-frame operations are finite sets: `unique`, `merge` and sorting are not modelled; order is not
part of any claim. GO.db's parent map is the given edge list; `goDbTermAspects` is the given function `aspectOf`.
That GO.db's `GO*ANCESTOR` maps equal this closure on the real graph is KIT-173 T2's executed check, not claimed.
-/

@[expose] public section

namespace MultiScholaRProofs.GoClosure

variable {T G L Asp : Type} [DecidableEq T] [DecidableEq G] [DecidableEq L] [DecidableEq Asp]

/-- GO.db's parent table of one aspect (`goParentTable`): (child, parent, relation label) rows, and the top term
`all`. -/
structure GoGraph (T L : Type) where
  edges : List (T × T × L)
  all : T

/-- One step up the graph over the allowed labels `Ls`: an edge `(c, p, l)` with `l ∈ Ls` and `p ≠ all` (`:79`). -/
def parentRel (Γ : GoGraph T L) (Ls : Finset L) (c p : T) : Prop :=
  ∃ l, (c, p, l) ∈ Γ.edges ∧ l ∈ Ls ∧ p ≠ Γ.all

/-- The ancestor relation: the reflexive-transitive closure of `parentRel`. No DAG assumption: cycles and
self-loops are allowed. -/
def anc (Γ : GoGraph T L) (Ls : Finset L) : T → T → Prop :=
  Relation.ReflTransGen (parentRel Γ Ls)

/-- The edges `goAncestorPairs` follows (`:79`), as (child, parent) pairs. -/
def keptEdges (Γ : GoGraph T L) (Ls : Finset L) : List (T × T) :=
  (Γ.edges.filter fun e => decide (e.2.2 ∈ Ls) && decide (e.2.1 ≠ Γ.all)).map fun e => (e.1, e.2.1)

/-- One pass of the `while` loop (`:82-90`): the frontier's (term, node) pairs joined to the edges whose child is
the node give (term, parent) pairs; those not yet seen are the new frontier and are added to the seen set. -/
def bfsStep (E : List (T × T)) (st : Finset (T × T) × Finset (T × T)) : Finset (T × T) × Finset (T × T) :=
  let found : Finset (T × T) :=
    st.2.biUnion fun x => ((E.filter fun e => decide (e.1 = x.2)).map fun e => (x.1, e.2)).toFinset
  let fresh := found \ st.1
  (st.1 ∪ fresh, fresh)

/-- The loop, run for at most `fuel` passes; it stops early when the frontier is empty, as R's `while` does. -/
def bfsLoop (E : List (T × T)) : ℕ → Finset (T × T) × Finset (T × T) → Finset (T × T) × Finset (T × T)
  | 0, st => st
  | f + 1, st => if st.2 = ∅ then st else bfsLoop E f (bfsStep E st)

/-- Every term the search can reach: the start terms and every parent in the edges. -/
def nodes (E : List (T × T)) (terms : Finset T) : Finset T :=
  terms ∪ (E.map Prod.snd).toFinset

/-- The pass budget: one more than the number of (term, node) pairs, so that the loop, which adds at least one new
pair per pass that does not stop, always reaches an empty frontier (`bfs_terminates_stmt`). -/
def fuel (E : List (T × T)) (terms : Finset T) : ℕ :=
  terms.card * (nodes E terms).card + 1

/-- The loop's start (`:80-81`): each term with itself, as seen set and frontier. -/
def bfsStart (terms : Finset T) : Finset (T × T) × Finset (T × T) :=
  (terms.image fun t => (t, t), terms.image fun t => (t, t))

/-- The search's final state. -/
def bfsState (Γ : GoGraph T L) (Ls : Finset L) (terms : Finset T) : Finset (T × T) × Finset (T × T) :=
  bfsLoop (keptEdges Γ Ls) (fuel (keptEdges Γ Ls) terms) (bfsStart terms)

/-- `goAncestorPairs(terms, parents, relations)`: the seen set at the end, as (term, node) pairs. -/
def bfs (Γ : GoGraph T L) (Ls : Finset L) (terms : Finset T) : Finset (T × T) :=
  (bfsState Γ Ls terms).1

/-- The input filter of `propagateGoAnnotations` (`:128-131`): the (term, gene) pairs whose term GO.db holds in the
aspect `o`; a term GO.db does not hold (`aspectOf t = none`) or of another aspect is dropped. -/
def kept (aspectOf : T → Option Asp) (o : Asp) (A : Set (T × G)) : Set (T × G) :=
  {x | x ∈ A ∧ aspectOf x.1 = some o}

/-- Propagation as a closure image: `(a, g)` whenever some kept `(t, g)` has `a` as an ancestor of `t`. -/
def propagate (Γ : GoGraph T L) (Ls : Finset L) (aspectOf : T → Option Asp) (o : Asp) (A : Set (T × G)) :
    Set (T × G) :=
  {x | ∃ t, (t, x.2) ∈ kept aspectOf o A ∧ anc Γ Ls t x.1}

/-- `propagateGoAnnotations` as the code computes it (`:128-138`): the kept pairs, the search from their distinct
terms, the join on the term, the (node, gene) pairs. -/
def propagateFin (Γ : GoGraph T L) (Ls : Finset L) (aspectOf : T → Option Asp) (o : Asp) (A : Finset (T × G)) :
    Finset (T × G) :=
  let keptA := A.filter fun x => decide (aspectOf x.1 = some o)
  let pairs := bfs Γ Ls (keptA.image Prod.fst)
  ((keptA ×ˢ pairs).filter fun x => decide (x.1.1 = x.2.1)).image fun x => (x.2.2, x.1.2)

/-- `buildUniprotGoGeneSets` for one aspect (`:204-208`): the propagated pairs with the aspect root removed. -/
def build (Γ : GoGraph T L) (Ls : Finset L) (aspectOf : T → Option Asp) (o : Asp) (root : T)
    (A : Set (T × G)) : Set (T × G) :=
  {x | x ∈ propagate Γ Ls aspectOf o A ∧ x.1 ≠ root}

/-- (3a) The loop ends with an empty frontier within its pass budget (termination of R's `while`). No hypothesis. -/
def bfs_terminates_stmt : Prop :=
  ∀ (T L : Type) [DecidableEq T] [DecidableEq L] (Γ : GoGraph T L) (Ls : Finset L) (terms : Finset T),
    (bfsState Γ Ls terms).2 = ∅

/-- (3a) The search returns exactly the (term, ancestor) pairs of the start terms: `(t, a) ∈ bfs terms` iff `t` is
a start term and `a` an ancestor of `t` (reflexively). No hypothesis; cycles and self-loops allowed. -/
def bfs_correct_stmt : Prop :=
  ∀ (T L : Type) [DecidableEq T] [DecidableEq L] (Γ : GoGraph T L) (Ls : Finset L) (terms : Finset T) (t a : T),
    (t, a) ∈ bfs Γ Ls terms ↔ t ∈ terms ∧ anc Γ Ls t a

/-- (3b) The code's propagation equals the closure image. No hypothesis. -/
def propagate_fin_stmt : Prop :=
  ∀ (T G L Asp : Type) [DecidableEq T] [DecidableEq G] [DecidableEq L] [DecidableEq Asp] (Γ : GoGraph T L)
    (Ls : Finset L) (aspectOf : T → Option Asp) (o : Asp) (A : Finset (T × G)),
    ((propagateFin Γ Ls aspectOf o A : Finset (T × G)) : Set (T × G)) = propagate Γ Ls aspectOf o (A : Set (T × G))

/-- (3b) Idempotent: propagating propagated annotations adds nothing. No hypothesis (the input filter included). -/
def idempotent_stmt : Prop :=
  ∀ (T G L Asp : Type) [DecidableEq T] [DecidableEq G] [DecidableEq L] [DecidableEq Asp] (Γ : GoGraph T L)
    (Ls : Finset L) (aspectOf : T → Option Asp) (o : Asp) (A : Set (T × G)),
    propagate Γ Ls aspectOf o (propagate Γ Ls aspectOf o A) = propagate Γ Ls aspectOf o A

/-- (3b) Monotone in the annotations. Hypothesis `A ⊆ B`. -/
def monotone_stmt : Prop :=
  ∀ (T G L Asp : Type) [DecidableEq T] [DecidableEq G] [DecidableEq L] [DecidableEq Asp] (Γ : GoGraph T L)
    (Ls : Finset L) (aspectOf : T → Option Asp) (o : Asp) (A B : Set (T × G)),
    A ⊆ B → propagate Γ Ls aspectOf o A ⊆ propagate Γ Ls aspectOf o B

/-- (3b) Extensive on what is kept: every kept direct annotation stays. No hypothesis. -/
def extensive_kept_stmt : Prop :=
  ∀ (T G L Asp : Type) [DecidableEq T] [DecidableEq G] [DecidableEq L] [DecidableEq Asp] (Γ : GoGraph T L)
    (Ls : Finset L) (aspectOf : T → Option Asp) (o : Asp) (A : Set (T × G)),
    kept aspectOf o A ⊆ propagate Γ Ls aspectOf o A

/-- (3b) Not extensive on the input: a two-term instance where a direct annotation of a term outside the aspect is
dropped, so `A ⊄ propagate A`. -/
def not_extensive_stmt : Prop :=
  ∃ (Γ : GoGraph ℕ ℕ) (Ls : Finset ℕ) (aspectOf : ℕ → Option ℕ) (o : ℕ) (A : Set (ℕ × ℕ)),
    ¬ A ⊆ propagate Γ Ls aspectOf o A

/-- (3c) Labels: following more relation labels can only add pairs (the two-relation set inside the five-relation
set, `goPropagationRelations`). Hypothesis `Ls ⊆ Ls'`. -/
def labels_mono_stmt : Prop :=
  ∀ (T G L Asp : Type) [DecidableEq T] [DecidableEq G] [DecidableEq L] [DecidableEq Asp] (Γ : GoGraph T L)
    (Ls Ls' : Finset L) (aspectOf : T → Option Asp) (o : Asp) (A : Set (T × G)),
    Ls ⊆ Ls' → propagate Γ Ls aspectOf o A ⊆ propagate Γ Ls' aspectOf o A

/-- (3d) Roots removed: no pair of `build` has the root; for any other term, `build` and `propagate` agree. No
hypothesis. -/
def root_removed_stmt : Prop :=
  ∀ (T G L Asp : Type) [DecidableEq T] [DecidableEq G] [DecidableEq L] [DecidableEq Asp] (Γ : GoGraph T L)
    (Ls : Finset L) (aspectOf : T → Option Asp) (o : Asp) (root : T) (A : Set (T × G)),
    (∀ x ∈ build Γ Ls aspectOf o root A, x.1 ≠ root) ∧
    (∀ (t : T) (g : G), t ≠ root →
      ((t, g) ∈ build Γ Ls aspectOf o root A ↔ (t, g) ∈ propagate Γ Ls aspectOf o A))

/-- (3d) `all` is never a term of a propagated pair: it has no aspect and no followed edge leads to it.
Hypothesis `aspectOf Γ.all = none` (`goDbTermAspects` keeps BP, MF and CC only, `:52-56`). -/
def all_absent_stmt : Prop :=
  ∀ (T G L Asp : Type) [DecidableEq T] [DecidableEq G] [DecidableEq L] [DecidableEq Asp] (Γ : GoGraph T L)
    (Ls : Finset L) (aspectOf : T → Option Asp) (o : Asp) (A : Set (T × G)),
    aspectOf Γ.all = none → ∀ g : G, (Γ.all, g) ∉ propagate Γ Ls aspectOf o A

/-- (3d) Of `goAncestorPairs` alone: `(t, all)` is a pair only when `t = all` is itself a start term. No
hypothesis. -/
def bfs_all_stmt : Prop :=
  ∀ (T L : Type) [DecidableEq T] [DecidableEq L] (Γ : GoGraph T L) (Ls : Finset L) (terms : Finset T) (t : T),
    (t, Γ.all) ∈ bfs Γ Ls terms → t = Γ.all ∧ Γ.all ∈ terms

/-- (3e) Aspect: if no followed edge joins two aspects (KIT-173 2.2's table, "Edges to another aspect: none", read in
GO.db 3.23.1), every ancestor of a term is in its aspect. Hypothesis: every edge with an allowed label and a parent
other than `all` has parent and child of the same aspect. -/
def aspect_closed_stmt : Prop :=
  ∀ (T L Asp : Type) [DecidableEq T] [DecidableEq L] [DecidableEq Asp] (Γ : GoGraph T L) (Ls : Finset L)
    (aspectOf : T → Option Asp),
    (∀ c p l, (c, p, l) ∈ Γ.edges → l ∈ Ls → p ≠ Γ.all → aspectOf p = aspectOf c) →
    ∀ t a, anc Γ Ls t a → aspectOf a = aspectOf t

/-! ## Witnesses (section 5.1) -/

/-- Witness for `monotone_stmt`: the inclusion can be strict. -/
def monotone_witness_stmt : Prop :=
  ∃ (Γ : GoGraph ℕ ℕ) (Ls : Finset ℕ) (aspectOf : ℕ → Option ℕ) (o : ℕ) (A B : Set (ℕ × ℕ)),
    A ⊆ B ∧ propagate Γ Ls aspectOf o A ≠ propagate Γ Ls aspectOf o B

/-- Witness for `labels_mono_stmt`: a `regulates`-like label adds pairs. -/
def labels_mono_witness_stmt : Prop :=
  ∃ (Γ : GoGraph ℕ ℕ) (Ls Ls' : Finset ℕ) (aspectOf : ℕ → Option ℕ) (o : ℕ) (A : Set (ℕ × ℕ)),
    Ls ⊆ Ls' ∧ propagate Γ Ls aspectOf o A ≠ propagate Γ Ls' aspectOf o A

/-- Witness for `all_absent_stmt`: `all` has no aspect, an edge to `all` exists, and some pair is propagated. -/
def all_absent_witness_stmt : Prop :=
  ∃ (Γ : GoGraph ℕ ℕ) (Ls : Finset ℕ) (aspectOf : ℕ → Option ℕ) (o : ℕ) (A : Set (ℕ × ℕ)),
    aspectOf Γ.all = none ∧ (∃ c l, (c, Γ.all, l) ∈ Γ.edges ∧ l ∈ Ls) ∧ (propagate Γ Ls aspectOf o A).Nonempty

/-- Witness for `aspect_closed_stmt`: the hypothesis holds on a graph with a followed edge and a term of another
aspect. -/
def aspect_closed_witness_stmt : Prop :=
  ∃ (Γ : GoGraph ℕ ℕ) (Ls : Finset ℕ) (aspectOf : ℕ → Option ℕ),
    (∀ c p l, (c, p, l) ∈ Γ.edges → l ∈ Ls → p ≠ Γ.all → aspectOf p = aspectOf c) ∧
    (∃ c p l, (c, p, l) ∈ Γ.edges ∧ l ∈ Ls ∧ p ≠ Γ.all ∧ c ≠ p) ∧
    (∃ t u, aspectOf t ≠ aspectOf u)

end MultiScholaRProofs.GoClosure
