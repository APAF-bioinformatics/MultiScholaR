module

public import Mathlib.Data.List.Basic
public import Mathlib.Data.List.Dedup
public import Mathlib.Logic.Function.Basic

/-!
# Target 2: replay-rule selection (KIT-193 section 6.2)

Source: KIT-187 r1.3 section 7.3 and its code on `kit-187-filter-by-count` at e6c42d1d
(`R/utils_workflow_analysis_record.R:421-432` the builder's registry fallback, `:644-693` the marker table,
`workflowReplayRule`, `workflowReplayRuleOptions`, `workflowFixedInForce`, `:785` the per-step wrapper;
`R/utils_workflow_dia_state.R:729` `setWorkflowStepProvenance` applies `workflowFixedInForce` for every writer);
KIT-190 r1.2 section 7.4 (one more table row,
`ruv_control_grouping = list(option = "MultiScholaR.ruv_control_grouping", current = "kept", earlier = "dropped")`,
and the managed writer's removal of that marker from a new-design step, `R/utils_workflow_execution.R:3327-3332`
as KIT-190 cites it; read as ticket text, KIT-190's branch does not yet carry KIT-187's table).

Model. Names `N`, option names `O` and recorded values `V` are abstract types with decidable equality. Values are
compared after `decodeWorkflowRecordValue` (the encode-decode round trip of a token is an assumption, checked by the
paired test). An option state is `O → Option V`, `none` for unset. `withr::with_options` restores options at the
end of each step, error or not (modelled as per-step scope); a single process (KIT-190 7.4). The marker table is a
list constant, bound by its own digest (section 8.2). That every operation whose kernel reads a marker's option
also declares that marker is a property of the code, not of this model (KIT-187 2.4 and T8, KIT-190 T7 (iv), (v)).
That a new-design step's kernel does not read a marker's option is a property of the code (KIT-190 T7 (v)); no
statement here depends on it.
-/

@[expose] public section

namespace MultiScholaRProofs.ReplayRule

/-- The source of a recorded parameter entry. R records further sources for declared parameters; the rule reads only
whether the source is `fixed`, so they are represented by `args` and `declared`. -/
inductive Source
  | fixed
  | registry
  | args
  | declared
  deriving DecidableEq

variable {N O V : Type} [DecidableEq N] [DecidableEq O] [DecidableEq V]

/-- One parameter entry of a recorded step: `list(name, source, value, …)`. -/
structure Entry (N V : Type) where
  name : N
  source : Source
  value : V

/-- A recorded step: its parameter entries, in record order. -/
abbrev Step (N V : Type) := List (Entry N V)

/-- One row of `workflowReplayRuleMarkers`: the marker's name, the option that selects the rule inside the package,
its current token (the declaration's fixed value) and its earlier token. -/
structure Marker (N O V : Type) where
  name : N
  option : O
  current : V
  earlier : V

/-- A declaration's (or a provenance's) fixed constants: a named list, as (name, value) pairs in order. -/
abbrev Fixed (N V : Type) := List (N × V)

/-- An option state: the value each option holds, `none` when unset. -/
abbrev Opts (O V : Type) := O → Option V

/-- The marker table's shape, as the constant has it (two rows after KIT-190): distinct names and distinct options.
A property of the constant, bound by its digest (section 8.2). -/
def TableOK (tbl : List (Marker N O V)) : Prop :=
  (tbl.map Marker.name).Nodup ∧ (tbl.map Marker.option).Nodup

/-- Whether an entry records marker `μ` as made under the current rule: name `μ`, source `fixed`, value current
(the three `identical()` tests of `:659-660`). -/
def recordsCurrent (μ : Marker N O V) (e : Entry N V) : Bool :=
  decide (e.name = μ.name) && decide (e.source = Source.fixed) && decide (e.value = μ.current)

/-- `workflowReplayRule(step, marker)` as written (`:657-664`): the first entry, in record order, that records `μ`
as current returns the current token; with none, the earlier token. -/
def ruleLoop (step : Step N V) (μ : Marker N O V) : V :=
  match step.find? (recordsCurrent μ) with
  | some _ => μ.current
  | none => μ.earlier

/-- The rule a step was made with, stated existentially: current if some entry records `μ` as current, earlier
otherwise. `rule_first_match_stmt` states that it is `ruleLoop`. -/
def rule (step : Step N V) (μ : Marker N O V) : V :=
  if ∃ e ∈ step, recordsCurrent μ e = true then μ.current else μ.earlier

/-- `workflowFixedInForce(fixed)` (`:686-693`): each table marker whose option holds its earlier token is dropped;
an unset option or any other token keeps it. (R returns `NULL` when nothing is left; the empty list here.) -/
def inForce (tbl : List (Marker N O V)) (fixed : Fixed N V) (opts : Opts O V) : Fixed N V :=
  fixed.filter fun x => decide (¬ ∃ μ ∈ tbl, μ.name = x.1 ∧ opts μ.option = some μ.earlier)

/-- What a writer records as the provenance's fixed constants for a step whose declaration has fixed constants
`decl`: `setWorkflowStepProvenance` applies `workflowFixedInForce` (`R/utils_workflow_dia_state.R:729`) in process,
on replay (`runSavedOperation`) and in the managed legacy writer (`completeWorkflowAction`). -/
def compute (tbl : List (Marker N O V)) (decl : Fixed N V) (opts : Opts O V) : Fixed N V :=
  inForce tbl decl opts

/-- The managed writer for a new-design step (KIT-190 7.4): the names in `dropped` (`ruv_control_grouping`) are
removed from the registry declaration's fixed constants before they reach `setWorkflowStepProvenance`. -/
def computeManagedNewDesign (tbl : List (Marker N O V)) (dropped : List N) (decl : Fixed N V)
    (opts : Opts O V) : Fixed N V :=
  inForce tbl (decl.filter fun x => decide (x.1 ∉ dropped)) opts

/-- The record builder's fixed entries (`:421-432`): the provenance's fixed constants with source `fixed` when there
are any; otherwise the registry declaration's fixed constants minus every table marker, with source `registry`
(a marker is never back-filled). Only these entries can have source `fixed` or `registry`; `others` are the step's
other entries, whatever their place: `rule` does not depend on order. -/
def build (tbl : List (Marker N O V)) (others : Step N V) (prov decl : Fixed N V) : Step N V :=
  others ++
    if prov ≠ [] then prov.map fun x => ⟨x.1, Source.fixed, x.2⟩
    else (decl.filter fun x => decide (x.1 ∉ tbl.map Marker.name)).map fun x => ⟨x.1, Source.registry, x.2⟩

/-- The rows of the table that a declaration fixes, in the declaration's order without repeats
(`intersect(names(declaration$fixed), names(workflowReplayRuleMarkers))`, `:677`). -/
def markersFixed (tbl : List (Marker N O V)) (decl : Fixed N V) : List (Marker N O V) :=
  (decl.map Prod.fst).dedup.filterMap fun nm => tbl.find? fun μ => decide (μ.name = nm)

/-- `workflowReplayRuleOptions(step, declaration)` (`:675-680`): one option per table marker the declaration fixes,
set to the step's own rule. -/
def replayOptions (tbl : List (Marker N O V)) (step : Step N V) (decl : Fixed N V) : List (O × V) :=
  (markersFixed tbl decl).map fun μ => (μ.option, rule step μ)

/-- `withr::with_options(list, …)` over an ambient option state: the listed options set, every other one as it
was. R's `options[[key]] <- value` keeps the last value for a repeated key and `List.lookup` the first; under
`TableOK` no key repeats. -/
def withOpts (ambient : Opts O V) (set : List (O × V)) : Opts O V :=
  fun o => match set.lookup o with
    | some v => some v
    | none => ambient o

/-- The rule the helper applies under an option state: it reads `getOption(option, current)`, and only the earlier
token selects the earlier rule (KIT-187 7.3). -/
def helperRule (opts : Opts O V) (μ : Marker N O V) : V :=
  if opts μ.option = some μ.earlier then μ.earlier else μ.current

/-- `ruleLoop` (the R loop, first match) equals the existential `rule`. Light specification lemma. -/
def rule_first_match_stmt : Prop :=
  ∀ (N O V : Type) [DecidableEq N] [DecidableEq O] [DecidableEq V] (step : Step N V) (μ : Marker N O V),
    ruleLoop step μ = rule step μ

/-- (2a) Totality: the rule of any step for any marker is its current or its earlier token; and the replay options
name exactly the options of the table markers the declaration fixes. Hypothesis `TableOK tbl` (distinct names: a
marker name the declaration fixes finds its own row). -/
def totality_stmt : Prop :=
  ∀ (N O V : Type) [DecidableEq N] [DecidableEq O] [DecidableEq V] (tbl : List (Marker N O V))
    (step : Step N V) (decl : Fixed N V),
    TableOK tbl →
      (∀ μ : Marker N O V, rule step μ = μ.current ∨ rule step μ = μ.earlier) ∧
      (∀ o : O, o ∈ (replayOptions tbl step decl).map Prod.fst ↔
        ∃ μ ∈ tbl, μ.option = o ∧ μ.name ∈ decl.map Prod.fst)

/-- (2b) Faithfulness: for a declaration that fixes marker `μ` with its current token, and every option state `o`
(unset, or any token), the step every writer of `compute` records, through the builder, has the rule the helper
applies under `o`: earlier if `o` holds the earlier token, current otherwise.
Hypotheses: `TableOK tbl`; `μ ∈ tbl`; `(μ.name, μ.current) ∈ decl` (the declaration fixes the marker with its current
token, as the table's definition says); every entry of `others` has a source other than `fixed` and `registry`
(declared parameters and arguments; the builder writes `fixed`/`registry` entries only from the fixed constants). -/
def faithful_stmt : Prop :=
  ∀ (N O V : Type) [DecidableEq N] [DecidableEq O] [DecidableEq V] (tbl : List (Marker N O V))
    (μ : Marker N O V) (others : Step N V) (decl : Fixed N V) (o : Opts O V),
    TableOK tbl → μ ∈ tbl → (μ.name, μ.current) ∈ decl →
    (∀ e ∈ others, e.source ≠ Source.fixed ∧ e.source ≠ Source.registry) →
      rule (build tbl others (compute tbl decl o) decl) μ = helperRule o μ

/-- (2b, legacy steps) and (2e, new-design steps): a step whose provenance carries no fixed entry named `μ` (a step
made before the marker existed, or a new-design step) records the earlier rule, whether the builder lists the
provenance or falls back to the registry (whose entries have source `registry`, which the rule never counts).
Hypotheses: the provenance has no entry named `μ`; `others` as in `faithful_stmt`. -/
def without_marker_stmt : Prop :=
  ∀ (N O V : Type) [DecidableEq N] [DecidableEq O] [DecidableEq V] (tbl : List (Marker N O V))
    (μ : Marker N O V) (others : Step N V) (prov decl : Fixed N V),
    μ.name ∉ prov.map Prod.fst →
    (∀ e ∈ others, e.source ≠ Source.fixed ∧ e.source ≠ Source.registry) →
      rule (build tbl others prov decl) μ = μ.earlier

/-- (2e) The managed writer's new-design step records no marker named in `dropped`, so replays under the earlier
token, whatever the option state. Hypotheses: `μ.name ∈ dropped`; `others` as above. -/
def managed_new_design_stmt : Prop :=
  ∀ (N O V : Type) [DecidableEq N] [DecidableEq O] [DecidableEq V] (tbl : List (Marker N O V))
    (μ : Marker N O V) (dropped : List N) (others : Step N V) (decl : Fixed N V) (o : Opts O V),
    μ.name ∈ dropped →
    (∀ e ∈ others, e.source ≠ Source.fixed ∧ e.source ≠ Source.registry) →
      rule (build tbl others (computeManagedNewDesign tbl dropped decl o) decl) μ = μ.earlier

/-- (2c) Round trip: a step replayed under its replay options records a step whose rule equals the original's, for
every marker the declaration fixes with its current token, every ambient option state, and every writer of
`compute`. Hypotheses: `TableOK tbl`; `μ ∈ tbl`; `(μ.name, μ.current) ∈ decl`; `others` as in `faithful_stmt`. -/
def round_trip_stmt : Prop :=
  ∀ (N O V : Type) [DecidableEq N] [DecidableEq O] [DecidableEq V] (tbl : List (Marker N O V))
    (μ : Marker N O V) (step others : Step N V) (decl : Fixed N V) (ambient : Opts O V),
    TableOK tbl → μ ∈ tbl → (μ.name, μ.current) ∈ decl →
    (∀ e ∈ others, e.source ≠ Source.fixed ∧ e.source ≠ Source.registry) →
      rule (build tbl others (compute tbl decl (withOpts ambient (replayOptions tbl step decl))) decl) μ =
        rule step μ

/-- (2d) Locality: a step's rule for `μ` depends only on its entries named `μ` (so the two markers do not
interact). No hypothesis. -/
def locality_stmt : Prop :=
  ∀ (N O V : Type) [DecidableEq N] [DecidableEq O] [DecidableEq V] (step : Step N V) (μ : Marker N O V),
    rule step μ = rule (step.filter fun e => decide (e.name = μ.name)) μ

/-- (2d) A step whose declaration does not fix table marker `ν` sees, for `ν`'s option, the helper's default
(current), when that option is unset outside the replay (the option is internal and set nowhere else, KIT-187 7.3).
Hypotheses: `TableOK tbl`; `ν ∈ tbl`; `ν.name ∉ decl` names; the ambient option state leaves `ν.option` unset. -/
def unfixed_default_stmt : Prop :=
  ∀ (N O V : Type) [DecidableEq N] [DecidableEq O] [DecidableEq V] (tbl : List (Marker N O V))
    (ν : Marker N O V) (step : Step N V) (decl : Fixed N V) (ambient : Opts O V),
    TableOK tbl → ν ∈ tbl → ν.name ∉ decl.map Prod.fst → ambient ν.option = none →
      helperRule (withOpts ambient (replayOptions tbl step decl)) ν = ν.current

/-! ## Witnesses (section 5.1)

Each statement with hypotheses has a witness: its hypotheses hold together on an explicit instance with a table of
two rows (the table after KIT-190), distinct current and earlier tokens, and a declaration that fixes a marker. -/

/-- A non-degenerate marker table of two rows, with distinct tokens per row. -/
def TwoRowTable (tbl : List (Marker N O V)) : Prop :=
  TableOK tbl ∧ tbl.length = 2 ∧ ∀ μ ∈ tbl, μ.current ≠ μ.earlier

/-- Witness for `totality_stmt`, `faithful_stmt` and `round_trip_stmt` (their hypotheses together, with a
non-empty step and an `others` entry). -/
def faithful_witness_stmt : Prop :=
  ∃ (tbl : List (Marker String String String)) (μ : Marker String String String)
    (others : Step String String) (decl : Fixed String String),
    TwoRowTable tbl ∧ μ ∈ tbl ∧ (μ.name, μ.current) ∈ decl ∧ others ≠ [] ∧
    (∀ e ∈ others, e.source ≠ Source.fixed ∧ e.source ≠ Source.registry)

/-- Witness for `without_marker_stmt` (a provenance with an entry, not named `μ`). -/
def without_marker_witness_stmt : Prop :=
  ∃ (tbl : List (Marker String String String)) (μ : Marker String String String)
    (others : Step String String) (prov : Fixed String String),
    TwoRowTable tbl ∧ μ ∈ tbl ∧ prov ≠ [] ∧ μ.name ∉ prov.map Prod.fst ∧
    (∀ e ∈ others, e.source ≠ Source.fixed ∧ e.source ≠ Source.registry)

/-- Witness for `managed_new_design_stmt`. -/
def managed_new_design_witness_stmt : Prop :=
  ∃ (tbl : List (Marker String String String)) (μ : Marker String String String) (dropped : List String),
    TwoRowTable tbl ∧ μ ∈ tbl ∧ μ.name ∈ dropped

/-- Witness for `unfixed_default_stmt`. -/
def unfixed_default_witness_stmt : Prop :=
  ∃ (tbl : List (Marker String String String)) (ν : Marker String String String) (decl : Fixed String String),
    TwoRowTable tbl ∧ ν ∈ tbl ∧ decl ≠ [] ∧ ν.name ∉ decl.map Prod.fst

end MultiScholaRProofs.ReplayRule
