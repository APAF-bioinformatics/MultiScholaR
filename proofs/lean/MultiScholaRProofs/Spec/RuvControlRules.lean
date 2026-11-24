module

public import Mathlib.Algebra.Order.Field.Rat
public import Mathlib.Data.Fintype.Card
public import Mathlib.Data.List.Sort
public import Mathlib.Data.Rat.Floor

/-!
# Group R, R-N1: selection specifications (KIT-193 section 7.2)

Source: KIT-193 7.2 row R-N1 (audit U4); code read on `v0.9-staging` at 19997323:
* the sweep's control selection, `getNegCtrlProtAnovaHelper` (`R/func_prot_norm_optimization_helpers.R:463` onward):
  missing p-values set to 1 (`:565`), q-values (`:570`, `:572`), `filtered_list <- aov[aov > ruv_qval_cutoff]`
  (`:578`), `list_size <- ifelse(num_neg_ctrl > length(filtered_list), length(filtered_list), num_neg_ctrl)`
  (`:580`: the smaller of the two), `names(sort(filtered_list, decreasing = TRUE)[1:list_size])` (`:582`). With every
  q equal (VAL-82: the grouping lost, every p and q 1) the selection is positional: the first rows. This is the rule
  of every caller at 19997323. KIT-190 (VAL-110, VAL-126) gives lipidomics and metabolomics another rule (ranked by
  p-value, ties by feature identifier, features observed in one group left out), which this file does not model; the
  rule modelled here stays for proteomics, and for the replay of lipidomics and metabolomics records made before
  KIT-190.
* the first-pass controls of proteomics, `protRuvFirstPassControls` (`:1380`):
  `rank(first, ties.method = "first", na.last = "keep") > 0.25 * sum(is.finite(first))`, missing ranks `FALSE`; and of
  lipidomics and metabolomics, `omicsRuvFirstPassControls` (`R/func_omics_ruv_factors.R:54`, the same expression on
  the smallest p-value over the comparisons, among the features whose coefficients are estimable).
* the first-pass controls of transcriptomics, `trxRuvDesignControls` (`R/func_trx_ruv.R:296-297`):
  `genes[sort(order(first_p)[-seq_len(ceiling(0.25 * length(genes)))])]`: the first `⌈G/4⌉` of all `G` genes in
  `order()` of the p-values (missing last) are left out, a different rule (audit U4: an omics-parity question).

Model. P-values are rationals, a missing one is `none` (R's `NA`; `is.finite` is `FALSE` exactly there in the model).
R's `rank(ties.method = "first", na.last = "keep")` is the count model `rankFirst`; `order()` (missing last, ties in
index order) is `orderPos`. Two assumptions about R, checked by the paired test: `sort` and `order` on doubles are
stable here (method `radix`, which `auto` chooses below 2^31 values); and `x[1:n]` for `0 ≤ n ≤ length(x)` is
`rHead` (R's `1:0` is `c(1, 0)`, which selects the first element, so `list_size = 0` with a non-empty list selects
one control). Row names are distinct (the code returns `rownames(data_matrix) %in% control_genes`). A missing q-value
is outside the model (the code sets missing p-values to 1 before `qvalue`, `:565`).
-/

@[expose] public section

namespace MultiScholaRProofs.RuvControls

open Finset

/-- R's `x[1:n]` for `0 ≤ n ≤ length(x)`. -/
def rHead {α : Type} (n : ℕ) (l : List α) : List α := if n = 0 then l.take 1 else l.take n

/-- The sweep's controls (`:578-582`): the features with `q` above the cutoff, sorted by `q` decreasing (stable), the
first `min(num, length)`. -/
def sweepControls {p : ℕ} (q : Fin p → ℚ) (cutoff : ℚ) (num : ℕ) : List (Fin p) :=
  rHead (min num ((List.finRange p).filter fun i => decide (cutoff < q i)).length)
    (((List.finRange p).filter fun i => decide (cutoff < q i)).mergeSort fun i j => decide (q j ≤ q i))

/-- R-N1 (a) (kind S). With every `q` equal, the sweep selects by position: the first `min(num, #kept)` kept rows
(all rows when the common `q` is above the cutoff: VAL-82's positional selection). Hypothesis: every `q` equal. -/
def N1_stable_stmt : Prop :=
  ∀ (p : ℕ) (q : Fin p → ℚ) (cutoff : ℚ) (num : ℕ), (∀ i j, q i = q j) →
    sweepControls q cutoff num =
      rHead (min num ((List.finRange p).filter fun i => decide (cutoff < q i)).length)
        ((List.finRange p).filter fun i => decide (cutoff < q i))

/-- R-N1 (a'): ties in general. The sweep's order is "larger q first, then earlier row", so inside a block of tied
q-values the table's row order decides (VAL-126). No hypothesis. -/
def N1_ties_stmt : Prop :=
  ∀ (p : ℕ) (q : Fin p → ℚ) (cutoff : ℚ),
    (((List.finRange p).filter fun i => decide (cutoff < q i)).mergeSort fun i j => decide (q j ≤ q i)).Pairwise
      fun i j => q j < q i ∨ (q j = q i ∧ i < j)

/-- The number of p-values present. -/
def nPresent {F : ℕ} (pv : Fin F → Option ℚ) : ℕ := (univ.filter fun i => (pv i).isSome).card

/-- `j` comes before `i` among the present p-values: smaller p, or equal p and earlier row. -/
def before {F : ℕ} (pv : Fin F → Option ℚ) (j i : Fin F) : Prop :=
  match pv j, pv i with
  | some a, some b => a < b ∨ (a = b ∧ j < i)
  | _, _ => False

instance {F : ℕ} (pv : Fin F → Option ℚ) (j i : Fin F) : Decidable (before pv j i) := by
  unfold before
  split <;> infer_instance

/-- R's `rank(pv, ties.method = "first", na.last = "keep")` at a present p-value. -/
def rankFirst {F : ℕ} (pv : Fin F → Option ℚ) (i : Fin F) : ℕ := 1 + (univ.filter fun j => before pv j i).card

/-- The first-pass controls of proteomics, lipidomics and metabolomics: present and of rank above `F/4`. -/
def firstPassControls {F : ℕ} (pv : Fin F → Option ℚ) : Finset (Fin F) :=
  univ.filter fun i => (pv i).isSome ∧ (1 / 4 : ℚ) * nPresent pv < rankFirst pv i

/-- R-N1 (b) (kind S). The first-pass controls are exactly the present p-values but the `⌊F/4⌋` most significant
(in the order "smaller p, then earlier row"): there are `F - ⌊F/4⌋` of them, and each present p-value left out comes
before each control. No hypothesis. -/
def N1_first_pass_stmt : Prop :=
  ∀ (F : ℕ) (pv : Fin F → Option ℚ),
    (firstPassControls pv).card = nPresent pv - nPresent pv / 4 ∧
      ∀ i j, (pv i).isSome → i ∉ firstPassControls pv → j ∈ firstPassControls pv → before pv i j

/-- `j` comes before `i` in R's `order()`: missing last, smaller p first, ties in row order. -/
def orderBefore {G : ℕ} (pv : Fin G → Option ℚ) (j i : Fin G) : Prop :=
  ((pv j).isSome = true ∧ pv i = none) ∨ before pv j i ∨ (pv j = none ∧ pv i = none ∧ j < i)

instance {G : ℕ} (pv : Fin G → Option ℚ) (j i : Fin G) : Decidable (orderBefore pv j i) := by
  unfold orderBefore
  infer_instance

/-- The position of `i` in R's `order(pv)` (1-based). -/
def orderPos {G : ℕ} (pv : Fin G → Option ℚ) (i : Fin G) : ℕ := 1 + (univ.filter fun j => orderBefore pv j i).card

/-- The first-pass controls of transcriptomics: every gene after the first `⌈G/4⌉` in `order()`. -/
def trxFirstPassControls {G : ℕ} (pv : Fin G → Option ℚ) : Finset (Fin G) :=
  univ.filter fun i => ⌈(1 / 4 : ℚ) * G⌉₊ < orderPos pv i

/-- R-N1 (c) (kind S). Transcriptomics leaves out the first `⌈G/4⌉` of all `G` genes in `order()`: there are
`G - ⌈G/4⌉` controls, a gene without a p-value can be one, and each gene left out comes before each control. No
hypothesis. -/
def N1_trx_stmt : Prop :=
  ∀ (G : ℕ) (pv : Fin G → Option ℚ),
    (trxFirstPassControls pv).card = G - ⌈(1 / 4 : ℚ) * G⌉₊ ∧
      ∀ i j, i ∉ trxFirstPassControls pv → j ∈ trxFirstPassControls pv → orderBefore pv i j

/-- R-N1 (d): the two rules differ (audit U4) in two ways. With five genes all with a p-value the first leaves out
`⌊5/4⌋ = 1` and the second `⌈5/4⌉ = 2`. With four genes of which the first has no p-value, that gene is a
transcriptomics control and is not a first-pass control of the other omics. Closed statement. -/
def N1_rules_differ_stmt : Prop :=
  (let pv : Fin 5 → Option ℚ := fun i => some ((i : ℕ) : ℚ)
   firstPassControls pv ≠ trxFirstPassControls pv) ∧
  (let pv : Fin 4 → Option ℚ := fun i => if (i : ℕ) = 0 then none else some ((i : ℕ) : ℚ)
   (0 : Fin 4) ∈ trxFirstPassControls pv ∧ (0 : Fin 4) ∉ firstPassControls pv)

/-! ## Witnesses (section 5.1) -/

/-- Witness for `N1_stable_stmt`: every `q` equal and above the cutoff, more rows than selected. -/
def N1_stable_witness_stmt : Prop :=
  ∃ (p : ℕ) (q : Fin p → ℚ) (cutoff : ℚ) (num : ℕ), (∀ i j, q i = q j) ∧ (∀ i, cutoff < q i) ∧ 0 < num ∧ num < p

end MultiScholaRProofs.RuvControls
