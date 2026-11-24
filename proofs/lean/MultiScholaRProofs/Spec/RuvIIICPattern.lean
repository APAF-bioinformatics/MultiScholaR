module

public import MultiScholaRProofs.Spec.RuvIIIC

/-!
# Group R, R-J3: RUV-III-C's missing-value pattern (KIT-193 section 7.2; draft for its statements review)

Not frozen. Written after the statements gate as its register (section 4) allows: this file is new, has its own line
in `SPEC.sha256`, and goes through its own statements review before any proof (ticket 5.3 item 4). No other Spec
file changes.

Source: KIT-193 7.2 row R-J3 and 5.3 item 1; the branch conditions as the statements register (section 4) read them
in RUVIIIC 1.0.19 at the locked commit c534151f, `src/RUVIIIC_Varying.cpp`: every value of the target missing when
`min(residual dimensions, controls complete on the target's runs) < k` (`:168-171`); the closed-form branch when the
residual dimension equals `k` (`:173-184`; R-K1 (ii)); otherwise the general branch (`:229-303`), every value missing
when fewer than `k` eigenvalues converge (`:252-263`) or the eigen step throws (`:304` onward); a corrected target is
written on its observed runs and missing elsewhere (`:186-198`, `:286-301`). The callers drop the features (and runs)
left with no value (`R/func_prot_s4_norm_methods.R:264-265`; lipid and metabolite `R/func_lipid_s4_ruv_methods.R:713-715`),
and the post-RUV clean-up applies KIT-187's count rule to the corrected matrix (`R/utils_workflow_dia_products.R:489-492`).

Model. `RuvIIICBr` adds to `RuvIIIC` the residual dimension on a row set (`currentResidualDimensions`, abstract), a
flag for the eigen step (at least `k` eigenvalues converged and no exception; abstract), and the closed-form branch's
values (abstract here; R-K1 (ii) states them). Values are reals: a non-finite value (from the closed form's
unchecked inverse of a singular `Y_C Y_C'`, audit U3, or from the eigen-solver) is outside the model, as is the
percentile threshold the clean-up recomputes on the corrected values (not covered, ticket 7.2).
-/

@[expose] public section

namespace MultiScholaRProofs.Ruv

open Matrix Finset

/-- RUV-III-C with its branches. -/
structure RuvIIICBr (n p s k : ℕ) extends RuvIIIC n p s k where
  resDim : Finset (Fin n) → ℕ
  eigenOK : Finset (Fin n) → Bool
  closedForm : Fin p → Fin n → ℝ

namespace RuvIIICBr

variable {n p s k : ℕ}

/-- Target `j` is corrected: at least `k` residual dimensions and `k` controls complete on its runs, and, in the
general branch (more than `k` residual dimensions), the eigen step succeeded. -/
def correctable (r : RuvIIICBr n p s k) (j : Fin p) : Bool :=
  decide (k ≤ r.resDim (r.rows j)) && decide (k ≤ (r.ctlOn (r.rows j)).card) &&
    (decide (r.resDim (r.rows j) = k) || r.eigenOK (r.rows j))

/-- The written value of target `j` at run `i`: missing unless `j` is corrected and observed at `i`; the closed form
when the residual dimension is `k`, else the general branch's `y_j - W_j a_j`. -/
noncomputable def out (r : RuvIIICBr n p s k) (j : Fin p) (i : Fin n) : Option ℝ :=
  if r.correctable j = true ∧ r.obs i j = true then
    some (if r.resDim (r.rows j) = k then r.closedForm j i
      else r.Y i j - (r.WOn (r.rows j) * r.aOn (r.rows j)) i j)
  else none

end RuvIIICBr

/-- R-J3 (kind S). A corrected feature has a value exactly on its observed runs; a feature that is not corrected has
no value at any run (so the caller drops it). No hypothesis. -/
def J3_pattern_stmt : Prop :=
  ∀ (n p s k : ℕ) (r : RuvIIICBr n p s k) (j : Fin p),
    (r.correctable j = true → ∀ i, (r.out j i).isSome = r.obs i j) ∧
      (r.correctable j = false → ∀ i, r.out j i = none)

/-- R-J3, the clean-up's counts: for every feature the caller keeps (one with a value at some run) and every set of
runs `g` (a group, on the runs that remain), the number of missing values after RUV equals the number before. No
hypothesis beyond the feature being kept. -/
def J3_counts_stmt : Prop :=
  ∀ (n p s k : ℕ) (r : RuvIIICBr n p s k) (j : Fin p) (g : Finset (Fin n)),
    (∃ i, (r.out j i).isSome = true) →
      (g.filter fun i => r.out j i = none).card = (g.filter fun i => r.obs i j = false).card

/-- Witness for `J3_counts_stmt`: a corrected feature with a missing run (so a count is not trivially 0) and a
feature that is not corrected. -/
def J3_witness_stmt : Prop :=
  ∃ (n p s k : ℕ) (r : RuvIIICBr n p s k) (j j' : Fin p),
    r.correctable j = true ∧ r.correctable j' = false ∧ (∃ i, r.obs i j = false) ∧ ∃ i, r.obs i j = true

end MultiScholaRProofs.Ruv
