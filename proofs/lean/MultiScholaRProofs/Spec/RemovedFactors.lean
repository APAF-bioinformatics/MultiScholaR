module

public import MultiScholaRProofs.Spec.RuvAlgebra
public import Mathlib.Data.Matrix.ColumnRowPartitioned

/-!
# Group R, R-I1: KIT-184's removed factors (KIT-193 section 7.2; draft for its statements review)

Not frozen. Written after the statements gate as its register (section 4) allows: this file is new, has its own line
in `SPEC.sha256`, and goes through its own statements review before any proof (ticket 5.3 item 4). No other Spec
file changes.

Source: KIT-193 7.2 row R-I1 and 7.1 row I1; KIT-184 r1.3 5.2 (what is equal to what), 5.3 (averaging runs per
biological sample, `integrationAverageRuns()`, `R/utils_workflow_integration_model.R:647` on `v0.9-staging`), 5.4
(`integrationInputRemovedFactors`: `F` the features finite in every run in both matrices; `D` the difference
uncorrected minus corrected on `F`, each feature centred over the runs; its right singular vectors, numerical rank at
`1e-8 · s_1`; status `exact` when the rank is at most `k`). KIT-184's writer is not merged: the statements below are
its algebra over the RUV-III model of `Spec/RuvAlgebra.lean` (runs × features, so KIT-184's rows of `D` are columns
here and its row space is a column space here); the numerical rank rule and the writer's status are not stated and
wait for KIT-184's merged code (statements register, section 4).

Which corrected values: RUV-III's `Y - Wα` (L1 and T1's integration input, not restored) and the centred form
`Y - W_c α` (P5, and by R-C1 the restored values of P1 to P3 and L2); on RUV-III-C a feature observed in every run
is corrected with the one `W` of the complete features (R-J1), so the statements apply to the features `F` with that
`W`. Hypotheses, each with its reason: `Hac` (the code's refusal, `solve()` on a singular matrix); `HU` (`R_M U = U`,
as in R-G0); `Hone` (every run in exactly one replicate set, so `1 ∈ col(M)` and, with `HU`, `U'1 = 0`); for the
equality of spans, `rank α_F = k` (the removed directions are all seen on `F`); for the averaging, an averaging map
`A` (biological samples × runs) whose rows sum to 1 (`integrationAverageRuns()` averages each sample's runs).
-/

@[expose] public section

namespace MultiScholaRProofs.Ruv

open Matrix RuvIII

/-- R-I1 (a): `U'W_c = I`, so the centred factors have rank `k`. Hypotheses `Hac`, `HU`, `Hone`. -/
def I1_UtWc_stmt : Prop :=
  ∀ (n p s k : ℕ) (r : RuvIII n p s k), Hac r → HU r → Hone r → r.Uᵀ * r.Wc = 1 ∧ r.Wc.rank = k

/-- R-I1 (b): the difference uncorrected minus corrected, each feature centred over the runs, is `W_c α`, for the
corrected values `Y - Wα` and for `Y - W_c α` alike. No hypothesis: meant for the junk inverse too. -/
def I1_centred_diff_stmt : Prop :=
  ∀ (n p s k : ℕ) (r : RuvIII n p s k),
    centre (r.Y - r.Ystar) = r.Wc * r.alpha ∧ centre (r.Y - r.adjusted) = r.Wc * r.alpha

/-- R-I1 (c): on features `F`, the span of the centred difference `W_c α_F` (KIT-184's row space of `D`; its right
singular vectors when the status is `exact`) lies in `col(W_c)`, with equality when `rank α_F = k`. The inclusion
has no hypothesis; the equality has `rank α_F = k`. -/
def I1_span_stmt : Prop :=
  ∀ (n p s k : ℕ) (r : RuvIII n p s k) (F : Finset (Fin p)),
    LinearMap.range (r.Wc * r.alpha.submatrix id (fun j : F => (j : Fin p))).mulVecLin ≤
        LinearMap.range r.Wc.mulVecLin ∧
      ((r.alpha.submatrix id (fun j : F => (j : Fin p))).rank = k →
        LinearMap.range (r.Wc * r.alpha.submatrix id (fun j : F => (j : Fin p))).mulVecLin =
          LinearMap.range r.Wc.mulVecLin)

/-- R-I1 (d): averaging runs per biological sample (a linear map `A`, rows summing to 1) keeps the removed part of
every feature in `span[1, A W_c]` (KIT-184 5.3's `W-bar`), for both forms of the corrected values. Hypothesis: the
rows of `A` sum to 1. -/
def I1_average_stmt : Prop :=
  ∀ (n p s k d : ℕ) (r : RuvIII n p s k) (A : Matrix (Fin d) (Fin n) ℝ),
    A *ᵥ (fun _ => (1 : ℝ)) = (fun _ => 1) → ∀ j : Fin p,
      A *ᵥ (fun i => (r.Y - r.Ystar) i j) ∈
          LinearMap.range (fromCols (Matrix.of fun (_ : Fin d) (_ : Fin 1) => (1 : ℝ)) (A * r.Wc)).mulVecLin ∧
        A *ᵥ (fun i => (r.Y - r.adjusted) i j) ∈
          LinearMap.range (fromCols (Matrix.of fun (_ : Fin d) (_ : Fin 1) => (1 : ℝ)) (A * r.Wc)).mulVecLin

/-! ## Witnesses (section 5.1) -/

/-- Witness for `I1_UtWc_stmt` and the equality of `I1_span_stmt`: every hypothesis together, `0 < k`, a set of
features on which `α` has rank `k`, and at least one more run than sets. -/
def I1_witness_stmt : Prop :=
  ∃ (n p s k : ℕ) (r : RuvIII n p s k) (F : Finset (Fin p)),
    0 < k ∧ s < n ∧ Hac r ∧ HU r ∧ Hone r ∧ (r.alpha.submatrix id (fun j : F => (j : Fin p))).rank = k

/-- Witness for `I1_average_stmt`: an averaging map with fewer samples than runs. -/
def I1_average_witness_stmt : Prop :=
  ∃ (d n : ℕ) (A : Matrix (Fin d) (Fin n) ℝ), 0 < d ∧ d < n ∧ A *ᵥ (fun _ => (1 : ℝ)) = fun _ => 1

end MultiScholaRProofs.Ruv
