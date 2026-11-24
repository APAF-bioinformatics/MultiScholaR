module

public import MultiScholaRProofs.Spec.QvalueBh
import Mathlib.Tactic

/-!
# Target 4 (KIT-193 section 6.4): a counterexample to a frozen statement

`bh_sorted_eq_stmt` (R's sort-and-cumulative-minimum transcription `bhSorted` equals the closed form `bh`, for every
`p`) is false as frozen: it carries no hypothesis on the sign of `p`, and with negative tied values the later member of
a tie group gets a smaller cumulative minimum than the closed form. Instance: `K = ℚ`, `m = 2`, `p = (-1, -1)`;
whatever order the sort gives the two tied indices, the second has `bhSorted = -2` and `bh = -1`. Under ticket 5.4 the
item stops here and goes back to the statements gate; p-values are never negative, and the statement with
`∀ j, 0 ≤ p j` is the candidate for the re-review (not proved here).
-/

public section

namespace MultiScholaRProofs.QvalueBh

open Finset

/-- The counter-instance: two tied p-values of `-1`. -/
def tiedNeg : Fin 2 → ℚ := fun _ => -1

theorem bh_tiedNeg (i : Fin 2) : bh tiedNeg i = -1 := by
  have hcnt : cntLe tiedNeg (-1) = 2 := by
    simp [cntLe, tiedNeg]
  have hterm : ∀ j, bhTerm tiedNeg j = -1 := by
    intro j
    simp only [bhTerm, tiedNeg, hcnt]
    norm_num
  have hfold : ((univ.filter fun j => tiedNeg i ≤ tiedNeg j).fold min (bhTerm tiedNeg i) (bhTerm tiedNeg)) = -1 := by
    apply le_antisymm
    · exact (Finset.fold_min_le _).2 (Or.inl (le_of_eq (hterm i)))
    · exact (Finset.le_fold_min _).2 ⟨le_of_eq (hterm i).symm, fun j _ => le_of_eq (hterm j).symm⟩
  unfold bh
  rw [hfold]
  norm_num

/-- `bh_sorted_eq_stmt` is false (a result of the proof attempt; ticket 5.4). -/
theorem bh_sorted_eq_counterexample : ¬ bh_sorted_eq_stmt := by
  intro h
  have heq := h ℚ 2 tiedNeg
  have hperm := List.mergeSort_perm (List.finRange 2) (fun a b => decide (tiedNeg b ≤ tiedNeg a))
  have hlen : ((List.finRange 2).mergeSort fun a b => decide (tiedNeg b ≤ tiedNeg a)).length = 2 := by
    rw [hperm.length_eq]
    simp
  obtain ⟨a, b, hab_eq⟩ :
      ∃ a b : Fin 2, ((List.finRange 2).mergeSort fun a b => decide (tiedNeg b ≤ tiedNeg a)) = [a, b] := by
    match hm : ((List.finRange 2).mergeSort fun a b => decide (tiedNeg b ≤ tiedNeg a)), hlen with
    | [a, b], _ => exact ⟨a, b, rfl⟩
  have hnd : [a, b].Nodup := by
    rw [← hab_eq]
    exact hperm.nodup_iff.2 (List.nodup_finRange 2)
  have hab : a ≠ b := by
    intro hab
    subst hab
    simp at hnd
  have hval : bhSorted tiedNeg b = -2 := by
    unfold bhSorted
    rw [hab_eq]
    have hidx : [a, b].idxOf b = 1 := by
      simp [List.idxOf_cons_ne _ hab]
    rw [hidx]
    simp [cumMin, tiedNeg]
    norm_num
  have := congrFun heq b
  rw [hval, bh_tiedNeg] at this
  norm_num at this

end MultiScholaRProofs.QvalueBh
