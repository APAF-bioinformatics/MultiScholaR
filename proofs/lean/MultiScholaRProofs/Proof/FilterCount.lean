module

public import MultiScholaRProofs.Spec.FilterCount
import Mathlib.Tactic

/-!
# Proofs for target 1 (KIT-193 section 6.1)

Theorems whose types are the frozen statements of `Spec/FilterCount.lean`. The exactness within the lcm bound rests
on one integer fact: when `k/n > a/n*`, the difference `k n* - a n` is a positive multiple of `gcd(n, n*)`, so it is at
least `gcd(n, n*)`, and `gcd · lcm = n n*` with `lcm < 100000` leaves the rounding slack of `P` no room.
-/

public section

namespace MultiScholaRProofs.FilterCount

theorem roundedOK_bounds {P : ℤ} {s : ℚ} (h : RoundedOK P s) :
    100000 * s - 1 / 2 ≤ (P : ℚ) ∧ (P : ℚ) ≤ 100000 * s + 1 / 2 := by
  unfold RoundedOK at h
  rw [abs_le] at h
  constructor <;> linarith [h.1, h.2]

theorem failsRule_iff_rat (P : ℤ) (n k : ℕ) :
    failsRule P n k ↔ (n : ℚ) * (2 * P + 1) < 200000 * k := by
  unfold failsRule
  constructor
  · intro h
    exact_mod_cast h
  · intro h
    exact_mod_cast h

/-- (6.1 (a)) never stricter. -/
theorem never_stricter : never_stricter_stmt := by
  intro s P n k hn hP hfail
  obtain ⟨hlo, -⟩ := roundedOK_bounds hP
  rw [failsRule_iff_rat] at hfail
  unfold failsIntended
  have hn' : (0 : ℚ) < n := by exact_mod_cast hn
  rw [gt_iff_lt, lt_div_iff₀ hn']
  nlinarith

/-- The exactness core, for a share `a / n*` with any numerator `a`. -/
theorem exact_core (a nstar n k : ℕ) (P : ℤ) (hns : 0 < nstar) (hn : 0 < n)
    (hP : RoundedOK P ((a : ℚ) / nstar)) (hl : Nat.lcm n nstar < 100000) :
    failsRule P n k ↔ failsIntended ((a : ℚ) / nstar) n k := by
  constructor
  · exact never_stricter _ P n k hn hP
  · intro hint
    have hns' : (0 : ℚ) < nstar := by exact_mod_cast hns
    have hn' : (0 : ℚ) < n := by exact_mod_cast hn
    -- the intended rule in integers: a n < k n*
    have hcross : a * n < k * nstar := by
      unfold failsIntended at hint
      rw [gt_iff_lt, div_lt_div_iff₀ hns' hn'] at hint
      exact_mod_cast hint
    -- k n* - a n is a positive multiple of g = gcd(n, n*)
    set g := Nat.gcd n nstar with hg
    have hgpos : 0 < g := Nat.gcd_pos_of_pos_left _ hn
    have hdvd : g ∣ k * nstar - a * n :=
      Nat.dvd_sub (dvd_mul_of_dvd_right (Nat.gcd_dvd_right n nstar) k)
        (dvd_mul_of_dvd_right (Nat.gcd_dvd_left n nstar) a)
    have hge : g ≤ k * nstar - a * n := Nat.le_of_dvd (by omega) hdvd
    have hgl : g * Nat.lcm n nstar = n * nstar := Nat.gcd_mul_lcm n nstar
    -- 2 n n* < 200000 (k n* - a n), in ℕ
    have hkey : 2 * (n * nstar) < 200000 * (k * nstar - a * n) := by
      rw [← hgl]
      have : g * Nat.lcm n nstar < g * 100000 := Nat.mul_lt_mul_of_pos_left hl hgpos
      nlinarith
    have hkeyQ : (2 : ℚ) * (n * nstar) < 200000 * ((k : ℚ) * nstar - a * n) := by
      have h1 : ((k * nstar - a * n : ℕ) : ℚ) = (k : ℚ) * nstar - a * n := by
        rw [Nat.cast_sub hcross.le]
        push_cast
        ring
      have h2 : ((2 * (n * nstar) : ℕ) : ℚ) < ((200000 * (k * nstar - a * n) : ℕ) : ℚ) := by
        exact_mod_cast hkey
      push_cast at h2
      rw [h1] at h2
      linarith
    obtain ⟨-, hhi⟩ := roundedOK_bounds hP
    rw [failsRule_iff_rat]
    -- n* (2P + 1) ≤ 200000 a + 2 n*
    have hPn : (nstar : ℚ) * (2 * P + 1) ≤ 200000 * a + 2 * nstar := by
      have : (nstar : ℚ) * P ≤ nstar * (100000 * ((a : ℚ) / nstar) + 1 / 2) :=
        mul_le_mul_of_nonneg_left hhi hns'.le
      rw [mul_add, ← mul_assoc, mul_comm (nstar : ℚ) 100000, mul_assoc, mul_div_cancel₀ _ hns'.ne'] at this
      linarith
    have hmul : (n : ℚ) * nstar * (2 * P + 1) < 200000 * k * nstar := by
      have := mul_le_mul_of_nonneg_left hPn hn'.le
      nlinarith
    have hmul' : (n : ℚ) * (2 * P + 1) * nstar < 200000 * k * nstar := by linarith
    exact lt_of_mul_lt_mul_right hmul' hns'.le

/-- (6.1 (b)) exact within the lcm bound. -/
theorem exact_lcm : exact_lcm_stmt := by
  intro nstar m n k P hns hn hP hl
  exact exact_core _ nstar n k P hns hn hP hl

/-- (6.1 (b)) every design with groups of at most 316. -/
theorem exact_316 : exact_316_stmt := by
  intro nstar m n k P hns hn hP hn316 hns316
  apply exact_core _ nstar n k P hns hn hP
  have h1 : Nat.lcm n nstar ≤ n * nstar := Nat.le_of_dvd (Nat.mul_pos hn hns) (Nat.lcm_dvd_mul n nstar)
  have h2 : n * nstar ≤ 316 * 316 := Nat.mul_le_mul hn316 hns316
  omega

/-- (6.1 (c)) balanced designs. -/
theorem balanced : balanced_stmt := by
  intro n m k P hn hlt hP
  have hcore := exact_core (n - min n m) n n k P hn hn hP (by rw [Nat.lcm_self]; exact hlt)
  rw [hcore]
  unfold failsIntended
  have hn' : (0 : ℚ) < n := by exact_mod_cast hn
  rw [gt_iff_lt, div_lt_div_iff_of_pos_right hn']
  exact_mod_cast Iff.rfl

/-- (6.1 (d)) the share grows with group size. -/
theorem share_mono : share_mono_stmt := by
  intro m a b ha hab
  have hb : 0 < b := lt_of_lt_of_le ha hab
  have ha' : (0 : ℚ) < a := by exact_mod_cast ha
  have hb' : (0 : ℚ) < b := by exact_mod_cast hb
  unfold allowedShare
  rw [div_le_div_iff₀ ha' hb']
  by_cases hma : a ≤ m
  · have : a - min a m = 0 := by omega
    rw [this]
    simp only [Nat.cast_zero, zero_mul]
    positivity
  · have hmin_a : min a m = m := by omega
    have hmin_b : min b m = m := by omega
    rw [hmin_a, hmin_b, Nat.cast_sub (by omega), Nat.cast_sub (by omega)]
    have : (m : ℚ) * a ≤ m * b := by
      apply mul_le_mul_of_nonneg_left _ (Nat.cast_nonneg m)
      exact_mod_cast hab
    nlinarith

/-- (6.1 (f)) the peptide path's floor form. -/
theorem peptide_floor : peptide_floor_stmt := by
  intro P n k
  unfold failsRule
  generalize (n : ℤ) * (2 * P + 1) = X
  omega

/-- (6.1 (g)) beyond the bound, with one admissible rounding. -/
theorem counterexample : counterexample_stmt := by
  unfold counterexample_stmt RoundedOK failsIntended failsRule allowedShare
  norm_num [abs_le]

/-- (6.1 (g)) beyond the bound, whatever the rounding. -/
theorem counterexample_any_rounding : counterexample_any_rounding_stmt := by
  intro P hP
  obtain ⟨hlo, hhi⟩ := roundedOK_bounds hP
  have hs : allowedShare 321 320 = 1 / 321 := by
    unfold allowedShare
    norm_num
  rw [hs] at hlo hhi
  have h1 : (311 : ℚ) < P := by linarith
  have h2 : (P : ℚ) < 313 := by linarith
  have h1' : (311 : ℤ) < P := by exact_mod_cast h1
  have h2' : P < (313 : ℤ) := by exact_mod_cast h2
  have hP312 : P = 312 := by omega
  subst hP312
  refine ⟨?_, ?_⟩
  · rw [hs]
    unfold failsIntended
    norm_num
  · unfold failsRule
    norm_num

/-! Witnesses -/

theorem never_stricter_witness : never_stricter_witness_stmt :=
  ⟨1 / 2, 50000, 2, 2, by norm_num, by norm_num, by norm_num, by unfold RoundedOK; norm_num,
    by unfold failsRule; norm_num⟩

theorem roundedOK_3_2 : RoundedOK 33333 (allowedShare 3 2) := by
  unfold RoundedOK allowedShare
  norm_num [abs_le]

theorem exact_lcm_witness : exact_lcm_witness_stmt :=
  ⟨3, 2, 3, 33333, by norm_num, by norm_num, by norm_num, roundedOK_3_2, by norm_num⟩

theorem exact_316_witness : exact_316_witness_stmt :=
  ⟨3, 2, 3, 33333, by norm_num, by norm_num, by norm_num, roundedOK_3_2, by norm_num, by norm_num⟩

theorem balanced_witness : balanced_witness_stmt :=
  ⟨3, 2, 33333, by norm_num, by norm_num, by norm_num, roundedOK_3_2⟩

theorem share_mono_witness : share_mono_witness_stmt :=
  ⟨1, 2, 3, by norm_num, by norm_num, by norm_num⟩

theorem counterexample_any_rounding_witness : counterexample_any_rounding_witness_stmt :=
  ⟨312, by unfold RoundedOK allowedShare; norm_num [abs_le]⟩

end MultiScholaRProofs.FilterCount
