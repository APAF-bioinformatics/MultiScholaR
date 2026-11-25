# KIT-254 correction pass R1d: predictions, written before the request

Commit: staging c4505c7f plus `tests/testthat/test-kit254-omics-correlation-skip.R` from 44e3e910 only (branch
`kit-254-r1d`, never merged): no fix, no statement code. The new file, one job. Purpose: E7 failing first at the
constants of stop rule 2 (AUDIT-KIT-254-impl I01 step 2), and E17 failing first after its instrument fix.

- Every E row fails at least one assertion at base's value, as in R1 (45604366935) for E1 to E14, now on staging.
- E7 (lipid with RUV, t_hi 0.5, nothing removed): fails on its state and record assertions (E1's history and
  components, E3's threshold line, E4's record, E5's source); its design assertions pass on both (stop rule 2).
- E8, E15, E16, E17, E18: fail at "no statement"; E17's and E18's preconditions pass.
- No G row (G1 to G13) and no precondition fails (I01's precondition passes: the constants now say nothing removed).
