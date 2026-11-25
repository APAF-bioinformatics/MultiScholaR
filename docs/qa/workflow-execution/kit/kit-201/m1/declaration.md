# KIT-201 M1 declaration (committed alone, before M1 runs)

Fixed by ticket r1.7 section 5.2 (D8, VAL-127 (3)); nothing below is chosen by hand.

- **Script:** `docs/qa/workflow-execution/kit/kit-201/m1-sim.R`, SHA-256
  `19bf9b14ae8c2c7141438f0f2a8125fac5abaa071ecb58848d17432ea8abbd55` (as committed at ed56c7a6; unchanged since).
- **Family:** `m1/family.tsv`, SHA-256 `834e29df811fdcc660a05b9224611bb619ec99aebf32c1a88fe73a6b6b5d61a0`: K = 990
  gated comparisons (3 kernels x 18 cells; FDP, and the null rate overall, per third of the true level and, in VD20
  cells, among features with a missing value).
- **Precision:** `m1/precision.tsv`, SHA-256 `81419fe9362d959de5859866b48c74dc8283a9b00f2f0a103743724a62c3b185`, the
  formula's output from the blind R0 (`m1/r0-sd.tsv`, SHA-256
  `8ce7f5d3bd1dd3d879297004e8e84a2f5c3622fd653a36c59657c5e53be70afb`): SE* = 0.001712; R per cell
  max(200, ceil((1.5 SD_max / SE*)^2)), at most 1,000; projected M1 cost 11.63 core-hours (cap 180.0, not binding),
  3 jobs; 982 comparisons gated, 8 reported, not gated (projected SE above SE* at R = 1,000).
- **R0:** run 45573413958, base dd769671, head 65da7eab: every lifecycle check passes (`m1/r0-checks.tsv`).
- **Seeds:** data set r of cell c: 201000000 + 10000 c + r; M1 main r = 1..R, calibration r = 5001..5000 + R;
  R0 used r = 9001..9005 (no data set shared). Bootstrap: 2,000 resamples, seed 201999.
- **Sides:** base = staging dd769671 (KIT-196, KIT-200 and KIT-177 merged); head = this commit, whose `R/` tree is
  65da7eab's (the fix, merged with dd769671).
- **Rule:** margin 0.01; FDP: pass when the upper 95 % limit of head - base is at most 0.01, fails materially when its
  point estimate exceeds 0.01; null rate: pass when the upper limit of head - base or of head - 0.05 is at most 0.01,
  fails materially when both point estimates exceed 0.01; otherwise imprecise. Per kernel: no worse (every gated
  comparison passes), worse (one fails materially), inconclusive (none fails, one imprecise). Worse or inconclusive
  stops that kernel for the owner; M1 is never re-run after a classification.
- **Untested rows:** a feature without a p-value is not called, counted per data set and left out of every denominator.
