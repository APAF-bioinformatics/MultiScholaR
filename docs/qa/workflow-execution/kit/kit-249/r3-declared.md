# KIT-249 correction-pass run (R2 again): declared before it is sent

- **Purpose.** The ticket's one correction pass (stop rule 6), on the tree merged with staging 63143aa3 (merge train 5):
  the implementation audit's guard G1b (I03) and a re-check of R2's verdict on the current staging. It decides nothing
  about the stored captures (parked row 36, the owner's); it shows whether the moves are exactly I01's.
- **Commits.** Head `kit-249` 144d2c5d (the fix, G1b, merged with staging 63143aa3). Base `kit-249-r3-base` 968804de
  (144d2c5d with `R/func_prot_da_model_stats.R` as staging 63143aa3: the tests-only tree with G1b, without the fix).
- **Files.** The 77 of R2, re-taken at 144d2c5d: section 9's grep gives the same 70 (KIT-239's, KIT-214's, KIT-226's and
  KIT-241's new test files do not match it), + 7 named; heaviest first; 6 shards per side; R CMD check.
- **Ceiling.** 12 x 70 + 2 x 90 = 1,020 job-minutes (68.00 core-hours), the ticket's correction row.
- **Predictions.**
  1. `test-kit249-prot-da-qvalue-fallback.R` `fixed`: base 88 passed / 54 failed (R2's 82/54 plus G1b's 6, which pass
     on both sides), head 142 / 0.
  2. `test-workflow-saved-da.R` `message-changed`, base 171 -> head 192 failures, moving in exactly the leaves of
     AUDIT-KIT-249-impl-20261125 I01 (the 21 info lines of `r2-45601895895/saved-da-new.txt` newly failing; the 87
     changed messages in the leaves of I01's table; `limpa_one_contrast` unchanged).
  3. Every other file `same`.
  4. R CMD check identical on both sides (nothing only on head).
- **What follows.** All four hold: R2's verdict stands on the current staging, with G1b; the ticket waits only for the
  owner's word on row 36. Any other move: stop and report (stop rule 3); no further run.
