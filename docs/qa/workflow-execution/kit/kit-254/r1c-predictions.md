# KIT-254 correction pass R1c: predictions, written before the request

Commit 993c454b (kit-254 with staging c4505c7f merged, the fix, no replay statement code; the correction pass's tests).
The new file only, one job.

- **I01 precondition** (`Apply(t_hi) precondition`): lipid and metab pass (RUV skipped; KIT-224 changes the legacy RUV
  correction only, and these fixtures run no duplicate resolution, so KIT-214 does not apply). Lipid with RUV: **not
  predicted** (KIT-224 changes the RUV-corrected table): if it fails, I01's step 2 follows (re-derive t_hi, t_lo by G0's
  rule on the merged tree, a declared test-only edit of the two constants, E7 failing first again).
- E1 to E14 pass (the fix is present); G1 to G13 pass.
- Fail, at "no statement" (the statement code is absent): E8 (lipid, metab: the statement's type, its sample text, its
  threshold, its warning); E15 and E16 (lipid, metab: the second statement's text, its warning); E17 (lipid, metab);
  E18 (the correlation statement, its warning). E18's precondition (KIT-222's statement through its consumers' check)
  passes: KIT-222 is on this tree.
