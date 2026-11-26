# KIT-262 outcome through M2 (M3 not started)

Scope of this pass: step 0, M1, M2 and the skip rule's verdict. Branch `kit-262` from `origin/v0.9-staging` 15d86bd4. No package code, NEWS or version change; M3 is not started (needs KIT-242 Part B on staging and the coordinator's go). Files: `step0.md`, `m1-declaration.md`, `m1/`, `lifecycle.md`, `declaration.md`, `m2-continuation.md`, `m2/` (counted summaries; per-percentage rows gzipped), `results/` (scorer output), `runs.tsv`.

## Observed

**M1** (run 45712775309, b35b17d0; unhooked sweep read; hooked sweep identical, E1' pass). Verdict **confirmed**, lipidomics and metabolomics identical. Assay two: 2 %: k = 2, c = 0; 3 %: k = 2, c = 1; 4 %: k = 2, c = 1; 5 %: k = 3, c = 1; 6 %: k = 5, c = 2; every percentage R5, no usable percentage under P or C1; C2 (curve cut at the cap 3) is usable at 6 % only (k_b = 2, c = 2). Assay one: ok at 2 to 4 % (k = 2). The trace equals KIT-190's O1 print line for line (22 of 22). Predictions: R5 everywhere with c in 0 to 9, identical omics, assay one has an ok row: held. The ticket's C1 sentence (at least one percentage usable) was contradicted, as the declaration anticipated from the O1 print. E1, E2, E3, E7 pass.

**M2** (L2: runs 45714964374 and 45721854918, all checks pass after one instrument rerun; counted: runs 45727478832, 45727493119, 45727503364, 45727516470 and the continuation 45750756599; scoring 45758847377). All 198 cells x 40 data sets present, 0 error rows; E1, E2, E3, E4, E7 pass in all 7,920 data sets, E1' in the declared sample (r = 7, 27), E8 20 of 20. Two of the four jobs hit their 190-minute timeout (planned seconds came from the lifecycle, counted sweeps ran 1.3 to 1.6 times slower), no cell partial; the 17 missing cells were run once in the continuation (`m2-continuation.md`). Tables: `results/m2-cells.tsv`, `m2-causes.tsv`, `m2-gain.tsv`, `m2-resid-df.tsv`, `m2-stress.tsv`, `m2-verdict.tsv`; pooled over g and F in `results/m2-pooled.md`.

- Gain rule (z = 3.452, lower limit >= 0.05; 7 of 40): **a (C1 recovers RUV) shown in 17 of 180 cells; b (C2 against P) in 19; e (C2 over C1) in 1.** Largest C1-not-P: 36 of 40 (cell 114: n = 48, 2 groups, F = 1,000, MCAR 5 %). The cells with a: MCAR 20 % at n = 12 (F = 300); MCAR 5 to 10 % at n = 24; MCAR 5 % at n = 48; MCAR 1 % at n = 96 (F = 300); AD 20 % at n = 48 (2 groups) and n = 96 (both group counts), AD 10 % at n = 96 (2 groups, F = 300). No cell with n <= 6 and none with 24 samples under abundance-dependent missingness.
- **Skip rule (O-3 (2)): not applied; M3 is not skipped** (a and b both true). The stress cell by the declared selection: **cell 79** (n = 24, 2 groups, F = 300, MCAR 10 %; C1-not-P 34 of 40); runner-up cell 77 (MCAR 5 %, 8 of 40); the AD cells and MCAR 20 % have 0.
- Complete-controls cause (the owner's question): see the pooled table. Under abundance-dependent missingness a P-failure is rare below n = 96 (n = 96: 30 % at AD 20 %, 7 % at AD 10 %) and in every such data set C1 recovers RUV (cause "rescued by C1", no completeness-only failure). Under MCAR, failures concentrate where the complete-feature share collapses: P without RUV in 35 % (n = 12, 20 %), 66 % (n = 24, 10 %), 100 % (n = 24, 20 %), 84 % (n = 48, 5 %), 100 % (n = 48, 10 and 20 %), 41 % (n = 96, 1 %), 100 % (n = 96, 5 to 20 %). There C1 recovers most of the lost data sets at n = 12 to 48 up to MCAR 5 to 10 % (no RUV under C1: 6 %, 8 %, 0 %, 8 %) and at n = 96 MCAR 1 % (0 %), but not at the cells above (C1 still without RUV in 96 % at n = 24 / MCAR 20 %, 98 to 100 % at n = 48 / MCAR 10 to 20 % and n = 96 / MCAR 5 to 20 %); those remaining failures are all "completeness only" (c below min(k, cap) at every percentage). Row level: R5 rows are mostly cap-side (c at or above min(k, cap)) except in the heavy-MCAR cells, where completeness dominates (for example n = 24, MCAR 20 %: 3,065 of 3,120 R5 rows).
- Residual degrees of freedom (`results/m2-resid-df.tsv`): k used = cap in most data sets with RUV (P: for example 605 of 719 at n = 6, 2 groups); n - g - k used = 0 for none at n >= 6 and 2 groups; at n = 6 with 3 groups, 1 residual df in every data set with RUV and a k used of 2; at 2 x 2 (n = 4, O-4) k = 2 leaves 0 of 2 in 410 of 719 data sets with RUV (all three rules).
- Changed choice where P already runs RUV: C1's percentage differs from P's in a large minority of those data sets (for example 135 of 159 at n = 48, MCAR 1 %; 88 of 150 at n = 24, MCAR 5 %; about 25 to 65 % under abundance-dependent missingness at n = 24 to 96).

## Inferred (not tested)

M2 describes the search's behaviour on simulated lipid data with one unwanted factor, no between-sample normalisation, MCAR and one abundance-dependent form; it says nothing about error rates (M3) or real data. Whether C1's changed controls are harmful is exactly what M3 would test.

## Under-specified items and how they were resolved (declaration.md)

E4 for abundance-dependent missingness (expected share at the solved lod, realised share recorded); error rows and incomplete cells; complete controls from R5 reasons only; the skip rule read as "neither a nor b" (b = C2 against P, so a C2-only gain would also prevent a skip; it did not decide here, a holds).

## Not checked

M3 and everything it contains; E9, E10, E11, E12, L3, the re-estimation of R; the diff over the sweep, RUV-step and DA files between the bindings and the staging head (staging moved to 176f9d3f, ledger only, at branching); real-data missingness; metabolomics beyond E8 (identical in 50 data sets and M1); the ledger row (not written: it needs the outcome).

## Cost (core-hours = job-minutes x 4 / 60)

M1 0.88; L2 5.56 (2.99 + 2.57); M2 51.33 (13.43 + 8.03 + 13.45 + 11.86 + 4.56 for the continuation); scoring 0.80. Total 58.57 against the ticket's expected 88.81 for everything (42.99 if M3 is skipped) and ceiling 221.63; M2 alone above its expected 39.33 and inside its ceiling 58.60.
