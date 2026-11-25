# KIT-247 second correction (owner's VAL-277; AUDIT-KIT-247-impl-20261126 I01, I02, I07, G3), declared before sending

This is the last correction for KIT-247. Merged origin/v0.9-staging abfc99ac (0f9196a5): NEWS and ledger appends kept
both, staging's first (duplicate context rows removed); nothing else conflicted. Tests only: fbf729ef. Fix: 40948dd8.
Not in scope (not decided by the owner, parked row 47 (2)): whether a formula-only or an organism/FASTA-only change
counts as a change; the confirmation's wording. Behaviour and texts left as they are (wording.md).

## R1'' (tests only fbf729ef on the merged tree, without the fix; `test-kit247-design-resave.R`; ceiling 100 job-minutes)

Predictions for the new rows:
- E7 prot_lfq (confirm, cancel): PASSES (the gate and departure 3 already hold on the import; I01 is about errors).
- E7 error prot_lfq: FAILS (the refusal after the confirmation escapes the observer: the lane errors, no "Error during
  import" shown).
- E7 lipid_canonical (confirm, cancel) and E7 error lipid_canonical: PASS (the lipid import's commit is inside its
  handler; the refusal is before the bus).
- I02 (lipid): FAILS (the second save of the same design is judged "unchanged": "Design saved; nothing changed.", no
  cut).
- I07 prot_lfq and lipid_canonical: PASS (by construction since the first fix).
- G3: PASSES (guard).
- Every earlier row of the file: passes as at R2' (117).

## R2'' (base origin/v0.9-staging abfc99ac, head the commit carrying this file; 4 + 4 shards; check; ceiling 900 job-minutes)

Files: r2pp-files.txt (76 = R2''s 69 plus the audit's I12 scope: test-module-ci-trx-design.R, test-module-ci-trx-qc.R,
test-module-ci-lipid-norm.R, test-module-ci-metab-norm.R, test-kit243-summary-current-outputs.R,
test-kit244-limpa-pass-through.R, test-kit256-tech-rep-selection-cleared.R). Predictions: every file present on base
`same` (identical failure messages where it fails on base); test-kit247-design-resave.R, test-kit247-import-config.R and
test-kit247-resave-no-session.R only on head and passing; check findings identical, nothing only on head.

Budget: 90.8 of 148.67 core-hours used before this correction.
