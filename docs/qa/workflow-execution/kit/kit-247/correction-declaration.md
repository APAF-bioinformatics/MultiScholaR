# KIT-247 correction pass (the one allowed; coordinator's ruling after R2 45612382879), declared before sending

Merged origin/v0.9-staging 42846c0e into kit-247 (b88752e8): NEWS and ledger appends kept both, staging's first (the
ledger's duplicated VAL-230 to VAL-233 context rows removed); no other conflict.

## R1' (tests only, e6871939; ceiling 30 job-minutes)

Files: test-kit247-resave-no-session.R, test-workflow-native-session.R. Predictions:
- N1 FAILS (compareOmicsDesignResult returns "changed" when the bus holds no earlier design).
- N2 FAILS (no classed refusal: the gate errors on session$ns, class not omics_design_resave_unconfirmed).
- N3 PASSES (with a session the gate asks; unchanged behaviour).
- test-workflow-native-session.R: the KIT-026 test errors as in R2 ("attempt to apply non-function"); every other
  test of that file passes.

## R2' (base origin/v0.9-staging 42846c0e, the merged staging, so the comparison holds only this ticket's change; head the commit carrying this file; 4 + 4 shards, check; ceiling 740 job-minutes)

Files: r2p-files.txt (69 = R2's 66, plus test-kit247-resave-no-session.R and the two files staging gained that match
the grep: test-kit224-ruv-feature-levels.R, test-kit251-replicate-filter-one-click.R). Predictions:
- test-kit247-design-resave.R, test-kit247-import-config.R, test-kit247-resave-no-session.R: only on head, all pass.
- test-workflow-native-session.R: same (KIT-026 passes again).
- Every other file same. (The R1 base 526ec2aa sat on e0747fc1; R2' compares against the merged staging instead, so a
  difference cannot come from staging's trains 5 to 7.)
- Check: nothing only on head.
