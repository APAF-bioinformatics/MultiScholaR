# KIT-247 second correction, part 2 (coordinator's rulings on I07; the owner's VAL-288), declared before sending

Merged origin/v0.9-staging 12701142 (32c0ec57): NEWS and ledger appends kept both, staging's first; nothing else
conflicted. Tests only: 97fd133f (one commit for the I07 recording, I07b, E15b and VAL-288's rows).

## R1''' (tests only 97fd133f, on the merged tree with the I01/I02 fix 40948dd8 but WITHOUT the VAL-288 fix;
## test-kit247-design-resave.R; ceiling 110 job-minutes)

Predictions:
- I07 prot_lfq: passes. I07 lipid_canonical: fails as in R1'' (expectations unchanged); its KIT247-LOG lines record the
  notifications, the history after each press and whether the run committed.
- I07b (lipid, Quality Control re-applied): passes.
- E15b: passes (the run after the cut commits a state, or shows an error or warning).
- F288a prot_lfq and lipid_canonical: fail (a formula-only re-save is judged unchanged: no question).
- F288b: fails (an organism-only import is not asked about as a limited change: no question, or a full cut).
- F288c: errors (omicsDesignResaveLimitedText does not exist).
- Every other row passes (the E7 rows and I02 now pass: the I01/I02 fix is in).

## Classification of I07 on lipidomics (the coordinator's rule, written before the run)

- If the recorded run shows the application REFUSING, or not completing, the normalisation from the design state with
  a message that tells the user what to do: I07's sequence was not one a user can complete; I07 is recorded as a
  test-design error (its failure stays in the records as it happened) and I07b carries the assertion.
- If the run completes and gives a DA on the OLD groups, or commits nothing with NO message, or the record contradicts
  itself: a PRODUCT result (the original symptom, or a silent failure, on lipidomics): stop and report, no fix.
- If I07b fails: a product result too: stop and report.

## E15's gap

E15 asserts only that no QC state follows the design state, so a run that commits nothing passes it. E15b asserts that
the run after the cut either commits at least one state or shows an error or warning notification; never neither.

## The comparison (after the VAL-288 fix): base staging 12701142, head the commit carrying the fix and this record;
## the I12 scope (r2pp-files.txt, 76 files) + check; 4 + 4 shards; ceiling 900 job-minutes

Predictions: every file present on base same; the three kit247 files only on head and passing (all rows, I07 lipid
excepted if it is classified a test-design error: then it fails on head exactly as in R1'''); check identical.
