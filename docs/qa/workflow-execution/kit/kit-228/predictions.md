# KIT-228 predictions (committed before the run they concern)

## G0 (inventory, ticket section 5): done by reading, no run requested

`omicSummaryW1Frozen` (`tests/testthat/test-omic-summary-module.R`): every case has one Save step, or (`trx_main`) the
Save before the RUVIII export (the gate, on a fresh folder) followed by the export and a successful Save; the `*_no_label`
steps never reach the Save; `omicSummaryW1Start` leaves `integration/` empty. No case holds `integration/<omic>/` at the
start of a Save that ends without a new input. GUI-entry lanes (`tests/testdata/gui-entry/*.json`): `summary_refusal`
precedes `summary_save` in all ten and lists no `integration/` artifact (ticket section 6, re-read at this head: the
lane files were not edited since). Prediction: no listed case; stop rule 2 does not fire. R2 (base against head on the
same files) is the check that nothing moved; the G0 run (2.0 core-hours expected) is not requested: a by-reading inventory
the comparison run then tests is cheaper and stronger than an instrumented base run. Recorded as a departure.

## R1 (tests-only commit 96f32e91, R/ equal to staging 0bdf7213; one file, one job; ceiling 60 min = 4.0 core-hours)

Purpose: show each new failing-first row failing on the old code for its stated reason (stop rule 1). Decides: whether
the fix may be written. The E rows fail at the folder's presence and the missing notice; G rows and the class-lost checks
pass; counts of passing expectations are not predicted. Predicted failing expectations (22) by test:

| Test | Failing on the old code | Passing |
|---|---|---|
| R01 | 1: the message names "kept" with no previous input | the "could not be published at" match, the second half |
| E4, G1, G3, G4 (proteomics) | E4: folder present (1), notice absent (1); G3: no warning naming the folder (1), no T4 (1) = 4 | first Saves, G1, fallback file, G4 |
| E1 (trx) | folder present (1), T1 absent (1) = 2 | preconditions, `failed` class, G1 |
| E2, G4 (trx gate) | folder present (1), T2 absent (1) = 2 | gate warning, "saved successfully", G4 |
| E5, E6 (trx managed) | E5: present (1), notice absent (1); E6: present (1) = 3 | class lost, `done`/`failed` counts |
| E3, E5, E6, G1, G2, G4 (lipid) | E3: 2; E5: 2; E6: 1 = 5 | class-lost, G1, G2 text, G4 |
| the same (metab) | 5 | as lipid |

G3 is skipped only when the runner user can write into a read-only folder (not predicted on a runner). If an instrument
defect shows (a precondition fails, a fixture helper errors before a result), it is repaired under the instrument-defect
rule and sent again; an E row passing on the old code stops the work (stop rule 1).

## Comparison run R2 (after the fix; base = staging 0bdf7213, head = the fix commit)

Files: the new file; the 27-file grep set of ticket section 6 (R2) at this head plus the 14 module-level Summary files by
name; KIT-212's, KIT-216's, KIT-219's (when it exists on staging; its branch is not mine), KIT-184's new test files;
`test-kit240-integration-default.R` (the merge check after the retired identity script); `check: true`. Predicted:
every file `same` except `only-on-head` for the new file (all rows pass on head, G3 included); no stored reference moves;
no `new-failure`; `check` identical. Observation to record: base-side counts that already differ from the last record
(checkpoint findings of train 35) are "unchanged" when identical on both sides.

## R2b (instrument repair of R2; declared before sending)

R2's only red item: `test-kit247-design-resave.R` hit the file timeout I set (3600 s) on base and on head, so it gave no result on either side (an instrument defect of my request: no result seen, no expectation touched). It matters to KIT-228: its E23 asserts a refused Save removes the folder AND the top-level file, and E21 the design cut. Request R2b: that one file, head only (a7e3a582), `file_timeouts` 150 min, ceiling 160 job-min = 10.7 core-h; running total before it 37.0 of 92.0. Prediction: all expectations pass (the start removal keeps the folder-and-top-level-file removal); any failure there is a result, to be judged against base by a second request only if needed.
R2b attempt 1 (45794764741): infrastructure failure at library setup (apt 404 for libtiff, 16 min, no test started); rerun of the failed job only, instrument rerun 1 of 3, no change to request, ceiling unchanged (running total 37.0 + 0.3 core-h before it)
