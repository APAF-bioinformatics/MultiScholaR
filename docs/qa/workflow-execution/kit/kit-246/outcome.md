# KIT-246: outcome

Ticket `tickets/kit/KIT-246.md` r1.2 (step 0: AUDIT-KIT-246-plan-20261125 Re-check 1, R1-01 to R1-05, applied
with the register's texts). Owner's wording: VAL-251 (2026-11-25, the Model tab's help) and VAL-293 (2026-11-26,
the Methods sentence: option (a), the help's audited passage in full, without the effect size).

**State: fix done, R2 as expected; ready for the implementation audit.**

## Branches and commits

- `kit-246` (from 3eeaf100): `2d81e890` tests only; `13e81329` R1 records. Kept as R1's evidence.
- `kit-246b` (from `origin/v0.9-staging` 127011429, after KIT-239 and KIT-240 merged; the branch to merge):
  - `a86922c0` tests only (2d81e890 carried over; the conflict with KIT-239's appended blocks resolved by keeping
    both, no existing block edited);
  - `1ee26efa` R1 records carried over;
  - `c74f2a95` the fix, E2 and E5, NEWS, ledger INT-31 and VAL-327;
  - then the R2 records (this file, `wording.md`, the tables).

## The change

- `R/utils_integration_gui.R:34-45` (new, `@noRd`): `integrationScalingEvidenceText()`, the evidence passage
  written once.
- `R/utils_integration_gui.R:291-296`: `help_scaling` = its first sentence + `integrationScalingEvidenceText()`.
  Its text is unchanged (string comparison against 127011429 before the commit).
- `R/utils_integration_record.R:164-167`: `methods_scaling_feature_zscore` = "each scaled to mean 0 and standard
  deviation 1 after the selection (which used the unscaled values)." + `integrationScalingEvidenceText()`.
- So the Methods sentence and the help cannot drift: both end with the one passage, and E5 asserts it. Old and new
  sentences in full: `wording.md`.
- Not changed: `integrationRecordAnalysis()`, `integrationMethodsText()`, `methods_scaling_none`, the record's
  fields and schema, any computed value. `man/`, `NAMESPACE`, `DESCRIPTION`, `inst/` untouched.

## Runs (`runs.tsv`)

| Step | Package / base | Run | Result |
|---|---|---|---|
| R1 (refused) | 2d81e890 / none | 45584833302 | request committed blank; plan refused it in 7 s; no test run, no result |
| R1 | 2d81e890 / none | 45584885011 | as expected: 139 passed, 6 failed (E1 only) |
| R2 | c74f2a95 / a86922c0 | 45675980031 | as expected: 10 files `same`, the record file `fixed`; check identical |

Declared before R2:
- Purpose: the fix against the tests-only commit. Binds c74f2a95 and a86922c0.
- 11 files, one shard per side, `check: true`.
- Ceiling 280 job-minutes (2 x 80 test, 2 x 60 check).
- Decides A2 and stop rules 3 and 4.

Scope, recounted at the branch point with `grep -lE "integrationMethodsText|integrationStudyParameterLines|
integrationRecordTexts|methods_scaling|help_scaling" tests/testthat/test-*.R`:
- The ticket's eight files (`test-integration-workbook-api.R` is still present).
- KIT-239's reader of `help_scaling` (`test-integration-results.R`).
- KIT-240's reader of the paragraph (`test-kit240-integration-default.R`).
- `test-integration-gui.R`, because `R/utils_integration_gui.R` joined the files (section 3 item 2 (ii)).
- Left out as in section 10: `-enrichment-tests.R`, `-enrichment-options.R`, `-outputs-mofa.R`, `-record-mofa.R`
  (they read named keys of other sentences only).
- No MOFA lane (section 9).

## Observed

R1 (run 45584885011, at 3eeaf100 + tests; `r1-*`):
- `test-integration-record.R` 139 passed, 6 failed, 0 errors.
- All six failures are E1's three assertions (lines 434-436: "weak factors", "no practical difference", the old
  sentence), each for the paragraph as built and as read back from JSON, actual TRUE.
- E1's "contains the key text" check and G1 passed; no other block failed.

R2 (run 45675980031; `r2-compare.tsv`, `r2-tests.tsv`, `r2-tests-base.tsv`, `r2-failures-base.txt`,
`r2-check-compare.txt`, `r2-meta.txt`):
- Base (tests-only a86922c0, on staging 127011429): `test-integration-record.R` 173 passed, 6 failed: E1 at lines
  550-552, twice each. E1 fails first on the new base too.
- Head (c74f2a95): `test-integration-record.R` 200 passed, 0 failed, verdict `fixed`. The KIT-246 block went from
  6 to 33 passed, and no other block's count changed (173 - 6 = 200 - 33 = 167).
- The other 10 files: `same` (identical counts on both sides, no failures). Skips, identical on both sides:
  removed-factors 3, enrichment-mofa 3, kit240-integration-default 1.
- R CMD check: 0 ERROR, 3 WARNING, 3 NOTE on both; 6 findings identical; none only on head, none only on base.

## Inferred

- **E1:** every record built with the fix prints the owner's sentence, as built and from JSON.
- **G1:** a record with `scaling = "none"` keeps its own clause.
- **E2:** the key holds the decided numbers and terms, and none of "weak", "no practical difference", "too small",
  "there was none" or "rose by".
- **E5:** the help and the Methods sentence end with the same passage and share every term E5 lists.
- The count change in `test-integration-record.R` is E1's six now-passing assertions plus E2's and E5's 21
  expectations: 6 + 21 = 27 more passes, so the block went from 6 to 33. The ticket's A2 expected the verdict
  "count-changed". The runner reports "fixed" because the base side carries E1's six failing-first failures. That is
  the declared outcome under the runner's name for it.
- No stored reference moved (stop rule 4 not reached).

## Acceptance criteria

| # | State | Evidence |
|---|---|---|
| A1 | met | R1 45584885011: E1 fails by assertion at its base value; G1 passes; nothing else fails. Shown again on the staging base in R2 (base side) |
| A2 | met | R2 45675980031: the new block passes (E1, G1, E2, E5); 10 files `same`; `test-integration-record.R` changed only in the KIT-246 block; check findings only on head: 0 |
| A3 | open | the next checkpoint's full suite (shared) |
| A4 | met (numbering at merge) | NEWS (`### Fixes`); ledger INT-31 (after INT-30) and VAL-327 (after VAL-294) citing VAL-251 and VAL-293; this file with "For the change report"; `wording.md` |

## Budget

| Run | Job-minutes | Core-hours |
|---|---|---|
| R1 refused | 0.1 | 0.01 |
| R1 | 24.0 | 1.60 |
| R2 | 115.4 | 7.69 |
| **Total** | **139.5** | **9.30** |

Expected 10.53, ceiling 53.33 core-hours. No correction, instrument rerun or wording run was used.

## Findings (ticket section 12; not changed here)

- F1: KIT-151 d's `z-results.md` proposes the old wording; INT-31 says not to use it.
- F2: `git grep -n -i "weak" c74f2a95 -- R/utils_integration_*.R`: no "weak" remains in the scaling texts. The
  remaining hits are in the enrichment texts ("factors weaker than the main effect", "Weaker factors", and internal
  cell names such as `F_weak`; `R/utils_integration_gui.R:421`, `:772-1063`; `R/utils_integration_record.R:468`,
  `:472`; `R/utils_integration_enrichment.R:928`, `:1061`, `:1063`). They describe the enrichment check's factors
  relative to the main effect, from other evidence; left for their owners.
- F3, F4: as the ticket states (selection and z-score in `k151d.R` tied to the defaults by settings; "variance
  explained" on z-scored features).
- The comment above KIT-239's L7 in `tests/testthat/test-integration-results.R` still says "draft pending the owner";
  the owner has decided (VAL-251). Left as it is (not this ticket's file).

## What a user with an existing project sees

- Reports, `study_parameters.txt` files and workbooks already saved are files on disk and do not change.
- When the project's record is saved again with this version (Record tab Save, or `integrationWriteRecord()`), the
  new report, `study_parameters.txt` and Record tab carry the corrected sentence.
- A replay prints no paragraph; any paragraph built afterwards carries the corrected sentence.
- The saved record names the package version and commit that wrote it.

## Not checked

- No R was run locally; everything above is from the runners' tables.
- The rendered report's text with the new sentence (the MOFA lane renders it; not needed, section 9: no assertion
  needs a real model).
- Per-block counts outside the KIT-246 block (taken from the file totals and the runner's per-test note).
- The full suite (A3, at the next checkpoint).
- "About a fifth" is the declared 19 % (the realised share is not in KIT-151 d's committed files).

## For the change report

- **What was wrong, who was affected.** Every integration record made with the default per-feature scaling (INT-24)
  printed this in its draft Methods (report, `study_parameters.txt`, Record tab): "in simulation, with few samples
  and missing values this scaling improved the recovery of weak factors, and otherwise made no practical difference."
  - "Weak factors" was never defined.
  - "Otherwise made no practical difference" was not shown. At 20 samples without missing values the result was
    inconclusive: 29 against 26 of 30 recovered, 3 data sets only with scaling and none only without, exact McNemar
    p = 0.25. At 60 samples both arms recovered 30 of 30, which is a ceiling.
  - Anyone copying the draft into a paper would have stated more than the evidence. No computed number was affected.
- **Legacy or new:** new to 0.9 (the integration and its Methods paragraph do not exist on `main`; KIT-142, KIT-151).
- **Old and new sentence:** `wording.md`.
- **Code at base 127011429:**
  - `R/utils_integration_record.R:164-167`: the entry, printed through `integrationMethodsText()` at `:1035`.
  - `R/utils_integration_gui.R:281-287`: `help_scaling`, already corrected by KIT-239.
- **Code at fix c74f2a95:**
  - `R/utils_integration_gui.R:34-45`: `integrationScalingEvidenceText()`.
  - `R/utils_integration_gui.R:291-296`: `help_scaling`.
  - `R/utils_integration_record.R:164-167`: the entry.
- **What changed in results:** nothing computed; one printed sentence.
- **Evidence:** KIT-151 d (run 44786925438, commit 25784597 on `run-mofa-evidence-k151d`):
  - 30 paired simulated data sets per cell; four data sets of 20,000, 2,000, 200 and 100 features.
  - The pattern F2 is on a tenth of the features of the two smallest, loading ±1.5 against noise sd 1.5 at the
    trend's centre.
  - Missing values (when on) depend on abundance, about 19 % as declared.
  - N = 20 and 60; recovery means |r| ≥ 0.9; one generator and one training seed.
  - Limits: synthetic data only; nothing between 20 and 60 samples.
  - Ledger: INT-31 (correction of INT-24's evidence), VAL-251, VAL-293, VAL-327. Runs: R1 45584885011, R2
    45675980031.

## Follow-up (VAL-305), 2026-11-26: predictions (written before the change and before any run)

Bounded item: put the four phrases of the text the owner was shown ("made for this software", "(two groups of 10)",
"this scaling", "of four data sets") into the one shared passage `integrationScalingEvidenceText()`; both the Model
tab's help and the Methods sentence change with it. Branch kit-246c from staging 176f9d3ff. The passage's last
sentence (VAL-251) is unchanged.

Tests-only commit bfc57ce7: one added block of assertions in the KIT-246 block (E5) of
`tests/testthat/test-integration-record.R`: 4 phrases x 2 expectations (the shared passage holds the phrase; the help
and the Methods sentence both hold it) + 1 (the passage's last sentence is unchanged) = 9 expectations.

Predictions for the run (base = tests-only commit, head = the change):
- `test-integration-record.R`: base 201 passed, 8 failed (the eight phrase expectations; the last-sentence one
  passes); head 209 passed, 0 failed; verdict `fixed`. Staging has the file at 200 passed, 0 failed before the 9 added.
- The other ten files: `same` (reading of the term lists of E2, E5, L7 and the whole-table scans: none of
  "weak", "every computer|any computer|on every machine", "same total variance|unit variance", "too small",
  "there was none", "no practical difference", "rose by" occurs in the four phrases).
- R CMD check: findings identical on base and head.

### Predictions restated for the run's scope (coordinator, option 1, 2026-11-26; new commit, before the run)

The run is the ticket's conditional wording line: `test-integration-record.R`, `test-integration-results.R`,
`test-integration-gui.R`, `test-integration-model.R`; one shard per side; no R CMD check in this run; base bfc57ce7,
head 10df1d31; env MULTISCHOLAR_TRX_MANIFEST_NOTE=1; ceiling 60 job-minutes (4.00 core-hours).
- `test-integration-record.R`: base 201 passed / 8 failed, head 209 passed / 0 failed; verdict `fixed`.
- `test-integration-results.R`, `test-integration-gui.R`, `test-integration-model.R`: `same`.
- No check in this run. The full 11-file scope with R CMD check rides the merge train (coordinator's ruling,
  2026-11-26). The first predictions above that name ten `same` files and an identical check are superseded by this
  scope.

## Follow-up (VAL-305), 2026-11-26: result of the wording run

Run 45712506511 (request kit246c-wording-10df1d31, request commit d2599e58; package 10df1d31 against tests-only
bfc57ce7; `test-integration-record.R`, `-results.R`, `-gui.R`, `-model.R`; one shard per side; no check; env
MULTISCHOLAR_TRX_MANIFEST_NOTE=1). Tables: `r3-compare.tsv`, `r3-tests.tsv`, `r3-tests-base.tsv`,
`r3-compare-failures.txt`, `r3-meta.txt`.

Observed, against the predictions:

| File | Predicted | Observed |
|---|---|---|
| test-integration-record.R | base 201 / 8, head 209 / 0, fixed | base 201 passed, 8 failed; head 209 passed, 0 failed; `fixed` |
| test-integration-results.R | same | 258 / 0 on both; `same` |
| test-integration-gui.R | same | 405 / 0 on both; `same` |
| test-integration-model.R | same | 265 / 0 on both; `same` |

- The eight base failures are the phrase expectations: each of the four phrases, for the shared passage (4) and for
  the help and the Methods sentence together (4); each actual FALSE (`r3-compare-failures.txt`). The expectation that
  the passage's last sentence (VAL-251) is unchanged passed on both sides.
- Cost: 54.6 job-minutes (head 22.6, base 31.1, plan 0.4, results 0.5; GitHub job timings) = 3.64 core-hours, under
  the 4.00 ceiling (60 job-minutes). One run, no correction used.
- All predictions met. The 11-file scope with R CMD check was not run here; it rides the merge train (coordinator's
  ruling, 2026-11-26).

Inferred: no other reader of the passage moves, from the three unchanged files and the search of step 0.

Not checked: R CMD check on this change; the other seven files of KIT-246's scope; the rendered Word report.
