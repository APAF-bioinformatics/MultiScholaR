# KIT-184: outcome

Ticket `tickets/kit/KIT-184.md` revision 2.0. Branch `kit-184-ruv-factors-in-integration`, from `v0.9-staging`
67d38df3, merged with staging at 8cb536cc for the full suite. All computation on GitHub-hosted runners; no R on the
laptop. Owner's decisions: methods ledger VAL-53, VAL-58, VAL-62 (the design), VAL-150 (the declaration), VAL-160 (the
results: the adjustment ships off by default with a stated limit).

**State (2026-11-25): ready for the auditor's re-check. The implementation audit's repairs are made, the control
held, `tests/testdata/gui-entry/prot_dia_limpa.json` is re-frozen, and the one further run the owner allowed has
returned as declared (VAL-178; results at the end).** The audit (`tickets/kit/audits/AUDIT-KIT-184-impl-20261125.md`) confirmed the method and
the shipped path and found three majors in the texts and the record (K184-I01 to I03), minors (I04 to I09) and one
stored reference that moves (I10). The sections from "What was built" to "Not checked" are the record as it stood at
the audit, with the audit's record repairs applied; "After the implementation audit" at the end has what was done
since, the control, the declared run and its results.

## What was built (observed in the code and the runs)

1. **The removed factors, written with the integration input.** Where an omic applied RUV to a data set, saving its
   Summary writes `removed_<data set>.parquet`: an orthonormal basis of the difference between the normalised and the
   corrected matrix over the features measured in every run, checked against the correction's k (`exact`,
   `approximate` with its share, or `not_determined` with a reason). The reader checks it; the import carries it in
   the inputs state. An input without a correction is written and read as before.
2. **The adjusted test.** For the tests that have a design (`camera`, `camera_floor`), the setting
   `adjust_ruv_factors` ("Adjust for removed RUV factors") puts those factors, averaged per biological sample as the
   values are, into the test's design (intercept, factor, removed factors), on the values before RUV. It is **off
   unless chosen, on every model** (VAL-160). Every result of such a test carries the label of its form; the adjusted
   form has its own key on the head (`...|camera|adjusted`); the replay repeats the saved form; an entry stored before
   this ticket replays without the adjustment. With the setting off nothing changes for any model: KIT-183's test,
   identity and stored results.
3. **The simulation** (declared, audited, run once: `declaration.md`, `results.md`, `s/`), and its numbers and limit in
   the package's texts (`wording.md`).

## Runs, in order (`runs.tsv` has every run with its minutes)

| Run | Commit | Observed |
|---|---|---|
| R1a 45321879760 | a5229552 (tests only) | every test of new behaviour fails, as declared |
| R1 45325536775 (with R0 and X) | 939a9d10 against a5229552 | the new file passes; check identical on both sides; base-store job passes; four mutations caught; gate-off identical; the writer's cost (`r1/`) |
| R1b 45332543296 (VAL-54) | 621a32b5 | the two K4 fixtures given missing values: both files pass; gate G2 (`r1b/`) |
| M 45332530206, B 45332726284, B24 45336831513 | 621a32b5, 83f745ee, 5e5df9b9 | K13 on real MOFA2; the new file on the installed bundle on three platforms |
| L1 45333818847 | 21053957 | lifecycle; stopped at the seconds-per-replicate rule; not scored blind (disclosed) |
| L1b 45348572248 | 9e525fd9 | blind lifecycle extension; every item of its list held (`l1b/`) |
| **S 45353248796** | 81ee3052 (the freeze commit) | the counted simulation: 2,000 of 2,000 replicates; outcome R1 with R7 (`results.md`) |
| R2 45388084875 (with X) | 74ce8305 against a613b9cb (tests only) | after VAL-160: the restated expectations fail on the tests-only commit (79 in 10 blocks) and pass on the change; check identical; four mutations caught (`r2/`) |
| M2 45398052678, B2 45398091058, B24b 45401782250 | 5ba85c1f and the lanes | K13 with the adjustment asked for by name; the changed file on the installed bundle on three platforms: 493 passed, 0 failed on each |
| **F 45402062487** | 0a511795 (branch merged with staging 8cb536cc) against staging | the one full suite with check: one stored reference moved (below); otherwise no file worse; check identical (`f/`) |
| G 45463933843 | 5a9cd699 against staging | gate-off pages identical at the final head |

Two requests were refused before any job ran and were sent again under VAL-54: the evidence lane's first push
(45333717225, an unquoted colon in the workflow's name) and R2's first request (45387864719, a requester line over
300 characters).

## Acceptance criteria

| # | Criterion | State | Evidence |
|---|---|---|---|
| A1 | Gates G0 to G5; G1 exact with rank equal to k on every route's fixture | met | `feasibility.md`; R1, R1b, L1 (G3); VAL-94 (G4) |
| A2 | K1, K2 | met | R1; F (the new file passes at the final head) |
| A3 | K3, with the base store | met | R1's and R2's base-store job (`r2/guards.tsv`, variant `head`) |
| A4 | K4 per route | met, two details not built (below) | R1, R1b; `feasibility.md` |
| A5 | K5 | met | R2, F |
| A6 | K6 | met | R2, F |
| A7 | K7, as restated in ticket 26.4: covariates only when the setting is switched on; runs without covariates keep their identity | met | R2; base-store job |
| A8 | K8: refusals with their classes and words; no test run in another's place | met (three sentences, the same on every model, after VAL-160) | R2 |
| A9 | The setting, as restated: off by default on every model; each choice gives the design, the label and the sentences; `TRUE` refused where no RUV was applied or no usable factors exist; both forms kept on the head; the replay repeats the saved form; an entry stored before this ticket replays without the adjustment; the test of corrected values without factors reachable by no setting | met | R2; `r2/mutation-verdict.txt` |
| A10 | K10, K11, K12: record, Methods, tab, workbook | met | R2 |
| A11 | K14, K15 | met | R2 |
| A12 | MOFA2 lane, K13 | met | M2 (`m2/`) |
| A13 | Simulation: L1's checks; declaration committed alone before the counted run; outcome by 6.5; identity check; determinism; texts filled | met but for determinism, which is **recorded as not met as declared** by the owner's decision (VAL-160 (3)): 7 of 12 tables identical as text, five differing only in real-number digits (at most 6.4e-13) between two processors. Outcome R1 with R7; identity check at most 6.2e-14; texts filled, no placeholder left | `declaration.md`, `results.md`, `s/`, `wording.md`; the results audit |
| A14 | Gate off | met, at R1's head and again at the final head | `r1/gate-off-verdict.txt`, `r2/gate-off-verdict-5a9cd699.txt` |
| A15 | Cross-platform record; the writer's memory and time | met: the new file passes on the installed bundle on ubuntu-22.04, windows-2025 and ubuntu-24.04 before and after VAL-160; the writer 9.6 s and 740 MiB at 20,000 x 500 (ceilings 60 s, 1 GiB) | `cross-platform.tsv`; `r1/removed-factors-cost.tsv` |
| A16 | Failing first; the guards fail on a broken head | met, twice (the first implementation; the VAL-160 change) | R1a; R2's base side; `r1/x-verdict.txt`, `r2/mutation-verdict.txt` |
| A17 | Full suite base against head with check: no new failure; check findings unchanged | **not met as written: one new failure, which is the ticket's own intended change meeting a stored reference** (`test-gui-entry-proteomics.R`, lane `prot_dia_limpa`, 4 expectations: the reference lists the integration input's files and lacks `removed_proteomics.parquet`). Everything else: 455 files same, this ticket's files pass, check identical | F (`f/summary.md`, `f/moved-references.txt`) |
| A18 | Ledger rows | met: three rows as `VAL-<merge>` (implementation; G1; the simulation and what was built from it, citing VAL-150 and VAL-160) | `docs/methods-decision-ledger.md` |

## Open for the coordinator before merge

1. **The moved reference** `tests/testdata/gui-entry/prot_dia_limpa.json` (one artifact entry in each of four steps).
   Ticket 8.5 foresees a reference of this kind moving and its re-freeze "on merge"; it names the kind, not the file
   (when the ticket was written no GUI-entry reference listed a Parquet file; KIT-199 re-froze this one on staging
   since). The implementer's brief lets me re-record a reference only where the ticket names the file, so I did not.
2. **F's cost**: 1,268.4 job-minutes against the row's ceiling of 1,080 (declared before sending as possibly up to
   1,200). By ticket section 11 it is not repeated or split, and its times go to the owner.
3. **Numbering the three ledger rows** at merge.

## Core-hours

| Row (ticket section 11) | Ceiling, job-minutes | Used |
|---|---|---|
| L1 (with L1b) | 75 | 41.4 |
| S (with the score job) | 1,415 | 1,043.0 |
| S' and re-scores | 275 | 0 |
| R0, R1 (with R1a, R1b), X | 525 | 308.0 |
| R2 (with its mutation script and the gate-off run) | 450 | 282.7 |
| M | 170 | 32.3 |
| B | 280 | 101.6 |
| F | 1,080 | **1,268.4** |
| Ticket | 4,550 (303.33 core-hours) | 3,077.4 (205.2 core-hours) |

The repository's own push check (Module CI) ran on every push of these branches, outside the ticket's budget.

## What changed after the owner's decision VAL-160, and what had assumed "on"

The first implementation made the adjustment the default for a model trained on RUV-corrected values (VAL-53, VAL-62).
After the simulation the owner decided it ships off by default. Changed, tests first (a613b9cb, then 74ce8305,
5ba85c1f):

- The default: `FALSE` in the parameter specification, the resolver, the tab's control and `integrationEnrichFactor()`.
  The head key: `|adjusted` for the adjusted form on any model; no `|unadjusted` (no released entry has that suffix:
  the first implementation never merged). The refusals for missing factors: three sentences, raised only when the
  setting is on. NEWS: one entry under New features; the Scientific changes entry withdrawn, since no default result
  changes.
- Tests that had assumed "on" and were restated: K5, K6, K8 (both blocks), K10, K11 (both blocks), K12, K14, K15,
  K7 (ii) with K16, K7 (iii) with K9 (a), the cross-platform block, and K13 in the MOFA2 lane's file. The
  cross-platform block was missed in the tests-only commit: R2 found it (4 failures: it compared the form without the
  adjustment with a reference of the adjusted form); fixed in 5ba85c1f by naming the form, no expected value changed.
- Scripts that had assumed "on": `base-store-job.R`'s mutations m2 (now: the default switched on) and m3 (now: the
  replay forgets the saved form; the old m3 would change nothing under a default that is off). Both are caught (R2).
  `simulate-k184.R` is the frozen evidence and is not changed: for metabolomics on a learned factor its arm 3 asked the
  product for its default and checked that the product then declared the adjusted form, so it would stop on a build
  whose default is off. The counted evidence does not depend on the default (`results.md`: every arm's form was fixed
  by the request or checked by the script).

## For the change report (plain words)

- **What the adjustment is.** When an omic's RUV correction has removed unwanted variation from a data set, the
  integration's enrichment test can take the removed unwanted factors into its model, so that a feature counts only
  for what it shares with the factor beyond them. It is an option the user switches on; without it every model is
  tested as before. Code: `integrationInputRemovedFactors()` and the writer and reader in
  `R/utils_workflow_integration_input.R`; `integrationResolveEnrichmentCovariates()`,
  `integrationEnrichmentRemovedFactors()`, `integrationEnrichmentTestValues()` and `integrationEnrichmentKey()` in
  `R/utils_integration_enrichment.R`; the control in `R/mod_integration_enrichment.R`; the texts in
  `R/utils_integration_gui.R` and `R/utils_integration_record.R`.
- **What was measured** (one simulated design, 2,000 replicates, 12 biological samples, RUV-III; run 45353248796).
  With the adjustment, sets with no link to what was tested were called in at most 5.4 % of cases (upper limit 6.4 %,
  against the limit of 7.5 %) in all 33 null classes of both deciding cells; every called true set was in the right
  direction; and more true sets were found than without it (sets of 10: 45 to 48 % against 40 %; 30 to 34 % against
  17 % where unwanted variation followed the groups). Testing the corrected values directly, which was never built,
  called measurably more.
- **The limit.** The adjustment does not protect sets of features that share unwanted variation when that variation
  follows the groups. There such sets were called in 29.6 % (genes), 11.2 and 13.3 % (proteins), 8.1 % (metabolites)
  with the adjustment, against 2.0 %, 7.3 and 11.2 %, 10.1 % without it and about 2 % had the true factors been known.
  The excess comes from the factors being estimated, rises with set size, and is confined to such sets. No real data.
- **What changed from the pipeline before.** Nothing by default. New: one more file with a saved integration input
  where RUV was applied; a line on every correlation-adjusted result saying whether it was adjusted; the option.

## Hand-offs

- **The follow-up measurement (VAL-160 (2); a candidate, not this ticket).** What it needs, from the results audit: per
  set, the correlation `camera` estimates and the common shift, written out; the dependence on how strongly the
  unwanted factor follows the groups and on the loadings; a candidate remedy evaluated before any claim. To know: the
  existing overlap warning (adjusted share of the factor explained by the removed factors at least 0.5) would not
  fire in the simulation's confounded cell (median share 0.25 to 0.29). At a build whose default is off,
  `simulate-k184.R` line 977 `designTests("3", adjusted, NULL, TRUE)` must read
  `designTests("3", adjusted, TRUE, TRUE)`; a run so changed is a new declared run with a new script hash, not a
  reproduction of S, which is reproduced at package 81ee3052. The per-replicate tables of run 45353248796 exist only as its artifacts
  (about 86 MB compressed, kept 90 days from 2026-11-24).
- **KIT-193's pair R-I1** (the removed factors' span) waits for this ticket's merge. It binds to
  `integrationInputRemovedFactors()` (`R/utils_workflow_integration_input.R`) and
  `integrationEnrichmentRemovedFactors()` (`R/utils_integration_enrichment.R`); the tests K1, K2, K5 and K6 of
  `tests/testthat/test-integration-removed-factors.R` hold the cases. The default being off does not touch it.
- **The text owed to KIT-181** (`tickets/kit/audits/OWED-TEXT-KIT-181-from-KIT-184-20261123.md`) **is stale after
  VAL-160**: it says the resolver's default follows the model's input and names the suffixes of both forms. I do not
  edit that file. What it should say now: "A Gene Ontology run of a test that has a design
  (`INTEGRATION_ENRICHMENT_DESIGN_TESTS`: `camera`, `camera_floor`, `camera_pr`; `camera_pr` is in the constant but is
  not offered, so in practice `camera` and `camera_floor`) follows KIT-184's setting `adjust_ruv_factors` exactly as KEGG and Reactome do (VAL-62, VAL-160): the request carries
  it; absent, it is off whatever the model was trained on; `values` is built with the covariates
  `integrationEnrichmentRemovedFactors()` returns when it is on and with `covariates = NULL` otherwise; the entry keeps
  the resolved setting; the head key gains `|adjusted` for the adjusted form and nothing otherwise; the result carries
  the label and the sentences of KIT-184 (`kit-184/wording.md`). The permutation and the pre-ranked test refuse `TRUE`.
  Whichever of KIT-181 and KIT-184 merges second adds the test: on a model trained on RUV-corrected values a Gene
  Ontology run of `camera` has `intercept, factor` and the label 'no' when the setting is absent, and `ruv_W` columns
  only when the request set it to `TRUE`; an entry replays in the form it was saved in. A Gene Ontology entry with the
  adjustment brings the record's clause `methods_enrich_limits_ruv` as KEGG and Reactome entries do."
  Checked against the code at the merge with staging (d7928ac9, `R/utils_integration_enrichment.R`): line 862 is
  `INTEGRATION_ENRICHMENT_DESIGN_TESTS <- c("camera", "camera_floor", "camera_pr")`; line 853,
  `INTEGRATION_ENRICHMENT_TESTS_OFFERED`, is `camera`, `camera_floor`, `permutation`, `gsea`, without `camera_pr`; a
  request that names `camera_pr` is refused with `integration_enrichment_test_not_offered` (line 1236) before the
  setting is read, so no `camera_pr` result can be made with or without the adjustment. The constant, the key's suffix
  rule and the record's texts cover it should it ever be offered. An earlier version of this hand-off named only
  `camera` and `camera_floor`; the text above now agrees with the one in `tickets/kit/KIT-181.md` r1.5.
- **The Manuscripts deliverable is not written by me.** It needs: the Methods fragment and the limits clause
  (`methods_enrich_limits_ruv` and the fragments of `wording.md`, after the owner's wording pass); the simulation's
  summary figures from `s/` (the deciding rates with their intervals, `ke1-checks.tsv`; the limit per data set and
  arm, `ke4-checks.tsv` and `error-rates.tsv`; by set size, `error-rates-by-size.tsv`; power, `power.tsv`; what the
  removed factors caught, `caught.tsv`; the feature-level rates, `feature-rates.tsv`); the statements of ticket 15.2 as
  restated in its 26.5; the lifecycle disclosure sentence; "one simulated design, no real data".
- **The owner's wording pass**: every user-facing text of this ticket is a draft, verbatim in `wording.md` (32 keys,
  the setting's line, the documentation of `integrationEnrichFactor()`, the starter workbook's paragraph, NEWS).
- **The stored reference** `tests/testdata/gui-entry/prot_dia_limpa.json`: re-frozen by this ticket after its control
  (VAL-178 (2); below). Three references that fail on staging already will list `removed_transcriptomics.parquet`
  when they are re-frozen under their own tickets: `tests/testdata/gui-entry/trx_screen.json` (13 expectations fail on
  staging) and the captures `trx_main` and `trx_managed` of the summary-module tests (`f/moved-references.txt`).

## Where the tests are (differs from ticket section 10 in location only)

K1 to K3, K5 to K12 and K14 to K16, and the placeholder test, are in the one new file
`tests/testthat/test-integration-removed-factors.R` (with its own fixtures); the ticket places K6, K7, K14 and K15 in
KIT-183's file and K10 to K12 in the record, tab and workbook files. K4 is in the per-route writer files through one
helper, `expectIntegrationRemovedFactors()` (`tests/testthat/helper-integration-input.R`). K13 is in the MOFA2 lane's
file.

Existing tests restated because the ticket changes their behaviour by design: `test-integration-enrichment-tests.R`
(KIT-183 T10: a result of a test with a design names its form; T11: the RUV note beside the selector stays with the
tests that have no design); `test-integration-input-matrix.R` (the file list of an input with a corrected matrix);
`test-integration-record.R` and `test-integration-gui.R` (a fixture whose corrected matrix equals the normalised one
now raises the classed warning); `test-workflow-integration-input.R` (the stand-in test after the correlation filter
expects `not_determined`).

K4, not built as the ticket words it: the lipidomics declared design "with and without a missing value in the
fixture" has one fixture only (without), and the principal angle to the DA's stored factors is printed for it; the
transcriptomics test prints no angle to the state's RUVSeq factors and asserts instead that the written factors span
RUV-III's own W; **classic proteomics with a declared design and metabolomics with a declared design have no K4
fixture and were not run** (G1 holds on every RUV route that has a writer fixture; the writer is route-agnostic and
both use RUV-III-C as G2 showed, so exactness on those two is inferred, not shown; implementation audit K184-I07).

Other departures recorded in the ticket (20.1): 6.6's "residual degrees of freedom are 16" is 15 on
metabolomics_GCMS; the workbook does not send the setting with `test = "gsea"` or `"permutation"` (a refusal listed
under `skipped` would be the stricter reading).

## Inferred, not observed

- That Module CI's failure on the branch's pushes (job "module-ci impacted (transcriptomics)") is not this ticket's:
  its failing expectations are the corpus manifest's SHA-256 recorded on the laptop, which fail identically on base
  and head in R1 and R2.
- That the classic proteomics route's larger residual (3e-09 against 1e-12 elsewhere) is the rounding of its level
  restoration.
- The mechanism behind the limit beyond what the paired arms show (`results.md`, "K-E4 in cell B").
- That no other stored reference will move: F ran every test file, so any reference that passes on staging and lists
  an integration input's files with RUV applied would have shown; references that already fail on staging are not
  covered by that.

## Not checked

- The tab in a browser (the module is driven with `shiny::testServer` and `renderTags`).
- Real data; any generator, sample size or correction other than the simulation's.
- The distance between the DA's stored factors and the removed ones on a declared design with missing values.
- That the CRAN archive's ruv 0.9.7.1 tarball is the one the lock's library was installed from (G0).
- The package at the merged head on Windows beyond the one new test file (the bundle lane ran that file only; the full
  suite ran on ubuntu-22.04).
- The MOFA2 lane at the merged head (M2 ran at 5ba85c1f, before the merge with staging; the merge changed no file of
  the integration enrichment).
- Cell C of the simulation (not run: KIT-183's outcome is O1).
- A second reading of the texts by the owner: all are drafts.

## After the implementation audit (2026-11-25; owner's decision VAL-178)

VAL-178: (1) the claim narrowed to what was measured, with the 9.4 % and 8.6 % figures printed, the limit printed in
the record whenever any entry is adjusted, an approximate set of factors stated with its share; (2)
`tests/testdata/gui-entry/prot_dia_limpa.json` re-frozen, that file only, after the control; (3) one further bounded
run beyond the spent correction pass, inside the ticket's unchanged ceiling; (4) VAL-160 governs this ticket's
proteomics and metabolomics figures.

### Commits

| Commit | What |
|---|---|
| 5352a498 | tests only, for K184-I01 to I06 (nothing under `R/`) |
| 33facbb5 | the fix: texts, the approximate fragment, the clause printed for any adjusted entry; the record repairs I07 to I09; `wording.md`; `gui-entry-control.R` |
| d7928ac9 | the merge with `v0.9-staging` 32766889 (one conflict, the methods ledger: staging's rows to VAL-182 first, then this ticket's three `VAL-<merge>` rows; NEWS merged without conflict; nothing else conflicted) |

### Tests first: each changed or new expectation, what it gives before the fix (5352a498) and after (33facbb5)

All in `tests/testthat/test-integration-removed-factors.R`.

| Finding | Expectation | Before the fix | After |
|---|---|---|---|
| I01 | K10: the clause holds "of the kinds the check decided on (random sets, correlated blocks of features,", "the first three also against the model's main factor), were called at the 5 % level in at most 5.4 %", the most variable features' sets "up to 9.4 % (6.8 % without the adjustment, 9.8 % had the true unwanted factors been known)", the planted sets "up to 8.6 % (up to 7.0 % ..., up to 8.5 % ...)" | absent: 4 failures | present |
| I01 | K10: the clause does not hold "any set against a covariate" | present: fails | absent |
| I01 | K11: the help holds "(at most 5.4 %) for the kinds of sets with no link to what was tested that the check decided on; ... (up to 9.4 %, 6.8 % without it)" and "(up to 8.6 %, up to 7.0 % without it)" | absent: 2 failures | present |
| I01 | the placeholder block: the starter workbook's paragraph holds the same two phrases | absent: 2 failures | present |
| I02 | K10, a model trained before RUV, record with the adjusted entry: the limits paragraph contains `methods_enrich_limits_ruv` and "in 8.1 to 29.6 % with the adjustment" | absent: 2 failures (the head asserted the absence) | present |
| I02 | K10, the same model, record with only the unadjusted entry: no clause, Methods as before | passes | passes |
| I02 | K10: the clause opens "The correlation-adjusted tests can add to their model the unwanted factors the RUV correction removed from each data set" | opens "For a model trained on RUV-corrected values": fails | passes |
| I03 | new block "K10 (K184-I03)": `integrationMethodsCovariates()` on an entry with approximate factors (share 0.9) gives "; approximate: 90.0 % of what the correction removed lies within these factors (runs were removed after the correction), so they do not match the corrected values exactly", directly after "8 residual degrees of freedom)" and before the words for a model trained before RUV | absent: 3 failures | present |
| I03 | the same block: an exact entry and an unadjusted one hold no "approximate" | passes | passes |
| I03 | K10: the clause holds "for a model trained on RUV-corrected values and features measured in every run where no run was removed after the correction, this is identical to fitting that model" | absent: fails | present |
| I04 | K11: the note holds "the regression is exactly that of the same model (the factor and these 2 unwanted factors) fitted to the corrected values, with 8 residual degrees of freedom; the test's variance trend reads ..." and not "this gives exactly the result" | 2 failures | passes |
| I05 | K11 and the starter: "(with as many factors removed as were planted; with four removed where two were planted it found fewer)" | absent: 2 failures | present |
| I06 (a), (b) | K10: the clause holds the proteins' two compositions and the simulated correction, and not "11.2 to 13.3 %" | 4 failures | passes |
| I06 (c) | "yes (2 factor(s))." (K11; K10's label and ranking cell), "2 removed unwanted factor(s)" (K8's refusal, K10's ranking text, the design sentence), "4 removed factor(s)" (K15) | "factors": 6 failures | passes |
| I06 (d) | the starter holds "In that check it did not protect pathways ..." | "It does not protect": fails | passes |
| I06 (f) | K11: the four texts that state the limit give its numbers; no tab text says "(see Methods)" | 5 failures | passes |

### The control before the re-freeze (as the audit register gives it, section 5; VAL-178 (2))

It runs on a runner (`gui-entry-control.R`, in the pattern of `kit-199/control.R` for the same lane):

1. A base dump of the lane `prot_dia_limpa` at the staging commit the merge takes (32766889). It must be byte-equal
   to the stored file; this shows the dumper reproduces the stored reference on that runner.
2. Two head dumps at the merge commit (d7928ac9), identical to each other.
3. The head dump against the stored file must differ in exactly four added keys,
   `"integration/proteomics/removed_proteomics.parquet": "<present>"` in `summary_save`, `summary_copy`,
   `summary_report` and `summary_export`, and in nothing else.

Then the head dump is committed alone as the new reference, and the lane must pass in the run. If the control shows
anything else: stop and report; the file is not re-frozen and nothing is widened.

One thing about "byte-equal to the stored file", read in the helper (`tests/testthat/helper-gui-entry.R`,
`guiEntryExpectLane()`): a frozen file holds, beside its steps, the commit, the date and R's version of its freeze, so
no later dump can equal the whole file byte for byte. The control compares the `steps` member, written alone by the
helper's own writer, by SHA-256, as KIT-199's control did. The whole-file hashes are written too (`sha256.tsv`).

### The one further run, declared before it is sent

- **Purpose.** The last run of this ticket (VAL-178 (3)): the tests-first evidence for I01 to I06, the regression at
  the merge, the control, and, only if the control's three verdicts hold, the lane passing on the re-frozen file. An
  unfavourable valid result is reported, not repaired.
- **Queue.** Checked with `gh run list` before each send (recorded in `lifecycle.md`).

| Part | Request | Commits | What it decides | Expected; ceiling (job-minutes) |
|---|---|---|---|---|
| a, tests first | `test-integration-removed-factors.R`, head only | 5352a498 (tests only) | that the restated expectations fail before the fix at the values of the table above | 15; 45 |
| b, the merge | 29 files (R2's 25 and the audit's `test-kit206-prot-norm-export-guard.R`, `test-omic-summary-module.R`, `test-omics-summary-core.R`, `test-gui-entry-proteomics.R`), 3 shards per side, R CMD check on both sides; scripts: `gui-entry-control.R 32766889`, then `base-store-job.R a5229552 mutations` | head d7928ac9 (the merge) against base 32766889 (staging) | A2 to A11 and the repairs at the merge; nothing worse than staging but for the known four expectations of lane `prot_dia_limpa` (the reference is re-frozen only after the control); check unchanged; the control's three verdicts; the base-store blocks and the four mutations | 330; 600 |
| c, the re-frozen lane (sent only if the control's three verdicts hold) | `test-gui-entry-proteomics.R` | head = the commit that holds the head dump of `prot_dia_limpa.json` alone, against base 32766889 | that the lane passes on the re-frozen file: the file's verdict `same` as staging (189 passed, 4 failed in lane `prot_dia`, which fails on staging already) | 90; 180 |
| **Total** | | | | **435 (29.0 core-hours); 825 (55.0 core-hours)** |

- **Budget.** 205.16 of 303.33 core-hours used before this run; with the ceiling 260.2. No row of ticket section 11
  holds this run: it is the further run of VAL-178 (3), inside the ticket's unchanged ceiling.
- **Not in it, and why.** The MOFA2 lane and the bundle lane are not repeated (the audit's regression scope does not
  ask for them: K13's assertions and the cross-platform block read nothing the repairs change); no full suite (F at
  8cb536cc stands for the rest, staging's code since being disjoint from this ticket's, audit section 8).
- **Mechanics.** A request refused before a job runs, or a script that cannot run (no verdict written), is an
  instrument failure (VAL-54): fixed and sent again inside the ceiling. A verdict or a test result is a result.

### The run's results (2026-11-25; `runs.tsv`, `lifecycle.md`, `r3/`)

| Part | Run | Result |
|---|---|---|
| a, tests first | 45474659283 at 5352a498 | 36 expectations fail in 6 blocks, 0 errors, as declared (K10 15, K11 12, the starter's paragraph 4, the new I03 block 3, K8 1, K15 1) |
| b, the merge | 45474665699, d7928ac9 against staging 32766889 | the new file 526 passed, 0 failed, 3 skipped (the base-store blocks, which pass in the base-store job); 19 files `same`, 6 `count-changed` (this ticket's, all passing), 2 `message-changed` (the transcriptomics captures failing on staging already); `test-gui-entry-proteomics.R` worse by the four expectations of lane `prot_dia_limpa`, as declared; check identical on both sides; the four mutations caught |
| the control | script 1 of part b | **base dump equal to the stored file: TRUE. Two head dumps identical: TRUE. Head differing from the stored file in exactly the four added entries: TRUE** (`r3/control/verdict.tsv`, `differences.tsv`, `sha256.tsv`) |
| the re-freeze | commit 0c73ee7c | `tests/testdata/gui-entry/prot_dia_limpa.json` alone: the runner's head dump as written |
| c, the re-frozen lane | 45490115851, 0c73ee7c against staging | `same` as staging: 189 passed, 4 failed on both sides (lane `prot_dia`, failing on staging already); lane `prot_dia_limpa` passes |

Cost: 401.0 job-minutes (26.73 core-hours); declared 435 expected, 825 ceiling. Ticket total 3,478.4 job-minutes =
231.89 core-hours of 303.33.

### The audit's findings

| Finding | State | Evidence |
|---|---|---|
| I01 (major): the claim wider than what was checked | closed (VAL-178 (1)); the words are drafts for the owner | 33facbb5; parts a and b; `wording.md` |
| I02 (major): the limit missing from the record of a model trained before RUV with an adjusted entry | closed | 33facbb5; K10 fails in part a and passes in part b |
| I03 (major): an approximate set of factors not stated in the record | closed | 33facbb5; the new block fails in part a and passes in part b |
| I04, I05 (minor): "exactly the result"; "more true terms" without scope | closed; drafts for the owner | 33facbb5; `wording.md` |
| I06 (minor): wording items | applied as recommended, each listed in `wording.md` as "owner sets final wording"; (c) only in the label and the count phrases (the running sentences still read "the 1 unwanted factors" at one factor) | `wording.md` |
| I07 (minor): records | closed for ledger rows 1 and 2, `feasibility.md` and this file; VAL-93's superseded clause is the coordinator's to note at the merge | 33facbb5 |
| I08 (minor): `trx_screen.json` left out | closed | `f/moved-references.txt`; Hand-offs above |
| I09 (minor): the line of the frozen script | closed | `results.md`; Hand-offs above |
| I10 (owner): the re-freeze | decided (VAL-178 (2)) and done after the control | 0c73ee7c; `r3/control/`; part c |
| I11 (owner): VAL-98 against VAL-160 | decided (VAL-178 (4)): VAL-160 governs this ticket's figures; recorded in `results.md` | VAL-178 |
| I12 (record): F over its row's ceiling | recorded; the times are with the owner | "Core-hours" above |

### Acceptance after this run

A1 to A16 and A18 as in the table above, now at the merge with staging 32766889 (part b) for the criteria the repairs
touch (A8 to A10). **A17** (the full suite, no new failure): F (45402062487, at staging 8cb536cc) found one new
failure, the reference that the owner has had re-frozen; part c shows it closed. F was not repeated: the audit's
regression scope at the newer staging head ran in part b, and staging's code since F is disjoint from this ticket's
(audit section 8).

### Not checked by this run

- The MOFA2 lane and the bundle lane (three platforms) at the final head: last run at 5ba85c1f and its lane commits,
  before these text repairs. K13 and the cross-platform block read nothing the repairs change (inferred, not run).
- A full suite on the final head.
- A record of a real run with an approximate set of factors (the test builds the entry from a unit fixture).
- The control script was written for this run and ran once; nothing but its own checks (the frozen file's steps equal
  the head dump's; the checkout's reference untouched) and my reading of its outputs vouches for it.
