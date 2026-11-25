# KIT-256 outcome: the technical-replicate selection is emptied after each successful Assign

Ticket KIT-256 r1.2 (owner's decision VAL-249). Branch `kit-256` from `v0.9-staging` 63143aa3.

| Commit | What |
|---|---|
| 88f60b00 | the browser script `browser-selection.R` only (R/ as staging 63143aa3): G0's binding |
| e65ad78a | G0 records |
| b45c8779 | tests only: `tests/testthat/test-kit256-tech-rep-selection-cleared.R`; the browser script's notification recorder repaired (below) |
| 37451255 | R1 records |
| 57238dcd | the fix, `R/mod_omic_design_builder.R` only |
| (records commit) | GB and R2 records, NEWS, ledger row `VAL-<merge>`, `wording.md`, this file |

## Runs (runs.tsv)

| Run | Commit | Purpose | Job-min | Verdict |
|---|---|---|---|---|
| 45617875300 | 88f60b00 | G0, base browser values | 12 | stop rule 1 not met: the box keeps the runs |
| 45620278822 | b45c8779 | R1, failing-first | 12 | E1, E2, E3, G3 fail at the base values; G1, G2, G4 pass |
| 45622591472 | 57238dcd vs b45c8779 | R2, 18 files and check | 142 | every predicted verdict held |
| 45622587312 | 57238dcd | GB, head browser values | 8 | 24 of 24 predictions held; stop rule 2 not met |

Total 174 job-minutes = 11.6 core-hours (expected 25.33, ceiling 110.0). No correction pass and no instrument rerun.

## Observed

**G0 (base, run 45617875300; `g0-base-values.tsv`, `g0/g0-compare.tsv`).** Phase 0 passed: Chrome
`/usr/bin/google-chrome`, the app answered, AppDriver attached by URL, the app process's `warn` was 1, the box offered
S1 to S4 (shinytest2 0.5.1, chromote 0.5.1, shiny 1.14.0). After Assign of {S1, S2} (`lowest`) the summary showed one
biological replicate {S1, S2} and the box read `["S1","S2"]` (step 3). The box then held S1 to S4 (step 4); Assign
gave one biological replicate of four runs, and the box kept them. A refused Assign ({S1}) left the box `["S1"]` and the
summary unchanged. The saved design had one key `A_1={S1,S2,S3,S4}`. App log: 0 warnings, 0 errors. 18 of 19
predictions held (the G0 gate has 19: 14 shared with GB and 5 of its own); the one that did not is the
instrument note below.

**R1 (run 45620278822; `r1/`).** 14 passed, 6 failed, 0 errors. Each failure is an assertion at the ticket's base value:
E1 lipidomics and transcriptomics, last message `c("S1","S2")`; E2 proteomics `c("S1","S2")`; E3 lipidomics and
transcriptomics, the box after the first Assign `S1 S2` and one key of S1 to S4; G3 `c("S1","S2")`. G1 (too few runs;
two groups), G2 and G4 passed.

**R2 (run 45622591472; `r2/compare.tsv`, `r2/check-compare.txt`).** 17 files `same`. The pre-existing failures are
identical on both sides: `test-gui-entry-proteomics.R` 4, `test-gui-entry-transcriptomics.R` 13,
`test-module-ci-trx-design.R` 1 error, `test-workflow-ruv-pearson-diagnostic.R` 2 errors. The new file is `fixed`
(base 14 passed 6 failed, head 20 passed 0 failed). R CMD check: 0 errors, 3 warnings, 3 notes on both sides; findings
identical (none only on head). The regression list was recounted at the branch point 63143aa3: grep 1 of section 10
gives the same 13 files, plus the 4 GUI-entry lanes.

**GB (the fix, run 45622587312; `gb-head-values.tsv`, `gb/gb-compare.tsv`).** Phase 0 passed as in G0. Step 3 the box
read `[]`; step 4 sent {S3, S4}; Assign gave two biological replicates, {S1, S2} and {S3, S4}, and the box read `[]`; the
refused Assign ({S1}) showed the warning "Please select at least two samples.", left the summary and kept the box
`["S1"]`; the saved design had two keys `A_1={S1,S2} A_2={S3,S4}`. Second session (no-change path): {S1, S2} marked,
then {S1, S2} again in the same mode, Assign showed "Assigned 2 samples as technical replicates", the summary was
unchanged and the box read `[]`. App log: 0 warnings, 0 errors.

**Instrument note (VAL-54).** In G0 the notification recorder read each notification's text when Shiny inserted the
element, before Shiny filled it, so `step5.acted` recorded `shiny-notification :: ×` and its prediction (`warning.*at
least two`) did not hold. The click's effect itself was observed (a notification element appeared; the summary and the
box were unchanged). No acceptance value depends on that text. The recorder was changed in b45c8779 to read the class
and text when the notes are read; no expected value, pattern or stop rule changed. G0 was not re-run; GB ran with the
repaired recorder and read the texts.

## Inferred

- The flush-order reading of the ticket (section 2) holds in a real browser: on base the selection observer's re-send
  keeps the box filled; with the flag the box is empty after the Assign (GB step 3, step 4), and the direct clear empties
  it after an Assign that leaves the design identical (GB's second session, in the browser). Under testServer, G3 as
  first written read the last message of the whole session, which the first Assign's clear already supplied, so it
  would have passed with the direct clear (`R/mod_omic_design_builder.R:667`) removed (AUDIT-KIT-256-impl-20261126
  I01). G3 now reads from a mark set just before its second Assign (83d6ffd7), so it fails on the old code and guards
  the direct clear; R3 below.
- Every omic and route uses the same server and input id; proteomics was observed at module level (E2), not in a
  browser.

## Acceptance criteria

| # | Criterion | Met | Evidence |
|---|---|---|---|
| A1 | G0 on base: step 3 `["S1","S2"]`; step 6 one key of four runs | met | `g0-base-values.tsv`; run 45617875300 |
| A2 | E1, E2, E3, G3 fail at G0's values; G1, G2, G4 pass; all pass on the fix | met | R1 45620278822; R2 new file head 20/0 |
| A3 | R2: every file `same` except the new file; pre-existing failures identical | met | `r2/compare.tsv`, run 45622591472 |
| A4 | GB: step 3 `[]`, step 4's union {S3, S4}, two keys; no-change Assign leaves `[]`; refused Assign keeps the box | met | `gb-head-values.tsv`; run 45622587312 |
| A5 | NEWS, ledger, `outcome.md` with "For the change report" | met | `NEWS.md` (`### Fixes`), ledger `VAL-<merge>`, this file, `wording.md` |

## Implementation summary and recorded divergence

- `R/mod_omic_design_builder.R`: the state list gains `tech_rep_selection_clear = shiny::reactiveVal(FALSE)` (:436);
  the selection observer, after its isolated read, sends `character(0)` for `tech_rep_samples` and resets the flag
  when it is set (:532-537); the Assign observer, on `ok`, sets the flag only if the design changes, sets the design,
  then sends `selected = character(0)` (:660-667). A refused Assign does neither.
- Divergence from `conventions/R-shiny.md:309-310` (populate reactive values in the function that renders dependent
  content): the decision to clear is made in the Assign observer and carried to the selection observer by one
  `reactiveVal`, read under `isolate` (no dependency, no loop). The single-observer alternative (filtering marked runs
  out of the kept selection) was rejected in the ticket (section 3), and `session$onFlushed` is unproven.
- The tests copy KIT-222's `kit222ProtOwner` and `kit222ProtBuilder` as `kit256ProtOwner`/`kit256ProtBuilder`, and
  KIT-222's clicks (`kit222BuilderClicks`) split at its technical-replicate clicks into `kit256SetupClicks` and
  `kit256SaveClicks` (the setup and the save are what this file needs; the clicks in between are its own rows).

## Not checked

- Proteomics' own UI and the DIA-limpa builder in a browser (same server and input id; E2 at module level only).
- Metabolomics in the new file (the shared server under the metabolomics configuration; the GUI-entry metabolomics
  lane ran `same`).
- Windows.
- How often users built a sample in two Assigns, or merged two samples this way (VAL-249: not known).
- The help sentence is a draft only (`wording.md`); not in the interface.
- Head with the direct clear (`:667`) removed was not run; G3 failing there is inferred from R3's base side.

## Unrelated finding recorded (I05; ticket F1, no change here)

`selected_runs` (Assign Metadata) carries over the same way: nothing clears it after Assign Metadata
(`R/mod_omic_design_builder.R:531`, `:544-546`; the observer `:650-656`). Adding the next group's runs to a leftover
selection and pressing Assign moves the earlier runs to the new group too, numbered again from `replicate_start`, and
also clears those runs' technical-replicate marks (`R/func_omics_design_shared.R:622-623`). It is visible to the user
(design table, groups, contrast choices, technical-replicate summary), not a silent merge, and legacy on `main`
(`main:R/mod_lipid_design_builder.R:594`, `:607`). Recorded only; no change without a ruling
(AUDIT-KIT-256-impl-20261126 I05).

## Merge-train declaration (I06)

- Branch state: `kit-256` merged with `v0.9-staging` 42846c0e (3ae09360; NEWS and ledger appends kept both, staging's
  first; nothing else conflicted). Paths changed against staging: about 30, the longest under 100 bytes.
- Scope: the 17 files of ticket section 10 and `test-kit256-tech-rep-selection-cleared.R`.
  - The new file is `only-on-head`, 20 passed 0 failed.
  - The 17 are predicted `same`, base against head; check identical.
  - Do not assert 189/4 for `test-gui-entry-proteomics.R` (unmeasured since KIT-250 landed): read it from the run.
    `test-gui-entry-transcriptomics.R` 16/13 carries train 6's changed messages (VAL-236) on both sides.
- If KIT-247 lands first, in the same train or before: add `test-kit247-design-resave.R`, `test-kit247-import-config.R`
  and `test-kit247-resave-no-session.R`, predicted at KIT-247's own head values. If KIT-256 lands first, KIT-247's
  train runs the KIT-256 file.
- Conflicts at the train as in train 7: the ledger keeps staging's rows, then KIT-256's row as `VAL-<n>` (this file's
  `VAL-<merge>` replaced to match); NEWS keeps every entry under the one unreleased heading.
- At landing (I04): the KIT-256 row's status cell is in the train form with `<date>`, `<n>`, `<id>` to fill; VAL-249's
  status cell becomes "decided; implemented by KIT-256 (`VAL-<n>`), merged into `v0.9-staging` <date> by merge train
  <n>".

## For the change report

### Main part: legacy (proteomics, lipidomics, metabolomics)

**What was wrong.** In the design builder, after 'Assign Technical Replicates' the "Select Technical Replicates" box
kept the runs just marked. A user who then added the next biological sample's runs to the box and pressed Assign
marked all of them as technical replicates of one biological sample. Nothing refused it: KIT-222's refusals allow
re-selecting a marked sample's runs together with new runs (the documented way to add runs to a sample). Who was
affected: any user who marked technical replicates for more than one sample in a row without first clearing the box,
on any omic. How often is not known.

**Legacy, confirmed on `main` 436cdbe8 with `git show` (main not checked out).** The lipid builder re-sends the kept
selection after every design change (`main:R/mod_lipid_design_builder.R:586-610`: `selected_tech <-
input$tech_rep_samples` under `isolate`, then `updateSelectizeInput(session, "tech_rep_samples", choices = available_runs,
selected = selected_tech)`), and its Assign observer sets the design without touching the selection (`:913-962`).
Metabolomics has the same lines (`main:R/mod_metab_design_builder.R:586-610`, `:913-962`). Proteomics the same
(`main:R/mod_prot_design_builder.R:455-475`; Assign from `:789`). On `main` the edit also gives the four runs one number.

**Code references.** Base 63143aa3: `R/mod_omic_design_builder.R:520-541` (the selection observer re-sends the isolated
selection, `:530`, `:539`), `:650-656` (the Assign observer; the selection not touched). Fix 57238dcd:
`R/mod_omic_design_builder.R:436`, `:532-537`, `:660-667`.

**The merge, traced with G0's values.** Group A, runs S1 to S4 numbered 1 to 4. Assign {S1, S2} (`lowest`): one
biological replicate {S1, S2}; the box still shows `["S1","S2"]`. The user adds S3 and S4 and presses Assign: one
biological replicate of four runs, saved as one key `A_1={S1,S2,S3,S4}` (run 45617875300).

**What changed in results.** After the fix the box is empty after each successful Assign (run 45622587312: the same
clicks give `A_1={S1,S2}`, `A_2={S3,S4}`). For an analysis whose design was built through the leftover selection, a
corrected design changes results: transcriptomics no longer sums the counts of two samples; lipidomics, metabolomics and
proteomics no longer block or RUV-adjust them as one sample; the integration no longer averages them. Designs built
without this sequence are unchanged (the edit is not changed). The reverse direction remains possible by user action:
to add runs to a sample already marked, all of that sample's runs must be selected again with the new ones; runs
selected alone become a separate biological sample (technical replicates analysed as separate samples,
pseudo-replication; p-values too small). A user who built one sample in two Assigns got one sample before (the leftover
selection did it) and gets two after unless they re-select. NEWS states this; a help sentence is a draft (`wording.md`).

**Can an affected design be recognised?** No, not from the saved design. The merged key holds technical-replicate
numbers 1 to n with no repeat, and all its runs are marked, so KIT-222's save checks D1 and D2 see nothing. A user must
compare the design's technical-replicate summary with the laboratory's sample sheet. Recompute: correct the design and
run every step again.

**Evidence.** G0 45617875300 (base browser, 4-run one-group lipid table); R1 45620278822 (module level, lipidomics,
transcriptomics and proteomics `dia_long`, 8-run two-group design); R2 45622591472 (17 regression files and check,
base against head); GB 45622587312 (the fix in the browser, incl. a refused Assign and an Assign that leaves the design
unchanged). Limits: one omic in the browser (lipidomics); Linux, Chrome; synthetic designs only.

### Appendix: code new to 0.9 (transcriptomics)

The transcriptomics design tab uses the same shared builder (`omicGuiTrxDesign`, `wide_assay_list`); `main` has no
transcriptomics. Its behaviour before the fix was the same as above (R1: E1 and E3 transcriptomics fail on base at the
same values), and the fix applies to it identically (R2: the new file passes on head). Consequence before the fix: the
counts of two biological samples were summed into one. The same remarks on recognition and recompute apply.
