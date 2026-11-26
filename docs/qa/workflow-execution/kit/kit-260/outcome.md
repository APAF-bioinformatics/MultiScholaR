# KIT-260: outcome

Ticket r1.3 (owner's decisions of 2026-11-25, methods ledger VAL-241: OD-1 (c), OD-2 (b), OD-3 (a), OD-4 (a), stored
expectations yes to all three). Branch `kit-260` from `v0.9-staging` e0747fc198817596e3e66397a87e06b2f869c6a7.

| Commit | What |
|---|---|
| 42e87abb | G0's records brought over from branch `kit-260-g0` (records only) |
| 63e29c04 | tests only: `tests/testthat/test-kit260-itsd-default-pattern.R` and the approved stored expectations; nothing under `R/` |
| 288abf49 | the fix and NEWS (`R/` and `inst/` identical to the P0 scratch commit 406d8229 on branch `kit-260-p0`) |
| (this) | records, wording drafts, ledger row `VAL-<merge>` |

Runs: `runs.tsv`; tables in `g0-45592957172/`, `p0-45603194585/`, `p0b-45612339270/`, `r1-45618785756/`,
`r2-45621200497/`. G0 in detail: `g0-outcome.md`. Wording as implemented: `wording.md`.

## What was changed (observed in the diff)

- **OD-1 (c)**: both MS-DIAL defaults are `^IS_|_d[0-9]+$|ISTD|^IS-|\(d[0-9]+\)` (`R/func_metab_import_detection.R:210`,
  `R/func_lipid_import_detection.R:235` at 288abf49; base `:207`, `:234`). No other format's default changed.
- **OD-2 (b)**: the ITSD declaration's kernel returns, beside its checkpoint, the rows used per assay (identifier and
  annotation, read from the object before the step by the identifiers the S4 method stored) and the identifier column
  (`omicsItsdRowsUsed`, `R/utils_workflow_omics_declarations.R:338`, called `:437`; engineering E-1). The tab logs one
  line per assay naming each row "annotation (identifier column identifier)" with the count, and shows one notification
  (`omicsItsdReportResult`, `:408`; called in the managed route `R/utils_workflow_omics_actions.R:362` and in process
  `:487`). The record's ITSD step outcome gains `itsd_features_per_assay` (`:1290`; E-2), read from the step's saved
  product as the RUV outcomes are. The import box's help text says where the pattern is matched, that it is
  case-sensitive and that matches are removed (`R/mod_omic_import_feature_table.R:182`, text in `omicGuiItsdTexts()`).
  Removal stays on. The S4 ITSD methods are not changed; no state's values or `args` change.
- **OD-3 (a)**: where a run's pattern-mode matches include a name containing `w/o MS2:` or a name matched only through
  an unanchored `IS-` alternative (`omicsItsdLegacyPatternMatches`, `:315`), the step's outcome carries `itsd_note`
  (`:1293`); the study-parameters file writes it as `* Note:` in the ITSD block when the analysis record carries it
  (`R/func_general_filemgmt_study_parameter_writer.R:299`, `:500`; `omicsItsdRecordNote`, `:389`); both reports append it
  to the ITSD paragraph; the integration Sample Matching tab states it where marked rows are such names. No replay
  marker (the stored pattern governs). The `IS-` case counts a name that contains 'IS-' after its first character and
  that no other alternative of the stored pattern matches; a name that starts with 'IS-' (a label the corrected default
  also matches, as `^IS-`) is not counted. Read literally, 'matched only through the unanchored IS-' would count it too
  (E8's object: 10 rows per assay instead of 9) (I08). Correction pass (I02 (a)): when the ITSD step's product names no
  rows (a product saved by a build before KIT-260), the record reads them from the stored states
  (`omicsItsdRowsFromStates`: the ITSD state's stored identifiers, their annotations from the state the step read, two
  columns a page at a time, the note's rule with that state's stored pattern), so an older project reopened and
  summarised without re-running normalisation gets the note in the record, the study-parameters file and the reports.
- **OD-4 (a)**: INT-3's rule unchanged. The Sample Matching tab names, per data set, the marked rows "display name
  (identifier)" (`integrationInternalStandardNames`, `R/utils_workflow_integration_model.R:659`; facts in
  `R/utils_integration_gui.R`); the model product and the record's model outcome carry `internal_standards_marked`
  only when a data set has marked rows (`:1398`, `:1824`), so a model without them is saved and recorded as before.

## Runs and what each showed (observed)

- **G0** (run 45592957172, at 16112917 + script): D0 to D5 and D7 as predicted; D6 recorded. See `g0-outcome.md`.
- **P0** (run 45603194585, head only, scratch 406d8229, 129 files): 14 files with failures. Five of them fail
  identically, by count and by every message, in merge train 4's run 45586240983 (same staging code, base and head):
  `test-omic-summary-module.R` (28), `test-omics-summary-core.R` (25), `test-workflow-kit-harness.R` (1),
  `test-workflow-operation-registry.R` (2), `test-workflow-omics-binding-steps.R` (1). For the 8 files not in that run,
  P0's base side (run 45612339270, staging e0747fc1) showed `test-prot-norm-core.R` (306),
  `test-workflow-omics-declarations.R` (23) and `test-workflow-omics-owner.R` (1) failing identically on base.
  `test-integration-gui.R` passes on both sides of the train's run and failed once in P0. **The moving expectations,
  by file:**

  | File | Key | Movement | Approval (VAL-241) |
  |---|---|---|---|
  | `test-metab-00-import-detection-characterization.R` | `defaults$is_pattern` | the decided string | (i) |
  | `test-metab-import-detection-direct-shared.R` | `msdial_defaults$is_pattern` | the decided string | (i) |
  | `test-workflow-omics-import.R` | `omicsImportFrozen_i2_lipid` msdial, no_annotation, custom, malformed; `omicsImportFrozen_i2_metab` msdial, custom | the six digests the ticket named; lipidsearch and the trx previews did not move | (ii) |
  | `test-omic-norm-module.R` | `omicNormFrozen` lipid/metab × manual/automatic: `run$log` (+2 lines, one per assay), `run$notes` (+1 notification), `outputs` of run, skip_correlation, reset (manual) and run, final_plot, final_plot_failing (automatic) | 20 expectations; `states`, `provenance`, `normdata` and every other key unchanged | (iii), captures under OD-2 (b) |
  | `test-omics-norm-core.R` | the error run's `notes[[1L]]` | the ITSD notification now comes first; the error note is `notes[[2L]]` | (iii), captures under OD-2 (b) |
  | `test-integration-gui.R` | KIT-139 facts | one fact added after the standards count: "lipidomics_LCMS_Pos: internal standards marked: itsd1 (ITSD1)." | (iii) as read by the implementer: a capture P0 lists, from OD-4 (a)'s naming (the ticket's 8.1 did not foresee it; flagged for the audit) |

  The module-CI norm captures, the GUI-entry lanes and the record captures the ticket thought might move did not move.
  The KIT-073 and KIT-094 fixtures were not exercised by a moving expectation.
- **R1** (run 45618785756, base = tests-only 63e29c04, head = fix 288abf49, the new file): base 33 passed, 81 failed,
  1 error; head 117 passed. On base each E row fails at the ticket's base value: E1 TRUE TRUE; E2 TRUE TRUE (also on the
  custom preview); E3/E4 counts 13 per assay with the 9 returned rows absent and the 5 decided standards present; E5 13
  marked, then the OD-4 naming function absent (the one error); E6 no log line, no notification, no record outcome
  (both routes); E7 the help text absent; E8 (ii) no `itsd_note`, no `* Note:` line, no report note (E8 (i), 13 and 4
  rows, passes on base); E9 the old string. G2, G3 and G1's four intended labels pass on base; G1's `(d7)` element and
  G4 fail on base, as the ticket's table states for them. **Stop rule 1 not triggered.**
- **R2** (run 45621200497, base = tests-only, head = fix, 130 files, 7 shards per side, R CMD check both sides): 123
  `same`, 7 `fixed` (the new file and exactly the six files of the table above); no `new-failure`, no
  `message-changed`, no `count-changed`. Check: ERRORs 0, WARNINGs 3, NOTEs 3 on both sides, findings identical (6),
  none only on head. **Stop rule 2 not triggered.**

R2's file list (recounted at e0747fc1): the ticket's grep, 106 files of 491; its four named files (none in the grep); 19 further lipidomics,
metabolomics and integration files that exercise the other code r1.3's file list touches (the study-parameters writer,
the reports, the integration facts and model record, the tab pipelines: e.g. `test-integration-gui.R`,
`test-integration-record.R`, `test-kit235-lipid-metab-report-parser.R`, `test-omic-norm-module.R`,
`test-workflow-analysis-record.R`, the two GUI-entry lanes and e2e files); and the new file: 106 + 4 + 19 + 1 = 130 (24 beyond the grep). Proteomics- and
transcriptomics-only matches of the wider grep were left out.

## Acceptance criteria

| # | Criterion | Evidence | State |
|---|---|---|---|
| A1 | G0 observes D0, D0b, D1 to D7 as predicted; C14's record cited | `g0-outcome.md`, run 45592957172 | met (D6 recorded; its 'none' route had no DA table) |
| A2 | G-R recorded if run | not run (optional, by instruction) | not applicable |
| A3 | Tests-only commit: every E row fails at G0's base value; G rows pass; all pass on the fix | R1 45618785756; R2 45621200497 | met, with G1's `(d7)` element and G4 failing on base as their own table rows state |
| A4 | P0 lists the moving expectations; R2 as 8.2; approved files differ only in their listed expectations | P0 45603194585 and 45612339270; R2 | met; `test-integration-gui.R` is a capture the ticket's 8.1 did not list (above) |
| A5 | R3 (the merge train's full suite) | the train's | open (the coordinator's) |
| A6 | NEWS, ledger row, `outcome.md` with "For the change report"; `man/`, `NAMESPACE`, `DESCRIPTION` unchanged | `NEWS.md` (Scientific changes), `docs/methods-decision-ledger.md` row `VAL-<merge>`, this file; `git diff e0747fc1 -- man NAMESPACE DESCRIPTION` empty | met |

## Budget

Job-minutes: G0 18.0, P0 236.2, P0 base side 54.0, R1 22.3, R2 447.5, blank request 0.2: 778.2 = **51.9 core-hours**
(expected 56.67 under OD-2 (b) plus the train's share; ceiling 201.33). No correction pass, no instrument rerun.

## Observed and inferred

Observed: every run result above; the diff. Inferred: that `test-integration-gui.R`'s movement falls under the owner's
"(iii) the captures P0 lists"; that the record viewer shows `itsd_features_per_assay` and `itsd_note` (it prints every
step outcome, `inst/reports/shared/analysis_record_SI.rmd:55-57`; not rendered here).

## For KIT-242's implementer

- KIT-242 adds a per-assay ITSD count line to the tab and `itsd_counts_per_assay` to the record. KIT-260 already logs,
  per assay, "ITSD, assay <a>: <n> rows used as internal standards, then removed: <rows>" and records
  `itsd_features_per_assay`, from `omicsItsdReportResult` and `omicsRouteRecordOutcomes` in
  `R/utils_workflow_omics_declarations.R` and the two calls in `R/utils_workflow_omics_actions.R` (managed and in
  process, right after the ITSD step). Merge the count into that line or keep both; do not log the count twice.
- The ITSD kernel now returns `itsd_used`, `itsd_id_column` and, when present, `itsd_legacy_rows` beside `final`. A
  skipped assay has no entry there; KIT-242's "matched nothing" notice and record belong beside them.
- `test-omic-norm-module.R`'s `omicNormFrozen` and `test-omics-norm-core.R`'s error-run notes now carry KIT-260's lines;
  KIT-242's own log or notification changes will move the same keys again.
- **One name, two meanings (I04).** Both branches define `omicsItsdReportResult`: KIT-260 in
  `R/utils_workflow_omics_declarations.R` (`(value, addLog, notify)`, reading `value$itsd_used`), KIT-242 in
  `R/utils_workflow_omics_actions.R` (`(value, addLog, warn)`, reading `value$itsd`). The definitions are in different
  files and hunks, so a merge keeps both without a conflict, and the Collate order (declarations before actions) lets
  the later file's definition win silently when the namespace loads: KIT-260's call sites would run KIT-242's function.
  Whichever ticket merges second renames one of them and keeps a single per-assay log line (the coordinator's ruling,
  2026-11-26; nothing renamed now).
- **The empty pattern box (I09).** With 'Apply ITSD normalization' ticked (the default) an empty pattern box becomes NA
  in the object's slot and pattern-mode ITSD stops with R's "missing value where TRUE/FALSE needed"
  (`R/func_metab_s4_transformation_methods.R:340`, lipid `R/func_lipid_s4_normalization_methods.R:111`; KIT-260's ticket
  2.1 (4)). KIT-260's new import help ends "Leave it empty for none." (wording for the owner's pass, `wording.md`); the
  error is KIT-242's.
- **Merge order (I03).** The coordinator's ruling of 2026-11-26: KIT-260 may merge before KIT-242 if it is ready first;
  it waits for the owner's by-file answer on `test-integration-gui.R` (parked row 42) in either order.

## The correction pass (AUDIT-KIT-260-impl-20261126; one pass, tests first)

- Commits: 7af2d51e tests (I02 (a) per consumer, metab and lipid: record, study parameters, report; I05; I10: lipid
  study parameters and report, the integration-tab note with flagged rows, the model record's marked rows through a
  stand-in trainer); eca09bc5 the I02 (a) fix (`omicsItsdRowsFromStates`); b71b2d62 a coverage test of E6's managed arm
  through the worker's own delivery (`omicsTabServer` with `omicsManagedChainMocks`); then these records; then the merge
  of staging 42846c0e.
- **I05 (the ruling's test).** An older-build product is made by running ITSD with `omicsItsdRowsUsed` returning
  nothing; the record is written as that build wrote it; the replay runs this version's kernel (the ITSD product gains
  the rows-used fields). The test requires no mismatch from the replay harness (`kitReplayJourney`, which judges every
  state and every product by sha256) and identical assays in every replayed state. Its result is in the correction
  pass's runs below; if the harness reports the ITSD product, the coordinator rules (a replay marker, as KIT-242's
  OD-4 (c)); none is added here.
- **I06, I07, I08, I09, I11 to I13** applied as the register's texts (ledger row in five cells with escaped pipes and the
  open by-file question; the R2 file count; `wording.md`).
- **I14** (the DA join failure with a numeric identifier after between-sample 'None'): recorded in `g0-outcome.md`;
  lipidomics likely alike, legacy likely; candidate row K260-G0-D6 is the coordinator's. **I15** the old default string
  survives in `tools/profiling/kit_omics_corpus.R:10` (a comment; that folder is not edited by rule) and in
  `tools/refactor/staging/wave1_*` (historic copies, built out). **I16** the exported `integrationMatchSamples()` returns
  the new names line in its facts and two more elements in its `preview` (additive); the S4 methods list an assay's
  standards before its factor is computed, so "then removed" can be false where the factor fails (moot once KIT-242's
  Part A stops on such a failure).

## Not checked

- The real-data share of MS1-only annotations (G-R not run).
- The report rendered as a document (E8 checks the parser and the Methods paragraph's text, not the rendered file).
- The record viewer rendered (`renderWorkflowAnalysisRecord`).
- The integration model with MOFA2 itself (the correction pass's test trains with a stand-in trainer and reads
  `internal_standards_marked` back from the product and the record outcomes).
- The application in a browser.
- R3 (the merge train's full suite).

## The correction pass: results (runs 45635735448 and 45635740954)

- **R1b** (the new file; base 7af2d51e, head e83c8011): on base every I02 test fails at 'no note' (metab and lipid; the
  record's `itsd_note` and `itsd_features_per_assay`, the study-parameters `* Note:` line, the report's note and
  Methods paragraph); on head they pass, the twins carry no note. The I10 tests (lipid study parameters and report, the
  integration-tab note, the model record's marked rows through the stand-in trainer) and E6's managed arm through the
  worker's delivery pass on head. Head: 150 passed, 1 failed.
- **I05: the stop condition of the coordinator's ruling.** On head the replay of an older-build record (its ITSD product
  without the rows-used fields) reports six mismatches from the replay harness (`kitReplayJourney`, every product
  judged by sha256): `metab_itsd_norm/products/metab_itsd_normalisation_science`, and on `metab_log2` and
  `metab_normalized` the ITSD product they carry and their own products, which chain to it
  (`metab_log2_science`, `metab_between_sample_normalisation_science`). Every replayed state's assays are identical to
  the origin's (the test's value expectations pass). So the replay comparison counts KIT-260's additive product fields
  as a difference. By the ruling the work stops here: a replay marker (as KIT-242's OD-4 (c)) or another standard for
  additive display fields is the coordinator's decision; none is added.
- **R2b** (130 files, base staging 42846c0e, head e83c8011, check both sides): 128 `same`; `test-omics-norm-core.R`
  `count-changed` (one more passing expectation, predicted); the new file `only-on-head` with the one I05 failure; no
  `new-failure`, no `message-changed`; check findings identical (6 each side), nothing only on head. Predictions:
  `predictions-correction.md` (all held; I05's mismatch was the likely outcome stated there).
- Budget: the correction pass used 23.6 + 467.3 = 490.9 job-minutes = 32.7 core-hours (R2b's ceiling 73.33); the
  ticket's total is now 778.2 + 490.9 = 1,269.1 job-minutes = 84.6 core-hours (ceiling 201.33).

## The second correction (VAL-276; the owner's decisions VAL-275 and VAL-276, staging abfc99ac)

Why: the first correction's I05 test showed that the replay harness counts KIT-260's additive ITSD product fields as a
difference when an older record is replayed (run 45635735448). VAL-276 allowed one narrow correction: those fields are
written only under a replay marker, the one KIT-242 uses, so one marker covers both tickets. VAL-275 allowed the one fact
line in `test-integration-gui.R` and confirmed the 20 captures of `test-omic-norm-module.R` and the 1 key of
`test-omics-norm-core.R`. No correction follows this one.

- a649d90b, tests only: I05 now builds its older record as a build before KIT-260 made it (the ITSD step under the
  marker's earlier token, so no marker recorded and none of the fields); its mirror: a record made at this version
  replays exactly and its replayed ITSD product carries the rows used. I05 was the failing-first evidence (it fails on
  the head before the fix, run 45635735448).
- 933671b5, the merge of origin/v0.9-staging abfc99ac (fetched first): the ledger's append collision kept both,
  staging's first; NEWS and every file under `R/` and `tests/` merged without conflict. The merge was made before the
  fix so that the marker row follows staging's rows in `workflowReplayRuleMarkers` without a conflict.
- 35496462, the fix: `workflowReplayRuleMarkers$itsd_outcome` (`option = "MultiScholaR.itsd_outcome"`,
  `current = "per_assay"`, `earlier = "applied_flag"`, the same name, option and tokens as KIT-242's row on
  `origin/kit-242` 2776b76c); the ITSD declaration fixes `itsd_outcome = "per_assay"` (as KIT-242's does); the kernel
  writes `itsd_used`, `itsd_id_column` and `itsd_legacy_rows` only when the option does not hold the earlier token. The
  record's `itsd_features_per_assay` and `itsd_note` are record outcomes, not part of any compared product, and stay as
  they are (for an older record they are read from the stored states, I02 (a)).
- KIT-242's semantics carry these fields: under its earlier token the ITSD step writes "no per-assay fields", which is
  exactly the product of a build before KIT-260. When KIT-242 merges, both rows and both `fixed` entries are identical
  (KIT-242 must keep one row and one entry).
- **Stored references this correction may move beyond VAL-275 (predicted, not approved):** the ITSD step now records a
  fixed constant, so the `provenance` digests of `omicNormFrozen` in `test-omic-norm-module.R` move wherever ITSD runs
  (the lipid and metabolomics manual and automatic captures), and any capture that pins the ITSD step's recorded
  parameters or step provenance may move (KIT-224's marker moved 12 such digests, re-pinned by the owner's word by file).
  The comparison run lists them; none is re-pinned without the owner's word.
- **Merge order:** KIT-260 merges with or after KIT-242 (the coordinator's instruction of 2026-11-26 with VAL-276).

### Result of the comparison run (run 45665921270; stopped)

- `test-kit260-itsd-default-pattern.R` (only on head): 156 passed, 0 failed. **I05 now passes** (the older record replays
  with no mismatch and identical assays) and its mirror passes (a new record replays exactly; its replayed ITSD product
  carries the rows used).
- 122 files `same`; `test-omics-norm-core.R` `count-changed` (predicted); check findings identical (6 each side).
- **Expectations outside VAL-275 that move (stop; none re-pinned):**

  | File | Test | Expectation | Why |
  |---|---|---|---|
  | `test-omic-norm-module.R` | (b) lipid and metab, manual and automatic | 12 `provenance` digests of `omicNormFrozen` (manual: run, skip_correlation, reset; automatic: run, final_plot, final_plot_failing; both omics) | the ITSD step records the fixed marker `itsd_outcome` (predicted) |
  | `test-workflow-lipid-declarations.R` | KIT-031: the lipidomics declarations are registered with the units, merges and statistics of the step map | `declarations$lipid_itsd_normalisation$fixed` identical to `list(method = "ITSD")` | the declaration fixes the marker |
  | `test-workflow-metab-declarations.R` | KIT-032: the metabolomics declarations … | `declarations$metab_itsd_normalisation$fixed` identical to `list(method = "ITSD")` | the same |
  | `test-integration-gui.R` | KIT-139: internal standards, features seen once and few samples are stated with their counts | `expect_match(facts[[6L]], "proteomics: 1 features measured in fewer than 3 of the 5 biological samples")`, added by KIT-240 on staging | the fact line VAL-275 allowed comes before it, so that fact is now `facts[[7L]]`; a merge interplay, no textual conflict |

- `test-workflow-omics-binding-steps.R` (KIT-094) and `test-workflow-omics-declarations.R` (KIT-073) fail on both sides
  (laptop fixtures) and their messages change: the ITSD and later steps record one more parameter (the marker) and the
  inputs' digests differ. These are the fixtures VAL-241 (iii) defers to the end-of-project re-record (VAL-181, VAL-189).
- When KIT-242 merges, its own marker makes the same `provenance`, `fixed` and KIT-073/KIT-094 movements; one by-file
  answer can cover both tickets.
- Budget: 484.4 job-minutes (32.3 core-hours, within this run's 73.33); the ticket's total is 1,753.5 job-minutes =
  116.9 core-hours (ceiling 201.33).

### The coordinator's rulings on run 45665921270 (2026-11-26), and where things stand

1. **`test-integration-gui.R`, KIT-240's `facts[[6L]]`**: test code reading a position, a declared test edit
   (the coordinator's). 513ae741 (tests only) makes the KIT-139 test find the fact by its text:
   `expect_match(facts, "proteomics: 1 features measured in fewer than 3 of the 5 biological samples", fixed = TRUE,
   all = FALSE)`, so a later added fact cannot move it.
2. **Open with the owner (row 50), put once together with KIT-242's same movements; not edited here.** Stored references
   and registration pins moved by the replay marker `itsd_outcome` (old values on staging abfc99ac, new values from
   `r2c-45665921270/failures.txt`):

   | File | Key | Old (staging) | New (kit-260 be377464) |
   |---|---|---|---|
   | `tests/testthat/test-omic-norm-module.R` | `omicNormFrozen` lipid manual run `provenance` | 1ab8a05076a605820dc0d14e237215f9abfaa251d04c7076ca03825877e57dc2 | 5fb69439ca6f755c340b20d38d5ffec9c5a2061349a7f8c5f198afe98924079f |
   | | lipid manual skip_correlation `provenance` | 0b68624dba919ae6e6852528493ba338d672d75d235b8437c2895e8c6f724bf5 | 6c1bc031bc84407f9885c22dd5d36832fbbdda02c0a96bfb9e1ed6e0a4415f0b |
   | | lipid manual reset `provenance` | 6f03ec2fff1c90a3a1728422205f6ebe8d929a0e97a02b42dc0e053b026e5942 | d675df2b63edd2a834774dd9f4c2261228d6aa35460ae436cda34c966da40c28 |
   | | metab manual run `provenance` | 24204497ce33beb525005c4869fd60e6f560a7c1f6eb2293cc928049e1410ba2 | 9424d6af9de6e2ec5323a33ccb1a2f5fd96390ad2597f0734d6d368dbd7aa300 |
   | | metab manual skip_correlation `provenance` | 667cc9562bd90c598ef04d93d5fe30d6ad1818b89c02f50d86dcfec344a11371 | a586d6828e1ba97bcd7b3f4fe9161c1480589c1e1b9466d96e843b44950bfc87 |
   | | metab manual reset `provenance` | 279b01b63aa09b94e9ccaa8e085154e7b281eb0fd48361548782ceead23456e3 | 29b24751a67690c1827837ad097d4f5fd5435fa7a9a6e02c57362f505afc2c45 |
   | | lipid automatic run, final_plot, final_plot_failing `provenance` (3 keys, one value) | 024a3e8b9582da5526ed27f5f48f518ccda26236c0eeb37c4ddd94ad294e955a | 846c6b2c3242fdcbf84e906b76f7026a87afb5e23dad2c0f4379c586c43445cd |
   | | metab automatic run, final_plot, final_plot_failing `provenance` (3 keys, one value) | 6858408bb39aa5c6a558ccd712b51e5ad5ce9aca6ddf7abd28832b96d50a8848 | 96cbc0f23bdb70283b02662cab7d48cb60b4156c44dde8b167700d06a426e64f |
   | `tests/testthat/test-workflow-lipid-declarations.R:47` | `declarations$lipid_itsd_normalisation$fixed` | `list(method = "ITSD")` | `list(method = "ITSD", itsd_outcome = "per_assay")` |
   | `tests/testthat/test-workflow-metab-declarations.R:47` | `declarations$metab_itsd_normalisation$fixed` | `list(method = "ITSD")` | `list(method = "ITSD", itsd_outcome = "per_assay")` |

3. **The KIT-073 and KIT-094 laptop fixtures** (`test-workflow-omics-declarations.R`, `test-workflow-omics-binding-steps.R`)
   fail on both sides already and their messages change (one more recorded parameter per ITSD and later step); they wait
   for the end-of-project re-record, as VAL-241 (iii) says (VAL-181, VAL-189).
4. **No further comparison run now. The next run of this branch is the merge train's.** Its predictions for these
   files, under either owner answer on row 50:

   | File | Owner allows (re-pinned with the new values before the train) | Left for the re-record (not re-pinned) |
   |---|---|---|
   | `test-omic-norm-module.R` | `same` (head passes its new values; base staging passes its own) | 12 `provenance` expectations fail on head only: `new-failure` |
   | `test-workflow-lipid-declarations.R`, `test-workflow-metab-declarations.R` | `same` | 1 expectation each fails on head only: `new-failure` |
   | `test-integration-gui.R` | `same` (the KIT-139 test passes: VAL-275's line and the fact found by text) | `same` |
   | `test-workflow-omics-declarations.R`, `test-workflow-omics-binding-steps.R` | `message-changed` (fail on both; deferred, VAL-181/189) | `message-changed` (the same) |
   | `test-omics-norm-core.R` | `count-changed` (one more passing expectation) | `count-changed` |
   | `test-kit260-itsd-default-pattern.R` | `only-on-head`, all pass | `only-on-head`, all pass |
   | every other file of the 130 | `same` | `same` |

   R CMD check: findings identical on both sides.
5. **Merge order:** with or after KIT-242. One `itsd_outcome` row in `workflowReplayRuleMarkers` and one `fixed` entry on
   each ITSD declaration (KIT-242 declares the same); the ticket merging second renames one of the two
   `omicsItsdReportResult` definitions and keeps one per-assay log line.

### Closure check (AUDIT-KIT-260-impl-20261126, "Closure check (2026-11-26)"): what was done, and the hand-off

- **7c900dfe**, the merge of origin/v0.9-staging bbbfdc1a (fetched first): the ledger's append collision kept both,
  staging's rows first (VAL-280 to VAL-297); everything under `R/`, `inst/` and `tests/` merged without conflict.
- **c545f24e**, tests only (declared test edits, the coordinator's; no expected value changes):
  - CC01: both chunk lists of `test-kit260-itsd-default-pattern.R` load "Hierarchical Parser", "Sample Analysis
    Functions" and "Method Section": since KIT-236 the reports' `FormatMethodDetails` calls `ReportKnown()`, defined in
    "Sample Analysis Functions" (as KIT-235's test loads it).
  - CC03: I05 asserts directly that the older-form origin's ITSD product holds only `final` (the gate itself).
  - CC07: the comment at `test-integration-gui.R:447` cites VAL-275; the ledger row's authority cell names VAL-241,
    VAL-275 and VAL-276.
- **No run was made** for these edits; the next run of this branch is the merge train's.
- **CC03, recorded.** No tests-first run was made for the second correction, although VAL-276 named one: its
  failing-first evidence is I05 as it stood before the fix (it failed on head in run 45635735448) and the comparison run
  45665921270 in which I05 and its mirror pass. I05's replay comparison cannot by itself tell a gated kernel from an
  ungated one, because its origin and its replay come from one build: both sides would carry the same fields. What
  states the gate is I05's record-outcome assertion (with the state fallback stubbed), the mirror's marker assertion,
  and now the direct assertion that the older-form product holds only `final`.
- **CC06, a stated assumption.** `omicsItsdRowsFromStates` reads the annotations from `record$previous_state`, taken as
  the state the ITSD step read (the head, on the application's routes). Were it another state, the same identifiers'
  annotations would be read from it; the annotation column does not change across these steps (INFERRED).
- **CC02, the order (the coordinator's ruling).** KIT-260 and KIT-242 land in ONE merge train, KIT-242 first within it,
  so that no staging commit holds one ticket without the other and every record replays exactly. KIT-260 must not land
  before KIT-242 (the recorded token "per_assay" would change meaning when KIT-242 lands), and not in a later train
  (records made in between would replay with product mismatches once KIT-260 lands). The branch rests until KIT-242
  Part B exists.
- **CC05, checks on the merged tree** (each must give 1):
  - `git grep -c "itsd_outcome = list(" R/utils_workflow_analysis_record.R` (exactly one marker row);
  - `git grep -c "omicsItsdReportResult <- function" R/` summed over files (one reporting function, after the rename of
    I04), and one per-assay log line in the tab;
  - one `fixed = list(method = "ITSD", itsd_outcome = "per_assay")` on the ITSD declaration.
  Row 50's option (b) (the 12 `provenance` digests left failing until the re-record) would hide any later provenance move
  in those 12 captures until then.
- **CC04, the train's predictions.** The predictions above (ruling 4 of run 45665921270) are for KIT-260 alone on
  staging abfc99ac; they must be rewritten for the combined head (KIT-190, then KIT-242, then KIT-260). In that head
  KIT-242's count line moves again the `log` and `outputs` keys of `omicNormFrozen` that VAL-275 approved at KIT-260's
  values, the two `fixed` pins and the 12 `provenance` keys take the combined tree's values (KIT-190 re-dumps the same
  12 lines), and `test-omics-norm-core.R`'s notes may move with KIT-242's notices. The by-file question for the combined
  values is the coordinator's to put once, with row 50. (`predictions-second.md` had predicted `test-integration-gui.R`
  `same`; it failed on KIT-240's position, handled by 513ae741.)

## KIT-260 on KIT-242 (2026-11-26): one merge, the comparison against origin/kit-242 (run 45750259637)

- **State before:** `kit-260` rested at 014fd530 (origin equal); no run was out; every run was recorded (the last, R2c
  45665921270).
- **666ace2e: one merge of origin/kit-242 8162ea52** (KIT-242 Part B on staging aa119308). Conflicts: NEWS and the ledger
  (appends, kit-242's first); `R/utils_workflow_analysis_record.R` (the marker list: kit-242's rows, exactly one
  `itsd_outcome` row, ending `..., correlation_exclusion, normalisation_fallback, protein_normalisation_fallback,
  itsd_outcome, ruv_fallback`); `R/utils_workflow_omics_declarations.R` (one `fixed` entry; the kernel: KIT-242's value
  with KIT-260's rows used added to `value$itsd` as `rows_used`, `id_column`, `legacy_rows`, under KIT-242's own
  condition, which holds only under the marker's current token; the record outcomes: KIT-242's per-assay fields, then
  KIT-260's rows and note read from `itsd$rows_used`, else from the stored states, then KIT-242's between-sample block);
  `R/utils_workflow_omics_actions.R` (the in-process call site: KIT-242's `itsd_run`). Nothing else under `R/`
  conflicted.
- **The agreed interface:** one `omicsItsdReportResult(value, addLog, warn)` (KIT-242's), which now also logs KIT-260's
  per-assay rows line (text unchanged: "ITSD, assay A1: 5 rows used as internal standards, then removed: ..."); KIT-260's
  own reporter is removed; its "rows used" notification (type message) is sent at the two call sites from
  `omicsItsdRowsNotice(value)`. CC05 on the merged tree: `itsd_outcome = list(` 1, `omicsItsdReportResult <- function` 1.
- **No test or expectation was changed in the merge.** Failing-first by value stands from R1 45618785756 and R1b
  45635735448.
- **Result (observed), against `predictions-on-kit242.md`:** 127 `same`; `test-omics-norm-core.R` `count-changed` and
  `test-omic-norm-module.R` `message-changed` (583/56 on both sides), as predicted; `test-prot-norm-core.R` 1547/307,
  `test-gui-entry-lipidomics.R` 19/16, `test-gui-entry-metabolomics.R` 57/48 on both sides, as staging's known counts;
  the KIT-073/KIT-094 fixtures, the pins, the declarations' `fixed` pins and `test-integration-gui.R` `same`; R CMD check
  findings identical (6 each), nothing only on head.
- **Not as predicted:** the new file is `only-on-head` with 156 passed and **1 failed**: the mirror test's last
  assertion reads `product$itsd_used$A1$annotation`; under the agreed interface the rows used travel in
  `product$itsd$rows_used`, so the path is NULL. The value is carried (the record's `itsd_features_per_assay` and the
  tab's log, which read `itsd$rows_used`, pass in the same file). The test predates the interface; by the
  coordinator's rule it is handed back, not edited. The one-line repair, for the coordinator's ruling:
  `expect_identical(product$itsd$rows_used$A1$annotation, KIT260_STANDARDS_HEAD)`.
- **Stored values that moved, by file and line (not edited; the ITSD marker's digests are covered by VAL-303 and
  VAL-352):** `tests/testthat/test-omic-norm-module.R`, `omicNormFrozen`, keys failing on both sides whose message
  changes because the merged tab adds KIT-260's per-assay rows lines and notification to KIT-242's lines: lipid manual
  run `log` (:1016) and `outputs` of run (:1050), skip_correlation (:1102), reset (:1139); metab manual run `log`
  (:1158) and `outputs` of run (:1192), skip_correlation (:1244), reset (:1281); lipid automatic run `log` (:1300) and
  `outputs` of run (:1333), final_plot (:1371), final_plot_failing (:1388); metab automatic run `log` (:1430) and
  `outputs` of run (:1463), final_plot (:1501), final_plot_failing (:1518). These are keys VAL-275 approved at KIT-260's
  stand-alone values; their combined values (KIT-242 and KIT-260) are a by-file question for the coordinator to put once
  with row 50 (CC04).
- **Budget:** 498.8 job-minutes; running total 2,252.3 of 3,020 (150.2 of 201.33 core-hours). A further comparison of
  this size (about 500) fits the remainder (767.7); only one correction pass remains.

## The ticket's one correction pass (2026-11-26), and the 16 keys for the owner

**The correction pass, spent.** The coordinator's ruling, quoted: "the one failing assertion
(`product$itsd_used$A1$annotation`) asserts a product path that the agreed interface moved (your fields live in
`value$itsd`, by my instruction); the expected value is unchanged and the same value is asserted through the record
outcome and the tab log in the same file. The repair you propose is allowed: `expect_identical(product$itsd$rows_used$A1$annotation,
KIT260_STANDARDS_HEAD)`, nothing else. It is recorded as the ticket's ONE correction pass (a declared edit of an
expectation's path after a result was seen), not as an instrument repair." Applied in f59a2af7 (tests only, that line
alone). No correction pass remains.

**The 16 keys of `tests/testthat/test-omic-norm-module.R` (`omicNormFrozen`)** whose stored values VAL-275 approved at
KIT-260's stand-alone values and which on the combined head (KIT-242 + KIT-260) carry both tickets' lines. "Today" is
the stored value on origin/kit-242 8162ea52 (staging's, not yet re-pinned by KIT-242); "KIT-260 re-pin" is commit
63e29c04 (and the merged tree, which keeps it); "combined head" is what kit-260 0c7fbdf7 produced in run 45750259637
(from its `failures.txt`). Not re-pinned further.

| Key (line on the merged tree) | Today (origin/kit-242) | KIT-260 re-pin (63e29c04) | Combined head (run 45750259637) |
|---|---|---|---|
| `lipid manual run log` (:1016) | the 0.9 lines (20 entries; no internal-standard line) | today's plus KIT-260's two per-assay lines before "ITSD normalization complete" (22 entries) | "ITSD normalization complete", then KIT-242's "Internal standards used per assay: LCMS_Pos 1; GCMS 1", then KIT-260's two per-assay lines; and KIT-242's "RUV-III not applied: no assay could be corrected (see the warnings above); the data continue without RUV" in place of "RUV-III correction applied" (23 entries) |
| `lipid manual run outputs` (:1050) | a8a2cb84e449727a31e0dc2112ec3cf9faa4a4fb81bf734ea0b31536ec285c6c | d3cc98f4ee7167dee0764a30300f87b96c7f394c3fb464d70cce6590934892d5 | 6b3257380985053f717deb20ad965f7f196ccf5943582871e36bdd4f5017b8d2 |
| `lipid manual skip_correlation outputs` (:1102) | 701a6ac1d87ea09b693aac28f20a465e43f8c30008b3c50d175c8b49582187f8 | e1058a016b6290767718154fea3003b12384442225ca601b158afa38b1a02ffd | 0100859bb3a50c0122d87f9629c23467ee1bf1ef02c582a430e536f6f481fc96 |
| `lipid manual reset outputs` (:1139) | 3a4f085c66ea18873418f52fa26eb53d66b7b88bcad7ecf7f680d0c6cb248b6f | 1c3646ad814d6ae58bf628301f6aeae3b1ee8f3a751ddca396a3682a1114a0d1 | 257b7fcacf230304c4bf38116b9857241d66e94b604d5caf770b7f725d41306f |
| `metab manual run log` (:1158) | the 0.9 lines (20 entries; no internal-standard line) | today's plus KIT-260's two per-assay lines before "ITSD normalization complete" (22 entries) | "ITSD normalization complete", then KIT-242's "Internal standards used per assay: LCMS_Pos 1; GCMS 1", then KIT-260's two per-assay lines; and KIT-242's "RUV-III not applied: no assay could be corrected (see the warnings above); the data continue without RUV" in place of "RUV-III correction applied" (23 entries) |
| `metab manual run outputs` (:1192) | a8a2cb84e449727a31e0dc2112ec3cf9faa4a4fb81bf734ea0b31536ec285c6c | 45624e015319acf091d92ea85574964e93ebe9e350ab00b2e473cfd4bf939690 | eb58e3d9a7b494bd7ee959f5b6c5caf9fb4b3207b4e23531174a80042037b63e |
| `metab manual skip_correlation outputs` (:1244) | 701a6ac1d87ea09b693aac28f20a465e43f8c30008b3c50d175c8b49582187f8 | 9d473dbca0ede8af25475634752283310652012facc2cf169780a2303873f12b | 433336fdb161af9ba079ecc28bc249a9e34ce9524511425608a0f5abe9e38796 |
| `metab manual reset outputs` (:1281) | 3a4f085c66ea18873418f52fa26eb53d66b7b88bcad7ecf7f680d0c6cb248b6f | 962bce5d0ee602e1b52a4204061fa1396a7922e70baf4d6aa93d5df95d6fe7e7 | be2771ee2999e73d3ee7214bfe81e6425604f4486e997d60d6f683d5e4da4f4e |
| `lipid automatic run log` (:1300) | the 0.9 lines (18 entries; no internal-standard line) | today's plus KIT-260's two per-assay lines before "ITSD normalization complete" (20 entries) | "ITSD normalization complete", then KIT-242's "Internal standards used per assay: LCMS_Pos 2; GCMS 2", then KIT-260's two per-assay lines; and KIT-242's "RUV-III not applied: no assay could be corrected (see the warnings above); the data continue without RUV" in place of "RUV-III correction applied" (21 entries) |
| `lipid automatic run outputs` (:1333) | 2ae8151f2d2a90bababe866ccf367b11b51981e4f62d6fd807352bb70c0649b9 | b024bcd03404063dad170d139517809a5d00a445de94f5ca3949225003e3a1f3 | be600cf55f440c4ae3bbde44503cf741275fe7a469e3376e6b303493e7007477 |
| `lipid automatic final_plot outputs` (:1371) | 2ae8151f2d2a90bababe866ccf367b11b51981e4f62d6fd807352bb70c0649b9 | b024bcd03404063dad170d139517809a5d00a445de94f5ca3949225003e3a1f3 | be600cf55f440c4ae3bbde44503cf741275fe7a469e3376e6b303493e7007477 |
| `lipid automatic final_plot_failing outputs` (:1388) | 2ae8151f2d2a90bababe866ccf367b11b51981e4f62d6fd807352bb70c0649b9 | b024bcd03404063dad170d139517809a5d00a445de94f5ca3949225003e3a1f3 | be600cf55f440c4ae3bbde44503cf741275fe7a469e3376e6b303493e7007477 |
| `metab automatic run log` (:1430) | the 0.9 lines (18 entries; no internal-standard line) | today's plus KIT-260's two per-assay lines before "ITSD normalization complete" (20 entries) | "ITSD normalization complete", then KIT-242's "Internal standards used per assay: LCMS_Pos 2; GCMS 2", then KIT-260's two per-assay lines; and KIT-242's "RUV-III not applied: no assay could be corrected (see the warnings above); the data continue without RUV" in place of "RUV-III correction applied" (21 entries) |
| `metab automatic run outputs` (:1463) | 2ae8151f2d2a90bababe866ccf367b11b51981e4f62d6fd807352bb70c0649b9 | 3704d8a5534af4a8d928fccf30ee451efa65ba801e30321ead9441c572df5018 | 295c409a7898fb8180ab1f22ae0ed6493df785b0421a9a68f5057984970ce642 |
| `metab automatic final_plot outputs` (:1501) | 2ae8151f2d2a90bababe866ccf367b11b51981e4f62d6fd807352bb70c0649b9 | 3704d8a5534af4a8d928fccf30ee451efa65ba801e30321ead9441c572df5018 | 295c409a7898fb8180ab1f22ae0ed6493df785b0421a9a68f5057984970ce642 |
| `metab automatic final_plot_failing outputs` (:1518) | 2ae8151f2d2a90bababe866ccf367b11b51981e4f62d6fd807352bb70c0649b9 | 3704d8a5534af4a8d928fccf30ee451efa65ba801e30321ead9441c572df5018 | 295c409a7898fb8180ab1f22ae0ed6493df785b0421a9a68f5057984970ce642 |

The combined head's `log` and `outputs` values differ from KIT-260's re-pin by KIT-242's lines (its count line and its
RUV outcome line), so the owner's by-file question for these 16 keys is the combined values above, put once with row 50.

## For the change report

**What was wrong, and who was affected.** On `main` and in 0.9.0, every metabolomics table and every lipidomics table
except LipidSearch was read with the MS-DIAL parser, which pre-filled the import's internal-standard pattern with
`^IS_|_d[0-9]+$|ISTD|IS-|w/o MS2:`. MS-DIAL 4 and 5.0 to 5.1 prefix 'w/o MS2:' to every feature annotated from MS1
without an accepted MS/MS match (ordinary, MSI level 3 features); `IS-`, unanchored and case-sensitive, also matched
names written in capitals such as 'CIS-ACONITIC ACID'. With ITSD ticked (the default) and no standard picked by hand,
these features set each sample's normalisation factor and were then removed from the data, the differential analysis,
the reports and (new in 0.9, even with ITSD unticked) the integration; the study-parameters file and the report listed
them as the standards used (as integer `Alignment ID`s for MS-DIAL). Affected: lipidomics and metabolomics users who
imported any file except a LipidSearch one with the proposed pattern and ran ITSD without picking standards, or ran the
integration. On the default route (log2 then cyclic loess) with complete data the factor is removed by the next step,
so the lasting effect was the removal (and its knock-on on RUV-III-C's controls and k, the trend prior and the FDR
family); with 'None' or 'Scale' the factor also shifted fold changes; with missing values a small residual remains
under cyclic loess.

**Code.** Base e0747fc1: `R/func_metab_import_detection.R:207`, `R/func_lipid_import_detection.R:234` (the string);
legacy `main:R/func_metab_import.R:393`, `main:R/func_lipid_import.R:421`, pre-filled at `main:R/mod_metab_import.R:507-508`,
`main:R/mod_lipid_import.R:517-518`, since 45160aa5 (2025-12-17) and ce640728 (2026-01-18); tags v0.5.0-rc1, v0.9.0.
Fix 288abf49: the string at `R/func_metab_import_detection.R:210`, `R/func_lipid_import_detection.R:235`; the display
and record changes listed under "What was changed".

**What changed in results.** New imports of an MS-DIAL 4 or 5.0 to 5.1 table no longer use or remove MS1-only
features, nor capitalised names containing 'IS-'; labelled standards written `(d7)` are now used. On the section-4
objects (two assays, 60 rows each) the rows used fall from 13 to 5 per assay and the 9 returned rows are normalised like
every other row (R1 E3/E4). Where no labelled standard remains in an assay, no ITSD is applied to it (KIT-242 shows and
records that). Projects imported before keep their stored pattern and replay unchanged (G0 D7; R1/R2 E8 (i)); where
their pattern matched such features, the record, the study-parameters file, the report and the integration tab say so.

**Evidence.** MS-DIAL's source as cited in the ticket (section 1); G0 run 45592957172 (D0 to D7: the base pattern's 13
matches of 20 names, every import layout's proposal, the declaration and worker routes, the study-parameters listing,
the integration marking, the QC tab, D6's statistics as recorded, D7's storage and replay); P0 45603194585 and
45612339270; R1 45618785756; R2 45621200497. Limits: synthetic objects only; the real-data share of MS1-only
annotations is not measured; D6's 'None' route could not be observed at the base commit (the DA step refused integer
identifiers after method 'none'; recorded in `g0-outcome.md` as an adjacent discovery).
