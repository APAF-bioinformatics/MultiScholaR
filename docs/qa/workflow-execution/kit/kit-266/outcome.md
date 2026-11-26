# KIT-266 outcome (branch kit-266, from origin/v0.9-staging 7b4ae7da)

Work item (AGENTS.md, stated before editing): the DIA-with-limpa route builds Figure 1 (VAL-308; answers VAL-316; ceiling
VAL-333). Ticket KIT-266 r1.5. Files: the ticket's section 12 list. Checks: R1 (tests-only, failing first), R2 (base against
head, 55 files and R CMD check), R3 (render), M1 (controller memory). Exit: sections 10 and 13 of the ticket. One correction
pass; instrument defects per VAL-54 (at most three re-sends across R1 and M1).

## Step 0: the ticket's lines read against the merged tree (7b4ae7da)

Read: AGENTS.md, conventions (4), IMPLEMENTER-BRIEF, KIT-266 r1.5, AUDIT-KIT-266-plan-20261126 (Re-check 2 in full), ledger rows
VAL-308, 316, 333, 334, 342, 343, 346, 348, 349.

G0's reuse condition (ticket section 11): `git diff 894b9cc8 7b4ae7da` over the files G0 observed, by hunk. No observed hunk
changed: the tab entry (`R/mod_prot_norm_workflow_helpers.R` 374-458), the gate and the aesthetics path
(`R/mod_prot_norm_qc_support_helpers.R` 116-155), the pre generator and source rule (`R/mod_prot_norm_qc_generation_helpers.R`
1-160), the pass-through (`R/utils_workflow_dia_norm.R`, now 745-769), the Reset (`R/mod_prot_norm_session_helpers.R`), Copy's
composite entries and the harness (`tests/testthat/helper-gui-entry.R`) are unchanged in the observed hunks. G0's premise
stands; G0 is not repeated.

What moved (lines now, at 7b4ae7da; the ticket's are at 176f9d3f / origin/kit-259b):

| Item | Ticket | Now |
|---|---|---|
| KIT-259 merged (train 23): `protNormFigureTitle(filename, adjusted)` | `origin/kit-259b:338-445` | `R/utils_workflow_dia_norm.R:339`; `writeProtNormFigurePlaceholder` :360; `drawProtNormQcFigure` :396; `writeProtNormCompositePlaceholders(ruvMode, qcDir, omicType)` :424 (no `resolveInputsFn` yet); `withProtNormFigureNotice` :464 |
| resolver / builder / pass-through | `:617`, `:673`, `:702-725` | :617, :673, :745-769 (+43 above them); `diaLimpaPassThroughSource` :777; `removeProtNormOwnFiles` :850 |
| KIT-271 merged (train 26): `resolveProtNormCorrelationTableNames` | `:742-768` | :822; Copy's limpa branch `R/func_general_filemgmt_results_summary_helpers.R:386-400` |
| composite call and list call in the publication | `:1107`, `:1110` (kit-259b) | `R/utils_workflow_dia_products.R:1107-1110` |
| `computeWorkflowDiaQcFiguresAction`, `diaQcFiguresParameters`, twin `saveWorkflowDiaCompositeFigure`, `writeDIACompositeFigureFiles` | `:2090`, `:2068`, `:2177-2221`, `:2232` | :2095, :2073, :2182-2226, :2237 |
| module helpers: covariates `finish` | `:100-101` | :100-101 (the observer's call is now wrapped in `withProtNormFigureNotice`, :107-113); summary text :180; `dispatchProtNormPassThrough` :268-283; `dispatchProtNormRuvCovariates` :293-325 (`publish` wrapped :300-310); `runProtNormNormalizationWorkflow` limpa branch :504-533; Apply :685-704; `dispatchProtNormCorrelation` :733-765; Skip :785; `dispatchProtNormNormalization` `preQc` wrapper :656-658 |
| limpa template | `:586-631`, `:1429-1436`, `:1476`, `:1486` | Result Section :577-640, Main results call :1429-1436, `fig.cap` :1476, `QC_metrics_placeholders` chunk :1486-1517 (`describe()` maps `ruv_corrected_cancor` only :1502-1503) |

Nothing contradicts the ticket's design. Ticket predictions restated on top of KIT-271:
- `prot_dia_limpa.json` moves as the ticket's section 8 says, plus KIT-271's 13 lines (the table is named
  `normalised_results_cln_with_replicates.tsv` / `Publication_tables/normalised_results.tsv`); not edited here.
- `test-gui-entry-proteomics.R` is 92/101 on staging, with the four file-name key lines already failing at the ten steps from
  `correlation` on (VAL-349 (c)). R2's prediction: head failures = base failures plus the limpa steps from `norm` on that pass on
  base (they gain the composite and `pre_norm_*` keys), compared with the four file-name key lines masked, as merge train 26 did.
- `test-prot-norm-core.R` stays 1547/307 on base (stored values already failing); its declared edits add base-only failures
  (section 8 item 6) and no head-only failure.

## Coordinator's notice (KIT-272 on merge train 31; staging ec27f69e)

KIT-266 is the second to merge. To apply when R2's predictions are restated, after KIT-272 has landed on staging (merge staging
once then, and record what conflicted in `DIANN_limpa_report.rmd`; this branch edits `FormatResults`, the Main results call, the
`QC_metrics` caption and the placeholder chunk; KIT-272 edits the density/PCA bullet lines):
1. G3 (`test-kit266-limpa-combined-qc-figure.R`, declared edit, AUDIT-KIT-272-impl I07): base-side substitution of the bullet and the
   PCA bullet, 2 x 3 states per template, each with an exactly-one-match assertion, before the comparison with `kit236/base`.
2. R2 also runs `test-kit272-prot-box-panel-text.R`.
3. R3's V0 reference text (kit-236 render) holds the old bullets: restated from KIT-272's text; KIT-266's text keeps "the same
   scores that panels a-c plot as points" true for the figure it builds.
Budget: one more test file per side in R2 (about one test-minute); the declared figures (expected 48.03, ceiling 163.46
core-hours) are unchanged. No current run is disturbed.

## R1c (run 45745220661, head 449720a7, new file only; an extra declared head check, charged to the ticket) against predictions-r1c.md

183 passed, 2 failed, 0 errors, 34 test-minutes (42.4 job-minutes incl. setup, 2.83 core-hours). Prediction (all rows pass; E12 (c) among the
rows at risk): confirmed for E1-E21 and G1-G6 except E12 (c). Read: E16 (managed pass-through), E17, E18, E19 (fallback logged), E20 (A111-01
order, final figure 1950 px wide, 16 actions), E21 pass: through the in-process double the worker's files are found saved by the controller (the
sibling save, its identity and its root work).
E12 (c) failed on its two later assertions: after a second covariates press with the same inputs the marker and the list remain. Cause by
reading and the log (the marker WAS removed by `removeDiaLimpaRuvMarkers`, line 47 of the log): `publishDIAActionOutputs` replays an earlier
publication with the same identity (product, colour, shape), copying its staged files, markers included; the figure is not drawn again, so
placeholder image, marker and list stay consistent with each other. The test's premise was wrong, not the code. Corrected in the tests-only
base (no expected value changes: marker and list gone after a run that draws the figure): the second press follows a change of the
aesthetics, which changes the identity. This is a test-premise correction before R2 (R2's base side carries it); not a re-send.

## M1's failed path control: diagnosis (before any re-send; coordinator's rule (a))

R1c shows E16 to E21 pass on the head, so the controller finds the worker's files saved when the worker role runs through the real entry point
in process. The product side (sibling save, identity, root, acceptance, index) is therefore not the cause. The one difference in M1 between
the controller and the child is the package: `m1.R:30` (first send) loaded the controller with `pkgload::load_all` (source tree) while the
child ran the installed package; the identities carry the software fingerprint of the store's context. Not proven by reading which field
differed. Repair (script only): the controller now runs the installed package (library + the test helpers sourced into an environment under
its namespace), and a failed path control now prints which of the three files are found saved, the layout string and the software
fingerprint, so a further failure is self-explaining. M1 is re-sent as instrument re-send 2 of 3 (VAL-54), declaration restated in
m1/declaration.md.

## Deviations accepted by the coordinator
(1) The pass-through composite is the last action of the same runner run, not chained by `done` (the KIT-244 runner double has no `done`
argument): same actions as the ticket (three). (2) The wrapper deletes an earlier list of replaced panels before writing the new one.
(3) D1 caption: the ticket's draft named left, centre and right columns; the pass-through figure has two, so the caption now says "the first two
columns are before the sample-correlation filter, and a third column, when the figure has one, is after it" (test row: it names no left, centre
or right column). (4) The cyclic loess centre title has a 32-character first line (ticket limit 30); R3's lane uses method none and does not render it.

## M1 (A5a): the controller's memory at the limpa_composite delivery

First send 45743367152 (path control failed), second 45750228451 (script died in helper sourcing, no measurement), third (valid) 45751962444, bound to 495af273
(source tree in controller and child). Observed: the controller process's rise from the reset to the end of the delivery is 0.007 GiB at the
two-column build and 0.010 GiB at the three-column build (VmRSS before 817,008 and 848,216 kB), against KIT-122's 0.3 GiB applied to the moved event; the
in-process build in a fresh process reads 0.540 and 0.914 GiB (positive control). Controls: reset (VmHWM equals VmRSS after "5" to
clear_refs), path (files written from the store at both points, no fallback logged, PNGs 1500 x 3300 and 1950 x 3900). Cause of the first
failed path control (inferred from the third send): the controller and the child ran different code (a `load_all` tree against the
installed package); with the same source tree in both the worker's files are found saved. Instrument re-sends used: 3 of 3 (R1b, M1b, M1c).

## R2 (run 45750231427): base 5e3de9c0 (staging 36e3180b + the branch's tests) against head 41035e75, 57 files, check on both sides

Verdicts: 50 same; fixed: test-kit266 (base 71 passed / 101 failed / 5 skipped, head 210 / 0), test-kit259b (base 5 failed, head 0), test-kit227 (base 1,
head 0), test-prot-norm-worker-figures (base 3, head 0): each as predicted. test-prot-norm-core.R message-changed (base 1531 passed / 327 failed / 1 error, head 1555 / 307 / 0):
head 307 = the staging stored-value failures; the 20 failures and 1 error that are base-only are all in the two declared blocks (W3 limpa: 34 on base, 16 on
head; A111-01 limpa: 2 failures and 1 error, the index [[14]] on base); no other test differs. test-gui-entry-proteomics.R message-changed, 92/101 on both: every head-only line is
a declared leaf of the ticket's section 8 (composite artifact keys at 11 steps and in the Copy steps, pre_norm_* at norm, the eight post_filter outputs at norm and
correlation, dia-qc store counts, the nine ruv_covariates_summary lines); the four file-name key lines are KIT-271's and equal on both sides. Check: identical (3 WARNINGs, 3 NOTEs, 0 ERRORs on
both). test-omic-norm-module.R same; every other file same.
UNPREDICTED: test-kit272-prot-box-panel-text.R, 373/0 on base, 372/1 on head (new-failure): T7 (DIANN_limpa) pins the limpa placeholder sentence "...panel d
(principal component scores by group after log2 transformation)...", while this ticket's design (item 8c, D4) gives the limpa template its own stage words ("of DPC-Quant's
protein estimates"). A consequence of the two tickets' merge, not a defect of the figure. Stop rule 3 applies; the one correction pass is used for it: a declared edit of
that single expectation (limpa: the stage words of this ticket; the other three templates unchanged), re-run on that file only, base and head.

## R2c (run 45759491318) and R3b (run 45759493931): the correction pass

R2c: test-kit272-prot-box-panel-text.R, base ec6e282b (R2 base plus the declared T7 limpa edit) against head 0daa4812: `fixed`, base 372 passed / 1 failed (T7 DIANN_limpa, the old
words "after log2 transformation"), head 373 / 0; as predicted. 39.9 job-minutes (2.66 core-hours).
R3b: the same 19 checks as R3 (run 45750234015), all as predicted, now with the two limpa PNGs uploaded and opened. Three-column figure (1950 x 3900): titles "DPC-Quant estimates (log2)" (2 lines),
"After the peptide count filter (no between-sample normalisation applied)" (3 lines) and "RUV-III-adjusted (for quality control only; not the data analysed)" (3 lines) all fit their columns, none clipped (the
three-line titles come close to the top edge of the strip but are whole); panels a to l and m) present, m) centred. Two-column figure (1500 x 3300): the first two titles fit, panels a to h.
The cyclic loess title (first line 32 characters) is not rendered by this lane and is unobserved. The panel d to f box plots match KIT-272's description. 30.5 job-minutes (2.03).

## Budget (core-hours = job-minutes x 4 / 60)

G0 1.79; R1 1.85; R1b 1.38; M1 0.89; R1c 2.83; M1b 0.75; R3 1.64; M1c 1.45; R2 38.52; R2c 2.66; R3b 2.03 = 55.79 against the ceiling 163.46 (expected 48.03; over the expected
by the extra declared runs: R1c 2.83, M1 sends 3.09 instead of 2.0, R3b 2.03, R2c 2.66, because R2's expected 40.57 came in at 38.52).

## Acceptance criteria

A1 met (R1, R1b). A2 met (G0). A3 met: every E and G row passes on the head (R1c 183/2 with the E12 (c) premise corrected, R2 head 210/0). A4 met with one unpredicted verdict (kit272 T7), repaired by a declared edit (R2c).
A5a met (M1c). A5 met (R3, R3b: 19 of 19; drawings 3, 3, 3 and 2; caption; sizes; titles fit). A6 met (E13 and R3 V0/V1/V2). A7: NEWS, ledger row VAL-<merge>, this file; section 8's lines below.

## Stored values that move (not edited; for the end-of-project re-record)

tests/testdata/gui-entry/prot_dia_limpa.json: the limpa lane's artifacts gain proteomics_composite_QC_figure.{pdf,png,rds} at 11 steps (norm, correlation, ruv_covariates, session_export, da_load, da,
summary_refusal, summary_save, summary_copy, summary_report, summary_export) and the two Copy files at the three Copy steps; pre_norm_* appear at norm; the eight *_post_filter_assay1 outputs at norm and
correlation turn from <table> to <plot>; the dia-qc-image/product/state store counts move; the nine ruv_covariates_summary lines change (D5); observed by R2 (45750231427) in test-gui-entry-proteomics.R
(92/101 both sides, messages only). Declared test-literal edits (not stored values): see the ledger row.

## For the change report

Wrong: on DIA with limpa (0.9.0 and staging) the report's Results paragraph described the combined quality-control figure (panels a to l) and the report contained no figure; Copy said the files were copied while its log listed
the figure as missing; the Normalisation tab's first column of QC figures was drawn from a later table, the RUV-III-corrected matrix, or not at all, labelled "Pre-normalisation". Affected: every user of the limpa route.
Code: the builder was called at v0.9.0:R/utils_workflow_dia_products.R:867 and :1900 only (classic); the pre figures' source rule at v0.9.0:R/mod_prot_norm_qc_generation_helpers.R:8-11; now the limpa closures call the builder
(R/mod_prot_norm_workflow_helpers.R: pass_through closure, runProtNormRuvCovariatesObserver finish, runProtNormApplyCorrelationObserver, runProtNormSkipCorrelationObserver; R/utils_workflow_dia_norm.R: buildDiaLimpaCombinedQcFigure; worker:
R/utils_workflow_dia_products.R saveWorkflowDiaLimpaCompositeFigure). Results unchanged. Evidence: runs in runs.tsv; PNG sizes 1500x3300 and 1950x3900; M1 rises 0.007 and 0.010 GiB. On main the single normalisation module built the figure for limpa projects.

## Not checked

The cyclic loess image title (32-character first line); a real worker on Windows or with systemd (M1 used a child Rscript, the production command shape without the unit); sessions saved before this change beyond the labelled-note path; the DA tab at a re-applied correlation_filtered head (VAL-334); the
G2b fallback beyond the shape tested; Copy's log text beyond the assertion on 'Source file not found'; test-omics-module-globals.R (not in scope); R CMD check beyond the identical findings; the real-data behaviour.
