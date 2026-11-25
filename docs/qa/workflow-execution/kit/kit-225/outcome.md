# KIT-225 outcome (r1.2): proteomics "Reset to Pre-Normalization" returns to the table the normalisation started from

Branch `kit-225` from `v0.9-staging` cf8af245. Commits: e3bb2b1e (G0 script), 398b9af9 (tests only: new file, the
declared 05b edit, G0 instrument repair), 17643d77 (fix, NEWS, ledger, wording drafts). Runs: `runs.tsv`.
State: ready for the implementation audit; R3 (full suite) not run, by instruction (own or shared with KIT-224).

## Observed

- **G0** (run 45506246515, then VAL-54 instrument rerun 1 of 3, run 45508719013; `g0-base-values.tsv`, attempt 1 in
  `g0/`). Base Reset selected: E1, E2 `protein_s4_initial` (P0; the second run's source P0); E3 `precursor_rollup`;
  E4 `imputed`; E5 `protein_replicate_filtered` holding the previous design's table (P4); E6, E6b, E6c
  `protein_s4_initial` (Quality Control cut back with no normalisation run); E7 no error, `protein_replicate_filtered`;
  G1, G2 history, head, provenance, source and tab status as before the run. Attempt 1 read E4 to E6c's state before
  the lazily evaluated Reset ran (their `previous_state` was valid); repaired with `force(reset)`, nothing else changed.
- **E3 A1 on base (stop rule 3).** After the Reset to `precursor_rollup`: the pass-through completed and designated the
  roll-up table; correlation Apply refused (`workflow_operation_input`: "pearson-pairs" expects
  ProteinQuantitativeData, the source is PeptideQuantitativeData). Fresh copy: pass-through completed; correlation Skip
  committed `correlation_filtered` and then stopped (`saved_entity_unsupported`: "table read scope requires a Parquet
  descriptor"); RUV covariates refused (`saved_entity_outputs_pending`). No numbers came out.
- **R1** (run 45508723674, tests-only 398b9af9): new file 24 passed, 32 failed, 0 errors, each failure an assertion at
  G0's value; G1, G2 passed; 05b failed only at `calls$revert$normData`.
- **R2** (run 45510988952, 17643d77 against 398b9af9; 30 files; `r2/`): 28 `same`, new file and 05b `fixed` (56/0,
  440/0); `test-kit129-run-path.R` (6) and `test-prot-norm-core.R` (306) fail identically on both sides; check 0
  errors, 3 warnings, 3 notes on both, nothing only on head.

## Acceptance criteria

| # | Status | Evidence |
|---|---|---|
| A1 | met | `g0-base-values.tsv` (G0-2), E3 A1 rows above |
| A2 | met | R1: every E row and the 05b expectation fail at G0's value; G rows pass; R2 head: all pass |
| A3 | met | R2: every file `same` except the new file and 05b (`fixed`); check: nothing only on head |
| A4 | not run | R3 left to the coordinator (own, or shared with KIT-224 under E-4) |
| A5 | met, R3 pending | NEWS (Fixes), ledger `VAL-<merge>`; `man/`, `NAMESPACE`, `DESCRIPTION` unchanged; `wording.md` |

Guards: G3 (KIT-182 T5, T6) and G4 (05d history-read error) are `same` in R2 (in `test-prot-norm-input-current.R`,
`test-prot-05d-...`); G5 (05b no-manager test; `test-omic-norm-module.R` KIT-122 F6, F10) `same`.

## Inferred (not executed)

- In the application the manager is always a reference owner, so the refusal and the no-op are reached as in the tests.
- The managed path records the same source (`normalization_input_reference`), so the same selection follows.
- DIA with limpa, 'Run Normalization & RUV' pressed again while the correlation filter's result is current: the
  pass-through keeps its earlier input only at `ruv_covariates` or `ruv_corrected` (`R/utils_workflow_dia_norm.R:604`,
  `:628`), so the recorded starting table becomes the `correlation_filtered` state itself. A Reset then selects that
  state (the recorded source, as VAL-203 decides), changes no state, clears the tab, and its notice names '(the
  peptide quality-control steps)', which is not true. Base returned `precursor_rollup`. Not tested. Candidate (not
  this ticket): add `correlation_filtered` to the pass-through's kept-input list at both lines; that also changes the
  table a correlation filter re-applied after such a press reads (a results change for that sequence, the owner's).
- E3 on base (stop rule 3) shows a defect outside this ticket that stays reachable after the fix: on DIA with limpa,
  when the pass-through designates a peptide table (on base after every Reset; before and after the fix whenever 'Run
  Normalization & RUV' is pressed before DPC Quantification, which the tab allows with a warning, VAL-44), the
  correlation Apply refuses with a class message, the Skip commits the peptide table as `correlation_filtered` (U2)
  and stops in its publication, which writes the protein table (`R/utils_workflow_dia_products.R:1935`), and the RUV
  covariates then refuse with 'correlation outputs are pending; retry that action', which retrying cannot cure. No
  numbers result. Reset (now returning the designated peptide state) or a Quality Control Revert is the way out. The
  fix removes only the Reset route into it. Candidate (with U2): the pass-through refuses a peptide-class designation
  in words.
- Design saved again after a normalisation, with no new run (audit I08, read only): a name re-committed by Save Design
  keeps its earlier history position (`R/utils_workflow_state.R:114-117`), so Reset returns the recorded (previous
  design's) table, as base does by name, or refuses with the true reason when the source's own name was re-committed.
  A rarer order (a normalisation from the design table, a Save Design, new QC names appended after `normalized`, a
  run, then Reset) can leave `normalized` before the head, so a second Reset refuses instead of reporting nothing to
  undo. Whether Save Design can follow a normalisation in the application is not verified; no repair proposed.

## Candidates for the owner (not this ticket)

- I01 (audit AUDIT-KIT-225-impl-20261125): the limpa pass-through keeps its earlier input at `correlation_filtered` as
  at the RUV states (`R/utils_workflow_dia_norm.R:604`, `:628`); a results change for that sequence, the owner's to
  decide.
- I02 (with U2): the pass-through refuses a peptide-class designation in words, and the class messages are reworded;
  no numbers change.

## Not checked

- The running application (no GUI drive); the notice and error texts on screen.
- Sessions saved before KIT-182 and legacy RDS projects (the restore's name fill, F1).
- Whether "Save Design" is reachable after normalisation in the application (E5 is at the function level).
- R3 (full suite).
- Not exercised by any test (audit I09; traced by reading only): RUV-mode runs (every new test uses RUV "skip"), the
  managed path, a reopened project, selection through the `qc_sources$pre_normalization` fallback, and the refusal's
  text (E7 asserts the class only).

## Hand-off to KIT-223

- T4 (a) also applies to `buildProtNormResetNotificationMessage` (`R/mod_prot_norm_session_helpers.R:467-471` after
  KIT-225), which carries its own copy of the fallback for a state without a title: the same return rule (an untitled
  state other than `initial` returns its own name). `test-kit225-prot-norm-reset-source.R` joins KIT-223's R2. Note
  I01: the Reset notice can also reach the fallback with `correlation_filtered` on DIA with limpa.

## For the change report

- **What was wrong (legacy, on `main`).** The proteomics Reset went back to the newest of ten fixed state names, not
  to the table the normalisation started from (`main:R/mod_prot_norm.R:2559-2640`; the list at `:2569`, the choice at
  `:2571`), and the next normalisation reads the current table (`main:R/mod_prot_norm.R:1318-1319`). Label-free
  (MaxQuant, FragPipe, Spectronaut: every import other than TMT and DIA-NN, `main:R/mod_prot_import.R:699-705`) and
  TMT users who normalised before the Protein Replicate Filter and then pressed Reset were returned to the design
  table; the next normalisation, RUV, DA and enrichment ran on unfiltered data with no warning, while the report's
  Methods still described the filters. DIA-NN (iq) users who normalised before the replicate filter were returned to a
  peptide table, losing the protein roll-up and protein Quality Control (on the development branch the next run was
  then refused; what it did on `main` was not checked). Reset pressed before any normalisation also cut Quality
  Control back. After a design saved again it could return the previous design's table (by reading; whether Save
  Design can follow a normalisation in the application is not verified).
- **On the development branch only.** The DIA with limpa route was added during development (KIT-041, KIT-046) and
  does not exist on `main`. On it, every Reset after the precursor roll-up undid five Quality Control steps (the later
  steps then refused; no numbers). Fixed by the same change, before release.
- **Code.** Legacy: `main:R/mod_prot_norm.R:2559-2640`. Development base cf8af245:
  `R/mod_prot_norm_session_helpers.R:352-392` (revert by name list), `:443-454` (notice), `:467-470` (call). Fix
  17643d77: `:364-400` (revert by recorded source identity; no change without a normalisation-tab state; refusal
  `prot_norm_reset_source_unknown`), `:457-474` (notice), `:486-490` (call passes `normData`).
- **What changed in results.** Only for the affected users: the table after Reset, hence every later number.
  Unchanged where the source was the newest listed name (QC complete on DIA, label-free, TMT: G1, G2, KIT-182 T5).
- **Evidence.** G0 45506246515 / 45508719013, R1 45508723674, R2 45510988952; synthetic "small" recipe fixtures;
  function-level drives, not the GUI.
