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

## Not checked

- The running application (no GUI drive); the notice and error texts on screen.
- Sessions saved before KIT-182 and legacy RDS projects (the restore's name fill, F1).
- Whether "Save Design" is reachable after normalisation in the application (E5 is at the function level).
- R3 (full suite).

## For the change report

- **What was wrong.** The proteomics Reset went back to the newest of ten fixed state names, not to the table the
  normalisation started from. Label-free (MaxQuant, FragPipe, Spectronaut) and TMT users who normalised before the
  Protein Replicate Filter and then pressed Reset were returned to the design table; the next normalisation, RUV, DA
  and enrichment ran on unfiltered data with no warning, while the report's Methods still described the filters.
  DIA with limpa users lost five Quality Control steps on every Reset (later steps then refused; no numbers). DIA (iq)
  users normalised before the replicate filter lost the protein roll-up and protein QC (the next run was refused).
  Reset before any normalisation also cut Quality Control back; after a design saved again it could return the
  previous design's table.
- **Code.** Base cf8af245: `R/mod_prot_norm_session_helpers.R:352-392` (revert by name list), `:443-454` (notice),
  `:467-470` (call). Fix 17643d77: `:364-404` (revert by recorded source identity; no-op without a normalisation-tab
  state; refusal `prot_norm_reset_source_unknown`), `:457-474` (notice), `:486-490` (call passes `normData`).
- **What changed in results.** Only for the affected users: the table after Reset, hence every later number.
  Unchanged where the source was the newest listed name (QC complete on DIA, label-free, TMT: G1, G2, KIT-182 T5).
- **Evidence.** G0 45506246515 / 45508719013, R1 45508723674, R2 45510988952; synthetic "small" recipe fixtures;
  function-level drives, not the GUI.
