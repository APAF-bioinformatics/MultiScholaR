# KIT-027: the analysis record as a derived view

Status 2026-10-17: **complete.** Install `kit_027` (`run_workflow_kit.sh install`): `browser-smoke.tXkjVS`, installed `source_commit` `fb1aee89…` verified, the record functions present; bindings regenerated (control `kit_036` unchanged). Control `kit_036` (unchanged; no counted member).

## What changed (one mechanism: the record is built on demand from head metadata)

- **Builder, writer, validator** (`R/utils_workflow_analysis_record.R`):
  - `buildWorkflowAnalysisRecord(manager)` reads the head and each state's metadata record (`readDIAStateMetadata`); no assay is hydrated. The `initial` placeholder is read directly. A state stored as a whole legacy object is listed without its arguments and reported as incomplete.
  - `writeWorkflowAnalysisRecord(manager, path)` writes the JSON.
  - `validateWorkflowAnalysisRecord(record)` returns `valid`, `errors` and `incomplete`.
- **Schema `multischolar.analysis_record.v1`**, as the ticket lists it. Additions: per step `defaults_sha256`, `description` and `stored_as_object`; per parameter `source` (`args`, `config`, `fixed`, `missing`); per input `entity_id` (the GUI design); per external resource `release_kind`; record-level `defaults_sha256`, `defaults_changed` and `limits`.
- **Parameters from both places (rev 1.2).** Each declared name (KIT-025 `parameter_names`) is read from `slots$args[[section]]`, else the state's `config`, and its source is named. The accession cleanup's `fasta_sha256` is read from `config`. A `legacy` step lists its arguments whole.
- **Untruncated values.** Atomic vectors are stored in full with type, names, levels and class; data frames as `{rows, columns, sha256}`; lists recursed. Doubles are written with `digits = I(17)`: the ticket's `digits = NA` writes 15 significant digits, so 0.1 + 0.2 does not restore `identical()`.
- **Provenance through a `defaults` map, never by name.** `workflowOperationDefaults` is attached to each registry entry as `defaults`; registration refuses a map naming an undeclared parameter. A mapping is `list(section, key)` (provenance `default` or `user_set` against the shipped `config.ini`), `"fixed"`, `"input"` (a value from the data or an earlier step: the inclusion list, contrasts, `fasta_sha256`), or `"no_default"` (a setting `config.ini` does not ship, such as the protein intensity filter's group minimums and the aggregation method). The last two categories are added to the ticket's `default`/`user_set`/`fixed`/`unknown`, because these parameters have no shipped key. `unknown` remains only for a step whose operation has no registry entry.
- **Package metadata.** The new file is added to `DESCRIPTION` `Collate` (outside the allowlist; required by the install, which refused the first attempt at `6a91ffb4`: "files missing from Collate"; the failed install directory was removed and nothing bound it).
- **Summary page and export.** `completeProtSummaryWorkflowArgsSave` writes `analysis_record.json` after `study_parameters.txt` for an on-disk owner; a failure is logged and does not stop the summary. `collectDIANativeGuiFiles` exports the file when present. `kitFingerprintFiles` excludes it.
- **Verifier check on the FASTA identity (rev 1.2).** `kitCanonical` no longer drops `fasta_sha256` by name. `kitFingerprintState` removes it from the top level of a state's `config` only, before the summary, and returns it beside the fingerprint. `kitFingerprintStore` sets `fasta_matches_input` against the head's `fasta` input. `kitReproduceMismatches` reports `<step>/fasta_sha256` when they differ. The oracles carry no such field and are unchanged.
- **Stated limit (rev 1.2).** `limits$external_resources`: services queried inside a worker are not recorded. The validator always lists "external resources: controller-process queries only" under `incomplete`.

## UniProt release header (rev 1.2 item 8)

One live request on 2026-10-17 to `https://rest.uniprot.org/uniprotkb/search` (the endpoint and form `directUniprotDownload` uses; query P04637, organism 9606): status 200. The response carries `x-uniprot-release: 2026_03` and `x-uniprot-release-date: 02-September-2026`. The header KIT-026 reads exists under that name, so no package change is needed.

## Evidence

| Check | Result |
|---|---|
| In-process check on the screen fixture ([record-check.R](record-check.R) → [record-check-screen.json](record-check-screen.json), record [analysis_record-screen.json](analysis_record-screen.json)), from the source tree. The KIT-026 check flow: the KIT-008 project reflinked, the corpus inputs recorded, IQ and the accession cleanup rerun in process; then the record built twice. | Record valid; steps = `state_history` (11 after the rerun truncates the history); two builds byte-identical; `full_s4_reads` unchanged; 918 atomic argument values restored `identical()` from the JSON (0 mismatches). The accession step lists `fasta_sha256` with `source = "config"`, equal to the record's `fasta` input. Inputs `iq_output`, `peptide_input`, `diann_report`, `fasta`. Incomplete: the nine pre-KIT-025 `legacy` steps and the external-resources limit. Verifier: `fasta_matches_input` TRUE, no mismatch on the two recomputed states. 18 fingerprint comparisons identical. |
| `test-workflow-analysis-record.R` (new): steps = history with operations from provenance, every `step_label` non-empty; `qvalue_threshold = 0.05` is `user_set`, an untouched key `default`; `proteotypic_only` read from config; `fasta_sha256` source `config` equal to the input; full-length round trip (the 40-entry inclusion list); byte-identical rebuild; `full_s4_reads` unchanged; a legacy store reported incomplete, not invalid; the summary writes the record after `study_parameters.txt` | 41 pass |
| `test-workflow-operation-registry.R` (new: every declared parameter of every registered operation is mapped) | 219 pass |
| `test-workflow-kit-harness.R` (new: nested `fasta_sha256` no longer dropped; the verifier flags a FASTA differing from its input) | 335 pass |
| Summary tests: `test-module-ci-prot-summary.R`, `test-prot-12a-…`, `test-prot-12-summary-module-contracts.R` | 110, 62, 446 pass |
| `test-workflow-operations.R`, `test-workflow-dia-products.R`, `test-workflow-native-session.R` | 69, 176, 521 pass (2 skip) |
| Full suite at `fb1aee89` from a clean worktree ([logs](check-logs/)) | 344 files, 29,877 pass, 0 fail, 140 skip; the one error is the PERF-101 runner artefact; `R CMD check` 3 warnings, 3 notes (baseline). A first suite start at `6a91ffb4` was stopped when the install showed the Collate omission. |

## Claim limits

The screen fixture predates KIT-025, so most of its steps are `legacy`, and its journey has no enrichment, so its record lists no external resources. Those sections are proven by the unit tests and by the two steps rerun on this build. The record describes; it does not compute. Replay is KIT-028.
