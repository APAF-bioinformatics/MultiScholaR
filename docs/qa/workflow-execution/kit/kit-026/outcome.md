# KIT-026: inputs and external resources in the head

Status 2026-10-17: **complete.** Control `kit_036` (unchanged; no counted member).

## What changed (one mechanism: the head names the analysis's external inputs)

- **Inputs by role.** `file_dependencies` entries carry `role`. The two existing execution entries are `iq_output` and `peptide_input`. `recordDIAWorkflowInputs(manager, paths, entities)` (`R/utils_workflow_dia_state.R`) adds file inputs (name, bytes, sha256) and saved inputs (entity id):
  - `diann_report` and `fasta`: the import call site names its two paths by role; `completeDIAInputFiles` hashes them after parsing, once the leases are verified unchanged.
  - `design_builder`: on the GUI-built path, the design builder writes it into the same head as `design_initial`, with the entity id of `design_initial$data_reference`.
  - `design_file` and `fasta`: the directory import (`runProtDesignImportConfirmationFlow`) records them from its files.
- **Inputs are not archive members.** A recorded input has no `source`, so `resolveDIANativeDependencies` leaves it out of the export: the report and the FASTA are named by hash, not copied. Entries written before roles existed open with `role = NA`.
- **FASTA identity on the accession cleanup.** `fasta_sha256` is a declared parameter (declaration, registry entry, validator). When the caller gives none, it defaults to the head's `fasta` input, else `NA`. It is recorded in the state's `config` (revision 1.2: args carry forward into every later state; config does not).
- **External resources.** `recordWorkflowExternalResource(manager, service, release, release_kind, query_date, query_summary)` appends to `snapshot$external_resources`; `query_date` is ISO 8601. Workers never publish the head, so they record nothing. Callers:
  - Enrichment product save: each reported service, g:Profiler `database_version` or clusterProfiler `package_version`.
  - STRING submit: `string-db`, `12.0`, `api_version`.
  - UniProt fetch: the first successful response's `x-uniprot-release` header, `release_kind = "header"`, NA when absent. Whether the live service sends that header is unverified (no live call made).
- **Harness.** The canonical form omits `fasta_sha256` (revision 1.2). The oracles are unchanged and the accession state's tables are still compared.

## Evidence

| Check | Result |
|---|---|
| In-process check on the screen fixture ([inputs-check.R](inputs-check.R) → [inputs-check-screen.json](inputs-check-screen.json)), from the source tree: KIT-008 member project reflinked; every state fingerprinted; store reopened on this build; screen corpus report and FASTA recorded by role; IQ rollup rerun in process from `imputed`; accession cleanup rerun with the bound state's parameters | `diann_report` and `fasta` sha256 equal to `sha256sum` and the bindings' corpus hashes, bytes equal; the rerun IQ commit writes `iq_output` and `peptide_input` with sha256 equal to `sha256sum` of their files; the accession cleanup was applied and its `config$fasta_sha256` equals the head's `fasta` input (not in `slots$args`); provenance parameter names `delimiter`, `aggregation_method`, `fasta_sha256` all present in config; pre-KIT-026 entries reopened as `role = NA`; 18 fingerprint comparisons identical (16 states after recording, the two recomputed states) |
| `test-workflow-native-session.R` (new: roles, entity and file design inputs, hashes, no archive copy of recorded inputs, export then `openDIANativeSession` restores `file_dependencies` and `external_resources`, `role = NA` backfill, workers record nothing) | 521 pass, 2 skip |
| `test-workflow-dia-products.R` (enrichment on the test double: one `gprofiler2` / `gprofiler2_test_double` / `database_version` entry, ISO date) | 176 pass |
| `test-workflow-prot-declarations.R` (`fasta_sha256` NA without a recorded FASTA, the head's hash with one) | 119 pass |
| `test-multiomics-enrich-shared.R` (STRING submit under an owner records `12.0`, `api_version`) | 91 pass |
| `test-prot-annotation-uniprot-direct-shared.R` (one `uniprot` entry for three chunks, the header value) | 20 pass |
| Enrichment module contracts: `test-prot-11-enrichment-module-contracts.R`, `test-prot-11a-…`, `test-module-ci-prot-enrich.R`, `test-e2e-enrichment-doubles.R` | 1253, 143, 92 (1 skip), 34 pass |
| `test-workflow-operation-registry.R`, `test-workflow-operations.R`, `test-workflow-kit-harness.R` | 140, 69, 332 pass |
| Full suite at `7909897b` from a clean worktree ([logs](check-logs/)) | 343 files, 29,754 pass, 0 fail, 140 skip; the one error is the PERF-101 runner artefact; `R CMD check` 3 warnings, 3 notes (baseline). A first suite attempt wrote nowhere (relative output path; runner setup error, no test ran) and was restarted. |

Install `kit_026` (`run_workflow_kit.sh install`): `browser-smoke.qyJkTq`, installed `source_commit` `7909897b126e…` verified, the registry declares `fasta_sha256`; bindings regenerated (control `kit_036` unchanged).

## Cost of hashing the report

The full corpus report (`report.tsv`, 1,223,346,269 bytes) hashes in 9.6 s in process (`digest`, the helper's no-log path), with process RSS 57 MB. `sha256sum` takes 1.4 s on the same file. The import's memory is unchanged. Time is reported, not judged.

## Claim limits

The design builder entry is unit-tested; the screen fixture's design was built before this build, so the check does not rebuild it. A STRING or UniProt query inside a worker is not recorded. The UniProt release header is unverified against the live service.
