# Remediation record: analysis-record series audit (2026-10-14)

Closes the findings in `AUDIT-20261014.md` (which is not edited). Authority: user "fix everything" (2026-10-14). Every closure is a ticket revision 1.1 or a plan/manuscript edit; no code was changed and nothing was executed. The user's decision on AR-14: the DIA imputation method (limpa or otherwise) is not yet chosen, so the manuscript stops naming one rather than describing either.

| ID | Closure | Where |
|---|---|---|
| AR-01 | Software identity is `build` (`workflowExecutionBuildIdentity()`) + `version` + `source_commit` from an install-written `MultiScholaR/SourceCommit` DESCRIPTION field (the runner already records the sha at `run_workflow_checkpoint_proof.sh:459`), `RemoteSha` only as a fallback; acceptance rewritten. | KIT-025 items 2, acceptance; `files_to_modify` gains the runner |
| AR-02 | Provenance anchored on every `state_history` writer, enumerated by grep in Simulation Setup; known sites listed; `publishDIAStateHead` named as not the anchor. | KIT-025 item 1 |
| AR-03 | Variant derived at the `completeWorkflowAction` caller from `job$request` and in `runSavedOperation` from `parameters$variant`. | KIT-025 item 1 |
| AR-04 | Seed capture removed; `seed = NA` for every DIA step; a no-RNG source test over registered kernels; declared-seed rule in the plan. | KIT-025 item 5; plan science section |
| AR-05 | Per-entry `defaults` map (parameter → `config.ini` section/key) on `workflowOperationEntry`; provenance computed only through it; registry test that every declared parameter is mapped or `fixed`. | KIT-027 item 4, acceptance; `files_to_modify` gains `utils_workflow_operations.R`, registry test |
| AR-06 | Replay scope defined: import and design are verified preconditions (entity id + input hash), registry operations are replayed in `state_history` order; the browser driver is not reused; `kitReplayJourney` added as an in-process harness kind; state count corrected (17 names, 16 distinct); DA UI-args commit replayed through `DA-contrast-output`. | KIT-028 items 1-2, setup, acceptance |
| AR-07 | `renderWorkflowAnalysisRecord` calls `rmarkdown::render` on `system.file("reports","shared",...)` with `params`; `RenderReport` untouched and removed from the allowlist. | KIT-029 item 1 |
| AR-08 | Report hashed during import adoption (cost recorded); GUI-built design recorded by `design_initial`'s entity id (`design_builder`), `design_file` only on the directory-import path; FASTA hashed at import and declared as `fasta_sha256` on the accession cleanup (its oracle recomputed under the generator; KIT-015 binding not reused). | KIT-026 item 1, acceptance; allowlist widened |
| AR-09 | `release_kind` added (`database_version`, `package_version`, `api_version`, `header`); `query_date` always; UniProt header read once, recorded only if present and flagged unverified. | KIT-026 item 2, acceptance |
| AR-10 | Acceptance judges `full_s4_reads` only; `dimensions` from `sizes[[state]]`; lineage from `previous_state`. | KIT-027 item 2, acceptance |
| AR-11 | The DOCX/XLSX/JSON bundle added to `kitFingerprintFiles`' exclusion list in the same change; acceptance that the file fingerprint is unchanged by the bundle. | KIT-029 item 3, acceptance |
| AR-12 | Real-render tests that skip without pandoc, run and recorded on the host with the pandoc version; the mocked tests named as insufficient. | KIT-029 setup/acceptance; KIT-033 setup |
| AR-13 | `step_label` falls back to the state name; the five `NULL`-label entries named; acceptance that every step has a label. | KIT-027 item 2, acceptance |
| AR-14 | Manuscript: the Methods no longer name an imputation method for the DIA path (user: undecided, limpa or not); a TODO marks the decision; the Tyanova citation removed from that sentence. KIT-029's template describes steps from the registry and the record, never from prose. | `Manuscripts/MultiScholaR/main.tex`; KIT-029 item 1 |
| AR-15 | Claim limit stated at the series end and in the manuscript plan: DIA path only until KIT-031/032/040 declare their steps. | KIT-033 handoff; `notes/manuscript-plan.md` |
| AR-16 | Environment table keyed by `build`. | KIT-025 item 2 |
| AR-17 | "18 screen states" replaced by `state_history` membership in every ticket. | KIT-025/026/027/028 |
| AR-18 | Edited-record criterion names the peptide q-value threshold (0.01 → 0.001) with the expected post-filter row count recorded from a direct kernel call before replay. | KIT-028 setup, acceptance |
| AR-19 | Plan rule notes the replay member kind is not in protocol v1.4 and is a user decision when KIT-028 lands. | plan science section; KIT-033 handoff |

Open after remediation: none blocking. Two items are recorded as unverified until implementation: whether the UniProt REST API sends a release header (AR-09) and the measured cost of hashing a proof-size report during adoption (AR-08). Both are reported in the respective outcome files, not assumed.
