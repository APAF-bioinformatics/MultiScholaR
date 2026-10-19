# KIT-020: LFQ/TMT imports on the on-disk state manager

Status 2026-10-20: **completed; implemented at `e61aefd0`** (ticket revision 2.1); the three protein-level fixtures run on the owner from import through DA, their records validate and replay exactly, and a native export reopens with its type; tests, check and install recorded below.

## What changed (one mechanism: the state owner at import for the protein-level formats)

- `initializeDIANativeImport` (`R/utils_workflow_dia_sessions.R`) attaches the owner for every format `workflowTypeForImportFormat` knows (`diann` → DIA, `pd_tmt` → TMT, `maxquant`/`fragpipe`/`spectronaut` → LFQ); an unknown format never attaches; an existing owner accepts a re-import of its own type, judged by the store's workflow identity (the manager's own type is the fresh default `LFQ` until the import sets it, so a comparison on it refused the first DIA-NN import on a pre-built owner in the GUI lifetime tests; found by the dependents run, repaired before commit), and refuses a change of type; the store is `<FORMAT>-native` with the type as its workflow id (`DIA-IQ-native`/`DIA` unchanged). The Limpa refusal and the Linux gate are untouched (KIT-046, KIT-044).
- `adapter$setWorkflowType` (`R/utils_workflow_dia_state.R`) accepts `DIA`, `LFQ`, `TMT`.
- The native session manifest carries the head's type and a route by type (`workflowRoutesByType`: DIA → `iq`, LFQ/TMT → `protein`); `openDIANativeSession` accepts the three (DIA manifests unchanged).
- The import's owned path (`applyProtImportResultToWorkflow`, the import observer's `deferFasta`) covers every owned format: the table is owned by the store, the import metadata recorded (`recordDIAImportMetadata` works on the protein-level readers' long tables: `Run`, `Intensity`/`Abundance`), the FASTA processed after the table is released. The roles `<format>_report`/`fasta` were already recorded for every format (verified).
- The design step: `adapter$saveState` already commits an S4 as a component checkpoint; the LFQ/TMT branch of `buildProtDesignStateCheckpoint` now carries the KIT-038 provenance (`design_builder`, fixed pivot and transform). The fallback S4-creation tab (`runProteinS4CreationStep`) reads the owned table through `getDIAOwnedPayload` and records the same provenance.
- The two `protein_s4_initial` probes in `mod_prot_qc.R` read the state reference on an owner instead of hydrating.
- Replay accepts the four report roles.

## Evidence

- Owner check ([owner-check.R](owner-check.R), [owner-check.json](owner-check.json), [log](owner-check.log)), the `prot_lfq` MaxQuant and FragPipe fixtures and the `prot_tmt` Proteome Discoverer fixture, each on a fresh owner: store `MAXQUANT-native`/`FRAGPIPE-native`/`PD_TMT-native` with workflow id `LFQ`/`LFQ`/`TMT`; `protein_s4_initial` committed as a component checkpoint whose hydrated object is identical (table, design, slots and args) to the `ProteinQuantitativeData` built in memory from the same table; provenance `design_builder` with the fixed transform; inputs recorded by role with sha256; the protein intensity, duplicate and replicate declarations, normalisation (RUV skipped), correlation and DA (one contrast, 4 rows) run on the owner through the registry handlers and the reference entry points; the record validates with no legacy step and `workflow` LFQ/TMT; **replay: 7 states, 0 mismatches** on each; a native export reopens through `openDIANativeSession` with the reopened owner typed `LFQ`/`TMT` at `correlation_filtered` (a DIA/iq-only manifest check would have refused it).
- DIA route unchanged: the KIT-038 record check rerun on the kit-journey-028 project on this tree ([record-check-kit020-dia-unchanged.json](../kit-038/record-check-kit020-dia-unchanged.json)): the record validates, replay **16 states, 0 mismatches**; the same store, head provenance and fingerprints as the KIT-038 run.
- Tests (proof library, `devtools::load_all`): `test-workflow-native-session.R` 547 pass (the format test extended; one new KIT-020 test), `test-workflow-checkpoint-dia-state.R` 236 (the refusal case now `metabolomics_standard`), `test-prot-01-import.R` 15, `test-prot-01c-import-module-contracts.R` 410, `test-prot-01d-import-module-characterization.R` 79, `test-prot-04-design.R` 1571, `test-workflow-analysis-record.R` 177, `test-module-ci-prot-import.R` 90, `test-prot-02o-qc-protein-s4-shared.R` 30; `test-e2e-proteomics-tmt-lfq.R` skips (no corpus). Dependents grep (`initializeDIANativeImport`, `setWorkflowType`, `buildProtDesignStateCheckpoint`, `recordDIAImportMetadata`, `stageDIANativeSession`, `protein_s4_initial`): 13 files, all green after the gate repair (`test-workflow-gui-memory-lifetimes.R` 252 pass; it had 7 failures from the first gate).

## Deviations and claim limits

- The in-memory equivalence on the protein-level fixtures is proven at the design state (identical objects) and for the downstream steps by the shared declarations, which are the same code over `ProteinQuantitativeData` as on the DIA route (KIT-005 oracles); the full in-memory cross-check on an LFQ corpus is KIT-021's. The owner check runs RUV in `skip` mode because a five-protein fixture has too few proteins for the negative-control sweep.
- The accession cleanup is not run in the owner check (it needs the FASTA table the GUI's import produces); it is the same declaration the DIA route runs after IQ and is exercised by KIT-021's journey.
- The TMT reader names runs `126_WT_1`; the design must use those names (the GUI's design builder builds on the imported runs, as the check does).
- The fallback S4-creation tab records `transform = "none"` (it does not log-transform), the design builder `log2(x + 1)`.

## Check and install

- `R CMD check --no-manual --no-tests` at `e61aefd0` (clean worktree, `_R_CHECK_FORCE_SUGGESTS_=false`, proof library on `R_LIBS`): 3 WARNINGs, 3 NOTEs, the baseline ([log](check-e61aefd0/r-cmd-check.log)).
- Install `kit_020` = `core/browser-smoke.Zc4Uyb` (source commit `e61aefd0` verified on the installed `DESCRIPTION`); bindings regenerated. The control stays `kit_036` (no counted member).

## Hand-off

KIT-042 adopts the owner and the typed manifest on every workflow type; KIT-043 renders the record on LFQ/TMT projects; KIT-021 takes the protein-level step list (import → design → intensity → duplicates → replicate → normalisation → correlation → DA) and the in-memory cross-check on a real LFQ corpus. Rollback: `git revert e61aefd0`; LFQ/TMT native sessions written by this build are refused by the previous build's `openDIANativeSession` (its manifest check was `DIA`/`iq` only).
