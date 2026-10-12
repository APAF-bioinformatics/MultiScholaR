# KIT-012: holders to references

Status 2026-10-12: **in progress; stopped before any code change or counted run for a scope decision** (below). Authority: `tickets/kit/00-KIT-PLAN.md` (run ledger KIT-012: a journey pair, control the KIT-008 build, candidate; the reclamation default decided before launch); protocol v1.3. Control: `kit_011` (`core/browser-smoke.AFaqC1`); the pair's control is `kit_008` (`core/browser-smoke.otE3Ck`) as the ledger names it. Predecessor cases: PERF-125, PERF-126, PERF-107 holders; audit findings B26, B27, B28.

## Pre-step done: the KIT-005 oracles under the current canonical form

The KIT-005 `journey_<step>` oracle fingerprints in `bindings.json` were the in-run fingerprints of control-1 (kit_001, 2026-10-10) under the harness's canonical form of that day. The form has since dropped `declaration` (KIT-011 correction 1) and `core_utilisation`, the args-embedded support copies and the step description (KIT-008 correction 2), so a record-less journey candidate verified against those oracles would mismatch on bookkeeping on every state, as the KIT-008 member did against the KIT-011 screen member before its correction 2 (hand-off recorded in the KIT-008 outcome and the ledger).

The generator now recomputes the state fingerprints from control-1's immutable store (`kitFingerprintStore`, the same function the verifier uses) under the current form and overlays them on the in-run record as the verifier overlays a candidate's: class, rows, components, support and the args summary from the store, the state record's products recomputed, the per-step assets the driver attached kept from the in-run fingerprint (`bindings$journey$baseline$oracle_canonical_form`). Result: the args summary changed on all 17 states with a record (every record carries fields the current form drops) and one product (`normalization_science`, the KIT-011 correction 1 field); every component and support hash is unchanged. Bindings regenerated; harness tests 241 green.

The decision for a ticket's journey pair takes the ticket's predeclared `rule`: `whole_slice_anon` (KIT-011, default) or `resident_end` (KIT-012: `resident_end_bytes` ≤ 0.95 × the control's, whole-slice anonymous peak ≤ 1.02 ×, the per-step guard from q-value onward, `invalid` when the control has no resident-end reading). The KIT-012 binding (`tickets$012 = {journey, control kit_008, candidate kit_012, record FALSE, rule resident_end}`) is written by the generator once a `kit_012` build is bound. Test added (`test-workflow-kit-harness.R`).

## Holder survey on the measured path (acceptance 1)

The ticket's worklist is the KIT-001 holder inventory (2026-10-09), written when the GUI's default owner was the in-memory `WorkflowState`. Since KIT-011 the DIA-NN import attaches the on-disk owner with reference states, and every reference-aware holder resolves to a reference on that path. Re-read under the default (`usesDIAReferenceStates(workflow_data$state_manager)` TRUE):

| Holder (KIT-001 worklist) | On the DIA default | Where |
|---|---|---|
| `da_data$current_s4_object` | checkpoint reference (`readDIAStateGuiMetadata`, `resolveDIANativeSessionView`) | `mod_prot_da_handlers_init.R`, `mod_prot_da_handlers_load_session.R` |
| `da_data$da_results_list`, `workflow_data$da_analysis_results_list` | the DA manifest (`dia.da_results.v1`: references to the per-contrast tables and products; 98 kB at 300 samples), one object on both fields; the `theObject` list exists only on the in-memory branch | `utils_workflow_dia_da.R` (`runDIAReferenceDaAnalysis`) |
| `enrichmentData$current_s4_object`, `da_results_data` | `manifest$source` (reference) and the manifest | `mod_prot_enrich_registration_helpers.R` (both observers take the manifest branch) |
| `createProtNormReactiveState` fields (`normalized_protein_obj`, `ruv_normalized_obj`, `correlation_filtered_obj`) | references through `setDIAOwnedPayload` / `saveDIAOwnedValue` (an S4 is checkpointed and the reference stored); the observer restores references with `getDIAStateReference` | `mod_prot_norm_*_helpers.R`, `utils_workflow_dia_products.R` |
| `workflow_data$ruv_normalised_for_da_analysis_obj` | reference (`getDIAStateReference` / `saveDIAOwnedValue`) | `mod_prot_norm_correlation_helpers.R` |
| design builder `resultRv`, `dataClnReactive` | `savedReactiveVal` references | `mod_prot_design_builder_server_helpers.R` |
| import `data_tbl`, `data_cln` | table references (`setDIAOwnedPayload`) | import helpers |
| peptide and protein QC modules | no S4 holder: `.savePeptideQcState`/`saveState` commit the state and the local object leaves scope; the modules keep one plot grob each | `mod_prot_qc_*.R` |
| hydrated working copies | the adapter's working environment is cleared at every flush and at every scope exit | `attachDIAStateSession`, `withDIAStateScope` |

The in-memory branches (`theObject` in the DA list, `getState` copies in the enrichment and session helpers) remain for LFQ/TMT, Limpa and non-Linux sessions, which keep the in-memory `WorkflowState` (KIT-011 claim limit); they are not on the measured DIA path. Grep for the worklist names in `R/mod_prot_*.R`, `R/mod_proteomics.R`, `R/utils_workflow_dia_*.R` is the line list this table was read from.

Consistent with the KIT-011 candidate's trace: the slice's anonymous memory at each step's completion does not accumulate (design 2.12 GiB, replicate 1.88, protein_replicate 2.28, correlation 2.24, session export 2.26 GiB); a holder kept across steps would raise it monotonically. The end reading (3.02 GiB) is 0.76 GiB above the pre-DA level, so the DA step owns the question.

## What the DA step leaves resident (B27 predeclaration evidence)

[da-resident-probe.R](da-resident-probe.R) → [da-resident-probe.json](da-resident-probe.json) (not counted): the DA tab's production path replayed in process on a reflink copy of the KIT-011 candidate's project (300 samples, `correlation_filtered`): native Open of the journey's session export, `runDIAReferenceDaAnalysis`, the results-table page, the static volcano, the Glimma volcano, the heatmap (NULL: no protein under the threshold in the synthetic corpus), then every module holder released. At each boundary: the process's anonymous resident memory as sampled, after `gc(full = TRUE)` and after the allocator released its free pages (`workflowAllocatorRelease`), the live R heap, and each holder's kind and size.

| Boundary | Anonymous (sampled) | After allocator release | Live R heap |
|---|---|---|---|
| session loaded (33 s) | 0.60 GiB | 0.58 GiB | 0.37 GiB |
| DA analysis (94 s) | 1.40 GiB | 1.13 GiB | 0.73 GiB |
| table page, static volcano, Glimma volcano, heatmap | 1.15 to 1.24 GiB | 1.12 to 1.19 GiB | 0.72 GiB |
| every module holder set to NULL | 1.13 GiB | 1.13 GiB | 0.72 GiB |

Holders at the end of the DA step: `current_s4_object` a 5.6 kB reference, the manifest 98 kB, `contrasts_tbl` 1 kB, `design_matrix` 80 kB, `config_list` 26 kB; the static volcano ggplot 17.5 MB and the Glimma widget 7.6 MB are render values the module does not keep. Releasing every holder changes nothing: the DA step adds about 0.55 GiB of anonymous memory, of which about 0.36 GiB is live R heap that no module holder owns.

[da-namespace-check.R](da-namespace-check.R) → [da-namespace-check.json](da-namespace-check.json): loading the DA step's packages in a fresh process. limma, openxlsx and htmlwidgets add nothing (already loaded); **Glimma loads 20 namespaces (DESeq2, edgeR, SummarizedExperiment, GenomicRanges, IRanges, S4Vectors, DelayedArray, Biobase, BiocParallel, …) and adds 0.44 GiB anonymous memory and 0.30 GiB live R heap** (4.5 million R nodes: code, S4 class and method tables). The DIA journey's DA step writes the interactive volcano HTML through Glimma (`writeDIAInteractiveVolcanoOutputs`) and renders it in the tab, so every DIA journey pays this once, in the application process, for the life of the session. The remainder of the DA step's leftover (about 0.2 GiB) is allocator and DuckDB working memory that `gc` and a page release do not return.

A walk of every loaded namespace, the store's caches, the session's user data, Shiny's globals and the options after the DA step ([da-resident-probe-walk.json](da-resident-probe-walk.json)) finds no other owner; the walk itself forces the lazy-loaded data of the newly loaded namespaces (about 0.35 GiB more), which is why the clean record is the run without it.

## Scope decision needed before any code or counted run

1. **Acceptance 1 is met by KIT-011's work**: on the DIA default no module holds an S4 across steps (table above). There is nothing left for this ticket's mechanism to convert on the measured path.
2. **Acceptance 2 cannot be reached by this mechanism**: the required reduction is ≥ 5 % of `resident_end_bytes` (about 150 MB of the 3.02 GiB the KIT-011 candidate ended at). The end-of-journey resident memory above the pre-DA level is the Glimma namespace tree (0.44 GiB) plus allocator and DuckDB working memory, not holders. A `kit_012` build with no holder change would measure at the KIT-005 noise (control-1 vs control-2 resident end 0.983) and the pair would be `rejected` on `required_reduction` by construction.
3. **Reclamation (B26)** is not reachable in production: `reclaim_enabled` is set FALSE at every adapter creation and no code sets it TRUE; `reclaimWorkflowMemory` is a gc plus an allocator page release, which the probe shows returns about 0.2 GiB after DA. Enabling it would be a second mechanism and would not touch the namespace cost. Recommended default for this ticket: unchanged (off).

Options for the user:

- **A (recommended).** Close KIT-012 on acceptance 1 with this record (holders are references on the DIA default; the KIT-001 worklist's in-memory branches stay for the in-memory owner), withdraw acceptance 2 as not measurable by this mechanism, spend no counted member (ledger amendment: the two journey members return to the budget), keep the harness pre-step and the `resident_end` rule for KIT-017's closure pair, and record the Glimma namespace cost as a discovery for a later ticket (interactive volcano written and rendered outside the application process, or Glimma loaded only on demand of the interactive view; KIT-016 owns the DA display products).
- **B.** Run the pair as declared (control `kit_008`, candidate a `kit_012` build with no holder change beyond this record), accepting the predictable `rejected` as the recorded unfavourable result.
- **C.** Retarget KIT-012's one mechanism to the measured owner (the DA step's Glimma publication out of the application process) and judge it with the pair under the `resident_end` rule; that is a new mechanism and a ticket revision, and the expected reduction (about 0.4 GiB of 3.0, 13 %) would be predeclared from the namespace check.

No code was changed in `R/`; nothing was launched.
