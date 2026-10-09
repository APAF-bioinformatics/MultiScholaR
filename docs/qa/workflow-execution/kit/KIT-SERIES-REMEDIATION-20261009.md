# KIT series remediation record (2026-10-09)

Register: `KIT-SERIES-AUDIT-20261009.md` (unchanged). Authority: user, "fix all blockers" (2026-10-09). Scope of this pass: the ten blockers only. Should-fix and note items stay open.

| Finding | State | Closure evidence |
|---|---|---|
| A1 peptide intensity filter unit | closed | `global_statistics` pre-pass added to the declaration contract (KIT-004 item 1); KIT-009 declares the intensity filter's dataset-wide quantile as the pre-pass; plan/README unit tables moved it to feature group after the pre-pass; KIT-014 wording updated. |
| A2 protein missing-value filter unit | closed | Plan/README unit tables and KIT-015 move the `removeRowsWithMissingValuesPercent` family to whole-matrix; KIT-010 proves the row partitioner on the protein replicate filter and declares q-value (also closes A3's contradiction for q-value: row, KIT-010; KIT-014 depends on KIT-010). |
| B1 CORE-pinned harness files | closed | KIT-001 and KIT-005 file lists contain only new harness files (`workflow_kit.R`, `workflow_kit_driver.R`, `workflow_kit_coordinator.R`, `run_workflow_kit.sh`); regression check now verifies the four pinned files' sha256 against `core/bindings.json` and `core-018-floor --validate-only`. |
| B2 populations/protocol circularity | closed | KIT-001 proves snapshot-and-fingerprint on the named screen corpus and freezes protocol v1; KIT-005's counted control journey produces every post-rollup population and oracle, its second member reproduces them; plan ledger and protocol §3/§5 updated. |
| B3 off-ledger oracle runs | closed | No in-process full-size oracle runs remain; oracles come from the counted control journey (KIT-005 item 2; protocol §3). |
| B4 runner identity scheme | closed | KIT-004 item 2: the runner calls `diaComponentScientificIdentity` verbatim with existing operation names and factors; normalisation and RUV stay one declaration; rollback note corrected; `R/mod_prot_norm_ruv_helpers.R` and DESCRIPTION added. |
| B18 no worker in unmanaged launches | closed | KIT-008 item 1 states the dispatch rule (worker when `session$userData$workflow_execution` is present, in-process otherwise) and tests both. |
| B19 KIT-008 before the adapter | closed | Index and dependencies reordered: KIT-011 (deps KIT-004, KIT-007) precedes KIT-008 (dep KIT-011); KIT-012 depends on KIT-008; controls re-pointed (KIT-011 vs KIT-007 build, KIT-012 vs KIT-008 build). |
| B20 member cannot exercise the GUI | closed | KIT-008's member is a `gui`-mode managed launch driven through the rollup tab (shinytest2 against `managed_url`), 10 min; legacy-vs-P14 page identity checked before launch. |
| B29 entry fixture moved to a helper | closed | KIT-013 item 4 keeps a generic production entry hook injected by the managed launch; file list adds `utils_workflow_execution.R`, the KIT driver, `test-perf107-entry.R`; the non-existent test file names are replaced. |

Validation after repair: all 26 KIT tickets pass the schema and OpenAI validators; the selector returns KIT-001; stage-1 counted total recounted to 34.

## Should-fix and note pass (same day; user: "continue fixing")

| Finding | State | Closure evidence |
|---|---|---|
| A3 q-value unit | closed | Row on `partitionByRows` (KIT-010); KIT-009 no longer claims it; KIT-014 depends on KIT-010; plan/README tables agree. |
| A4 accession cleanup | closed | KIT-015: DIA branch = store, other branch = whole-matrix; `aa_seq_tbl_final` is a KIT-013 global pair. |
| A5 duplicate resolution | closed | KIT-015: verbatim extraction of the inline dplyr into a kernel with an exactness proof; feature group by `protein_id_column`. |
| A6 IQ side product | closed | KIT-015: PDF disabled per batch or declared a product; column-set comparison. |
| A7 per-component merge | closed | `merge` rules in the declaration (KIT-004 item 1); KIT-009 and KIT-014 exactness on every component; protocol and README state it. |
| A8 oracle scope | closed | KIT-005 and protocol: every component, support entry, summary, product; KIT-004 adds publication-file hashes and a holder-field snapshot; zero-row safe fingerprints. |
| A9 KIT-040 library sizes | closed | Declared `global_statistics` pre-pass for `filterByExpr`/CPM; TMM whole-matrix. |
| A10 enrichment identity | closed | KIT-016: service/database version in the identity or non-cacheable; exactness on the test double. |
| A11 Limpa | closed | Removed from KIT-015 and KIT-017; owned by KIT-050 once a Limpa population exists. |
| A12 closure rule | closed | Per-pair rule with KIT-005's baseline pair as the noise estimate and a per-step 1.02 guard (KIT-017, protocol, plan). |
| A13 journey parameters | closed | KIT-001 freezes `diagnostic_settings = FALSE`, `ruv_pearson` and per-step settings in the bindings. |
| A14 state-manager mode | closed | KIT-001 freezes adapter mode for journeys and records it as a claim limit. |
| B5 real kernels | closed | KIT-004 names the kernels, splits the two GUI adapters, adds `mod_prot_norm_ruv_helpers.R` and a callr test. |
| B6 DESCRIPTION | closed | Added to KIT-004, KIT-006 and KIT-015. |
| B7 metadata gate | closed | KIT-002: generic `readComponentMetadata`; `readDIAStateMetadata` stays a DIA wrapper with the same error class. |
| B8 hydrate defaults | closed | KIT-002 keeps the protein-only defaults; any change is a later tested behaviour change. |
| B9 declaration field | closed | Top-level only, never in `spec$slots`; identity-digest test. |
| B10 slot kinds | closed | KIT-003: `nested_object` and `opaque` kinds; serialize-hash exactness; nested component layout. |
| B11 DA output class | closed | The existing `dia.da_results.v1` manifest is the DA output; KIT-012 depends on KIT-003. |
| B12 filter oracle | closed | KIT-005: `verifyKitMember` verifies pages, summary, support and design hashes; screen check first; "existing operations" until KIT-007. |
| B13 sample-filter mode | closed | KIT-006 names `output_mode = "physical"`; production default unchanged (KIT-008); physical-vs-native screen check. |
| B14 partitioner shape | closed | `partitionBySamples(scope, declaration, stage, consume, cap, phase_labels)`. |
| B15 estimator | closed | Estimator per declaration; P14 batch count equals the control's; reduced-cap test forces > 3 batches. |
| B16 registry schema | closed | KIT-007 entry carries handlers, step label, byte limit, resolver, variants, diagnostic flag; grep lists each site. |
| B17 registry in workers | closed | Registration at `.onLoad`; worker resolution test; sample filter declared explicitly (KIT-006). |
| B21 KIT-008 scope | closed | Kernel switch declared and checked; gui-member validation and fixture observer handled via KIT-013's hook; real test file names; `runPerf107GuiLifecycle` left as retired history. |
| B22 group cap and oracle | closed | Batches from the group index; oracle is the control build's output; sample section routed to `partitionBySamples`. |
| B23 worker plumbing | closed | KIT-007's entry schema (peptide and protein completion handlers); KIT-010 adds `finishDIAProteinAction`. |
| B24 projection and caps | closed | KIT-010: projection unit test on the 300-sample protein population; codec column caps lifted for declared wide assays. |
| B25 KIT-011 scope | closed | DIA-NN/IQ on Linux; Limpa and non-Linux a claim limit; file list completed; real test file names. |
| B26 reclamation | closed | `reclaimDIACompletedMemory` stays reachable; default decided before KIT-012's pair (ledger row). |
| B27 measurability | closed | KIT-001: driver-emitted uptime marks, coordinator trace, teardown mark; journeys under a `gui`-mode managed launch; KIT-012 predeclares the expected resident share. |
| B28 KIT-012 bundle | closed | KIT-012 is holders only; DA products to KIT-016; report tables to KIT-018; rollback note corrected. |
| B30 global scope | closed | KIT-013 bounded to proteomics write/read pairs with every file listed. |
| B31 validator direction | closed | KIT-013 reproduces the failed reopen first, then records the direction. |
| B32 KIT-016 bundle | closed | Split into KIT-016 and KIT-018; `ingestDIANativeLegacySession` named; Limpa/non-Linux handling; KIT-018 depends on KIT-013. |
| B33 journey order | closed | KIT-001 production order; real source path; cleanup in the step map; KIT-009 measures imputation. |
| B34 phase marks | closed | Driver-emitted uptime-clocked marks, no `R/` hook. |
| B35 ceiling | closed | Set from the measured screen duration, never above 90 min. |
| B36 ledger total | closed | 34, corrections capped at one per ticket. |
| B37 control definition | closed | Baseline-pair member kind; `core_018` labelled by source and install; controls defined per ticket; KIT-021 oracle note. |
| B38 gate timing | closed | Confirmations for proof steps before KIT-004 and KIT-006; the gate applies to proof steps. |
| B39 regex count | closed | Moot: KIT selectors live in new files. |
| C1 governance | closed | Approval row for the abstractions; confirmations in the tracked step map. |
| C2 case ownership | closed | KIT-017/021 partial evidence; KIT-050 owns the cases with a non-synthetic corpus; default activations recorded as scoped. |
| C3 admission point | closed | KIT-050 names registry dispatch as the single admission point. |
| C4 stub scope | closed | Split requirement and "declarations only" rule in KIT-020/030/032/040; KIT-031's unit gap pre-recorded. |
| C5 KIT-040 deps | closed | Depends on KIT-031; rebase onto the stage-3 closure. |
| C6 history wording | closed | KIT-011 states entities are immutable and the history list is cut on revert. |

Validation after this pass: all 27 KIT tickets pass both validators; the selector returns KIT-001; index order 001, 002, 003, 004, 005, 006, 007, 009, 010, 011, 013, 008, 012, 014, 015, 016, 018, 017, then stages 2–4 and the gates. No finding remains open.

