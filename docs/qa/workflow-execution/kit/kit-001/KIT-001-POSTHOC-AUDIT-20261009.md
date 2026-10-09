# KIT-001 post-hoc audit (2026-10-09)

Mode: `post_hoc`, audit only. No audited artifact was changed; repairs, if any, go to a separate remediation record. Scope: the completed KIT-001 execution at `1ef2b8b3` against protocol v1.1 (`../protocol.md`, sha256 `8998db8f…`), the plan run ledger (`tickets/kit/00-KIT-PLAN.md`, KIT-001 row and decision rows of 2026-10-09) and the ticket's own acceptance criteria.

## Authority and inventory (mechanical)

| Check | Result |
|---|---|
| `audit_inventory.py` on `tickets/kit` KIT-001 with the protocol and step-map hashes as authorities | 0 findings, exit 0 |
| Harness file hashes in `bindings.json` versus the working tree (6 files) | all match |
| Protocol and step-map hashes versus `bindings.json` | match |
| Four CORE-pinned files versus `core/bindings.json` | all match; `docs/qa/workflow-execution/core/` unchanged since `cf0c8383` |
| Member binding digests (both stages) recomputed from the current bindings | identical to the receipts |
| Build identity of both members | `kit_001` (`browser-smoke.y8rDOi`, source `9075b361`) |
| Launch snapshots versus the launch commits (`ea731185` complete, `3512f84d` fail) | harness file hashes identical |
| Installed control carries the decisions | `resolveProtRuvReplicateBlocks` present; the RUV method uses it; `require_complete_design_runs` present; no `ceiling(quantile` in either intensity filter |
| Exact regressions re-run under a source load | kit harness 54, launch record 16, RUV structure 15, intensity batches 5, peptide methods 119, CORE harness 131: all pass |

## Execution facts (from the receipts and the member store)

- Complete member: 18 steps, 269 s, whole-slice peak 2.76 GiB, anonymous peak 2.15 GiB (normalisation/RUV), resident end 1.86 GiB, sampling gap 70 ms, zero swap, slice stopped, host profile `low-power`, CPU busy 4 %. Journey record `kit-journey-v1`, native sessions on, no qualification. Store: 16 states, 31 products; the sample filter removed `calibration_run_001` (cutoff 3,000) leaving 7 samples; 2,500 proteins after IQ, 2,395 after RUV.
- Failure member: import to intensity complete, guard injected, `injected_guard_stop` by the application monitor, slice stopped, head equal to the intensity-filter snapshot, zero swap, 87 s.
- Decisions: `lifecycle_complete_observed`, `lifecycle_failure_observed`.

## Findings

| ID | Severity | Observed | Affected claim and impact | Smallest sufficient repair | Blocked state |
|---|---|---|---|---|---|
| PA-01 | medium | The complete member's in-run fingerprinting failed (C stack on the RUV plot product); fingerprints were recomputed post hoc from the immutable store by `attempt4-post-hoc-fingerprints.R`; the verifier was rerun from the working tree (not a frozen stage snapshot) and `verifier_completion$accepted` plus `decision.rds` were written by an ad-hoc script recorded only in the outcome and the session, not in a tracked script. In-run snapshots (head, state, entity id, rows) are the driver's for all 18 steps. | The acceptance "snapshot-and-fingerprint mechanism proven on the screen corpus" is proven for snapshots in-run and for fingerprints post hoc. The decision stands with that claim limit (recorded in `outcome.md`). | Record the verify/decide commands in the correction script (successor version) so the correction is reproducible; require KIT-005's lifecycle members to fingerprint in-run. | none for KIT-002; KIT-005 acceptance must include in-run fingerprints |
| PA-02 | medium | `journey$correction$harness_commit` records `3512f84d`, but the script and driver that produced the accepted record include the state-attach fix committed afterwards at `dfd6aacd` (uncommitted when the correction ran). | Provenance of the correction points at the wrong revision. | Successor note in the correction record naming `dfd6aacd` as the code that ran (the committed script is byte-identical to what ran). | none |
| PA-03 | medium | In the failure member the guard was injected and the sample step was marked `start` at the same uptime; the step failed 0.27 s later on disconnection, before any sample-filter work. The head check therefore proves stability across a stop *between* steps, not during a step's write. | The protocol's lifecycle requirement (completion and forced failure through the changed launch path) is met; the stronger reading, "head unchanged after a failed step", is weaker than CORE's mid-aggregate fault. | State the claim scope in the outcome; design KIT-005's operation-member failure as a mid-step fault (as CORE-004's `after_aggregate`). | none; KIT-005 design |
| PA-04 | low | Four launches of the complete member against a ledger allowing one correction; attempts 3 and 4 were user-authorised and each cause is recorded. The ledger row says "one recorded install"; three install runs happened (`sLJ4aC` failed under the kit root for missing exact-sum artifacts, `uvS55B` superseded by the native-path fix, `y8rDOi` bound); the outcome names only `y8rDOi`. | Ledger accounting incomplete. | Add the install sequence to the outcome's launch record. | none |
| PA-05 | low | `kitJourneySteps()` declares `import` leaves `raw_data_s4`; the store shows import leaves `initial` and design creates `raw_data_s4`. The journey record flags `state_as_expected = FALSE` for `import` (and trivially for `session_export`). | Step table labelling; populations are attributed correctly after the `dfd6aacd` fix. | Set the `import` step's state to `initial`. | none |
| PA-06 | low | On screen corpus v2 every design `replicates` value is 1, so the technical-replicate blocks coincide with the groups (3 blocks, 7 samples). The record shows `ruv_replicate_structure = technical_replicates` (under the `ruvIII_C_VaryingFn` args section; the `ruvIII_C_Varying` section holds a stale `ruv_number_k = 5` from the search), but the journey cannot distinguish the new structure from the old one. The full corpus (`replicates = ceiling(i/9)`) does. | The RUV decision is installed and recorded but not exercised by this proof. | None for KIT-001; KIT-005's control journey is the first exercise. Note the `…Fn` section naming for KIT-004's identity work. | none |
| PA-07 | info | The protein intensity threshold is not recorded in the state args (helper-internal); the ceiling removal is verified from the installed code and from the recorded peptide threshold (782261.4, non-integer). | None. | KIT-015 may record the protein threshold in args for exactness evidence. | none |
| PA-08 | low | The ceiling was removed at five proteomics sites, including `removePeptidesWithMissingValuesPercent` (`func_prot_qc_peptide_methods.R:213`) and `func_prot_s4_qc_filter_methods.R:629`, which the user's decision did not name explicitly (the two flagged filters); the extension was disclosed in the outcome. | Scope of a science change slightly wider than the words of the decision. | User confirmation that the five-site scope is intended (otherwise revert the two unnamed sites and reinstall). | decision (user) |
| PA-09 | info | Host profile `low-power` for both members; recorded, not a refusal (protocol §1). Time is reported only. | None. | — | none |
| PA-10 | low | Ticket text cites 51 harness assertions (now 54); the outcome says the failure member ran "with the corrected harness", which held the masking fix but not the attach fix (irrelevant to its decision). | Documentation precision. | Correct both sentences. | none |

## Evidence layers

- Mechanical: hashes, bindings, pins, schema validation: pass.
- Semantic: the receipts' fields match the protocol's member kind (journey member metrics present; selectors as §4); decisions computed by `kitLifecycleDecision` on the amended receipt.
- Empirical: both members ran under the managed slice with valid measurements; one correction applied post hoc (PA-01).
- Promotion: none claimed; ticket status is workflow metadata.

## States

| State | Value |
|---|---|
| implementation_verified | yes (tests, hashes, installed build inspection) |
| execution_complete | yes (both lifecycle members) |
| decision_supported | yes, with the claim limits PA-01 and PA-03 |
| promotion_authorized | not applicable |
| claim_scope | lifecycle proof on screen corpus v2 under the `kit-journey-v1` record; per-step snapshots in-run, fingerprints post hoc; the failure is a between-step stop; RUV structure not distinguishable on this corpus |

## Next action

No blocker for KIT-002. PA-08 needs the user's confirmation. PA-01 to PA-05 and PA-10 are bounded documentation and harness corrections for a separate remediation record if the user asks.
