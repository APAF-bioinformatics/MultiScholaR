# KIT-011 pre-execution audit (2026-10-11)

Mode: `audit_only` (no repairs; nothing launched). Scope: the KIT-011 implementation commit `65f0ef2e` (R: the native default, the removed capacity trial and rollup arm, revert by reference, the saved-reactive alias; harness: the ticket journey stage, its decision and the oracle verification; the generator; the test updates; `outcome.md` with the declared pair) against protocol v1.3, the plan ledger (KIT-011: two journey members) and the user's KIT-005 amendment.

Mechanical inventory: `audit_inventory.py` exit 0, 0 findings; harness hashes match the tree; working tree clean.

## What was verified

| Check | Result |
|---|---|
| R diff (128 changed lines) | every hunk belongs to one of: `isDIANativeDefaultEnabled` constant, removal of `initializeDIACapacitySession` / `computeDIACapacityRollup` / `prepareDIACapacityRollupInput` / `getDIARollupEnvelope` and their call sites, the `sql_rollup` branch of `runPeptideRollupApplyStep`, the launch-record option line, the protein-S4 revert, the `savedReactiveVal` alias, the attach scope label |
| Remaining callers of the default flag | `initializeDIANativeImport` (platform and format limits kept), `ingestDIANativeLegacySession` (adopts legacy DIA sessions on Linux by default), `mod_prot_qc_protein_rollup.R` (the route-mismatch refusal is now default behaviour) |
| `reference_states` | set explicitly TRUE by the import and the legacy adoption; no derivation from the session format remains |
| Reclamation | `reclaim_enabled` defaults FALSE at adapter creation, as in KIT-005's native members (B26 stays a KIT-012 decision) |
| Harness | `kit-journey-011` expands to `control-1` (kit_006, no record) and `candidate-1` (kit_011, no record); the coordinator writes no `gui-member.rds` for either and checks `is.null(config$journey)`; the driver snapshots a store when the build creates one; the verifier recomputes the candidate's fingerprints and compares them to the 17 bound oracles with the same comparator KIT-005's control-2 used; `kitJourneyTicketDecision` tested on synthetic receipts (guarding control accepted, missing reduction rejected, worse step rejected, guarding candidate invalid, complete control with ratios) |
| Tests | 19 affected files green except the known installed-runtime failure in `test-workflow-native-session.R`, identical on the committed baseline |

## Findings

| ID | Severity | Layer | Observed | Impact | Smallest repair | Blocked |
|---|---|---|---|---|---|---|
| PX-01 | medium | decision design | The per-step not-worse guard (≤ 1.02 × control on every step the control completed) is applied to import and design, where the two paths run the same code and KIT-005 measured the native path above the default on import (2.9 vs 2.7 GiB, ratio 1.07; design 3.3 vs 3.5). A single noisy import step would reject a candidate whose data steps are 2 to 4 times lower. | A `rejected` decision for a reason unrelated to the mechanism. | User decision before launch: keep the rule as predeclared, or apply the per-step guard from the first saved state onward (q-value) and report import and design. One line in `kitJourneyTicketDecision` plus a ledger row. | the decision rule is the user's |
| PX-02 | medium | exactness assumption | The candidate is verified against the KIT-005 oracles produced on `kit_001`. Since then KIT-006 (partitioner, proven exact on `full_p14` pages and on the native filter record) and KIT-007 (registry) changed the build, and KIT-011 removed the rollup-input preparation; no whole-journey cross-build fingerprint check has been run on any later build. A benign cross-build difference in one state's summary would make the 65-minute candidate `invalid` at verification, which a verifier correction cannot repair. | Risk of an uninformative long run. | Before the pair: one screen lifecycle member on `kit_011` (`kit-lifecycle-complete`, about 8 minutes, the KIT-001 screen directories archived) and `kitReproduceMismatches` of its 18 populations against the KIT-001 screen member's (`kit-lifecycle-complete/member/journey-result.rds`, produced on `kit_001`, 18 fingerprints present). This also satisfies the protocol's "small lifecycle test before a long run" for a changed launch path, which the production default now is. A guarded run: it needs the user's go. | launch of the pair, by recommendation |
| PX-03 | low | ticket metadata | `files_to_modify` lists files not changed (`mod_prot_import_*`, `utils_workflow_dia_qc_rollup.R`, `workflow_kit_driver.R`) and omits changed ones (`mod_prot_qc_peptide_protein.R`, `mod_prot_qc_protein_s4.R`, `workflow_kit.R`, the coordinator and runner, the generator, the 01c/01d import tests, `test-workflow-saved-qc-rollup.R`, `test-workflow-kit-harness.R`). | Allowlist discipline. | Ticket revision at closeout. | none |
| PX-04 | info | behaviour change | Two behaviours gated on the flag become default: a legacy DIA session opened on Linux is adopted into the store; a DIA-NN project whose import recorded the Limpa route refuses an IQ rollup until reset (and vice versa). Both are the behaviours the gated path already had. | Users of the GUI see them without an option. | State both in the outcome's claim limits. | none |
| PX-05 | info | compatibility | `attachDIAStateSession`'s default scope label is `dia-native-session`; the PERF-era driver and proof runner still name `full-journey-capacity-v1` for their own receipts. The label is written only when a receipt directory is given. | None. | None. | none |
| PX-06 | info | storage | Older stores carry `rollup_input` in their head snapshot; nothing reads it now, and `publishDIAStateHead` passes unknown snapshot fields through. | None observed; a reopen of a KIT-005 store was not re-tested after the removal. | Covered by PX-02's screen member only partly (fresh store); a reopen of the KIT-005 control-1 store through `kitReadHead` is a ten-second check. | none |

## Closeout states

- `implementation_verified`: true for the code paths exercised by tests; the production default end to end on an installed build is unverified until the PX-02 screen member.
- `execution_complete`: false (install and pair declared).
- `decision_supported`: not yet; the rule carries the PX-01 risk.
- `promotion_authorized`: false.
- `claim_scope` so far: the on-disk path is the DIA-NN/IQ default on Linux in the source tree; every DIA holder the import writes is a saved reference; no memory claim.

Recommended order: decide PX-01; run the PX-02 screen member and its comparison; then the install and the pair.
