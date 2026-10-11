# KIT-006 post-hoc audit (2026-10-11)

Mode: `post_hoc`, audit only (no repairs; no member rerun). Scope: the KIT-006 execution as closed at `53ec213d`: the install `kit_006` (`core/browser-smoke.qBR5tm`), the installed preflight (`core/browser-smoke.095yCi`), stage `kit-006` attempt 1 (`8d8efccc`), rerun 1 of the two sample-filter members, correction 1 (verifier only), the `accepted` decision, the promotion of `kit_006` to control, and the records (`outcome.md`, `kit-006-stage.json`, the ledger, the ticket) against protocol v1.3 and the plan ledger. The pre-execution audit (`KIT-006-PREEXEC-AUDIT-20261011.md`) covered the implementation; its findings were fixed before launch.

Mechanical inventory: `audit_inventory.py` exit 0, 0 findings; harness hashes match the tree; working tree clean; the selector returns KIT-007.

## What was verified

| Check | Result |
|---|---|
| Member receipts against the raw `memory-trace.tsv` (samples, slice peak, anonymous peak, swap, max gap) | control-1 1,314 / 6.607 / 5.254 GiB / 0 / 80 ms; candidate-1 1,284 / 6.547 / 5.232 / 0 / 70 ms; control-2 517 / 1.584 / 1.332 / 0 / 60 ms; candidate-2 565 / 1.578 / 1.324 / 0 / 70 ms: all as recorded; the two refused attempt 1 members (42 and 46 samples, 0.35 GiB) preserved as `*.attempt1` |
| Pair ratios | P14 0.991 (slice) / 0.996 (anonymous); sample filter 0.996 / 0.994; both ≤ 1.02 |
| Exactness | control-1 and candidate-1: 1,200 P14 pages each against the P14 oracle, `exact_sum_sha256` equal to each build's own bound identity (`49798e98…`, `5a7c1153…`); control-2 and candidate-2: 615 pages, 6,150,000 rows, 225 samples, design, summary and support identical to CORE's native record; output mode `physical` on both |
| Batch rules | P14 3 batches of 100 on both builds (from the verification records); sample filter 3 on the control (fixed policy) and 3 on the candidate (adaptive, 4.07 GiB estimate under a 14.3 GiB cap) |
| Eligibility | every member complete, verified, zero swap, gap ≤ 250 ms, elapsed ≤ 600 s, peak below the early stop, slice stopped |
| Rerun bookkeeping | `stage.rds` `reruns[[1]]` records the two cases, the rerun bindings hash, the rerun harness hashes (`source.rerun1`, `workflow_kit.R` `6e15b6f6…` = the file at `8d8efccc`), host at start and build; the attempt 1 and rerun 1 decisions kept in `decisions_before_rerun` (both `invalid`); case bindings identical between attempt and rerun (asserted by the coordinator) |
| Correction 1 | verifier only; `verification.harness-failed.rds` kept beside the accepted `verification.rds`; `verify-correction-1` run directory present; one correction, within the ledger's allowance |
| Authority | the rerun was authorised by the user ("go" on the stated amendment) and recorded in the ledger before launch; the promotion follows protocol section 2 (latest accepted build) |
| Install | source commit `62fb9aec` = the tree at launch; empty runtime patch; preflight 174 checks, 0 failed, five exact P14 cases |

## Findings

| ID | Severity | Layer | Observed | Impact | Smallest repair | Blocked |
|---|---|---|---|---|---|---|
| PX-01 | low | provenance | The correction 1 record in `stage.rds` names `harness_commit 8d8efccc`, but its `harness_sha256` (`a4afa9f7…`) is the `workflow_kit.R` with the verifier fallback, which was uncommitted when the correction ran and was committed in `53ec213d`; the file at `8d8efccc` hashes `6e15b6f6…` (the rerun snapshot). The hash is the truthful identity; the commit label is wrong. | A reader following the commit label finds the pre-fix verifier. | One sentence in the outcome (the receipt's hash stands); future correction scripts should record the commit after committing or label it "uncommitted". | none |
| PX-02 | low | generator | `bindings.json` has `control = kit_006` but `journey.pair.control = kit_001`: the generator builds the journey pair before the block that promotes `kit_006`. For KIT-005's stage the historical value is right; a later journey pair (KIT-011) generated the same way would bind the wrong control unless the ticket binds its own pair. | A stale default for the next journey pair. | Build the pair after the promotion, or derive the pair control from `b$control` at the end; regenerate. | none |
| PX-03 | info | decision record | `decision$pairs[["1"]]$batches` is NULL for the P14 pair because the attempt 1 receipts predate the `science_execution` field; the batch counts are in the verification records and the summary reads them from there. Stated in the outcome's claim limits. | None. | None. | none |
| PX-04 | info | budget | Counted members used: four, as declared. Two additional slice launches (attempt 1's sample-filter members) were refused in 3 to 4 seconds before any work and are recorded in the ledger; the harness defect behind them was the kit's own case-parameter typing, caught on the first operation with a vector parameter. | None; recorded. | None. | none |
| PX-05 | info | host hygiene | Two user units remain listed after the stage (inactive, as in the KIT-005 audit). | None. | `systemctl --user reset-failed`. | none |

## Closeout states

- `implementation_verified`: true (pre-execution audit plus the member verifications on both builds).
- `execution_complete`: true (install, preflight, four members, one authorised rerun of two members, one verifier correction).
- `decision_supported`: true (`accepted`: both pairs exact and not worse; eligibility rules hold on every member; the rerun and correction are recorded and within authority).
- `promotion_authorized`: true by protocol section 2 (`kit_006` is the latest accepted build; bound as the control, recorded in the ledger and the bindings).
- `claim_scope`: the complete-sample partitioner under the P14 rollup and the physical sample filter is exact on `full_p14` and `full_primary` and not worse in memory than the previous control; no reduction is claimed; operation members run far below the cap (under 7 GiB and 2 GiB), so this is a mechanism and exactness proof.

Next action: none blocking. PX-01 and PX-02 are a sentence and a generator reorder, suitable for the next ticket's closeout or a short fix pass.
