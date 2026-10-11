# KIT-011 post-hoc audit (2026-10-12)

Mode: `post_hoc`, audit only (no repairs; no member rerun). Scope: the KIT-011 execution as closed at `5ade1884`: the install `kit_011` (`core/browser-smoke.AFaqC1`, source `b3d44bdb`), the screen lifecycle member on `kit_011` and its cross-build check, stage `kit-journey-011` (control `kit_006`, candidate `kit_011`, no launch record), correction 1 (verifier only), the `accepted` decision, the promotion of `kit_011` to control, and the records (`outcome.md`, `kit-journey-011-stage.json`, the ledger, the ticket) against protocol v1.3, the user's rules of 2026-10-11 (per-step guard from q-value; whole-slice against the control's observed peak) and the plan ledger.

Mechanical inventory: `audit_inventory.py` exit 0, 0 findings; harness hashes match the tree; working tree clean.

## What was verified

| Check | Result |
|---|---|
| Receipts against the raw traces (samples, slice peak, anonymous peak, swap, max gap, gaps > 250 ms) | control-1 12,567 / 15.538 / 15.398 GiB / 0 / 170 ms / 0; candidate-1 49,867 / 7.979 / 5.508 / 0 / 170 ms / 0; screen member 8,335 / 2.527 / 2.099 / 0 / 80 ms / 0: all as recorded |
| Control's guard | `guard-stop.txt` = `memory_guard`; slice peak 15.538 GiB ≥ the 15.5 GiB early stop; 8 completed steps with finite peaks |
| Decision arithmetic | whole-slice ratio 5.508 / 15.398 = 0.3577 as recorded; guarded steps q-value to replicate 0.23 to 0.48; import 0.98 and design 0.81 reported, not judged, as decided before launch (ledger row of 2026-10-11 precedes the launch file at 15:57 AEDT) |
| Candidate verification | 18 populations with in-run fingerprints; recomputed from the immutable store; class, rows, components, support and summary identical in-run vs recomputed on every state; 0 oracle mismatches against the 17 KIT-005 oracles; the four operation identity products reported; `canonical_form_differences` = normalisation_ruv only |
| Correction 1 | verifier only; `verification.harness-failed.rds` kept (`normalisation_ruv/products/normalization_science`); the launch decision (`invalid`) kept in `decisions_before_correction`; the recorded harness sha256 (`1cb0220a…`) equals the committed `workflow_kit.R`; no rerun directory |
| Authority | the install, the screen member and the pair were authorised ("let's go with recommendations"); one launch of the pair; the promotion follows protocol section 2; `tickets$011$control` keeps `kit_006` as the control it ran against |
| Host | 21.3 GiB available at launch, no managed unit, low-power profile |

## Findings

| ID | Severity | Layer | Observed | Impact | Smallest repair | Blocked |
|---|---|---|---|---|---|---|
| PY-01 | low | verifier design | Correction 1 bundles two verifier changes: the canonical form drops the KIT-002 `declaration` field, and products are now judged on the recomputed state-record entries overlaid on the in-run per-step assets (the first change exposed that the recomputed state fingerprints carry only the record's own products, not the per-step assets the driver attaches from the head). The integrity rule between in-run and recomputed fingerprints is therefore strict on class, rows, components, support and summary and reports (rather than refuses) product differences. | The per-step assets are judged from the in-run record (whose hashes equal the oracle's) rather than recomputed from the store. | None required for this decision; a later ticket may recompute per-step assets from the head's product history. | none |
| PY-02 | low | claim scope | The record-less verifier branch fingerprints the candidate's 528 publication files but compares them with nothing (the KIT-005 file comparison exists only between members of the same build). The outcome's claim limits do not say so. | The exactness claim covers states, support, summaries and content-addressed products, not publication files. | One sentence in the outcome's claim limits. | none |
| PY-03 | info | provenance | The correction ran on the uncommitted verifier and the recorded hash equals the file committed in `5ade1884`; the KIT-006 audit's commit-label defect did not recur (the record carries the hash, not a commit). | None. | None. | none |
| PY-04 | info | evidence shape | Acceptance item 2 ("the legacy v1 archive still opens") rests on `test-workflow-native-session.R` and the legacy ingest tests, not on a run; the candidate's own session export step completed on the default path. | Stated in the ticket. | None. | none |
| PY-05 | info | bindings | The promotion sets `journey.pair.control` to `kit_011`; the KIT-005 baseline block still names `kit_001` as the control it ran under. | None. | None. | none |

## Closeout states

- `implementation_verified`: true (pre-execution audit; screen member on the installed build; candidate journey on the installed build).
- `execution_complete`: true (install, screen member, two journey members, one verifier correction).
- `decision_supported`: true (`accepted`; every eligibility rule holds; the reduction is judged against a lower bound of the control's peak, which can only understate it; the per-step guard holds on every judged step).
- `promotion_authorized`: true by protocol section 2 (`kit_011` is the latest accepted build; bound as the control, recorded in the ledger and the bindings).
- `claim_scope`: on Linux with DIA-NN and the IQ route, the GUI as users run it now completes the 300-sample journey at 5.5 GiB anonymous and 3.0 GiB resident where the shipping build stopped at the 15.5 GiB guard after 8 steps; every state equals the KIT-005 oracles; publication files and the operation identity products are not judged; Limpa, other formats and non-Linux are unchanged and unmeasured; times are reported only.

Next action: none blocking. PY-02 is a sentence for the next closeout.
