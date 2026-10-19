# Proposed amendment: KIT measurement protocol v1.5 (2026-10-19, pending the user's decision)

`protocol.md` (v1.4) is frozen; its SHA-256 is bound in `bindings.json`, and any change is a new version needing a user decision. This file is the proposal; nothing in it is in force until the user records the decision in the plan's decision record, after which the text below is appended to `protocol.md` as version 1.5 and the bindings are regenerated. KIT-021 and KIT-053 run no counted member before that.

Why: the series audit (`KIT-SERIES-AUDIT-20261019-stage2.md`, A-25, A-26, C-24, S-10) found that the LFQ/TMT and `DIA_LIMPA` tickets name member kinds, oracle sources and selectors that v1.4 defines only for the DIA journey on the 300-sample corpus against the KIT-005 oracles.

## Text proposed for §4 and §5

1. **Oracle source for a route or corpus without a KIT-005 oracle (§5).** Where a member's route (LFQ, TMT, `DIA_LIMPA`) or corpus (`lfq_300`, `tmt_fixture`) has no KIT-005 oracle, the per-step oracle is an in-process run of the **bound build** on the **bound corpus**, executed read-only and recorded in the ticket's evidence before the member is launched; the member's states are verified against it. Such an oracle proves that the owner path reproduces the in-process path of the same kernels; it is not a cross-build claim. The DIA journey keeps the KIT-005 oracles.
2. **Single journey member, whole journey (§4).** The single journey member (KIT-018) may be bound with no stop step, in which case it runs the whole journey of its route; metrics as the journey member (whole-slice anonymous peak, per-step peaks, `resident_end_bytes`), completion and exactness judged, no ratio and no control. Used by KIT-021 (LFQ, TMT), KIT-053 (`DIA_LIMPA`) and KIT-048 (DIA after the retirement).
3. **Selectors (§4).** A single journey member is `kit:journey-0NN<suffix>:candidate:1`, where the optional suffix distinguishes several single journeys under one ticket (`kit:journey-021:candidate:1`, `kit:journey-021t:candidate:1`); each ticket key binds one `kind` and one `ceiling_seconds`. The operation-member form `kit:0NN:(control|candidate):[12]` is never used for a journey.
4. **Recorded runs without a managed slice (§1, §4).** A run launched without a managed slice (for example with the GNU tools or `systemctl` removed from `PATH`, KIT-044) is not a member of any kind, carries no measurement claim, and is recorded in the ticket's evidence, never counted.
5. **Route step maps (§5).** Each route has its own step map file hashed in the bindings (`kit-001/step-map.json` for DIA, `step-map-protein.json`, `step-map-limpa.json`); a step is declared only when its row carries no `flag`.

No decision rule, threshold, workload, ceiling or limit of v1.4 changes.
