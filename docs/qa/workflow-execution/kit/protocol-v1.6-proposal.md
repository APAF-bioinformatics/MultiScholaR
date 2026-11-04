# Proposed amendment: KIT measurement protocol v1.6 (2026-11-01) — SUPERSEDED: in force as `protocol.md` section 8 since 2026-11-04, with the user's change to the per-step guard at worker steps

`protocol.md` (v1.5) is frozen; its SHA-256 is bound in `bindings.json`, and any change is a new version needing a user decision. This file is the proposal. Nothing in it is in force until the user approves the text and KIT-081 appends it to `protocol.md` as section 8, before any stage-3b counted member (KIT-080). The proposal reuses the protocol's existing thresholds (0.95 reduction, 1.02 not-worse); it adds no new threshold.

Why: the stage-3b audit (`STAGE3B-AUDIT-20261101.md`) found the following, and the user's "fix all please" (2026-11-01) took the audit's recommended repairs:
- **S3B-19:** the managed slice holds the controller, the worker and the monitor, so the protocol's whole-slice metrics cannot show what KIT-078 is for (the application process holds no matrix during a step).
- **S3B-26:** the operation-member selectors take no suffix and allow two pairs per key, so KIT-080's per-step operation pairs cannot be bound.
- **S3B-27:** the ticket journey pair is defined on the 300-sample DIA corpus with the KIT-005 oracles and the latest accepted control, which a lipid/metab pair on a real corpus is not.

## Text proposed for §4, §5 and §2

1. **Stage-3b control (§2).** For the lipidomics and metabolomics journey pairs of stage 3b, the control is `kit_032`, the stage-3 closure build: the last build before the stage-3b mechanisms, with the on-disk owner and declarations of stage 3. The "latest accepted KIT build" rule is unchanged for proteomics.
2. **Pair corpus and oracles (§4, §5).** A ticket journey pair may run on a real corpus registered in the bindings by sha256 and aggregate dimensions only (the KIT-054 data rule). Each member is verified against an in-process oracle of its own build on that corpus (§7.1). The decision judges the candidate's exactness. The control's exactness is reported.
3. **Rule `app_process` (§4).** Measured from the controller service's cgroup inside the managed slice, sampled with the slice by the KIT-081 sampler at the same 50 ms cadence and gap rules. The candidate's controller anonymous peak must be ≤ 0.95 × the control's controller anonymous peak, at each step whose compute moves to the worker in the candidate (normalisation, correlation, DA; listed in the binding). The control has no worker, so its controller peak is its application's peak. A pair may bind `app_process` together with `whole_slice_anon`. Each rule's decision is recorded verbatim, neither is retried, and the ticket states beforehand which claim each rule decides.
4. **Operation-member selectors (§4, §7.3).** An operation stage may carry a suffix, `kit-0NN[a-z]`, with selectors `kit:0NN[a-z]:(control|candidate):[12]`, so one ticket may bind several operation keys of at most two pairs each. Journey selectors are unchanged.
5. **Operation ceiling (§4).** A suffixed operation key binds its `ceiling_seconds` (at most 1,800 s). Unsuffixed operation keys keep 600 s.
6. **Operation members on a list assay (§4).** For a class whose assay is a named list of tables, the operation member's input component is the list. It is identified by the digest of its sorted element entity ids, which is the scientific identity's rule since KIT-031, and fingerprinted per element.

## What does not change

Every existing rule, threshold, guard, gap tolerance, ceiling and selector for the proteomics members; the run ledger's authority; the correction-rerun rule; the KIT-005 oracles for the DIA journey.
