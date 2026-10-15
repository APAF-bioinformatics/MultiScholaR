# KIT-017: DIA stage closure

Status 2026-10-15: **counted pair run; stage `invalid`, awaiting the user's decision on the candidate.** The control member is complete and verified; the candidate failed at the DA step on a defect found only at full size (remediation record CP-02); the ledger's two members are spent. Ticket revision 1.1: control `kit_018` (`core/browser-smoke.yXPrcl`, `6a7304e7`), candidate `kit_017` (`core/browser-smoke.sOtQzu`, `c9635938` = the `kit_024` source plus the screen-run repairs), rule `resident_end`.

## What the closure measures and why

The pair is a protocol v1.4 ticket journey pair: each build one journey member as users run the GUI (no launch record) on the 300-sample corpus, in one session, the candidate's states verified against the KIT-005 oracles. The bound rule is `resident_end`: the one memory mechanism between the two builds is KIT-019's DA publication through the worker, whose claim is end-of-journey resident memory (candidate ≤ 0.95 × control), with the whole-slice anonymous peak not worse than 1.02 × and the per-step guard from the first saved state. The stage's whole-slice reduction is reported against the production-default baselines (KIT-005 `default-1` on `kit_001`, KIT-011's `kit_006` control; both guard readings, so lower bounds), as protocol section 4 says; no mechanism between `kit_018` and `kit_024` claims a whole-slice reduction, so the pair does not judge one ([audit](KIT-017-PREEXEC-AUDIT-20261015.md) PE-02).

## Before the counted launch

- [Pre-execution audit](KIT-017-PREEXEC-AUDIT-20261015.md) (ten findings) and its [remediation record](KIT-017-PREEXEC-REMEDIATION-20261015.md), which also carries the three screen lifecycle attempts. PE-01: the bindings committed with KIT-034 hashed the journey driver at its pre-KIT-034 content (regenerated). PE-03: neither build had been observed through a managed journey past design, so the pair ran on the screen corpus first (the harness pair branch takes a bound corpus key, the KIT-018 pattern).
- Screen attempt 1 failed both members and exposed two plain code paths only a live browser reaches: the server-side paged tables (KIT-024) resolved the state owner inside Shiny's request handler, which carries no reactive domain (SR-01), and the application monitor died reading the worker registration file in the window after the controller unlinked it, taking the application with it (SR-02). Attempt 2, on the repaired candidate, hung at the same point: the served page was a named data frame, so its rows serialised as JSON objects where the widget reads arrays, and DataTables' native error dialog blocked the renderer (SR-04). The managed driver now records and dismisses native dialogs and fails the next command by name (SR-05). Attempt 3 completed and verified on both builds (control 561 s, candidate 563 s; [decision](screen-lifecycle-attempt3-decision.rds)).
- Installs: `kit_017` first as `core/browser-smoke.UpIfdA` (`9293655c`, superseded), then `core/browser-smoke.sOtQzu` (`c9635938`); both recorded, neither counted. The control build keeps SR-02's race; a recurrence there is a demonstrated instrumentation failure (one correction rerun of that member).
- Predeclared before the outcome was visible (remediation record, attempt 3): the DA publication worker is a second R process inside the measured slice (about 0.5 GiB on the screen corpus, where the DA step read 1.28 × and the whole slice 1.14 × while resident end read 0.825), so the DA step, which is the journey's whole-slice peak on 300 samples, and possibly the whole slice may exceed the frozen 1.02 bound. The rule is not changed; a rejection on that ground closes the stage as measured-not-accepted with the shortfall stated.

## Result

Stage `kit-journey-017` (tracked summary [kit-journey-017-stage.json](kit-journey-017-stage.json); 300-sample corpus, no launch record, 14,400 s ceiling; launched 12:04, ended 15:09):

| Member | Build | Outcome | Elapsed | Slice peak | Anonymous peak | Resident end | Max gap | Verified |
|---|---|---|---|---|---|---|---|---|
| control-1 | `kit_018` | **complete, all 18 steps** | 4,777 s | 7.27 GiB | 5.74 GiB | 2.87 GiB | 150 ms | every state equal to the KIT-005 oracles; the rendered volcano page reported as an additional product (correction 1, verifier only) |
| candidate-1 | `kit_017` (`sOtQzu`) | **failed at the DA step** (`da_publication`: unsupported completion context) | 4,830 s to the failure | 7.33 GiB | 5.69 GiB | 2.64 GiB at failure | 160 ms | not verified (failed journey) |

Decision after correction 1: `invalid` (`candidate-1 verifier_accepted`). No memory claim is made from this pair. Per step through correlation the two builds read within 1.03 of each other (control / candidate anonymous peaks in GiB: q-value 4.94 / 4.79, rollup 5.74 / 5.69, IQ 4.50 / 4.33, normalisation 3.00 / 3.00); the DA step, where the predeclared worker cost would show, was not reached by the candidate (control 5.29 GiB). These are reported, not judged.

Against the production-default baselines the control's whole-slice anonymous peak is 0.37 × KIT-011's `kit_006` guard reading (15.40 GiB, a lower bound) and the journey completes where the in-memory default could not; this is the stage's standing whole-slice evidence, unchanged by this pair.

The candidate's defect (CP-02) is repaired in the working tree (`workflowManagedCompletionContext`); a relaunch needs an install and a candidate member the ledger does not hold.

## Regression (standing rule: the whole suite once per ticket, classified at the previous commit)

| Commit | Suite (343 files) | `R CMD check` (tarball from a clean worktree under `/home`, proof library on `R_LIBS`, suggests not forced) |
|---|---|---|
| `c439ce88` (binding, harness only) | 29,627 pass, 0 fail, 1 error, 139 skip | 3 warnings, 3 notes ([logs](check-logs/)) |
| `9293655c` (round-one repairs) | 29,631 pass, 0 fail, 1 error, 139 skip | 3 warnings, 3 notes |
| `c9635938` (round-two repairs, the candidate's source) | 29,632 pass, 0 fail, 1 error, 140 skip | 3 warnings, 3 notes |

The one error in every run is `test-workflow-perf101-diagnostic.R` "request source snapshots are self-contained and detect missing or changed bytes", which passes alone (20 pass) and in KIT-034's run: the regression script builds the tarball in the same worktree while that test snapshots source bytes. An artefact of the runner, not of the code; recorded, the runner left as is. KIT-002's check record had 4 warnings and 3 notes; the code/documentation-mismatch warning is gone and no warning or note is new.

## Kit completeness checklist

[kit-checklist.md](kit-checklist.md): every kit function exists with no DIA logic in its body (a naming debt recorded: the storage primitives keep the `DIA` prefix), all 16 data steps are registered declarations, `partitionByPlex` is exempt until KIT-020, the three kernel adapters are named, the quarantine state and the overtaken-tests record are stated.

## Hand-off

Stage 2 and stage 3 inputs are in the checklist's section G; the technical-replicate declaration (KIT-035) and the analysis-record series (KIT-025 to KIT-029, KIT-033) follow this closure.
